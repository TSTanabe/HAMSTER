#!/usr/bin/python

import statistics
import os
import sqlite3
import csv
import multiprocessing
import pickle
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from pprint import pprint




def group_gene_cluster_statistic(options):
    """
    Berechnet und speichert (oder lädt) die Statistiken zu Gen-Clustern.
    
    Args:
        options: Konfigurationsoptionen, die Datenbank- und Verzeichnisinformationen enthalten.
    
    Returns:
        tuple: 
            - filtered_stats_dict (dict)
            - grouped_keywords (dict)
            - clustered_excluded_keywords (dict)
    """
    cache_dir = os.path.join(options.result_files_directory, "pkl_cache") 
    os.makedirs(cache_dir, exist_ok=True)
    # Load cache if available
    filtered_stats_dict = load_cache(cache_dir, "filtered_stats.pkl")
    domain_score_limits = load_cache(cache_dir, "domain_score_limits.pkl")
    grouped_keywords = load_cache(cache_dir, "grouped_keywords.pkl")
    distant_keywords = load_cache(cache_dir, "distant_keywords.pkl")
    clustered_excluded_keywords = load_cache(cache_dir, "clustered_excluded_keywords.pkl")
    query_score_dict = load_cache(cache_dir, "query_score_dict.pkl")  # Try loading from cache first
    
    
    if domain_score_limits and grouped_keywords and clustered_excluded_keywords:
        print("Loaded existing CSB grouping from cache.")
        return domain_score_limits, filtered_stats_dict, grouped_keywords, clustered_excluded_keywords

    # Schritt 1: Berechnung der Keyword-Statistiken
    if not filtered_stats_dict or not query_score_dict:
        print("Computing keyword statistics")
        stats_dict = get_keyword_statistics_parallel(options.database_directory, options.cores)

        # Schritt 2: Extrahiere das höchste Bitscore pro Domain in QUERY
        print("Extracting highest bitscores per domain")
        query_score_dict = get_highest_bitscores_for_genome(options.database_directory, "QUERY")
        save_cache(cache_dir, "query_score_dict.pkl", query_score_dict)

        # Schritt 3: Entferne CSBs, deren Domains unter 80% der Query-Referenz liegen
        print("Filtering out low-quality CSBs")
        filtered_stats_dict = filter_out_low_quality_csb(stats_dict, query_score_dict, 0.2, min(10,options.min_seqs))
        
        save_cache(cache_dir, "filtered_stats.pkl", filtered_stats_dict)
    
    # Schritt 4: Speichere gefilterte Statistiken als TSV
    save_stats_to_tsv(filtered_stats_dict, options.Csb_directory)

    # Schritt 5: Gruppiere Keywords pro Domain mit einer maximalen Abweichung von 30%
    if not grouped_keywords or not distant_keywords:
        print("Grouping keywords by domain")

        grouped_keywords, distant_keywords = group_keywords_by_domain_extended(filtered_stats_dict, query_score_dict, 0.30)
        save_cache(cache_dir, "grouped_keywords.pkl", grouped_keywords)
        save_cache(cache_dir, "distant_keywords.pkl", distant_keywords)

    # Schritt 6: Clustere ausgeschlossene Keywords basierend auf statistischer Ähnlichkeit
    if not clustered_excluded_keywords:
        print("Clustering excluded keywords by similarity")

        clustered_excluded_keywords = group_excluded_keywords_by_similarity(filtered_stats_dict, distant_keywords)
        save_cache(cache_dir, "clustered_excluded_keywords.pkl", clustered_excluded_keywords)

    # Schritt 7: Oberes und unteres Score-Limit der gruppierten CSBs ohne Outliers
    if not domain_score_limits:
        print("Computing score limits for grouped CSBs")
 
        domain_score_limits = compute_score_limits(filtered_stats_dict, grouped_keywords)
        save_cache(cache_dir, "domain_score_limits.pkl", domain_score_limits)

    return domain_score_limits, filtered_stats_dict, grouped_keywords, clustered_excluded_keywords








#################################################################################
# Fetch for each csb for each domain all the possible scores. Calculate the min,max,mean,median and std_dev and print them
    

def wrapped_fetch_keyword_scores(database, chunk, progress_counter, lock):
    """Fetches keyword scores and updates progress."""
    result = fetch_keyword_scores(database, chunk)
    with lock:
        progress_counter.value += 1
        print(f"Progress: {progress_counter.value} genecluster statistics completed", end="\r")
    return result

def get_keyword_statistics_parallel(database, num_workers=4):
    """
    Parallelized version of keyword statistics calculation using multiprocessing.
    Shows a progress update when workers complete their tasks.

    Args:
        database (str): Path to SQLite database.
        num_workers (int): Number of parallel workers.

    Returns:
        dict: { keyword: { domain: { 'n': int, 'min': float, 'max': float, 
                                     'mean': float, 'median': float, 'std_dev': float } } }
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("SELECT DISTINCT keyword FROM Keywords;")
        keywords = [row[0] for row in cur.fetchall()]

    # Split keywords into chunks for parallel processing
    chunk_size = max(1, len(keywords) // num_workers)  # Avoid division by zero
    keyword_chunks = [keywords[i:i + chunk_size] for i in range(0, len(keywords), chunk_size)]

    # Progress tracking
    manager = multiprocessing.Manager()
    progress_counter = manager.Value("i", 0)  # Shared counter for completed workers
    lock = manager.Lock()  # Lock to prevent race conditions

    # Use multiprocessing pool with progress tracking
    with multiprocessing.Pool(num_workers) as pool:
        raw_results = pool.starmap(wrapped_fetch_keyword_scores, 
                                   [(database, chunk, progress_counter, lock) for chunk in keyword_chunks])

    # Merge raw scores
    merged_raw_scores = {}
    for raw_result in raw_results:
        for keyword, domain_dict in raw_result.items():
            if keyword not in merged_raw_scores:
                merged_raw_scores[keyword] = {}
            for domain, scores in domain_dict.items():
                if domain not in merged_raw_scores[keyword]:
                    merged_raw_scores[keyword][domain] = []
                merged_raw_scores[keyword][domain].extend(scores)

    # Compute statistics
    return compute_statistics_for_keywords(merged_raw_scores)


def fetch_keyword_scores(database, keyword_list):
    """
    Fetches raw domain scores for a given list of keywords.

    Args:
        database (str): Path to SQLite database.
        keyword_list (list): List of keywords to process.

    Returns:
        dict: { keyword: { domain: [scores] } }
    """
    raw_scores = {}

    with sqlite3.connect(database) as con:
        cur = con.cursor()
        
        query = """
        SELECT k.keyword, d.domain, d.score
        FROM Keywords k
        JOIN Domains d ON k.clusterID = (SELECT clusterID FROM Proteins WHERE proteinID = d.proteinID)
        WHERE k.keyword IN ({})
        AND d.domain IS NOT NULL 
        AND d.score IS NOT NULL;
        """.format(",".join(["?"] * len(keyword_list)))
        
        cur.execute(query, keyword_list)
        
        for keyword, domain, score in cur.fetchall():
            if keyword not in raw_scores:
                raw_scores[keyword] = {}
            if domain not in raw_scores[keyword]:
                raw_scores[keyword][domain] = []
            raw_scores[keyword][domain].append(score)

    return raw_scores


def compute_statistics_for_keywords(raw_scores):
    """
    Computes statistical values (count, min, max, mean, median, std_dev) for the given raw scores.

    Args:
        raw_scores (dict): { keyword: { domain: [bitscores] } }

    Returns:
        dict: { keyword: { domain: { 'n': int, 'min': float, 'max': float, 
                                     'mean': float, 'median': float, 'std_dev': float } } }
    """
    stats_dict = {}

    for keyword, domain_dict in raw_scores.items():
        stats_dict[keyword] = {}
        for domain, scores in domain_dict.items():
            if scores:
                stats_dict[keyword][domain] = {
                    "n": len(scores),
                    "min": min(scores),
                    "max": max(scores),
                    "mean": round(statistics.mean(scores), 2),
                    "median": round(statistics.median(scores), 2),
                    "std_dev": round(statistics.stdev(scores), 2) if len(scores) > 1 else 0.0,
                }

    return stats_dict


#########################################################################################



def save_stats_to_tsv(stats_dict, output_directory, filename="bitscore_statistics.tsv"):
    """
    Saves the computed bitscore statistics to a TSV file.

    Args:
        stats_dict (dict): { keyword: { domain: { 'n': int, 'min': float, 'max': float, 
                                                  'mean': float, 'median': float, 'std_dev': float } } }
        output_directory (str): Path to the directory where the file should be saved.
        filename (str): Name of the output file (default: "bitscore_statistics.tsv").

    Returns:
        str: Full path of the saved file.
    """
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Define the full file path
    file_path = os.path.join(output_directory, filename)

    # Write to TSV file
    with open(file_path, mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")

        # Write header
        writer.writerow(["Keyword", "Domain", "n", "Min", "Max", "Mean", "Median", "Std_Dev"])

        # Write data
        for keyword, domain_dict in stats_dict.items():
            for domain, stats in domain_dict.items():
                if stats:  # Skip empty entries
                    writer.writerow([
                        keyword, domain, stats["n"], stats["min"], stats["max"],
                        stats["mean"], stats["median"], stats["std_dev"]
                    ])
    
    return file_path

###########################################################################################################
# Rountines to combine specific csbs that have similarly high bitscore hits to the reference sequences

def get_highest_bitscores_for_genome(database, genome_id):
    """
    Retrieves all domains for a given genomeID and stores only the highest bitscore for each domain.

    Args:
        database (str): Path to the SQLite database.
        genome_id (str): The genomeID for which to retrieve domain information.

    Returns:
        dict: { domain: highest_bitscore }
    """
    domain_max_bitscores = {}

    with sqlite3.connect(database) as con:
        cur = con.cursor()

        query = """
        SELECT d.domain, MIN(d.score) as max_score
        FROM Domains d
        JOIN Proteins p ON d.proteinID = p.proteinID
        WHERE p.genomeID = ?
        GROUP BY d.domain;
        """

        cur.execute(query, (genome_id,))
        
        for domain, max_score in cur.fetchall():
            domain_max_bitscores[domain] = max_score

    return domain_max_bitscores

def filter_out_low_quality_csb(stats_dict, query_score_dict, threshold=0.2, min_occurrences=10):
    """
    Filters out CSBs where all domains have a median score less than the given threshold of the query reference.

    Args:
        stats_dict (dict): { keyword: { domain: { 'n': int, 'min': float, 'max': float, 
                                                  'mean': float, 'median': float, 'std_dev': float } } }
        query_score_dict (dict): { domain: max_bitscore }
        threshold (float): The percentage threshold of the query reference score.

    Returns:
        dict: Filtered CSBs where at least one domain meets the threshold.
    """
    filtered_csb = {}

    # Precompute threshold scores for each domain
    domain_thresholds = {domain: query_score * threshold for domain, query_score in query_score_dict.items()}

    for keyword, domain_dict in stats_dict.items():
        csb_should_be_removed = True  # Assume CSB will be removed unless a domain passes
        max_occurrences = max((stats["n"] for stats in domain_dict.values() if "n" in stats), default=0)

        # Skip CSB if its highest occurrence count is below the minimum required
        if max_occurrences < min_occurrences:
            continue
    
        for domain, stats in domain_dict.items():
            if domain in domain_thresholds and stats:
                # Check if at least one domain meets the threshold
                if stats["median"] >= domain_thresholds[domain]:
                    csb_should_be_removed = False
                    break  # No need to check other domains for this CSB

        if not csb_should_be_removed:
            filtered_csb[keyword] = domain_dict  # Keep CSB if at least one domain passed

    return filtered_csb

def group_keywords_by_domain(stats_dict, query_score_dict,acceptable_deviation=0.30,score_type='max'):
    """
    Compares median bitscores from stats_dict with QUERY values and groups keywords by domain 
    if the deviation is ≤ 30%. Also provides ungrouped keywords.

    Args:
        stats_dict (dict): { keyword: { domain: { 'n': int, 'min': float, 'max': float, 
                                                  'mean': float, 'median': float, 'std_dev': float } } }
        query_score_dict (dict): { domain: max_bitscore }

    Returns:
        tuple: 
            - grouped_domains (dict): { domain: [grouped_keywords] }
            - excluded_domains (dict): { domain: [excluded_keywords] }
    """
    grouped_domains = {}
    excluded_domains = {}

    for keyword, domain_dict in stats_dict.items():
        for domain, stats in domain_dict.items():
            if stats and domain in query_score_dict:
                query_score = query_score_dict[domain]  # Max score from QUERY
                median_score = stats[score_type]  # Median score for this keyword-domain

                # Calculate deviation
                deviation = abs(median_score - query_score) / query_score

                # Group keyword under domain if deviation is ≤ 30%
                if deviation <= acceptable_deviation:
                    if domain not in grouped_domains:
                        grouped_domains[domain] = [[]]
                    grouped_domains[domain][0].append(keyword)
                else:
                    if domain not in excluded_domains:
                        excluded_domains[domain] = []
                    excluded_domains[domain].append(keyword)

    return grouped_domains, excluded_domains

def group_keywords_by_domain_extended(stats_dict, query_score_dict, acceptable_deviation=0.30,score_type='max'):
    """
    Erweitert die Gruppierung von Keywords pro Domain, indem alle Proteine innerhalb eines Genclusters berücksichtigt werden,
    wenn eine der Domänen eines Genclusters das Kriterium erfüllt. Das Keyword dieses Genclusters wird dann zu allen
    Domänen im Ausgabecluster unter grouped hinzugefügt, die in diesem Gencluster enthalten sind.
    
    Args:
        stats_dict (dict): { gene_cluster: { domain: { 'n': int, 'min': float, 'max': float, 
                                                       'mean': float, 'median': float, 'std_dev': float } } }
        query_score_dict (dict): { domain: max_bitscore }
        acceptable_deviation (float): Maximal zulässige Abweichung für die Gruppierung.
    
    Returns:
        tuple: 
            - grouped_domains (dict): { domain: [grouped_keywords] }
            - excluded_domains (dict): { domain: [excluded_keywords] }
    """
    grouped_domains = {}
    excluded_domains = {}
    
    for gene_cluster, domain_dict in stats_dict.items():
        eligible_domains = set()
        for domain, stats in domain_dict.items():
            if stats and domain in query_score_dict:
                query_score = query_score_dict[domain]  # Max score from QUERY
                median_score = stats[score_type]  # Median score für dieses gene_cluster-domain Paar

                # Berechnung der Abweichung
                deviation = abs(median_score - query_score) / query_score

                if deviation <= acceptable_deviation:
                    eligible_domains.add(domain)

        if eligible_domains:
            for domain in domain_dict.keys():  # Füge das Keyword für alle Domänen im Gencluster hinzu
                if domain not in grouped_domains:
                    grouped_domains[domain] = [[]]
                if gene_cluster not in grouped_domains[domain][0]:
                    grouped_domains[domain][0].append(gene_cluster)
        else:
            for domain in domain_dict.keys():
                if domain not in excluded_domains:
                    excluded_domains[domain] = []
                if gene_cluster not in excluded_domains[domain]:
                    excluded_domains[domain].append(gene_cluster)
    
    return grouped_domains, excluded_domains
    
    

def group_excluded_keywords_by_similarity(stats_dict, excluded_domains, threshold=0.2, num_clusters=None):
    """
    Groups excluded keywords into clusters based on statistical similarity, preserving domain information.

    Args:
        stats_dict (dict): { keyword: { domain: { stats } } }
        excluded_domains (dict): { domain: [excluded_keywords] }
        threshold (float): Clustering distance threshold.
        num_clusters (int): If provided, force a fixed number of clusters.

    Returns:
        dict: { domain: [[grouped_keywords]] }  # List of lists per domain
    """
    keyword_vectors = []
    keyword_info = []  # Stores (keyword, domain) pairs
    
    # Extract statistics for excluded keywords
    for domain, keywords in excluded_domains.items():
        for keyword in keywords:
            if keyword in stats_dict and domain in stats_dict[keyword]:
                stats = stats_dict[keyword][domain]
                keyword_info.append((keyword, domain))
                keyword_vectors.append([
                    stats["median"],
                    stats["mean"],
                    stats["std_dev"]
                ])

    if not keyword_vectors:
        print("Warning: There were no keyword statistics for distinct keywords")
        return {}

    keyword_vectors = np.array(keyword_vectors)

    # Use Agglomerative Clustering for hierarchical grouping
    clustering = AgglomerativeClustering(
        n_clusters=num_clusters, metric="euclidean", linkage="ward",
        distance_threshold=threshold if num_clusters is None else None
    )
    labels = clustering.fit_predict(keyword_vectors)
    clustered_keywords = {}

    # Map cluster labels to grouped keywords
    cluster_mapping = {}
    for idx, label in enumerate(labels):
        keyword, domain = keyword_info[idx]

        if domain not in clustered_keywords:
            clustered_keywords[domain] = []  # Initialize domain with an empty list

        if label not in cluster_mapping:
            cluster_mapping[label] = len(clustered_keywords[domain])  # Assign a cluster index
            clustered_keywords[domain].append([])  # Create a new list for this cluster

        # Append the keyword to the correct cluster list

        try:
            clustered_keywords[domain][cluster_mapping[label]].append(keyword)
        except IndexError:
            print(f"Warning: Skipping domain {domain} with label {label} with keyword {keyword} due to high dissimilarity.")
            #print("Clustered Keywords")
            #pprint(clustered_keywords)
            #print("Mapping")
            #pprint(cluster_mapping)

            continue  # Move to the next keyword instead of failing

    return clustered_keywords

####################################################################################

def compute_score_limits(filtered_stats_dict, grouped_keywords_dict):
    """
    Compute upper and lower score limits per domain based on the grouped CSBs
    and their statistical properties from filtered_stats_dict.

    Args:
        filtered_stats_dict (dict): {csb: {domain: {'mean': X, 'std_dev': Y, ...}}}
        grouped_keywords_dict (dict): {domain: [[csb1, csb2, ...]]}

    Returns:
        dict: {domain: {'lower_limit': min_value, 'upper_limit': max_value}}
    """
    domain_limits = {}

    for domain, csb_groups in grouped_keywords_dict.items():
        lower_limits = []
        upper_limits = []

        # Flatten list of lists to get all CSBs for this domain
        csbs = [csb for group in csb_groups for csb in group]

        for csb in csbs:
            if csb in filtered_stats_dict and domain in filtered_stats_dict[csb]:
                stats = filtered_stats_dict[csb][domain]
                #mean = stats.get("mean", 0)
                #std_dev = stats.get("std_dev", 0)
                minimum = stats.get("min", 0)
                maximum = stats.get("max", 0)
                lower_limits.append(minimum)
                upper_limits.append(maximum)

        # Store limits only if valid statistics exist
        if lower_limits and upper_limits:
            domain_limits[domain] = {
                "lower_limit": min(lower_limits),
                "upper_limit": max(upper_limits)
            }

    return domain_limits



####################################################################################
# serilize dictionaries

    
def save_cache(directory, filename, data):
    """
    Speichert Daten in einer Datei im Pickle-Format.
    """
    with open(os.path.join(directory, filename), "wb") as f:
        pickle.dump(data, f)

def load_cache(directory, filename):
    """
    Lädt zwischengespeicherte Daten aus einer Datei, falls vorhanden.
    """
    file_path = os.path.join(directory, filename)
    if os.path.exists(file_path):
        with open(file_path, "rb") as f:
            return pickle.load(f)
    return None
    
    
    
    
    
    
    






