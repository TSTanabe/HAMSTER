#!/usr/bin/python

import statistics
import os
import sqlite3
import csv
import multiprocessing
import numpy as np

from typing import Dict, Any, List, Tuple, Optional
from sklearn.cluster import AgglomerativeClustering

from src.core import myUtil
from src.core.logging import get_logger

logger = get_logger(__name__)


def group_gene_cluster_statistic(options: Any):
    """
    Backward-compatible API:
        Returns domain_score_limits, filtered_stats_dict, grouped_keywords, clustered_excluded_keywords

    Internally runs:
        (1) compute_cluster_stats()
        (2) apply_cluster_selection()
    """
    filtered_stats_dict, query_score_dict = compute_cluster_stats(options)
    domain_score_limits, grouped_keywords = apply_cluster_selection(
        options, filtered_stats_dict, query_score_dict
    )

    return domain_score_limits, filtered_stats_dict, grouped_keywords


def compute_cluster_stats(options: Any) -> Tuple[Dict, Dict]:
    """
    Phase 1: Compute the statistical basis for cluster selection.
      - Per keyword × domain hitscore statistics
      - QUERY reference scores for each domain
    Saves intermediate results to cache.

    Returns:
        filtered_stats_dict (dict): {keyword: {domain: {'n','min','max','mean','median','std_dev'}}}
        query_score_dict    (dict): {domain: best_query_score}
    """
    # Load caches if available
    filtered_stats_dict = myUtil.load_cache(options, "stat_filtered_stats.pkl")
    query_score_dict = myUtil.load_cache(options, "stat_query_score_dict.pkl")

    # Compute from scratch if not cached
    if not filtered_stats_dict or not query_score_dict:
        logger.info("Computing hitscore range per protein per collinear syntenic block")
        stats_dict = get_keyword_statistics_parallel(
            options.database_directory, options.cores
        )

        logger.info("Extracting highest bitscores per hit for query-selfblast")
        query_score_dict = get_highest_bitscores_for_genome(
            options.database_directory, "QUERY"
        )
        myUtil.save_cache(options, "stat_query_score_dict.pkl", query_score_dict)

        logger.info(
            f"Filtering out CSBs where all hits are below exclude_csb_hitscore {options.low_hitscore_csb_cutoff}"
        )
        filtered_stats_dict = filter_out_low_quality_csb(
            stats_dict,
            query_score_dict,
            options.low_hitscore_csb_cutoff,
            min(10, options.min_seqs),
            options.csb_name_prefix,
        )

        if not filtered_stats_dict:
            logger.warning(
                "No CSBs passed the strict low-quality filter. "
                "Falling back to best-per-domain CSB selection."
            )
            filtered_stats_dict = rescue_best_csb_per_domain(
                stats_dict,
                csb_name_prefix=options.csb_name_prefix,
                score_type="max",  # oder "median", wenn dir das lieber ist
            )

        logger.info(
            f"Filtering out CSBs where a hit is on the exclusion list {options.exclude_csb_proteins}"
        )
        filtered_stats_dict = filter_out_csb_with_protein_types(
            filtered_stats_dict, options.exclude_csb_proteins, options.csb_name_prefix
        )

        myUtil.save_cache(options, "stat_filtered_stats.pkl", filtered_stats_dict)

    return filtered_stats_dict, query_score_dict


def rescue_best_csb_per_domain(
    stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    csb_name_prefix: str = "csb-",
    score_type: str = "max",
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Fallback-Routine, falls der strenge Filter keine CSBs übrig lässt.

    Idee:
        Für jede Domain (Protein-Typ) wird der CSB mit dem höchsten Hitscore
        (standardmäßig stats[score_type] = 'max') ausgewählt.
        Alle ausgewählten CSBs werden zurückgegeben.

    Args:
        stats_dict:
            {keyword: {domain: {'n','min','max','mean','median','std_dev'}}}
        csb_name_prefix:
            Optional: Nur CSBs berücksichtigen, deren Keyword mit diesem Prefix beginnt.
        score_type:
            Welcher Statistik-Wert als "Hitscore" gewertet wird ('max', 'median', ...).

    Returns:
        dict:
            Gefiltertes stats_dict, das nur die „best-per-domain“-CSBs enthält.
            {keyword: {domain: {...}}}
    """

    if not stats_dict:
        return {}

    # Für jede Domain: bester CSB und zugehöriger Score
    # domain -> (best_keyword, best_score)
    best_csb_for_domain: Dict[str, Tuple[str, float]] = {}

    for keyword, domain_dict in stats_dict.items():
        # Optional: Nur CSBs mit Prefix betrachten
        if csb_name_prefix and not keyword.startswith(csb_name_prefix):
            continue

        for domain, stats in domain_dict.items():
            if not stats:
                continue
            val = stats.get(score_type)
            if val is None:
                continue

            current = best_csb_for_domain.get(domain)
            if current is None or val > current[1]:
                best_csb_for_domain[domain] = (keyword, val)

    if not best_csb_for_domain:
        # Es gibt nichts Sinnvolles zu retten
        return {}

    # Jetzt alle CSBs einsammeln, die für mindestens eine Domain "beste" sind
    rescued: Dict[str, Dict[str, Dict[str, float]]] = {}
    for domain, (keyword, _) in best_csb_for_domain.items():
        # komplette Domain-Statistik dieses CSB übernehmen
        rescued[keyword] = stats_dict[keyword]

    return rescued


def apply_cluster_selection(
    options: Any, filtered_stats_dict: Dict, query_score_dict: Dict
):
    """
    Phase 2: Apply selection rules based on precomputed statistics.
      - Save stats to TSV
      - Group keywords to domains (extended mode)
      - Optionally cluster 'excluded' CSBs
      - Compute domain score limits
    Saves results to cache.

    Returns:
        domain_score_limits         (dict)
        grouped_keywords            (dict): {domain: [[keywords]]} (extended)
        clustered_excluded_keywords (dict)
    """
    # Load caches if available
    domain_score_limits = myUtil.load_cache(options, "stat_domain_score_limits.pkl")
    grouped_keywords = myUtil.load_cache(options, "stat_grouped_keywords.pkl")

    # Always save TSV for filtered statistics
    save_stats_to_tsv(filtered_stats_dict, options.Csb_directory)

    # Extended grouping: if a CSB passes, its keyword is added to all domains in the CSB
    if not grouped_keywords:
        logger.info(
            "Selecting CSB patterns that encode for highly similar query homologs"
        )

        grouped_keywords = group_keywords_by_domain_passers(
            filtered_stats_dict, query_score_dict, options.low_hitscore_csb_cutoff
        )
        # The following routine includes all domains from a csb were a single domain passes
        # grouped_keywords = group_keywords_by_domain_extended(
        #    filtered_stats_dict,
        #    query_score_dict,
        #    options.low_hitscore_csb_cutoff
        # )
        myUtil.save_cache(options, "stat_grouped_keywords.pkl", grouped_keywords)

    # Compute score limits per domain based on grouped keywords
    if not domain_score_limits:
        logger.info(
            "Computing for each protein hit score range across all selected CSBs"
        )
        domain_score_limits = compute_score_limits(
            filtered_stats_dict, grouped_keywords
        )
        myUtil.save_cache(options, "stat_domain_score_limits.pkl", domain_score_limits)

    return domain_score_limits, grouped_keywords


#################################################################################
# Fetch for each csb for each domain all the possible scores. Calculate the min,max,mean,median and std_dev and print them


def wrapped_fetch_keyword_scores(
    database: str,
    chunk: List[str],
    progress_counter: multiprocessing.Value,
    lock: multiprocessing.Lock,
) -> Dict[str, Dict[str, List[float]]]:
    """
    Fetches keyword scores and updates progress for multiprocessing.

    Args:
        database (str): Path to SQLite database.
        chunk (list): List of keywords for this worker.
        progress_counter (multiprocessing.Value): Shared counter for progress.
        lock (multiprocessing.Lock): For progress.

    Returns:
        dict: {keyword: {domain: [scores]}}
    """

    result = fetch_keyword_scores(database, chunk)
    with lock:
        progress_counter.value += 1
        print(
            f"Progress: {progress_counter.value} cluster statistics completed", end="\r"
        )
    return result


def get_keyword_statistics_parallel(
    database: str, num_workers: int = 4
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Parallelized version of keyword statistics calculation.

    Args:
        database (str): Path to SQLite database.
        num_workers (int): Number of parallel workers.

    Returns:
        dict: {keyword: {domain: {'n', 'min', 'max', 'mean', 'median', 'std_dev'}}}
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("SELECT DISTINCT keyword FROM Keywords;")
        keywords = [row[0] for row in cur.fetchall()]

    # Split keywords into chunks for parallel processing
    chunk_size = max(1, len(keywords) // num_workers)  # Avoid division by zero
    keyword_chunks = [
        keywords[i : i + chunk_size] for i in range(0, len(keywords), chunk_size)
    ]

    # Progress tracking
    manager = multiprocessing.Manager()
    progress_counter = manager.Value("i", 0)  # Shared counter for completed workers
    lock = manager.Lock()  # Lock to prevent race conditions

    # Use multiprocessing pool with progress tracking
    with multiprocessing.Pool(num_workers) as pool:
        raw_results = pool.starmap(
            wrapped_fetch_keyword_scores,
            [(database, chunk, progress_counter, lock) for chunk in keyword_chunks],
        )

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


def fetch_keyword_scores(
    database: str, keyword_list: List[str]
) -> Dict[str, Dict[str, List[float]]]:
    """
    Fetches raw domain scores for a list of keywords.

    Args:
        database (str): Path to SQLite DB.
        keyword_list (list): List of keywords.

    Returns:
        dict: {keyword: {domain: [scores]}}
    """
    raw_scores: Dict[str, Dict[str, List[float]]] = {}
    if not keyword_list:
        return raw_scores

    with sqlite3.connect(database) as con:
        cur = con.cursor()
        query = f"""
            SELECT k.keyword, d.domain, d.score
            FROM Keywords k
            JOIN Domains d ON k.clusterID = (SELECT clusterID FROM Proteins WHERE proteinID = d.proteinID)
            WHERE k.keyword IN ({",".join(["?"] * len(keyword_list))})
            AND d.domain IS NOT NULL
            AND d.score IS NOT NULL;
        """
        cur.execute(query, keyword_list)
        for keyword, domain, score in cur.fetchall():
            raw_scores.setdefault(keyword, {}).setdefault(domain, []).append(score)
    return raw_scores


def compute_statistics_for_keywords(
    raw_scores: Dict[str, Dict[str, List[float]]],
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Computes count, min, max, mean, median, std_dev for each domain per keyword.

    Args:
        raw_scores (dict): {keyword: {domain: [scores]}}

    Returns:
        dict: {keyword: {domain: {'n', 'min', 'max', 'mean', 'median', 'std_dev'}}}
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
                    "std_dev": round(statistics.stdev(scores), 2)
                    if len(scores) > 1
                    else 0.0,
                }

    return stats_dict


#########################################################################################


def save_stats_to_tsv(
    stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    output_directory: str,
    filename: str = "bitscore_statistics.tsv",
) -> str:
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
        writer.writerow(
            ["Keyword", "Domain", "n", "Min", "Max", "Mean", "Median", "Std_Dev"]
        )

        # Write data
        for keyword, domain_dict in stats_dict.items():
            for domain, stats in domain_dict.items():
                if stats:  # Skip empty entries
                    writer.writerow(
                        [
                            keyword,
                            domain,
                            stats["n"],
                            stats["min"],
                            stats["max"],
                            stats["mean"],
                            stats["median"],
                            stats["std_dev"],
                        ]
                    )

    return file_path


###########################################################################################################
# Rountines to combine specific csbs that have similarly high bitscore hits to the reference sequences


def get_highest_bitscores_for_genome(database: str, genome_id: str) -> Dict[str, float]:
    """
    Gets the highest bitscore for each domain for a given genome.

    Args:
        database (str): Path to SQLite DB.
        genome_id (str): GenomeID.

    Returns:
        dict: {domain: highest_bitscore}
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


def filter_out_low_quality_csb(
    stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    query_score_dict: Dict[str, float],
    threshold: float = 0.2,
    min_occurrences: int = 10,
    csb_name_prefix: str = "csb-",
) -> Dict[str, Dict[str, Dict[str, float]]]:
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
    domain_thresholds = {
        domain: query_score * (1 - threshold)
        for domain, query_score in query_score_dict.items()
    }

    for keyword, domain_dict in stats_dict.items():
        csb_should_be_removed = (
            True  # Assume CSB will be removed unless a domain passes
        )
        max_occurrences = max(
            (stats["n"] for stats in domain_dict.values() if "n" in stats), default=0
        )

        # Skip CSB if its highest occurrence count is below the minimum required
        if max_occurrences < min_occurrences:
            continue

        for domain, stats in domain_dict.items():
            if domain in domain_thresholds and stats:
                # Check if at least one domain meets the threshold
                if stats["median"] >= domain_thresholds[domain]:
                    csb_should_be_removed = False
                    break  # No need to check other domains for this CSB

        if not (csb_should_be_removed and keyword.startswith(csb_name_prefix)):
            filtered_csb[keyword] = (
                domain_dict  # Keep CSB if at least one domain passed
            )

    return filtered_csb


def filter_out_csb_with_protein_types(
    stats_dict, exclude_types: list, csb_name_prefix: str
) -> dict:
    """
    Removes all CSBs (keywords) from stats_dict that contain any domain listed in exclude_types.

    Args:
        stats_dict (dict): {keyword: {domain: {...}}}
        exclude_types (list): List of domains/types to exclude.

    Returns:
        dict: Filtered stats_dict.
    """
    filtered = {}
    exclude_set = set(exclude_types)
    for csb_keyword, domain_dict in stats_dict.items():
        domains_in_csb = set(domain_dict.keys())
        intersection = domains_in_csb & exclude_set

        # intersection means excluded proteins were found
        if intersection and csb_keyword.startswith(csb_name_prefix):
            logger.debug(f"Removed csb {csb_keyword} due a hit for {intersection}")
        elif intersection and not csb_keyword.startswith(csb_name_prefix):
            logger.debug(
                f"Kept pattern file defined csb {csb_keyword} with hit for {intersection}"
            )
            filtered[csb_keyword] = domain_dict
        else:
            filtered[csb_keyword] = domain_dict

    return filtered


def group_keywords_by_domain_extended_deprecated(
    stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    query_score_dict: Dict[str, float],
    acceptable_deviation: float = 0.30,
    score_type: str = "max",
) -> Tuple[Dict[str, List[List[str]]], Dict[str, List[str]]]:
    """
    DEPRECATED BECAUSE OPTIMIZED ROUTINE WAS WRITTEN DIRECTLY AFTER THIS ROUTINE
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
    """
    grouped_domains = {}

    for gene_cluster, domain_dict in stats_dict.items():
        eligible_domains = set()
        for domain, stats in domain_dict.items():
            if stats and domain in query_score_dict:
                query_score = query_score_dict[domain]  # Max score from QUERY
                median_score = stats[
                    score_type
                ]  # Median score für dieses gene_cluster-domain Paar

                # Berechnung der Abweichung
                deviation = abs(median_score - query_score) / query_score

                if deviation <= (1 - acceptable_deviation):
                    eligible_domains.add(domain)

        if eligible_domains:
            for domain in (
                domain_dict.keys()
            ):  # Füge das Keyword für alle Domänen im Gencluster hinzu
                if domain not in grouped_domains:
                    grouped_domains[domain] = [[]]
                if gene_cluster not in grouped_domains[domain][0]:
                    grouped_domains[domain][0].append(gene_cluster)

    return grouped_domains


def group_keywords_by_domain_extended(
    stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    query_score_dict: Dict[str, float],
    acceptable_deviation: float = 0.30,
    score_type: str = "max",
) -> Dict[str, List[List[str]]]:
    """
    Nimmt ein CSB (gene_cluster), prüft: erreicht *irgendeine* Domain im CSB den
    Schwellenwert (relativ zum QUERY-Score)? Wenn ja, wird das Keyword (gene_cluster)
    bei *allen* Domains des CSB eingetragen.

    Rückgabe:
        { domain: [[keyword1, keyword2, ...]] }  # genau EINE Gruppe je Domain
    """
    if not stats_dict or not query_score_dict:
        return {}

    # Precompute Schwellen je Domain (z. B. 70% des Query-Referenzscores bei 0.30 Abweichung)
    thresholds = {
        d: s * (1 - acceptable_deviation) for d, s in query_score_dict.items()
    }
    if not thresholds:
        return {}

    grouped: Dict[str, List[List[str]]] = {}

    for gene_cluster, domain_stats in stats_dict.items():
        # Prüfen, ob *irgendeine* Domain in diesem CSB die Schwelle erfüllt
        passes = False
        for d, st in domain_stats.items():
            thr = thresholds.get(d)
            if thr is None or not st:
                continue
            val = st.get(score_type)
            if val is not None and val >= thr:
                passes = True
                break

        if not passes:
            continue

        # Wenn ja: Keyword bei *allen* Domains dieses CSB eintragen
        for d in domain_stats.keys():
            bucket = grouped.setdefault(d, [[]])[0]
            bucket.append(gene_cluster)

    return grouped


def group_keywords_by_domain_passers(
    stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    query_score_dict: Dict[str, float],
    acceptable_deviation: float = 0.30,
    score_type: str = "max",
) -> Dict[str, List[List[str]]]:
    """
    Pro Domain werden nur jene Keywords (CSBs) gesammelt, bei denen
    diese Domain selbst den Schwellenwert erreicht.

    Rückgabeformat kompatibel zu downstream:
        { domain: [[keyword1, keyword2, ...]] }  # genau EINE Gruppe je Domain
    """
    if not stats_dict or not query_score_dict:
        return {}

    # Precompute: domain-spezifische Schwellen (z. B. 70% des Query-Scores bei acceptable_deviation=0.30)
    thresholds = {d: s * acceptable_deviation for d, s in query_score_dict.items()}
    if not thresholds:
        return {}

    grouped: Dict[str, List[List[str]]] = {}

    # Iteriere über CSBs (gene_cluster = Keyword)
    for gene_cluster, domain_stats in stats_dict.items():
        # Finde alle Domains in diesem CSB, die selbst >= cutoff sind
        # (nur Domains berücksichtigen, für die es einen Query-Referenzscore gibt)
        eligible_domains = []
        for domain, stats in domain_stats.items():
            thr = thresholds.get(domain)
            if thr is None or not stats:
                continue
            val = stats.get(score_type)
            if val is None:
                continue
            if val >= thr:
                eligible_domains.append(domain)

        # Trage das Keyword NUR bei den passierenden Domains ein
        if not eligible_domains:
            continue
        for domain in eligible_domains:
            # genau EINE Gruppe pro Domain
            bucket = grouped.setdefault(domain, [[]])[0]
            bucket.append(gene_cluster)

    return grouped


def group_excluded_keywords_by_similarity(
    stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    excluded_domains: Dict[str, List[str]],
    threshold: float = 0.2,
    num_clusters: Optional[int] = None,
) -> Dict[str, List[List[str]]]:
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
                keyword_vectors.append(
                    [stats["median"], stats["mean"], stats["std_dev"]]
                )

    if not keyword_vectors:
        logger.warning("There were no keyword statistics for distinct keywords")
        return {}

    keyword_vectors = np.array(keyword_vectors)

    # Use Agglomerative Clustering for hierarchical grouping
    clustering = AgglomerativeClustering(
        n_clusters=num_clusters,
        metric="euclidean",
        linkage="ward",
        distance_threshold=threshold if num_clusters is None else None,
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
            cluster_mapping[label] = len(
                clustered_keywords[domain]
            )  # Assign a cluster index
            clustered_keywords[domain].append([])  # Create a new list for this cluster

        # Append the keyword to the correct cluster list

        try:
            clustered_keywords[domain][cluster_mapping[label]].append(keyword)
        except IndexError:
            logger.warning(
                f"Skipping {domain} with label {label} and keyword {keyword} due to high dissimilarity."
            )
            continue

    return clustered_keywords


####################################################################################


def compute_score_limits(
    filtered_stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    grouped_keywords_dict: Dict[str, List[List[str]]],
) -> Dict[str, Dict[str, float]]:
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
                # mean = stats.get("mean", 0)
                # std_dev = stats.get("std_dev", 0)
                minimum = stats.get("min", 0)
                maximum = stats.get("max", 0)
                lower_limits.append(minimum)
                upper_limits.append(maximum)

        # Store limits only if valid statistics exist
        if lower_limits and upper_limits:
            domain_limits[domain] = {
                "lower_limit": min(lower_limits),
                "upper_limit": max(upper_limits),
            }

    return domain_limits
