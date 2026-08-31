#!/usr/bin/python

import csv
import multiprocessing
import os
import sqlite3
import statistics
from typing import Any, Dict, List, Optional, Tuple

from src.core import myUtil
from src.core.logging import get_logger

logger = get_logger(__name__)


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


def get_keyword_statistics(
    database: str,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Calculate SQL-native score statistics for every keyword × domain pair.

    Returns:
        {
            keyword: {
                domain: {
                    "n": int,
                    "min": float,
                    "max": float,
                    "avg": float,
                }
            }
        }
    """

    stats_dict = {}

    query = """
        SELECT
            k.keyword,
            d.domain,
            COUNT(*) AS n,

            MIN(d.score) AS min_score,
            MAX(d.score) AS max_score,
            AVG(d.score) AS avg_score,

            MIN(d.blast_score_ratio) AS min_bsr,
            MAX(d.blast_score_ratio) AS max_bsr,
            AVG(d.blast_score_ratio) AS avg_bsr

        FROM Keywords AS k
        JOIN Proteins AS p
          ON p.clusterID = k.clusterID
        JOIN Domains AS d
          ON d.proteinID = p.proteinID

        WHERE d.domain IS NOT NULL
          AND d.score IS NOT NULL
          AND d.blast_score_ratio IS NOT NULL

        GROUP BY
            k.keyword,
            d.domain;
    """

    with sqlite3.connect(database, timeout=120.0) as con:
        cur = con.cursor()

        cur.execute("PRAGMA query_only=TRUE;")
        cur.execute("PRAGMA cache_size=-262144;")
        cur.execute("PRAGMA mmap_size=2147483648;")
        cur.execute("PRAGMA automatic_index=ON;")

        cur.execute(query)

        for (
            keyword,
            domain,
            n,
            min_score,
            max_score,
            avg_score,
            min_bsr,
            max_bsr,
            avg_bsr,
        ) in cur:
            stats_dict.setdefault(keyword, {})[domain] = {
                "n": n,
                "min": min_score,
                "max": max_score,
                "mean": round(avg_score, 2),
                "bsr_min": min_bsr,
                "bsr_max": max_bsr,
                "bsr_mean": float(avg_bsr),
            }

    return stats_dict


def _filter_out_low_quality_csb(
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
        domain: query_score * threshold  # This was previously noted as 1 - threshold
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
                if stats["max"] >= domain_thresholds[domain]:
                    csb_should_be_removed = False
                    break  # No need to check other domains for this CSB

        if not (csb_should_be_removed and keyword.startswith(csb_name_prefix)):
            filtered_csb[keyword] = (
                domain_dict  # Keep CSB if at least one domain passed
            )

    return filtered_csb


####################################################################################


def _filter_out_csb_with_protein_types(
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


def compute_cluster_stats(config: Any) -> Tuple[Dict, Dict]:
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
    filtered_stats_dict = myUtil.load_cache(config, "stat_filtered_stats.pkl")
    query_score_dict = myUtil.load_cache(config, "stat_query_score_dict.pkl")
    if filtered_stats_dict and query_score_dict:
        return filtered_stats_dict, query_score_dict

    # Get the selfblast scores to get maximum scores
    logger.info("Extracting highest bitscores per hit for query-selfblast")
    query_score_dict = get_highest_bitscores_for_genome(
        config.database_directory, "QUERY"
    )
    myUtil.save_cache(config, "stat_query_score_dict.pkl", query_score_dict)

    # Compute from scratch if not cached
    logger.info("Computing hitscore range per protein per collinear syntenic block")
    stats_dict = get_keyword_statistics(
        config.database_directory  # , options.cores
    )
    logger.info(
        f"Filtering out CSBs where all hits are below exclude_csb_hitscore {config.low_hitscore_csb_cutoff}"
    )
    filtered_stats_dict = _filter_out_low_quality_csb(
        stats_dict=stats_dict,
        query_score_dict=query_score_dict,
        threshold=config.low_hitscore_csb_cutoff,
        min_occurrences=min(10, config.min_seqs),
        csb_name_prefix=config.csb_name_prefix,
    )

    if not filtered_stats_dict:
        logger.warning(
            "No CSBs passed the strict low-quality filter. "
            "Falling back to best-per-domain CSB selection."
        )
        # TODO identify the domains that have no csb. try to find very abundant big clusters for these
        # not only if there is no csb found but also in case distant genomes from the query were used
        # Dafür braucht es eine routine die checkd welche domänen gar keinen above threshold hit in einem csb haben
        filtered_stats_dict = rescue_best_csb_per_domain(
            stats_dict,
            csb_name_prefix=config.csb_name_prefix,
            score_type="max",  # oder "median", wenn dir das lieber ist
        )

    if config.exclude_csb_proteins:
        logger.info(
            f"Filtering hits and csb on the exclusion list {config.exclude_csb_proteins}"
        )
        filtered_stats_dict = _filter_out_csb_with_protein_types(
            filtered_stats_dict, config.exclude_csb_proteins, config.csb_name_prefix
        )
    myUtil.save_cache(config, "stat_filtered_stats.pkl", filtered_stats_dict)

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
                        ]
                    )

    return file_path


###########################################################################################################


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
                logger.debug(
                    f"For gene cluster pattern {gene_cluster} added {domain}: {val} >= {thr}"
                )

        # Trage das Keyword NUR bei den passierenden Domains ein
        if not eligible_domains:
            continue
        for domain in eligible_domains:
            # genau EINE Gruppe pro Domain
            bucket = grouped.setdefault(domain, [[]])[0]
            bucket.append(gene_cluster)

    return grouped


####################################################################################


def compute_score_limits(
    filtered_stats_dict: Dict[str, Dict[str, Dict[str, float]]],
    grouped_keywords_dict: Dict[str, List[List[str]]],
) -> Dict[str, Dict[str, float]]:
    domain_limits = {}

    for domain, csb_groups in grouped_keywords_dict.items():
        # --------------------------------------------------------------
        # Bitscore statistics
        # --------------------------------------------------------------

        lower_limits = []
        upper_limits = []

        weighted_sum = 0.0
        total_n = 0

        # --------------------------------------------------------------
        # BSR statistics
        # --------------------------------------------------------------

        bsr_lower_limits = []
        bsr_upper_limits = []

        bsr_weighted_sum = 0.0
        bsr_total_n = 0

        # --------------------------------------------------------------

        csbs = [csb for group in csb_groups for csb in group]

        for csb in csbs:
            if csb not in filtered_stats_dict or domain not in filtered_stats_dict[csb]:
                continue

            stats = filtered_stats_dict[csb][domain]

            n = stats.get("n", 0)

            # ----------------------------------------------------------
            # Bitscore
            # ----------------------------------------------------------

            minimum = stats.get("min")
            maximum = stats.get("max")
            mean = stats.get("mean")

            if minimum is not None:
                lower_limits.append(float(minimum))

            if maximum is not None:
                upper_limits.append(float(maximum))

            if mean is not None and n:
                weighted_sum += float(mean) * int(n)

                total_n += int(n)

            # ----------------------------------------------------------
            # Blast Score Ratio
            # ----------------------------------------------------------

            bsr_minimum = stats.get("bsr_min")
            bsr_maximum = stats.get("bsr_max")
            bsr_mean = stats.get("bsr_mean")

            if bsr_minimum is not None:
                bsr_lower_limits.append(float(bsr_minimum))

            if bsr_maximum is not None:
                bsr_upper_limits.append(float(bsr_maximum))

            if bsr_mean is not None and n:
                bsr_weighted_sum += float(bsr_mean) * int(n)

                bsr_total_n += int(n)

        # --------------------------------------------------------------
        # Final limits
        # --------------------------------------------------------------

        if lower_limits and upper_limits:
            average = weighted_sum / total_n if total_n > 0 else None

            bsr_average = bsr_weighted_sum / bsr_total_n if bsr_total_n > 0 else None

            domain_limits[domain] = {
                # Bitscore
                "lower_limit": min(lower_limits),
                "average": average,
                "upper_limit": max(upper_limits),
                # Blast Score Ratio
                "bsr_lower_limit": (
                    min(bsr_lower_limits) if bsr_lower_limits else None
                ),
                "bsr_average": bsr_average,
                "bsr_upper_limit": (
                    max(bsr_upper_limits) if bsr_upper_limits else None
                ),
            }

    return domain_limits


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

        if not grouped_keywords:
            logger.warning(
                f"No gene cluster encoded hits with at a bitscore of at least {options.low_hitscore_csb_cutoff} of the query. Consider reducing the threshold with the -exclude_csb_score option."
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
