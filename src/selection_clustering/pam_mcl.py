#!/usr/bin/python


import sqlite3
from pathlib import Path
from typing import Dict, Set, Any, Tuple

import pandas as pd

from src.selection_defragmentation import pam_mx_algorithm, predictor
from src.selection_seed import csb_proteins_selection
from src.core import myUtil

from src.core.logging import get_logger

logger = get_logger(__name__)


def _select_ref_seq_mcl_sequences(
    mcl_file: str, reference_sequences: Set[str]
) -> Set[str]:
    """
    Returns all protein IDs from MCL clusters that contain at least one reference sequence.

    Args:
        mcl_file (str): Path to the MCL output file.
        domain (str): Domain label.
        reference_sequences (set): Set of known reference protein IDs.

    Returns:
        set: All protein IDs from clusters with at least one reference protein.

    Example:
        select_ref_seq_mcl_sequences("a_mcl_clusters.txt", "A", {"ref1","ref2"})
    """

    selected: set[str] = set()

    with open(mcl_file, "r") as handle:
        for line in handle:
            if not line.strip():
                continue

            cluster = set(line.split())

            if cluster & reference_sequences:
                selected.update(cluster)

    return selected


def select_seqs_with_truncated_csb_vicinity(
    database_path: str,
    protein_ids: Set[str],
    common_domains: Set[frozenset],
    chunk_size: int = 900,
) -> Set[str]:
    """
    Selects proteins whose cluster neighborhood only contains domains from the given sets.

    Args:
        database_path (str): SQLite DB path.
        protein_ids (set): Input protein IDs.
        common_domains (set): Set of acceptable domain sets (frozenset per neighborhood).
        chunk_size (int): SQL chunk size.

    Returns:
        set: ProteinIDs with valid neighborhood only.

    Example:
        select_seqs_with_truncated_csb_vicinity("db.sqlite", {"p1"}, {frozenset({"D1", "D2"})})
    """

    selected = set()

    def chunked(iterable, size):
        items = list(iterable)

        for start in range(0, len(items), size):
            yield items[start : start + size]

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        # Step 1: Fetch clusterIDs for proteinIDs in chunks
        pid_to_cluster = {}
        for chunk in chunked(protein_ids, chunk_size):
            placeholders = ",".join("?" for _ in chunk)
            cur.execute(
                f"""
                SELECT proteinID, clusterID 
                FROM Proteins 
                WHERE proteinID IN ({placeholders})
                  AND clusterID IS NOT NULL
            """,
                tuple(chunk),
            )
            pid_to_cluster.update(cur.fetchall())

        if not pid_to_cluster:
            return selected

        cluster_ids = set(pid_to_cluster.values())

        # Step 2: Fetch all proteinIDs from those clusters in chunks
        cluster_to_proteins = {}
        for chunk in chunked(cluster_ids, chunk_size):
            placeholders = ",".join("?" for _ in chunk)
            cur.execute(
                f"""
                SELECT clusterID, proteinID 
                FROM Proteins 
                WHERE clusterID IN ({placeholders})
            """,
                tuple(chunk),
            )
            for clusterID, proteinID in cur.fetchall():
                cluster_to_proteins.setdefault(clusterID, set()).add(proteinID)

        # Step 3: Fetch all domains for all involved proteinIDs in chunks
        all_neighbor_proteins = {p for ps in cluster_to_proteins.values() for p in ps}
        protein_to_domains = {}
        for chunk in chunked(all_neighbor_proteins, chunk_size):
            placeholders = ",".join("?" for _ in chunk)
            cur.execute(
                f"""
                SELECT proteinID, domain 
                FROM Domains 
                WHERE proteinID IN ({placeholders})
            """,
                tuple(chunk),
            )
            for proteinID, domain in cur.fetchall():
                protein_to_domains.setdefault(proteinID, set()).add(domain)

        # Step 4: Check each input protein for domain validity in its cluster
        for pid in protein_ids:
            cluster = pid_to_cluster.get(pid)
            if not cluster:
                continue
            neighbor_proteins = cluster_to_proteins.get(cluster, set())
            domains_in_cluster = set()
            for neighbor in neighbor_proteins:
                domains_in_cluster.update(protein_to_domains.get(neighbor, set()))

            # Vergleich: domains_in_cluster ⊆ common_domains[i]
            if any(
                domains_in_cluster.issubset(domain_set) for domain_set in common_domains
            ):
                selected.add(pid)

    return selected


#############################################################################################################


def select_gene_cluster_vicinity_domains(
    db_path: str, hit_ids: Set[str]
) -> Tuple[Set[frozenset], Dict]:
    """
    For a set of protein IDs, extracts their cluster neighborhoods (as frozensets of domains).

    Args:
        db_path (str): SQLite DB path.
        hit_ids (set): Protein IDs.

    Returns:
        tuple:
            - set of frozensets (unique gene clusters as sets of domains)
            - neighbors dict (for possible debugging)
    """

    neighbors, _ = csb_proteins_selection.fetch_neighbouring_genes_with_domains(
        db_path, hit_ids
    )

    # Collapse the neighbourhoods by similarity

    unique_clusters = set()

    for cluster in neighbors.values():
        # Without genomic context no entry
        if cluster == [["singleton"]]:
            continue
        # Gene cluster in correct order
        flattened = frozenset(domain for gene in cluster for domain in gene)
        unique_clusters.add(flattened)

    return unique_clusters, neighbors


#############################################################################################################


def log_all_mcl_cluster_statistics(
    reference_dict: Dict[str, Set[str]],
    cluster_proteins_dict: Dict[str, Set[str]],
    grouped_3_dict: Dict[str, Set[str]],
    grouped_4_dict: Dict[str, Set[str]],
) -> None:
    """
    Logs selection statistics for all domains, reporting merged set composition.

    Args:
        reference_dict (dict): {domain: set(reference seqs)}
        cluster_proteins_dict (dict): {domain: set(cluster hits)}
        grouped_3_dict (dict): {domain: set(by CSB)}
        grouped_4_dict (dict): {domain: set(by PAM)}

    Returns:
        None. Logs detailed stats.
    """
    for domain in reference_dict:
        reference_sequences = reference_dict.get(domain, set())
        mcl_cluster_protID_set = cluster_proteins_dict.get(domain, set())
        group3_set = grouped_3_dict.get(domain, set())
        group4_set = grouped_4_dict.get(domain, set())

        total_clusters_with_refs = len(mcl_cluster_protID_set)
        total_reference_count = len(reference_sequences)
        true_new_hits = mcl_cluster_protID_set - reference_sequences
        total_new_hits = len(true_new_hits)

        # Final merged set
        final_merged_set = reference_sequences | group3_set | group4_set

        # Detaillierte Zusammensetzung
        ref_part = final_merged_set & reference_sequences
        grp3_grp4_part = (
            final_merged_set & group3_set & group4_set - reference_sequences
        )
        grp3_part = final_merged_set & group3_set - reference_sequences - group4_set
        grp4_part = final_merged_set & group4_set - reference_sequences - group3_set

        logger.debug(f"\n[{domain}]")
        logger.debug(
            f"  Reference sequences:                 {len(reference_sequences)}"
        )
        logger.debug(
            f"  Total sequences in ref-mcl-clusters: {total_clusters_with_refs}"
        )
        logger.debug(f"  New candidate sequences:             {total_new_hits}")
        logger.debug(f"  Selected by (CSB):                   {len(group3_set)}")
        logger.debug(f"  Selected by (PAM):                   {len(group4_set)}")
        logger.debug(f"  Final set:                           {len(final_merged_set)}")
        logger.debug(f"     ├─ from references:               {len(ref_part)}")
        logger.debug(f"     ├─ from (CSB) & (PAM):            {len(grp3_grp4_part)}")
        logger.debug(f"     ├─ from (CSB) only:               {len(grp3_part)}")
        logger.debug(f"     └─ from (PAM) only:               {len(grp4_part)}")


def _select_truncated_csb_hits(
    config: Any,
    clustering_results: Dict[str, str],
    basis_grouped: Dict[str, Set[str]],
) -> tuple[
    Dict[str, Set[str]],
    Dict[str, Set[str]],
]:
    truncated_csb_hits = myUtil.load_cache(
        config,
        "mcl_truncated_csb_hits.pkl",
    )

    if truncated_csb_hits:
        return truncated_csb_hits, {}

    truncated_csb_hits = {}
    cluster_candidates = {}

    total = len(clustering_results)

    for index, (domain, cluster_file) in enumerate(
        clustering_results.items(),
        start=1,
    ):
        reference_sequences = basis_grouped.get(domain, set())

        if not reference_sequences:
            logger.debug(
                "No reference sequences found for domain '%s' - skipping",
                domain,
            )
            continue

        logger.debug(
            "[%d/%d] Processing domain %s",
            index,
            total,
            domain,
        )

        common_gene_vicinity = _get_common_gene_vicinity(
            config=config,
            domain=domain,
            reference_sequences=reference_sequences,
        )

        cluster_protein_ids = _select_ref_seq_mcl_sequences(
            mcl_file=cluster_file,
            reference_sequences=reference_sequences,
        )

        cluster_candidates[domain] = cluster_protein_ids

        candidate_ids = cluster_protein_ids - reference_sequences

        truncated_csb_hits[domain] = select_seqs_with_truncated_csb_vicinity(
            config.database_directory,
            candidate_ids,
            common_gene_vicinity,
        )

    myUtil.save_cache(
        config,
        "mcl_truncated_csb_hits.pkl",
        truncated_csb_hits,
    )

    return truncated_csb_hits, cluster_candidates


def _get_common_gene_vicinity(
    config: Any,
    domain: str,
    reference_sequences: Set[str],
) -> Set[frozenset]:
    cache_name = f"mcl_common_gene_vicinity_{domain}.pkl"

    common_gene_vicinity = myUtil.load_cache(
        config,
        cache_name,
    )

    if common_gene_vicinity:
        return common_gene_vicinity

    common_gene_vicinity, neighbors = select_gene_cluster_vicinity_domains(
        config.database_directory,
        reference_sequences,
    )

    myUtil.save_cache(
        config,
        cache_name,
        common_gene_vicinity,
    )

    myUtil.save_cache(
        config,
        f"mcl_gene_vicinity_dict_{domain}.pkl",
        neighbors,
    )

    return common_gene_vicinity


def _select_pam_hits(
    config: Any,
    basis_seed_sequences: Dict[str, Set[str]],
    basis_score_limit: Dict[str, Dict[str, float]],
) -> Dict[str, Set[str]]:
    pam_hits = myUtil.load_cache(
        config,
        "mcl_PAM_plausible_hits.pkl",
    )

    if pam_hits:
        return pam_hits

    logger.info(
        "Selecting sequences with plausible genomic co-occurrence"
    )

    pam_hits = predictor.predictor_training_calibration_application(
        config=config,
        basis_seed_sequences=basis_seed_sequences,
        basis_score_limit=basis_score_limit,
        probability_cutoff=config.pam_threshold,
        support_models_name="grp3_support_models.pkl",
    )

    myUtil.save_cache(
        config,
        "mcl_PAM_plausible_hits.pkl",
        pam_hits,
    )

    return pam_hits


def select_hits_by_pam_csb_mcl(
    config: Any,
    clustering_results: Dict[str, str],
    basis_seed_sequences: Dict[str, Set[str]],
    basis_score_limit: Dict[str, Dict[str, float]],
) -> Dict[str, Set[str]]:

    cached_result = myUtil.load_cache(
        config,
        "mcl_PAM_csb_merged_hits.pkl",
    )
    if cached_result:
        return cached_result

    truncated_csb_hits, cluster_candidates = _select_truncated_csb_hits(
        config=config,
        clustering_results=clustering_results,
        basis_grouped=basis_seed_sequences,
    )

    pam_seed_sequences = csb_proteins_selection.merge_protein_sets(
        basis_seed_sequences,
        truncated_csb_hits,
    )

    pam_hits = _select_pam_hits(
        config=config,
        basis_seed_sequences=pam_seed_sequences,
        basis_score_limit=basis_score_limit,
    )

    log_all_mcl_cluster_statistics(
        reference_dict=basis_seed_sequences,
        cluster_proteins_dict=cluster_candidates,
        grouped_3_dict=truncated_csb_hits,
        grouped_4_dict=pam_hits,
    )

    merged_grouped = csb_proteins_selection.merge_protein_sets(
        basis_seed_sequences,
        truncated_csb_hits,
        pam_hits,
    )

    myUtil.save_cache(
        config,
        "mcl_PAM_csb_merged_hits.pkl",
        merged_grouped,
    )

    return merged_grouped