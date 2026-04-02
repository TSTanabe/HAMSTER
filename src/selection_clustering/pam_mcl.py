#!/usr/bin/python


import sqlite3
from pathlib import Path
from typing import Dict, Set, Any, Tuple

import pandas as pd

from src.selection_defragmentation import pam_mx_algorithm, pam_defragmentation
from src.selection_seed import csb_proteins_selection
from src.core import myUtil

from src.core.logging import get_logger

logger = get_logger(__name__)


def select_hits_by_pam_csb_mcl(
    options: Any,
    mcl_clustering_results_dict: Dict[str, str],
    basis_grouped: Dict[str, Set[str]],
) -> Dict[str, Set[str]]:
    """
    ### Main routine to the module ###

    Get all Proteins from MCL clusters with a reference sequence

    Get the csb vicinity from all non-ref seqs => definiert durch die keys in grouped. Diese kommen schließlich aus Genclustern

    Get the csb vinitiy from all ref-seqs => Müssen aus der Datenbank gelesen werden

    Select all non-ref seqs that have a partial ref-seq genomic vicinity into the extended_refseqs

    Wie bisher auch get the sequences that are plausible by liberal selection into the extended_refseqs



    Args:
        options (Any): Config/options object.
        mcl_clustering_results_dict (dict): {domain: mcl_cluster_file_path}
        basis_grouped (dict): {domain: set(reference protein IDs)}

    Returns:
        dict: {domain: set(protein IDs)} merged final set.
    """

    grouped_3_dict = myUtil.load_cache(options, "mcl_truncated_csb_hits.pkl")
    grouped_4_dict = myUtil.load_cache(options, "mcl_PAM_plausible_hits.pkl")
    extended_grouped = myUtil.load_cache(options, "mcl_PAM_csb_merged_hits.pkl")

    if extended_grouped:
        return extended_grouped
    else:
        extended_grouped = {}

    # Process reference_dict to extract the actual domain names
    processed_reference_dict = {
        key.split("_", 1)[-1]: value for key, value in basis_grouped.items()
    }

    total = len(mcl_clustering_results_dict)

    if not grouped_3_dict:
        grouped_3_dict = {}
        logger.info("Selecting sequences from mcl clusters with truncated csb pattern")
        for i, (domain, mcl_file) in enumerate(mcl_clustering_results_dict.items(), 1):
            p = Path(mcl_file)
            short_path = Path(*p.parts[-3:])  # letzte 2 Ordner + Datei
            logger.debug(f"[{i}/{total}] Processing: {domain} in file {short_path}")

            # Get reference sequences for the domain (if exists in processed reference dict)
            reference_sequences = processed_reference_dict.get(domain, set())
            if not reference_sequences:
                logger.warning(
                    f"No reference sequences found for domain '{domain}' - skipping"
                )
                continue

            # Get the common gene cluster vicinity (presence without order or doublication)
            common_gene_vicinity = myUtil.load_cache(
                options, f"mcl_common_gene_vicinity_{domain}.pkl"
            )
            if not common_gene_vicinity:
                common_gene_vicinity, neighbors_dict = (
                    select_gene_cluster_vicinity_domains(
                        options.database_directory, reference_sequences
                    )
                )  # Neighbors dict ist proteinID => vicinity
                myUtil.save_cache(
                    options,
                    f"mcl_common_gene_vicinity_{domain}.pkl",
                    common_gene_vicinity,
                )
                myUtil.save_cache(
                    options, f"mcl_gene_vicinity_dict_{domain}.pkl", neighbors_dict
                )

            mcl_cluster_protID_set = select_ref_seq_mcl_sequences(
                mcl_file, domain, reference_sequences
            )

            # Save for later use in the PAM calculation
            extended_grouped[domain] = mcl_cluster_protID_set

            # only the new proteinIDs without the refseq are processed for common csb
            new_proteinID_set = mcl_cluster_protID_set - reference_sequences

            # Select the proteins with a truncated csb, that fits the other csbs
            filtered_proteinIDs = select_seqs_with_truncated_csb_vicinity(
                options.database_directory, new_proteinID_set, common_gene_vicinity
            )
            grouped_3_dict[domain] = filtered_proteinIDs
    myUtil.save_cache(options, "mcl_truncated_csb_hits.pkl", grouped_3_dict)


    logger.info(
        f"Selecting sequences from with plausible co-occurence in the genome"
    )
    # Select the proteins with plausible PAM
    grouped_3_dict = myUtil.merge_grouped_refseq_dicts_simple(
        grouped_3_dict, basis_grouped
    ) # Full representation of all sequences after the group 3 addition
    if not grouped_4_dict:
        grouped_4_dict = pam_defragmentation.pam_genome_defragmentation_hit_finder(
            options=options, basis_grouped=grouped_3_dict, plausability_cutoff=options.pam_threshold, support_models_name="grp3_support_models.pkl"
        )

    # Print statistics on selection to the terminal
    logger.debug("Extended reference sequence datasets")
    log_all_mcl_cluster_statistics(
        processed_reference_dict, extended_grouped, grouped_3_dict, grouped_4_dict
    )

    merged_dict = myUtil.merge_grouped_refseq_dicts_simple(
        grouped_3_dict, grouped_4_dict
    )

    merged_dict = myUtil.merge_grouped_refseq_dicts_simple(
        processed_reference_dict, merged_dict
    )


    myUtil.save_cache(options, "mcl_PAM_plausible_hits.pkl", grouped_4_dict)
    myUtil.save_cache(options, "mcl_PAM_csb_merged_hits.pkl", merged_dict)

    return merged_dict


def select_ref_seq_mcl_sequences(
    mcl_file: str, domain: str, reference_sequences: Set[str]
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

    def parse_mcl_clusters(path):
        # Read MCL file and return a list of clusters (each cluster is a list of protein IDs)
        with open(path, "r") as f:
            return [line.strip().split() for line in f if line.strip()]

    clusters = parse_mcl_clusters(mcl_file)

    all_protein_ids_with_refseqs = set()

    for cluster in clusters:
        cluster_set = set(cluster)

        # If the cluster contains at least one reference protein, add all proteins from that cluster
        if cluster_set & reference_sequences:
            all_protein_ids_with_refseqs.update(cluster_set)

    # print(f"[{domain}] Clusters with reference sequences found: {len(all_protein_ids_with_refseqs)} proteins total")

    return all_protein_ids_with_refseqs


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
        for i in range(0, len(iterable), size):
            yield list(iterable)[i : i + size]

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
