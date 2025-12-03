#!/usr/bin/python
import sqlite3
from collections import defaultdict
from typing import Any, Dict, Set

from src.core import myUtil
from src.selection_defragmentation import (
    seq_clustering,
    protein_mcl,
    pam_defragmentation,
)
from src.selection_seed import (
    csb_proteins_selection,
    singleton_finder,
    singleton_finder2,
)

from src.core.logging import get_logger

logger = get_logger(__name__)


##########################################################################################################
######################### Extend grouped reference proteins with similar csb #############################
##########################################################################################################


def extend_merged_grouped_by_csb_similarity(
    options: Any, grouped: Dict[str, Set[str]]
) -> Dict[str, Set[str]]:
    """
    Main routine: Extends grouped protein sets by including proteins from highly similar CSB patterns.

    Args:
        options: Options/config object, must have csb_output_file, jaccard, sqlite_chunks, etc.
        grouped (dict): {domain: set(proteinIDs)} (starting protein sets).

    Returns:
        dict: Updated grouped dict with new proteins added.

    Example:
        {'ABC': {'p1', 'p2'}} → {'ABC': {'p1', 'p2', 'p3', 'p4'}}
    """
    # Main routine for proteins from highly similar csb to the basic csb

    # Find keywords that are in jaccard distance to the csb of reference sequences in grouped
    protein_to_new_keywords_dict = select_similar_csb_patterns_per_protein(
        options, grouped, options.jaccard
    )

    # Integrate proteins with the new keywords into merged_grouped dataset
    extended_grouped = integrate_csb_variants_into_merged_grouped(
        options, grouped, protein_to_new_keywords_dict, options.sqlite_chunks
    )

    return extended_grouped


def select_similar_csb_patterns_per_protein(
    options: Any, merged_grouped: Dict[str, Set[str]], jaccard_threshold: float = 0.7
) -> Dict[str, Set[str]]:
    """
    For each domain, finds all CSB patterns with high Jaccard similarity to those associated with the proteins in that domain.

    Args:
        options: Options/config object.
        merged_grouped (dict): {domain: set(proteinIDs)}
        jaccard_threshold (float): Similarity threshold.

    Returns:
        dict: {domain: set(similar CSB keywords)}

    Example:
        {'ABC': {'p1', 'p2'}} → {'ABC': {'csb1', 'csb2'}}
    """

    # Stores the similar CSB keywords for each domain
    jaccard_included_patterns = {}  # { domain : set(csb_keywords) }

    # Load all CSB patterns from csb_output_file
    csb_dictionary = csb_proteins_selection.parse_csb_file_to_dict(
        options.csb_output_file
    )

    # Process each domain
    for domain, protein_ids in merged_grouped.items():
        logger.debug(f"Processing domain {domain}")

        # Fetch all keywords associated with proteins of this domain
        all_keywords = fetch_keywords_for_proteins(
            options.database_directory, protein_ids
        )

        if not all_keywords:
            logger.warning(f"No keywords found for domain {domain}")
            continue

        # Build union of patterns from these keywords
        domain_pattern_union = set().union(
            *[csb_dictionary[k] for k in all_keywords if k in csb_dictionary]
        )

        logger.debug(
            f"Domain {domain}: {len(domain_pattern_union)} unique CSB encode reference sequences."
        )

        # Compare against all CSB patterns and collect similar ones
        similar_csb_keywords = set()

        for csb_key, csb_pattern in csb_dictionary.items():
            intersection = len(domain_pattern_union & csb_pattern)
            union = len(domain_pattern_union | csb_pattern)
            similarity = intersection / union if union > 0 else 0.0

            if similarity >= jaccard_threshold:
                similar_csb_keywords.add(csb_key)

        logger.debug(
            f"Domain {domain}: {len(similar_csb_keywords)} similar CSB patterns (Jaccard >= {jaccard_threshold})"
        )

        # Store the result for this domain
        jaccard_included_patterns[domain] = similar_csb_keywords

    # Return dictionary: domain → set of similar CSB keywords
    return jaccard_included_patterns


def fetch_keywords_for_proteins(
    database_path: str, protein_ids: Set[str], chunk_size: int = 999
) -> Set[str]:
    """
    Fetches all CSB keywords (cluster identifiers) for a set of protein IDs.

    Args:
        database_path (str): Path to SQLite database.
        protein_ids (set): Protein IDs to fetch keywords for.
        chunk_size (int): DB chunk size.

    Returns:
        set: All unique keywords.

    Example:
        {'p1','p2'} → {'csb1', 'csb2'}
    """

    protein_to_keywords = defaultdict(set)

    if not protein_ids:
        logger.warning("No proteinIDs provided to fetch_keywords_for_proteins.")
        return set()

    conn = sqlite3.connect(database_path)
    cur = conn.cursor()

    protein_ids = list(protein_ids)
    total = len(protein_ids)
    # print(f"[INFO] Fetching keywords for {total} proteins")

    for start in range(0, total, chunk_size):
        end = start + chunk_size
        chunk = protein_ids[start:end]

        protein_placeholders = ",".join(["?"] * len(chunk))
        query = f"""
        SELECT Proteins.proteinID, Keywords.keyword
        FROM Proteins
        LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
        WHERE Proteins.proteinID IN ({protein_placeholders})
        """

        cur.execute(query, chunk)
        rows = cur.fetchall()

        for protein_id, keyword in rows:
            if keyword:
                protein_to_keywords[protein_id].add(keyword)

        # print(f"[INFO] Processed proteins {start+1} to {min(end, total)} of {total}")

    conn.close()

    logger.debug(f"Retrieved keywords for {len(protein_to_keywords)} proteins.")

    # Union of all keywords associated with these reference sequences
    all_keywords = set().union(*protein_to_keywords.values())

    return all_keywords


def integrate_csb_variants_into_merged_grouped(
    options: Any,
    merged_grouped: Dict[str, Set[str]],
    domain_to_new_keywords_dict: Dict[str, Set[str]],
    chunk_size: int = 999,
) -> Dict[str, Set[str]]:
    """
    For each domain, integrates all proteins whose CSB keyword matches those in domain_to_new_keywords_dict.

    Args:
        options: Options/config object (needs database_directory)
        merged_grouped (dict): {domain: set(proteinIDs)} to be expanded
        domain_to_new_keywords_dict (dict): {domain: set(keywords) to add}
        chunk_size (int): DB chunk size.

    Returns:
        dict: Updated merged_grouped.

    Example:
        ({'ABC': {...}}, {'ABC': {'csb5'}}) → {'ABC': {..., 'p55', 'p56'}}
    """

    logger.debug("Integration of added CSB proteins to grouped dataset")

    conn = sqlite3.connect(options.database_directory)
    cur = conn.cursor()

    total_proteins_added = 0

    for domain, new_keywords in domain_to_new_keywords_dict.items():
        if not new_keywords:
            logger.debug(f"Domain {domain}: No new keywords to integrate.")
            continue

        logger.debug(
            f"Domain {domain}: Integrating sequences from {len(new_keywords)} new keywords."
        )

        new_keywords_list = list(new_keywords)
        domain_proteins_before = len(merged_grouped.get(domain, set()))
        proteins_added_this_domain = 0

        for start in range(0, len(new_keywords_list), chunk_size):
            end = start + chunk_size
            chunk = new_keywords_list[start:end]

            placeholders = ",".join(["?"] * len(chunk))
            query = f"""
            SELECT Proteins.proteinID
            FROM Proteins
            LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
            LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
            WHERE Domains.domain = ? AND Keywords.keyword IN ({placeholders})
            """

            params = [domain] + chunk
            cur.execute(query, params)
            rows = cur.fetchall()

            for (protein_id,) in rows:
                if protein_id not in merged_grouped[domain]:
                    merged_grouped[domain].add(protein_id)
                    proteins_added_this_domain += 1

        domain_proteins_after = len(merged_grouped[domain])
        logger.info(
            f"Domain {domain}: Added {proteins_added_this_domain} new proteins with matching synteny to reference. New total: {len(merged_grouped[domain])}"
        )

        total_proteins_added += proteins_added_this_domain

    conn.close()

    logger.info(
        f"Integration completed: {total_proteins_added} proteins added across all domains."
    )

    return merged_grouped


def pam_defragmentation_stage(options) -> object | None:
    """
    Find additional plausible hits based on presence absence patterns. This should include hits
    from fragmented assemblies or split csb

    Prepare the presence absence matrix for grp0 and train logistic regression on the matrix
    With the trained matrix exclude each column and predict presence.
    For predicted presences select from the genome the best hit

    Also add csb that are below jaccard distance threshold from the grp0 csb

    Output: are the grp1 fasta files

    Finds additional plausible hits based on presence/absence patterns (grp1 dataset).

    Args:
        options (Options): Pipeline options

    Output:
        - Updates options.grouped for further analysis.
        - Writes grp1 FASTA files.
    """
    basis_grouped = (
        options.grouped
        if hasattr(options, "grouped")
        else myUtil.load_cache(options, "basis_merged_grouped.pkl")
    )
    basis_score_limit_dict = (
        options.score_limit_dict
        if hasattr(options, "score_limit_dict")
        else myUtil.load_cache(options, "basis_merged_score.pkl")
    )

    # Load precomputed grp1 results if available
    grp1_merged_dict = myUtil.load_cache(options, "grp1_merged_grouped.pkl")
    if grp1_merged_dict:
        options.grouped = grp1_merged_dict
        return grp1_merged_dict

    # If precomputed grp1 is not available then compute the grp1

    # Adds protein sequences from csb that are below jaccard distance threshold distance to grp0 csb
    added_similar_csb_proteins = extend_merged_grouped_by_csb_similarity(
        options, basis_grouped
    )

    # Adds potential hits by presence absence matrix
    added_pam_propability_proteins = (
        pam_defragmentation.pam_genome_defragmentation_hit_finder(
            options, basis_grouped, basis_score_limit_dict
        )
    )

    # Merge the added proteins, for same key in both sets sum up the sets
    merged_grouped = csb_proteins_selection.merge_grouped_protein_ids(
        added_similar_csb_proteins, added_pam_propability_proteins
    )

    # Calculate the score limits for the reference sequences
    score_limit_dict = csb_proteins_selection.generate_score_limit_dict_from_grouped(
        options.database_directory, merged_grouped
    )

    # Write fasta files with the reference sequences and similar sequences within the score cutoff range of the reference seqs for the linclustering
    csb_proteins_selection.fetch_protein_family_sequences(
        options, options.fasta_initial_hit_directory, score_limit_dict, merged_grouped
    )

    # Cluster sequences at 90 % identity and 70 % coverage to select highly similar proteins without context
    linclust_mcl_format_output_files_dict = seq_clustering.run_mmseqs_linclust_lowlevel(
        options.fasta_initial_hit_directory, "0.9", "0.7", options.cores
    )  # seq identity=> float 0.9 und min aln length => float 0.7

    # Add the clustered hits to the reference sequence sets
    # _linclust_mcl_format.txt select from these files in fasta_initial_hit_directory
    mcl_extended_grouped, mcl_cutoffs = protein_mcl.select_hits_by_csb_mcl(
        options, linclust_mcl_format_output_files_dict, merged_grouped, 0.0, 0.0001
    )  # low cutoffs for closely related protein clusters

    # Save computed grp1 datasets
    myUtil.save_cache(options, "grp1_merged_grouped.pkl", mcl_extended_grouped)
    myUtil.save_cache(options, "grp1_merged_score_limits.pkl", score_limit_dict)

    # Print the grp0 csb and singletons to fasta
    csb_proteins_selection.fetch_training_data_to_fasta(options, merged_grouped, "grp1")

    # Result dictionary is stores in options.grouped, overwriting the grp0 with grp1 key_domain pairs
    options.grouped = merged_grouped

    return


"""
grp1 has the 90 % identity sequences to basic set
proteins with same synteny but not collinearity
similar presence absence pattern
"""
