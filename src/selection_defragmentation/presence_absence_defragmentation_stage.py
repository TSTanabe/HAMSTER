#!/usr/bin/python
import sqlite3
from collections import defaultdict
from typing import Any, Dict, Set

from src.core import myUtil
from src.selection_defragmentation import (
    seq_clustering,
    protein_mcl,
    pam_defragmentation
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

    # Find keywords that are in jaccard distance to the csb of reference sequences in grouped. Load if possible
    protein_to_new_keywords_dict = myUtil.load_cache(options, "grp1_protein_to_key.pkl")
    if not protein_to_new_keywords_dict:
        protein_to_new_keywords_dict = select_similar_csb_patterns_per_protein( # can be loaded from here as a single thread version
            options, grouped, options.jaccard
        )
        myUtil.save_cache(options, "grp1_protein_to_key.pkl", protein_to_new_keywords_dict)

    # Integrate proteins with the new keywords into merged_grouped dataset. Load if possible
    extended_grouped = myUtil.load_cache(options, "grp1_extended_grouped.pkl")
    if not extended_grouped:
        extended_grouped = integrate_csb_variants_into_merged_grouped(
            options, grouped, protein_to_new_keywords_dict, options.sqlite_chunks
        )
        myUtil.save_cache(options, "grp1_extended_grouped.pkl", extended_grouped)

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

        domain_pattern_union_len = len(domain_pattern_union)

        for csb_key, csb_pattern in csb_dictionary.items():

            csb_pattern_len = len(csb_pattern)

            # Upper bound pruning (sehr effektiv)
            max_possible = min(domain_pattern_union_len, csb_pattern_len) / max(domain_pattern_union_len, csb_pattern_len) if max(domain_pattern_union_len, csb_pattern_len) > 0 else 0.0
            if max_possible < jaccard_threshold:
                continue

            # Intersection ohne temporäres Set
            if domain_pattern_union_len < csb_pattern_len:
                intersection = sum(1 for x in domain_pattern_union if x in csb_pattern)
            else:
                intersection = sum(1 for x in csb_pattern if x in domain_pattern_union)

            if intersection == 0:
                if jaccard_threshold > 0.0:
                    continue

            union = domain_pattern_union_len + csb_pattern_len - intersection
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
    database_path: str, protein_ids: Set[str], chunk_size: int = 50000
) -> Set[str]:
    """
    Fetches all CSB keywords for a set of protein IDs.

    Speed-up:
      - Uses a TEMP table (tmp_protein_ids) instead of IN(...) chunking
      - Single join query with DISTINCT keywords
      - Avoids building protein->keywords dict when only union is needed

    Notes:
      - chunk_size here controls batch size for executemany inserts into TEMP, not SQL placeholders.
    """
    if not protein_ids:
        logger.warning("No proteinIDs provided to fetch_keywords_for_proteins.")
        return set()

    all_keywords: Set[str] = set()

    with sqlite3.connect(database_path, timeout=120.0) as con:
        #logger.info(f"Fetching {len(all_keywords)} keywords for {len(protein_ids)} proteins.")
        cur = con.cursor()

        # Pragmas: TEMP in memory; speed tuning
        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")      # ~256 MiB
        cur.execute("PRAGMA mmap_size=2147483648;")    # 2 GiB
        cur.execute("PRAGMA automatic_index=ON;")

        # TEMP table once
        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_protein_ids (proteinID TEXT PRIMARY KEY);"
        )

        # Fill TEMP table in batches (fast + low peak memory)
        ids = list(protein_ids)
        cur.execute("DELETE FROM tmp_protein_ids;")

        for start in range(0, len(ids), chunk_size):
            batch = ids[start : start + chunk_size]
            cur.executemany(
                "INSERT OR IGNORE INTO tmp_protein_ids(proteinID) VALUES (?);",
                ((pid,) for pid in batch),
            )

        # Protect persistent DB now (TEMP already filled)
        cur.execute("PRAGMA query_only=TRUE;")

        # One query: get unique keywords for these proteins
        cur.execute(
            """
            SELECT DISTINCT k.keyword
            FROM tmp_protein_ids t
            JOIN Proteins p ON p.proteinID = t.proteinID
            LEFT JOIN Keywords k ON k.clusterID = p.clusterID
            WHERE k.keyword IS NOT NULL
            """
        )

        all_keywords = {row[0] for row in cur.fetchall()}

    logger.debug(f"Retrieved {len(all_keywords)} keywords for {len(protein_ids)} proteins.")
    return all_keywords


def integrate_csb_variants_into_merged_grouped(
    options: Any,
    merged_grouped: Dict[str, Set[str]],
    domain_to_new_keywords_dict: Dict[str, Set[str]],
    chunk_size: int = 50000,  # now used for TEMP insert batching, not SQL IN()
) -> Dict[str, Set[str]]:
    """
    For each domain, integrates all proteins whose CSB keyword matches those in
    domain_to_new_keywords_dict, using TEMP tables for speed (no IN-chunking).

    chunk_size controls batching for executemany inserts into TEMP tables.
    """
    logger.info("Integration of added CSB proteins to grouped dataset")

    if not domain_to_new_keywords_dict:
        return merged_grouped

    with sqlite3.connect(options.database_directory, timeout=120.0) as con:
        cur = con.cursor()

        # Pragmas: TEMP in memory; read workload tuning
        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")      # ~256 MiB
        cur.execute("PRAGMA mmap_size=2147483648;")    # 2 GiB
        cur.execute("PRAGMA automatic_index=ON;")

        # TEMP table once
        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_keywords (keyword TEXT PRIMARY KEY);"
        )

        total_proteins_added = 0

        for domain, new_keywords in domain_to_new_keywords_dict.items():
            if not new_keywords:
                logger.debug(f"Domain {domain}: No new keywords to integrate.")
                continue

            # Ensure target set exists
            merged_grouped.setdefault(domain, set())

            logger.debug(
                f"Domain {domain}: Integrating sequences from {len(new_keywords)} new keywords."
            )

            before = len(merged_grouped[domain])

            # Fill tmp_keywords for this domain (batched)
            cur.execute("DELETE FROM tmp_keywords;")
            kw_list = [k for k in new_keywords if k]
            for start in range(0, len(kw_list), chunk_size):
                batch = kw_list[start : start + chunk_size]
                cur.executemany(
                    "INSERT OR IGNORE INTO tmp_keywords(keyword) VALUES (?);",
                    ((k,) for k in batch),
                )

            # Protect persistent DB after TEMP is filled (optional but safe)
            cur.execute("PRAGMA query_only=TRUE;")

            # One query: fetch proteinIDs that (a) belong to clusters annotated with these keywords
            # and (b) have the given domain in Domains table.
            #
            # Note: JOIN Domains is required to enforce domain membership.
            cur.execute(
                """
                SELECT DISTINCT p.proteinID
                FROM tmp_keywords tk
                JOIN Keywords k
                  ON k.keyword = tk.keyword
                JOIN Proteins p
                  ON p.clusterID = k.clusterID
                JOIN Domains d
                  ON d.proteinID = p.proteinID
                 AND d.domain    = ?
                """,
                (domain,),
            )

            added = 0
            for (protein_id,) in cur:
                if protein_id not in merged_grouped[domain]:
                    merged_grouped[domain].add(protein_id)
                    added += 1

            # allow next TEMP refill (some builds restrict TEMP writes under query_only)
            cur.execute("PRAGMA query_only=FALSE;")

            after = len(merged_grouped[domain])
            logger.info(
                f"Domain {domain}: Added {added} new proteins with matching synteny to reference. "
                f"New total: {after} (was {before})"
            )
            total_proteins_added += added

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
