#!/usr/bin/python
import sqlite3
import sys
from typing import Any, Dict, Set

from src.core import myUtil
from src.selection_defragmentation import (
    seq_clustering,
    protein_mcl,
    predictor,
    identical_synteny,
)
from src.selection_seed import (
    csb_proteins_selection,
    fetch_seed_proteins,
)

from src.core.logging import get_logger

logger = get_logger(__name__)


##########################################################################################################
######################### Extend grouped reference proteins with similar csb #############################
##########################################################################################################


def generate_score_limit_dict_from_grouped(
    database_path: str,
    grouped_dict: Dict[str, Set[str]],
    default_lower: int = 100,
    default_upper: int = 2000,
    chunk_size: int = 999,
) -> Dict[str, Dict[str, float]]:
    """
    Generate a score limit dictionary for each domain from grouped proteinIDs.

    Args:
        database_path (str): SQLite DB path.
        grouped_dict (dict): {domain: set(proteinIDs)}
        default_lower (int): Fallback if no data.
        default_upper (int): Fallback if no data.
        chunk_size (int): Max number of SQL params.

    Returns:
        dict: {domain: {'lower_limit': X, 'upper_limit': Y}}
    """
    result = {}

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        for domain, protein_ids in grouped_dict.items():
            if not protein_ids:
                result[domain] = {
                    "lower_limit": default_lower,
                    "upper_limit": default_upper,
                }
                continue

            min_score = float("inf")
            max_score = float("-inf")
            protein_id_list = list(protein_ids)

            for i in range(0, len(protein_id_list), chunk_size):
                chunk = protein_id_list[i : i + chunk_size]
                placeholders = ",".join("?" for _ in chunk)
                query = f"""
                    SELECT MIN(score), MAX(score)
                    FROM Domains
                    WHERE domain = ? AND proteinID IN ({placeholders});
                """
                cur.execute(query, (domain, *chunk))
                row = cur.fetchone()
                if row:
                    chunk_min, chunk_max = row
                    if chunk_min is not None:
                        min_score = min(min_score, chunk_min)
                    if chunk_max is not None:
                        max_score = max(max_score, chunk_max)

            if min_score == float("inf") or max_score == float("-inf"):
                # No valid score data found
                result[domain] = {
                    "lower_limit": default_lower,
                    "upper_limit": default_upper,
                }
            else:
                result[domain] = {"lower_limit": min_score, "upper_limit": max_score}

    return result


def pam_defragmentation_stage(options) -> object | None:
    """
    Find additional plausible hits based on presence absence patterns. This should include hits
    from fragmented assemblies or split csb

    Prepare the presence absence matrix for grp0 and plausibility model on selected matrices
    Then iterate all genomes and for each protein define if presence is expected and add if above expectation
    threshold
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
    # Load precomputed grp1 results if available
    grp1_merged_dict = myUtil.load_cache(options, "grp1_merged_grouped.pkl")
    if grp1_merged_dict:
        options.grouped = grp1_merged_dict
        return grp1_merged_dict

    basis_grouped = (
        options.grouped
        if hasattr(options, "grouped")
        else myUtil.load_cache(options, "basis_merged_grouped.pkl")
    )
    basis_score_limits = (
        options.score_limit_dict
        if hasattr(options, "score_limit_dict")
        else myUtil.load_cache(options, "basis_merged_score.pkl")
    )

    # Adds potential hits by presence absence matrix
    predictor_probability_proteins = predictor.predictor_training_calibration_application(
            config=options,
            basis_seed_sequences=basis_grouped,
            basis_score_limit=basis_score_limits,
            probability_cutoff=options.pam_threshold,
            support_models_name="grp2_predictor_models.pkl",
        )

    # Adds protein sequences from csb that are below jaccard distance threshold distance to grp0 csb
    syntenic_proteins = identical_synteny.extend_merged_grouped_by_csb_similarity(
        options, basis_grouped
    )

    # Merge the added proteins, for same key in both sets sum up the sets
    merged_grouped = csb_proteins_selection.merge_protein_sets(
        syntenic_proteins, predictor_probability_proteins
    )

    # Calculate the score limits for the reference sequences
    score_limit_dict = generate_score_limit_dict_from_grouped(
        options.database_directory, merged_grouped
    )

    ## Clustering at 90 % identitity
    # Write fasta files with the reference sequences and similar sequences within the score cutoff range of the reference seqs for the linclustering
    csb_proteins_selection.fetch_protein_family_sequences(
        config=options, directory=options.fasta_initial_hit_directory, score_limit_dict=score_limit_dict, domain_to_proteinID=merged_grouped
    )

    # Cluster sequences at 90 % identity and 70 % coverage to select highly similar proteins without context
    linclust_mcl_format_output_files_dict = seq_clustering.run_mmseqs_linclust_lowlevel(
        directory=options.fasta_initial_hit_directory, min_seq_id= 0.9, min_aln_len= 0.7, cores=options.cores
    )  # seq identity=> float 0.9 und min aln length => float 0.7

    # Add the clustered hits to the reference sequence sets
    # _linclust_mcl_format.txt select from these files in fasta_initial_hit_directory
    mcl_extended_grouped, mcl_cutoffs = protein_mcl.select_hits_by_csb_mcl(
        config=options, mcl_output_dict=linclust_mcl_format_output_files_dict, reference_dict=merged_grouped, density_threshold=0.0, reference_threshold=0.0001
    )  # low cutoffs for closely related protein clusters

    ## Storage
    # Save computed grp1 datasets
    myUtil.save_cache(options, "grp1_merged_grouped.pkl", mcl_extended_grouped)
    myUtil.save_cache(options, "grp1_merged_score_limits.pkl", score_limit_dict)

    # Print the grp0 csb and singletons to fasta
    csb_proteins_selection.fetch_training_data_to_fasta(options, merged_grouped, "ds2")

    # Result dictionary is stores in options.grouped, overwriting the grp0 with grp1 key_domain pairs
    options.grouped = merged_grouped

    return


"""
grp1 has the 90 % identity sequences to basic set
proteins with same synteny but not collinearity
similar presence absence pattern with 0.8 plausability score
"""
