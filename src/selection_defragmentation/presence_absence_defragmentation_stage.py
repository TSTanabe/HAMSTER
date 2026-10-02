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
    csb_type_statistic,
)

from src.core.logging import get_logger

logger = get_logger(__name__)


##########################################################################################################
######################### Extend grouped reference proteins with similar csb #############################
##########################################################################################################

def pam_defragmentation_stage(config) -> object | None:
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
        config (Options): Pipeline options

    Output:
        - Updates options.grouped for further analysis.
        - Writes grp1 FASTA files.
    """
    # Load precomputed grp1 results if available
    grp1_merged_dict = myUtil.load_cache(config, "grp1_merged_grouped.pkl")
    if grp1_merged_dict:
        config.grouped = grp1_merged_dict
        return grp1_merged_dict

    basis_grouped = (
        config.grouped
        if hasattr(config, "grouped")
        else myUtil.load_cache(config, "basis_merged_grouped.pkl")
    )
    basis_score_limits = (
        config.score_limit_dict
        if hasattr(config, "score_limit_dict")
        else myUtil.load_cache(config, "basis_merged_score.pkl")
    )

    # Adds potential hits by presence absence matrix
    predictor_probability_proteins = predictor.predictor_training_calibration_application(
            config=config,
            basis_seed_sequences=basis_grouped,
            basis_score_limit=basis_score_limits,
            probability_cutoff=config.pam_threshold,
            support_models_name="grp2_predictor_models.pkl",
        )

    # Adds protein sequences from csb that are below jaccard distance threshold distance to grp0 csb
    syntenic_proteins = identical_synteny.extend_merged_grouped_by_csb_similarity(
        config, basis_grouped
    )

    # Merge the added proteins, for same key in both sets sum up the sets
    merged_grouped = csb_proteins_selection.merge_protein_sets(
        syntenic_proteins, predictor_probability_proteins
    )

    # Calculate the score limits for the reference sequences
    score_limit_dict = csb_type_statistic.generate_score_limits_from_seed_dict(
        config.database_directory, merged_grouped
    )

    ## Clustering at 90 % identitity
    # Write fasta files with the reference sequences and similar sequences within the score cutoff range of the reference seqs for the linclustering
    csb_proteins_selection.fetch_protein_family_sequences(
        config=config, directory=config.fasta_initial_hit_directory, score_limit_dict=score_limit_dict, domain_to_proteinID=merged_grouped
    )

    # Cluster sequences at 90 % identity and 70 % coverage to select highly similar proteins without context
    linclust_mcl_format_output_files_dict = seq_clustering.run_mmseqs_linclust_lowlevel(
        directory=config.fasta_initial_hit_directory, min_seq_id= 0.9, min_aln_len= 0.7, cores=config.cores
    )  # seq identity=> float 0.9 und min aln length => float 0.7

    # Add the clustered hits to the reference sequence sets
    # _linclust_mcl_format.txt select from these files in fasta_initial_hit_directory
    mcl_extended_grouped, mcl_cutoffs = protein_mcl.select_hits_by_csb_mcl(
        config=config, mcl_output_dict=linclust_mcl_format_output_files_dict, reference_dict=merged_grouped, density_threshold=0.0, reference_threshold=0.0001
    )  # low cutoffs for closely related protein clusters

    ## Storage
    # Save computed grp1 datasets
    myUtil.save_cache(config, "grp1_merged_grouped.pkl", mcl_extended_grouped)
    myUtil.save_cache(config, "grp1_merged_score_limits.pkl", score_limit_dict)

    # Print the grp0 csb and singletons to fasta
    csb_proteins_selection.fetch_training_data_to_fasta(config, merged_grouped, "ds2")

    # Result dictionary is stores in options.grouped, overwriting the grp0 with grp1 key_domain pairs
    config.grouped = merged_grouped

    return


"""
grp1 has the 90 % identity sequences to basic set
proteins with same synteny but not collinearity
similar presence absence pattern with 0.8 plausability score
"""
