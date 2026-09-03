from typing import Any

from src.core import myUtil
from src.selection_clustering import pam_mcl
from src.selection_defragmentation import seq_clustering, protein_mcl
from src.selection_seed import csb_proteins_selection

from src.core.logging import get_logger

logger = get_logger(__name__)


def _get_grp1_grouped(config: Any) -> dict[str, set[str]]:
    if hasattr(config, "grouped"):
        return config.grouped

    return myUtil.load_cache(
        config,
        "grp1_merged_grouped.pkl",
    )


def _load_linclust_results(config: Any) -> dict:
    clustering_results = myUtil.load_cache(
        config,
        "linclust_clustering_results.pkl",
    )

    clustering_results = protein_mcl.validate_mcl_cluster_paths(
        clustering_results,
        config.result_files_directory,
    )

    myUtil.save_cache(
        config,
        "linclust_clustering_results.pkl",
        clustering_results,
        overwrite=True,
    )

    return clustering_results


def mcl_family_clustering_sequences(config: Any) -> None:
    """
    This routine prepares the sequence clustering via linclust

    Formerly the MCL clustering algorithm was used, but this requires a lot of computational time
    due to the all vs all blast


    Prepare protein family file, including all hits above 25% identity.
    All vs. all diamond blast (Can this be accellerated like in proteinortho?)

    With the trained matrix exclude each column and predict presence.
    For predicted presences select from the genome the best hit

    Output: are the grp2 fasta files

    Prepares sequence clustering via linclust (replaces MCL).

    Args:
        config (Options): Pipeline options

    Output:
        - FASTA files for protein families.
        - Clustering result files for further analysis.
    """

    # Load grp1 datasets, that includes basis + proteins with similar csb and presence absence patterns
    grouped = _get_grp1_grouped(config)
    
    score_limit_dict = (
        config.score_limit_dict
        if hasattr(config, "score_limit_dict")
        else myUtil.load_cache(config, "grp1_merged_score_limits.pkl")
    )

    logger.info("Prepare protein sequence identity clustering")
    csb_proteins_selection.fetch_protein_family_sequences(
        config, config.phylogeny_directory, score_limit_dict, grouped
    )

    # Cluster sequences with linclust at 40 % identitiy
    linclust_mcl_format_output_files_dict = seq_clustering.run_mmseqs_linclust_lowlevel(
        directory=config.phylogeny_directory, min_seq_id=config.mcl_min_seq_id, min_aln_len=0.7, cores=config.cores
    )

    myUtil.save_cache(
        config,
        "linclust_clustering_results.pkl",
        linclust_mcl_format_output_files_dict,
    )

    config.mcl_clustering_results_dict = linclust_mcl_format_output_files_dict

    return


def mcl_select_grp2_clusters(config) -> dict:
    """
    Selects MCL clusters with sufficient fraction of reference sequences (grp2).

    Args:
        config (Options): Pipeline options

    Output:
        - grp2 FASTA files written to disk.
        - Returns mcl_extended_grouped dictionary.

    Returns:
        mcl_extended_grouped: dict[str, set[str]]
    """

    # Load grp1 reference sets
    grouped = _get_grp1_grouped(config)

    # Load and validate clustering results
    clustering_results = _load_linclust_results(config)

    logger.info(
        "Generating dataset 3: Selecting MCL clusters with sufficient number of reference sequences"
    )
    mcl_extended_grouped, mcl_cutoffs = protein_mcl.select_hits_by_csb_mcl(
        config,
        clustering_results,
        grouped,
        config.mcl_density_thrs,
        config.mcl_reference_thrs,
    )

    myUtil.save_cache(config, "mcl_grp2_cluster_selection_cutoffs.pkl", mcl_cutoffs)
    myUtil.save_cache(config, "grp2_merged_grouped.pkl", mcl_extended_grouped)

    csb_proteins_selection.fetch_training_data_to_fasta(
        config, mcl_extended_grouped, "ds3"
    )

    return mcl_extended_grouped


def mcl_select_grp3_clusters(config, grouped) -> dict:
    """
    Extends grp2 by PAM model, produces grp3.

    Args:
        config (Options): Pipeline options
        mcl_extended_grouped_grp2 (dict): Output from mcl_select_grp2_clusters

    Output:
        - grp3 FASTA files written to disk.
        - Returns merged_grouped (grp3) dictionary.

    Returns:
        merged_grouped: dict[str, set[str]]
    """
    logger.info(
        "Generating grp3: Extended MCL cluster selection by csb and presence plausibility"
    )

    clustering_results = _load_linclust_results(config)

    score_limit_dict = (
        config.score_limit_dict
        if hasattr(config, "score_limit_dict")
        else myUtil.load_cache(config, "grp1_merged_score_limits.pkl")
    )

    # Extend references via PAM model
    regrouped = pam_mcl.select_hits_by_pam_csb_mcl(
        config=config,
        clustering_results=clustering_results,
        basis_seed_sequences=grouped,
        basis_score_limit=score_limit_dict
    )
    myUtil.save_cache(config, "grp3_selection_ref_seqs.pkl", regrouped)

    pam_mcl_extended_grouped, mcl_cutoffs = protein_mcl.select_hits_by_csb_mcl(
        config,
        clustering_results,
        regrouped,
        config.mcl_density_thrs,
        config.mcl_reference_thrs,
    )
    myUtil.save_cache(config, "mcl_grp3_cluster_selection_cutoffs.pkl", mcl_cutoffs)

    # Further extend via high coverage, low threshold
    mcl_extended_grouped_final, _ = protein_mcl.select_hits_by_csb_mcl(
        config, clustering_results, pam_mcl_extended_grouped, 0.0, 0.0001
    )

    # Merge grp2 + extended grp3
    merged_grouped = csb_proteins_selection.merge_protein_sets(
        mcl_extended_grouped_final, pam_mcl_extended_grouped
    )

    myUtil.save_cache(config, "grp3_merged_grouped.pkl", merged_grouped)

    # Export fasta for grp3
    csb_proteins_selection.fetch_training_data_to_fasta(config, merged_grouped, "ds4")

    return merged_grouped
