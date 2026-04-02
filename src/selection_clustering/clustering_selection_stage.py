from src.core import myUtil
from src.selection_clustering import pam_mcl
from src.selection_defragmentation import seq_clustering, protein_mcl
from src.selection_seed import csb_proteins_selection

from src.core.logging import get_logger

logger = get_logger(__name__)


def mcl_family_clustering_sequences(options) -> None:
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
        options (Options): Pipeline options

    Output:
        - FASTA files for protein families.
        - Clustering result files for further analysis.
    """

    # Load grp1 datasets, that includes basis + proteins with similar csb and presence absence patterns
    grouped = (
        options.grouped
        if hasattr(options, "grouped")
        else myUtil.load_cache(options, "grp1_merged_grouped.pkl")
    )
    score_limit_dict = (
        options.score_limit_dict
        if hasattr(options, "score_limit_dict")
        else myUtil.load_cache(options, "grp1_merged_score_limits.pkl")
    )

    logger.info("Prepare protein sequence identity clustering")
    csb_proteins_selection.fetch_protein_family_sequences(
        options, options.phylogeny_directory, score_limit_dict, grouped
    )

    # Cluster sequences with linclust at 40 % identitiy
    linclust_mcl_format_output_files_dict = seq_clustering.run_mmseqs_linclust_lowlevel(
        options.phylogeny_directory, options.mcl_min_seq_id, "0.7", options.cores
    )

    myUtil.save_cache(
        options,
        "linclust_clustering_results.pkl",
        linclust_mcl_format_output_files_dict,
    )

    ## This block is replaced by linclust.
    # Linclust is much faster than mcl clustering with similar performance
    # Markov chain clustering for grp1 fasta files
    # if not options.csb_mcl_clustering:
    #    return

    # print("\n[INFO] Calculating Markov Chain Clustering")
    # Clustering sequences with MCL algorithm
    # mcl_clustering_results_dict = Csb_mcl.csb_mcl_datasets(options,grouped) # markov chain clustering grouped training data

    # myUtil.save_cache(options, 'mcl_clustering_results.pkl', mcl_clustering_results_dict)

    options.mcl_clustering_results_dict = linclust_mcl_format_output_files_dict

    return


def mcl_select_grp2_clusters(options) -> dict:
    """
    Selects MCL clusters with sufficient fraction of reference sequences (grp2).

    Args:
        options (Options): Pipeline options

    Output:
        - grp2 FASTA files written to disk.
        - Returns mcl_extended_grouped dictionary.

    Returns:
        mcl_extended_grouped: dict[str, set[str]]
    """

    # Load grp1 reference sets
    grouped = (
        options.grouped
        if hasattr(options, "grouped")
        else myUtil.load_cache(options, "grp1_merged_grouped.pkl")
    )

    # Load and validate clustering results
    mcl_clustering_results_dict = myUtil.load_cache(
        options, "linclust_clustering_results.pkl"
    )
    mcl_clustering_results_dict = protein_mcl.validate_mcl_cluster_paths(
        mcl_clustering_results_dict, options.result_files_directory
    )
    myUtil.save_cache(
        options,
        "linclust_clustering_results.pkl",
        mcl_clustering_results_dict,
        overwrite=True,
    )

    logger.info(
        "Generating grp2: Selecting MCL clusters with sufficient reference fraction"
    )
    mcl_extended_grouped, mcl_cutoffs = protein_mcl.select_hits_by_csb_mcl(
        options,
        mcl_clustering_results_dict,
        grouped,
        options.mcl_density_thrs,
        options.mcl_reference_thrs,
    )

    myUtil.save_cache(options, "mcl_grp2_cluster_selection_cutoffs.pkl", mcl_cutoffs)
    myUtil.save_cache(options, "grp2_merged_grouped.pkl", mcl_extended_grouped)

    csb_proteins_selection.fetch_training_data_to_fasta(
        options, mcl_extended_grouped, "grp2"
    )

    return mcl_extended_grouped


def mcl_select_grp3_clusters(options, mcl_extended_grouped_grp2: dict) -> dict:
    """
    Extends grp2 by PAM model, produces grp3.

    Args:
        options (Options): Pipeline options
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

    grouped = (
        options.grouped
        if hasattr(options, "grouped")
        else myUtil.load_cache(options, "grp1_merged_grouped.pkl")
    )
    mcl_clustering_results_dict = myUtil.load_cache(
        options, "linclust_clustering_results.pkl"
    )
    mcl_clustering_results_dict = protein_mcl.validate_mcl_cluster_paths(
        mcl_clustering_results_dict, options.result_files_directory
    )
    myUtil.save_cache(
        options,
        "linclust_clustering_results.pkl",
        mcl_clustering_results_dict,
        overwrite=True,
    )

    # Extend references via PAM model
    regrouped = pam_mcl.select_hits_by_pam_csb_mcl(
        options, mcl_clustering_results_dict, grouped
    )
    myUtil.save_cache(options, "grp3_selection_ref_seqs.pkl", regrouped)

    pam_mcl_extended_grouped, mcl_cutoffs = protein_mcl.select_hits_by_csb_mcl(
        options,
        mcl_clustering_results_dict,
        regrouped,
        options.mcl_density_thrs,
        options.mcl_reference_thrs,
    )
    myUtil.save_cache(options, "mcl_grp3_cluster_selection_cutoffs.pkl", mcl_cutoffs)

    # Further extend via high coverage, low threshold
    mcl_extended_grouped_final, _ = protein_mcl.select_hits_by_csb_mcl(
        options, mcl_clustering_results_dict, pam_mcl_extended_grouped, 0.0, 0.0001
    )

    # Merge grp2 + extended grp3
    merged_grouped = csb_proteins_selection.merge_grouped_protein_ids(
        mcl_extended_grouped_final, pam_mcl_extended_grouped
    )

    myUtil.save_cache(options, "grp3_merged_grouped.pkl", merged_grouped)

    # Export fasta for grp3
    csb_proteins_selection.fetch_training_data_to_fasta(options, merged_grouped, "grp3")

    return merged_grouped
