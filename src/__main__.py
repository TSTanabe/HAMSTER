#!/usr/bin/python

import os
import sys

from src.cli import cli

from src.fasta_preparation import fasta_preparation_stage

from src.core import logging, queue
from src.core import project
from src.core import myUtil

from src.search import initial_search

from src.csb import csb_prediction_stage
from src.selection_clustering import clustering_selection_stage
from src.selection_seed import basis_selection_stage, csb_proteins_selection

from src.dataset_testing import validation, alignment, Reports_printing, Reports, Reports_plotting

from src.selection_defragmentation import (
    presence_absence_defragmentation_stage,
    protein_mcl,
)

# Initilize the logger and location
from src.core.logging import get_logger

logger = get_logger(__name__)

if getattr(sys, "frozen", False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__))
    )


def model_alignment(options) -> None:
    """
    Runs sequence alignments for all clusters.

    Args:
        options (Options): Pipeline options

    Input Example:
        options.fasta_output_directory: str = "./output/aln"

    Output:
        - Alignments written to disk (e.g., .fasta_aln files)
    """

    alignment.initial_alignments(options, options.fasta_output_directory)

    return


def cross_validation(options) -> None:
    """
    Performs cross-validation on HMMs and protein sets.

    Args:
        options (Options): Pipeline options

    Output:
        - HMM files, validation reports.
        - Updates options.sequence_faa_file and targeted sequence sets.

    Example Output:
        options.sequence_faa_file: str = "./crossval/sequences.faa"
        options.targeted_sequence_faa_file_dict: dict[str, str] = {...}
    """

    logger.info("Generating hidden Markov models")

    validation.create_hmms_from_msas(
        options.fasta_output_directory,
        options.Hidden_markov_model_directory,
        "fasta_aln",
        "hmm",
        options.cores,
    )

    myUtil.move_HMMs(
        options.fasta_output_directory, options.Hidden_markov_model_directory, "hmm"
    )  # move the hmms to the Hidden markov model folder

    csb_proteins_selection.fetch_domains_superfamily_to_fasta(
        options, options.cross_validation_directory
    )  # Target file for each HMM excluding seqs already below threshold

    csb_proteins_selection.fetch_all_proteins(
        options.database_directory,
        options.cross_validation_directory + "/sequences.faa",
    )  # Backup target file if something fails

    logger.info("Generating validation sequence sets")

    options.sequence_faa_file = (
        options.cross_validation_directory + "/sequences.faa"
    )  # File with all sequences to be searched

    options.targeted_sequence_faa_file_dict = validation.get_target_sets(
        options.cross_validation_directory
    )

    validation.initial_self_recognition_validation(options)

    if options.cross_validation_activated:
        validation.parallel_cross_validation(options)

    return


def report_cv_performance(options) -> None:
    """
    Generates and saves all validation/cutoff performance reports.

    Args:
        options (Options): Pipeline options

    Output:
        - Performance plots, cutoff files, summary reports.
        - Various reports saved to options.Hidden_markov_model_directory

    Example Output:
        Files: cv_cutoff_performance.txt, cv_strict_cutoffs.txt, plots, etc.
    """
    logger.info(
        f"Saving the cutoffs and performance reports from initial calculation to {options.Hidden_markov_model_directory}"
    )
    mcl_clustering_results_dict = myUtil.load_cache(
        options, "linclust_clustering_results.pkl"
    )

    mcl_clustering_results_dict = protein_mcl.validate_mcl_cluster_paths(
        mcl_clustering_results_dict, options.result_files_directory
    )  # Check for path existence

    myUtil.save_cache(
        options,
        "mcl_clustering_results.pkl",
        mcl_clustering_results_dict,
        overwrite=True,
    )

    Reports_printing.process_initial_validations(
        options,
        options.result_files_directory,
        options.fasta_alignment_directory,
        options.database_directory,
    )

    # External R scripts
    try:
        logger.info("Plotting confusion matrices per genomic context")
        Reports_plotting.process_initial_plotting(options)
    except Exception as e:
        logger.error("An error occurred during the R plotting: %s", str(e))
    # cross validation
    logger.info("Starting cross-validation")
    logger.info(
        f"Saving the cutoffs and performance reports from the cross-validation to {options.Hidden_markov_model_directory}"
    )

    options.reports = Reports.parse_matrices_to_report(
        options.cross_validation_directory, "_cv_matrices.txt"
    )
    Reports.create_performance_file(
        options,
        options.reports,
        options.Hidden_markov_model_directory,
        "/cv_cutoff_performance.txt",
    )
    Reports.concatenate_cv_cutoff_files(
        options.cross_validation_directory,
        "_cv_thresholds.txt",
        options.Hidden_markov_model_directory + "/cv_strict_cutoffs.txt",
    )


def main(args: list = None) -> None:
    """
    Main entrypoint for the HAMSTER pipeline.

    Args:
        args (list, optional): sys.argv[1:] or custom list.

    Steps:
        1. Prepare directories
        2. Prepare/queue genome files
        3. Run search
        4. CSB prediction
        5. Basis FASTA
        6. PAM defragmentation
        7. Linclust/MCL
        8. Alignment
        9. Cross-validation
       10. Reporting

    Output:
        None
    """
    options = cli.parse_arguments(args)

    # 1
    project.prepare_result_space(options)
    log_file = os.path.join(options.result_files_directory, "execution_logfile.txt")
    logging.setup_logging(getattr(options, "verbose", 0), log_file)

    # 2
    if options.stage <= 2 and options.end >= 2:
        myUtil.print_header(
            "\n 2. Prokaryotic gene recognition and translation via prodigal"
        )
        fasta_preparation_stage.fasta_preparation(options)

    # 3
    if options.stage <= 3 and options.end >= 3:
        myUtil.print_header("\n 3. Searching for homologoues sequences")
        queue.queue_protein_annotation_inputs(options)
        initial_search.initial_search(options)

    # 4
    # csb naming
    if options.stage <= 4 and options.end >= 4:
        myUtil.print_header("\n 4. Searching for collinear syntenic blocks")
        csb_prediction_stage.csb_finder(options)

    # 5
    if options.stage <= 5 and options.end >= 5:
        myUtil.print_header(
            "\n 5. Selecting homologs to reference seqs by similarity and synteny"
        )
        basis_selection_stage.basis_sequence_fasta(options)  # grp0 dataset

    # 6
    if options.stage <= 6 and options.end >= 6:
        myUtil.print_header(
            "\n 6. Sequence similarity clustering for protein sequences"
        )
        presence_absence_defragmentation_stage.pam_defragmentation_stage(
            options
        )  # grp1 dataset
        #
        # include seqs with 90% ident or same synteny or same presence absence

    # 7
    if options.stage <= 7 and options.end >= 7:
        myUtil.print_header(
            "\n 7. Selecting protein similarity clusters with reference sequences"
        )
        clustering_selection_stage.mcl_family_clustering_sequences(options)
        mcl_extended_grouped_grp2 = clustering_selection_stage.mcl_select_grp2_clusters(
            options
        )
        mcl_extended_grouped_grp3 = clustering_selection_stage.mcl_select_grp3_clusters(
            options, mcl_extended_grouped_grp2
        )
    # 8
    if options.stage <= 8 and options.end >= 8:
        myUtil.print_header("\n 8. Aligment of sequence datasets")
        model_alignment(options)

    # 9
    # make cross validation files
    # Validation.CV(options.fasta_output_directory,options.cores)
    if options.stage <= 9 and options.end >= 9:
        myUtil.print_header("\n 9. Performing validation procedure")
        cross_validation(options)

    # 10
    # mach eine weitere HMMsearch mit dem vollen model auf alle seqs und prüfe die treffer
    if options.stage <= 10 and options.end >= 10:
        myUtil.print_header("\n 10. Writing dataset performance reports")
        report_cv_performance(options)


if __name__ == "__main__":
    args = sys.argv[1:]

    main(args)
