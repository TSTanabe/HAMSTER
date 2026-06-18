#!/usr/bin/python

import csv
import os
import sys
import tarfile
from datetime import datetime

from src.core.logging import get_logger

logger = get_logger(__name__)

PROJECT_DIRS = {
    "fasta_initial_hit_directory": "Sequence_Clusters70",
    "fasta_output_directory": "Sequences",
    "fasta_alignment_directory": "Initial_validation",
    "Hidden_markov_model_directory": "Hidden_markov_models",
    "cross_validation_directory": "Cross_validation",
    "phylogeny_directory": "Sequence_Clusters40",
    "Csb_directory": "Collinear_syntenic_blocks",
}

def prepare_directory_structure(directory: str) -> None:
    """
    Sets up the `bin` and `src` directories within the specified directory.
    Unpacks any `.tar.gz` files in `bin` and deletes the archives.
    Creates an empty `patterns` file in `src` if it does not exist.

    Args:
        directory (str): Path to the directory to set up.
    """
    # Define the paths for `bin` and `src` directories
    bin_dir = os.path.join(directory, "bin")
    src_dir = os.path.join(directory, "src")

    # Ensure `bin` and `src` directories exist
    os.makedirs(bin_dir, exist_ok=True)
    os.makedirs(src_dir, exist_ok=True)

    # Process `.tar.xz` files in `bin` directory
    for file_name in os.listdir(bin_dir):
        if file_name.endswith(".tar.xz"):
            file_path = os.path.join(bin_dir, file_name)

            # Unpack the .tar.xz file
            with tarfile.open(file_path, "r:xz") as tar:
                tar.extractall(bin_dir)

            # Delete the original .tar.xz file
            os.remove(file_path)

    # Create an empty `patterns` file in `src` if it does not exist
    patterns_file = os.path.join(src_dir, "Patterns")
    if not os.path.exists(patterns_file):
        with open(patterns_file, "w") as f:
            pass  # Creates an empty file


def _set_project_paths(options, res: str) -> None:
    """
    Set all canonical project paths from one central directory mapping.
    """

    res = os.path.abspath(os.path.normpath(res))
    options.result_files_directory = res

    options.database_directory = os.path.join(res, "database.db")

    for attr, dirname in PROJECT_DIRS.items():
        setattr(options, attr, os.path.join(res, dirname))

    options.filtered_blast_table = os.path.join(res, "filtered_blast_results_table")
    options.divergent_output_file = os.path.join(res, "div_output_file.faa")

    options.csb_output_file = os.path.join(options.Csb_directory, "Csb_output.txt")
    options.gene_clusters_file = os.path.join(
        options.Csb_directory, "All_gene_clusters.txt"
    )
    options.data_computed_Instances_json = os.path.join(
        options.Csb_directory, "csb_instances.json"
    )

    options.report_output_file = os.path.join(res, "Report.txt")
    options.thresholds_output_file = os.path.join(res, "Thresholds.txt")


def _create_project_dirs(options) -> None:
    """
    Create all project subdirectories defined in PROJECT_DIRS.
    """

    for attr in PROJECT_DIRS:
        path = getattr(options, attr)
        try:
            os.makedirs(path, exist_ok=True)
        except Exception as e:
            logger.error(f"Could not create directory {path}: {e}")
            sys.exit(1)


def prepare_result_space(options, project: str = "project") -> None:
    """
    Creates and sets up the results directory for a project, including all needed subdirectories.
    Updates multiple attributes of options.
    """

    if getattr(options, "location", None):
        options.location = os.path.abspath(
            os.path.normpath(os.path.expanduser(options.location))
        )

    if getattr(options, "result_files_directory", None):
        options.result_files_directory = os.path.abspath(
            os.path.normpath(os.path.expanduser(options.result_files_directory))
        )

    default_dir = os.path.join(options.location, "results")

    fasta_and_query_provided = (
        getattr(options, "fasta_file_directory", None)
        and getattr(options, "query_file", None)
        and os.path.isdir(options.fasta_file_directory)
        and os.path.isfile(options.query_file)
    )

    if options.new_project:
        if not os.path.isdir(default_dir):
            try:
                os.makedirs(default_dir, exist_ok=True)
                logger.info(f"Created results directory: {default_dir}")
            except Exception as e:
                logger.error(f"No writing rights in default results directory. {e}")
                sys.exit(1)

        options.result_files_directory = create_project(default_dir, project)
        write_options_to_tsv(options, options.result_files_directory)

    elif fasta_and_query_provided:
        if not os.path.isdir(options.result_files_directory):
            try:
                os.makedirs(options.result_files_directory, exist_ok=True)
                logger.info(
                    f"Created results directory: {options.result_files_directory}"
                )
            except Exception as e:
                logger.error(f"No writing rights in results directory. {e}")
                sys.exit(1)

        options.result_files_directory = create_project(
            options.result_files_directory, project
        )
        options.new_project = True
        write_options_to_tsv(options, options.result_files_directory)

    elif not is_project_folder(options):
        if not os.path.isdir(options.result_files_directory):
            try:
                os.makedirs(options.result_files_directory, exist_ok=True)
                logger.info(f"Created results dir: {options.result_files_directory}")
            except Exception as e:
                logger.error(
                    f"No writing rights for directory {options.result_files_directory}\n {e}"
                )
                sys.exit(1)

        options.result_files_directory = create_project(
            options.result_files_directory, project
        )
        options.new_project = True
        write_options_to_tsv(options, options.result_files_directory)

    _set_project_paths(options, options.result_files_directory)
    _create_project_dirs(options)

    if not options.new_project and options.stage < 3:
        logger.warning("Existing project directory detected. Setting start stage to 3.")
        options.stage = 3


def create_project(directory: str, projectname: str = "project") -> str:
    """
    Creates a new timestamped project subdirectory in the given directory.

    Args:
        directory (str): Parent directory
        projectname (str): Project subfolder name (default 'project')

    Returns:
        str: Full path to the created project directory
    """

    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")  # z. B. "2025-04-16_14-53-21"
    directory = os.path.join(directory, f"{timestamp}_{projectname}")

    try:
        os.mkdir(directory)
    except Exception as e:
        logger.error(f"No writing rights. Project was not created.\n {e}")
        sys.exit(1)

    return directory


def is_project_folder(options) -> bool:
    """
    Checks if a result_files_directory contains a valid project structure.
    If a database is found, adjusts options accordingly.
    """

    if options.database_directory and os.path.isfile(options.database_directory):
        options.result_files_directory = os.path.dirname(options.database_directory)

    try:
        database_path = os.path.join(options.result_files_directory, "database.db")

        if not os.path.isfile(database_path):
            logger.warning(
                "No database file found in expected directory. Searching recursively..."
            )

            found_db = None
            for root, dirs, files in os.walk(options.result_files_directory):
                if "database.db" in files:
                    found_db = os.path.join(root, "database.db")
                    break

            if found_db:
                logger.info(f"Found database at {found_db}")
                options.result_files_directory = os.path.dirname(found_db)
            else:
                logger.warning(
                    "No database found at all. Existing project was not found."
                )
                return False

        for subdir in PROJECT_DIRS.values():
            path = os.path.join(options.result_files_directory, subdir)
            if not os.path.isdir(path):
                logger.info(f"Missing directory {subdir} created")
                os.makedirs(path, exist_ok=True)

        self_query_path = os.path.join(options.result_files_directory, "self_blast.faa")
        if os.path.isfile(self_query_path):
            options.self_query = self_query_path
        else:
            logger.warning(
                "No internal query file (self_blast.faa) found. Will need to create it later."
            )

        if options.glob_table is None:
            for file_name in os.listdir(options.result_files_directory):
                if file_name.startswith("filtered_"):
                    options.glob_table = os.path.join(
                        options.result_files_directory, file_name
                    )
                    break

    except Exception as e:
        logger.error(
            f"An error occurred while checking for existing project folder: {e}"
        )
        sys.exit(1)

    return True


def any_process_args_provided(args, default_values: dict) -> bool:
    """
    Checks if any CLI/process arguments are different from their default values.

    Args:
        args: The argument object (e.g., argparse.Namespace)
        default_values (dict): Mapping of argname -> default

    Returns:
        bool: True if any arg has a value different from default, else False.
    """
    for arg, default in default_values.items():
        if getattr(args, arg) != default:
            # print(f"{arg} was not {default} but {getattr(args, arg)}")
            return True
    return False


def write_options_to_tsv(
    options, output_directory: str, filename: str = "parameters_summary.tsv"
) -> None:
    """
    Writes all parameters from an options object as a TSV file.

    Args:
        options: The argument/options object (argparse.Namespace or custom class)
        output_directory (str): Output directory for the file
        filename (str): Name of output file (default: parameters_summary.tsv)

    Output:
        TSV file listing all parameters and their values.

    Example Output:
        /mydir/parameters_summary.tsv

    """
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    output_path = os.path.join(output_directory, filename)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Parameter", "Value"])
        for key, value in vars(options).items():
            writer.writerow([key, value])

    logger.info(f"Parameters saved: {output_path}")
