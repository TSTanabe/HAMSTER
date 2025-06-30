#!/usr/bin/python

import csv
import os
import sys
import tarfile
from datetime import datetime

from . import myUtil

logger = myUtil.logger

def prepare_directory_structure(directory: str) -> None:
    """
    Sets up the `bin` and `src` directories within the specified directory.
    Unpacks any `.tar.gz` files in `bin` and deletes the archives.
    Creates an empty `patterns` file in `src` if it does not exist.
    
    Args:
        directory (str): Path to the directory to set up.
    """
    # Define the paths for `bin` and `src` directories
    bin_dir = os.path.join(directory, 'bin')
    src_dir = os.path.join(directory, 'src')

    # Ensure `bin` and `src` directories exist
    os.makedirs(bin_dir, exist_ok=True)
    os.makedirs(src_dir, exist_ok=True)
    
    # Process `.tar.xz` files in `bin` directory
    for file_name in os.listdir(bin_dir):
        if file_name.endswith('.tar.xz'):
            file_path = os.path.join(bin_dir, file_name)
            
            # Unpack the .tar.xz file
            with tarfile.open(file_path, 'r:xz') as tar:
                tar.extractall(bin_dir)
            
            # Delete the original .tar.xz file
            os.remove(file_path)
                
    # Create an empty `patterns` file in `src` if it does not exist
    patterns_file = os.path.join(src_dir, 'Patterns')
    if not os.path.exists(patterns_file):
        with open(patterns_file, 'w') as f:
            pass  # Creates an empty file


def prepare_result_space(options, project: str = "project") -> None:
    """
    Creates and sets up the results directory for a project, including all needed subdirectories.

    Args:
        options (Options): Has .location, .result_files_directory, etc.
        project (str): Project name suffix (default: 'project')

    Output:
        Updates multiple attributes of options (see docstring).

    Example Input:
        options.location = '/home/user/myproject'
        options.result_files_directory = '/home/user/myproject/results'
    """
    
    now = datetime.now()
    timestamp = str(datetime.timestamp(now))
    

    # 2. Ensure the results directory 
    default_dir = os.path.join(options.location, "results")

    # Use standard directory if given
    if options.new_project:
        # Make a new project in default directory
        if not os.path.isdir(default_dir):
            try:
                os.mkdir(default_dir)
                logger.info(f"Created results directory: {default_dir}")
            except Exception as e:
                 logger.error(f"No writing rights in default results directory. {e}")
                 sys.exit(1)
                 
        options.result_files_directory = create_project(default_dir, project)
        write_options_to_tsv(options, options.result_files_directory)

    # User-defined directory: check for project
    elif not isProjectFolder(options):
        if not os.path.isdir(options.result_files_directory):
            try:
                os.mkdir(options.result_files_directory)
                logger.info(f"Created results dir: {options.result_files_directory}")
            except Exception as e:
                logger.error(f"No writing rights for directory {options.result_files_directory}\n {e}")
                sys.exit(1)
        
        options.result_files_directory = create_project(options.result_files_directory, project)
        options.new_project = True
        write_options_to_tsv(options, options.result_files_directory)

        


    # Project structure, applies in both cases
    options.database_directory = options.result_files_directory + "/database.db"
    options.fasta_initial_hit_directory = options.result_files_directory + "/Hit_list"
    options.fasta_output_directory = options.result_files_directory + "/Sequences"
    options.fasta_alignment_directory = options.result_files_directory + "/Initial_validation"
    options.Hidden_markov_model_directory = options.result_files_directory + "/Hidden_markov_models"
    options.cross_validation_directory = options.result_files_directory + "/Cross_validation"
    options.phylogeny_directory = options.result_files_directory + "/Protein_Phylogeny"
    options.Csb_directory = options.result_files_directory + "/Collinear_syntenic_blocks"
    options.filtered_blast_table = options.result_files_directory + "/filtered_blast_results_table"


    options.divergent_output_file = options.result_files_directory + "/div_output_file.faa"
    options.csb_output_file = options.Csb_directory + "/Csb_output.txt"
    options.gene_clusters_file = options.Csb_directory + "/All_gene_clusters.txt"
    options.data_computed_Instances_json = options.Csb_directory + "/csb_instances.json"
    options.report_output_file = options.result_files_directory + "/Report.txt"
    options.thresholds_output_file = options.result_files_directory + "/Thresholds.txt"

    # Create required directories
    for path in [
        options.fasta_initial_hit_directory,
        options.fasta_output_directory,
        options.fasta_alignment_directory,
        options.Hidden_markov_model_directory,
        options.cross_validation_directory,
        options.phylogeny_directory,
        options.Csb_directory
    ]:
        if not os.path.exists(path):
            os.mkdir(path)
    
    # 5. If using a pre-existing project, force pipeline to start at stage 3
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
    

def isProjectFolder(options) -> bool:
    """
    Checks if a result_files_directory contains a valid project structure.
    If a database is found, adjusts options accordingly.

    Args:
        options (Options): Has .result_files_directory etc.

    Returns:
        bool: True if project folder found or created, else False.
    """

    #Database may be provided, then consider this as a result directory
    if options.database_directory and os.path.isfile(options.database_directory):
        options.result_files_directory = os.path.dirname(options.database_directory)
    
    try:
        database_path = os.path.join(options.result_files_directory, "database.db")

        if not os.path.isfile(database_path):
            logger.warning("No database file found in expected directory. Searching recursively...")

            # Recursive search for database.db
            found_db = None
            for root, dirs, files in os.walk(options.result_files_directory):
                if "database.db" in files:
                    found_db = os.path.join(root, "database.db")
                    break

            if found_db:
                logger.info(f"Found database at {found_db}")
                options.result_files_directory = os.path.dirname(found_db)
            else:
                logger.warning("No database found at all. Creating new project folder structure.")

        # Now ensure required directories exist
        required_dirs = [
            "Sequences",
            "Initial_validation",
            "Hit_list",
            "Hidden_markov_models",
            "Cross_validation",
            "Collinear_syntenic_blocks"
        ]

        for subdir in required_dirs:
            path = os.path.join(options.result_files_directory, subdir)
            if not os.path.isdir(path):
                logger.info(f"Missing directory {subdir} created")
                os.makedirs(path, exist_ok=True)


        self_query_path = os.path.join(options.result_files_directory, "self_blast.faa")
        if os.path.isfile(self_query_path):
            options.self_query = self_query_path
        else:
            logger.warning("No internal query file (self_blast.faa) found. Will need to create it later.")


        # Find blast table if not already defined
        if options.glob_table is None:
            for file_name in os.listdir(options.result_files_directory):
                if file_name.startswith("filtered_"):
                    options.glob_table = os.path.join(options.result_files_directory, file_name)
                    break

    except Exception as e:
        logger.error(f"An error occurred while checking project folder: {e}")
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
            #print(f"{arg} was not {default} but {getattr(args, arg)}")
            return True
    return False
    

def write_options_to_tsv(options, output_directory: str, filename: str = "parameters_summary.tsv") -> None:
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
     






