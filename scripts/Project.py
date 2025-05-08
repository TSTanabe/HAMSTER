#!/usr/bin/python

import csv
import os
import tarfile
from datetime import datetime

def prepare_directory_structure(directory):
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


def prepare_result_space(options,project="project"):
    
    logfile = 0
    now = datetime.now()
    timestamp = str(datetime.timestamp(now))
    
    if isProjectFolder(options): #check if existing project folder was given
        
        #print(options.result_files_directory)
        #print(options.fasta_initial_hit_directory)
        
        #Define directories from the existing project
        options.database_directory = options.result_files_directory+"/database.db"
        options.fasta_initial_hit_directory = options.result_files_directory+"/Hit_list"
        options.fasta_output_directory = options.result_files_directory+"/Sequences"
        options.fasta_alignment_directory = options.result_files_directory+"/Initial_validation"        
        options.Hidden_markov_model_directory = options.result_files_directory+"/Hidden_markov_models"
        options.cross_validation_directory = options.result_files_directory+"/Cross_validation"
        options.phylogeny_directory = options.result_files_directory+"/Protein_Phylogeny"
        options.Csb_directory = options.result_files_directory+"/Collinear_syntenic_blocks"
        
        
        #Define files from the existing project
        options.divergent_output_file = options.result_files_directory+"/div_output_file.faa"
        options.csb_output_file = options.Csb_directory+"/Csb_output.txt"
        options.gene_clusters_file = options.Csb_directory+"/All_gene_clusters.txt"
        
        
        options.report_output_file = options.result_files_directory+"/Report.txt" #for CV
        options.thresholds_output_file = options.result_files_directory+"/Thresholds.txt" #for CV





    else: #if new project is created
        if not os.path.isdir(options.result_files_directory):
            try:
                os.mkdir(options.result_files_directory)
            except:
                raise Exception(f"\nERROR: No writing rights.")
        #Creates an new project and overwrites the result file directory with the project folder. All
        #other subfolders and files are stored in this project folder
        options.result_files_directory = create_project(options.result_files_directory,project)        
        write_options_to_tsv(options, options.result_files_directory)
        
        options.database_directory = options.result_files_directory+"/database.db"
        

        options.fasta_initial_hit_directory = options.result_files_directory+"/Hit_list"
        os.mkdir(options.fasta_initial_hit_directory)
        options.fasta_output_directory = options.result_files_directory+"/Sequences"
        os.mkdir(options.fasta_output_directory)
        options.fasta_alignment_directory = options.result_files_directory+"/Initial_validation"
        os.mkdir(options.fasta_alignment_directory)
        options.Hidden_markov_model_directory = options.result_files_directory+"/Hidden_markov_models"
        os.mkdir(options.Hidden_markov_model_directory)
        options.cross_validation_directory = options.result_files_directory+"/Cross_validation"
        os.mkdir(options.cross_validation_directory)
        options.phylogeny_directory = options.result_files_directory+"/Protein_Phylogeny"
        os.mkdir(options.phylogeny_directory)
        options.Csb_directory = options.result_files_directory+"/Collinear_syntenic_blocks"
        os.mkdir(options.Csb_directory)
        
        
        options.divergent_output_file = options.result_files_directory+"/div_output_file.faa"
        options.csb_output_file = options.Csb_directory+"/Csb_output.txt"
        options.gene_clusters_file = options.Csb_directory+"/All_gene_clusters.txt"
        options.data_computed_Instances_json = options.Csb_directory+"/csb_instances.json"

        options.report_output_file = options.result_files_directory+"/Report.txt"
        options.thresholds_output_file = options.result_files_directory+"/Thresholds.txt"
        


def create_project(directory, projectname="project"):
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")  # z. B. "2025-04-16_14-53-21"
    directory = os.path.join(directory, f"{timestamp}_{projectname}")
    
    try:
        os.mkdir(directory)
    except Exception:
        raise Exception("\nERROR: No writing rights.")
    
    return directory
    

def isProjectFolder(options):

    #Database may be provided, then consider this as a result directory
    if options.database_directory and os.path.isfile(options.database_directory):
        options.result_files_directory = os.path.dirname(options.database_directory)
    
    try:
        database_path = os.path.join(options.result_files_directory, "database.db")

        if not os.path.isfile(database_path):
            print("No database file found in expected directory. Searching recursively...")

            # Recursive search for database.db
            found_db = None
            for root, dirs, files in os.walk(options.result_files_directory):
                if "database.db" in files:
                    found_db = os.path.join(root, "database.db")
                    break

            if found_db:
                print(f"Found database at {found_db}")
                options.result_files_directory = os.path.dirname(found_db)
            else:
                print("No database found at all. Creating new project folder structure.")

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
                print(f"Missing directory {subdir} created")
                os.makedirs(path, exist_ok=True)


        self_query_path = os.path.join(options.result_files_directory, "self_blast.faa")
        if os.path.isfile(self_query_path):
            options.self_query = self_query_path
        else:
            print("No internal query file (self_blast.faa) found. Will need to create later.")

        # Find blast table if not already defined
        if options.glob_table is None:
            for file_name in os.listdir(options.result_files_directory):
                if file_name.startswith("filtered_"):
                    options.glob_table = os.path.join(options.result_files_directory, file_name)
                    break

    except Exception as e:
        raise Exception(f"An error occurred while checking project folder: {e}")

    return 1
    
    
    
def any_process_args_provided(args, default_values):
    for arg, default in default_values.items():
        if getattr(args, arg) != default:
            #print(f"{arg} was not {default} but {getattr(args, arg)}")
            return True
    return False
    

def write_options_to_tsv(options, output_directory, filename="parameters_summary.tsv"):
    """
    Schreibt alle Parameter aus einem argparse-Objekt als TSV-Datei.

    Args:
        options (argparse.Namespace): Dein Argument-Objekt.
        output_directory (str): Zielverzeichnis für die Datei.
        filename (str): Name der Ausgabedatei (Standard: parameters_summary.tsv).
    """
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    output_path = os.path.join(output_directory, filename)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Parameter", "Value"])
        for key, value in vars(options).items():
            writer.writerow([key, value])
    
    print(f"Parameters saved: {output_path}")   
     






