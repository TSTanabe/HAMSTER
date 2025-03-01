#!/usr/bin/python



import csv
import heapq
import os
import sqlite3
import subprocess

from collections import defaultdict
from multiprocessing import Pool, Manager
import numpy as np




def parallel_cross_validation(options):
    #####
    # MAIN VALIDATION ROUTINE
    #####
    print("Initilizing cross-validation")
    
    #Prepare shared data
    all_seq_number = count_fasta_headers(options.sequence_faa_file)
    alignment_files = [os.path.join(options.fasta_output_directory, f) for f in os.listdir(options.fasta_output_directory) if f.endswith(".fasta_aln")]
    files_number = len(alignment_files)
    
    print("Cutoff and performance validation without folds")

    #Prepare task-specific arguments
    #TODO only include the ones that are in the included and exlude the ones that are excluded by options
    args_list = [(alignment_file, options, all_seq_number, files_number) for alignment_file in alignment_files]

    # Initial validation of full HMMs and assigning cutoffs
    with Pool(processes=options.cores) as pool:
        pool.starmap(initial_process_folds, args_list)
    
    if options.cross_validation_deactivated: # if deactivated skip the cross validation
        return
        
    print(f"Filter HMMs by performance MCC >{options.MCC_threshold} and prepare for cross-validation")
    
    
    # Filter out the HMMs with a performance below threshold
    initial_HMM_validation_MCC_dict = process_MCC_files(options.fasta_alignment_directory,"_MCC.txt")
    initial_HMM_validation_MCC_dict = filter_dict_by_value(initial_HMM_validation_MCC_dict,options.MCC_threshold)
    
    alignment_files = filter_filepaths(alignment_files, initial_HMM_validation_MCC_dict)
    files_number = len(alignment_files)
    
    args_list = [(alignment_file, options, all_seq_number, files_number) for alignment_file in alignment_files]

    print("Cross-validation initilized")
    #Makes the cross validation and assigns the cutoff values 
    with Pool(processes=options.cores) as pool:
        
        #starts the workers
        pool.starmap(process_cross_folds,args_list)
    print(f"Processed all validations.")

    return

    
def create_hmms_from_msas(directory, target_dir, ending="fasta_aln", extension="hmm_cv", cores=1):
    # Ensure the target directory exists
    os.makedirs(target_dir, exist_ok=True)

    # Get a list of all MSA files with the specified ending in the directory
    msa_files = [f for f in os.listdir(directory) if f.endswith(ending)]

    # Run hmmbuild on each MSA file
    for msa_file in msa_files:
        msa_path = os.path.join(directory, msa_file)
        hmm_file = os.path.join(target_dir, msa_file.replace(ending, extension))

        # Skip if the HMM file already exists in the target directory
        if os.path.isfile(hmm_file):
            print(f"HMM {hmm_file} already exists. Skipping...")
            continue

        # Construct the hmmbuild command
        hmmbuild_cmd = f"hmmbuild --amino --cpu {cores} {hmm_file} {msa_path}"
        
        # Run the hmmbuild command
        print(f"Running hmmbuild for {msa_file} -> {hmm_file}")
        subprocess.run(hmmbuild_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def count_fasta_headers(fasta_file):
    with open(fasta_file, "r") as f:
        return sum(1 for line in f if line.startswith(">"))

def get_target_sets(directory):
    """
    Finds all files in the given directory that start with 'superfamily_' and end with '.faa'.
    Stores the file paths in a dictionary with keys as the filename without the prefix and extension.
    
    Parameters:
        directory (str): The directory to search for files.

    Returns:
        dict: A dictionary with keys as modified filenames and values as file paths.
    """
    result = {}
    
    for file_name in os.listdir(directory):
        if file_name.startswith("superfamily_") and file_name.endswith(".faa"):
            # Remove 'superfamily_' prefix and '.faa' extension to form the key
            key = file_name[len("superfamily_"):-len(".faa")]
            # Construct the full path and store in the dictionary
            result[key] = os.path.join(directory, file_name)
    
    return result
    
####################################################################################
####################  Initial validation procedure start  ##########################
####################################################################################

def initial_process_folds(alignment_file, options, all_sequence_number, files_number):
    """
    Runs in parallel to test the full HMM against a single test target.
    Calculates the confusion matrix and returns:
    - HMM name
    - MCC
    - Cutoffs
    - Confusion matrix from validation
    """
    
    # Unpacking arguments
    #alignment_file, options, all_sequence_number, files_number = args_tuple

    # Define directory and file paths
    initial_validation_directory = options.fasta_alignment_directory 
    hv_subfolder_name = os.path.splitext(os.path.basename(alignment_file))[0]
    hv_directory = os.path.join(initial_validation_directory, hv_subfolder_name)
    
    hmm = os.path.join(options.Hidden_markov_model_directory, f"{hv_subfolder_name}.hmm")
    
    # Check if the HMM file exists
    if not os.path.isfile(hmm):
        return
    
    os.makedirs(hv_directory, exist_ok=True)
    
    if os.path.exists(os.path.join(hv_directory, f"{hv_subfolder_name}_matrices.tsv")):
        return None  # This was already finished return to next HMM validation
    
    print(f"Initial validation of HMM {hv_subfolder_name}")    
    
    # Determine target file
    name = hv_subfolder_name.split('_')[-1]
    sequence_faa_file = options.targeted_sequence_faa_file_dict.get(name, options.sequence_faa_file)

    # Define the output report path
    report_name = os.path.splitext(os.path.basename(alignment_file))[0]
    output_report = os.path.join(hv_directory, f"{report_name}.report")
    
    # Run the HMMsearch
    HMMsearch(hmm, sequence_faa_file, output_report, 1)

    # Check if the report file was generated    
    hv_report = output_report

    if not os.path.exists(hv_report):
        return None  # Skip if no report was created

        
    # Extract true positives from the alignment file
    TP_seqIDs = extract_record_ids_from_alignment(alignment_file)
    TN = all_sequence_number - len(TP_seqIDs)
        
    # Calculate cutoff thresholds and MCC
    optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix, hit_id_distribution, all_matrices = cutoffs(TP_seqIDs, TN, hv_report)
    print(f"Finished HMM {hv_subfolder_name} cutoffs: {optimized_cutoff} {trusted_cutoff} {noise_cutoff} and MCC: {MCC}")
    
    #Write down the matrices
    save_matrices_to_tsv(all_matrices, os.path.join(hv_directory, f"{hv_subfolder_name}_matrices.tsv"))
        
    #Calculate the optimum        
    strict_report = calculate_performance_and_sum_matrices([[optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix]])  # Returns a tuple
        
    # Write down the matrices and cutoffs in the CV subfolder
    with open(os.path.join(hv_directory, f"{hv_subfolder_name}_thresholds.txt"), 'w') as writer:
        writer.write(f"{hv_subfolder_name}\t{strict_report[2]}\t{strict_report[3]}\t{strict_report[4]}\n")
        
    with open(os.path.join(hv_directory, f"{hv_subfolder_name}_MCC.txt"), 'w') as writer:
        fold_matrix = strict_report[1]
        writer.write(f"{hv_subfolder_name}\t{MCC}\t{matrix}\n")
    
    
    #More informative output of positive and negative hits
    if not options.hit_report_deactivated:
        hit_id_reports(hv_directory, options.database_directory, [hit_id_distribution])
    
    return None



def process_MCC_files(directory, file_ending):
    """
    Recursively finds all files with a specified ending, concatenates their content, 
    splits into the first and last column, and filters the entries to return a dictionary.

    Parameters:
        directory (str): Path to the root directory to search.
        file_ending (str): File extension to look for (e.g., ".txt").

    Returns:
        dict: A dictionary with the first column (name) as keys and the last column (float) as values.
    """
    result_dict = {}

    # Recursively find all files with the specified ending
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_ending):
                file_path = os.path.join(root, file)

                # Process each file
                with open(file_path, 'r') as f:
                    for line in f:
                        # Split line into columns, assuming tab-separated
                        columns = line.strip().split('\t')
                        if len(columns) >= 2:  # Ensure there are at least two columns
                            name = columns[0].strip()  # First column (name)
                            try:
                                value = float(columns[1].strip())  # Last column (float)
                                result_dict[name] = value  # Add to dictionary
                            except ValueError:
                                # Skip lines where the last column is not a float
                                print(f"Skipping line due to invalid float value: {line.strip()}")

    return result_dict

def filter_dict_by_value(d, threshold):
    """
    Removes all key-value pairs from a dictionary where the values are below a given threshold.

    Parameters:
    d (dict): The dictionary to filter.
    threshold (float): The threshold value. Key-value pairs with values below this will be removed.

    Returns:
    dict: A new dictionary with only the key-value pairs meeting the threshold criteria.
    """
    return {key: value for key, value in d.items() if isinstance(value, (int, float)) and value >= threshold}

    
def filter_filepaths(filepaths, name_dict):
    """
    Removes file paths from a list if the filename (without extension) is not in the keys of the given dictionary.

    Parameters:
    filepaths (list): A list of file paths to filter.
    name_dict (dict): A dictionary whose keys will be used to determine whether to keep a file path.

    Returns:
    list: A filtered list of file paths.
    """
    valid_names = set(name_dict.keys())
    filtered_filepaths = [
        filepath for filepath in filepaths
        if os.path.splitext(os.path.basename(filepath))[0] in valid_names
    ]
    return filtered_filepaths
    
####################################################################################
####################  Cross Validation procedure start  ############################
####################################################################################

def process_cross_folds(alignment_file, options, all_sequence_number, files_number):
    # This process is running in parallel
    # Prepare the folds and HMMs
    #alignment_file, options, all_sequence_number, files_number  = args_tuple
    # Generate the cross-validation directory
    cv_directory, cv_subfolder_name = setup_cross_validation_directory(alignment_file, options.cross_validation_directory)


    print(f"Cross-validation of HMM {cv_subfolder_name}")
    # Check if the task was already completed
    if os.path.exists(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_matrices.txt")):
        return None


    create_cross_validation_sets(alignment_file, cv_directory)
    create_hmms_from_msas(cv_directory, cv_directory, "cv")
    
    # Retrieve the HMMs and target sequence files
    cv_hmms = get_cv_hmms(cv_directory)
    sequence_faa_file = get_target_sequence_file(cv_subfolder_name, options)
    
    # Run HMMsearch for each cross-validation model
    run_hmmsearch_for_cv_hmms(cv_hmms, sequence_faa_file)
    
    # Retrieve reprots and calculate cutoffs
    cv_reports = get_cv_reports(cv_directory)
    process_cutoffs_and_save_results(alignment_file, cv_reports, all_sequence_number, cv_directory, cv_subfolder_name, options.database_directory)

    # Cleanup temporary files
    cleanup_temp_files(cv_directory)
    
    
    return None

    
### HELPER ROUTINES FOR PROCESS CROSS FOLDS
    
def setup_cross_validation_directory(alignment_file, cross_validation_directory):
    """Generates the cross-validation directory and returns its path."""
    cv_subfolder_name = os.path.splitext(os.path.basename(alignment_file))[0]
    cv_directory = os.path.join(cross_validation_directory, cv_subfolder_name)
    return cv_directory, cv_subfolder_name   
    
def get_cv_hmms(cv_directory):
    """Returns a list of all HMM files used in cross-validation."""
    return [os.path.join(cv_directory, f) for f in os.listdir(cv_directory) if f.endswith(".hmm_cv")]

def get_target_sequence_file(cv_subfolder_name, options):
    """Determines the correct sequence file based on the alignment file name."""
    name = cv_subfolder_name.split('_')[-1]
    return options.targeted_sequence_faa_file_dict.get(name, options.sequence_faa_file)

def run_hmmsearch_for_cv_hmms(cv_hmms, sequence_faa_file):
    """Executes HMMsearch for each cross-validation HMM model."""
    for hmm in cv_hmms:
        HMMsearch(hmm, sequence_faa_file, hmm + ".report", 1)

def get_cv_reports(cv_directory):
    """Retrieves all HMMsearch reports from the cross-validation directory."""
    return [os.path.join(cv_directory, f) for f in os.listdir(cv_directory) if f.endswith(".report")]
            

def process_cutoffs_and_save_results(alignment_file, cv_reports, all_sequence_number, cv_directory, cv_subfolder_name, database_directory):
    """Processes cutoffs from cross-validation reports, saves results, and generates hit reports."""
    # Get true positive sequence IDs
    TP_seqIDs = extract_record_ids_from_alignment(alignment_file)
    TN = all_sequence_number - len(TP_seqIDs)

    # Calculate cutoffs for each report
    model_values = []
    model_hit_distribution = []

    for report in cv_reports:
        optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix, hit_id_distribution, all_matrices = cutoffs(TP_seqIDs, TN, report)
        model_hit_distribution.append(hit_id_distribution)  # Save distribution
        model_values.append([optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix])

    # Determine final performance metrics
    strict_report = calculate_performance_and_sum_matrices(model_values)

    # Save results to files
    save_thresholds(cv_directory, cv_subfolder_name, strict_report)
    save_matrices(cv_directory, cv_subfolder_name, strict_report[1])

    # Generate detailed hit reports
    # hit_id_reports(cv_directory, database_directory, model_hit_distribution)

def save_thresholds(cv_directory, cv_subfolder_name, strict_report):
    """Writes the calculated cutoff thresholds to a file."""
    with open(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_thresholds.txt"), 'w') as writer:
        writer.write(f"{cv_subfolder_name}\t{strict_report[2]}\t{strict_report[3]}\t{strict_report[4]}\n")

def save_matrices(cv_directory, cv_subfolder_name, fold_matrix):
    """Writes the confusion matrices to a file."""
    with open(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_matrices.txt"), 'w') as writer:
        writer.write('\n'.join(['\t'.join(map(str, row)) for row in fold_matrix]))

def cleanup_temp_files(cv_directory):
    """Removes temporary training data files from the cross-validation directory."""
    for file in os.listdir(cv_directory):
        if file.startswith("training_data_"):
            os.remove(os.path.join(cv_directory, file))
            
            
            
## CROSS FOLD CREATION            
            
def create_cross_validation_sets(alignment_file, output_directory, ending=".cv", num_folds=5):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # **Change `sequences` to a dictionary**
    sequences = {}
    record_ids = []

    with open(alignment_file, "r") as f:
        record_id = None
        sequence_lines = []
        
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if record_id is not None:
                    sequences[record_id] = "".join(sequence_lines)
                    sequence_lines = []
                
                record_id = line[1:]  # Remove '>'
                record_ids.append(record_id)
            else:
                sequence_lines.append(line)
        
        if record_id is not None:
            sequences[record_id] = "".join(sequence_lines)

    num_sequences = len(sequences)
    fold_size = num_sequences // num_folds

    # **Rewriting fold generation using dictionary**
    record_id_list = list(sequences.keys())  # Get record IDs in a list

    for fold in range(num_folds):
        start_idx = fold * fold_size
        end_idx = (fold + 1) * fold_size
        
        train_record_ids = record_id_list[:start_idx] + record_id_list[end_idx:]  # Remove fold from training set
        train_file = os.path.join(output_directory, f"training_data_{fold}{ending}")

        with open(train_file, "w") as train_f:
            for record_id in train_record_ids:
                train_f.write(f">{record_id}\n{sequences[record_id]}\n")

                

def extract_record_ids_from_alignment(alignment_file):
    record_ids = set()

    # Open the alignment file and iterate over each line
    with open(alignment_file, 'r') as file:
        for line in file:
            # Check if the line starts with '>', indicating a new record ID
            if line.startswith('>'):
                # Extract the record ID by stripping the '>' and splitting by the first space
                record_id = line[1:].strip().split()[0]
                # Add the record ID to the set
                record_ids.add(record_id)

    return record_ids


def HMMsearch(library, input_file, output_file, cores=1):
    
    if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
        return
        
    # Construct the hmmsearch command
    hmmsearch_cmd = f"hmmsearch --incdomT 10 --noali --domtblout {output_file} --cpu {cores} {library} {input_file}"

    # Run the hmmsearch command using subprocess
    status = subprocess.call(hmmsearch_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if status > 0:
        print(f"ERROR: {hmmsearch_cmd}\nDied with exit code: {status}")
        # You may raise an exception here or handle the error as needed

def calculateMetric(metric, TP, FP, FN, TN):
    if metric == "MCC":
        numerator = TP * TN - FP * FN
        denominator = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
        mcc = numerator / denominator if denominator != 0 else 0
        return mcc
    else:
        raise ValueError("Unsupported metric. Only 'MCC' is supported.")


def cutoffs(true_positives,true_negatives,report_filepath):
    """
    Calculate the cutoffs and the confusion matrix. Sort the proteinIDs into the categories TP,FP,FN,TN
    """
    
    # Initialize confusion matrix variables and identifier data structures
    sum_TP = len(true_positives) #true_positives is a set, sum_TP is the total number of TP
    TP, FP, FN, TN = 0, 0, 0, true_negatives
    
    trusted_cutoff = float("-inf")
    noise_cutoff = 10
    optimized_cutoff = 0
    MCC = 0
    matrix = []
    all_matrices = []

    true_positives_set = set(true_positives)
    processed_TPs = set()  # Speichert bereits verarbeitete TPs

    # Track hit occurrences and their lowest scores
    true_positive_dict = defaultdict(float)
    false_positive_dict = defaultdict(float)
    false_negative_dict = {key: 0 for key in true_positives}
        
    with open(report_filepath, 'r') as file:
        
        for line in file:
            if line.startswith('#'):
                continue

            # Split the line by any whitespace, treating consecutive spaces as a single separator
            fields = line.strip().split()
            hit_id, bitscore = fields[0], float(fields[13])
            
            # Process true positives                        
            if hit_id in true_positives_set:
                if hit_id in processed_TPs:
                    continue #skip doublicate hit
                
                TP += 1
                true_positive_dict[hit_id] = bitscore # collect the positive hit identifiers set
                processed_TPs.add(hit_id)  # Ensure TPs are only counted once

                FN = sum_TP - TP #wenn tp die hit_ids speichert dann könnte man hier direkt die FN hit IDs rausfiltern über vergleiche von sets
                
                # Store last TP as trusted bitscore                
                if FP == 0:
                    trusted_cutoff = bitscore
                
                # Update false negative dict
                if hit_id in false_negative_dict: #track the bitscores, TPs will be removed from here when better MCC is found
                    false_negative_dict[hit_id] = bitscore
                    
            # Process false positives
            else:
                if hit_id in false_positive_dict:
                    continue #skip doublicate hit
                
                FP += 1
                false_positive_dict[hit_id] = bitscore # collect the false positive hit identifiers
                TN = true_negatives - FP
            
            # Calculate matthews correlation coefficient            
            new_MCC = calculateMetric("MCC",TP,FP,FN,TN)
            #print(f"{TP},{FP},{FN},{TN} with hit_id {hit_id} as next positive hit")
            #print(new_MCC, " > ", MCC)
            
            # Store the confusion matrix     
            all_matrices.append([TP, FP, FN, TN, bitscore, new_MCC])
            
            # Update MCC if improved
            if new_MCC > MCC:
                MCC = new_MCC
                optimized_cutoff = bitscore
                matrix = [TP,FP,FN,TN]
                
                # Remove all TPs at this point from the false_negative_dict thus scores can be recorded for all hits but only those in the
                # below the threshold will be returned              
                false_negative_dict = {k: v for k, v in false_negative_dict.items() if k not in processed_TPs}
                
                if FN != len(false_negative_dict): # This should never happen
                    print(f"ERROR: Asynchrone false negative count and dictionary in the cutoff routine {FN}/{len(false_negative_dict)}")
            
            # Stop if all TPs are found
            if TP == sum_TP:
                noise_cutoff = bitscore
                break #MCC there is no better MCC without more TP
                
            if noise_cutoff > bitscore:
                noise_cutoff = bitscore
    
    # Filter false positives above cutoff
    false_positive_dict_filtered = {key: value for key, value in false_positive_dict.items() if value >= optimized_cutoff}
    
    # Prepare hit id distribution for return
    hit_id_distribution = [true_positive_dict, false_positive_dict_filtered, false_negative_dict] #distribution of hit ids at optimum MCC
    
    return optimized_cutoff,trusted_cutoff,noise_cutoff,MCC,matrix, hit_id_distribution, all_matrices


def calculate_performance_and_sum_matrices(cross_validation_folds):
    """
    Computes the best performance metrics from multiple cross-validation folds.

    Args:
        cross_validation_folds (list of tuples): Each tuple contains:
            (optimized_cutoff, trusted_cutoff, noise_cutoff, mcc, confusion_matrix)

    Returns:
        tuple: (summed_confusion_matrix, all_folds_matrices, best_optimized_cutoff, 
                highest_trusted_cutoff, lowest_noise_cutoff)
    """

    # Initilize
    highest_trusted_cutoff = float("-inf")
    lowest_noise_cutoff = float("inf")
    best_mcc = float("-inf")
    best_optimized_cutoff = None
    sum_matrix = [0, 0, 0, 0]
    fold_matrices = []
    
    for optimized_cutoff, trusted_cutoff, noise_cutoff, mcc, matrix in cross_validation_folds:
        fold_matrices.append(matrix)
        
        # Update highest trusted cutoff
        highest_trusted_cutoff = max(highest_trusted_cutoff, trusted_cutoff)
        lowest_noise_cutoff = min(lowest_noise_cutoff, noise_cutoff)

        # Update best MCC and corresponding optimized cutoff
        if mcc > best_mcc or (mcc == best_mcc and optimized_cutoff < best_optimized_cutoff):
            best_mcc = mcc
            best_optimized_cutoff = optimized_cutoff

    # Ensure no value is infinite or negative
    highest_trusted_cutoff = max(highest_trusted_cutoff, 10)
    lowest_noise_cutoff = max(lowest_noise_cutoff, 10)
    best_optimized_cutoff = best_optimized_cutoff if best_optimized_cutoff and best_optimized_cutoff > 0 else 10
    
    return (sum_matrix, fold_matrices, best_optimized_cutoff, highest_trusted_cutoff, lowest_noise_cutoff)



#######################################################################
##################    Hit report routines      ########################
#######################################################################



def write_report_to_file(cv_directory, dict_label, sorted_data):
    """
    Writes the combined data (key, score, count, genomic neighborhood) to a file.
    
    Args:
        cv_directory (str): Path to the directory where the report file will be saved.
        dict_label (str): Label indicating the type of report (TP, FP, FN).
        sorted_data (list): A list of tuples containing (key, score, count, genomic neighborhood).
    
    Returns:
        None
    """
    # Create the file path
    report_file = os.path.join(cv_directory, f"{dict_label}_hits.hit_report")
    
    # Write the report data to the file
    with open(report_file, 'w') as f:
        # Write the header
        f.write("Key\tScore\tHits\tGenomic_Neighborhood\n")
        
        # Write each row of the sorted data
        for key, score, count, neighborhood in sorted_data:
            f.write(f"{key}\t{score}\t{count}\t{neighborhood}\n")



def write_report_to_iTolBinary(cv_directory, dict_label, sorted_data):
    # Create the file path
    report_file = os.path.join(cv_directory, f"{dict_label}_iTolBinary.hit_report")
    
    # Prepare the header for the iTOL binary dataset
    header = """DATASET_BINARY
SEPARATOR COMMA
DATASET_LABEL,Full Binary Dataset
COLOR,#ff0000
FIELD_SHAPES,1
FIELD_LABELS,ID
FIELD_COLORS,#ff0000
DATA
"""
    
    # Open the output file for writing
    with open(report_file, 'w') as outfile:
        outfile.write(header)
        for key, score, count, neighborhood in sorted_data:
            outfile.write(f"{key},1\n")



def hit_id_reports(cv_directory, database, data):
    """
    Generates a report with the true positive, false positive, and false negative hits, 
    and how often they were found through the cross-validation folds.
    
    Args:
        cv_directory (str): Path to the cross-validation directory.
        database (str): Path to the database file.
        data (list): A list of lists, each containing TP_dict, FP_dict, FN_dict.
    
    Returns:
        None
    """
    # Initialize empty dictionaries to hold key occurrence counts for TP, FP, and FN
    key_counts = {label: defaultdict(int) for label in ['TP', 'FP', 'FN']}
    
    # lowest hit score for each proteinID
    hit_scores = {label: defaultdict(lambda: float("inf")) for label in ['TP', 'FP', 'FN']}

    # Process all TP, FP, FN in a single pass
    for tp_dict, fp_dict, fn_dict in data:
        for dict_label, hit_dict in zip(['TP', 'FP', 'FN'], [tp_dict, fp_dict, fn_dict]):
            for hit, score in hit_dict.items():
                key_counts[dict_label][hit] += 1
                hit_scores[dict_label][hit] = min(hit_scores[dict_label][hit], score)  # Track lowest score

    # Batch-fetch gene vicinity information for all keys in one DB call per category
    gene_vicinity = {label: fetch_neighbouring_genes_with_domains(database, list(key_counts[label].keys()))
                     for label in ['TP', 'FP', 'FN']}

    # Process each category once
    for dict_label in ['TP', 'FP', 'FN']:
        score_dict = hit_scores[dict_label]
        sorted_data = sorted(
            ((key, score, key_counts[dict_label][key], gene_vicinity[dict_label].get(key, "No data")) 
             for key, score in score_dict.items()),
            key=lambda x: x[1],  # Sort by lowest score
            reverse=True
        )
        
            # Write the sorted data to a file
        write_report_to_file(cv_directory, dict_label, sorted_data)
        write_report_to_iTolBinary(cv_directory, dict_label, sorted_data)


    # At this point, the key counts and gene vicinity data for TP, FP, and FN have been processed
    return 
    



def fetch_neighbouring_genes_with_domains(database, protein_ids):
    """
    Fetches neighboring genes for a list of proteinIDs based on clusterID and gene order,
    and returns the neighboring genes with their domains. For proteinIDs without clusterID,
    only the entry for the proteinID itself is returned.

    Args:
        database (str): Pathway to the database file.
        protein_ids (list): List of proteinIDs to fetch neighbors for.

    Returns:
        dict: A dictionary where each key is a proteinID, and the value is a list of 
              neighboring genes in the same cluster, sorted by gene start position, or
              just the entry for the proteinID if it has no cluster.
              Each list entry contains the domains and the proteinID in the format 
              '_domain_proteinID'.
    """
    if not protein_ids:
        return {}
    
    chunk_size = 900    
    
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cluster_ids = set()
        
        # Step 1: Fetch all relevant proteinID, clusterIDs (chunked)
        for i in range(0, len(protein_ids), chunk_size):
            chunk = protein_ids[i:i + chunk_size]
            placeholders = ','.join('?' * len(chunk))

            cur.execute(f"SELECT DISTINCT clusterID FROM Proteins WHERE proteinID IN ({placeholders})", chunk)
            cluster_ids.update(row[0] for row in cur.fetchall() if row[0] is not None)
            
            
            
        # Step 2: Fetch all proteins that belong to these clusters (chunked)
        protein_results = []
        for i in range(0, len(cluster_ids), chunk_size):
            chunk = list(cluster_ids)[i:i + chunk_size]  # Convert set to list for slicing
            placeholders = ','.join('?' * len(chunk))

            query = f"""
            SELECT p.proteinID, p.clusterID, p.genomeID, p.start, COALESCE(d.domain, 'no_domain')
            FROM Proteins p
            LEFT JOIN Domains d ON p.proteinID = d.proteinID
            WHERE p.clusterID IN ({placeholders})
            ORDER BY p.clusterID, p.start
            """
            cur.execute(query, chunk)
            protein_results.extend(cur.fetchall())


    # Step 3: Organize results efficiently
    cluster_dict = defaultdict(list)
    protein_cluster_map = {}

    for protein_id, cluster_id, genome_id, start, domain in protein_results:
        domain_entry = f"{domain}_{protein_id}"
        if cluster_id:
            cluster_dict[cluster_id].append((start, domain_entry))
            protein_cluster_map[protein_id] = cluster_id
        else:
            protein_cluster_map[protein_id] = None  # No cluster assigned

    # Step 4: Sort by start position
    for cluster_id in cluster_dict:
        cluster_dict[cluster_id].sort()  # Sort by tuples first value: `start` position


    # Step 5: Construct final output
    neighbors_dict = {}
    for protein_id in protein_ids:
        cluster_id = protein_cluster_map.get(protein_id)
        if cluster_id:
            neighbors_dict[protein_id] = [protein for _, protein in cluster_dict[cluster_id]]
        else:
            neighbors_dict[protein_id] = [f"singleton_{protein_id}"]

    return neighbors_dict


def save_matrices_to_tsv(all_matrices, output_filepath):
    """
    Saves all confusion matrices to a tab-separated (.tsv) file.

    Args:
        all_matrices (list): List of confusion matrices, each in the format [TP, FP, FN, TN].
        output_filepath (str): Path to the output TSV file.

    Returns:
        str: The path to the saved file.
    """
    # Define headers
    headers = ["TP", "FP", "FN", "TN", "cutoff", "MCC"]
    
    # Write to TSV file
    with open(output_filepath, mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        
        # Write header row
        writer.writerow(headers)
        
        # Write confusion matrices
        for matrix in all_matrices:
            writer.writerow(matrix)

    return output_filepath

    
