#!/usr/bin/python



import csv
import heapq
import os
import sqlite3
import subprocess
import pickle
from . import myUtil

from collections import defaultdict
from multiprocessing import Pool, Manager
import numpy as np




def initial_self_recognition_validation(options):
    #####
    # MAIN VALIDATION ROUTINE
    #####
    print("[INFO] Initilize validation of training dataset recognition")
    
    # Prepare shared data
    all_seq_number = count_fasta_headers(options.sequence_faa_file) # Number of all TP + TN
    alignment_files = [os.path.join(options.fasta_output_directory, f) for f in os.listdir(options.fasta_output_directory) if f.endswith(".fasta_aln")]
    files_number = len(alignment_files)
    
    print("[INFO] Calculate cutoffs and performance on training dataset without folds")

    # Prepare task-specific arguments
    args_list = [(options, alignment_file, all_seq_number, files_number) for alignment_file in alignment_files]

    # Initial validation of full HMMs and assigning cutoffs
    with Pool(processes=options.cores) as pool:
        results = pool.starmap(initial_training_data_performance, args_list)
    
    
    results_sorted = sorted([res for res in results if res is not None], key=lambda x: x[0])

    output_file = os.path.join(options.Hidden_markov_model_directory, "_ini_cutoffs.txt")
    with open(output_file, "w") as out:
        out.write("hmm_protein_name\toptimized_cutoff\ttrusted_cutoff\tnoise_cutoff\tMCC\n")
        for res in results_sorted:
            out.write("\t".join(str(x) for x in res) + "\n")
    return
    

####################################################################################
####################  Initial validation procedure start  ##########################
####################################################################################

def initial_training_data_performance(options, alignment_file, all_sequence_number, files_number):
    """
    Runs in parallel to test the full HMM against a single test target.
    Calculates the confusion matrix and returns:
    - HMM name
    - MCC
    - Cutoffs
    - Confusion matrix from validation
    """
    
    # Define directory and protein name
    initial_validation_directory = options.fasta_alignment_directory
    hmm_protein_name = os.path.splitext(os.path.basename(alignment_file))[0]
    protein_name = hmm_protein_name.split('_',1)[-1]
    
    load = os.path.join(initial_validation_directory, f"{hmm_protein_name}.ini_cutoffs_pkl")

    # Skip if results are already present
    data = myUtil.load_cache(options, hmm_protein_name, load)
    if data is not None:
        hmm_protein_name, optimized_cutoff, trusted_cutoff, noise_cutoff, best_MCC, matrix = data
        return hmm_protein_name, optimized_cutoff, trusted_cutoff, noise_cutoff, best_MCC, matrix
    
    # Define HMM for the test
    if not (hmm := os.path.join(options.Hidden_markov_model_directory, f"{hmm_protein_name}.hmm")) or not os.path.isfile(hmm):
        print(f"[WARN] HMM file not found for {hmm_protein_name}: {hmm}")
        return
    # Define target file    
    if not (sequence_faa_file := options.targeted_sequence_faa_file_dict.get(protein_name, options.sequence_faa_file)):
        print(f"[WARN] No target sequence file found for {protein_name}")
        return
    
    print(f"[INFO] Initial validation of HMM {hmm_protein_name}")    
    
    # Define the output report path
    output_report = os.path.join(initial_validation_directory, f"{hmm_protein_name}.hmmreport")
    
    # Run the HMMsearch
    HMMsearch(hmm, sequence_faa_file, output_report, 1)

    # Extract true positives from the alignment file
    TP_seqIDs = extract_record_ids_from_alignment(alignment_file)
    TN = all_sequence_number - len(TP_seqIDs)
        
    # Calculate cutoff thresholds and MCC
    optimized_cutoff, trusted_cutoff, noise_cutoff, report, best_MCC, matrix = cutoffs(TP_seqIDs, TN, output_report)
    
    # Save the labelling at the optimal position for the bigger report
    myUtil.save_cache(options, f"{hmm_protein_name}.ini_performance_pkl", report, initial_validation_directory)
    myUtil.save_cache(options, f"{hmm_protein_name}.ini_cutoffs_pkl", (hmm_protein_name, optimized_cutoff, trusted_cutoff, noise_cutoff, best_MCC, matrix), initial_validation_directory)
    
    print(f"[INFO] Finished HMM {hmm_protein_name} optimized cutoff: {optimized_cutoff} trusted cutoff: {trusted_cutoff} noise cutoff: {noise_cutoff} and MCC: {best_MCC}")
    
    return hmm_protein_name, optimized_cutoff, trusted_cutoff, noise_cutoff, best_MCC




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
####################  Cross Validation procedure start  ############################
####################################################################################

    
def parallel_cross_validation(options):    

        
    all_seq_number = count_fasta_headers(options.sequence_faa_file) # Number of all TP + TN
    alignment_files = [os.path.join(options.fasta_output_directory, f) for f in os.listdir(options.fasta_output_directory) if f.endswith(".fasta_aln")]
    files_number = len(alignment_files)

    args_list = [(alignment_file, options, all_seq_number, files_number) for alignment_file in alignment_files]

    print("Cross-validation initilized")
    #Makes the cross validation and assigns the cutoff values 
    with Pool(processes=options.cores) as pool:
        
        #starts the workers
        pool.starmap(process_cross_folds,args_list)

    return






    
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

    for report in cv_reports:
        optimized_cutoff, trusted_cutoff, noise_cutoff, report, MCC, matrix = cutoffs(TP_seqIDs, TN, report)
        model_values.append([optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix])

    # Determine final performance metrics
    strict_report = calculate_performance_and_sum_matrices(model_values)

    # Save results to files
    save_thresholds(cv_directory, cv_subfolder_name, strict_report)
    save_matrices(cv_directory, cv_subfolder_name, strict_report[1])

    # Generate detailed hit reports

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







####################################################################################
############################  Common subroutines  ##################################
####################################################################################


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
        raise ValueError("[ERROR] Unsupported metric. Only 'MCC' is supported.")




def cutoffs(true_positives, true_negatives, report_filepath):
    """
    Calculate the optimal cutoff based on MCC, sort proteinIDs into TP, FP, FN, TN.
    Returns:
        optimized_cutoff, trusted_cutoff, noise_cutoff, best_MCC, best_matrix, hit_id_distribution, all_matrices
    """

    sum_TP = len(true_positives)
    TP, FP, FN, TN = 0, 0, 0, true_negatives
    trusted_cutoff = float("-inf")
    noise_cutoff = float("inf")  # start high for minimum search!
    optimized_cutoff = None      # None to signal "not found yet"
    best_MCC = float("-inf")
    best_matrix = []
    true_positives_set = set(true_positives)
    processed_IDs = set()

    # Hit-Score Dictionaries
    hit_report = {}

    with open(report_filepath, 'r') as file:
        for line in file:
            if line.startswith('#'): continue
            fields = line.strip().split()
            if len(fields) < 8:
                print(f"[WARN] Malformed line (skipped): {line.strip()}")
                continue

            hit_id, bitscore = fields[0], float(fields[7])
            
            if hit_id in processed_IDs:
                continue
            
            processed_IDs.add(hit_id)

            
            # Increase counters in confusion matrix
            if hit_id in true_positives_set:
                
                TP += 1
                FN = sum_TP - TP
                if FP == 0:  # first TP (highest bitscore, da absteigend sortiert)
                    trusted_cutoff = bitscore

                true_value = 'TP' # This hit_id is actually a true positive

            else:

                FP += 1
                TN = true_negatives - FP
                
                true_value = 'TN' # This hit_id is actually a true positive
                
            # Calculate MCC
            current_MCC = calculateMetric("MCC", TP, FP, FN, TN)
            #print(f"[DEBUG] TP={TP} FP={FP} FN={FN} TN={TN} MCC={current_MCC:.3f} bitscore={bitscore:.2f} hit_id={hit_id}")

            # Get best MCC
            if current_MCC > best_MCC:
                best_MCC = current_MCC
                optimized_cutoff = bitscore
                best_matrix = [TP, FP, FN, TN]
            

            # Save hit_id for report
            hit_report[hit_id] = {
                "bitscore": bitscore,
                "MCC": current_MCC,
                "TP": TP,
                "FP": FP,
                "FN": FN,
                "TN": TN,
                "assignment": 'None',
                "true_value": true_value,
            }
            
            if TP == sum_TP:
                noise_cutoff = bitscore
                break
            if bitscore < noise_cutoff:
                noise_cutoff = bitscore




    # Sanity-Check: Warn if no optimal cutoff was found
    if optimized_cutoff is None:
        print('[WARN] No optimized cutoff found')
        optimized_cutoff = trusted_cutoff

    # Second pass: assign based on optimized_cutoff
    for hit_id, info in hit_report.items():
        bs = info['bitscore']
        true_val = info['true_value']

        # Assignment at optimized cutoff
        if bs >= optimized_cutoff:
            info['assignment_optimized_score'] = 'TP' if true_val == 'TP' else 'FP'
        else:
            info['assignment_optimized_score'] = 'FN' if true_val == 'TP' else 'TN'

        # Assignment at trusted cutoff
        if bs >= trusted_cutoff:
            info['assignment_trusted_cutoff'] = 'TP' if true_val == 'TP' else 'FP'
        else:
            info['assignment_trusted_cutoff'] = 'FN' if true_val == 'TP' else 'TN'

        # Assignment at noise cutoff
        if bs >= noise_cutoff:
            info['assignment_noise_cutoff'] = 'TP' if true_val == 'TP' else 'FP'
        else:
            info['assignment_noise_cutoff'] = 'FN' if true_val == 'TP' else 'TN'


    return optimized_cutoff, trusted_cutoff, noise_cutoff, hit_report, best_MCC, best_matrix





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
        if safe_compare(mcc, best_mcc, optimized_cutoff, best_optimized_cutoff):
            best_mcc = mcc
            best_optimized_cutoff = optimized_cutoff

    # Ensure no value is infinite or negative
    highest_trusted_cutoff = max(highest_trusted_cutoff, 10)
    lowest_noise_cutoff = max(lowest_noise_cutoff, 10)
    best_optimized_cutoff = best_optimized_cutoff if best_optimized_cutoff and best_optimized_cutoff > 0 else 10
    
    return (sum_matrix, fold_matrices, best_optimized_cutoff, highest_trusted_cutoff, lowest_noise_cutoff)

def safe_compare(mcc, best_mcc, optimized_cutoff, best_optimized_cutoff):
    if any(x is None for x in [mcc, best_mcc, optimized_cutoff, best_optimized_cutoff]):
        return False
    return mcc > best_mcc or (mcc == best_mcc and optimized_cutoff < best_optimized_cutoff)


def count_fasta_headers(fasta_file):
    with open(fasta_file, "r") as f:
        return sum(1 for line in f if line.startswith(">"))


    
