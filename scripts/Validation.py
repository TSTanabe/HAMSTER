#!/usr/bin/python

import os
import subprocess
import sqlite3

from multiprocessing import Pool, Manager


def create_hmms_from_msas(directory, ending="fasta_aln", extension="hmm_cv", cores=1):
    # Get a list of all MSA files with the specified ending in the directory
    msa_files = [f for f in os.listdir(directory) if f.endswith(ending)]

    # Run hmmbuild on each MSA file
    for msa_file in msa_files:
        msa_path = os.path.join(directory, msa_file)
        hmm_file = os.path.join(directory, msa_file.replace(ending, extension))

        # Construct the hmmbuild command
        hmmbuild_cmd = f"hmmbuild --amino --cpu {cores} {hmm_file} {msa_path}"

        # Run the hmmbuild command, redirecting stdout and stderr to null to suppress output
        subprocess.run(hmmbuild_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

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

def initial_process_folds(args_tuple):
    # This process is running in parallel
    # It test the full HMM against a single test target and calculates the confusion matrix
    # It returns the name of the HMM + the MCC + cutoffs + confusion matrix from the validation of the full model
    
    # Prepare the folds and HMMs
    alignment_file, options, index = args_tuple
    
    # Use global variables initialized by the worker
    all_sequence_number = global_all_seq_number
    limit = global_files_number-1
    
    
        

    initial_validation_directory = options.fasta_alignment_directory #Uses the currently unused fasta_alignment_directory
    
    # Generate the cross-validation directory
    hv_subfolder_name = os.path.splitext(os.path.basename(alignment_file))[0]
    hv_directory = os.path.join(initial_validation_directory, hv_subfolder_name)
    os.makedirs(hv_directory, exist_ok=True)
    
    print(f"Initial validation of HMM {index} of {limit} file {hv_subfolder_name}")
    if os.path.exists(os.path.join(hv_directory, f"{hv_subfolder_name}_hv_matrices.txt")):
        return None  # This was already finished return to next HMM validation
    
    
    # Get the HMM for the initial search
    hmm = options.Hidden_markov_model_directory+f"/{hv_subfolder_name}.hmm"
    
    # Get the target file
    name = hv_subfolder_name.split('_')[-1]
    
    # Use the key from the dictionary if it exists; otherwise, use the fallback
    sequence_faa_file = options.targeted_sequence_faa_file_dict.get(name, options.sequence_faa_file)

    # Define the outputfile
    report_name = os.path.splitext(os.path.basename(alignment_file))[0]
    output = f"{hv_directory}/{report_name}.report"
    
    # Search with the HMMs
    HMMsearch(hmm, sequence_faa_file, output, 1)
    
    # Locate the single report
    hv_report = next((os.path.join(hv_directory, filename) 
                      for filename in os.listdir(hv_directory) 
                      if filename.endswith(".report")), None)

    # Check if the report exists
    if hv_report:
        # Make a report for the cutoffs
        model_values = []
        model_hit_distribution = []
        
        # Get the TP seq IDs and save in a set   
        TP_seqIDs = extract_record_ids_from_alignment(alignment_file)
        TN = all_sequence_number - len(TP_seqIDs)
        
        # Calculate strict values for cutoffs
        optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix, hit_id_distribution = cutoffs(TP_seqIDs, TN, hv_report)
        model_hit_distribution.append(hit_id_distribution)  # Store hit distribution
        model_values.append([optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix])

    
        strict_report = calculate_performance_and_sum_matrices(model_values)  # Returns a tuple
        
        # Write down the matrices and cutoffs in the CV subfolder
        with open(os.path.join(hv_directory, f"{hv_subfolder_name}_thresholds.txt"), 'w') as writer:
            writer.write(f"{hv_subfolder_name}\t{strict_report[2]}\t{strict_report[3]}\t{strict_report[4]}\n")
        
        with open(os.path.join(hv_directory, f"{hv_subfolder_name}_MCC.txt"), 'w') as writer:
            fold_matrix = strict_report[1]
            writer.write(f"{hv_subfolder_name}\t{MCC}\n")
    
    
        #More informative output of positive and negative hits
        hit_id_reports(hv_directory, options.database_directory, model_hit_distribution)
    
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
                                value = float(columns[-1].strip())  # Last column (float)
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


def init_validation_worker(all_seq_number, alignment_files, files_number):
    global global_all_seq_number
    global global_alignment_files
    global global_files_number
    
    # Assign values to global variables
    global_all_seq_number = all_seq_number
    global_alignment_files = alignment_files
    global_files_number = files_number
    
    
    
def parallel_cross_validation(options):
    #Main cross validation routine
    print("Initilizing cross-validation")
    
    #Prepare shared data
    all_seq_number = fetch_db_seq_count(options.database_directory)
    alignment_dict = get_alignment_files(options.fasta_output_directory)  #aligned_seqs => unaligned_seqs filepaths for keys and values
    alignment_files = list(alignment_dict.keys())
    files_number = len(alignment_files)
    print("Cutoff and performance validation without folds")
    #Prepare task-specific arguments
    args_list = [(alignment_file, options, index) for index, alignment_file in enumerate(alignment_files)]

    # Make the inital validation of full hmms and assign cutoffs
    with Pool(processes=options.cores, initializer=init_validation_worker, 
              initargs=(all_seq_number, alignment_files, files_number)) as pool:
        
        #starts the workers
        pool.map(initial_process_folds,args_list)

    print("Filter by overall performance and prepare for cross-validation")
    # Filter out the HMMs with a performance below threshold
    initial_HMM_validation_MCC_dict = process_MCC_files(options.fasta_alignment_directory,"_MCC.txt")
    initial_HMM_validation_MCC_dict = filter_dict_by_value(initial_HMM_validation_MCC_dict,options.MCC_threshold)
    
    alignment_files = filter_filepaths(alignment_files, initial_HMM_validation_MCC_dict)
    files_number = len(alignment_files)
    
    args_list = [(alignment_file, options, index) for index, alignment_file in enumerate(alignment_files)]

    print("Cross-validation initilized")
    #Makes the cross validation and assigns the cutoff values 
    with Pool(processes=options.cores, initializer=init_validation_worker, 
              initargs=(all_seq_number, alignment_files, files_number)) as pool:
        
        #starts the workers
        pool.map(process_cross_folds,args_list)
    print(f"Processed all validations.")



    return
    

    
    
def process_cross_folds(args_tuple):
    # This process is running in parallel
    # Prepare the folds and HMMs
    alignment_file, options, index = args_tuple
    
    # Use global variables initialized by the worker
    all_sequence_number = global_all_seq_number
    limit = global_files_number-1
    
    
        

    cross_validation_directory = options.cross_validation_directory
    
    # Generate the cross-validation directory
    cv_subfolder_name = os.path.splitext(os.path.basename(alignment_file))[0]
    cv_directory = os.path.join(cross_validation_directory, cv_subfolder_name)
    
    print(f"Cross-validation of HMM {index} of {limit} file {cv_subfolder_name}")
    if os.path.exists(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_matrices.txt")):
        return None  # This was already finished return to next HMM validation
    
    create_cross_validation_sets(alignment_file, cv_directory)  # make CV folds in subfolder    
    
    # Generate the cross-validation folds for the HMMs
    create_hmms_from_msas(cv_directory, "cv")
    cv_hmms = [os.path.join(cv_directory, filename) for filename in os.listdir(cv_directory) if filename.endswith(".hmm_cv")]
    
    #Get the target file
    name = cv_subfolder_name.split('_')[-1]
    sequence_faa_file = options.targeted_sequence_faa_file_dict.get(name, options.sequence_faa_file)
    
    
    
    # Search with the HMMs
    for hmm in cv_hmms:  # Perform the hmmsearch for each fold
        HMMsearch(hmm, sequence_faa_file, hmm + ".report", 1)
    
    cv_reports = [os.path.join(cv_directory, filename) for filename in os.listdir(cv_directory) if filename.endswith(".report")]
    
    # Make a report for the cutoffs
    model_values = []
    model_hit_distribution = []
    
    # Get the TP seq IDs and save in a set   
    TP_seqIDs = extract_record_ids_from_alignment(alignment_file)
    TN = all_sequence_number - len(TP_seqIDs)
    
    for report in cv_reports:  # Read each of the reports from hmmsearch
        # Calculate strict values for cutoffs, only exact matches with csb and sequence cluster count as TP
        optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix, hit_id_distribution = cutoffs(TP_seqIDs, TN, report)
        model_hit_distribution.append(hit_id_distribution) # each hit_distribution is a list [TP_dict,FP_dict,FN_dict]
        model_values.append([optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix])
    
    strict_report = calculate_performance_and_sum_matrices(model_values)  # Returns a tuple
    
    # Write down the matrices and cutoffs in the CV subfolder
    with open(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_thresholds.txt"), 'w') as writer:
        writer.write(f"{cv_subfolder_name}\t{strict_report[2]}\t{strict_report[3]}\t{strict_report[4]}\n")
    
    with open(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_matrices.txt"), 'w') as writer:
        fold_matrix = strict_report[1]
        writer.write('\n'.join(['\t'.join(map(str, row)) for row in fold_matrix]))
    
    
    #More informative output of positive and negative hits
    hit_id_reports(cv_directory, options.database_directory, model_hit_distribution)
    
    return None

    
    
    
    
    
            


def get_alignment_files(directory):
    alignment_faa_dict = {}

    for file_name in os.listdir(directory):
        if file_name.endswith(".fasta_aln"):
            alignment_file = os.path.join(directory, file_name)
            faa_file = os.path.join(directory, file_name.replace(".fasta_aln", ".faa"))

            # Check if the corresponding .faa file exists
            if os.path.isfile(faa_file):
                alignment_faa_dict[alignment_file] = faa_file

    return alignment_faa_dict #aligned_seqs => unaligned_seqs



def create_cross_validation_sets(alignment_file, output_directory, ending=".cv", num_folds=5):
    # Check if the output directory exists, if not, create it
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Read the alignment file manually
    sequences = []
    record_ids = []
    with open(alignment_file, "r") as f:
        record_id = None
        sequence_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if record_id is not None:
                    sequences.append("".join(sequence_lines))
                    sequence_lines = []
                record_id = line[1:]  # Remove the ">"
                record_ids.append(record_id)
            else:
                sequence_lines.append(line)
        if record_id is not None:
            sequences.append("".join(sequence_lines))

    num_sequences = len(sequences)
    fold_size = num_sequences // num_folds

    # Create cross-validation folds
    for fold in range(num_folds):
        start_idx = fold * fold_size
        end_idx = (fold + 1) * fold_size
        train_sequences = sequences[:start_idx] + sequences[end_idx:]
        train_record_ids = record_ids[:start_idx] + record_ids[end_idx:]
        train_file = os.path.join(output_directory, f"training_data_{fold}{ending}")

        with open(train_file, "w") as train_f:
            for record_id, sequence in zip(train_record_ids, train_sequences):
                train_f.write(f">{record_id}\n{sequence}\n")

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


def count_sequences_in_fasta(options):
    """
    CURRENTLY NOT IN USE
    Counts the number of sequences in a FASTA file. If the specified FASTA file does not exist,
    it will look for 'sequences.faa' in the cv_directory and use it if available.
    
    Args:
        fasta_file (str): Path to the FASTA file.
        cv_directory (str): Path to the directory where 'sequences.faa' may exist.
    
    Returns:
        int: Number of sequences in the FASTA file.
    """
    fasta_file = options.sequence_faa_file
    cv_directory = options.cross_validation_directory

    # Check if the provided fasta_file exists
    if fasta_file is None or not os.path.exists(fasta_file):
        # If the fasta_file does not exist, check if 'sequences.faa' exists in cv_directory
        alt_fasta_file = os.path.join(cv_directory, 'sequences.faa')
        if os.path.exists(alt_fasta_file):
            print(f"FASTA file not found, using '{alt_fasta_file}' instead.")
            fasta_file = alt_fasta_file
            options.sequence_faa_file = fasta_file
        else:
            raise FileNotFoundError(f"Neither {fasta_file} nor 'sequences.faa' found in {cv_directory}")
    
    # Count the number of sequences in the selected FASTA file
    num_sequences = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                num_sequences += 1

    return num_sequences
    
def HMMsearch(library, input_file, output_file, cores=1):

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
#vom höchsten an
#true positives shall be a datastructure with the seqIDs
#true negatives is an int
    sum_TP = len(true_positives) #true_positives is a set, sum_TP is the total number of TP
    true_positive_dict = {}
    false_positive_dict = {}
    false_negative_dict = {key: 0 for key in true_positives}
    
    TP = 0
    FP = 0
    FN = 0
    TN = true_negatives # is an integer
    
    trusted_cutoff = float("-inf")
    noise_cutoff = 10
    optimized_cutoff = 0
    MCC = 0
    matrix = []
    

    with open(report_filepath, 'r') as file:
        for line in file:
            # Skip lines that start with '#'
            if line.startswith('#'):
                continue

            # Split the line by any whitespace, treating consecutive spaces as a single separator
            fields = line.strip().split()

            # Extract the relevant fields (assuming 'hit.id' is in the first column and 'bitscore' in a specific column)
            hit_id = fields[0]
            bitscore = float(fields[13])  # domain for bitscore
            
                        
            if hit_id in true_positives:
                if hit_id in true_positive_dict:
                    continue #skip if hits are double because of domains
                    
                TP = TP + 1
                true_positive_dict[hit_id] = bitscore # collect the positive hit identifiers set
                if FP == 0:
                    trusted_cutoff = bitscore
                FN = sum_TP - TP #wenn tp die hit_ids speichert dann könnte man hier direkt die FN hit IDs rausfiltern über vergleiche von sets
                
                
                if hit_id in false_negative_dict: #track the bitscores, TPs will be removed from here when better MCC is found
                    false_negative_dict[hit_id] = bitscore
            else:
                if hit_id in false_positive_dict:
                    continue #skip if hits are double because of domains
                FP = FP + 1
                false_positive_dict[hit_id] = bitscore # collect the false positive hit identifiers
          
                TN = true_negatives - FP
            
            new_MCC = calculateMetric("MCC",TP,FP,FN,TN)
            #print(f"{TP},{FP},{FN},{TN} with hit_id {hit_id} as next positive hit")
            #print(new_MCC, " > ", MCC)
            
            if new_MCC > MCC:
                MCC = new_MCC
                optimized_cutoff = bitscore
                matrix = [TP,FP,FN,TN]
                
                # Remove all TPs at this point from the false_negative_dict thus scores can be recorded for all hits but only those in the
                # below the threshold will be returned              
                false_negative_dict = {k: v for k, v in false_negative_dict.items() if k not in set(true_positive_dict)}
                
                if not FN == len(false_negative_dict.keys()): # This should never happen
                    print("\n ERROR False negative dict in the cutoff routine")
                    print(FN)
                    print(false_negative_dict)
                    print("\n\n\n")
                    print(true_positives)
                    print("\n\n")
                    print(set(true_positive_dict))
                    print("\n\n\n")


            if TP == sum_TP:
                noise_cutoff = bitscore
                break #MCC there is no better MCC without more TP
            if noise_cutoff > bitscore:
                noise_cutoff = bitscore
    
    false_positive_dict_filtered = {key: value for key, value in false_positive_dict.items() if value >= optimized_cutoff}
    hit_id_distribution = [true_positive_dict.copy(), false_positive_dict_filtered.copy(), false_negative_dict.copy()] #distribution of hit ids at optimum MCC   
    return optimized_cutoff,trusted_cutoff,noise_cutoff,MCC,matrix, hit_id_distribution


def failed_hits(true_positives, report_filepath, cutoff):
    TP = set(true_positives)  # Copy the set to avoid modifying the original
    FP = {}
    
    with open(report_filepath, 'r') as file:
        for line in file:
            # Skip lines that start with '#'
            if line.startswith('#'):
                continue

            # Split the line into fields by whitespace
            fields = line.strip().split()

            # Assuming the ID is in the first column and bitscore in a specific column
            hit_id = fields[0]
            bitscore = float(fields[13])  # Adjust this index based on actual bitscore column

            # If the bitscore is below the cutoff, stop processing
            if bitscore < cutoff:
                return FP, TP

            # Process the hit
            if hit_id in TP:
                TP.discard(hit_id)
            else:
                # Count how many times each hit_id is added as FP
                if hit_id in FP:
                    FP[hit_id] += 1
                else:
                    FP[hit_id] = 1

    return FP, TP  # Return the false positives with their counts and the remaining true positives (FN)



                
#from all fold. decide which is the highest trusted cutoff, the lowest noise cutoff and the optimized cutoff with the best MCC, if equal the lower cutoff. MCC is then the value for the model performance under these circumstances
def calculate_performance_and_sum_matrices(cross_validation_folds):
    highest_trusted_cutoff = float("-inf")
    lowest_noise_cutoff = float("inf")
    best_mcc = float("-inf")
    best_optimized_cutoff = None
    sum_matrix = [0, 0, 0, 0]
    fold_matrices = list()
    for fold_data in cross_validation_folds:
        optimized_cutoff, trusted_cutoff, noise_cutoff, mcc, matrix = fold_data
        
        #List all matrices
        fold_matrices.append(matrix)
        
        # Update highest trusted cutoff
        highest_trusted_cutoff = max(highest_trusted_cutoff, trusted_cutoff)

        # Update lowest noise cutoff
        #print(lowest_noise_cutoff," verglichen mit ", noise_cutoff)
        lowest_noise_cutoff = min(lowest_noise_cutoff, noise_cutoff)

        # Update best MCC and corresponding optimized cutoff
        if mcc > best_mcc or (mcc == best_mcc and optimized_cutoff < best_optimized_cutoff):
            best_mcc = mcc
            best_optimized_cutoff = optimized_cutoff

    # Ensure no value is infinite or negative
    if highest_trusted_cutoff == float("-inf") or highest_trusted_cutoff < 0:
        highest_trusted_cutoff = 10
    if lowest_noise_cutoff == float("inf") or lowest_noise_cutoff < 0:
        lowest_noise_cutoff = 10
    if best_optimized_cutoff is None or best_optimized_cutoff < 0:
        best_optimized_cutoff = 10

    
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
    key_counts = {
        'TP': {},
        'FP': {},
        'FN': {}
    }
    
    # lowest hit score for each proteinID
    hit_scores = {
        'TP': {},
        'FP': {},
        'FN': {}
    }

    # Iterate over the list of lists
    for inner_list in data:
        # Process TP, FP, FN dicts in a loop to avoid redundancy
        for i, dict_label in enumerate(['TP', 'FP', 'FN']):
            hit_dictionary = inner_list[i]
            #print(f"\t{dict_label}\n")
            # Update the overall count for TP, FP, or FN
            for key, count in hit_dictionary.items():
                
                if key in key_counts[dict_label]:
                    key_counts[dict_label][key] += 1
                    #print(f"+1 {key}")
                else:
                    key_counts[dict_label][key] = 1
                    #print(f"Initialized {key}")
            #print(key_counts)            
            #Now organize the hit score values
            for hit, score in hit_dictionary.items():
                if hit in hit_scores[dict_label]:
                    if hit_scores[dict_label][hit] > score:
                        hit_scores[dict_label][hit] = score
                else:
                    hit_scores[dict_label][hit] = score
                
    

    #print(hit_scores)       
    # Fetch and print gene vicinity information for the current dictionary
    for dict_label in ['TP', 'FP', 'FN']:
        gene_vicinity_dictionary = fetch_neighbouring_genes_with_domains(database, list(key_counts[dict_label].keys()))
        
        # Prepare data for the report
        sorted_data = []
                
        # Sort the keys based on their score (from score_dict) and generate the report rows
        score_dict = hit_scores[dict_label]
        for key in sorted(score_dict.keys(), key=lambda k: score_dict[k], reverse=True):
            score = score_dict[key]  # The score from the inner_list[i] (TP_dict, FP_dict, FN_dict)
            count = key_counts[dict_label].get(key, 0)  # The count of hits for this key
            neighborhood = gene_vicinity_dictionary.get(key, "No data")  # The genomic neighborhood
                    
            # Append a tuple (key, score, count, neighborhood) to the sorted data
            sorted_data.append((key, score, count, neighborhood))
                
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
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Fetch all proteins and their domains in the same clusters as the input protein_ids
        query = f"""
        SELECT 
            p.proteinID, p.clusterID, p.genomeID, p.start, d.domain
        FROM Proteins p
        LEFT JOIN Domains d ON p.proteinID = d.proteinID
        WHERE p.proteinID IN ({','.join('?' * len(protein_ids))})
           OR p.clusterID IN (
                SELECT DISTINCT clusterID FROM Proteins WHERE proteinID IN ({','.join('?' * len(protein_ids))})
           )
        ORDER BY p.clusterID, p.start
        """
        
        # Execute the query with the protein_ids list
        cur.execute(query, protein_ids + protein_ids)  # Pass protein_ids twice for both conditions
        protein_results = cur.fetchall()
    
    # Organize results by cluster and by proteinID
    cluster_dict = {}
    protein_cluster_map = {}  # Map each protein to its cluster

    for protein_id, cluster_id, genome_id, start, domain in protein_results:
        if cluster_id is not None:  # Only process if a cluster exists
            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = []
            # If domain is None, set it to 'no_domain'
            domain_str = domain if domain else "no_domain"
            cluster_dict[cluster_id].append((protein_id, start, f"{domain_str}_{protein_id}"))
            protein_cluster_map[protein_id] = cluster_id
        else:
            # Handle proteins without clusters separately in the result later
            protein_cluster_map[protein_id] = None

    # Prepare the final dictionary with proteinID as the key and neighboring genes sorted by start
    neighbors_dict = {}
    for protein_id in protein_ids:
        cluster_id = protein_cluster_map.get(protein_id)
        
        # Handle proteins with a cluster
        if cluster_id:
            # Find all proteins in the same cluster
            cluster_proteins = cluster_dict[cluster_id]
            # Sort proteins by start position
            sorted_proteins = sorted(cluster_proteins, key=lambda x: x[1])
            # Create the neighbor list including the current protein's own entry
            neighbors = [protein for pid, _, protein in sorted_proteins]
            neighbors_dict[protein_id] = neighbors
        else:
            # Handle proteins without a cluster (just return the protein itself)
            neighbors_dict[protein_id] = [f"no_domain_{protein_id}"]

    return neighbors_dict





#######################################################################
##################    Count TP + TN in DB      ########################
#######################################################################



def fetch_db_seq_count(database):
    """
    Fetches the count of unique protein sequences from the SQLite database.

    Args:
        database (str): Path to the SQLite database.

    Returns:
        int: The number of unique sequences in the database (excluding 'QUERY' genomeID).
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        
        # Query to count distinct proteins excluding genomeID 'QUERY'
        query = """
                SELECT COUNT(DISTINCT Proteins.proteinID)
                FROM Proteins
                WHERE NOT Proteins.genomeID = ?
            """
        cur.execute(query, ('QUERY',))
        seq_count = cur.fetchone()[0]

    return seq_count











    
    
    
    
    
    
    
    
    
    
    
