#!/usr/bin/python

import os
import subprocess
from multiprocessing import Pool, Manager
from Bio import SeqIO


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


def parallel_cross_validation(options):    
    All_seq_number = count_sequences_in_fasta(options.sequence_faa_file)
    
    alignment_dict = get_alignment_files(options.fasta_output_directory)  #aligned_seqs => unaligned_seqs filepaths for keys and values
    alignment_files = list(alignment_dict.keys())
    files_number = len(alignment_files)
    
    
    #prepare the arguments
    args_list = []
    for index,alignment_file in enumerate(alignment_files):
        
        args_list.append([alignment_file, All_seq_number, options, index, files_number])
        
    
    #Makes the cross validation and assigns the cutoff values 
    with Pool(processes=options.cores) as pool:
        
        #starts the workers
        pool.map(process_cross_folds,args_list)
        

    print(f"Processed all validations.")
    return
    

def process_cross_folds(args_tuple):
    # This process is running in parallel
    # Prepare the folds and HMMs
    
    alignment_file, all_sequence_number, options, index, limit = args_tuple
    sequence_faa_file = options.sequence_faa_file
    cross_validation_directory = options.cross_validation_directory
    
    
    # Generate the cross-validation directory
    cv_subfolder_name = os.path.splitext(os.path.basename(alignment_file))[0]
    cv_directory = os.path.join(cross_validation_directory, cv_subfolder_name)
    
    create_cross_validation_sets(alignment_file, cv_directory)  # make CV folds in subfolder
    
    if os.path.exists(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_matrices.txt")):
        return None  # This was already finished
    
    print(f"Validation of HMM {index} of {limit} file {cv_subfolder_name}")
    
    # Get the TP seq IDs and save in a set   
    TP_seqIDs = extract_record_ids_from_alignment(alignment_file)
    TN = all_sequence_number - len(TP_seqIDs)
    
    # Generate the cross-validation folds for the HMMs
    create_hmms_from_msas(cv_directory, "cv")
    cv_hmms = [os.path.join(cv_directory, filename) for filename in os.listdir(cv_directory) if filename.endswith(".hmm_cv")]
    
    # Search with the HMMs
    for hmm in cv_hmms:  # Perform the hmmsearch for each fold
        HMMsearch(hmm, sequence_faa_file, hmm + ".report", 1)
    
    cv_reports = [os.path.join(cv_directory, filename) for filename in os.listdir(cv_directory) if filename.endswith(".report")]
    
    # Make a report for the cutoffs
    model_values = []
    for report in cv_reports:  # Read each of the reports from hmmsearch
        # Calculate strict values for cutoffs, only exact matches with csb and sequence cluster count as TP
        optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix = cutoffs(TP_seqIDs, TN, report)
        model_values.append([optimized_cutoff, trusted_cutoff, noise_cutoff, MCC, matrix])
    
    strict_report = calculate_performance_and_sum_matrices(model_values)  # Returns a tuple
    
    # Write down the matrices and cutoffs in the CV subfolder
    with open(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_thresholds.txt"), 'w') as writer:
        writer.write(f"{cv_subfolder_name}\t{strict_report[2]}\t{strict_report[3]}\t{strict_report[4]}\n")
    
    with open(os.path.join(cv_directory, f"{cv_subfolder_name}_cv_matrices.txt"), 'w') as writer:
        fold_matrix = strict_report[1]
        writer.write('\n'.join(['\t'.join(map(str, row)) for row in fold_matrix]))
    
    # Makes a list of all FP and FN hits, forms an average over all CV folds, and stores them in the cv_directory for later analysis
    sum_FN = set()
    sum_FP = dict()
    for report in cv_reports:
        FP, FN = failed_hits(TP_seqIDs, report, strict_report[2])
        sum_FN.update(FN)
        sum_FP.update(FP)
    
    with open(os.path.join(cv_directory, f"sum_FN_{cv_subfolder_name}.txt"), "w") as file_FN:
        for protID,hits in sum_FN.items():
            it = protID.replace("**", "\t")
            file_FN.write(f"{it}\t{hits}\n")
    
    with open(os.path.join(cv_directory, f"sum_FP_{cv_subfolder_name}.txt"), "w") as file_FP:
        for item in sum_FP:
            it = item.replace("**", "\t")
            file_FP.write(f"{it}\n")
    
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

    # Read the alignment file
    sequences = []
    record_ids = []
    for record in SeqIO.parse(alignment_file, "fasta"):
        sequences.append(str(record.seq))
        record_ids.append(record.id)

    num_sequences = len(sequences)
    fold_size = num_sequences // num_folds

    # Create custom cross-validation folds
    for fold in range(num_folds):
        start_idx = fold * fold_size
        end_idx = (fold + 1) * fold_size
        train_sequences = sequences[:start_idx] + sequences[end_idx:]
        train_record_ids = record_ids[:start_idx] + record_ids[end_idx:]
        train_file = os.path.join(output_directory, f"training_data_{fold}{ending}")

        with open(train_file, "w") as train_f:
            for record_id, sequence in zip(train_record_ids, train_sequences):
                train_f.write(f">{record_id}\n{sequence}\n")
                

def deprecated_concat_files_with_extension(directory, extension, output_file):
    file_list = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(extension)]
    if not file_list:
        print(f"No files with the '{extension}' extension found in {directory}")
        return None
    
    with open(output_file, 'w') as output:
        
        for filename in file_list:
            #file_path = os.path.join(directory, filename)
            with open(filename, 'r') as f:
                output.write(f.read())

    return output_file

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


def count_sequences_in_fasta(fasta_file):
    num_sequences = 0
    
    # Open the FASTA file and iterate over each line
    with open(fasta_file, 'r') as file:
        for line in file:
            # Check if the line starts with '>'
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
    TP = 0
    sum_TP = len(true_positives) #true_positives is a set
    FP = 0
    FN = 0
    TN = true_negatives
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
                TP = TP + 1
                if FP == 0:
                    trusted_cutoff = bitscore
            else:
                FP = FP + 1
                
            FN = sum_TP - TP
            TN = true_negatives - FP
            
            new_MCC = calculateMetric("MCC",TP,FP,FN,TN)
            #print(f"{TP},{FP},{FN},{TN}")
            #print(new_MCC, " > ", MCC)
            
            if new_MCC > MCC:
                MCC = new_MCC
                optimized_cutoff = bitscore
                matrix = [TP,FP,FN,TN]
            if TP == sum_TP:
                noise_cutoff = bitscore
                break #MCC there is no better MCC without more TP
            if noise_cutoff > bitscore:
                noise_cutoff = bitscore
   
    return optimized_cutoff,trusted_cutoff,noise_cutoff,MCC,matrix
    

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


    
    return (sum_matrix, fold_matrices, best_optimized_cutoff, highest_trusted_cutoff, lowest_noise_cutoff)

 


    
    
    




























    
    
    
    
    
    
    
    
    
    
    
