#!/usr/bin/python

import os
import shutil
import subprocess
from . import myUtil
from sklearn.model_selection import KFold
from Bio import AlignIO
from Bio import SearchIO
from Bio import SeqIO


class ReportObject:
    def __init__(self, matrix = None, fold_matrices = None, optimized_cutoff = None, trusted_cutoff = None, noise_cutoff = None):
        self.matrix = matrix
        self.fold_matrices = fold_matrices
        self.optimized_cutoff = optimized_cutoff
        self.trusted_cutoff = trusted_cutoff 
        self.noise_cutoff = noise_cutoff #TODO if infinite than 10
        


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
                
def concatenate_fasta_files(fasta_dict, output_directory): #TODO delete? 021223 keine anwendung dafür gefunden
    
    fasta_files = list(fasta_dict.values()) # values are filepaths to fasta files
    output_file = os.path.join(output_directory, f"Sequence.faa") # should be a persistent dir

    # Concatenate fasta files using the operating system command
    try:
        if len(fasta_files) > 0:
            if len(fasta_files) == 1:
                # If there is only one file, use shutil.copyfile directly
                shutil.copyfile(fasta_files[0], output_file)
            else:
                # Use operating system command to concatenate multiple files
                if os.name == 'posix':  # Linux/MacOS
                    subprocess.run(["cat"] + fasta_files, stdout=open(output_file, 'wb'))
                elif os.name == 'nt':  # Windows
                    subprocess.run(["type"] + fasta_files, stdout=open(output_file, 'wb'))

        else:
            print("No fasta files to concatenate.")

    except Exception as e:
        print("Error:", e)
    
    return output_file

def concat_files_with_extension(directory, extension, output_file):
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
    record_ids = {}

    # Read the alignment file and get the record IDs
    alignment = AlignIO.read(alignment_file, "fasta")
    for record in alignment:
        record_ids[record.id] = record.seq

    return record_ids

# build hmm and search the true positive file parse the report file and compute the confusion matrix for each threshold. => makes roc curve. find optimal threshold and write down the matrix



def create_hmms_from_msas(directory, ending="cv", extension="hmm_cv", cores=1):
    # Get a list of all MSA files with .cv extension in the directory
    
    msa_files = [f for f in os.listdir(directory) if f.endswith(ending)]

    # Run hmmbuild on each MSA file
    for msa_file in msa_files:
        msa_path = os.path.join(directory, msa_file)
        hmm_file = os.path.join(directory, msa_file.replace(ending, extension))

        # Run the hmmbuild command using subprocess
        hmmbuild_cmd = f"hmmbuild --amino --cpu {cores} {hmm_file} {msa_path}"
        #subprocess.run(hmmbuild_cmd, shell=True)
        with open(os.devnull, 'w') as null_output:
            subprocess.run(hmmbuild_cmd, shell=True, stdout=null_output, stderr=subprocess.STDOUT)

#dann sind jetzt im directory die .cv files und die dazugehörigen .msa files


def HMMsearch(library, input_file, output_file, cores=1):

    # Construct the hmmsearch command
    hmmsearch_cmd = f"hmmsearch -E 0.0001 --cpu {cores} {library} {input_file} > {output_file}"

    # Run the hmmsearch command using subprocess
    status = subprocess.call(hmmsearch_cmd, shell=True)

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


def Cutoffs(true_positives,true_negatives,report_filepath,relaxed=0):
#vom höchsten an
#true positives shall be a datastructure with the seqIDs
#true negatives is an int
    TP = 0
    sum_TP = len(true_positives)
    FP = 0
    FN = 0
    TN = true_negatives - sum_TP
    trusted_cutoff = float("-inf")
    noise_cutoff = float("inf")
    optimized_cutoff = 0
    MCC = 0
    matrix = []
    
    parts = true_positives[0].split('_')
    current_seq_cluster_class = '_'.join(parts[-2:])
    
    for hmmer_qresult in SearchIO.parse(report_filepath,"hmmer3-text"):
        
        for hit in hmmer_qresult:
            
            if hit.id in true_positives:
                TP = TP + 1
                if FP == 0:
                    trusted_cutoff = hit.bitscore
            elif relaxed:
                #check for sequence cluster similarity
                parts = hit.id.split('_')
                seq_cluster_class = '_'.join(parts[-2:])
                if not seq_cluster_class == current_seq_cluster_class:
                    FP = FP + 1
            else:
                FP = FP + 1
                
            FN = sum_TP - TP
            TN = true_negatives - FP
            
            new_MCC = calculateMetric("MCC",TP,FP,FN,TN)
            #print(f"{TP},{FP},{FN},{TN}")
            #print(new_MCC, " > ", MCC)
            if new_MCC > MCC:
                MCC = new_MCC
                optimized_cutoff = hit.bitscore
                matrix = [TP,FP,FN,TN]
            if TP == sum_TP:
                noise_cutoff = hit.bitscore
                break #MCC there is no better MCC without more TP
            if noise_cutoff > hit.bitscore:
                noise_cutoff = hit.bitscore
   
    return optimized_cutoff,trusted_cutoff,noise_cutoff,MCC,matrix

                
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

        # Sum up the matrix values
        for i in range(len(matrix)):
            sum_matrix[i] += matrix[i]

    
    reportObj = ReportObject(sum_matrix, fold_matrices, best_optimized_cutoff, highest_trusted_cutoff, lowest_noise_cutoff)
    print(sum_matrix, fold_matrices, best_optimized_cutoff, highest_trusted_cutoff, lowest_noise_cutoff) #TODO Delete
    return reportObj

 
def count_sequences_in_fasta(fasta_file):
    num_sequences = 0
    
    # Iterate over the records in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        num_sequences += 1
    
    return num_sequences




def process_cross_folds(args_tuple):
    #prepare the folds and HMMs
    alignment_file, TN, sequence_faa_file,options = args_tuple
    f = myUtil.get_filename_without_extension(alignment_file) #for the name of the subfolder
    create_cross_validation_sets(alignment_file,options.Cross_validation_directory+f"/{f}") #make CV folds in subfolder
    TP_seqIDs = extract_record_ids_from_alignment(alignment_file) # get the TP seq IDs
    create_hmms_from_msas(options.Cross_validation_directory+f"/{f}","cv") # create the cross validatin fold
        
    cv_hmms = [os.path.join(options.Cross_validation_directory+f"/{f}", filename) for filename in os.listdir(options.Cross_validation_directory+f"/{f}") if filename.endswith(f".hmm_cv")]
    
    #search the HMMs
    for hmm in cv_hmms: # perform the hmmsearch for each fold
        HMMsearch(hmm,sequence_faa_file,hmm+".report",1)
    cv_reports = [os.path.join(options.Cross_validation_directory+f"/{f}", filename) for filename in os.listdir(options.Cross_validation_directory+f"/{f}") if filename.endswith(f".report")]
    
    #make a report for the cutoffs
    model_values = []
    relaxed_model_values = []
    for report in cv_reports: #read each of the reports from hmmsearch
        #Calculate strict values for cutoffs only exact matches with csb and sequence cluster count as TP
        optimized_cutoff,trusted_cutoff,noise_cutoff,MCC,matrix = Cutoffs(TP_seqIDs,TN,report)
        model_values.append([optimized_cutoff,trusted_cutoff,noise_cutoff,MCC,matrix])
 
        #Calculate relaxed values for cutofs matches with same sequence cluster but different csb do not count as FP
        optimized_cutoff,trusted_cutoff,noise_cutoff,MCC,matrix = Cutoffs(TP_seqIDs,TN,report,1)
        relaxed_model_values.append([optimized_cutoff,trusted_cutoff,noise_cutoff,MCC,matrix])
        
    report = calculate_performance_and_sum_matrices(model_values) #returns a report object

    relaxed_report = calculate_performance_and_sum_matrices(relaxed_model_values) #returns a report object


    
    return f, (report,relaxed_report)
    
    
    
    
    
    
    
    
def convert_HMMreport_to_sets(Filepath,Thresholds,extension=".tmp"):
    output_file_path = os.path.splitext(Filepath)[0] + extension   
    with open(output_file_path, 'a') as output_file:
    
    
        for hmmer_qresult in SearchIO.parse(Filepath,"hmmer3-text"):
            query = hmmer_qresult.id    # Name of HMM without description
            true_positives = set()
            false_positives = set()
            false_negatives = set()
            
            if query in Thresholds:
                threshold = Thresholds[query] # specific threshold
            
            
            for hit in hmmer_qresult:
                if threshold < hit.bitscore:
                    if hit.desc == query:
                        # do something when bitscore is above threshold and desc matches query
                        true_positives.append(hit.bitscore," ",hit.id)
                    else:
                        # do something when bitscore is above threshold and desc doesn't match query
                        false_positives.append(hit.bitscore," ",hit.id)
                elif threshold > hit.bitscore:
                    false_negatives.append(hit.bitscore," ",hit.id)

            #make a report of this as a tmp file
            combined_list = list(zip(true_positives, false_positives, false_negatives))

            # Write the combined sets to a text file for later use
            for row in combined_list:
                output_file.write(threshold," ",query,"\t")
                output_file.write("\t".join(map(str, row)) + "\n")

    return output_file_path
    
    
    
    
    
    
    
    
    
    
