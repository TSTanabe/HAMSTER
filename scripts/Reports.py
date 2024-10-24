#!/usr/bin/python

import os
import shutil

class ReportObject:
    def __init__(self, matrix = None, fold_matrices = None, optimized_cutoff = None, trusted_cutoff = None, noise_cutoff = None):
        self.matrix = matrix # is an array []
        self.fold_matrices = fold_matrices # is a list
        self.optimized_cutoff = optimized_cutoff
        self.trusted_cutoff = trusted_cutoff 
        self.noise_cutoff = noise_cutoff
        
        self.accuracy = 0
        self.f1score = 0
        self.mcc = 0        

def move_HMMs(input_folder,output_folder,file_extension):
    print(f"Saving HMMs in the  directory: {output_folder}")
    for datafile in os.listdir(input_folder):
        if datafile.endswith(file_extension):
            source = os.path.join(input_folder, datafile)
            target = os.path.join(output_folder, datafile)
            shutil.move(source,target)

def concatenate_cv_cutoff_files(directory, file_extension, output_file):
    with open(output_file, 'w') as outfile:
        for root, dirs, files in os.walk(directory):
            for file_path in files:
                if file_path.endswith(file_extension):
                    file_path = os.path.join(root, file_path)
                    with open(file_path, 'r') as infile:
                        outfile.write(infile.read())
    return output_file
            
def create_cutoff_file(options,cutoff_dict,directory,filename = "/thresholds.txt"):
    file_path = directory + filename
    with open(file_path, 'w') as file:
        sorted_keys = sorted(cutoff_dict.keys())
        for key in sorted_keys:
            report = cutoff_dict[key]
            file.write(key+"\t"+str(report.optimized_cutoff)+"\t"+str(report.trusted_cutoff)+"\t"+str(report.noise_cutoff)+"\n")        
    
    return file_path            

def parse_matrices_to_report(directory,file_extension):
    report_dict = dict()  # List to store folder names

    for root, dirs, files in os.walk(directory):
        for file_path in files:
            file_path = os.path.join(root, file_path)
            if file_path.endswith(file_extension):
                folder_name = os.path.basename(root)
                    
                data = []  # List to store lists of four values
                column_sums = []
                with open(file_path, 'r') as infile:

                    for line in infile:     #Make the fold matrices
                        # Split each line into a list of four values and strip the newline
                        values = [int(value.strip()) for value in line.strip().split('\t')]
                        
                        # Ensure that there are exactly four values in the line
                        if len(values) == 4:
                            data.append(values)
                        else:
                            print(f"Warning: Skipping line with incorrect number of values: {line}")
                        
                            #sum up the fold matrices
                    columns = list(zip(*data))
                    column_sums = [sum(map(int, column)) for column in columns]     # Calculate the sum for each column
                report = ReportObject(column_sums,data)
                report_dict[folder_name] = report
                
                       
    return report_dict # dictionary key=> model, value= report object
            
def create_performance_file(options,performance_dict,directory,filename = "/performance.txt"):
    file_path = directory + filename
    
    for key,report in performance_dict.items():
        report.accuracy = calculate_balanced_accuracy(report.fold_matrices)
        report.f1score = calculate_f1_score(report.matrix)
        report.mcc = calculate_mcc(report.matrix)        
    
    with open(file_path, 'w') as file:
        file.write("HMM\tbalanced accuracy\tF1 score\tMCC\tmatrix[TP,FP,FN,TN]\t\n")
        sorted_items = sorted(performance_dict.items(), key=lambda x: (x[0].split('_')[-1], x[1].mcc), reverse=True)
  
        for key, report in sorted_items:
            file.write(key+"\t"+str(report.accuracy)+"\t"+str(report.f1score)+"\t"+str(report.mcc)+"\t"+str(report.matrix)+"\n")        

    return file_path
            
            
def calculate_balanced_accuracy(confusion_matrices):
    balanced_accuracies = []
    
    for matrix in confusion_matrices:
        TP, FP, FN, TN = matrix
        sensitivity = TP / (TP + FN)
        specificity = TN / (TN + FP)
        balanced_accuracy = (sensitivity + specificity) / 2
        balanced_accuracies.append(balanced_accuracy)
    
    # Calculate the average balanced accuracy
    average_balanced_accuracy = sum(balanced_accuracies) / len(balanced_accuracies)
    
    return average_balanced_accuracy

def calculate_f1_score(confusion_matrix):
    TP, FP, FN, TN = confusion_matrix
    
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    
    f1_score = 2 * (precision * recall) / (precision + recall)
    
    return f1_score

def calculate_mcc(confusion_matrix):
    TP, FP, FN, TN = confusion_matrix

    numerator = (TP * TN) - (FP * FN)
    denominator = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5

    if denominator == 0:
        mcc = 0  # To handle cases where the denominator is zero
    else:
        mcc = numerator / denominator

    return mcc






def read_grouped_csb(filepath):
    # Initialize an empty dictionary to store the data
    data = {}

    # Read the tab-separated file
    with open(filepath, 'r') as f:
        for line in f:
            # Split each line by tabs
            elements = line.strip().split('\t')
            if elements:
                key = elements[0]
                values = elements[1:]
                data[key] = values

    return data




def csbs_report(input_string,csbs):
    key = "csb_"+str(input_string.split('csb', 1))
    result = csbs.get(key, "")
    return str(result)


def print_report(filepath):

    input_file_path = 'sets_output.txt'
    output_file_path = 'output_result.txt'

    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:

        for line in input_file:

            
            values = line.strip().split('\t')
            set1, set2, set3 = map(set, values[1:3])

            # Iterate through elements in each set
            for current_set in [set2,set3]:
                # Initialize a defaultdict to store the count of last parts
                last_part_count = defaultdict(int)
                # Initialize a dictionary to store sets with last parts as keys
                sets_with_last_parts = dict()

                for element in current_set:
                    # Split the element by '**'
                    split_parts = element.split('**')

                    # Get the last part
                    last_part = split_parts[-1]

                    # Update the count of last parts
                    last_part_count[last_part] += 1

                    # Store the leading part in an anonymous set in the dictionary
                    if last_part not in sets_with_last_parts:
                        sets_with_last_parts[last_part] = {split_parts[0]}
                    else:
                        sets_with_last_parts[last_part].add(split_parts[0])

                # Write the results to the output file to the output file
                output_file.write("FALSE POSITIVE HITS:\n")
                for last_part, count in last_part_count.items():
                    output_file.write(f"\t{last_part} was found {count} times\n")
                    output_file.write(f"protein IDs\t {last_part}: {leading_parts_set}")


