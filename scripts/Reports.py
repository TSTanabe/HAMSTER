#!/usr/bin/python

import os
import shutil

        

def move_HMMs(input_folder,output_folder,file_extension):

    for datafile in os.listdir(input_folder):
        if datafile.endswith(file_extension):
            source = os.path.join(input_folder, datafile)
            target = os.path.join(output_folder, datafile)
            shutil.move(source,target)
            
def create_cutoff_file(options,cutoff_dict,directory,filename = "/thresholds.txt"):
    file_path = directory + filename
    with open(file_path, 'w') as file:
        sorted_keys = sorted(cutoff_dict.keys())
        for key in sorted_keys:
            report = cutoff_dict[key]
            file.write(key+"\t"+str(report.optimized_cutoff)+"\t"+str(report.trusted_cutoff)+"\t"+str(report.noise_cutoff)+"\n")        
    
    return file_path            
            
def create_performance_file(options,performance_dict,directory,filename = "/performance.txt"):
    file_path = directory + filename
    with open(file_path, 'w') as file:
        file.write("HMM\tbalanced accuracy\tF1 score\tMCC\tmatrix\tCsb\n")
        sorted_keys = sorted(performance_dict.keys())
        
        csbs = read_grouped_csb(options.csb_output_file)
        
        for key in sorted_keys:
            report = performance_dict[key]
            csb = csbs_report(key,csbs)
            
            #welche werte brauchen wir?
            #balanced accuracy
            accuracy = calculate_balanced_accuracy(report.fold_matrices)
            #F1 score
            f1score = calculate_f1_score(report.matrix)
            #MCC
            mcc = calculate_mcc(report.matrix)
            file.write(key+"\t"+str(accuracy)+"\t"+str(f1score)+"\t"+str(mcc)+"\t"+str(report.matrix)+"\t"+csb+"\n")        

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
#    # Print the results
#    print("Last Part Count:")
#    for last_part, count in last_part_count.items():
#        print(f"{last_part}: {count}")

#    print("\nSets with Last Parts:")
#    for last_part, leading_parts_set in sets_with_last_parts.items():


#with :
#output_file.write(line)















#make a new textfile to write to in the same location as the filepath

#for each line in the filepath
#split the line by the tab separation
#columns 1 2 and 3 are printed sets revert them into sets
#for each set iterate the elements
#each element split at the **
#the last part of the split may occur in other elements of the set. count how often this last part occurs in the intire set at hand and store the string of the last part in a hast as key and the number of occurence as value
#for the leading part from the split, store this string in a set. the set shall be an anonymous set that is itself the value of a dict. the corresponding key shall be the last part of the split


