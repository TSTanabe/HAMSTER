#!/usr/bin/python

import os
import shutil
import sqlite3
import pickle
from collections import defaultdict
from multiprocessing import Pool, cpu_count


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


def concat_and_sort_files(input_dir, file_extension, output_dir, output_filename="concatenated_sorted.tsv"):
    """
    Recursively concatenates all files with the given extension in a directory and its subdirectories,
    sorts the resulting file by the first and second columns, and saves it to the output directory.

    Parameters:
        input_dir (str): Directory to search for files.
        file_extension (str): File extension to filter (e.g., ".txt").
        output_dir (str): Directory to save the concatenated sorted file.
        output_filename (str): Name of the output file. Default is 'concatenated_sorted.tsv'.

    Returns:
        str: Path to the resulting output file.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Collect all matching files recursively
    files_to_concat = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(file_extension):
                files_to_concat.append(os.path.join(root, file))

    if not files_to_concat:
        raise ValueError(f"No files with extension '{file_extension}' found in directory: {input_dir}")

    # Read all lines from the matching files
    all_lines = []
    for file in files_to_concat:
        with open(file, 'r') as f:
            all_lines.extend(f.readlines())

    # Sort lines by the first and second columns
    sorted_lines = sorted(all_lines, key=lambda line: tuple(line.split("\t")[:2]))

    # Write the sorted lines to the output file
    output_file_path = os.path.join(output_dir, output_filename)
    with open(output_file_path, 'w') as f:
        f.writelines(sorted_lines)

    return output_file_path
    
    
    












#######################################################################
##################    Hit report routines      ########################
#######################################################################


def load_and_process_hit_distributions(directory, database_directory):
    """
    Recursively searches for serialized hit_id_distribution objects in a given directory, 
    loads them, and executes `hit_id_reports` in the same directory where the `.pkl` file was found.

    Args:
    - directory (str): The root directory where serialized `.pkl` files are stored.
    - database_directory (str): Path to the database directory for `hit_id_reports`.

    Returns:
    - hit_distributions (dict): A dictionary where keys are HMM subfolder names 
                                and values are their corresponding hit_id_distribution data.
    """

    hit_distributions = {}

    # Walk through all subdirectories and files in the given directory
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith("_hit_id_distribution.pkl"):
                file_path = os.path.join(root, file)
                
                # Extract the exact directory where the `.pkl` file is located
                hmm_directory = os.path.dirname(file_path)  

                # Extract the HMM subfolder name from the filename
                hv_subfolder_name = os.path.basename(file).replace("_hit_id_distribution.pkl", "")

                # Load the serialized data
                with open(file_path, 'rb') as f:
                    try:
                        hit_id_distribution = pickle.load(f)
                        hit_distributions[hv_subfolder_name] = hit_id_distribution
                        print(f"Successfully loaded {hv_subfolder_name} from {file_path}")

                        # Execute hit_id_reports using the exact directory where the `.pkl` was found
                        hit_id_reports(hmm_directory, database_directory, [hit_id_distribution])
                        print(f"Processed hit_id_reports for {hv_subfolder_name} in {hmm_directory}")

                    except Exception as e:
                        print(f"Error: Failed to load {file_path}: {e}")

    return hit_distributions





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
    Generates a hit report and includes neighborhood data.
    Also saves sorted domain compositions separately.
    """
    key_counts = {label: defaultdict(int) for label in ['TP', 'FP', 'FN']}
    hit_scores = {label: defaultdict(lambda: float("inf")) for label in ['TP', 'FP', 'FN']}

    for tp_dict, fp_dict, fn_dict in data:
        for dict_label, hit_dict in zip(['TP', 'FP', 'FN'], [tp_dict, fp_dict, fn_dict]):
            for hit, score in hit_dict.items():
                key_counts[dict_label][hit] += 1
                hit_scores[dict_label][hit] = min(hit_scores[dict_label][hit], score)

    # Fetch gene neighborhood information and sorted compositions
    gene_vicinity, sorted_compositions = fetch_neighbouring_genes_with_domains(database, list(set.union(*[set(d.keys()) for d in key_counts.values()])))

    # Write sorted compositions report
    write_sorted_compositions_report(cv_directory, sorted_compositions)

    for dict_label in ['TP', 'FP', 'FN']:
        score_dict = hit_scores[dict_label]
        sorted_data = sorted(
            ((key, score, key_counts[dict_label][key], gene_vicinity.get(key, [["No data"]])) 
             for key, score in score_dict.items()),
            key=lambda x: x[1], reverse=True
        )

        write_report_to_file(cv_directory, dict_label, sorted_data)
        write_report_to_iTolBinary(cv_directory, dict_label, sorted_data)

    



def fetch_neighbouring_genes_with_domains_deprecated(database, protein_ids):
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




##############################################################################
######################## Fetch neighbourhood routines ########################
##############################################################################

def fetch_cluster_ids(db_path, protein_ids):
    """
    Fetch cluster IDs for a batch of protein IDs.
    """
    with sqlite3.connect(db_path) as con:
        cur = con.cursor()
        cur.execute(f"SELECT proteinID, clusterID FROM Proteins WHERE proteinID IN ({','.join(['?']*len(protein_ids))})", tuple(protein_ids))
        return {row[0]: row[1] for row in cur.fetchall()}

def fetch_cluster_proteins(db_path, cluster_ids):
    """
    Fetch all proteins belonging to a set of cluster IDs.
    """
    with sqlite3.connect(db_path) as con:
        cur = con.cursor()
        cur.execute(f"SELECT proteinID, clusterID FROM Proteins WHERE clusterID IN ({','.join(['?']*len(cluster_ids))})", tuple(cluster_ids))
        cluster_to_proteins = defaultdict(set)
        for protein_id, cluster_id in cur.fetchall():
            cluster_to_proteins[cluster_id].add(protein_id)
        return cluster_to_proteins

def fetch_protein_domains(db_path, protein_ids):
    """
    Fetch domains for a batch of proteins.
    """
    with sqlite3.connect(db_path) as con:
        cur = con.cursor()
        cur.execute(f"SELECT proteinID, domain FROM Domains WHERE proteinID IN ({','.join(['?']*len(protein_ids))})", tuple(protein_ids))
        protein_to_domains = defaultdict(list)
        for protein_id, domain in cur.fetchall():
            protein_to_domains[protein_id].append(domain)
        return protein_to_domains

def fetch_gene_start_positions(db_path, protein_ids):
    """
    Fetches the start positions of proteins in the genome.
    """
    with sqlite3.connect(db_path) as con:
        cur = con.cursor()
        cur.execute(f"""
            SELECT proteinID, start FROM Proteins
            WHERE proteinID IN ({','.join(['?']*len(protein_ids))})
            """, tuple(protein_ids))
        
        return {row[0]: row[1] for row in cur.fetchall()}

def write_sorted_compositions_report(directory, sorted_compositions):
    """
    Saves the sorted domain compositions to a file.
    """
    report_path = os.path.join(directory, "sorted_domain_compositions.txt")
    with open(report_path, "w") as f:
        for composition, count in sorted_compositions:
            f.write(f"Cluster Occurrences: {count}\n")
            for domains in composition:
                f.write("; ".join(domains) + "\n")
            f.write("-" * 40 + "\n")

    print(f"Sorted domain compositions saved: {report_path}")
        
def fetch_neighbouring_genes_with_domains(database, protein_ids, num_workers=8):
    """
    Fetches neighboring genes and their domain compositions for given protein IDs.
    
    Returns:
    - neighbors_dict (dict): Ordered gene neighborhood information (for hit_id_reports).
    - sorted_compositions (list): Unique domain compositions sorted by occurrence.
    """
    # 1. Parallel Fetch: Get cluster IDs
    protein_id_chunks = [list(protein_ids)[i::num_workers] for i in range(num_workers)]
    with Pool(num_workers) as pool:
        results = pool.starmap(fetch_cluster_ids, [(database, chunk) for chunk in protein_id_chunks])
    
    protein_to_cluster = {k: v for result in results for k, v in result.items()}

    # 2. Parallel Fetch: Get all proteins for clusters
    cluster_ids = set(protein_to_cluster.values())
    cluster_id_chunks = [list(cluster_ids)[i::num_workers] for i in range(num_workers)]
    with Pool(num_workers) as pool:
        results = pool.starmap(fetch_cluster_proteins, [(database, chunk) for chunk in cluster_id_chunks])

    cluster_to_proteins = defaultdict(set)
    for result in results:
        for k, v in result.items():
            cluster_to_proteins[k].update(v)

    # 3. Parallel Fetch: Get protein domains
    all_protein_ids = {p for proteins in cluster_to_proteins.values() for p in proteins}
    protein_id_chunks = [list(all_protein_ids)[i::num_workers] for i in range(num_workers)]
    with Pool(num_workers) as pool:
        results = pool.starmap(fetch_protein_domains, [(database, chunk) for chunk in protein_id_chunks])

    protein_to_domains = defaultdict(list)
    for result in results:
        for k, v in result.items():
            protein_to_domains[k].extend(v)

    # 4. Fetch start positions (order genes by their genomic location)
    with Pool(num_workers) as pool:
        results = pool.starmap(fetch_gene_start_positions, [(database, chunk) for chunk in protein_id_chunks])

    protein_to_start = {k: v for result in results for k, v in result.items()}

    # 5. Construct the neighborhood dictionary for `hit_id_reports`
    neighbors_dict = {}
    for protein_id in protein_ids:
        cluster_id = protein_to_cluster.get(protein_id)
        if cluster_id and cluster_id in cluster_to_proteins:
            # Create a list of tuples (proteinID, start_position, domains)
            protein_list = [
                (p, protein_to_start.get(p, float("inf")), protein_to_domains.get(p, ["no_neighbours"]))
                for p in cluster_to_proteins[cluster_id]
            ]

            # Sort by start position
            protein_list.sort(key=lambda x: x[1])

            # Store as a single sorted list of domains
            neighbors_dict[protein_id] = [domains for _, _, domains in protein_list]
        else:
            neighbors_dict[protein_id] = [["singleton"]]

    # 6. Construct and sort unique domain compositions
    cluster_compositions = defaultdict(int)
    for cluster_id, protein_list in cluster_to_proteins.items():
        domain_composition = sorted({tuple(sorted(protein_to_domains.get(protein, ["no_neighbours"]))) for protein in protein_list})
        domain_tuple = tuple(domain_composition)
        cluster_compositions[domain_tuple] += 1

    sorted_compositions = sorted(cluster_compositions.items(), key=lambda x: x[1], reverse=True)

    return neighbors_dict, sorted_compositions

    
            
