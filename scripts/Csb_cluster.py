#!/usr/bin/python


import os
from . import Csb_Mp_Algorithm
from scipy.spatial import distance
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from itertools import combinations
from collections import defaultdict

# For the clustering of csbs by jaccard and agglomerativeClustering
def csb_prediction(options):
    #Finds collinear syntenic blocks with the csb finder algorithm using a printed representation of the clusters.
    options.redundant,options.non_redundant = dereplicate(options.gene_clusters_file) #returns two filepaths, dereplicates identical gene clusters
    
    options.redundancy_hash = create_redundancy_hash(options.redundant) # value is an integer, number of existing replicates
    gene_clusters = create_gene_cluster_hash(options.non_redundant)
    extend_redundancy_hash(options.non_redundant,options.redundancy_hash) # for all which do not have a redundant gene cluster
    #modified CsbfinderS algorithm
    computed_Instances_dict = Csb_Mp_Algorithm.csb_finderS_matchpoint_algorithm(options.redundancy_hash,gene_clusters,options.insertions,options.occurence) #k insertions und q occurences müssen über die optionen festgelegt werden

    # Combine reverse csbs
    reverse_pairs = find_reverse_pairs(computed_Instances_dict)
    computed_Instances_dict = merge_sets_for_pairs(computed_Instances_dict, reverse_pairs)
    
    options.computed_Instances_dict = computed_Instances_dict
    
def csb_jaccard(options):
    #convert the keys in computed_instances_dict into a list
    computed_Instances_key_list = csb_Instance_key_list(options.computed_Instances_dict, options.min_csb_size)
    cluster_dict = dict()
    if len(computed_Instances_key_list) > 1:
        matrix = calculate_similarity_matrix_jaccard(computed_Instances_key_list) #matrix of similarity and the corresponding clusterID for each row and column as names
        cluster_dict = hierachy_clustering(matrix,options.jaccard) # 0.2 means that 80 % have to be the same genes
    elif len(computed_Instances_key_list) == 1:
        cluster_dict[0] = [0]
    else:
        return
    #Sorts the csb keywords to the geneclusters based on the present csb
    csb_gene_cluster_dict, grouped_csb_tuples = csb_index_to_gene_clusterID(cluster_dict,computed_Instances_key_list,options.computed_Instances_dict,options.csb_name_prefix,options.csb_name_suffix)
    #adds the replicates again for saving
    replicates(csb_gene_cluster_dict,options.redundancy_hash,options.redundant)
    
    write_grouped_csb(options.csb_output_file,grouped_csb_tuples)

    return csb_gene_cluster_dict
    








def print_cluster(filepath, cluster_dict):
    with open(filepath, 'a') as writer:
        for cluster in cluster_dict:
            writer.write(cluster.clusterID)  # Write cluster ID (assuming cluster.clusterID is a string)
            domains = '\t'.join(cluster.domains)  # Join the domains list into a tab-separated string
            writer.write("\t" + domains + "\n")  # Write cluster domains followed by a newline

    return filepath



########################################################################################
################ Secondary subroutines #################################################
################ Subroutines used here #################################################
########################################################################################

def dereplicate(filepath):
    # Dictionary to store lines based on variable columns
    line_dict = {}

    # List to store non-redundant lines
    non_redundant_lines = []

    # Read the file and check for duplicates
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            identifier, *variables = line.split('\t')
            variable_container = tuple(variables)
            if variable_container in line_dict:
                line_dict[variable_container].append(identifier)
            elif reversed(variable_container) in line_dict:
                line_dict[variable_container].append(identifier)
            else:
                line_dict[variable_container] = [identifier]
                non_redundant_lines.append(line)

    # Print the list of redundant identifiers
    redundant_identifiers = []
    for identifiers in line_dict.values():
        if len(identifiers) > 1:
            redundant_identifiers.extend(identifiers)

    # Write redundant identifiers to a file
    redundant_file = os.path.join(os.path.dirname(filepath), 'redundant.txt')
    with open(redundant_file, 'w') as f:
        for identifiers in line_dict.values():
            if len(identifiers) > 1:
                f.write('\t'.join(identifiers)+'\n')
                

    # Write non-redundant lines to a file
    non_redundant_file = os.path.join(os.path.dirname(filepath), 'non-redundant.txt')
    with open(non_redundant_file, 'w') as f:
        f.write('\n'.join(non_redundant_lines))

    return redundant_file, non_redundant_file


def create_redundancy_hash(file_path):
    result_dict = {}  # Initialize an empty dictionary

    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')

            if columns:  # Check if there are any columns in the line
                key = columns[0]  # The value from the first column is the key
                value = len(columns)  # Number of columns is the value
                
                # Store the key and value in the dictionary
                result_dict[key] = value
    return result_dict

def create_gene_cluster_hash(file_path):
    result_dict = {}  # Initialize an empty dictionary

    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')

            if columns:  # Check if there are any columns in the line
                key = columns[0]  # The value from the first column is the key
                values = columns[1:]  # Values from the other columns as a list

                # Store the key and values in the dictionary
                result_dict[key] = values
                

    return result_dict

def extend_redundancy_hash(filepath,redundancy_hash):
    with open(filepath, 'r') as file:
        for line in file:
            # Split the line into columns using tab as the delimiter
            columns = line.strip().split('\t')

            if columns:  # Check if there are any columns in the line
                key = columns[0]  # The value from the first column is the key
                if not key in redundancy_hash:
                    redundancy_hash[key] = 1
    return redundancy_hash
    
    
def csb_Instance_key_list(Instance_dict,threshold):
    listing = list()
    for k in Instance_dict.keys():
        if len(k)>=threshold:
            listing.append(k)

    return listing

def write_grouped_csb(filepath, data):
    # Writes down which csb keyword belongs to which csb in a tab-separated file
    with open(filepath, 'w') as f:
        # Iterate through the dictionary items
        for key, value in data.items():
            # Write the key to the file, followed by a tab
            f.write(f"{key}\t")
            # Iterate through the list of tuples
            for item in value:
                # Write each tuple as a tab-separated string
                f.write(f"{item}\t")
            # Add a newline after each entry
            f.write('\n')            
###########################################################################
# Added 19.10.2024 to specifically cluster the reverse csb

def find_reverse_pairs(my_dict):
    reverse_pairs = []
    seen = set()  # Um Duplikate zu vermeiden
    
    for key in my_dict.keys():
        reversed_key = key[::-1]  # Erzeuge das Revers-Tupel
        
        # Überprüfe, ob das Revers in den Keys ist und das Paar noch nicht überprüft wurde
        if reversed_key in my_dict and reversed_key not in seen:
            reverse_pairs.append((key, reversed_key))
            seen.add(key)  # Füge das Paar zu "gesehenen" hinzu
            seen.add(reversed_key)
    
    return reverse_pairs

def merge_sets_for_pairs(my_dict, reverse_pairs):
    """
    Merges sets of reverse key pairs and keeps only one of the keys in the dictionary.
    
    Args:
        my_dict (dict): Original dictionary with tuples as keys and sets as values.
        reverse_pairs (list): List of tuple pairs where the second tuple is the reverse of the first.
    
    Returns:
        dict: The updated dictionary with merged sets and only one key per reverse pair.
    """
    # Iteriere über die Reverse-Paare und führe die Sets zusammen
    for key1, key2 in reverse_pairs:
        if key1 in my_dict and key2 in my_dict:
            # Führe die Sets von key1 und key2 zusammen
            my_dict[key1] = my_dict[key1].union(my_dict[key2])
            # Entferne key2 aus dem Dictionary
            del my_dict[key2]
    
    return my_dict

###########################################################################
# Read the gene clusters from the file
def jaccard(set1,set2):
    if not isinstance(set1, set) or not isinstance(set2,set):
        try:
            set1 = set(set1)
            set2 = set(set2)
        except Exception as e:
            print("An error occurred while converting to sets:", e)
            return 0
            
    intersection = len(set1.intersection(set2))
    union = len(set1) + len(set2) - intersection
    jaccard_similarity = intersection / union if union != 0 else 0

    return jaccard_similarity

def calculate_similarity_matrix_jaccard_deprectaed(cluster_columns):
    #cluster_columns needs to be a list of sets, not tuple because jaccard relies on the set datatype
    #
    #THIS ROUTINE NEEDS NUMPY. wich i tried to remove
    
    # Determine the number of clusters
    num_clusters = len(cluster_columns) if cluster_columns else 0


    # Calculate Jaccard similarity matrix
    jaccard_generator = (jaccard(row1, row2) for row1, row2 in combinations(cluster_columns, r=2))
    flattened_matrix = np.fromiter(jaccard_generator, dtype=np.float64)

    # since flattened_matrix is the flattened upper triangle of the matrix
    # we need to expand it.
    similarity_matrix = distance.squareform(flattened_matrix)
    similarity_matrix += np.identity(len(cluster_columns))# replacing zeros with ones at the diagonal. 
    
    print(similarity_matrix)
    sys.exit()
    return similarity_matrix

def calculate_similarity_matrix_jaccard(cluster_columns):
    #cluster_columns needs to be a list of sets, not tuple because jaccard relies on the set datatype
    
    # Convert sets to a 2D binary matrix
    unique_elements = sorted(set().union(*cluster_columns))  # Find all unique elements across the sets
    binary_matrix = [[1 if element in cluster else 0 for element in unique_elements] for cluster in cluster_columns]
    
    # Compute the Jaccard distance (1 - Jaccard similarity)
    jaccard_distances = distance.pdist(binary_matrix, metric='jaccard')
    
    # Convert the distances to a full similarity matrix using squareform
    similarity_matrix = 1 - distance.squareform(jaccard_distances)
    
    # Add 1s to the diagonal to represent self-similarity
    for i in range(len(cluster_columns)):
        similarity_matrix[i, i] = 1.0


    return similarity_matrix

def hierachy_clustering(similarity_matrix,threshold):

    avg_dissim_threshold = threshold

    # Compute the dissimilarity matrix
    dissimilarity_matrix = 1 - similarity_matrix
    #print(dissimilarity_matrix)
    # Perform Agglomerative Clustering with average linkage
    clustering = AgglomerativeClustering(n_clusters=None, metric='precomputed', linkage='average', distance_threshold=avg_dissim_threshold)
    
    # Fit the clustering model and obtain the cluster labels
    cluster_labels = clustering.fit_predict(dissimilarity_matrix)

    # Create a dictionary to store the clusters
    clusters = defaultdict(list)
    for i, label in enumerate(cluster_labels):
        clusters[label].append(i)

    #returns a dictionary, key is a number identifying the cluster (eine ganze zahl). value is a list with numbers. these numbers
    #are the cluster labels. this dictionary is pushed forward to the subroutine Clusters.keyword_dictionary
    return clusters

def csb_index_to_gene_clusterID(cluster_dict,computed_Instances_key_list,computed_Instances_dict,prefix="csb*",suffix="*"):
    result_dict = dict() # csbID to geneclusters
    result_dict2 = dict() # csbID to csb tuples
    
    for key,indices in cluster_dict.items():
        tuples = [computed_Instances_key_list[i] for i in indices] #list of tuples corresponding to the indices from jaccard clustering
        result_dict2[prefix+str(key)+suffix] = tuples
        for i in tuples:
            for e in computed_Instances_dict[i]: #returns all clusterIDs of the csb 'i'
                if prefix+str(key)+suffix in result_dict.keys():
                    result_dict[prefix+str(key)+suffix].append(e)
                else:
                    result_dict[prefix+str(key)+suffix] = [e]
    return result_dict,result_dict2
    

def replicates(csb_gene_cluster_dict, redundancy_hash, filepath_redundant):
    redundant_dict = dict()
    # Read cluster IDs from the redundant file
    with open(filepath_redundant, 'r') as file:
        for line in file:
            clusterIDs = line.strip().split('\t')
            first_clusterID = clusterIDs[0]
            redundant_dict[first_clusterID] = clusterIDs[1:]
    # Iterate through the csb_gene_cluster_dict and merge with redundant_dict
    for key, lst in csb_gene_cluster_dict.items():
        new_list = lst.copy()  # Create a copy to avoid modifying the original list
        for cluster_id in lst:
            if cluster_id in redundant_dict.keys():
                new_list.extend(redundant_dict[cluster_id])
        csb_gene_cluster_dict[key] = new_list
    return csb_gene_cluster_dict  # Return the updated dictionary

######################################################################
################ For proteinID clustering ############################
################ in Csb_proteins          ############################
######################################################################

def apply_dbscan_on_similarity_matrix(similarity_matrix, eps=0.5, min_samples=2):
    """
    Applies DBSCAN clustering on a Jaccard similarity matrix.

    Args:
        similarity_matrix (numpy.ndarray): A square matrix of Jaccard similarities between sets.
        eps (float): The maximum distance between two points to be considered as neighbors.
        min_samples (int): The minimum number of points required to form a cluster.

    Returns:
        numpy.ndarray: An array of cluster labels, where -1 indicates noise.
    """
    # Step 1: Convert similarity matrix to distance matrix (1 - similarity)
    distance_matrix = 1 - similarity_matrix

    # Step 2: Apply DBSCAN with the precomputed distance matrix
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    labels = dbscan.fit_predict(distance_matrix)

    return labels    



