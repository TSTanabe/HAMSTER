#!/usr/bin/python

import sqlite3
import os

from . import Csb_cluster

#first get the csb with their identifiers. make sets of appearence as value to key name of the csb
#do not compare the csb appearences as they jaccard itself should have managed it.
def csb_proteins_datasets(options):
    
    csb_dictionary = parse_csb_file_to_dict(options.csb_output_file) #dictionary with cbs_name => csb items
    pattern_dictionary = parse_csb_file_to_dict(options.patterns_file)  # Fetch the ones that are in the pattern file
    csb_dictionary = {**csb_dictionary, **pattern_dictionary}
    options.csb_dictionary = csb_dictionary # save the patterns for later use
    
    
    
    #Fetch for each csb id all the domains in the csb that are query domains
    dictionary = fetch_proteinIDs_dict(options.database_directory,csb_dictionary,options.min_seqs)
    #dictionary is: dict[(keyword, domain)] => set(proteinIDs)
    dictionary = remove_non_query_clusters(options.database_directory, dictionary) #delete all that are not in accordance with query
    
    #Remove domains that are excluded by user options
    dictionary = filter_dictionary_by_inclusion_domains(dictionary, options.include_list)
    dictionary = filter_dictionary_by_excluding_domains(dictionary, options.exclude_list)

    
    #Grouping routines to reduce redundancy in training datasets
    grouped = group_proteinID_sets(dictionary) #removed doublicates #key: frozenset of proteinIDs value: list of tuples with (keyword,domain) pairs
    
    merged, singles = merge_similar_groups(grouped,options.dbscan_epsilon,"m")
    #print("\nMerged\n\n")
    #print(merged)
    #print("\nSingles\n\n")
    #print(singles)
    options.TP_merged = merged
    options.TP_singles = singles
    print_grouping_report("Merged csb groupes that include highly similar or identical proteins groups",merged,options.Csb_directory+"/Merged_csb_groups_1")
    print_grouping_report("Unique csb groupes that include distinct proteins groups",merged,options.Csb_directory+"/Singles_csb_groups_1")


def training_data_fasta(options):
    
    training_datasets = {**options.TP_merged, **options.TP_singles, **options.TP_monophyla, **options.superfamily}
    
    dicts = [options.TP_merged, options.TP_singles, options.TP_monophyla]
    training_datasets = {}
    
    for d in dicts:
        for key,value in d.items():
            if key in training_datasets:
                training_datasets[key] += value
            else:
                training_datasets[key] = value
    
    merged, singles = merge_similar_groups(training_datasets, options.dbscan_epsilon,"p")
    print_grouping_report("Merged csb groupes that include highly similar or identical proteins groups",merged,options.Csb_directory+"/Merged_phylogenetic_groups_2")
    print_grouping_report("Unique csb groupes that include distinct proteins groups",merged,options.Csb_directory+"/Singles_csb_groups_2")

    #print(len(options.TP_merged))
    #print(len(options.TP_singles))
    #print(len(options.TP_monophyla))
    #print(len(merged))
    #print(len(singles))

    
    fetch_protein_to_fasta(options,merged)
    fetch_protein_to_fasta(options,singles)
    
    
def parse_csb_file_to_dict(file_path):
    data_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Split the line by tabs
            parts = line.strip().split('\t')

            # The first part is the identifier (key)
            identifier = parts[0]

            # Initialize a set to store unique values from the tuples
            value_set = set()

            # Iterate over the remaining parts (tuples)
            for tuple_str in parts[1:]:
                if tuple_str:  # Check if the part is not empty
                    # Manually parse the tuple string
                    # Remove parentheses and split by commas
                    tuple_str = tuple_str.strip("()")  # Remove parentheses
                    values = tuple_str.split(", ")  # Split by comma and space
                    
                    # Remove quotes around the values and add them to the set
                    for value in values:
                        value_set.add(value.strip("'"))

            # Add the identifier and the set of values to the dictionary
            data_dict[identifier] = value_set

    return data_dict    


def fetch_proteinIDs_dict(database, csb_dictionary, min_seqs):
    """
    Fetches the results from the SQLite database based on domain and keyword, stores it in a dictionary.

    Args:
        database: Name of the SQLite database.
        csb_dictionary: Dictionary where the key is 'csb' and the value is a set of domains.
        min_seqs: Minimum number of sequences required to include in the result.

    Returns:
        A dictionary where the key is a tuple (keyword, domain) and the value is a set of protein IDs.
    """
    csb_proteinIDs_dict = {}  # (keyword, domain) => set of proteinIDs
    
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Iterate over the keywords and domains in the dictionary
        for keyword, domains in csb_dictionary.items():
            if keyword == 'default':
                continue
            
            # Prepare a query using `LIKE` for pattern matching
            query_conditions = ' OR '.join(['Domains.domain LIKE ?'] * len(domains))
            query = f"""
                SELECT DISTINCT Proteins.proteinID, Domains.domain
                FROM Proteins
                INNER JOIN Domains ON Proteins.proteinID = Domains.proteinID
                INNER JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
                WHERE ({query_conditions}) AND Keywords.keyword = ?;
            """

            # Add wildcard patterns for LIKE (e.g., '%Hydrolase%')
            like_domains = [f"%{domain}%" for domain in domains]

            # Execute the query with the domains (with wildcards) and keyword as parameters
            cur.execute(query, (*like_domains, keyword))
            
            rows = cur.fetchall()
            #print(f"{query} {like_domains} {keyword}" + str(len(rows)))
            # Group results by (keyword, domain) and apply filtering based on min_seqs
            for proteinID, domain in rows:
                key = (keyword, domain)
                if key not in csb_proteinIDs_dict:
                    csb_proteinIDs_dict[key] = set()
                csb_proteinIDs_dict[key].add(proteinID)

            # Apply filtering based on min_seqs
            for key in list(csb_proteinIDs_dict.keys()):
                if len(csb_proteinIDs_dict[key]) <= min_seqs:
                    del csb_proteinIDs_dict[key]

    return csb_proteinIDs_dict  # dict[(keyword, domain)] => set(proteinIDs)


######################################################################################################


def fetch_query_clusters(database, dictionary):
    """
    First selects the protein domain types that were assigned by the selfblast.
    Then removes every protein from the dictionary that is not equal to the
    domain types from the selfblast.
    
    Args:
        database (str): Path to the SQLite database.
        dictionary (dict): A dictionary where the keys are tuples, and the second
                           element in each tuple represents a protein domain.
    
    Returns:
        dict: A dictionary with only the entries where the domain matches one
              of the domains from the query in the database.
    """
    query_domains = set()  # Set to hold the selected domains

    
    # Connect to the SQLite database
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        # Query to get distinct domains from the selfblast associated with genomeID 'QUERY'
        query = """
            SELECT DISTINCT Domains.domain
            FROM Proteins
            LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
            WHERE Proteins.genomeID = ?;
        """
        # Execute the query with 'QUERY' as the genomeID
        cur.execute(query, ('QUERY',))
        rows = cur.fetchall()
        
        # Add each domain to the query_domains set
        for row in rows:
            query_domains.add(row[0])
    return query_domains

    
    

def remove_non_query(dictionary,query_domains):    
    #dictionary has the csb => set of domains
    #query domains
    # Loop through the dictionary and filter based on the domain in the key
    selected_dictionary = dict()  # Dictionary to hold filtered results
    for key, present_genes in dictionary.items():
        common_genes = present_genes.intersection(query_domains)
        if common_genes:
            selected_dictionary[key] = common_genes
    
    return selected_dictionary   

def remove_non_query_clusters(database, dictionary):
    """
    First selects the protein domain types that were assigned by the selfblast.
    Then removes every protein from the dictionary that is not equal to the
    domain types from the selfblast.
    
    Args:
        database (str): Path to the SQLite database.
        dictionary (dict): A dictionary where the keys are tuples, and the second
                           element in each tuple represents a protein domain.
    
    Returns:
        dict: A dictionary with only the entries where the domain matches one
              of the domains from the query in the database.
    """
    query_domains = set()  # Set to hold the selected domains
    selected_dictionary = dict()  # Dictionary to hold filtered results
    
    # Connect to the SQLite database
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        # Query to get distinct domains from the selfblast associated with genomeID 'QUERY'
        query = """
            SELECT DISTINCT Domains.domain
            FROM Proteins
            LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
            WHERE Proteins.genomeID = ?;
        """
        # Execute the query with 'QUERY' as the genomeID
        cur.execute(query, ('QUERY',))
        rows = cur.fetchall()
        
        # Add each domain to the query_domains set
        for row in rows:
            query_domains.add(row[0])
    # Loop through the dictionary and filter based on the domain in the key
    for key, value in dictionary.items():
        if key[1] in query_domains:
            # If the domain in the key is part of query_domains, copy it to selected_dictionary
            selected_dictionary[key] = value
    
    return selected_dictionary

###################################################################################################

def group_proteinID_sets(csb_proteinIDs_dict):
    """
    Groups keyword-domain pairs by their sets of protein IDs.

    Args:
        csb_proteinIDs_dict: Dictionary returned by fetch_proteinIDs_dict, where the keys are (keyword, domain) tuples and the values are sets of protein IDs.

    Returns:
        A dictionary where the keys are frozensets of protein IDs and the values are lists of keyword-domain pairs that share the same set.
    """
    # Dictionary to store groups of keyword-domain pairs by their sets
    grouped_sets = {}

    # Iterate over the csb_proteinIDs_dict and group by proteinID sets
    for (keyword, domain), proteinIDs in csb_proteinIDs_dict.items():
        # Convert the set to a frozenset to use it as a key in the dictionary
        proteinID_frozenset = frozenset(proteinIDs)
        
        # Add the keyword-domain pair to the corresponding group
        if proteinID_frozenset not in grouped_sets:
            grouped_sets[proteinID_frozenset] = []
        grouped_sets[proteinID_frozenset].append((keyword, domain))
    
    # This loop is for debugging
    #for proteinID_set, keyword_domain_pairs in grouped_sets.items():
    #    print(f"Group with protein IDs: {proteinID_set}")
    #    print("Keyword-Domain Pairs:", keyword_domain_pairs)
    #    print()
    
    return grouped_sets #key: frozenset of proteinIDs value: list of tuples with (keyword,domain) pairs


def merge_similar_groups(grouped_sets, epsilon, merge_extension="m"):
    #epsilon is the acceptable dissimilarity to merge points
    
    # Step 1: Extract the frozen sets (keys) from the dictionary
    cluster_columns = list(grouped_sets.keys())
    cluster_columns.sort(key=lambda fs: tuple(sorted(fs))) #sort the list to have consistent input order between runs

    # Step 2: Calculate the Jaccard similarity matrix for the frozen sets
    similarity_matrix = Csb_cluster.calculate_similarity_matrix_jaccard(cluster_columns)
    
    # Step 3: Apply DBSCAN to get the cluster labels based on the similarity matrix
    labels = Csb_cluster.apply_dbscan_on_similarity_matrix(similarity_matrix, eps=epsilon, min_samples=2)
    
    # Step 4: Create two dictionaries: one for merged groups and one for noise points
    merged_groups = {}
    noise_points = {}
    
    # Step 5: Group sets by their DBSCAN label
    visited_labels = set()
    for idx, label in enumerate(labels):  # Loop through unique labels
        merged_set = set()  # Will contain the union of the sets
        merged_list = []  # Will contain the merged list of tuples
        
        if label == -1:
            # Handle noise points (label -1)
            modified_identifiers = [("s" + item[0], item[1]) for item in grouped_sets[cluster_columns[idx]]]
            noise_points[cluster_columns[idx]] = modified_identifiers
        elif not label in visited_labels:
            # For non-noise points that were not already visited, merge sets with the same label
        
        
            for idx2, set_label in enumerate(labels):
                if set_label == label:
                    # Merge the sets with the same label
                    merged_set = merged_set.union(cluster_columns[idx2])  # Union of frozen sets
                
                    # Extend merged_list, adding merge_extension to each identifier in the list
                    modified_identifiers = [(merge_extension + item[0], item[1]) for item in grouped_sets[cluster_columns[idx2]]]
                    merged_list.extend(modified_identifiers)  # Merge with modified identifiers
                    visited_labels.add(set_label)

            # Store the merged set (converted to frozenset to make it a valid dictionary key)
            merged_groups[frozenset(merged_set)] = merged_list
    
    # Return both the merged groups and the noise points
    return merged_groups, noise_points

def print_grouping_report(description, group, output):
    """
    Writes a description followed by a linewise print of the group dictionary values to the output file,
    separating lines by a blank line.

    :param description: A string containing the description to be printed.
    :param group: A dictionary where the values are lists or other iterables.
    :param output: The file object to write the report to.
    """
    with open(output, 'w') as file:
        # Write the description to the file
        file.write(f"{description}\n\n")

        # Iterate over the group dictionary
        for key, values in group.items():
            # Write the key and values to the file
            file.write(f"Group: {key}\n")
            for value in values:
                file.write(f"{value}\n")
            file.write("\n")  # Add a blank line between groups





###############################################################################
#################### Protein to fasta operations ##############################
###############################################################################

def fetch_protein_to_fasta(options, grouped, prefix=""):
    """
    Fetches protein sequences from the Proteins table for protein IDs in the grouped sets,
    and writes the sequences directly to a FASTA file.

    Args:
        options: Object containing configuration options, including the database directory.
        grouped: Dictionary where the keys are frozensets of protein IDs and the values are lists of keyword-domain pairs.

    Returns:
        None, but writes the sequences to a FASTA file.
    """
    with sqlite3.connect(options.database_directory) as con:
        cursor = con.cursor()

        for proteinID_frozenset, keyword_domain_pairs in grouped.items():
            # Convert frozenset to a tuple for the SQL IN clause
            proteinID_tuple = tuple(proteinID_frozenset)
            
            #Do not fetch if the minimal number of sequences is not reached
            if len(proteinID_tuple) < options.min_seqs:
                print(f"Info: Less than {options.min_seqs} with {keyword_domain_pairs}. Skipping this set.")
                continue
                
            # Fetch protein sequences for the given protein IDs
            query = """
                SELECT proteinID, sequence
                FROM Proteins
                WHERE proteinID IN ({})
            """.format(','.join(['?'] * len(proteinID_tuple)))

            cursor.execute(query, proteinID_tuple)
            rows = cursor.fetchall()

            # Write the sequences directly to the FASTA file
            keyword, domain = keyword_domain_pairs[0]
            filename = f"{prefix}{keyword}_{domain}.faa"
            output_fasta_path = os.path.join(options.fasta_output_directory, filename)
            with open(output_fasta_path, 'w') as fasta_file:
                for proteinID, sequence in rows:
                    # Create a header with keyword-domain information
                    header = f">{proteinID}"
                    
                    # Write the header and sequence to the file
                    fasta_file.write(f"{header}\n")
                    fasta_file.write(f"{sequence}\n")

    
    







def fetch_all_proteins(database, filepath):
    """
    Fetches the results from the SQLite database based on domain and keyword, and stores them in a FASTA file.

    Args:
        database: Name of the SQLite database.
        results_dict: Dictionary containing domain as keys and keyword as values.
        output_directory: Directory path to store the FASTA file.
        
        Output is a files with dir/sequence.faa
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        
        #get all proteins with id and key output shall be "proteinId csb_proteintype" 
        query = """
                SELECT DISTINCT Proteins.proteinID,Proteins.Sequence
                FROM Proteins WHERE NOT Proteins.genomeID = ?
            """
        cur.execute(query, ('QUERY',))
        rows = cur.fetchall()
        with open(filepath,'w') as fasta_file:
            for row in rows:
                fasta_file.write(f">{row[0]}\n")
                fasta_file.write(f"{row[1]}\n")

    return filepath
    
    
def filter_dictionary_by_inclusion_domains(input_dict, include_list):
    """
    Filters a dictionary to retain only entries with keys that match domains in a given list.
    If the domain list is None or empty, the original dictionary is returned.

    Args:
        input_dict (dict): The dictionary with keys as (keyword, domain) tuples and values as sets of proteinIDs.
        domain_list (list): A list of domains to retain. If None or empty, no filtering is performed.

    Returns:
        dict: A filtered dictionary with only the specified domains retained, or the original dictionary.
    """
    if not include_list:  # Check if domain_list is None or empty
        return input_dict

    # Use dictionary comprehension to filter the dictionary
    filtered_dict = {
        key: value
        for key, value in input_dict.items()
        if key[1] in include_list  # Check if the domain is in the provided list
    }
    return filtered_dict

def filter_dictionary_by_excluding_domains(input_dict, exclude_list):
    """
    Filters a dictionary to exclude entries with keys that match domains in a given list.
    If the exclude list is None or empty, the original dictionary is returned.

    Args:
        input_dict (dict): The dictionary with keys as (keyword, domain) tuples and values as sets of proteinIDs.
        exclude_list (list): A list of domains to exclude. If None or empty, no filtering is performed.

    Returns:
        dict: A filtered dictionary with specified domains excluded, or the original dictionary.
    """
    if not exclude_list:  # Check if exclude_list is None or empty
        return input_dict

    # Use dictionary comprehension to exclude the specified domains
    filtered_dict = {
        key: value
        for key, value in input_dict.items()
        if key[1] not in exclude_list  # Exclude keys where domain is in the exclude list
    }
    return filtered_dict
