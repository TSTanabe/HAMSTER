#!/usr/bin/python

import sqlite3
import os


#first get the csb with their identifiers. make sets of appearence as value to key name of the csb
#do not compare the csb appearences as they jaccard itself should have managed it.
def csb_proteins_fasta(options):
    
    
    csb_dictionary = parse_csb_file_to_dict(options.csb_output_file)
    dictionary = fetch_proteinIDs_dict(options.database_directory,csb_dictionary,options.min_seqs) #Fetches for each csb id all the domains in the csb, regardless of clustering: 
    #dictionary is: dict[(keyword, domain)] => set(proteinIDs)
    dictionary = remove_non_query_clusters(options.database_directory, dictionary) #delete all that are not in accordance with query
    
    #Grouping routines to reduce redundancy in training datasets
    grouped = group_proteinID_sets(dictionary) #removed doublicates
    subset = find_subsets(grouped)
    merged = merge_subsets(subset, grouped, options.csb_overlap_factor, options.Csb_directory) #merged core subsets that are above the overlap_factor

    #fetch the grouped protein sets
    if options.csb_distinct_grouping:
        fetch_protein_to_fasta(options,grouped)
    
    
    #fetch the merged subsets
    fetch_protein_to_fasta(options,merged,"Merged_")
    
    options.TP_grouped = grouped
    options.TP_merged = merged
    
    #generate target fasta file    
    options.sequence_faa_file = fetch_all_proteins(options.database_directory, options.cross_validation_directory+"/sequences.faa") #File with all sequences to be searched
    
    
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


def find_subsets(grouped_sets):
    """
    Routine just for the information of the user
    Finds and stores information about subset relationships between sets in grouped_sets.

    Args:
        grouped_sets: Dictionary where keys are frozensets of protein IDs and values are lists of (keyword, domain) tuples.

    Returns:
        A dictionary where the keys are frozensets and the values are lists of frozensets that are supersets of the key.
    """
    subset_info = {}  # To store subset relationships

    # Convert the grouped_sets dictionary into a list of tuples for easy comparison
    grouped_list = list(grouped_sets.keys())

    # Compare each set against every other set
    for i, smaller_set in enumerate(grouped_list):
        for j, larger_set in enumerate(grouped_list):
            if i != j and smaller_set < larger_set:  # Check if smaller_set is a proper subset of larger_set
                if smaller_set not in subset_info:
                    subset_info[smaller_set] = []
                subset_info[smaller_set].append(larger_set)
    

    return subset_info

def merge_subsets(subset_info, grouped_sets, grouping_threshold, csb_directory):
    reversed_merged_subsets = {}
    output_file_path = os.path.join(csb_directory, "Csb_sequence_subsets.txt")

    with open(output_file_path, 'w') as output_file:
        # Debugging loop for showing the results
        for smaller_set, larger_sets in subset_info.items():
            key_smaller = grouped_sets[smaller_set]
            core_len = len(smaller_set)

            output_file.write(f"Set {key_smaller} {core_len}\n is a subset of the following sets:\n")
            var_len = 0
            for larger_set in larger_sets:
                key_larger = grouped_sets[larger_set]
                var_len += len(larger_set) - core_len
                output_file.write(f" {key_larger} \n")
            
            output_file.write(f"Core length is {core_len} and var length is {var_len}\n\n")
            
            total_len = core_len + var_len
            if core_len / total_len > grouping_threshold:
                # Merge the large and small sets
                merge_set = set(smaller_set.copy())
                for larger_set in larger_sets:
                    merge_set.update(larger_set)
                
                output_file.write(f"Merged the sets core {key_smaller} to the larger ones: \n")
                
                merge_set_frozenset = frozenset(merge_set)
                
                # If the frozenset already exists, append the new key_smaller[0]
                if merge_set_frozenset in reversed_merged_subsets:
                    reversed_merged_subsets[merge_set_frozenset].append(key_smaller[0])
                else:
                    # Otherwise, create a new entry with the frozenset as the key
                    reversed_merged_subsets[merge_set_frozenset] = [key_smaller[0]]

    return reversed_merged_subsets

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

    
    
    
def deprecated_fetch_proteins_to_fasta(database, protein, proteinID_sets_csb_dict, output_directory):
    """
    Fetches the results from the SQLite database based on domain and keyword, and stores them in a FASTA file.

    Args:
        database: Name of the SQLite database.
        output_directory: Directory path to store the FASTA file.
        
    """
    with sqlite3.connect(database) as con:
        cursor = con.cursor()
        
        for proteinIDs,csbs in proteinID_sets_csb_dict.items():
            placeholders_csbs = ','.join(['?']*len(csbs))
            placeholders_proteinIDs = ','.join(['?']*len(proteinIDs))
            query = f"""
                SELECT DISTINCT Proteins.proteinID,Keywords.keyword,Domains.domain,Proteins.Sequence
                FROM Proteins
                LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
                LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
                WHERE Domains.domain = ? AND Keywords.keyword IN ({placeholders_csbs}) AND Proteins.proteinID IN ({placeholders_proteinIDs})
            """
            
            args = list()
            args.append(protein) #domain type binding variable
            
            for csb in csbs: #csb identifier binding variables
                args.append(csb)
            
            for proteinID in proteinIDs: #protein identifier binding variables
                args.append(proteinID)
              
            cursor.execute(query, tuple(args))
            
            results = cursor.fetchall()
            
            csb = '-'.join(sorted(csbs))
            filename = f"{csb}_{protein}.faa"
            filepath = os.path.join(output_directory, filename)
            
            if len(csbs) == 1:
                with open(filepath, 'w') as fasta_file:
                    for row in results:
                        fasta_file.write(f">{row[0]}**{row[1]}_{row[2]}\n")
                        fasta_file.write(f"{row[3]}\n")
            else:
                written_proteinIDs = set() #proteinsIDs that are already written
                with open(filepath, 'w') as fasta_file:
                    for row in results:
                        if not row[1] in written_proteinIDs:
                            fasta_file.write(f">{row[0]}**{row[1]}_{row[2]}\n")
                            fasta_file.write(f"{row[3]}\n")
                            written_proteinIDs.add(row[1])
            
    return #fasta files have >proteinID**key_domain













def deprecated_fetch_domains_per_key_to_fasta(database, pairs_set, output_directory, min_seqs=2):
    """
    Fetches the results from the SQLite database based on domain and keyword, and stores them in a FASTA file.

    Args:
        database: Name of the SQLite database.
        results_dict: Dictionary containing domain as keys and keyword as values.
        output_directory: Directory path to store the FASTA file.
        
        Output are files with {domain}_{keyword}.faa
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Iterate over the domain and keyword dictionary
        for tup in pairs_set:
            query = """
                SELECT DISTINCT Proteins.proteinID,Keywords.keyword,Domains.domain,Proteins.Sequence
                FROM Proteins
                LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
                LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
                WHERE domain = ? AND keyword = ?
            """
            
            if tup[1] == 'default':
                continue
            cur.execute(query, tup)

                # Fetch the results
            rows = cur.fetchall()
            if len(rows) > min_seqs:
                # Create the filename for the FASTA file
                filename = f"{tup[1]}_{tup[0]}.faa"
                protein_type = f"{tup[1]}_{tup[0]}"
                filepath = os.path.join(output_directory, filename)
                # Write the results to the FASTA file
                with open(filepath, 'w') as fasta_file:
                    for row in rows:
                        # Write the data in FASTA format
                        # Modify this part according to your table structure and data
                        fasta_file.write(f">{row[0]}**{row[1]}_{row[2]}\n")
                        fasta_file.write(f"{row[3]}\n")

    return







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
    
    



