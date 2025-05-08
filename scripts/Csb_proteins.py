#!/usr/bin/python
import sqlite3
import os
import pickle
import multiprocessing
from multiprocessing import Pool, Manager, Value, Lock

from . import Csb_cluster

#first get the csb with their identifiers. make sets of appearence as value to key name of the csb
#do not compare the csb appearences as they jaccard itself should have managed it.
def csb_proteins_datasets(options, sglr_dict):
    
    #Get the domain types as set per csb and predefined pattern
    print("Reading collinear syntenic block files")
    csb_dictionary = parse_csb_file_to_dict(options.csb_output_file) #dictionary with cbs_name => csb items
    pattern_dictionary = parse_csb_file_to_dict(options.patterns_file)  # Fetch the ones that are in the pattern file
    csb_dictionary = {**csb_dictionary, **pattern_dictionary}
    options.csb_dictionary = csb_dictionary # save the patterns for later use

    # Try to retrieve saved dictionary
    cache_dir = os.path.join(options.result_files_directory, "pkl_cache") 
    os.makedirs(cache_dir, exist_ok=True)
    
    # Load existing data
    dictionary = load_cache(cache_dir, "csb_protein_dataset.pkl")
    
    if dictionary:
        print("Loaded existing protein dataset")
        dictionary = filter_dictionary_by_inclusion_domains(dictionary, options.include_list)
        dictionary = filter_dictionary_by_excluding_domains(dictionary, options.exclude_list)
        return dictionary    

    #Fetch for each csb id all the domains in the csb that are query domains
    #dictionary is: dict[(keyword, domain)] => set(proteinIDs)
    
    if not options.sglr: #if not sglr is true, default sglr is false
        csbs_to_remove = {csb for csb_list in sglr_dict.values() for sublist in csb_list for csb in sublist}
        csb_dictionary = {csb: domains for csb, domains in csb_dictionary.items() if csb not in csbs_to_remove}

    print("Reading proteinIDs from database")
    dictionary = fetch_proteinIDs_dict_multiprocessing(options.database_directory,csb_dictionary,options.min_seqs,options.cores)

    dictionary = remove_non_query_clusters(options.database_directory, dictionary) #delete all that are not in accordance with query
    
    save_cache(cache_dir, "csb_protein_dataset.pkl", dictionary)
    
    #Remove domains that are excluded by user options
    dictionary = filter_dictionary_by_inclusion_domains(dictionary, options.include_list)
    
    dictionary = filter_dictionary_by_excluding_domains(dictionary, options.exclude_list)

    return dictionary

################################################################################################

def csb_proteins_datasets_combine(keyword_lists, csb_proteins_dict, category):
    """
    Kombiniert Protein-IDs basierend auf den gruppierten und ausgeschlossenen Keywords.
    
    Args:
        keyword_lists (dict): Dictionary mit gruppierten Keywords als Listen von Listen.
        csb_proteins_dict (dict): Das Rückgabedictionary aus csb_proteins_datasets, 
                                  das (keyword, domain) als Schlüssel und Protein-IDs als Werte enthält.
        category (str): Kategorie der Keywords (z. B. "grouped" oder "excluded").
    
    Returns:
        dict: Ein Dictionary mit kombinierten Protein-IDs pro Schlüsselgruppe.
    """
    combined_protein_sets = {}
    
    for domain, keyword_groups in keyword_lists.items():
        for i, keyword_group in enumerate(keyword_groups):
            protein_set = set()
            for keyword in keyword_group:
                key = (keyword, domain)
                if key in csb_proteins_dict:
                    protein_set.update(csb_proteins_dict[key])
            if protein_set:  # Nur hinzufügen, wenn nicht leer
                combined_protein_sets[f"{category}{i}_{domain}"] = protein_set
    
    return combined_protein_sets


def add_query_ids_to_proteinIDset(combined_protein_sets, database_path):
    """
    Fetches proteinIDs from the SQLite database where the domain matches and genomeID is 'QUERY',
    then adds them to the corresponding protein_set.

    Args:
        combined_protein_sets (dict): Dictionary containing protein sets per domain.
        database_path (str): Path to the SQLite database.

    Returns:
        dict: Updated combined_protein_sets with additional proteinIDs from the database.
    """
    # Connect to the database
    with sqlite3.connect(database_path, timeout=120.0) as conn:
        cursor = conn.cursor()

        # **Schritt 1: Finde alle proteinIDs mit genomeID = 'QUERY'**
        cursor.execute(
            "SELECT proteinID FROM Proteins WHERE genomeID = 'QUERY'"
        )
        query_protein_ids = {row[0] for row in cursor.fetchall()}  # Set für schnellere Suche

        if not query_protein_ids:
            return combined_protein_sets  # Falls leer, sofort zurückgeben

        for key in combined_protein_sets:
            # Extract domain from key (assuming format 'categoryX_domain')
            parts = key.split("_")
            if len(parts) < 2:
                continue  # Skip malformed keys

            domain = "_".join(parts[1:])  # Reconstruct domain if it contains '_'

            # **Schritt 2: Hole proteinIDs aus Domains, aber nur, wenn sie in query_protein_ids sind**
            cursor.execute(
                f"""
                SELECT proteinID FROM Domains 
                WHERE domain = ? AND proteinID IN ({','.join(['?'] * len(query_protein_ids))})
                """,
                (domain, *query_protein_ids)
            )

            # Add fetched proteinIDs to the protein set
            protein_ids = {row[0] for row in cursor.fetchall()}
            combined_protein_sets[key].update(protein_ids)

    return combined_protein_sets



################################################################################################


def csb_granular_datasets(options, dictionary):
    #Currently not in use. Granular datasets have usually issues with false positives
    #Grouping routines to reduce redundancy in training datasets
    grouped = group_proteinID_sets(dictionary) #removed doublicates #key: frozenset of proteinIDs value: list of tuples with (keyword,domain) pairs
    
    merged, singles = merge_similar_groups(grouped,options.dbscan_epsilon,"m") # reduce redundancy by merging sets that have identical proteinIDs (like the extension of an existing csb)
    #print("\nMerged\n\n")
    #print(merged)
    #print("\nSingles\n\n")
    #print(singles)
    options.TP_merged = merged
    options.TP_singles = singles
    print_grouping_report("Merged csb groupes that include highly similar or identical proteins groups",merged,options.Csb_directory+"/Merged_csb_groups_1")
    print_grouping_report("Unique csb groupes that include distinct proteins groups",singles,options.Csb_directory+"/Singles_csb_groups_1")


################################################################################################

def decorate_training_data(options, score_limit_dict, grouped):

    # Training datasets with additional sequences
    score_limit_dict = filter_existing_faa_files(score_limit_dict, options.phylogeny_directory) # Do not fetch again for existing files
    decorated_grouped_dict = fetch_protein_ids_parallel(options.database_directory, score_limit_dict, options.cores) # get the proteinIDs within the score limits for each domain, new keys are domain only
    decorated_grouped_dict = merge_grouped_protein_ids(decorated_grouped_dict, grouped)
    fetch_seqs_to_fasta_parallel(options.database_directory, decorated_grouped_dict, options.phylogeny_directory, options.min_seqs, options.max_seqs, options.cores)
    
    return
################################################################################################




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
    print_grouping_report("Unique csb groupes that include distinct proteins groups",singles,options.Csb_directory+"/Singles_csb_groups_2")

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



def process_keyword_domains(args):
    """
    Optimierte Funktion zur Verarbeitung eines Keywords und seiner Domains mit Chunks für große `IN`-Abfragen.
    """
    database, keyword, domains, min_seqs = args
    result = {}
    chunk_size = 900  # Sicher unter dem SQLite-Limit von 999 bleiben

    with sqlite3.connect(database, timeout=120.0) as con:
        cur = con.cursor()
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA synchronous=NORMAL;')
        con.execute('PRAGMA temp_store=MEMORY;')
        con.execute('PRAGMA cache_size=-25000;')  # ca. 100MB Cache
        
        # Schritt 1: Alle Proteins.proteinID für das Keyword abrufen
        cur.execute("""
            SELECT DISTINCT Proteins.proteinID
            FROM Proteins
            INNER JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
            WHERE Keywords.keyword = ?;
        """, (keyword,))
        
        protein_ids = {row[0] for row in cur.fetchall()}  # Set für schnelleres Nachschlagen

        if not protein_ids:
            return result  # Falls keine Treffer, direkt zurückgeben

        # Schritt 2: Passende Protein-IDs mit den gewünschten Domains abrufen (in Chunks)
        query_conditions = ' OR '.join(['Domains.domain = ?'] * len(domains))

        for i in range(0, len(protein_ids), chunk_size):
            chunk = list(protein_ids)[i:i + chunk_size]  # Nimm max. 900 Protein-IDs
            
            query = f"""
                SELECT DISTINCT Proteins.proteinID, Domains.domain
                FROM Proteins
                INNER JOIN Domains ON Proteins.proteinID = Domains.proteinID
                WHERE ({query_conditions}) AND Proteins.proteinID IN ({','.join(['?'] * len(chunk))});
            """

            cur.execute(query, (*domains, *chunk))
            rows = cur.fetchall()

            # Ergebnisse in Dictionary speichern
            for proteinID, domain in rows:
                key = (keyword, domain)
                if key not in result:
                    result[key] = set()
                result[key].add(proteinID)

    return result




def fetch_proteinIDs_dict_multiprocessing(database, csb_dictionary, min_seqs, num_workers=4):
    """
    Optimized version of fetch_proteinIDs_dict using multiprocessing to
    parallelize the processing of csb_dictionary items, with ordered task printing.
    """
    # Prepare the arguments for the worker function
    tasks = [
        (database, keyword, domains, min_seqs)
        for keyword, domains in csb_dictionary.items()
        if keyword != 'default'
    ]

    total_tasks = len(tasks)

    # Use multiprocessing to process the items in parallel
    with Pool(processes=num_workers) as pool:
        results = []
        for i, result in enumerate(pool.imap(process_keyword_domains, tasks), start=1):
            # Print task start in order
            print(f"{i}/{total_tasks}", end="\r")
            results.append(result)

    # Combine results from all workers
    combined_dict = {}
    for result in results:
        for key, value in result.items():
            if key not in combined_dict:
                combined_dict[key] = set()
            combined_dict[key].update(value)

    print("\nCompleted training data selection")
    return combined_dict

    


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
    with sqlite3.connect(database, timeout=120.0) as con:
        cur = con.cursor()
        con.execute('PRAGMA journal_mode=WAL;')
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
    with sqlite3.connect(database, timeout=120.0) as con:
        cur = con.cursor()
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA synchronous=NORMAL;')
        con.execute('PRAGMA temp_store=MEMORY;')
        con.execute('PRAGMA cache_size=-25000;')  # ca. 100MB Cache

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

def fetch_seqs_to_fasta_parallel(database, dataset_dict, output_directory, min_seq, max_seq, cores=4, chunk_size=990):
    """
    Forks off the fetching of sequences for each domain using multiprocessing.

    Args:
        database (str): Path to the SQLite database.
        dataset_dict (dict): { domain: set(proteinIDs) }
        output_directory (str): Directory where FASTA files will be stored.
        min_seq (int): Minimum number of sequences required for processing.
        max_seq (int): Maximum number of sequences allowed for processing.
        cores (int): Number of parallel processes.
        chunk_size (int): Number of protein IDs to process per query batch.

    Returns:
        None
    """
    os.makedirs(output_directory, exist_ok=True)  # Ensure the output directory exists

    tasks = []
    
    for domain, protein_ids in dataset_dict.items():
        num_sequences = len(protein_ids)
        output_file = os.path.join(output_directory, f"{domain}.faa")

        # Check limits and file existence
        if num_sequences < min_seq:
            print(f"WARNING: Domain '{domain}' skipped (too few sequences: {num_sequences} < {min_seq})")
            continue  # Skip this domain

        if num_sequences > max_seq:
            print(f"WARNING: Domain '{domain}' skipped (too many sequences: {num_sequences} > {max_seq})")
            continue  # Skip this domain

        if os.path.exists(output_file):
            print(f"Skipping '{domain}', FASTA file already exists: {output_file}")
            continue  # Skip existing files

        # If all checks pass, add to tasks
        tasks.append((database, domain, protein_ids, output_directory, chunk_size))

    if not tasks:
        return  # Exit if no tasks remain

    # Use multiprocessing to run fetch_seq_to_fasta in parallel
    with multiprocessing.Pool(processes=cores) as pool:
        pool.starmap(fetch_seq_to_fasta, tasks)


        
        
def fetch_seq_to_fasta(database, domain, protein_ids, output_directory, chunk_size=990):
    """
    Fetch protein sequences from the database for a specific domain and save them into a FASTA file.
    """
    fasta_file_path = os.path.join(output_directory, f"{domain}.faa")

    if not protein_ids or os.path.exists(fasta_file_path):
        return  # Skip empty domains or if the file already exists

    with sqlite3.connect(database, timeout=120.0) as con:
        cur = con.cursor()
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA synchronous=NORMAL;')
        con.execute('PRAGMA temp_store=MEMORY;')
        con.execute('PRAGMA cache_size=-25000;')  # ca. 100MB Cache
        
        with open(fasta_file_path, "w") as fasta_file:
            protein_id_list = list(protein_ids)
            for i in range(0, len(protein_id_list), chunk_size):
                chunk = protein_id_list[i:i + chunk_size]

                if not chunk:
                    continue  # Ensure chunk is not empty
                # Convert all protein IDs to strings before passing to SQL
                query = f"""
                    SELECT proteinID, sequence FROM Proteins
                    WHERE proteinID IN ({','.join(['?'] * len(chunk))});
                """
                
                try:
                    cur.execute(query, tuple(str(protein_id) for protein_id in chunk))
                    rows = cur.fetchall()

                    # Write sequences to the domain-specific FASTA file
                    for protein_id, sequence in rows:
                        fasta_file.write(f">{protein_id}\n{sequence}\n")
                except sqlite3.InterfaceError as e:
                    print(f"SQL Error in domain {domain}: {e}")  # Debugging message

    print(f"FASTA file saved: {fasta_file_path}")




def fetch_protein_ids_for_domain(database, domain, lower_limit, upper_limit):
    """
    Fetch all protein IDs for a single domain that fall within the specified score limits.

    Args:
        database (str): Path to the SQLite database.
        domain (str): The domain name.
        lower_limit (float): Lower score limit.
        upper_limit (float): Upper score limit.

    Returns:
        tuple: (domain, set(proteinIDs))
    """
    protein_ids = set()

    with sqlite3.connect(database, timeout=120.0) as con:
        cur = con.cursor()
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA synchronous=NORMAL;')
        con.execute('PRAGMA temp_store=MEMORY;')
        con.execute('PRAGMA cache_size=-25000;')  # ca. 100MB Cache
        
        # Query to fetch protein IDs within score limits for this domain
        query = """
            SELECT DISTINCT proteinID
            FROM Domains
            WHERE domain = ?
            AND score BETWEEN ? AND ?;
        """
        cur.execute(query, (domain, lower_limit, upper_limit))
        rows = cur.fetchall()

        # Store results in a set
        protein_ids = {row[0] for row in rows}

    return domain, protein_ids


def filter_existing_faa_files(domain_dict, directory):
    """
    Entfernt Einträge aus `domain_dict`, falls bereits eine .faa-Datei 
    für die jeweilige Domain im `directory` existiert.

    Args:
        domain_dict (dict): Dictionary mit Domains als Keys.
        directory (str): Verzeichnis, in dem nach bestehenden .faa-Dateien gesucht wird.

    Returns:
        dict: Gefiltertes Dictionary ohne Domains, die bereits .faa-Dateien haben.
    """
    # Hole eine Liste aller existierenden .faa-Dateien im Verzeichnis
    existing_files = {f for f in os.listdir(directory) if f.endswith(".faa")}

    # Bereinige das Dictionary: Entferne Domains mit vorhandener .faa-Datei
    filtered_dict = {domain: values for domain, values in domain_dict.items() 
                     if f"{domain}.faa" not in existing_files}

    return filtered_dict
    
    
def fetch_protein_ids_parallel(database, score_limit_dict, cores):
    """
    Fetch all protein IDs per domain in parallel using multiprocessing.

    Args:
        database (str): Path to the SQLite database.
        score_limit_dict (dict): { domain: {'lower_limit': X, 'upper_limit': Y} }

    Returns:
        dict: { domain: set(proteinIDs) }
    """
    tasks = [
        (database, domain, limits["lower_limit"], limits["upper_limit"])
        for domain, limits in score_limit_dict.items()
    ]

    # Use multiprocessing to run fetch_protein_ids_for_domain in parallel
    with multiprocessing.Pool(processes=cores) as pool:
        results = pool.starmap(fetch_protein_ids_for_domain, tasks)

    # Combine results into a dictionary
    domain_protein_ids = {domain: protein_ids for domain, protein_ids in results}
    #print(domain_protein_ids.keys())
    return domain_protein_ids


def merge_grouped_protein_ids(protein_ids_by_domain, grouped_dict):

    # Merge the grp0 seqIDs and the decorate seqIDs into a new dictionary which has only domain without prefix
    # this value set is used for the fetching of the sequences to a file for the further classification steps
    merged = protein_ids_by_domain.copy()
    for grouped_key, protein_ids in grouped_dict.items():
        domain = grouped_key.replace("grp0_", "", 1)
        if domain in merged:
            merged[domain].update(protein_ids)
        else:
            merged[domain] = set(protein_ids)

    return {k: v for k, v in merged.items() if not k.startswith("grp0_")}











def fetch_protein_to_fasta(options, grouped, prefix=""):
    #Currently not in use
    """
    Fetches protein sequences from the Proteins table for protein IDs in the grouped sets,
    and writes the sequences directly to a FASTA file.

    Args:
        options: Object containing configuration options, including the database directory.
        grouped: Dictionary where the keys are frozensets of protein IDs and the values are lists of keyword-domain pairs.

    Returns:
        None, but writes the sequences to a FASTA file.
    """
    with sqlite3.connect(options.database_directory, timeout=120.0) as con:
        cursor = con.cursor()
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA synchronous=NORMAL;')
        con.execute('PRAGMA temp_store=MEMORY;')
        con.execute('PRAGMA cache_size=-25000;')  # ca. 100MB Cache

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

    
    
def fetch_domains_superfamily_to_fasta(options, directory):
    """
    Fetches all proteins associated with each domain from the database and writes them to FASTA files.

    Parameters:
        options (object): Contains the database connection and other options.
        directory (str): Target directory for the FASTA files.

    Returns:
        dict: Mapping of domain names to the generated FASTA file paths.
    """
    output_files = {}

    with sqlite3.connect(options.database_directory, timeout=120.0) as con:
        cur = con.cursor()
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA synchronous=NORMAL;')
        con.execute('PRAGMA temp_store=MEMORY;')
        con.execute('PRAGMA cache_size=-25000;')  # ca. 100MB Cache



        # Retrieve all distinct domains from the database
        query = "SELECT DISTINCT domain FROM Domains"
        cur.execute(query)
        domains = [row[0] for row in cur.fetchall()]

        for domain in domains:
            # Create the output filename for the domain
            name = domain.split('_')[-1]
            filename = f"superfamily_{name}.faa"
            output_fasta_path = os.path.join(directory, filename)

            # Skip if the file already exists
            if os.path.isfile(output_fasta_path):
                print(f"Skipping. File {output_fasta_path} already exists.")
                continue
            else:
                print(f"{output_fasta_path} did not exist")

            # Fetch all proteins associated with this domain
            query = '''
                SELECT DISTINCT P.proteinID, P.sequence
                FROM Proteins P
                INNER JOIN Domains D ON P.proteinID = D.proteinID
                WHERE D.domain = ?
            '''
            cur.execute(query, (domain,))
            proteins = cur.fetchall()

            # Write all sequences to the FASTA file
            with open(output_fasta_path, 'w') as fasta_file:
                for proteinID, sequence in proteins:
                    fasta_file.write(f'>{proteinID}\n{sequence}\n')

            # Add the file to the output dictionary
            output_files[name] = output_fasta_path

    return output_files






def fetch_all_proteins(database, filepath):
    """
    Fetches the results from the SQLite database based on domain and keyword, and stores them in a FASTA file.

    Args:
        database: Name of the SQLite database.
        results_dict: Dictionary containing domain as keys and keyword as values.
        output_directory: Directory path to store the FASTA file.
        
        Output is a files with dir/sequence.faa
    """
    if os.path.isfile(filepath):
        return filepath
    with sqlite3.connect(database, timeout=120.0) as con:
        cur = con.cursor()
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA journal_mode=WAL;')
        con.execute('PRAGMA synchronous=NORMAL;')
        con.execute('PRAGMA temp_store=MEMORY;')
        con.execute('PRAGMA cache_size=-25000;')  # ca. 100MB Cache
        
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
    
    
    
    
##################################################################################################
    
    
def save_cache(directory, filename, data):
    """
    Speichert Daten in einer Datei im Pickle-Format.
    """
    with open(os.path.join(directory, filename), "wb") as f:
        pickle.dump(data, f)

def load_cache(directory, filename):
    """
    Lädt zwischengespeicherte Daten aus einer Datei, falls vorhanden.
    """
    file_path = os.path.join(directory, filename)
    if os.path.exists(file_path):
        print("Loading from cache")
        with open(file_path, "rb") as f:
            return pickle.load(f)
    return None
    
    
    

    
