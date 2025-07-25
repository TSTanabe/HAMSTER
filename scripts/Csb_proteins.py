#!/usr/bin/python
import sqlite3
import os
import multiprocessing
import time
import threading
from multiprocessing import Pool, Manager, Value, Lock
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
import itertools
from typing import Dict, Set, Tuple, List, Any

from . import Csb_statistic
from . import myUtil

logger = myUtil.logger

#first get the csb with their identifiers. make sets of appearence as value to key name of the csb
#do not compare the csb appearences as they jaccard itself should have managed it.

def prepare_csb_grouped_training_proteins(options: Any) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Set[str]]]:
    """
    Prepares grouped protein sequences for training by analyzing gene clusters (CSBs)
    and extracting associated proteins with conserved genomic context.

    This function:
    1. Computes score limits and cluster keyword groupings.
    2. Extracts training protein sequences.
    3. Loads or computes grouped training sets from cache.
    4. Writes grouped sequences to FASTA files if needed.

    Args:
        options (argparse.Namespace): The runtime configuration.

    Returns:
        tuple:
            - grp_score_limit_dict (dict): Domain-wise score limits.
            - grouped (dict): Mapping from group protein label → set of protein IDs.
    """

    # Step 1: Try loading cached grouped training data
    grouped = myUtil.load_cache(options, 'grp0_training_proteinIDs.pkl') # Name is defined by fetch_to_fasta routine
    grp_score_limit_dict = myUtil.load_cache(options, 'grp0_score_limit_dict.pkl')

    if grouped and grp_score_limit_dict:
        return grp_score_limit_dict, grouped
        
    # Step 2: Compute score limits and keyword clusters
    grp_score_limit_dict, _, grouped_keywords_dict, clustered_excluded_keywords_dict = Csb_statistic.group_gene_cluster_statistic(options)
    
    logger.info("Collecting sequences for training datasets with similar csb")
    csb_proteins_dict = csb_proteins_datasets(options, grouped_keywords_dict)

    logger.info("Processing highly similar homologs with specific genomic context")

    # Step 3: Export one fasta per protein
    grouped = csb_proteins_datasets_combine(grouped_keywords_dict, csb_proteins_dict)
    grouped = add_query_ids_to_proteinIDset(grouped, options.database_directory)

    # Step 4: Save in pkl cache
    myUtil.save_cache(options, 'grp0_training_proteinIDs.pkl', grouped)
    myUtil.save_cache(options, 'grp0_score_limit_dict.pkl', grp_score_limit_dict)
    	
    return grp_score_limit_dict, grouped
    
    
    
def csb_proteins_datasets(options: Any, grouped_keywords_dict: Dict) -> Dict[Tuple[str, str], Set[str]]:
    """
    Fetches all proteinIDs for every (keyword, domain) from the csb clusters.
    First get all csb from the CSB output file

    Args:
        options (Namespace): Configuration.
        grouped_keywords_dict (dict): dict: { domain: [[grouped_keywords]] }

    Returns:
        dict: {(keyword, domain): set(proteinIDs)}. Example: {("ClusterA", "ABC_trans"): {"prot1", "prot2"}}
    """
 
    # Load existing data
    dictionary = myUtil.load_cache(options, "csb_protein_dataset.pkl")
    
    if dictionary:
        logger.debug("Loaded existing reference protein sequence dataset from cache: csb_protein_dataset.pkl")
        dictionary = filter_dictionary_by_inclusion_domains(dictionary, options.include_list)
        dictionary = filter_dictionary_by_excluding_domains(dictionary, options.exclude_list)
        return dictionary  
        
    ## If not precomputed then compute the csb to select and fetch die proteinIDs from the dataset         
    # 1. Get the domain types as set by the selection statistic module for the csb finder algorithm
    csb_dictionary = parse_csb_file_to_dict(options.csb_output_file) # dictionary with cbs_name => csb items

    # Determine the csb that shall be kept.
    csbs_to_keep = set()
    for domain, group_lists in grouped_keywords_dict.items():
        for group in group_lists:
            csbs_to_keep.update(group)

    # Now, filter the dictionary to keep only the wanted CSBs
    csb_dictionary = {csb: val for csb, val in csb_dictionary.items() if csb in csbs_to_keep}

    
    # 2. add the static csb from the pattern file
    pattern_dictionary = parse_csb_file_to_dict(options.patterns_file)  # Fetch the ones that are in the pattern file
    csb_dictionary = {**csb_dictionary, **pattern_dictionary}
    options.csb_dictionary = csb_dictionary # save the patterns for later use

    #Fetch for each csb id all the domains in the csb that are query domains
    #dictionary is: dict[(keyword, domain)] => set(proteinIDs)
    myUtil.save_cache(options, "csb_selected_csb_to_domain_dataset.pkl", csb_dictionary)

    logger.info("Fetching the protein sequence identifiers from local database")
    dictionary = fetch_proteinIDs_dict_multiprocessing(options.database_directory,csb_dictionary,options.min_seqs,options.cores)
    dictionary = remove_non_query_clusters(options.database_directory, dictionary) #delete all that are not in accordance with query

    myUtil.save_cache(options, "csb_protein_dataset.pkl", dictionary)
    
    # Remove domains that are excluded by user options
    dictionary = filter_dictionary_by_inclusion_domains(dictionary, options.include_list)
    dictionary = filter_dictionary_by_excluding_domains(dictionary, options.exclude_list)

    return dictionary

def generate_score_limit_dict_from_grouped(
    database_path: str,
    grouped_dict: Dict[str, Set[str]],
    default_lower: int = 100,
    default_upper: int = 2000,
    chunk_size: int = 999
) -> Dict[str, Dict[str, float]]:
    """
    Generate a score limit dictionary for each domain from grouped proteinIDs.

    Args:
        database_path (str): SQLite DB path.
        grouped_dict (dict): {domain: set(proteinIDs)}
        default_lower (int): Fallback if no data.
        default_upper (int): Fallback if no data.
        chunk_size (int): Max number of SQL params.

    Returns:
        dict: {domain: {'lower_limit': X, 'upper_limit': Y}}
    """
    result = {}

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        for domain, protein_ids in grouped_dict.items():
            if not protein_ids:
                result[domain] = {'lower_limit': default_lower, 'upper_limit': default_upper}
                continue

            min_score = float('inf')
            max_score = float('-inf')
            protein_id_list = list(protein_ids)

            for i in range(0, len(protein_id_list), chunk_size):
                chunk = protein_id_list[i:i + chunk_size]
                placeholders = ",".join("?" for _ in chunk)
                query = f"""
                    SELECT MIN(score), MAX(score)
                    FROM Domains
                    WHERE domain = ? AND proteinID IN ({placeholders});
                """
                cur.execute(query, (domain, *chunk))
                row = cur.fetchone()
                if row:
                    chunk_min, chunk_max = row
                    if chunk_min is not None:
                        min_score = min(min_score, chunk_min)
                    if chunk_max is not None:
                        max_score = max(max_score, chunk_max)

            if min_score == float('inf') or max_score == float('-inf'):
                # No valid score data found
                result[domain] = {'lower_limit': default_lower, 'upper_limit': default_upper}
            else:
                result[domain] = {'lower_limit': min_score, 'upper_limit': max_score}

    return result

################################################################################################

def csb_proteins_datasets_combine(
    keyword_lists: List[str],
    csb_proteins_dict: Dict[str, Set[str]]
) -> Dict[str, Set[str]]:
    """
    Combines the grouped keywords dict and csb proteins dict to create final grouped training sets.

    Args:
        grouped_keywords_dict (dict): {group: [keyword1, keyword2, ...]}
        csb_proteins_dict (dict): {keyword: set(proteinIDs)}

    Returns:
        dict: {group: set(proteinIDs)}

    Example:
        grouped_keywords_dict = {'ABC': ['a','b']}
        csb_proteins_dict = {'a': {'prot1'}, 'b': {'prot2'}}
        Output: {'ABC': {'prot1','prot2'}}
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
                combined_protein_sets[domain] = protein_set
    
    return combined_protein_sets


def add_query_ids_to_proteinIDset(combined_protein_sets: Dict[str, Set[str]], database_path: str) -> Dict[str, Set[str]]:
    """
    Adds any proteinIDs from the query cluster to each protein group set.

    Args:
        grouped (dict): {group: set(proteinIDs)}
        database_path (str): Path to SQLite DB.

    Returns:
        dict: {group: set(proteinIDs)} (now includes query IDs)
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

            # **Schritt 2: Hole proteinIDs aus Domains, aber nur, wenn sie in query_protein_ids sind**
            cursor.execute(
                f"""
                SELECT proteinID FROM Domains 
                WHERE domain = ? AND proteinID IN ({','.join(['?'] * len(query_protein_ids))})
                """,
                (key, *query_protein_ids)
            )

            # Add fetched proteinIDs to the protein set
            protein_ids = {row[0] for row in cursor.fetchall()}
            combined_protein_sets[key].update(protein_ids)

    return combined_protein_sets



################################################################################################

def fetch_protein_family_sequences(options, directory, score_limit_dict, grouped):

    # Training datasets with additional sequences
    score_limit_dict = filter_existing_faa_files(score_limit_dict, directory) # Do not fetch again for existing files
    decorated_grouped_dict = fetch_protein_ids_parallel(options.database_directory, score_limit_dict, options.cores, options.max_seqs) # get the proteinIDs within the score limits for each domain, new keys are domain only
    decorated_grouped_dict = merge_grouped_protein_ids(decorated_grouped_dict, grouped)
    fetch_seqs_to_fasta_parallel(options.database_directory, decorated_grouped_dict, directory, options.min_seqs, options.max_seqs, options.cores)
    
    return
################################################################################################

    
def parse_csb_file_to_dict(filepath: str) -> Dict[str, List[str]]:
    """
    Parses a CSB output or patterns file into a dict.

    Args:
        filepath (str): Path to CSB or patterns file.

    Returns:
        dict: {key: [domains, ...]}

    Example:
        "csb1\tABC\tDEF\n" → {'csb1': ['ABC','DEF']}
    """
    data_dict = {}

    if filepath is None or not os.path.isfile(filepath):
        return data_dict

    with open(filepath, 'r') as file:
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



def process_keyword_domains(args: Tuple[str, str, List[str], int]) -> Dict[Tuple[str, str], Set[str]]:
    """
    Efficiently process a keyword and its domains, using chunked SQL queries for large IN clauses.

    Args:
        args (tuple): (database_path, keyword, domains, min_seqs)
            - database_path (str): Path to SQLite DB.
            - keyword (str): Cluster keyword.
            - domains (list): List of domains.
            - min_seqs (int): Minimum number of sequences.

    Returns:
        dict: {(keyword, domain): set(proteinIDs)}

    Example:
        ('db.sqlite', 'kword', ['ABC', 'DEF'], 5) 
        → {('kword','ABC'): {'prot1','prot2'}}
    """
    database, keyword, domains, min_seqs = args
    result = {}
    chunk_size = 900  # Sicher unter dem SQLite-Limit von 999 bleiben
    logger.debug(f"Selecting proteinIDs for {keyword} {domains}")
    
    with sqlite3.connect(database, timeout=120.0) as con:
        cur = con.cursor()
        
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




def fetch_proteinIDs_dict_multiprocessing(
    database_path: str,
    csb_dictionary: Dict[str, List[str]],
    min_seqs: int,
    cores: int
) -> Dict[Tuple[str, str], Set[str]]:
    """
    Multiprocessing wrapper for process_keyword_domains.

    Args:
        database_path (str): SQLite DB path.
        csb_dictionary (dict): {keyword: [domains]}
        min_seqs (int): Minimum sequence count to keep.
        cores (int): Number of parallel processes.

    Returns:
        dict: {(keyword, domain): set(proteinIDs)}

    Example Output:
        {('A','DEF'): {'prot1', ...}}
    """
    # Prepare the arguments for the worker function
    tasks = [
        (database_path, keyword, domains, min_seqs)
        for keyword, domains in csb_dictionary.items()
        if keyword != 'default'
    ]

    total_tasks = len(tasks)

    # Use multiprocessing to process the items in parallel
    with Pool(processes=cores) as pool:
        results = []
        for i, result in enumerate(pool.imap(process_keyword_domains, tasks), start=1):
            results.append(result)

    # Combine results from all workers
    combined_dict = {}
    for result in results:
        for key, value in result.items():
            if key not in combined_dict:
                combined_dict[key] = set()
            combined_dict[key].update(value)

    return combined_dict

    


######################################################################################################

    
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





###############################################################################
#################### Protein to fasta operations ##############################
###############################################################################

def fetch_training_data_to_fasta(options: Any, grouped: Dict[str, Set[str]], prefix: str) -> None:
    """
    Writes training FASTA files per group/domain, using a prefix.

    Args:
        options (Namespace): Configuration/options.
        grouped (dict): {group: set(proteinIDs)}
        prefix (str): Prefix for each group name (output FASTA: {prefix}_group.faa)

    Returns:
        None (writes FASTA files).
    """

    extended_grouped_prefixed = {f"{prefix}_{key}": value for key, value in grouped.items()} # Extend dictionary with a prefix
    
    fetch_seqs_to_fasta_parallel(
        options.database_directory,
        extended_grouped_prefixed,
        options.fasta_output_directory,
        min_seq=options.min_seqs,
        max_seq=options.max_seqs,
        cores=options.cores,
        hardcap=options.hardcap
    )

def fetch_seqs_to_fasta_parallel(
    database: str,
    dataset_dict: Dict[str, Set[str]],
    output_directory: str,
    min_seq: int,
    max_seq: int,
    cores: int = 4,
    chunk_size: int = 990,
    hardcap: int = 5000
) -> None:
    """
    Write each protein family in dataset_dict to a FASTA file using multiprocessing.

    Args:
        database (str): SQLite DB.
        dataset_dict (dict): {domain: set(proteinIDs)}
        output_directory (str): Directory for FASTA output.
        min_seq (int): Minimum seqs required.
        max_seq (int): Maximum seqs allowed.
        cores (int): Multiprocessing processes.
        chunk_size (int): SQL batch.
        hardcap (int): Absolute max allowed (protect from runaway).

    Returns:
        None (writes FASTA).
    """
    os.makedirs(output_directory, exist_ok=True)  # Ensure the output directory exists

    tasks = []
    
    for domain, protein_ids in dataset_dict.items():
        num_sequences = len(protein_ids)
        output_file = os.path.join(output_directory, f"{domain}.faa")

        # Check limits and file existence
        if num_sequences < min_seq:
            logger.warn(f"'{domain}' was skipped due to too few sequences: {num_sequences} < {min_seq}")
            continue  # Skip this domain

        if num_sequences > (max_seq + hardcap):
            logger.warn(f"'{domain}' was skipped due to too many sequences {num_sequences} > {max_seq+hardcap})")
            continue  # Skip this domain

        if os.path.exists(output_file):
            logger.debug(f"Skipped '{domain}'. FASTA file already exists: {output_file}")
            continue  # Skip existing files

        # If all checks pass, add to tasks
        tasks.append((database, domain, protein_ids, output_directory, chunk_size))

    if not tasks:
        return  # Exit if no tasks remain

    # Use multiprocessing to run fetch_seq_to_fasta in parallel
    with multiprocessing.Pool(processes=cores) as pool:
        pool.starmap(fetch_seq_to_fasta, tasks)



def connect_readonly(db_path):
    # helper routine for the fetch seq to fasta database connection
    uri = f"file:{db_path}?mode=ro&immutable=1"
    return sqlite3.connect(uri, uri=True, check_same_thread=False)

def fetch_seq_to_fasta(
    database: str,
    domain: str,
    protein_ids: Set[str],
    output_directory: str,
    chunk_size: int = 990
) -> None:
    """
    Fetch sequences for a domain and write to a FASTA file.

    Args:
        database (str): SQLite DB.
        domain (str): Domain/family name.
        protein_ids (set): Protein IDs.
        output_directory (str): Where to save.
        chunk_size (int): SQL chunk size.

    Returns:
        None (writes FASTA).
    """
    fasta_file_path = os.path.join(output_directory, f"{domain}.faa")

    if not protein_ids or os.path.exists(fasta_file_path):
        return  # Skip empty domains or if the file already exists


    with connect_readonly(database) as con:
        cur = con.cursor()
        
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
                    logger.error(f"SQL Error in domain {domain}: {e}")  # Debugging message

    logger.debug(f"FASTA file saved: {fasta_file_path}")




def fetch_protein_ids_for_domain(
    database: str,
    domain: str,
    lower_limit: float,
    upper_limit: float,
    max_count: int = 50000
) -> Tuple[str, Set[str]]:
    """
    Efficiently fetches the top N proteinIDs for a domain within a score range, including all ties.

    Args:
        database (str): SQLite DB.
        domain (str): Domain.
        lower_limit (float): Lower threshold.
        upper_limit (float): Upper threshold.
        max_count (int): Top N (includes ties at boundary).

    Returns:
        tuple: (domain, set(proteinIDs))

    Example Output:
        ('ABC', {'prot1','prot2'})
    """

    with sqlite3.connect(database) as con:
        cur = con.cursor()

        query = """
            SELECT proteinID, score
            FROM Domains
            WHERE domain = ?
              AND score BETWEEN ? AND ?
            ORDER BY score DESC;
        """
        cur.execute(query, (domain, lower_limit, upper_limit))
        rows = cur.fetchall()

    protein_ids = set()
    count = 0
    last_score = None
    tie_count = 0

    for protein_id, score in rows:
        if count < max_count:
            protein_ids.add(protein_id)
            last_score = score
            count += 1
        elif score == last_score:
            protein_ids.add(protein_id)
            tie_count += 1
        else:
            break

    if tie_count > 0:
        logger.warning(f"max_count {max_count} reached for {domain}, but {tie_count} more proteins with same score ({last_score}) included.")

    return domain, protein_ids



def filter_existing_faa_files(domain_dict: Dict[str, Any], directory: str) -> Dict[str, Any]:
    """
    Removes entries from domain_dict if a .faa file already exists.

    Args:
        domain_dict (dict): {domain: ...}
        directory (str): Directory to check.

    Returns:
        dict: Only domains with no .faa in directory.

    Example:
        {'A':{...}}, 'dir/' (and 'dir/A.faa' exists) → {}
    """
    # Hole eine Liste aller existierenden .faa-Dateien im Verzeichnis
    existing_files = {f for f in os.listdir(directory) if f.endswith(".faa")}

    # Bereinige das Dictionary: Entferne Domains mit vorhandener .faa-Datei
    filtered_dict = {domain: values for domain, values in domain_dict.items() 
                     if f"{domain}.faa" not in existing_files}

    return filtered_dict
    
    
def fetch_protein_ids_parallel(database, score_limit_dict, cores, max_seqs=50000):
    """
    Fetch all protein IDs per domain in parallel using multiprocessing.

    Args:
        database (str): Path to the SQLite database.
        score_limit_dict (dict): { domain: {'lower_limit': X, 'upper_limit': Y} }

    Returns:
        dict: { domain: set(proteinIDs) }
    """
    tasks = [
        (database, domain, limits["lower_limit"], limits["upper_limit"], max_seqs)
        for domain, limits in score_limit_dict.items()
    ]

    # Use multiprocessing to run fetch_protein_ids_for_domain in parallel. This routine fetches seqs up to a limiter and then all seqs with the same score
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
                logger.debug(f"File {output_fasta_path} already exists - skipping.")
                continue

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
    
    

def filter_dictionary_by_inclusion_domains(dictionary: Dict, include_list: List[str]) -> Dict:
    """
    Filters the dictionary to include only domains listed in include_list.

    Args:
        dictionary (dict): {domain: ...}
        include_list (list): List of allowed domains.

    Returns:
        dict: Filtered dictionary.

    Example:
        dictionary = {'ABC': {...}, 'DEF': {...}}
        include_list = ['ABC']
        Output: {'ABC': {...}}
    """
    if not include_list:
        return dictionary
    return {k: v for k, v in dictionary.items() if k in include_list}


def filter_dictionary_by_excluding_domains(dictionary: Dict, exclude_list: List[str]) -> Dict:
    """
    Filters out domains listed in exclude_list.

    Args:
        dictionary (dict): {domain: ...}
        exclude_list (list): List of forbidden domains.

    Returns:
        dict: Filtered dictionary.

    Example:
        dictionary = {'ABC': {...}, 'DEF': {...}}
        exclude_list = ['DEF']
        Output: {'ABC': {...}}
    """
    if not exclude_list:
        return dictionary
    return {k: v for k, v in dictionary.items() if k not in exclude_list}

def remove_non_query_clusters(database_path: str, dictionary: Dict) -> Dict:
    """
    Remove all clusters that only contain query protein IDs (no 'real' proteins).

    Args:
        database_path (str): Path to SQLite DB.
        dictionary (dict): {domain: set(proteinIDs)}

    Returns:
        dict: Filtered dictionary.
    """
    try:
        with sqlite3.connect(database_path) as con:
            cur = con.cursor()
            cur.execute("SELECT proteinID FROM Proteins WHERE genomeID != 'QUERY'")
            non_query_proteins = {row[0] for row in cur.fetchall()}
        return {k: v & non_query_proteins for k, v in dictionary.items() if v & non_query_proteins}
    except Exception as e:
        logger.warning(f"Could not filter non-query clusters: {e}")
        return dictionary
    
    
    

##########################################################################################
###################### Get genomic context per protein ###################################
##########################################################################################

_thread_local = threading.local()

def _get_thread_cursor(database_path: str) -> sqlite3.Cursor:
    """
    Opens or reuses a read-only SQLite cursor in a thread, with proper PRAGMAs set.

    Args:
        database_path (str): Path to SQLite DB.

    Returns:
        sqlite3.Cursor: SQLite cursor, thread-local, always read-only.

    Example:
        cur = _get_thread_cursor("db.sqlite")
    """
    if getattr(_thread_local, "cur", None) is None:
        # 1) Open connection as immutable, no lockfile, no writes
        uri = f"file:{database_path}?mode=ro&immutable=1&nolockfile=1"
        con = sqlite3.connect(uri, uri=True, check_same_thread=False)
        con.row_factory = sqlite3.Row
        cur = con.cursor()
        cur.executescript("""
            PRAGMA query_only = TRUE;
            PRAGMA journal_mode = OFF;
            PRAGMA temp_store = MEMORY;
            PRAGMA defer_foreign_keys = TRUE;
        """)
        _thread_local.con = con
        _thread_local.cur = cur
    return _thread_local.cur

def _fetch_neighbours_chunk(args: Tuple[str, List[str]]) -> Tuple[Dict[str, List[List[str]]], Dict[Tuple, int]]:
    """
    Worker for a chunk of protein IDs: fetch cluster, neighbours, and their domains.

    Args:
        args (tuple): (db_path, protein_ids)

    Returns:
        tuple:
          - neighbours: {proteinID: [ [domain1,...], ... ] } (domains for all neighbours in cluster)
          - comp_counts: {composition tuple: count}

    Example:
        ("db.sqlite", ["p1","p2"]) ->
          ({"p1":[["ABC","DEF"],["DEF"]], ...}, {("ABC", "DEF"): 2, ...})
    """
    db_path, protein_ids = args
    neighbours: Dict[str, List[List[str]]] = {}
    comp_counts: DefaultDict[Tuple, int] = defaultdict(int)
    if not protein_ids:
        return {}, {}

    t0 = time.perf_counter()
    cur = _get_thread_cursor(db_path)

    # 1) Map hit → clusterID
    ph1 = ",".join("?" for _ in protein_ids)
    cur.execute(
        f"SELECT proteinID, clusterID FROM Proteins WHERE proteinID IN ({ph1})",
        protein_ids
    )
    hit2cid = {r["proteinID"]: r["clusterID"] for r in cur.fetchall()}
    if not hit2cid:
        dt = time.perf_counter() - t0
        logger.debug(f"[TIMING] chunk {len(protein_ids)} (no clusters) in {dt:.3f}s")
        return {pid: [["singleton"]] for pid in protein_ids}, {}

    # 2) clusterID → all neighbour proteins (sorted)
    cids = list(set(hit2cid.values()))
    ph2 = ",".join("?" for _ in cids)
    cur.execute(
        f"SELECT proteinID, clusterID, start FROM Proteins WHERE clusterID IN ({ph2})",
        cids
    )
    nbr_rows = list(cur.fetchall())
    nbr_rows.sort(key=lambda r: (r["clusterID"], r["start"]))

    # 3) neighbour protein → domains
    nids = [r["proteinID"] for r in nbr_rows]
    ph3 = ",".join("?" for _ in nids)
    cur.execute(
        f"SELECT proteinID, domain FROM Domains WHERE proteinID IN ({ph3})",
        nids
    )
    dom_rows = cur.fetchall()

    dt = time.perf_counter() - t0
    logger.debug(f"[TIMING] chunk {len(protein_ids)} → {len(nids)} neighbours in {dt:.3f}s")

    # Aggregation
    domains_by_nid: DefaultDict[str, List[str]] = defaultdict(list)
    for r in dom_rows:
        domains_by_nid[r["proteinID"]].append(r["domain"] or "no_neighbours")

    cluster_to_nids: DefaultDict[Any, List[str]] = defaultdict(list)
    for r in nbr_rows:
        cluster_to_nids[r["clusterID"]].append(r["proteinID"])

    # For each input hit: all its neighbours' domains
    neighbours = {
        hit: [
            domains_by_nid.get(nid, ["singleton"])
            for nid in cluster_to_nids.get(cid, [])
        ] or [["singleton"]]
        for hit, cid in hit2cid.items()
    }

    # Count unique compositions per cluster
    for cid, ids in cluster_to_nids.items():
        comp = tuple(
            sorted(tuple(sorted(domains_by_nid.get(nid, ["singleton"])))
                   for nid in ids)
        )
        comp_counts[comp] += 1

    return neighbours, comp_counts

def fetch_neighbouring_genes_with_domains(
    db_path: str,
    protein_ids: Set[str],
    chunk_size: int = 999,
    threads: int = 8
) -> Tuple[Dict[str, List[List[str]]], List[Tuple[Tuple, int]]]:
    """
    Threaded: process protein_ids in chunks and get all their cluster neighbours' domain arrays.

    Args:
        db_path (str): Path to SQLite DB.
        protein_ids (set): Input protein IDs to analyze.
        chunk_size (int): How many IDs per thread.
        threads (int): Thread pool size.

    Returns:
        tuple:
            - full_neighbours: {proteinID: [ [domain,...], ... ] }
            - sorted_compositions: [ (composition tuple, count), ... ] (sorted by count)

    Example:
        full_neigh, comp = fetch_neighbouring_genes_with_domains("db.sqlite", {"p1","p2"}, 500, 8)
    """
    protein_list = list(protein_ids)
    if not protein_ids:
        return {}, []

    chunks = [
        protein_list[i : i + chunk_size]
        for i in range(0, len(protein_ids), chunk_size)
    ]
    args = [(db_path, chunk) for chunk in chunks]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        partials = list(executor.map(_fetch_neighbours_chunk, args))

    full_neighbours: Dict[str, List[List[str]]] = {}
    total_counts: DefaultDict[Tuple, int] = defaultdict(int)
    for neigh_dict, comp_dict in partials:
        full_neighbours.update(neigh_dict)
        for comp, cnt in comp_dict.items():
            total_counts[comp] += cnt

    sorted_compositions = sorted(
        total_counts.items(), key=lambda x: x[1], reverse=True
    )
    return full_neighbours, sorted_compositions
    
    
    
##########################################################################################################
######################### Extend grouped reference proteins with similar csb #############################
##########################################################################################################

def extend_merged_grouped_by_csb_similarity(
    options: Any,
    grouped: Dict[str, Set[str]]
) -> Dict[str, Set[str]]:
    """
    Main routine: Extends grouped protein sets by including proteins from highly similar CSB patterns.

    Args:
        options: Options/config object, must have csb_output_file, jaccard, sqlite_chunks, etc.
        grouped (dict): {domain: set(proteinIDs)} (starting protein sets).

    Returns:
        dict: Updated grouped dict with new proteins added.

    Example:
        {'ABC': {'p1', 'p2'}} → {'ABC': {'p1', 'p2', 'p3', 'p4'}}
    """    
    # Main routine for proteins from highly similar csb to the basic csb
    
    # Find keywords that are in jaccard distance to the csb of reference sequences in grouped
    protein_to_new_keywords_dict = select_similar_csb_patterns_per_protein(options, grouped, options.jaccard)
    
    # Integrate proteins with the new keywords into merged_grouped dataset
    extended_grouped = integrate_csb_variants_into_merged_grouped(options, grouped, protein_to_new_keywords_dict, options.sqlite_chunks)

    return extended_grouped







def select_similar_csb_patterns_per_protein(
    options: Any,
    merged_grouped: Dict[str, Set[str]],
    jaccard_threshold: float = 0.7
) -> Dict[str, Set[str]]:
    """
    For each domain, finds all CSB patterns with high Jaccard similarity to those associated with the proteins in that domain.

    Args:
        options: Options/config object.
        merged_grouped (dict): {domain: set(proteinIDs)}
        jaccard_threshold (float): Similarity threshold.

    Returns:
        dict: {domain: set(similar CSB keywords)}

    Example:
        {'ABC': {'p1', 'p2'}} → {'ABC': {'csb1', 'csb2'}}
    """

    # Stores the similar CSB keywords for each domain
    jaccard_included_patterns = {}  # { domain : set(csb_keywords) }

    # Load all CSB patterns from csb_output_file
    csb_dictionary = parse_csb_file_to_dict(options.csb_output_file)

    # Process each domain
    for domain, protein_ids in merged_grouped.items():
        logger.debug(f"Processing domain {domain}")

        # Fetch all keywords associated with proteins of this domain
        all_keywords = fetch_keywords_for_proteins(options.database_directory, protein_ids)

        if not all_keywords:
            logger.warning(f"No keywords found for domain {domain}")
            continue

        # Build union of patterns from these keywords
        domain_pattern_union = set().union(
            *[csb_dictionary[k] for k in all_keywords if k in csb_dictionary]
        )

        logger.debug(f"Domain {domain}: {len(domain_pattern_union)} unique CSB encode reference sequences.")

        # Compare against all CSB patterns and collect similar ones
        similar_csb_keywords = set()

        for csb_key, csb_pattern in csb_dictionary.items():
            intersection = len(domain_pattern_union & csb_pattern)
            union = len(domain_pattern_union | csb_pattern)
            similarity = intersection / union if union > 0 else 0.0

            if similarity >= jaccard_threshold:
                similar_csb_keywords.add(csb_key)

        logger.debug(f"Domain {domain}: {len(similar_csb_keywords)} similar CSB patterns (Jaccard >= {jaccard_threshold})")
        
        # Store the result for this domain
        jaccard_included_patterns[domain] = similar_csb_keywords

    # Return dictionary: domain → set of similar CSB keywords
    return jaccard_included_patterns



def fetch_keywords_for_proteins(
    database_path: str,
    protein_ids: Set[str],
    chunk_size: int = 999
) -> Set[str]:
    """
    Fetches all CSB keywords (cluster identifiers) for a set of protein IDs.

    Args:
        database_path (str): Path to SQLite database.
        protein_ids (set): Protein IDs to fetch keywords for.
        chunk_size (int): DB chunk size.

    Returns:
        set: All unique keywords.

    Example:
        {'p1','p2'} → {'csb1', 'csb2'}
    """

    protein_to_keywords = defaultdict(set)

    if not protein_ids:
        logger.warning("No proteinIDs provided to fetch_keywords_for_proteins.")
        return set()

    conn = sqlite3.connect(database_path)
    cur = conn.cursor()

    protein_ids = list(protein_ids)
    total = len(protein_ids)
    #print(f"[INFO] Fetching keywords for {total} proteins")

    for start in range(0, total, chunk_size):
        end = start + chunk_size
        chunk = protein_ids[start:end]

        protein_placeholders = ','.join(['?'] * len(chunk))
        query = f"""
        SELECT Proteins.proteinID, Keywords.keyword
        FROM Proteins
        LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
        WHERE Proteins.proteinID IN ({protein_placeholders})
        """

        cur.execute(query, chunk)
        rows = cur.fetchall()

        for protein_id, keyword in rows:
            if keyword:
                protein_to_keywords[protein_id].add(keyword)

        #print(f"[INFO] Processed proteins {start+1} to {min(end, total)} of {total}")

    conn.close()

    logger.debug(f"Retrieved keywords for {len(protein_to_keywords)} proteins.")

    # Union of all keywords associated with these reference sequences
    all_keywords = set().union(*protein_to_keywords.values())

    return all_keywords



def integrate_csb_variants_into_merged_grouped(
    options: Any,
    merged_grouped: Dict[str, Set[str]],
    domain_to_new_keywords_dict: Dict[str, Set[str]],
    chunk_size: int = 999
) -> Dict[str, Set[str]]:
    """
    For each domain, integrates all proteins whose CSB keyword matches those in domain_to_new_keywords_dict.

    Args:
        options: Options/config object (needs database_directory)
        merged_grouped (dict): {domain: set(proteinIDs)} to be expanded
        domain_to_new_keywords_dict (dict): {domain: set(keywords) to add}
        chunk_size (int): DB chunk size.

    Returns:
        dict: Updated merged_grouped.

    Example:
        ({'ABC': {...}}, {'ABC': {'csb5'}}) → {'ABC': {..., 'p55', 'p56'}}
    """

    logger.debug("Integration of added CSB proteins to grouped dataset")


    conn = sqlite3.connect(options.database_directory)
    cur = conn.cursor()

    total_proteins_added = 0

    for domain, new_keywords in domain_to_new_keywords_dict.items():
        if not new_keywords:
           logger.debug(f"Domain {domain}: No new keywords to integrate.")
           continue

        logger.debug(f"Domain {domain}: Integrating sequences from {len(new_keywords)} new keywords.")

        new_keywords_list = list(new_keywords)
        domain_proteins_before = len(merged_grouped.get(domain, set()))
        proteins_added_this_domain = 0

        for start in range(0, len(new_keywords_list), chunk_size):
            end = start + chunk_size
            chunk = new_keywords_list[start:end]

            placeholders = ','.join(['?'] * len(chunk))
            query = f"""
            SELECT Proteins.proteinID
            FROM Proteins
            LEFT JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
            LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
            WHERE Domains.domain = ? AND Keywords.keyword IN ({placeholders})
            """

            params = [domain] + chunk
            cur.execute(query, params)
            rows = cur.fetchall()

            for (protein_id,) in rows:
                if protein_id not in merged_grouped[domain]:
                    merged_grouped[domain].add(protein_id)
                    proteins_added_this_domain += 1

        domain_proteins_after = len(merged_grouped[domain])
        logger.info(f"Domain {domain}: Added {proteins_added_this_domain} new proteins with matching synteny to reference. New total: {len(merged_grouped[domain])}")

        total_proteins_added += proteins_added_this_domain

    conn.close()

    logger.info(f"Integration completed: {total_proteins_added} proteins added across all domains.")

    return merged_grouped
