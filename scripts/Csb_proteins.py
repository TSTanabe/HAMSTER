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
from . import Csb_statistic
from . import myUtil

#first get the csb with their identifiers. make sets of appearence as value to key name of the csb
#do not compare the csb appearences as they jaccard itself should have managed it.

def prepare_csb_grouped_training_proteins(options):
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
            - grouped (dict): Mapping from group label → set of protein IDs.
    """

    # Step 1: Compute score limits and keyword clusters
    grp_score_limit_dict, _, grouped_keywords_dict, clustered_excluded_keywords_dict = Csb_statistic.group_gene_cluster_statistic(options)

    print("[INFO] Collecting sequences for training datasets with similar csb")
    csb_proteins_dict = csb_proteins_datasets(options, clustered_excluded_keywords_dict)

    print("[INFO] Processing highly similar homologs with specific genomic context")

    # Step 2: Try loading cached grouped training data
    grouped = myUtil.load_cache(options, 'grp_training_proteinIDs.pkl')

    # Step 3: If cache is missing, recompute and export FASTAs
    if not grouped:
        grouped = csb_proteins_datasets_combine(grouped_keywords_dict, csb_proteins_dict)
        grouped = add_query_ids_to_proteinIDset(grouped, options.database_directory)
        
        
        extended_grouped_prefixed = {f"grp0_{key}": value for key, value in grouped.items()}
        fetch_seqs_to_fasta_parallel(
            options.database_directory,
            extended_grouped_prefixed,
            options.fasta_output_directory,
            min_seq=options.min_seqs,
            max_seq=options.max_seqs,
            cores=options.cores
        )
        myUtil.save_cache(options, 'grp_training_proteinIDs.pkl', grouped)

    return grp_score_limit_dict, grouped
    
    
    
def csb_proteins_datasets(options, sglr_dict):
    
    # Get the domain types as set per csb and predefined pattern
    print("[INFO] Initilize collinear syntenic block sorting")
    csb_dictionary = parse_csb_file_to_dict(options.csb_output_file) # dictionary with cbs_name => csb items
    pattern_dictionary = parse_csb_file_to_dict(options.patterns_file)  # Fetch the ones that are in the pattern file
    csb_dictionary = {**csb_dictionary, **pattern_dictionary}
    options.csb_dictionary = csb_dictionary # save the patterns for later use

    # Load existing data
    dictionary = myUtil.load_cache(options, "csb_protein_dataset.pkl")
    
    if dictionary:
        print("[LOAD] Loaded existing reference protein sequence dataset")
        dictionary = filter_dictionary_by_inclusion_domains(dictionary, options.include_list)
        dictionary = filter_dictionary_by_excluding_domains(dictionary, options.exclude_list)
        return dictionary    

    #Fetch for each csb id all the domains in the csb that are query domains
    #dictionary is: dict[(keyword, domain)] => set(proteinIDs)
    
    csbs_to_remove = {csb for csb_list in sglr_dict.values() for sublist in csb_list for csb in sublist}
    csb_dictionary = {csb: domains for csb, domains in csb_dictionary.items() if csb not in csbs_to_remove}
    
    myUtil.save_cache(options, "csb_selected_csb_to_domain_dataset.pkl", csb_dictionary)

    print("[INFO] Collecting protein sequence identifiers from local database")
    dictionary = fetch_proteinIDs_dict_multiprocessing(options.database_directory,csb_dictionary,options.min_seqs,options.cores)

    dictionary = remove_non_query_clusters(options.database_directory, dictionary) #delete all that are not in accordance with query
    
    myUtil.save_cache(options, "csb_protein_dataset.pkl", dictionary)
    
    # Remove domains that are excluded by user options
    dictionary = filter_dictionary_by_inclusion_domains(dictionary, options.include_list)
    
    dictionary = filter_dictionary_by_excluding_domains(dictionary, options.exclude_list)

    return dictionary

################################################################################################

def csb_proteins_datasets_combine(keyword_lists, csb_proteins_dict):
    """
    Kombiniert Protein-IDs basierend auf den gruppierten und ausgeschlossenen Keywords.
    
    Args:
        keyword_lists (dict): Dictionary mit gruppierten Keywords als Listen von Listen.
        csb_proteins_dict (dict): Das Rückgabedictionary aus csb_proteins_datasets, 
                                  das (keyword, domain) als Schlüssel und Protein-IDs als Werte enthält.
    
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
                combined_protein_sets[domain] = protein_set
    
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

def decorate_training_data(options, score_limit_dict, grouped):

    # Training datasets with additional sequences
    score_limit_dict = filter_existing_faa_files(score_limit_dict, options.phylogeny_directory) # Do not fetch again for existing files
    decorated_grouped_dict = fetch_protein_ids_parallel(options.database_directory, score_limit_dict, options.cores, options.max_seqs) # get the proteinIDs within the score limits for each domain, new keys are domain only
    decorated_grouped_dict = merge_grouped_protein_ids(decorated_grouped_dict, grouped)
    fetch_seqs_to_fasta_parallel(options.database_directory, decorated_grouped_dict, options.phylogeny_directory, options.min_seqs, options.max_seqs, options.cores)
    
    return
################################################################################################

    
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
            print(f"  {i}/{total_tasks}", end="\r")
            results.append(result)

    # Combine results from all workers
    combined_dict = {}
    for result in results:
        for key, value in result.items():
            if key not in combined_dict:
                combined_dict[key] = set()
            combined_dict[key].update(value)

    print("\n[INFO] Completed training data selection")
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
            print(f"[WARN] Sequences for '{domain}' skipped (too few sequences: {num_sequences} < {min_seq})")
            continue  # Skip this domain

        if num_sequences > (max_seq + 5000):
            print(f"[WARN] Sequences for '{domain}' skipped (too many sequences: {num_sequences} > {max_seq+5000})")
            continue  # Skip this domain

        if os.path.exists(output_file):
            print(f"[SKIP] '{domain}', FASTA file already exists: {output_file}")
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
                    print(f"[ERROR] SQL Error in domain {domain}: {e}")  # Debugging message

    print(f"[INFO] FASTA file saved: {fasta_file_path}")




def fetch_protein_ids_for_domain(database, domain, lower_limit, upper_limit, max_count=50000):
    """
    Efficiently fetch top protein IDs for a domain, including all with same score at cutoff,
    with only one sorted query. Prints a warning if max_count is exceeded due to score ties.

    Args:
        database (str): Path to SQLite database.
        domain (str): Domain name.
        lower_limit (float): Lower score threshold.
        upper_limit (float): Upper score threshold.
        max_count (int): Minimum number of top scoring entries to fetch (ties allowed above).

    Returns:
        tuple: (domain, set(proteinIDs))
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
        print(f"[WARN] max_count of {max_count} reached for {domain}, but {tie_count} additional protein(s) "
              f"with the same score ({last_score}) were included.")

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
                print(f"[SKIP] File {output_fasta_path} already exists.")
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
    
    
    

##########################################################################################
###################### Get genomic context per protein ###################################
##########################################################################################

_thread_local = threading.local()


def _get_thread_cursor(db_path: str):
    """
    Returns a thread-local Cursor. On the first call in each thread:
     - Opens the DB in truly read-only+immutable mode
     - Disables any attempt to create journals or temporary files

    Subsequent calls from the same thread return the same Cursor.
    """
    if getattr(_thread_local, "cur", None) is None:
        # 1) Open connection with mode=ro & immutable=1
        #    &nolockfile forces SQLite to never try to create a rollback journal
        uri = f"file:{db_path}?mode=ro&immutable=1&nolockfile=1"
        con = sqlite3.connect(uri, uri=True, check_same_thread=False)

        # 2) We still set row_factory so fetch returns Row objects
        con.row_factory = sqlite3.Row
        cur = con.cursor()

        # 3) Explicitly disable any journaling or temp file usage
        #    - journal_mode = OFF means no rollback or WAL is ever used
        #    - temp_store = MEMORY means any temporary table/index is purely in RAM
        #    - query_only = TRUE forbids any accidental writes
        cur.executescript("""
            PRAGMA query_only = TRUE;
            PRAGMA journal_mode = OFF;
            PRAGMA temp_store = MEMORY;
            PRAGMA defer_foreign_keys = TRUE;
        """)

        # 4) Save to thread-local
        _thread_local.con = con
        _thread_local.cur = cur

    return _thread_local.cur
def _fetch_neighbours_chunk(args):
    """
    Worker für einen Chunk, mit gecachedem Cursor und Timing.
    """
    db_path, protein_ids = args
    if not protein_ids:
        return {}, {}

    # Timer starten
    t0 = time.perf_counter()

    # Cursor holen (erstellt Connection + PRAGMAs einmal pro Thread)
    cur = _get_thread_cursor(db_path)

    # --- 1) hit → clusterID ---
    ph1 = ",".join("?" for _ in protein_ids)
    cur.execute(
        f"SELECT proteinID, clusterID FROM Proteins WHERE proteinID IN ({ph1})",
        protein_ids
    )
    hit2cid = {r["proteinID"]: r["clusterID"] for r in cur.fetchall()}

    if not hit2cid:
        dt = time.perf_counter() - t0
        print(f"[TIMING] chunk {len(protein_ids)} (no clusters) in {dt:.3f}s")
        return {pid:[["singleton"]] for pid in protein_ids}, {}

    # --- 2) clusterID → neighbours (+ start) ---
    cids = list(set(hit2cid.values()))
    ph2 = ",".join("?" for _ in cids)
    cur.execute(
        f"SELECT proteinID, clusterID, start "
        f"FROM Proteins WHERE clusterID IN ({ph2})",
        cids
    )
    nbr_rows = list(cur.fetchall())
    nbr_rows.sort(key=lambda r: (r["clusterID"], r["start"]))

    # --- 3) neighbour → domains ---
    nids = [r["proteinID"] for r in nbr_rows]
    ph3 = ",".join("?" for _ in nids)
    cur.execute(
        f"SELECT proteinID, domain FROM Domains WHERE proteinID IN ({ph3})",
        nids
    )
    dom_rows = cur.fetchall()

    # Timer stoppen
    dt = time.perf_counter() - t0
    print(f"[TIMING] chunk {len(protein_ids)} → "
          f"{len(nbr_rows)} neighbours in {dt:.3f}s")

    # --- Rest wie gehabt (Aggregation) ---
    domains_by_nid = defaultdict(list)
    for r in dom_rows:
        domains_by_nid[r["proteinID"]].append(r["domain"] or "no_neighbours")

    cluster_to_nids = defaultdict(list)
    for r in nbr_rows:
        cluster_to_nids[r["clusterID"]].append(r["proteinID"])

    neighbours = {
        hit: [
            domains_by_nid.get(nid, ["singleton"])
            for nid in cluster_to_nids.get(cid, [])
        ] or [["singleton"]]
        for hit, cid in hit2cid.items()
    }

    comp_counts = defaultdict(int)
    for cid, ids in cluster_to_nids.items():
        comp = tuple(
            sorted(tuple(sorted(domains_by_nid.get(nid, ["singleton"])))
                   for nid in ids)
        )
        comp_counts[comp] += 1

    return neighbours, comp_counts


def fetch_neighbouring_genes_with_domains(db_path,
                                                  protein_ids,
                                                  chunk_size: int = 999,
                                                  threads: int = 8):
    """
    Threaded wrapper: chunk_size große Listen parallel in Threads verarbeiten,
    Ergebnisse zusammenführen und sortierte Kompositionen zurückgeben.
    """
    # Ensure it's a list so we can index/slice
    protein_list = list(protein_ids)

    if not protein_ids:
        return {}, []

    # 1) Chunking
    chunks = [
        protein_list[i : i + chunk_size]
        for i in range(0, len(protein_ids), chunk_size)
    ]
    args = [(db_path, chunk) for chunk in chunks]

    # 2) Parallel mit ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=threads) as executor:
        partials = list(executor.map(_fetch_neighbours_chunk, args))
    
    # 3) Merge results
    full_neighbours = {}
    total_counts = defaultdict(int)
    for neigh_dict, comp_dict in partials:
        full_neighbours.update(neigh_dict)
        for comp, cnt in comp_dict.items():
            total_counts[comp] += cnt
    
    # 4) Sort compositions
    sorted_compositions = sorted(
        total_counts.items(), key=lambda x: x[1], reverse=True
    )

    return full_neighbours, sorted_compositions
    
    
    
    
##########################################################################################################
######################### Extend grouped reference proteins with similar csb #############################
##########################################################################################################

def extend_merged_grouped_by_csb_similarity(options, grouped):
    
    # Main routine for addition of similar proteins from similar csb to the grouped datasets
    
    # Find keywords that are in jaccard distance to the csb of reference sequences in grouped
    protein_to_new_keywords_dict = select_similar_csb_patterns_per_protein(options, grouped, options.jaccard)
    
    # Integrate proteins with the new keywords into merged_grouped dataset
    extended_grouped = integrate_csb_variants_into_merged_grouped(options, grouped, protein_to_new_keywords_dict, options.sqlite_chunks)

    return extended_grouped







def select_similar_csb_patterns_per_protein(options, merged_grouped, jaccard_threshold=0.7):
    """
    Extend merged_grouped by finding CSB patterns with high similarity to the patterns 
    associated with the current proteins in each domain.

    Args:
        options: Options object containing paths and settings.
        merged_grouped (dict): Current merged_grouped dictionary { domain: set(proteinIDs) }.
        jaccard_threshold (float): Jaccard similarity threshold.

    Returns:
        dict: jaccard_included_patterns { domain: set(similar CSB keywords) }.
    """

    # Stores the similar CSB keywords for each domain
    jaccard_included_patterns = {}  # { domain : set(csb_keywords) }

    # Load all CSB patterns from csb_output_file
    csb_dictionary = parse_csb_file_to_dict(options.csb_output_file)

    # Process each domain
    for domain, protein_ids in merged_grouped.items():
        print(f"[INFO] Processing domain {domain}")

        # Fetch all keywords associated with proteins of this domain
        all_keywords = fetch_keywords_for_proteins(options.database_directory, protein_ids)

        if not all_keywords:
            print(f"[WARN] No keywords found for domain {domain}")
            continue

        # Build union of patterns from these keywords
        domain_pattern_union = set().union(
            *[csb_dictionary[k] for k in all_keywords if k in csb_dictionary]
        )

        print(f"[INFO] Domain {domain} - Pattern size: {len(domain_pattern_union)}")

        # Compare against all CSB patterns and collect similar ones
        similar_csb_keywords = set()

        for csb_key, csb_pattern in csb_dictionary.items():
            intersection = len(domain_pattern_union & csb_pattern)
            union = len(domain_pattern_union | csb_pattern)
            similarity = intersection / union if union > 0 else 0.0

            if similarity >= jaccard_threshold:
                similar_csb_keywords.add(csb_key)

        print(f"[INFO] Domain {domain} - Found {len(similar_csb_keywords)} similar CSB patterns (Jaccard >= {jaccard_threshold})")

        # Store the result for this domain
        jaccard_included_patterns[domain] = similar_csb_keywords

    # Return dictionary: domain → set of similar CSB keywords
    return jaccard_included_patterns



def fetch_keywords_for_proteins(database_path, protein_ids, chunk_size=999):
    """
    Fetch keywords for a list of proteinIDs from the database (chunked).

    Args:
        database_path (str): Path to the SQLite database.
        protein_ids (iterable): List or set of proteinIDs.
        chunk_size (int): Maximum size of chunks (default 999).

    Returns:
        set: Union of all keywords associated with these proteinIDs.
    """

    protein_to_keywords = defaultdict(set)

    if not protein_ids:
        print("[WARN] No proteinIDs provided to fetch_keywords_for_proteins.")
        return set()

    conn = sqlite3.connect(database_path)
    cur = conn.cursor()

    protein_ids = list(protein_ids)
    total = len(protein_ids)
    print(f"[INFO] Fetching keywords for {total} proteins (chunk size {chunk_size})")

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

        print(f"[INFO] Processed proteins {start+1} to {min(end, total)} of {total}")

    conn.close()

    print(f"[INFO] Retrieved keywords for {len(protein_to_keywords)} proteins.")

    # Union of all keywords associated with these reference sequences
    all_keywords = set().union(*protein_to_keywords.values())

    return all_keywords



def integrate_csb_variants_into_merged_grouped(options, merged_grouped, domain_to_new_keywords_dict, chunk_size=999):
    """
    Integrates proteins matching new CSB keywords into merged_grouped.

    Args:
        options: Options object with paths and settings.
        merged_grouped (dict): Current merged_grouped { domain : set(proteinIDs) }.
        domain_to_new_keywords_dict (dict): { domain : set(new_keywords) }.

    Returns:
        updated merged_grouped (dict).
    """

    print(f"[INFO] Starting integration of new keywords into merged_grouped...")

    conn = sqlite3.connect(options.database_directory)
    cur = conn.cursor()

    total_proteins_added = 0

    for domain, new_keywords in domain_to_new_keywords_dict.items():
        if not new_keywords:
            print(f"[INFO] Domain {domain}: No new keywords to integrate.")
            continue

        print(f"[INFO] Domain {domain}: Integrating sequences from {len(new_keywords)} new keywords")

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
        print(f"[INFO] Added {proteins_added_this_domain} new {domain} based on similar synteny. New total: {domain_proteins_after}")

        total_proteins_added += proteins_added_this_domain

    conn.close()

    print(f"[INFO] Integration complete: {total_proteins_added} proteins added across all domains.")

    return merged_grouped
