#!/usr/bin/python
import sqlite3
import os

from typing import Dict, Set, Tuple, List, Any, DefaultDict

from src.core import myUtil

logger = myUtil.logger


def parse_csb_file_to_dict(filepath: str) -> Dict[str, Set[str]]:
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

    with open(filepath, "r") as file:
        for line in file:
            # Split the line by tabs
            parts = line.strip().split("\t")

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


def filter_dictionary_by_inclusion_domains(
    dictionary: Dict, include_list: List[str]
) -> Dict:
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


def filter_dictionary_by_excluding_domains(
    dictionary: Dict, exclude_list: List[str]
) -> Dict:
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


def _prepare_keyword_domain_tasks_temp(
    cur: sqlite3.Cursor, csb_dictionary: Dict[str, List[str]]
) -> int:
    """
    TEMP table tmp_tasks(keyword TEXT, domain TEXT) befüllen.
    Analog zu db_fetch_protein.py: TEMP table + executemany, damit keine 999-Variable-Limits.
    Returns: number of inserted task pairs.
    """
    cur.execute(
        "CREATE TEMP TABLE IF NOT EXISTS tmp_tasks (keyword TEXT, domain TEXT);"
    )
    cur.execute("DELETE FROM tmp_tasks;")

    pairs = []
    for keyword, domains in csb_dictionary.items():
        if keyword == "default":
            continue
        if not domains:
            continue
        for d in domains:
            if d:
                pairs.append((keyword, d))

    if not pairs:
        return 0

    # Optional: entdoppeln (kann DB auch selber, aber spart Insertarbeit)
    pairs = list(dict.fromkeys(pairs))

    cur.executemany(
        "INSERT INTO tmp_tasks(keyword, domain) VALUES (?, ?);",
        pairs,
    )
    return len(pairs)


def fetch_proteinIDs_dict(
    database_path: str,
    csb_dictionary: Dict[str, List[str]],
    min_seqs: int,
    cores: int,  # bleibt im Signature, wird aber nicht genutzt
) -> Dict[Tuple[str, str], Set[str]]:
    """
    FAST bulk-fetch: holt alle ProteinIDs für (keyword, domain) mit einem TEMP-table Join.

    Ersetzt den alten multiprocessing+chunked-IN Ansatz komplett.
    cores ist absichtlich ungenutzt: SQLite profitiert hier i.d.R. mehr von 1 Bulk-Query
    als von vielen parallelen Readern.

    Returns:
        {(keyword, domain): set(proteinIDs)}
    """
    out: Dict[Tuple[str, str], Set[str]] = {}

    # Normal connection (no mode=ro), because TEMP TABLE needs write capability.
    # We will protect the main DB with PRAGMA query_only=TRUE after TEMP setup.
    with sqlite3.connect(database_path, timeout=120.0) as con:
        cur = con.cursor()

        # --- Performance pragmas (safe) ---
        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")  # ~256 MiB
        cur.execute("PRAGMA mmap_size=2147483648;")  # 2 GiB
        cur.execute("PRAGMA automatic_index=ON;")

        # Build + fill TEMP tasks table
        n_tasks = _prepare_keyword_domain_tasks_temp(cur, csb_dictionary)
        if n_tasks == 0:
            return out

        # Now lock down the main database: no writes to persistent tables possible.
        cur.execute("PRAGMA query_only=TRUE;")

        # --- Single bulk join query ---
        cur.execute(
            """
            SELECT
                t.keyword   AS keyword,
                t.domain    AS domain,
                p.proteinID AS proteinID
            FROM tmp_tasks t
            JOIN Keywords k
              ON k.keyword = t.keyword
            JOIN Proteins p
              ON p.clusterID = k.clusterID
            JOIN Domains d
              ON d.proteinID = p.proteinID
             AND d.domain    = t.domain
            """
        )

        for keyword, domain, proteinID in cur:
            key = (keyword, domain)
            s = out.get(key)
            if s is None:
                out[key] = {proteinID}
            else:
                s.add(proteinID)

    # Apply min_seqs filter (keep only (keyword,domain) with enough proteins)
    if min_seqs and min_seqs > 1:
        out = {k: v for k, v in out.items() if len(v) >= min_seqs}

    return out


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
        con.execute("PRAGMA journal_mode=WAL;")
        con.execute("PRAGMA journal_mode=WAL;")
        con.execute("PRAGMA synchronous=NORMAL;")
        con.execute("PRAGMA temp_store=MEMORY;")
        con.execute("PRAGMA cache_size=-25000;")  # ca. 100MB Cache

        # Query to get distinct domains from the selfblast associated with genomeID 'QUERY'
        query = """
            SELECT DISTINCT Domains.domain
            FROM Proteins
            LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
            WHERE Proteins.genomeID = ?;
        """
        # Execute the query with 'QUERY' as the genomeID
        cur.execute(query, ("QUERY",))
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


def csb_proteins_datasets(
    options: Any, grouped_keywords_dict: Dict
) -> Dict[Tuple[str, str], Set[str]]:
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
        logger.debug(
            "Loaded existing reference protein sequence dataset from cache: csb_protein_dataset.pkl"
        )
        dictionary = filter_dictionary_by_inclusion_domains(
            dictionary, options.include_list
        )
        dictionary = filter_dictionary_by_excluding_domains(
            dictionary, options.exclude_list
        )
        return dictionary

    ## If not precomputed then compute the csb to select and fetch die proteinIDs from the dataset
    # 1. Get the domain types as set by the selection statistic module for the csb finder algorithm
    csb_dictionary = parse_csb_file_to_dict(
        options.csb_output_file
    )  # dictionary with cbs_name => csb items

    # Determine the csb that shall be kept.
    csbs_to_keep = set()
    for domain, group_lists in grouped_keywords_dict.items():
        for group in group_lists:
            csbs_to_keep.update(group)

    # Now, filter the dictionary to keep only the wanted CSBs
    csb_dictionary = {
        csb: val for csb, val in csb_dictionary.items() if csb in csbs_to_keep
    }

    # 2. add the static csb from the pattern file
    pattern_dictionary = parse_csb_file_to_dict(
        options.patterns_file
    )  # Fetch the ones that are in the pattern file
    csb_dictionary = {**csb_dictionary, **pattern_dictionary}
    options.csb_dictionary = csb_dictionary  # save the patterns for later use

    # Fetch for each csb id all the domains in the csb that are query domains
    # dictionary is: dict[(keyword, domain)] => set(proteinIDs)
    myUtil.save_cache(options, "csb_selected_csb_to_domain_dataset.pkl", csb_dictionary)

    logger.info("Fetching the protein sequence identifiers from local database")
    dictionary = fetch_proteinIDs_dict(
        options.database_directory, csb_dictionary, options.min_seqs, options.cores
    )
    dictionary = remove_non_query_clusters(
        options.database_directory, dictionary
    )  # delete all that are not in accordance with query

    myUtil.save_cache(options, "csb_protein_dataset.pkl", dictionary)

    # Remove domains that are excluded by user options
    dictionary = filter_dictionary_by_inclusion_domains(
        dictionary, options.include_list
    )
    dictionary = filter_dictionary_by_excluding_domains(
        dictionary, options.exclude_list
    )

    return dictionary
