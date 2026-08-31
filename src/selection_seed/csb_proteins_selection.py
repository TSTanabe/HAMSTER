#!/usr/bin/python
import random
import sqlite3
import os
import multiprocessing

import time
import threading

from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from typing import Dict, Set, Tuple, List, Any, DefaultDict

from src.selection_seed import csb_type_statistic
from src.selection_seed import fetch_seed_proteins
from src.core import myUtil


logger = myUtil.logger


# first get the csb with their identifiers. make sets of appearence as value to key name of the csb
# do not compare the csb appearences as they jaccard itself should have managed it.
def _csb_proteins_datasets_combine(
    keyword_lists: Dict[str, str], csb_proteins_dict: Dict[tuple[str, str], set[str]]
) -> Dict[str, Set[str]]:
    """
    Combines the grouped keywords dict and csb proteins dict to create final grouped training sets.

    Args:
        keyword_lists (dict): {group: [keyword1, keyword2, ...]}
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


def prepare_csb_grouped_seed_proteins(
    config: Any,
) -> Tuple[
    Dict[str, Dict[str, float]],
    Dict[str, Set[str]],
]:
    """
    Prepares grouped protein sequences for training by analyzing gene clusters (CSBs)
    and extracting associated proteins with conserved genomic context.

    This function:
    1. Computes score limits and cluster keyword groupings.
    2. Extracts training protein sequences.
    3. Loads or computes grouped training sets from cache.
    4. Writes grouped sequences to FASTA files if needed.

    Args:
        config (argparse.Namespace): The runtime configuration.

    Returns:
        tuple:
            - grp_score_limit_dict (dict): Domain-wise score limits.
            - grouped (dict): Mapping from group protein label → set of protein IDs.
    """

    # Step 1: Try loading cached grouped training data
    grouped = myUtil.load_cache(
        config, "grp0_training_proteinIDs.pkl"
    )  # Name is defined by fetch_to_fasta routine
    grp_score_limit_dict = myUtil.load_cache(config, "grp0_score_limit_dict.pkl")
    if grouped and grp_score_limit_dict:
        return grp_score_limit_dict, grouped

    # Step 2: Compute score limits and keyword clusters
    # seed_grouped_keywords_dict => domain: [[csb_keyword1, csb_keyword2, ...]]
    grp_score_limit_dict, _, seed_grouped_keywords_dict = (
        csb_type_statistic.group_gene_cluster_statistic(config)
    )

    logger.info(
        f"Collecting sequences for training datasets with similar csb for "
        f"{len(seed_grouped_keywords_dict)} proteins: "
        f"{', '.join(sorted(seed_grouped_keywords_dict.keys()))}"
    )

    csb_proteins_dict: dict[tuple[str, str], set[str]] = (
        fetch_seed_proteins.csb_proteins_datasets(config, seed_grouped_keywords_dict)
    )

    logger.info("Processing highly similar homologs with specific genomic context")

    # Step 3: Export one fasta per protein
    grouped = _csb_proteins_datasets_combine(
        seed_grouped_keywords_dict, csb_proteins_dict
    )
    # grouped = _add_query_ids_to_proteinIDset(grouped, options.database_directory)

    logger.info(
        f"Found {len(grouped)} proteins types with syntenic gene cluster patterns: {', '.join(sorted(grouped.keys()))}"
    )

    # Step 4: Save in pkl cache
    myUtil.save_cache(config, "grp0_training_proteinIDs.pkl", grouped)
    myUtil.save_cache(config, "grp0_score_limit_dict.pkl", grp_score_limit_dict)

    return grp_score_limit_dict, grouped


################################################################################################


def process_keyword_domains(
    args: Tuple[str, str, List[str], int],
) -> Dict[Tuple[str, str], Set[str]]:
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
        cur.execute(
            """
            SELECT DISTINCT Proteins.proteinID
            FROM Proteins
            INNER JOIN Keywords ON Proteins.clusterID = Keywords.clusterID
            WHERE Keywords.keyword = ?;
        """,
            (keyword,),
        )

        protein_ids = {
            row[0] for row in cur.fetchall()
        }  # Set für schnelleres Nachschlagen

        if not protein_ids:
            return result  # Falls keine Treffer, direkt zurückgeben

        # Schritt 2: Passende Protein-IDs mit den gewünschten Domains abrufen (in Chunks)
        query_conditions = " OR ".join(["Domains.domain = ?"] * len(domains))

        for i in range(0, len(protein_ids), chunk_size):
            chunk = list(protein_ids)[i : i + chunk_size]  # Nimm max. 900 Protein-IDs

            query = f"""
                SELECT DISTINCT Proteins.proteinID, Domains.domain
                FROM Proteins
                INNER JOIN Domains ON Proteins.proteinID = Domains.proteinID
                WHERE ({query_conditions}) AND Proteins.proteinID IN ({",".join(["?"] * len(chunk))});
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


######################################################################################################


###############################################################################
#################### Protein to fasta operations ##############################
###############################################################################


def _connect_readonly(db_path):
    # helper routine for the fetch seq to fasta database connection
    uri = f"file:{db_path}?mode=ro&immutable=1"
    return sqlite3.connect(uri, uri=True, check_same_thread=False)


def _fetch_seq_to_fasta(
    database: str,
    domain: str,
    protein_ids: Set[str],
    output_directory: str,
    chunk_size: int = 990,
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

    with _connect_readonly(database) as con:
        cur = con.cursor()

        with open(fasta_file_path, "w") as fasta_file:
            protein_id_list = list(protein_ids)
            for i in range(0, len(protein_id_list), chunk_size):
                chunk = protein_id_list[i : i + chunk_size]

                if not chunk:
                    continue  # Ensure chunk is not empty
                # Convert all protein IDs to strings before passing to SQL
                query = f"""
                    SELECT proteinID, sequence FROM Proteins
                    WHERE proteinID IN ({",".join(["?"] * len(chunk))});
                """

                try:
                    cur.execute(query, tuple(str(protein_id) for protein_id in chunk))
                    rows = cur.fetchall()

                    # Write sequences to the domain-specific FASTA file
                    for protein_id, sequence in rows:
                        fasta_file.write(f">{protein_id}\n{sequence}\n")
                except sqlite3.InterfaceError as e:
                    logger.error(
                        f"SQL Error in domain {domain}: {e}"
                    )  # Debugging message

    logger.debug(f"FASTA file saved: {fasta_file_path}")


def _fetch_seqs_to_fasta_parallel(
    database: str,
    dataset_dict: Dict[str, Set[str]],
    output_directory: str,
    min_seq: int,
    max_seq: int,
    cores: int = 4,
    chunk_size: int = 990,
    hardcap: int = 5000,
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
            logger.warn(
                f"'{domain}' was skipped due to too few sequences: {num_sequences} < {min_seq}"
            )
            continue  # Skip this domain

        if num_sequences > (max_seq + hardcap):
            logger.warning(
                f"'{domain}' has too many sequences ({num_sequences}), "
                f"randomly subsampling to {max_seq}"
            )
            protein_ids = set(random.sample(sorted(protein_ids), max_seq))

        if os.path.exists(output_file):
            logger.debug(
                f"Skipped '{domain}'. FASTA file already exists: {output_file}"
            )
            continue  # Skip existing files

        # If all checks pass, add to tasks
        tasks.append((database, domain, protein_ids, output_directory, chunk_size))

    if not tasks:
        return  # Exit if no tasks remain
    logger.info(
        "Fetching amino acid sequences from local database. Might take some minutes ..."
    )
    # Use multiprocessing to run fetch_seq_to_fasta in parallel
    with multiprocessing.Pool(processes=cores) as pool:
        pool.starmap(_fetch_seq_to_fasta, tasks)


def fetch_training_data_to_fasta(
    options: Any, grouped: Dict[str, Set[str]], prefix: str
) -> None:
    """
    Writes training FASTA files per group/domain, using a prefix.

    Args:
        options (Namespace): Configuration/options.
        grouped (dict): {group: set(proteinIDs)}
        prefix (str): Prefix for each group name (output FASTA: {prefix}_group.faa)

    Returns:
        None (writes FASTA files).
    """

    extended_grouped_prefixed = {
        f"{prefix}_{key}": value for key, value in grouped.items()
    }  # Extend dictionary with a prefix

    _fetch_seqs_to_fasta_parallel(
        options.database_directory,
        extended_grouped_prefixed,
        options.fasta_output_directory,
        min_seq=options.min_seqs,
        max_seq=options.max_seqs,
        cores=options.cores,
        hardcap=options.hardcap,
    )


def fetch_protein_family_sequences(options, directory, score_limit_dict, grouped):
    # Training datasets with additional sequences
    score_limit_dict = filter_existing_faa_files(
        score_limit_dict, directory
    )  # Do not fetch again for existing files
    decorated_grouped_dict = fetch_protein_ids_parallel(
        options.database_directory, score_limit_dict, options.cores, options.max_seqs
    )  # get the proteinIDs within the score limits for each domain, new keys are domain only
    decorated_grouped_dict = merge_grouped_protein_ids(decorated_grouped_dict, grouped)
    _fetch_seqs_to_fasta_parallel(
        options.database_directory,
        decorated_grouped_dict,
        directory,
        options.min_seqs,
        options.max_seqs,
        options.cores,
    )

    return


def fetch_protein_ids_for_domain(
    database: str,
    domain: str,
    lower_limit: float,
    upper_limit: float,
    max_count: int = 50000,
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
        logger.warning(
            f"max_count {max_count} reached for {domain}, but {tie_count} more proteins with same bitscore ({last_score}) included."
        )

    return domain, protein_ids


def filter_existing_faa_files(
    domain_dict: Dict[str, Any], directory: str
) -> Dict[str, Any]:
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
    filtered_dict = {
        domain: values
        for domain, values in domain_dict.items()
        if f"{domain}.faa" not in existing_files
    }

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


def merge_protein_sets(
    *grouped_dicts: Dict[str, Set[str]],
) -> Dict[str, Set[str]]:
    merged: Dict[str, Set[str]] = {}

    for grouped in grouped_dicts:
        for domain, protein_ids in grouped.items():
            merged.setdefault(domain, set()).update(protein_ids)

    return merged


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
        con.execute("PRAGMA journal_mode=WAL;")
        con.execute("PRAGMA synchronous=NORMAL;")
        con.execute("PRAGMA temp_store=MEMORY;")
        con.execute("PRAGMA cache_size=-25000;")  # ca. 100MB Cache

        # Retrieve all distinct domains from the database
        query = "SELECT DISTINCT domain FROM Domains"
        cur.execute(query)
        domains = [row[0] for row in cur.fetchall()]

        for domain in domains:
            # Create the output filename for the domain
            name = domain.split("_")[-1]
            filename = f"superfamily_{name}.faa"
            output_fasta_path = os.path.join(directory, filename)

            # Skip if the file already exists
            if os.path.isfile(output_fasta_path):
                logger.debug(f"File {output_fasta_path} already exists - skipping.")
                continue

            # Fetch all proteins associated with this domain
            query = """
                SELECT DISTINCT P.proteinID, P.sequence
                FROM Proteins P
                INNER JOIN Domains D ON P.proteinID = D.proteinID
                WHERE D.domain = ?
            """
            cur.execute(query, (domain,))
            proteins = cur.fetchall()

            # Write all sequences to the FASTA file
            with open(output_fasta_path, "w") as fasta_file:
                for proteinID, sequence in proteins:
                    fasta_file.write(f">{proteinID}\n{sequence}\n")

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
        con.execute("PRAGMA journal_mode=WAL;")
        con.execute("PRAGMA synchronous=NORMAL;")
        con.execute("PRAGMA temp_store=MEMORY;")
        con.execute("PRAGMA cache_size=-25000;")  # ca. 100MB Cache

        # get all proteins with id and key output shall be "proteinId csb_proteintype"
        query = """
                SELECT DISTINCT Proteins.proteinID,Proteins.Sequence
                FROM Proteins WHERE NOT Proteins.genomeID = ?
            """
        cur.execute(query, ("QUERY",))
        rows = cur.fetchall()
        with open(filepath, "w") as fasta_file:
            for row in rows:
                fasta_file.write(f">{row[0]}\n")
                fasta_file.write(f"{row[1]}\n")

    return filepath


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


def _fetch_neighbours_chunk(
    args: Tuple[str, List[str]],
) -> Tuple[Dict[str, List[List[str]]], Dict[Tuple, int]]:
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
        protein_ids,
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
        cids,
    )
    nbr_rows = list(cur.fetchall())
    nbr_rows.sort(key=lambda r: (r["clusterID"], r["start"]))

    # 3) neighbour protein → domains
    nids = [r["proteinID"] for r in nbr_rows]
    ph3 = ",".join("?" for _ in nids)
    cur.execute(
        f"SELECT proteinID, domain FROM Domains WHERE proteinID IN ({ph3})", nids
    )
    dom_rows = cur.fetchall()

    dt = time.perf_counter() - t0
    logger.debug(
        f"[TIMING] chunk {len(protein_ids)} → {len(nids)} neighbours in {dt:.3f}s"
    )

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
        ]
        or [["singleton"]]
        for hit, cid in hit2cid.items()
    }

    # Count unique compositions per cluster
    for cid, ids in cluster_to_nids.items():
        comp = tuple(
            sorted(tuple(sorted(domains_by_nid.get(nid, ["singleton"]))) for nid in ids)
        )
        comp_counts[comp] += 1

    return neighbours, comp_counts


def fetch_neighbouring_genes_with_domains(
    db_path: str, protein_ids: Set[str], chunk_size: int = 999, threads: int = 8
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
        protein_list[i : i + chunk_size] for i in range(0, len(protein_ids), chunk_size)
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

    sorted_compositions = sorted(total_counts.items(), key=lambda x: x[1], reverse=True)
    return full_neighbours, sorted_compositions
