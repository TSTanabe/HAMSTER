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
) -> Dict[str, Set[str]]:
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
    grouped_seed_sequences = myUtil.load_cache(
        config, "grp0_training_proteinIDs.pkl"
    )  # Name is defined by fetch_to_fasta routine
    if grouped_seed_sequences:
        return grouped_seed_sequences

    # Step 2: Compute score limits and keyword clusters
    # seed_grouped_keywords_dict => domain: [[csb_keyword1, csb_keyword2, ...]]
    seed_grouped_keywords_dict = csb_type_statistic.group_gene_cluster_statistic(config)

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
    grouped_seed_sequences = _csb_proteins_datasets_combine(
        seed_grouped_keywords_dict, csb_proteins_dict
    )
    # grouped = _add_query_ids_to_proteinIDset(grouped, options.database_directory)

    logger.info(
        f"Found {len(grouped_seed_sequences)} proteins types with syntenic gene cluster patterns: {', '.join(sorted(grouped_seed_sequences.keys()))}"
    )

    # Step 4: Save in pkl cache
    myUtil.save_cache(config, "grp0_training_proteinIDs.pkl", grouped_seed_sequences)

    return grouped_seed_sequences


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
) -> None:
    output_file = os.path.join(output_directory, f"{domain}.faa")

    if not protein_ids or os.path.exists(output_file):
        return

    uri = f"file:{database}?mode=ro"

    with sqlite3.connect(uri, uri=True) as con:
        con.execute("PRAGMA temp_store=FILE;")

        con.execute("""
            CREATE TEMP TABLE tmp_ids (
                proteinID TEXT PRIMARY KEY
            ) WITHOUT ROWID;
        """)

        con.executemany(
            "INSERT OR IGNORE INTO tmp_ids VALUES (?);",
            ((str(pid),) for pid in protein_ids),
        )

        tmp_file = output_file + ".tmp"

        with open(tmp_file, "w") as fasta:
            cur = con.execute("""
                SELECT p.proteinID, p.sequence
                FROM Proteins AS p
                JOIN tmp_ids AS t
                  ON t.proteinID = p.proteinID
                ORDER BY p.proteinID;
            """)

            for protein_id, sequence in cur:
                fasta.write(f">{protein_id}\n{sequence}\n")

        os.replace(tmp_file, output_file)


def _fetch_seq_to_fasta_deprecated(
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


def _fetch_seqs_to_fasta(
    database: str,
    dataset_dict: Dict[str, Set[str]],
    output_directory: str,
    min_seq: int,
    max_seq: int,
    hardcap: int | None = None,
) -> None:
    os.makedirs(output_directory, exist_ok=True)

    for domain, protein_ids in dataset_dict.items():
        num_sequences = len(protein_ids)
        output_file = os.path.join(output_directory, f"{domain}.faa")

        if num_sequences < min_seq:
            logger.warning(f"{domain:<15} | skipped: {num_sequences} < {min_seq}")
            continue

        if os.path.exists(output_file):
            continue

        if hardcap is not None and num_sequences > max_seq + hardcap:
            protein_ids = set(sorted(protein_ids)[:max_seq])
            logger.warning(f"{domain:<15} | reduced {num_sequences} -> {max_seq}")

        _fetch_seq_to_fasta(
            database,
            domain,
            protein_ids,
            output_directory,
        )


def fetch_training_data_to_fasta(
    options: Any,
    grouped: Dict[str, Set[str]],
    prefix: str,
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
    }

    _fetch_seqs_to_fasta(
        options.database_directory,
        extended_grouped_prefixed,
        options.fasta_output_directory,
        min_seq=options.min_seqs,
        max_seq=options.max_seqs,
        hardcap=options.hardcap,
    )


###############################################################################
#################### Get proteinIDs from database #############################
###############################################################################


def _filter_existing_faa_files(
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


def _fetch_protein_ids_for_domain(
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
            f"Maximum sequence number ({max_count}) reached for {domain}. {tie_count} sequences with the same bitscore ({last_score}) included."
        )

    return domain, protein_ids


def _fetch_protein_ids_parallel(database, score_limit_dict, cores, max_seqs=50000):
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
        results = pool.starmap(_fetch_protein_ids_for_domain, tasks)

    # Combine results into a dictionary
    domain_protein_ids = {domain: protein_ids for domain, protein_ids in results}

    return domain_protein_ids


def _stream_protein_family_to_fasta(
    con: sqlite3.Connection,
    domain: str,
    lower_limit: float,
    upper_limit: float,
    seed_ids: set[str],
    output_file: str,
    min_seqs: int,
    max_seqs: int,
    batch_size: int = 2000,
) -> None:
    con.execute("DELETE FROM tmp_selected;")

    # Always retain seed/reference proteins.
    con.executemany(
        "INSERT OR IGNORE INTO tmp_selected(proteinID) VALUES (?);",
        ((str(pid),) for pid in seed_ids),
    )

    seed_count = con.execute("SELECT COUNT(*) FROM tmp_selected;").fetchone()[0]

    remaining = max(0, max_seqs - seed_count)

    if seed_count > max_seqs:
        logger.warning(
            f"{domain:<12} | {seed_count} seed sequences exceed "
            f"max_seqs={max_seqs}; retaining all seeds"
        )

    # ---------------------------------------------------------
    # Stream candidates ordered by score.
    # No fetchall().
    # ---------------------------------------------------------

    if remaining > 0:
        read_cur = con.cursor()
        insert_cur = con.cursor()

        read_cur.execute(
            """
            SELECT proteinID, score
            FROM Domains
            WHERE domain = ?
              AND score BETWEEN ? AND ?
            ORDER BY score DESC, proteinID ASC;
            """,
            (domain, lower_limit, upper_limit),
        )

        added = 0

        while added < remaining:
            n_to_fetch = min(batch_size, remaining - added)
            rows = read_cur.fetchmany(n_to_fetch)

            if not rows:
                break

            before = con.total_changes

            insert_cur.executemany(
                """
                INSERT OR IGNORE INTO tmp_selected(proteinID)
                VALUES (?);
                """,
                ((protein_id,) for protein_id, _ in rows),
            )

            added += con.total_changes - before

    selected_count = con.execute("SELECT COUNT(*) FROM tmp_selected;").fetchone()[0]

    if selected_count < min_seqs:
        logger.warning(
            f"{domain:<12} | skipped: {selected_count} < {min_seqs} sequences"
        )
        return

    # ---------------------------------------------------------
    # Stream selected sequences directly to FASTA.
    # ---------------------------------------------------------

    tmp_file = output_file + ".tmp"
    written = 0

    try:
        with open(tmp_file, "w") as fasta:
            cur = con.execute(
                """
                SELECT p.proteinID, p.sequence
                FROM Proteins AS p
                JOIN tmp_selected AS s
                  ON s.proteinID = p.proteinID
                ORDER BY p.proteinID;
                """
            )

            for protein_id, sequence in cur:
                fasta.write(f">{protein_id}\n{sequence}\n")
                written += 1

        os.replace(tmp_file, output_file)

    except Exception:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        raise

    logger.info(
        f"{domain:<12} | seeds: {seed_count:>6} | "
        f"selected: {selected_count:>7} | written: {written:>7}"
    )


def fetch_protein_family_sequences(
    config: Any,
    directory: str,
    score_limit_dict: dict[str, dict[str, float]],
    domain_to_proteinID: dict[str, set[str]],
) -> None:
    """Stream protein-family sequences directly from SQLite to FASTA."""

    os.makedirs(directory, exist_ok=True)

    uri = f"file:{config.database_directory}?mode=ro"
    with sqlite3.connect(uri, uri=True, timeout=120.0) as con:
        con.execute("PRAGMA temp_store=FILE;")
        con.execute("PRAGMA cache_size=-32768;")

        con.execute("""
            CREATE TEMP TABLE tmp_selected (
                proteinID TEXT PRIMARY KEY
            ) WITHOUT ROWID;
        """)

        for domain, limits in score_limit_dict.items():
            output_file = os.path.join(directory, f"{domain}.faa")

            if os.path.exists(output_file):
                logger.debug(f"{domain:<12} | FASTA exists - skipping")
                continue

            _stream_protein_family_to_fasta(
                con=con,
                domain=domain,
                lower_limit=limits["lower_limit"],
                upper_limit=limits["upper_limit"],
                seed_ids=domain_to_proteinID.get(domain, set()),
                output_file=output_file,
                min_seqs=config.min_seqs,
                max_seqs=config.max_seqs,
            )


def fetch_protein_family_sequences_deprecated(
    config: Any,
    directory: str,
    score_limit_dict: dict[str, dict[str, float]],
    domain_to_proteinID: dict[str, set[str]],
) -> None:
    """Fetch protein-family sequences within score limits and write them to FASTA."""

    # Do not fetch again for existing files.
    score_limit_dict = _filter_existing_faa_files(
        domain_dict=score_limit_dict,
        directory=directory,
    )

    # Get protein IDs within the score limits for each domain.

    decorated_grouped_dict: dict[str, set[str]] = _fetch_protein_ids_parallel(
        database=config.database_directory,
        score_limit_dict=score_limit_dict,
        cores=4,
        max_seqs=config.max_seqs,
    )
    logger.info("Merging fetched sequence identifiers and seed sequence identifiers")
    decorated_grouped_dict = merge_protein_sets(
        decorated_grouped_dict,
        domain_to_proteinID,
    )

    _fetch_seqs_to_fasta(
        database=config.database_directory,
        dataset_dict=decorated_grouped_dict,
        output_directory=directory,
        min_seq=config.min_seqs,
        max_seq=config.max_seqs,
        hardcap=None,
        cores=4,
    )


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
    db_path: str,
    protein_ids: Set[str],
) -> Tuple[Dict[str, List[List[str]]], List[Tuple[Tuple, int]]]:
    if not protein_ids:
        return {}, []

    uri = f"file:{db_path}?mode=ro"
    with sqlite3.connect(uri, uri=True) as con:
        con.row_factory = sqlite3.Row
        con.execute("PRAGMA temp_store=FILE;")
        con.execute("PRAGMA cache_size=-32768;")

        # Input proteins
        con.execute("""
            CREATE TEMP TABLE tmp_hits (
                proteinID TEXT PRIMARY KEY
            ) WITHOUT ROWID;
        """)

        con.executemany(
            "INSERT OR IGNORE INTO tmp_hits VALUES (?);",
            ((str(pid),) for pid in protein_ids),
        )

        # Unique clusters represented by the input proteins
        con.execute("""
            CREATE TEMP TABLE tmp_clusters (
                clusterID TEXT PRIMARY KEY
            ) WITHOUT ROWID;
        """)

        con.execute("""
            INSERT OR IGNORE INTO tmp_clusters(clusterID)
            SELECT p.clusterID
            FROM Proteins p
            JOIN tmp_hits h ON h.proteinID = p.proteinID
            WHERE p.clusterID IS NOT NULL;
        """)

        # Input protein -> cluster
        hit2cid = {
            row["proteinID"]: row["clusterID"]
            for row in con.execute("""
                SELECT h.proteinID, p.clusterID
                FROM tmp_hits h
                LEFT JOIN Proteins p
                    ON p.proteinID = h.proteinID;
            """)
        }

        # All proteins from these clusters, fetched only once
        cluster_to_nids = defaultdict(list)

        for row in con.execute("""
            SELECT p.clusterID, p.proteinID
            FROM Proteins p
            JOIN tmp_clusters c
                ON c.clusterID = p.clusterID
            ORDER BY p.clusterID, p.start, p.proteinID;
        """):
            cluster_to_nids[row["clusterID"]].append(row["proteinID"])

        # Make lookup table for all neighbour proteins
        con.execute("""
            CREATE TEMP TABLE tmp_neighbors (
                proteinID TEXT PRIMARY KEY
            ) WITHOUT ROWID;
        """)

        con.execute("""
            INSERT OR IGNORE INTO tmp_neighbors(proteinID)
            SELECT p.proteinID
            FROM Proteins p
            JOIN tmp_clusters c
                ON c.clusterID = p.clusterID;
        """)

        # Fetch all domains once
        domains_by_nid = defaultdict(list)

        for row in con.execute("""
            SELECT d.proteinID, d.domain
            FROM Domains d
            JOIN tmp_neighbors n
                ON n.proteinID = d.proteinID
            ORDER BY d.proteinID, d.domain;
        """):
            domains_by_nid[row["proteinID"]].append(row["domain"] or "no_neighbours")

    # Build result for individual hits
    neighbours = {
        hit: [
            domains_by_nid.get(nid, ["singleton"])
            for nid in cluster_to_nids.get(cid, [])
        ]
        or [["singleton"]]
        for hit, cid in hit2cid.items()
    }

    # Count each cluster only once
    comp_counts = defaultdict(int)

    for cid, ids in cluster_to_nids.items():
        comp = tuple(
            sorted(tuple(sorted(domains_by_nid.get(nid, ["singleton"]))) for nid in ids)
        )
        comp_counts[comp] += 1

    sorted_compositions = sorted(
        comp_counts.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    return neighbours, sorted_compositions
