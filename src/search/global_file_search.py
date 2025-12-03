#!/usr/bin/python
import gzip
import os
import csv
import shutil
import traceback
from multiprocessing import Manager, Pool
from typing import List, Set, Dict, Any, Tuple

from src.core.logging import get_logger
from src.parse_reports import csb_finder, parse_reports
from src.db import database
from src.search import query_selfblast_search
from src.search import diamond_search

logger = get_logger(__name__)


def unpackgz(path: str) -> str:
    """
    Decompresses a .gz file if not already extracted.

    Args:
        path (str): Path to .gz file.

    Returns:
        str: Path to decompressed file.
    """
    if not path.endswith(".gz"):
        return path
    file = path[:-3]
    if os.path.exists(file):
        return file
    with gzip.open(path, "rb") as f_in:
        with open(file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return file


def filter_blast_table(
    output_file: str,
    blast_file: str,
    evalue_cutoff: float,
    score_cutoff: float,
    coverage_cutoff: float,
    identity_cutoff: float,
    bsr_cutoff: float,
    sequence_lengths: Dict[str, float],
    selfblast_scores: Dict[str, float],
    buffer_size: int = 10000,
) -> str:
    """
    Filters a BLASTP table by e-value, score, coverage, identity, and BSR.

    Args:
        output_file: Path for output
        blast_file: Path to input
        evalue_cutoff: Max e-value
        score_cutoff: Min bitscore
        coverage_cutoff: Min coverage (0-1)
        identity_cutoff: Min identity (0-1)
        bsr_cutoff: Min Blast Score Ratio
        sequence_lengths: dict of qseqid -> seq len
        selfblast_scores: dict of qseqid -> self-hit bitscore
        buffer_size: Write buffer

    Returns:
        str: Output file path

    Output Example:
        "results/filtered_glob.faa.diamond.tab"
    """

    # Determine the output file path
    if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
        logger.info(
            f"Filtered hit results file already exists and is non-empty: {output_file}"
        )
        return output_file

    # Open input and output files using CSV reader/writer
    with open(blast_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        buffer = []  # Buffer to store valid rows

        # Process each row efficiently
        for row in reader:
            try:
                sseqid, qseqid, evalue, bitscore, sstart, send, pident = row[:7]

                # Convert necessary values only when needed
                evalue = float(evalue)
                bitscore = float(bitscore)
                pident = float(pident)
                sstart, send = int(sstart), int(send)

                # Fetch precomputed query length and self-blast score
                query_length = sequence_lengths.get(qseqid)
                selfblast_score = selfblast_scores.get(qseqid)

                # Pre-filter based on missing values
                if not query_length or not selfblast_score:
                    buffer.append(row)
                else:
                    # Compute alignment length and coverage
                    alignment_length = abs(send - sstart) + 1
                    coverage = alignment_length / query_length

                    # Compute Blast Score Ratio (BSR)
                    bsr = bitscore / selfblast_score

                    # Apply all filtering criteria
                    if (
                        evalue <= evalue_cutoff
                        and bitscore >= score_cutoff
                        and coverage >= coverage_cutoff
                        and pident >= identity_cutoff
                        and bsr >= bsr_cutoff
                    ):
                        row.append(f"{bsr:.3f}")
                        buffer.append(row)

                # **Write buffer to disk when it reaches buffer_size**
                if len(buffer) >= buffer_size:
                    writer.writerows(buffer)
                    buffer.clear()  # Reset buffer

            except ValueError as ve:
                logger.warning(f"Skipping malformed row: {row} (ValueError: {ve})")
                continue  # Skip invalid rows gracefully

        # Final flush: Write any remaining data in the buffer
        if buffer:
            writer.writerows(buffer)

    return output_file


def collect_genome_ids(report_path: str, div: str = "___") -> Set[str]:
    """
    Collect genome IDs from a TSV file where the hit identifier is in the first column.

    Args:
        report_path (str): Path to the TSV file.
        div (str): Divider in proteinID to split out genomeID.

    Returns:
        set: Unique genome IDs

    Example Output:
        {"GCF_000001405.39", ...}
    """
    genome_ids = set()

    with open(report_path, "r") as tsvfile:
        for line in tsvfile:
            line = line.strip()  # Remove leading/trailing whitespace
            if line:  # Ensure the line is not empty
                columns = line.split("\t")
                key = columns[0].split(div)[0]
                genome_ids.add(key)
    return genome_ids


def filter_genome_ids_with_existing_files(genome_ids, faa_files, gff_files):
    """
    Entfernt alle genomeIDs, für die entweder FAA oder GFF nicht als existierende Datei vorliegt.

    Args:
        genome_ids (Iterable[str]): Zu prüfende genomeIDs.
        faa_files (dict[str, str]): Map genomeID → faa-Dateipfad.
        gff_files (dict[str, str]): Map genomeID → gff-Dateipfad.

    Returns:
        set[str]: Gefilterte genomeIDs mit existierenden FAA- und GFF-Files.
    """
    valid_ids = set()
    missing_ids = []

    for gid in genome_ids:
        faa = faa_files.get(gid)
        gff = gff_files.get(gid)
        if faa and gff and os.path.isfile(faa) and os.path.isfile(gff):
            valid_ids.add(gid)
        else:
            missing_ids.append(gid)

    logger.info(
        f"Filtered genomeIDs: {len(valid_ids)} valid, {len(missing_ids)} with missing files (FAA/GFF)"
    )
    if missing_ids:
        logger.debug(
            f"GenomeIDs removed due to missing files: {', '.join(missing_ids)}"
        )

    return valid_ids


def fix_sseqid_prefix_suffix_uniqueness(filtered_file: str, options) -> str:
    """
    Checks if sseqid prefix matches queued_genomes and that (suffix, col2) pairs are unique.
    If not, prepends the prefix and replaces the file in place.
    """
    logger.info(
        "Validating genomeIDs and proteinIDs in the filtered blast result table"
    )

    from collections import Counter

    temp_file = filtered_file + ".tmp"
    rows = []
    prefixes = []
    suffixes = []
    col2_vals = []
    queued_genomes = options.queued_genomes  # set of genomeIDs

    with open(filtered_file, "r") as infile:
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            sseqid = row[0]
            # Split at the FIRST occurrence of ___
            parts = sseqid.split("___", 1)
            if len(parts) == 2:
                prefix, suffix = parts
            else:
                prefix = sseqid
                suffix = ""
            prefixes.append(prefix)
            suffixes.append(suffix)
            col2_vals.append(row[1] if len(row) > 1 else "")
            rows.append(row)

    # Schritt 1: Prefix-Check
    prefix_in_queued = [p for p in prefixes if p in queued_genomes]
    if not prefix_in_queued:
        logger.error(
            "No prefix in any sseqid matches a genomeID in the provided fasta files"
        )
        raise ValueError("No prefix in sseqid matches a genomeID in queued_genomes")

    # Schritt 2: Suffix/Col2-Pair-Check
    pair_counter = Counter(zip(suffixes, col2_vals))
    duplicated_pairs = [
        pair for pair, count in pair_counter.items() if count > 1 and pair[0]
    ]
    if duplicated_pairs:
        logger.warning(f"Duplicated proteinID/query pairs found: {duplicated_pairs}")

    has_duplicate_pair = bool(duplicated_pairs)
    if not has_duplicate_pair:
        logger.info(
            "All genomeID, proteinID and query combination generate unique identifiers."
        )
        return filtered_file

    # Korrektur: Prepende Prefix und '___' zu sseqid, dann Datei ersetzen
    with open(temp_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        for row, prefix in zip(rows, prefixes):
            row[0] = f"{prefix}___{row[0]}"
            writer.writerow(row)
    os.replace(temp_file, filtered_file)
    logger.info(
        f"Sequence identifier entries fixed with genomeID prefix in {filtered_file}"
    )
    return filtered_file


def split_genome_ids_into_batches(
    genome_ids_set: List[str], num_batches: int
) -> List[List[str]]:
    """
    Splits a list of genomeIDs into nearly equal-sized batches.

    Args:
        genome_ids_set (list): The list of genomeIDs to split.
        num_batches (int): Number of batches.

    Returns:
        list of lists: Each a batch of genomeIDs

    Example Output:
        [["GCF1", "GCF2"], ["GCF3", "GCF4"], ...]
    """

    # Calculate the approximate size of each batch
    batch_size = len(genome_ids_set) // num_batches
    remainder = len(genome_ids_set) % num_batches

    # Create the batches, distributing the remainder across the first few batches
    batches = []
    start = 0
    for i in range(num_batches):
        # Distribute the remainder across the first 'remainder' batches
        end = start + batch_size + (1 if i < remainder else 0)
        batches.append(genome_ids_set[start:end])
        start = end

    return batches


def init_worker(
    diamond_blast_results_table: str,
    score_threshold_diction: dict,
    csb_patterns_diction: dict,
    csb_pattern_names: list,
) -> None:
    """
    Sets up global variables for each worker of the parsing process.
    """
    global global_diamond_blast_results_table
    global global_score_threshold_diction
    global global_csb_patterns_diction
    global global_csb_pattern_names

    # Assign the arguments to the global variables in worker processes
    global_diamond_blast_results_table = diamond_blast_results_table
    global_score_threshold_diction = score_threshold_diction
    global_csb_patterns_diction = csb_patterns_diction
    global_csb_pattern_names = csb_pattern_names


def process_writer(queue: Any, options: Any, counter: Any) -> None:
    """
    Handles output of search and writes to database.

    Args:
        queue: Multiprocessing queue with parsed results
        options: Options object
        counter: Shared counter

    Output Example:
        Writes proteins, clusters, genomeIDs to database and updates files.
    """
    genome_id_batch = set()
    protein_batch = {}
    cluster_batch = {}
    batch_size = options.glob_chunks
    batch_counter = 0

    while True:
        tup = queue.get()
        if tup is None:
            # Process any remaining data in the batch
            if protein_batch and cluster_batch:
                submit_batches(protein_batch, cluster_batch, genome_id_batch, options)
            break

        else:
            counter.value += 1
            print(f"[INFO] Processed {counter.value} genomes ", end="\r")  #

        genome_id, protein_dict, cluster_dict = tup

        # Concatenate the data
        genome_id_batch.add(genome_id)
        protein_batch.update(protein_dict)
        cluster_batch.update(cluster_dict)
        batch_counter += 1

        # If batch size is reached, process the batch
        if batch_counter >= batch_size:
            submit_batches(protein_batch, cluster_batch, genome_id_batch, options)
            protein_batch = {}
            cluster_batch = {}
            batch_counter = 0

    return


def submit_batches(
    protein_batch: Dict, cluster_batch: Dict, genome_id_set: Set[str], options: Any
) -> None:
    """
    Submits batch results to the database and writes files.

    Args:
        protein_batch: All proteins to write (dict)
        cluster_batch: All clusters to write (dict)
        genome_id_set: All genome IDs in batch (set)
        options: Options object
    """

    # Insert into the database
    database.insert_database_genome_ids(options.database_directory, genome_id_set)
    database.insert_database_proteins(options.database_directory, protein_batch)
    database.insert_database_clusters(options.database_directory, cluster_batch)

    # Append to the gene clusters file
    with open(options.gene_clusters_file, "a") as file:
        for clusterID, cluster in cluster_batch.items():
            domains = cluster.types
            file.write(clusterID + "\t" + "\t".join(domains) + "\n")


def process_parallel_bulk_parse_batch(args: Tuple[Any, Any, Any, Any]) -> None:
    """
    Processes a batch of genomeIDs in parallel (bulk parsing).

    Args:
        args: (data_queue, genomeID_batch, options, counter)

    Input Example:
        genomeID_batch = ["GCF1", "GCF2"]
    """
    data_queue, genomeID_batch, options, counter = args

    # Get global variables for read only
    diamond_blast_results_table = global_diamond_blast_results_table
    score_threshold_diction = global_score_threshold_diction
    csb_patterns_diction = global_csb_patterns_diction
    csb_pattern_names = global_csb_pattern_names

    # Process each genomeID in the batch
    for genomeID in genomeID_batch:
        process_single_genome(
            data_queue,
            genomeID,
            options,
            diamond_blast_results_table,
            score_threshold_diction,
            csb_patterns_diction,
            csb_pattern_names,
            counter,
        )


def process_single_genome(
    data_queue: Any,
    genomeID: str,
    options: Any,
    report: str,
    score_threshold_diction: Dict,
    csb_patterns_diction: Dict,
    csb_pattern_names: List,
    counter: Any,
) -> None:
    """
    Processes a single genomeID, extracts hits, filters, stores the results.

    Args:
        data_queue: Multiprocessing queue for results
        genomeID (str): Genome ID to process
        options: Options object
        report (str): Path to BLAST result table
        score_threshold_diction (dict): Score thresholds
        csb_patterns_diction (dict): CSB patterns
        csb_pattern_names (list): CSB names
        counter: Shared multiprocessing.Value counter
    """

    faa_file = ""

    try:
        # Unpack required files if necessary
        faa_file = unpackgz(options.faa_files[genomeID])
        gff_file = unpackgz(options.gff_files[genomeID])

        # Parse BLAST report for protein hits
        protein_dict = parse_reports.parse_bulk_blastreport_genomize(
            genomeID, report, score_threshold_diction, options.thrs_score
        )

        # Remove multidomain proteins if they are not allowed
        if not options.multidomain_allowed:
            protein_dict = remove_multi_domain_proteins(protein_dict)

        # if options sagt, dass es einzelgenome sind, dann die genomID aus den keys und den protein identifiern entfernen
        if not options.glob_gff:
            protein_dict = parse_reports.clean_dict_keys_and_protein_ids(
                protein_dict, genomeID
            )

        # Complete the hit information
        parse_reports.parse_gff_file(gff_file, protein_dict)
        parse_reports.get_protein_sequence(faa_file, protein_dict)

        # Find the syntenic regions and insert to database
        cluster_dict = csb_finder.find_syntenic_blocks(
            genomeID, protein_dict, options.nucleotide_range
        )

        # Store parsed data in multiprocessing queue
        data_queue.put((genomeID, protein_dict, cluster_dict))

    except Exception as e:
        logger.warning(f"Skipped {faa_file} due to an error - {e}")
        logger.debug(traceback.format_exc())
        return

    return


def remove_multi_domain_proteins(input_dict):
    # returns only the key:protein pairs that have not more than 1 domain
    return {
        key: protein
        for key, protein in input_dict.items()
        if len(protein.domains.values()) <= 1
    }


def initial_glob_search(options: Any) -> None:
    """
    This function initializes and performs a DIAMOND BLASTp search,
    filters the results, processes genome hits in parallel, and writes parsed data.

    Steps:
    1. Perform self BLAST query to establish baseline.
    2. If a BLAST table is not provided, run DIAMOND BLASTp.
    3. Filter BLAST results based on e-value, sequence identity, and score thresholds.
    4. Extract unique genome IDs from filtered hits and insert them into the database.
    5. Batch process the genome hits in parallel using multiprocessing.

    Args:
        options: Options object with paths and search params.

    Example Input:
        options.glob_faa: str = "glob.faa"
        options.query_file: str = "query.faa"
        options.cores: int = 8

    Output:
        Writes filtered results to options.filtered_blast_table,
        database gets all hits & clusters inserted.
    """

    # Step 1: Perform self BLAST query to establish reference scores for blast score ratio
    logger.info("Initilize DIAMOND BLASTp self-blast against query")
    selfblast_scores_dict, query_length_dict = query_selfblast_search.self_blast_query(
        options
    )

    # Step 2: Run DIAMOND BLASTp if no precomputed BLAST table is provided or exists
    blast_results_table = options.glob_table
    if not blast_results_table:
        logger.info("Initilize DIAMOND BLASTp against target")
        blast_results_table = diamond_search.diamond_search(
            options.glob_faa,
            options.query_file,
            options.cores,
            options.evalue,
            options.searchcoverage,
            options.minseqid,
            options.diamond_report_hits_limit,
            options.alignment_mode,
        )
    else:
        logger.info(f"Using DIAMOND BLASTp result table {blast_results_table}")

        # Step 3: Filter BLAST results based on e-value, sequence identity, and score thresholds
    logger.info(
        "Filtering raw DIAMOND BLASTp results by score, e-value, identity and blast score ratio parameters"
    )

    # Get baseline BLAST scores

    blast_results_table = filter_blast_table(
        output_file=options.filtered_blast_table,
        blast_file=blast_results_table,
        evalue_cutoff=options.evalue,
        score_cutoff=options.thrs_score,
        coverage_cutoff=options.searchcoverage,
        identity_cutoff=options.minseqid,
        bsr_cutoff=options.thrs_bsr,
        sequence_lengths=query_length_dict,
        selfblast_scores=selfblast_scores_dict,
    )

    # Store the filtered table in the options object
    options.glob_table = blast_results_table

    # Step 4: Extract and store unique genome IDs
    genome_ids_set = collect_genome_ids(
        blast_results_table
    )  # returns a set of all genomeIDs
    logger.info(f"Found hits in {len(genome_ids_set)} genomes in {blast_results_table}")

    # Remove genomeIDs without faa and gff file. It is possible that a blast result table has more/other IDs as the provided files
    genome_ids_set = filter_genome_ids_with_existing_files(
        genome_ids_set, options.faa_files, options.gff_files
    )

    database.insert_database_genome_ids(
        options.database_directory, genome_ids_set
    )  # Insert genomeIDs into database

    # Step 5 Check if the filterd fasta file has correct genomeIDs ___ proteinID combinations with unique ids
    # options.glob_table = fix_sseqid_prefix_suffix_uniqueness(options.glob_table, options)

    # Step 6: Process genome hits in parallel
    logger.info("Parsing hits from filtered results table")
    genome_id_batches = split_genome_ids_into_batches(
        list(genome_ids_set), options.cores - 1
    )  # batch sets of genomeIDs for each worker

    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value("i", 0)

    csb_patterns_diction, csb_pattern_names = (
        {},
        {},
    )  # csb_finder.make_pattern_dict(options.patterns_file)
    score_threshold_diction = {}  # empty because all have to be parsed

    with Pool(
        processes=options.cores,
        initializer=init_worker,
        initargs=(
            blast_results_table,
            score_threshold_diction,
            csb_patterns_diction,
            csb_pattern_names,
        ),
    ) as pool:
        p_writer = pool.apply_async(process_writer, (data_queue, options, counter))

        pool.map(
            process_parallel_bulk_parse_batch,
            [(data_queue, batch, options, counter) for batch in genome_id_batches],
        )

        # Signal the end of data to the process_writer
        for _ in range(options.cores):
            data_queue.put(None)

        p_writer.get()
    logger.info("Finished parsing BLASTp results")
    return
