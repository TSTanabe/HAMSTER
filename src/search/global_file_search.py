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
from src.parse_reports.parse_reports import Protein
from src.search import query_selfblast_search, blastp_n_filter
from src.fasta_preparation.materialize import materialize_pair_gz_next_to_input
from src.search import diamond_search

logger = get_logger(__name__)


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


def init_worker(
    blast_results_table,
    score_threshold_diction,
    csb_patterns_diction,
    csb_pattern_names,
    faa_files,
    gff_files,
    multidomain_allowed,
    glob_gff,
    nucleotide_range,
    thrs_score,
):
    global global_diamond_blast_results_table
    global global_score_threshold_diction
    global global_csb_patterns_diction
    global global_csb_pattern_names
    global global_faa_files
    global global_gff_files
    global global_multidomain_allowed
    global global_glob_gff
    global global_nucleotide_range
    global global_thrs_score

    global_thrs_score = thrs_score
    global_diamond_blast_results_table = blast_results_table
    global_score_threshold_diction = score_threshold_diction
    global_csb_patterns_diction = csb_patterns_diction
    global_csb_pattern_names = csb_pattern_names
    global_faa_files = faa_files
    global_gff_files = gff_files
    global_multidomain_allowed = multidomain_allowed
    global_glob_gff = glob_gff
    global_nucleotide_range = nucleotide_range


def remove_multi_domain_proteins(input_dict):
    # returns only the key:protein pairs that have not more than 1 domain
    return {
        key: protein
        for key, protein in input_dict.items()
        if len(protein.domains.values()) <= 1
    }

def process_single_genome_imap(genomeID: str):
    """
    Worker für pool.imap_unordered.
    Rückgabe ist exakt das, was bisher in die Queue ging:
      (genomeID, protein_dict, cluster_dict)
    Bei Fehler: None (damit der Main-Prozess einfach skippen kann)
    """
    # Globals aus init_worker()
    report = global_diamond_blast_results_table
    score_threshold_diction = global_score_threshold_diction

    try:
        faa_in = global_faa_files[genomeID]
        gff_in = global_gff_files[genomeID]

        with materialize_pair_gz_next_to_input(faa_in, gff_in) as (faa_file, gff_file):
           # Parse BLAST report for protein hits
            protein_dict = parse_reports.parse_bulk_blastreport_genomize(
                genomeID, report, score_threshold_diction, global_thrs_score
            )

            # Remove multidomain proteins if they are not allowed
            if not global_multidomain_allowed:
                protein_dict = remove_multi_domain_proteins(protein_dict)

            # Einzelgenom-Modus: genomeID aus Keys/ProteinIDs entfernen
            if not global_glob_gff:
                protein_dict = parse_reports.clean_dict_keys_and_protein_ids(
                    protein_dict, genomeID
                )

            # Complete hit information
            parse_reports.parse_gff_file(gff_file, protein_dict)
            parse_reports.get_protein_sequence(faa_file, protein_dict)

            # Syntenie/Cluster
            cluster_dict = csb_finder.find_syntenic_blocks(
                genomeID, protein_dict, global_nucleotide_range
            )

            return genomeID, protein_dict, cluster_dict

    except Exception as e:
        logger.warning(f"Skipped {faa_file} due to an error - {e}")
        logger.debug(traceback.format_exc())
        return None


def initial_glob_search(config: Any) -> None:
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
        config
    )

    # Step 2+3: Run DIAMOND BLASTp if no precomputed BLAST table is provided or exists
    blast_results_table = blastp_n_filter.run_and_filter_diamond_blastp(
        config,
        query_length_dict=query_length_dict,
        selfblast_scores_dict=selfblast_scores_dict,
    )

    # Step 4: Extract and store unique genome IDs
    genome_ids_set = collect_genome_ids(
        blast_results_table
    )  # returns a set of all genomeIDs
    logger.info(f"Found hits in {len(genome_ids_set)} genomes in {blast_results_table}")

    # Remove genomeIDs without faa and gff file. It is possible that a blast result table has more/other IDs as the provided files
    genome_ids_set = filter_genome_ids_with_existing_files(
        genome_ids_set, config.faa_files, config.gff_files
    )

    database.insert_database_genome_ids(
        config.database_directory, genome_ids_set
    )  # Insert genomeIDs into database

    # Step 5 Check if the filterd fasta file has correct genomeIDs ___ proteinID combinations with unique ids
    # options.glob_table = fix_sseqid_prefix_suffix_uniqueness(options.glob_table, options)

    # Step 6: Process genome hits in parallel
    logger.info("Parsing hits from filtered results table")
    genome_id_list = list(genome_ids_set)

    n_genomes: int = len(genome_id_list)
    worker_processes = max(1, (config.cores - 1))
    worker_processes = min(worker_processes, n_genomes)

    # Batch buffers im Main
    protein_batch: dict[str, Protein] = {}
    cluster_batch: dict = {}
    batch_size: int = config.glob_chunks
    batch_counter = 0

    # Fortschritt
    genomes_done = 0
    log_step = max(1, n_genomes // 100)

    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value("i", 0)

    csb_patterns_diction, csb_pattern_names = (
        {},
        {},
    )  # csb_finder.make_pattern_dict(options.patterns_file)
    score_threshold_diction = {}  # empty because all have to be parsed

    with Pool(
        processes=worker_processes,
        initializer=init_worker,
        initargs=(
            blast_results_table,
            score_threshold_diction,
            csb_patterns_diction,
            csb_pattern_names,
            config.faa_files,
            config.gff_files,
            config.multidomain_allowed,
            config.glob_gff,
            config.nucleotide_range,
            config.thrs_score,
        ),
    ) as pool:
        it = pool.imap_unordered(process_single_genome_imap, genome_id_list, chunksize=1)
        for genomeID, protein_dict, cluster_dict in it:
            genomes_done += 1
            if (genomes_done % log_step == 0) or (genomes_done == n_genomes):
                pct = (genomes_done * 100) // max(1, n_genomes)
                logger.info(
                    f"[Progress] {genomes_done}/{n_genomes} ({pct}%) genomes processed"
                )

            # Batch sammeln
            protein_batch.update(protein_dict)
            cluster_batch.update(cluster_dict)
            batch_counter += 1
            if batch_counter >= batch_size:
                database.insert_database_proteins(
                    config.database_directory, protein_batch
                )
                database.insert_database_clusters(
                    config.database_directory, cluster_batch
                )
                protein_batch.clear()
                cluster_batch.clear()
                batch_counter = 0
                with open(config.gene_clusters_file, "a") as file:
                    for clusterID, cluster in cluster_batch.items():
                        domains = cluster.types
                        file.write(clusterID + "\t" + "\t".join(domains) + "\n")

        if protein_batch or cluster_batch:
            database.insert_database_proteins(config.database_directory, protein_batch)
            database.insert_database_clusters(config.database_directory, cluster_batch)
            with open(config.gene_clusters_file, "a") as file:
                for clusterID, cluster in cluster_batch.items():
                    domains = cluster.types
                    file.write(clusterID + "\t" + "\t".join(domains) + "\n")

    logger.info("Finished parsing BLASTp results")
    return
