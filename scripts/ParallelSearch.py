#!/usr/bin/python

import os
import csv
import re
import traceback
import subprocess
import multiprocessing
from multiprocessing import Manager, Pool
from typing import List, Set, Dict, Any, Tuple, Optional

from . import myUtil
from . import Csb_finder
from . import Database
from . import ParseReports
from . import Alignment

logger = myUtil.logger





################################################################################################        
################################## Genomize search #############################################
################################################################################################


def initial_genomize_search(options: Any) -> None:
    """
    Initializes and performs a DIAMOND BLASTp search per genome (genomize mode).
    Uses multiprocessing for parallel search and collects all hits and clusters.

    Args:
        options: The options object with all search parameters and queue settings.

    Example Input:
        options.queued_genomes: set[str] = {"GCF_...", ...}
        options.faa_files: dict[str, str] genomeID => filepath
        options.gff_files: dict[str, str] genomeID => filepath
        options.cores: int = 8

    Output:
        Populates options.database_directory with parsed and clustered search results.
    """
     
    score_threshold_diction = Search.makeThresholdDict(options.score_threshold_file, options.threshold_type)
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
    self_blast_query(options)
    
    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value('i', 0)
   
    process_instances = []
    
    with Pool(processes = options.cores) as pool:
    
        p_writer = pool.apply_async(process_writer, (data_queue, options))
        
        pool.map(process_parallel_search, [(data_queue, genomeID, options, counter, score_threshold_diction, csb_patterns_diction,csb_pattern_names) for genomeID in options.queued_genomes])

        # Signal the end of data to the process_writer
        for _ in range(options.cores):
            data_queue.put(None)
        
        p_writer.get()
    logger.info(f"Finished {counter} DIAMOND BLASTp searches")
    logger.info(f"Saved BLASTp hits to database {options.database_directory}")
    logger.info(f"Individual hit lists to {options.fasta_initial_hit_directory}")

    return    

def process_parallel_search(args_tuple: Tuple) -> None:
    """
    Worker for parallelized DIAMOND search per genome.

    Args:
        args_tuple: Tuple containing (data_queue, genomeID, options, counter, ... dicts)

    Input Example:
        genomeID: "GCF_000001405.39"
        options: Options object

    Output:
        Puts (genomeID, protein_dict, cluster_dict) to multiprocessing queue.
    """
    queue,genomeID,options,counter,score_threshold_diction, csb_patterns_diction,csb_pattern_names = args_tuple
    with counter.get_lock():
        counter.value += 1
        print(f"[INFO] Searching assembly ({counter.value}/{len(options.queued_genomes)})", end='\r', flush=True)
    cluster_dict = {}
    try:
        faa_file = myUtil.unpackgz(options.faa_files[genomeID])
        gff_file = myUtil.unpackgz(options.gff_files[genomeID])
        cores = 1 # worker process has only one core
        report = DiamondSearch(options.glob_faa, options.query_file, cores, options.evalue, options.searchcoverage, options.minseqid, options.diamond_report_hits_limit)
        
        #Parse the hits
        protein_dict = parse_bulk_blastreport_genomize(genomeID,report,score_threshold_diction,options.thrs_score)
        
        #Complete the hit information
        ParseReports.parseGFFfile(gff_file,protein_dict)
        ParseReports.getProteinSequence(faa_file,protein_dict)
        
        #Find the syntenic regions and insert to database
        cluster_dict = Csb_finder.find_syntenicblocks(genomeID,protein_dict,options.nucleotide_range)
        Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_dict,1)

        data_queue.put((genomeID,protein_dict,cluster_dict))
    except Exception as e:
        logger.warning(f"Skipped {faa_file} due to an error - {e}")
        logger.debug(traceback.format_exc())
        return

    return

    
    
    
     

    
################################################################################################        
############################################ Glob search #######################################
################################################################################################

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
    self_blast_report = self_blast_query(options)

    # Step 2: Run DIAMOND BLASTp if no precomputed BLAST table is provided or exists
    blast_results_table = options.glob_table
    if not blast_results_table:
        logger.info("Initilize DIAMOND BLASTp against target")
        blast_results_table = DiamondSearch(
            options.glob_faa, options.query_file, options.cores, 
            options.evalue, options.searchcoverage, options.minseqid, 
            options.diamond_report_hits_limit, options.alignment_mode,
        )
    
    # Step 3: Filter BLAST results based on e-value, sequence identity, and score thresholds
    logger.info("Filtering raw DIAMOND BLASTp results by score, e-value, identity and blast score ratio parameters")
       
    query_length_dict = get_sequence_legth(options.self_query) # Retrieve sequence lengths
    selfblast_scores_dict = get_sequence_hits_scores(self_blast_report) # Get baseline BLAST scores
    
    blast_results_table = filter_blast_table(
        options.filtered_blast_table,
        blast_results_table, options.evalue, options.thrs_score, 
        options.searchcoverage, options.minseqid, options.thrs_bsr, 
        query_length_dict, selfblast_scores_dict
    )
    
    # Store the filtered table in the options object
    options.glob_table = blast_results_table

    # Step 4: Extract and store unique genome IDs
    genomeIDs_set = collect_genomeIDs(blast_results_table) #returns a set of all genomeIDs
    logger.info(f"Found hits in {len(genomeIDs_set)} genomes")
    
    Database.insert_database_genomeIDs(options.database_directory, genomeIDs_set) # Insert genomeIDs into database
    
    # Step 5: Process genome hits in parallel
    logger.info("Parsing hits from filtered results table")
    genomeID_batches = split_genomeIDs_into_batches(list(genomeIDs_set), options.cores-1) # batch sets of genomeIDs for each worker
    
    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value('i', 0)
    
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
    score_threshold_diction = {} #empty because all have to be parsed
    
    with Pool(processes = options.cores, initializer=init_worker, 
              initargs=(blast_results_table, score_threshold_diction, csb_patterns_diction, csb_pattern_names)) as pool:
    
        p_writer = pool.apply_async(process_writer, (data_queue, options, counter))
        
        pool.map(process_parallel_bulk_parse_batch, [(data_queue, batch, options, counter) for batch in genomeID_batches])

        # Signal the end of data to the process_writer
        for _ in range(options.cores):
            data_queue.put(None)
        
        p_writer.get()
    logger.info("Finished parsing BLASTp results")
    return


# Initializer function to set up global variables for each worker of the parsing process
def init_worker(diamond_blast_results_table: str, score_threshold_diction: dict, csb_patterns_diction: dict, csb_pattern_names: list) -> None:
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
        
        

    


####### Search routine #############

def DiamondSearch(path: str, query_fasta: str, cores: int, evalue: float, coverage: float, minseqid: float, diamond_report_hits_limit: int, alignment_mode: int = 2, sensitivity: str = "ultra-sensitive") -> str:
    """
    Runs DIAMOND BLASTp search on a FASTA file.

    Args:
        path (str): Path to the FASTA file (target db)
        query_fasta (str): Query FASTA
        cores (int): Number of cores
        evalue (float): E-value cutoff
        coverage (float): Min coverage
        minseqid (float): Min sequence identity
        diamond_report_hits_limit (int): Report hits limit
        alignment_mode (int): DIAMOND alignment mode
        sensitivity (str): Sensitivity level

    Returns:
        str: Path to output .diamond.tab file

    Example Output:
        "results/glob.faa.diamond.tab"
    """

    diamond = myUtil.find_executable("diamond")
    target_db_name = f"{path}.dmnd"
    os.system(f'{diamond} makedb --quiet --in {path} -d {target_db_name} --threads {cores} 1>/dev/null 0>/dev/null')
    
    output_results_tab = f"{path}.diamond.tab"
    #{hit_id}\t{query_id}\t{e_value}\t{score}\t{bias}\t{hsp_start}\t{hsp_end}\t{description}
    
    logger.debug(f'{diamond} blastp --quiet --{sensitivity} -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} --outfmt 6 sseqid qseqid evalue bitscore sstart send pident 1>/dev/null 0>/dev/null')
    os.system(f'{diamond} blastp --quiet --{sensitivity} -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} -k {diamond_report_hits_limit} --outfmt 6 sseqid qseqid evalue bitscore sstart send pident 1>/dev/null 0>/dev/null')
    #output format hit query evalue score identity alifrom alito
    return output_results_tab

    

def collect_genomeIDs(report_path: str, div: str = '___') -> Set[str]:
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
    genomeIDs = set()
    
    with open(report_path, 'r') as tsvfile:
        for line in tsvfile:
            line = line.strip()  # Remove leading/trailing whitespace
            if line:  # Ensure the line is not empty
                columns = line.split('\t')
                key = columns[0].split(div)[0]
                genomeIDs.add(key)
    return genomeIDs
    
def split_genomeIDs_into_batches(genomeIDs_list: List[str], num_batches: int) -> List[List[str]]:
    """
    Splits a list of genomeIDs into nearly equal-sized batches.

    Args:
        genomeIDs_list (list): The list of genomeIDs to split.
        num_batches (int): Number of batches.

    Returns:
        list of lists: Each a batch of genomeIDs

    Example Output:
        [["GCF1", "GCF2"], ["GCF3", "GCF4"], ...]
    """
    
    # Calculate the approximate size of each batch
    batch_size = len(genomeIDs_list) // num_batches
    remainder = len(genomeIDs_list) % num_batches

    # Create the batches, distributing the remainder across the first few batches
    batches = []
    start = 0
    for i in range(num_batches):
        # Distribute the remainder across the first 'remainder' batches
        end = start + batch_size + (1 if i < remainder else 0)
        batches.append(genomeIDs_list[start:end])
        start = end

    return batches
    
############################################################################
################# Bulk parsing routines ####################################
############################################################################

def process_parallel_bulk_parse_batch(args: Tuple[Any, Any, Any, Any]) -> None:
    """
    Processes a batch of genomeIDs in parallel (bulk parsing).

    Args:
        args: (data_queue, genomeID_batch, options, counter)

    Input Example:
        genomeID_batch = ["GCF1", "GCF2"]
    """
    data_queue, genomeID_batch, options, counter = args
    
    #Get global variables for read only
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
            counter
        )
   
def process_single_genome(
    data_queue: Any,
    genomeID: str,
    options: Any,
    report: str,
    score_threshold_diction: Dict,
    csb_patterns_diction: Dict,
    csb_pattern_names: List,
    counter: Any
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
    protein_dict = {}
    cluster_dict = {}
    
    try:
        # Unpack required files if necessary
        faa_file = myUtil.unpackgz(options.faa_files[genomeID])
        gff_file = myUtil.unpackgz(options.gff_files[genomeID])
        
        #Parse BLAST report for protein hits
        protein_dict = parse_bulk_blastreport_genomize(genomeID,report,score_threshold_diction,options.thrs_score)
        
        #Remove multidomain proteins if they are not allowed
        if not options.multidomain_allowed:
            protein_dict = remove_multi_domain_proteins(protein_dict)
        
        #if options sagt, dass es einzelgenome sind, dann die genomID aus den keys und den protein identifiern entfernen
        if not options.glob_gff:
            protein_dict = clean_dict_keys_and_protein_ids(protein_dict, genomeID)
        
        #Complete the hit information
        ParseReports.parseGFFfile(gff_file,protein_dict)
        ParseReports.getProteinSequence(faa_file,protein_dict)
        
        #Find the syntenic regions and insert to database
        cluster_dict = Csb_finder.find_syntenicblocks(genomeID,protein_dict,options.nucleotide_range)
        Csb_finder.name_syntenicblocks(csb_patterns_diction,csb_pattern_names,cluster_dict,1)
        
        # Store parsed data in multiprocessing queue
        data_queue.put((genomeID,protein_dict,cluster_dict))
        
    except Exception as e:
        logger.warning(f"Skipped {faa_file} due to an error - {e}")
        logger.debug(traceback.format_exc())
        return


    return


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
    genomeID_batch = set()
    protein_batch = {}
    cluster_batch = {}
    batch_size = options.glob_chunks
    batch_counter = 0

    while True:
        tup = queue.get()
        if tup is None:
            # Process any remaining data in the batch
            if protein_batch and cluster_batch:
                submit_batches(protein_batch, cluster_batch, genomeID_batch, options)
            break
        
        else:
            counter.value += 1
            print(f"[INFO] Processed {counter.value} genomes ", end="\r")#
        
        genomeID, protein_dict, cluster_dict = tup

        # Concatenate the data
        genomeID_batch.add(genomeID)
        protein_batch.update(protein_dict)
        cluster_batch.update(cluster_dict)
        batch_counter += 1
        
        # If batch size is reached, process the batch
        if batch_counter >= batch_size:
            submit_batches(protein_batch, cluster_batch, genomeID_batch, options)
            protein_batch = {}
            cluster_batch = {}
            batch_counter = 0

    return

def submit_batches(protein_batch: Dict, cluster_batch: Dict, genomeID_set: Set[str], options: Any) -> None:
    """
    Submits batch results to the database and writes files.

    Args:
        protein_batch: All proteins to write (dict)
        cluster_batch: All clusters to write (dict)
        genomeID_set: All genome IDs in batch (set)
        options: Options object
    """
    
    # Insert into the database
    Database.insert_database_genomeIDs(options.database_directory, genomeID_set)
    Database.insert_database_proteins(options.database_directory, protein_batch)
    Database.insert_database_clusters(options.database_directory, cluster_batch)
    
    # Append to the gene clusters file
    with open(options.gene_clusters_file, "a") as file:
        for clusterID, cluster in cluster_batch.items():
            domains = cluster.types
            file.write(clusterID + '\t' + '\t'.join(domains) + '\n')


def write_query_hit_sequence_fasta(directory, protein_dict):
    
    sorted_proteins = sorted(protein_dict.values(), key=lambda x: x.get_domains())
    
    current_file_path = None
    current_file = None

    try:
        for protein in sorted_proteins:
            genomeID = protein.genomeID
            proteinID = f"{genomeID}-{protein.proteinID}" #same concatenation is done for the database, do not alter
            domains = protein.get_domains()
            file_path = f"{directory}/Query_{domains}.faa"
            
            # If the file path changes, close the current file and open a new one
            if file_path != current_file_path:
                if current_file:
                    current_file.close()
                current_file_path = file_path
                current_file = open(file_path, "a")
            
            # Write the protein data to the file
            current_file.write(f">{proteinID}\n")
            current_file.write(f"{protein.protein_sequence}\n")
    finally:
        # Ensure the last file is closed properly
        if current_file:
            current_file.close()

def parse_bulk_blastreport_genomize(genome_id: str, filepath: str, thresholds: Dict, cut_score: int = 10) -> Dict[str, Any]:
    """
    Parses a BLAST report for a specific genomeID and extracts protein domain hits.

    Args:
        genome_id (str): The genome ID to filter hits.
        filepath (str): Path to BLAST report.
        thresholds (dict): Score thresholds (unused here).
        cut_score (int): Score threshold for filtering.

    Returns:
        dict: proteinID -> Protein object (ParseReports.Protein)

    Output Example:
        {"WP_012345678": ParseReports.Protein(...), ...}
    """

    protein_dict = {}

    try:
        # Efficiently filter relevant lines with grep
        with subprocess.Popen(['grep', genome_id, filepath], stdout=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                columns = line.strip().split('\t')  # Assuming tab-separated format

                if len(columns) < 6:  # Ensure enough columns exist
                    continue

                try:
                    # Expected BLAST format: hit_id, query_id, e_value, bitscore, hsp_start, hsp_end
                    hit_protein_id = columns[0]  # Protein identifier in BLAST
                    query_id = columns[1]  # Query sequence ID
                    hit_bitscore = int(float(columns[3]))  # Bitscore as integer
                    hsp_start = int(columns[4])  # Start position
                    hsp_end = int(columns[5])  # End position
                    hsp_ident = int(float(columns[6]))
                    try:
                        hsp_bsr = float(columns[7])
                    except IndexError:
                        hsp_bsr = 1.0  # oder ein anderer sinnvoller Default-Wert


                    # If protein already exists, add domain info
                    if hit_protein_id in protein_dict:
                        protein_dict[hit_protein_id].add_domain(query_id, hsp_start, hsp_end, hit_bitscore)
                    else:
                        # Create new protein object
                        protein_dict[hit_protein_id] = ParseReports.Protein(
                            hit_protein_id, query_id, hsp_start, hsp_end, hit_bitscore, genome_id, hsp_ident, hsp_bsr
                        )

                except ValueError as ve:
                    logger.warning(f"Skipped malformed line in {filepath}: {line.strip()} (ValueError: {ve})")
    except Exception as e:
        logger.error(f"Failed to parse {filepath} for genome {genome_id}: {e}")
        logger.debug(traceback.format_exc())

    return protein_dict
    




def remove_multi_domain_proteins(input_dict):
    #returns only the key:protein pairs that have not more than 1 domain
    return {
        key: protein
        for key, protein in input_dict.items()
        if len(protein.domains.values()) <= 1
    }


def clean_dict_keys_and_protein_ids(input_dict: Dict[str, Any], genomeID: str) -> Dict[str, Any]:
    """
    Removes the genomeID prefix from dictionary keys and proteinID fields.

    Args:
        input_dict: dict of proteinID -> Protein
        genomeID: genomeID prefix to remove

    Returns:
        dict: updated protein dict
    """
    prefix = genomeID + '___'
    updated_dict = {}
    
    for key, protein in input_dict.items():
        # Entferne das Präfix von jedem Key, wenn es vorhanden ist
        new_key = key[len(prefix):] if key.startswith(prefix) else key
        
        # Entferne das Präfix von proteinID, falls es vorhanden ist
        if hasattr(protein, 'proteinID') and protein.proteinID.startswith(prefix):
            protein.proteinID = protein.proteinID[len(prefix):]
        
        # Füge das geänderte Key-Value-Paar zum neuen Dictionary hinzu
        updated_dict[new_key] = protein
    
    return updated_dict





################################################################################################        
################################## Filter blast hit report #####################################
################################################################################################


def get_sequence_legth(file_path: str) -> Dict[str, float]:
    """
    Gets sequence length per query name from FASTA (averaged if multiple copies).

    Args:
        file_path (str): Path to FASTA

    Returns:
        dict: query_name -> average sequence length

    Output Example:
        {"QueryA": 237.0, "QueryB": 180.0}
    """
    
    sequence_data = {}  # Dictionary to store lengths per identifier
    sequence_counts = {}  # Dictionary to track counts of sequences per identifier
    
    with open(file_path, 'r') as fasta_file:
        current_id = None
        current_sequence = []
        
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                # Process previous sequence if any
                if current_id:
                    sequence_length = len("".join(current_sequence))
                    if current_id in sequence_data:
                        sequence_data[current_id] += sequence_length
                        sequence_counts[current_id] += 1
                    else:
                        sequence_data[current_id] = sequence_length
                        sequence_counts[current_id] = 1
                
                # Extract the identifier, considering only the part before "___"
                current_id = line[1:].split('___')[1]
                current_sequence = []  # Reset sequence for the new identifier
            else:
                current_sequence.append(line)
        
        # Process the last sequence in the file
        if current_id:
            sequence_length = len("".join(current_sequence))
            if current_id in sequence_data:
                sequence_data[current_id] += sequence_length
                sequence_counts[current_id] += 1
            else:
                sequence_data[current_id] = sequence_length
                sequence_counts[current_id] = 1
    
    # Calculate average length for identifiers with multiple sequences
    averaged_data = {}
    for key in sequence_data:
        averaged_data[key] = sequence_data[key] / sequence_counts[key]
    
    return averaged_data
    
    
def get_sequence_hits_scores(blast_file: str) -> Dict[str, float]:
    """
    Generates a dict of self-blast scores from a BLAST table file.

    Args:
        blast_file (str): Path to BLASTP table file

    Returns:
        dict: qseqid -> highest bitscore

    Output Example:
        {"Q12345": 180.0, ...}
    """
    selfblast_scores = {}

    with open(blast_file, 'r') as infile:
        for line in infile:
            row = line.strip().split('\t')
            sseqid, qseqid, evalue, bitscore, sstart, send, pident = row

            # Convert bitscore to float for comparison
            bitscore = float(bitscore)

            # Update the dictionary with the highest bitscore for self-hits
            if qseqid in selfblast_scores:
                selfblast_scores[qseqid] = max(selfblast_scores[qseqid], bitscore)
            else:
                selfblast_scores[qseqid] = bitscore

    return selfblast_scores
    
        
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
    buffer_size: int = 10000
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
    #output_file = os.path.join(os.path.dirname(blast_file), f"filtered_{os.path.basename(blast_file)}")
    
    # Open input and output files using CSV reader/writer
    with open(blast_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

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
                        evalue <= evalue_cutoff and
                        bitscore >= score_cutoff and
                        coverage >= coverage_cutoff and
                        pident >= identity_cutoff and
                        bsr >= bsr_cutoff
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

    
################################################################################################        
##################################### Selfblast query ##########################################
################################################################################################

def self_blast_query(options: Any) -> str:
    """
    Performs a self-BLAST of the query file.

    Args:
        options: Options object with .self_query, .query_file, etc.

    Output:
        Writes to database, returns BLAST report path.

    Example Output:
        "results/self_query.diamond.tab"
    """
    #Selfblast the query file
    report = DiamondSearch(options.self_query, options.query_file, options.cores, options.evalue, 100, 100, options.diamond_report_hits_limit) #Selfblast, coverage and identity have to be 100 % or weakly similar domains may occur
    protein_dict = parse_bulk_blastreport_genomize("QUERY",report,{},10) #Selfblast should not have any cutoff score
    
    ParseReports.getProteinSequence(options.self_query,protein_dict) #Get the protein Sequences
    protein_dict = clean_dict_keys_and_protein_ids(protein_dict, "QUERY")
    
    Database.insert_database_genomeIDs(options.database_directory, {"QUERY"})
    Database.insert_database_proteins(options.database_directory, protein_dict)
    
    return report



                        
