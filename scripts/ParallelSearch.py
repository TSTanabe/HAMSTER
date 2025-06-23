#!/usr/bin/python

import os
import csv
import re

from . import myUtil
from . import Csb_finder
from . import Database
from . import ParseReports
from . import Alignment

import traceback
import subprocess
import multiprocessing
from multiprocessing import Manager, Pool







################################################################################################        
################################## Genomize search #############################################
################################################################################################


def initial_genomize_search(options):
    
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
    print("[INFO] Finished DIAMOND BLASTp search")
    print(f"[SAVE] BLASTp hits to database {options.database_directory}")
    print(f"[SAVE] individual hit lists to {options.fasta_initial_hit_directory}")
    return    

def process_parallel_search(args_tuple):

    queue,genomeID,options,counter,score_threshold_diction, csb_patterns_diction,csb_pattern_names = args_tuple
    counter.value += 1
    print(f"[INFO] Searching assembly ({counter.value}/{len(options.queued_genomes)})", end="\r")
    cluster_dict
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
        error_message = f"\nError: occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"[WARN] Skipped {faa_file} due to an error - {error_message}")
        print(f"Traceback details:\n{traceback_details}")
        return

    return

    
    
    
     

    
################################################################################################        
############################################ Glob search #######################################
################################################################################################

def initial_glob_search(options):
    """
    This function initializes and performs a DIAMOND BLASTp search,
    filters the results, processes genome hits in parallel, and writes parsed data.
    
    Steps:
    1. Perform self BLAST query to establish baseline.
    2. If a BLAST table is not provided, run DIAMOND BLASTp.
    3. Filter BLAST results based on e-value, sequence identity, and score thresholds.
    4. Extract unique genome IDs from filtered hits and insert them into the database.
    5. Batch process the genome hits in parallel using multiprocessing.
    """
    
    print(f"[INFO] Initilize DIAMOND BLASTp self-blast against query")
    
    # Step 1: Perform self BLAST query to establish reference scores for blast score ratio
    self_blast_report = self_blast_query(options)

    # Step 2: Run DIAMOND BLASTp if no precomputed BLAST table is provided or exists
    blast_results_table = options.glob_table
    if not blast_results_table:
        print(f"[INFO] Initilize DIAMOND BLASTp against target")
        blast_results_table = DiamondSearch(
            options.glob_faa, options.query_file, options.cores, 
            options.evalue, options.searchcoverage, options.minseqid, 
            options.diamond_report_hits_limit, options.alignment_mode,
        )
    
    # Step 3: Filter BLAST results based on e-value, sequence identity, and score thresholds
    print("[INFO] Filtering raw BLASTp results by score, e-value, identity and blast score ratio parameters")   
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
    print(f"[INFO] Found hits in {len(genomeIDs_set)} genomes")
    
    Database.insert_database_genomeIDs(options.database_directory, genomeIDs_set) # Insert genomeIDs into database
    
    # Step 5: Process genome hits in parallel
    print("[INFO] Parsing hits from filtered results table per genome")
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
    print("[INFO] Finished parsing BLASTp results")
    return


# Initializer function to set up global variables for each worker of the parsing process
def init_worker(diamond_blast_results_table, score_threshold_diction, csb_patterns_diction, csb_pattern_names):
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

def DiamondSearch(path, query_fasta, cores, evalue, coverage, minseqid, diamond_report_hits_limit, alignment_mode=2, sensitivity = "ultra-sensitive"):
    # path ist die Assembly FASTA-Datei
    # query_fasta ist die statische FASTA-Datei mit den Abfragesequenzen
    
    #Alignment mode is propably not needed anywhere, but I like args to be consistent with mmseqs
    

    diamond = myUtil.find_executable("diamond")
    # Erstellen der Diamond-Datenbank
    target_db_name = f"{path}.dmnd"
    os.system(f'{diamond} makedb --quiet --in {path} -d {target_db_name} --threads {cores} 1>/dev/null 0>/dev/null')
    
    # Suchergebnisse-Dateien
    output_results_tab = f"{path}.diamond.tab"

    # Durchführen der Diamond-Suche
    #{hit_id}\t{query_id}\t{e_value}\t{score}\t{bias}\t{hsp_start}\t{hsp_end}\t{description}
    print(f'[INFO] {diamond} blastp --quiet --{sensitivity} -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} --outfmt 6 sseqid qseqid evalue bitscore sstart send pident 1>/dev/null 0>/dev/null')
    os.system(f'{diamond} blastp --quiet --{sensitivity} -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} -k {diamond_report_hits_limit} --outfmt 6 sseqid qseqid evalue bitscore sstart send pident 1>/dev/null 0>/dev/null')
    #output format hit query evalue score identity alifrom alito
    return output_results_tab

    
def collect_genomeIDs(report_path, div='___'):
    """
    Collect genome IDs from a TSV file where the hit identifier is in the first column.

    Parameters:
    report_path (str): Path to the TSV file.
    div (str): Divider used when the proteinID includes the genomeID.

    Returns:
    set: A set of unique genome IDs.
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
    
def split_genomeIDs_into_batches(genomeIDs_list, num_batches):
    """
    Splits a list of genomeIDs into nearly equal-sized batches without using math.ceil().
    
    Args:
        genomeIDs_list (list): The list of genomeIDs to be split.
        num_batches (int): The number of batches to create.
    
    Returns:
        list: A list containing `num_batches` sublists, each containing a portion of the genomeIDs.
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

def process_parallel_bulk_parse_batch(args):
    """
    This function processes a batch of genomeIDs in parallel.
    
    Args:
        args: Tuple containing data_queue, genomeID_batch, and other required parameters.
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
   
def process_single_genome(data_queue,genomeID,options,report,score_threshold_diction, csb_patterns_diction,csb_pattern_names,counter):
    """
    Processes a single genomeID, extracts hits, filters, and stores the results.

    Args:
        data_queue: Multiprocessing queue for passing results.
        genomeID (str): The genome ID being processed.
        options (object): Configuration options.
        report (str): DIAMOND BLAST results table.
        score_threshold_diction (dict): Score thresholds.
        csb_patterns_diction (dict): CSB pattern dictionary.
        csb_pattern_names (list): CSB pattern names.
        counter (multiprocessing.Value): Shared counter for progress tracking.
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
        error_message = f"\nError: occurred: {str(e)}"
        traceback_details = traceback.format_exc()
        print(f"\t[WARN] Skipped {faa_file} due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")
        return

    return


def process_writer(queue, options, counter):
    # This routine handles the output of the search and writes it into the database
    # It gets input from multiple workers as the database connection to sqlite is unique
    
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

def submit_batches(protein_batch, cluster_batch, genomeID_set, options):
    # Submit the batches to the database and write to the file

    # Insert into the database
    Database.insert_database_genomeIDs(options.database_directory, genomeID_set)
    Database.insert_database_proteins(options.database_directory, protein_batch)
    Database.insert_database_clusters(options.database_directory, cluster_batch)
    
    #write initial hits fasta
    #write_query_hit_sequence_fasta(options.fasta_initial_hit_directory, protein_batch)
    
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

def parse_bulk_blastreport_genomize(genome_id, filepath, thresholds, cut_score=10):
    """
    Parses a BLAST report for a specific genomeID and extracts protein domain hits.

    Args:
        genome_id (str): The genome ID to filter hits.
        filepath (str): Path to the BLAST report file.
        thresholds (dict): Threshold values (currently unused in the function).
        cut_score (int): Score threshold for filtering hits.

    Returns:
        dict: A dictionary mapping protein IDs to ParseReports.Protein objects.
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
                    print(f"\t[SKIP] Skipped malformed line in {filepath}: {line.strip()} (ValueError: {ve})")
                    continue

    except Exception as e:
        print(f"\n[ERROR] Failed to parse {filepath} for genome {genome_id}")
        print(f"Exception: {str(e)}")
        print(f"Traceback:\n{traceback.format_exc()}")

    return protein_dict
    




def remove_multi_domain_proteins(input_dict):
    #returns only the key:protein pairs that have not more than 1 domain
    return {
        key: protein
        for key, protein in input_dict.items()
        if len(protein.domains.values()) <= 1
    }

def clean_dict_keys_and_protein_ids(input_dict, genomeID):
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


def get_sequence_legth(file_path):
    #get the sequence length per query name without the numbering. If multiple queries have the same
    #name, then take the average
    
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
    
    
def get_sequence_hits_scores(blast_file):
    """
    Generates a dictionary of self-blast scores from a BLAST table file.

    :param blast_file: Path to the BLASTP table file.
    :return: Dictionary with qseqid as keys and the highest self-hit bitscore as values.
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
    
        
def filter_blast_table(output_file, blast_file, evalue_cutoff, score_cutoff, coverage_cutoff, identity_cutoff, bsr_cutoff, sequence_lengths, selfblast_scores, buffer_size=10000):
    """
    Filters a BLASTP table based on given criteria including Blast Score Ratio (BSR),
    and writes the results to a new file using buffered writing to optimize RAM usage.

    :param output_file
    :param blast_file: Path to the BLASTP table file.
    :param evalue_cutoff: Maximum allowable e-value.
    :param score_cutoff: Minimum required bit score.
    :param coverage_cutoff: Minimum required coverage (percentage as a decimal, e.g., 0.8 for 80%).
    :param identity_cutoff: Minimum required percentage identity (pident as a decimal, e.g., 0.9 for 90%).
    :param bsr_cutoff: Minimum required Blast Score Ratio (BSR).
    :param sequence_lengths: Dictionary of sequence lengths for qseqid.
    :param selfblast_scores: Dictionary of self-blast scores for each query (qseqid).
    :param buffer_size: Number of rows to store in memory before writing to disk.
    :return: Path to the filtered output file.
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
                print(f"[WARN] Skipping malformed row: {row} (ValueError: {ve})")
                continue  # Skip invalid rows gracefully

        # Final flush: Write any remaining data in the buffer
        if buffer:
            writer.writerows(buffer)

    return output_file  

    
################################################################################################        
##################################### Selfblast query ##########################################
################################################################################################

def self_blast_query(options):
    #Selfblast the query file
    report = DiamondSearch(options.self_query, options.query_file, options.cores, options.evalue, 100, 100, options.diamond_report_hits_limit) #Selfblast, coverage and identity have to be 100 % or weakly similar domains may occur
    protein_dict = parse_bulk_blastreport_genomize("QUERY",report,{},10) #Selfblast should not have any cutoff score
    
    ParseReports.getProteinSequence(options.self_query,protein_dict) #Get the protein Sequences
    protein_dict = clean_dict_keys_and_protein_ids(protein_dict, "QUERY")
    
    Database.insert_database_genomeIDs(options.database_directory, {"QUERY"})
    Database.insert_database_proteins(options.database_directory, protein_dict)
    #write_query_hit_sequence_fasta(options.fasta_initial_hit_directory, protein_dict)
    
    return report



                        
