#!/usr/bin/python
import os
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
    print("\nFinished searching")
    print(f"Saved hits to database {options.database_directory}")
    print(f"Saved individual hit lists to {options.fasta_initial_hit_directory}")
    return    

def process_parallel_search(args_tuple):

    queue,genomeID,options,counter,score_threshold_diction, csb_patterns_diction,csb_pattern_names = args_tuple
    counter.value += 1
    print(f"Searching assembly ({counter.value}/{len(options.queued_genomes)})", end="\r")
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
        print(f"\tWARNING: Skipped {faa_file} due to an error - {error_message}")
        print(f"\tTraceback details:\n{traceback_details}")
        return

    return

    
    
    
     

    
################################################################################################        
############################################ Glob search #######################################
################################################################################################

def initial_glob_search(options):
    print(f"Starting diamond blastp search")
    blast_results_table = options.glob_table
    if not blast_results_table: # If a blast table is not provided make the blastp
        blast_results_table = DiamondSearch(options.glob_faa, options.query_file, options.cores, options.evalue, options.searchcoverage, options.minseqid, options.diamond_report_hits_limit, options.alignment_mode) #Diamond search the target fasta file
    
    genomeIDs_set = collect_genomeIDs(blast_results_table) #returns a set of all genomeIDs
    print(f"Found hits in {len(genomeIDs_set)} genomes")
    Database.insert_database_genomeIDs(options.database_directory, genomeIDs_set) # Insert genomeIDs into database
    
    genomeID_batches = split_genomeIDs_into_batches(list(genomeIDs_set), options.cores-1) # batch sets of genomeIDs for each worker
    
    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value('i', 0)
    csb_patterns_diction,csb_pattern_names = Csb_finder.makePatternDict(options.patterns_file)
    score_threshold_diction = {} #empty because all have to be parsed
    
    with Pool(processes = options.cores*4, initializer=init_worker, 
              initargs=(blast_results_table, score_threshold_diction, csb_patterns_diction, csb_pattern_names)) as pool:
    
        p_writer = pool.apply_async(process_writer, (data_queue, options, counter))
        
        pool.map(process_parallel_bulk_parse_batch, [(data_queue, batch, options, counter) for batch in genomeID_batches])

        # Signal the end of data to the process_writer
        for _ in range(options.cores):
            data_queue.put(None)
        
        p_writer.get()
    print("Finished parsing reports")
    return


# Initializer function to set up global variables for each worker
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

def DiamondSearch(path, query_fasta, cores, evalue, coverage, minseqid, diamond_report_hits_limit, alignment_mode=2):
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
    print(f'{diamond} blastp --quiet --ultra-sensitive -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} --outfmt 6 sseqid qseqid evalue bitscore sstart send pident 1>/dev/null 0>/dev/null')
    os.system(f'{diamond} blastp --quiet --ultra-sensitive -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} -k {diamond_report_hits_limit} --outfmt 6 sseqid qseqid evalue bitscore sstart send pident 1>/dev/null 0>/dev/null')
    #output format hit query evalue score identity alifrom alito
    return output_results_tab

def MMseqsSearch(path, query, cores, evalue, coverage, minseqid, alignment_mode=2):
    #query is the static fasta file with the query sequences
    #path is the assembly fasta file
    
    mmseqs = Alignment.find_executable("mmseqs")
    
    target_db_name=f"{path}.targetDB"
    os.system(f'{mmseqs} createdb {path} {target_db_name} 1>/dev/null')
    tmp = f"{path}.tmp"
    
    query_db_name=f"{query}.queryDB" #can be done during argument parsing if possible
    os.system(f'{mmseqs} createdb {query} {query_db_name} 1>/dev/null')
    
    output_results=f"{path}.alndb"
    output_results_tab=f"{path}.tab"

    os.system(f'{mmseqs} search {query_db_name} {target_db_name} {output_results} {tmp} --threads {cores} --alignment-mode {alignment_mode} --min-seq-id {minseqid} -e {evalue} -c {coverage} 1>/dev/null') #attention cores is used here --min-seq-id {coverage}
    os.system(f'{mmseqs} convertalis {query_db_name} {target_db_name} {output_results} {output_results_tab} --threads {cores} --format-output "target,query,evalue,bits,tstart,tend,pident" 1>/dev/null')
    
    #os.system(f'mmseqs rmdb {query_db_name}')
    #os.system(f'mmseqs rmdb {target_db_name}')
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
    
        
######## Parsing routines ############
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
        process_parallel_bulk_parse((data_queue, genomeID, options, diamond_blast_results_table, score_threshold_diction, csb_patterns_diction, csb_pattern_names, counter))

   
def process_parallel_bulk_parse(args_tuple):
    
    data_queue,genomeID,options,report,score_threshold_diction, csb_patterns_diction,csb_pattern_names,counter = args_tuple
    #counter.value +=1
    #print(f"Processed {counter.value} genomes by worker", end="\r")#
    protein_dict = dict()
    cluster_dict = dict()
    try:
        faa_file = myUtil.unpackgz(options.faa_files[genomeID])
        gff_file = myUtil.unpackgz(options.gff_files[genomeID])
        
        #Parse the hits
        protein_dict = parse_bulk_blastreport_genomize(genomeID,report,score_threshold_diction,options.thrs_score)
        
        #if options sagt, dass es einzelgenome sind, dann die genomID aus den keys und den protein identifiern entfernen
        if not options.glob_gff:
            protein_dict = clean_dict_keys_and_protein_ids(protein_dict, genomeID)
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
        print(f"\tWARNING: Skipped {faa_file} due to an error - {error_message}")
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
            print(f"Processed {counter.value} genomes ", end="\r")#
        
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
    writeQueryHitsSequenceFasta(options.fasta_initial_hit_directory, protein_batch)
    
    # Append to the gene clusters file
    with open(options.gene_clusters_file, "a") as file:
        for clusterID, cluster in cluster_batch.items():
            domains = cluster.get_domains()
            file.write(clusterID + '\t' + '\t'.join(domains) + '\n')


def writeQueryHitsSequenceFasta(directory, protein_dict):
    # Sort proteins by domains
    #Fetch operation from the database would be more stable but possibly slower
    sorted_proteins = sorted(protein_dict.values(), key=lambda x: x.get_domains())
    
    current_file_path = None
    current_file = None

    try:
        for protein in sorted_proteins:
            genomeID = protein.genomeID
            proteinID = f"{genomeID}-{protein.proteinID}" #same concatenation is done for the database, do not alter
            domains = protein.get_domains()
            file_path = f"{directory}/Query_{domains}.hit_list"
            
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

def parse_bulk_blastreport_genomize(genomeID,Filepath,Thresholds,cut_score=10):
    #
    # Parse a hit table with sseqid qseqid evalue bitscore sstart send pident
    # First use grep to get all lines with the genomeID. Requisite: GENOMEID has to be part of the Hit identifier
    #
    
    
    protein_dict = {}

    result = subprocess.run(['grep', genomeID, Filepath], stdout=subprocess.PIPE, text=True)
    
    lines = result.stdout.splitlines()  # Split output into lines

    for line in lines:
        columns = line.split('\t')  # Assuming columns are tab-separated
        if columns:
            try:
                #{hit_id}\t{query_id}\t{e_value}\t{score}\t{bias}\t{hsp_start}\t{hsp_end}\t{description}
                #the hit identifier in the fasta file must have a genomeID___proteinID$ format, only the first ___ will be recognized
                
                #Bei Query___proteinID$
                
                hit_proteinID = columns[0]
                query = columns[1]
                hit_bitscore = int(float(columns[3]))
                hsp_start = int(float(columns[4]))
                hsp_end = int(float(columns[5]))
                
                if hit_proteinID in protein_dict:
                    protein = protein_dict[hit_proteinID]
                    protein.add_domain(query,hsp_start,hsp_end,hit_bitscore)
                else:
                    protein_dict[hit_proteinID] = ParseReports.Protein(hit_proteinID,query,hsp_start,hsp_end,hit_bitscore,genomeID)
            except Exception as e:
                error_message = f"\nError occurred: {str(e)}"
                traceback_details = traceback.format_exc()
                print(f"\tWARNING: Skipped {Filepath} due to an error - {error_message}")
                print(f"\tTraceback details:\n{traceback_details}")
                continue
                
    return protein_dict


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
##################################### Selfblast query ##########################################
################################################################################################

def self_blast_query(options):
    #Selfblast the query file
    report = DiamondSearch(options.self_query, options.query_file, options.cores, options.evalue, 100, 100, options.diamond_report_hits_limit) #Selfblast, coverage and identity have to be 100 % or weakly similar domains may occur
    protein_dict = parse_bulk_blastreport_genomize("QUERY",report,{},10) #Selfblast should not have any cutoff score
    print(protein_dict.keys())
    ParseReports.getProteinSequence(options.self_query,protein_dict) #Get the protein Sequences
    protein_dict = clean_dict_keys_and_protein_ids(protein_dict, "QUERY")
    print(protein_dict.keys())
    Database.insert_database_genomeIDs(options.database_directory, {"QUERY"})
    Database.insert_database_proteins(options.database_directory, protein_dict)
    writeQueryHitsSequenceFasta(options.fasta_initial_hit_directory, protein_dict)
    
    
