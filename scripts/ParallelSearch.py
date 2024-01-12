#!/usr/bin/python
from . import myUtil
from . import Search
from . import ParseReports
from . import Csb_finder
from . import Database
import multiprocessing
from multiprocessing import Process, Manager, Pool, Semaphore
#import re
#import os



def process_parallel_search(args_tuple):
    
    queue,genomeID,options,counter = args_tuple
    counter.value += 1
    print(f"Processing assembly {counter.value} of {len(options.queued_genomes)}")
    try:
        query = options.query_file
        faa_file = myUtil.unpackgz(options.faa_files[genomeID])
        gff_file = myUtil.unpackgz(options.gff_files[genomeID])
        search_results = Search.MMseqsSearch(faa_file,options.query_file,options) #TODO könnte probleme geben wenn mehrere mmseqs auf die selbe query zugreifen wollen
        protein_dict = ParseReports.parseBlastreporttab(search_results,dict()) # returns the hits from a single
        vicinity_genes_dict = ParseReports.parseGFFfile(gff_file,protein_dict,options.gene_vicinity_range) # adds the upstream and downstream genes to a separate dict and gff features to protein_dict
        ParseReports.getProteinSequence(faa_file,protein_dict)
        ParseReports.getProteinSequence(faa_file,vicinity_genes_dict)
        
        combined_protein_diction = dict()
        combined_protein_diction.update(vicinity_genes_dict)
        combined_protein_diction.update(protein_dict)
        cluster_diction = Csb_finder.find_syntenicblocks(genomeID,combined_protein_diction,options.distance)
        queue.put((genomeID,protein_dict,vicinity_genes_dict,combined_protein_diction,cluster_diction))
    except Exception as e:
        error_message = f"\nError occurred: {str(e)}"
        print(f"\t WARNING: Skipped {faa_file} due to an error - {error_message}")
        myUtil.unlink(faa_file)
        myUtil.unlink(gff_file)
        return

    
    myUtil.unlink(faa_file)
    myUtil.unlink(gff_file)
    myUtil.unlink_mmseq_target_tmps(faa_file)
    
    
     
def process_writer(queue,options):
    
    while True:
        tup = queue.get()
        if tup is None :
            break

        genomeID,protein_dict,vicinity_genes_dict,combined_protein_diction,cluster_diction = tup
        ParseReports.writeProteinSequenceFasta(vicinity_genes_dict,genomeID,options.divergent_output_file)     # write the vicinity genes down for later clustering
        ParseReports.writeQueryHitsSequenceFasta(protein_dict,genomeID,options.fasta_initial_hit_directory)
        
        Database.insert_database_genomeID(options.database_directory,genomeID)
        Database.insert_database_protein(options.database_directory,genomeID,combined_protein_diction)
        Database.insert_database_cluster(options.database_directory,genomeID,cluster_diction)
        #this might be slow because the con is made three times with the database
    return
        
def multi_search_process(options):
    
    
    
    query = options.query_file
    query_db_name=f"{query}.queryDB" #can be done during argument parsing if possible
    myUtil.command(f'mmseqs createdb {query} {query_db_name} 1>/dev/null')
    
    manager = Manager()
    data_queue = manager.Queue()
    counter = manager.Value('i', 0)
    
    
    #p2 = multiprocessing.Process(target=process_writer, args=(data_queue,options,))
    #p2.start()
    
    process_instances = []
    
    with Pool(processes = options.cores) as pool:
    
        p_writer = pool.apply_async(process_writer, (data_queue, options))

        pool.map(process_parallel_search, [(data_queue, genomeID, options, counter) for genomeID in options.queued_genomes])

        # Signal the end of data to the process_writer
        for _ in range(options.cores):
            data_queue.put(None)
        print("Waiting for writer")
        p_writer.get()
    print("Finished searching")
    return    
        
     
     
     
     
     
     
        
def initial_search(options):
    
    query = options.query_file
    query_db_name=f"{query}.queryDB" #can be done during argument parsing if possible
    myUtil.command(f'mmseqs createdb {query} {query_db_name} 1>/dev/null')
    
    for index,genomeID in enumerate(options.queued_genomes):
        now = datetime.now()
        print(f"{now} Processing assembly {index+1} of {len(options.queued_genomes)}") #print control sequence output
        #Prepare names
        try:
            faa_file = myUtil.unpackgz(options.faa_files[genomeID])
            gff_file = myUtil.unpackgz(options.gff_files[genomeID]) 
            search_results = Search.MMseqsSearch(faa_file,options.query_file,options)
            protein_dict = ParseReports.parseBlastreporttab(search_results,dict()) # returns the hits from a single genome TODO tmp debug TODO dict shall be a Threshold dictionary
            vicinity_genes_dict = ParseReports.parseGFFfile(gff_file,protein_dict,options.gene_vicinity_range) # adds the upstream and downstream genes to a separate dict and gff features to protein_dict
            ParseReports.getProteinSequence(faa_file,protein_dict)
            ParseReports.getProteinSequence(faa_file,vicinity_genes_dict)
            ParseReports.writeProteinSequenceFasta(vicinity_genes_dict,genomeID,options.divergent_output_file)     # write the vicinity genes down for later clustering
            ParseReports.writeQueryHitsSequenceFasta(protein_dict,genomeID,options.fasta_initial_hit_directory)
            
            #find the syntenic regions and insert to database
            combined_protein_diction = {**vicinity_genes_dict, **protein_dict}
            cluster_diction = Csb_finder.find_syntenicblocks(genomeID,combined_protein_diction,options.distance)
            Database.insert_database_genomeID(options.database_directory,genomeID)
            Database.insert_database_protein(options.database_directory,genomeID,combined_protein_diction)
            Database.insert_database_cluster(options.database_directory,genomeID,cluster_diction)
        except Exception as e:
            error_message = f"\nError occurred: {str(e)}"
            print(f"\t WARNING: Skipped {faa_file} due to an error - {error_message}")
            continue
        myUtil.unlink(faa_file)
        myUtil.unlink(gff_file)
        #myUtil.unlink_mmseq_target_tmps(faa_file)
    #myUtil.unlink_mmseq_query_tmps(options.query_file)
