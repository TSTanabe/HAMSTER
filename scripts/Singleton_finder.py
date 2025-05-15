#!/usr/bin/python

import os
from . import Csb_proteins
from . import myUtil

    
def extract_protein_ids_from_fasta(fasta_path):
    protein_ids = set()
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                full_id = line[1:].split()[0]
                parts = full_id.split('___')
                protein_id = parts[1] if len(parts) > 1 else full_id
                protein_ids.add(protein_id)
    return protein_ids
    
def extract_domain_names_from_directory(directory_path):
    """
    Extrahiert eindeutige Domänennamen aus Dateinamen in einem Verzeichnis.
    Entfernt Dateiendung und alles bis einschließlich dem ersten '_'.
    """
    domain_names = set()
    
    for filename in os.listdir(directory_path):
        if os.path.isfile(os.path.join(directory_path, filename)):
            name_no_ext = os.path.splitext(filename)[0]  # z.B. "Query_ABC123"
            parts = name_no_ext.split('_', 1)  # nur am ersten "_" splitten
            domain = parts[1] if len(parts) > 1 else name_no_ext
            domain_names.add(domain)
    
    return domain_names
    
def get_min_bitscore_for_query(report_path, query_id, blast_score_ratio = 0.9):
    """
    Liest einen DIAMOND BLAST Report ein und gibt den kleinsten Bitscore für eine gegebene Query-ID zurück.
    
    Args:
        report_path (str): Pfad zur Diamond-Tabelle (.tab).
        query_id (str): Die qseqid (Query-ID), nach der gesucht werden soll.
    
    Returns:
        float: Der kleinste Bitscore für die Query, oder None falls nicht gefunden.
    """
    bitscores = []
    min_cutoff = 100
    max_cutoff = 100
    with open(report_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            hit_id, qid, evalue, bitscore = parts[:4]
            hit_id = hit_id.split('___')[1]
            if hit_id == query_id and qid == query_id:
                try:
                    bitscores.append(float(bitscore))
                except ValueError:
                    continue  # ungültiger bitscore

    if len(bitscores) == 1:
        max_cutoff = bitscores[0]*(2-blast_score_ratio)
        min_cutoff = bitscores[0]*blast_score_ratio #if only a single sequence is given than 0.9 blast score ratio is accepted
    else:
        max_cutoff = max(bitscores)
        min_cutoff = min(bitscores)
    return min_cutoff, max_cutoff



#### Main routine of this module
def singleton_reference_sequences(options):

    # Load cache if available
    domain_score_limits = myUtil.load_cache(options, "sng_domain_score_limits.pkl")
    singleton_reference_seqs_dict = myUtil.load_cache(options, "sng_reference_seqs_dict.pkl")
    
    if domain_score_limits and singleton_reference_seqs_dict:
        print("Loaded existing reference sequences for genes without conserved genomic context")
        return domain_score_limits,singleton_reference_seqs_dict
        

    # If cache is not available calculate score limits and define reference sequences
    
    # Get the domains as a set that are present in the query file but not in training seq data
    query_names = extract_protein_ids_from_fasta(options.self_query)
    training_set_domains = extract_domain_names_from_directory(options.fasta_output_directory)
    singletons = query_names - training_set_domains
    print(f"Proteins without recognized genomic context {singletons}")
    
    # Get the bitscores from the Query
    domain_score_limits = {}
    output_results_tab = f"{options.self_query}.diamond.tab"
    for singleton in singletons:
        min_cutoff,max_cutoff = get_min_bitscore_for_query(output_results_tab,singleton) #TODO the last argument, blastscore ratio, is yet to be included to hits routine
        domain_score_limits[singleton] = {
            "lower_limit": min_cutoff,
            "upper_limit": max_cutoff
        }
    
    myUtil.save_cache(options, "sng_domain_score_limits.pkl", domain_score_limits)
    
    # fetch proteins with the singleton domain and above thrs score to fasta
    singleton_reference_seqs_dict = Csb_proteins.fetch_protein_ids_parallel(options.database_directory,domain_score_limits, options.cores)
    singleton_reference_seqs_dict = {f"sng0_{k}": v for k, v in singleton_reference_seqs_dict.items()} # Add prefix
    myUtil.save_cache(options, "sng_reference_seqs_dict.pkl", singleton_reference_seqs_dict)
    Csb_proteins.fetch_seqs_to_fasta_parallel(options.database_directory, singleton_reference_seqs_dict, options.fasta_output_directory, 5, options.max_seqs, options.cores)
    
    return domain_score_limits, singleton_reference_seqs_dict








    
