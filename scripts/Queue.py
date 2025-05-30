#!/usr/bin/python

import os
from . import Database
from . import myUtil

     

def queue_files(options):
    """
    Args:
        directory   directory for the result files to be stored
        finished    set with genome identifiers already processed
        options     current options object
    Operation:
        collect all zip protein fasta files and unzipped protein fasta files
        collect all corresponding gff files and check for existence of this file
        if both files present
        get the genome identifiers
    """
    print("\n[INFO] Filling the queue with faa files --", end="\r")
    genomeID_queue = set()
    faa_files = {}
    gff_files = {}

    pairs = find_faa_gff_pairs(options.fasta_file_directory)
    
    for faa_file,gff_file in pairs:
        genomeID = myUtil.getGenomeID(faa_file)
        genomeID_queue.add(genomeID)
        faa_files[genomeID] = faa_file
        gff_files[genomeID] = gff_file
        
    # compare two sets
    find_missing_genomes(genomeID_queue, options.fasta_file_directory)
    
    options.queued_genomes = genomeID_queue
    options.faa_files = faa_files
    options.gff_files = gff_files
    print("[INFO] Filling the queue with faa files -- ok")
    print(f"[INFO] Queued {len(options.queued_genomes)} faa/gff pairs")
    return
    

def compare_with_existing_database(options,genomeIDs):
    
    #compare with database leaving out all genomes already searched previously
    if options.redo_search: #redo the search for all
        options.finished_genomes = {}
        for genomeID in genomeIDs:
            options.finished_genomes[genomeID] = 1 #required for extending the protein dict in csb assignment with all old results
    else: #leave out all genomeIDs in the database
        genomeIDs = Database.fetch_genomeIDs(options.database_directory)
        for genomeID in genomeIDs:
            if genomeID in options.faa_files.keys():
                print(f"[INFO] Found assembly {genomeID} in database leaving out {myUtil.getFileName(options.faa_files[genomeID])}")
                del options.faa_files[genomeID]
                del options.gff_files[genomeID]
                options.queued_genomes.remove(genomeID)
    
    print(f"[INFO] Queued {len(options.queued_genomes)} for processing")
    if len(options.queued_genomes) == 0:
        print("[ERROR] There were 0 genomes queued. Use -redo_search option if genomes are already present in database")
    
    return
    

def find_faa_gff_pairs(directory):
    """
    Find pairs of files with the same name but different extensions (.faa/.faa.gz and .gff/.gff.gz)
    in the given directory and its subdirectories.

    Args:
        directory (str): The directory to search for file pairs.

    Returns:
        list: A list of tuples, each containing the paths to a paired .faa and .gff file.
    """
    # Dictionary to store files with the same basename
    files_dict = {}

    # Traverse the directory and its subdirectories
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            
            # Check for .faa or .faa.gz files
            if file.endswith('.faa'): #or file.endswith('.faa.gz'):
                basename = file.replace('.faa', '').replace('.gz', '')
                if basename not in files_dict:
                    files_dict[basename] = {}
                files_dict[basename]['faa'] = file_path
            
            # Check for .gff or .gff.gz files
            elif file.endswith('.gff'): # or file.endswith('.gff.gz'):
                basename = file.replace('.gff', '').replace('.gz', '')
                if basename not in files_dict:
                    files_dict[basename] = {}
                files_dict[basename]['gff'] = file_path

    # Find and store pairs of .faa and .gff files
    pairs = []
    for basename, file_paths in files_dict.items():
        if 'faa' in file_paths and 'gff' in file_paths:
            pairs.append((file_paths['faa'], file_paths['gff']))
    return pairs



def find_missing_genomes(genomeIDs, faa_file_directory):
    """Find .faa and .faa.gz files in the directory whose genome IDs are not in the provided list."""
    
    def list_faa_files(directory):
        """List all .faa and .faa.gz files in the directory."""
        return [f for f in os.listdir(directory) if f.endswith('.faa') or f.endswith('.faa.gz')]

    def extract_genomeID_from_faa(filename):
        """Extract genome ID from filename using myUtil.get_genomeID."""
        return myUtil.getGenomeID(filename)

    missing_files = []
    all_faa_files = list_faa_files(faa_file_directory)
    
    for faa_file in all_faa_files:
        genomeID = extract_genomeID_from_faa(faa_file)
        if genomeID not in genomeIDs:
            missing_files.append(faa_file)
    
    return missing_files



#TODO Wird die routine gebraucht?
def query_names(options,query_file):
    
    with open(query_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                header = line.strip()
                options.query_names.append(header)








