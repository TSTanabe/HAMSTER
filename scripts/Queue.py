#!/usr/bin/python

import os
from typing import List, Tuple, Set

from . import Database
from . import myUtil

     
logger = myUtil.logger

def queue_files(options) -> None:
    """
    Fills the options object with genome IDs, .faa and .gff file mappings.
    
    Args:
        options: options object with at least .fasta_file_directory, will be filled with:
            .queued_genomes (set[str])
            .faa_files (dict[str, str])
            .gff_files (dict[str, str])
        
    Operation:
        - Collect all zipped/unzipped protein fasta files and corresponding gff files.
        - Queue only if both files present, by genome identifier.
        
    Output Example:
        options.queued_genomes = {'GCF_000001405.39', ...}
        options.faa_files = {'GCF_000001405.39': '/dir/xxx.faa', ...}
        options.gff_files = {'GCF_000001405.39': '/dir/xxx.gff', ...}
    """
    
    logger.info("Filling the queue with faa files to be processed")
    genomeID_queue = set()
    faa_files = {}
    gff_files = {}

    pairs = find_faa_gff_pairs(options.fasta_file_directory)
    
    for faa_file,gff_file in pairs:
        genomeID = myUtil.getGenomeID(faa_file)
        genomeID_queue.add(genomeID)
        faa_files[genomeID] = faa_file
        gff_files[genomeID] = gff_file
        
    # compare two sets (find missing)
    find_missing_genomes(genomeID_queue, options.fasta_file_directory)
    
    options.queued_genomes = genomeID_queue
    options.faa_files = faa_files
    options.gff_files = gff_files
    
    logger.info(f"Queued {len(options.queued_genomes)} faa/gff pairs")
    
    return
    

def find_faa_gff_pairs(directory: str) -> List[Tuple[str, str]]:
    """
    Find pairs of files with the same name but different extensions (.faa/.faa.gz and .gff/.gff.gz)
    in the given directory and its subdirectories.

    Args:
        directory (str): The directory to search for file pairs.

    Returns:
        list of tuple: Each containing the paths to a paired .faa and .gff file.

    Output Example:
        [('/path/xxx.faa', '/path/xxx.gff'), ...]
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



def find_missing_genomes(genomeIDs: Set[str], faa_file_directory: str) -> List[str]:
    """
    Find .faa files in the directory whose genome IDs are not in the provided list.

    Args:
        genomeIDs (set): Set of genome IDs
        faa_file_directory (str): Directory to search

    Returns:
        List of missing .faa file names (not present in genomeIDs)

    Output Example:
        ['GCF_000001405.39.faa', ...]
    """
        
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







