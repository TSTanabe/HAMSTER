#!/usr/bin/python

import os
import sys
import csv
import re
import multiprocessing

from . import myUtil

from Bio import SeqIO

def parallel_translation(directory,cores):
    """
    18.5.24
        Args:  
            directory   fasta file containing directory
            
        Uses prodigal to translate all nucleotide fasta files with fna or fna.gz or .fasta ending to faa
        Warning: if the directory path includes parentheses function prodigal is not working
    """
    zipFnaFiles = myUtil.compareFileLists(directory,".fna.gz",".faa.gz") 
    FnaFiles = myUtil.compareFileLists(directory,".fna",".faa")
    fastaFiles = myUtil.getAllFiles(directory,".fasta")
    NucleotideFastaFiles = zipFnaFiles + FnaFiles + fastaFiles
    print(f"Found {len(NucleotideFastaFiles)} assemblies in nucleotide or ambigous format for prodigal")
    
    manager  = multiprocessing.Manager()
    counter = manager.Value('i',0)
    lock = manager.Lock()
    length = len(NucleotideFastaFiles)
    with multiprocessing.Pool(processes=cores) as pool:
        args_list = [(fasta, length ,counter,lock) for fasta in NucleotideFastaFiles]
        pool.map(translate_fasta, args_list)

    return



def translate_fasta(args):
    fasta,length, counter, lock = args

    #unpack if required
    if os.path.splitext(fasta)[-1] == ".gz":
        fasta = myUtil.unpackgz(fasta)

    #Run prodigal
    output = os.path.splitext(fasta)[0]
    faa = output + ".faa"
    
    prodigal = myUtil.find_executable("prodigal")
    
    string = f"{prodigal} -a {faa} -i {fasta} >/dev/null 2>&1"
    try:
        os.system(string)
    except Exception as e:
        print(f"\tWARNING: Could not translate {fasta} - {e}")
        
    with lock:
        counter.value += 1
        print(f"\rProcessing assembly {counter.value} of {length}", end ='',flush=True)        
    return




############################################################################
############### Parallel Transcription #####################################
############################################################################


def parallel_transcription(directory,cores):
    """
    8.10.22
        Args:  
            directory   fasta file containing directory
            
        Transcribe for all faa files gff3 files
        Secure a packed and unpacked version is present
        Unlink unpacked versions afterwards
        Warning: if the directory path includes parentheses function prodigal is not working
    """
            
    

    gzfaaFiles = myUtil.getAllFiles(directory,".faa.gz")
    gzgffFiles = myUtil.getAllFiles(directory,".gff.gz")
    print(f"Found {len(gzfaaFiles)} zipped faa files")
    print(f"Found {len(gzgffFiles)} zipped gff files")

    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(unpacker,gzgffFiles)
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(unpacker,gzfaaFiles)   
    
    faaFiles = myUtil.getAllFiles(directory,".faa")
    gffFiles = myUtil.getAllFiles(directory,".gff")   
    print(f"Found {len(gffFiles)} gff files")
    print(f"Found {len(faaFiles)} faa files")    
        	
    FaaFiles = myUtil.compareFileLists(directory,".faa",".gff")
    print(f"Found {len(FaaFiles)} protein fasta files without gff")
    
    manager = multiprocessing.Manager()
    counter = manager.Value('i',0)
    lock = manager.Lock()
    length = len(FaaFiles)
    
    with multiprocessing.Pool(processes=cores) as pool:
        args_list = [(fasta, length, counter, lock) for fasta in FaaFiles]
        pool.map(transcripe_fasta, args_list)
    
    print("\nFinished faa and gff file preparation")
    return


def transcripe_fasta(args):
    fasta, length, counter, lock = args
    gff = ""
            
    if check_prodigal_format(fasta):        
        gff = prodigalFaaToGff(fasta)
    
        with lock:
            counter.value += 1
            print(f"\rProcessing file {counter.value} of {length}", end='', flush=True)
  
    return

def check_prodigal_format(File):
    
    with open(File, "r") as reader:
        for line in reader.readlines():
            if line[0] == ">":
                line = line[1:]
                ar = line.split("#")
                if len(ar) == 5:
                    return 1
                else:
                    return 0
    return Gff

def prodigalFaaToGff(filepath):
#Translates faa files from prodigal to normal gff3 formated files. returns the gff3 file name    
    dir_path = os.path.dirname(filepath)
    # Extract the filename with extensions
    filename_with_ext = os.path.basename(filepath)
    
    if filename_with_ext.endswith('.gz'):
        filename_with_ext = os.path.splitext(filename_with_ext)[0]

    # Remove the .faa extension if present
    if filename_with_ext.endswith('.faa'):
        filename_without_ext = os.path.splitext(filename_with_ext)[0]
    else:
        filename_without_ext = filename_with_ext
    
    Gff = dir_path+"/"+filename_without_ext+'.gff'
    
    writer = open(Gff,"w")
    
    with open(filepath, "r") as reader:
        genomeID = myUtil.getGenomeID(filepath)
        for line in reader.readlines():
            if line[0] == ">":
                try:
                    line = line[1:]
                    ar = line.split("#")
                    #print(ar)
                    contig = re.split("\_{1}\d+\W+$",ar[0])
                    #print(contig)
                    strand = '+' if ar[3] == ' 1 ' else '-'
                    writer.write(contig[0]+"\tprodigal\tcds\t"+ar[1]+"\t"+ar[2]+"\t0.0\t"+strand+"\t0\tID=cds-"+ar[0]+";Genome="+genomeID+"\n")
                except Exception as e:
                    print(f"Error: Missformated header\n {line} - {e}")
    writer.close()
    return Gff


def packer(file):
    myUtil.packgz(file)
def unpacker(file):
    myUtil.unpackgz(file)

######################################################################
################# Glob files creation routines #######################
######################################################################

################# Concatenate subroutines     
def create_glob_file(options):
    # Check if a glob fasta file and deconcatenated files were already provided
    if not options.glob_faa is None and os.path.isfile(options.glob_faa) and os.path.isdir(options.fasta_file_directory):
        return

    # Define the output file path for the glob.faa file
    options.glob_faa = os.path.join(options.result_files_directory, "glob.faa")
    
    with open(options.glob_faa, 'w') as outfile:
        for genomeID in options.queued_genomes:
            faa_file = options.faa_files[genomeID]  # Get the FASTA file for the current genomeID

            # Open each faa_file explicitly and close it after processing
            infile = open(faa_file, 'r')
            try:
                sequence = ""
                header = None

                for line in infile:
                    line = line.strip()
                    if line.startswith(">"):
                        # Write the previous sequence if it exists
                        if header and sequence:
                            outfile.write(f"{header}\n{sequence}\n")
                        
                        # Start a new record
                        original_id = line[1:].split()[0]  # Extract original ID up to the first whitespace
                        header = f">{genomeID}___{original_id}"  # Create the modified header
                        sequence = ""  # Reset sequence for the new record
                    else:
                        sequence += line  # Append sequence lines
                
                # Write the last record if there was one
                if header and sequence:
                    outfile.write(f"{header}\n{sequence}\n")
            finally:
                infile.close()  # Ensure the file is closed after reading



def create_selfquery_file(options):
    # Construct the paths for the output files
    options.self_query = os.path.join(options.result_files_directory, "self_blast.faa")
    options.self_seqs = os.path.join(options.result_files_directory, "self_seqs.faa")
    
    # Dictionary to keep track of the number of times each original_id is encountered
    id_count = {}
    
    # Open both output files for writing
    with open(options.self_query, 'w') as outfile_query, open(options.self_seqs, 'w') as outfile_seqs:
        
        # Open the input query file and process it line by line
        with open(options.query_file, 'r') as infile:
            sequence = ""
            header = None

            for line in infile:
                line = line.strip()
                if line.startswith(">"):
                    # If a header is found, process the previous record
                    if header and sequence:
                        # Extract the original ID (up to the first whitespace)
                        original_id = header[1:].split()[0]
                        
                        # Increment the counter for this ID
                        id_count[original_id] = id_count.get(original_id, 0) + 1
                        count = id_count[original_id]

                        # Create modified headers for both output files
                        query_header = f">QUERY___{original_id}___{count}"
                        seqs_header = f">{original_id}___{count}"

                        # Write to self_blast.faa with the QUERY___ prefix
                        outfile_query.write(f"{query_header}\n{sequence}\n")
                        # Write to self_seqs.faa without the QUERY___ prefix
                        outfile_seqs.write(f"{seqs_header}\n{sequence}\n")
                    
                    # Start a new record
                    header = line
                    sequence = ""  # Reset sequence for the new record
                else:
                    sequence += line  # Append sequence lines

            # Write the last record if there was one
            if header and sequence:
                original_id = header[1:].split()[0]
                id_count[original_id] = id_count.get(original_id, 0) + 1
                count = id_count[original_id]

                query_header = f">QUERY___{original_id}___{count}"
                seqs_header = f">{original_id}___{count}"

                outfile_query.write(f"{query_header}\n{sequence}\n")
                outfile_seqs.write(f"{seqs_header}\n{sequence}\n")


############### Deconcat subroutines

def deconcat(options):

    options.fasta_file_directory = options.result_files_directory+"/deconcatenated_fasta_files"
    deconcatenate_faa(options.glob_faa,options.fasta_file_directory)
    deconcatenate_gff(options.glob_gff,options.fasta_file_directory)
    
def deconcatenate_faa(input_filepath, output_directory):
    """
    Deconcatenates a protein FASTA file into genome-wise FASTA files,
    removing everything in the header up to the first '___'.
    
    Args:
        input_filepath (str): Path to the concatenated FASTA file.
        output_directory (str): Directory to save the genome-wise FASTA files.
    """
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    current_genome_id = None
    output_handle = None

    def close_current_handle():
        if output_handle:
            output_handle.close()

    # Read the concatenated FASTA file
    with open(input_filepath, "r") as infile:
        for line in infile:
            if line.startswith(">"):
                #Extract the genomeID from the identifier line (header)
                identifier = line[1:].strip() # Remove '>' and any surrounding spaces
                genome_id = identifier.split('___')[0]
                
                if genome_id != current_genome_id:
                    close_current_handle()
                    current_genome_id = genome_id
                    output_filepath = os.path.join(output_directory, f"{genome_id}.faa")
                    output_handle = open (output_filepath, "a") 
                output_handle.write(line)
            else:
                output_handle.write(line)
    close_current_handle()
    print(f"Deconcatenation of faa complete. Files saved to {output_directory}")

  
 

def deconcatenate_gff(input_filepath, output_directory):
    """
    Deconcatenates a protein GFF file into genome-wise GFF files.
    
    Args:
        input_filepath (str): Path to the concatenated GFF file.
        output_directory (str): Directory to save the genome-wise GFF files.
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    current_genome_id = None
    output_handle = None

    # Pre-compile the regular expression for better performance
    genome_id_regex = re.compile(r'ID=(.*?)___')

    # Open the input file and process it line by line
    with open(input_filepath, "r") as reader:
        for row in reader:
            # Search for genome ID in the current row
            match = genome_id_regex.search(row)
            if match:
                genomeID = match.group(1)
                # If the genomeID changes, switch to a new file
                if genomeID != current_genome_id:
                    # Close the previous file handle if open
                    if output_handle:
                        output_handle.close()

                    # Update the current genome ID and open a new output file for writing
                    current_genome_id = genomeID
                    output_filepath = os.path.join(output_directory, f"{genomeID}.gff")
                    output_handle = open(output_filepath, "a")

                # Write the current row to the appropriate genome file
                output_handle.write(row)

    # Close the last output handle if open
    if output_handle:
        output_handle.close()

    print(f"Deconcatenation of GFF complete. Files saved to {output_directory}")
