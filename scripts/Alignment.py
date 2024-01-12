import os
import re
from Bio import SeqIO
import subprocess


def align_tcs(options,):
    #input files need extension .faa
    #t coffee in path
    #mafft available -> possible to give these as executables? TODO
    
    
    current_dir = os.getcwd()
    os.chdir(options.fasta_output_directory)
    directory_to_clean = options.fasta_output_directory

    # Überprüfen, ob das Verzeichnis existiert
    if os.path.exists(directory_to_clean):
        for filename in os.listdir(directory_to_clean):
            if filename.endswith(".fasta_aln") or filename.endswith(f"tcs_residue_filter{options.tcs_filter}_fasta") or filename.endswith(".dnd"):
                file_path = os.path.join(directory_to_clean, filename)
                os.remove(file_path)
               
    
    infiles = [os.path.join(options.fasta_output_directory, filename) for filename in os.listdir(options.fasta_output_directory) if filename.endswith(".faa")]
    length = len(infiles)
    for index,infile in enumerate(infiles):
        
        name, ext = os.path.splitext(os.path.basename(infile))
        print(f"\tMcoffee Alignment for {name} {index+1} of {length}")

        if ext != ".faa":
            continue

        # Alignment with t-coffee
        # Filter with TCS 4 to improve phylogenetic tree Mprobcons_msa
        # Mpcma_msa Mclustalw_msa Mdialigntx_msa Mpoa_msa Mmuscle_msa Mprobcons_msa Mt_coffee_msa Mclustalw_msa 
        
        output = f"tcs_residue_filter{options.tcs_filter}_fasta"
        mcoffee = f"t_coffee -in={infile} Mmafft_msa Mclustalw_msa -output={output} -maxnseq={options.alnmax_seqs} -maxlen={options.alnmax_len} -case=upper -seqnos=off -outorder=input -run_name=Mcoffee-{name} -multi_core=[templates,jobs,relax,msa] -thread={options.cores} -quiet=stdout 1>/dev/null"
        status = subprocess.run(mcoffee, shell=True).returncode
        if status > 0:
            print(f"\nERROR:\t{mcoffee}\nDied with exit code:\t{status}\n")
            continue

        # Remove 90% gap columns with t-coffee
        infile = f"Mcoffee-{name}.{output}"
        #infile = output
        output = f"{name}.fasta_aln"
        
        mcoffee = f"t_coffee -other_pg seq_reformat -in={infile}  -action +rm_gap {options.gap_col_remove} -output fasta_aln > {output} "
        status = subprocess.run(mcoffee, shell=True).returncode
        if status > 0:
            print(f"\nERROR:\t{mcoffee}\nDied with exit code:\t{status}\n")
            continue

    os.chdir(current_dir)


