import os
import sys
import shutil
import subprocess
from . import myUtil



def align_tcs(options,):
    #input files need extension .faa
    #t coffee in path
    
    
    
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
            sys.exit("Terminating the process")

        # Remove 90% gap columns with t-coffee
        infile = f"Mcoffee-{name}.{output}"
        #infile = output
        output = f"{name}.fasta_aln"
        
        mcoffee = f"t_coffee -other_pg seq_reformat -in={infile}  -action +rm_gap {options.gap_col_remove} -output fasta_aln > {output} "
        status = subprocess.run(mcoffee, shell=True).returncode
        if status > 0:
            print(f"\nERROR:\t{mcoffee}\nDied with exit code:\t{status}\n")
            sys.exit("Terminating the process")

    os.chdir(current_dir)

####################################################
##########   MAFFT Alignments routines    ##########
####################################################

#concats .fasta_aln files in a directory. if none are found
#it aligns with mafft and removes gaps with trimal at 95 % 
            
            
def initial_alignments(options):

    directory_to_clean = options.fasta_output_directory

    # Clean up preexisting alignments
    if os.path.exists(directory_to_clean):
        for filename in os.listdir(directory_to_clean):
            if filename.endswith(".fasta_aln"):
                file_path = os.path.join(directory_to_clean, filename)
                os.remove(file_path)
                
                
    filepaths = get_alignments(options)
    

    # Trim the concatenated alignment using trimAl
    for fasta_file in filepaths:
        base_name = os.path.splitext(os.path.basename(fasta_file))[0] # Extract the filename without the .faa extension
        output_fasta = os.path.join(options.fasta_output_directory, f"{base_name}.fasta_aln")
        remove_gaps_with_trimal(fasta_file, output_fasta, options.gap_remove_threshold)
    
    alignment_files = myUtil.getAllFiles(options.fasta_output_directory,".fasta_aln")
    return alignment_files
    
def get_alignments(options):
    """
    02.11.22
        Args:
            filepaths with the files to be processed
            directory for the concatenated seqs file
        Alignment files should be concated if header is the same            
    """

    fasta_files = myUtil.getAllFiles(options.fasta_output_directory,".faa") 
    for fasta_file in fasta_files:
        input_dir = os.path.dirname(fasta_file)
        base_name = os.path.splitext(os.path.basename(fasta_file))[0] # Extract the filename without the .faa extension
        output_fasta = os.path.join(input_dir, f"{base_name}.faa_aln") # Create the output filename with .aln extension in the same directory as the input file

        #Align with default mafft
        align_fasta_with_mafft(fasta_file, output_fasta, options.cores)
    
    alignment_files = myUtil.getAllFiles(options.fasta_output_directory,".faa_aln")            
    
    return alignment_files

def get_executable_dir():
    """
    Get the directory of the current executable or script.
    This works whether the script is compiled or run directly as a Python script.
    """
    if getattr(sys, 'frozen', False):
        # If the program is compiled, sys.frozen is True, and sys.executable gives the path to the executable
        return os.path.dirname(sys.executable)
    else:
        # If running as a script, __file__ gives the path to the script
        return os.path.dirname(os.path.abspath(__file__))

def find_executable(executable):
    """
    Find the MAFFT executable.
    First check in the system's PATH, then in the local ./bin directory relative to the executable/script.
    Returns the path to the MAFFT executable.
    """
    # Check if MAFFT is in the system's PATH
    executable_path = shutil.which(f"{executable}")
    if executable_path:
        return executable_path
    
    # If not found, check in the local ./bin directory relative to the executable or script
    executable_dir = get_executable_dir()
    local_executable_path = os.path.join(executable_dir, "bin", f"{executable}")
    if os.path.isfile(local_executable_path) and os.access(local_executable_path, os.X_OK):
        return local_executable_path
    
    raise FileNotFoundError(f"{executable} executable not found in system PATH or local bin directory.")

def align_fasta_with_mafft(input_fasta, output_fasta, cores=2):
    """
    Aligns a given .faa FASTA file using MAFFT.
    
    Args:
        input_fasta: Path to the input FASTA (.faa) file.
        output_fasta: Path to the output aligned FASTA file.
    """
    # Find MAFFT executable
    mafft = find_executable("mafft")  # Ensure this function finds the MAFFT executable

    # Run MAFFT alignment
    try:
        with open(output_fasta, "w") as output_file:
            # Pass the file object to stdout and stderr
            subprocess.run([mafft, "--thread", str(cores), "--auto", input_fasta], stdout=output_file, stderr=subprocess.PIPE, check=True)
        print(f"Alignment complete: {input_fasta} -> {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during MAFFT alignment: {e.stderr.decode('utf-8')}")
        raise
    except FileNotFoundError:
        print("MAFFT executable not found. Please make sure MAFFT is installed and in your system's PATH.")



def remove_gaps_with_trimal(input_fasta, output_alignment, gap_threshold=0.95):
    """
    Remove columns with gaps using trimAl based on the specified gap threshold.

    Args:
        input_fasta: Path to the input FASTA file (aligned).
        output_alignment: Path to the output trimmed alignment file.
        gap_threshold: Proportion of gaps allowed in a column (default: 0.95).
    """
    trimal = find_executable("trimal")
    try:
        # Run the trimAl command with the gap threshold
        subprocess.run([
            trimal, 
            "-in", input_fasta, 
            "-out", output_alignment, 
            "-gt", str(gap_threshold),
            "-keepheader"
        ], check=True)
        #print(f"Trimming complete: {input_fasta} -> {output_alignment}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during trimming with trimAl: {e.stderr.decode('utf-8')}")
        raise

