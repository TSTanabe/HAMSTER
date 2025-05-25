#!/usr/bin/python
import os
import sys
import shutil
import subprocess
import multiprocessing

from . import myUtil


####################################################
##########   MAFFT Alignments routines    ##########
####################################################

#concats .fasta_aln files in a directory. if none are found
#it aligns with mafft and removes gaps with trimal at 95 % 
            
            
def initial_alignments(options, fasta_output_directory):

    filepaths = make_alignments(options, fasta_output_directory)
    

    # Trim the concatenated alignment using trimAl
    for fasta_file in filepaths:
        base_name = os.path.splitext(os.path.basename(fasta_file))[0] # Extract the filename without the .faa extension
        output_fasta = os.path.join(fasta_output_directory, f"{base_name}.fasta_aln")
        
        if os.path.exists(output_fasta):
            continue
            
        remove_gaps_with_trimal(fasta_file, output_fasta, options.gap_remove_threshold)
    
    alignment_files = myUtil.getAllFiles(fasta_output_directory,".fasta_aln")
    return alignment_files
    
def make_alignments(options,fasta_output_directory):
    """
    02.11.22
        Args:
            filepaths with the files to be processed
            directory for the concatenated seqs file
        Alignment files should be concated if header is the same            
    """

    fasta_files = myUtil.getAllFiles(fasta_output_directory,".faa") 
    for fasta_file in fasta_files:
        input_dir = os.path.dirname(fasta_file)
        base_name = os.path.splitext(os.path.basename(fasta_file))[0] # Extract the filename without the .faa extension
        output_fasta = os.path.join(input_dir, f"{base_name}.faa_aln") # Create the output filename with .aln extension in the same directory as the input file

        if os.path.isfile(output_fasta):
            print(f"[SKIP] {output_fasta} already exists")
            continue
        #Align with default mafft
        align_fasta_with_mafft(fasta_file, output_fasta, options.cores)
    
    alignment_files = myUtil.getAllFiles(fasta_output_directory,".faa_aln")            
    
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
        print(f"[INFO] Alignment complete: {output_fasta}")
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode('utf-8') if e.stderr else "Unknown error"
        print(f"[ERROR] occurred during MAFFT alignment: {error_message}")
        print(f"[SKIP] file {input_fasta} due to the error.")

    except FileNotFoundError:
        print("[ERROR] MAFFT executable not found. Please make sure MAFFT is installed and in your system's PATH.")



def remove_gaps_with_trimal(input_fasta, output_alignment, gap_threshold=0.05):
    """
    Remove columns with gaps using trimAl based on the specified gap threshold.

    Args:
        input_fasta: Path to the input FASTA file (aligned).
        output_alignment: Path to the output trimmed alignment file.
        gap_threshold: Proportion of gaps allowed in a column (default: 0.95).
    """
    trimal = find_executable("trimal")
    print(f"[INFO] Trimming {input_fasta} with {gap_threshold} gap threshold")
    try:
        # Run the trimAl command with the gap threshold
        subprocess.run([
            trimal, 
            "-in", input_fasta, 
            "-out", output_alignment, 
            "-gt", str(gap_threshold),
            "-keepseqs",
            "-keepheader"
        ], check=True)
        #print(f"Trimming complete: {input_fasta} -> {output_alignment}")
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode('utf-8') if e.stderr else "Unknown error"
        print(f"[ERROR] occurred during trimming with trimAl: {error_message}")
        print(f"[SKIP] file {input_fasta} due to the error.")



def run_fasttree_on_alignment(alignment_file, cores=2):
    """
    Aligns a given .faa FASTA file using MAFFT.
    
    Args:
        alignment_file: Path to the input FASTA (.fasta_aln) file.
    """
    # Find MAFFT executable
    veryfasttree = find_executable("VeryFastTree-avx2")  # Ensure this function finds the MAFFT executable

    # Replace the original suffix with '.tree'
    base_name = os.path.splitext(alignment_file)[0]  # Get the file name without the extension
    output_tree_file = base_name + ".tree"  # Add the .tree extension

    if os.path.exists(output_tree_file):
        print(f"Tree file already exists for {alignment_file}, skipping FastTree.")
        return output_tree_file  # Return the existing tree file


    # Run fasttree alignment
    try:
        with open(output_tree_file, "w") as f:
            subprocess.run(
                [veryfasttree, "-nosupport", "-threads", str(cores), "-ext", "AVX2", alignment_file],
                stdout=f,         # Direct stdout to the output file
                stderr=subprocess.PIPE,  # Keep stderr for error handling
                check=True
            )
        print(f"VeryFastTree-avx2 complete: {alignment_file}")
        return output_tree_file
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode('utf-8') if e.stderr else "Unknown error"
        print(f"Error occurred during tVeryFastTree-avx2: {error_message}")
        print(f"Skipping file {alignment_file} due to the error.")
    except FileNotFoundError:
        print("VeryFastTree-avx2 executable not found. Please make sure VeryFastTree-avx2 is installed and in your system's PATH or in the local bin directory.")





def calculate_phylogeny_parallel(options, alignment_files):
    """
    Parallelizes the execution of FastTree on multiple alignment files using multiprocessing.

    Args:
        alignment_files (list): A list of paths to alignment files.
        options: An object containing options, specifically the number of cores available (options.cores).

    Returns:
        dict: A dictionary where the key is the filename (without extension) and the value is the path to the tree file.
    """
    # Create a multiprocessing pool with the specified number of cores
    with multiprocessing.Pool(processes=options.cores) as pool:
        # Distribute the workload among the cores and collect the results (tree file paths)
        tree_files = pool.map(run_fasttree_on_alignment, alignment_files)

    # Create a dictionary to map the filenames (without extension) to the tree file paths
    tree_dict = {}
    for tree_file in tree_files:
        if tree_file is not None:
            filename = os.path.splitext(os.path.basename(tree_file))[0]
            tree_dict[filename] = tree_file
    
    return tree_dict




















