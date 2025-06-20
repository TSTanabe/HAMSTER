#!/usr/bin/python

import csv
import os
import numpy as np
from . import Csb_proteins
from . import myUtil




def csb_mcl_datasets(options, reference_dict):

    # Diamond self blast
    blast_results_dict = run_diamond_self_blast(options, options.phylogeny_directory) # domain => result_file
    
    # Markov Chain Clustering
    mcl_results_dict = run_mcl_for_all_domains(blast_results_dict, options.phylogeny_directory, options.mcl_inflation)
    
    mcl_results_dict = {f.replace("_mcl_clusters.txt", ""): os.path.join(options.phylogeny_directory, f) for f in os.listdir(options.phylogeny_directory) if f.endswith("_mcl_clusters.txt")}

    return mcl_results_dict





##################################################################

##################################################################

def run_diamond_self_blast(options, directory):
    """
    Runs a DIAMOND BLASTp self-search for all .faa files in a directory.
    If the output file already exists, the DIAMOND search is skipped.

    Parameters:
    - options: Settings for the DIAMOND search (e.g., number of threads, memory limits).
    - directory: Directory containing the .faa files.
    """

    # List all .faa files in the directory
    faa_files = [f for f in os.listdir(directory) if f.endswith(".faa")]
    # Dictionary to store output file paths
    output_files_dict = {}
    
    for faa_file in faa_files:
        faa_path = os.path.join(directory, faa_file)

        # Generate the output file name
        input_prefix = os.path.splitext(faa_file)[0]  # Removes ".faa"
        output_file_path = os.path.join(directory, f"{input_prefix}_diamond_selfblast.tsv")
        mcl_output_file = os.path.join(directory, f"{input_prefix}_mcl_clusters.txt")

        # Debugging output: Check if files exist before processing
        if os.path.exists(mcl_output_file):
            print(f"[SKIP] Found existing MCL results for {input_prefix}")
            continue

        if os.path.exists(output_file_path):
            print(f"[SKIP] DIAMOND BLAST results already exist: {output_file_path}")
            output_files_dict[input_prefix] = output_file_path
            continue

        print(f"[INFO] Starting DIAMOND self-BLASTp for {faa_file}")

        # Run DIAMOND BLASTp self-search
        blast_results_path = DiamondSearch(
            faa_path, faa_path, options.cores, 
            options.mcl_evalue, options.mcl_searchcoverage, options.mcl_minseqid, 
            options.mcl_hit_limit, options.alignment_mode, options.mcl_sensitivity
        )

        # Debugging: Check if DIAMOND output was actually created
        if not os.path.exists(blast_results_path):
            print(f"[ERROR] DIAMOND output file {blast_results_path} not found after self-BLASTp for MCL")
            continue

        # If DIAMOND produces an output file, rename it to match the input file prefix
        os.rename(blast_results_path, output_file_path)
        output_files_dict[input_prefix] = output_file_path

    return output_files_dict



def run_mcl_clustering(mcl_input_file, mcl_output_file, mcl_inflation, domain='unknown'):
    """
    Converts a Diamond BLAST output file into an MCL-compatible format
    and runs MCL clustering.

    Args:
    - diamond_tab_file (str): Path to the Diamond BLAST output (tabular format).
    - mcl_input_file (str): Name of the MCL input file.
    - mcl_output_file (str): Name of the MCL output file.
    """
    
    if os.path.isfile(mcl_output_file):
        print(f"[SKIP] Found existing MCL clustering for {domain}")
        return mcl_output_file
    
    # Step 2: Run MCL clustering
    print(f"[INFO] Starting MCL clustering for {domain}")
    mcl = myUtil.find_executable("mcl")
    os.system(f"{mcl} {mcl_input_file} --abc -I {mcl_inflation} -o {mcl_output_file} > /dev/null 2>&1")


    # Step 3: Display a preview of the MCL clusters
    return mcl_output_file


def run_mcl_for_all_domains(output_files_dict, mcl_directory, mcl_inflation):
    """
    Iterates over a dictionary of Diamond BLAST output files and runs MCL clustering for each.

    Args:
    - output_files_dict (dict): Dictionary where keys are domains and values are Diamond BLAST output file paths.
    - mcl_directory (str): Path to the directory where MCL output files should be stored.

    Returns:
    - mcl_output_dict (dict): Dictionary where keys are domains and values are the corresponding MCL output file paths.
    """

    # Ensure the output directory exists
    if not os.path.exists(mcl_directory):
        os.makedirs(mcl_directory)

    # Dictionary to store MCL output file paths
    mcl_output_dict = {}

    # Iterate over all Diamond BLAST output files
    for domain, diamond_tab_file in output_files_dict.items():

        # Define input and output file paths for MCL
        mcl_output_file = os.path.join(mcl_directory, f"{domain}_mcl_clusters.txt")

        # Run MCL clustering
        mcl_output_path = run_mcl_clustering(diamond_tab_file, mcl_output_file, mcl_inflation, domain)

        # Store the MCL output file in the dictionary
        if mcl_output_path:
            mcl_output_dict[domain] = mcl_output_path

    return mcl_output_dict


##################################################################################################################
#### Main for the mcl cluster iteration ####

def select_hits_by_csb_mcl(options, mcl_output_dict, reference_dict, density_threshold=None, reference_threshold=None):
    """
    Iterates over multiple MCL output files, extracts reference sequences, 
    and processes them using process_single_mcl_file.
    Selection is based on F1 score optimized density and reference threshold

    Args:
    - mcl_output_dict (dict): Dictionary where keys are domains and values are MCL output file paths. These are either grp0_ or without prefix
    - reference_dict (dict): Dictionary where keys need to be split ('_'), and values are sets of reference sequence IDs. keys are grp0_domain
    - density_threshold (float): Minimum reference density required.

    Returns:
    - all_clusters (dict): Merged dictionary of all high-density clusters.
    """
    all_clusters = {}
    all_cutoffs = {}
    # Process reference_dict to extract the actual domain names
    processed_reference_dict = {
        key.split("_", 1)[-1]: value for key, value in reference_dict.items()
    }

    for domain, mcl_file in mcl_output_dict.items():
        print(f"\n[INFO] Selecting hits by csb for {domain}")

        # Get reference sequences for the domain (if exists in processed reference dict)
        reference_sequences = processed_reference_dict.get(domain, set())
        if not reference_sequences:
            print(f"[WARN] No reference sequences found for {domain}")
            continue
        
        
        local_density_thrs = density_threshold
        local_reference_thrs = reference_threshold
        
        # Step 1 Calculate the optimal density and reference thresholds
        local_density_thrs, local_reference_thrs = calculate_optimized_mcl_threshold(mcl_file, domain, reference_sequences,fixed_density_threshold=local_density_thrs, fixed_reference_threshold=local_reference_thrs)
        
        # Get fallback values if thresholds were not provided
        if local_density_thrs is None:
            print(f"[WARN] MCL density threshold was not calculated for {domain}, fallback 0.1")
            local_density_thrs = 0.1
        if local_reference_thrs is None:
            print(f"[WARN] MCL reference threshold was not calculated for {domain}, fallback 0.001")
            local_reference_thrs = 0.01
        
        
        #print(f"[INFO] {domain} with local_density_thrs {local_density_thrs} and reference thrs {local_reference_thrs}")
        
        # Step 1 Process the MCL file for this domain
        mcl_domain_clusters_dict = process_single_mcl_file(mcl_file, domain, reference_sequences, local_density_thrs, local_reference_thrs)
        
        
        # Step 2 Combine all individual cluster sets into one merged set
        combined_set = set().union(*mcl_domain_clusters_dict.values()).union(reference_sequences)

        all_clusters[domain] = combined_set
        all_cutoffs[domain] = {'density_threshold': local_density_thrs, 'reference_threshold': local_reference_thrs}
        
    return all_clusters, all_cutoffs


def calculate_optimized_mcl_threshold(mcl_file, domain, reference_sequences, 
                                      fixed_density_threshold=None, 
                                      fixed_reference_threshold=None):
    """
    Calculates the optimal (density, reference) thresholds based on F1 score.

    Args:
        mcl_file (str): Path to MCL output file.
        domain (str): Domain name (used for logging/debugging).
        reference_sequences (set): Set of known reference sequence IDs.
        fixed_density_threshold (float or None): Optional fixed density threshold.
        fixed_reference_threshold (float or None): Optional fixed reference threshold.

    Returns:
        (float, float): Tuple of best (density_threshold, reference_threshold)
    """
    
    def parse_mcl_clusters(path):
        with open(path, 'r') as f:
            return [line.strip().split() for line in f if line.strip()]
    
    clusters = parse_mcl_clusters(mcl_file)
    n_ref_total = len(reference_sequences)

    if n_ref_total == 0:
        print("[WARN] No reference sequence provided")
        return fixed_density_threshold,fixed_reference_threshold

    cluster_metrics = []
    for cluster in clusters:
        cluster_set = set(cluster)
        cluster_size = len(cluster_set)
        ref_count = len(cluster_set & reference_sequences)
        non_ref_count = cluster_size - ref_count

        if cluster_size == 0 or ref_count == 0:
            #print("[WARN] No reference sequence in cluster found")
            continue

        ref_ratio_in_cluster = ref_count / cluster_size
        ref_ratio_of_total = ref_count / n_ref_total

        cluster_metrics.append({
            'ref_ratio_in_cluster': ref_ratio_in_cluster,
            'ref_ratio_of_total': ref_ratio_of_total,
            'ref_count': ref_count,
            'non_ref_count': non_ref_count
        })

    if not cluster_metrics:
        return fixed_density_threshold,fixed_reference_threshold

    # Threshold grids
    x_bins = [fixed_reference_threshold] if fixed_reference_threshold is not None else np.round(np.arange(0.001, 1.01, 0.005), 3)
    y_bins = [fixed_density_threshold] if fixed_density_threshold is not None else np.round(np.arange(0.01, 1.01, 0.005), 3)

    best_score = 0.0
    best_coords = (None, None)

    for x_thresh in x_bins:
        for y_thresh in y_bins:
            selected = [
                c for c in cluster_metrics
                if c['ref_ratio_of_total'] >= x_thresh and c['ref_ratio_in_cluster'] >= y_thresh
            ]

            TP = sum(c['ref_count'] for c in selected)
            FP = sum(c['non_ref_count'] for c in selected)
            FN = n_ref_total - TP

            if TP + FP == 0 or TP + FN == 0:
                continue

            precision = TP / (TP + FP)
            recall = TP / (TP + FN)
            f1 = 2 * (precision * recall) / (precision + recall)

            if f1 >= best_score:
                best_score = f1
                best_coords = (x_thresh, y_thresh)


    print(f"[INFO] {domain} optimal MCL cluster inclusion thresholds: density = {best_coords[1]}, reference = {best_coords[0]}, F1 = {best_score:.3f}")
    return best_coords[1], best_coords[0]  # Return in order: (density, reference)



def process_single_mcl_file(mcl_file, domain, reference_sequences, density_threshold, reference_fraction_threshold):
    """
    Processes a single MCL output file and extracts high-density clusters.

    Returns:
    - cluster_dict (dict): Dictionary where keys are "mcl_X_domain" and values are sets of sequence IDs.
    """

    cluster_dict = {}
    cluster_index = 1

    # Read MCL file and process clusters
    with open(mcl_file, "r") as f:
        for line in f:
            seqs = set(line.strip().split())

            ref_count = len(seqs & reference_sequences)
            cluster_size = len(seqs)
            ref_density = ref_count / cluster_size if cluster_size > 0 else 0
            ref_fraction = ref_count / len(reference_sequences) if reference_sequences else 0

            if ref_density >= density_threshold and ref_fraction >= reference_fraction_threshold:
                cluster_dict[f"mcl{cluster_index}_{domain}"] = seqs
                cluster_index += 1

    # Berechne alle enthaltenen Sequenzen in den akzeptierten Clustern
    all_cluster_seqs = set()
    for seqs in cluster_dict.values():
        all_cluster_seqs.update(seqs)

    ref_in_clusters = all_cluster_seqs & reference_sequences
    new_candidates = all_cluster_seqs - reference_sequences

    print(f"  Reference sequences provided:    {len(reference_sequences)}")
    print(f"  Total sequences in kept clusters:{len(all_cluster_seqs)}")
    print(f"     ├─ from references:           {len(ref_in_clusters)}")
    print(f"     └─ newly selected sequences:  {len(new_candidates)}")

    return cluster_dict



        

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

    os.system(f'{diamond} blastp --quiet --{sensitivity} -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} -k {diamond_report_hits_limit} --outfmt 6 sseqid qseqid bitscore 1>/dev/null 0>/dev/null')
    
    #output format hit query evalue score identity alifrom alito
    
    return output_results_tab

def validate_mcl_cluster_paths(path_dict, result_files_directory):
    """
    Validates a dictionary of {key: path_to_mcl_file}. For each path:
    - If the file exists: keep it.
    - If not: try looking in result_files_directory/Protein_Phylogeny/.
    - If that works, update the path.
    - If neither exists, remove the entry.

    Args:
        path_dict (dict): Dictionary with keys (e.g., domain names) and values as file paths.
        result_files_directory (str): Root directory to search fallback file under 'Protein_Phylogeny'.

    Returns:
        dict: Filtered and updated dictionary with only valid paths.

    Raises:
        FileNotFoundError: If no valid MCL files are found.
    """
    validated = {}

    for key, path in path_dict.items():
        if os.path.isfile(path):
            validated[key] = path
        else:
            fallback_path = os.path.join(result_files_directory, "Protein_Phylogeny", os.path.basename(path))
            if os.path.isfile(fallback_path):
                validated[key] = fallback_path
            else:
                print(f"[WARN] No valid MCL file for {key}: '{path}' or fallback '{fallback_path}'")

    if not validated:
        raise FileNotFoundError("No valid MCL cluster files found in any specified or fallback location.")

    return validated









