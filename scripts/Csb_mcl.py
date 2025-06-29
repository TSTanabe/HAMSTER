#!/usr/bin/python


import csv
import os
import numpy as np
from typing import Dict, Any, Set, Tuple, Optional

from . import Csb_proteins
from . import myUtil

logger = myUtil.logger





def csb_mcl_datasets(options: Any, reference_dict: Dict[str, Set[str]]) -> Dict[str, str]:
    """
    Runs DIAMOND self-blast and MCL clustering for all .faa files in the given directory.

    Args:
        options (Any): Settings object, must contain phylogeny_directory and mcl_inflation.
        reference_dict (dict): {domain: set(reference protein IDs)}.

    Returns:
        dict: {domain: mcl_cluster_file_path}
    """
    
    # Diamond self blast
    blast_results_dict = run_diamond_self_blast(options, options.phylogeny_directory) # domain => result_file
    
    # Markov Chain Clustering
    mcl_results_dict = run_mcl_for_all_domains(blast_results_dict, options.phylogeny_directory, options.mcl_inflation)
    
    mcl_results_dict = {f.replace("_mcl_clusters.txt", ""): os.path.join(options.phylogeny_directory, f) for f in os.listdir(options.phylogeny_directory) if f.endswith("_mcl_clusters.txt")}

    return mcl_results_dict


##################################################################

def run_diamond_self_blast(options: Any, directory: str) -> Dict[str, str]:
    """
    Runs DIAMOND BLASTp self-search for all .faa files in the directory.
    Skips if output already exists.

    Args:
        options: DIAMOND settings.
        directory (str): Input directory with .faa files.

    Returns:
        dict: {domain: path_to_diamond_selfblast.tsv}
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
            logger.debug(f"Found existing MCL results for {input_prefix}")
            continue

        if os.path.exists(output_file_path):
            logger.debug(f"DIAMOND BLAST results already exist: {output_file_path}")
            output_files_dict[input_prefix] = output_file_path
            continue

        logger.debug(f"Starting DIAMOND self-BLASTp for {faa_file}")

        # Run DIAMOND BLASTp self-search
        blast_results_path = DiamondSearch(
            faa_path, faa_path, options.cores, 
            options.mcl_evalue, options.mcl_searchcoverage, options.mcl_minseqid, 
            options.mcl_hit_limit, options.alignment_mode, options.mcl_sensitivity
        )

        # Debugging: Check if DIAMOND output was actually created
        if not os.path.exists(blast_results_path):
            logger.error(f"DIAMOND output file {blast_results_path} not found after self-BLASTp for MCL")
            continue

        # If DIAMOND produces an output file, rename it to match the input file prefix
        os.rename(blast_results_path, output_file_path)
        output_files_dict[input_prefix] = output_file_path

    return output_files_dict



def run_mcl_clustering(
    mcl_input_file: str,
    mcl_output_file: str,
    mcl_inflation: float,
    domain: str = "unknown"
) -> Optional[str]:
    """
    Runs MCL clustering on a DIAMOND tabular output file.

    Args:
        mcl_input_file (str): Path to DIAMOND BLAST output.
        mcl_output_file (str): Output file path for MCL clusters.
        mcl_inflation (float): MCL inflation parameter.
        domain (str): Domain label (for logging).

    Returns:
        str: Path to MCL output file if created, else None.
    """
    
    if os.path.isfile(mcl_output_file):
        logger.debug(f"Found existing MCL clustering for {domain}")
        return mcl_output_file
    
    # Step 2: Run MCL clustering
    logger.debug(f"Starting MCL clustering for {domain}")
    mcl = myUtil.find_executable("mcl")
    os.system(f"{mcl} {mcl_input_file} --abc -I {mcl_inflation} -o {mcl_output_file} > /dev/null 2>&1")


    # Step 3: Display a preview of the MCL clusters
    return mcl_output_file


def run_mcl_for_all_domains(
    output_files_dict: Dict[str, str],
    mcl_directory: str,
    mcl_inflation: float
) -> Dict[str, str]:
    """
    Runs MCL clustering for all DIAMOND output files.

    Args:
        output_files_dict (dict): {domain: path_to_diamond_selfblast.tsv}
        mcl_directory (str): Directory for output.
        mcl_inflation (float): MCL inflation value.

    Returns:
        dict: {domain: path_to_mcl_clusters.txt}
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

def select_hits_by_csb_mcl(
    options: Any,
    mcl_output_dict: Dict[str, str],
    reference_dict: Dict[str, Set[str]],
    density_threshold: Optional[float] = None,
    reference_threshold: Optional[float] = None
) -> Tuple[Dict[str, Set[str]], Dict[str, Dict[str, float]]]:
    """
    Selects clusters from MCL output based on density and reference thresholds.

    Args:
        options: Not used here but passed for consistency.
        mcl_output_dict (dict): {domain: mcl_cluster_file_path}
        reference_dict (dict): {domain: set(reference sequence IDs)}
        density_threshold (float): Minimum reference density required.
        reference_threshold (float): Minimum reference coverage required.

    Returns:
        tuple:
            - all_clusters (dict): {domain: set(selected protein IDs)}
            - all_cutoffs (dict): {domain: {'density_threshold': x, 'reference_threshold': y}}
    """
    
    all_clusters = {}
    all_cutoffs = {}
    # Process reference_dict to extract the actual domain names
    processed_reference_dict = {
        key.split("_", 1)[-1]: value for key, value in reference_dict.items()
    }

    for domain, mcl_file in mcl_output_dict.items():
        logger.debug(f"Selecting hits by sequence clustering for {domain}")

        # Get reference sequences for the domain (if exists in processed reference dict)
        reference_sequences = processed_reference_dict.get(domain, set())
        if not reference_sequences:
            logger.warning(f"No reference sequences found for {domain}")
            continue
        
        
        local_density_thrs = density_threshold
        local_reference_thrs = reference_threshold
        
        # Step 1 Calculate the optimal density and reference thresholds
        local_density_thrs, local_reference_thrs = calculate_optimized_mcl_threshold(mcl_file, domain, reference_sequences,fixed_density_threshold=local_density_thrs, fixed_reference_threshold=local_reference_thrs)
        
        # Get fallback values if thresholds were not provided
        if local_density_thrs is None:
            logger.warning(f"MCL density threshold was not calculated for {domain}, fallback 0.1")
            local_density_thrs = 0.1
        if local_reference_thrs is None:
            logger.warning(f"MCL reference threshold was not calculated for {domain}, fallback 0.001")
            local_reference_thrs = 0.01
        
        
        #print(f"[INFO] {domain} with local_density_thrs {local_density_thrs} and reference thrs {local_reference_thrs}")
        
        # Step 1 Process the MCL file for this domain
        mcl_domain_clusters_dict = process_single_mcl_file(mcl_file, domain, reference_sequences, local_density_thrs, local_reference_thrs)
        
        
        # Step 2 Combine all individual cluster sets into one merged set
        combined_set = set().union(*mcl_domain_clusters_dict.values()).union(reference_sequences)

        all_clusters[domain] = combined_set
        all_cutoffs[domain] = {'density_threshold': local_density_thrs, 'reference_threshold': local_reference_thrs}
        
    return all_clusters, all_cutoffs


def calculate_optimized_mcl_threshold(
    mcl_file: str,
    domain: str,
    reference_sequences: Set[str],
    fixed_density_threshold: Optional[float] = None,
    fixed_reference_threshold: Optional[float] = None
) -> Tuple[Optional[float], Optional[float]]:
    """
    Finds optimal (density, reference) thresholds for MCL cluster selection based on F1 score.

    Args:
        mcl_file (str): MCL clusters file.
        domain (str): Domain for log output.
        reference_sequences (set): Known references.
        fixed_density_threshold (float): Optional override.
        fixed_reference_threshold (float): Optional override.

    Returns:
        (density_threshold, reference_threshold)
    """
    
    def parse_mcl_clusters(path):
        with open(path, 'r') as f:
            return [line.strip().split() for line in f if line.strip()]
    
    clusters = parse_mcl_clusters(mcl_file)
    n_ref_total = len(reference_sequences)

    if n_ref_total == 0:
        logger.warning("No reference sequence provided")
        return fixed_density_threshold,fixed_reference_threshold

    cluster_metrics = []
    for cluster in clusters:
        cluster_set = set(cluster)
        cluster_size = len(cluster_set)
        ref_ids = cluster_set & reference_sequences
        ref_count = len(ref_ids)
        non_ref_count = cluster_size - ref_count


        if cluster_size == 0 or ref_count == 0:
            #print("[WARN] No reference sequence in cluster found")
            continue

        ref_ratio_in_cluster = ref_count / cluster_size
        ref_ratio_of_total = ref_count / n_ref_total

        cluster_metrics.append({
            'ref_ratio_in_cluster': ref_count / cluster_size,
            'ref_ratio_of_total': len(ref_ids) / n_ref_total,
            'ref_count': ref_count,
            'non_ref_count': non_ref_count,
            'ref_ids': ref_ids  # neue Info
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

            TP_ids = set().union(*(c['ref_ids'] for c in selected))
            TP = len(TP_ids)
            FP = sum(c['non_ref_count'] for c in selected)
            FN = n_ref_total - TP

            if TP + FP == 0 or TP + FN == 0:
                continue

            precision = TP / (TP + FP)
            recall = TP / (TP + FN)

            # Sanity check
            if not (0 <= precision <= 1):
                logger.error(f"Invalid precision: {precision} (TP={TP}, FP={FP})")
                continue
            if not (0 <= recall <= 1):
                logger.error(f"Invalid recall: {recall} (TP={TP}, FN={FN})")
                continue

            f1 = 2 * (precision * recall) / (precision + recall)

            if f1 > 1:
                logger.error(f"F1 > 1: TP={TP}, FP={FP}, FN={FN}, Precision={precision}, Recall={recall}, F1={f1}")
                continue

            if f1 >= best_score:
                best_score = f1
                best_coords = (x_thresh, y_thresh)


    logger.debug(f"{domain} optimal MCL cluster inclusion thresholds: density = {best_coords[1]}, reference = {best_coords[0]}, F1 = {best_score:.3f}")
    return best_coords[1], best_coords[0]  # Return in order: (density, reference)


def process_single_mcl_file(
    mcl_file: str,
    domain: str,
    reference_sequences: Set[str],
    density_threshold: float,
    reference_fraction_threshold: float
) -> Dict[str, Set[str]]:
    """
    Processes a single MCL output file and extracts clusters above threshold.

    Args:
        mcl_file (str): Path to cluster file.
        domain (str): Domain label.
        reference_sequences (set): Known references.
        density_threshold (float): Min ratio of reference in cluster.
        reference_fraction_threshold (float): Min coverage of total reference.

    Returns:
        dict: {cluster_label: set(proteinIDs)}
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

    logger.debug(f"  Reference sequences provided:    {len(reference_sequences)}")
    logger.debug(f"  Total sequences in kept clusters:{len(all_cluster_seqs)}")
    logger.debug(f"     ├─ from references:           {len(ref_in_clusters)}")
    logger.debug(f"     └─ newly selected sequences:  {len(new_candidates)}")

    return cluster_dict



        

def DiamondSearch(
    path: str,
    query_fasta: str,
    cores: int,
    evalue: float,
    coverage: float,
    minseqid: float,
    diamond_report_hits_limit: int,
    alignment_mode: int = 2,
    sensitivity: str = "ultra-sensitive"
) -> str:
    """
    Runs DIAMOND BLASTP on a FASTA against itself, output in tab format.

    Args:
        path (str): FASTA file path.
        query_fasta (str): Query FASTA.
        cores (int): Number of threads.
        evalue (float): E-value cutoff.
        coverage (float): Coverage threshold.
        minseqid (float): Minimum sequence identity.
        diamond_report_hits_limit (int): Report limit.
        alignment_mode (int): Not used, for consistency.
        sensitivity (str): Sensitivity flag.

    Returns:
        str: Path to DIAMOND output file (.diamond.tab)
    """    

    diamond = myUtil.find_executable("diamond")
    # Erstellen der Diamond-Datenbank
    target_db_name = f"{path}.dmnd"
    os.system(f'{diamond} makedb --quiet --in {path} -d {target_db_name} --threads {cores} 1>/dev/null 0>/dev/null')
    
    # Suchergebnisse-Dateien
    output_results_tab = f"{path}.diamond.tab"

    os.system(f'{diamond} blastp --quiet --{sensitivity} -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} -k {diamond_report_hits_limit} --outfmt 6 sseqid qseqid bitscore 1>/dev/null 0>/dev/null')
    
    #output format hit query evalue score identity alifrom alito
    
    return output_results_tab

def validate_mcl_cluster_paths(
    path_dict: Dict[str, str],
    result_files_directory: str,
    target_dir: str = "Protein_Phylogeny"
) -> Dict[str, str]:
    """
    Validates a dictionary of {key: path_to_mcl_file}. For each path:
    - If the file exists: keep it.
    - If not: try looking in result_files_directory/Protein_Phylogeny/.
    - If that works, update the path.
    - If neither exists, remove the entry.

    Args:
        path_dict (dict): {key: mcl_path}
        result_files_directory (str): Root directory to check for fallback.
        target_dir (str): Fallback subdir.

    Returns:
        dict: Only valid MCL file paths by key.
    """
    validated = {}

    for key, path in path_dict.items():
        if os.path.isfile(path):
            validated[key] = path
        else:
            fallback_path = os.path.join(result_files_directory, target_dir, os.path.basename(path))
            if os.path.isfile(fallback_path):
                validated[key] = fallback_path
            else:
                logger.warning(f"No valid MCL file for {key}: '{path}' or fallback '{fallback_path}'")

    if not validated:
        logger.error("No valid MCL cluster files found in any specified or fallback location.")
        return {}

    return validated









