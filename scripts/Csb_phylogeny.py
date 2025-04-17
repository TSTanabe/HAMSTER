#!/usr/bin/python

import sqlite3
import os
import subprocess
from multiprocessing import Pool
from collections import defaultdict
from . import Alignment

from Bio import Phylo

#first get the csb with their identifiers. make sets of appearence as value to key name of the csb
#do not compare the csb appearences as they jaccard itself should have managed it.


def csb_phylogeny_datasets(options, training_datasets):

    # Align the protein family
    family_alignment_files = Alignment.initial_alignments(options, options.phylogeny_directory)
    
    # Calculate phylogeny for each alignment file
    tree_files = Alignment.calculate_phylogeny_parallel(options, family_alignment_files) #creates .tree files in the phylogeny directory

    # Find the Last Common Ancestors (LCA) for each potential training set and write the fasta file
    monophylums = get_last_common_ancestor_fasta(options,training_datasets,tree_files,'')
    
    return monophylums

    
#################################################################################
############  Finding lca in phylogenetic trees  ################################
#################################################################################




def find_lca_and_monophyly(protein_ids, tree_file):
    """
    Given a phylogenetic tree and a set of protein identifiers, 
    this function finds the last common ancestor (LCA) of the given proteins,
    checks if the clade is monophyletic, and returns relevant information.
    
    Args:
    - protein_ids: A list of protein identifiers (as present in the tree tips).
    - tree_file: Path to the Newick file containing the tree.
    
    Returns:
    - LCA: The last common ancestor of the given proteins.
    - is_monophylum: True if the clade is monophyletic, False otherwise.
    - included_identifiers: If not monophyletic, returns all included protein identifiers in the LCA clade.
    """
    # Load the phylogenetic tree
    tree = Phylo.read(tree_file, 'newick')
    
    # Get all terminal names (tree tip identifiers)
    tree_tip_names = {tip.name for tip in tree.get_terminals()}
    
    # Filter protein_ids to include only those present in the tree
    filtered_protein_ids = [pid for pid in protein_ids if pid in tree_tip_names]
    
    # If no valid protein IDs remain after filtering, return early
    if not filtered_protein_ids:
        print(f"Warning: None of the provided protein IDs are in the tree.")
        return None, None

    # Find the LCA of the filtered protein identifiers
    try:
        lca = tree.common_ancestor(filtered_protein_ids)
    except Exception as e:
        print(f"Error: finding LCA in tree file: {tree_file}")
        print(f"Original error: {e}")
        return None, None
    
    # Check if the clade is monophyletic
    is_monophylum = tree.is_monophyletic(filtered_protein_ids)
    
    if is_monophylum:
        return lca, set()  # LCA and no additional identifiers needed
    else:
        # Get all the identifiers in the LCA clade if it's not monophyletic
        included_identifiers = {tip.name for tip in lca.get_terminals()}
        return lca, included_identifiers


def get_lca_worker_task(args):
    """
    Worker function to process one proteinID_frozenset and its associated keyword_domain_pairs.
    Args:
        args: A tuple containing (proteinID_frozenset, keyword_domain_pairs, trees_dict, tree_prefix, key_prefix, csb_dictionary)
    Returns:
        tuple: The resulting monophyly entry and its associated filtered monophylum.
    """
    proteinID_set, domain, trees_dict, tree_prefix = args

    try:
        # Initialize monophyly
        monophyly = {}
        
        # Attempt to get the tree
        tree = trees_dict.get(domain)
        if not tree:
            print(f"Skipping {domain}. Tree is not available")
            return None  # Skip if tree is missing

        # Get the monophyletic group
        lca, clade_identifier_set = find_lca_and_monophyly(proteinID_set, tree)

        if clade_identifier_set:
            monophyly["lca_"+domain] = clade_identifier_set
        else:
            return {},{} #if no addition was found and set is already monophylum, return emtpy dicts

        
        return monophyly

    except Exception as e:
        print(f"Error: processing {proteinID_set}: {e}")
        return {},{}
        
        
        
def get_last_common_ancestor_fasta(options, grouped, trees_dict, tree_prefix):
    """
    Main function to process grouped data in parallel using multiprocessing.
    """

    # Prepare arguments for the worker function
    worker_args = [
        (proteinID_set, grouped_key.replace("grp0_", "", 1), trees_dict, tree_prefix)
        for grouped_key, proteinID_set in grouped.items()
    ]


    # Use multiprocessing to fork the outer loop
    monophyly_results = []
    with multiprocessing.Pool(processes=options.cores) as pool:
        results = pool.map(get_lca_worker_task, worker_args)

    # Collect results from workers
    monophyly = {}
    for result in results:
        if result:
            worker_monophyly = result
            monophyly.update(worker_monophyly)

    return monophyly


#########################################################################################################        
##################################### Hits sorted by query hit ##########################################
#########################################################################################################


def get_domain_superfamily_proteinIDs(blast_table, domain, query_length, tolerance=0.5):
    """
    Filters protein IDs based on length constraints relative to the query length.

    Parameters:
        blast_table (str): Path to the BLAST table file.
        domain (str): Domain to grep in the BLAST table.
        query_length (float): The length of the query.
        tolerance (float): Allowed length deviation as a fraction (default is 0.5, i.e., ±50%).

    Returns:
        set: A set of protein IDs that satisfy the length constraint.
    """
    proteinIDs = set()

    # Run the grep command to find lines containing the domain
    result = subprocess.run(['grep', domain, blast_table], stdout=subprocess.PIPE, text=True)
    lines = result.stdout.splitlines()  # Split output into lines

    # Length bounds
    lower_bound = (1 - tolerance) * query_length
    upper_bound = (1 + tolerance) * query_length

    for line in lines:
        columns = line.split('\t')  # Assuming columns are tab-separated
        if columns:
            try:
                # Extract relevant columns
                hit_proteinID = columns[0]
                hsp_start = int(float(columns[4]))
                hsp_end = int(float(columns[5]))

                # Calculate hit length
                hit_length = hsp_end - hsp_start

                # Check if hit length is within bounds
                if lower_bound <= hit_length <= upper_bound:
                    proteinIDs.add(hit_proteinID)

            except (ValueError, IndexError):
                # Skip lines with parsing errors
                print(f"Skipping line due to parsing error: {line}")
                continue

    return proteinIDs

    
def get_sequence_legth(file_path):
    #get the sequence length per query name without the numbering. If multiple queries have the same
    #name, then take the average
    
    sequence_data = {}  # Dictionary to store lengths per identifier
    sequence_counts = {}  # Dictionary to track counts of sequences per identifier
    
    with open(file_path, 'r') as fasta_file:
        current_id = None
        current_sequence = []
        
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                # Process previous sequence if any
                if current_id:
                    sequence_length = len("".join(current_sequence))
                    if current_id in sequence_data:
                        sequence_data[current_id] += sequence_length
                        sequence_counts[current_id] += 1
                    else:
                        sequence_data[current_id] = sequence_length
                        sequence_counts[current_id] = 1
                
                # Extract the identifier, considering only the part before "___"
                current_id = line[1:].split('___')[1]
                current_sequence = []  # Reset sequence for the new identifier
            else:
                current_sequence.append(line)
        
        # Process the last sequence in the file
        if current_id:
            sequence_length = len("".join(current_sequence))
            if current_id in sequence_data:
                sequence_data[current_id] += sequence_length
                sequence_counts[current_id] += 1
            else:
                sequence_data[current_id] = sequence_length
                sequence_counts[current_id] = 1
    
    # Calculate average length for identifiers with multiple sequences
    averaged_data = {}
    for key in sequence_data:
        averaged_data[key] = sequence_data[key] / sequence_counts[key]
    
    return averaged_data   



#########################################################################################################
#################      Routine for the sequences below PAM regression        ############################
#########################################################################################################

def _analyze_single_tree(args):
    domain, tree_path, unplausible_proteins, grp1_proteins, max_branch_length, max_neighbour_distance = args
    close_grp1_hits = set()
    close_external_hits = set()

    try:
        tree = Phylo.read(tree_path, "newick")
        terminals = tree.get_terminals()
        name_to_terminal = {t.name: t for t in terminals}
        print(f"Analysing placement of hits for {domain} in the phylogenetic tree")
    except Exception:
        return domain, close_grp1_hits, close_external_hits

    for pid in unplausible_proteins:
        if pid not in name_to_terminal:
            continue
        clade = name_to_terminal[pid]

        # skip if branch too long
        if clade.branch_length and clade.branch_length > max_branch_length:
            continue

        path = tree.get_path(clade)
        if len(path) < 2:
            continue

        parent = path[-2]
        sibling_clades = [c for c in parent.clades if c != clade]
        found = False

        for sibling in sibling_clades:
            sib_terms = sibling.get_terminals()
            for neighbor in sib_terms:
                nname = neighbor.name
                if not nname or nname == pid:
                    continue
                if nname in unplausible_proteins:
                    continue
                if nname in grp1_proteins:
                    try:
                        dist = tree.distance(name_to_terminal[pid], name_to_terminal[nname])
                        if dist <= max_neighbour_distance:
                            close_grp1_hits.add(pid)
                    except Exception as e:
                        print(f"Distance error between {pid} and {nname}: {e}")
                    found = True
                    break
                else:
                    close_external_hits.add(pid)
                    found = True
                    break
            if found:
                break

    return domain, close_grp1_hits, close_external_hits

def analyze_unplausible_proteins_in_trees_parallel(tree_dir, unplausible_dict, grp1_dict, max_branch_length=0.8, max_neighbour_distance=0.15, cores=4):
    dict_close_grp1 = defaultdict(set)
    dict_close_external = defaultdict(set)
    print(f"Acceptable distance to a reference sequence {max_neighbour_distance}")
    args_list = []
    for domain, unplausible_proteins in unplausible_dict.items():
        tree_path = os.path.join(tree_dir, domain + ".tree")
        if os.path.isfile(tree_path):
            grp1_proteins = grp1_dict.get(domain, set())
            args_list.append((domain, tree_path, unplausible_proteins, grp1_proteins, max_branch_length, max_neighbour_distance))

    with Pool(processes=cores) as pool:
        results = pool.map(_analyze_single_tree, args_list)

    for domain, close_grp1_hits, close_external_hits in results:
        dict_close_grp1[domain].update(close_grp1_hits)
        dict_close_external[domain].update(close_external_hits)
    
    return dict_close_grp1, dict_close_external

