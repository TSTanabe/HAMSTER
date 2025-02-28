#!/usr/bin/python

import sqlite3
import os
import subprocess
import multiprocessing
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

def csb_phylogeny_target_sets(options):
    target_files = fetch_domains_superfamily_to_fasta(options, options.cross_validation_directory)
    return target_files #dictionary with domain => target file with TP and TN
        
def get_domain_key_list_pairs(input_dict, output_dict=None):
    # If no output_dict is provided, initialize a new empty dictionary
    if output_dict is None:
        output_dict = {}

    # Iterate over the input_dict values, which contain lists of tuples
    for tuples_list in input_dict.values():
        for key, domain in tuples_list:
            # Check if the domain exists in the output dictionary
            if domain not in output_dict:
                # If the domain is not present, create a new list for it
                output_dict[domain] = []
            # Append the key to the list of keys for the domain
            output_dict[domain].append(key)
    return output_dict

    
def fetch_protein_type_to_fasta(options, domain_keyword_dict):
#Fetches the TP proteinsIDs of one type regardless of the csb to combine csbs with a mixed distribution across the tree


    with sqlite3.connect(options.database_directory) as con:
        cur = con.cursor()

        # Iterate over the domain-keyword dictionary
        for domain, keywords in domain_keyword_dict.items():
            print(domain)
            modified_keywords = [keyword[1:] for keyword in keywords] #remove the marker character
            print(modified_keywords)
            # Use a parameterized query to check for the domain and the keywords
            query = '''
                SELECT DISTINCT P.proteinID, P.sequence
                FROM Proteins P
                JOIN Domains D ON P.proteinID = D.proteinID
                JOIN Keywords K ON P.clusterID = K.clusterID
                WHERE D.domain = ? AND K.keyword IN ({})
            '''.format(','.join('?' * len(keywords)))
            # Execute the query with the domain and the list of keywords
            cur.execute(query, (domain, *modified_keywords))

            proteins = cur.fetchall()

            # Write the results to the domain-specific FASTA file
            # Create the filename for this domain's output file
            #this leaves only the protein type not the cluster number
            name = domain.split('_')[-1]
            filename = f"lca_{name}.faa"
            output_fasta_path = os.path.join(options.phylogeny_directory, filename)
            
            # Use a set to track already written proteinIDs
            written_proteins = set()
            if os.path.exists(output_fasta_path):
                # Open the existing file and extract protein IDs
                with open(output_fasta_path, 'r') as fasta_file:
                    for line in fasta_file:
                        if line.startswith('>'):
                            # Extract protein ID from the FASTA header line (assuming ID is everything after '>')
                            protein_id = line[1:].strip().split()[0]  # Assuming ID is the first part after '>'
                            written_proteins.add(protein_id)

            #Append with the fetched sequences
            with open(output_fasta_path, 'a') as fasta_file:
                for proteinID, sequence in proteins:
                    if proteinID not in written_proteins:
                        fasta_file.write(f'>{proteinID}\n{sequence}\n')
                        written_proteins.add(proteinID)




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
        print(f"Error finding LCA in tree file: {tree_file}")
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






def fetch_neighbouring_domains_deprecated(database, protein_ids):
    """
    Fetches neighboring domains for a list of proteinIDs based on clusterID and gene order,
    and returns a dictionary where each proteinID is a key, and the value is a set of 
    neighboring domains (including the domain of the protein itself).

    Args:
        database (str): Pathway to the database file.
        protein_ids (list): List of proteinIDs to fetch neighbors for.

    Returns:
        dict: A dictionary where each key is a proteinID, and the value is a set of 
              domains of neighboring genes in the same cluster or just the domain of the 
              proteinID if it has no cluster.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()

        # Fetch all proteins and their domains in the same clusters as the input protein_ids
        query = f"""
        SELECT 
            p.proteinID, p.clusterID, d.domain
        FROM Proteins p
        LEFT JOIN Domains d ON p.proteinID = d.proteinID
        WHERE p.proteinID IN ({','.join('?' * len(protein_ids))})
           OR p.clusterID IN (
                SELECT DISTINCT clusterID FROM Proteins WHERE proteinID IN ({','.join('?' * len(protein_ids))})
           )
        """
        
        # Execute the query with the protein_ids list
        cur.execute(query, protein_ids + protein_ids)  # Pass protein_ids twice for both conditions
        protein_results = cur.fetchall()
    
    # Organize results by cluster and by proteinID
    cluster_dict = {}
    protein_cluster_map = {}  # Map each protein to its cluster
    domain_map = {}           # Map each protein to its domain

    for protein_id, cluster_id, domain in protein_results:
        domain_str = domain if domain else "no_domain"  # Handle missing domain
        domain_map[protein_id] = domain_str

        if cluster_id is not None:  # Only process if a cluster exists
            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = set()
            cluster_dict[cluster_id].add(domain_str)
            protein_cluster_map[protein_id] = cluster_id
        else:
            # Handle proteins without clusters separately in the result later
            protein_cluster_map[protein_id] = None

    # Prepare the final dictionary with proteinID as the key and neighboring domains as a set
    neighbors_dict = {}
    for protein_id in protein_ids:
        cluster_id = protein_cluster_map.get(protein_id)
        
        # Handle proteins with a cluster
        if cluster_id:
            neighbors_dict[protein_id] = cluster_dict[cluster_id]
        else:
            # Handle proteins without a cluster (just return the domain of the protein itself)
            neighbors_dict[protein_id] = {domain_map[protein_id]}

    return neighbors_dict

    

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
            print(f" Tree for {domain} is missing")
            return None  # Skip if tree is missing

        # Get the monophyletic group
        lca, clade_identifier_set = find_lca_and_monophyly(proteinID_set, tree)

        if clade_identifier_set:
            monophyly[domain] = clade_identifier_set
        else:
            return {},{} #if no addition was found and set is already monophylum, return emtpy dicts

        
        return monophyly

    except Exception as e:
        print(f"Error processing {proteinID_set}: {e}")
        return {},{}
        
        
        
def get_last_common_ancestor_fasta(options, grouped, trees_dict, tree_prefix):
    """
    Main function to process grouped data in parallel using multiprocessing.
    """

    # Prepare arguments for the worker function
    worker_args = [
        (proteinID_set, grouped_key.replace("grpd0_", "", 1), trees_dict, tree_prefix)
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



def fetch_protein_superfamily_to_fasta(options, blast_table, domain_keyword_dict, directory, max_seqs, deviation=0.5):
    output_files = {}

    with sqlite3.connect(options.database_directory) as con:
        cur = con.cursor()

        # Retrieve the SQLite variable limit dynamically
        cur.execute("PRAGMA max_variable_number;")
        result = cur.fetchone()
        sqlimit = result[0] if result else 999  # Fallback to 999

        for domain, query_length in domain_keyword_dict.items():  # Outer loop
            # Write the results to the domain-specific FASTA file
            name = domain.split('_')[-1]
            filename = f"superfamily_{name}.faa"
            output_fasta_path = os.path.join(directory, filename)
            
            #Skip if file exists
            if os.path.isfile(output_fasta_path):
                print(f" Skipping. File {output_fasta_path} already exists.")
                continue
            else:
                print(f" {output_fasta_path} did not exist")
            
            
            # Get the proteinIDs set
            proteinIDs = get_domain_superfamily_proteinIDs(blast_table, domain, query_length, deviation)
            proteinIDs = {string.strip() for string in proteinIDs}

            # Adjust the identifiers from glob to match the database
            if options.glob_gff:
                proteinIDs = {f"{string.split('___', 1)[0]}-{string}" for string in proteinIDs}
            else:
                proteinIDs = {string.replace('___', '-', 1) for string in proteinIDs}

            # Convert proteinIDs set to a list for slicing
            proteinIDs = list(proteinIDs)

            # Flag to indicate if the inner loop was broken
            skip_domain = False

            # Prepare to fetch in chunks based on the SQLite limit
            proteins = []
            total_sequences = 0  # Counter to track the number of sequences fetched
            for i in range(0, len(proteinIDs), sqlimit):  # Inner loop
                chunk = proteinIDs[i:i + sqlimit]
                
                # Use a parameterized query for each chunk
                query = '''
                    SELECT DISTINCT P.proteinID, P.sequence
                    FROM Proteins P
                    WHERE P.proteinID IN ({})
                '''.format(','.join('?' * len(chunk)))
                
                cur.execute(query, tuple(chunk))
                fetched = cur.fetchall()

                proteins.extend(fetched)
                total_sequences += len(fetched)

                # Stop fetching if the sequence limit is reached, never stop if max_seqs is equal 0
                if max_seqs > 0 and total_sequences >= max_seqs:
                    proteins = proteins[:max_seqs]  # Trim to max allowed sequences
                    skip_domain = True  # Set flag to skip further processing
                    break  # Exit the inner loop

            # Skip processing and file writing for this domain if limit was reached
            if skip_domain:
                print(f"Warning: Skipping domain {domain} because max sequence limit was reached.")
                continue  # Skip the current iteration of the outer loop

            # Use a set to track already written proteinIDs
            written_proteins = set()
            if os.path.exists(output_fasta_path):
                with open(output_fasta_path, 'r') as fasta_file:
                    for line in fasta_file:
                        if line.startswith('>'):
                            protein_id = line[1:].strip().split()[0]
                            written_proteins.add(protein_id)

            # Append the fetched sequences
            with open(output_fasta_path, 'a') as fasta_file:  # 'w' to overwrite the file
                for proteinID, sequence in proteins:
                    if proteinID not in written_proteins:
                        fasta_file.write(f'>{proteinID}\n{sequence}\n')
                        written_proteins.add(proteinID)

            # Add file to the output dictionary
            output_files[name] = output_fasta_path

    return output_files

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


