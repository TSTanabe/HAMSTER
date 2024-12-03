#!/usr/bin/python

import sqlite3
import os
import subprocess

from . import Alignment

from Bio import Phylo

#first get the csb with their identifiers. make sets of appearence as value to key name of the csb
#do not compare the csb appearences as they jaccard itself should have managed it.


def csb_phylogeny_datasets(options):

    # options.phylogeny_directory #to save the protein family trees
    # options.TP_merged = merged #already merged because of highly similar training data key: set => value tuple(keyword,domain)
    # options.TP_singles = singles #not merged because unique dataset data key: set => value tuple(keyword,domain)

    #collect for each domain the relevant csb keywords
    training_datasets = {**options.TP_merged, **options.TP_singles}
    
    domain_family_dict = get_domain_key_list_pairs(training_datasets)

    #fetch the proteins for each 
    fetch_protein_type_to_fasta(options, domain_family_dict)
    
    query_length_dict = get_sequence_legth(options.self_query)
    fetch_protein_superfamily_to_fasta(options, options.glob_table, query_length_dict, options.phylogeny_directory) #for the phylogeny and as TN set

    #align the protein family
    family_alignment_files = Alignment.initial_alignments(options, options.phylogeny_directory)
    
    #calculate phylogeny for each alignment file
    tree_files = Alignment.calculate_phylogeny_parallel(options, family_alignment_files) #creates .tree files in the phylogeny directory
    
    #Find the Last Common Ancestors (LCA) for each potential training set and write the fasta file
    monophylums = get_last_common_ancestor_fasta(options,training_datasets,tree_files,'lca','l')
    
    superfamily_monophylums = get_last_common_ancestor_fasta(options,training_datasets,tree_files,'superfamily','u')
    
    options.TP_monophyla = monophylums
    options.superfamily = superfamily_monophylums

def csb_phylogeny_target_sets(options):
    training_datasets = {**options.TP_merged, **options.TP_singles}
    domain_family_dict = get_domain_key_list_pairs(training_datasets)
    query_length_dict = get_sequence_legth(options.self_query)
    target_files = fetch_protein_superfamily_to_fasta(options, options.glob_table, query_length_dict, options.cross_validation_directory,1) #for the phylogeny and as TN set
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


def get_last_common_ancestor_fasta(options, grouped, trees_dict, tree_prefix, key_prefix):
    monophyly = {}
    for proteinID_frozenset, keyword_domain_pairs in grouped.items():
        try:
            # Define domain, key, and tree
            keyword, domain = keyword_domain_pairs[0]
            domain = domain.split('_')[-1]
            tree_key = tree_prefix + '_' + domain
            
            # Attempt to get the tree from the dictionary
            tree = trees_dict[tree_key]  # Throws KeyError if tree_key is missing
            
            # Get the monophyletic group
            lca, clade_identifier_set = find_lca_and_monophyly(proteinID_frozenset, tree)
            
            # Check distance to the original group
            if clade_identifier_set:
                monophyly[frozenset(clade_identifier_set)] = [(key_prefix + keyword, domain)]
        
        except KeyError:
            # Handle the missing tree
            print(f"Warning: Tree for key '{tree_key}' not found. Skipping.")
            continue  # Skip this iteration and proceed with the next one

    return monophyly 




#########################################################################################################        
##################################### Hits sorted by query hit ##########################################
#########################################################################################################



def fetch_protein_superfamily_to_fasta(options, blast_table, domain_keyword_dict, directory, deviation=0.5):
    output_files = {}
    max_seqs = options.max_seqs  # Maximum allowed sequences

    with sqlite3.connect(options.database_directory) as con:
        cur = con.cursor()

        # Retrieve the SQLite variable limit dynamically
        cur.execute("PRAGMA max_variable_number;")
        result = cur.fetchone()
        sqlimit = result[0] if result else 999  # Fallback to 999

        for domain, query_length in domain_keyword_dict.items():  # Outer loop
            # Flag to indicate if the inner loop was broken
            skip_domain = False

            # Get the proteinIDs set
            proteinIDs = get_domain_superfamily_proteinIDs(blast_table, domain, query_length, deviation)
            proteinIDs = {string.strip() for string in proteinIDs}

            # Adjust the identifiers from glob to match the database
            if options.glob_gff:
                proteinIDs = {f"{string.split('___', 1)[0]}-{string}" for string in proteinIDs}
            else:
                proteinIDs = {string.replace('___', '-', 1) for string in proteinIDs}

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

                # Stop fetching if the sequence limit is reached
                if total_sequences >= max_seqs:
                    proteins = proteins[:max_seqs]  # Trim to max allowed sequences
                    skip_domain = True  # Set flag to skip further processing
                    break  # Exit the inner loop

            # Skip processing and file writing for this domain if limit was reached
            if skip_domain:
                print(f"Warning: Skipping domain {domain} because max sequence limit was reached.")
                continue  # Skip the current iteration of the outer loop

            # Write the results to the domain-specific FASTA file
            name = domain.split('_')[-1]
            filename = f"superfamily_{name}.faa"
            output_fasta_path = os.path.join(directory, filename)

            # Use a set to track already written proteinIDs
            written_proteins = set()
            if os.path.exists(output_fasta_path):
                with open(output_fasta_path, 'r') as fasta_file:
                    for line in fasta_file:
                        if line.startswith('>'):
                            protein_id = line[1:].strip().split()[0]
                            written_proteins.add(protein_id)

            # Append the fetched sequences
            with open(output_fasta_path, 'w') as fasta_file:  # 'w' to overwrite the file
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


