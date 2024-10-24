#!/usr/bin/python

import sqlite3
import os

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
    #domain_family_dict = get_domain_key_list_pairs(options.TP_singles, domain_family_dict)

    #fetch the proteins for each 
    fetch_protein_family_to_fasta(options, domain_family_dict)
    
    #align the protein family
    family_alignment_files = Alignment.initial_alignments(options, options.phylogeny_directory)
    
    #calculate phylogeny for each alignment file
    tree_files = Alignment.calculate_phylogeny_parallel(options, family_alignment_files) #creates .tree files in the phylogeny directory
    
    #Find the Last Common Ancestors (LCA) for each potential training set and write the fasta file
    monophylums = get_last_common_ancestor_fasta(options,training_datasets,tree_files)
    #get_last_common_ancestor_fasta(options,options.TP_singles,tree_files)
    
    options.TP_monophyla = monophylums

    
        
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

    
def fetch_protein_family_to_fasta(options, domain_keyword_dict):

    with sqlite3.connect(options.database_directory) as con:
        cur = con.cursor()

        # Iterate over the domain-keyword dictionary
        for domain, keywords in domain_keyword_dict.items():
            
            # Use a parameterized query to check for the domain and the keywords
            query = '''
                SELECT DISTINCT P.proteinID, P.sequence
                FROM Proteins P
                JOIN Domains D ON P.proteinID = D.proteinID
                JOIN Clusters C ON P.clusterID = C.clusterID
                JOIN Keywords K ON C.clusterID = K.clusterID
                WHERE D.domain = ? AND K.keyword IN ({})
            '''.format(','.join('?' * len(keywords)))

            # Execute the query with the domain and the list of keywords
            cur.execute(query, (domain, *keywords))

            proteins = cur.fetchall()

            # Use a set to track already written proteinIDs
            written_proteins = set()

            # Write the results to the domain-specific FASTA file
            # Create the filename for this domain's output file
            filename = f"{domain}.faa"
            output_fasta_path = os.path.join(options.phylogeny_directory, filename)
            with open(output_fasta_path, 'w') as fasta_file:
                for proteinID, sequence in proteins:
                    if proteinID not in written_proteins:
                        # Add protein to the file and mark it as written
                        fasta_file.write(f'>{proteinID}\n{sequence}\n')
                        written_proteins.add(proteinID)





def find_lca_and_monophyly(protein_ids, tree_file):
    """
    Given a phylogenetic tree and a set of protein identifiers, 
    this function finds the last common ancestor (LCA) of the given proteins,
    checks if the clade is monophyletic, and returns relevant information.
    
    Args:
    - tree: A Bio.Phylo tree object.
    - protein_ids: A list of protein identifiers (as present in the tree tips).
    
    Returns:
    - LCA: The last common ancestor of the given proteins.
    - is_monophylum: True if the clade is monophyletic, False otherwise.
    - included_identifiers: If not monophyletic, returns all included protein identifiers in the LCA clade.
    """
    tree = Phylo.read(tree_file, 'newick')
    
    # Find the LCA of the given protein identifiers
    lca = tree.common_ancestor(protein_ids)
    
    # Check if the clade is monophyletic
    is_monophylum = tree.is_monophyletic(protein_ids)
    
    if is_monophylum:
        return lca, set()  # LCA, monophylum status, and no need for included identifiers
    else:
        # Get all the identifiers in the LCA clade if it's not monophyletic
        included_identifiers = {tip.name for tip in lca.get_terminals()}
        return lca, included_identifiers


def get_last_common_ancestor_fasta(options,grouped,trees_dict):
    monophyly = {}
    for proteinID_frozenset, keyword_domain_pairs in grouped.items():
        
        #define domain, key and tree
        keyword, domain = keyword_domain_pairs[0]
        tree = trees_dict[domain]
        
        #get the monophyletic group
        lca, clade_identifier_set = find_lca_and_monophyly(proteinID_frozenset, tree)
        
        #check distance to the original group
        
        
        if clade_identifier_set:
            monophyly[frozenset(clade_identifier_set)] = [(keyword,domain)]

    return monophyly    







