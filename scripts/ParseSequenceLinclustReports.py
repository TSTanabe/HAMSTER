#!/usr/bin/python
import os

from . import Database
from . import myUtil

def cluster_sequences(options):
    #Groups sequences via linclust and updates the database
    
    #print(f"MMseqs linclust with --threads {options.cores} --min-seq-id {options.cminseqid} --alignment-mode {options.alignment_mode} -e {options.evalue} -c {options.clustercoverage}")    


    #Get all the hits for the sequences in the query file
    prefix = "Query_"
    
    initial_hit_files = [os.path.join(options.fasta_initial_hit_directory, f) for f in os.listdir(options.fasta_initial_hit_directory) if f.startswith(prefix)]
    limit = len(initial_hit_files)
    
    for index,fasta_file in enumerate(initial_hit_files):
        if not os.path.getsize(fasta_file) == 0:
            base_filename = os.path.basename(fasta_file)
            print(f"Sequence similarity clustering for {base_filename} ({index+1}/{limit})      ",end="\r")
            similarity_cluster_file = MMseqsLinclust(fasta_file,options)
            similarity_cluster_groups = parse_linclust_find_connected_groups(similarity_cluster_file) #a list of sets
            #print(similarity_cluster_groups)
            #get the name of the file
            fasta_file_name, extension = os.path.splitext(fasta_file)
            domain_type = fasta_file_name.split(prefix, 1)[1]
            
            #assigns the protein names: number_queryname
            iterator = 0
            for group in similarity_cluster_groups:
                iterator = iterator + 1
                similarity_cluster = str(iterator) + "_" + domain_type
                Database.update_domain(options.database_directory,group,domain_type,similarity_cluster)


def MMseqsLinclust(query,options):    
    
    mmseqs = myUtil.find_executable("mmseqs")
    #query_db_name="${query%.*}.queryDB" #can be done during argument parsing if possible
    query_db_name=f"{query}.queryDB"
    os.system(f'{mmseqs} createdb {query} {query_db_name} 1>/dev/null')

    output_results=f"{query}.output"
    output_results2=f"{query}_cluster.tsv"

    os.system(f'{mmseqs} cluster {query_db_name} {output_results} tmp --remove-tmp-files --threads {options.cores} --min-seq-id {options.cminseqid} --alignment-mode {options.alignment_mode} -e {options.evalue} -c {options.clustercoverage} 1>/dev/null') #-c is the minimal coverage
    os.system(f'{mmseqs} createtsv {query_db_name} {query_db_name} {output_results} {output_results2} --threads {options.cores} 1>/dev/null')

    #os.system(f'{mmseqs} rmdb {query_db_name}')
    return output_results2

######
#Parser for linclust results
######

def parse_linclust_find_connected_groups(table_file):
    # Build the graph
    graph = {}
    with open(table_file, 'r') as file:
        for line in file:
            identifier1, identifier2 = line.strip().split('\t')

            # Skip line if both columns are identical
            if identifier1 == identifier2:
                continue
            
            if identifier1 not in graph:
                graph[identifier1] = set()
            if identifier2 not in graph:
                graph[identifier2] = set()
            graph[identifier1].add(identifier2)
            graph[identifier2].add(identifier1)

    # Perform DFS to find connected groups
    visited = set()
    groups = []
    for identifier in graph:
        if identifier not in visited:
            group = set()
            dfs(graph, visited, identifier, group)
            groups.append(group)

    return groups

            
def dfs(graph, visited, identifier, group):
    visited.add(identifier)
    if not isinstance(group, set):
        group = set()  # Initialize group as a set if it's not already
    group.add(identifier)
    for neighbor in graph[identifier]:
        if neighbor not in visited:
            dfs(graph, visited, neighbor, group)




