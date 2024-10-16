#!/usr/bin/python
import os
from . import myUtil

def HMMsearch(Path,HMMLibrary,cores = 1):
    #Path to faa File, HMMLibrary for the HMMs and cores number of cores used by HMMER3
    #29.8.22

    Output = myUtil.getReportName(Path)
    #print("hmmsearch -E 0.0001 --cpu "+str(cores)+" "+HMMLibrary+" "+Path+">"+Output)
    os.system(f'hmmsearch -E 0.0001 --cpu {str(cores)} {HMMLibrary} {Path}>{Output}')
    return Output


def MMseqsSearch(path,query_db_name,options,cores=1):
    #query is the static fasta file with the query sequences
    #path is the assembly fasta file
    evalue = options.evalue
    coverage = options.searchcoverage
    minseqid = options.minseqid
    
    
    target_db_name=f"{path}.targetDB"
    os.system(f'mmseqs createdb {path} {target_db_name} 1>/dev/null')
    tmp = f"{path}.tmp"
    
    output_results=f"{path}.alndb"
    output_results_tab=f"{path}.tab"

    os.system(f'mmseqs search {query_db_name} {target_db_name} {output_results} {tmp} --threads {cores} --alignment-mode {options.alignment_mode} --min-seq-id {options.minseqid} -e {options.evalue} -c {options.searchcoverage} 1>/dev/null') #attention cores is used here --min-seq-id {coverage}
    os.system(f'mmseqs convertalis {query_db_name} {target_db_name} {output_results} {output_results_tab} --format-mode 0 --threads {cores} 1>/dev/null')
    
    
    return output_results_tab

def MMseqsLinclust(query,options):    

    #query_db_name="${query%.*}.queryDB" #can be done during argument parsing if possible
    query_db_name=f"{query}.queryDB"
    os.system(f'mmseqs createdb {query} {query_db_name} 1>/dev/null')

    output_results=f"{query}.output"
    output_results2=f"{query}_cluster.tsv"

    os.system(f'mmseqs cluster {query_db_name} {output_results} tmp --threads {options.cores} --min-seq-id {options.cminseqid} --alignment-mode {options.alignment_mode} -e {options.evalue} -c {options.clustercoverage} 1>/dev/null') #-c is the minimal coverage
    os.system(f'mmseqs createtsv {query_db_name} {query_db_name} {output_results} {output_results2} --threads {options.cores} 1>/dev/null')

    return output_results2





def convert_diamond_to_mmseqs_format(diamond_tab_file, output_file):
    with open(diamond_tab_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            columns = line.strip().split('\t')
            query, target = columns[0], columns[1]
            percent_identity = float(columns[2])
            alignment_length = int(columns[3])
            mismatches = int(columns[4])
            gapopens = int(columns[5])
            q_start, q_end = int(columns[6]), int(columns[7])
            t_start, t_end = int(columns[8]), int(columns[9])
            evalue = float(columns[10])
            bits = float(columns[11])
            
            nident = int(percent_identity / 100 * alignment_length)
            q_len = q_end - q_start + 1
            t_len = t_end - t_start + 1
            qcov = alignment_length / q_len
            tcov = alignment_length / t_len
            
            output_line = f"{query}\t{target}\t{q_start}\t{q_end}\t{t_start}\t{t_end}\t{q_len}\t{t_len}\t{alignment_length}\t{nident}\t{mismatches}\t{gapopens}\t{qcov:.2f}\t{tcov:.2f}\t{evalue}\t{bits}\n"
            outfile.write(output_line)










