#!/usr/bin/python
from . import myUtil
#import re
#import os


def makeThresholdDict(File,threshold_type=1):
    """
        30.3.23
        Threshold_type 1 => optimized
        Threshold_type 2 => noise cutoff
        Threshold_type 3 => trusted_cutoff
    """
    Thresholds = {}
    with open(File, "r") as reader:
        for line in reader.readlines():
            #print(line)
            #lines = "key1=value1;key2=value2;key3=value3"
            l = line.split("\t")

            try:
                Thresholds[l[0]] = float(l[threshold_type])
            except:
                print("Error in cutoff line: "+line)
        #for k, v in Thresholds.items():
        #   print(k, v)
    return Thresholds

def HMMsearch(Path,HMMLibrary,cores = 1):
    #Path to faa File, HMMLibrary for the HMMs and cores number of cores used by HMMER3
    #29.8.22

    Output = myUtil.getReportName(Path)
    #print("hmmsearch -E 0.0001 --cpu "+str(cores)+" "+HMMLibrary+" "+Path+">"+Output)
    myUtil.command(f'hmmsearch -E 0.0001 --cpu {str(cores)} {HMMLibrary} {Path}>{Output}')
    return Output

def Diamondsearch(path,query,cores,evalue):
    
    database_name="${path%.*}_database"
    diamond_results="${path%.*}_diamondresults"
    myUtil.command(f'diamond makedb --in {path} --db {database_name}')
    myUtil.command(f'diamond blastp -d {database_name} -q {query} -o {diamond_results} -p {cores} -e {evalue}')
    unlink(database_name)

    return diamond_results

def MMseqsSearch(path,query,options):
    #query is the static fasta file with the query sequences
    #path is the assembly fasta file
    cores = options.cores
    evalue = options.evalue
    coverage = options.searchcoverage
    minseqid = options.minseqid
    
    query_db_name=f"{query}.queryDB" #can be done during argument parsing if possible
    #myUtil.command(f'mmseqs createdb {query} {query_db_name} 1>/dev/null')
    
    target_db_name=f"{path}.targetDB"
    myUtil.command(f'mmseqs createdb {path} {target_db_name} 1>/dev/null')
    tmp = f"{path}.tmp"
    
    output_results=f"{path}.alndb"
    output_results2=f"{path}.tab"

    myUtil.command(f'mmseqs search {query_db_name} {target_db_name} {output_results} {tmp} --threads {options.cores} --alignment-mode {options.alignment_mode} --min-seq-id {options.minseqid} -e {options.evalue} -c {options.searchcoverage} 1>/dev/null') #attention cores is used here --min-seq-id {coverage}
    myUtil.command(f'mmseqs convertalis {query_db_name} {target_db_name} {output_results} {output_results2} --format-mode 0 --threads {options.cores} 1>/dev/null')

    return output_results2

def MMseqsLinclust(query,options):    

    #query_db_name="${query%.*}.queryDB" #can be done during argument parsing if possible
    query_db_name=f"{query}.queryDB"
    myUtil.command(f'mmseqs createdb {query} {query_db_name} 1>/dev/null')

    output_results=f"{query}.output"
    output_results2=f"{query}_cluster.tsv"

    myUtil.command(f'mmseqs cluster {query_db_name} {output_results} tmp --threads {options.cores} --min-seq-id {options.cminseqid} --alignment-mode {options.alignment_mode} -e {options.evalue} -c {options.clustercoverage} 1>/dev/null') #-c is the minimal coverage
    myUtil.command(f'mmseqs createtsv {query_db_name} {query_db_name} {output_results} {output_results2} --threads {options.cores} 1>/dev/null')

    return output_results2


#def getGenomeID(Path):
    #should return genomeID according to standard genome identifiers from common source
    #29.8.22
    #9.10.22 deprecated, moved to my Util
#    File = myUtil.getFileName(Path)
#    File = myUtil.removeExtension2(File)
#    return File








