#!/usr/bin/python

import os
import sys
import argparse

from . import Project
from . import Queue

from . import Database
from . import Translation
from . import ParallelSearch
from . import ParseSequenceLinclustReports

from . import Csb_cluster
from . import Csb_proteins

from . import Alignment

from . import Validation
from . import Reports
from . import myUtil


if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))



class Options:
    def __init__(self):
    

        
        #variables used by other subroutines, possibly they can be moved somewhere else
        #in order of appearence
        self.queued_genomes = set()
        self.faa_files = {}
        self.gff_files = {}
        
        
        #csb prediction dereplication
        self.redundant = 0
        self.non_redundant = 0
        self.redundancy_hash = dict()
        
        self.csb_overlap_factor = 0.75 #merge subsets with at least this overlap
        self.gap_remove_threshold = 0.95 #remove columns with gaps with trimal

        self.sequence_faa_file = None #filepath to all sequences fasta file
        self.reports = dict()
        self.relaxed_reports = dict()
        self.standard_cutoff_report_file = None
        self.standard_performance_report_file = None
        self.relaxed_cutoff_report_file = None
        self.relaxed_performance_report_file = None
        
        self.csb_name_prefix = "csb-" #prefix of clusterIDs determined by csb finder algorithm
        self.csb_name_suffix = "_" #suffix of clusterIDs determined by csb finder algorithm
        
        self.self_query = None


def parse_arguments(arguments):
    formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=96,width =300)
    
    parser = argparse.ArgumentParser(formatter_class=argparse.HelpFormatter, description = "HAMSTER version 0.0.6 \nSyntax: HAMSTER [OPTIONS]",epilog = "")
    
    
    parser.add_argument('-f', dest='fasta_file_directory', type=myUtil.dir_path, default = None, help='Directory of the target fasta files')
    parser.add_argument('-q', dest='query_file', type=myUtil.file_path, default = None, help='Query sequences fasta file')

    resources = parser.add_argument_group("Optional parameters")
    resources.add_argument('-s', dest='stage', type=int, default = 0, choices= [0,1,2,3,4,5,6,7,8,9],help='Start at stage')
    resources.add_argument('-c', dest='cores', type=int, default = 2, help='Number of CPUs')
    resources.add_argument('-t', dest='taxonomy_file', type=myUtil.file_path, default = None, help='Taxonomy csv file')
    resources.add_argument('-r', dest='result_files_directory', type=myUtil.dir_path, default = __location__+"/results", help='Directory for the result files/results from a previous run') # project folder TODO print the directory that is used for reconfirmation with the user
    resources.add_argument('-db',dest='database_directory', type=myUtil.file_path, metavar='<filepath>', help='Filepath to existing database')
    
    resources.add_argument('-glob_faa', dest='glob_faa', type=myUtil.file_path, default=None, help='Concatenated fasta file')
    resources.add_argument('-glob_gff', dest='glob_gff', type=myUtil.file_path, default=None, help='Concatenated gff file')
    resources.add_argument('-glob_chunks', dest='glob_chunks', type=int, default=3000, metavar='<int>',help='Chunk size for parsing results from glob before entering into database')
    resources.add_argument('-no_glob', dest='glob_search', action='store_false', help='Do not concatenated fasta file for search')
    
        
    search = parser.add_argument_group("Optional search parameters and sequence clustering parameters for mmseqs2")
    search.add_argument('-evalue', dest='evalue', type=float, default = 0.1, help='E-value cutoff [0,inf]')
    search.add_argument('-thrs_score', dest='thrs_score', type=int, default = 10, help='Score cutoff [0,inf]')
    search.add_argument('--min-seq-id',dest='minseqid',type=float, default=0.000, help='Sequence search matches above this sequence identity [0.0,1.0]')
    search.add_argument('--search-coverage', dest='searchcoverage', type=float, default=0.000, help='Min. coverage used for searching')
    
    #Linclust parameters
    search.add_argument('--alignment-mode',dest='alignment_mode',type=int, default=2, choices=[0,1,2,3,4], help='mmseqs2 linclust search alignment mode')
    search.add_argument('--cluster-coverage', dest='clustercoverage', type=float, default = 0.000, help='mmseqs2 linclust min. coverage used for clustering sequences')
    search.add_argument('--cluster:min-seq-id',dest='cminseqid',type=float, default=0.000, help='mmseqs2 search list matches above this sequence identity [0.0,1.0]')


    
    genecluster = parser.add_argument_group("Optional gene cluster parameters")
    genecluster.add_argument('--distance', dest='nucleotide_range', type=int, default = 3500, help='Max. nucleotide distance between synthenic genes')
    genecluster.add_argument('-p', dest= 'patterns_file' , type=myUtil.file_path, default=__location__+"/src/Patterns", metavar='<filepath>', help='Filepath to patterns file')
    
    csb = parser.add_argument_group("Optional collinear synthenic block parameters")
    csb.add_argument('-insertions', dest='insertions', type=int,default = 2, help='Max. insertions in a csb')
    csb.add_argument('-occurence', dest='occurence', type=int,default = 2, help='Min. number of csb occurs at least times')
    csb.add_argument('-min_csb_size', dest='min_csb_size', type=int,default = 2, help='Min. csb size before recognized as csb')
    csb.add_argument('-jaccard', dest='jaccard', type=float,default = 0.4, help='Acceptable dissimilarity in jaccard clustering. 0.2 means that 80 percent have to be the same genes')
    csb.add_argument('-csb_overlap', dest='csb_overlap_factor', type=float,default = 0.75, help='Merge if sequences from two csb is identical above this threshold')
    
    csb.add_argument('-csb_distinct', dest='csb_distinct_grouping', action='store_true', help='Set to merge only if sequences from two csb are distinct.')
    
    

    alignment = parser.add_argument_group("Optional alignment parameters")
    alignment.add_argument('--min_seqs', dest='min_seqs', type=int, default = 5, help='Min. number of required sequences for the alignment')
    alignment.add_argument('--gap_col_remove', dest='gap_remove_threshold', type=int, default = 95, help='[0,100] remove alignment columns with percent gaps')
    
    #hhfilter = parser.add_argument_group("Optional HHfilter parameters")
    #hhfilter.add_argument('--hhfilter_id', dest='hhfilter_id', type=int, default = 90, choices = [], help='[0,100]  maximum pairwise percent sequence identity def=90')
    #hhfilter.add_argument('--hhfilter_diff', dest='hhfilter_diff', type=int, default = 0, help='[0,inf[  filter MSA by selecting most diverse set of sequences, keeping at least this many seqs in each MSA block of length 50 def=0 ')
    #hhfilter.add_argument('--hhfilter_cov', dest='hhfilter_cov', type=int, default = 0, choices = [], help='[0,100]  minimum percent coverage with query def=0')
    #hhfilter.add_argument('--hhfilter_qid', dest='hhfilter_qid', type=int, default = 0, choices = [], help='[0,100]  minimum percent sequence identity with query def=0')
    #hhfilter.add_argument('--hhfilter_maxseq', dest='hhfilter_maxseq', type=int, default = 65535, help='<int>  max number of input rows def=65535') 
    #hhfilter.add_argument('--hhfilter_maxres', dest='hhfilter_maxres', type=int, default = 20001, help='<int>  max number of HMM columns def=20001') 





    options = Options()
    parser.parse_args(namespace=options)
    
    if options.fasta_file_directory is None or options.query_file is None:  
        sys.exit("ERROR: Missing query file or target genomes. Please use -f and -q arguments")
    
    return options
    

def fasta_preparation(options):
    #Unpacks and translates fasta
    Translation.parallel_translation(options.fasta_file_directory, options.cores)
    Translation.parallel_transcription(options.fasta_file_directory, options.cores)
    
    if not options.glob_search:
        print("Not concatenating protein fasta files")
        Queue.queue_files(options)
        return
    
    #concats fasta if wished
    if options.glob_faa and options.glob_gff:
        print(f"Preparing separate files from {options.glob_faa}")
        #create separated files
        Translation.deconcat(options)
        Queue.queue_files(options)
    else:
        #concat to to globfile
        Queue.queue_files(options)
        print(f"Generating glob file")
        Translation.create_glob_file(options) #fasta_file_directory, options.cores, concat the files with the genomeIdentifier+ ___ + proteinIdentifier   
    Translation.create_selfquery_file(options)

def initial_search(options):
    #Writes a database with the protein hits and their gene clusters for later use.
    #Writes fasta files for each hit for linclust    
    if os.path.isfile(options.database_directory):
        genomeIDs = Database.fetch_genomeIDs(options.database_directory)
        Queue.compare_with_existing_database(options,genomeIDs)
    else:
        Database.create_database(options.database_directory)


    if options.glob_search:
        #Search the glob file
        ParallelSearch.initial_glob_search(options)
    else:
        #Process all files separately
        ParallelSearch.initial_genomize_search(options)
    ParallelSearch.self_blast_query(options)
    
    
    
def cluster_sequences(options):
    #Groups sequences via linclust and updates the database with grouped identifiers
    
    print(f"MMseqs linclust with --threads {options.cores} --min-seq-id {options.cminseqid} --alignment-mode {options.alignment_mode} -e {options.evalue} -c {options.clustercoverage}")    
    ParseSequenceLinclustReports.cluster_sequences(options)


def csb_finder(options):

    Csb_cluster.csb_prediction(options)
    csb_gene_cluster_dict = Csb_cluster.csb_jaccard(options)
    #print(csb_gene_cluster_dict)
    Database.delete_keywords_from_csb(options.database_directory, options) #remove keys with options.csb_name_prefix options.csb_name_suffix to avoid old keyword interference
    Database.update_keywords(options.database_directory,csb_gene_cluster_dict) #assigns the names of the keywords to the clusters



def generate_csb_sequence_fasta(options):
    #prepares the sequence fasta files for the alignments
    
    Csb_proteins.csb_proteins_fasta(options)
    
    return
    

def model_alignment(options):
    Alignment.initial_alignments(options)


def cross_validation(options):

    Validation.create_hmms_from_msas(options.fasta_output_directory,"fasta_aln","hmm",options.cores) #create the full hmms for later use
    Reports.move_HMMs(options.fasta_output_directory,options.Hidden_markov_model_directory,"hmm") #move the hmms to the Hidden markov model folder
    Validation.parallel_cross_validation(options)


       
def report_cv_performance(options):
    print(f"Saving the HMMs to {options.Hidden_markov_model_directory}")
    options.reports = Reports.parse_matrices_to_report(options.cross_validation_directory,"_cv_matrices.txt")
    
    print(f"Saving the performance values to {options.Hidden_markov_model_directory}")
    options.standard_performance_report_file = Reports.create_performance_file(options,options.reports,options.Hidden_markov_model_directory,"/cutoff_performance.txt")
    
    print(f"Saving the cutoff values to {options.Hidden_markov_model_directory}")
    options.strict_cutoff_report_file = Reports.concatenate_cv_cutoff_files(options.cross_validation_directory, "_cv_thresholds.txt", options.Hidden_markov_model_directory+"/strict_cutoffs.txt")





    
    
            



    

    
    
  
def main(args=None):
#1
    options = parse_arguments(args)
    myUtil.print_header("\nPreparing space for the results")
    Project.prepare_result_space(options)
    
    
    if options.stage <= 1:
        myUtil.print_header("\nProkaryotic gene recognition and translation via prodigal")
        fasta_preparation(options)
        

    
#2    
    if options.stage <= 2:
        myUtil.print_header("\nSearching for homologoues sequences")
        initial_search(options)
        
#3    
    #cluster the adjacent protein sequences, which are not annotated yet into groups called pcs combined with a number
    #cluster with mmseqs the divergent_output_file
    #Info: this resulted in an error when only a single genome was used. linclust did this error and dumped the core. the tsv file was not created and subsequent errors occured therefore

    if options.stage <= 3:
        myUtil.print_header("\nClustering sequences by similarity")
        cluster_sequences(options)
    
#4    
    #csb naming
    if options.stage <= 4:
        myUtil.print_header("\nSearching for collinear syntenic blocks")
        csb_finder(options)
#5
    if options.stage <= 5:
        generate_csb_sequence_fasta(options)   

#6    
    #align
    if options.stage <= 6:
        myUtil.print_header("\nAligning sequences")
        model_alignment(options)
        
    
#7        
    #make cross validation files
    #Validation.CV(options.fasta_output_directory,options.cores)
    if options.stage <= 7:
        myUtil.print_header("\nPerforming cross valdation procedure")
        cross_validation(options)
        #report_cv_performance(options)
        

#8  
    #mach eine weitere HMMsearch mit dem vollen model auf alle seqs und prüfe die treffer
    #die FP Sequenzen ausgeben
    #die FN Sequenzen ausgeben
    if options.stage <= 8:
        report_cv_performance(options)
    
    
    

if __name__ == "__main__":
    args = sys.argv[1:]
    
    main(args)

