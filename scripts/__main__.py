#!/usr/bin/python

import os
import sys
import argparse
import pprint
#from datetime import datetime

from . import Project
from . import Queue

from . import Database
from . import Translation
from . import ParallelSearch
from . import ParseSequenceLinclustReports

from . import Csb_cluster
from . import Csb_proteins
from . import Csb_phylogeny
from . import Csb_statistic
from . import Csb_mcl

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
        

        self.sequence_faa_file = None #dictionary to the target files for the validation
        self.reports = dict()
        self.standard_cutoff_report_file = None
        self.standard_performance_report_file = None
        
        self.csb_name_prefix = "csb-" #prefix of clusterIDs determined by csb finder algorithm
        self.csb_name_suffix = "_" #suffix of clusterIDs determined by csb finder algorithm
        
        self.self_query = None #Query fasta like concatenated fasta
        self.self_seqs = None #Query sequences

        self.MCC_threshold = 0.8 #TODO move this to options, below this threshold the HMMs are not validated to save runtime
        
        self.deconcat_flag = 0 #Memorize if a deconcatenation was done
        self.glob_flag = 0 #Memorize if a concatenation was done
        
        
def parse_arguments(arguments):
    formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=96,width =300)
    
    parser = argparse.ArgumentParser(formatter_class=argparse.HelpFormatter, description = "HAMSTER version 0.0.8 \nSyntax: HAMSTER [OPTIONS]",epilog = "")
    
    
    parser.add_argument('-f', dest='fasta_file_directory', type=myUtil.dir_path, default = None, metavar = '<directory>', help='Directory of the target fasta files')
    parser.add_argument('-q', dest='query_file', type=myUtil.file_path, default = None, metavar = '<filepath>', help='Query sequences fasta file')

    resources = parser.add_argument_group("Optional parameters")
    resources.add_argument('-s', dest='stage', type=int, default = 0, choices= [0,1,2,3,4,5,6,7,8,9], metavar='<int>', help='Start at stage')
    resources.add_argument('-x', dest='end', type=int, default = 9, choices= [0,1,2,3,4,5,6,7,8,9], metavar='<int>', help='End at stage')
    resources.add_argument('-c', dest='cores', type=int, default = 2, metavar='<int>', help='Number of CPUs')
    resources.add_argument('-t', dest='taxonomy_file', type=myUtil.file_path, default = None, metavar = '<filepath>', help='Taxonomy csv file')
    resources.add_argument('-r', dest='result_files_directory', type=myUtil.dir_path, default = __location__+"/results", metavar = '<directory>', help='Directory for the result files/results from a previous run') # project folder TODO print the directory that is used for reconfirmation with the user
    resources.add_argument('-db',dest='database_directory', type=myUtil.file_path, metavar='<filepath>', help='Filepath to existing database')
   
    
    resources.add_argument('-glob_faa', dest='glob_faa', type=myUtil.file_path, default=None, metavar = '<filepath>', help='Predefined concatenated fasta file')
    resources.add_argument('-glob_gff', dest='glob_gff', type=myUtil.file_path, default=None, metavar = '<filepath>', help='Concatenated gff file')
    resources.add_argument('-glob_blast_table', dest='glob_table', type=myUtil.file_path, default=None, metavar = '<filepath>', help='Concatenated blast result table')
    resources.add_argument('-glob_chunks', dest='glob_chunks', type=int, default=3000, metavar='<int>', help='Chunk size for parsing results from glob before entering into database')
    resources.add_argument('-no_glob', dest='glob_search', action='store_false', help='Do not concatenated fasta file for search')
    resources.add_argument('-cv-off', dest='cross_validation_deactivated', action='store_true', help='Skip cross-validation step')
    resources.add_argument('-rep-off', dest='hit_report_deactivated', action='store_true', help='Skip hit report')
    
    
    search = parser.add_argument_group("Optional search parameters for diamond")
    search.add_argument('-evalue', dest='evalue', type=float, default = 1e-5, metavar = '<float>', help='E-value cutoff [0,inf]. Default: 1e-5')
    search.add_argument('-thrs_score', dest='thrs_score', type=int, default = 100, metavar = '<int>', help='Score cutoff [0,inf]. Default: 100')
    search.add_argument('-min-seq-id',dest='minseqid',type=float, default=25, metavar = '<float>', help='Sequence search matches above this sequence identity [0,100.0]. Default: 25')
    search.add_argument('-search-coverage', dest='searchcoverage', type=float, default=0.6, metavar = '<float>', help='Min. coverage used for searching [0.0,1.0]. Default: 0.6')
    search.add_argument('-alignment-mode',dest='alignment_mode',type=int, default=2, metavar='<int>', choices=[0,1,2], help='DIAMOND BLASTp alignment mode. Default: 2')
    search.add_argument('-memory',dest='diamond_memory_limit',type=int, default=8, metavar='<int>', help='DIAMOND BLASTp RAM limit')
    search.add_argument('-blast-score-ratio', dest='thrs_bsr', type=float, default=0.0, metavar = '<float>', help='Blast score ratio for hits [0.0,1.0]. Default: 0.0')
    search.add_argument('-allow_multidomain', dest='multidomain_allowed', action='store_true', help='Allow multiple query hits for each sequence')
    search.add_argument('-reports_hit', dest='diamond_report_hits_limit', type=int, default=0, metavar = '<int>', help='Limit to this number of top hits per query. 0 = no limit')



    #Cluster parameters
    #protein_cluster = parser.add_argument_group("Sequence clustering parameters for mmseqs2")
    #protein_cluster.add_argument('-cluster-active', dest='protein_cluster_active', action='store_true', help='Cluster initial blastp hits with mmseqs2 cluster')
    #protein_cluster.add_argument('-alignment-mode',dest='alignment_mode',type=int, default=2, metavar='<int>', choices=[0,1,2,3,4], help='mmseqs2 cluster search alignment mode')
    #protein_cluster.add_argument('-cluster-coverage', dest='clustercoverage', type=float, default = 0.800, metavar='<float>', help='mmseqs2 cluster min. coverage used for clustering sequences')
    #protein_cluster.add_argument('-cluster-min-seq-id',dest='cminseqid',type=float, default=0.000, metavar='<float>', help='mmseqs2 search list matches above this sequence identity [0.0,1.0]')


    
    genecluster = parser.add_argument_group("Optional gene cluster prediction parameters")
    genecluster.add_argument('-distance', dest='nucleotide_range', type=int, default = 3500, metavar='<int>', help='Max. nucleotide distance between synthenic genes. Default: 3500')
    genecluster.add_argument('-p', dest= 'patterns_file' , type=myUtil.file_path, default=__location__+"/src/Patterns", metavar='<filepath>', help='Filepath to patterns file')
    
    csb = parser.add_argument_group("Optional collinear synthenic block parameters")
    csb.add_argument('-insertions', dest='insertions', type=int,default = 2, metavar='<int>', help='Max. insertions in a csb. Default: 2')
    csb.add_argument('-occurence', dest='occurence', type=int,default = 10, metavar='<int>', help='Min. number of csb occurs at least times. Default: 10')
    csb.add_argument('-min_csb_size', dest='min_csb_size', type=int,default = 4, metavar='<int>', help='Min. csb size before recognized as csb. Default: 4')
    csb.add_argument('-jaccard', dest='jaccard', type=float,default = 0.4, metavar='<float>', help='Acceptable dissimilarity in jaccard clustering. 0.2 means that 80 percent have to be the same genes. Default: 0.4')
    csb.add_argument('-csb_overlap', dest='csb_overlap_factor', type=float, default = 0.75, metavar='<float>', help='Merge if sequences from two csb is identical above this threshold. Default: 0.75')
    
    csb.add_argument('-no_phylogeny', dest='csb_distinct_grouping', action='store_false', help='Skip phylogenetic supported training dataset clustering')
    csb.add_argument('-no_mcl', dest='csb_mcl_clustering', action='store_false', help='Skip markov chain clustering')
    csb.add_argument('-scan_eps', dest='dbscan_epsilon', type=float,default = 0.3, metavar='<float>', help='Acceptable dissimilarity for protein training datasets to be clustered. Default: 0.3') #Wir das noch genutzt?
    csb.add_argument('-distant_homologs', dest='sglr', action='store_true', help='Include alignments for distantly related proteins with conserved genomic vicinity')

    mcl_search = parser.add_argument_group("Optional Markov Chain Clustering parameters")
    mcl_search.add_argument('-mcl_evalue', dest='mcl_evalue', type=float, default = 1e-10, metavar = '<float>', help='MCL matrix e-value cutoff [0,inf]. Default: 1e-10')
    mcl_search.add_argument('-mcl_min-seq-id',dest='mcl_minseqid',type=float, default=25, metavar = '<float>', help='MCL matrix sequence identity cutoff [0,100.0]. Default: 25')
    mcl_search.add_argument('-mcl_search-coverage', dest='mcl_searchcoverage', type=float, default=0.6, metavar = '<float>', help='MCL matrix min. coverage [0.0,1.0]. Default: 0.6')
    mcl_search.add_argument('-mcl_hit_limit', dest='mcl_hit_limit', type=int, default=500, metavar = '<int>', help='MCL maximum number of edges between sequences. Default: 500')
    mcl_search.add_argument('-mcl_inflation', dest='mcl_inflation', type=float, default=2.0, metavar = '<float>', help='MCL inflation factor for granularity control. Default: 2.0')
    mcl_search.add_argument('-mcl_sensitivity', dest='mcl_sensitivity', type=str, choices=["fast", "more-sensitive", "sensitive", "very-sensitive", "ultra-sensitive"], default="more-sensitive", metavar='<sensitivity>', help="Set DIAMOND sensitivity for MCL clustering. Choices: fast, more-sensitive, sensitive, very-sensitive, ultra-sensitive. Default: more-sensitive")
    mcl_search.add_argument('-mcl_density_thrs', dest='mcl_density_thrs', type=float, default=0.01, metavar = '<float>', help='Required proportion of reference sequences in the total number of sequences in the MCL cluster to label it as true positive. Default: 0.01')
    mcl_search.add_argument('-mcl_reference_thrs', dest='mcl_reference_thrs', type=float, default=0.1, metavar = '<float>', help='Required proportion of reference sequences from the total reference sequences in the MCL cluster to label it as true positive. Default: 0.1')
    
    #            options.mcl_hit_limit, options.alignment_mode, options.mcl_sensitivity
    alignment = parser.add_argument_group("Optional alignment parameters")
    alignment.add_argument('-min_seqs', dest='min_seqs', type=int, default = 5, metavar='<int>', help='Min. number of required sequences for the alignment. Default: 5')
    alignment.add_argument('-max_seqs', dest='max_seqs', type=int, default = 100000, metavar='<int>', help='Max. number of sequences that are aligned. Default: 100 000')
    alignment.add_argument('-gap_col_remove', dest='gap_remove_threshold', type=float, default = 0.05, metavar='<float>', help='[0,1] remove alignment columns with only percent amino acids. Default: 0.05')
    alignment.add_argument('-include_domains', dest='include_list',nargs='+', default=[], metavar='<list>', help='List of domains, separated by spaces, that are specifically included')
    alignment.add_argument('-exclude_domains', dest='exclude_list', nargs='+', default=[], metavar='<list>', help='List of domains, separated by spaces that are specifically excluded')
    
    
    if len(arguments) == 0: # Print help if no arguments were provided
        parser.print_help()
        sys.exit("No arguments were provided.")

    options = Options()
    parser.parse_args(namespace=options)
    
    validate_options(options)
    
    return options
    
def validate_options(options):
    # Check if a database file exists
    database_provided = options.database_directory and os.path.isfile(options.database_directory)
    
    # Check if a result directory exists
    result_directory_provided = options.result_files_directory and os.path.isdir(options.result_files_directory)
    
    # Check if both fasta file directory and query file exist
    fasta_and_query_provided = (
        options.fasta_file_directory and 
        options.query_file and 
        os.path.isdir(options.fasta_file_directory) and 
        os.path.isfile(options.query_file)
    )
    
    # Only proceed if one of the conditions is met
    if not (database_provided or result_directory_provided or fasta_and_query_provided):
        sys.exit("ERROR: Please provide either a database file, a result directory, or both a fasta file directory and a query file using the -f and -q arguments.")
    

def fasta_preparation(options):
    
   
    if not options.glob_search:
        print("Not concatenating protein fasta files")
        Queue.queue_files(options)
        return
    
    #Comcatenate the fasta files, in case of a glob fasta file is provided deconcat it
    if options.glob_faa and options.glob_gff:
        print(f"Preparing separate files from {options.glob_faa}")
        #create separated files
        Translation.deconcat(options)
        Queue.queue_files(options)
        options.deconcat_flag = 1
    else:

        #Unpacks and translates fasta
        Translation.parallel_translation(options.fasta_file_directory, options.cores)
        Translation.parallel_transcription(options.fasta_file_directory, options.cores)

        #concat to to globfile if search results are not already provided
        Queue.queue_files(options)
        if not options.glob_table:
            print(f"Generating glob faa file")
            Translation.create_glob_file(options) #fasta_file_directory, options.cores, concat the files with the genomeIdentifier+ ___ + proteinIdentifier
            options.glob_flag = 1   
    Translation.create_selfquery_file(options)

def initial_search(options):
    #Writes a database with the protein hits and their gene clusters for later use.
    #Writes fasta files for each hit for linclust    

    if not os.path.isfile(options.database_directory):
        Database.create_database(options.database_directory)

    if options.glob_search:
        #Search the glob file
        ParallelSearch.initial_glob_search(options)
    else:
        #Process all files separately
        ParallelSearch.initial_genomize_search(options)

    if options.deconcat_flag: #deconcatenation was done but is not needed anymore
        myUtil.remove_directory(options.fasta_file_directory)#remove the deconcatenated files
    elif options.glob_flag:
        os.remove(options.glob_faa)
    
def cluster_sequences(options):
    # Currently not in use
    if not options.protein_cluster_active: # skip if not activated to cluster the sequences
        return
    
    #Groups sequences via linclust and updates the database with grouped identifiers
    Database.index_database(options.database_directory)
    
    print(f"MMseqs cluster with --threads {options.cores} --min-seq-id {options.cminseqid} --alignment-mode {options.alignment_mode} -e {options.evalue} -c {options.clustercoverage}")    
    ParseSequenceLinclustReports.cluster_sequences(options)


def csb_finder(options):

    Csb_cluster.csb_prediction(options)
    csb_gene_cluster_dict = Csb_cluster.csb_jaccard(options)

    Database.index_database(options.database_directory)
    Database.delete_keywords_from_csb(options.database_directory, options) #remove keys with options.csb_name_prefix options.csb_name_suffix to avoid old keyword interference
    Database.update_keywords(options.database_directory,csb_gene_cluster_dict) #assigns the names of the keywords to the clusters



def generate_csb_sequence_fasta(options):
    #prepares the sequence fasta files for the alignments
    Database.index_database(options.database_directory)
    
    score_limit_dict, filtered_stats_dict, grouped_keywords_dict, clustered_excluded_keywords_dict = Csb_statistic.group_gene_cluster_statistic(options) #keywords are in list of lists with the
    print("Fetching sequences for training datasets")
    csb_proteins_dict = Csb_proteins.csb_proteins_datasets(options, clustered_excluded_keywords_dict) # groups training data

    # Collect proteins from singular distantly related keywords
    if options.sglr:
        print("Processing dissimilar homologs with specific genomic context")
        singular = Csb_proteins.csb_proteins_datasets_combine(clustered_excluded_keywords_dict, csb_proteins_dict, "sglr")
        Csb_proteins.fetch_seqs_to_fasta_parallel(options.database_directory, singular, options.fasta_output_directory, options.min_seqs, options.max_seqs, options.cores)
        singular = {}
    
    
    # Collect proteins from grouped keywords
    print("Processing highly similar homologs with specific genomic context")
    grouped = Csb_proteins.csb_proteins_datasets_combine(grouped_keywords_dict, csb_proteins_dict, "grpd")
    grouped = Csb_proteins.add_query_ids_to_proteinIDset(grouped, options.database_directory)
    Csb_proteins.fetch_seqs_to_fasta_parallel(options.database_directory, grouped, options.fasta_output_directory, 10, options.max_seqs, options.cores)
	
    # Like before per TPs per csb, erscheint mir nicht sinnvoll, weil bisher waren diese ergebnisse immer etwas schlechter als die plcsb
    #Csb_proteins.csb_granular_datasets(options,csb_proteins_dict)
    if options.csb_distinct_grouping or options.csb_mcl_clustering:
        print("Including homologs without genomic context based on protein sequence phylogeny.")
        Csb_proteins.decorate_training_data(options, score_limit_dict, grouped)
        #decorated_grouped_dict = Csb_proteins.fetch_protein_ids_parallel(options.database_directory, score_limit_dict, options.cores)
        #decorated_grouped_dict = Csb_proteins.merge_grouped_protein_ids(decorated_grouped_dict, grouped)
        #Csb_proteins.fetch_seqs_to_fasta_parallel(options.database_directory, decorated_grouped_dict, options.phylogeny_directory, 20, options.max_seqs, options.cores)
        
        # Phylogeny clustering
        if options.csb_distinct_grouping:
            # Decorate grouped high scoring csb with similarly high seqs from the same phylogenetic clade
            print("Calculating phylogeny")
            decorated_grouped_dict = Csb_phylogeny.csb_phylogeny_datasets(options, grouped) # phylogenetic grouped training data
            Csb_proteins.fetch_seqs_to_fasta_parallel(options.database_directory, decorated_grouped_dict, options.fasta_output_directory, options.min_seqs, options.max_seqs, options.cores)
    
        # Markov chain clustering 
        if options.csb_mcl_clustering:
            # Eingabe daten grpd0
            Csb_mcl.csb_mcl_datasets(options,grouped) # markov chain clustering grouped training data


        
def model_alignment(options):
    Alignment.initial_alignments(options, options.fasta_output_directory)


def cross_validation(options):
    print("Initialize training data subsamples")
    Validation.create_hmms_from_msas(options.fasta_output_directory, options.Hidden_markov_model_directory, "fasta_aln","hmm",options.cores) #create the full hmms for later use
    Reports.move_HMMs(options.fasta_output_directory,options.Hidden_markov_model_directory,"hmm") #move the hmms to the Hidden markov model folder
    
    Csb_proteins.fetch_domains_superfamily_to_fasta(options, options.cross_validation_directory) #Target file for each HMM excluding seqs already below threshold
    Csb_proteins.fetch_all_proteins(options.database_directory, options.cross_validation_directory+"/sequences.faa") #Backup target file if something fails

    print("Initialize target sequence sets")
    options.sequence_faa_file = options.cross_validation_directory+"/sequences.faa" #File with all sequences to be searched
    options.targeted_sequence_faa_file_dict = Validation.get_target_sets(options.cross_validation_directory)
    print(f"From {options.cross_validation_directory} the following target files were selected {options.targeted_sequence_faa_file_dict}")
    
    Validation.parallel_cross_validation(options)


       
def report_cv_performance(options):
    #Initial validation
    print(f"Saving the cutoffs and performance reports from initial calculation to {options.Hidden_markov_model_directory}")
    Reports.concat_and_sort_files(options.fasta_alignment_directory, '_MCC.txt', options.Hidden_markov_model_directory, "_ini_performance_matrices.txt")
    Reports.concat_and_sort_files(options.fasta_alignment_directory, '_thresholds.txt', options.Hidden_markov_model_directory, "_ini_cutoffs.txt")
    

    #cross validation
    print(f"Saving the cutoffs and performance reports from the cross-validatio to {options.Hidden_markov_model_directory}")
    options.reports = Reports.parse_matrices_to_report(options.cross_validation_directory,"_cv_matrices.txt")
    
    options.standard_performance_report_file = Reports.create_performance_file(options,options.reports,options.Hidden_markov_model_directory,"/cutoff_performance.txt")
    
    options.strict_cutoff_report_file = Reports.concatenate_cv_cutoff_files(options.cross_validation_directory, "_cv_thresholds.txt", options.Hidden_markov_model_directory+"/strict_cutoffs.txt")

    Reports.load_and_process_hit_distributions(options.fasta_alignment_directory, options.database_directory)




   
    
            



    

    
    
  
def main(args=None):
#1
    myUtil.print_header("\n 1. Preparing space for the results")
    Project.prepare_directory_structure(__location__)
    options = parse_arguments(args)
    Project.prepare_result_space(options)
    
    
    if options.stage <= 2 and options.end >= 2:
        myUtil.print_header("\n 2. Prokaryotic gene recognition and translation via prodigal")
        fasta_preparation(options)
        

    
#2    
    if options.stage <= 3 and options.end >= 3:
        myUtil.print_header("\n 3. Searching for homologoues sequences")
        initial_search(options)
        
#3    
    #cluster the adjacent protein sequences, which are not annotated yet into groups called pcs combined with a number
    #cluster with mmseqs the divergent_output_file
    #Info: this resulted in an error when only a single genome was used. linclust did this error and dumped the core. the tsv file was not created and subsequent errors occured therefore

    #if options.stage <= 3 and options.end >= 3:
    #    myUtil.print_header("\n 3. Clustering sequences by similarity")
    #    cluster_sequences(options)
    
#4    
    #csb naming
    if options.stage <= 4 and options.end >= 4:
        myUtil.print_header("\n 4. Searching for collinear syntenic blocks")
        csb_finder(options)
#5
    if options.stage <= 5 and options.end >= 5:
        myUtil.print_header("\n 5. Preparing training data fasta files")
        generate_csb_sequence_fasta(options)   

#6    
    #align
    if options.stage <= 6 and options.end >= 6:
        myUtil.print_header("\n 6. Aligning sequences")
        model_alignment(options)
        
    
#7        
    #make cross validation files
    #Validation.CV(options.fasta_output_directory,options.cores)
    if options.stage <= 7 and options.end >= 7:
        myUtil.print_header("\n 7. Performing cross valdation procedure")
        cross_validation(options)

#8  
    #mach eine weitere HMMsearch mit dem vollen model auf alle seqs und prüfe die treffer
    if options.stage <= 8 and options.end >= 8:
        report_cv_performance(options)
    
    
    

if __name__ == "__main__":
    args = sys.argv[1:]
    
    main(args)

