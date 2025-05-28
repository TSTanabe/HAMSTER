#!/usr/bin/python

import os
import sys
import argparse

from . import Project
from . import Queue

from . import Database
from . import Translation
from . import ParallelSearch

from . import Csb_cluster
from . import Csb_proteins
from . import Csb_statistic
from . import Csb_mcl

from . import Alignment
from . import Validation
from . import myUtil

from . import Pam_Singleton_finder
from . import Pam_defragmentation
from . import Pam_mcl

from . import Reports_plotting
#from . import Reports_printing
from . import Reports

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
        
        self.csb_name_prefix = "csb-" #prefix of clusterIDs determined by csb finder algorithm
        self.csb_name_suffix = "_" #suffix of clusterIDs determined by csb finder algorithm
        
        self.self_query = None #Query fasta like concatenated fasta
        self.self_seqs = None #Query sequences

        self.deconcat_flag = 0 #Memorize if a deconcatenation was done
        self.glob_flag = 0 #Memorize if a concatenation was done
        
        
def parse_arguments(arguments):
    formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=96,width =300)
    
    #parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description = "HAMSTER version 0.0.8 \nSyntax: HAMSTER [OPTIONS]",epilog = "")
    parser = argparse.ArgumentParser(
    description="HAMSTER",
    epilog="",
    usage=argparse.SUPPRESS,   # <— ganz ausschalten
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    essential = parser.add_argument_group("Provide either (-q & -f) or -db")
    essential.add_argument('-f', dest='fasta_file_directory', type=myUtil.dir_path, default = __location__, metavar = '<directory>', help='Directory of the target fasta files, need -q option')
    essential.add_argument('-q', dest='query_file', type=myUtil.file_path, default = None, metavar = '<filepath>', help='Query sequences fasta file, need -f option')
    essential.add_argument('-db',dest='database_directory', type=myUtil.file_path, metavar='<filepath>', help='Existing HAMSTER sqlite database')
    

    resources = parser.add_argument_group("Optional parameters")
    resources.add_argument('-s', dest='stage', type=int, default = 0, choices= [0,1,2,3,4,5,6,7,8,9,10], metavar='<int>', help='Start pipeline at stage')
    resources.add_argument('-x', dest='end', type=int, default = 10, choices= [0,1,2,3,4,5,6,7,8,9,10], metavar='<int>', help='End pipeline at stage')
    resources.add_argument('-c', dest='cores', type=int, default = 2, metavar='<int>', help='Use number of CPUs')
    resources.add_argument('-t', dest='taxonomy_file', type=myUtil.file_path, default = None, metavar = '<filepath>', help='Taxonomy csv file fro input assemblies')
    resources.add_argument('-r', dest='result_files_directory', type=myUtil.dir_path, default = __location__+"/results", metavar = '<directory>', help='Resultfile directory (of previous run)')
    resources.add_argument('-glob_off', dest='glob_search', action='store_false', help='Do not concatenated fasta file for initial diamond search')
    resources.add_argument('-cv_off', dest='cross_validation_deactivated', action='store_true', help='Skip cross-validation step')
   
    resources2 = parser.add_argument_group("Optional glob assembly file parameters")
    resources2.add_argument('-glob_faa', dest='glob_faa', type=myUtil.file_path, default=None, metavar = '<filepath>', help='Concatenated multi-assembly fasta file')
    resources2.add_argument('-glob_gff', dest='glob_gff', type=myUtil.file_path, default=None, metavar = '<filepath>', help='Concatenated multi-assembly gff file')
    resources2.add_argument('-glob_blast_table', dest='glob_table', type=myUtil.file_path, default=None, metavar = '<filepath>', help='Concatenated multi-assembly blast result table')
    resources2.add_argument('-glob_chunks', dest='glob_chunks', type=int, default=3000, metavar='<int>', help='Chunk size for parsing glob results')
    
    
    
    search = parser.add_argument_group("Optional search parameters for diamond")
    search.add_argument('-evalue', dest='evalue', type=float, default = 1e-5, metavar = '<float>', help='E-value cutoff [0,inf]')
    search.add_argument('-thrs_score', dest='thrs_score', type=int, default = 100, metavar = '<int>', help='Score cutoff [0,inf]')
    search.add_argument('-min-seq-id',dest='minseqid',type=float, default=25, metavar = '<float>', help='Sequence search matches above this sequence identity [0,100.0]')
    search.add_argument('-search-coverage', dest='searchcoverage', type=float, default=0.6, metavar = '<float>', help='Min. coverage used for searching [0.0,1.0]')
    search.add_argument('-alignment-mode',dest='alignment_mode',type=int, default=2, metavar='<int>', choices=[0,1,2], help='DIAMOND BLASTp alignment mode')
    search.add_argument('-memory',dest='diamond_memory_limit',type=int, default=8, metavar='<int>', help='DIAMOND BLASTp RAM limit')
    search.add_argument('-blast-score-ratio', dest='thrs_bsr', type=float, default=0.0, metavar = '<float>', help='Blast score ratio for hits [0.0,1.0]')
    search.add_argument('-allow_multidomain', dest='multidomain_allowed', action='store_true', help='Allow multiple query hits for each sequence')
    search.add_argument('-reports_hit', dest='diamond_report_hits_limit', type=int, default=0, metavar = '<int>', help='Limit to this number of top hits per query')



    #Cluster parameters
    #protein_cluster = parser.add_argument_group("Sequence clustering parameters for mmseqs2")
    #protein_cluster.add_argument('-cluster-active', dest='protein_cluster_active', action='store_true', help='Cluster initial blastp hits with mmseqs2 cluster')
    #protein_cluster.add_argument('-alignment-mode',dest='alignment_mode',type=int, default=2, metavar='<int>', choices=[0,1,2,3,4], help='mmseqs2 cluster search alignment mode')
    #protein_cluster.add_argument('-cluster-coverage', dest='clustercoverage', type=float, default = 0.800, metavar='<float>', help='mmseqs2 cluster min. coverage used for clustering sequences')
    #protein_cluster.add_argument('-cluster-min-seq-id',dest='cminseqid',type=float, default=0.000, metavar='<float>', help='mmseqs2 search list matches above this sequence identity [0.0,1.0]')


    
    genecluster = parser.add_argument_group("Optional gene cluster prediction parameters")
    genecluster.add_argument('-distance', dest='nucleotide_range', type=int, default = 3500, metavar='<int>', help='Max. nucleotide distance between synthenic genes')
    genecluster.add_argument('-p', dest= 'patterns_file' , type=myUtil.file_path, default=__location__+"/src/Patterns", metavar='<filepath>', help='Filepath to predefined tab separated patterns file')
    
    csb = parser.add_argument_group("Optional collinear synthenic block parameters")
    csb.add_argument('-insertions', dest='insertions', type=int,default = 2, metavar='<int>', help='Max. insertions in a csb')
    csb.add_argument('-occurence', dest='occurence', type=int,default = 10, metavar='<int>', help='Min. number of csb occurs at least time')
    csb.add_argument('-min_csb_size', dest='min_csb_size', type=int,default = 4, metavar='<int>', help='Min. csb size before recognized as csb')
    csb.add_argument('-jaccard', dest='jaccard', type=float,default = 0.0, metavar='<float>', help='Acceptable dissimilarity in jaccard clustering. 0.2 means that 80 percent have to be the same genes')
    csb.add_argument('-csb_overlap', dest='csb_overlap_factor', type=float, default = 0.75, metavar='<float>', help='Merge if sequences from two csb is identical above this threshold')
    
    #csb.add_argument('-no_phylogeny', dest='csb_distinct_grouping', action='store_false', help='Skip phylogenetic supported training dataset clustering')
    csb.add_argument('-exclude_csb_hitscore', dest='low_hitscore_csb_cutoff', type=float,default = 0.3, metavar='<float>', help='Exclude csb with all hits below this deviation from the query self-hits')
    csb.add_argument('-group_csb_hitscore', dest='group_hitscore_csb_cutoff', type=float,default = 0.3, metavar='<float>', help='Group hits from csb above this deviation from the query self-hit')
    csb.add_argument('-csb_singleton_correlation', dest='csb_singleton_correlation', type=float,default = 0.3, metavar='<float>', help='Required correlation of singletons to existing csb')

    pam_search = parser.add_argument_group("Optional presence/absence matrix parameters")
    pam_search.add_argument('-mx_thrs', dest='pam_threshold', type=float, default=0.3, metavar = '<float>', help='Presence/absence matrix co-occurence significance parameter [0.0,1.0]')
    pam_search.add_argument('-mx_bsr', dest='pam_bsr_threshold', type=float, default=0.6, metavar = '<float>', help='Blast score ratio inclusion cutoff from PAM prediction [0.0,1.0]')
    
    #pam_search.add_argument('-mx_long_branch_thrs', dest='pam_long_branch_thres', type=float, default=0.5, metavar = '<float>', help='Do not consider hits with this branch length in the phylogenetic tree. Default: 0.5')
    #pam_search.add_argument('-mx_max_phylo_distance', dest='pam_phylogenetic_distance', type=float, default=0.05, metavar = '<float>', help='Max. phylogenetic distance to reference sequence. Default: 0.05')

    mcl_search = parser.add_argument_group("Optional Markov Chain Clustering parameters")
    mcl_search.add_argument('-mcl_off', dest='csb_mcl_clustering', action='store_false', help='Skip markov chain clustering')
    mcl_search.add_argument('-mcl_evalue', dest='mcl_evalue', type=float, default = 1e-10, metavar = '<float>', help='MCL matrix e-value cutoff [0,inf]')
    mcl_search.add_argument('-mcl_min-seq-id',dest='mcl_minseqid',type=float, default=50, metavar = '<float>', help='MCL matrix sequence identity cutoff [0,100.0]')
    mcl_search.add_argument('-mcl_search-coverage', dest='mcl_searchcoverage', type=float, default=0.6, metavar = '<float>', help='MCL matrix min. coverage [0.0,1.0]')
    mcl_search.add_argument('-mcl_hit_limit', dest='mcl_hit_limit', type=int, default=100, metavar = '<int>', help='MCL maximum number of edges between sequences')
    mcl_search.add_argument('-mcl_inflation', dest='mcl_inflation', type=float, default=2.0, metavar = '<float>', help='MCL inflation factor for granularity control')
    mcl_search.add_argument('-mcl_sensitivity', dest='mcl_sensitivity', type=str, choices=["fast", "more-sensitive", "sensitive", "very-sensitive", "ultra-sensitive"], default="sensitive", metavar='<sensitivity>', help="Set DIAMOND sensitivity for MCL clustering. Choices: fast, more-sensitive, sensitive, very-sensitive, ultra-sensitive")
    mcl_search.add_argument('-mcl_density_thrs', dest='mcl_density_thrs', type=float, default=None, metavar = '<float>', help='Required proportion of reference sequences in the total number of sequences in the MCL cluster to label it as true positive [0.0,1.0]')
    mcl_search.add_argument('-mcl_reference_thrs', dest='mcl_reference_thrs', type=float, default=None, metavar = '<float>', help='Required proportion of reference sequences from the total reference sequences in the MCL cluster to label it as true positive [0.0,1.0]')
    
    #This is how the MCL reference thresholds work
    #ref_count = len(seqs.intersection(reference_sequences)) # absolute number of reference sequences in the cluster
    #cluster_size = len(seqs) # absoulte number of sequences in the cluster
    #ref_density = ref_count / cluster_size if cluster_size > 0 else 0  # Fraction of reference sequences in the total number of sequencse in the current mcl cluster
    #ref_fraction = ref_count / len(reference_sequences) if len(reference_sequences) > 0 else 0  # Fraction of total reference sequences in this cluster
    #ref_density >= density_threshold and ref_fraction >= reference_fraction_threshold:
    
    
    alignment = parser.add_argument_group("Optional alignment parameters")
    alignment.add_argument('-min_seqs', dest='min_seqs', type=int, default = 5, metavar='<int>', help='Min. number of required sequences for the alignment')
    alignment.add_argument('-max_seqs', dest='max_seqs', type=int, default = 100000, metavar='<int>', help='Max. number of sequences that are aligned')
    alignment.add_argument('-gap_col_remove', dest='gap_remove_threshold', type=float, default = 0.05, metavar='<float>', help='[0,1] remove alignment columns with only percent amino acids')
    alignment.add_argument('-include_domains', dest='include_list',nargs='+', default=[], metavar='<list>', help='Space separated protein list, that are specifically included')
    alignment.add_argument('-exclude_domains', dest='exclude_list', nargs='+', default=[], metavar='<list>', help='Space separated protein list, that are specifically excluded')
    
    
    if len(arguments) == 0: # Print help if no arguments were provided
        parser.print_help()
        sys.exit("No arguments were provided.")

    options = Options()
    parser.parse_args(namespace=options)
    options.location = __location__
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
        print("[INFO] Not concatenating protein fasta files")
        Queue.queue_files(options)
        return
    
    # Concatenate the fasta files, in case of a glob fasta file is provided deconcat it
    if options.glob_faa and options.glob_gff:
        print(f"[INFO] Preparing separate files from {options.glob_faa}")
        #create separated files
        Translation.deconcat(options)
        Queue.queue_files(options)
        options.deconcat_flag = 1
    else:

        # Unpacks and translates fasta
        Translation.parallel_translation(options.fasta_file_directory, options.cores)
        Translation.parallel_transcription(options.fasta_file_directory, options.cores)

        # concat to to globfile if search results are not already provided
        Queue.queue_files(options)
        if not options.glob_table:
            print("[INFO] Generating concatenated glob faa file", end="\r")
            Translation.create_glob_file(options) #fasta_file_directory, options.cores, concat the files with the genomeIdentifier+ ___ + proteinIdentifier
            print("[INFO] Generating concatenated glob faa file -- ok")
            options.glob_flag = 1   
    
    # Create self query file
    Translation.create_selfquery_file(options)
    
    return

def initial_search(options):
    # Writes a database with the protein hits and their gene clusters for later use.
    # Writes fasta files for each hit for linclust    

    if not os.path.isfile(options.database_directory):
        Database.create_database(options.database_directory)

    if options.glob_search:
        # Search the glob file
        ParallelSearch.initial_glob_search(options)
    else:
        # Process all files separately
        ParallelSearch.initial_genomize_search(options)

    if options.deconcat_flag: # deconcatenation was done but is not needed anymore
        myUtil.remove_directory(options.fasta_file_directory)# remove the deconcatenated files
    elif options.glob_flag:
        os.remove(options.glob_faa)
    

def csb_finder(options):

    Csb_cluster.csb_prediction(options)
    csb_gene_cluster_dict = Csb_cluster.csb_jaccard(options)

    Database.index_database(options.database_directory)
    Database.delete_keywords_from_csb(options.database_directory, options) # remove keys with options.csb_name_prefix options.csb_name_suffix to avoid old keyword interference
    Database.update_keywords(options.database_directory,csb_gene_cluster_dict) # assigns the names of the keywords to the clusters



def basis_sequence_fasta(options):
    """
        prepare the sequence fasta and identifier lists. Seqs are derived from
        named gene clusters and csb finder gene clusters that encode a protein sequence
        highly similar to the input reference sequence.
        
        Input: Database, options
        Output are the grp0 fasta files
    """
    Database.index_database(options.database_directory)
    grp_score_limit_dict, grouped = Csb_proteins.prepare_csb_grouped_training_proteins(options)
    
    ### Collect singletons without any conserved csb
    print("\n\n[INFO] Collecting highly similar homologs from query hits without any conserved genomic context")
    sng_score_limit_dict, sng_ref_seqs_dict = Pam_Singleton_finder.singleton_reference_finder(options, grouped) #Very strict reference seqs with 0.95 blast score ratio

    # Merge groups and limits from csb and sng
    merged_score_limit_dict = {**grp_score_limit_dict, **sng_score_limit_dict}
    merged_grouped = {**grouped, **sng_ref_seqs_dict}


    # Save the merged groups
    myUtil.save_cache(options, 'basis_merged_grouped.pkl', merged_grouped)
    myUtil.save_cache(options, 'basis_merged_score.pkl', merged_score_limit_dict)
    
    # Update options object with the fetched proteinID groups and score limits
    options.grouped = merged_grouped
    options.score_limit_dict = merged_score_limit_dict
    
    return

def pam_defragmentation(options):
    """
        Try to find additional hits based on presence absence patterns
        
        Prepare the presence absence matrix for grp0 and train logistic regression on the matrix        
        With the trained matrix exclude each column and predict presence.
        For predicted presences select from the genome the best hit

        Output: are the grp1 fasta files
    """
    basis_grouped = options.grouped if hasattr(options, 'grouped') else myUtil.load_cache(options,'basis_merged_grouped.pkl')
    basis_score_limit_dict = options.score_limit_dict if hasattr(options, 'score_limit_dict') else myUtil.load_cache(options, 'basis_merged_score.pkl')
    
    merged_grouped, merged_score_limit_dict = Pam_defragmentation.pam_genome_defragmentation_hit_finder(options, basis_grouped, basis_score_limit_dict)
    
    # Save computed grp1 datasets
    myUtil.save_cache(options, 'grp1_merged_grouped.pkl', merged_grouped)
    myUtil.save_cache(options, 'grp1_merged_score_limits.pkl', merged_score_limit_dict)
    
    # Result dictionary is stores in options.grouped, overwriting the grp0 with grp1 key_domain pairs
    options.grouped = merged_grouped
    options.score_limit_dict = merged_score_limit_dict
    
    return





def mcl_family_clustering_sequences(options):
    """
        Try to find additional hits based on Markov Chain Clustering
        
        Prepare protein family file, including all hits above 25% identity.
        All vs. all diamond blast (Can this be accellerated like in orthofinder?)
        
        With the trained matrix exclude each column and predict presence.
        For predicted presences select from the genome the best hit

        Output: are the grp1 fasta files
    """
    grouped = options.grouped if hasattr(options, 'grouped') else myUtil.load_cache(options,'grp1_merged_grouped.pkl')
    score_limit_dict = options.score_limit_dict if hasattr(options, 'score_limit_dict') else myUtil.load_cache(options, 'grp1_merged_score_limits.pkl')

    print("\n[INFO] Including homologs without genomic context based on protein sequence clustering")
    print("[INFO] Generating protein family fasta files")
    Csb_proteins.decorate_training_data(options, score_limit_dict, grouped)
        
    # Markov chain clustering for grp1 fasta files
    if not options.csb_mcl_clustering:
        return
    
    print("\n[INFO] Calculating Markov Chain Clustering")            
    mcl_clustering_results_dict = Csb_mcl.csb_mcl_datasets(options,grouped) # markov chain clustering grouped training data

    myUtil.save_cache(options, 'mcl_clustering_results.pkl', mcl_clustering_results_dict)

    options.mcl_clustering_results_dict = mcl_clustering_results_dict
    
    
    
def mcl_decorate_training_sequences(options):

    grouped = options.grouped if hasattr(options, 'grouped') else myUtil.load_cache(options,'grp1_merged_grouped.pkl')
    score_limit_dict = options.score_limit_dict if hasattr(options, 'score_limit_dict') else myUtil.load_cache(options, 'grp1_merged_score_limits.pkl')
    mcl_clustering_results_dict = options.mcl_clusterin_results_dict if hasattr(options, 'mcl_clustering_results_dict') else myUtil.load_cache(options, 'mcl_clustering_results.pkl')
    
    mcl_clustering_results_dict = Csb_mcl.validate_mcl_cluster_paths(mcl_clustering_results_dict, options.result_files_directory) # Check for path existence
    
    
    # MCL cluster analysis, cluster selection with F1 optimized density for basic grouped sequences
    print("\n[INFO] Generating grp1: Selecting MCL clusters with sufficient fraction of reference sequences with conserved genomic vicinity")
    extended_grouped, mcl_cutoffs = Csb_mcl.select_hits_by_csb_mcl(options, mcl_clustering_results_dict, grouped, options.mcl_density_thrs, options.mcl_reference_thrs)
    extended_grouped_prefixed = {f"grp1_{key}": value for key, value in extended_grouped.items()}

    # Save the mcl cluster selection cutoffs    
    myUtil.save_cache(options, 'mcl_grp1_cluster_selection_cutoffs.pkl', mcl_cutoffs)
    
    # Write down the extended
    Csb_proteins.fetch_seqs_to_fasta_parallel(
            options.database_directory,
            extended_grouped_prefixed,
            options.fasta_output_directory,
            options.min_seqs,
            options.max_seqs,
            options.cores
        )
    
    print("\n[INFO] Generating grp2: Extended MCL cluster selection by csb and presence plausibility")
    # MCL cluster analysis with pam and csb selection
    regrouped = Pam_mcl.select_hits_by_pam_csb_mcl(options, mcl_clustering_results_dict, grouped)
    
    # Save the regrouped reference seqs used to select inclusion of the grp2 mcl clusters
    myUtil.save_cache(options, 'grp2_selection_ref_seqs.pkl', regrouped)
    
    # MCL cluster analysis, cluster selection with F1 optimized density for extended grouped
    extended_grouped, mcl_cutoffs = Csb_mcl.select_hits_by_csb_mcl(options, mcl_clustering_results_dict, regrouped, options.mcl_density_thrs, options.mcl_reference_thrs)
    extended_grouped_prefixed = {f"grp2_{key}": value for key, value in extended_grouped.items()}
    
    # Save the mcl cluster selection cutoffs
    myUtil.save_cache(options, 'mcl_grp2_cluster_selection_cutoffs.pkl', mcl_cutoffs)
    
    # Write down the extended
    Csb_proteins.fetch_seqs_to_fasta_parallel(
            options.database_directory,
            extended_grouped_prefixed,
            options.fasta_output_directory,
            options.min_seqs,
            options.max_seqs,
            options.cores
        )

def mcl_decorate_training_sequences2(options):
    
    print("\n[INFO] Generating grp3: Extended MCL cluster selection by all sequences with similar genomic vicinity")

    # Get the 















           
    # Phylogeny clustering for lca fasta files
    #if options.csb_distinct_grouping:
        
    #    # Decorate grouped high scoring csb with similarly high seqs from the same phylogenetic clade
    #    print("\nCalculating phylogeny for each protein family")
        
    #    decorated_grouped_dict = Csb_phylogeny.csb_phylogeny_datasets(options, grouped) # phylogenetic grouped training data
        
    #    Csb_proteins.fetch_seqs_to_fasta_parallel(options.database_directory, decorated_grouped_dict, options.fasta_output_directory, options.min_seqs, options.max_seqs, options.cores)
    





        
def model_alignment(options):
    
    Alignment.initial_alignments(options, options.fasta_output_directory)

    return

def cross_validation(options):

    print("Initialize training data subsamples")

    Validation.create_hmms_from_msas(options.fasta_output_directory, options.Hidden_markov_model_directory, "fasta_aln","hmm",options.cores) #create the full hmms for later use

    Reports.move_HMMs(options.fasta_output_directory,options.Hidden_markov_model_directory,"hmm") #move the hmms to the Hidden markov model folder
    
    Csb_proteins.fetch_domains_superfamily_to_fasta(options, options.cross_validation_directory) #Target file for each HMM excluding seqs already below threshold

    Csb_proteins.fetch_all_proteins(options.database_directory, options.cross_validation_directory+"/sequences.faa") #Backup target file if something fails


    print("Initialize target sequence sets")
    
    options.sequence_faa_file = options.cross_validation_directory+"/sequences.faa" #File with all sequences to be searched
    
    options.targeted_sequence_faa_file_dict = Validation.get_target_sets(options.cross_validation_directory)

    Validation.initial_self_recognition_validation(options)
    
    if options.cross_validation_deactivated:
        return
    Validation.parallel_cross_validation(options)

    return
       
def report_cv_performance(options):
    
    #Initial validation
    print(f"Saving the cutoffs and performance reports from initial calculation to {options.Hidden_markov_model_directory}")
    
    Reports_plotting.process_initial_validations(options, options.result_files_directory, options.fasta_alignment_directory, options.database_directory)
    
    #cross validation
    print(f"Saving the cutoffs and performance reports from the cross-validatio to {options.Hidden_markov_model_directory}")
    
    options.reports = Reports.parse_matrices_to_report(options.cross_validation_directory,"_cv_matrices.txt")
    
    Reports.create_performance_file(options,options.reports,options.Hidden_markov_model_directory,"/cv_cutoff_performance.txt")
    
    Reports.concatenate_cv_cutoff_files(options.cross_validation_directory, "_cv_thresholds.txt", options.Hidden_markov_model_directory+"/cv_strict_cutoffs.txt")





   
    
            



    

    
    
  
def main(args=None):
#1
    Project.prepare_directory_structure(__location__)
    options = parse_arguments(args)
    myUtil.print_header("\n 1. Preparing space for the results")
    Project.prepare_result_space(options)
    
#2    
    if options.stage <= 2 and options.end >= 2:
        myUtil.print_header("\n 2. Prokaryotic gene recognition and translation via prodigal")
        fasta_preparation(options)
    
#3    
    if options.stage <= 3 and options.end >= 3:
        myUtil.print_header("\n 3. Searching for homologoues sequences")
        initial_search(options)
        
#4    
    #csb naming
    if options.stage <= 4 and options.end >= 4:
        myUtil.print_header("\n 4. Searching for collinear syntenic blocks")
        csb_finder(options)

#5
    if options.stage <= 5 and options.end >= 5:
        myUtil.print_header("\n 5. Preparing training data fasta files")
        basis_sequence_fasta(options)
        pam_defragmentation(options)

#6
    if options.stage <= 6 and options.end >= 6:
        myUtil.print_header("\n 6. Markov chain clustering for protein family sequences")
        mcl_family_clustering_sequences(options)
        
#7   
    if options.stage <= 7 and options.end >= 7:
        myUtil.print_header("\n 7. Selecting MCL clusters with reference sequences")
        mcl_decorate_training_sequences(options)

#8    
    #align
    if options.stage <= 8 and options.end >= 8:
        myUtil.print_header("\n 8. Aligning sequences")
        model_alignment(options)
        
#9        
    #make cross validation files
    #Validation.CV(options.fasta_output_directory,options.cores)
    if options.stage <= 9 and options.end >= 9:
        myUtil.print_header("\n 9. Performing cross valdation procedure")
        cross_validation(options)

#10  
    #mach eine weitere HMMsearch mit dem vollen model auf alle seqs und prüfe die treffer
    if options.stage <= 10 and options.end >= 10:
        myUtil.print_header("\n 10. Writing dataset performance reports")
        report_cv_performance(options)
    
    
    

if __name__ == "__main__":
    args = sys.argv[1:]
    
    main(args)

