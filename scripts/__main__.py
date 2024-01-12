#!/usr/bin/python

import os
import sys
import argparse
import msgpack
from datetime import datetime
from multiprocessing import Pool, cpu_count

from . import Algorithm
from . import Alignment
from . import Arguments
from . import Csb_finder
from . import Cluster
from . import Validation


from . import Database
from . import myUtil
from . import ParseReports
from . import Project
from . import Search
from . import Translation
from . import Reports
from . import ParallelSearch




#get location of script or executable
#for the output module report 

#TODO deal with query sequences that do not have a csb but are query proteins
#TODO die bezeichnung für die cross validation von ** trennzeichen erweiterin sodass auch die beschreibung vorhanden ist.
#TODO bei CV weitere werte einfügen, sollte es sich bei FP und FN um proteine des gleichen clusters aber anderen csb handeln so sollte das nicht automatisch zu punktabzug führen, vllt zwei matrizen ausgeben

#TODO result folder in options.result_directory auf einen ordner ändern der nicht innerhalb der scripts liegt
if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
print(__location__)

class Options:
    def __init__(self):
        self.stage = 0
        self.fasta_file_directory = None
        self.query_file = None
        self.cores = 4
        self.evalue = 0.1

        self.query_names = [] #the names of the sequences used as querys
        self.result_files_directory = __location__+"/results"
        self.database_directory = None
        self.fasta_initial_hit_directory = __location__+"/results/Hit_lists"
        self.divergent_output_file = __location__+"/results/div_output_file.faa"
        self.fasta_output_directory = __location__+"/results/Sequences"
        self.fasta_alignment_directory = __location__+"/results/Alignments"
        self.Hidden_markov_model_directory = __location__+"/results/Hidden_markov_models"
        self.Cross_validation_directory = __location__+"/results/Cross_validation"
        self.report_output_file = __location__+"/Report.tsv"
        self.thresholds_output_file = __location__+"/Thresholds.tsv"
        self.csb_output_file = __location__+"/Csb.tsv"
        self.gene_clusters_file = __location__+"/All_gene_clusters.tsv"
        self.data_computed_Instances_json = __location__+"/csb_instances.json"
        
        
        self.redundancy_hash = dict()
        self.redundant = 0
        self.non_redundant = 0
        
        self.queued_genomes = {}
        self.faa_files = {}
        self.gff_files = {}
        self.taxonomyfile = ""
        
        self.redo_search = 0

        self.distance = 3500        
        
        self.computed_Instances_dict = 0
        
        
        self.alignment_mode = 2
        self.searchcoverage = 0.000
        self.clustercoverage = 0.100
        self.minseqid = 0.000
        self.gene_vicinity_range = 0
        
        self.min_seqs = 1 #minimal number of sequences per protein domain and key to be saved to fasta in step 5
        self.alnmax_seqs = 2000
        self.alnmax_len = 3000
        self.tcs_filter = 4
        self.gap_col_remove = 90

        self.sequence_faa_file = None #filepath to all sequences fasta file
        self.reports = dict()
        self.relaxed_reports = dict()
        self.standard_cutoff_report_file = None
        self.standard_performance_report_file = None
        self.relaxed_cutoff_report_file = None
        self.relaxed_performance_report_file = None
        

def parse_arguments(arguments):
    formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=72,width =200)
    
    parser = argparse.ArgumentParser(formatter_class=formatter, description = "",epilog = "")
    parser.add_argument('-s', dest='stage', type=int, default = 0, choices= [0,1,2,3,4,5,6,7,8,9],help='Start at stage')
    parser.add_argument('-c', dest='cores', type=int, help='Available CPUs')
    resources = parser.add_argument_group("Search resources")
    resources.add_argument('-f', dest='fasta_file_directory', type=myUtil.dir_path, help='Directory of the fasta files')
    resources.add_argument('-r', dest='result_files_directory', type=myUtil.dir_path, help='Directory for the result files/results from a previous run') # project folder
    resources.add_argument('-q', dest='query_file', type=myUtil.file_path, help='Query sequences fasta file')
    
    search = parser.add_argument_group("Search parameters and sequence clustering parameters")
    search.add_argument('-evalue', dest='evalue', type=float, default = 0.1, help='E-value cutoff [0,inf]')
    search.add_argument('--alignment-mode',dest='alignment_mode',type=int,default=2, choices=[0,1,2,3,4], help='mmseqs2 search alignment mode')
    search.add_argument('--min-seq-id',dest='minseqid',type=float,default=0.000, help='mmseqs2 search list matches above this sequence identity [0.0,1.0]')
    search.add_argument('--search-coverage', dest='searchcoverage', type=float,default=0.000, help='Min. coverage used for searching')
    search.add_argument('--gene-vicinity-range', dest='gene_vicinity_range', type=int,default=0, help='Number of genes up- and downstream of a hit, taken into account for gene vicinity')
    search.add_argument('--cluster-coverage', dest='clustercoverage', type=float,default = 0.000, help='Min. coverage used for clustering sequences')
    search.add_argument('--cluster:min-seq-id',dest='cminseqid',type=float,default=0.000, help='mmseqs2 search list matches above this sequence identity [0.0,1.0]')
    
    genecluster = parser.add_argument_group("Gene cluster parameters")
    genecluster.add_argument('--distance', dest='distance', type=int,default = 3500, help='Max. distance between synthenic genes')
    
    csb = parser.add_argument_group("Collinear synthenic block parameters")
    csb.add_argument('-insertions', dest='insertions', type=int,default = 2, help='Max. insertions in a csb')
    csb.add_argument('-occurence', dest='occurence', type=int,default = 2, help='Min. number of csb occurs at least times')
    csb.add_argument('-jaccard', dest='jaccard', type=float,default = 0.2, help='Acceptable dissimilarity in jaccard clustering. 0.2 means that 80 % have to be the same genes')
    csb.add_argument('-min_csb_size', dest='min_csb_size', type=int,default = 2, help='Min. csb occurence before recognized as csb')
    
    alignment = parser.add_argument_group("Alignment parameters")
    alignment.add_argument('--min_seqs', dest='min_seqs', type=int, help='Description of min_seqs')
    alignment.add_argument('--max_seqs', dest='alnmax_seqs', type=int, default = 2000, help='Max. number of sequences that are accepted for alignment')
    alignment.add_argument('--max_len', dest='alnmax_len', type=int, default = 3000, help='Max. length of sequences that are accepted for alignment')
    alignment.add_argument('--tcs_level', dest='tcs_filter', type=int, default = 4, choices = [1,2,3,4,5,6,7,8,9], help='M-coffee tcs filter level [1-9]')
    alignment.add_argument('--gap_col_remove', dest='gap_col_remove', type=int, default = 90, help=' Remove alignment columns with % gaps [0,100]')

    options = Options()
    parser.parse_args(namespace=options)

    return options
    
#parser.add_argument('-db', dest='database_directory', type=myUtil.file_path, help='Description of database_directory')
    

    

def initial_search(options):
    print(f"MMseqs search with -alignment-mode {options.alignment_mode} --min-seq-id {options.minseqid} -e {options.evalue} -c {options.searchcoverage}")
    Project.queue_files(options)
    if options.cores == 1:
        ParallelSearch.initial_search(options)
    elif options.cores >= 1:
        ParallelSearch.multi_search_process(options)
    
    
    
    
def cluster_sequences(options):
    #groups sequences via linclust
    #TODO clean after clustering the tmp files.    
    Database.index_database(options.database_directory)
    if os.path.isfile(options.divergent_output_file) and not os.path.getsize(options.divergent_output_file) == 0 and options.gene_vicinity_range > 0:
        similarity_cluster_file = Search.MMseqsLinclust(options.divergent_output_file,options.cores)
        similarity_cluster_groups = ParseReports.parse_linclust_find_connected_groups(similarity_cluster_file) #a list of sets
        #print(similarity_cluster_groups)
        Database.update_database_protein_type(options.database_directory,similarity_cluster_groups)

    prefix = "Query_"
    matching_files = [os.path.join(options.fasta_initial_hit_directory, f) for f in os.listdir(options.fasta_initial_hit_directory) if f.startswith(prefix)]
    #matching_files = [f for f in os.listdir(options.fasta_initial_hit_directory) if f.startswith(prefix)]
    
    print(f"MMseqs linclust with --threads {options.cores} --min-seq-id {options.cminseqid} --alignment-mode {options.alignment_mode} -e {options.evalue} -c {options.clustercoverage}")
    for fasta_file in matching_files:
        if not os.path.getsize(fasta_file) == 0:
            similarity_cluster_file = Search.MMseqsLinclust(fasta_file,options)
            similarity_cluster_groups = ParseReports.parse_linclust_find_connected_groups(similarity_cluster_file) #a list of sets
            #print(similarity_cluster_groups)
            #get the name of the file
            fasta_file_name, extension = os.path.splitext(fasta_file)
            domain_type = fasta_file_name.split(prefix, 1)[1]
            
            #assigns the protein names: number_queryname
            Database.update_database_protein_type(options.database_directory,similarity_cluster_groups,domain_type,domain_type)



def csb_finder(options):
    genomeIDs = Database.fetch_genomeIDs(options.database_directory)
     
    with open(options.gene_clusters_file, "w") as file: # write down the gene cluster from each genome into one file
        for genomeID in genomeIDs:
            cluster_diction = Database.fetch_cluster_dict(options.database_directory, genomeID)
            for clusterID, cluster in cluster_diction.items():
                domains = cluster.get_domains()
                file.write(clusterID + '\t' + '\t'.join(domains) + '\n')
    

    options.redundant,options.non_redundant = Cluster.dereplicate(options.gene_clusters_file) #returns two filepaths, dereplicates identical gene clusters
    
    
    options.redundancy_hash = Cluster.create_redundancy_hash(options.redundant) # value is an integer, number of existing replicates
    gene_clusters = Cluster.create_gene_cluster_hash(options.non_redundant)
    Cluster.extend_redundancy_hash(options.non_redundant,options.redundancy_hash)

    #modified CsbfinderS algorithm
    computed_Instances_dict = Algorithm.csb_finderS_matchpoint_algorithm(options.redundancy_hash,gene_clusters,options.insertions,options.occurence) #k und q müssen über die optionen festgelegt werden
    options.computed_Instances_dict = computed_Instances_dict
    
    data_serializable = {key: list(value) for key, value in computed_Instances_dict.copy().items()}

    serialized_data = msgpack.packb(data_serializable,use_bin_type=True)
    with open(f"{options.result_files_directory}/csb_instances.json",'wb') as file:
        file.write(serialized_data)
    

    
    
def csb_jaccard(options):
    #convert the keys in computed_instances_dict into a list
    if not options.computed_Instances_dict and options.data_computed_Instances_json:
        with open(f"{options.data_computed_Instances_json}",'rb') as file:
            serialized_data = file.read()
        options.computed_Instances_dict = msgpack.unpackb(serialized_data,raw=False)
        options.computed_Instances_dict = {key: set(value) for key, value in options.computed_Instances_dict.items()}
    elif not options.computed_Instances_dict and not options.data_computed_Instances_json:
        print(f"ERROR: No {options.data_computed_Instances_json} found")
        sys.exit()
    computed_Instances_key_list = Cluster.csb_Instance_key_list(options.computed_Instances_dict, options.min_csb_size)
    
    print("\n computed_Instances_key_list: ",computed_Instances_key_list)
    cluster_dict = dict()
    if len(computed_Instances_key_list) > 1:
        matrix = Cluster.calculate_similarity_matrix_jaccard(computed_Instances_key_list) #matrix of similarity and the corresponding clusterID for each row and column as names
        cluster_dict = Cluster.hierachy_clustering(matrix,options.jaccard) # 0.2 means that 80 % have to be the same genes
    elif len(computed_Instances_key_list) == 1:
        cluster_dict[0] = [0]
    else:
        return
    #Sorts the csb keywords to the geneclusters based on the present csb
    csb_gene_cluster_dict, grouped_csb_tuples = Cluster.csb_index_to_gene_clusterID(cluster_dict,computed_Instances_key_list,options.computed_Instances_dict)
    
    #adds the replicates again for saving
    Cluster.replicates(csb_gene_cluster_dict,options.redundancy_hash,options.redundant)
    Cluster.write_grouped_csb(options.csb_output_file,grouped_csb_tuples)
    Database.update_keywords(options.database_directory,csb_gene_cluster_dict) #assigns the names of the keywords to the clusters



def generate_csb_sequence_fasta(options):
    keyword_domain_tuples = Database.fetch_unique_keyword_domain_pairs(options.database_directory)
    Database.fetch_domains_per_key_to_fasta(options.database_directory, keyword_domain_tuples, options.fasta_output_directory, options.min_seqs) #Output are files with {domain}_{keyword}.faa



def cross_validation(options):
    
    #File with all sequences to be searched
    #sequence_faa_file = Validation.concat_files_with_extension(options.fasta_initial_hit_directory,"hit_list",options.Cross_validation_directory+"/sequences.faa") # all seqs
    options.sequence_faa_file = Database.fetch_all_proteins(options.database_directory, options.Cross_validation_directory+"/sequences.faa")
    
    
    TN = Validation.count_sequences_in_fasta(options.sequence_faa_file)
    alignment_dict = Validation.get_alignment_files(options.fasta_output_directory) #aligned_seqs => unaligned_seqs filepaths for keys and values
    alignment_files = list(alignment_dict.keys())
    
    with Pool(processes=options.cores) as pool:
        results = pool.map(Validation.process_cross_folds, [(alignment_file, TN, options.sequence_faa_file,options) for alignment_file in alignment_files])

    for f, reports in results:
        options.reports[f] = report[0] #returns ReportObject(sum_matrix, fold_matrices, best_optimized_cutoff, highest_trusted_cutoff, lowest_noise_cutoff)
        options.relaxed_reports[f] = report[1]








       
def report_performance(options):
    #soll die performance von jedem nehmen die werte berechnen und zusammen mit den cutoffs aufschreiben.
    #cutoffs gesondert aufschreiben für die spätere verwendung von einer HMM lib
    
    #das ausgabepaket schnüren:
    #HMMs in voll die erstellt wurden
    Validation.create_hmms_from_msas(options.fasta_output_directory,"fasta_aln","hmm") #create the full hmms for later use
    Reports.move_HMMs(options.fasta_output_directory,options.Hidden_markov_model_directory,"hmm") #move the hmms to the Hidden markov model folder
    
    options.standard_cutoff_report_file = Reports.create_cutoff_file(options,options.reports,options.Hidden_markov_model_directory)
    options.standard_performance_report_file = Reports.create_performance_file(options,options.reports,options.Hidden_markov_model_directory)

    options.relaxed_cutoff_report_file = Reports.create_cutoff_file(options,options.relaxed_reports,options.Hidden_markov_model_directory)
    options.relaxed_performance_report_file = Reports.create_performance_file(options,options.relaxed_reports,options.Hidden_markov_model_directory)
    
    
def interference_performance(options):

    #make the HMMsearch for all sequences
    concat_HMMlib = Validation.concat_files_with_extension(options.Hidden_markov_model_directory, ".hmm", options.Hidden_markov_model_directory+"concat_HMMlib.cat")
    output_file = options.Hidden_markov_model_directory+"all_seqs_search_results"
    Validation.HMMsearch(concat_HMMlib,options.sequence_faa_file, output_file, cores=1)

    Thresholds = Search.makeThresholdDict(options.cutoff_file,1) #standard optimized cutoff
    standard_cutoff_restext = Validation.convert_HMMreport_to_sets(output_file,Thresholds,".std") #returns a textfile
    
    Thresholds = Search.makeThresholdDict(options.relaxed_cutoff_file,1) #relaxed optimized cutoff
    relaxed_cutoff_restext = Validation.convert_HMMreport_to_sets(output_file,Thresholds,".rel") #returns a textfile
    
    #the textfiles are ordered by query name, then TP,FP,FN
    
    
  
def main(args=None):
#1
    options = parse_arguments(args)
    Project.prepare_result_space(options)
    if options.stage <= 1:
        myUtil.print_header("\nProkaryotic gene recognition and translation via prodigal")
        Translation.translation(options.fasta_file_directory)      #fine
        Translation.transcription(options.fasta_file_directory)    #fine
    

    
#2    
    #get the homologues via mmseqs
    if options.stage <= 2:
        initial_search(options)
    
#3    
    #cluster the adjacent protein sequences, which are not annotated yet into groups called pcs combined with a number
    #cluster with mmseqs the divergent_output_file
    #Info: this resulted in an error when only a single genome was used. linclust did this error and dumped the core. the tsv file was not created and subsequent errors occured therefore

    if options.stage <= 3:
        cluster_sequences(options)
    
#4    
    #csb naming
    ##per genome   
    if options.stage <= 4:
        csb_finder(options)
        csb_jaccard(options)
#5
    if options.stage <= 5:
        generate_csb_sequence_fasta(options)   
    
#6    
    #align
    if options.stage <= 6:
        Alignment.align_tcs(options)
    
#7        
    #make cross validation files
    #Validation.CV(options.fasta_output_directory,options.cores)
    if options.stage <= 7:
        cross_validation(options)
        report_performance(options)
        

#8  #TODO interferenz zwischen HMMs messen
    #mach eine weitere HMMsearch mit dem vollen model auf alle seqs und prüfe die treffer
    #die FP Sequenzen ausgeben
    #die FN Sequenzen ausgeben
    if options.stage <= 8:
        interference_performance(options)
    
    
    

if __name__ == "__main__":
    args = sys.argv[1:]
    
    main(args) #calls the main method of __main__

