#!/usr/bin/python

import os
import sys
import shutil
from datetime import datetime
from . import Database
from . import myUtil

def create_project(directory,projectname="project"):
    now = datetime.now()
    timestamp = datetime.timestamp(now)
    directory = directory + "/" + str(timestamp) + "_" + projectname
    try:
        os.mkdir(directory)
    except:
        sys.exit("\nERROR: Could not create result directory. No writing rights.")
        raise Exception(f"\nERROR: No writing rights.")
    return directory
    
def create_logfile(directory):
    logfile = directory +"/logfile"
    #print(logfile)
    with open (logfile,"w") as writer:
        writer.write("Date\tgenomeID\tLibrary\tdatabase\tfile\n")
    
    return logfile

def create_genomize_dir(directory):
    os.mkdir(directory+"/per_genome")
    return  directory+"/per_genome"

def read_logfile(filepath):
    used_files = set()
    #f = open(filepath,"r")
    #print(f.readlines())
    with open(filepath,"r") as reader:
        line = reader.readline()  # skip header
        for line in reader.readlines():
            line = line.replace("\n","")
            array = line.split("\t")
            #used_files[array[1]] = 1    # genomeID into array
            used_files.add(array[1])
    return used_files

def write_logfile_entry(genomeID,options):
    array = [str(datetime.now()),genomeID,myUtil.getFileName(options.library),options.database_directory,options.faa_files[genomeID]]
    string = "\t".join(array)+"\n"
    with open (options.logfile,"a") as writer:
        writer.write(string)
    writer.close()
    return

def create_taxonomyfile(directory):
    taxon = directory +"/taxonomy"
    #print(logfile)
    with open (taxon,"w") as writer:
                       #filename superkingdom clade phylum class order family genus species strain taxid biosample bioproject genbank refseq completeness contamination typestrain
                        writer.write("genomeID\tsuperkingdom\tclade\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\ttaxid\tbiosample\tbioproject\tgenbank\trefseq\tcompleteness\tcontamination\ttypestrain\n")
    return taxon

def isProjectFolder(options):
    if not os.path.isfile(options.result_files_directory+"/database.db"):
        print("Miss db")
        return 0
    if not os.path.isdir(options.result_files_directory+"/Sequences"):
        print("Miss seq")
        return 0
    if not os.path.isdir(options.result_files_directory+"/Alignments"):
        print("Miss aln")
        return 0
    if not os.path.isdir(options.result_files_directory+"/Hit_list"):
        print("Miss hit list")
        return 0
    if not os.path.isdir(options.result_files_directory+"/Hidden_markov_models"):
        print("Miss HMM")
        return 0
    if not os.path.isdir(options.result_files_directory+"/Cross_validation"):
        print("Miss CV")
        return 0
    return 1
    
    
    

    
def prepare_result_space(options,project="project"):
    """
        6.10.22
        Args:
            directory   directory for the result files to be stored
        Return:
            path        directory path for the result files
            continue    boolean if old project should be continues
        Options taken are
        Flowchart for this routine
            Checks if user input directory is already present
             if not 
                try to create this one
             if yes 
                check if database is already present
                 if no
                    create new one
                 if yes
                    continue writing on this existing database
            when database is created follow with the other report files        
            
    """
    myUtil.print_header("\nPreparing space for the results")
    logfile = 0
    now = datetime.now()
    timestamp = str(datetime.timestamp(now))
    #feststellen ob es sich um einen project folder handelt
    if isProjectFolder(options):
        
        print(options.result_files_directory)
        options.database_directory = options.result_files_directory+"/database.db"
        options.fasta_initial_hit_directory = options.result_files_directory+"/Hit_list"
        print(options.fasta_initial_hit_directory)
        options.divergent_output_file = options.result_files_directory+"/div_output_file.faa"
        options.fasta_output_directory = options.result_files_directory+"/Sequences"
        options.fasta_alignment_directory = options.result_files_directory+"/Alignments"
        
        options.Hidden_markov_model_directory = options.result_files_directory+"/Hidden_markov_models"
        options.Cross_validation_directory = options.result_files_directory+"/Cross_validation"
        print(options.Cross_validation_directory)
        options.divergent_output_file = options.result_files_directory+"/div_output_file.faa"
        options.report_output_file = options.result_files_directory+"/_Report.txt"
        options.thresholds_output_file = options.result_files_directory+"/_Thresholds.txt"
        options.csb_output_file = options.result_files_directory+"/Csb_output.txt"
        options.gene_clusters_file = options.result_files_directory+"/_All_gene_clusters.txt"
        if os.path.isfile(options.result_files_directory+"/csb_instances.json"):
            options.data_computed_Instances_json = options.result_files_directory+"/csb_instances.json"
        else:
            options.data_computed_Instances_json = 0

    else:
        if not os.path.isdir(options.result_files_directory):
            try:
                os.mkdir(options.result_files_directory)
            except:
                sys.exit("\nERROR: Could not create result directory. Missing writing rights.")
                raise Exception(f"\nERROR: No writing rights.")
    
        options.result_files_directory = create_project(options.result_files_directory,project)        

        options.database_directory = options.result_files_directory+"/database.db"
        Database.create_database(options.database_directory)

        options.fasta_initial_hit_directory = options.result_files_directory+"/Hit_list"
        os.mkdir(options.fasta_initial_hit_directory)
        options.fasta_output_directory = options.result_files_directory+"/Sequences"
        os.mkdir(options.fasta_output_directory)
        options.fasta_alignment_directory = options.result_files_directory+"/Alignments"
        os.mkdir(options.fasta_alignment_directory)
        options.Hidden_markov_model_directory = options.result_files_directory+"/Hidden_markov_models"
        os.mkdir(options.Hidden_markov_model_directory)
        options.Cross_validation_directory = options.result_files_directory+"/Cross_validation"
        os.mkdir(options.Cross_validation_directory)

        options.divergent_output_file = options.result_files_directory+"/div_output_file.faa"
        options.report_output_file = options.result_files_directory+"/"+timestamp+"_Report.txt"
        options.thresholds_output_file = options.result_files_directory+"/"+timestamp+"_Thresholds.txt"
        options.csb_output_file = options.result_files_directory+"/Csb_output.txt"
        options.gene_clusters_file = options.result_files_directory+"/"+timestamp+"_All_gene_clusters.txt"
        options.data_computed_Instances_json = options.result_files_directory+"/csb_instances.json"
   
    
          
def queue_files(options):
    """
        9.10.22
        Args:
            directory   directory for the result files to be stored
            finished    set with genome identifiers already processed
            options     current options object
        Operation:
            collect all zip protein fasta files and unzipped protein fasta files
            collect all corresponding gff files and check for existance of this file
            if both files present
            get the genome identifiers
    """
    print("\nFilling the queue with faa files --")
    queue = set()
    faa_files = {}
    gff_files = {}

    #For all zipped .faa files find corresponding zipped or unzipped gff. Fill hashes with zipped format only
    all_zip_FaaFiles = myUtil.getAllFiles(options.fasta_file_directory,".faa.gz")
    for zip_FaaFile in all_zip_FaaFiles:
        faa_file = myUtil.removeExtension(zip_FaaFile)
        file_name = myUtil.removeExtension(faa_file)
        #print(file_name+".gff.gz")
        if os.path.isfile(file_name+".gff.gz"):
            gff_file = file_name+".gff.gz"
            genomeID = myUtil.getGenomeID(faa_file)
            queue.add(genomeID) #GenomeIDs that are queued
            faa_files[genomeID] = zip_FaaFile
            gff_files[genomeID] = gff_file
        elif os.path.isfile(file_name+".gff"):
            gff_file = myUtil.packgz(file_name+".gff")
            myUtil.unlink(file_name+".gff")
            genomeID = myUtil.getGenomeID(faa_file)
            queue.add(genomeID)
            faa_files[genomeID] = zip_FaaFile
            gff_files[genomeID] = gff_file
        else:
            print(f"\tWARNING: No corresponding .gff file found for {zip_FaaFile}")
    
    #For all unzipped .faa files find corresponding zipped or unzipped gff. Fill hashes with zipped format only
    all_FaaFiles = myUtil.getAllFiles(options.fasta_file_directory,".faa")
    for FaaFile in all_FaaFiles:
        zip_FaaFile = myUtil.packgz(FaaFile)
        myUtil.unlink(FaaFile)
        faa_file = myUtil.removeExtension(zip_FaaFile)
        file_name = myUtil.removeExtension(faa_file)
        if os.path.isfile(file_name+".gff.gz"):
            gff_file = file_name+".gff.gz"
            genomeID = myUtil.getGenomeID(faa_file)
            queue.add(genomeID)
            faa_files[genomeID] = zip_FaaFile
            gff_files[genomeID] = gff_file
        elif os.path.isfile(file_name+".gff"):
            gff_file = myUtil.packgz(file_name+".gff")
            myUtil.unlink(file_name+".gff")
            genomeID = myUtil.getGenomeID(faa_file)
            queue.add(genomeID)
            faa_files[genomeID] = zip_FaaFile
            gff_files[genomeID] = gff_file
        else:
            print(f"\tWARNING: No corresponding .gff file found for {zip_FaaFile}")

    #compare two sets
    options.queued_genomes = queue
    options.faa_files = faa_files
    options.gff_files = gff_files
    
    
    #compare with database leaving out all genomes already searched previously
    if options.redo_search: #redo the search for all
        genomeIDs = Database.fetch_genomeIDs(options.database_directory)
        options.finished_genomes = {}
        for genomeID in genomeIDs:
            options.finished_genomes[genomeID] = 1 #requires for extending the protein dict in csb assignment with all old results
    else: #leave out all genomeIDs in the database
        genomeIDs = Database.fetch_genomeIDs(options.database_directory)
        for genomeID in genomeIDs:
            if genomeID in faa_files.keys():
                print(f"\tFound assembly {genomeID} in database leaving out {myUtil.getFileName(faa_files[genomeID])}")
                del faa_files[genomeID]
                del gff_files[genomeID]
                options.queued_genomes.remove(genomeID)
    
    print(f"Queued {len(options.queued_genomes)} for processing")
    if len(options.queued_genomes) == 0:
        print("There were 0 genomes queued. Use -redo_search option if genomes are already present in database")
    print("\nFilling the queue with faa files -- ok")




def query_names(options,query_file):
    
    with open(query_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                header = line.strip()
                options.query_names.append(header)








