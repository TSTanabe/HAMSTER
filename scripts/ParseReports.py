#!/usr/bin/python
from Bio import SearchIO
from Bio import SeqIO
from collections import deque
from . import myUtil
import re


#TODO e value in protein speichern // erstmal auslassen
#TODO check validity of input in class Protein

class Protein:
    """
    The class Protein organizes protein domains. When constructed the firt domain has to be added and assigned to a proteinID. The proteinID and the list of domains are accessible from outside. Also the coordinates, scores and HMM names are accessible as strings separated by "-". When a domain is added after the construction it is checked for overlapping sequence coordinates. If coordinates overlap in any way the novel domain has to have a higher score than all overlapped domains. If new domain has a lower score than any previously added domain new domain is not added.
    This follows the assumption that the HMM with highest domain is normally assigned to the protein. Here the additional information of other domains is added if it does not interfere with this assumption
    
    Args:
        protein_ID - unique string as identifier
        HMM - protein type designation
        start - start coordinate of the protein type in the AA sequence
        end - end coordinate of the protein type in the AA sequence
        score - bitscore, propability of the dignated protein type
    """


    def __init__(self,proteinID,HMM,start=0,end=0,score=1): 
    #TODO check validity of input or exception throw
        #Protein attributes
        self.proteinID = proteinID
        self.genomeID = ""
        self.protein_sequence = ""
        
        #Gene attributes
        self.gene_contig = ""
        self.gene_start = 0
        self.gene_end = 0
        self.gene_strand = "."
        self.gene_locustag = ""
        
        #Cluster attributes
        self.cluster_ID = "" #reziprok mit Cluster class
        self.keywords = {} #reziprok mit Cluster class
        
        self.domains = {}   # dictionary start coordinate => Domain object
        self.add_domain(HMM,start,end,score)
        

    ##### Getter ####
    def get_proteinID(self):
        return self.proteinID
        
    def get_genomeID(self):
        return self.genomeID
            
    def get_domains(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key].get_HMM())
        return '-'.join(listing)
    
    def get_domains_dict(self):
    #return dict
        return self.domains
            
    def get_domain_listing(self):
    #return list
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key])
        return listing
                        
    def get_domain_coordinates(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(f"{self.domains[key].get_start()}:{self.domains[key].get_end()}")
        return '-'.join(listing)

    def get_domain_scores(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(f"{self.domains[key].get_score()}")
        return '-'.join(listing)   

    def get_domain_count(self):
        return len(self.domains)
        
    def get_protein_string(self):
        #3.9.22 representation of the whole protein in one line
        a = self.proteinID
        b = self.get_domains()
        c = self.get_domain_scores()
        d = self.gene_contig
        e = self.gene_start
        f = self.gene_end
        #d = self.get_protein_sequence()
        string = f"{a} {b} {c} {d} {e} {f}"
        return string
        
    def get_protein_list(self):
        #3.9.22 representation of the whole protein in one line
        listing = []
        listing.append(self.proteinID)
        listing.append(self.get_domains())
        listing.append(str (self.get_domain_scores()))
        listing.append(str (self.get_domain_coordinates()))
        listing.append(self.gene_contig)
        listing.append(str (self.gene_start))
        listing.append(str (self.gene_end))
        listing.append(self.gene_strand)
        listing.append(self.gene_locustag)
        #d = self.get_protein_sequence()
        #string = f"{a} {b} {c} {d} {e} {f}"
        return listing
            
    def get_gene_contig(self):
        return self.gene_contig
        
    def get_gene_start(self):
        return self.gene_start
        
    def get_gene_end(self):
        return self.gene_end
        
    def get_gene_strand(self):
        return self.gene_strand
        
    def get_gene_locustag(self):
        return self.gene_locustag
        
    def get_protein_sequence(self):
        return self.protein_sequence
        
    def get_clusterID(self):
        return self.cluster_ID

    def get_sequence(self):
        return str(self.protein_sequence)
        
    def get_length(self):
        return len(self.protein_sequence)
    ##### Setter #####
    def set_genomeID(self,string):
        self.genomeID = string
        return
        
    def set_gene_contig(self,string):
        self.gene_contig = string
        return
        
    def set_gene_start(self,integer):
        self.gene_start = int(integer)
        return
        
    def set_gene_end(self,integer):
        self.gene_end = int(integer)
        return
        
    def set_gene_strand(self,string):
        self.gene_strand = string
        return
        
    def set_gene_locustag(self,string):
        self.gene_locustag = string
        return
        
    def set_protein_sequence(self,string):
        self.protein_sequence = string
        return
        
    def set_clusterID(self,string):
        self.cluster_ID = string
        return
    
    def check_domain_overlap(self,new_start, new_end,\
    current_start,current_end):
    #2.9.22

        if (current_start <= new_start and new_start <= current_end)\
        or (current_start <= new_end and new_end <= current_end):
        #start oder endpunkt innerhalb er grenzen
            return 1

        elif (new_start <= current_start and current_end <= new_end)\
        or (current_start <= new_start and new_end <= current_end):
        #start und end innerhalb der grenzen oder alte domäne innerhalb der neuen
            return 1
            
        elif (new_start <= current_start and new_end <= current_start)\
        or (new_start >= current_end and new_end >= current_end):
        #start und end kleiner als self start oder start und end größer als self end dann
        #domäne außerhalb der alten domäne und adden egal welcher score
            return 0
        return None
        

    def add_domain(self,HMM,start,end,score):
        """
        2.9.22
        Adds a domain to an existing protein. Only if there is no overlap with a 
        currently existing domain or of the overlapping domain scores higher than
        any existing overlapped domain. Overlapped domains with minor scores are deleted
        
        Args:
            HMM - protein type designation
            start - start coordinate of the protein type in the AA sequence
            end - end coordinate of the protein type in the AA sequence
            score - bitscore, propability of the dignated protein type
        Return: 
            0 - no domain was added
            1 - domain was added
        """

        del_domains = [] # start coordinates/keys of domains to be replace
        for domain in self.domains.values():
            if self.check_domain_overlap(start,end,domain.get_start(),domain.get_end()):
                if domain.get_score() < score:
                    del_domains.append(domain.get_start())
                else:
                    return 0


        
        for key in del_domains:
            self.domains.pop(key)
        self.domains.update({start:Domain(HMM,start,end,score)}) # if loop complete
        
        return 1
        
        
        


class Domain:
#2.9.22
    def __init__(self,HMM,start,end,score):
        self.HMM = HMM
        self.start = int(start)
        self.end = int(end)
        self.score = int(score)
    
    def __hash__(self):
        return hash((self.HMM, self.start, self.end, self.score))
    
    def __eq__(self, other):
        if isinstance(other, Domain):
            return self.HMM == other.HMM and self.start == other.start and self.end == other.end and self.score == other.score
        return False
    
    def get_HMM(self):
        return self.HMM
    def get_start(self):
        return self.start
    def get_end(self):
        return self.end
    def get_score(self):
        return self.score
    
#   Parsing subroutines
#------------------------------------------------------------

def parseGFFfile(Filepath,protein_dict,n=0):
    """
    3.9.22
    
    Adds the general genomic features to Protein Objects in a dictionary
    
    Args:
        Filepath - GFF3 formatted file
        protein_dict - Dictionary with key proteinID and value Protein Objects
        n - number of upstream and downstream genes extending the original hit vicinity
    Return:
        protein_dict (even though possibly not necessary)
    """
    queue = deque(maxlen=n)
    vicinity_protein_data = {}
    last_match = 0
    
    with open(Filepath,"r") as reader:
        for line in reader.readlines():
            if line.startswith("#"): #skip comments
                continue
            match = re.search('ID=(cds-){0,1}(\S+?)\W{0,1};',line)
            match = match.group(2) #find the ID of the protein in the line
            if match in protein_dict:
                #print(protein_dict[match].get_protein())
                #string teilen und dem 
                #protein hinzufügen
                protein = protein_dict[match]
                gff = line.split("\t")
                protein.set_gene_contig(gff[0])
                protein.set_gene_start(gff[3])
                protein.set_gene_end(gff[4])
                protein.set_gene_strand(gff[6])
                locustag = getLocustag(line)
                protein.set_gene_locustag(locustag)
                #add the information for the hit to the dictionary

                #add genes from queue to the data dict
                #Add upstream dictionary
                counter = 0
                while queue:
                    protein = queue.popleft() # a list [proteinID, contig, start, end, strand]
                    vicinity_protein_data[protein.proteinID] = protein
                    
                    counter = counter +1
                
                last_match = n
                    
            elif last_match:
                #Add the downstream elements to dictionary
                last_match = last_match-1
                gff = line.split("\t")
                protein = Protein(match,"default")
                protein.set_gene_contig(gff[0])
                protein.set_gene_start(gff[3])
                protein.set_gene_end(gff[4])
                protein.set_gene_strand(gff[6])
                vicinity_protein_data[match] = protein
                
            else:
                gff = line.split("\t")
                protein = Protein(match,"default")
                protein.set_gene_contig(gff[0])
                protein.set_gene_start(gff[3])
                protein.set_gene_end(gff[4])
                protein.set_gene_strand(gff[6])

                queue.append(protein) #queue is a dqueue with length argument n. default 0
                
    return vicinity_protein_data


    
def getProteinSequence(Filepath,protein_dict):
    """
    3.9.22
    Adds the protein Sequence to the Protein Objects in a dictionary. ProteinIDs of the dictionary have
    to match the header of the .faa file
    
    Args:
        Filepath - fasta formatted amino acid sequence containing file
        protein_dict - Dictionary with key proteinID and value Protein Objects
    Return:
        protein_dict (even though possibly not necessary)    
    """
    with open(Filepath) as reader:
        for record in SeqIO.parse(reader,"fasta"):
            if record.id in protein_dict:
                protein = protein_dict[record.id]
                protein.set_protein_sequence(record.seq)
                
    return protein_dict


def writeProteinSequenceFasta(vicinity_genes_dict,genomeID,output_file):
    with open(output_file,"a") as writer:
        for protein in vicinity_genes_dict.values():
            proteinID = f"{genomeID}-{protein.proteinID}"
            writer.write(f">{proteinID}\n")
            writer.write(f"{protein.protein_sequence}\n")


def getLocustag(string):
    match = re.search('locus\_tag=(\S*?)[\n|;]',string)
    if match:
        return match.group(1)
    else:
        return ""

def writeQueryHitsSequenceFasta(protein_dict,genomeID,directory):

    for protein in protein_dict.values():
        
        proteinID = f"{genomeID}-{protein.proteinID}"
        
        domains = protein.get_domains()
        file_path = directory+"/Query_"+domains+".hit_list"
        #print(proteinID," ",domains)
        File = open(file_path, "a")
        try:
            File.write(f">{proteinID}\n")
            File.write(f"{protein.protein_sequence}\n")
        finally:
            # Close the file to ensure proper handling
            File.close()

def parseBlastreportXML(Filepath,Thresholds):
    
    protein_dict = {}
    threshold = 10

        
    

    for blast_record in SearchIO.parse(Filepath,"blast-xml"):
        query = blast_record.id    # Name of query
        if query in Thresholds:
            threshold = Thresholds[query] # specific threshold
        else:
            threshold = 10
                
        for hit in blast_record.hits:
            if threshold<hit.hsps[0].bitscore:
                # Get the protein identifier of the hit
                hit_proteinID = hit.id

                # Get the hit score
                hit_bitscore = hit.hsps[0].bitscore

                # Get the start and stop coordinates of the matching sequence
                hsp_start = hit.hsps[0].hit_start
                hsp_end = hit.hsps[0].hit_end
                #print(f"HMM: {query} ProteinID:{hit_proteinID} Hitscore:{hit_bitscore} Start:{hsp_start} End:{hsp_end}") #Debugging Zeile
                if hit_proteinID in protein_dict:
                    protein = protein_dict[hit_proteinID]
                    
                    protein.add_domain(query,hsp_start,hsp_end,hit_bitscore)
                else:
                    protein_dict[hit_proteinID] = Protein(hit_proteinID,query,hsp_start,hsp_end,hit_bitscore)
                    
    return protein_dict   
    
                            
def parseBlastreporttab(Filepath,Thresholds):
    
    protein_dict = {}
    threshold = 10


    for blast_record in SearchIO.parse(Filepath,"blast-tab"):
        query = blast_record.id    # Name of query
        if query in Thresholds:
            threshold = Thresholds[query] # specific threshold
        else:
            threshold = 10
                
        for hit in blast_record.hits:
            if threshold<hit.hsps[0].bitscore:
                # Get the protein identifier of the hit
                hit_proteinID = hit.id

                # Get the hit score
                hit_bitscore = hit.hsps[0].bitscore

                # Get the start and stop coordinates of the matching sequence
                hsp_start = hit.hsps[0].hit_start
                hsp_end = hit.hsps[0].hit_end
                #print(f"HMM: {query} ProteinID:{hit_proteinID} Hitscore:{hit_bitscore} Start:{hsp_start} End:{hsp_end}") #Debugging Zeile
                if hit_proteinID in protein_dict:
                    protein = protein_dict[hit_proteinID]
                    
                    protein.add_domain(query,hsp_start,hsp_end,hit_bitscore)
                else:
                    protein_dict[hit_proteinID] = Protein(hit_proteinID,query,hsp_start,hsp_end,hit_bitscore)
                    
    return protein_dict                          


#Parser for linclust

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






