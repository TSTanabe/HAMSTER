#!/usr/bin/python
import re
import subprocess
from typing import Dict, Optional
from . import Database
from . import myUtil

logger = myUtil.logger


class Protein:
    """
    The class Protein organizes protein domains. When constructed the firt domain has to be added and assigned to a proteinID. The proteinID and the list of domains are accessible from outside. Also the coordinates, scores and HMM names are accessible as strings separated by "-". When a domain is added after the construction it is checked for overlapping sequence coordinates. If coordinates overlap in any way the novel domain has to have a higher score than all overlapped domains. If new domain has a lower score than any previously added domain new domain is not added.
    This follows the assumption that the HMM with highest domain is normally assigned to the protein. Here the additional information of other domains is added if it does not interfere with this assumption
    
    Organizes protein domains and related attributes for a protein. 
    New domains are only added if they do not overlap with higher-scoring existing domains.
    All coordinates, scores, and HMM names are accessible as dash-separated strings.

    Args:
        proteinID (str): Unique identifier.
        HMM (str): Domain name.
        start (int): Start coordinate.
        end (int): End coordinate.
        score (float): Domain bitscore.
        genomeID (str): Genome identifier (optional).
        ident (int): Percent identity (default: 25).
        bsr (float): Blast score ratio (default: 1.0).

    Example:
        p = Protein("Prot1", "HMM_A", 10, 120, 40)
        p.add_domain("HMM_B", 130, 200, 30)
    """


    def __init__(
        self,
        proteinID: str,
        HMM: str,
        start: int = 0,
        end: int = 0,
        score: float = 1,
        genomeID: str = "",
        ident: int = 25,
        bsr: float = 1.0
    ):
        self.proteinID: str = proteinID
        self.genomeID: str = genomeID
        self.protein_sequence: str = ""
        self.gene_contig: str = ""
        self.gene_start: int = 0
        self.gene_end: int = 0
        self.gene_strand: str = "."
        self.gene_locustag: str = ""
        self.clusterID: str = ""
        self.keywords: Dict = {}
        self.domains: Dict[int, Domain] = {}  # start coordinate → Domain object
        self.add_domain(HMM, start, end, score, ident, bsr)
        

    ##### Getter ####
            
    def get_domains(self):
    #return string
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key].get_HMM())
        return '-'.join(listing)
    
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

    def get_sequence(self):
        return str(self.protein_sequence)
        
    ##### Setter #####

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
        


    def add_domain(
        self,
        HMM: str,
        start: int,
        end: int,
        score: float,
        ident: int = 25,
        bsr: float = 1.0
    ) -> int:
        """
        Adds a domain to the protein only if it does not overlap
        with a higher-scoring existing domain. If overlap exists with lower-scoring
        domain, that domain is removed.

        Returns:
            int: 1 if domain added, 0 if not added.
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
        self.domains.update({start:Domain(HMM,start,end,score,ident,bsr)}) # if loop complete
        
        return 1
        
        
        


class Domain:
#2.9.22
    """
    Stores domain information (HMM name, coordinates, score, identity, bsr).

    Args:
        HMM (str): Domain name.
        start (int): Start coord.
        end (int): End coord.
        score (float): Bitscore.
        ident (int): Percent identity.
        bsr (float): Blast score ratio.
    """
    def __init__(self, HMM: str, start: int, end: int, score: float, ident: int = 1, bsr: float = 1.0):
        self.HMM: str = HMM
        self.start: int = int(start)
        self.end: int = int(end)
        self.score: float = float(score)
        self.identity: int = int(ident)
        self.bsr: float = float(bsr)
    
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




def parseGFFfile(
    filepath: str, protein_dict: Dict[str, Protein]
) -> Dict[str, Protein]:
    """
    3.9.22
    Adds GFF attributes to each Protein object in the dictionary.

    Args:
        filepath (str): Path to GFF3 file.
        protein_dict (dict): {proteinID: Protein object}

    Returns:
        dict: Updated protein_dict.
    """
    locustag_pattern = re.compile(r'locus_tag=(\S*?)(?:[;\s]|$)')
    geneID_pattern = re.compile(r'ID=(cds-)?(\S+?)(?:[;\s]|$)')
    
    grep_pattern = "|".join(protein_dict.keys())
    try:
        grep_process = subprocess.Popen(['grep', '-E', grep_pattern, filepath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = grep_process.communicate()

        if stderr:
            logger.error("Grep process:", stderr)
            return protein_dict
        
        for line in stdout.decode('utf-8').split('\n'):
            if not line:
                continue
            gff = line.split("\t")
            match = geneID_pattern.search(gff[-1])
            if not match:
                continue
            match = match.group(2)
            if match in protein_dict:
                # Add protein
                protein = protein_dict[match]
                
                protein.gene_contig = str(gff[0])
                protein.gene_start = int(gff[3])
                protein.gene_end = int(gff[4])
                protein.gene_strand = str(gff[6])
                locustag = getLocustag(locustag_pattern, line)
                protein.gene_locustag = str(locustag)
    except Exception as e:
        logger.warning(f"Error occurred while parsing GFF: {str(e)}")
        return protein_dict
    return protein_dict


def getProteinSequence(Filepath, protein_dict):
    """
    3.9.22
    Fügt die Proteinsequenz zu den Proteinobjekten in einem Dictionary hinzu. Die Protein-IDs im Dictionary
    müssen mit den Headern der .faa-Datei übereinstimmen.

    Args:
        Filepath: Pfad zur Datei im FASTA-Format mit Aminosäuresequenzen
        protein_dict: Dictionary mit Protein-ID als Schlüssel und Protein-Objekten als Werte
    Return:
        protein_dict (obwohl möglicherweise nicht notwendig)
    """
    
    reader = None
    try:
        reader = open(Filepath, "r")
        sequence = ""
        header = None
        save_sequence = False  # Indikator, ob die Sequenz gespeichert werden muss
        
        for line in reader:
            line = line.strip()
            
            if line.startswith(">"):
                # Speichere die bisherige Sequenz, falls notwendig
                if header and save_sequence and sequence:
                    protein = protein_dict[header]
                    protein.protein_sequence = sequence

                # Extrahiere die Protein-ID aus dem Header (bis zum ersten Leerzeichen)
                header = line[1:].split()[0]
                
                # Prüfe, ob diese Protein-ID im Dictionary ist
                if header in protein_dict:
                    save_sequence = True  # Nur dann sammeln wir die Sequenz
                    sequence = ""  # Setze die Sequenz zurück für das neue Protein
                else:
                    save_sequence = False  # Ignoriere Sequenzen, die nicht im Dictionary sind

            elif save_sequence:
                # Füge die Sequenzzeile nur hinzu, wenn die ID im Dictionary ist
                sequence += line
        
        # Verarbeite die letzte Sequenz, falls nötig
        if header and save_sequence and sequence:
            protein = protein_dict[header]
            protein.protein_sequence = sequence
    
    except IOError:
        logger.error(f"Cannot open {filepath}")
    
    finally:
        if reader is not None:
            reader.close()  # Datei explizit schließen
    
    return protein_dict


def getLocustag(locustag_pattern: re.Pattern, string: str) -> str:
    """
    Extracts locus_tag from a string using a regex pattern.

    Args:
        locustag_pattern (re.Pattern): Regex pattern for locus_tag.
        string (str): Line to search.

    Returns:
        str: locus_tag or empty string if not found.
    """
    match = locustag_pattern.search(string)
    return match.group(1) if match else ""

