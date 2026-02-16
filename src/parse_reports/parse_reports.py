#!/usr/bin/python
import re
import subprocess
import traceback
from typing import Dict, Set, Any, Optional

from src.core.logging import get_logger

logger = get_logger(__name__)


class Protein:
    """
    The class Protein organizes protein domains. When constructed the firt domain has to be added and assigned
    to a proteinID. The proteinID and the list of domains are accessible from outside. Also the coordinates, scores
    and HMM names are accessible as strings separated by "-". When a domain is added after the construction it is
    checked for overlapping sequence coordinates. If coordinates overlap in any way the novel domain has to have
    a higher score than all overlapped domains. If new domain has a lower score than any previously added domain
    new domain is not added.
    This follows the assumption that the HMM with highest domain is normally assigned to the protein.
    Here the additional information of other domains is added if it does not interfere with this assumption

    Organizes protein domains and related attributes for a protein.
    New domains are only added if they do not overlap with higher-scoring existing domains.
    All coordinates, scores, and HMM names are accessible as dash-separated strings.

    Args:
        protein_id (str): Unique identifier.
        hmm (str): Domain name.
        start (int): Start coordinate.
        end (int): End coordinate.
        score (float): Domain bitscore.
        genome_id (str): Genome identifier (optional).
        ident (int): Percent identity (default: 25).
        bsr (float): Blast score ratio (default: 1.0).

    Example:
        p = Protein("Prot1", "HMM_A", 10, 120, 40)
        p.add_domain("HMM_B", 130, 200, 30)
    """

    def __init__(
        self,
        protein_id: str,
        hmm: str,
        start: int = 0,
        end: int = 0,
        score: float = 1,
        genome_id: str = "",
        ident: int = 25,
        bsr: float = 1.0,
    ):
        self.proteinID: str = protein_id
        self.genomeID: str = genome_id
        self.protein_sequence: str = ""
        self.gene_contig: str = ""
        self.gene_start: int = 0
        self.gene_end: int = 0
        self.gene_strand: str = "."
        self.gene_locustag: str = ""
        self.clusterID: str = ""
        self.keywords: Dict = {}
        self.domains: Dict[int, Domain] = {}  # start coordinate → Domain object
        self.deleted_domains: Dict[
            str, Domain
        ] = {}  # Domain objects that have been removed
        self.add_domain(hmm, start, end, score, ident, bsr)
        self.selection_comment: Set[str] = set()  # Trusted cutoff Flag or cooccurrence
        self.alternative_hit: str = ""
        self.valid_hit = True

    ##### Getter ####

    def get_domains(self):
        # return string
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key].get_domain())
        return "-".join(listing)

    def get_domains_dict(self):
        # return dict
        return self.domains

    def get_domain_listing(self):
        # return list
        listing = []
        for key in sorted(self.domains):
            listing.append(self.domains[key])
        return listing

    def get_domain_set(self):
        domains = set()
        for v in self.domains.values():
            domains.add(v.get_domain())
        return domains

    def get_domain_coordinates(self):
        # return string
        listing = []
        for key in sorted(self.domains):
            listing.append(
                f"{self.domains[key].get_start()}:{self.domains[key].get_end()}"
            )
        return "-".join(listing)

    def get_domain_scores(self):
        # return string
        listing = []
        for key in sorted(self.domains):
            listing.append(f"{self.domains[key].get_score()}")
        return "-".join(listing)

    def get_domain_count(self):
        return len(self.domains)

    def get_protein_list(self):
        # 3.9.22 representation of the whole protein in one line
        listing = [
            self.proteinID,
            self.get_domains(),
            str(self.get_domain_scores()),
            str(self.get_domain_coordinates()),
            self.gene_contig,
            str(self.gene_start),
            str(self.gene_end),
            self.gene_strand,
            self.gene_locustag,
        ]
        # d = self.get_protein_sequence()
        # string = f"{a} {b} {c} {d} {e} {f}"
        return listing

    def get_sequence(self):
        return str(self.protein_sequence)

    def get_selection_comment_csv(self, sep: str = ",") -> str:
        """
        Return the selection_comment as a sep-separated string.
        Accepts both set[str] (preferred) and str (fallback).
        """
        val = self.selection_comment
        if not val:
            return ""
        if isinstance(val, str):
            return val  # bereits CSV-String
        # erwarteter Fall: Menge/Tokens
        return sep.join(sorted(val))

    ##### Setter #####

    @staticmethod
    def check_domain_overlap(new_start, new_end, current_start, current_end):
        # 2.9.22

        if (current_start <= new_start <= current_end) or (
            current_start <= new_end <= current_end
        ):
            # start oder endpunkt innerhalb er grenzen
            return 1

        elif (new_start <= current_start and current_end <= new_end) or (
            current_start <= new_start and new_end <= current_end
        ):
            # start und end innerhalb der grenzen oder alte domäne innerhalb der neuen
            return 1

        elif (new_start <= current_start and new_end <= current_start) or (
            new_start >= current_end and new_end >= current_end
        ):
            # start und end kleiner als self start oder start und end größer als self end dann
            # domäne außerhalb der alten domäne und adden egal welcher score
            return 0
        return None

    def add_selection_comment(self, comment: str, sep: str = ",") -> None:
        """
        Add one or multiple comment tokens to this protein.
        - Accepts a single token ("Tc") or a CSV string ("Tc,Coo").
        """
        if not comment:
            return
        tokens = [t.strip() for t in str(comment).split(sep) if t.strip()]
        self.selection_comment.update(tokens)

    def add_domain(
        self,
        hmm: str,
        start: int,
        end: int,
        score: float,
        ident: int = 25,
        bsr: float = 1.0,
        *,
        force: bool = False,
    ) -> int:
        """
        Adds a domain to the protein.

        Standardverhalten (force=False):
          - Wenn Überlappung mit existierender Domäne vorliegt:
              * Entferne überlappende Domänen mit geringerem Score
              * Brich ab (return 0), wenn eine überlappende Domäne >= Score hat.
          - Ansonsten einfügen (return 1).

        Force-Modus (force=True):
          - Ignoriere die Score-Vergleiche bei Überlappung.
          - Entferne alle überlappenden Domänen und füge die neue ein (return 1).

        Returns:
            1 wenn hinzugefügt, 0 wenn nicht hinzugefügt.
        """
        del_domains = []  # start-Koordinaten der zu entfernenden Domänen

        for domain in self.domains.values():
            if self.check_domain_overlap(
                start, end, domain.get_start(), domain.get_end()
            ):
                if force:
                    # im Force-Modus: jede überlappende Domäne räumen
                    del_domains.append(domain.get_start())
                else:
                    # Standard: nur schwächere Domänen räumen, sonst abbrechen
                    if domain.get_score() < score:
                        del_domains.append(domain.get_start())
                    else:
                        return 0

        # überlappende domänen verschieben in deleted_domains
        for key in del_domains:
            dom = self.domains.pop(key, None)
            if dom is not None:
                key = dom.domain
                self.deleted_domains[key] = dom

        # neue Domäne eintragen (Schlüssel = start)
        self.domains[start] = Domain(hmm, start, end, score, ident, bsr)

        return 1


class Domain:
    # 2.9.22
    """
    Stores domain information (HMM name, coordinates, score, identity, bsr).

    Args:
        domain (str): Domain name.
        start (int): Start coord.
        end (int): End coord.
        score (float): Bitscore.
        ident (int): Percent identity.
        bsr (float): Blast score ratio.
    """

    def __init__(
        self,
        domain: str,
        start: int,
        end: int,
        score: float,
        ident: int = 1,
        bsr: float = 1.0,
    ):
        self.domain: str = domain
        self.start: int = int(start)
        self.end: int = int(end)
        self.score: float = float(score)
        self.identity: int = int(ident)
        self.bsr: float = float(bsr)

    def __hash__(self):
        return hash((self.domain, self.start, self.end, self.score))

    def __eq__(self, other):
        if isinstance(other, Domain):
            return (
                self.domain == other.domain
                and self.start == other.start
                and self.end == other.end
                and self.score == other.score
            )
        return False

    def get_domain(self):
        return self.domain

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_score(self):
        return self.score


#########################################
########   Parsing subroutines ##########
#########################################


def parse_gff_file(
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
    locustag_pattern = re.compile(r"locus_tag=(\S*?)(?:[;\s]|$)")
    gene_id_pattern = re.compile(r"ID=(cds-)?(\S+?)(?:[;\s]|$)")

    grep_pattern = "|".join(protein_dict.keys())
    try:
        grep_process = subprocess.Popen(
            ["grep", "-E", grep_pattern, filepath],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = grep_process.communicate()

        if stderr:
            logger.error("Grep process:", stderr)
            return protein_dict

        for line in stdout.decode("utf-8").split("\n"):
            if not line:
                continue
            gff = line.split("\t")
            match = gene_id_pattern.search(gff[-1])
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
                locustag = get_locustag(locustag_pattern, line)
                protein.gene_locustag = str(locustag)
    except Exception as e:
        logger.warning(f"Error occurred while parsing GFF: {str(e)}")
        return protein_dict
    return protein_dict


def get_protein_sequence(filepath, protein_dict):
    """
    Add protein sequences from a FASTA file to Protein objects in a dictionary.

    Iterates through a FASTA file of amino acid sequences and sets the `.protein_sequence`
    attribute for each Protein object in `protein_dict`, keyed by their protein ID (header).
    Only protein IDs present in `protein_dict` are processed.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file containing amino acid sequences.
    protein_dict : dict
        Dictionary with protein IDs as keys and Protein objects as values.

    Returns
    -------
    dict
        The updated protein_dict with sequences added to the Protein objects.

    Examples
    --------
    >>> protein_dict = {'WP_0001': Protein(...), 'WP_0002': Protein(...)}
    >>> get_protein_sequence('proteins.faa', protein_dict)
    {'WP_0001': <Protein with sequence>, 'WP_0002': <Protein with sequence>}
    """
    reader = None
    try:
        reader = open(filepath, "r")
        sequence = ""
        header = None
        save_sequence = False

        for line in reader:  # type: str
            line = line.strip()
            if line.startswith(">"):
                if header and save_sequence and sequence:
                    # Save sequence for previous protein
                    protein = protein_dict[header]
                    protein.protein_sequence = sequence
                # Parse header up to first whitespace
                header = line[1:].split()[0]
                if header in protein_dict:
                    save_sequence = True
                    sequence = ""
                else:
                    save_sequence = False
            elif save_sequence:
                sequence += line

        # Handle the last protein
        if header and save_sequence and sequence:
            protein = protein_dict[header]
            protein.protein_sequence = sequence

    except IOError as e:
        logger.error(f"Cannot open {filepath}: {e}")
    finally:
        if reader is not None:
            reader.close()

    return protein_dict


def get_locustag(locustag_pattern: re.Pattern, string: str) -> str:
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


#######################################################################################################################
######### Blast report parsing routine
def clean_dict_keys_and_protein_ids(
    input_dict: Dict[str, Any], genomeID: str
) -> Dict[str, Any]:
    """
    Removes the genomeID prefix from dictionary keys and proteinID fields.

    Args:
        input_dict: dict of proteinID -> Protein
        genomeID: genomeID prefix to remove

    Returns:
        dict: updated protein dict
    """
    prefix = genomeID + "___"
    updated_dict = {}

    for key, protein in input_dict.items():
        # Entferne das Präfix von jedem Key, wenn es vorhanden ist
        new_key = key[len(prefix) :] if key.startswith(prefix) else key

        # Entferne das Präfix von proteinID, falls es vorhanden ist
        if hasattr(protein, "proteinID") and protein.proteinID.startswith(prefix):
            protein.proteinID = protein.proteinID[len(prefix) :]

        # Füge das geänderte Key-Value-Paar zum neuen Dictionary hinzu
        updated_dict[new_key] = protein

    return updated_dict


def strip_query_suffix(query_id: str) -> str:
    """
    Entfernt eine numerische Suffix-Nummerierung der Form '___<int>'.
    Beispiel:
        'SmoC___1'  -> 'SmoC'
        'SmoC___abc' -> 'SmoC___abc' (wird nicht abgeschnitten)
    """
    if "___" in query_id:
        base, suffix = query_id.rsplit("___", 1)
        if suffix.isdigit():
            return base
    return query_id


def parse_bulk_blastreport_genomize(
    genome_id: str, filepath: str, thresholds: Dict, cut_score: int = 10
) -> Dict[str, Any]:
    """
    Parses a BLAST report for a specific genomeID and extracts protein domain hits.

    Args:
        genome_id (str): The genome ID to filter hits.
        filepath (str): Path to BLAST report.
        thresholds (dict): Score thresholds (unused here).
        cut_score (int): Score threshold for filtering.

    Returns:
        dict: proteinID -> Protein object (ParseReports.Protein)

    Output Example:
        {"WP_012345678": ParseReports.Protein(...), ...}
    """

    protein_dict = {}

    try:
        # Efficiently filter relevant lines with grep
        with subprocess.Popen(
            ["grep", genome_id, filepath], stdout=subprocess.PIPE, text=True
        ) as proc:
            for line in proc.stdout:
                columns = line.strip().split("\t")  # Assuming tab-separated format

                if len(columns) < 6:  # Ensure enough columns exist
                    continue

                try:
                    # Expected BLAST format: hit_id, query_id, e_value, bitscore, hsp_start, hsp_end
                    hit_protein_id = columns[0]  # Protein identifier in BLAST
                    raw_query_id = columns[1]  # Query sequence ID
                    query_id = strip_query_suffix(raw_query_id)
                    hit_bitscore = int(float(columns[3]))  # Bitscore as integer
                    hsp_start = int(columns[4])  # Start position
                    hsp_end = int(columns[5])  # End position
                    hsp_ident = int(float(columns[6]))
                    try:
                        hsp_bsr = float(columns[7])
                    except IndexError:
                        hsp_bsr = 1.0  # oder ein anderer sinnvoller Default-Wert

                    # If protein already exists, add domain info
                    if hit_protein_id in protein_dict:
                        protein_dict[hit_protein_id].add_domain(
                            hmm=query_id,
                            start=hsp_start,
                            end=hsp_end,
                            score=hit_bitscore,
                            ident=hsp_ident,
                            bsr=hsp_bsr,
                        )
                    else:
                        # Create new protein object
                        protein_dict[hit_protein_id] = Protein(
                            protein_id=hit_protein_id,
                            hmm=query_id,
                            start=hsp_start,
                            end=hsp_end,
                            score=hit_bitscore,
                            genome_id=genome_id,
                            ident=hsp_ident,
                            bsr=hsp_bsr,
                        )

                except ValueError as ve:
                    logger.warning(
                        f"Skipped malformed line in {filepath}: {line.strip()} (ValueError: {ve})"
                    )
    except Exception as e:
        logger.error(f"Failed to parse {filepath} for genome {genome_id}: {e}")
        logger.debug(traceback.format_exc())

    return protein_dict


def parse_bulk_blastreport_consecutive(
    genome_id: str,
    filepath: str,
) -> Dict[str, Any]:
    """
    Parses a BLAST report for a specific genomeID and extracts protein domain hits.

    Args:
        genome_id (str): The genome ID to store in the Protein object.
        filepath (str): Path to BLAST report.
        thresholds (dict): Score thresholds (unused here).
        cut_score (int): Score threshold for filtering (optional, derzeit nicht genutzt).

    Returns:
        dict: proteinID -> Protein object (ParseReports.Protein)

    Output Example:
        {"WP_012345678": ParseReports.Protein(...), ...}
    """

    protein_dict: Dict[str, Protein] = {}

    try:
        with open(filepath, "r") as infile:
            for line in infile:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                columns = line.split("\t")  # Assuming tab-separated format

                # Mindestens bis hsp_end vorhanden? (0..5)
                if len(columns) < 6:
                    continue

                try:
                    hit_protein_id = columns[0]  # Protein identifier in BLAST
                    raw_query_id = columns[1]  # Query sequence ID
                    query_id = strip_query_suffix(raw_query_id)
                    hit_bitscore = int(float(columns[3]))  # Bitscore as integer
                    hsp_start = int(columns[4])
                    hsp_end = int(columns[5])
                    hsp_ident = int(float(columns[6]))
                    try:
                        hsp_bsr = float(columns[7])
                    except (IndexError, ValueError):
                        hsp_bsr = 1.0  # sinnvoller Default, wenn keine BSR-Spalte

                    # Wenn Protein schon existiert → Domain anhängen
                    if hit_protein_id in protein_dict:
                        protein_dict[hit_protein_id].add_domain(
                            hmm=query_id,
                            start=hsp_start,
                            end=hsp_end,
                            score=hit_bitscore,
                            ident=hsp_ident,
                            bsr=hsp_bsr,
                        )
                    else:
                        # Neues Protein-Objekt erzeugen
                        protein_dict[hit_protein_id] = Protein(
                            protein_id=hit_protein_id,
                            hmm=query_id,
                            start=hsp_start,
                            end=hsp_end,
                            score=hit_bitscore,
                            genome_id=genome_id,
                            ident=hsp_ident,
                            bsr=hsp_bsr,
                        )

                except ValueError as ve:
                    logger.warning(
                        f"Skipped malformed line in {filepath}: {line} "
                        f"(ValueError: {ve})"
                    )

    except Exception as e:
        logger.error(f"Failed to parse {filepath} for genome {genome_id}: {e}")
        logger.debug(traceback.format_exc())

    return protein_dict



def hit_passes_thresholds(
    *,
    raw_query_id: str,
    bitscore: float,
    evalue: float,
    sstart: int,
    send: int,
    pident: float,
    evalue_cutoff: float,
    score_cutoff: float,
    coverage_cutoff: float,
    identity_cutoff: float,
    sequence_lengths: Optional[Dict[str, float]] = None,
) -> bool:
    """
    Decide whether a hit should be kept based on threshold logic:

        Keep if:
            (A) evalue <= evalue_cutoff AND bitscore >= score_cutoff AND coverage >= coverage_cutoff
         OR (B) pident >= identity_cutoff AND coverage >= coverage_cutoff

    Coverage requires `sequence_lengths`. If missing or query length unknown -> False.
    """
    if not sequence_lengths:
        return False

    qlen = sequence_lengths.get(raw_query_id)
    if not qlen:
        return False

    alignment_length = abs(send - sstart) + 1
    coverage = alignment_length / float(qlen)

    cond_a = (
        evalue <= evalue_cutoff
        and bitscore >= score_cutoff
        and coverage >= coverage_cutoff
    )
    cond_b = (
        pident >= identity_cutoff
        and coverage >= coverage_cutoff
    )
    return cond_a or cond_b


def compute_bsr(
    *,
    raw_query_id: str,
    bitscore: float,
    selfblast_scores: Optional[Dict[str, float]] = None,
    default: float = 1.0,
) -> float:
    """
    Compute Blast Score Ratio (BSR) = bitscore / selfblast_score.

    If `selfblast_scores` is None or does not contain the query id -> return `default`.
    """
    if not selfblast_scores:
        return default

    selfscore = selfblast_scores.get(raw_query_id)
    if not selfscore:
        return default

    try:
        return bitscore / float(selfscore)
    except Exception:
        return default


def parse_filter_single_blastreport(
    genome_id: str,
    filepath: str,
    *,
    # threshold params for filtering:
    evalue_cutoff: float,
    score_cutoff: float,
    coverage_cutoff: float,
    identity_cutoff: float,
    # optional dicts:
    sequence_lengths: Optional[Dict[str, float]] = None,
    selfblast_scores: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """
    Parses a BLAST/DIAMOND report for a specific genomeID and extracts protein domain hits.

    Expected tab columns (at least 7):
        0 sseqid   (hit protein id)
        1 qseqid   (query id)
        2 evalue
        3 bitscore
        4 sstart
        5 send
        6 pident
        7 optional precomputed BSR (ignored if selfblast_scores provided and contains query id)

    Filtering:
        A hit is only added to the Protein dict if hit_passes_thresholds(...) returns True.

    Returns:
        dict: hit_protein_id -> Protein object
    """
    protein_dict: Dict[str, Protein] = {}

    try:
        with open(filepath, "r") as infile:
            for line in infile:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                columns = line.split("\t")
                if len(columns) < 7:
                    continue

                try:
                    hit_protein_id = columns[0]
                    raw_query_id = columns[1]
                    query_id = strip_query_suffix(raw_query_id)

                    evalue = float(columns[2])
                    bitscore = float(columns[3])
                    sstart = int(columns[4])
                    send = int(columns[5])
                    pident = float(columns[6])

                    # threshold filter (skip early)
                    if not hit_passes_thresholds(
                        raw_query_id=raw_query_id,
                        bitscore=bitscore,
                        evalue=evalue,
                        sstart=sstart,
                        send=send,
                        pident=pident,
                        evalue_cutoff=evalue_cutoff,
                        score_cutoff=score_cutoff,
                        coverage_cutoff=coverage_cutoff,
                        identity_cutoff=identity_cutoff,
                        sequence_lengths=sequence_lengths,
                    ):
                        continue

                    # BSR computation (kept separate)
                    bsr = compute_bsr(
                        raw_query_id=raw_query_id,
                        bitscore=bitscore,
                        selfblast_scores=selfblast_scores,
                        default=1.0,
                    )

                    # If you prefer to use an existing BSR column when selfblast_scores missing:
                    if (not selfblast_scores) and len(columns) > 7:
                        try:
                            bsr = float(columns[7])
                        except ValueError:
                            pass

                    hit_bitscore = int(bitscore)
                    hsp_ident = int(pident)

                    if hit_protein_id in protein_dict:
                        protein_dict[hit_protein_id].add_domain(
                            hmm=query_id,
                            start=sstart,
                            end=send,
                            score=hit_bitscore,
                            ident=hsp_ident,
                            bsr=bsr,
                        )
                    else:
                        protein_dict[hit_protein_id] = Protein(
                            protein_id=hit_protein_id,
                            hmm=query_id,
                            start=sstart,
                            end=send,
                            score=hit_bitscore,
                            genome_id=genome_id,
                            ident=hsp_ident,
                            bsr=bsr,
                        )

                except ValueError as ve:
                    logger.warning(
                        f"Skipped malformed line in {filepath}: {line} (ValueError: {ve})"
                    )

    except Exception as e:
        logger.error(f"Failed to parse {filepath} for genome {genome_id}: {e}")
        logger.debug(traceback.format_exc())

    return protein_dict