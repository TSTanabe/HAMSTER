#!/usr/bin/python
import os
import sys
import traceback
from typing import Dict, Any

from src.search import diamond_search
from src.parse_reports import parse_reports

from src.db import database

from src.core.logging import get_logger

logger = get_logger(__name__)


def make_numbered_query_fasta(input_query: str, result_dir: str) -> str:
    """
    Nimmt eine Query-FASTA und erzeugt eine neue FASTA,
    in der jede Sequenz einen eindeutigen Header ID___<laufende Zahl> hat.

    Beispiel:
        >SmoC         -> >SmoC___1
        >SmoC         -> >SmoC___2
    """

    output_path = os.path.join(result_dir, "query_numbered.faa")
    if os.path.isfile(output_path) and os.path.getsize(output_path) > 0:
        return output_path

    if input_query is None or os.path.getsize(input_query) == 0:
        logger.error("Missing query file. Please use '-q' argument")
        sys.exit()

    counter = 0
    with open(input_query, "r") as infile, open(output_path, "w") as outfile:
        header = None
        seq = []

        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                # alten Record schreiben
                if header and seq:
                    counter += 1
                    base_id = header[1:].split()[0]
                    new_id = f"{base_id}___{counter}"
                    outfile.write(f">{new_id}\n{''.join(seq)}\n")
                header = line
                seq = []
            else:
                seq.append(line)

        # letzten Record schreiben
        if header and seq:
            counter += 1
            base_id = header[1:].split()[0]
            new_id = f"{base_id}___{counter}"
            outfile.write(f">{new_id}\n{''.join(seq)}\n")

    return output_path


def _clean_dict_keys_and_protein_ids(
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


def _strip_query_suffix(query_id: str) -> str:
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


def _parse_self_blastreport(
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

    protein_dict: Dict[str, parse_reports.Protein] = {}

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
                    query_id = _strip_query_suffix(raw_query_id)
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
                        protein_dict[hit_protein_id] = parse_reports.Protein(
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


def _get_sequence_hits_scores(blast_file: str) -> Dict[str, float]:
    """
    Generates a dict of self-blast scores from a BLAST table file.

    Args:
        blast_file (str): Path to BLASTP table file

    Returns:
        dict: qseqid -> highest bitscore

    Output Example:
        {"Q12345": 180.0, ...}
    """
    selfblast_scores = {}

    with open(blast_file, "r") as infile:
        for line in infile:
            row = line.strip().split("\t")
            sseqid, qseqid, evalue, bitscore, sstart, send, pident = row

            # Convert bitscore to float for comparison
            bitscore = float(bitscore)

            # Update the dictionary with the highest bitscore for self-hits
            if qseqid in selfblast_scores:
                selfblast_scores[qseqid] = max(selfblast_scores[qseqid], bitscore)
            else:
                selfblast_scores[qseqid] = bitscore

    return selfblast_scores


def _get_sequence_legth(file_path: str) -> Dict[str, float]:
    """
    Gibt die Sequenzlänge pro Query-ID zurück.

    Bei durchnummerierten Queries ist jede ID eindeutig
    (z.B. SmoC___1, SmoC___2)
    """
    lengths: Dict[str, float] = {}

    with open(file_path, "r") as fasta_file:
        current_id = None
        current_seq = []

        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = len("".join(current_seq))
                # komplette ID bis zum ersten Leerzeichen
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            lengths[current_id] = len("".join(current_seq))

    return lengths


def self_blast_query(config: Any) -> tuple[Dict[str, float], Dict[str, float]]:
    """
    Self-BLAST der (bereits durchnummerierten) Query-FASTA gegen sich selbst.

    Liefert:
      - selfblast_scores_dict: qseqid -> self-Bitscore
      - query_length_dict:     qseqid -> Sequenzlänge
    """

    query_fasta = config.query_file  # das ist jetzt query_numbered.faa

    # Self-BLAST: Query und Target sind gleich
    report = diamond_search.diamond_search(
        query_fasta,  # path = Target-DB
        query_fasta,  # query_fasta
        config.cores,
        config.evalue,
        1.0,  # coverage 100%
        100.0,  # minseqid 100% (oder 1.0, je nach Skala)
        config.diamond_report_hits_limit,
        config.alignment_mode,
    )

    ## These are for the return
    # Self-BLAST-Scores pro qseqid
    selfblast_scores_dict = _get_sequence_hits_scores(report)
    # Sequenzlängen pro qseqid (muss an neues Header-Schema angepasst werden)
    query_length_dict = _get_sequence_legth(query_fasta)

    ## Save in the database
    protein_dict = _parse_self_blastreport(genome_id="QUERY", filepath=report)

    parse_reports.get_protein_sequence(query_fasta, protein_dict)
    protein_dict = _clean_dict_keys_and_protein_ids(protein_dict, "QUERY")
    database.insert_database_genome_ids(config.database_directory, genome_ids={"QUERY"})

    database.insert_database_proteins(config.database_directory, protein_dict)
    return selfblast_scores_dict, query_length_dict
