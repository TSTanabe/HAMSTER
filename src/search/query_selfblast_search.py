#!/usr/bin/python
import os
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
    import os

    output_path = os.path.join(result_dir, "query_numbered.faa")

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


def create_selfquery_file(
    query_file: str, result_dir: str, prefix: str = "QUERY"
) -> tuple[str, str]:
    """
    Create two FASTA files derived from a query FASTA:

    - self_blast.faa : headers get prefix + unique counter
    - self_seqs.faa  : same headers but without prefix

    Parameters
    ----------
    query_file : str
        Path to input FASTA with query sequences.
    result_dir : str
        Output directory.
    prefix : str, optional
        Prefix for self_blast file (default "QUERY___").

    Returns
    -------
    (self_blast_path, self_seqs_path)
    """

    # Output paths
    self_query_path = os.path.join(result_dir, "self_query.faa")
    self_target_path = os.path.join(result_dir, "self_target.faa")

    id_count = {}

    with (
        open(self_query_path, "w") as outfile_query,
        open(self_target_path, "w") as outfile_seqs,
        open(query_file, "r") as infile,
    ):
        header = None
        sequence = ""

        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header and sequence:
                    original_id = header[1:].split()[0]
                    id_count[original_id] = id_count.get(original_id, 0) + 1
                    count = id_count[original_id]

                    outfile_query.write(
                        f">{prefix}___{original_id}___{count}\n{sequence}\n"
                    )
                    outfile_seqs.write(f">{original_id}___{count}\n{sequence}\n")

                header = line
                sequence = ""
            else:
                sequence += line

        # write last record
        if header and sequence:
            original_id = header[1:].split()[0]
            id_count[original_id] = id_count.get(original_id, 0) + 1
            count = id_count[original_id]

            outfile_query.write(f">{prefix}{original_id}___{count}\n{sequence}\n")
            outfile_seqs.write(f">{original_id}___{count}\n{sequence}\n")

    return self_query_path, self_target_path


def get_sequence_hits_scores(blast_file: str) -> Dict[str, float]:
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


def get_sequence_legth(file_path: str) -> Dict[str, float]:
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


def self_blast_query(options: Any) -> tuple[Dict[str, float], Dict[str, float]]:
    """
    Self-BLAST der (bereits durchnummerierten) Query-FASTA gegen sich selbst.

    Liefert:
      - selfblast_scores_dict: qseqid -> self-Bitscore
      - query_length_dict:     qseqid -> Sequenzlänge
    """

    query_fasta = options.query_file  # das ist jetzt query_numbered.faa

    # Self-BLAST: Query und Target sind gleich
    report = diamond_search.diamond_search(
        query_fasta,  # path = Target-DB
        query_fasta,  # query_fasta
        options.cores,
        options.evalue,
        1.0,  # coverage 100%
        100.0,  # minseqid 100% (oder 1.0, je nach Skala)
        options.diamond_report_hits_limit,
        options.alignment_mode,
    )

    ## These are for the return
    # Self-BLAST-Scores pro qseqid
    selfblast_scores_dict = get_sequence_hits_scores(report)
    # Sequenzlängen pro qseqid (muss an neues Header-Schema angepasst werden)
    query_length_dict = get_sequence_legth(query_fasta)

    ## Save in the database
    protein_dict = parse_reports.parse_bulk_blastreport_consecutive(
        genome_id="QUERY", filepath=report
    )
    print(protein_dict)
    parse_reports.get_protein_sequence(query_fasta, protein_dict)
    protein_dict = parse_reports.clean_dict_keys_and_protein_ids(protein_dict, "QUERY")
    database.insert_database_genome_ids(
        options.database_directory, genome_ids={"QUERY"}
    )
    print(protein_dict)
    database.insert_database_proteins(options.database_directory, protein_dict)
    return selfblast_scores_dict, query_length_dict
