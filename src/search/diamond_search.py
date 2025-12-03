#!/usr/bin/python
import os
from src.core import myUtil
from src.core import logging

logger = logging.get_logger(__name__)


def diamond_search(
    path: str,
    query_fasta: str,
    cores: int,
    evalue: float,
    coverage: float,
    minseqid: float,
    diamond_report_hits_limit: int,
    alignment_mode: int = 2,
    sensitivity: str = "ultra-sensitive",
) -> str:
    """
    Runs DIAMOND BLASTp search on a FASTA file.

    Args:
        path (str): Path to the FASTA file (target db)
        query_fasta (str): Query FASTA
        cores (int): Number of cores
        evalue (float): E-value cutoff
        coverage (float): Min coverage
        minseqid (float): Min sequence identity
        diamond_report_hits_limit (int): Report hits limit
        alignment_mode (int): DIAMOND alignment mode
        sensitivity (str): Sensitivity level

    Returns:
        str: Path to output .diamond.tab file

    Example Output:
        "results/glob.faa.diamond.tab"
    """

    diamond = myUtil.find_executable("diamond")
    target_db_name = f"{path}.dmnd"
    logger.debug(
        f"{diamond} makedb --quiet --in {path} -d {target_db_name} --threads {cores} 1>/dev/null 0>/dev/null"
    )
    os.system(
        f"{diamond} makedb --quiet --in {path} -d {target_db_name} --threads {cores} 1>/dev/null 0>/dev/null"
    )

    output_results_tab = f"{path}.diamond.tab"
    # {hit_id}\t{query_id}\t{e_value}\t{score}\t{bias}\t{hsp_start}\t{hsp_end}\t{description}

    logger.debug(
        f"{diamond} blastp --quiet --{sensitivity} -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} --outfmt 6 sseqid qseqid evalue bitscore sstart send pident 1>/dev/null 0>/dev/null"
    )
    os.system(
        f"{diamond} blastp --quiet --{sensitivity} -d {target_db_name} -q {query_fasta} -o {output_results_tab} --threads {cores} -e {evalue} -k {diamond_report_hits_limit} --outfmt 6 sseqid qseqid evalue bitscore sstart send pident 1>/dev/null 0>/dev/null"
    )
    # output format hit query evalue score identity alifrom alito
    return output_results_tab
