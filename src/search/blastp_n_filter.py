#!/usr/bin/python
import os
import csv

from typing import Dict
from src.search import diamond_search
from src.core.logging import get_logger
logger = get_logger(__name__)


def filter_blast_table(
    output_file: str,
    blast_file: str,
    evalue_cutoff: float,
    score_cutoff: float,
    coverage_cutoff: float,
    identity_cutoff: float,
    bsr_cutoff: float,
    sequence_lengths: Dict[str, float],
    selfblast_scores: Dict[str, float],
    buffer_size: int = 10000,
) -> str:
    """
    Filters a BLASTP table by e-value, score, coverage, identity, and BSR.

    Args:
        output_file: Path for output
        blast_file: Path to input
        evalue_cutoff: Max e-value
        score_cutoff: Min bitscore
        coverage_cutoff: Min coverage (0-1)
        identity_cutoff: Min identity (0-1)
        bsr_cutoff: Min Blast Score Ratio
        sequence_lengths: dict of qseqid -> seq len
        selfblast_scores: dict of qseqid -> self-hit bitscore
        buffer_size: Write buffer

    Returns:
        str: Output file path

    Output Example:
        "results/filtered_glob.faa.diamond.tab"
    """

    # Determine the output file path
    if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
        logger.info(
            f"Filtered hit results file already exists and is non-empty: {output_file}"
        )
        return output_file

    # Open input and output files using CSV reader/writer
    with open(blast_file, "r") as infile, open(output_file, "w", newline="") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        buffer = []  # Buffer to store valid rows

        # Process each row efficiently
        for row in reader:
            try:
                sseqid, qseqid, evalue, bitscore, sstart, send, pident = row[:7]

                # Convert necessary values only when needed
                evalue = float(evalue)
                bitscore = float(bitscore)
                pident = float(pident)
                sstart, send = int(sstart), int(send)

                # Fetch precomputed query length and self-blast score
                query_length = sequence_lengths.get(qseqid)
                selfblast_score = selfblast_scores.get(qseqid)

                # Pre-filter based on missing values
                if not query_length or not selfblast_score:
                    buffer.append(row)
                else:
                    # Compute alignment length and coverage
                    alignment_length = abs(send - sstart) + 1
                    coverage = alignment_length / query_length

                    # Compute Blast Score Ratio (BSR)
                    bsr = bitscore / selfblast_score

                    # Apply all filtering criteria
                    if ( # If (bitscore or evalue) and coverage ok
                        evalue <= evalue_cutoff
                        and bitscore >= score_cutoff
                        and coverage >= coverage_cutoff
                    ):
                        row.append(f"{bsr:.3f}")
                        buffer.append(row)
                    elif ( # If identity and coverage ok
                        coverage >= coverage_cutoff
                        and pident >= identity_cutoff
                    ):
                        row.append(f"{bsr:.3f}")
                        buffer.append(row)
                # **Write buffer to disk when it reaches buffer_size**
                if len(buffer) >= buffer_size:
                    writer.writerows(buffer)
                    buffer.clear()  # Reset buffer

            except ValueError as ve:
                logger.warning(f"Skipping malformed row: {row} (ValueError: {ve})")
                continue  # Skip invalid rows gracefully

        # Final flush: Write any remaining data in the buffer
        if buffer:
            writer.writerows(buffer)

    return output_file


def run_and_filter_diamond_blastp(
    options,
    *,
    query_length_dict: Dict[str, float],
    selfblast_scores_dict: Dict[str, float],
) -> str:
    """
    Ensure a DIAMOND BLASTp table exists (run it if needed), then filter it and
    store the filtered path in options.glob_table.

    Returns:
        Path to the filtered blast table.
    """

    # Step 2: obtain raw DIAMOND table (either precomputed or freshly generated)
    blast_results_table = options.glob_table
    if not blast_results_table:
        logger.info("Initialize DIAMOND BLASTp against target")
        blast_results_table = diamond_search.diamond_search(
            options.glob_faa,
            options.query_file,
            options.cores,
            options.evalue,
            options.searchcoverage,
            options.minseqid,
            options.diamond_report_hits_limit,
            options.alignment_mode,
        )
    else:
        logger.info(f"Using DIAMOND BLASTp result table {blast_results_table}")

    # Step 3: filter raw table
    logger.info(
        "Filtering raw DIAMOND BLASTp results by score, e-value, identity, coverage and BSR"
    )
    filtered_table = filter_blast_table(
        output_file=options.filtered_blast_table,
        blast_file=blast_results_table,
        evalue_cutoff=options.evalue,
        score_cutoff=options.thrs_score,
        coverage_cutoff=options.searchcoverage,
        identity_cutoff=options.minseqid,
        bsr_cutoff=options.thrs_bsr,
        sequence_lengths=query_length_dict,
        selfblast_scores=selfblast_scores_dict,
    )

    options.glob_table = filtered_table
    return filtered_table