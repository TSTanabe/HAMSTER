#!/usr/bin/python
import sys

from src.db import database
from src.core import myUtil
from src.selection_seed import (
    csb_proteins_selection,
    singleton_finder2,
)
from src.core.logging import get_logger

logger = get_logger(__name__)


def basis_sequence_fasta(options) -> None:
    """
        prepare the sequence fasta and identifier lists. Seqs are derived from
        named gene clusters and csb finder gene clusters that encode a protein sequence
        highly similar to the input reference sequence. The hit score for the inclusion of
        sequences has to be 30 % of the best scoring query, which is similar to a blast score ration
        of 0.7 and the conserved sequence identity of 70 %

        Input: Database, options
        Output are the grp0 fasta files

        Example:
        A score limiter dictionary
        singleton_score_limits[domain] = {
                    "lower_limit": min(bitscores),
                    "upper_limit": max(bitscores)
                }
    Args:
        options (Options): Pipeline options

    Output:
        - FASTA files for the 'grp0' (basis) dataset.
        - Updates options.grouped and options.score_limit_dict with results.

    Example Output:
        options.grouped: dict[str, set[str]] = {'domainA': {'seq1', 'seq2', ...}, ...}
        options.score_limit_dict: dict[str, dict[str, float]] = {'domainA': {'lower_limit': 10.0, 'upper_limit': 200.0}}
    """
    logger.debug("Indexing local database")
    database.index_database(options.database_directory)

    ### Collect the sequences from csb where at least one query is encoded
    grp_score_limit_dict, grouped = (
        csb_proteins_selection.prepare_csb_grouped_training_proteins(options)
    )

    ### Collect singletons without any conserved csb
    logger.info(
        "Collecting highly similar homologs from query hits without any conserved genomic context"
    )
    sng_score_limit_dict, sng_ref_seqs_dict = (
        singleton_finder2.singleton_reference_finder(options)
    )

    # Merge groups and limits from csb and sng
    if grouped or sng_ref_seqs_dict:
        merged_score_limit_dict = {**grp_score_limit_dict, **sng_score_limit_dict}
        merged_grouped = {**grouped, **sng_ref_seqs_dict}
    else:
        logger.error(
            "There were no proteins selected for basic training. Consider increasing the input data or lowering selection thresholds - stopping execution"
        )
        sys.exit()

    # Print the grp0 csb and singletons to fasta
    csb_proteins_selection.fetch_training_data_to_fasta(options, merged_grouped, "grp0")

    # Save the merged groups
    myUtil.save_cache(options, "basis_merged_grouped.pkl", merged_grouped)
    myUtil.save_cache(options, "basis_merged_score.pkl", merged_score_limit_dict)

    # Update options object with the fetched proteinID groups and score limits
    options.grouped = merged_grouped
    options.score_limit_dict = merged_score_limit_dict

    return
