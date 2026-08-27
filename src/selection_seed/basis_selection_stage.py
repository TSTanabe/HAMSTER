#!/usr/bin/python
import sqlite3
import sys
from typing import Dict, Set

from src.db import database
from src.core import myUtil
from src.selection_seed import (
    csb_proteins_selection,
    singleton_finder2,
)
from src.core.logging import get_logger

logger = get_logger(__name__)

def _add_query_ids_to_proteinIDset(
    combined_protein_sets: Dict[str, Set[str]], database_path: str
) -> Dict[str, Set[str]]:
    """
    Adds any proteinIDs from the query cluster to each protein group set.

    Args:

        database_path (str): Path to SQLite DB.

    Returns:
        dict: {group: set(proteinIDs)} (now includes query IDs)
    """
    # Connect to the database
    with sqlite3.connect(database_path, timeout=120.0) as conn:
        cursor = conn.cursor()

        # **Schritt 1: Finde alle proteinIDs mit genomeID = 'QUERY'**
        cursor.execute("SELECT proteinID FROM Proteins WHERE genomeID = 'QUERY'")
        query_protein_ids = {
            row[0] for row in cursor.fetchall()
        }  # Set für schnellere Suche

        if not query_protein_ids:
            return combined_protein_sets  # Falls leer, sofort zurückgeben

        for key in combined_protein_sets:
            # **Schritt 2: Hole proteinIDs aus Domains, aber nur, wenn sie in query_protein_ids sind**
            cursor.execute(
                f"""
                SELECT proteinID FROM Domains 
                WHERE domain = ? AND proteinID IN ({",".join(["?"] * len(query_protein_ids))})
                """,
                (key, *query_protein_ids),
            )

            # Add fetched proteinIDs to the protein set
            protein_ids = {row[0] for row in cursor.fetchall()}
            combined_protein_sets[key].update(protein_ids)

    return combined_protein_sets


def basis_sequence_fasta(config) -> None:
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
        config (Options): Pipeline options

    Output:
        - FASTA files for the 'grp0' (basis) dataset.
        - Updates options.grouped and options.score_limit_dict with results.

    Example Output:
        options.grouped: dict[str, set[str]] = {'domainA': {'seq1', 'seq2', ...}, ...}
        options.score_limit_dict: dict[str, dict[str, float]] = {'domainA': {'lower_limit': 10.0, 'upper_limit': 200.0}}
    """
    logger.debug("Indexing local database")
    database.index_database(config.database_directory)

    ### Collect the sequences from csb where at least one query is encoded
    grp_score_limit_dict, grouped = (
        csb_proteins_selection.prepare_csb_grouped_seed_proteins(config)
    )
    config.grouped = grouped

    ### Collect singletons without any conserved csb
    logger.info(
        "Collecting highly similar homologs from query hits without any conserved genomic context"
    )
    sng_score_limit_dict, sng_ref_seqs_dict = (
        singleton_finder2.prepare_singleton_seed_proteins(config)
    )

    # Merge groups and limits from csb and sng and add queries
    merged_score_limit_dict = {**grp_score_limit_dict, **sng_score_limit_dict}
    merged_basis_seed_proteins_dict = csb_proteins_selection.merge_protein_sets(grouped, sng_ref_seqs_dict)

    if not merged_basis_seed_proteins_dict:
        logger.error(
            "There were no proteins selected for basic training. Consider increasing the input data or lowering selection thresholds - stopping execution"
        )
    merged_basis_seed_proteins_dict = _add_query_ids_to_proteinIDset(
        merged_basis_seed_proteins_dict, config.database_directory
    )

    # Print the grp0 csb and singletons to fasta
    csb_proteins_selection.fetch_training_data_to_fasta(config, merged_basis_seed_proteins_dict, "ds1")

    # Save the merged groups
    myUtil.save_cache(config, "basis_merged_grouped.pkl", merged_basis_seed_proteins_dict)
    myUtil.save_cache(config, "basis_merged_score.pkl", merged_score_limit_dict)

    # Update options object with the fetched proteinID groups and score limits
    config.grouped = merged_basis_seed_proteins_dict
    config.score_limit_dict = merged_score_limit_dict

    return
