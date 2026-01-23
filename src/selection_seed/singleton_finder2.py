#!/usr/bin/python

import copy
import os
import sqlite3
import pandas as pd
from collections import defaultdict
from typing import Dict, Any, Set, Tuple

from src.selection_defragmentation import pam_mx_algorithm
from src.core import myUtil

logger = myUtil.logger

"""
I want a new selection process for singletons. First the ones without a context must be identified

For this task take all hits with >80 % identity to the query and see if there are some examples without clusterID

Then take these genomes and take all proteins that have an identity of 80% or higher (adjustable by input)

Consider these as co occurring patterns. Select proteins that have the same co occurence pattern and at least 70 % 
identity to query. Do this for each protein.
"""


########################################
# New co-occurrence based singleton selection
########################################


def _find_context_free_high_identity_hits(
    database_path: str,
    identity_cutoff: float = 80.0,
) -> dict[str, set[str]]:
    """
    Find for each domain the genomes where identity >= cutoff
    AND clusterID IS NULL
    AND genomeID is not QUERY.

    Returns:
        dict:
            {
                domain: {genomeID1, genomeID2, ...},
                ...
            }
    """
    context_free: dict[str, set[str]] = defaultdict(set)

    query = """
        SELECT DISTINCT d.domain, p.genomeID
        FROM Domains d
        JOIN Proteins p ON p.proteinID = d.proteinID
        WHERE d.identity >= ?
          AND p.clusterID IS NULL
          AND p.genomeID != 'QUERY'
    """

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()
        for domain, genome_id in cur.execute(query, (identity_cutoff,)):
            context_free[domain].add(genome_id)

    return context_free


def _fetch_high_identity_domain_intersection(
    database_path: str,
    domain_to_genomes: Dict[str, Set[str]],
    min_identity_cutoff: float = 70.0,
) -> Dict[str, Set[str]]:
    """
    Für jede gegebene Domain und deren zugehörige Genome:
      - hole alle Domaintypen mit identity >= min_identity_cutoff in diesen Genomen
      - berechne die Schnittmenge der Domain-Sets über die Genome dieser Domain

    Args:
        database_path: Pfad zur SQLite-Datenbank.
        domain_to_genomes: dict {seed_domain: set(genomeIDs)},
                           z.B. Ausgabe von _find_context_free_high_identity_hits().
        min_identity_cutoff: minimale Identity in Prozent.

    Returns:
        dict:
            {
                seed_domain: {domain1, domain2, ...},  # Domains, die in ALLEN zugehörigen Genomen
                                                     # mit identity >= cutoff vorkommen
                ...
            }
    """

    # Ergebnis: für jede Seed-Domain die Schnittmenge der Co-Occurrence-Domains
    intersections_per_seed: Dict[str, Set[str]] = {}

    if not domain_to_genomes:
        return intersections_per_seed

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()
        chunk_size = 900  # wegen SQLite-Placeholder-Limit

        for seed_domain, genomes in domain_to_genomes.items():
            genomes_list = [g for g in genomes if g]
            if not genomes_list:
                intersections_per_seed[seed_domain] = set()
                continue

            # pro Genom: welche Domains erfüllen den Identity-Cutoff?
            domains_per_genome: Dict[str, Set[str]] = defaultdict(set)

            # DB in Chunks abfragen
            for i in range(0, len(genomes_list), chunk_size):
                chunk = genomes_list[i : i + chunk_size]
                placeholders = ",".join(["?"] * len(chunk))

                query = f"""
                    SELECT DISTINCT d.domain, p.genomeID
                    FROM Domains d
                    JOIN Proteins p ON p.proteinID = d.proteinID
                    WHERE d.identity >= ?
                      AND p.genomeID IN ({placeholders})
                      AND p.genomeID != 'QUERY'
                """
                cur.execute(query, (min_identity_cutoff, *chunk))
                for domain, genome_id in cur.fetchall():
                    domains_per_genome[genome_id].add(domain)

            # Schnittmenge über alle Genome dieser Seed-Domain bilden
            intersection: Set[str] | None = None
            for genome_id in genomes_list:
                doms = domains_per_genome.get(genome_id, set())
                if intersection is None:
                    intersection = set(doms)
                else:
                    intersection &= doms

                # Early exit: wenn leer, kann nichts mehr dazu kommen
                if not intersection:
                    break

            intersections_per_seed[seed_domain] = intersection or set()

    return intersections_per_seed


def select_singleton_refs_by_domain_pattern(
    database_path: str,
    seed_to_pattern_domains: Dict[str, Set[str]],
    min_identity_cutoff: float = 70.0,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Set[str]]]:
    limits_dict: Dict[str, Dict[str, float]] = {}
    sng_reference_seq_dict: Dict[str, Set[str]] = defaultdict(set)

    if not seed_to_pattern_domains:
        return limits_dict, sng_reference_seq_dict

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        # TEMP table once
        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_pattern_domains (domain TEXT PRIMARY KEY)"
        )

        for seed_domain, pattern_domains in seed_to_pattern_domains.items():
            if not pattern_domains:
                continue

            # fill temp table
            cur.execute("DELETE FROM tmp_pattern_domains")
            cur.executemany(
                "INSERT OR IGNORE INTO tmp_pattern_domains(domain) VALUES (?);",
                ((d,) for d in pattern_domains),
            )

            # 1) Genome that have ALL pattern domains
            query_genomes = """
                SELECT p.genomeID
                FROM Domains d
                JOIN Proteins p ON p.proteinID = d.proteinID
                JOIN tmp_pattern_domains t ON t.domain = d.domain
                WHERE d.identity >= ?
                  AND p.genomeID != 'QUERY'
                GROUP BY p.genomeID
                HAVING COUNT(DISTINCT d.domain) = (SELECT COUNT(*) FROM tmp_pattern_domains)
            """

            cur.execute(query_genomes, (min_identity_cutoff,))
            genomes = [row[0] for row in cur.fetchall()]

            if not genomes:
                continue

            # 2) fetch ALL seed_domain hits in these genomes
            all_proteinIDs = set()
            all_bitscores = []

            chunk_size = 900
            for i in range(0, len(genomes), chunk_size):
                chunk = genomes[i : i + chunk_size]
                placeholders = ",".join(["?"] * len(chunk))

                query_hits = f"""
                    SELECT d.proteinID, d.score
                    FROM Domains d
                    JOIN Proteins p ON p.proteinID = d.proteinID
                    WHERE d.domain = ?
                      AND d.identity >= ?
                      AND p.genomeID IN ({placeholders})
                      AND p.genomeID != 'QUERY'
                """
                cur.execute(query_hits, (seed_domain, min_identity_cutoff, *chunk))
                for protein_id, score in cur.fetchall():
                    all_proteinIDs.add(protein_id)
                    try:
                        all_bitscores.append(float(score))
                    except:
                        pass

            if not all_proteinIDs or not all_bitscores:
                continue

            sng_reference_seq_dict[seed_domain] = all_proteinIDs
            limits_dict[seed_domain] = {
                "lower_limit": min(all_bitscores),
                "upper_limit": max(all_bitscores),
            }

    return limits_dict, sng_reference_seq_dict


#### Main routine of this module
def singleton_reference_finder(
    options: Any,
) -> (
    tuple[object , object]
    | tuple[dict[str, dict[str, float]], dict[str, set[str]]]
):
    """
    Main routine: finds and predicts reference singletons for each protein/domain.

    Args:
        options (Any): Configuration/options.

    Returns:
        tuple:
            - domain_score_limits (dict): {domain: {'lower_limit', 'upper_limit'}}
            - singleton_reference_seqs_dict (dict): {domain: set(proteinIDs)}
    """

    high_identity_cutoff = getattr(options, "singleton_identity_cutoff", 70.0)
    pattern_identity_cutoff = getattr(
        options, "singleton_pattern_identity_cutoff", 70.0
    )

    # 1) Find genomes & domains with context-free high-identity hits
    context_free_domains_dict = _find_context_free_high_identity_hits(
        options.database_directory,
        identity_cutoff=high_identity_cutoff,
    )

    logger.info(f"Found {len(context_free_domains_dict)} context-free hits with {high_identity_cutoff} percent identity")

    if not context_free_domains_dict:
        logger.warning(
            "No context-free high-identity hits found – stopping singleton selection."
        )
        return {}, {}

    # 2) For these genomes, fetch all hits with identity >= pattern_identity_cutoff
    domain_presence_intersection_pattern = _fetch_high_identity_domain_intersection(
        database_path=options.database_directory,
        domain_to_genomes=context_free_domains_dict,
        min_identity_cutoff=pattern_identity_cutoff,  # identity cutoff for co occurring pattern domains
    )

    if not domain_presence_intersection_pattern:
        logger.warning(
            f"No presence/absence pattern detected for hits with {high_identity_cutoff} percent identity – stopping singleton selection."
        )
        return {}, {}

    # 3) Select co-occurrence-based singleton candidates
    domain_score_limits, singleton_reference_seqs_dict = (
        select_singleton_refs_by_domain_pattern(
            database_path=options.database_directory,
            seed_to_pattern_domains=domain_presence_intersection_pattern,
            min_identity_cutoff=pattern_identity_cutoff,
        )
    )

    myUtil.save_cache(
        options, "sng0_training_proteinIDs.pkl", singleton_reference_seqs_dict
    )
    myUtil.save_cache(
        options, "sng0_training_proteinIDs_limits.pkl", domain_score_limits
    )

    return domain_score_limits, singleton_reference_seqs_dict
