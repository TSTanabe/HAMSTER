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

    Uses TEMP table to avoid repeated Python-side bookkeeping and can be extended
    later to reuse tmp_cf across subsequent steps.
    """
    context_free: dict[str, set[str]] = defaultdict(set)

    with sqlite3.connect(database_path, timeout=120.0) as con:
        cur = con.cursor()

        # TEMP is a write -> allow it, then lock down persistent DB afterwards if you want
        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")   # ~256 MiB
        cur.execute("PRAGMA mmap_size=2147483648;") # 2 GiB
        cur.execute("PRAGMA automatic_index=ON;")

        cur.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_cf (domain TEXT, genomeID TEXT);")
        cur.execute("DELETE FROM tmp_cf;")

        # Fill TEMP table
        cur.execute(
            """
            INSERT INTO tmp_cf(domain, genomeID)
            SELECT DISTINCT d.domain, p.genomeID
            FROM Domains d
            JOIN Proteins p ON p.proteinID = d.proteinID
            WHERE d.identity >= ?
              AND p.clusterID IS NULL
              AND p.genomeID != 'QUERY'
            """,
            (identity_cutoff,),
        )

        # Now (optional) protect persistent schema (TEMP already populated)
        cur.execute("PRAGMA query_only=TRUE;")

        # Read back
        for domain, genome_id in cur.execute("SELECT domain, genomeID FROM tmp_cf;"):
            context_free[domain].add(genome_id)

    return context_free



def _fetch_high_identity_domain_intersection(
    database_path: str,
    domain_to_genomes: Dict[str, Set[str]],
    min_identity_cutoff: float = 70.0,
) -> Dict[str, Set[str]]:
    """
    For each seed_domain and its genomes:
      - compute the intersection of domains with identity >= cutoff across those genomes
    using TEMP tables and a single GROUP BY/HAVING query per seed_domain.

    This replaces the chunked IN() queries and Python-side set intersection.
    """
    intersections_per_seed: Dict[str, Set[str]] = {}
    if not domain_to_genomes:
        return intersections_per_seed

    with sqlite3.connect(database_path, timeout=120.0) as con:
        cur = con.cursor()

        # TEMP allowed; then lock down main schema after filling temp each iteration if desired
        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")   # ~256 MiB
        cur.execute("PRAGMA mmap_size=2147483648;") # 2 GiB
        cur.execute("PRAGMA automatic_index=ON;")

        # Create TEMP table once
        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_seed_genomes (genomeID TEXT PRIMARY KEY);"
        )

        for seed_domain, genomes in domain_to_genomes.items():
            logger.info(f"Fetching data for {seed_domain} in {len(genomes)} genomes")
            genomes_list = [g for g in genomes if g and g != "QUERY"]
            if not genomes_list:
                intersections_per_seed[seed_domain] = set()
                continue

            # Fill tmp_seed_genomes for this seed
            cur.execute("DELETE FROM tmp_seed_genomes;")
            cur.executemany(
                "INSERT OR IGNORE INTO tmp_seed_genomes(genomeID) VALUES (?);",
                ((g,) for g in genomes_list),
            )

            # Optional: lock down persistent db after TEMP is ready
            cur.execute("PRAGMA query_only=TRUE;")

            # N genomes in the current seed set
            cur.execute("SELECT COUNT(*) FROM tmp_seed_genomes;")
            n_genomes = cur.fetchone()[0]
            if n_genomes == 0:
                intersections_per_seed[seed_domain] = set()
                cur.execute("PRAGMA query_only=FALSE;")
                continue

            # Intersection: domains that appear in ALL genomes (identity cutoff)
            # Key trick: GROUP BY domain + HAVING COUNT(DISTINCT genomeID) = n_genomes
            cur.execute(
                """
                SELECT d.domain
                FROM Domains d
                JOIN Proteins p ON p.proteinID = d.proteinID
                JOIN tmp_seed_genomes g ON g.genomeID = p.genomeID
                WHERE d.identity >= ?
                  AND p.genomeID != 'QUERY'
                GROUP BY d.domain
                HAVING COUNT(DISTINCT p.genomeID) = ?
                """,
                (min_identity_cutoff, n_genomes),
            )

            intersections_per_seed[seed_domain] = {row[0] for row in cur.fetchall()}

            # Allow next TEMP refill (some SQLite builds block TEMP writes under query_only)
            cur.execute("PRAGMA query_only=FALSE;")

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

    with sqlite3.connect(database_path, timeout=120.0) as con:
        cur = con.cursor()

        # Pragmas (TEMP muss erlaubt sein)
        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")      # ~256 MiB
        cur.execute("PRAGMA mmap_size=2147483648;")    # 2 GiB
        cur.execute("PRAGMA automatic_index=ON;")

        # TEMP tables einmalig
        cur.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_pattern_domains (domain TEXT PRIMARY KEY);")
        cur.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_genomes (genomeID TEXT PRIMARY KEY);")

        for seed_domain, pattern_domains in seed_to_pattern_domains.items():
            if not pattern_domains:
                continue

            logger.info(f"Fetching data for {seed_domain} with {len(pattern_domains)} pattern domains")

            # --- tmp_pattern_domains füllen ---
            cur.execute("DELETE FROM tmp_pattern_domains;")
            cur.executemany(
                "INSERT OR IGNORE INTO tmp_pattern_domains(domain) VALUES (?);",
                ((d,) for d in pattern_domains if d),
            )

            # --- passende Genomes direkt in tmp_genomes materialisieren ---
            cur.execute("DELETE FROM tmp_genomes;")
            cur.execute(
                """
                INSERT OR IGNORE INTO tmp_genomes(genomeID)
                SELECT p.genomeID
                FROM Domains d
                JOIN Proteins p ON p.proteinID = d.proteinID
                JOIN tmp_pattern_domains t ON t.domain = d.domain
                WHERE d.identity >= ?
                  AND p.genomeID != 'QUERY'
                GROUP BY p.genomeID
                HAVING COUNT(DISTINCT d.domain) = (SELECT COUNT(*) FROM tmp_pattern_domains)
                """,
                (min_identity_cutoff,),
            )

            # Schnell abbrechen, wenn keine Genomes
            cur.execute("SELECT COUNT(*) FROM tmp_genomes;")
            n_genomes = cur.fetchone()[0]
            if n_genomes == 0:
                continue

            # --- Seed-domain Hits in diesen Genomes: 1 Query, kein IN-chunking ---
            # 1) Limits (MIN/MAX) direkt in SQL (kein Python-List-Aufbau)
            cur.execute(
                """
                SELECT MIN(d.score), MAX(d.score)
                FROM Domains d
                JOIN Proteins p ON p.proteinID = d.proteinID
                JOIN tmp_genomes g ON g.genomeID = p.genomeID
                WHERE d.domain = ?
                  AND d.identity >= ?
                  AND p.genomeID != 'QUERY'
                """,
                (seed_domain, min_identity_cutoff),
            )
            row = cur.fetchone()
            if not row or row[0] is None or row[1] is None:
                continue
            lower, upper = float(row[0]), float(row[1])

            # 2) ProteinIDs holen (DISTINCT)
            cur.execute(
                """
                SELECT DISTINCT d.proteinID
                FROM Domains d
                JOIN Proteins p ON p.proteinID = d.proteinID
                JOIN tmp_genomes g ON g.genomeID = p.genomeID
                WHERE d.domain = ?
                  AND d.identity >= ?
                  AND p.genomeID != 'QUERY'
                """,
                (seed_domain, min_identity_cutoff),
            )
            protein_ids = {r[0] for r in cur.fetchall()}
            if not protein_ids:
                continue

            sng_reference_seq_dict[seed_domain] = protein_ids
            limits_dict[seed_domain] = {"lower_limit": lower, "upper_limit": upper}

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
    logger.info("Fetching co-occurence-based singleton candidates")
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
