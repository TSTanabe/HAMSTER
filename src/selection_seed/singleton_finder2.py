#!/usr/bin/python

import copy
import math
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
        cur.execute("PRAGMA cache_size=-262144;")  # ~256 MiB
        cur.execute("PRAGMA mmap_size=2147483648;")  # 2 GiB
        cur.execute("PRAGMA automatic_index=ON;")

        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_cf (domain TEXT, genomeID TEXT);"
        )
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

def _add_missing_query_domains_to_context_free_dict(
    database_path: str,
    already_grouped_domains: set[str],
    blast_score_ratio_cutoff: float = 0.7,
) -> dict[str, set[str]]:
    """
    Find all QUERY domains not already covered by CSB-derived training data
    and collect non-QUERY genomes containing sufficiently similar homologs.

    Similarity is evaluated using the stored Blast Score Ratio (BSR).

    Args:
        database_path:
            Path to SQLite database.

        already_grouped_domains:
            Domains already represented by CSB-derived training data.

        blast_score_ratio_cutoff:
            Minimum Blast Score Ratio required for a homolog.

    Returns:
        {
            domain: {genomeID1, genomeID2, ...}
        }
    """

    missing_domains_dict: dict[str, set[str]] = defaultdict(set)

    with sqlite3.connect(database_path, timeout=120.0) as con:
        cur = con.cursor()

        # --------------------------------------------------------------
        # Temporary lookup table for already covered domains
        # --------------------------------------------------------------

        cur.execute(
            """
            CREATE TEMP TABLE IF NOT EXISTS tmp_grouped_domains (
                domain TEXT PRIMARY KEY
            ) WITHOUT ROWID;
            """
        )

        cur.execute("DELETE FROM tmp_grouped_domains;")

        grouped_domains = set(already_grouped_domains or set())

        if grouped_domains:
            cur.executemany(
                """
                INSERT OR IGNORE INTO tmp_grouped_domains(domain)
                VALUES (?);
                """,
                ((domain,) for domain in grouped_domains),
            )

        # --------------------------------------------------------------
        # Find missing QUERY domains and their homolog-containing genomes
        # directly in one SQL query.
        # --------------------------------------------------------------

        cur.execute(
            """
            WITH missing_query_domains AS (
                SELECT DISTINCT d.domain
                FROM Domains AS d
                JOIN Proteins AS p
                  ON p.proteinID = d.proteinID
                LEFT JOIN tmp_grouped_domains AS g
                  ON g.domain = d.domain
                WHERE p.genomeID = 'QUERY'
                  AND d.domain IS NOT NULL
                  AND g.domain IS NULL
            )
            SELECT DISTINCT
                m.domain,
                p.genomeID
            FROM missing_query_domains AS m
            JOIN Domains AS d
              ON d.domain = m.domain
            JOIN Proteins AS p
              ON p.proteinID = d.proteinID
            WHERE p.genomeID != 'QUERY'
              AND d.blast_score_ratio >= ?;
            """,
            (blast_score_ratio_cutoff,),
        )

        for domain, genome_id in cur:
            if domain and genome_id:
                missing_domains_dict[domain].add(genome_id)

    logger.info(
        f"Found {len(missing_domains_dict)} QUERY domains not covered by "
        f"CSB training with homologs at BSR >= {blast_score_ratio_cutoff}"
    )

    return missing_domains_dict


def _fetch_conserved_co_occurrence_pattern(
    database_path: str,
    domain_to_genomes: Dict[str, Set[str]],
    min_bsr_cutoff: float = 0.7,
    min_presence_fraction: float = 0.90,
) -> Dict[str, Set[str]]:
    """
    For each seed domain, identify domains that occur in at least a defined
    fraction of the associated genomes.

    Example:
        min_presence_fraction = 0.90

        A domain occurring in >= 90 % of the seed genomes is considered
        part of the conserved co-occurrence pattern.
    """

    patterns_per_seed: Dict[str, Set[str]] = {}

    if not domain_to_genomes:
        return patterns_per_seed

    with sqlite3.connect(database_path, timeout=120.0) as con:
        cur = con.cursor()

        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")
        cur.execute("PRAGMA mmap_size=2147483648;")
        cur.execute("PRAGMA automatic_index=ON;")

        cur.execute(
            """
            CREATE TEMP TABLE IF NOT EXISTS tmp_seed_genomes (
                genomeID TEXT PRIMARY KEY
            );
            """
        )

        for seed_domain, genomes in domain_to_genomes.items():

            genomes_list = [
                g for g in genomes
                if g and g != "QUERY"
            ]

            if not genomes_list:
                patterns_per_seed[seed_domain] = set()
                continue

            cur.execute("DELETE FROM tmp_seed_genomes;")

            cur.executemany(
                """
                INSERT OR IGNORE INTO tmp_seed_genomes(genomeID)
                VALUES (?);
                """,
                ((g,) for g in genomes_list),
            )

            cur.execute(
                "SELECT COUNT(*) FROM tmp_seed_genomes;"
            )
            n_genomes = cur.fetchone()[0]

            if n_genomes == 0:
                patterns_per_seed[seed_domain] = set()
                continue

            # Minimum number of genomes in which a domain must occur.
            min_genomes = math.ceil(
                n_genomes * min_presence_fraction
            )

            cur.execute(
                """
                SELECT d.domain
                FROM Domains AS d
                JOIN Proteins AS p
                  ON p.proteinID = d.proteinID
                JOIN tmp_seed_genomes AS g
                  ON g.genomeID = p.genomeID
                WHERE d.blast_score_ratio >= ?
                  AND p.genomeID != 'QUERY'
                GROUP BY d.domain
                HAVING COUNT(DISTINCT p.genomeID) >= ?;
                """,
                (
                    min_bsr_cutoff,
                    min_genomes,
                ),
            )

            patterns_per_seed[seed_domain] = {
                row[0]
                for row in cur
                if row[0]
            }

            logger.debug(
                f"Seed {seed_domain}: "
                f"{n_genomes} genomes, "
                f"pattern domains must occur in >= "
                f"{min_genomes} genomes "
                f"({min_presence_fraction:.0%})."
            )

    return patterns_per_seed


def select_singleton_refs_by_domain_pattern(
    database_path: str,
    seed_to_pattern_domains: Dict[str, Set[str]],
    min_bsr_cutoff: float = 70.0,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Set[str]]]:
    limits_dict: Dict[str, Dict[str, float]] = {}
    sng_reference_seq_dict: Dict[str, Set[str]] = defaultdict(set)

    if not seed_to_pattern_domains:
        return limits_dict, sng_reference_seq_dict

    with sqlite3.connect(database_path, timeout=120.0) as con:
        cur = con.cursor()

        # Pragmas (TEMP muss erlaubt sein)
        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")  # ~256 MiB
        cur.execute("PRAGMA mmap_size=2147483648;")  # 2 GiB
        cur.execute("PRAGMA automatic_index=ON;")

        # TEMP tables einmalig
        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_pattern_domains (domain TEXT PRIMARY KEY);"
        )
        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_genomes (genomeID TEXT PRIMARY KEY);"
        )

        for seed_domain, pattern_domains in seed_to_pattern_domains.items():
            if not pattern_domains:
                continue

            logger.debug(
                f"Fetching data for {seed_domain} with {len(pattern_domains)} pattern domains"
            )

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
                (min_bsr_cutoff,),
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
                (seed_domain, min_bsr_cutoff),
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
                (seed_domain, min_bsr_cutoff),
            )
            protein_ids = {r[0] for r in cur.fetchall()}
            if not protein_ids:
                continue

            sng_reference_seq_dict[seed_domain] = protein_ids
            limits_dict[seed_domain] = {"lower_limit": lower, "upper_limit": upper}

    return limits_dict, sng_reference_seq_dict


def _add_bsr_fallback_for_domains_without_pattern(
    database_path: str,
    context_free_domains_dict: Dict[str, Set[str]],
    domain_presence_intersection_pattern: Dict[str, Set[str]],
    singleton_reference_seqs_dict: Dict[str, Set[str]],
    domain_score_limits: Dict[str, Dict[str, float]],
    bsr_cutoff: float,
) -> Tuple[
    Dict[str, Dict[str, float]],
    Dict[str, Set[str]],
]:
    """
    Add BSR-based fallback reference sequences for seed domains for which
    no informative co-occurrence pattern could be identified.

    A pattern is considered non-informative if:
        1. no pattern exists for the seed domain,
        2. the pattern is empty, or
        3. the pattern contains only the seed domain itself.

    For these domains, all non-QUERY hits with blast_score_ratio >=
    bsr_cutoff are added to the singleton reference sequence set.

    Score limits are calculated from the qualifying hits and merged with
    existing limits if necessary.

    Args:
        database_path:
            Path to SQLite database.

        context_free_domains_dict:
            Seed domain -> genomes used for singleton/pattern detection.

        domain_presence_intersection_pattern:
            Seed domain -> conserved co-occurring domains.

        singleton_reference_seqs_dict:
            Existing singleton training sets:
            {domain: set(proteinIDs)}

        domain_score_limits:
            Existing score limits:
            {
                domain: {
                    "lower_limit": float,
                    "upper_limit": float
                }
            }

        bsr_cutoff:
            Minimum blast score ratio for fallback hits.

    Returns:
        Updated:
            domain_score_limits,
            singleton_reference_seqs_dict
    """

    # --------------------------------------------------------------
    # 1. Identify domains without an informative co-occurrence pattern
    # --------------------------------------------------------------

    fallback_domains = set()

    for seed_domain in context_free_domains_dict:

        pattern = domain_presence_intersection_pattern.get(seed_domain)

        if not pattern:
            fallback_domains.add(seed_domain)
            continue

        # Pattern consists only of the seed itself
        if set(pattern) == {seed_domain}:
            fallback_domains.add(seed_domain)

    if not fallback_domains:
        logger.info(
            "All singleton seed domains have informative co-occurrence patterns."
        )
        return domain_score_limits, singleton_reference_seqs_dict

    logger.info(
        f"Using BSR-only fallback for {len(fallback_domains)} domains "
        f"without informative co-occurrence pattern: "
        f"{', '.join(sorted(fallback_domains))}"
    )

    # --------------------------------------------------------------
    # 2. Fetch all qualifying hits for these domains
    # --------------------------------------------------------------

    with sqlite3.connect(database_path, timeout=120.0) as con:
        cur = con.cursor()

        cur.execute(
            """
            CREATE TEMP TABLE IF NOT EXISTS tmp_fallback_domains (
                domain TEXT PRIMARY KEY
            ) WITHOUT ROWID;
            """
        )

        cur.execute("DELETE FROM tmp_fallback_domains;")

        cur.executemany(
            """
            INSERT OR IGNORE INTO tmp_fallback_domains(domain)
            VALUES (?);
            """,
            ((domain,) for domain in fallback_domains),
        )

        # ----------------------------------------------------------
        # Protein IDs
        # ----------------------------------------------------------

        cur.execute(
            """
            SELECT DISTINCT
                d.domain,
                d.proteinID
            FROM Domains AS d
            JOIN Proteins AS p
              ON p.proteinID = d.proteinID
            JOIN tmp_fallback_domains AS f
              ON f.domain = d.domain
            WHERE d.blast_score_ratio >= ?
              AND p.genomeID != 'QUERY';
            """,
            (bsr_cutoff,),
        )

        added_counts = defaultdict(int)

        for domain, protein_id in cur:

            protein_set = singleton_reference_seqs_dict.setdefault(
                domain,
                set(),
            )

            old_size = len(protein_set)
            protein_set.add(protein_id)

            if len(protein_set) > old_size:
                added_counts[domain] += 1

        # ----------------------------------------------------------
        # 3. Calculate MIN/MAX score for exactly the same BSR-selected
        #    fallback hits
        # ----------------------------------------------------------

        cur.execute(
            """
            SELECT
                d.domain,
                MIN(d.score),
                MAX(d.score)
            FROM Domains AS d
            JOIN Proteins AS p
              ON p.proteinID = d.proteinID
            JOIN tmp_fallback_domains AS f
              ON f.domain = d.domain
            WHERE d.blast_score_ratio >= ?
              AND p.genomeID != 'QUERY'
              AND d.score IS NOT NULL
            GROUP BY d.domain;
            """,
            (bsr_cutoff,),
        )

        for domain, lower, upper in cur:

            if lower is None or upper is None:
                continue

            lower = float(lower)
            upper = float(upper)

            # Domain may already have sequences from another selection step.
            # In this case widen the existing limits.
            if domain in domain_score_limits:

                existing = domain_score_limits[domain]

                domain_score_limits[domain] = {
                    "lower_limit": min(
                        existing["lower_limit"],
                        lower,
                    ),
                    "upper_limit": max(
                        existing["upper_limit"],
                        upper,
                    ),
                }

            else:
                domain_score_limits[domain] = {
                    "lower_limit": lower,
                    "upper_limit": upper,
                }

    for domain in sorted(fallback_domains):
        logger.debug(
            f"BSR fallback for {domain}: "
            f"added {added_counts.get(domain, 0)} proteins "
            f"with BSR >= {bsr_cutoff}"
        )

    return domain_score_limits, singleton_reference_seqs_dict


#### Main routine of this module
def prepare_singleton_seed_proteins(
    options: Any,
) -> tuple[object, object] | tuple[dict[str, dict[str, float]], dict[str, set[str]]]:
    """
    Main routine: finds and predicts reference singletons for each protein/domain.

    Args:
        options (Any): Configuration/options.

    Returns:
        tuple:
            - domain_score_limits (dict): {domain: {'lower_limit', 'upper_limit'}}
            - singleton_reference_seqs_dict (dict): {domain: set(proteinIDs)}

    These two values define the domains that are taken to form the co-occurrence pattern:
    high_identitiy_cutoff is to detect the genomes where the singleton is likely correct
    cooccurence bsr cutoff is to define hits that are candidates for the pattern
    in the genomes where singleton is present

    Then hits that are cooccurring in at least 90 % of the genomes with this high homology singleton are taken
    """

    high_bsr_cutoff = getattr(options, "singleton_identity_cutoff", 70.0)
    cooccurence_bsr_cutoff = getattr(
        options, "singleton_identity_cutoff", 30.0
    )


    # 1) Add QUERY domains that were not selected by CSB grouping.
    # Genomes that are considered still need a seq with high identity
    context_free_domains_dict = _add_missing_query_domains_to_context_free_dict(
        database_path=options.database_directory,
        already_grouped_domains=options.grouped,
        blast_score_ratio_cutoff=options.low_hitscore_csb_cutoff,
    )

    logger.info(
        f"Found {len(context_free_domains_dict)} proteins without syntenic gene cluster and blast score ratio >= {high_bsr_cutoff}"
    )

    if not context_free_domains_dict:
        logger.warning(
            "No context-free high-identity hits found – stopping singleton selection."
        )
        return {}, {}

    # 2) For these genomes, determine coccurrence patterns of hits with bsr >= cooccurence_bsr_cutoff
    domain_presence_intersection_pattern = _fetch_conserved_co_occurrence_pattern(
        database_path=options.database_directory,
        domain_to_genomes=context_free_domains_dict,
        min_bsr_cutoff=cooccurence_bsr_cutoff,  # identity cutoff for co occurring pattern domains
        min_presence_fraction=cooccurence_bsr_cutoff,
    )

    # 3) Select seed hits for singletons from genomes with this cooccurence. Seed hit needs score >= high_bsr/2
    logger.info("Fetching co-occurence-based singleton candidates")
    domain_score_limits, singleton_reference_seqs_dict = (
        select_singleton_refs_by_domain_pattern(
            database_path=options.database_directory,
            seed_to_pattern_domains=domain_presence_intersection_pattern,
            min_bsr_cutoff=cooccurence_bsr_cutoff*0.5,
        )
    )

    # 4) Fallback for all that have no cooccurence pattern add everything above high bsr score
    domain_score_limits, singleton_reference_seqs_dict = (
        _add_bsr_fallback_for_domains_without_pattern(
            database_path=options.database_directory,
            context_free_domains_dict=context_free_domains_dict,
            domain_presence_intersection_pattern=domain_presence_intersection_pattern,
            singleton_reference_seqs_dict=singleton_reference_seqs_dict,
            domain_score_limits=domain_score_limits,
            bsr_cutoff=high_bsr_cutoff,
        )
    )

    myUtil.save_cache(
        options, "sng0_training_proteinIDs.pkl", singleton_reference_seqs_dict
    )
    myUtil.save_cache(
        options, "sng0_training_proteinIDs_limits.pkl", domain_score_limits
    )

    return domain_score_limits, singleton_reference_seqs_dict
