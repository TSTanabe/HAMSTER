#!/usr/bin/python

import sqlite3
from typing import Any, Dict, Set

from src.core import myUtil
from src.core.logging import get_logger
from src.selection_seed import fetch_seed_proteins


logger = get_logger(__name__)


def extend_merged_grouped_by_csb_similarity(
    options: Any,
    grouped: Dict[str, Set[str]],
) -> Dict[str, Set[str]]:
    """
    Extend grouped protein sets by including proteins from
    highly similar CSB patterns.
    """

    protein_to_new_keywords_dict = myUtil.load_cache(
        options,
        "grp1_protein_to_key.pkl",
    )

    if not protein_to_new_keywords_dict:
        protein_to_new_keywords_dict = select_similar_csb_patterns_per_protein(
            options,
            grouped,
            options.jaccard,
        )

        myUtil.save_cache(
            options,
            "grp1_protein_to_key.pkl",
            protein_to_new_keywords_dict,
        )

    extended_grouped = myUtil.load_cache(
        options,
        "grp1_extended_grouped.pkl",
    )

    if not extended_grouped:
        extended_grouped = integrate_csb_variants_into_merged_grouped(
            options,
            grouped,
            protein_to_new_keywords_dict,
            options.sqlite_chunks,
        )

        myUtil.save_cache(
            options,
            "grp1_extended_grouped.pkl",
            extended_grouped,
        )

    return extended_grouped


def select_similar_csb_patterns_per_protein(
    options: Any,
    merged_grouped: Dict[str, Set[str]],
    jaccard_threshold: float = 0.7,
) -> Dict[str, Set[str]]:
    """
    For each domain, find CSB patterns with sufficient Jaccard
    similarity to CSBs containing its reference proteins.

    Returns
    -------
    Dict[str, Set[str]]
        domain -> similar CSB keywords
    """

    jaccard_included_patterns: Dict[str, Set[str]] = {}

    csb_dictionary = fetch_seed_proteins.parse_csb_file_to_dict(
        options.csb_output_file
    )

    for domain, protein_ids in merged_grouped.items():
        logger.debug(f"Processing domain {domain}")

        all_keywords = fetch_keywords_for_proteins(
            options.database_directory,
            protein_ids,
        )

        if not all_keywords:
            logger.warning(f"No keywords found for domain {domain}")
            continue

        domain_pattern_union = set().union(
            *(csb_dictionary[keyword] for keyword in all_keywords if keyword in csb_dictionary)
        )

        logger.debug(
            f"Domain {domain}: {len(domain_pattern_union)} unique CSB "
            f"elements encode reference sequences."
        )

        similar_csb_keywords: Set[str] = set()
        domain_pattern_union_len = len(domain_pattern_union)

        for csb_key, csb_pattern in csb_dictionary.items():
            csb_pattern_len = len(csb_pattern)

            denominator = max(
                domain_pattern_union_len,
                csb_pattern_len,
            )

            max_possible = (
                min(domain_pattern_union_len, csb_pattern_len) / denominator
                if denominator > 0
                else 0.0
            )

            if max_possible < jaccard_threshold:
                continue

            if domain_pattern_union_len < csb_pattern_len:
                intersection = sum(
                    1 for element in domain_pattern_union if element in csb_pattern
                )
            else:
                intersection = sum(
                    1 for element in csb_pattern if element in domain_pattern_union
                )

            if intersection == 0 and jaccard_threshold > 0.0:
                continue

            union = (
                domain_pattern_union_len
                + csb_pattern_len
                - intersection
            )

            similarity = (
                intersection / union
                if union > 0
                else 0.0
            )

            if similarity >= jaccard_threshold:
                similar_csb_keywords.add(csb_key)

        logger.debug(
            f"Domain {domain}: {len(similar_csb_keywords)} similar CSB patterns "
            f"(Jaccard >= {jaccard_threshold})"
        )

        jaccard_included_patterns[domain] = similar_csb_keywords

    return jaccard_included_patterns


def fetch_keywords_for_proteins(
    database_path: str,
    protein_ids: Set[str],
    chunk_size: int = 50000,
) -> Set[str]:
    """
    Fetch all CSB keywords associated with a set of protein IDs.
    """

    if not protein_ids:
        logger.warning(
            "No proteinIDs provided to fetch_keywords_for_proteins."
        )
        return set()

    with sqlite3.connect(
        database_path,
        timeout=120.0,
    ) as con:
        cur = con.cursor()

        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")
        cur.execute("PRAGMA mmap_size=2147483648;")
        cur.execute("PRAGMA automatic_index=ON;")

        cur.execute(
            """
            CREATE TEMP TABLE IF NOT EXISTS tmp_protein_ids (
                proteinID TEXT PRIMARY KEY
            );
            """
        )

        cur.execute(
            "DELETE FROM tmp_protein_ids;"
        )

        protein_id_list = list(protein_ids)

        for start in range(
            0,
            len(protein_id_list),
            chunk_size,
        ):
            batch = protein_id_list[
                start : start + chunk_size
            ]

            cur.executemany(
                """
                INSERT OR IGNORE INTO tmp_protein_ids(proteinID)
                VALUES (?);
                """,
                ((protein_id,) for protein_id in batch),
            )

        cur.execute(
            "PRAGMA query_only=TRUE;"
        )

        cur.execute(
            """
            SELECT DISTINCT k.keyword
            FROM tmp_protein_ids t

            JOIN Proteins p
              ON p.proteinID = t.proteinID

            LEFT JOIN Keywords k
              ON k.clusterID = p.clusterID

            WHERE k.keyword IS NOT NULL
            """
        )

        all_keywords = {
            keyword
            for (keyword,) in cur
        }

    logger.debug(
        f"Retrieved {len(all_keywords)} keywords for "
        f"{len(protein_ids)} proteins."
    )

    return all_keywords


def integrate_csb_variants_into_merged_grouped(
    options: Any,
    merged_grouped: Dict[str, Set[str]],
    domain_to_new_keywords_dict: Dict[str, Set[str]],
    chunk_size: int = 50000,
) -> Dict[str, Set[str]]:
    """
    Add proteins belonging to similar CSB patterns to the
    corresponding grouped domain sets.
    """

    logger.info(
        "Integration of added CSB proteins to grouped dataset"
    )

    if not domain_to_new_keywords_dict:
        return merged_grouped

    with sqlite3.connect(
        options.database_directory,
        timeout=120.0,
    ) as con:
        cur = con.cursor()

        cur.execute("PRAGMA temp_store=MEMORY;")
        cur.execute("PRAGMA cache_size=-262144;")
        cur.execute("PRAGMA mmap_size=2147483648;")
        cur.execute("PRAGMA automatic_index=ON;")

        cur.execute(
            """
            CREATE TEMP TABLE IF NOT EXISTS tmp_keywords (
                keyword TEXT PRIMARY KEY
            );
            """
        )

        total_proteins_added = 0

        for domain, new_keywords in domain_to_new_keywords_dict.items():
            if not new_keywords:
                logger.debug(
                    f"Domain {domain}: No new keywords to integrate."
                )
                continue

            merged_grouped.setdefault(
                domain,
                set(),
            )

            logger.debug(
                f"Domain {domain}: Integrating sequences from "
                f"{len(new_keywords)} new keywords."
            )

            before = len(
                merged_grouped[domain]
            )

            cur.execute(
                "DELETE FROM tmp_keywords;"
            )

            keyword_list = [
                keyword
                for keyword in new_keywords
                if keyword
            ]

            for start in range(
                0,
                len(keyword_list),
                chunk_size,
            ):
                batch = keyword_list[
                    start : start + chunk_size
                ]

                cur.executemany(
                    """
                    INSERT OR IGNORE INTO tmp_keywords(keyword)
                    VALUES (?);
                    """,
                    ((keyword,) for keyword in batch),
                )

            cur.execute(
                "PRAGMA query_only=TRUE;"
            )

            cur.execute(
                """
                SELECT DISTINCT p.proteinID

                FROM tmp_keywords tk

                JOIN Keywords k
                  ON k.keyword = tk.keyword

                JOIN Proteins p
                  ON p.clusterID = k.clusterID

                JOIN Domains d
                  ON d.proteinID = p.proteinID
                 AND d.domain = ?
                """,
                (domain,),
            )

            added = 0

            for (protein_id,) in cur:
                if protein_id in merged_grouped[domain]:
                    continue

                merged_grouped[domain].add(
                    protein_id
                )
                added += 1

            cur.execute(
                "PRAGMA query_only=FALSE;"
            )

            after = len(
                merged_grouped[domain]
            )

            logger.info(
                f"Domain {domain}: Added {added} new proteins with matching "
                f"synteny to reference. New total: {after} (was {before})"
            )

            total_proteins_added += added

    logger.info(
        f"Integration completed: {total_proteins_added} proteins "
        f"added across all domains."
    )

    return merged_grouped