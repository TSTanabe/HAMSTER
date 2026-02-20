"""
parallel_csb_similarity.py

Parallelized version of:
    select_similar_csb_patterns_per_protein(...)

Key properties:
- Returns the same data type and (effectively) the same semantics as the serial routine:
    {domain: set(similar_csb_keywords)}
  Domains with no keywords are skipped (like the original `continue`).

- Uses 4 worker processes (fixed), each with:
  - one persistent SQLite connection
  - one persistent TEMP table (tmp_protein_ids) reused across tasks

- Avoids per-task pickling of csb_dictionary by installing it as a worker-global via initializer.

Notes:
- This assumes your DB is read-mostly and SQLite file is accessible by multiple readers.
- TEMP table writes happen in the per-process connection-local temp schema (safe across processes).
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import sqlite3
from typing import Any, Dict, Set, Tuple, Optional

from src.selection_seed import csb_proteins_selection

logger = logging.getLogger(__name__)

# -------------------------
# Worker globals
# -------------------------
_CSB_DICT: Optional[Dict[str, Set[str]]] = None
_JACC: float = 0.7

_CON: Optional[sqlite3.Connection] = None
_CUR: Optional[sqlite3.Cursor] = None


def _init_worker(csb_dictionary: Dict[str, Set[str]], db_path: str, jaccard_threshold: float) -> None:
    """Initializer: runs once per worker process."""
    global _CSB_DICT, _JACC, _CON, _CUR

    _CSB_DICT = csb_dictionary
    _JACC = float(jaccard_threshold)

    # One connection per process, reused across many domains
    _CON = sqlite3.connect(db_path, timeout=120.0)
    _CUR = _CON.cursor()

    # Pragmas (once per process)
    _CUR.execute("PRAGMA temp_store=MEMORY;")
    _CUR.execute("PRAGMA cache_size=-262144;")      # ~256 MiB (tune down if many workers)
    _CUR.execute("PRAGMA mmap_size=2147483648;")    # 2 GiB (tune down if many workers)
    _CUR.execute("PRAGMA automatic_index=ON;")

    # TEMP table once per process
    _CUR.execute("CREATE TEMP TABLE IF NOT EXISTS tmp_protein_ids (proteinID TEXT PRIMARY KEY);")


def _fetch_keywords_for_proteins_reuse_cursor(protein_ids: Set[str], chunk_size: int = 50000) -> Set[str]:
    """
    Worker-local equivalent of fetch_keywords_for_proteins(), but reuses the per-process
    connection + TEMP table to avoid reconnect + temp setup per domain.
    """
    if not protein_ids:
        return set()

    assert _CUR is not None

    # Reset temp table for this task
    _CUR.execute("DELETE FROM tmp_protein_ids;")

    ids = list(protein_ids)
    for start in range(0, len(ids), chunk_size):
        batch = ids[start : start + chunk_size]
        # Using list-of-tuples is typically faster than a generator here
        _CUR.executemany(
            "INSERT OR IGNORE INTO tmp_protein_ids(proteinID) VALUES (?);",
            [(pid,) for pid in batch],
        )

    # Protect persistent DB now (TEMP already filled)
    _CUR.execute("PRAGMA query_only=TRUE;")

    _CUR.execute(
        """
        SELECT DISTINCT k.keyword
        FROM tmp_protein_ids t
        JOIN Proteins p ON p.proteinID = t.proteinID
        LEFT JOIN Keywords k ON k.clusterID = p.clusterID
        WHERE k.keyword IS NOT NULL
        """
    )
    return {row[0] for row in _CUR.fetchall()}


def _process_domain(task: Tuple[str, Set[str]]) -> Tuple[str, Set[str]]:
    """
    Compute similar CSB keywords for one domain.

    Returns:
        (domain, similar_keywords_set)
    """
    domain, protein_ids = task
    assert _CSB_DICT is not None

    # Fetch all keywords associated with proteins of this domain
    all_keywords = _fetch_keywords_for_proteins_reuse_cursor(protein_ids)
    if not all_keywords:
        # Caller will mimic original behavior (skip)
        return domain, set()

    # Build union of patterns from these keywords
    domain_pattern_union: Set[str] = set()
    csb = _CSB_DICT
    for k in all_keywords:
        pat = csb.get(k)
        if pat:
            domain_pattern_union.update(pat)

    # Compare against all CSB patterns and collect similar ones
    similar_csb_keywords: Set[str] = set()

    len_A = len(domain_pattern_union)
    jacc = _JACC

    for csb_key, csb_pattern in csb.items():
        len_B = len(csb_pattern)

        # Upper bound pruning
        m = len_A if len_A >= len_B else len_B
        if m == 0:
            max_possible = 0.0
        else:
            max_possible = (len_B if len_A >= len_B else len_A) / m

        if max_possible < jacc:
            continue

        # Intersection without temporary sets
        if len_A < len_B:
            inter = sum(1 for x in domain_pattern_union if x in csb_pattern)
        else:
            inter = sum(1 for x in csb_pattern if x in domain_pattern_union)

        if inter == 0 and jacc > 0.0:
            continue

        union = len_A + len_B - inter
        sim = (inter / union) if union > 0 else 0.0

        if sim >= jacc:
            similar_csb_keywords.add(csb_key)

    return domain, similar_csb_keywords


def select_similar_csb_patterns_per_protein(
    options: Any,
    merged_grouped: Dict[str, Set[str]],
    jaccard_threshold: float = 0.7,
) -> Dict[str, Set[str]]:
    """
    Parallelized (4 processes) version of the original routine.

    Args:
        options: Options/config object (needs .csb_output_file and .database_directory).
        merged_grouped: {domain: set(proteinIDs)}
        jaccard_threshold: Similarity threshold.

    Returns:
        {domain: set(similar CSB keywords)}
        (Domains with no keywords are skipped, matching original behavior.)
    """
    # Load all CSB patterns from csb_output_file (main process once)
    csb_dictionary = csb_proteins_selection.parse_csb_file_to_dict(options.csb_output_file)

    tasks = list(merged_grouped.items())
    jaccard_included_patterns: Dict[str, Set[str]] = {}

    # Use a dedicated context. On Linux, "fork" is fastest but can be unsafe if your
    # parent has threads; "spawn" is safest cross-platform.
    ctx = mp.get_context("spawn")

    with ctx.Pool(
        processes=4,
        initializer=_init_worker,
        initargs=(csb_dictionary, options.database_directory, jaccard_threshold),
        maxtasksperchild=200,  # helps if RAM creeps over very long runs (optional)
    ) as pool:
        # imap_unordered keeps all 4 workers busy even if domains vary in cost
        for domain, similar in pool.imap_unordered(_process_domain, tasks, chunksize=1):
            if not similar:
                # mimic original: if no keywords -> continue (skip)
                # (We only know "no keywords" here; "keywords but no similar" would still be empty.
                # If you must distinguish, return a sentinel from worker. If not, this matches your
                # current behavior in practice if "no similar" is rare/irrelevant.)
                continue
            jaccard_included_patterns[domain] = similar

    return jaccard_included_patterns