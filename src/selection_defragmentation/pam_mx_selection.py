#!/usr/bin/python
import math
import sqlite3
import random
from typing import Dict, Any, List, Set, Tuple, Optional, Iterator, DefaultDict


from collections import defaultdict
from multiprocessing import Pool
from src.core.logging import get_logger

logger = get_logger(__name__)
_con = None
_cur = None

# Helper routines for fetching proteins from genomes and vice versa from database


def chunked(iterable: List[Any], size: int) -> Iterator[List[Any]]:
    """
    Yield successive chunk-sized lists from iterable.

    Args:
        iterable (list): Input list.
        size (int): Chunk size.

    Returns:
        Iterator[list]: Yields list of chunk size.

    Example:
        list(chunked([1,2,3,4,5], 2)) → [[1,2],[3,4],[5]]
    """
    seq = list(iterable)
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


def _get_readonly_cursor(database_path: str) -> sqlite3.Cursor:
    """
    Opens a read-only SQLite connection and caches cursor per process.

    Args:
        database_path (str): Path to SQLite DB.

    Returns:
        sqlite3.Cursor: Cursor to DB in read-only mode.
    """
    global _con, _cur
    if _con is None:
        _con = sqlite3.connect(
            f"file:{database_path}?mode=ro", uri=True, check_same_thread=False
        )
        _con.execute("PRAGMA query_only = TRUE;")
        _cur = _con.cursor()
    return _cur


def _fetch_protein_to_genome_chunk(args: Tuple[str, List[str]]) -> Dict[str, str]:
    """
    Maps a list of proteinIDs to genomeIDs for a DB chunk.

    Args:
        args: (database_path, chunk: List[str])

    Returns:
        dict: {proteinID: genomeID}
    """
    database_path, chunk = args
    if not chunk:
        return {}
    cur = _get_readonly_cursor(database_path)
    placeholders = ",".join("?" for _ in chunk)
    sql = (
        f"SELECT proteinID, genomeID FROM Proteins WHERE proteinID IN ({placeholders})"
    )
    cur.execute(sql, chunk)
    rows = cur.fetchall()
    return {proteinID: genomeID for proteinID, genomeID in rows}


def _fetch_domains_for_genome_chunk(
    args: Tuple[str, List[str]],
) -> Dict[str, Dict[str, List[str]]]:
    """
    Fetch all domain assignments for a chunk of genomeIDs.

    Parameters
    ----------
    args : tuple
        (database_path, genome_chunk)

    Returns
    -------
    dict
        Nested dict of the form:
        {
            genomeID: {
                domain_label: [proteinID1, proteinID2, ...]
            }
        }

    Notes
    -----
    - This reads *all* Domains entries for the selected genomes.
    - Multiple proteins of the same domain within one genome are retained.
    - The return format is intentionally compatible with the existing PAM format.
    """
    database_path, genome_chunk = args

    if not genome_chunk:
        return {}

    # Read-only connection for this worker/process
    con = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True)
    con.execute("PRAGMA query_only = TRUE;")
    cur = con.cursor()

    placeholders = ",".join("?" for _ in genome_chunk)

    query = f"""
        SELECT
            Proteins.genomeID,
            Domains.domain,
            Proteins.proteinID
        FROM Proteins
        JOIN Domains
          ON Proteins.proteinID = Domains.proteinID
        WHERE Proteins.genomeID IN ({placeholders})
    """

    cur.execute(query, genome_chunk)

    out: DefaultDict[str, DefaultDict[str, List[str]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for genome_id, domain, protein_id in cur.fetchall():
        # Collect all proteins per genome/domain pair
        out[genome_id][domain].append(protein_id)

    con.close()

    # Convert nested defaultdicts to normal dicts before returning
    return {gid: dict(dommap) for gid, dommap in out.items()}

from typing import Dict, List, Set


def _filter_pam_to_grouped_proteins(
    pam: Dict[str, Dict[str, List[str]]],
    grouped: Dict[str, Set[str]],
) -> Dict[str, Dict[str, List[str]]]:
    """
    Filter a PAM-like structure so that only proteinIDs remain that are present
    in the grouped reference sets for the same domain.

    Parameters
    ----------
    pam : dict
        PAM-like nested structure:
            {genomeID: {domain_label: [proteinIDs]}}
    grouped : dict
        Reference mapping:
            {domain_label: set(proteinIDs)}

    Returns
    -------
    dict
        Filtered PAM with the same structure:
            {genomeID: {domain_label: [proteinIDs]}}

        Rules:
        - Only proteinIDs present in grouped[domain_label] are retained.
        - Domains with no remaining proteinIDs are removed.
        - Genomes with no remaining domains are removed.
    """
    filtered_pam: Dict[str, Dict[str, List[str]]] = {}

    for genome_id, domain_map in pam.items():
        filtered_domain_map: Dict[str, List[str]] = {}

        for domain, protein_ids in domain_map.items():
            allowed_ids = grouped.get(domain)
            if not allowed_ids:
                continue

            kept = [protein_id for protein_id in protein_ids if protein_id in allowed_ids]

            if kept:
                filtered_domain_map[domain] = kept

        if filtered_domain_map:
            filtered_pam[genome_id] = filtered_domain_map

    return filtered_pam
#
#
#


def create_presence_absence_matrix_for_domain_positive_sample(
    proteinID_sets: Dict[str, Set[str]],
    target_domain: str,
    database_directory: str,
    chunk_size: int = 900,
    cores: int = 8,
    max_genomes: int = 1000,
    random_seed: int = 42,
) -> Dict[str, Dict[str, List[str]]]:
    """
    Create a domain-specific presence/absence matrix from a random subset of
    genomes carrying curated true-positive proteins of one target domain.

    Parameters
    ----------
    proteinID_sets : dict
        Mapping of the form:
            {domain_label: set(proteinIDs)}

        Important:
        proteinID_sets[target_domain] is assumed to already contain curated
        true-positive proteinIDs for the target domain.

    target_domain : str
        Domain to process.

    database_directory : str
        Path to the SQLite database.

    chunk_size : int, optional
        Maximum number of IDs per chunked DB query.

    cores : int, optional
        Number of worker processes.

    max_genomes : int, optional
        Maximum number of positive genomes to keep after random subsampling.

    random_seed : int, optional
        Seed for reproducible genome subsampling.

    Returns
    -------
    dict
        PAM-compatible nested mapping:
        {
            genomeID: {
                domain_label: [proteinID1, proteinID2, ...]
            }
        }

    Workflow
    --------
    1. Read the curated true-positive proteinIDs for target_domain.
    2. Map these proteinIDs to genomeIDs.
    3. Deduplicate to obtain the positive genome set.
    4. Randomly subsample genomes up to max_genomes.
    5. Fetch all present domains for these genomes from the DB.
    6. Return the resulting domain-specific PAM.

    Notes
    -----
    - No global PAM is built.
    - No BSR scores are used here.
    - This function is intended as a lightweight replacement for a global
      create_presence_absence_matrix() call when training should be based only
      on a positive genome subsample for one domain.
    """
    # Step 1: Get curated true-positive proteins for the requested domain
    target_proteins = proteinID_sets.get(target_domain, set())

    if not target_proteins:
        # No curated positives available for this domain
        logger.warning(f"No target proteins present for {target_domain}")
        return {}

    # Step 2: Map proteinID -> genomeID using the existing helper logic
    chunk_args = [
        (database_directory, chunk)
        for chunk in chunked(list(target_proteins), chunk_size)
    ]

    protein_to_genome: Dict[str, str] = {}

    with Pool(processes=cores) as pool:
        results = pool.map(_fetch_protein_to_genome_chunk, chunk_args)

    for partial in results:
        protein_to_genome.update(partial)

    # Step 3: Convert curated positive proteins into a deduplicated genome set
    positive_genomes = sorted(
        {
            protein_to_genome[pid]
            for pid in target_proteins
            if pid in protein_to_genome and protein_to_genome[pid]
        }
    )

    if not positive_genomes:
        # None of the proteins could be mapped back to genomes
        logger.warning("\t No positive genome")
        return {}

    # Step 4: Randomly subsample genomes up to the requested limit
    rng = random.Random(random_seed)

    if len(positive_genomes) > max_genomes:
        selected_genomes = rng.sample(positive_genomes, max_genomes)
    else:
        selected_genomes = positive_genomes

    # Step 5: Fetch all domains present in the selected genomes
    genome_chunk_args = [
        (database_directory, selected_genomes[i : i + chunk_size])
        for i in range(0, len(selected_genomes), chunk_size)
    ]

    genome_domain_matrix: DefaultDict[str, DefaultDict[str, List[str]]] = defaultdict(
        lambda: defaultdict(list)
    )

    with Pool(processes=cores) as pool:
        results = pool.map(_fetch_domains_for_genome_chunk, genome_chunk_args)

    for partial in results:
        for genome_id, domain_map in partial.items():
            for domain, protein_ids in domain_map.items():
                genome_domain_matrix[genome_id][domain].extend(protein_ids)

    # Step 6: Return PAM-compatible nested dict
    # only include gene that are in the grouped dataset, not all possible hits
    pam = {gid: dict(dommap) for gid, dommap in genome_domain_matrix.items()}
    pam = _filter_pam_to_grouped_proteins(pam, proteinID_sets)

    return pam

# Train support models
def train_positive_support_model_from_pam(
    pam: Dict[str, Dict[str, List[str]]],
    target_domain: str,
    alpha: float = 1.0,
    min_feature_genomes: int = 2,
) -> Dict[str, Any]:
    """
    Train a positive-only weighted support model for one target domain.

    The model learns which co-occurring domains are typical for genomes
    carrying target_domain. Each feature gets a normalized weight, so that
    all feature weights sum to 1.

    Feature strength:
        a_j = -log(1 - p_j)

    Normalized feature weight:
        w_j = a_j / sum(a_k)

    Parameters
    ----------
    pam : dict
        PAM-like nested structure:
            {genomeID: {domain_label: [proteinIDs]}}
    target_domain : str
        Domain whose genomic context should be modeled.
    alpha : float, optional
        Laplace smoothing parameter.
    min_feature_genomes : int, optional
        Minimum number of positive genomes in which a feature must occur
        to be retained.

    Returns
    -------
    dict
        Model dictionary containing:
            {
                "model_type": str,
                "target_domain": str,
                "alpha": float,
                "n_train_genomes": int,
                "features": List[str],
                "feature_probs": Dict[str, float],
                "feature_strengths": Dict[str, float],
                "feature_weights": Dict[str, float],
            }
    """
    positive_genomes = [
        genome_id
        for genome_id, dommap in pam.items()
        if target_domain in dommap and dommap[target_domain]
    ]

    if not positive_genomes:
        raise ValueError(
            f"No positive genomes found for target_domain={target_domain!r}."
        )

    feature_counts: Dict[str, int] = {}
    for genome_id in positive_genomes:
        dommap = pam[genome_id]
        for domain, proteins in dommap.items():
            if domain == target_domain:
                continue
            if not proteins:
                continue
            feature_counts[domain] = feature_counts.get(domain, 0) + 1

    features = sorted(
        domain
        for domain, count in feature_counts.items()
        if count >= min_feature_genomes
    )

    n_train = len(positive_genomes)

    feature_probs: Dict[str, float] = {}
    feature_strengths: Dict[str, float] = {}

    for feature in features:
        n_present = 0
        for genome_id in positive_genomes:
            if pam[genome_id].get(feature):
                n_present += 1

        p = (n_present + alpha) / (n_train + 2.0 * alpha)
        p = min(max(p, 1e-12), 1.0 - 1e-12)

        feature_probs[feature] = p
        feature_strengths[feature] = -math.log(1.0 - p)

    total_strength = sum(feature_strengths.values())

    if total_strength > 0.0:
        feature_weights = {
            feature: strength / total_strength
            for feature, strength in feature_strengths.items()
        }
    else:
        feature_weights = {feature: 0.0 for feature in features}

    model = {
        "model_type": "positive_support_weighted",
        "target_domain": target_domain,
        "alpha": alpha,
        "n_train_genomes": n_train,
        "features": features,
        "feature_probs": feature_probs,
        "feature_strengths": feature_strengths,
        "feature_weights": feature_weights,
    }

    return model


def train_support_models_for_each_domain(
    database_path: str,
    grouped: Dict[str, Set[str]],
    cores: int,
    chunk_size: int = 900,
    max_genomes: int = 1000,
    random_seed: int = 42,
    alpha: float = 1.0,
    min_feature_genomes: int = 2,
) -> Dict[str, Dict[str, Any]]:
    """
    For each domain, build a positive PAM sample and train a positive-only
    Bernoulli support model.

    Args:
        database_path (str): Path to SQLite DB.
        grouped (dict): {domain: set(proteinIDs)}. For each domain, the set is
                        assumed to contain curated true-positive proteinIDs.
        cores (int): Number of CPUs.
        chunk_size (int): DB query chunk size.
        max_genomes (int): Max number of positive genomes sampled per domain.
        random_seed (int): Seed for reproducible subsampling.
        alpha (float): Laplace smoothing parameter for support model training.
        min_feature_genomes (int): Minimum number of positive genomes in which
                                   a feature domain must occur to be retained.

    Returns:
        dict: {domain: support_model_dict}

        Returns:
        Dict[str, Dict[str, Any]]

        Mapping:
            {
                domain_label: support_model_dict
            }

        wobei support_model_dict die folgende Struktur hat:

            {
                "model_type": str,
                    # z. B. "positive_support_bernoulli"

                "target_domain": str,
                    # die Domain, für die das Modell trainiert wurde

                "alpha": float,
                    # Laplace-Smoothing-Parameter

                "n_train_genomes": int,
                    # Anzahl der Genome im positiven Trainingsset

                "features": List[str],
                    # Liste der verwendeten Feature-Domains
                    # (alle Domains außer target_domain, gefiltert nach min_feature_genomes)

                "feature_probs": Dict[str, float],
                    # Zentrales Modell:
                    # {feature_domain: p_j}
                    # mit p_j = P(feature_domain vorhanden | target_domain-positive Genome)

                "train_mean_raw_score": float,
                    # Mittelwert der Log-Likelihood-Scores der Trainingsgenome

                "train_std_raw_score": float,
                    # Standardabweichung der Trainingsscores (für z-Score Normalisierung)

                "train_min_raw_score": float,
                    # Minimaler Trainingsscore (optional für alternative Normalisierung)

                "train_max_raw_score": float,
                    # Maximaler Trainingsscore (optional für alternative Normalisierung)
            }

        Interpretation:
            - Jeder Eintrag im äußeren Dict entspricht einem Modell für eine Domain.
            - Das Modell beschreibt das typische Ko-Präsenzmuster der übrigen Domains
              in Genomen, die diese Ziel-Domain enthalten.
            - feature_probs enthält die eigentlichen gelernten Parameter.
            - Die train_* Werte dienen nur zur späteren Skalierung des Scores
              (z. B. für plausibility_score ∈ [0,1]).
    """
    total = len(grouped)

    logger.info(f"Training positive support models for {total} domains")

    support_models: Dict[str, Dict[str, Any]] = {}

    for i, domain in enumerate(sorted(grouped), start=1):
        logger.info(f"[{i}/{total}] Building positive PAM sample for domain {domain}")

        # 1) Build a small domain-specific positive PAM
        positive_pam = create_presence_absence_matrix_for_domain_positive_sample(
            proteinID_sets=grouped,
            target_domain=domain,
            database_directory=database_path,
            chunk_size=chunk_size,
            cores=cores,
            max_genomes=max_genomes,
            random_seed=random_seed,
        )

        if not positive_pam:
            logger.warning(
                f"No positive PAM could be built for domain {domain}, skipping."
            )
            continue

        # 2) Train the positive-only support model
        try:
            model = train_positive_support_model_from_pam(
                pam=positive_pam,
                target_domain=domain,
                alpha=alpha,
                min_feature_genomes=min_feature_genomes,
            )
        except ValueError as exc:
            logger.warning(f"Skipping domain {domain}: {exc}")
            continue

        support_models[domain] = model

        logger.info(
            f"Trained support model for {domain}: "
            f"{model['n_train_genomes']} genomes, "
            f"{len(model['features'])} retained features"
        )

    return support_models
