import math
import sqlite3
from collections import defaultdict
from typing import Dict, Set, List, Tuple, Any, DefaultDict

from src.core.logging import get_logger

logger = get_logger(__name__)
# Global variables for the databse read only mode
_con = None
_cur = None


def _fetch_all_genome_ids(database_path: str) -> List[str]:
    """
    Fetch all distinct genomeIDs from the Proteins table.

    Parameters
    ----------
    database_path : str
        Path to SQLite DB.

    Returns
    -------
    List[str]
        Sorted list of all genomeIDs present in the DB.
    """
    con = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True)
    con.execute("PRAGMA query_only = TRUE;")
    cur = con.cursor()

    cur.execute("SELECT DISTINCT genomeID FROM Proteins WHERE genomeID IS NOT NULL")
    genome_ids = [row[0] for row in cur.fetchall() if row[0]]

    con.close()
    genome_ids.sort()
    return genome_ids


def _fetch_best_domain_hits_for_genome_chunk(
    database_path: str,
    genome_chunk: List[str],
    score_field: str = "score",
) -> Dict[str, Dict[str, Tuple[str, float]]]:
    """
    Fetch all domain hits for a chunk of genomes and keep only the highest-scoring
    protein per genome/domain pair.

    Parameters
    ----------
    database_path : str
        Path to SQLite DB.
    genome_chunk : List[str]
        GenomeIDs to fetch.
    score_field : str, optional
        Which score field to use for selecting the best hit per domain.
        Allowed: "score" or "blast_score_ratio".

    Returns
    -------
    Dict[str, Dict[str, Tuple[str, float]]]
        Nested mapping:
            {
                genomeID: {
                    domain_label: (best_proteinID, best_score)
                }
            }

    Notes
    -----
    - The query joins Proteins and Domains on proteinID, consistent with the
      current PAM-building code.
    - If multiple proteins of the same domain exist in one genome, only the one
      with the highest selected score is retained.
    """
    if not genome_chunk:
        return {}

    if score_field not in {"score", "blast_score_ratio"}:
        raise ValueError("score_field must be 'score' or 'blast_score_ratio'")

    con = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True)
    con.execute("PRAGMA query_only = TRUE;")
    cur = con.cursor()

    placeholders = ",".join("?" for _ in genome_chunk)

    query = f"""
        SELECT
            Proteins.genomeID,
            Domains.domain,
            Proteins.proteinID,
            Domains.{score_field}
        FROM Proteins
        JOIN Domains
          ON Proteins.proteinID = Domains.proteinID
        WHERE Proteins.genomeID IN ({placeholders})
          AND Proteins.genomeID IS NOT NULL
          AND Domains.domain IS NOT NULL
    """

    cur.execute(query, genome_chunk)

    best_hits: DefaultDict[str, Dict[str, Tuple[str, float]]] = defaultdict(dict)

    for genome_id, domain, protein_id, score in cur.fetchall():
        if score is None:
            continue

        current = best_hits[genome_id].get(domain)
        if current is None or score > current[1]:
            best_hits[genome_id][domain] = (protein_id, float(score))

    con.close()
    return dict(best_hits)


def score_genome_against_support_model(
    model: Dict[str, Any],
    genome_pam: Dict[str, List[str]],
    require_target_domain: bool = True,
    return_feature_contributions: bool = False,
) -> Dict[str, Any]:
    """
    Score one genome against a weighted positive-only support model.

    Scoring logic
    -------------
    - Each feature has a normalized weight w_j, where all feature weights sum to 1.
    - If a feature is present, its weight contributes to positive evidence.
    - If a feature is missing, its weight contributes to negative evidence.
    - Therefore:
          positive_evidence + negative_evidence = 1
      (up to numerical rounding)

    Returned scores
    ---------------
    plausibility_balance:
        A directed evidence score in [-1, 1]

            +1  -> all weighted evidence is present
             0  -> positive and negative evidence are equally strong
            -1  -> all weighted evidence is missing

        Computed as:
            (positive_evidence - negative_evidence) / (positive_evidence + negative_evidence)

    plausibility_score:
        Same quantity rescaled to [0, 1]

            1.0 -> perfectly plausible
            0.5 -> ambiguous / balanced evidence
            0.0 -> strongly implausible

        Computed as:
            (plausibility_balance + 1.0) / 2.0

    Parameters
    ----------
    model : dict
        Model returned by train_positive_support_model_from_pam().
        Expected keys:
            - "target_domain"
            - "features"
            - "feature_probs"
            - "feature_weights"
    genome_pam : dict
        Domain mapping for a single genome:
            {domain_label: [proteinIDs]}
    require_target_domain : bool, optional
        If True, raise an error when the genome does not contain the model's
        target_domain.
    return_feature_contributions : bool, optional
        If True, also return per-feature contribution details.

    Returns
    -------
    dict
        {
            "target_domain": str,
            "has_target_domain": bool,
            "n_features_used": int,
            "plausibility_balance": float,   # [-1, 1]
            "plausibility_score": float,     # [0, 1]
            "positive_evidence": float,
            "negative_evidence": float,
            "present_feature_count": int,
            "missing_feature_count": int,
            "feature_contributions": Optional[Dict[str, Dict[str, Any]]],
        }
    """
    target_domain = model["target_domain"]
    features = model["features"]
    feature_probs = model["feature_probs"]
    feature_weights = model["feature_weights"]

    has_target = bool(genome_pam.get(target_domain))

    if require_target_domain and not has_target:
        raise ValueError(
            f"Genome does not contain target_domain={target_domain!r}, "
            f"but require_target_domain=True."
        )

    if not features:
        result = {
            "target_domain": target_domain,
            "has_target_domain": has_target,
            "n_features_used": 0,
            "plausibility_balance": 0.0,
            "plausibility_score": 0.5,
            "positive_evidence": 0.0,
            "negative_evidence": 0.0,
            "present_feature_count": 0,
            "missing_feature_count": 0,
        }
        if return_feature_contributions:
            result["feature_contributions"] = {}
        return result

    positive_evidence = 0.0
    negative_evidence = 0.0
    present_count = 0
    missing_count = 0
    feature_contributions = {}

    for feature in features:
        p = feature_probs[feature]
        weight = feature_weights[feature]
        present = bool(genome_pam.get(feature))

        if present:
            positive_evidence += weight
            present_count += 1
            if return_feature_contributions:
                feature_contributions[feature] = {
                    "present": 1,
                    "p_feature_given_positive": float(p),
                    "weight": float(weight),
                    "contribution_type": "positive",
                }
        else:
            negative_evidence += weight
            missing_count += 1
            if return_feature_contributions:
                feature_contributions[feature] = {
                    "present": 0,
                    "p_feature_given_positive": float(p),
                    "weight": float(weight),
                    "contribution_type": "negative",
                }

    total_evidence = positive_evidence + negative_evidence

    if total_evidence <= 0.0:
        plausibility_balance = 0.0
    else:
        plausibility_balance = (positive_evidence - negative_evidence) / total_evidence

    plausibility_score = (plausibility_balance + 1.0) / 2.0

    result = {
        "target_domain": target_domain,
        "has_target_domain": has_target,
        "n_features_used": len(features),
        "plausibility_balance": float(plausibility_balance),
        # "plausibility_score": float(plausibility_score),
        "positive_evidence": float(positive_evidence),
        "negative_evidence": float(negative_evidence),
        "present_feature_count": int(present_count),
        "missing_feature_count": int(missing_count),
    }

    if return_feature_contributions:
        result["feature_contributions"] = feature_contributions

    return result


def collect_plausible_domain_hits_from_support_models(
    support_models: Dict[str, Dict[str, Any]],
    database_path: str,
    chunk_size: int = 900,
    plausibility_cutoff: float = 0.8,
    score_field: str = "score",
) -> Dict[str, Set[str]]:
    """
    Iterate over all genomes in the database, build a genome-wise reduced PAM
    (one best protein per domain), score each present domain against its support
    model, and collect plausible proteinIDs.

    Parameters
    ----------
    support_models : Dict[str, Dict[str, Any]]
        Mapping:
            {domain_label: support_model_dict}

        Each support_model_dict is assumed to come from
        train_positive_support_model_from_pam(...).

    database_path : str
        Path to SQLite DB.

    chunk_size : int, optional
        Number of genomeIDs processed per DB chunk.

    plausibility_cutoff : float, optional
        Minimum plausibility_score required to retain a domain hit.

    score_field : str, optional
        Score column used to choose the best protein per genome/domain.
        Allowed: "score" or "blast_score_ratio".

    Returns
    -------
    Dict[str, Set[str]]
        Final mapping:
            {
                domain_label: {proteinID1, proteinID2, ...}
            }

        Only proteinIDs whose domain occurrence is considered plausible are returned.

    Workflow
    --------
    1. Fetch all genomeIDs from the DB.
    2. Process genomes in chunks.
    3. For each genome:
       - fetch all domain hits
       - keep only the highest-scoring protein for each domain
       - build a single-genome PAM-like dict: {domain: [proteinID]}
       - for each present domain that has a support model:
         score its plausibility in the context of the full genome
       - if plausibility_score >= cutoff:
         keep the proteinID for that domain
    4. Aggregate results as {domain: set(proteinIDs)}

    Notes
    -----
    - The current global PAM code stores genome-domain-protein relationships in
      the form {genomeID: {domain: [proteinIDs]}}. This routine builds the same
      kind of per-genome information, just reduced to one best protein per domain. :contentReference[oaicite:3]{index=3}
    - The score used here is a support/plausibility score, not a calibrated posterior
      probability.
    """
    logger.info(
        f"Scoring all genomes in DB against {len(support_models)} support models"
    )

    all_genomes = _fetch_all_genome_ids(database_path)
    total = len(all_genomes)

    if total == 0:
        logger.warning("No genomes found in database.")
        return {}

    plausible_hits: DefaultDict[str, Set[str]] = defaultdict(set)

    for start in range(0, total, chunk_size):
        end = min(start + chunk_size, total)
        genome_chunk = all_genomes[start:end]

        logger.info(
            f"Scoring genomes chunk {start + 1}-{end} / {total} "
            f"({(end / total) * 100:.1f}%)"
        )

        # genome -> domain -> (best_proteinID, best_score)
        best_hits_by_genome = _fetch_best_domain_hits_for_genome_chunk(
            database_path=database_path,
            genome_chunk=genome_chunk,
            score_field=score_field,
        )

        for genome_id, domain_map in best_hits_by_genome.items():
            if not domain_map:
                continue

            # Convert to the single-genome PAM format expected by the scoring routine:
            # {domain: [proteinID]}
            genome_pam = {
                domain: [protein_id]
                for domain, (protein_id, _score) in domain_map.items()
            }
            # print("Genome presence domain:")
            # print(genome_pam)
            # Score each *present* domain separately, if a support model exists
            for domain, (protein_id, _score) in domain_map.items():
                model = support_models.get(domain)
                if model is None:
                    continue

                try:
                    result = score_genome_against_support_model(
                        model=model,
                        genome_pam=genome_pam,
                        require_target_domain=True,
                        return_feature_contributions=False,
                    )
                except ValueError:
                    # Should normally not happen, because domain is present in genome_pam
                    continue

                if result["plausibility_balance"] >= plausibility_cutoff:
                    plausible_hits[domain].add(protein_id)
                # else:
                #    print("Non plausible")
                #    print(result)
    return dict(plausible_hits)
