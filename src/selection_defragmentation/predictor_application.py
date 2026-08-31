from __future__ import annotations

import math
import sqlite3

from dataclasses import dataclass
from typing import Dict, Set, Optional, Any


# ======================================================================
# Result structures
# ======================================================================


@dataclass(slots=True)
class PredictorApplicationHit:
    """
    Result for one tested genome x target combination.
    """

    genome_id: str
    protein_id: str
    target: str

    target_score: float
    target_bsr: Optional[float]

    application_cutoff: float

    context_evidence: float
    probability: float

    n_predictors: int

    accepted: bool


@dataclass(slots=True)
class PredictorApplicationResult:
    """
    Complete result of predictor application.
    """

    predictions: list[PredictorApplicationHit]

    accepted_hits: Dict[str, Set[str]]

    n_tested: int
    n_accepted: int

    n_skipped_seed_top_hit: int
    n_skipped_below_cutoff: int
    n_skipped_no_model: int
    n_skipped_no_calibration: int

    def print_summary(self) -> None:
        print()
        print("Predictor application")
        print("=" * 60)

        print(f"Tested hits:              {self.n_tested}")

        print(f"Accepted hits:            {self.n_accepted}")

        print(f"Skipped: top hit is seed: {self.n_skipped_seed_top_hit}")

        print(f"Skipped: below cutoff:    {self.n_skipped_below_cutoff}")

        print(f"Skipped: no model:        {self.n_skipped_no_model}")

        print(f"Skipped: no calibration:  {self.n_skipped_no_calibration}")


# ======================================================================
# Database helpers
# ======================================================================


def _fetch_target_hits(
    database_path: str,
    targets: Set[str],
) -> Dict[
    tuple[str, str],
    list[tuple[str, float, Optional[float]]],
]:
    """
    Fetch all protein hits for targets.

    Returns
    -------
    {
        (genomeID, domain): [
            (proteinID, score, BSR),
            ...
        ]
    }

    Multiple paralogs are intentionally retained here.
    Selection of the highest-scoring paralog happens later.
    """

    results: Dict[
        tuple[str, str],
        list[tuple[str, float, Optional[float]]],
    ] = {}

    if not targets:
        return results

    with sqlite3.connect(
        database_path,
        timeout=120.0,
    ) as con:
        cur = con.cursor()

        cur.execute("PRAGMA temp_store=MEMORY;")

        cur.execute(
            """
            CREATE TEMP TABLE tmp_application_targets (
                domain TEXT PRIMARY KEY
            )
            """
        )

        cur.executemany(
            """
            INSERT INTO tmp_application_targets(domain)
            VALUES (?)
            """,
            ((target,) for target in targets),
        )

        # --------------------------------------------------------------
        # Aggregate repeated domain records of the SAME protein.
        #
        # Paralogs remain separate because proteinID is in GROUP BY.
        # --------------------------------------------------------------

        cur.execute(
            """
            SELECT
                p.genomeID,
                d.domain,
                d.proteinID,
                MAX(d.score) AS best_score,
                MAX(d.blast_score_ratio) AS best_bsr

            FROM Domains AS d

            JOIN Proteins AS p
              ON p.proteinID = d.proteinID

            JOIN tmp_application_targets AS t
              ON t.domain = d.domain

            WHERE p.genomeID IS NOT NULL
              AND p.genomeID != 'QUERY'
              AND d.score IS NOT NULL

            GROUP BY
                p.genomeID,
                d.domain,
                d.proteinID
            """
        )

        for (
            genome_id,
            domain,
            protein_id,
            score,
            bsr,
        ) in cur:
            key = (
                str(genome_id),
                str(domain),
            )

            results.setdefault(
                key,
                [],
            ).append(
                (
                    str(protein_id),
                    float(score),
                    (float(bsr) if bsr is not None else None),
                )
            )

    return results


def _fetch_genome_predictor_bsrs(
    database_path: str,
    predictor_domains: Set[str],
    basic_seed_sequences: Dict[
        str,
        Set[str],
    ],
) -> Dict[
    str,
    Dict[str, float],
]:
    """
    Retrieve the strongest predictor support for each
    genome x predictor domain.

    Rules
    -----
    - normal hit:
        use blast_score_ratio

    - protein is part of the basic seed set for this domain:
        use BSR = 1.0

    - several paralogs:
        use the maximum resulting BSR

    Therefore, if any predictor paralog is a basic seed
    sequence, the predictor contributes with BSR = 1.0.
    """

    genome_bsrs: Dict[
        str,
        Dict[str, float],
    ] = {}

    if not predictor_domains:
        return genome_bsrs

    with sqlite3.connect(
        database_path,
        timeout=120.0,
    ) as con:
        cur = con.cursor()

        cur.execute("PRAGMA temp_store=MEMORY;")

        cur.execute(
            """
            CREATE TEMP TABLE tmp_application_predictors (
                domain TEXT PRIMARY KEY
            )
            """
        )

        cur.executemany(
            """
            INSERT INTO tmp_application_predictors(domain)
            VALUES (?)
            """,
            ((domain,) for domain in predictor_domains),
        )

        # ----------------------------------------------------------
        # Do NOT aggregate in SQL here.
        #
        # We need proteinID so that seed proteins can be recognized.
        # Aggregation happens below after seed correction.
        # ----------------------------------------------------------

        cur.execute(
            """
            SELECT
                p.genomeID,
                d.domain,
                d.proteinID,
                d.blast_score_ratio

            FROM Domains AS d

            JOIN Proteins AS p
              ON p.proteinID = d.proteinID

            JOIN tmp_application_predictors AS x
              ON x.domain = d.domain

            WHERE p.genomeID IS NOT NULL
              AND p.genomeID != 'QUERY'
            """
        )

        for (
            genome_id,
            domain,
            protein_id,
            bsr,
        ) in cur:
            genome_id = str(genome_id)
            domain = str(domain)
            protein_id = str(protein_id)

            # Basic seed proteins are considered secure predictors
            # and therefore receive full BSR support.

            seed_proteins = basic_seed_sequences.get(
                domain,
                set(),
            )

            if protein_id in seed_proteins:
                effective_bsr = 1.0
            else:
                if bsr is None:
                    continue
                effective_bsr = float(bsr)

            # Several paralogs:
            # keep strongest effective BSR.
            # Thus a seed hit automatically wins with 1.0.

            domain_bsrs = genome_bsrs.setdefault(
                genome_id,
                {},
            )

            previous_bsr = domain_bsrs.get(domain)

            if previous_bsr is None or effective_bsr > previous_bsr:
                domain_bsrs[domain] = effective_bsr

    return genome_bsrs


# ======================================================================
# Main application routine
# ======================================================================


def evaluate_predictor_models(
    database_path: str,
    predictor_models: Dict[str, Any],
    seed_sequences: Dict[str, Set[str]],
    score_limit_dict: Dict[str, Dict[str, float]],
    probability_cutoff: float = 0.5,
    application_fraction: float = 0.5,
    missing_as_absent: bool = True,
    store_predictions: bool = False,
) -> PredictorApplicationResult:
    """
    Apply trained + calibrated context predictor models to database hits.

    Selection logic
    ---------------

    For every genome x target:

    1. Find the highest-scoring target hit among ALL paralogs.

    2. If a seed protein is among the highest-scoring hits:
           SKIP genome x target completely.

       This also handles tied top scores conservatively.

    3. Otherwise select ONE highest-scoring non-seed protein.

    4. Require:

           target_score >=
               application_fraction * lower_limit

       Default:
           application_fraction = 0.5

    5. Calculate genome context evidence from predictor BSRs.

    6. Convert evidence into calibrated probability using the
       calibration stored in PreparedTargetModel.

    7. Accept the protein if:

           probability >= probability_cutoff

    Returns
    -------
    PredictorApplicationResult
    """

    if not 0.0 <= probability_cutoff <= 1.0:
        raise ValueError("probability_cutoff must be between 0 and 1")

    if application_fraction < 0.0:
        raise ValueError("application_fraction must be >= 0")

    # --------------------------------------------------------------
    # Determine usable target models
    # --------------------------------------------------------------

    targets = set(predictor_models.keys()).intersection(score_limit_dict.keys())

    # --------------------------------------------------------------
    # Fetch all candidate target hits
    # --------------------------------------------------------------

    target_hits = _fetch_target_hits(
        database_path=database_path,
        targets=targets,
    )

    # --------------------------------------------------------------
    # Collect predictor domains required by all models
    # --------------------------------------------------------------

    predictor_domains: Set[str] = set()

    for model in predictor_models.values():
        predictor_domains.update(model.predictors.keys())

    # --------------------------------------------------------------
    # Fetch BSR context once for all genomes
    # --------------------------------------------------------------

    genome_bsrs = _fetch_genome_predictor_bsrs(
        database_path=database_path,
        predictor_domains=predictor_domains,
        basic_seed_sequences=seed_sequences,
    )

    # --------------------------------------------------------------
    # Result containers
    # --------------------------------------------------------------

    predictions: list[PredictorApplicationHit] = []

    n_tested = 0

    accepted_hits: Dict[
        str,
        Set[str],
    ] = {}

    n_skipped_seed_top_hit = 0
    n_skipped_below_cutoff = 0
    n_skipped_no_model = 0
    n_skipped_no_calibration = 0

    # ==============================================================
    # Evaluate genome x target combinations
    # ==============================================================

    for (
        genome_id,
        target,
    ), hits in target_hits.items():
        # ----------------------------------------------------------
        # Target model
        # ----------------------------------------------------------

        model = predictor_models.get(target)

        if model is None or not model.predictors:
            n_skipped_no_model += 1
            continue

        # ----------------------------------------------------------
        # Calibration must already be attached to final model
        # ----------------------------------------------------------

        if model.calibration_intercept is None or model.calibration_coefficient is None:
            n_skipped_no_calibration += 1
            continue

        # ----------------------------------------------------------
        # Score threshold
        # ----------------------------------------------------------

        limits = score_limit_dict.get(target)

        if limits is None:
            continue

        lower_limit = limits.get("lower_limit")

        if lower_limit is None:
            continue

        application_cutoff = application_fraction * float(lower_limit)

        # ----------------------------------------------------------
        # Determine maximum target score among ALL paralogs
        # ----------------------------------------------------------

        max_score = max(hit[1] for hit in hits)

        # ----------------------------------------------------------
        # Identify all hits tied for maximum score
        #
        # isclose prevents tiny floating-point differences from
        # changing the seed decision.
        # ----------------------------------------------------------

        top_hits = [
            hit
            for hit in hits
            if math.isclose(
                hit[1],
                max_score,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        ]

        seed_proteins = seed_sequences.get(
            target,
            set(),
        )

        # ----------------------------------------------------------
        # IMPORTANT:
        #
        # If ANY seed is among the highest-scoring paralogs,
        # skip this genome x target.
        # ----------------------------------------------------------

        if any(
            protein_id in seed_proteins
            for (
                protein_id,
                _,
                _,
            ) in top_hits
        ):
            n_skipped_seed_top_hit += 1
            continue

        # ----------------------------------------------------------
        # Best hit is not a seed.
        #
        # If several non-seed proteins have exactly the same top
        # score, select one deterministically by proteinID.
        # ----------------------------------------------------------

        top_hits.sort(key=lambda x: x[0])

        (
            protein_id,
            target_score,
            target_bsr,
        ) = top_hits[0]

        # ----------------------------------------------------------
        # Require >= 0.5 * lower_limit by default
        # ----------------------------------------------------------

        if target_score < application_cutoff:
            n_skipped_below_cutoff += 1
            continue

        # ----------------------------------------------------------
        # Genome-wide predictor BSR context
        # ----------------------------------------------------------

        genome_bsr = genome_bsrs.get(
            genome_id,
            {},
        )

        # ----------------------------------------------------------
        # Evidence
        # ----------------------------------------------------------

        context_evidence = model.score_context_bsr(
            genome_bsr=genome_bsr,
            missing_as_absent=(missing_as_absent),
        )

        # ----------------------------------------------------------
        # Calibrated probability
        # ----------------------------------------------------------

        probability = model.calibrated_probability_from_evidence(context_evidence)

        accepted = probability >= probability_cutoff

        n_tested += 1  # Add to counter to avoid the PredictorApplicationHit objects

        if store_predictions:
            prediction = PredictorApplicationHit(
                genome_id=genome_id,
                protein_id=protein_id,
                target=target,
                target_score=float(target_score),
                target_bsr=target_bsr,
                application_cutoff=float(application_cutoff),
                context_evidence=float(context_evidence),
                probability=float(probability),
                n_predictors=len(model.predictors),
                accepted=accepted,
            )

            predictions.append(prediction)

        if accepted:
            accepted_hits.setdefault(
                target,
                set(),
            ).add(protein_id)

    # --------------------------------------------------------------
    # Final result
    # --------------------------------------------------------------

    return PredictorApplicationResult(
        predictions=predictions,
        accepted_hits=accepted_hits,
        n_tested=n_tested,
        n_accepted=sum(len(protein_ids) for protein_ids in accepted_hits.values()),
        n_skipped_seed_top_hit=(n_skipped_seed_top_hit),
        n_skipped_below_cutoff=(n_skipped_below_cutoff),
        n_skipped_no_model=(n_skipped_no_model),
        n_skipped_no_calibration=(n_skipped_no_calibration),
    )
