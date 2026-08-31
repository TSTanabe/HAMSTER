from __future__ import annotations

import math
import sqlite3
from array import array
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Set, Iterable, Any

import numpy as np

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold

from src.selection_defragmentation import mutual_information


# ======================================================================
# Data structures
# ======================================================================


@dataclass(slots=True)
class CalibrationSample:
    """
    One out-of-fold prediction used for calibration.
    """

    target: str
    genome_id: str
    fold: int

    truth: int
    evidence: float

    n_predictors: int

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class TargetCalibration:
    """
    Logistic calibration for one target.

    P(target present | context) =
        sigmoid(intercept + coefficient * context_evidence)
    """

    target: str

    intercept: float
    coefficient: float

    n_samples: int
    n_positive: int
    n_negative: int

    evidence_mean: float
    evidence_std: float
    evidence_min: float
    evidence_max: float

    def probability(
        self,
        evidence: float,
    ) -> float:
        """
        Convert context evidence into calibrated probability.
        """

        z = self.intercept + self.coefficient * float(evidence)

        # numerically stable sigmoid
        if z >= 0:
            return 1.0 / (1.0 + math.exp(-z))

        exp_z = math.exp(z)

        return exp_z / (1.0 + exp_z)

    def evidence_at_probability(
        self,
        probability: float,
    ) -> float:
        """
        Evidence score corresponding to a selected probability.

        Example:
            evidence_at_probability(0.5)
            evidence_at_probability(0.9)
            evidence_at_probability(0.95)
        """

        if not 0.0 < probability < 1.0:
            raise ValueError("probability must be between 0 and 1")

        if abs(self.coefficient) < 1e-12:
            raise ValueError(
                f"Calibration coefficient for {self.target!r} is approximately zero."
            )

        log_odds = math.log(probability / (1.0 - probability))

        return (log_odds - self.intercept) / self.coefficient

    def to_dict(self) -> Dict[str, Any]:
        result = asdict(self)

        if self.coefficient > 0:
            result["evidence_50"] = self.evidence_at_probability(0.50)

            result["evidence_90"] = self.evidence_at_probability(0.90)

            result["evidence_95"] = self.evidence_at_probability(0.95)

        else:
            result["evidence_50"] = None
            result["evidence_90"] = None
            result["evidence_95"] = None

        return result


@dataclass(slots=True)
class CalibrationResult:
    """
    Result of cross-validated predictor calibration.
    """

    n_folds: int
    n_genomes: int

    calibrations: Dict[
        str,
        TargetCalibration,
    ]

    samples: Dict[
        str,
        List[CalibrationSample],
    ]

    skipped_targets: Dict[
        str,
        str,
    ]

    def probability(
        self,
        target: str,
        evidence: float,
    ) -> float:
        calibration = self.calibrations.get(target)

        if calibration is None:
            raise KeyError(f"No calibration available for target {target!r}")

        return calibration.probability(evidence)

    def print_summary(
        self,
    ) -> None:
        print()
        print("Predictor calibration")
        print("=" * 80)

        print(f"Genomes: {self.n_genomes}")

        print(f"Folds: {self.n_folds}")

        print(f"Calibrated targets: {len(self.calibrations)}")

        print(f"Skipped targets: {len(self.skipped_targets)}")

        print()

        header = (
            f"{'Target':<20}"
            f"{'n':>8}"
            f"{'pos':>8}"
            f"{'neg':>8}"
            f"{'coef':>12}"
            f"{'E50':>12}"
            f"{'E90':>12}"
            f"{'E95':>12}"
        )

        print(header)
        print("-" * len(header))

        for target in sorted(self.calibrations):
            calibration = self.calibrations[target]

            if calibration.coefficient > 0:
                e50 = calibration.evidence_at_probability(0.50)

                e90 = calibration.evidence_at_probability(0.90)

                e95 = calibration.evidence_at_probability(0.95)

            else:
                e50 = float("nan")
                e90 = float("nan")
                e95 = float("nan")

            print(
                f"{target:<20}"
                f"{calibration.n_samples:>8}"
                f"{calibration.n_positive:>8}"
                f"{calibration.n_negative:>8}"
                f"{calibration.coefficient:>12.4f}"
                f"{e50:>12.3f}"
                f"{e90:>12.3f}"
                f"{e95:>12.3f}"
            )


# ======================================================================
# Database helpers
# ======================================================================


def _get_genome_ids(
    database_path: str,
) -> List[str]:
    """
    Return all genomes except QUERY.
    """

    with sqlite3.connect(
        database_path,
        timeout=120.0,
    ) as con:
        cur = con.cursor()

        cur.execute(
            """
            SELECT DISTINCT genomeID
            FROM Proteins
            WHERE genomeID IS NOT NULL
              AND genomeID != 'QUERY'
            """
        )

        genomes = sorted(genome_id for (genome_id,) in cur.fetchall())

    return genomes


def _fetch_genome_domain_values(
    database_path: str,
    genome_ids: Set[str],
    score_field: str = "score",
) -> tuple[
    Dict[str, Dict[str, float]],
    Dict[str, Dict[str, float]],
]:
    """
    Retrieve for each selected genome/domain:

        1. maximum target score
        2. maximum blast score ratio

    Returns
    -------
    genome_scores
        {
            genomeID: {
                domain: best_score
            }
        }

    genome_bsrs
        {
            genomeID: {
                domain: best_bsr
            }
        }
    """

    if score_field not in {
        "score",
        "blast_score_ratio",
    }:
        raise ValueError("score_field must be 'score' or 'blast_score_ratio'")

    genome_scores: Dict[
        str,
        Dict[str, float],
    ] = {genome_id: {} for genome_id in genome_ids}

    genome_bsrs: Dict[
        str,
        Dict[str, float],
    ] = {genome_id: {} for genome_id in genome_ids}

    if not genome_ids:
        return (
            genome_scores,
            genome_bsrs,
        )

    with sqlite3.connect(
        database_path,
        timeout=120.0,
    ) as con:
        cur = con.cursor()

        cur.execute("PRAGMA temp_store=MEMORY;")

        cur.execute(
            """
            CREATE TEMP TABLE tmp_calibration_genomes (
                genomeID TEXT PRIMARY KEY
            )
            """
        )

        cur.executemany(
            """
            INSERT INTO tmp_calibration_genomes(
                genomeID
            )
            VALUES (?)
            """,
            ((genome_id,) for genome_id in genome_ids),
        )

        query = f"""
            SELECT
                p.genomeID,
                d.domain,

                MAX(d.{score_field}) AS best_score,
                MAX(d.blast_score_ratio) AS best_bsr

            FROM Domains AS d

            JOIN Proteins AS p
              ON p.proteinID = d.proteinID

            JOIN tmp_calibration_genomes AS g
              ON g.genomeID = p.genomeID

            WHERE d.domain IS NOT NULL

            GROUP BY
                p.genomeID,
                d.domain
        """

        cur.execute(query)

        for (
            genome_id,
            domain,
            best_score,
            best_bsr,
        ) in cur:
            if best_score is not None:
                genome_scores[genome_id][domain] = float(best_score)

            if best_bsr is not None:
                genome_bsrs[genome_id][domain] = float(best_bsr)

    return (
        genome_scores,
        genome_bsrs,
    )


# ======================================================================
# Target truth
# ======================================================================


def _classify_target_truth(
    best_score: Optional[float],
    threshold: mutual_information.DomainThreshold,
) -> Optional[int]:
    """
    Convert target score into secure calibration truth.

    PRESENT:
        score >= average_cutoff

    ABSENT:
        no hit
        OR score < negative_cutoff

    UNKNOWN:
        between negative and average cutoff

    Returns
    -------
    1
        secure present

    0
        secure absent

    None
        ambiguous / unknown
    """

    # No target hit at all
    if best_score is None:
        return 0

    score = float(best_score)

    if score >= threshold.average_cutoff:
        return 1

    if score < threshold.negative_cutoff:
        return 0

    return None


# ======================================================================
# Logistic calibration
# ======================================================================


def _fit_target_calibration(
    target: str,
    evidence_values: array,
    truth_values: bytearray,
    min_samples: int = 10,
    min_positive: int = 3,
    min_negative: int = 3,
    calibration_C: float = 1.0,
) -> tuple[
    Optional[TargetCalibration],
    Optional[str],
]:
    """
    Fit:

        P(target present)
            =
        sigmoid(
            intercept
            + coefficient * evidence
        )

    Evidence is standardized during fitting for numerical stability,
    but the returned intercept/coefficient operate directly on the
    original evidence scale.
    """
    n_samples = len(evidence_values)

    if n_samples != len(truth_values):
        raise RuntimeError(
            f"Calibration evidence/truth size "
            f"mismatch for {target}: "
            f"{n_samples} vs "
            f"{len(truth_values)}"
        )

    if n_samples < min_samples:
        return (
            None,
            (f"only {n_samples} usable calibration samples"),
        )

    n_positive = sum(truth_values)

    n_negative = n_samples - n_positive

    if n_positive < min_positive:
        return (
            None,
            (f"only {n_positive} positive calibration samples"),
        )

    if n_negative < min_negative:
        return (
            None,
            (f"only {n_negative} negative calibration samples"),
        )

    evidence = np.frombuffer(
        evidence_values,
        dtype=np.float64,
    )

    truth = np.frombuffer(
        truth_values,
        dtype=np.uint8,
    )

    evidence_mean = float(evidence.mean())

    evidence_std = float(evidence.std())

    if evidence_std < 1e-12:
        return (
            None,
            "context evidence has no variation",
        )

    # --------------------------------------------------------------
    # Standardize only during logistic fitting
    # --------------------------------------------------------------

    x = (evidence - evidence_mean) / evidence_std

    classifier = LogisticRegression(
        C=calibration_C,
        solver="lbfgs",
        max_iter=2000,
        class_weight=None,
    )

    classifier.fit(
        x.reshape(-1, 1),
        truth,
    )

    standardized_coefficient = float(classifier.coef_[0][0])

    standardized_intercept = float(classifier.intercept_[0])

    # --------------------------------------------------------------
    # Convert back to raw evidence scale:
    #
    # a_std + b_std * ((E - mean) / std)
    #
    # =
    #
    # (a_std - b_std * mean/std)
    # +
    # (b_std/std) * E
    # --------------------------------------------------------------

    coefficient = standardized_coefficient / evidence_std

    intercept = (
        standardized_intercept - standardized_coefficient * evidence_mean / evidence_std
    )

    calibration = TargetCalibration(
        target=target,
        intercept=intercept,
        coefficient=coefficient,
        n_samples=n_samples,
        n_positive=n_positive,
        n_negative=n_negative,
        evidence_mean=evidence_mean,
        evidence_std=evidence_std,
        evidence_min=float(evidence.min()),
        evidence_max=float(evidence.max()),
    )

    return (
        calibration,
        None,
    )


# ======================================================================
# Main calibration routine
# ======================================================================


def calibrate_predictor_models(
    database_path: str,
    score_limit_dict: Dict[
        str,
        Dict[str, float],
    ],
    # --------------------------------------------------------------
    # Cross validation
    # --------------------------------------------------------------
    n_folds: int = 5,
    random_state: int = 42,
    # Optional subset of targets
    targets: Optional[Iterable[str]] = None,
    # --------------------------------------------------------------
    # Pairwise training
    # --------------------------------------------------------------
    negative_fraction: float = 0.5,
    score_field: str = "score",
    # additional arguments for
    # calculate_pairwise_domain_cooccurrence()
    #
    # useful if you later add e.g.
    # bsr_negative_fraction
    pairwise_extra_kwargs: Optional[Dict[str, Any]] = None,
    # --------------------------------------------------------------
    # Predictor selection
    # --------------------------------------------------------------
    max_fdr: Optional[float] = 0.05,
    min_log2_odds_ratio: Optional[float] = 2.0,
    min_both: int = 5,
    positive_only: bool = True,
    # --------------------------------------------------------------
    # Predictor redundancy
    # --------------------------------------------------------------
    redundancy_max_fdr: Optional[float] = 0.05,
    redundancy_min_phi: float = 0.3,
    alpha: float = 0.5,
    # --------------------------------------------------------------
    # Application
    # --------------------------------------------------------------
    missing_as_absent: bool = True,
    # --------------------------------------------------------------
    # Calibration requirements
    # --------------------------------------------------------------
    min_calibration_samples: int = 10,
    min_positive: int = 3,
    min_negative: int = 3,
    calibration_C: float = 1.0,
    # Store individual OOF samples for debugging.
    # False is strongly recommended for large production datasets.
    store_samples: bool = False,
) -> CalibrationResult:
    """
    Generate out-of-fold context-evidence predictions and calibrate
    them to target-presence probabilities.

    Workflow
    --------

    1. Split genomes into n_folds.

    2. For every fold:

        a. exclude held-out genomes

        b. calculate pairwise associations
           using only training genomes

        c. build predictor models
           using only training genomes

        d. retrieve BSR values for held-out genomes

        e. determine secure target truth:
               present / absent / unknown

        f. predict context evidence
           for held-out genomes

    3. Pool all out-of-fold evidence values.

    4. Fit target-specific logistic calibration:

           P(target present | evidence)
               =
           sigmoid(a + b * evidence)

    Notes
    -----
    The final predictor models used in production should still be
    trained separately on ALL genomes after calibration.

    The fold-specific models generated here exist only to produce
    unbiased out-of-fold evidence values for calibration.
    """

    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")

    if pairwise_extra_kwargs is None:
        pairwise_extra_kwargs = {}

    # Protect arguments controlled directly by this function
    reserved_pairwise_args = {
        "database_path",
        "score_limit_dict",
        "negative_fraction",
        "score_field",
        "exclude_genomes",
    }

    invalid_extra_args = reserved_pairwise_args.intersection(pairwise_extra_kwargs)

    if invalid_extra_args:
        raise ValueError(
            "pairwise_extra_kwargs contains "
            "reserved arguments: "
            f"{sorted(invalid_extra_args)}"
        )

    # --------------------------------------------------------------
    # Genome list
    # --------------------------------------------------------------

    genome_ids = _get_genome_ids(database_path)

    n_genomes = len(genome_ids)

    if n_genomes < n_folds:
        raise ValueError(f"Cannot create {n_folds} folds from only {n_genomes} genomes")

    if targets is not None:
        selected_targets = set(targets)

    else:
        selected_targets = None

    # --------------------------------------------------------------
    # Cross-validation folds
    # --------------------------------------------------------------

    fold_splitter = KFold(
        n_splits=n_folds,
        shuffle=True,
        random_state=random_state,
    )

    # For less RAM
    evidence_values: Dict[
        str,
        array,
    ] = {}

    # For less RAM
    truth_values: Dict[
        str,
        bytearray,
    ] = {}

    # For more RAM and debugging
    samples: Dict[
        str,
        List[CalibrationSample],
    ] = {}

    genome_array = np.asarray(
        genome_ids,
        dtype=object,
    )

    # ==============================================================
    # Fold loop
    # ==============================================================

    for (
        fold_index,
        (
            train_indices,
            test_indices,
        ),
    ) in enumerate(
        fold_splitter.split(genome_array),
        start=1,
    ):
        test_genomes = {str(genome_array[index]) for index in test_indices}

        # IMPORTANT:
        # calculate_pairwise_domain_cooccurrence()
        # only automatically excludes QUERY when
        # exclude_genomes is None.
        #
        # Therefore explicitly add QUERY here.
        excluded_genomes = set(test_genomes) | {"QUERY"}

        # ----------------------------------------------------------
        # 1. Train pairwise model on training folds
        # ----------------------------------------------------------

        fold_statistics = mutual_information.calculate_pairwise_domain_cooccurrence(
            database_path=database_path,
            score_limit_dict=score_limit_dict,
            negative_fraction=negative_fraction,
            score_field=score_field,
            exclude_genomes=excluded_genomes,
            **pairwise_extra_kwargs,
        )

        # ----------------------------------------------------------
        # 2. Build predictor models from training folds
        # ----------------------------------------------------------

        fold_models = fold_statistics.prepare_predictor_models(
            max_fdr=max_fdr,
            min_log2_odds_ratio=(min_log2_odds_ratio),
            min_both=min_both,
            positive_only=positive_only,
            redundancy_max_fdr=(redundancy_max_fdr),
            redundancy_min_phi=(redundancy_min_phi),
            alpha=alpha,
        )

        # ----------------------------------------------------------
        # 3. Read scores + BSRs for held-out genomes
        # ----------------------------------------------------------

        (
            genome_scores,
            genome_bsrs,
        ) = _fetch_genome_domain_values(
            database_path=database_path,
            genome_ids=test_genomes,
            score_field=score_field,
        )

        # ----------------------------------------------------------
        # 4. Apply fold models to held-out genomes
        # ----------------------------------------------------------

        for (
            target,
            model,
        ) in fold_models.items():
            if selected_targets is not None and target not in selected_targets:
                continue

            if not model.predictors:
                continue

            threshold = fold_statistics.thresholds.get(target)

            if threshold is None:
                continue

            # ----------------------------------------------------------
            # Compact target-specific buffers
            # ----------------------------------------------------------

            target_evidence = evidence_values.setdefault(
                target,
                array("d"),
            )

            target_truth = truth_values.setdefault(
                target,
                bytearray(),
            )

            for genome_id in test_genomes:
                best_target_score = genome_scores[genome_id].get(target)

                truth = _classify_target_truth(
                    best_score=best_target_score,
                    threshold=threshold,
                )

                if truth is None:
                    continue

                genome_bsr = genome_bsrs[genome_id]

                evidence = model.score_context_bsr(
                    genome_bsr=genome_bsr,
                    missing_as_absent=(missing_as_absent),
                )

                # ------------------------------------------------------
                # Compact data used for actual calibration
                # ------------------------------------------------------

                target_evidence.append(float(evidence))

                target_truth.append(int(truth))

                # ------------------------------------------------------
                # Optional full debugging information
                # ------------------------------------------------------

                if store_samples:
                    samples.setdefault(
                        target,
                        [],
                    ).append(
                        CalibrationSample(
                            target=target,
                            genome_id=genome_id,
                            fold=fold_index,
                            truth=truth,
                            evidence=float(evidence),
                            n_predictors=len(model.predictors),
                        )
                    )
    # ==============================================================
    # Fit target-specific calibration models
    # ==============================================================

    calibrations: Dict[
        str,
        TargetCalibration,
    ] = {}

    skipped_targets: Dict[
        str,
        str,
    ] = {}

    candidate_targets = (
        set(score_limit_dict) if selected_targets is None else selected_targets
    )

    for target in sorted(candidate_targets):
        target_evidence = evidence_values.get(target)

        target_truth = truth_values.get(target)

        if target_evidence is None:
            target_evidence = array("d")

        if target_truth is None:
            target_truth = bytearray()

        (
            calibration,
            reason,
        ) = _fit_target_calibration(
            target=target,
            evidence_values=(target_evidence),
            truth_values=(target_truth),
            min_samples=(min_calibration_samples),
            min_positive=min_positive,
            min_negative=min_negative,
            calibration_C=(calibration_C),
        )

        if calibration is None:
            skipped_targets[target] = reason if reason is not None else "unknown reason"

            continue

        calibrations[target] = calibration

    return CalibrationResult(
        n_folds=n_folds,
        n_genomes=n_genomes,
        calibrations=calibrations,
        samples=samples,
        skipped_targets=(skipped_targets),
    )


# ======================================================================
# Convenience function for final application
# ======================================================================


def predict_calibrated_probability(
    predictor_model: mutual_information.PreparedTargetModel,
    calibration: TargetCalibration,
    genome_bsr: Dict[str, float],
    missing_as_absent: bool = True,
) -> tuple[float, float]:
    """
    Apply final predictor model + calibration to a new genome.

    Returns
    -------
    evidence
    probability
    """

    if predictor_model.target != calibration.target:
        raise ValueError(
            "Predictor model and calibration "
            "refer to different targets: "
            f"{predictor_model.target!r} vs "
            f"{calibration.target!r}"
        )

    evidence = predictor_model.score_context_bsr(
        genome_bsr=genome_bsr,
        missing_as_absent=(missing_as_absent),
    )

    probability = calibration.probability(evidence)

    return (
        evidence,
        probability,
    )


def attach_calibrations_to_models(
    predictor_models: Dict[
        str,
        mutual_information.PreparedTargetModel,
    ],
    calibration_results: CalibrationResult,
) -> None:
    """
    Attach target-specific calibration parameters
    to final predictor models in place.
    """

    for (
        target,
        model,
    ) in predictor_models.items():
        calibration = calibration_results.calibrations.get(target)

        if calibration is None:
            continue

        model.calibration_intercept = calibration.intercept

        model.calibration_coefficient = calibration.coefficient
