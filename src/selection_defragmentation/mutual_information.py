from __future__ import annotations

import math
import random
import sqlite3

from dataclasses import dataclass, asdict
from itertools import combinations
from typing import Dict, Set, Tuple, Any, Optional

from scipy.stats import fisher_exact


# ----------------------------------------------------------------------
# Data structures
# ----------------------------------------------------------------------


@dataclass(slots=True)
class DomainThreshold:
    minimum_cutoff: float
    average_cutoff: float
    negative_cutoff: float


@dataclass(slots=True)
class BSRThreshold:
    lower_limit: float
    average: float
    upper_limit: float
    negative_cutoff: float


@dataclass(slots=True)
class DomainState:
    """
    Genome classification for one domain.

    absent genomes are not stored explicitly:
        absent = all_genomes - present - unknown
    """

    present: Set[int]
    unknown: Set[int]


@dataclass(slots=True)
class PairwiseAssociation:
    domain_a: str
    domain_b: str

    n_total_genomes: int
    n_tested: int
    n_excluded_unknown: int

    n11: int
    n10: int
    n01: int
    n00: int

    n_a_present: int
    n_b_present: int

    p_b_given_a: Optional[float]
    p_a_given_b: Optional[float]

    expected_both: Optional[float]
    fold_enrichment: Optional[float]

    jaccard: Optional[float]
    mutual_information: Optional[float]
    phi: Optional[float]

    odds_ratio: Optional[float]
    log2_odds_ratio: Optional[float]

    fisher_p: float
    fdr_q: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def format_summary(self) -> str:
        """Return a human-readable summary of this domain-pair association."""

        def fmt(value, digits=3):
            if value is None:
                return "NA"
            return f"{value:.{digits}f}"

        def fmt_p(value):
            if value is None:
                return "NA"
            if value == 0:
                return "<1e-300"
            if value < 0.001:
                return f"{value:.2e}"
            return f"{value:.4f}"

        return (
            f"{self.domain_a} ↔ {self.domain_b}\n"
            f"  genomes: tested={self.n_tested:,} / total={self.n_total_genomes:,} "
            f"(excluded unknown={self.n_excluded_unknown:,})\n"
            f"  contingency: both={self.n11:,}, "
            f"{self.domain_a}-only={self.n10:,}, "
            f"{self.domain_b}-only={self.n01:,}, "
            f"neither={self.n00:,}\n"
            f"  P({self.domain_b}|{self.domain_a}) = {fmt(self.p_b_given_a)}\n"
            f"  P({self.domain_a}|{self.domain_b}) = {fmt(self.p_a_given_b)}\n"
            f"  mutual information = {fmt(self.mutual_information, 4)} bits\n"
            f"  log2 odds ratio    = {fmt(self.log2_odds_ratio)}\n"
            f"  fold enrichment    = {fmt(self.fold_enrichment)}\n"
            f"  phi                = {fmt(self.phi)}\n"
            f"  Jaccard            = {fmt(self.jaccard)}\n"
            f"  Fisher p           = {fmt_p(self.fisher_p)}\n"
            f"  FDR q              = {fmt_p(self.fdr_q)}"
        )

    def print_summary(self) -> None:
        print(self.format_summary())


@dataclass(slots=True)
class DomainPredictor:
    """
    Directed interpretation of a pairwise association:

        predictor -> target
    """

    target: str
    predictor: str

    n_tested: int
    n_both: int

    p_target: float
    p_target_given_predictor: float
    p_target_without_predictor: float

    lift: float

    log2_odds_ratio: Optional[float]
    mutual_information: Optional[float]
    phi: Optional[float]

    fisher_p: float
    fdr_q: Optional[float]


@dataclass(slots=True)
class PreparedPredictor:
    """
    Predictor with precomputed likelihood evidence and
    redundancy correction.
    """

    target: str
    predictor: str

    n_tested: int
    n_both: int
    fdr_q: Optional[float]

    # probabilities used for likelihood calculation
    p_predictor_given_target: float
    p_predictor_given_not_target: float

    # raw evidence
    log_lr_present: float
    log_lr_absent: float

    # predictor-predictor redundancy
    redundancy_load: float
    weight: float

    # final evidence to use during genome scoring
    weighted_log_lr_present: float
    weighted_log_lr_absent: float

    # BSR -> presence transformation. Instead presence absence low bsr limit > score = absent
    bsr_lower_limit: float
    bsr_average: float
    bsr_upper_limit: float
    bsr_negative_cutoff: float

    def bsr_to_soft_presence(
        self,
        bsr: float,
    ) -> float:
        """
        Convert BSR into continuous presence support [0, 1].

        BSR <= negative cutoff
            -> 0.0

        BSR >= average BSR of trusted/reference hits
            -> 1.0

        Between both boundaries
            -> linear interpolation.
        """

        bsr = max(
            0.0,
            min(1.0, float(bsr)),
        )

        lower = self.bsr_negative_cutoff
        upper = self.bsr_average

        if bsr <= lower:
            return 0.0

        if bsr >= upper:
            return 1.0

        if upper <= lower:
            return 1.0

        return (bsr - lower) / (upper - lower)

    def evidence_from_bsr(
        self,
        bsr: float,
    ) -> float:
        soft_presence = self.bsr_to_soft_presence(bsr)

        return (
            soft_presence * self.weighted_log_lr_present
            + (1.0 - soft_presence) * self.weighted_log_lr_absent
        )


@dataclass(slots=True)
class PreparedTargetModel:
    target: str

    predictors: Dict[
        str,
        PreparedPredictor,
    ]

    calibration_intercept: Optional[float] = None
    calibration_coefficient: Optional[float] = None

    def score_context_bsr(
        self,
        genome_bsr: Dict[str, float],
        missing_as_absent: bool = True,
    ) -> float:
        evidence = 0.0

        for (
            predictor_name,
            predictor,
        ) in self.predictors.items():
            bsr = genome_bsr.get(predictor_name)

            if bsr is None:
                if missing_as_absent:
                    bsr = 0.0
                else:
                    continue

            evidence += predictor.evidence_from_bsr(bsr)

        return evidence

    def calibrated_probability_from_evidence(
        self,
        evidence: float,
    ) -> float:
        """
        Convert context evidence into calibrated
        target-presence probability.
        """

        if self.calibration_intercept is None or self.calibration_coefficient is None:
            raise ValueError(f"No calibration available for target {self.target!r}")

        z = self.calibration_intercept + self.calibration_coefficient * float(evidence)

        # numerically stable sigmoid
        if z >= 0:
            return 1.0 / (1.0 + math.exp(-z))

        exp_z = math.exp(z)

        return exp_z / (1.0 + exp_z)

    def probability_from_bsr(
        self,
        genome_bsr: Dict[str, float],
        missing_as_absent: bool = True,
    ) -> float:
        """
        Directly calculate calibrated target-presence
        probability from genome BSRs.
        """

        evidence = self.score_context_bsr(
            genome_bsr=genome_bsr,
            missing_as_absent=missing_as_absent,
        )

        return self.calibrated_probability_from_evidence(evidence)


# ======================================================================
# Main result object
# ======================================================================


@dataclass(slots=True)
class PairwiseCooccurrenceResult:
    n_genomes: int

    thresholds: Dict[
        str,
        DomainThreshold,
    ]

    pairs: Dict[
        Tuple[str, str],
        PairwiseAssociation,
    ]

    bsr_thresholds: Dict[
        str,
        BSRThreshold,
    ]

    # ------------------------------------------------------------------
    # Internal helper:
    # get pair independent of A/B storage orientation
    # ------------------------------------------------------------------

    def _get_pair(
        self,
        domain_a: str,
        domain_b: str,
    ) -> Optional[PairwiseAssociation]:
        pair = self.pairs.get((domain_a, domain_b))

        if pair is not None:
            return pair

        return self.pairs.get((domain_b, domain_a))

    # ------------------------------------------------------------------
    # Orient contingency table:
    #
    #              predictor+
    #              predictor-
    #
    # target+      n11       n_target_only
    # target-      n_pred_only   n00
    # ------------------------------------------------------------------

    def _get_oriented_counts(
        self,
        target: str,
        predictor: str,
    ):
        stats = self._get_pair(
            target,
            predictor,
        )

        if stats is None:
            return None

        # Stored orientation already:
        #
        # target = A
        # predictor = B
        if stats.domain_a == target and stats.domain_b == predictor:
            return (
                stats,
                stats.n11,
                stats.n10,  # target only
                stats.n01,  # predictor only
                stats.n00,
            )

        # Reverse orientation:
        #
        # target = B
        # predictor = A
        return (
            stats,
            stats.n11,
            stats.n01,  # target only
            stats.n10,  # predictor only
            stats.n00,
        )

    # ==================================================================
    # Obtain directed predictors for one target
    # ==================================================================

    def get_predictors(
        self,
        target: str,
        max_fdr: Optional[float] = 0.05,
        min_log2_odds_ratio: Optional[float] = 2.0,
        min_both: int = 3,
        positive_only: bool = True,
    ) -> list[DomainPredictor]:
        """
        Select useful predictors for one target domain.

        A predictor is retained if:

            1. sufficient shared observations exist (n11 >= min_both)

            2. if positive_only=True:
                   P(target | predictor) > P(target | no predictor)

            3. AND at least one of the following is true:
                   FDR q <= max_fdr
                   OR
                   log2 odds ratio >= min_log2_odds_ratio

        This prevents strong predictors from being discarded solely because
        statistical significance is limited by a small number of evaluable
        genomes.
        """

        predictors = []

        for stats in self.pairs.values():
            # ----------------------------------------------------------
            # Identify predictor for this target
            # ----------------------------------------------------------

            if stats.domain_a == target:
                predictor = stats.domain_b

            elif stats.domain_b == target:
                predictor = stats.domain_a

            else:
                continue

            oriented = self._get_oriented_counts(
                target,
                predictor,
            )

            if oriented is None:
                continue

            (
                stats,
                n11,
                n_target_only,
                n_predictor_only,
                n00,
            ) = oriented

            # ----------------------------------------------------------
            # Minimum support
            # ----------------------------------------------------------

            if n11 < min_both:
                continue

            if stats.n_tested <= 0:
                continue

            # ----------------------------------------------------------
            # Baseline P(target)
            # ----------------------------------------------------------

            p_target = (n11 + n_target_only) / stats.n_tested

            # ----------------------------------------------------------
            # P(target | predictor)
            # ----------------------------------------------------------

            predictor_present_n = n11 + n_predictor_only

            if predictor_present_n > 0:
                p_target_given_predictor = n11 / predictor_present_n
            else:
                continue

            # ----------------------------------------------------------
            # P(target | predictor absent)
            # ----------------------------------------------------------

            predictor_absent_n = n_target_only + n00

            if predictor_absent_n > 0:
                p_target_without_predictor = n_target_only / predictor_absent_n
            else:
                p_target_without_predictor = 0.0

            # ----------------------------------------------------------
            # Require positive predictive relationship
            # ----------------------------------------------------------

            if positive_only:
                if p_target_given_predictor <= p_target_without_predictor:
                    continue

                # additionally reject explicitly negative OR
                if stats.log2_odds_ratio is not None and stats.log2_odds_ratio <= 0:
                    continue

            # ----------------------------------------------------------
            # Statistical significance criterion
            # ----------------------------------------------------------

            passes_fdr = False

            if max_fdr is None:
                passes_fdr = True

            elif stats.fdr_q is not None:
                passes_fdr = stats.fdr_q <= max_fdr

            # ----------------------------------------------------------
            # Effect-strength criterion
            # ----------------------------------------------------------

            passes_effect_size = False

            if min_log2_odds_ratio is None:
                passes_effect_size = True

            elif stats.log2_odds_ratio is not None:
                passes_effect_size = stats.log2_odds_ratio >= min_log2_odds_ratio

            # ----------------------------------------------------------
            # Predictor must have either:
            #
            # good statistical support
            # OR
            # strong predictive effect
            # ----------------------------------------------------------

            if not (passes_fdr or passes_effect_size):
                continue

            # ----------------------------------------------------------
            # Lift
            # ----------------------------------------------------------

            lift = p_target_given_predictor / p_target if p_target > 0 else 0.0

            predictors.append(
                DomainPredictor(
                    target=target,
                    predictor=predictor,
                    n_tested=stats.n_tested,
                    n_both=n11,
                    p_target=p_target,
                    p_target_given_predictor=(p_target_given_predictor),
                    p_target_without_predictor=(p_target_without_predictor),
                    lift=lift,
                    log2_odds_ratio=(stats.log2_odds_ratio),
                    mutual_information=(stats.mutual_information),
                    phi=stats.phi,
                    fisher_p=stats.fisher_p,
                    fdr_q=stats.fdr_q,
                )
            )

        # --------------------------------------------------------------
        # Sort strongest predictors first
        # --------------------------------------------------------------

        predictors.sort(
            key=lambda x: (
                x.p_target_given_predictor - x.p_target_without_predictor,
                x.log2_odds_ratio if x.log2_odds_ratio is not None else float("-inf"),
                x.n_both,
            ),
            reverse=True,
        )

        return predictors

    # ==================================================================
    # Human-readable predictor printout
    # ==================================================================

    def print_predictors(
        self,
        target: str,
        max_fdr: Optional[float] = None,
        min_both: int = 0,
        positive_only: bool = False,
        max_predictors: Optional[int] = None,
    ) -> None:
        predictors = self.get_predictors(
            target=target,
            max_fdr=max_fdr,
            min_both=min_both,
            positive_only=positive_only,
        )

        if max_predictors is not None:
            predictors = predictors[:max_predictors]

        print(f"\nPredictors for {target}\n{'=' * (15 + len(target))}")

        print(
            f"{'Predictor':<20}"
            f"{'P(T|P)':>10}"
            f"{'P(T|!P)':>10}"
            f"{'Lift':>10}"
            f"{'log2OR':>10}"
            f"{'n11':>8}"
            f"{'FDR q':>12}"
        )

        print("-" * 80)

        for predictor in predictors:
            log_or = (
                f"{predictor.log2_odds_ratio:.3f}"
                if predictor.log2_odds_ratio is not None
                else "NA"
            )

            fdr = f"{predictor.fdr_q:.2e}" if predictor.fdr_q is not None else "NA"

            print(
                f"{predictor.predictor:<20}"
                f"{predictor.p_target_given_predictor:>10.3f}"
                f"{predictor.p_target_without_predictor:>10.3f}"
                f"{predictor.lift:>10.3f}"
                f"{log_or:>10}"
                f"{predictor.n_both:>8}"
                f"{fdr:>12}"
            )

    # ==================================================================
    # Prepare all target models
    # ==================================================================

    def prepare_predictor_models(
        self,
        max_fdr: Optional[float] = 0.05,
        min_log2_odds_ratio: Optional[float] = 2.0,
        min_both: int = 3,
        positive_only: bool = True,
        redundancy_max_fdr: Optional[float] = 0.05,
        redundancy_min_phi: float = 0.3,
        alpha: float = 0.5,
    ) -> Dict[str, PreparedTargetModel]:
        """
        Build precomputed predictor models for every target domain.

        Steps
        -----
        1. Select statistically useful predictors for target.
        2. Calculate present/absent likelihood ratios.
        3. Determine redundancy between predictors.
        4. Downweight redundant predictors.
        5. Store final weighted evidence.
        """

        models: Dict[
            str,
            PreparedTargetModel,
        ] = {}

        # All domains for which score thresholds exist
        targets = sorted(self.thresholds.keys())

        for target in targets:
            selected = self.get_predictors(
                target=target,
                max_fdr=max_fdr,
                min_log2_odds_ratio=min_log2_odds_ratio,
                min_both=min_both,
                positive_only=positive_only,
            )

            if not selected:
                models[target] = PreparedTargetModel(
                    target=target,
                    predictors={},
                )

                continue

            predictor_names = {p.predictor for p in selected}

            # ----------------------------------------------------------
            # Calculate redundancy load for each predictor
            # ----------------------------------------------------------

            redundancy_loads: Dict[str, float] = {
                predictor: 0.0 for predictor in predictor_names
            }

            predictor_list = sorted(predictor_names)

            for i in range(len(predictor_list)):
                predictor_a = predictor_list[i]

                for j in range(i + 1, len(predictor_list)):
                    predictor_b = predictor_list[j]

                    pair_stats = self._get_pair(
                        predictor_a,
                        predictor_b,
                    )

                    if pair_stats is None:
                        continue

                    # Require significant predictor-predictor relationship
                    if redundancy_max_fdr is not None:
                        if (
                            pair_stats.fdr_q is None
                            or pair_stats.fdr_q > redundancy_max_fdr
                        ):
                            continue

                    phi = pair_stats.phi

                    if phi is None:
                        continue

                    # Only positive co-occurrence counts as redundancy
                    if phi <= redundancy_min_phi:
                        continue

                    redundancy = max(
                        0.0,
                        phi,
                    )

                    redundancy_loads[predictor_a] += redundancy

                    redundancy_loads[predictor_b] += redundancy

            # ----------------------------------------------------------
            # Build prepared predictor objects
            # ----------------------------------------------------------

            prepared: Dict[
                str,
                PreparedPredictor,
            ] = {}

            for predictor_info in selected:
                predictor = predictor_info.predictor

                bsr_threshold = self.bsr_thresholds.get(predictor)

                if bsr_threshold is None:
                    # Cannot use this predictor for continuous BSR scoring
                    continue

                oriented = self._get_oriented_counts(
                    target,
                    predictor,
                )

                if oriented is None:
                    continue

                (
                    pair_stats,
                    n11,
                    n_target_only,
                    n_predictor_only,
                    n00,
                ) = oriented

                # ------------------------------------------------------
                # P(predictor | target)
                #
                # pseudocount prevents 0 / infinity
                # ------------------------------------------------------

                p_predictor_given_target = (n11 + alpha) / (
                    n11 + n_target_only + 2.0 * alpha
                )

                # ------------------------------------------------------
                # P(predictor | NOT target)
                # ------------------------------------------------------

                p_predictor_given_not_target = (n_predictor_only + alpha) / (
                    n_predictor_only + n00 + 2.0 * alpha
                )

                # ------------------------------------------------------
                # PRESENT likelihood ratio
                # ------------------------------------------------------

                lr_present = p_predictor_given_target / p_predictor_given_not_target

                log_lr_present = math.log(lr_present)

                # ------------------------------------------------------
                # ABSENT likelihood ratio
                # ------------------------------------------------------

                p_absent_given_target = 1.0 - p_predictor_given_target

                p_absent_given_not_target = 1.0 - p_predictor_given_not_target

                lr_absent = p_absent_given_target / p_absent_given_not_target

                log_lr_absent = math.log(lr_absent)

                # ------------------------------------------------------
                # Redundancy weighting
                # ------------------------------------------------------

                redundancy_load = redundancy_loads[predictor]

                weight = 1.0 / (1.0 + redundancy_load)

                prepared[predictor] = PreparedPredictor(
                    target=target,
                    predictor=predictor,
                    n_tested=pair_stats.n_tested,
                    n_both=n11,
                    fdr_q=pair_stats.fdr_q,
                    p_predictor_given_target=(p_predictor_given_target),
                    p_predictor_given_not_target=(p_predictor_given_not_target),
                    log_lr_present=log_lr_present,
                    log_lr_absent=log_lr_absent,
                    redundancy_load=redundancy_load,
                    weight=weight,
                    weighted_log_lr_present=(log_lr_present * weight),
                    weighted_log_lr_absent=(log_lr_absent * weight),
                    # BSR thresholds
                    bsr_lower_limit=(bsr_threshold.lower_limit),
                    bsr_average=(bsr_threshold.average),
                    bsr_upper_limit=(bsr_threshold.upper_limit),
                    bsr_negative_cutoff=(bsr_threshold.negative_cutoff),
                )

            models[target] = PreparedTargetModel(
                target=target,
                predictors=prepared,
            )

        return models

    def print_summary(
        self,
        sort_by: str = "mutual_information",
        descending: bool = True,
        max_pairs: Optional[int] = None,
        max_fdr: Optional[float] = None,
        min_n11: int = 0,
    ) -> None:
        """
        Print pairwise co-occurrence results in readable form.

        Parameters
        ----------
        sort_by
            Attribute of PairwiseAssociation used for sorting.

            Examples:
                "mutual_information"
                "log2_odds_ratio"
                "fold_enrichment"
                "phi"
                "fdr_q"

        descending
            Sort high -> low if True.

        max_pairs
            Maximum number of pairs to print.

        max_fdr
            If set, show only pairs with FDR q <= max_fdr.

        min_n11
            Minimum number of genomes where both domains are present.
        """

        results = list(self.pairs.values())

        # --------------------------------------------------------------
        # Filter
        # --------------------------------------------------------------

        filtered_results = []

        for stats in results:
            if stats.n11 < min_n11:
                continue

            if max_fdr is not None:
                if stats.fdr_q is None or stats.fdr_q > max_fdr:
                    continue

            filtered_results.append(stats)

        results = filtered_results

        # --------------------------------------------------------------
        # Validate sorting attribute
        # --------------------------------------------------------------

        if results and not hasattr(
            results[0],
            sort_by,
        ):
            raise ValueError(f"Unknown PairwiseAssociation attribute: {sort_by!r}")

        # --------------------------------------------------------------
        # Sort
        # --------------------------------------------------------------

        def sort_value(stats):
            value = getattr(
                stats,
                sort_by,
            )

            if value is None:
                if descending:
                    return float("-inf")

                return float("inf")

            return value

        results.sort(
            key=sort_value,
            reverse=descending,
        )

        if max_pairs is not None:
            results = results[:max_pairs]

        # --------------------------------------------------------------
        # Header
        # --------------------------------------------------------------

        print()
        print("Pairwise domain co-occurrence analysis")
        print("=" * 40)

        print(f"Total genomes: {self.n_genomes:,}")

        print(f"Total domain pairs: {len(self.pairs):,}")

        print(f"Shown pairs: {len(results):,}")

        print(f"Sorted by: {sort_by}")

        if max_fdr is not None:
            print(f"Maximum FDR q: {max_fdr}")

        if min_n11 > 0:
            print(f"Minimum shared occurrences: {min_n11}")

        print()

        # --------------------------------------------------------------
        # Individual pair summaries
        # --------------------------------------------------------------

        for i, stats in enumerate(
            results,
            start=1,
        ):
            print(f"[{i}] {stats.format_summary()}")

            print("-" * 70)


# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------


def _get_threshold_value(
    values: Dict[str, float],
    candidate_keys: Tuple[str, ...],
) -> Optional[float]:
    for key in candidate_keys:
        value = values.get(key)
        if value is not None:
            return float(value)

    return None


def _prepare_thresholds(
    score_limit_dict: Dict[str, Dict[str, float]],
    negative_fraction: float = 0.5,
) -> Dict[str, DomainThreshold]:
    """
    Convert score_limit_dict into explicit present / unknown / absent thresholds.

    Supported aliases:

    minimum:
        minimum_cutoff
        min
        lower_limit

    average:
        average_cutoff
        average
        mean
        mean_score
    """

    thresholds: Dict[str, DomainThreshold] = {}

    for domain, values in score_limit_dict.items():
        minimum = _get_threshold_value(
            values,
            (
                "minimum_cutoff",
                "min",
                "lower_limit",
            ),
        )

        average = _get_threshold_value(
            values,
            (
                "average_cutoff",
                "average",
                "mean",
                "mean_score",
            ),
        )

        if minimum is None:
            raise ValueError(f"No minimum cutoff found for domain {domain!r}")

        if average is None:
            print(score_limit_dict)
            raise ValueError(
                f"No average cutoff found for domain {domain!r}. "
                f"The score_limit_dict must contain a mean/average score."
            )

        if average < minimum:
            raise ValueError(
                f"Invalid score limits for {domain}: "
                f"average={average} < minimum={minimum}"
            )

        thresholds[domain] = DomainThreshold(
            minimum_cutoff=minimum,
            average_cutoff=average,
            negative_cutoff=minimum * negative_fraction,
        )

    return thresholds


def _mutual_information(
    n11: int,
    n10: int,
    n01: int,
    n00: int,
) -> Optional[float]:
    """
    Mutual information in bits.
    """

    total = n11 + n10 + n01 + n00

    if total == 0:
        return None

    cells = (
        (1, 1, n11),
        (1, 0, n10),
        (0, 1, n01),
        (0, 0, n00),
    )

    row_counts = {
        1: n11 + n10,
        0: n01 + n00,
    }

    col_counts = {
        1: n11 + n01,
        0: n10 + n00,
    }

    mi = 0.0

    for a, b, count in cells:
        if count == 0:
            continue

        p_ab = count / total
        p_a = row_counts[a] / total
        p_b = col_counts[b] / total

        mi += p_ab * math.log2(p_ab / (p_a * p_b))

    return mi


def _phi_coefficient(
    n11: int,
    n10: int,
    n01: int,
    n00: int,
) -> Optional[float]:
    """
    Phi coefficient / Matthews correlation for binary variables.
    """

    denominator = (n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00)

    if denominator <= 0:
        return None

    numerator = n11 * n00 - n10 * n01

    return numerator / math.sqrt(denominator)


def _benjamini_hochberg(
    p_values: list[float],
) -> list[float]:
    """
    Benjamini-Hochberg FDR correction.
    """

    n = len(p_values)

    if n == 0:
        return []

    order = sorted(
        range(n),
        key=lambda i: p_values[i],
    )

    q_values = [1.0] * n

    previous_q = 1.0

    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1

        q = p_values[original_index] * n / rank

        q = min(
            q,
            previous_q,
            1.0,
        )

        q_values[original_index] = q
        previous_q = q

    return q_values


def _prepare_bsr_thresholds(
    score_limit_dict: Dict[str, Dict[str, float]],
    negative_fraction: float = 0.5,
) -> Dict[str, BSRThreshold]:
    thresholds = {}

    for domain, values in score_limit_dict.items():
        lower = values.get("bsr_lower_limit")
        average = values.get("bsr_average")
        upper = values.get("bsr_upper_limit")

        if lower is None or average is None or upper is None:
            continue

        lower = float(lower)
        average = float(average)
        upper = float(upper)

        thresholds[domain] = BSRThreshold(
            lower_limit=lower,
            average=average,
            upper_limit=upper,
            negative_cutoff=(negative_fraction * lower),
        )

    return thresholds


def select_training_genomes(
    database_path: str,
    max_genomes: int = 10000,
    random_state: int = 42,
) -> Set[str]:
    """
    Select at most max_genomes non-QUERY genomes from the database.

    Selection is random but reproducible for a given random_state.
    """

    if max_genomes < 1:
        raise ValueError("max_genomes must be >= 1")

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

        genome_ids = sorted(str(genome_id) for (genome_id,) in cur.fetchall())

    if len(genome_ids) <= max_genomes:
        return set(genome_ids)

    rng = random.Random(random_state)

    return set(
        rng.sample(
            genome_ids,
            max_genomes,
        )
    )


# ----------------------------------------------------------------------
# Main routine
# ----------------------------------------------------------------------


def calculate_pairwise_domain_cooccurrence(
    database_path: str,
    score_limit_dict: Dict[str, Dict[str, float]],
    negative_fraction: float = 0.5,
    score_field: str = "score",
    exclude_genomes: Optional[Set[str]] = None,
    fisher_alternative: str = "greater",
) -> PairwiseCooccurrenceResult:
    """
    Calculate pairwise genome-level co-occurrence statistics for all domains.

    Classification per domain/genome
    --------------------------------

    best_score >= average_cutoff
        -> PRESENT

    best_score < negative_fraction * minimum_cutoff
        -> ABSENT

    otherwise
        -> UNKNOWN

    Genomes with UNKNOWN classification for either member of a pair are
    excluded only from that specific pairwise comparison.

    No detected hit is treated as ABSENT.

    Parameters
    ----------
    database_path
        SQLite database containing Proteins and Domains tables.

    score_limit_dict
        Domain-specific thresholds.

        Expected structure, e.g.:

        {
            "DsrA": {
                "min": 120.0,
                "mean": 220.0,
            },
            "DsrB": {
                "min": 100.0,
                "mean": 190.0,
            },
        }

        'lower_limit' can also be used instead of 'min'.

    negative_fraction
        Fraction of minimum cutoff used as upper limit for ABSENT.
        Default = 0.5.

    score_field
        Domains column used for classification.
        Currently "score" or "blast_score_ratio".

    exclude_genomes
        GenomeIDs excluded completely.
        Default: {"QUERY"}.

    fisher_alternative
        "greater":
            test specifically for positive co-occurrence.

        "two-sided":
            test for any dependency.

    Returns
    -------
    PairwiseCooccurrenceResult

        result.pairs[("DsrA", "DsrB")]

        returns a PairwiseAssociation object.
    """

    if score_field not in {
        "score",
        "blast_score_ratio",
    }:
        raise ValueError("score_field must be 'score' or 'blast_score_ratio'")

    if fisher_alternative not in {
        "greater",
        "less",
        "two-sided",
    }:
        raise ValueError("Invalid fisher_alternative")

    if exclude_genomes is None:
        exclude_genomes = {"QUERY"}

    thresholds = _prepare_thresholds(
        score_limit_dict=score_limit_dict,
        negative_fraction=negative_fraction,
    )

    bsr_thresholds = _prepare_bsr_thresholds(
        score_limit_dict=score_limit_dict,
        negative_fraction=negative_fraction,
    )

    domains = sorted(thresholds)

    states: Dict[str, DomainState] = {
        domain: DomainState(
            present=set(),
            unknown=set(),
        )
        for domain in domains
    }

    # --------------------------------------------------------------
    # Fetch genomes + best score per genome/domain
    # --------------------------------------------------------------

    with sqlite3.connect(
        database_path,
        timeout=120.0,
    ) as con:
        cur = con.cursor()

        cur.execute("PRAGMA temp_store=MEMORY;")

        cur.execute(
            """
            SELECT DISTINCT genomeID
            FROM Proteins
            WHERE genomeID IS NOT NULL
            """
        )

        genome_ids = sorted(
            genome_id
            for (genome_id,) in cur.fetchall()
            if genome_id not in exclude_genomes
        )

        genome_to_index = {
            genome_id: index for index, genome_id in enumerate(genome_ids)
        }

        # TEMP table avoids huge IN(...) statements
        cur.execute(
            """
            CREATE TEMP TABLE tmp_target_domains (
                domain TEXT PRIMARY KEY
            )
            """
        )

        cur.executemany(
            """
            INSERT INTO tmp_target_domains(domain)
            VALUES (?)
            """,
            ((domain,) for domain in domains),
        )

        query = f"""
            SELECT
                p.genomeID,
                d.domain,
                MAX(d.{score_field}) AS best_score

            FROM Domains AS d

            JOIN Proteins AS p
              ON p.proteinID = d.proteinID

            JOIN tmp_target_domains AS t
              ON t.domain = d.domain

            WHERE p.genomeID IS NOT NULL
              AND d.{score_field} IS NOT NULL

            GROUP BY
                p.genomeID,
                d.domain
        """

        cur.execute(query)

        for genome_id, domain, best_score in cur:
            if genome_id in exclude_genomes:
                continue

            genome_index = genome_to_index.get(genome_id)

            if genome_index is None:
                continue

            threshold = thresholds[domain]
            score = float(best_score)

            # Positive
            if score >= threshold.average_cutoff:
                states[domain].present.add(genome_index)

            # Negative
            elif score < threshold.negative_cutoff:
                # absent is stored implicitly
                pass

            # Gray zone
            else:
                states[domain].unknown.add(genome_index)

    # --------------------------------------------------------------
    # Pairwise statistics
    # --------------------------------------------------------------

    n_total_genomes = len(genome_ids)

    pair_results: Dict[
        Tuple[str, str],
        PairwiseAssociation,
    ] = {}

    pair_keys: list[Tuple[str, str]] = []
    p_values: list[float] = []

    for domain_a, domain_b in combinations(
        domains,
        2,
    ):
        state_a = states[domain_a]
        state_b = states[domain_b]

        pa = state_a.present
        pb = state_b.present

        ua = state_a.unknown
        ub = state_b.unknown

        # ----------------------------------------------------------
        # Number excluded because one or both domains are UNKNOWN
        # ----------------------------------------------------------

        n_unknown_both = len(ua.intersection(ub))

        n_excluded_unknown = len(ua) + len(ub) - n_unknown_both

        n_tested = n_total_genomes - n_excluded_unknown

        # ----------------------------------------------------------
        # 2x2 table
        #
        # Avoid explicitly constructing the huge valid-genome set.
        # ----------------------------------------------------------

        n11 = len(pa.intersection(pb))

        # A present but B unknown cannot be used
        n_a_present_valid = len(pa) - len(pa.intersection(ub))

        # B present but A unknown cannot be used
        n_b_present_valid = len(pb) - len(pb.intersection(ua))

        n10 = n_a_present_valid - n11

        n01 = n_b_present_valid - n11

        n00 = n_tested - n11 - n10 - n01

        # Numerical sanity check
        if (
            min(
                n11,
                n10,
                n01,
                n00,
            )
            < 0
        ):
            raise RuntimeError(
                f"Invalid contingency table for "
                f"{domain_a} / {domain_b}: "
                f"{n11}, {n10}, {n01}, {n00}"
            )

        # ----------------------------------------------------------
        # Conditional probabilities
        # ----------------------------------------------------------

        n_a_present = n11 + n10
        n_b_present = n11 + n01

        p_b_given_a = n11 / n_a_present if n_a_present > 0 else None

        p_a_given_b = n11 / n_b_present if n_b_present > 0 else None

        # ----------------------------------------------------------
        # Expected co-occurrence under independence
        # ----------------------------------------------------------

        if n_tested > 0:
            expected_both = n_a_present * n_b_present / n_tested

        else:
            expected_both = None

        if expected_both is not None and expected_both > 0:
            fold_enrichment = n11 / expected_both

        else:
            fold_enrichment = None

        # ----------------------------------------------------------
        # Jaccard
        # ----------------------------------------------------------

        jaccard_denominator = n11 + n10 + n01

        jaccard = n11 / jaccard_denominator if jaccard_denominator > 0 else None

        # ----------------------------------------------------------
        # Mutual information
        # ----------------------------------------------------------

        mutual_information = _mutual_information(
            n11,
            n10,
            n01,
            n00,
        )

        # ----------------------------------------------------------
        # Phi / Matthews correlation
        # ----------------------------------------------------------

        phi = _phi_coefficient(
            n11,
            n10,
            n01,
            n00,
        )

        # ----------------------------------------------------------
        # Odds ratio
        #
        # Haldane-Anscombe correction makes OR finite even if a cell = 0.
        # ----------------------------------------------------------

        c11 = n11 + 0.5
        c10 = n10 + 0.5
        c01 = n01 + 0.5
        c00 = n00 + 0.5

        odds_ratio = c11 * c00 / (c10 * c01)

        log2_odds_ratio = math.log2(odds_ratio)

        # ----------------------------------------------------------
        # Fisher exact test
        # ----------------------------------------------------------

        if n_tested > 0:
            _, fisher_p = fisher_exact(
                [
                    [n11, n10],
                    [n01, n00],
                ],
                alternative=fisher_alternative,
            )

            fisher_p = float(fisher_p)

        else:
            fisher_p = 1.0

        result = PairwiseAssociation(
            domain_a=domain_a,
            domain_b=domain_b,
            n_total_genomes=n_total_genomes,
            n_tested=n_tested,
            n_excluded_unknown=n_excluded_unknown,
            n11=n11,
            n10=n10,
            n01=n01,
            n00=n00,
            n_a_present=n_a_present,
            n_b_present=n_b_present,
            p_b_given_a=p_b_given_a,
            p_a_given_b=p_a_given_b,
            expected_both=expected_both,
            fold_enrichment=fold_enrichment,
            jaccard=jaccard,
            mutual_information=mutual_information,
            phi=phi,
            odds_ratio=odds_ratio,
            log2_odds_ratio=log2_odds_ratio,
            fisher_p=fisher_p,
        )

        key = (
            domain_a,
            domain_b,
        )

        pair_results[key] = result
        pair_keys.append(key)
        p_values.append(fisher_p)

    # --------------------------------------------------------------
    # Multiple-testing correction
    # --------------------------------------------------------------

    q_values = _benjamini_hochberg(p_values)

    for key, q_value in zip(
        pair_keys,
        q_values,
    ):
        pair_results[key].fdr_q = q_value

    return PairwiseCooccurrenceResult(
        n_genomes=n_total_genomes,
        thresholds=thresholds,
        bsr_thresholds=bsr_thresholds,
        pairs=pair_results,
    )
