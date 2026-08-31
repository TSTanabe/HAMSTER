#!/usr/bin/python

from typing import Dict, Set, Any

from src.selection_defragmentation import (
    mutual_information,
    predictor_calibration,
    predictor_application,
)

from src.core import myUtil
from src.core.logging import get_logger


logger = get_logger(__name__)


def _predictor_models_are_calibrated(
    predictor_models: Dict[str, Any],
) -> bool:
    if not predictor_models:
        return False

    for model in predictor_models.values():
        if not model.predictors:
            continue

        if model.calibration_intercept is None or model.calibration_coefficient is None:
            return False

    return True


def get_or_train_predictor_models(
    config: Any,
    basis_score_limit: Dict[
        str,
        Dict[str, float],
    ],
    support_models_name: str,
) -> Dict[str, Any]:
    predictor_models = myUtil.load_cache(
        config,
        support_models_name,
    )

    if predictor_models is not None and _predictor_models_are_calibrated(
        predictor_models
    ):
        logger.info(f"Loaded predictor models from {support_models_name}")

        return predictor_models

    # --------------------------------------------------------------
    # Pairwise statistics
    # --------------------------------------------------------------

    statistics = mutual_information.calculate_pairwise_domain_cooccurrence(
        database_path=(config.database_directory),
        score_limit_dict=(basis_score_limit),
        negative_fraction=0.5,
    )

    # --------------------------------------------------------------
    # Final predictor models
    # --------------------------------------------------------------

    predictor_models = statistics.prepare_predictor_models(
        max_fdr=0.05,
        min_log2_odds_ratio=2.0,
        min_both=1,
        positive_only=True,
        redundancy_max_fdr=0.05,
        redundancy_min_phi=0.3,
        alpha=0.5,
    )

    # --------------------------------------------------------------
    # OOF calibration
    # --------------------------------------------------------------

    calibration_results = predictor_calibration.calibrate_predictor_models(
        database_path=(config.database_directory),
        score_limit_dict=(basis_score_limit),
        n_folds=3,
        negative_fraction=0.5,
        max_fdr=0.05,
        min_log2_odds_ratio=2.0,
        min_both=1,
        positive_only=True,
        redundancy_max_fdr=0.05,
        redundancy_min_phi=0.3,
        alpha=0.5,
        min_calibration_samples=10,
        min_positive=3,
        min_negative=3,
        store_samples=False,
    )

    predictor_calibration.attach_calibrations_to_models(
        predictor_models=(predictor_models),
        calibration_results=(calibration_results),
    )

    # --------------------------------------------------------------
    # Persist fully trained models
    # --------------------------------------------------------------

    myUtil.save_cache(
        config,
        support_models_name,
        predictor_models,
    )

    return predictor_models


def apply_predictor_models(
    config: Any,
    predictor_models: Dict[str, Any],
    basis_seed_sequences: Dict[
        str,
        Set[str],
    ],
    basis_score_limit: Dict[
        str,
        Dict[str, float],
    ],
    probability_cutoff: float,
) -> Dict[str, Set[str]]:
    application_result = predictor_application.evaluate_predictor_models(
        database_path=(config.database_directory),
        predictor_models=(predictor_models),
        seed_sequences=(basis_seed_sequences),
        score_limit_dict=(basis_score_limit),
        probability_cutoff=(probability_cutoff),
        application_fraction=0.5,
        missing_as_absent=True,
        store_predictions=False,
    )

    application_result.print_summary()

    plausible_hits: Dict[
        str,
        Set[str],
    ] = {
        str(domain): set(protein_ids)
        for (
            domain,
            protein_ids,
        ) in (application_result.accepted_hits.items())
    }

    return plausible_hits


def predictor_training_calibration_application(
    config: Any,
    basis_seed_sequences: Dict[
        str,
        Set[str],
    ],
    basis_score_limit: Dict[
        str,
        Dict[str, float],
    ],
    probability_cutoff: float, # probability for presence at which new hits are accepted as TP
    support_models_name: str,
) -> Dict[str, Set[str]]:
    predictor_models = get_or_train_predictor_models(
        config=config,
        basis_score_limit=(basis_score_limit),
        support_models_name=(support_models_name),
    )

    plausible_hits = apply_predictor_models(
        config=config,
        predictor_models=(predictor_models),
        basis_seed_sequences=(basis_seed_sequences),
        basis_score_limit=(basis_score_limit),
        probability_cutoff=(probability_cutoff),
    )

    return plausible_hits
