#!/usr/bin/python


from typing import Dict, Set, Any

from src.selection_defragmentation import (
    pam_mx_selection,
    pam_mx_calculate_plausability,
)
from src.core import myUtil

from src.core.logging import get_logger

logger = get_logger(__name__)


#### Main routine of this module
def pam_genome_defragmentation_hit_finder(
    options: Any, basis_grouped: Dict[str, Set[str]], plausability_cutoff: float, support_models_name: str
) -> Dict[str, Set[str]]:
    """
    Main routine: For each domain, finds additional genomes/proteins via PAM-based logistic regression, and merges with the basis set.

    Args:
        options (Any): Config/options object.
        basis_grouped (dict): {domain: set(proteinIDs)} basic set.
        basis_score_limit_dict (dict): (not used in logic here, for compatibility).

    Returns:
        dict: {domain: set(proteinIDs)}, merged set.
    """

    support_models = myUtil.load_cache(options, support_models_name)
    if support_models is None:
        support_models = pam_mx_selection.train_support_models_for_each_domain(
            database_path=options.database_directory,
            grouped=basis_grouped,
            cores=options.cores,
            chunk_size=900,
            max_genomes=1000,
            random_seed=42,
            alpha=1.0,
            min_feature_genomes=2,
        )
        #myUtil.save_cache(options, support_models_name, support_models)


    plausible_hits = (
        pam_mx_calculate_plausability.collect_plausible_domain_hits_from_support_models(
            support_models=support_models,
            database_path=options.database_directory,
            chunk_size=900,
            plausibility_cutoff=plausability_cutoff,
            score_field="score",
        )
    )

    # Merge the new proteinIDs to the basic set to generate the grp1 refseq dataset
    merged_dict = myUtil.merge_grouped_refseq_dicts_simple(
        plausible_hits, basis_grouped
    )

    report_added_only_counts(
        merged_dict, plausible_hits, basis_grouped
    )  # print in terminal the number of added sequences

    myUtil.save_cache(options, "grp1_pam_defragmented_dict.pkl", merged_dict)
    return merged_dict


def report_added_only_counts(
    grp1_merged_dict, grp1_added_reference_seq_dict, grp1_basis_grouped_dict
):
    for key in grp1_merged_dict:
        added_set = grp1_added_reference_seq_dict.get(key, set())
        basis_set = grp1_basis_grouped_dict.get(key, set())

        added_only = added_set - basis_set
        logger.info(
            f"{len(added_only)} {key} sequences were added due to the presence propability"
        )
