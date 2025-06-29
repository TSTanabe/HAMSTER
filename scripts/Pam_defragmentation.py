#!/usr/bin/python


import copy
import os
import sqlite3
from typing import Dict, Set, Any, Tuple
import pandas as pd

from . import Pam_Mx_algorithm
from . import myUtil

logger = myUtil.logger

    

#### Main routine of this module
def pam_genome_defragmentation_hit_finder(
    options: Any,
    basis_grouped: Dict[str, Set[str]],
    basis_score_limit_dict: Dict
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
    
    # returns dict: {domain: prediction_series (genomeID → score), only new genomes}
    predicted_genomes = predict_reference_seqs_for_each_domain(options.database_directory, basis_grouped, options.cores, chunk_size=900)
    
    # Collect die zusätzlichen proteinIDs aus den entsprechenden Genomen der DB
    pam_added_reference_seq_dict = Pam_Mx_algorithm.collect_predicted_singleton_hits_from_db_parallel(predicted_genomes, options.database_directory, options.pam_threshold, options.pam_bsr_threshold)
    
    # Merge the new proteinIDs to the basic set to generate the grp1 refseq dataset
    grp1_merged_dict = myUtil.merge_grouped_refseq_dicts_simple(pam_added_reference_seq_dict, basis_grouped)

    report_added_only_counts(grp1_merged_dict, pam_added_reference_seq_dict, basis_grouped) # print in terminal the number of added sequences
    
    return grp1_merged_dict


#######################


def predict_reference_seqs_for_each_domain(
    database_path: str,
    grouped: Dict[str, Set[str]],
    cores: int,
    chunk_size: int = 900
) -> Dict[str, pd.Series]:
    """
    For each domain, train a logistic regression model and predict probability of presence for genomes NOT already in grouped.

    Args:
        database_path (str): Path to SQLite DB.
        grouped (dict): {domain: set(proteinIDs)}.
        cores (int): Number of CPUs.
        chunk_size (int): DB query chunk size.

    Returns:
        dict: {domain: pd.Series (genomeID → probability)} for new genomes only.
    """
    logger.info(f"Training logistic regression {len(grouped.values())} for proteins")
    predictions_all = {}

    # 1. PAM berechnen
    global_pam = Pam_Mx_algorithm.create_presence_absence_matrix(
        grouped,
        database_directory=database_path,
        output="pam_with_all_domains",
        chunk_size=chunk_size,
        cores=cores
    )

    # 2. BSR-Scores laden
    bsr_hit_scores = Pam_Mx_algorithm.fetch_bsr_scores(database_path, grouped, chunk_size=chunk_size)
    
    # 3. Modelle trainieren
    models, _ = Pam_Mx_algorithm.train_logistic_from_pam_with_scores(global_pam, bsr_hit_scores, cores=cores)
    
    # 4. Vorhersagen je Domain
    for domain in grouped:
        logger.debug(f"Fitting logistic regression model for {domain}")

        
        model = models.get(domain)
        if model is None:
            logger.warning(f"No model found for domain {domain}, skipping.")
            continue
    
        # 4.1 PAM ohne domain
        test_pam = copy.deepcopy(global_pam)
        for genome_id in test_pam:
            test_pam[genome_id].pop(domain, None)

        # 4.2 DataFrame
        genomes = sorted(test_pam.keys())
        features = sorted({d for doms in test_pam.values() for d in doms})
        df = pd.DataFrame(index=genomes, columns=features)

        for genome_id in genomes:
            for feat in features:
                df.at[genome_id, feat] = ",".join(test_pam[genome_id].get(feat, [])) if feat in test_pam[genome_id] else ""

        # 4.3 Matrix
        X = Pam_Mx_algorithm.build_presence_score_matrix(df, bsr_hit_scores)

        # 4.4 Spalten angleichen
        X = X.reindex(columns=model.feature_names_in_, fill_value=0)

        # 4.5 Vorhersage nur für Genome ohne bisherige Präsenz, keine hits doppelt holen wegen PAM
        raw_preds = pd.Series(model.predict_proba(X)[:, 1], index=X.index)
        novel_genomes = [gid for gid in raw_preds.index if domain not in global_pam.get(gid, {})]
        preds = raw_preds.loc[novel_genomes]

        if not preds.empty:
            predictions_all[domain] = preds

    return predictions_all



################ Helper routines and print routines ###################

def merge_dicts_of_sets(
    dict1: Dict[str, Set[str]],
    dict2: Dict[str, Set[str]]
) -> Dict[str, Set[str]]:
    """
    Merges two dicts of sets by keywise union.

    Args:
        dict1, dict2: {key: set(values)}

    Returns:
        dict: {key: union of all values}

    Example:
        {'A':{'x'}}, {'A':{'y'}, 'B':{'z'}} → {'A':{'x','y'}, 'B':{'z'}}
    """
    merged = {}

    all_keys = set(dict1) | set(dict2)  # Vereinigung aller Keys
    for key in all_keys:
        merged[key] = dict1.get(key, set()) | dict2.get(key, set())

    return merged

def report_added_only_counts(grp1_merged_dict, grp1_added_reference_seq_dict, grp1_basis_grouped_dict):
    for key in grp1_merged_dict:
        added_set = grp1_added_reference_seq_dict.get(key, set())
        basis_set = grp1_basis_grouped_dict.get(key, set())

        added_only = added_set - basis_set
        logger.info(f"{len(added_only)} {key} sequences were added due to the presence propability")



    
