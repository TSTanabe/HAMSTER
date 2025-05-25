#!/usr/bin/python

import copy
import os
import subprocess
import sqlite3
import traceback
import pandas as pd
from collections import defaultdict

from . import Pam_Mx_algorithm
from . import myUtil

    

#### Main routine of this module
def pam_genome_defragmentation_hit_finder(options, basis_grouped, basis_score_limit_dict):
    
    # Load precomputed grp1 results if available
    grp1_merged_dict = myUtil.load_cache(options,'grp1_merged_grouped.pkl')
    if grp1_merged_dict:
        return grp1_merged_dict
    
    # Eine Genomliste zurückbekommen
    predicted_genomes = predict_reference_seqs_for_each_domain(options.database_directory, basis_grouped, options.cores, chunk_size=900)
    
    # Collect die zusätzlichen proteinIDs aus den entsprechenden Genomen der DB
    pam_added_score_limit_dict, pam_added_reference_seq_dict = Pam_Mx_algorithm.collect_predicted_singleton_hits_from_db_parallel(predicted_genomes, options.database_directory, options.pam_threshold, options.pam_bsr_threshold)
    
    # Merge the new proteinIDs to the basic set to generate the grp1 refseq dataset
    grp1_merged_dict = myUtil.merge_grouped_refseq_dicts_simple(pam_added_reference_seq_dict, basis_grouped)
    grp1_merged_score_dict = myUtil.merge_score_limits(pam_added_score_limit_dict, basis_score_limit_dict)

    report_added_only_counts(grp1_merged_dict, pam_added_reference_seq_dict, basis_grouped) # print in terminal the number of added sequences
    
    return grp1_merged_dict, grp1_merged_score_dict


#######################

def predict_reference_seqs_for_each_domain(database_path, grouped, cores, chunk_size=900):
    """
    Für jede Domäne wird ein Modell trainiert und nur neue Präsenz-Vorhersagen gemacht (also für Genome,
    in denen die Domäne laut PAM nicht vorkommt).
    
    Args:
        database_path (str): Pfad zur SQLite-Datenbank.
        grouped (dict): Basis-Domänenstruktur {domain_label: set(proteinIDs)}.
        cores (int): Anzahl CPUs.
        chunk_size (int): DB-Abfragegröße.

    Returns:
        dict: {domain: prediction_series (genomeID → score), nur neue Genome}
    """

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
        print(f"[INFO] Processing prediction for domain {domain}")
        
        model = models.get(domain)
        if model is None:
            print(f"[WARN] No model found for domain {domain}, skipping.")
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

def merge_dicts_of_sets(dict1, dict2):
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
        print(f"[INFO] {len(added_only)} {key} sequences were added due to the presence absence matrix")



    
