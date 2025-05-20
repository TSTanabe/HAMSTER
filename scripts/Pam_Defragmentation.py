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
def pam_defragmentation_finder(options, basis_grouped, basis_score_limit_dict):
    
    # Load precomputed grp1 results if available
    grp1_merged_dict = myUtil.load_cache(options,'grp1_merged_grouped.pkl')
    if grp1_merged_dict:
        return grp1_merged_dict
    
    # Eine Genomliste zurückbekommen
    predicted_genomes = predict_reference_seqs_for_each_domain(database_path, grouped, cores, chunk_size=900)
    
    # Collect die zusätzlichen proteinIDs aus den entsprechenden Genomen der DB
    pam_added_score_limit_dict, pam_added_reference_seq_dict = collect_predicted_singleton_hits_from_db_parallel(predictions_all, options.database_directory, plausible_cutoff = 0.02, bsr_cutoff=0.5) #TODO diese routine könnte auch in das allgemeine PAM Modul, weil sie doch sehr allgemein ist

    # Merge the pam added seqs with the basis under the grp1 prefix key    
    grp1_added_reference_seq_dict = {f"grp1_{k}": v for k, v in pam_added_reference_seq_dict.items()}
    grp1_added_score_limits_dict = {f"grp1_{k}": v for k, v in pam_added_score_limit_dict.items()}
    
    grp1_basis_grouped_dict = {k.replace("grp0_", "grp1_"): v for k, v in basis_grouped.items()}
    grp1_score_limits_dict = {k.replace("grp0_", "grp1_"): v for k, v in basis_score_limit_dict.items()}

    grp1_merged_dict = merge_dicts_of_sets(grp1_added_reference_seq_dict,grp1_basis_grouped_dict)
    grp1_merged_score_dict = merge_dicts_of_sets(grp1_added_score_limits_dict,grp1_score_limits_dict)

    report_added_only_counts(grp1_merged_dict, grp1_added_reference_seq_dict, grp1_basis_grouped_dict) # print in terminal the number of added sequences

    myUtil.save_cache(options, 'grp1_merged_grouped.pkl', grp1_merged_dict)
    myUtil.save_cache(options, 'grp1_merged_score_limits.pkl', grp1_merged_score_dict)
    
    return grp1_merged_dict


#######################

def predict_reference_seqs_for_each_domain(database_path, grouped, cores, chunk_size=900):
    """
    Für jedes Singleton wird ein Modell auf einem kombinierten PAM (grouped + singleton_domain) trainiert.
    Danach wird das Modell verwendet, um Vorhersagen auf grouped zu treffen.
    
    Args:
        database_path (str): Pfad zur SQLite-Datenbank.
        grouped (dict): Basis-Domänenstruktur {domain_label: set(proteinIDs)}.
        singleton (dict): {singleton_domain_label: set(proteinIDs)}.
        cores (int): Anzahl CPUs.
        chunk_size (int): DB-Abfragegröße.

    Returns:
        dict: {singleton_domain: prediction_series (genomeID → score)}
    """

    predictions_all = {}

    # 1. Presence/Absence-Matrix erzeugen
    global_pam = Pam_Mx_algorithm.create_presence_absence_matrix(
        grouped_plus,
        database_directory=database_path,
        output=f"pam_with_{sng_domain}",
        chunk_size=chunk_size,
        cores=cores
    )

    # 2. Load BSR-Scores
    bsr_hit_scores = Pam_Mx_algorithm.fetch_bsr_scores(database_path, grouped, chunk_size=chunk_size)
    
    # 3. Train the logistic regression models
    models, _ = Pam_Mx_algorithm.train_logistic_from_pam_with_scores(global_pam, bsr_hit_scores, cores=cores)
    
    # 4.
    for domain in grouped:
        print(f"[INFO] Processing presence absence prediction for domain {domain}")
        
        model = models.get(domain)
        
        if model is None:
            print(f"[WARN] No logistic regression model found for domain {domain}, skipping.")
            continue
    
        # 4.1 Remove the domain column from the PAM
        grouped_minus_domain = {d: grouped[d] for d in grouped if d != domain}
        # TODO kann man das nicht auch aus der PAM einfach rausstreichen anstatt immer neue PAMs ohne die eine spalte zu erstellen?
        test_pam = Pam_Mx_algorithm.create_presence_absence_matrix(
            grouped_minus_domain,
            database_directory=database_path,
            output=f"pam_without_{domain}",
            chunk_size=chunk_size,
            cores=cores
        )
        # 4.2 Feature-DataFrame bauen
        genomes = sorted(test_pam.keys())
        features = sorted({d for doms in test_pam.values() for d in doms})
        df = pd.DataFrame(index=genomes, columns=features)

        for genome_id in genomes:
            for feat in features:
                df.at[genome_id, feat] = ",".join(test_pam[genome_id].get(feat, [])) if feat in test_pam[genome_id] else ""

        # 4.3 Feature-Matrix für Modell
        X = Pam_Mx_algorithm.build_presence_score_matrix(df, bsr_hit_scores)

        # 4.4 Sicherstellen, dass alle Modell-Features vorhanden sind
        X = X.reindex(columns=model.feature_names_in_, fill_value=0)

        # 4.5 Prediction
        preds = pd.Series(model.predict_proba(X)[:, 1], index=X.index)
        predictions_all[domain] = preds

    return predictions_all



########### TODO TODO TODO TODO TODO ####################

#move this processing to the PAM_MX_ALGORITHM.py
def process_domain_hits(domain, genome_list, database_path, bsr_threshold, return_hits, return_limits):
    best_hits = {}
    con = sqlite3.connect(database_path)
    cur = con.cursor()

    for i in range(0, len(genome_list), 900):
        chunk = genome_list[i:i + 900]
        placeholders = ','.join(['?'] * len(chunk))
        query = f"""
            SELECT genomeID, Proteins.proteinID, score, blast_score_ratio
            FROM Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID
            WHERE blast_score_ratio > ?
              AND domain = ?
              AND genomeID IN ({placeholders})
        """
        cur.execute(query, (bsr_threshold, domain, *chunk))
        for genome_id, proteinID, score, bsr in cur.fetchall():
            if genome_id not in best_hits or bsr > best_hits[genome_id][0]:
                best_hits[genome_id] = (bsr, proteinID, score)
    
    con.close()

    # Aggregate result
    proteinIDs = set()
    bitscores = []
    for genome_id, (bsr, proteinID, score) in best_hits.items():
        proteinIDs.add(proteinID)
        bitscores.append(score)

    if proteinIDs:
        return_hits[domain] = proteinIDs
    if bitscores:
        return_limits[domain] = {"lower_limit": min(bitscores), "upper_limit": max(bitscores)}

def collect_predicted_singleton_hits_from_db_parallel(predictions_all, database_path, plausible_cutoff=0.6, bsr_cutoff=0.5):
    from multiprocessing import Process, Manager

    singleton_hits = defaultdict(set)
    singleton_score_limits = {}

    domain_to_genomes = defaultdict(list)
    for singleton, prediction_series in predictions_all.items():
        domain = singleton.split('_', 1)[1]
        plausible_genomes = prediction_series[prediction_series > plausible_cutoff].index
        domain_to_genomes[domain].extend(plausible_genomes)

    manager = Manager()
    return_hits = manager.dict()
    return_limits = manager.dict()

    jobs = []
    for domain, genome_list in domain_to_genomes.items():
        p = Process(target=process_domain_hits, args=(
            domain, genome_list, database_path, bsr_threshold, return_hits, return_limits
        ))
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()

    # Merge Manager dicts to regular dicts
    singleton_hits = {d: set(p) for d, p in return_hits.items()}
    singleton_score_limits = dict(return_limits)

    return singleton_score_limits, singleton_hits


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
        print(f"[INFO] For {key}: {len(added_only)} sequences were added due to the presence absence matrix")



    
