#!/usr/bin/python

import os
import csv
import glob
import hashlib
import pickle
import pandas as pd

import numpy as np

from collections import Counter, defaultdict

from . import myUtil
from . import Reports

import matplotlib
matplotlib.use('Agg')    # <— ganz oben, vor pyplot!
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

def process_initial_validations(options,
                                report_directory: str,
                                init_val_dir: str,
                                db_path: str,
                                output_dir: str = None):
    """
    Processes all initial validation .pkl reports in a directory:
    - Loads reference sets for grp1 and grp2
    - Unpacks each .pkl report into a DataFrame
    - Fetches genetic neighborhood for each hit ID
    - Adds neighborhood info to the DataFrame
    - Saves enriched DataFrame as CSV in output_dir

    Parameters:
        init_val_dir (str): Directory containing '*.ini_performance_pkl' files
        db_path (str): Path to the SQLite database for neighborhood queries
        output_dir (str): Directory to save enriched CSVs. Defaults to init_val_dir/enriched_reports
        
    Info:
        Cannot be parallelized because the neighbourhood fetching is already a forked
        process.
    """
    # Prepare output directory
    if output_dir is None:
        output_dir = os.path.join(report_directory, 'Reports')
    os.makedirs(output_dir, exist_ok=True)

    # Load reference sets
    grp1_refs = myUtil.load_cache(options, 'grp1_merged_grouped.pkl')
    grp2_refs = myUtil.load_cache(options, 'grp2_selection_ref_seqs.pkl')
    
    # Load mcl cutoff sets
    mcl_grp1_cutoff = myUtil.load_cache(options, 'mcl_grp1_cluster_selection_cutoffs.pkl')
    mcl_grp2_cutoff = myUtil.load_cache(options, 'mcl_grp2_cluster_selection_cutoffs.pkl')
    
    # Precompute unified sets (safe if empty)
    all_grp1 = set().union(*grp1_refs.values()) if grp1_refs else set()
    all_grp2 = set().union(*grp2_refs.values()) if grp2_refs else set()

    # Find all initial performance pkl files
    pattern = os.path.join(init_val_dir, '*.ini_performance_pkl')
    pkl_files = glob.glob(pattern)

    # Find all mcl files
    mcl_results = myUtil.load_cache(options, 'mcl_clustering_results.pkl')
    
    # Collect cutoffs
    cutoff_collection = {}
    performance_collection = {}
    
    for pkl_file in pkl_files:
    
        protein = os.path.splitext(os.path.basename(pkl_file))[0]
        hmm_report_output_dir = os.path.join(output_dir, f"{protein}")
        os.makedirs(hmm_report_output_dir, exist_ok=True)
        
        cutoffs_file = os.path.join(init_val_dir, f"{protein}.ini_cutoffs_pkl")
        
        print(f"\n[INFO] Generating HMM training report for {protein}")
        
        #### Make the dataframe ###
        
        # Load report dict: hit_id -> metrics dict
        report_dict = myUtil.load_cache(options, "Report metrics", pkl_file)
        
        # Load reported cutoffs
        cutoffs = myUtil.load_cache(options, "Report", file_path=cutoffs_file)
        hmm_protein_name, optimized_cutoff, trusted_cutoff, noise_cutoff, best_MCC, best_matrix = cutoffs
        cutoffs_dict = {'trusted cutoff':trusted_cutoff, 'noise cutoff':noise_cutoff, 'optimized cutoff':optimized_cutoff}
        cutoff_collection[hmm_protein_name] = cutoffs_dict
        performance_collection[hmm_protein_name] = {'MCC': best_MCC, 'Matrix [TP,FP,FN,TN]': best_matrix}
        
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(report_dict, orient='index')
        df.index.name = 'hit_id'
        
        # Check validity of the dataframe
        required_columns = ['bitscore', 'MCC', 'TP', 'FP', 'FN', 'TN']

        if df.empty:
            print("[ERROR] Report for {protein} is missing.")
            continue
        else:

            # Check, ob alle erforderlichen Spalten vorhanden sind
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                print(f"[ERROR] {missing_columns} Report for {protein} is incomplete.")
                continue
        
        # Fetch neighborhood information
        hit_ids = df.index.tolist()
        neighbors, _ = Reports.fetch_neighbouring_genes_with_domains(db_path, hit_ids)

        # Map neighborhood to DataFrame
        df['neighborhood'] = df.index.map(lambda hid: neighbors.get(hid, [['no_data']]))

        # Flag grp1/grp2 membership
        df['in_grp1'] = df.index.isin(all_grp1)
        df['in_grp2'] = df.index.isin(all_grp2)

        # Save enriched DataFrame
        out_csv = os.path.join(hmm_report_output_dir, f"{protein}_enriched.txt")
        df.to_csv(out_csv, sep='\t', index=True)
        print(f"[SAVE] Saved validation report: {out_csv}")
    
        # ROC curve plotting
        roc_path = os.path.join(hmm_report_output_dir, f"{protein}_roc.txt")
        write_roc_from_metrics_to_tsv(df, score_col='bitscore', trusted_cutoff=trusted_cutoff, noise_cutoff = noise_cutoff, optimized_cutoff = optimized_cutoff, output_path=roc_path)

        # MCC curve per score
        mcc_path = os.path.join(hmm_report_output_dir, f"{protein}_mcc.txt")
        export_existing_mcc_curve(df, score_col="bitscore", mcc_col="MCC", output_path=mcc_path)
        
        # Confusion stats per gene neighbourhood

        for name,cutoff in cutoffs_dict.items():
            write_neighborhood_confusion_data(df, cutoff, name, output_dir = hmm_report_output_dir)

        # MCL selection plots
        protein_name = protein.split('_', 1)[-1]

        # Check if all required data for MCL plotting is present
        if (
            protein_name in mcl_results and
            protein_name in grp1_refs and
            protein_name in mcl_grp1_cutoff and
            protein_name in grp2_refs and
            protein_name in mcl_grp2_cutoff
        ):
            mcl_file = mcl_results[protein_name]
            
            grp1_ref_seqs_set = grp1_refs[protein_name]
            grp1_cuts = mcl_grp1_cutoff[protein_name]
            write_mcl_vs_references(
                mcl_file,
                "grp1_mcl_cluster_selection",
                grp1_ref_seqs_set,
                grp1_cuts['density_threshold'],
                grp1_cuts['reference_threshold'],
                output_dir=hmm_report_output_dir
            )
            print(f"[SAVE] MCL selection plot grp1")

            grp2_ref_seqs_set = grp2_refs[protein_name]
            grp2_cuts = mcl_grp2_cutoff[protein_name]
            write_mcl_vs_references(
                mcl_file,
                "grp2_mcl_cluster_selection",
                grp2_ref_seqs_set,
                grp2_cuts['density_threshold'],
                grp2_cuts['reference_threshold'],
                output_dir=hmm_report_output_dir
            )
            print(f"[SAVE] MCL selection plot grp2")

        else:
            print(f"[WARN] Skipping MCL selection plots for {protein_name}: missing clustering or cutoff data.")
        
    save_cutoffs_table(cutoff_collection, options.Hidden_markov_model_directory, 'cutoffs.txt')
    save_performance_table(performance_collection, options.Hidden_markov_model_directory, 'performance.txt')
    return


def write_roc_from_metrics_to_tsv(df: pd.DataFrame,
                          score_col: str = 'bitscore',
                          trusted_cutoff: float = None,
                          noise_cutoff: float = None,
                          optimized_cutoff: float = None,
                          output_path: str = 'roc_data.txt') -> None:
    """
    Writes ROC curve data (FPR, TPR, score) to a TSV file.
    Optionally adds labels for the closest points above specified cutoffs.

    Parameters:
        df (pd.DataFrame): Must contain:
            - score_col: continuous classifier score
            - columns: 'TP', 'FP', 'FN', 'TN'
        trusted_cutoff (float): Score threshold for 'Trusted'
        noise_cutoff (float): Score threshold for 'Noise'
        optimized_cutoff (float): Score threshold for 'Optimized'
        output_path (str): Path to save the TSV file
    """
    if os.path.isfile(output_path):
        return

    scores = df[score_col].values
    TP = df['TP'].values
    FP = df['FP'].values
    FN = df['FN'].values
    TN = df['TN'].values

    with np.errstate(divide='ignore', invalid='ignore'):
        tprs = TP / (TP + FN)
        fprs = FP / (FP + TN)
    tprs = np.nan_to_num(tprs)
    fprs = np.nan_to_num(fprs)

    unique_scores = np.unique(scores)[::-1]
    roc_points = []
    for s in unique_scores:
        idxs = np.where(scores == s)[0]
        roc_points.append((fprs[idxs].mean(), tprs[idxs].mean(), s))

    fpr_vals, tpr_vals, thresh_vals = zip(*roc_points)

    # Create label column
    labels = [''] * len(thresh_vals)

    def label_cutoff(cutoff, label):
        if cutoff is None:
            return
        cands = [s for s in unique_scores if s >= cutoff]
        s_sel = min(cands) if cands else max(unique_scores)
        i = thresh_vals.index(s_sel)
        labels[i] = label

    label_cutoff(trusted_cutoff, 'Trusted')
    label_cutoff(noise_cutoff, 'Noise')
    label_cutoff(optimized_cutoff, 'Optimized')

    out_df = pd.DataFrame({
        'FPR': fpr_vals,
        'TPR': tpr_vals,
        'Score': thresh_vals,
        'Label': labels
    })

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    out_df.to_csv(output_path, sep='\t', index=False)
    print(f"[SAVE] Saved ROC data to {output_path}")



def write_neighborhood_confusion_data(df: pd.DataFrame,
                                       cutoff: float,
                                       cutoff_name: str = 'cutoff',
                                       neighborhood_col: str = 'neighborhood',
                                       score_col: str = 'bitscore',
                                       true_label_col: str = 'true_value',
                                       top_n: int = None,
                                       output_dir: str = None) -> None:
    """
    Writes confusion matrix summary data per genomic neighborhood at a specific cutoff to a TSV file.

    Parameters:
        df (pd.DataFrame): Must contain score, true labels, and neighborhood columns.
        cutoff (float): Score threshold to classify predictions as positive (>=) or negative (<).
        cutoff_name (str): Descriptive name for the cutoff (used in titles/filenames).
        neighborhood_col (str): Column with neighborhood lists.
        score_col (str): Column with classifier scores.
        true_label_col (str): Column with true labels ('TP' or 'TN').
        top_n (int, optional): Number of top neighborhoods to include by total count. If None, include all.
        output_dir (str, optional): Directory to save the file. Defaults to cwd if None.
    """
    out_dir = output_dir or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    fname = f'neighborhood_confusion_{cutoff_name}.tsv'
    path = os.path.join(out_dir, fname)
    if os.path.isfile(path):
        return

    for col in [neighborhood_col, score_col, true_label_col]:
        if col not in df.columns:
            raise KeyError(f"DataFrame must contain column '{col}'")

    preds_pos = df[score_col] >= cutoff
    true_pos = df[true_label_col] == 'TP'
    df['_conf_label'] = np.where(preds_pos & true_pos, 'TP',
                          np.where(preds_pos & ~true_pos, 'FP',
                          np.where(~preds_pos & true_pos, 'FN', 'TN')))

    def canonical(nb):
        seq = tuple(tuple(x) for x in (nb or []))
        rev = tuple(reversed(seq))
        return seq if seq <= rev else rev

    df['_nb'] = df[neighborhood_col].map(canonical)

    counts = {nb: Counter(sub_df['_conf_label']) for nb, sub_df in df.groupby('_nb')}
    totals = {nb: sum(cnt.values()) for nb, cnt in counts.items()}

    selected = sorted(totals, key=totals.get, reverse=True)
    if top_n is not None:
        selected = selected[:top_n]

    out_rows = []
    for nb in selected:
        cluster_genes = [gene for tup in nb for gene in tup]
        out_rows.append({
            'Neighborhood': " ".join(cluster_genes),
            'TP': counts[nb].get('TP', 0),
            'FP': counts[nb].get('FP', 0),
            'FN': counts[nb].get('FN', 0),
            'TN': counts[nb].get('TN', 0),
            'Total': totals[nb]
        })

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(path, sep='\t', index=False)
    print(f"[SAVE] Saved neighborhood confusion data to {path}")

    df.drop(columns=['_conf_label','_nb'], inplace=True)
    
    
    
def write_mcl_vs_references(mcl_file: str,
                             name: str,
                             ref_ids: set,
                             ref_line_pct: float = 0.1,
                             pct_color_threshold: float = 0.05,
                             max_pct: int = 10,
                             output_dir: str = None) -> None:
    """
    Writes cluster size vs. number of reference sequences for a single MCL clusters file to a TSV file.

    Parameters:
        mcl_file (str): Path to the MCL cluster file.
        ref_ids (set): Set of reference sequence IDs.
        ref_line_pct (float): Fractional threshold for the reference line (e.g., 0.1 for 10%).
        pct_color_threshold (float): Fractional threshold for flagging low-ref points.
        max_pct (int): Maximum percent value for classification.
        output_dir (str, optional): Directory in which to save the TSV. Defaults to cwd.
    """
    # File name prefix and output location
    protein_type = os.path.splitext(os.path.basename(mcl_file))[0]
    out_dir = output_dir or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{name}_{protein_type}.tsv")
    if os.path.isfile(out_path):
        return

    # Read clusters
    clusters = parse_mcl_clusters(mcl_file)
    n_ref_total = len(ref_ids)
    if n_ref_total == 0:
        print(f"[WARN] No reference sequences for {protein_type}")
        return

    # Compute cluster statistics
    rows = []
    for cluster in clusters:
        cs = len(cluster)
        rc = len(set(cluster) & ref_ids)
        rf = rc / n_ref_total
        label = "gray" if rf < pct_color_threshold else f">={int(rf * 100)}%"
        rows.append({
            'ClusterSize': cs,
            'ReferenceCount': rc,
            'ReferenceFraction': round(rf, 5),
            'ClassLabel': label
        })

    # Save output TSV
    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f"[SAVE] Saved MCL-reference data to {out_path}")


def parse_mcl_clusters(mcl_file_path):
    with open(mcl_file_path, 'r') as f:
        return [line.strip().split() for line in f if line.strip()]



def export_existing_mcc_curve(df: pd.DataFrame, score_col: str, mcc_col: str, output_path: str) -> None:
    """
    Exports the given MCC values for each unique score threshold to a TSV file.

    Parameters:
        df (pd.DataFrame): Must contain columns score_col and mcc_col.
        score_col (str): Column containing classifier scores (e.g. bitscore).
        mcc_col (str): Column containing precomputed MCC values.
        output_path (str): Path to the TSV output file.
    """
    if os.path.isfile(output_path):
        return

    # Sort and drop duplicate scores, keeping the first occurrence
    sorted_df = df[[score_col, mcc_col]].dropna().drop_duplicates(subset=score_col)
    sorted_df = sorted_df.sort_values(by=score_col, ascending=False)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    sorted_df.to_csv(output_path, sep='\t', index=False)
    print(f"[SAVE] MCC curve data exported to {output_path}")




    
        
def save_cutoffs_table(cutoffs_collection, output_dir, filename='cutoffs.tsv'):
    """
    Saves the cutoffs_collection into a tab-separated file.

    Parameters:
        cutoffs_collection (dict): Mapping name -> dict with keys
                                   'optimized cutoff', 'trusted cutoff', 'noise cutoff'
        output_dir (str): Directory in which to save the TSV.
        filename (str): Name of the TSV file (default: 'cutoffs.tsv').
    """
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)

    # Define column order
    fieldnames = ['name', 'optimized', 'trusted', 'noise']

    with open(out_path, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile,
                                fieldnames=fieldnames,
                                delimiter='\t')
        # write header
        writer.writeheader()

        # write one line per entry
        for name, cuts in cutoffs_collection.items():
            writer.writerow({
                'name':      name,
                'optimized': cuts.get('optimized cutoff'),
                'trusted':   cuts.get('trusted cutoff'),
                'noise':     cuts.get('noise cutoff'),
            })

    print(f"[SAVE] Cutoffs table saved to {out_path}")
        

def save_performance_table(performance_collection, output_dir, filename='performance.tsv'):
    """
    Saves the performance_collection into a tab-separated file.

    Parameters:
        performance_collection (dict): Mapping name -> dict with keys
                                       'MCC' and 'Matrix [TP,FP,FN,TN]'
        output_dir (str): Directory in which to save the TSV.
        filename (str): Name of the TSV file (default: 'performance.tsv').
    """
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)

    # Define column order
    fieldnames = ['name', 'MCC', 'TP', 'FP', 'FN', 'TN']

    with open(out_path, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile,
                                fieldnames=fieldnames,
                                delimiter='\t')
        # Header schreiben
        writer.writeheader()

        # jede Zeile pro Eintrag
        for name, perf in performance_collection.items():
            matrix = perf.get('Matrix [TP,FP,FN,TN]', [])
            # Entpacke Matrix oder fülle mit None, falls unvollständig
            TP, FP, FN, TN = (matrix + [None]*4)[:4]
            writer.writerow({
                'name': name,
                'MCC':  perf.get('MCC'),
                'TP':   TP,
                'FP':   FP,
                'FN':   FN,
                'TN':   TN,
            })

    print(f"[SAVE] Performance table saved to {out_path}")
        
        
        
        
        
        
