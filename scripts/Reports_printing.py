#!/usr/bin/python


import os
import csv
import glob
import hashlib
import pickle
from typing import Dict, List, Set, Any, Optional

import pandas as pd
import numpy as np

from collections import Counter, defaultdict

from . import myUtil
from . import Reports
from . import Csb_proteins
from . import Pam_Mx_algorithm

logger = myUtil.logger


def process_initial_validations(
    options: Any,
    report_directory: str,
    init_val_dir: str,
    db_path: str,
    output_dir: Optional[str] = None
) -> None:
    """
    Main entry: loads performance and cutoff PKLs, writes summary tables and all per-protein reports.

    Args:
        options: Main options/config object.
        report_directory (str): Path to main reports directory.
        init_val_dir (str): Path to directory with initial validation PKLs.
        db_path (str): Path to SQLite DB.
        output_dir (Optional[str]): Where to save summary/output. If None, 'Reports' subdir of report_directory.

    Returns:
        None. Writes output files.
    """
    
    # Prepare directories and file lists
    output_dir = _prepare_output_dir(report_directory, output_dir)
    pkl_files = _find_pkl_files(init_val_dir)

    # Collect cutoffs and performance metrics
    cutoff_collection = _collect_cutoff_and_performance(pkl_files, init_val_dir, options)

    # Save the performance and cutoffs in the report directory
    save_cutoffs_table(cutoff_collection['cutoffs'], output_dir, filename='all_cutoffs.txt')
    save_performance_table(cutoff_collection['performance'], output_dir, filename='all_performance.txt')
    
    write_pkl_tsv_reports(options, db_path, cutoff_collection, output_dir, pkl_files)
    
    # TODO fix this routine
    #write_pam_tsv_report(options, output_dir)

    return
    
    # TODO MCL file writer needed
    # Load reference and MCL data
    #mcl_grp1_cutoff, mcl_grp2_cutoff, mcl_results = _load_mcl_cutoff_sets(options)

    # Save presence-absence matrices of all reference sequences
    
    # Process each report
    #for pkl_file in pkl_files:
    #    _handle_pkl_file(
    #        options, pkl_file, init_val_dir, db_path, output_dir,
    #        cutoff_collection, all_grp1, all_grp2,
    #        grp1_refs, grp2_refs,
    #        mcl_grp1_cutoff, mcl_grp2_cutoff, mcl_results
    #    )


####################################
# --- Directory & File Helpers --- #
####################################


def _prepare_output_dir(report_directory: str, output_dir: Optional[str]) -> str:

    if output_dir is None:
        output_dir = os.path.join(report_directory, 'Reports')
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def _find_pkl_files(init_val_dir: str) -> List[str]:
    pattern = os.path.join(init_val_dir, '*.ini_performance_pkl')
    return glob.glob(pattern)


# --- Data Loading & Collection ---

def _collect_cutoff_and_performance(
    pkl_files: List[str],
    init_val_dir: str,
    options: Any
) -> Dict[str, Dict]:

    cutoffs = {}
    performance = {}
    for pkl in pkl_files:
        name, cuts, perf = _extract_cutoff_and_performance(pkl, init_val_dir, options)
        cutoffs[name] = cuts
        performance[name] = perf
    # Optionally write summary files here
    return {'cutoffs': cutoffs, 'performance': performance}


def _extract_cutoff_and_performance(
    pkl_file: str,
    init_val_dir: str,
    options: Any
) -> (str, Dict, Dict):
    protein = os.path.splitext(os.path.basename(pkl_file))[0]
    cutoffs_file = os.path.join(init_val_dir, f"{protein}.ini_cutoffs_pkl")
    data = myUtil.load_cache(options, f"{protein} cutoffs", file_path=cutoffs_file)
    name, opt, trust, noise, best_MCC, best_matrix = data
    cuts = {'optimized cutoff': opt,
            'trusted cutoff': trust,
            'noise cutoff': noise}
    perf = {'MCC': best_MCC,
            'Matrix [TP,FP,FN,TN]': best_matrix}
    return name, cuts, perf





#################################################
####
#################################################

def write_pkl_tsv_reports(
    options: Any,
    db_path: str,
    collections: Dict,
    output_dir: str,
    pkl_files: List[str]
) -> None:
    """
    For each per-protein PKL, generates all tabular/TSV reports.

    Args:
        options: Main config/options.
        db_path (str): SQLite DB path.
        collections (dict): {'cutoffs':..., 'performance':...}
        output_dir (str): Directory for output.
        pkl_files (list): List of performance PKL files.

    Returns:
        None.
    """
    
    for pkl_file in pkl_files:
        protein = os.path.splitext(os.path.basename(pkl_file))[0]
        report_dir = os.path.join(output_dir, protein)
        os.makedirs(report_dir, exist_ok=True)

        if _skip_if_all_outputs_exist(protein, report_dir, collections):
            return

        df = _load_and_validate_report(options, pkl_file)

        if df is None:
            return

        df = _add_neighborhood_info(df, protein, options, db_path)

        _write_enriched_table(df, report_dir, protein)
        _generate_roc_and_mcc(df, report_dir, protein, collections)
        _generate_neighborhood_confusion(df, report_dir, collections, protein)




    return


def _load_and_validate_report(options: Any, pkl_file: str) -> Optional[pd.DataFrame]:

    report = myUtil.load_cache(options, 'Report metrics', pkl_file)
    df = pd.DataFrame.from_dict(report, orient='index')
    df.index.name = 'hit_id'

    if df.empty:
        logger.error(f"Missing data for pkl cache file {os.path.basename(pkl_file)}")
        return None

    required = ['bitscore', 'MCC', 'TP', 'FP', 'FN', 'TN']
    missing = [c for c in required if c not in df.columns]
    if missing:
        logger.error(f"Missing columns {missing} in {os.path.basename(pkl_file)}")
        return None

    return df


def _add_neighborhood_info(
    df: pd.DataFrame,
    protein: str,
    options: Any,
    db_path: str
) -> pd.DataFrame:

    domain = protein.split('_', 1)[-1]
    cache_name = f"mcl_gene_vicinity_dict_{domain}.pkl"
    neighbors = myUtil.load_cache(options, cache_name) or {}

    missing = [hid for hid in df.index if hid not in neighbors]
    if missing:
        new_nb, _ = Csb_proteins.fetch_neighbouring_genes_with_domains(db_path, missing, options.sqlite_chunks)
        neighbors.update(new_nb)
        myUtil.save_cache(options, cache_name, neighbors)

    df['neighborhood'] = df.index.map(lambda hid: neighbors.get(hid, [['singleton']]))
    return df
    

def _skip_if_all_outputs_exist(
    protein: str,
    report_dir: str,
    collections: Dict
) -> bool:

    enriched_path = os.path.join(report_dir, f"{protein}_enriched.txt")
    roc_file = os.path.join(report_dir, f"{protein}_roc.txt")
    mcc_file = os.path.join(report_dir, f"{protein}_mcc.txt")
    cut_labels = collections['cutoffs'].get(protein, {}).keys()
    confusion_files = [
        os.path.join(report_dir, f"neighborhood_confusion_{label.replace(' ', '_')}.tsv")
        for label in cut_labels
    ]
    all_expected = [enriched_path, roc_file, mcc_file] + confusion_files
    if all(os.path.exists(path) for path in all_expected):
        logger.debug(f"All outputs exist for {protein} - skipping.")
        return True
    return False


def _plot_roc_and_mcc(
    df: pd.DataFrame,
    report_dir: str,
    protein: str,
    cutoffs: Dict
) -> None:
    trusted = cutoffs.get('trusted cutoff')
    noise = cutoffs.get('noise cutoff')
    opt = cutoffs.get('optimized cutoff')

    roc_file = os.path.join(report_dir, f"{protein}_roc.txt")
    write_roc_from_metrics_to_tsv(
        df, score_col='bitscore',
        trusted_cutoff=trusted,
        noise_cutoff=noise,
        optimized_cutoff=opt,
        output_path=roc_file
    )
    logger.debug(f"Saved ROC: {roc_file}")

    mcc_file = os.path.join(report_dir, f"{protein}_mcc.txt")
    export_existing_mcc_curve(
        df, score_col='bitscore', mcc_col='MCC', output_path=mcc_file
    )
    logger.debug(f"Saved MCC: {mcc_file}")


def _write_neighborhood_confusions(df, report_dir, cutoffs):
    for label, cutoff in cutoffs.items():
        write_neighborhood_confusion_data(
            df, cutoff, label,
            output_dir=report_dir
        )
         
def _write_enriched_table(df, report_dir, protein):
    enriched_path = os.path.join(report_dir, f"{protein}_enriched.txt")
    df.to_csv(enriched_path, sep='\t')

def _generate_roc_and_mcc(df, report_dir, protein, collections):
    cut_dict = collections['cutoffs'].get(protein, {})
    _plot_roc_and_mcc(df, report_dir, protein, cut_dict)

def _generate_neighborhood_confusion(df, report_dir, collections, protein):
    cut_dict = collections['cutoffs'].get(protein, {})
    _write_neighborhood_confusions(df, report_dir, cut_dict)





####################################
# Save the presence absence matrix #
####################################

def write_pam_tsv_report(options, output_dir):
    # Load reference and MCL data
    grp0_refs = myUtil.load_cache(options, 'basis_merged_grouped.pkl')
    grp1_refs = myUtil.load_cache(options, 'grp1_merged_grouped.pkl')
    grp2_refs = myUtil.load_cache(options, 'grp2_merged_grouped.pkl')
    grp3_refs = myUtil.load_cache(options, 'grp3_merged_grouped.pkl')

    all0 = set().union(*grp0_refs.values()) if grp0_refs else set()
    all1 = set().union(*grp1_refs.values()) if grp1_refs else set()
    all2 = set().union(*grp2_refs.values()) if grp2_refs else set()
    all3 = set().union(*grp3_refs.values()) if grp3_refs else set()
    
    # Save presence-absence matrices of all reference sequences
    save_PAM_data(
        options.database_directory,
        all0,
        os.path.join(output_dir, 'pam_basis_reference_seqs.txt'),
        options.cores
    )
    save_PAM_data(
        options.database_directory,
        all1,
        os.path.join(output_dir, 'pam_grp1_reference_seqs.txt'),
        options.cores
    )
    save_PAM_data(
        options.database_directory,
        all2,
        os.path.join(output_dir, 'pam_grp2_reference_seqs.txt'),
        options.cores
    )
    save_PAM_data(
        options.database_directory,
        all3,
        os.path.join(output_dir, 'pam_grp3_reference_seqs.txt'),
        options.cores
    )


def save_PAM_data(database_path, grouped, output_path, cores):
    if os.path.isfile(output_path):
        return
    pam = Pam_Mx_algorithm.create_presence_absence_matrix(
            grouped,
            database_directory=database_path,
            output="pam_with_all_domains",
            chunk_size=900,
            cores=cores
        )        
    Pam_Mx_algorithm.write_pam_to_tsv(pam, output_path)
    return
















def _load_mcl_cutoff_sets(options):
    m1 = myUtil.load_cache(options, 'mcl_grp2_cluster_selection_cutoffs.pkl')
    m2 = myUtil.load_cache(options, 'mcl_grp3_cluster_selection_cutoffs.pkl')
    results = myUtil.load_cache(options, 'mcl_clustering_results.pkl')
    return m1, m2, results


# --- Per-File Processing ---
"""
		def _handle_pkl_file(
		    options, pkl_file, init_val_dir, db_path, output_dir,
		    collections, all_grp1, all_grp2,
		    grp1_refs, grp2_refs,
		    mcl1_cutoff, mcl2_cutoff, mcl_results
		):
		    protein = os.path.splitext(os.path.basename(pkl_file))[0]
		    report_dir = os.path.join(output_dir, protein)
		    os.makedirs(report_dir, exist_ok=True)

		    # Define expected output files
		    enriched_path = os.path.join(report_dir, f"{protein}_enriched.txt")
		    roc_file = os.path.join(report_dir, f"{protein}_roc.txt")
		    mcc_file = os.path.join(report_dir, f"{protein}_mcc.txt")
		    # neighborhood confusion files for each cutoff label
		    cut_labels = collections['cutoffs'].get(protein, {}).keys()
		    confusion_files = [
			os.path.join(report_dir, f"neighborhood_confusion_{label.replace(' ', '_')}.tsv")
			for label in cut_labels
		    ]


		    # Skip if all core and confusion outputs exist
		    all_expected = [enriched_path, roc_file, mcc_file] + confusion_files
		    if all(os.path.exists(path) for path in all_expected):
			print(f"[SKIP] All outputs exist for {protein}, skipping.")
			return
			

		    print(f"\n[INFO] Generating HMM report for {protein}")

		    # Load and validate basic metrics
		    df = _load_and_validate_report(options, pkl_file)
		    if df is None:
			return

		    # Enrich with neighborhood and group membership
		    df = _add_neighborhood_info(df, protein, options, db_path)
		    df['in_grp1'] = df.index.isin(all_grp1)
		    df['in_grp2'] = df.index.isin(all_grp2)

		    # Save enriched table
		    _save_enriched(df, report_dir, protein)

		    # Plot ROC and MCC curves
		    cut_dict = collections['cutoffs'].get(protein, {})
		    _plot_roc_and_mcc(df, report_dir, protein, cut_dict)

		    # Neighborhood confusion stats
		    _write_neighborhood_confusions(df, report_dir, cut_dict)

		    # MCL selection plots
		    _write_mcl_plots_if_available(
			protein, report_dir,
			grp1_refs, grp2_refs,
			mcl1_cutoff, mcl2_cutoff, mcl_results
		    )
"""

# --- Metric Loading & Validation ---




# --- Neighborhood Enrichment ---



# --- Saving & Plotting ---





def _write_mcl_plots_if_available(
    protein, report_dir,
    grp1_refs, grp2_refs,
    mcl1_cutoff, mcl2_cutoff, mcl_results
):
    key = protein.split('_', 1)[-1]
    if key in mcl_results:
        if key in grp1_refs and key in mcl1_cutoff:
            write_mcl_vs_references(
                mcl_results[key],
                'grp1_mcl_cluster_selection',
                grp1_refs[key],
                mcl1_cutoff[key]['density_threshold'],
                mcl1_cutoff[key]['reference_threshold'],
                output_dir=report_dir
            )
            logger.debug(f"Saved MCL grp1: {key}")
        if key in grp2_refs and key in mcl2_cutoff:
            write_mcl_vs_references(
                mcl_results[key],
                'grp2_mcl_cluster_selection',
                grp2_refs[key],
                mcl2_cutoff[key]['density_threshold'],
                mcl2_cutoff[key]['reference_threshold'],
                output_dir=report_dir
            )
            logger.debug(f"Saved MCL grp2: {key}")
    else:
        logger.warning(f"No MCL data for {key}, skipping plots.")




################################################################################################
#############  Subroutine for the cutoff and performance reports and mcl loading ###############
################################################################################################















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
    logger.debug(f"Saved ROC data to {output_path}")



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
    fname.replace(' ','_')
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
    logger.debug(f"Saved neighborhood confusion data to {path}")

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
        logger.warning(f"No reference sequences for {protein_type} in mcl vs reference printing")
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
    logger.debug(f"Saved MCL-reference data to {out_path}")


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
    logger.debug(f"Saved MCC curve data exported to {output_path}")



#########################################################################
### Save cutoff tables and performance tables in the report directory ###
#########################################################################
        
def save_cutoffs_table(
    cutoffs_collection: Dict[str, Dict[str, float]],
    output_dir: str,
    filename: str = 'cutoffs.tsv'
) -> None:
    """
    Saves the cutoffs_collection into a tab-separated file.

    Parameters:
        cutoffs_collection (dict): Mapping name -> dict with keys
                                   'optimized cutoff', 'trusted cutoff', 'noise cutoff'
        output_dir (str): Directory in which to save the TSV.
        filename (str): Name of the TSV file (default: 'cutoffs.tsv').
    Saves the cutoffs_collection into a tab-separated file.

    Args:
        cutoffs_collection (dict): name -> dict with keys 'optimized cutoff', 'trusted cutoff', 'noise cutoff'
        output_dir (str): Output directory.
        filename (str): Name of the TSV file.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)

    # Define column order
    fieldnames = ['name', 'optimized', 'trusted', 'noise']
    sorted_items = sorted(cutoffs_collection.items(), key=lambda kv: kv[0])
    with open(out_path, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile,
                                fieldnames=fieldnames,
                                delimiter='\t')
        # write header
        writer.writeheader()

        # write one line per entry
        for name, cuts in sorted_items:
            writer.writerow({
                'name':      name,
                'optimized': cuts.get('optimized cutoff'),
                'trusted':   cuts.get('trusted cutoff'),
                'noise':     cuts.get('noise cutoff'),
            })

    logger.info(f"Saved HMM score cutoffs table to {out_path}")
        

def save_performance_table(
    performance_collection: Dict[str, Dict[str, Any]],
    output_dir: str,
    filename: str = 'performance.tsv'
) -> None:
    """
    Saves the performance_collection into a tab-separated file, sorted by 'name'.

    Args:
        performance_collection (dict): name -> dict with keys 'MCC', 'Matrix [TP,FP,FN,TN]'
        output_dir (str): Directory for TSV file.
        filename (str): File name.

    Returns:
        None.
    """
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)

    # Definiere Spaltenreihenfolge
    fieldnames = ['name', 'MCC', 'TP', 'FP', 'FN', 'TN']

    # Sortiere alle Einträge alphabetisch nach dem Schlüssel (name)
    sorted_items = sorted(performance_collection.items(), key=lambda kv: kv[0])

    with open(out_path, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile,
                                fieldnames=fieldnames,
                                delimiter='\t')
        # Header schreiben
        writer.writeheader()

        # Schreibe jede Zeile in alphabetischer Reihenfolge (Name = erste Spalte)
        for name, perf in sorted_items:
            matrix = perf.get('Matrix [TP,FP,FN,TN]', [])
            TP, FP, FN, TN = (matrix + [None]*4)[:4]
            writer.writerow({
                'name': name,
                'MCC':  perf.get('MCC'),
                'TP':   TP,
                'FP':   FP,
                'FN':   FN,
                'TN':   TN,
            })

    logger.info(f"Saved HMM validation performance table to {out_path}")

        
        
        
        

