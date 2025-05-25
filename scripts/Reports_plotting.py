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
        report_dict = myUtil.load_cache(options, pkl_file)
        
        # Load reported cutoffs
        cutoffs = myUtil.load_cache(options, cutoffs_file)
        hmm_protein_name, optimized_cutoff, trusted_cutoff, noise_cutoff, best_MCC, best_matrix = cutoffs
        cutoffs_dict = {'trusted cutoff':trusted_cutoff, 'noise cutoff':noise_cutoff, 'optimized cutoff':optimized_cutoff}
        cutoff_collection[hmm_protein_name] = cutoffs_dict
        performance_collection[hmm_protein_name] = {'MCC': best_MCC, 'Matrix [TP,FP,FN,TN]': best_matrix}
        
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(report_dict, orient='index')
        df.index.name = 'hit_id'

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
        roc_path = os.path.join(hmm_report_output_dir, f"{protein}_roc.png")
        plot_roc_from_metrics(df, score_col='bitscore', trusted_cutoff=trusted_cutoff, noise_cutoff = noise_cutoff, optimized_cutoff = optimized_cutoff, output_path=roc_path)

        # Confusion stats per gene neighbourhood

        for name,cutoff in cutoffs_dict.items():
            plot_neighborhood_confusion_single(df, cutoff, name, output_dir = hmm_report_output_dir)

        # MCL selection plots
        protein_name = protein.split('_',1)[-1]
        mcl_file = mcl_results[protein_name]
        grp1_ref_seqs_set = grp1_refs[protein_name]
        grp1_cuts = mcl_grp1_cutoff[protein_name]
        
        plot_mcl_vs_references(mcl_file, "grp1_mcl_cluster_selection", grp1_ref_seqs_set, grp1_cuts['density_threshold'], grp1_cuts['reference_threshold'], output_dir = hmm_report_output_dir)
        print(f"[SAVE] MCL selection plot grp1")
        grp2_ref_seqs_set = grp2_refs[protein_name]
        grp2_cuts = mcl_grp2_cutoff[protein_name]
        
        plot_mcl_vs_references(mcl_file, "grp2_mcl_cluster_selection", grp2_ref_seqs_set, grp2_cuts['density_threshold'], grp2_cuts['reference_threshold'], output_dir = hmm_report_output_dir)
        print(f"[SAVE] MCL selection plot grp2")
        
    save_cutoffs_table(cutoff_collection, options.Hidden_markov_model_directory, 'cutoffs.txt')
    save_performance_table(performance_collection, options.Hidden_markov_model_directory, 'performance.txt')
    return



def plot_roc_from_metrics(df: pd.DataFrame,
                          score_col: str = 'bitscore',
                          trusted_cutoff: float = None,
                          noise_cutoff: float = None,
                          optimized_cutoff: float = None,
                          output_path: str = 'roc_curve.png') -> None:
    """
    Plots an ROC curve using precomputed confusion matrices in the DataFrame,
    and highlights points nearest above specified cutoffs with subtle markers.

    Parameters:
        df (pd.DataFrame): Must contain:
            - score_col: continuous classifier score
            - matrix_col: iterable of [TP, FP, FN, TN] per row
        trusted_cutoff (float): Score threshold for 'Trusted'
        noise_cutoff (float): Score threshold for 'Noise'
        optimized_cutoff (float): Score threshold for 'Optimized'
        output_path (str): Path to save the ROC plot
    """
    #
    if os.path.isfile(output_path):
        return
    
    # Extract scores and confusion matrix entries
    scores = df[score_col].values
    TP     = df['TP'].values
    FP     = df['FP'].values
    FN     = df['FN'].values
    TN     = df['TN'].values

    # Compute TPR and FPR
    with np.errstate(divide='ignore', invalid='ignore'):
        tprs = TP / (TP + FN)
        fprs = FP / (FP + TN)
    tprs = np.nan_to_num(tprs)
    fprs = np.nan_to_num(fprs)

    # Unique thresholds descending
    unique_scores = np.unique(scores)[::-1]
    # Map score -> corresponding average point
    roc_points = []
    for s in unique_scores:
        idxs = np.where(scores == s)[0]
        roc_points.append((fprs[idxs].mean(), tprs[idxs].mean(), s))
    fpr_vals, tpr_vals, thresh_vals = zip(*roc_points)

    # Plot base ROC in light gray
    plt.figure(figsize=(8, 6))
    plt.plot(fpr_vals, tpr_vals, color='lightgray', linewidth=2, label='ROC')
    plt.scatter(fpr_vals, tpr_vals,
                color='gray', edgecolor='black', s=40, alpha=0.8,
                label='_nolegend_')
    plt.plot([0, 1], [0, 1], linestyle='--', color='lightgray')

    # Helper to mark cutoff: smallest score >= cutoff
    def mark_cutoff(cutoff, label, color, marker):
        if cutoff is None:
            return
        # choose next higher threshold
        cands = [s for s in unique_scores if s >= cutoff]
        s_sel = min(cands) if cands else max(unique_scores)
        i = thresh_vals.index(s_sel)
        plt.scatter(fpr_vals[i], tpr_vals[i],
                    color=color, marker=marker,
                    edgecolor='black', linewidth=1,
                    s=80, alpha=1.0,
                    label=f"{label} ({s_sel})")

    # Subtle pastel colors for cutoffs
    mark_cutoff(trusted_cutoff,   'Trusted',   '#ADD8E6', 's')  # light blue square
    mark_cutoff(noise_cutoff,     'Noise',     '#FFB6C1', '^')  # light pink triangle
    mark_cutoff(optimized_cutoff, 'Optimized', '#90EE90', 'D')  # light green diamond

    # Labels and legend
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve with Highlighted Cutoffs')
    plt.legend(loc='lower right', frameon=True)
    plt.grid(True, color='lightgray', linestyle=':')
    plt.tight_layout()

    # Save plot
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"[SAVE] Saved ROC plot to {output_path}")





def plot_neighborhood_confusion_single(df: pd.DataFrame,
                                        cutoff: float,
                                        cutoff_name: str = 'cutoff',
                                        neighborhood_col: str = 'neighborhood',
                                        score_col: str = 'bitscore',
                                        true_label_col: str = 'true_value',
                                        top_n: int = None,
                                        output_dir: str = None) -> None:
    """
    Plots a stacked bar chart showing TP/FP/FN/TN counts per genomic neighborhood at a specific cutoff,
    grouping reverse complements together, sorting by cluster size, and coloring individual genes.

    Parameters:
        df (pd.DataFrame): Must contain score, true labels, and neighborhood columns.
        cutoff (float): Score threshold to classify predictions as positive (>=) or negative (<).
        cutoff_name (str): Descriptive name for the cutoff (used in titles/filenames).
        neighborhood_col (str): Column with neighborhood lists.
        score_col (str): Column with classifier scores.
        true_label_col (str): Column with true labels ('TP' or 'TN').
        top_n (int, optional): Number of top neighborhoods to include by total count. If None, include all.
        output_dir (str, optional): Directory to save the plot. Defaults to cwd if None.
    """
    # Define output and check if already exists
    out_dir = output_dir or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    fname = f'neighborhood_confusion_{cutoff_name}.png'
    path = os.path.join(out_dir, fname)
    if os.path.isfile(path):
        return

    # Validate inputs
    for col in [neighborhood_col, score_col, true_label_col]:
        if col not in df.columns:
            raise KeyError(f"DataFrame must contain column '{col}'")

    # Determine predicted classes
    preds_pos = df[score_col] >= cutoff
    true_pos = df[true_label_col] == 'TP'
    df['_conf_label'] = np.where(preds_pos & true_pos, 'TP',
                          np.where(preds_pos & ~true_pos, 'FP',
                          np.where(~preds_pos & true_pos, 'FN', 'TN')))

    # Canonicalize neighborhood: group reverse complements
    def canonical(nb):
        seq = tuple(tuple(x) for x in (nb or []))
        rev = tuple(reversed(seq))
        return seq if seq <= rev else rev

    df['_nb'] = df[neighborhood_col].map(canonical)

    # Count per neighborhood group and label
    counts = {nb: Counter(sub_df['_conf_label']) for nb, sub_df in df.groupby('_nb')}
    totals = {nb: sum(cnt.values()) for nb, cnt in counts.items()}

    # Select neighborhoods sorted by total descending or all
    selected = sorted(totals, key=totals.get, reverse=True)
    if top_n is not None:
        selected = selected[:top_n]

    # Build DataFrame for plotting
    plot_df = pd.DataFrame(
        [{ 'TP': counts[nb].get('TP', 0), 'FP': counts[nb].get('FP', 0),
           'FN': counts[nb].get('FN', 0), 'TN': counts[nb].get('TN', 0)}
         for nb in selected], index=selected
    )

    # Prepare clusters list of gene names
    clusters = [[gene for tup in nb for gene in tup] for nb in plot_df.index]
    # label_texts (space-separated gene names)
    label_texts = [" ".join(genes) for genes in clusters]

    # Plot stacked horizontal bars
    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df)*0.5)))
    fig.subplots_adjust(left=0.4, right=0.95)
    
    bottom = np.zeros(len(plot_df))
    cmap_conf = {'TN':'#d3d3d3','FP':'#ff9999','FN':'#ffcc99','TP':'#66b3ff'}
    for label in ['TN','FP','FN','TP']:
        vals = plot_df[label].values
        ax.barh(range(len(plot_df)), vals, left=bottom,
                color=cmap_conf[label], edgecolor='black', label=label)
        bottom += vals
        

    # Remove default y-tick labels
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(['']*len(plot_df))


    # Add individual gene names with background
    x_min, x_max = ax.get_xlim()
    x_span = x_max - x_min
    x_start = x_min - 2 * x_span
    for idx, genes in enumerate(clusters):
        x_offset = x_start
        for gene in genes:
            # consistent low-sat color per gene
            h = (int(hashlib.md5(gene.encode()).hexdigest(), 16) % 360) / 360.0
            rgb = mcolors.hsv_to_rgb((h, 0.3, 0.8))
            ax.text(x_offset, idx, gene,
                    va='center', ha='left',
                    color='black', backgroundcolor=rgb,
                    fontsize=9, clip_on=False)
            x_offset += len(gene) * 0.02 * x_span # Defines space between labels

    ax.invert_yaxis()
    ax.set_xlabel('Count')
    ax.set_title(f'Confusion matrix by genomic context at {cutoff_name}: hit score ≥ {cutoff}')
    ax.legend(loc='lower right')
    ax.grid(axis='x', linestyle='--', linewidth=0.5)

    # Save figure and close
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Clean up temporary columns
    df.drop(columns=['_conf_label','_nb'], inplace=True)
    
    
    
def plot_mcl_vs_references(mcl_file: str,
                            name: str,
                            ref_ids: set,
                            ref_line_pct: float = 0.1,
                            pct_color_threshold: float = 0.05,
                            max_pct: int = 10,
                            output_dir: str = None) -> None:
    """
    Plots cluster size vs. number of reference sequences for a single MCL clusters file.

    Parameters:
        mcl_file (str): Path to the MCL cluster file.
        ref_ids (set): Set of reference sequence IDs.
        ref_line_pct (float): Fractional threshold for the reference line (e.g., 0.1 for 10%).
        pct_color_threshold (float): Fractional threshold for coloring points; below this are gray.
        max_pct (int): Maximum percent value for the color scale (default 10).
        output_dir (str, optional): Directory in which to save the plot. Defaults to current working directory.
    """

    # Extract protein type from filename
    protein_type = os.path.splitext(os.path.basename(mcl_file))[0]

    # Output directory
    out_dir = output_dir or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{name}_{protein_type}.png")
    if os.path.isfile(out_path):
        return

    # Parse clusters
    clusters = parse_mcl_clusters(mcl_file)
    n_ref_total = len(ref_ids)
    if n_ref_total == 0:
        print(f"[WARN] No reference sequences for {protein_type}")
        return

    # Collect metrics
    cluster_sizes = []
    ref_counts = []
    ref_fractions = []
    for cluster in clusters:
        cs = len(cluster)
        rc = len(set(cluster) & ref_ids)
        cluster_sizes.append(cs)
        ref_counts.append(rc)
        ref_fractions.append(rc / n_ref_total)

    # Convert to arrays
    cs_arr = np.array(cluster_sizes)
    rc_arr = np.array(ref_counts)
    rf_arr = np.array(ref_fractions)

    # Mask for low-fraction points
    mask_low = rf_arr < pct_color_threshold
    pct_bins = np.clip((rf_arr[~mask_low] * 100).astype(int), 0, max_pct)
    cmap = plt.get_cmap("plasma", max_pct + 1)

    # Create plot
    plt.figure(figsize=(10, 6))
    # Gray for low
    plt.scatter(
        cs_arr[mask_low], rc_arr[mask_low],
        color='lightgray', s=80, edgecolor='k', alpha=0.75,
        label=f"<{int(pct_color_threshold*100)}% ref/all_ref cutoff"
    )
    # Colored for others
    scatter = plt.scatter(
        cs_arr[~mask_low], rc_arr[~mask_low],
        c=pct_bins, cmap=cmap, vmin=0, vmax=max_pct,
        s=80, edgecolor='k', alpha=0.75
    )

    # Colorbar
    cbar = plt.colorbar(scatter, ticks=range(0, max_pct + 1))
    cbar.set_label("Fraction of reference sequences from total number of reference sequences (%)")

    # Reference line
    x_vals = np.arange(1, max(cs_arr.max(), 200) + 1)
    plt.plot(
        x_vals, x_vals * ref_line_pct,
        linestyle='--', color='black', linewidth=1,
        label=f"{int(ref_line_pct*100)}% ref/total cutoff"
    )

    # Additional lines for orientation bei 25%, 50% und 75%
    for perc, shade in zip([0.25, 0.50, 0.75], ['darkgray', 'gray', 'lightgray']):
        plt.plot(
            x_vals, x_vals * perc,
            linestyle='--', color=shade, linewidth=1,
            label=f"{int(perc*100)}% ref/total"
        )

    # Labels, title, grid
    plt.xlabel("Cluster size (number of sequences)")
    plt.ylabel("Number of included reference sequences")
    plt.title(f"{protein_type}: MCL cluster selection")
    plt.grid(True)
    plt.legend(loc='upper right')

    # Save and close
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def parse_mcl_clusters(mcl_file_path):
    with open(mcl_file_path, 'r') as f:
        return [line.strip().split() for line in f if line.strip()]






    
        
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
        
        
        
        
        
        
