#!/usr/bin/python
import os
import sqlite3
import csv
import warnings

from . import Csb_proteins
from . import myUtil
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool


from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning

# Global variables for the databse read only mode
_con = None
_cur = None

# Routines for the logistic regression based on presence absence matrices

def train_logistic_from_pam_with_scores(pam, hit_scores, cores=4, target_domains=None):
    """
    Trains logistic regression models using both presence/absence and hit scores.
    Can be restricted to train only for selected target domains.

    Args:
        pam (dict): {genomeID: {domain_label: [proteinIDs]}}
        hit_scores (dict): {proteinID: score}
        cores (int): Number of parallel jobs for sklearn
        target_domains (list or set, optional): Domain labels to train models for. If None, train for all.

    Returns:
        models (dict): {domain_label: trained LogisticRegression model}
        coef_df (pd.DataFrame): Coefficient table (domain_label x features)
    """
    warnings.filterwarnings("ignore", category=ConvergenceWarning)

    genomes = sorted(pam.keys())
    domains = sorted({domain for dommap in pam.values() for domain in dommap})
    
    # Build binary matrix and score matrix
    binary_data = []
    score_data = []

    for genome in genomes:
        binary_row = {}
        score_row = {}
        for domain in domains:
            proteins = pam[genome].get(domain, [])
            binary_row[domain + "_presence"] = int(bool(proteins))

            # Compute average hit score of proteins, or 0 if none valid
            valid_scores = [hit_scores[p] for p in proteins if p in hit_scores]
            score_row[domain + "_score"] = sum(valid_scores) / len(valid_scores) if valid_scores else 0.0

        binary_data.append(binary_row)
        score_data.append(score_row)

    binary_matrix = pd.DataFrame(binary_data, index=genomes)
    score_matrix = pd.DataFrame(score_data, index=genomes)

    # Combine matrices
    combined_matrix = pd.concat([binary_matrix, score_matrix], axis=1)

    # Train logistic regression models
    models = {}
    coefficients_table = {}
    
    for domain in domains:
        if target_domains and domain not in target_domains:
            continue

        y = binary_matrix[domain + "_presence"]
        

        if y.nunique() <= 1:
            print(f"[SKIP] {domain} class distribution = {y.value_counts().to_dict()} not enough class variety")
            continue  # skip domains with only 1 class

        # Exclude self-domain features to prevent leakage
        feature_cols = [col for col in combined_matrix.columns if not col.startswith(domain + "_")]
        X = combined_matrix[feature_cols]

        model = LogisticRegression(max_iter=1000)
        model.fit(X, y)

        models[domain] = model
        coefs = pd.Series(model.coef_[0], index=X.columns)
        coefs["(Intercept)"] = model.intercept_[0]
        coefficients_table[domain] = coefs

    coef_df = pd.DataFrame(coefficients_table).T
    return models, coef_df



def fetch_bsr_scores(database_path, grouped, chunk_size=900):
    """
    Fetches blast score ratios (BSR) for all proteinIDs in grouped dict from the SQLite database.

    Args:
        database_path (str): Path to the SQLite database.
        grouped (dict): Dictionary of domain → set(proteinIDs).
        chunk_size (int): Number of proteinIDs per query batch.

    Returns:
        dict: {proteinID: BSR score} — Missing BSR values will not be included.
    """
    all_protein_ids = list({pid for pids in grouped.values() for pid in pids})
    bsr_scores = {}

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        for chunk in chunked(all_protein_ids, chunk_size):
            placeholders = ",".join(["?"] * len(chunk))
            query = f"""
                SELECT proteinID, blast_score_ratio 
                FROM Domains 
                WHERE proteinID IN ({placeholders}) AND blast_score_ratio IS NOT NULL
            """
            cur.execute(query, chunk)
            for proteinID, bsr in cur.fetchall():
                bsr_scores[proteinID] = bsr

    return bsr_scores



def build_presence_score_matrix(df, hit_scores):
    """
    Erstellt binäre Präsenz- und Score-Matrix aus einem Domain-DataFrame.

    Args:
        df (pd.DataFrame): DataFrame mit genomeID als Index, Domains als Spalten,
                           Zellinhalt = Komma-separierte Protein-IDs oder leer
        hit_scores (dict): Mapping {proteinID: float}

    Returns:
        pd.DataFrame: Kombination aus binärer Präsenzmatrix und Scorematrix
    """
    binary_matrix = pd.DataFrame(index=df.index)
    score_matrix = pd.DataFrame(index=df.index)

    def score_func(protein_str):
        if not protein_str or not isinstance(protein_str, str):
            return 0.0
        pids = protein_str.split(",")
        return max([hit_scores.get(pid, 0.0) for pid in pids if pid])

    for domain in df.columns:
        binary_matrix[domain + "_presence"] = df[domain].apply(lambda x: 1 if isinstance(x, str) and x != "" else 0)
        score_matrix[domain + "_score"] = df[domain].apply(score_func)

    combined_matrix = pd.concat([binary_matrix, score_matrix], axis=1)
    return combined_matrix







#########################################################################################
######################### Create presence absence matrix ################################
#########################################################################################


def create_presence_absence_matrix(proteinID_sets, database_directory, output, chunk_size=900, cores=8):
    """
    Creates a TSV file where each cell contains the proteinID(s) for a domain in a genome,
    using a grouped dictionary instead of reading FASTA files.

    Args:
        grouped (dict): {domain_label (e.g. "grp0_abc"): set(proteinIDs)}
        database_directory (str): Path to SQLite database.
        output_path (str): Path to the output TSV file.
        chunk_size (int): Max number of proteinIDs per DB query.
        cores (int): Number of parallel worker processes.
    """

    all_protein_ids = {pid for pids in proteinID_sets.values() for pid in pids}

    if not all_protein_ids:
        return

    # Step 1: Map proteinID → genomeID
    chunk_args = [
        (database_directory, chunk)
        for chunk in chunked(all_protein_ids, chunk_size)
    ]

    with Pool(processes=cores) as pool:
        results = pool.map(_fetch_protein_to_genome_chunk, chunk_args)

        protein_to_genome_dict = {}
        for partial in results:
            protein_to_genome_dict.update(partial)

        # Step 2: genomeID → domain → [proteinIDs]
        args_list = [
            (domain, protein_ids, protein_to_genome_dict)
            for domain, protein_ids in proteinID_sets.items()
        ]
 
    
        genome_domain_matrix = defaultdict(lambda: defaultdict(list))

        results = pool.map(_map_domain_proteins_to_genomes, args_list)

    for partial in results:
        for genome_id, domain_prot_pairs in partial.items():
            for domain, pid in domain_prot_pairs:
                genome_domain_matrix[genome_id][domain].append(pid)
    
    return genome_domain_matrix

def write_pam_to_tsv(pam, output_path):
    """
    Write the presence/absence matrix to a TSV file.

    Args:
        pam (dict): {genomeID: {domain_label: [proteinIDs]}}
        output_path (str): Path to TSV output file
    """
    # Automatisch alle Domänenlabels sammeln
    domain_set = set()
    for domain_map in pam.values():
        domain_set.update(domain_map.keys())
    domain_labels = sorted(domain_set)

    # Schreiben
    with open(output_path, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["genomeID"] + domain_labels)

        for genome_id in sorted(pam.keys()):
            row = [genome_id]
            for domain in domain_labels:
                prot_list = pam[genome_id].get(domain, [])
                row.append(",".join(prot_list) if prot_list else "")
            writer.writerow(row)
    
    
def _map_domain_proteins_to_genomes(args):
    domain, protein_ids, protein_to_genome = args
    result = defaultdict(list)
    for pid in protein_ids:
        genome_id = protein_to_genome.get(pid)
        if genome_id:
            result[genome_id].append((domain, pid))
            
    return result

def chunked(iterable, size):
    """Yield successive chunk-sized lists from iterable."""
    for i in range(0, len(iterable), size):
        yield list(iterable)[i:i + size]
        
#############################################################################
######## Read only database for the PAM mapping #############################
#############################################################################

def _get_readonly_cursor(database_path: str):
    """
    Öffnet eine nur-lesende SQLite-Verbindung (URI mode=ro) mit query_only,
    und cached Verbindung+Cursor pro Prozess.
    """
    global _con, _cur
    if _con is None:
        _con = sqlite3.connect(
            f"file:{database_path}?mode=ro",
            uri=True,
            check_same_thread=False
        )
        _con.execute("PRAGMA query_only = TRUE;")
        _cur = _con.cursor()
    return _cur

# -------------------------------------------------------------------
# 2) Chunk-weises Protein→Genome Mapping
# -------------------------------------------------------------------
def _fetch_protein_to_genome_chunk(args):
    """
    Args:
        args: Tuple(database_path, chunk: List[str])
    Returns:
        Dict[proteinID, genomeID]
    """
    database_path, chunk = args
    if not chunk:
        return {}
    cur = _get_readonly_cursor(database_path)
    placeholders = ",".join("?" for _ in chunk)
    sql = (
        f"SELECT proteinID, genomeID "
        f"FROM Proteins "
        f"WHERE proteinID IN ({placeholders})"
    )
    cur.execute(sql, chunk)
    rows = cur.fetchall()
    return {proteinID: genomeID for proteinID, genomeID in rows}



#########################################################################################
################## Train logistic model on presence absence matrix file #################
#########################################################################################


def process_domain_hits_pool(args):
    domain, genome_list, database_path, bsr_threshold = args
    best_hits = {}
    con = sqlite3.connect(database_path)
    cur = con.cursor()

    for i in range(0, len(genome_list), 999):
        chunk = genome_list[i:i + 999]
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

    return domain, proteinIDs


def collect_predicted_singleton_hits_from_db_parallel(predictions_all, database_path, plausible_cutoff=0.6, bsr_cutoff=0.5, max_parallel=4):
    domain_to_genomes = defaultdict(list)
    for singleton, prediction_series in predictions_all.items():
        domain = singleton.split('_', 1)[-1]
        plausible_genomes = prediction_series[prediction_series > plausible_cutoff].index
        domain_to_genomes[domain].extend(plausible_genomes)

    task_args = [(domain, genome_list, database_path, bsr_cutoff) for domain, genome_list in domain_to_genomes.items()]

    singleton_hits = {}

    with Pool(processes=max_parallel) as pool:
        results = pool.map(process_domain_hits_pool, task_args)

    for domain, proteinIDs in results:
        if proteinIDs:
            singleton_hits[domain] = proteinIDs

    return singleton_hits


















































