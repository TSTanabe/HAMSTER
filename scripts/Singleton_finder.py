#!/usr/bin/python

import copy
import os
import subprocess
import sqlite3
import traceback
import pandas as pd
from collections import defaultdict

from . import Pam_Mx_algorithm
from . import Csb_proteins
from . import myUtil

    

#### Main routine of this module
def singleton_reference_finder(options, grouped):
    
    singleton_reference_seqs_dict = myUtil.load_cache(options, 'sng_training_proteinIDs.pkl')
    domain_score_limits = myUtil.load_cache(options, 'sng_training_proteinIDs_limits.pkl')
    
    if singleton_reference_seqs_dict and domain_score_limit:
        return domain_score_limits, singleton_reference_seqs_dict

    domain_score_limits, singleton_reference_seqs_dict = get_singleton_reference_sequences(options)
    
    predictions_all = predict_singleton_reference_seqs_for_each_domain(options.database_directory, grouped, singleton_reference_seqs_dict, options.cores, chunk_size=900)
    print(predictions_all)
    limits_dict, sng_reference_seq_dict = collect_predicted_singleton_hits_from_db(predictions_all, options.database_directory, plausible_cutoff = 0.02)
    
    sng_reference_seq_dict = {f"grp0_{k}": v for k, v in sng_reference_seq_dict.items()}
    print(sng_reference_seq_dict)
    myUtil.save_cache(options, 'sng_training_proteinIDs.pkl', sng_reference_seq_dict)
    myUtil.save_cache(options, 'sng_training_proteinIDs_limits.pkl', limits_dict)

    # Step 3: If cache is missing, recompute and export FASTAs
    if sng_reference_seq_dict:
        Csb_proteins.fetch_seqs_to_fasta_parallel(
            options.database_directory,
            sng_reference_seq_dict,
            options.fasta_output_directory,
            min_seq=10,
            max_seq=options.max_seqs,
            cores=options.cores
        )
    

    return limits_dict, sng_reference_seq_dict

def get_singleton_reference_sequences(options):


    # Load cache if available
    domain_score_limits = myUtil.load_cache(options, "sng_domain_score_limits.pkl")
    singleton_reference_seqs_dict = myUtil.load_cache(options, "sng_reference_seqs_dict.pkl")
    
    if domain_score_limits and singleton_reference_seqs_dict:
        print("Loading existing reference sequences for genes without conserved genomic context")
        return domain_score_limits, singleton_reference_seqs_dict

    # Get domains present in query but missing from training sets
    query_names = extract_protein_ids_from_fasta(options.self_query)
    training_set_domains = extract_domain_names_from_directory(options.fasta_output_directory)
    singletons = query_names - training_set_domains
    print(f"Proteins without recognized genomic context {singletons}")

    # Step 3: Query database instead of BLAST file
    singleton_hits = defaultdict(set)
    domain_score_limits = {}


    query = """
        SELECT domain, proteinID, score
        FROM Domains
        WHERE domain = ? AND blast_score_ratio > ?
    """

    with sqlite3.connect(options.database_directory) as con:
        cur = con.cursor()
        for singleton in singletons:
            cur.execute(query, (singleton, 0.9)) # HERE IS THE HARDCODE FOR BLAST SCORE RATIO CUTOFF
            hits = cur.fetchall()

            bitscores = []
            for domain, proteinID, score in hits:
                singleton_hits[domain].add(proteinID)
                bitscores.append(score)

            if bitscores:
                domain_score_limits[singleton] = {
                    "lower_limit": min(bitscores),
                    "upper_limit": max(bitscores)
                }

    # Prefix keys and fetch sequences to fasta
    singleton_reference_seqs_dict = {f"sng0_{k}": v for k, v in singleton_hits.items()}
    Csb_proteins.fetch_seqs_to_fasta_parallel(
        options.database_directory,
        singleton_reference_seqs_dict,
        options.fasta_output_directory,
        2, options.max_seqs, options.cores
    )

    myUtil.save_cache(options, "sng_reference_seqs_dict.pkl", singleton_reference_seqs_dict)
    myUtil.save_cache(options, "sng_domain_score_limits.pkl", domain_score_limits)

    return domain_score_limits, singleton_reference_seqs_dict


#######################

def predict_singleton_reference_seqs_for_each_domain(database_path, grouped, singleton, cores, chunk_size=900):
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

    for sng_domain, sng_proteins in singleton.items():
        print(f"\n[INFO] Processing singleton domain: {sng_domain}")

        # 1. Kombiniertes grouped: Basis + aktuelles Singleton
        grouped_plus = copy.deepcopy(grouped)
        grouped_plus[sng_domain] = sng_proteins
        
        
        # 2. Presence/Absence-Matrix erzeugen
        pam = Pam_Mx_algorithm.create_presence_absence_matrix(
            grouped_plus,
            database_directory=database_path,
            output=f"pam_with_{sng_domain}",
            chunk_size=chunk_size,
            cores=cores
        )
        
        # 3. PAM-Filter: nur Genomes behalten, die Singleton enthalten
        genomeIDs = [genomeID for genomeID, domains in pam.items() if sng_domain in domains]
        filtered_pam = {genomeID: pam[genomeID] for genomeID in genomeIDs if genomeID in pam}
        filtered_pam = pam
        if not filtered_pam:
            print(f"[WARN] No genomes with presence of {sng_domain} found – skipping.")
            continue

        # 4. Hit-Scores für alle ProteinIDs laden
        bsr_hit_scores = Pam_Mx_algorithm.fetch_bsr_scores(database_path, grouped_plus, chunk_size=chunk_size)

        # 5. Modell trainieren auf gefilterter PAM
        models, _ = Pam_Mx_algorithm.train_logistic_from_pam_with_scores(filtered_pam, bsr_hit_scores, cores=cores, target_domains={sng_domain})

        if sng_domain not in models:
            print(models)
            print(f"[WARN] No model trained for {sng_domain} – skipping.")
            continue

        # 6. Testdaten: PAM nur für grouped (ohne das aktuelle Singleton)
        base_pam = Pam_Mx_algorithm.create_presence_absence_matrix(
            grouped,
            database_directory=database_path,
            output=f"pam_base_for_{sng_domain}",
            chunk_size=chunk_size,
            cores=cores
        )

        # 7. In DataFrame konvertieren für Vorhersage
        genomes = sorted(base_pam.keys())
        domains = sorted({d for m in base_pam.values() for d in m})
        df = pd.DataFrame(index=genomes, columns=domains)

        for genome_id in genomes:
            for domain in domains:
                df.at[genome_id, domain] = ",".join(base_pam[genome_id].get(domain, [])) if domain in base_pam[genome_id] else ""
        
        # 8. Feature-Matrix bauen über Hilfsfunktion
        X_test = Pam_Mx_algorithm.build_presence_score_matrix(df, bsr_hit_scores)
        
        # 9. Nur Modellrelevante Features
        model = models[sng_domain]
        X_test = X_test.reindex(columns=model.feature_names_in_, fill_value=0)

        # 10. Vorhersagen berechnen
        prediction_scores = pd.Series(model.predict_proba(X_test)[:, 1], index=df.index)
        predictions_all[sng_domain] = prediction_scores

    return predictions_all




def collect_predicted_singleton_hits_from_db(predictions_all, database_path, plausible_cutoff=0.6, bsr_threshold=0.5):

    singleton_hits = defaultdict(set)
    singleton_score_limits = {}

    # Schritt 1: Bereite die Struktur für Domains → plausible genomeIDs
    domain_to_genomes = defaultdict(list)
    for singleton, prediction_series in predictions_all.items():
        domain = singleton.replace('sng0_', '')
        plausible_genomes = prediction_series[prediction_series > plausible_cutoff].index
        domain_to_genomes[domain].extend(plausible_genomes)

    if not domain_to_genomes:
        return {}, {}

    # Schritt 2: Iteriere über Domains, hole Treffer genomeweise
    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        for domain, genome_list in domain_to_genomes.items():
            best_hits = {}
            for i in range(0, len(genome_list), 900):  # Max 999 placeholders in SQLite
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
                rows = cur.fetchall()

                for genome_id, proteinID, score, bsr in rows:
                    if genome_id not in best_hits or bsr > best_hits[genome_id][0]:
                        best_hits[genome_id] = (bsr, proteinID, score)

            # Speichere beste Treffer pro Domain
            bitscores = []
            for genome_id, (bsr, proteinID, score) in best_hits.items():
                singleton_hits[domain].add(proteinID)
                bitscores.append(score)

            if bitscores:
                singleton_score_limits[domain] = {
                    "lower_limit": min(bitscores),
                    "upper_limit": max(bitscores)
                }

    return singleton_score_limits, singleton_hits


################ Routines for first main routine singleton _reference_sequences ##########################
    
def extract_protein_ids_from_fasta(fasta_path):
    protein_ids = set()
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                full_id = line[1:].split()[0]
                parts = full_id.split('___')
                protein_id = parts[1] if len(parts) > 1 else full_id
                protein_ids.add(protein_id)
    return protein_ids
    
def extract_domain_names_from_directory(directory_path):
    """
    Extrahiert eindeutige Domänennamen aus Dateinamen in einem Verzeichnis.
    Entfernt Dateiendung und alles bis einschließlich dem ersten '_'.
    """
    domain_names = set()
    
    for filename in os.listdir(directory_path):
        if os.path.isfile(os.path.join(directory_path, filename)):
            name_no_ext = os.path.splitext(filename)[0]  # z.B. "Query_ABC123"
            parts = name_no_ext.split('_', 1)  # nur am ersten "_" splitten
            domain = parts[1] if len(parts) > 1 else name_no_ext
            domain_names.add(domain)
    
    return domain_names
    
def get_min_bitscore_for_query(report_path, query_id, blast_score_ratio = 0.9):
    """
    Liest einen DIAMOND BLAST Report ein und gibt den kleinsten Bitscore für eine gegebene Query-ID zurück.
    
    Args:
        report_path (str): Pfad zur Diamond-Tabelle (.tab).
        query_id (str): Die qseqid (Query-ID), nach der gesucht werden soll.
    
    Returns:
        float: Der kleinste Bitscore für die Query, oder None falls nicht gefunden.
    """
    bitscores = []
    min_cutoff = 100
    max_cutoff = 100
    with open(report_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            hit_id, qid, evalue, bitscore = parts[:4]
            hit_id = hit_id.split('___')[1]
            if hit_id == query_id and qid == query_id:
                try:
                    bitscores.append(float(bitscore))
                except ValueError:
                    continue  # ungültiger bitscore

    if len(bitscores) == 1:
        max_cutoff = bitscores[0]*(2-blast_score_ratio)
        min_cutoff = bitscores[0]*blast_score_ratio #if only a single sequence is given than 0.9 blast score ratio is accepted
    else:
        max_cutoff = max(bitscores)
        min_cutoff = min(bitscores)
    return min_cutoff, max_cutoff










    
