#!/usr/bin/python
import os
import sqlite3
import csv
import warnings

from . import Csb_proteins
from . import Csb_phylogeny

import pandas as pd
from collections import defaultdict
from multiprocessing import Pool


from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning


def create_presence_absence_matrix(faa_dir, database_directory, output_path, chunk_size=900, extensions=(".faa",), prefixes=("grp0", "sng0"), cores=4):
    """
    Creates a TSV file where each cell contains the proteinID(s) for a domain in a genome.
    Multiple proteinIDs per genome-domain pair are joined with commas. Empty if absent.

    Args:
        faa_dir (str): Directory with input FASTA files.
        database_directory (str): Path to SQLite database.
        output_path (str): Path to the output TSV file.
        chunk_size (int): Max number of proteinIDs per DB query.
        extensions (tuple): File extensions to include (e.g. (".faa", ".fa")).
        prefixes (tuple): Acceptable filename prefixes (e.g. ("grp0", "sng0")).
    """
    
    if os.path.isfile(output_path):
        return
    
    # Step 1: Find relevant FASTA files
    faa_files = [
        f for f in os.listdir(faa_dir)
        if f.endswith(tuple(extensions)) and f.startswith(tuple(prefixes))
    ]
    domain_to_proteins = {}
    all_protein_ids = set()
    
    if not faa_files:
        # There were no faa files found
        return
    
    for faa_file in faa_files:
        print(f"Fetching presence from faa_file {faa_file}")
        domain = os.path.splitext(faa_file)[0]#.replace("grp0_", "").replace("sng0_", "")
        file_path = os.path.join(faa_dir, faa_file)
        protein_ids = set()

        with open(file_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    prot_id = line[1:].strip().split()[0]
                    protein_ids.add(prot_id)
                    all_protein_ids.add(prot_id)

        domain_to_proteins[domain] = protein_ids


    ### Step 2: Map proteinID to genomeID
    # Erzeuge Chunk-Argumente
    chunk_args = [
        (database_directory, chunk)
        for chunk in chunked(all_protein_ids, chunk_size)
    ]

    # Parallel ausführen
    with Pool(processes=4) as pool:
        results = pool.map(_fetch_protein_to_genome_chunk, chunk_args)

    # Ergebnisse zusammenführen
    protein_to_genome = {}
    for partial in results:
        protein_to_genome.update(partial)




    ### Step 3: Build genomeID to domain to proteinID(s)
    # Prepare arguments per domain
    args_list = [
        (domain, protein_ids, protein_to_genome)
        for domain, protein_ids in domain_to_proteins.items()
    ]

    genome_domain_matrix = defaultdict(lambda: defaultdict(list))

    with Pool(processes=4) as pool:
        results = pool.map(_map_domain_proteins_to_genomes, args_list)

    # Merge results
    for partial in results:
        for genome_id, domain_prot_pairs in partial.items():
            for domain, pid in domain_prot_pairs:
                genome_domain_matrix[genome_id][domain].append(pid)

    ### Step 4: Write TSV
    all_domains = sorted(domain_to_proteins.keys())
    with open(output_path, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["genomeID"] + all_domains)

        for genome_id in sorted(genome_domain_matrix.keys()):
            row = [genome_id]
            for domain in all_domains:
                prot_list = genome_domain_matrix[genome_id].get(domain, [])
                row.append(",".join(prot_list) if prot_list else "")
            writer.writerow(row)
    return
    
    
def _map_domain_proteins_to_genomes(args):
    domain, protein_ids, protein_to_genome = args
    result = defaultdict(list)
    for pid in protein_ids:
        genome_id = protein_to_genome.get(pid)
        if genome_id:
            result[genome_id].append((domain, pid))
    return result

def _fetch_protein_to_genome_chunk(args):
    database_path, chunk = args
    local_result = {}
    with sqlite3.connect(database_path) as con:
        cur = con.cursor()
        placeholders = ",".join(["?"] * len(chunk))
        query = f"SELECT proteinID, genomeID FROM Proteins WHERE proteinID IN ({placeholders})"
        cur.execute(query, chunk)
        for proteinID, genomeID in cur.fetchall():
            local_result[proteinID] = genomeID
    return local_result




#########################################################################################
################## Routines for the correlation values ##################################
#########################################################################################

def compute_conditional_presence_correlations(matrix_path, output_path=None, prefix=None):
    """
    Computes asymmetric domain co-occurrence correlations.
    Assumes matrix contains comma-separated proteinIDs, or empty strings.
    """

    with open(matrix_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        domains = header[1:]  # skip genomeID

        matrix = []
        for row in reader:
            genome_id = row[0]
            values = {domain: row[i + 1].strip() for i, domain in enumerate(domains)}
            matrix.append((genome_id, values))

    # Select domains with given prefix
    if prefix:
        condition_domains = [d for d in domains if d.startswith(prefix)]
    else:
        condition_domains = domains[:]

    result = {}

    for cond_domain in condition_domains:
        co_counts = defaultdict(int)
        total_with_cond = 0

        for genome_id, values in matrix:
            if values.get(cond_domain):  # presence = non-empty string
                total_with_cond += 1
                for other_domain, v in values.items():
                    if other_domain != cond_domain and v:
                        co_counts[other_domain] += 1

        if total_with_cond > 0:
            result[cond_domain] = {
                d: round(count / total_with_cond, 3)
                for d, count in co_counts.items()
            }

    # Optional output
    if output_path:
        with open(output_path, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["condition_domain", "co_domain", "conditional_presence"])
            for cond_domain, co_dict in result.items():
                for other_domain, ratio in co_dict.items():
                    writer.writerow([cond_domain, other_domain, ratio])

    return result




def extract_high_correlation_partners(input_tsv, min_correlation=0.5, top_n=None, output_tsv=None):
    """
    Extracts domain combinations with high conditional correlation from a correlation TSV.

    Args:
        input_tsv (str): Path to the correlation TSV file.
        min_correlation (float): Minimum correlation value to consider.
        top_n (int, optional): If set, limits number of partners per domain.
        output_tsv (str, optional): If set, writes filtered results to this TSV file.

    Returns:
        dict: {condition_domain: [(co_domain, correlation), ...]}
    """
    result = defaultdict(list)

    # Read the TSV
    with open(input_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cond = row["condition_domain"]
            co = row["co_domain"]
            try:
                corr = float(row["conditional_presence"])
                if corr >= min_correlation:
                    result[cond].append((co, corr))
            except ValueError:
                continue  # skip malformed numbers

    # Sort and optionally limit results
    for cond_domain in result:
        result[cond_domain].sort(key=lambda x: x[1], reverse=True)
        if top_n:
            result[cond_domain] = result[cond_domain][:top_n]

    # Optional write to TSV
    if output_tsv:
        with open(output_tsv, "w", newline="") as out_file:
            writer = csv.writer(out_file, delimiter="\t")
            writer.writerow(["condition_domain", "co_domain", "correlation"])
            for cond, partners in result.items():
                for co_dom, corr in partners:
                    writer.writerow([cond, co_dom, corr])

    return result







#########################################
####### Subroutines for step 5 ##########
#########################################

def load_presence_absence_matrix(matrix_path):
    presence = {}

    with open(matrix_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        clean_fieldnames = [
            col if col == "genomeID" else col.replace("grp0_", "").replace("sng0_", "")
            for col in reader.fieldnames
        ]

        for row in reader:
            genome_id = row["genomeID"]
            domain_values = {
                clean_fieldnames[i]: row[col].strip()
                for i, col in enumerate(reader.fieldnames)
                if col != "genomeID"
            }
            presence[genome_id] = domain_values

    return presence





def chunked(iterable, size):
    """Yield successive chunk-sized lists from iterable."""
    for i in range(0, len(iterable), size):
        yield list(iterable)[i:i + size]


def _process_target_domain(args):
    """Optimized worker function for a single target domain with progress output."""
    database_path, target_domain, co_domains, presence_matrix = args
    result = set()

    print(f"Searching for occurence of {target_domain} with all proteins: {co_domains}")

    # Step 1: collect genomes that contain all co_domains
    genomes_with_all = set()

    for genome_id, domain_values in presence_matrix.items():
        if all(domain_values.get(domain) for domain in co_domains):
            genomes_with_all.add(genome_id)

    try:
        with sqlite3.connect(database_path) as con:
            cur = con.cursor()
            con.execute('PRAGMA journal_mode=WAL;')
            con.execute('PRAGMA synchronous=NORMAL;')
            con.execute('PRAGMA temp_store=MEMORY;')
            con.execute('PRAGMA cache_size=-25000;')  # ca. 100MB Cache
            
            # Step 2: Get best scoring proteins for target_domain in chunks
            genome_list = list(genomes_with_all)
            chunk_size = 990
            total_chunks = (len(genome_list) + chunk_size - 1) // chunk_size

            print(f"Fetching best hits for {target_domain} in {len(genome_list)} genomes (split into {total_chunks} chunks)")

            for idx, genome_chunk in enumerate(chunked(genome_list, chunk_size), 1):
                placeholders = ','.join(['?'] * len(genome_chunk))
                query = f"""
                    SELECT P.genomeID, D.proteinID, MAX(D.score)
                    FROM Domains D
                    JOIN Proteins P ON D.proteinID = P.proteinID
                    WHERE D.domain = ? AND P.genomeID IN ({placeholders})
                    GROUP BY P.genomeID
                """
                try:
                    cur.execute(query, (target_domain, *genome_chunk))
                    for genome_id, protein_id, max_score in cur.fetchall():
                        result.add(protein_id)
                except Exception as e:
                    print(f"ERROR: [{target_domain}] Chunk {idx}/{total_chunks} failed: {e}")
                
                print(f"Processed chunk {idx}/{total_chunks} for {target_domain}")

    except Exception as e:
        print(f"[{target_domain}] DB error: {e}")

    return (target_domain, result)


def find_best_proteins_parallel(database_path, correlation_dict, matrix_path, options):
    """
    Finds best-scoring proteins per domain using multiprocessing.
    Skips domains where a grp0_<domain>.faa file already exists in the output directory.
    
    Args:
        database_path (str): Path to the SQLite database.
        correlation_dict (dict): {domain: [co_domain1, ...]} correlations (no prefix).
        matrix_path (str): Path to presence/absence matrix TSV.
        options: An object with a .cores and .fasta_output_directory attribute.

    Returns:
        dict: {target_domain: [best_protein_ids]}
    """
    presence_matrix = load_presence_absence_matrix(matrix_path)

    task_args = []
    for target_domain, co_domains in correlation_dict.items():
        output_file = os.path.join(options.fasta_output_directory, f"grp0_{target_domain}.faa")
        if os.path.exists(output_file):
            print(f"[{target_domain}] Skipped: FASTA file already exists")
            continue
        task_args.append((database_path, target_domain, co_domains, presence_matrix))

    if not task_args:
        print("No domains to process – all target FASTA files already exist.")
        return {}

    with Pool(processes=options.cores) as pool:
        results = pool.map(_process_target_domain, task_args)

    return {
        domain: protein_ids
        for domain, protein_ids in results
        if protein_ids
    }

def compute_score_limit_dict(database_path, singleton_reference_seqs_dict):
    """
    Computes lower and upper score limits for each domain based on the
    protein scores in singleton_reference_seqs_dict.

    Args:
        database_path (str): Path to SQLite database.
        singleton_reference_seqs_dict (dict): {domain_label: [proteinIDs]}

    Returns:
        dict: {domain: {"lower_limit": float, "upper_limit": float}}
    """

    # Reverse mapping: proteinID → domain
    protein_to_domain = {}
    for domain_label, proteins in singleton_reference_seqs_dict.items():
        domain = domain_label.replace("grp0_", "").replace("sng0_", "")
        for pid in proteins:
            protein_to_domain[pid] = domain

    # Prepare result container
    domain_scores = defaultdict(list)

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()
        placeholders = ",".join(["?"] * len(protein_to_domain))
        query = f"""
            SELECT proteinID, domain, score
            FROM Domains
            WHERE proteinID IN ({placeholders})
        """
        cur.execute(query, list(protein_to_domain.keys()))

        for proteinID, domain, score in cur.fetchall():
            expected_domain = protein_to_domain.get(proteinID)
            if expected_domain == domain:
                domain_scores[domain].append(score)

    # Build score limits
    score_limit_dict = {}
    for domain, scores in domain_scores.items():
        if scores:
            score_limit_dict[domain] = {
                "lower_limit": float(min(scores)),
                "upper_limit": float(max(scores))
            }

    return score_limit_dict









def generate_singleton_reference_seqs(options):
    """
    ### Main routine to the module ###
    
    Wrapper to generate singleton reference sequences based on domain correlations.
    Produces a presence/absence matrix, computes co-occurrence correlations,
    identifies best-scoring proteins from correlated domains, and exports FASTA files.

    Args:
        options: An object with attributes:
            - faa_dir
            - database_directory
            - result_files_directory
            - phylogeny_directory
            - max_seqs
            - cores
    """

    # Define intermediate file paths
    matrix_filepath = os.path.join(options.result_files_directory, "grp0_presence_absence_matrix.tsv")
    correlation_filepath = os.path.join(options.result_files_directory, "sng0_protein_presence_correlations.tsv")

    # Step 1: Create presence/absence matrix
    create_presence_absence_matrix(
        faa_dir=options.fasta_output_directory,
        database_directory=options.database_directory,
        output_path=matrix_filepath,
        chunk_size=990,
        extensions=(".faa",),
        prefixes=("grp0", "sng0"),
        cores = options.cores
    )

    # Step 2: Compute domain co-occurrence correlations
    compute_conditional_presence_correlations(
        matrix_path=matrix_filepath,
        output_path=correlation_filepath,
        prefix='sng0'  # Use 'sng0' to limit to sng0 domains only
    )

    # Step 3: Extract highly correlated co-domains
    raw_correlation_dict = extract_high_correlation_partners(
        input_tsv=correlation_filepath,
        min_correlation=options.csb_singleton_correlation,
        top_n=None,
        output_tsv=correlation_filepath
    )

    # Step 4: Remove prefixes to normalize domain names
    correlation_dict = {
        domain.replace("grp0_", "").replace("sng0_", ""): [
            co.replace("grp0_", "").replace("sng0_", "") for co, _ in co_list
        ]
        for domain, co_list in raw_correlation_dict.items()
    }

    # Step 5: Identify genomes with correlated domains and extract top proteins
    print("Identification of genomes with co-occurrence of protein and gene cluster proteins")
    singleton_reference_seqs_dict = find_best_proteins_parallel(
        database_path=options.database_directory,
        correlation_dict=correlation_dict,
        matrix_path = matrix_filepath,
        options=options
    )

    # Step 6: Re-add grp0_ prefix for downstream labeling
    singleton_reference_seqs_dict = {
        f"grp0_{k}": v for k, v in singleton_reference_seqs_dict.items()
    }

    # Step 7: Export protein sequences to FASTA
    Csb_proteins.fetch_seqs_to_fasta_parallel(
        options.database_directory,
        singleton_reference_seqs_dict,
        options.fasta_output_directory,
        5, # min seqs
        options.max_seqs,
        options.cores
    )
    
    # Step 8: Calculate score limits    
    score_limit_dict = compute_score_limit_dict(
    database_path=options.database_directory,
    singleton_reference_seqs_dict=singleton_reference_seqs_dict
    )
    


    return score_limit_dict, singleton_reference_seqs_dict




















###############################################################################
####### Main routines for the regression of presence absence matrices #########
####### Step 7                                                        #########
###############################################################################



def train_logistic_models_from_tsv(tsv_path):
    """
    Liest eine TSV-Datei mit Präsenz-/Absenz-Matrix ein, entfernt das 'grp0_'-Prefix aus den Spaltennamen,
    wandelt die Werte in binär ("" = 0, alles andere = 1), und trainiert für jede Spalte ein
    logistische Regressionsmodell auf Basis der übrigen Spalten.

    Rückgabe:
        - Dictionary mit trainierten Modellen
        - DataFrame mit den Regressionsgewichten
    """
    # Einlesen
    df = pd.read_csv(tsv_path, sep="\t")

    # Spaltennamen bereinigen
    df.columns = [col.replace("grp0_", "").replace("sng0_","") for col in df.columns]

    # Binär-Matrix erzeugen
    binary_matrix = df.apply(lambda col: col.map(lambda x: 1 if isinstance(x, str) and x != "" else 0))


    # Ergebnisse
    models = {}
    coefficients_table = {}

    # Warnungen unterdrücken (z. B. falls kleine Datensätze nicht konvergieren)
    warnings.filterwarnings("ignore", category=ConvergenceWarning)

    for target_protein in binary_matrix.columns:
        X = binary_matrix.drop(columns=[target_protein])
        y = binary_matrix[target_protein]

        if y.nunique() > 1:
            model = LogisticRegression()
            model.fit(X, y)
            models[target_protein] = model
            coefs = pd.Series(model.coef_[0], index=X.columns)
            coefs["(Intercept)"] = model.intercept_[0]
            coefficients_table[target_protein] = coefs

    coef_df = pd.DataFrame(coefficients_table).T
    return models, coef_df


def check_new_presences_deprecated(base_tsv_path, extended_tsv_path, threshold=0.6):

    # Modelle auf Basisdaten trainieren
    models, coef_df = train_logistic_models_from_tsv(base_tsv_path)

    # Erweiterte Matrix laden und vorbereiten
    extended_df = pd.read_csv(extended_tsv_path, sep="\t")
    extended_df.columns = [col.replace("grp1_", "") for col in extended_df.columns]
    extended_df = extended_df.set_index("genomeID")
    extended_bin = extended_df.apply(lambda col: col.map(lambda x: 1 if isinstance(x, str) and x != "" else 0))

    # Basisdaten zur Differenzberechnung auch laden
    base_df = pd.read_csv(base_tsv_path, sep="\t")
    base_df.columns = [col.replace("grp0_", "") for col in base_df.columns]
    base_df = base_df.set_index("genomeID")
    base_bin = base_df.apply(lambda col: col.map(lambda x: 1 if isinstance(x, str) and x != "" else 0))
    
    # Ergebnisse
    non_plausible_hits = []

    # Warnungen unterdrücken (z. B. falls kleine Datensätze nicht konvergieren)
    warnings.filterwarnings("ignore", category=ConvergenceWarning)

    for genome_id in extended_bin.index:
        if genome_id not in base_bin.index:
            # Hits in new genomes are generally unplausible
            for protein in extended_bin.columns:
                if extended_bin.at[genome_id, protein] == 1:
                    non_plausible_hits.append((genome_id, protein, 0.0))
            continue

        for protein in extended_bin.columns:
            if protein not in models:
                continue  # kein Modell für konstante Spalte

            model = models[protein]
            was_0 = base_bin.at[genome_id, protein] == 0
            now_1 = extended_bin.at[genome_id, protein] == 1

            if was_0 and now_1:
                # Feature-Vektor ohne Zielspalte
                # Features holen aus grp0 für dieses Genom, ohne Zielprotein
                X_context = base_bin.loc[genome_id].drop(labels=[protein])

                # In DataFrame verwandeln und reindexieren auf trainierte Spalten
                expected_features = model.feature_names_in_
                X_aligned = pd.DataFrame([X_context], columns=base_bin.columns)
                X_aligned = X_aligned.reindex(columns=expected_features, fill_value=0)

                # Jetzt sicher und fehlerfrei
                score = model.predict_proba(X_aligned)[0, 1]


                if score < threshold:
                    non_plausible_hits.append((genome_id, protein, score))

    return non_plausible_hits

def check_new_presences_with_combined_context(base_tsv_path, extended_tsv_path, threshold=0.6):
    """
    Bewertet neue 1-Einträge in extended_tsv_path mit einem auf base_tsv_path trainierten Modell,
    wobei zusätzlich alle neuen Präsenzen im selben Genom berücksichtigt werden.

    Gibt eine Liste der unplausiblen neuen Einträge zurück (genomeID, protein, score).
    """

    # Modelle auf Basisdaten trainieren
    models, coef_df = train_logistic_models_from_tsv(base_tsv_path)

    # Daten laden
    extended_df = pd.read_csv(extended_tsv_path, sep="\t")
    extended_df.columns = [col.replace("grp1_", "") for col in extended_df.columns]
    extended_df = extended_df.set_index("genomeID")
    extended_bin = extended_df.apply(lambda col: col.map(lambda x: 1 if isinstance(x, str) and x != "" else 0))

    base_df = pd.read_csv(base_tsv_path, sep="\t")
    base_df.columns = [col.replace("grp0_", "").replace("sng0_","") for col in base_df.columns]
    base_df = base_df.set_index("genomeID")
    base_bin = base_df.apply(lambda col: col.map(lambda x: 1 if isinstance(x, str) and x != "" else 0))

    # Ergebnisse
    non_plausible_hits = []
    warnings.filterwarnings("ignore", category=ConvergenceWarning)

    for genome_id in extended_bin.index:
        if genome_id not in base_bin.index:
            # Komplett neues Genom – alles ist unplausibel
            for protein in extended_bin.columns:
                if extended_bin.at[genome_id, protein] == 1:
                    non_plausible_hits.append((genome_id, protein, 0.0))
            continue

        # Kontext von grp0
        base_context = base_bin.loc[genome_id].copy()

        # Zusätzlich neue 1en aus grp1 identifizieren
        new_hits = []
        for protein in extended_bin.columns:
            if protein not in base_bin.columns:
                continue  # if protein in extended but not in base continue
            was_0 = base_bin.at[genome_id, protein] == 0
            now_1 = extended_bin.at[genome_id, protein] == 1
            if was_0 and now_1:
                new_hits.append(protein)

        # Kombinierter Kontext: base + neue (außer dem Target selbst beim Test)
        for target in new_hits:
            if target not in models:
                continue
            model = models[target]

            # Erweitere Kontext: kopiere grp0 Präsenz, füge neue 1er hinzu
            combined_context = base_context.copy()
            for nh in new_hits:
                if nh != target:
                    combined_context[nh] = 1

            # Align auf Modellinput
            expected_features = model.feature_names_in_
            X_aligned = pd.DataFrame([combined_context], columns=base_bin.columns)
            X_aligned = X_aligned.reindex(columns=expected_features, fill_value=0)

            score = model.predict_proba(X_aligned)[0, 1]
            if score < threshold:
                non_plausible_hits.append((genome_id, target, score))

    return non_plausible_hits


def extract_protein_ids_for_hits(tsv_path, non_plausible_hits):
    """
    Extrahiert die Protein-IDs aus der Presence/Absence-Matrix für alle unplausiblen Treffer.

    Args:
        tsv_path (str): Pfad zur Presence/Absence-Matrix (.tsv).
        non_plausible_hits (list of tuples): [(genome_id, domain, score), ...]

    Returns:
        dict: {genome_id: {domain: [proteinIDs]}}
    """
    # Matrix laden, fehlende Zellen als leere Strings behandeln
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    
    # Spaltennamen bereinigen (entferne grp1_ oder grp0_)
    original_cols = df.columns.tolist()
    cleaned_cols = [
        col if col == "genomeID" else col.replace("grp0_", "").replace("grp1_", "")
        for col in original_cols
    ]
    df.columns = cleaned_cols

    # Index setzen
    df.set_index("genomeID", inplace=True)

    # Ergebnisse sammeln
    result = defaultdict(dict)
    for genome, domain, score in non_plausible_hits:
        if genome not in df.index:
            result[genome][domain] = []
            continue

        if domain not in df.columns:
            result[genome][domain] = []
            continue

        entry = df.at[genome, domain]
        proteins = entry.split(",") if entry else []
        result[genome][domain] = [p.strip() for p in proteins if p.strip()]

    return result

def invert_protein_mapping_by_domain(protein_mapping):
    """
    Wandelt das Mapping {genome: {domain: [proteinIDs]}} um in
    {domain: set(proteinIDs)} – unabhängig vom Genom.

    Args:
        protein_mapping (dict): Mapping mit Struktur {genome: {domain: [proteinIDs]}}

    Returns:
        dict: Mapping mit Struktur {domain: set(proteinIDs)}
    """
    domain_to_proteins = {}

    for genome, domain_dict in protein_mapping.items():
        for domain, proteins in domain_dict.items():
            if domain not in domain_to_proteins:
                domain_to_proteins[domain] = set()
            domain_to_proteins[domain].update(proteins)

    return domain_to_proteins

def extract_domain_to_proteins(tsv_path):
    """
    Liest eine Presence/Absence-TSV-Datei und erzeugt ein Dictionary:
    { domain_name: set(proteinIDs) }

    Annahmen:
    - Erste Spalte ist 'genomeID'
    - Alle anderen Spalten sind Domains
    - Zellen enthalten entweder leere Strings oder Kommaseparierte ProteinIDs
    """
    df = pd.read_csv(tsv_path, sep="\t")
    
    # Entferne Präfixe aus den Domain-Namen (falls vorhanden)
    df.columns = [col.replace("grp0_", "").replace("grp1_", "") for col in df.columns]

    domain_to_proteins = {}

    for domain in df.columns[1:]:  # Überspringe genomeID
        proteins = set()
        for entry in df[domain].dropna():
            if isinstance(entry, str) and entry.strip():
                proteins.update(entry.strip().split(","))
        domain_to_proteins[domain] = proteins

    return domain_to_proteins  

def write_plausible_fasta_sequences(faa_dir, unplausible_proteins_by_domain, output_dir):
    """
    Liest grp1_*.faa Dateien, filtert unplausible Protein-IDs und schreibt plausible Sequenzen
    in neue grp2_*.faa Dateien.

    Args:
        faa_dir (str): Verzeichnis mit den grp1_*.faa Dateien.
        unplausible_proteins_by_domain (dict): {domain: set(proteinIDs)}.
        output_dir (str): Zielverzeichnis für die grp2_*.faa Dateien.
    """
    for filename in os.listdir(faa_dir):
        if not filename.startswith("grp1_") or not filename.endswith(".faa"):
            continue

        domain = filename.replace("grp1_", "").replace(".faa", "")
        unplausible_ids = unplausible_proteins_by_domain.get(domain, set())

        input_path = os.path.join(faa_dir, filename)
        output_filename = filename.replace("grp1_", "grp2_")
        output_path = os.path.join(output_dir, output_filename)

        with open(input_path, "r") as infile, open(output_path, "w") as outfile:
            write_flag = False
            current_header = ""
            current_seq_lines = []

            for line in infile:
                if line.startswith(">"):
                    # Beim nächsten Header – prüfen, ob vorheriger Header geschrieben werden soll
                    if write_flag and current_header:
                        outfile.write(current_header)
                        outfile.writelines(current_seq_lines)

                    prot_id = line[1:].strip().split()[0]
                    write_flag = prot_id not in unplausible_ids
                    current_header = line if write_flag else ""
                    current_seq_lines = []

                else:
                    if write_flag:
                        current_seq_lines.append(line)

            # Letzte Sequenz schreiben (falls gültig)
            if write_flag and current_header:
                outfile.write(current_header)
                outfile.writelines(current_seq_lines)


def main_presence_absence_matrix_filter(options):
    """
    ### Main routine to the module ###
    
    Wrapper to generate singleton reference sequences based on domain correlations.
    Produces a presence/absence matrix, computes co-occurrence correlations,
    identifies best-scoring proteins from correlated domains, and exports FASTA files.

    Args:
        options: An object with attributes:
            - faa_dir
            - database_directory
            - result_files_directory
            - phylogeny_directory
            - max_seqs
            - cores
    """

    # Define intermediate file paths
    grp1_matrix_filepath = os.path.join(options.result_files_directory, "grp1_presence_absence_matrix.tsv")

    basis_matrix_filepath = os.path.join(options.result_files_directory, "basis_presence_absence_matrix.tsv")

    # Step 1: Create presence/absence matrix for grp0 as basis co-occurence
    create_presence_absence_matrix(
        faa_dir=options.fasta_output_directory,
        database_directory=options.database_directory,
        output_path=basis_matrix_filepath,
        chunk_size=990,
        extensions=(".faa",), # has to be a tuple
        prefixes=("grp0",), # has to be a tuple
        cores = options.cores
    )

    # Step 2: Create presence/absence matrix for grp1
    create_presence_absence_matrix(
        faa_dir=options.fasta_output_directory,
        database_directory=options.database_directory,
        output_path=grp1_matrix_filepath,
        chunk_size=990,
        extensions=(".faa",), # has to be a tuple
        prefixes=("grp1",), # has to be a tuple
        cores = options.cores
    )

    # Step 2.5 check if presence absence matrices are available
    basis_exists = os.path.isfile(basis_matrix_filepath)
    grp1_exists = os.path.isfile(grp1_matrix_filepath)

    basis_nonempty = os.path.getsize(basis_matrix_filepath) > 0 if basis_exists else False
    grp1_nonempty = os.path.getsize(grp1_matrix_filepath) > 0 if grp1_exists else False

    if basis_exists and basis_nonempty and grp1_exists and grp1_nonempty:
        print("Presence/absence matrices generated and non-empty")
    else:
        print("WARNING: Presence/absence matrices missing or empty. Cannot proceed.")

        if not basis_exists:
            print(f"Missing file: {basis_matrix_filepath}")
        elif not basis_nonempty:
            print(f"File is empty: {basis_matrix_filepath}")

        if not grp1_exists:
            print(f"Missing file: {grp1_matrix_filepath}")
        elif not grp1_nonempty:
            print(f"File is empty: {grp1_matrix_filepath}")

        return
    
    
    # Step 3: Learn regression from grp0 hit distribution and check grp1 hits
    #TODO option für liberal und konservativ
    #TODO regression lernen von sng0 mit grp0 und nicht nur von grp0 damit auch eine regression für diese einheiten besteht.
    # Conservative with new hits against base grp0 hits in the genome
    #grp1_non_plausible_hits = check_new_presences(basis_matrix_filepath,grp1_matrix_filepath) # Learn from base and check individually new hits against genome base hits
    #protein_mapping = extract_protein_ids_for_hits(grp1_matrix_filepath, grp1_non_plausible_hits)
 
    # Liberal with new hits against base grp0 hits in the genome + new hits in the genome
    cache_dir = os.path.join(options.result_files_directory, "pkl_cache")
    
    if os.path.exists(cache_dir+"/grp1_non_plausible_hits2.pkl") and os.path.exists(cache_dir+"/protein_mapping2.pkl"):
        grp1_non_plausible_hits2 = Csb_proteins.load_cache(cache_dir, "grp1_non_plausible_hits2.pkl")
        protein_mapping2 = Csb_proteins.load_cache(cache_dir, "protein_mapping2.pkl")
    else:
        print("\nCalculating logistic regression of grp0 presence/absence matrix")
        grp1_non_plausible_hits2 = check_new_presences_with_combined_context(
            basis_matrix_filepath,
            grp1_matrix_filepath,
            options.pam_threshold
        )
        protein_mapping2 = extract_protein_ids_for_hits(grp1_matrix_filepath, grp1_non_plausible_hits2)

        Csb_proteins.save_cache(cache_dir, "grp1_non_plausible_hits2.pkl", grp1_non_plausible_hits2)
        Csb_proteins.save_cache(cache_dir, "protein_mapping2.pkl", protein_mapping2)
    
    grp1_hits_dict = extract_domain_to_proteins(grp1_matrix_filepath)
    basis_hits_dict = extract_domain_to_proteins(basis_matrix_filepath)
    unplausible_hits_dict = invert_protein_mapping_by_domain(protein_mapping2)


    print("\nSelection of grp1 training sequence by computed regression of the presence/absence matrix of grp0")
    for domain in sorted(unplausible_hits_dict.keys()):
        total = len(grp1_hits_dict[domain]) - len(basis_hits_dict[domain]) # Number of all added proteinIDs by MCL
        excluded = len(unplausible_hits_dict.get(domain, set()))
        included = total - excluded
        
        print(f"Total {domain} candidates {total}\tIncluded: {included}\tRemoved: {excluded}")    
    
    # Step 4: Check unplausible hits in the phylogenetic tree
    # options für branch distances als zahlen in den dateinamen einfügen
    print("\nCross-check of included training sequences by phylogenetic placement")
    if os.path.exists(cache_dir+f"/Mx_included_hits_{str(options.pam_phylogenetic_distance)}.pkl"):
        include_hits_dict = Csb_proteins.load_cache(cache_dir, f"Mx_included_hits_{str(options.pam_phylogenetic_distance)}.pkl")
        
    else:
        grp1_hits_dict = extract_domain_to_proteins(grp1_matrix_filepath)
        include_hits_dict, exluded_hits_dict = Csb_phylogeny.analyze_unplausible_proteins_in_trees_parallel(options.phylogeny_directory, unplausible_hits_dict, grp1_hits_dict, options.pam_long_branch_thres, options.pam_phylogenetic_distance, options.cores)
    
        Csb_proteins.save_cache(cache_dir, f"Mx_included_hits_{str(options.pam_phylogenetic_distance)}.pkl", include_hits_dict) # proteinIDs in included hits set are below thrs distance to a plausible grp1 sequence

    # Delete the included hits from the unplausible dict
    for domain, included_proteins in include_hits_dict.items():
        if domain in unplausible_hits_dict:
            unplausible_hits_dict[domain] -= included_proteins
            # Clean up domains with empty sets
            if not unplausible_hits_dict[domain]:
                del unplausible_hits_dict[domain]
    
    # Print out the number of excluded sequences
    print("\nSelection of grp1 training sequence by computed regression of the presence/absence matrix and phylogenetic placement")
    for domain in sorted(unplausible_hits_dict.keys()):
        total = len(grp1_hits_dict[domain]) - len(basis_hits_dict[domain]) # Number of all added proteinIDs by MCL
        excluded = len(unplausible_hits_dict.get(domain, set()))
        included = total - excluded
        
        print(f"Total {domain} candidates {total}\tIncluded: {included}\tRemoved: {excluded}")

    # Step 5: Copy grp1 files without unplausible proteinIDs to the grp2 files
    print("\nWriting grp0 and included grp1 sequences to grp2 fasta files")
    write_plausible_fasta_sequences(options.fasta_output_directory, unplausible_hits_dict, options.fasta_output_directory)
    
    return
    
    

