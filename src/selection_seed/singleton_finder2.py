#!/usr/bin/python

import copy
import os
import sqlite3
import pandas as pd
from collections import defaultdict
from typing import Dict, Any, Set, Tuple

from src.selection_defragmentation import pam_mx_algorithm
from src.core import myUtil

logger = myUtil.logger

"""
I want a new selection process for singletons. First the ones without a context must be identified

For this task take all hits with >80 % identity to the query and see if there are some examples without clusterID

Then take these genomes and take all proteins that have an identity of 80% or higher (adjustable by input)

Consider these as co occurring patterns. Select proteins that have the same co occurence pattern and at least 70 % 
identity to query. Do this for each protein.
"""


########################################
# New co-occurrence based singleton selection
########################################


def _find_context_free_high_identity_hits(
    database_path: str,
    identity_cutoff: float = 80.0,
) -> dict[str, set[str]]:
    """
    Find for each domain the genomes where identity >= cutoff
    AND clusterID IS NULL
    AND genomeID is not QUERY.

    Returns:
        dict:
            {
                domain: {genomeID1, genomeID2, ...},
                ...
            }
    """
    context_free: dict[str, set[str]] = defaultdict(set)

    query = """
        SELECT DISTINCT d.domain, p.genomeID
        FROM Domains d
        JOIN Proteins p ON p.proteinID = d.proteinID
        WHERE d.identity >= ?
          AND p.clusterID IS NULL
          AND p.genomeID != 'QUERY'
    """

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()
        for domain, genome_id in cur.execute(query, (identity_cutoff,)):
            context_free[domain].add(genome_id)

    return context_free


def _fetch_high_identity_domain_intersection(
    database_path: str,
    domain_to_genomes: Dict[str, Set[str]],
    min_identity_cutoff: float = 70.0,
) -> Dict[str, Set[str]]:
    """
    Für jede gegebene Domain und deren zugehörige Genome:
      - hole alle Domaintypen mit identity >= min_identity_cutoff in diesen Genomen
      - berechne die Schnittmenge der Domain-Sets über die Genome dieser Domain

    Args:
        database_path: Pfad zur SQLite-Datenbank.
        domain_to_genomes: dict {seed_domain: set(genomeIDs)},
                           z.B. Ausgabe von _find_context_free_high_identity_hits().
        min_identity_cutoff: minimale Identity in Prozent.

    Returns:
        dict:
            {
                seed_domain: {domain1, domain2, ...},  # Domains, die in ALLEN zugehörigen Genomen
                                                     # mit identity >= cutoff vorkommen
                ...
            }
    """

    # Ergebnis: für jede Seed-Domain die Schnittmenge der Co-Occurrence-Domains
    intersections_per_seed: Dict[str, Set[str]] = {}

    if not domain_to_genomes:
        return intersections_per_seed

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()
        chunk_size = 900  # wegen SQLite-Placeholder-Limit

        for seed_domain, genomes in domain_to_genomes.items():
            genomes_list = [g for g in genomes if g]
            if not genomes_list:
                intersections_per_seed[seed_domain] = set()
                continue

            # pro Genom: welche Domains erfüllen den Identity-Cutoff?
            domains_per_genome: Dict[str, Set[str]] = defaultdict(set)

            # DB in Chunks abfragen
            for i in range(0, len(genomes_list), chunk_size):
                chunk = genomes_list[i : i + chunk_size]
                placeholders = ",".join(["?"] * len(chunk))

                query = f"""
                    SELECT DISTINCT d.domain, p.genomeID
                    FROM Domains d
                    JOIN Proteins p ON p.proteinID = d.proteinID
                    WHERE d.identity >= ?
                      AND p.genomeID IN ({placeholders})
                      AND p.genomeID != 'QUERY'
                """
                cur.execute(query, (min_identity_cutoff, *chunk))
                for domain, genome_id in cur.fetchall():
                    domains_per_genome[genome_id].add(domain)

            # Schnittmenge über alle Genome dieser Seed-Domain bilden
            intersection: Set[str] | None = None
            for genome_id in genomes_list:
                doms = domains_per_genome.get(genome_id, set())
                if intersection is None:
                    intersection = set(doms)
                else:
                    intersection &= doms

                # Early exit: wenn leer, kann nichts mehr dazu kommen
                if not intersection:
                    break

            intersections_per_seed[seed_domain] = intersection or set()

    return intersections_per_seed


def select_singleton_refs_by_domain_pattern(
    database_path: str,
    seed_to_pattern_domains: Dict[str, Set[str]],
    min_identity_cutoff: float = 70.0,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Set[str]]]:
    limits_dict: Dict[str, Dict[str, float]] = {}
    sng_reference_seq_dict: Dict[str, Set[str]] = defaultdict(set)

    if not seed_to_pattern_domains:
        return limits_dict, sng_reference_seq_dict

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        # TEMP table once
        cur.execute(
            "CREATE TEMP TABLE IF NOT EXISTS tmp_pattern_domains (domain TEXT PRIMARY KEY)"
        )

        for seed_domain, pattern_domains in seed_to_pattern_domains.items():
            if not pattern_domains:
                continue

            # fill temp table
            cur.execute("DELETE FROM tmp_pattern_domains")
            cur.executemany(
                "INSERT OR IGNORE INTO tmp_pattern_domains(domain) VALUES (?);",
                ((d,) for d in pattern_domains),
            )

            # 1) Genome that have ALL pattern domains
            query_genomes = """
                SELECT p.genomeID
                FROM Domains d
                JOIN Proteins p ON p.proteinID = d.proteinID
                JOIN tmp_pattern_domains t ON t.domain = d.domain
                WHERE d.identity >= ?
                  AND p.genomeID != 'QUERY'
                GROUP BY p.genomeID
                HAVING COUNT(DISTINCT d.domain) = (SELECT COUNT(*) FROM tmp_pattern_domains)
            """

            cur.execute(query_genomes, (min_identity_cutoff,))
            genomes = [row[0] for row in cur.fetchall()]

            if not genomes:
                continue

            # 2) fetch ALL seed_domain hits in these genomes
            all_proteinIDs = set()
            all_bitscores = []

            chunk_size = 900
            for i in range(0, len(genomes), chunk_size):
                chunk = genomes[i : i + chunk_size]
                placeholders = ",".join(["?"] * len(chunk))

                query_hits = f"""
                    SELECT d.proteinID, d.score
                    FROM Domains d
                    JOIN Proteins p ON p.proteinID = d.proteinID
                    WHERE d.domain = ?
                      AND d.identity >= ?
                      AND p.genomeID IN ({placeholders})
                      AND p.genomeID != 'QUERY'
                """
                cur.execute(query_hits, (seed_domain, min_identity_cutoff, *chunk))
                for protein_id, score in cur.fetchall():
                    all_proteinIDs.add(protein_id)
                    try:
                        all_bitscores.append(float(score))
                    except:
                        pass

            if not all_proteinIDs or not all_bitscores:
                continue

            sng_reference_seq_dict[seed_domain] = all_proteinIDs
            limits_dict[seed_domain] = {
                "lower_limit": min(all_bitscores),
                "upper_limit": max(all_bitscores),
            }

    return limits_dict, sng_reference_seq_dict


#### Main routine of this module
def singleton_reference_finder(
    options: Any,
) -> (
    tuple[object , object]
    | tuple[dict[str, dict[str, float]], dict[str, set[str]]]
):
    """
    Main routine: finds and predicts reference singletons for each protein/domain.

    Args:
        options (Any): Configuration/options.

    Returns:
        tuple:
            - domain_score_limits (dict): {domain: {'lower_limit', 'upper_limit'}}
            - singleton_reference_seqs_dict (dict): {domain: set(proteinIDs)}
    """

    high_identity_cutoff = getattr(options, "singleton_identity_cutoff", 80.0)
    pattern_identity_cutoff = getattr(
        options, "singleton_pattern_identity_cutoff", 70.0
    )

    # 1) Find genomes & domains with context-free high-identity hits
    context_free_domains_dict = _find_context_free_high_identity_hits(
        options.database_directory,
        identity_cutoff=high_identity_cutoff,
    )

    logger.info(f"Found {len(context_free_domains_dict)} context-free hits with {high_identity_cutoff} percent identity")

    if not context_free_domains_dict:
        logger.warning(
            "No context-free high-identity hits found – stopping singleton selection."
        )
        return {}, {}

    # 2) For these genomes, fetch all hits with identity >= pattern_identity_cutoff
    domain_presence_intersection_pattern = _fetch_high_identity_domain_intersection(
        database_path=options.database_directory,
        domain_to_genomes=context_free_domains_dict,
        min_identity_cutoff=pattern_identity_cutoff,  # identity cutoff for co occurring pattern domains
    )

    if not domain_presence_intersection_pattern:
        logger.warning(
            f"No presence/absence pattern detected for hits with {high_identity_cutoff} percent identity – stopping singleton selection."
        )
        return {}, {}

    # 3) Select co-occurrence-based singleton candidates
    domain_score_limits, singleton_reference_seqs_dict = (
        select_singleton_refs_by_domain_pattern(
            database_path=options.database_directory,
            seed_to_pattern_domains=domain_presence_intersection_pattern,
            min_identity_cutoff=pattern_identity_cutoff,
        )
    )

    myUtil.save_cache(
        options, "sng0_training_proteinIDs.pkl", singleton_reference_seqs_dict
    )
    myUtil.save_cache(
        options, "sng0_training_proteinIDs_limits.pkl", domain_score_limits
    )

    return domain_score_limits, singleton_reference_seqs_dict


#######################


def predict_singleton_reference_seqs_for_each_domain(
    database_path: str,
    grouped: Dict[str, Set[str]],
    singleton: Dict[str, Set[str]],
    cores: int,
    chunk_size: int = 900,
) -> Dict[str, pd.Series]:
    """
    For each singleton protein a linear regression model is generated. This is based on the combined presence absence matrix
    of singleton domains + the proteins from csb patterns where at least one reference query sequence occurs in the
    pattern. The logistic regression model is then used to predict the presence of the singleton in all genomes where a similar pam
    pattern exists.

    For each singleton protein, trains a logistic regression model and predicts probability of presence for all genomes.

    Args:
        database_path (str): Path to SQLite DB.
        grouped (dict): {domain: set(proteinIDs)}
        singleton (dict): {singleton_domain: set(proteinIDs)}
        cores (int): CPUs for parallelism.
        chunk_size (int): SQLite chunk size.

    Returns:
        dict: {singleton_domain: pd.Series(predicted plausibility for each genome)}

    Example:
        {'ABC': pd.Series(...), ...}
    """

    predictions_all = {}

    for sng_domain, sng_proteins in singleton.items():
        logger.info(f"Processing protein without conserved context {sng_domain}")

        # 1. Kombiniertes grouped: Basis + aktuelles Singleton
        grouped_plus = copy.deepcopy(grouped)
        grouped_plus[sng_domain] = sng_proteins

        # 2. Presence/Absence-Matrix erzeugen
        pam = Pam_Mx_algorithm.create_presence_absence_matrix(
            grouped_plus,
            database_directory=database_path,
            output=f"pam_with_{sng_domain}",
            chunk_size=chunk_size,
            cores=cores,
        )

        # 3. PAM-Filter: nur Genomes behalten, die Singleton enthalten
        genomeIDs = [
            genomeID for genomeID, domains in pam.items() if sng_domain in domains
        ]
        filtered_pam = {
            genomeID: pam[genomeID] for genomeID in genomeIDs if genomeID in pam
        }
        if not filtered_pam:
            logger.warning(f"No genomes with presence of {sng_domain} found – skipping")
            continue

        # 4. Hit-Scores für alle ProteinIDs laden
        bsr_hit_scores = Pam_Mx_algorithm.fetch_bsr_scores(
            database_path, grouped_plus, chunk_size=chunk_size
        )

        # 5. Modell trainieren auf gefilterter PAM
        models, _ = Pam_Mx_algorithm.train_logistic_from_pam_with_scores(
            filtered_pam, bsr_hit_scores, cores=cores, target_domains={sng_domain}
        )

        if sng_domain not in models:
            logger.warning(models)
            logger.warning(f"No model trained for {sng_domain} – skipping")
            continue

        # 6. Testdaten: PAM nur für grouped (ohne das aktuelle Singleton)
        base_pam = Pam_Mx_algorithm.create_presence_absence_matrix(
            grouped,
            database_directory=database_path,
            output=f"pam_base_for_{sng_domain}",
            chunk_size=chunk_size,
            cores=cores,
        )

        # 7. In DataFrame konvertieren für Vorhersage
        genomes = sorted(base_pam.keys())
        domains = sorted({d for m in base_pam.values() for d in m})
        df = pd.DataFrame(index=genomes, columns=domains)

        for genome_id in genomes:
            for domain in domains:
                df.at[genome_id, domain] = (
                    ",".join(base_pam[genome_id].get(domain, []))
                    if domain in base_pam[genome_id]
                    else ""
                )

        # 8. Feature-Matrix bauen über Hilfsfunktion
        X_test = Pam_Mx_algorithm.build_presence_score_matrix(df, bsr_hit_scores)

        # 9. Nur Modellrelevante Features
        model = models[sng_domain]
        X_test = X_test.reindex(columns=model.feature_names_in_, fill_value=0)

        # 10. Vorhersagen berechnen
        prediction_scores = pd.Series(model.predict_proba(X_test)[:, 1], index=df.index)
        predictions_all[sng_domain] = prediction_scores

    return predictions_all


def collect_predicted_singleton_hits_from_db(
    predictions_all: Dict[str, pd.Series],
    database_path: str,
    plausible_cutoff: float = 0.6,
    bsr_threshold: float = 0.5,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Set[str]]]:
    """
    Collects the best protein hit for each plausible singleton prediction, based on blast score ratio and bitscore.

    Args:
        predictions_all (dict): {singleton_domain: pd.Series(scores)}
        database_path (str): SQLite DB path.
        plausible_cutoff (float): Minimum prediction score.
        bsr_threshold (float): Minimum BLAST score ratio.

    Returns:
        tuple:
            - singleton_score_limits (dict): {domain: {'lower_limit', 'upper_limit'}}
            - singleton_hits (dict): {domain: set(proteinIDs)}

    Example:
        ({'ABC': {'lower_limit': 123, ...}}, {'ABC': {'prot1', ...}})
    """

    singleton_hits = defaultdict(set)
    singleton_score_limits = {}

    # Prepare the structure of the output dictionary domain => plausible genome genomeIDs > plausible cutoff
    domain_to_genomes = defaultdict(list)
    for singleton, prediction_series in predictions_all.items():
        domain = singleton.replace("sng0_", "")
        plausible_genomes = prediction_series[
            prediction_series > plausible_cutoff
        ].index
        domain_to_genomes[domain].extend(plausible_genomes)

    if not domain_to_genomes:
        return {}, {}

    # iterate over the proteins and get the best hits per genome
    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        for domain, genome_list in domain_to_genomes.items():
            best_hits = {}
            for i in range(0, len(genome_list), 900):  # Max 999 placeholders in SQLite
                chunk = genome_list[i : i + 900]
                placeholders = ",".join(["?"] * len(chunk))

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

            # Save best hit per protein
            bitscores = []
            for genome_id, (bsr, proteinID, score) in best_hits.items():
                singleton_hits[domain].add(proteinID)
                bitscores.append(score)

            if bitscores:
                singleton_score_limits[domain] = {
                    "lower_limit": min(bitscores),
                    "upper_limit": max(bitscores),
                }

    return singleton_score_limits, singleton_hits


################ Routines for first main routine singleton _reference_sequences ##########################


def extract_protein_ids_from_fasta(fasta_path: str) -> Set[str]:
    """
    Extracts all protein IDs from a FASTA file.

    Args:
        fasta_path (str): Path to the FASTA file.

    Returns:
        set: Set of protein IDs.

    Example:
        >prefix___ABC123 desc
        >DEF456
        Output: {"ABC123", "DEF456"}
    """

    protein_ids = set()
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                full_id = line[1:].split()[0]
                parts = full_id.split("___")
                protein_id = parts[1] if len(parts) > 1 else full_id
                protein_ids.add(protein_id)
    return protein_ids


def extract_domain_names_from_directory(directory_path: str) -> Set[str]:
    """
    Extracts all unique domain names from files in a directory (removes extension and prefix up to first "_").

    Args:
        directory_path (str): Path to directory.

    Returns:
        set: All unique domain names.
    """
    domain_names = set()

    for filename in os.listdir(directory_path):
        if os.path.isfile(os.path.join(directory_path, filename)):
            name_no_ext = os.path.splitext(filename)[0]  # z.B. "Query_ABC123"
            parts = name_no_ext.split("_", 1)  # nur am ersten "_" splitten
            domain = parts[1] if len(parts) > 1 else name_no_ext
            domain_names.add(domain)

    return domain_names


def get_min_bitscore_for_query(
    report_path: str, query_id: str, blast_score_ratio: float = 0.9
) -> Tuple[float, float]:
    """
    Reads a DIAMOND BLAST report and returns the min and max bitscore for a query ID.

    Args:
        report_path (str): Path to DIAMOND .tab file.
        query_id (str): Query ID to search for.
        blast_score_ratio (float): Used for single-hit fallback.

    Returns:
        tuple: (min_bitscore, max_bitscore)
    """

    bitscores = []
    min_cutoff = 100
    max_cutoff = 100
    with open(report_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            hit_id, qid, evalue, bitscore = parts[:4]
            hit_id = hit_id.split("___")[1]
            if hit_id == query_id and qid == query_id:
                try:
                    bitscores.append(float(bitscore))
                except ValueError:
                    continue  # ungültiger bitscore

    if len(bitscores) == 1:
        max_cutoff = bitscores[0] * (2 - blast_score_ratio)
        min_cutoff = (
            bitscores[0] * blast_score_ratio
        )  # if only a single sequence is given than 0.9 blast score ratio is accepted
    else:
        max_cutoff = max(bitscores)
        min_cutoff = min(bitscores)
    return min_cutoff, max_cutoff
