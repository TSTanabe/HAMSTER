#!/usr/bin/python
import os
import sqlite3
import csv
from . import Csb_proteins
from multiprocessing import Pool
from collections import defaultdict

def create_presence_absence_matrix(faa_dir, database_directory, output_path, chunk_size=500):
    """
    Generates a TSV presence/absence matrix of domains per genome without using pandas.
    It includes chunked database queries to handle large sets of protein IDs.

    Args:
        faa_dir (str): Directory containing .faa files (FASTA format with proteinIDs in headers).
        database_directory (str): Path to the SQLite database.
        output_path (str): Path to the output TSV file.
        chunk_size (int): Maximum number of proteinIDs per database query.
    """
    if os.path.isfile(output_path):
        return
    
    # 1. Collect all .faa files with the correct prefixes
    faa_files = [f for f in os.listdir(faa_dir)
                 if f.endswith(".faa") and (f.startswith("grp0") or f.startswith("sng0"))]

    domain_to_proteins = {}
    all_protein_ids = set()

    # 2. Parse protein IDs from each .faa file
    for faa_file in faa_files:
        domain = os.path.splitext(faa_file)[0]  # Use filename (without extension) as domain label
        file_path = os.path.join(faa_dir, faa_file)
        protein_ids = set()

        with open(file_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    prot_id = line[1:].strip().split()[0]
                    protein_ids.add(prot_id)
                    all_protein_ids.add(prot_id)

        domain_to_proteins[domain] = protein_ids

    # 3. Chunked query to map proteinID → genomeID from the database
    protein_to_genome = {}

    def chunked(iterable, size):
        it = list(iterable)
        for i in range(0, len(it), size):
            yield it[i:i + size]

    with sqlite3.connect(database_directory) as con:
        cur = con.cursor()
        for chunk in chunked(all_protein_ids, chunk_size):
            placeholders = ",".join(["?"] * len(chunk))
            query = f"SELECT proteinID, genomeID FROM Proteins WHERE proteinID IN ({placeholders})"
            cur.execute(query, chunk)
            for proteinID, genomeID in cur.fetchall():
                protein_to_genome[proteinID] = genomeID

    # 4. Build the presence/absence matrix: genomeID × domain
    genome_domain_matrix = defaultdict(dict)
    genome_ids = set()
    domain_list = sorted(domain_to_proteins.keys())

    for domain in domain_list:
        for prot_id in domain_to_proteins[domain]:
            genome_id = protein_to_genome.get(prot_id)
            if genome_id:
                genome_domain_matrix[genome_id][domain] = "1"
                genome_ids.add(genome_id)

    # 5. Write the matrix to TSV
    with open(output_path, "w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t")
        header = ["genomeID"] + domain_list
        writer.writerow(header)

        for genome_id in sorted(genome_ids):
            row = [genome_id]
            for domain in domain_list:
                row.append(genome_domain_matrix[genome_id].get(domain, "0"))
            writer.writerow(row)

    return


def compute_conditional_presence_correlations(matrix_path, output_path=None, prefix=None):
    """
    Computes asymmetric co-occurrence correlations between selected domains (by prefix)
    and all other domains, based on presence/absence data.

    For each domain matching the given prefix, it calculates the conditional probability
    that another domain is present, given that the selected domain is present.

    Args:
        matrix_path (str): Path to the TSV presence/absence matrix.
        output_path (str, optional): If given, results are written to this TSV file.
        prefix (str, optional): Domain prefix to filter (e.g., "sng0").
                                If None, all domains are used as condition sources.

    Returns:
        dict: {condition_domain: {co_domain: conditional_probability}}
    """
    # Read the TSV matrix
    with open(matrix_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        domains = header[1:]  # first column is genomeID

        matrix = []
        for row in reader:
            genome_id = row[0]
            values = list(map(int, row[1:]))
            matrix.append((genome_id, dict(zip(domains, values))))

    # Select condition domains (e.g., all sng0_* domains or all if no prefix)
    if prefix:
        condition_domains = [d for d in domains if d.startswith(prefix)]
    else:
        condition_domains = domains[:]

    result = {}

    for condition_domain in condition_domains:
        co_occurrence_counts = defaultdict(int)
        total_with_condition = 0

        for genome_id, domain_values in matrix:
            if domain_values.get(condition_domain) == 1:
                total_with_condition += 1
                for other_domain, presence in domain_values.items():
                    if other_domain != condition_domain:
                        co_occurrence_counts[other_domain] += presence

        if total_with_condition == 0:
            continue

        result[condition_domain] = {
            domain: round(count / total_with_condition, 3)
            for domain, count in co_occurrence_counts.items()
            if count > 0
        }

    # Optional output
    if output_path:
        with open(output_path, "w", newline="") as f_out:
            writer = csv.writer(f_out, delimiter="\t")
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



def chunked(iterable, size):
    """Yield successive chunk-sized lists from iterable."""
    for i in range(0, len(iterable), size):
        yield list(iterable)[i:i + size]


def _process_target_domain(args):
    """Optimized worker function for a single target domain with progress output."""
    database_path, target_domain, co_domains = args
    result = []
    print(f"[{target_domain}] Looking for genomes with all of: {co_domains}")
    
    try:
        with sqlite3.connect(database_path) as con:
            cur = con.cursor()

            # Step 1: Get genome sets for each co-domain (with progress print)
            # TODO dieser Schritt eins selektiert wohl noch zu viele Genome. Es sollen eigentlich nur die mit einem bitscore genommen werden, der innerhalb der score Grenzen liegt
            co_domain_genomes = []
            for idx, domain in enumerate(co_domains, 1):
                print(f"[{target_domain}] Fetching genomes for co-domain {idx}/{len(co_domains)}: {domain}")
                query = """
                    SELECT DISTINCT P.genomeID
                    FROM Domains D
                    JOIN Proteins P ON D.proteinID = P.proteinID
                    WHERE D.domain = ?
                """
                cur.execute(query, (domain,))
                genomes = {row[0] for row in cur.fetchall()}
                co_domain_genomes.append(genomes)

            if not co_domain_genomes:
                print(f"[{target_domain}] No genomes found for any co-domain.")
                return (target_domain, [])

            common_genomes = set.intersection(*co_domain_genomes)
            print(f"[{target_domain}] Found {len(common_genomes)} common genomes with all co-domains.")

            if not common_genomes:
                return (target_domain, [])
            
            # Step 2: Get best scoring proteins for target_domain in all common genomes
            genome_list = list(common_genomes)
            chunk_size = 900
            total_chunks = (len(genome_list) + chunk_size - 1) // chunk_size

            print(f"[{target_domain}] Fetching best proteins in {len(genome_list)} genomes (split into {total_chunks} chunks)")

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
                        result.append(protein_id)
                except Exception as e:
                    print(f"[{target_domain}] Chunk {idx}/{total_chunks} failed: {e}")
                
                print(f"[{target_domain}] Processed chunk {idx}/{total_chunks}")

    except Exception as e:
        print(f"[{target_domain}] Error: {e}")

    return (target_domain, result)


def find_best_proteins_parallel(database_path, correlation_dict, options):
    """
    Finds best-scoring proteins per domain using multiprocessing.

    Args:
        database_path (str): Path to the SQLite database.
        correlation_dict (dict): {domain: [co_domain1, ...]} correlations.
        options: An object or namespace with a `.cores` attribute.

    Returns:
        dict: {target_domain: [best_protein_ids]}
    """
    task_args = [
        (database_path, target_domain, co_domains)
        for target_domain, co_domains in correlation_dict.items()
    ]

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
    matrix_filepath = os.path.join(options.result_files_directory, "presence_absence_matrix.tsv")
    correlation_filepath = os.path.join(options.result_files_directory, "protein_presence_correlations.tsv")

    # Step 1: Create presence/absence matrix
    create_presence_absence_matrix(
        faa_dir=options.fasta_output_directory,
        database_directory=options.database_directory,
        output_path=matrix_filepath,
        chunk_size=900
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
        min_correlation=0.5,
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
    singleton_reference_seqs_dict = find_best_proteins_parallel(
        database_path=options.database_directory,
        correlation_dict=correlation_dict,
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
        5,
        options.max_seqs,
        options.cores
    )
    
    # Step 8: Calculate score limits    
    score_limit_dict = compute_score_limit_dict(
    database_path=options.database_directory,
    singleton_reference_seqs_dict=singleton_reference_seqs_dict
    )
    
    #TODO diese routine nimmt noch sehr viele proteinIDs auf im Schritt 5. Obwohl die korrelation wohl funktioniert, sind dabei wohl sehr viele Treffer die auch unterhalb des hit score 
    # Niveuas sind. Daher sollte d

    return score_limit_dict, singleton_reference_seqs_dict
