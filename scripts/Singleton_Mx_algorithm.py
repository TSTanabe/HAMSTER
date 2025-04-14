#!/usr/bin/python
import os
import sqlite3
import csv
from . import Csb_proteins
from multiprocessing import Pool
from collections import defaultdict

def create_presence_absence_matrix(faa_dir, database_directory, output_path, chunk_size=900, extensions=(".faa",), prefixes=("grp0", "sng0")):
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
    # Step 1: Find relevant FASTA files
    faa_files = [
        f for f in os.listdir(faa_dir)
        if f.endswith(tuple(extensions)) and f.startswith(tuple(prefixes))
    ]
    domain_to_proteins = {}
    all_protein_ids = set()
    print(faa_files)
    sys.exit()
    for faa_file in faa_files:
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

    # Step 2: Map proteinID → genomeID
    protein_to_genome = {}
    with sqlite3.connect(database_directory) as con:
        cur = con.cursor()

        def chunked(iterable, size):
            it = list(iterable)
            for i in range(0, len(it), size):
                yield it[i:i + size]

        for chunk in chunked(all_protein_ids, chunk_size):
            placeholders = ",".join(["?"] * len(chunk))
            query = f"""
                SELECT proteinID, genomeID FROM Proteins
                WHERE proteinID IN ({placeholders})
            """
            cur.execute(query, chunk)
            for proteinID, genomeID in cur.fetchall():
                protein_to_genome[proteinID] = genomeID

    # Step 3: Build genomeID → domain → proteinID(s)
    genome_domain_matrix = defaultdict(lambda: defaultdict(list))
    all_domains = sorted(domain_to_proteins.keys())

    for domain, protein_ids in domain_to_proteins.items():
        for prot_id in protein_ids:
            genome_id = protein_to_genome.get(prot_id)
            if genome_id:
                genome_domain_matrix[genome_id][domain].append(prot_id)

    # Step 4: Write TSV
    with open(output_path, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["genomeID"] + all_domains)

        for genome_id in sorted(genome_domain_matrix.keys()):
            row = [genome_id]
            for domain in all_domains:
                prot_list = genome_domain_matrix[genome_id].get(domain, [])
                row.append(",".join(prot_list) if prot_list else "")
            writer.writerow(row)

def create_presence_absence_matrix_deprecated(faa_dir, database_directory, output_path, chunk_size=500):
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

    print(f"[{target_domain}] Looking for genomes with all of: {co_domains}")

    # Step 1: collect genomes that contain all co_domains
    genomes_with_all = set()

    for genome_id, domain_values in presence_matrix.items():
        if all(domain_values.get(domain) for domain in co_domains):
            genomes_with_all.add(genome_id)

    try:
        with sqlite3.connect(database_path) as con:
            cur = con.cursor()

            # Step 2: Get best scoring proteins for target_domain in chunks
            genome_list = list(genomes_with_all)
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
                        result.add(protein_id)
                except Exception as e:
                    print(f"[{target_domain}] Chunk {idx}/{total_chunks} failed: {e}")
                
                print(f"[{target_domain}] Processed chunk {idx}/{total_chunks}")

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
    correlation_filepath = os.path.join(options.result_files_directory, "grp0_protein_presence_correlations.tsv")

    # Step 1: Create presence/absence matrix
    create_presence_absence_matrix(
        faa_dir=options.fasta_output_directory,
        database_directory=options.database_directory,
        output_path=matrix_filepath,
        chunk_size=990,
        extensions=(".faa",),
        prefixes=("grp0", "sng0")
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
        min_correlation=0.5, #TODO make this an option
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
    matrix_filepath = os.path.join(options.result_files_directory, "grp1_presence_absence_matrix.tsv")
    correlation_filepath = os.path.join(options.result_files_directory, "grp1_protein_presence_correlations.tsv")

    # Step 1: Create presence/absence matrix
    create_presence_absence_matrix(
        faa_dir=options.fasta_output_directory,
        database_directory=options.database_directory,
        output_path=matrix_filepath,
        chunk_size=990,
        extensions=(".faa",), # has to be a tuple
        prefixes=("grp1",) # has to be a tuple
    )

    # Step 2: Compute domain co-occurrence correlations
    compute_conditional_presence_correlations(
        matrix_path=matrix_filepath,
        output_path=correlation_filepath,
        prefix='grp1'  # Use 'grp1' to limit to grp1 domains only
    )

    # Step 3: Extract highly correlated co-domains
    raw_correlation_dict = extract_high_correlation_partners(
        input_tsv=correlation_filepath,
        min_correlation=0.5, #TODO make this an option
        top_n=None,
        output_tsv=correlation_filepath
    )
    
    #TODO hier weiter machen. Nummer 7 testen und die korrelationen in abhängigkeit von einander nutzen um in der presence absence matrix aufzuräumen
    return
    
    
    # Wie soll das danach ablaufen?
    # Meine eigene initiative wäre es...
    # markiere die domänen die nur in grp1 vorkommen
    # dann schau in der zeile, welche anderen domänen noch vorkommen
    # dann schau nach wie wahrscheinlich ist es, dass diese domänen vorkommen
    # addiere die wahrscheinlichkeiten zu einer wahrscheinlichkeits summe. 
    # Diese summe soll angeben, wie wahrscheinlich es ist, dass es ein zufall ist. 
    
    # Step 4: Remove prefixes to normalize domain names
    correlation_dict = {
        domain.replace("grp1_", ""): [
            co.replace("grp1_", "") for co, _ in co_list
        ]
        for domain, co_list in raw_correlation_dict.items()
    }

    # Step 5: Identify genomes with correlated domains and extract top proteins
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
