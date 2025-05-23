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


def select_hits_by_pam_csb_mcl(options, mcl_clustering_results_dict, basis_grouped):
    """
    ### Main routine to the module ###
    
    Get all Proteins from MCL clusters with a reference sequence
    
    Get the csb vicinity from all non-ref seqs => definiert durch die keys in grouped. Diese kommen schließlich aus Genclustern
    
    Get the csb vinitiy from all ref-seqs => Müssen aus der Datenbank gelesen werden
    
    Select all non-ref seqs that have a partial ref-seq genomic vicinity into the extended_refseqs
    
    Wie bisher auch get the sequences that are plausible by liberal selection into the extended_refseqs    
    
    

    Args:
        options: An object with attributes:
            - faa_dir
            - database_directory
            - result_files_directory
            - phylogeny_directory
            - max_seqs
            - cores
    """


    grouped_3_dict = myUtil.load_cache(options,'mcl_truncated_csb_hits.pkl')
    grouped_4_dict = myUtil.load_cache(options,'mcl_PAM_plausible_hits.pkl')
    extended_grouped = myUtil.load_cache(options,'mcl_PAM_csb_merged_hits.pkl')

    if extended_grouped:
        return extended_grouped
    else:
        extended_grouped = {}
    
    
    # Process reference_dict to extract the actual domain names
    processed_reference_dict = {
        key.split("_", 1)[-1]: value for key, value in basis_grouped.items()
    }

    common_domains = mcl_clustering_results_dict.keys()
    
    total = len(mcl_clustering_results_dict)
    
    if not grouped_3_dict:
        grouped_3_dict = {}
        print("Selection truncated csb")
        for i, (domain, mcl_file) in enumerate(mcl_clustering_results_dict.items(), 1):
            print(f"[{i}/{total}] Processing domain: {domain}")

            # Get reference sequences for the domain (if exists in processed reference dict)
            reference_sequences = processed_reference_dict.get(domain, set())
            if not reference_sequences:
                print(f"Warning: No reference sequences found for domain '{domain}', skipping.")
                continue
            
            mcl_cluster_protID_set = select_ref_seq_mcl_sequences(mcl_file, domain, reference_sequences)
            
            # Save for later use in the PAM calculation
            extended_grouped[domain] = mcl_cluster_protID_set
            
            # only the new proteinIDs without the refseq are processed for common csb
            new_proteinID_set = mcl_cluster_protID_set - reference_sequences


            # Select the proteins with a truncated csb, that fits the other csbs
            filtered_proteinIDs = select_seqs_with_truncated_csb_vicinity(options.database_directory, new_proteinID_set, common_domains)
            grouped_3_dict[domain] = filtered_proteinIDs

    
    print("Selection by PAM matrix")
    # Select the proteins with plausible PAM
    if not grouped_4_dict:
        grouped_4_dict = select_proteins_with_plausible_pam_from_mcl(options.database_directory, processed_reference_dict, extended_grouped, options.cores)
    
    # Print statistics on selection to the terminal
    log_all_mcl_cluster_statistics(processed_reference_dict, extended_grouped, grouped_3_dict, grouped_4_dict)
    
    merged_dict = myUtil.merge_grouped_refseq_dicts_simple(grouped_3_dict,grouped_4_dict)
    
    merged_dict = myUtil.merge_grouped_refseq_dicts_simple(processed_reference_dict,merged_dict)
    
    myUtil.save_cache(options, "mcl_truncated_csb_hits.pkl", grouped_3_dict)
    myUtil.save_cache(options, "mcl_PAM_plausible_hits.pkl", grouped_4_dict)
    myUtil.save_cache(options, "mcl_PAM_csb_merged_hits.pkl", merged_dict)
            
    return merged_dict

def select_ref_seq_mcl_sequences(mcl_file, domain, reference_sequences):
    """
    Returns all protein IDs from MCL clusters that contain at least one reference sequence.

    Args:
        mcl_file (str): Path to the MCL output file.
        domain (str): Domain name (used for logging/debugging).
        reference_sequences (set): Set of known reference protein IDs.

    Returns:
        set: All protein IDs from clusters that contain at least one reference sequence.
    """
    
    def parse_mcl_clusters(path):
        # Read MCL file and return a list of clusters (each cluster is a list of protein IDs)
        with open(path, 'r') as f:
            return [line.strip().split() for line in f if line.strip()]

    clusters = parse_mcl_clusters(mcl_file)

    all_protein_ids_with_refseqs = set()

    for cluster in clusters:
        cluster_set = set(cluster)

        # If the cluster contains at least one reference protein, add all proteins from that cluster
        if cluster_set & reference_sequences:
            all_protein_ids_with_refseqs.update(cluster_set)

    #print(f"[{domain}] Clusters with reference sequences found: {len(all_protein_ids_with_refseqs)} proteins total")
    
    
    return all_protein_ids_with_refseqs



def select_seqs_with_truncated_csb_vicinity(database_path, protein_ids, common_domains, chunk_size=900):
    """
    Select proteins whose neighbors (same clusterID) contain only domains from the common_domains set.
    
    Args:
        database_path (str): Path to the SQLite database.
        protein_ids (set): Set of proteinIDs to evaluate.
        common_domains (set): Set of acceptable domains.
        chunk_size (int): Max size for IN queries to avoid SQLite placeholder limits.
    
    Returns:
        set: ProteinIDs that have only valid domains in their cluster neighborhood.
    """
    selected = set()
    
    def chunked(iterable, size):
        for i in range(0, len(iterable), size):
            yield list(iterable)[i:i + size]

    with sqlite3.connect(database_path) as con:
        cur = con.cursor()

        # Step 1: Fetch clusterIDs for proteinIDs in chunks
        pid_to_cluster = {}
        for chunk in chunked(protein_ids, chunk_size):
            placeholders = ','.join('?' for _ in chunk)
            cur.execute(f"""
                SELECT proteinID, clusterID 
                FROM Proteins 
                WHERE proteinID IN ({placeholders})
                  AND clusterID IS NOT NULL
            """, tuple(chunk))
            pid_to_cluster.update(cur.fetchall())

        if not pid_to_cluster:
            return selected

        cluster_ids = set(pid_to_cluster.values())

        # Step 2: Fetch all proteinIDs from those clusters in chunks
        cluster_to_proteins = {}
        for chunk in chunked(cluster_ids, chunk_size):
            placeholders = ','.join('?' for _ in chunk)
            cur.execute(f"""
                SELECT clusterID, proteinID 
                FROM Proteins 
                WHERE clusterID IN ({placeholders})
            """, tuple(chunk))
            for clusterID, proteinID in cur.fetchall():
                cluster_to_proteins.setdefault(clusterID, set()).add(proteinID)

        # Step 3: Fetch all domains for all involved proteinIDs in chunks
        all_neighbor_proteins = {p for ps in cluster_to_proteins.values() for p in ps}
        protein_to_domains = {}
        for chunk in chunked(all_neighbor_proteins, chunk_size):
            placeholders = ','.join('?' for _ in chunk)
            cur.execute(f"""
                SELECT proteinID, domain 
                FROM Domains 
                WHERE proteinID IN ({placeholders})
            """, tuple(chunk))
            for proteinID, domain in cur.fetchall():
                protein_to_domains.setdefault(proteinID, set()).add(domain)

        # Step 4: Check each input protein for domain validity in its cluster
        for pid in protein_ids:
            cluster = pid_to_cluster.get(pid)
            if not cluster:
                continue
            neighbor_proteins = cluster_to_proteins.get(cluster, set())
            domains_in_cluster = set()
            for neighbor in neighbor_proteins:
                domains_in_cluster.update(protein_to_domains.get(neighbor, set()))

            if domains_in_cluster.issubset(common_domains):
                selected.add(pid)

    return selected


#############################################################################################################



def select_proteins_with_plausible_pam_from_mcl(database_path, basis_grouped, extended_grouped, cores, chunk_size=900):

    pam_plausible_all = {}

    # 1. Basis PAM berechnen
    global_pam = Pam_Mx_algorithm.create_presence_absence_matrix(
        basis_grouped,
        database_directory=database_path,
        output="basic_pam_with_all_domains",
        chunk_size=chunk_size,
        cores=cores
    )
    
    
    # 2. Extended PAM berechnen
    extended_pam = Pam_Mx_algorithm.create_presence_absence_matrix(
        extended_grouped,
        database_directory=database_path,
        output="extended_pam_with_all_domains",
        chunk_size=chunk_size,
        cores=cores
    )    
    
    
    # 3. Calculate logistic regression Models for basis PAM
    basis_bsr_hit_scores = Pam_Mx_algorithm.fetch_bsr_scores(database_path, basis_grouped, chunk_size=chunk_size)
    
    models, _ = Pam_Mx_algorithm.train_logistic_from_pam_with_scores(global_pam, basis_bsr_hit_scores, cores=cores)
    
    
    # 4. Calculate plausibility for the extended proteinIDs from the MCL
    extended_bsr_hit_scores = Pam_Mx_algorithm.fetch_bsr_scores(database_path, extended_grouped, chunk_size=chunk_size)
    
    plausible_genomes_df = validate_existing_domains_with_models(database_path, extended_pam, models, extended_bsr_hit_scores, cutoff=0.5) #cutoff muss noch festgelegt werden
    
    
    # 5. Extract the proteinIDs from the PAM for each genome with a plausible hit
    plausible_pam_proteinID_dict = extract_proteinIDs_from_plausible_predictions(plausible_genomes_df, extended_pam)
    
    return plausible_pam_proteinID_dict
    
    
    
    
    
def validate_existing_domains_with_models(database_path, extended_pam, models, bsr_hit_scores, cutoff=0.5):
    """
    Validates domain presence in the extended_pam using logistic models trained on basis_pam.
    Only domains that are already present in the extended_pam are evaluated.

    Args:
        database_path (str): Path to the SQLite database.
        extended_pam (dict): PAM structure {genomeID: {domain: [proteinIDs]}}.
        models (dict): Logistic regression models {domain: model}.
        bsr_hit_scores (dict): Score lookup {proteinID: BSR score}.
        cutoff (float): Probability threshold to consider a domain plausible.

    Returns:
        dict: {domain: pd.Series with plausible genomeIDs and probabilities}
    """


    # Prepare dataframe from extended PAM
    genomes = sorted(extended_pam.keys())
    features = sorted({d for doms in extended_pam.values() for d in doms})
    df = pd.DataFrame(index=genomes, columns=features)

    for genome_id in genomes:
        for feat in features:
            df.at[genome_id, feat] = ",".join(extended_pam[genome_id].get(feat, [])) if feat in extended_pam[genome_id] else ""

    X_all = Pam_Mx_algorithm.build_presence_score_matrix(df, bsr_hit_scores)

    plausible_results = {}

    for domain, model in models.items():
        if model is None:
            continue

        # Only evaluate genomes where the domain is present
        target_genomes = [gid for gid in extended_pam if domain in extended_pam[gid]]
        if not target_genomes:
            continue

        # Copy the PAM and remove the domain to prevent leakage
        X = X_all.loc[target_genomes]
        X = X.reindex(columns=model.feature_names_in_, fill_value=0)

        # Predict probabilities
        preds = pd.Series(model.predict_proba(X)[:, 1], index=X.index)

        # Filter by cutoff
        plausible_preds = preds[preds > cutoff]
        if not plausible_preds.empty:
            plausible_results[domain] = plausible_preds

    return plausible_results


def extract_proteinIDs_from_plausible_predictions(plausible_results, extended_pam):
    """
    Extracts proteinIDs from the extended PAM for genome/domain pairs that are marked as plausible.

    Args:
        plausible_results (dict): {domain: pd.Series of genomeID → probability}
        extended_pam (dict): {genomeID: {domain: [proteinIDs]}}

    Returns:
        dict: {domain: set(proteinIDs)} – same structure as "grouped"
    """
    output_grouped = {}

    for domain, genome_series in plausible_results.items():
        proteins = set()
        for genome_id in genome_series.index:
            dom_to_prots = extended_pam.get(genome_id, {})
            proteins.update(dom_to_prots.get(domain, []))
        if proteins:
            output_grouped[domain] = proteins

    return output_grouped



def dead_code():
    # Define intermediate file paths
    grp1_matrix_filepath = os.path.join(options.result_files_directory, "grp1_presence_absence_matrix.tsv")

    basis_matrix_filepath = os.path.join(options.result_files_directory, "basis_presence_absence_matrix.tsv")

    # Step 1: Create presence/absence matrix for grp0 as basis co-occurence
    if not os.path.isfile(basis_matrix_filepath):
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
    if not os.path.isfile(grp1_matrix_filepath):
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
    # Conservative with new hits against base grp0 hits in the genome
    #grp1_non_plausible_hits = check_new_presences(basis_matrix_filepath,grp1_matrix_filepath) # Learn from base and check individually new hits against genome base hits
    #protein_mapping = extract_protein_ids_for_hits(grp1_matrix_filepath, grp1_non_plausible_hits)
# 
    # Liberal with new hits against base grp0 hits in the genome + new hits in the genome
    grp1_non_plausible_hits = myUtil.load_cache(options, "grp1_non_plausible_hits.pkl")
    protein_mapping = myUtil.load_cache(options, "grp1_protein_mapping.pkl")
    if not grp1_non_plausible_hits or not protein_mapping:
        print("\nCalculating logistic regression of grp0 presence/absence matrix")
        grp1_non_plausible_hits = check_new_presences_with_combined_context(
            basis_matrix_filepath,
            grp1_matrix_filepath,
            options.pam_threshold
        )
        protein_mapping = extract_protein_ids_for_hits(grp1_matrix_filepath, grp1_non_plausible_hits)

        myUtil.save_cache(options, "grp1_non_plausible_hits.pkl", grp1_non_plausible_hits)
        myUtil.save_cache(options, "grp1_protein_mapping.pkl", protein_mapping)

    grp1_hits_dict = extract_domain_to_proteins(grp1_matrix_filepath)
    basis_hits_dict = extract_domain_to_proteins(basis_matrix_filepath)
    unplausible_hits_dict = invert_protein_mapping_by_domain(protein_mapping)


    print("\nSelection of grp1 training sequence by computed regression of the presence/absence matrix of grp0")
    summarize_domain_filtering(grp1_hits_dict, basis_hits_dict, unplausible_hits_dict)
    
    # Step 4: Check unplausible hits in the phylogenetic tree
    print("\nCross-check of included training sequences by phylogenetic placement")

    cache_name = f"grp2_included_hits_pd_{options.pam_phylogenetic_distance}.pkl"

    include_hits_dict = myUtil.load_cache(options, cache_name)
    if not include_hits_dict:
        grp1_hits_dict = extract_domain_to_proteins(grp1_matrix_filepath)

        include_hits_dict, exluded_hits_dict = Csb_phylogeny.analyze_unplausible_proteins_in_trees_parallel(
            options.phylogeny_directory,
            unplausible_hits_dict,
            grp1_hits_dict,
            options.pam_long_branch_thres,
            options.pam_phylogenetic_distance,
            options.cores
        )

        # proteinIDs in included hits set are below threshold distance to a plausible grp1 sequence
        myUtil.save_cache(options, cache_name, include_hits_dict)

    # Delete the included hits from the unplausible dict
    for domain, included_proteins in include_hits_dict.items():
        if domain in unplausible_hits_dict:
            unplausible_hits_dict[domain] -= included_proteins
            # Clean up domains with empty sets
            if not unplausible_hits_dict[domain]:
                del unplausible_hits_dict[domain]
    
    # Print out the number of excluded sequences
    print("\nSelection of grp1 training sequence by computed regression of the presence/absence matrix and phylogenetic placement")
    summarize_domain_filtering(grp1_hits_dict, basis_hits_dict, unplausible_hits_dict)

    # Step 5: Copy grp1 files without unplausible proteinIDs to the grp2 files
    print("\nWriting grp0 and included grp1 sequences to grp2 fasta files")
    write_plausible_fasta_sequences(options.fasta_output_directory, unplausible_hits_dict, options.fasta_output_directory)
    
    return



def log_all_mcl_cluster_statistics(reference_dict, cluster_proteins_dict, grouped_3_dict, grouped_4_dict):
    """
    Print selection statistics for all domains based on MCL clustering, including
    final merged set composition per domain.

    Args:
        reference_dict (dict): {domain: set of reference proteinIDs}
        cluster_proteins_dict (dict): {domain: set of all proteinIDs from clusters with refs}
        grouped_3_dict (dict): {domain: set of selected proteinIDs}
        grouped_4_dict (dict): {domain: set of selected proteinIDs}
    """
    for domain in reference_dict:
        reference_sequences = reference_dict.get(domain, set())
        mcl_cluster_protID_set = cluster_proteins_dict.get(domain, set())
        group3_set = grouped_3_dict.get(domain, set())
        group4_set = grouped_4_dict.get(domain, set())
        union_set = group3_set | group4_set

        total_clusters_with_refs = len(mcl_cluster_protID_set)
        total_reference_count = len(reference_sequences)
        total_new_hits = len(mcl_cluster_protID_set - reference_sequences)

        final_merged_set = reference_sequences | group3_set | group4_set
        ref_in_final = len(reference_sequences)
        grp3_in_final = len(group3_set - reference_sequences)
        grp4_in_final = len(group4_set - reference_sequences - group3_set)

        print(f"\n[{domain}]")
        print(f"  Reference sequences:             {total_reference_count}")
        print(f"  Total sequences in ref-clusters: {total_clusters_with_refs}")
        print(f"  New candidate sequences:         {total_new_hits}")
        print(f"  Selected by (CSB):               {len(group3_set)}")
        print(f"  Selected by (PAM):               {len(group4_set)}")
        print(f"  Final set:                       {len(final_merged_set)}")
        print(f"     ├─ from references:           {ref_in_final}")
        print(f"     ├─ from (CSB) only:           {grp3_in_final}")
        print(f"     └─ from (PAM) only:           {grp4_in_final}")




