#!/usr/bin/python

import csv
import os
import numpy as np
from . import Csb_proteins
from . import myUtil
import csv
from collections import defaultdict

import shutil
import subprocess

##################################################################
#### MCL tribe clustering with diamond blast alrgorithm ##########
##################################################################

    

def run_mmseqs_clustering(directory, min_seq_id, min_aln_len, cores): 
    faa_files = [f for f in os.listdir(directory) if f.endswith(".faa")]
    mmseqs = myUtil.find_executable("mmseqs")
    output_files_dict = {}

    for faa_file in faa_files:
        faa_path = os.path.join(directory, faa_file)
        input_prefix = os.path.splitext(faa_file)[0]

        # Output-Dateien explizit ins gleiche Verzeichnis legen
        output_prefix = os.path.join(directory, input_prefix)
        output_file_path = f"{output_prefix}_cluster.tsv"
        tmp_dir = os.path.join(directory, "tmp")

        if os.path.exists(output_file_path):
            print(f"[SKIP] Linclust results already exist: {output_file_path}")
            output_files_dict[input_prefix] = output_file_path
            continue

        # MMseqs2 aufrufen – nur cluster.tsv interessiert uns
        command = (
            f'{mmseqs} easy-linclust "{faa_path}" "{output_prefix}" "{tmp_dir}" '
            f'--min-seq-id {min_seq_id} --threads {cores} '
            f'1>/dev/null 2>/dev/null'
        )
        print(command)
        exit_code = os.system(command)

        if exit_code != 0:
            print(f"[ERROR] MMseqs2 failed on {faa_file} with exit code {exit_code}")
            continue

        # Cluster-Datei registrieren
        if os.path.exists(output_file_path):
            output_files_dict[input_prefix] = output_file_path
        else:
            print(f"[ERROR] Expected cluster file not found: {output_file_path}")

    return output_files_dict


def run_mmseqs_linclust_lowlevel(directory, min_seq_id, min_aln_len, cores):
    """
    Runs MMseqs2 clustering (low-level) for each .faa file in `directory`,
    producing only a *_cluster.tsv file (other intermediate files are deleted).

    Args:
        directory (str): Path to the directory with input .faa files.
        min_seq_id (float): Minimum sequence identity.
        min_aln_len (float): Minimum alignment length (fraction, e.g. 0.8).
        cores (int): Number of threads to use.

    Returns:
        dict: { input_basename (str) : path_to_cluster_tsv (str) }
    """
    faa_files = [f for f in os.listdir(directory) if f.endswith(".faa")]
    mmseqs = myUtil.find_executable("mmseqs")
    output_files = {}

    for faa_file in faa_files:
        faa_path = os.path.join(directory, faa_file)
        base = os.path.splitext(faa_file)[0]
        
        input_db     = os.path.join(directory, f"{base}_db")
        result_db    = os.path.join(directory, f"{base}_result")
        lincluster_tsv  = os.path.join(directory, f"{base}_lincluster.tsv")
        cluster_tsv  = os.path.join(directory, f"{base}_mcl_clusters.txt")
        tmp_dir      = os.path.join(directory, f"{base}_tmp")

        if os.path.exists(cluster_tsv):
            print(f"[SKIP] Already exists: {cluster_tsv}")
            output_files[base] = cluster_tsv
            continue

        print(f"[INFO] Clustering sequences for {base} at {float(min_seq_id)*100} % identity")

        try:
            os.makedirs(tmp_dir, exist_ok=True)

            # Step 1: createdb
            subprocess.run(
                [mmseqs, "createdb", faa_path, input_db],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            # Step 2: linclust
            subprocess.run(
                [
                    mmseqs, "linclust", input_db, result_db, tmp_dir,
                    "--min-seq-id", str(min_seq_id),
                    "--threads", str(cores)
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            # Step 3: createtsv
            subprocess.run(
                [mmseqs, "createtsv", input_db, input_db, result_db, lincluster_tsv],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            
            # Step 4: reformat cluster file to mcl format
            convert_linclust_to_mcl(lincluster_tsv, cluster_tsv)
            
            output_files[base] = cluster_tsv

        except subprocess.CalledProcessError as e:
            print(f"[ERROR] MMseqs2 failed for {faa_file}: {e}")
        #finally:
        #    # Cleanup intermediate files
        #    for f in [input_db, f"{input_db}.index", result_db, tmp_dir]:
        #        if os.path.isdir(f):
        #            shutil.rmtree(f, ignore_errors=True)
        #        elif os.path.exists(f):
        #            os.remove(f)

    return output_files

def convert_linclust_to_mcl(input_tsv, output_txt):
    """
    Converts a MMseqs linclust 2-column .tsv file (cluster_id, seq_id)
    to MCL-style format where each line is a space-separated cluster.
    """
    clusters = defaultdict(set)

    with open(input_tsv, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue  # Skip comment/header lines
            if len(row) < 2:
                print(f"[WARN] Malformed line (less than 2 columns): {row}")
                continue
            cluster_id, seq_id = row[:2]
            clusters[cluster_id].add(seq_id)

    with open(output_txt, 'w') as f:
        for members in clusters.values():
            f.write(' '.join(sorted(members)) + '\n')





