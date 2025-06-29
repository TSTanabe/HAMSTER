#!/usr/bin/python

import csv
import os
import shutil
import subprocess
from collections import defaultdict
from typing import Dict, Any

from . import Csb_proteins
from . import myUtil

logger = myUtil.logger

##################################################################
#### MCL tribe clustering with diamond blast alrgorithm ##########
##################################################################

    

def run_mmseqs_clustering(
    directory: str,
    min_seq_id: float,
    min_aln_len: float,
    cores: int
) -> Dict[str, str]:
    """
    Runs MMseqs2 'easy-linclust' on each .faa file in the given directory.

    Args:
        directory (str): Path with input .faa files.
        min_seq_id (float): Minimum sequence identity (e.g. 0.95).
        min_aln_len (float): Minimum alignment length (fraction, e.g. 0.8).
        cores (int): Number of threads.

    Returns:
        dict: {input_basename: output_cluster_tsv_path}

    Example:
        run_mmseqs_clustering("/data/fastas", 0.95, 0.8, 16)
    """
    logger.info(f"Clustering at proteins with mmseqs2 linclust") 
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
            logger.debug(f"Linclust results already exist for {output_prefix} - skipping: {output_file_path}")
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
            logger.error(f"MMseqs2 failed on {faa_file} with exit code {exit_code}")
            continue

        # Cluster-Datei registrieren
        if os.path.exists(output_file_path):
            output_files_dict[input_prefix] = output_file_path
        else:
            logger.error(f"Expected cluster file not found: {output_file_path}")

    return output_files_dict


def run_mmseqs_linclust_lowlevel(
    directory: str,
    min_seq_id: float,
    min_aln_len: float,
    cores: int
) -> Dict[str, str]:
    """
    Runs MMseqs2 clustering (low-level) for each .faa file in `directory`,
    producing only a *_cluster.tsv file (other intermediate files are deleted).

    Args:
        directory (str): Path to .faa files.
        min_seq_id (float): Minimum sequence identity (0.0-1.0).
        min_aln_len (float): Minimum alignment length (not used here).
        cores (int): Number of threads.

    Returns:
        dict: {input_basename: path_to_mcl_clusters_txt}

    Example:
        run_mmseqs_linclust_lowlevel("fastas/", 0.95, 0.8, 8)
    """
    logger.debug(f"Clustering sequences for at {float(min_seq_id)*100} % identity with mmseqs linclust")
    logger.debug(f"Skipping existing results")
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
            logger.debug(f"Results already exists: {cluster_tsv}")
            output_files[base] = cluster_tsv
            continue

        logger.debug(f"Clustering sequences for {base} at {float(min_seq_id)*100} % identity")

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
            logger.error(f"MMseqs2 failed for {faa_file}: {e}")
        #finally:
        #    # Cleanup intermediate files
        #    for f in [input_db, f"{input_db}.index", result_db, tmp_dir]:
        #        if os.path.isdir(f):
        #            shutil.rmtree(f, ignore_errors=True)
        #        elif os.path.exists(f):
        #            os.remove(f)

    return output_files

def convert_linclust_to_mcl(
    input_tsv: str,
    output_txt: str
) -> None:
    """
    Converts MMseqs linclust 2-column .tsv (cluster_id, seq_id)
    to MCL format (space-separated clusters, one per line).

    Args:
        input_tsv (str): Path to MMseqs cluster .tsv file.
        output_txt (str): Path for MCL-style .txt output.

    Returns:
        None.

    Example:
        convert_linclust_to_mcl("input.tsv", "clusters.txt")
    """
    
    clusters = defaultdict(set)

    with open(input_tsv, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue  # Skip comment/header lines
            if len(row) < 2:
                logger.warning(f"Malformed line in {input_tsv}: {row}")
                continue
            cluster_id, seq_id = row[:2]
            clusters[cluster_id].add(seq_id)

    with open(output_txt, 'w') as f:
        for members in clusters.values():
            f.write(' '.join(sorted(members)) + '\n')





