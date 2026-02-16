import gzip
import os
import shutil
import tempfile
import math
import csv
import traceback

from multiprocessing import Pool
from typing import Dict, Iterable, Tuple, Optional

from src.core.logging import get_logger
from src.db import database
from src.fasta_preparation.materialize import materialize_pair_gz_next_to_input
from src.parse_reports import parse_reports, csb_finder
from src.search import diamond_search, blastp_n_filter, query_selfblast_search

logger = get_logger(__name__)

# --- Global worker state (nur was wirklich gebraucht wird) ---
global_query_file: Optional[str] = None
global_evalue: Optional[float] = None
global_thrs_score: Optional[float] = None
global_coverage_cutoff: Optional[float] = None
global_identity_cutoff: Optional[float] = None
global_hits_limit: Optional[int] = None

global_query_length_dict: Optional[Dict[str, float]] = None
global_selfblast_scores_dict: Optional[Dict[str, float]] = None

global_faa_files = None
global_gff_files = None
global_nucleotide_range: Optional[int] = None



def init_worker(
    query_file: str,
    evalue_cutoff: float,
    score_cutoff: float,
    coverage_cutoff: float,
    identity_cutoff: float,
    hits_limit: int,
    query_length_dict: Dict[str, float],
    selfblast_scores_dict: Dict[str, float],
    faa_files,
    gff_files,
    nucleotide_range: int,
):
    global global_query_file, global_evalue, global_thrs_score
    global global_coverage_cutoff, global_identity_cutoff, global_hits_limit
    global global_query_length_dict, global_selfblast_scores_dict
    global global_faa_files, global_gff_files, global_nucleotide_range

    global_query_file = query_file
    global_evalue = evalue_cutoff
    global_thrs_score = score_cutoff
    global_coverage_cutoff = coverage_cutoff
    global_identity_cutoff = identity_cutoff
    global_hits_limit = hits_limit

    global_query_length_dict = query_length_dict
    global_selfblast_scores_dict = selfblast_scores_dict

    global_faa_files = faa_files
    global_gff_files = gff_files
    global_nucleotide_range = nucleotide_range



def process_single_genome(genomeID: str):
    """
    Worker: for one genomeID
      - materialize_pair_gz_next_to_input gives us plain faa + gff
      - ensure shared dmnd DB exists
      - diamond blastp -> raw tab
      - filter + BSR -> filtered tab
      - parse -> protein_dict
      - add GFF coords + sequences
      - CSB find -> cluster_dict
      - delete raw + filtered outputs
    Returns: (genomeID, protein_dict, cluster_dict)
    """
    try:
        faa_in = global_faa_files[genomeID]
        gff_in = global_gff_files[genomeID]

        # dein bestehender Contextmanager (entpackt / materialisiert bei Bedarf)
        with materialize_pair_gz_next_to_input(faa_in, gff_in) as (faa_file, gff_file):

            # 1) DIAMOND blastp
            # Ideal: diamond_search.diamond_search kann DB-Pfad akzeptieren.
            # Falls nicht, musst du diamond_search minimal erweitern.
            report = diamond_search.diamond_search(
                path=faa_file,
                query_fasta=global_query_file,
                cores=1,
                evalue=global_evalue,
                diamond_report_hits_limit=global_hits_limit,
            )

            # Parse and filter results
            protein_dict = parse_reports.parse_filter_single_blastreport(
                genome_id=genomeID,
                filepath=report,
                evalue_cutoff=global_evalue,
                score_cutoff=global_thrs_score,
                coverage_cutoff=global_coverage_cutoff,
                identity_cutoff=global_identity_cutoff,
                sequence_lengths=global_query_length_dict,
                selfblast_scores=global_selfblast_scores_dict,
            )

            # Get genomic features and protein sequences
            parse_reports.parse_gff_file(gff_file, protein_dict)
            parse_reports.get_protein_sequence(faa_file, protein_dict)

            # Find syntenic blocks
            cluster_dict = csb_finder.find_syntenic_blocks(
                genomeID, protein_dict, global_nucleotide_range
            )

            return genomeID, protein_dict, cluster_dict

    except Exception as e:
        logger.warning(f"Skipped {faa_file} due to an error - {e}")
        logger.debug(traceback.format_exc())
        return None


def _flush_batches_to_db_and_csb(
    config,
    *,
    genome_id_batch: set[str],
    protein_batch: dict,
    cluster_batch: dict,
):
    """
    Flush rule:
      1) INSERT genome IDs FIRST (as requested)
      2) INSERT proteins/clusters
      3) APPEND cluster types to gene_clusters_file
      4) CLEAR batches
    """
    if not genome_id_batch and not protein_batch and not cluster_batch:
        return

    # 1) genome IDs FIRST
    if genome_id_batch:
        database.insert_database_genome_ids(config.database_directory, genome_id_batch)

    # 2) proteins/clusters
    if protein_batch:
        database.insert_database_proteins(config.database_directory, protein_batch)

    if cluster_batch:
        database.insert_database_clusters(config.database_directory, cluster_batch)

        # 3) CSB output append (write BEFORE clear)
        with open(config.gene_clusters_file, "a") as f:
            for cluster_id, cluster in cluster_batch.items():
                f.write(cluster_id + "\t" + "\t".join(cluster.types) + "\n")

    # 4) clear
    genome_id_batch.clear()
    protein_batch.clear()
    cluster_batch.clear()

def run_consecutive_parallel_search(
    config,
):
    """
    Main-process routine (consecutive file search):
      - calls per-genome worker (DIAMOND -> filter -> parse -> CSB)
      - inserts genome IDs BEFORE each batch flush
      - performs batched inserts (proteins + clusters) and writes CSB output file
      - logs progress in ~1% steps
    """
    logger.info("Initilize DIAMOND BLASTp self-blast against query")
    selfblast_scores_dict, query_length_dict = query_selfblast_search.self_blast_query(
        config
    )

    # -------- genome list (FAA ∩ GFF) --------
    genome_ids = set(config.faa_files.keys()) & set(config.gff_files.keys())
    total = len(genome_ids)
    if total == 0:
        raise ValueError("No genomes to process: faa_files ∩ gff_files is empty.")

    # -------- worker outdir (per-genome tmp outputs) --------
    outdir = os.path.join(os.path.dirname(config.filtered_blast_table), "diamond_per_genome")
    os.makedirs(outdir, exist_ok=True)

    # -------- reset CSB file (optional; keep if you want fresh run) --------
    if os.path.exists(config.gene_clusters_file):
        os.remove(config.gene_clusters_file)

    # -------- batching --------
    batch_size = getattr(config, "glob_chunks", 50)  # interpreted as "hit genomes per flush"
    protein_batch: dict = {}
    cluster_batch: dict = {}

    # IMPORTANT: genome IDs must be inserted BEFORE flush
    genome_id_batch: set[str] = set()

    # -------- progress logging (1% steps) --------
    step_every = max(1, math.ceil(total / 100))
    next_step = step_every
    done = 0

    # -------- pool size --------
    worker_processes = max(1, min(getattr(config, "cores", 1), total))
    if worker_processes > 1:
        worker_processes = max(1, min(worker_processes - 1, total))  # keep 1 core for main/DB I/O

    logger.info(f"Starting consecutive parallel search: {total} genomes, workers={worker_processes}")

    # -------- Pool init args (EXPLICIT MAPPING) --------
    # init_worker signature (as previously defined):
    #
    # init_worker(
    #   query_file, evalue, searchcoverage, minseqid, hits_limit, alignment_mode, outdir,
    #   query_length_dict, selfblast_scores_dict,
    #   score_threshold_diction, faa_files, gff_files, multidomain_allowed, glob_gff,
    #   nucleotide_range, thrs_score
    # )
    #
    initargs = (
        config.query_file,  # query_file
        config.evalue,  # evalue_cutoff
        config.thrs_score,  # score_cutoff
        config.searchcoverage,  # coverage_cutoff
        config.minseqid,  # identity_cutoff
        config.diamond_report_hits_limit,  # hits_limit
        query_length_dict,  # query_length_dict
        selfblast_scores_dict,  # selfblast_scores_dict
        config.faa_files,  # faa_files
        config.gff_files,  # gff_files
        config.nucleotide_range,  # nucleotide_range
    )

    with Pool(
        processes=worker_processes,
        initializer=init_worker,
        initargs=initargs,
    ) as pool:

        for res in pool.imap_unordered(process_single_genome, genome_ids, chunksize=1):
            done += 1
            if done >= next_step or done == total:
                pct = int((done / total) * 100)
                logger.info(f"Progress: {done}/{total} genomes ({pct}%)")
                next_step += step_every

            if res is None:
                continue

            genome_id, protein_dict, cluster_dict = res
            if not protein_dict:
                continue

            # accumulate batches
            genome_id_batch.add(genome_id)
            protein_batch.update(protein_dict)
            cluster_batch.update(cluster_dict)

            # flush trigger: N hit-genomes accumulated
            if len(genome_id_batch) >= batch_size:
                _flush_batches_to_db_and_csb(
                    config,
                    genome_id_batch=genome_id_batch,
                    protein_batch=protein_batch,
                    cluster_batch=cluster_batch,
                )

    # final flush
    _flush_batches_to_db_and_csb(
        config,
        genome_id_batch=genome_id_batch,
        protein_batch=protein_batch,
        cluster_batch=cluster_batch,
    )

    logger.info(
        f"Consecutive parallel search finished. "
        f"Processed total={total}, successfully={len(genome_ids)}, "
        f"failed={len(genome_id_batch)}"
    )


