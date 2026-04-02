#!/usr/bin/python
import os
from src.db import database
from src.search import (
    global_file_search,
    create_glob_faa,
    query_selfblast_search,
    blastp_parallel,
)
from src.core.logging import get_logger

logger = get_logger(__name__)


def initial_search(config) -> None:
    """
    Executes the initial protein search using DIAMOND or BLAST.

    Args:
        config (Options): The main pipeline options

    Input Example:
        options.database_directory: str = "./results/my_db.sqlite"
        options.glob_search: bool = True
        options.glob_faa: str = "./output/all_genomes.faa"

    Output:
        Populates result tables in options.database_directory.
        Deletes deconcatenated files as needed.
    """

    # Writes a database with the protein hits and their gene clusters for later use.
    if not os.path.isfile(config.database_directory):
        logger.info("Created database")
        database.create_database(config.database_directory)

    config.query_file_original = config.query_file

    # Durchnummerierte Query erzeugen und als aktive Query verwenden
    config.query_file = query_selfblast_search.make_numbered_query_fasta(
        config.query_file_original,
        config.result_files_directory,
    )

    # header in glob file are genomeID___proteinID

    # global_file_search.initial_glob_search(config)

    blastp_parallel.run_consecutive_parallel_search(config)
    # else:
    #    # Diamond blastp all files separately
    #    parallel_search.initial_genomize_search(options)
