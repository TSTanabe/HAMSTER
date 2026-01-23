#!/usr/bin/python
import os
from src.db import database
from src.search import global_file_search, create_glob_faa, query_selfblast_search


def initial_search(options) -> None:
    """
    Executes the initial protein search using DIAMOND or BLAST.

    Args:
        options (Options): The main pipeline options

    Input Example:
        options.database_directory: str = "./results/my_db.sqlite"
        options.glob_search: bool = True
        options.glob_faa: str = "./output/all_genomes.faa"

    Output:
        Populates result tables in options.database_directory.
        Deletes deconcatenated files as needed.
    """

    # Writes a database with the protein hits and their gene clusters for later use.
    if not os.path.isfile(options.database_directory):
        database.create_database(options.database_directory)

    options.query_file_original = options.query_file

    # Durchnummerierte Query erzeugen und als aktive Query verwenden
    options.query_file = query_selfblast_search.make_numbered_query_fasta(
        options.query_file_original,
        options.result_files_directory,
    )

    # header in glob file are genomeID___proteinID
    create_glob_faa.create_glob_file(
        options
    )  # TODO make sure only single genomes can be provided
    global_file_search.initial_glob_search(options)
    os.remove(options.glob_faa)
    # else:
    #    # Diamond blastp all files separately
    #    parallel_search.initial_genomize_search(options)
