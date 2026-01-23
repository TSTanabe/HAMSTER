#!/usr/bin/env python3
"""
Definition des Options-Objekts und Hilfsfunktionen zur Validierung.
"""

import os
import sys

from src.core.logging import get_logger

logger = get_logger(__name__)

# Zentrale Basis-Location, identisch zu deiner bisherigen Logik
if getattr(sys, "frozen", False):
    BASE_DIR = os.path.split(sys.executable)[0]
else:
    # Pfad dieses Files (/HAMSTER/src/cli/options.py)
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

    # Projekt-Root ist zwei Ebenen oberhalb: /HAMSTER
    BASE_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", ".."))


class Options:
    """
    Holds all pipeline options and state variables.
    All attributes are set by CLI arguments or at runtime.
    """

    def __init__(self):
        # variables used by other subroutines (wie in deiner bisherigen Klasse)
        # in order of appearance

        # Queued genomes and file mappings
        self.queued_genomes = set()
        self.faa_files = {}
        self.gff_files = {}

        # csb prediction dereplication
        self.redundant = 0
        self.non_redundant = 0
        self.redundancy_hash = dict()

        # validation
        self.sequence_faa_file = (
            None  # dictionary to the target files for the validation
        )
        self.reports = dict()

        # csb naming
        self.csb_name_prefix = "csb-"
        self.csb_name_suffix = "_"

        # self query / concatenation flags
        self.self_query = None
        self.self_seqs = None
        self.deconcat_flag = 0
        self.glob_flag = 0

        # sqlite
        self.sqlite_chunks = 999  # chunks size for placeholders in sqlite fetch

        # plotting and script paths
        self.plotting_Rscripts = os.path.join(BASE_DIR, "src")

        # misc
        self.hardcap = 10000  # Allowed number of seqs to exceed the hardcap
        self.new_project = False  # make a new project

        # speichere Location auch im Objekt
        self.location = BASE_DIR

        # Platzhalter für andere Felder, die später via argparse gesetzt werden
        # (z.B. self.fasta_file_directory, self.result_files_directory, ...)


def validate_options(options: Options) -> None:
    """
    Validates input arguments (directories/files).

    Args:
        options (Options): The options object.

    Raises:
        SystemExit if no valid input provided.
    """
    # Check if a database file exists
    database_provided = bool(
        getattr(options, "database_directory", None)
        and os.path.isfile(options.database_directory)
    )

    # Check if a result directory exists
    result_directory_provided = bool(
        getattr(options, "result_files_directory", None)
        and os.path.isdir(options.result_files_directory)
    )

    # Check if both fasta file directory and query file exist
    fasta_and_query_provided = bool(
        getattr(options, "fasta_file_directory", None)
        and getattr(options, "query_file", None)
        and os.path.isdir(options.fasta_file_directory)
        and os.path.isfile(options.query_file)
    )

    # Only proceed if one of the conditions is met
    if not (database_provided or result_directory_provided or fasta_and_query_provided):
        logger.error(
            "Please provide either a fasta file directory and a query file, "
            "or an existing result directory from a previous run."
        )
        sys.exit(1)
