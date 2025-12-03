#!/usr/bin/env python3
"""
CLI-Definition: argparse, Parser-Aufbau und Argument-Parsing in ein Options-Objekt.
"""

from __future__ import annotations

import argparse
import sys

from .options import Options, validate_options, BASE_DIR


def _add_all_groups(parser: argparse.ArgumentParser, show_advanced: bool) -> None:
    """
    Definiert alle Argumentgruppen (essential + advanced) auf dem Parser.
    Diese Logik ist weitgehend 1:1 aus deiner bisherigen parse_arguments-Funktion
    übernommen, nur in eine Hilfsfunktion ausgelagert.
    """

    # --- Essential arguments ---
    essential = parser.add_argument_group("Essential parameters")
    essential.add_argument(
        "-f",
        dest="fasta_file_directory",
        type=dir_path,
        default=BASE_DIR,
        metavar="<directory>",
        help=(
            "Directory containing the target FASTA files to analyze (used with -q). "
            "Example: ./genomes/"
        ),
    )
    essential.add_argument(
        "-q",
        dest="query_file",
        type=file_path,
        default=None,
        metavar="<filepath>",
        help=(
            "FASTA file containing query protein or gene sequences (used with -f). "
            "Example: ./query.faa"
        ),
    )
    essential.add_argument(
        "-r",
        dest="result_files_directory",
        type=dir_path,
        default=os.path.join(BASE_DIR, "results"),
        metavar="<directory>",
        help=("Directory for output result files."),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        type=int,
        default=1,
        choices=[0, 1, 2],
        help="Set logging level: 0=WARNING, 1=INFO, 2=DEBUG",
    )
    parser.add_argument(
        "--help-all",
        action="store_true",
        help="Show all available options, including advanced parameters, and exit.",
    )

    # --- Advanced arguments (nur Helptext, wenn show_advanced True) ---
    def maybe(text: str) -> str:
        return text if show_advanced else argparse.SUPPRESS

    resources = parser.add_argument_group("Advanced parameters")
    resources.add_argument(
        "-s",
        dest="stage",
        type=int,
        default=0,
        choices=range(0, 11),
        metavar="<int>",
        help=maybe(
            "Pipeline stage to start execution from (0=full run, 1-10=partial)."
        ),
    )
    resources.add_argument(
        "-x",
        dest="end",
        type=int,
        default=10,
        choices=range(0, 11),
        metavar="<int>",
        help=maybe("Pipeline stage to stop after completion (inclusive)."),
    )
    resources.add_argument(
        "-c",
        dest="cores",
        type=int,
        default=2,
        metavar="<int>",
        help=maybe("Number of threads to use for parallel processing."),
    )
    resources.add_argument(
        "-t",
        dest="taxonomy_file",
        metavar="<filepath>",
        default=None,
        help=maybe("Path to taxonomy CSV file for input assemblies."),
    )
    resources.add_argument(
        "-db",
        dest="database_directory",
        metavar="<filepath>",
        default=None,
        help=maybe("Path to existing HAMSTER SQLite database."),
    )
    resources.add_argument(
        "-cv_off",
        dest="cross_validation_deactivated",
        action="store_true",
        help=maybe(
            "Disable cross-validation (recommended for faster validation only)."
        ),
    )

    resources2 = parser.add_argument_group("GlobDB file parameters (advanced)")
    resources2.add_argument(
        "-glob_faa",
        dest="glob_faa",
        metavar="<filepath>",
        default=None,
        help=maybe("Concatenated FASTA file with all input assemblies for speedup."),
    )
    resources2.add_argument(
        "-glob_gff",
        dest="glob_gff",
        metavar="<filepath>",
        default=None,
        help=maybe("Concatenated GFF annotation file for all assemblies."),
    )
    resources2.add_argument(
        "-glob_blast_table",
        dest="glob_table",
        metavar="<filepath>",
        default=None,
        help=maybe("Precomputed multi-assembly BLAST tabular result file."),
    )
    resources2.add_argument(
        "-glob_chunks",
        dest="glob_chunks",
        type=int,
        default=3000,
        metavar="<int>",
        help=maybe("Chunk size for batch parsing of large files."),
    )
    resources2.add_argument(
        "-glob_off",
        dest="glob_search",
        action="store_false",
        help=maybe(
            "Disable creation of concatenated files for the initial search step."
        ),
    )

    search = parser.add_argument_group("DIAMOND blastp search parameters (advanced)")
    search.add_argument(
        "-evalue",
        dest="evalue",
        type=float,
        default=1e-5,
        metavar="<float>",
        help=maybe("Maximum allowed E-value for BLASTp matches."),
    )
    search.add_argument(
        "-thrs_score",
        dest="thrs_score",
        type=int,
        default=100,
        metavar="<int>",
        help=maybe("Minimum required BLASTp alignment score."),
    )
    search.add_argument(
        "-min_seq_id",
        dest="minseqid",
        type=float,
        default=25,
        metavar="<float>",
        help=maybe("Minimum percentage sequence identity [%%] for BLASTp results."),
    )
    search.add_argument(
        "-search_coverage",
        dest="searchcoverage",
        type=float,
        default=0.6,
        metavar="<float>",
        help=maybe("Minimum fraction of query/target aligned (0.0–1.0)."),
    )
    search.add_argument(
        "-alignment_mode",
        dest="alignment_mode",
        type=int,
        default=2,
        choices=[0, 1, 2],
        metavar="<int>",
        help=maybe("DIAMOND BLASTp alignment mode."),
    )
    search.add_argument(
        "-blast_score_ratio",
        dest="thrs_bsr",
        type=float,
        default=0.0,
        metavar="<float>",
        help=maybe("Minimum BLAST score ratio threshold for hit inclusion."),
    )
    search.add_argument(
        "-allow_multidomain",
        dest="multidomain_allowed",
        action="store_true",
        help=maybe("Permit hits to multiple query domains per sequence."),
    )
    search.add_argument(
        "-reports_hit",
        dest="diamond_report_hits_limit",
        type=int,
        default=0,
        metavar="<int>",
        help=maybe("Limit number of BLASTp hits reported per query (0=all)."),
    )

    genecluster = parser.add_argument_group(
        "Gene cluster prediction parameters (advanced)"
    )
    genecluster.add_argument(
        "-distance",
        dest="nucleotide_range",
        type=int,
        default=3500,
        metavar="<int>",
        help=maybe("Maximum nucleotide distance between syntenic genes in a cluster."),
    )
    genecluster.add_argument(
        "-p",
        dest="patterns_file",
        metavar="<filepath>",
        default=os.path.join(BASE_DIR, "src", "Patterns"),
        help=maybe(
            "Tab-separated file with predefined syntenic gene cluster patterns."
        ),
    )

    csb = parser.add_argument_group(
        "Collinear syntenic block (csb) parameters (advanced)"
    )
    csb.add_argument(
        "-insertions",
        dest="insertions",
        type=int,
        default=2,
        metavar="<int>",
        help=maybe("Maximum insertions allowed between genes in a csb."),
    )
    csb.add_argument(
        "-occurence",
        dest="occurence",
        type=int,
        default=1,
        metavar="<int>",
        help=maybe("Minimum number of csb occurrences to be recognized."),
    )
    csb.add_argument(
        "-min_csb_size",
        dest="min_csb_size",
        type=int,
        default=3,
        metavar="<int>",
        help=maybe("Minimum number of genes in a csb."),
    )
    csb.add_argument(
        "-max_csb_size",
        dest="max_csb_size",
        type=int,
        default=40,
        metavar="<int>",
        help=maybe("Maximum number of genes in a csb."),
    )
    csb.add_argument(
        "-max_gene_repeats",
        dest="max_domain_repeats",
        type=int,
        default=2,
        metavar="<int>",
        help=maybe("Maximum number of repeated genes in a csb."),
    )
    csb.add_argument(
        "-jaccard",
        dest="jaccard",
        type=float,
        default=0.0,
        metavar="<float>",
        help=maybe("Maximum Jaccard dissimilarity allowed for csb clustering."),
    )

    csb_selection = parser.add_argument_group(
        "Base training data csb and protein selection parameters (advanced)"
    )
    csb_selection.add_argument(
        "-exclude_csb_score",
        dest="low_hitscore_csb_cutoff",
        type=float,
        default=0.8,
        metavar="<float>",
        help=maybe("Exclude csb with all hits below this BLAST score ratio."),
    )
    csb_selection.add_argument(
        "-exclude_csb_protein",
        dest="exclude_csb_proteins",
        nargs="+",
        default=[],
        metavar="<list>",
        help=maybe("Suppress csb with hits for these proteins."),
    )
    csb_selection.add_argument(
        "-redo_base_selection",
        dest="redo_base_selection",
        action="store_true",
        default=False,
        help=maybe("Redo all selection steps and overwrite all caches."),
    )

    pam_search = parser.add_argument_group(
        "Presence/absence matrix (pam) parameters (advanced)"
    )
    pam_search.add_argument(
        "-mx_thrs",
        dest="pam_threshold",
        type=float,
        default=0.3,
        metavar="<float>",
        help=maybe("Significance threshold for presence/absence matrix co-occurrence."),
    )
    pam_search.add_argument(
        "-mx_bsr",
        dest="pam_bsr_threshold",
        type=float,
        default=0.6,
        metavar="<float>",
        help=maybe("Minimum BLAST score ratio for pam prediction inclusion."),
    )

    mcl_search = parser.add_argument_group(
        "Protein sequence clustering parameters (advanced)"
    )
    mcl_search.add_argument(
        "-mcl_min_seq_id",
        dest="mcl_min_seq_id",
        type=float,
        default=0.4,
        metavar="<float>",
        help=maybe("Minimal sequence identity for clustering (0.0–1.0)."),
    )
    mcl_search.add_argument(
        "-mcl_density_thrs",
        dest="mcl_density_thrs",
        type=auto_float,
        default="auto",
        metavar="<float>",
        help=maybe("Required fraction of reference sequences in a cluster (0.0-1.0)."),
    )
    mcl_search.add_argument(
        "-mcl_reference_thrs",
        dest="mcl_reference_thrs",
        type=auto_float,
        default="auto",
        metavar="<float>",
        help=maybe(
            "Required fraction of all reference sequences found in a cluster (0.0-1.0)."
        ),
    )

    alignment = parser.add_argument_group("Alignment parameters (advanced)")
    alignment.add_argument(
        "-min_seqs",
        dest="min_seqs",
        type=int,
        default=3,
        metavar="<int>",
        help=maybe("Minimum number of sequences required for alignment."),
    )
    alignment.add_argument(
        "-max_seqs",
        dest="max_seqs",
        type=int,
        default=100000,
        metavar="<int>",
        help=maybe("Maximum number of sequences to align."),
    )
    alignment.add_argument(
        "-gap_col_remove",
        dest="gap_remove_threshold",
        type=float,
        default=0.05,
        metavar="<float>",
        help=maybe(
            "Remove alignment columns with gaps above this threshold (0.0–1.0)."
        ),
    )
    alignment.add_argument(
        "-include_domains",
        dest="include_list",
        nargs="+",
        default=[],
        metavar="<list>",
        help=maybe("Domains to specifically include (space-separated)."),
    )
    alignment.add_argument(
        "-exclude_domains",
        dest="exclude_list",
        nargs="+",
        default=[],
        metavar="<list>",
        help=maybe("Domains to specifically exclude (space-separated)."),
    )


def build_parser(show_advanced: bool) -> argparse.ArgumentParser:
    """
    Erstellt und konfiguriert den ArgumentParser.
    """
    parser = argparse.ArgumentParser(
        description="HAMSTER: Homolog and Synteny Mining Pipeline",
        epilog="Please cite:",
        usage=(
            "hamster.py -f <genomes_dir> -q <query.faa> [options]\n"
            "       hamster.py -r <results_dir> [options]\n"
            "       hamster.py --help-all"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_all_groups(parser, show_advanced=show_advanced)
    return parser


def parse_arguments(argv: list[str] | None = None) -> Options:
    """
    Parses CLI arguments into an Options object.

    - Wenn '--help-all' dabei ist, wird die komplette Hilfe gezeigt und beendet.
    - Wenn keine Argumente übergeben werden, wird die Standardhilfe gezeigt
      und mit Fehlermeldung beendet.

    Args:
        argv: Liste der CLI-Argumente (ohne Programmname).
              Wenn None, wird sys.argv[1:] verwendet.

    Returns:
        options: Fertig befülltes Options-Objekt.
    """
    if argv is None:
        argv = sys.argv[1:]

    show_all = "--help-all" in argv
    # '--help-all' nicht an argparse weiterreichen, da wir es selbst behandeln
    clean_args = [a for a in argv if a != "--help-all"]

    parser = build_parser(show_advanced=show_all)

    # --help-all: komplette Hilfe ausgeben und beenden
    if show_all:
        parser.print_help()
        sys.exit(0)

    # Keine Argumente -> kurze Hilfe + Hinweis
    if len(clean_args) == 0:
        parser.print_help()
        sys.exit("Please provide arguments. For help use -h or --help-all")

    # Options-Objekt instanziieren und Namespace direkt hineinparsen
    options = Options()
    parser.parse_args(clean_args, namespace=options)

    # Default: BASE_DIR als location setzen, falls noch nicht überschrieben
    options.location = getattr(options, "location", BASE_DIR)

    # new_project: wie bisher – Standard-Result-Dir => neues Projekt
    default_results_dir = os.path.join(BASE_DIR, "results")
    if getattr(options, "result_files_directory", None) == default_results_dir:
        options.new_project = True

    # 'auto'-Werte nachträglich anpassen (wie in deiner ursprünglichen Funktion)
    if getattr(options, "mcl_density_thrs", None) == "auto":
        options.mcl_density_thrs = None
    if getattr(options, "mcl_reference_thrs", None) == "auto":
        options.mcl_reference_thrs = None

    # Eingaben validieren
    validate_options(options)

    return options


import argparse
import os


def dir_path(path: str) -> str:
    """argparse type: existing directory."""
    if os.path.isdir(path):
        return path
    raise argparse.ArgumentTypeError(f"Directory does not exist: {path!r}")


def file_path(path: str) -> str:
    """argparse type: existing file."""
    if os.path.isfile(path):
        return path
    raise argparse.ArgumentTypeError(f"File does not exist: {path!r}")


def auto_float(value: str):
    """argparse type: float or 'auto'."""
    if value == "auto":
        return "auto"
    try:
        return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Value must be 'auto' or a float, got {value!r}"
        )
