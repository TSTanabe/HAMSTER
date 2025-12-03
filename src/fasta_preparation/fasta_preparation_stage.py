#!/usr/bin/python
from __future__ import annotations


from src.fasta_preparation.translation import (
    parallel_translation,
    parallel_transcription,
)
from src.core.queue import queue_fna_inputs, queue_faa_without_gff

from src.core.logging import get_logger

logger = get_logger(__name__)

"""
Stage: FASTA preparation.

Translates `.fna` files to `.faa`/`.gff` pairs using Prodigal
and transcribes additional FAA files if no GFF is present.
"""


def fasta_preparation(config) -> None:
    """Prepare FASTA files for the pipeline.

    - Collect `.fna` files and translate them to FAA/GFF using Prodigal.
    - Collect FAA files lacking GFF and transcribe them in parallel.

    Args:
        config: Configuration with `.cores` and input directories.

    Side Effects:
        Creates `.faa` and `.gff` files as required.
    """
    # print_header("FASTA preparation (gene calling / translation)", logger=logger)

    fna_files: dict[str, str] = queue_fna_inputs(config)
    parallel_translation(fna_files, config.cores)

    faa_files_without_gff: dict[str, str] = queue_faa_without_gff(config)
    parallel_transcription(faa_files_without_gff, config.cores)

    return
