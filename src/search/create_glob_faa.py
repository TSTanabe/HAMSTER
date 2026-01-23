import os
import gzip
from typing import TextIO

def open_fasta_maybe_gz(path: str) -> TextIO:
    """
    Open a FASTA file that may be plain text or gzip-compressed.
    Always returns a text-mode file handle.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")  # text mode
    else:
        return open(path, "r")


def create_glob_file(options) -> None:
    """
    Concatenate all queued genomes' .faa files into one glob.faa with unified headers.
    Args:
        options: Contains .queued_genomes, .faa_files, .result_files_directory, .glob_faa, .fasta_file_directory
    Output:
        options.glob_faa is written (and registered in options)
    """
    # Check if a glob fasta file and deconcatenated files were already provided
    if (
        not options.glob_faa is None
        and os.path.isfile(options.glob_faa)
        and os.path.isdir(options.fasta_file_directory)
    ):
        return

    # Define the output file path for the glob.faa file
    options.glob_faa = os.path.join(options.result_files_directory, "glob.faa")

    with open(options.glob_faa, "w") as outfile:
        for genomeID in options.queued_genomes:
            faa_file = options.faa_files[
                genomeID
            ]  # Get the FASTA file for the current genomeID

            # Open each faa_file explicitly and close it after processing
            infile = open_fasta_maybe_gz(faa_file)
            try:
                sequence = ""
                header = None

                for line in infile:
                    line = line.strip()
                    if line.startswith(">"):
                        # Write the previous sequence if it exists
                        if header and sequence:
                            outfile.write(f"{header}\n{sequence}\n")

                        # Start a new record
                        original_id = line[1:].split()[
                            0
                        ]  # Extract original ID up to the first whitespace
                        header = (
                            f">{genomeID}___{original_id}"  # Create the modified header
                        )
                        sequence = ""  # Reset sequence for the new record
                    else:
                        sequence += line  # Append sequence lines

                # Write the last record if there was one
                if header and sequence:
                    outfile.write(f"{header}\n{sequence}\n")
            finally:
                infile.close()  # Ensure the file is closed after reading
