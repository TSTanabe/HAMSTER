import argparse
import os
import pickle
import sys
import gzip
import shutil
import random

from datetime import datetime
from typing import Optional, Union, List, Dict, Set, Any
from pathlib import Path
from src.core.logging import get_logger

logger = get_logger(__name__)


def print_header(string: str, verbose: int = 0) -> None:
    """Prints a timestamped header if verbose is set."""
    if not verbose:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("\n" + string)
        print(current_time)
        print((len(current_time) + 3 + len(string)) * "-")


def getAllFiles(directory: str, ending: Union[str, int] = 0) -> List[str]:
    """
    Recursively retrieves files from a directory.

    Args:
        directory (str): Root directory to search.
        ending (str or int): Optional file suffix to filter by.

    Returns:
        List[str]: List of matching file paths.
    """
    matched = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            file = os.path.join(path, name)
            if ending == 0 or file.endswith(ending):
                matched.append(file)
    return matched


def get_genome_id(path: str) -> str:
    """Extracts genome ID from filename (before first dot)."""
    return os.path.basename(path).split(".")[0]


def find_executable(executable: str) -> str:
    """
    Locate an executable in PATH OR <project-root>/bin/
    """

    # 1) First check PATH
    path = shutil.which(executable)
    if path:
        return path

    # 2) Determine project root no matter where this is called from
    #    (this file is under HAMSTER/src/... so we climb up until HAMSTER)
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "bin").exists():
            bin_dir = parent / "bin"
            break
    else:
        raise FileNotFoundError("Could not locate project root containing bin/")

    local_path = bin_dir / executable

    if local_path.exists() and os.access(local_path, os.X_OK):
        return str(local_path)

    raise FileNotFoundError(
        f"{executable} executable not found in PATH or {local_path}"
    )


def save_cache(
    options: Any,
    name: str,
    data: object,
    redirect: Optional[str] = None,
    overwrite: bool = False,
) -> None:
    """
    Saves a Python object using pickle into a defined cache directory.

    Args:
        options: Options object with a .result_files_directory attribute.
        name (str): Filename for the pickle file.
        data (object): Any serializable Python object.
        redirect (str, optional): Alternate output directory.
        overwrite (bool): Whether to overwrite an existing file.
    """
    cache_dir = (
        redirect
        if redirect
        else os.path.join(options.result_files_directory, "pkl_cache")
    )
    os.makedirs(cache_dir, exist_ok=True)
    file_path = os.path.join(cache_dir, name)
    if os.path.exists(file_path) and not overwrite:
        logger.debug(f"Cache file exists and overwrite disabled: {file_path}")
        return
    with open(file_path, "wb") as f:
        pickle.dump(data, f)
    logger.debug(f"Saved cache file: {file_path}")


def load_cache(
    options: Any, name: str, file_path: Optional[str] = None
) -> Optional[object]:
    """
    Loads a Python object from a pickle file in the cache directory.

    Args:
        options: Options object with .result_files_directory.
        name (str): Name of the cache file.
        file_path (str, optional): Full path override.

    Returns:
        The loaded Python object or None if not found or on error.
    """
    cache_dir = os.path.join(options.result_files_directory, "pkl_cache")
    file_path = file_path or os.path.join(cache_dir, name)
    if os.path.exists(file_path):
        try:
            logger.debug(f"Loading cache file: {file_path}")
            with open(file_path, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            logger.error(f"Failed to load cache {file_path}: {e}")
            return None
    logger.debug(f"Precomputed cache file not found: {file_path}")
    return None


def merge_grouped_refseq_dicts_simple(
    grouped_3_dict: Dict[str, Set[str]], grouped_4_dict: Dict[str, Set[str]]
) -> Dict[str, Set[str]]:
    """
    Merges two protein dictionaries.

    Input example:
        grouped_3_dict = {"domainA": {"prot1", "prot2"}}
        grouped_4_dict = {"domainA": {"prot3"}}

    Output:
        {"domainA": {"prot1", "prot2", "prot3"}}
    """
    merged = {}
    all_domains = set(grouped_3_dict.keys()) | set(grouped_4_dict.keys())
    for domain in all_domains:
        merged[domain] = grouped_3_dict.get(domain, set()) | grouped_4_dict.get(
            domain, set()
        )
    return merged


def move_HMMs(input_folder, output_folder, file_extension):
    logger.info(f"Saving HMMs in the directory: {output_folder}")
    for datafile in os.listdir(input_folder):
        if datafile.endswith(file_extension):
            source = os.path.join(input_folder, datafile)
            target = os.path.join(output_folder, datafile)
            shutil.move(source, target)
