#!/usr/bin/python

import os
import pickle
import sys
import gzip
import shutil
import random
from datetime import datetime

def packgz(path):
# gzip file from path and return packed file name
    file = path+'.gz'
    with open(path, 'rb') as src, gzip.open(file, 'wb') as dst:
        dst.writelines(src)
    return file

def unpackgz(path):
    # Check if the file is a .gz file
    if not path.endswith('.gz'):
        return path
    
    # Determine the name of the unpacked file
    file = path[:-3]
    
    # Check if the unpacked file already exists
    if os.path.exists(file):
        return file  # Unpacked file already exists, no need to unpack
    
    # Unpack the .gz file
    with gzip.open(path, 'rb') as f_in:
        with open(file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    return file
    
def dir_path(string):
#check if path is valid dir
    if os.path.isdir(string):
        if string[-1] == "/":
            return string[:-1]
        else:
            return string
    else:
        sys.exit(f"\nERROR: {string} is not a valid directory")
        #raise Exception(f"\nERROR: {string} is not a valid directory")

def file_path(string):
#check if path is valid file
    if os.path.isfile(string):
        return string
    else:
        sys.exit(f"\nERROR: {string} is not a valid file")
        #raise Exception(f"\nERROR: {string} is not a valid directory")


def print_header(string,verbose=0):
    if not verbose:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("\n"+string)
        print(current_time)
        print((len(current_time) + 3 + len(string)) * "-")
        
        
        
def clean_empty_files(directory):

    # Loop over all files in the directory
    for filename in os.listdir(directory):
        # Check if the file is not a directory
        if not os.path.isdir(os.path.join(directory, filename)):
            # Check if the file is empty
            if os.path.getsize(os.path.join(directory, filename)) == 0:
                # Remove the file
                os.remove(os.path.join(directory, filename))
                #print(filename, 'was removed')

def generate_color(seed_int):
    random.seed(seed_int)
    color = '#{:06x}'.format(random.randint(0, 0xFFFFFF))
    return color        
        
def getAllFiles(directory, ending = 0):
#get all files of a directory with all subdirectories and return a list
    list = []
    for path, subdirs, files in os.walk(directory):
        for name in files:
            file = os.path.join(path, name)
            if ending == 0:
                list.append(file)
            elif file.endswith(ending):
                list.append(file)
            
    return list        
        
        
        
    
#### Next three routines are meant to work together
def compareFileLists(directory, ext1=0, ext2=0):
    """
    Returns a list of all files with extension 1 which have no equivalent with extension 2,
    considering only file names (ignoring directory paths).
    """
    if ext1 and ext2:
        # Get all files with the specified extensions
        Files1 = getAllFiles(directory, ext1)
        Files2 = getAllFiles(directory, ext2)
        
        # Normalize filenames by stripping directories and extensions
        CompareList1 = removeExtFromList(Files1, ext1)
        CompareList2 = removeExtFromList(Files2, ext2)
        
        # Find files present in ext1 but not in ext2
        Difference = set(CompareList1).difference(set(CompareList2))
        
        # Add back the extension to the differing file names
        listing = addExtToList(list(Difference), ext1)
        
        return listing

    return []

def getAllFiles(directory, ext):
    """
    Recursively retrieves all files with the specified extension in the given directory.
    """
    matched_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(ext):
                matched_files.append(os.path.join(root, file))
    return matched_files

def removeExtFromList(listing, ext):
    """
    Removes extension and directory paths from each element of a list.
    """
    return [os.path.splitext(os.path.basename(element))[0] for element in listing]

def addExtToList(listing, ext):
    """
    Adds an extension back to each element of a list.
    """
    return [element + ext for element in listing]



def getGenomeID(path):
    #return genomeID according to standard genome identifiers from common source. DO NOT USE '.' Where they do not belong in the filename!!!
    basename = os.path.basename(path)
    genomeID = basename.split('.')[0]
    return genomeID

def getReportName(path):
    file_name = os.path.splitext(path)[0]
    hmm_report = file_name+".HmmReport"
    return hmm_report


def taxonomy_lineage(array,trennzeichen):
    #04.11.22
    #concats the taxonomy lineage into one string. 
    #Returned string shall be added at the end of any string
    try:
        div = ''.join(trennzeichen)
        string = div.join(array)
        string = string.replace(" ","-")
        return str(string)
    except:
        return "NoTaxonomy"


def get_executable_dir():
    """
    Get the directory of the current executable or script.
    This works whether the script is compiled or run directly as a Python script.
    """
    if getattr(sys, 'frozen', False):
        # If the program is compiled, sys.frozen is True, and sys.executable gives the path to the executable
        return os.path.dirname(sys.executable)
    else:
        # If running as a script, __file__ gives the path to the script
        return os.path.dirname(os.path.abspath(__file__))

def find_executable(executable):
    """
    Find the any executable.
    First check in the system's PATH, then in the local ./bin directory relative to the executable/script.
    Returns the path to the executable.
    """
    # Check if MAFFT is in the system's PATH
    executable_path = shutil.which(f"{executable}")
    if executable_path:
        return executable_path
    
    # If not found, check in the local ./bin directory relative to the executable or script
    executable_dir = get_executable_dir()
    local_executable_path = os.path.join(executable_dir, "bin", f"{executable}")
    if os.path.isfile(local_executable_path) and os.access(local_executable_path, os.X_OK):
        return local_executable_path
    
    raise FileNotFoundError(f"{executable} executable not found in system PATH or local bin directory.")

def remove_directory(directory_path):
    """
    Removes a directory and all its contents.

    Args:
        directory_path (str): The path of the directory to remove.

    Returns:
        None
    """
    if os.path.exists(directory_path):  # Check if the directory exists
        # Iterate over all files and subdirectories and remove them
        for root, dirs, files in os.walk(directory_path, topdown=False):
            for file in files:
                file_path = os.path.join(root, file)
                os.remove(file_path)  # Remove file
            for dir in dirs:
                dir_path = os.path.join(root, dir)
                os.rmdir(dir_path)  # Remove empty subdirectory
        os.rmdir(directory_path)  # Remove the main directory itself
        print(f"[INFO] Removed directory and all its contents: {directory_path}")
    else:
        print(f"[INFO] Directory does not exist: {directory_path}")
        

def save_cache(options, name, data, redirect=None, overwrite=False):

    cache_dir = os.path.join(options.result_files_directory, "pkl_cache")
    if redirect:
        cache_dir = redirect

    os.makedirs(cache_dir, exist_ok=True)  # sicherer als os.system(mkdir ...)

    file_path = os.path.join(cache_dir, name)
    if os.path.exists(file_path) and overwrite == False:
        return  # Nicht überschreiben, wenn bereits vorhanden

    with open(file_path, "wb") as f:
        pickle.dump(data, f)


def load_cache(options, name, file_path = None):
    cache_dir = os.path.join(options.result_files_directory, "pkl_cache")

    if not file_path:
        file_path = os.path.join(cache_dir, name)
    
    if os.path.exists(file_path):
        print(f"[LOAD] Loading from cache: {name}")
        with open(file_path, "rb") as f:
            return pickle.load(f)
    else:
        print(file_path)
        return None
    
    
def merge_grouped_refseq_dicts_simple(grouped_3_dict, grouped_4_dict):
    """
    Merge grouped_3_dict and grouped_4_dict assuming they have the same domain keys.

    Args:
        grouped_3_dict (dict): {domain: set(proteinIDs)}
        grouped_4_dict (dict): {domain: set(proteinIDs)}

    Returns:
        dict: {domain: set(proteinIDs)} — unified non-redundant sets
    """
    merged = {}

    all_domains = set(grouped_3_dict.keys()) | set(grouped_4_dict.keys())

    for domain in all_domains:
        set_3 = grouped_3_dict.get(domain, set())
        set_4 = grouped_4_dict.get(domain, set())
        merged[domain] = set_3 | set_4  # union of both sets

    return merged
    
    
def merge_score_limits(dict1, dict2):
    """
    Merges two dictionaries of score limits by taking the minimum lower_limit
    and maximum upper_limit for each domain key.

    Args:
        dict1 (dict): {domain: {"lower_limit": int, "upper_limit": int}}
        dict2 (dict): {domain: {"lower_limit": int, "upper_limit": int}}

    Returns:
        dict: Merged score limits with min lower_limit and max upper_limit
    """
    merged = {}
    all_keys = set(dict1.keys()).union(dict2.keys())

    for key in all_keys:
        val1 = dict1.get(key, {})
        val2 = dict2.get(key, {})

        lower1 = val1.get("lower_limit", float('inf'))
        upper1 = val1.get("upper_limit", float('-inf'))
        lower2 = val2.get("lower_limit", float('inf'))
        upper2 = val2.get("upper_limit", float('-inf'))

        merged[key] = {
            "lower_limit": min(lower1, lower2),
            "upper_limit": max(upper1, upper2)
        }

    return merged

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
