# Generic draft
# HAMSTER: Homolog and Synteny Mining Pipeline

HAMSTER is a modular command-line pipeline for the high-throughput identification and analysis of homologous genes, gene clusters, and collinear syntenic blocks across multiple genome datasets. It provides automated, reproducible workflows for genome mining, synteny block detection, and protein sequence clustering with advanced visualization and reporting.

## Quick Start

**Typical usage for a new project:**
```bash
python hamster.py -f ./genomes -q queries.faa
```

**Resume from previous results:**
```bash
python hamster.py -r ./results --verbose 2
```

**Show all advanced options:**
```bash
python hamster.py --help-all
```

## Inputs

- **FASTA files**: Directory containing genome assemblies with .fna suffix or .faa with corresponding .gff files (required)
- **Query file**: FASTA file with protein sequences of interest usually encoded in a syntenic gene cluster (required)
- **Results directory**: (optional) Output directory. If an existing results folder is provided resume the analysis

## Outputs
For each new attempt to generate hidden Markov Models a new results folder is created that organizes the output files.
An example output folder structure with the most important output files is given below. The HMMs are located in the Hidden markov models
subdirectory and together with the cutoffs in the _ini_cutoffs.txt file. The training data for each HMM are located in the Sequences
subdirectory and a detailed report on each sequence in each training dataset in the Reports folder. HMMs, training data fasta files
and report files that belong together share the same basename.

### Output Folders

- **results/** 
  Main results directory (contains all project outputs)
    - **Sequences/**
      Protein sequences (FASTA) of assumed functional equivalent sequences from the analysis
    - **Hidden_markov_models/** 
      Generated profile Hidden Markov Models for each protein with each selection rule and corresponding cutoffs
    - **Reports/** 
      Generated detailed reports for each sequence set, including a list of all selected sequences with the genomic vicinity
    - **Collinear_syntenic_blocks/** 
      CSB (synteny block) files, cluster assignments, and instance summaries
    - **Protein_Phylogeny/** 
      Temporary files from the sequence sorting. May accelerate future runs
    - **pkl_cache/** 
      Cached files from the execution. May accelerate future runs
    - **Hit_list/** 
      Filtered BLAST hit tables for each genome/assembly
    - **Initial_validation/** 
      Reports from the HMM search against all hits
    - **Cross_validation/** 
      Cross-validation reports

### Main output files (in results/)

- **database.db** 
  SQLite database storing all annotation, cluster, and synteny data
- **filtered_blast_results_table** 
  Filtered and merged BLAST hits across all genomes
- **div_output_file.faa** 
  Non-redundant representative protein sequences (FASTA)
- **execution_logfile.txt** 
  Detailed logfile for the last run

### Key files in Collinear_syntenic_blocks/

- **Csb_output.txt** 
  Main collinear syntenic block patterns that were found
- **All_gene_clusters.txt** 
  All detected syntenic gene clusters across all genomes

### Key files in Reports/

- **grp[proteins]_enriched.txt** 
  Main report for each sequence in the set, including genomic vicinity and assignment during initial validation
- **all_cutoffs.txt** 
  Optimized, trusted and noise cutoff for all generated HMMs
- **all_performance.txt** 
  Performance of all HMMs during classification of the underlying training data
  
## Installation
HAMSTER is either available via the github directory or as a compiled binary file

1. **Clone this repository:**
    ```bash
    git clone https://github.com/TSTanabe/HAMSTER.git
    cd HAMSTER
    python hamster.py -f ./genomes -q queries.faa
    ```

2. **(Optional) Install R and required R packages** for PDF plotting (see `plotting_Rscripts/` for details).

## Dependencies
### Required libraries
- **Python** 3.8 or higher
- **Python packages**:
  - `numpy`
  - `pandas`
  - `scikit-learn`
  - `scipy`
- **R** (for optional plotting scripts)

### Required external tools for functionality 

  - **DIAMOND** — for fast protein BLAST-like searches (https://github.com/bbuchfink/diamond)
  - **Prodigal** — for gene prediction from prokaryotic genomes (https://github.com/hyattpd/Prodigal)
  - **mmseqs2** — for fast and sensitive protein sequence searching and clustering (https://github.com/soedinglab/MMseqs2)
  - **MAFFT** — for multiple sequence alignment (https://mafft.cbrc.jp/alignment/software/)
  - **trimAl** — for automated alignment trimming (http://trimal.cgenomics.org/downloads)

## Command-line Arguments

HAMSTER supports both essential and advanced configuration options.
- Use `--help` for a concise summary (recommended for most users)
- Use `--help-all` to see all available advanced options

```bash
python hamster.py --help
python hamster.py --help-all
```

## Citation

If you use HAMSTER in your research, please cite:

> [Add your publication/preprint/citation here]

---
