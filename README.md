# Generic draft
# HAMSTER: Homolog and Synteny Mining Pipeline

HAMSTER is a modular command-line pipeline for the high-throughput identification and analysis of homologous genes, gene clusters, and collinear syntenic blocks across multiple genome datasets. It provides automated, reproducible workflows for genome mining, synteny block detection, and protein sequence clustering with advanced visualization and reporting.

## Features

- Automated detection of homologs and synteny blocks across large numbers of genomes
- Flexible protein clustering using MCL, PAM, and CSB approaches
- Stepwise, checkpointed workflow—resume or start at any stage
- Supports both simple and highly configurable command-line operation
- Produces rich outputs: gene clusters, alignments, statistical reports, and PDF plots
- Modular structure for extension and integration

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

All outputs are organized into a results directory, including:

### Output Folders

- **results/** 
  Main results directory (contains all project outputs)
    - **Sequences/**
      Protein sequences (FASTA) of assumed functional equivalent sequences from the analysis
    - **Hidden_markov_models/** 
      Generated profile Hidden Markov Models for each protein
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

1. **Clone this repository:**
    ```bash
    git clone https://github.com/TSTanabe/HAMSTER.git
    cd HAMSTER
    ```

2. **Install Python dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3. **(Optional) Install R and required R packages** for PDF plotting (see `plotting_Rscripts/` for details).

## Dependencies

- **Python** 3.8 or higher
- **Python packages**:
  - `numpy`
  - `pandas`
  - `scikit-learn`
  - `scipy`
- **R** (for optional plotting scripts)
- External tools may be required for full functionality 

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
