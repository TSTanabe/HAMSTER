# HAMSTER: HMM generation with synteny support

HAMSTER is a modular command-line pipeline for the high-throughput identification of homologous genes with collinear syntenic blocks across multiple genome datasets. It provides an automated, reproducible workflow to sort sequences into functional equivalent groups based on genomic synteny block detection, presence propability and protein sequence clustering.

## Quick Start
### Usage of HAMSTER
**Typical usage for a new project:**
```bash
python HAMSTER.py -f ./genomes -q queries.faa
```

**Resume from previous results:**
```bash
python HAMSTER.py -r ./results --verbose 2
```

**Show all advanced options:**
```bash
python HAMSTER.py --help-all
```

### Installation with Conda
HAMSTER can be installed from the Git repository using the provided Conda environment file.
```bash
git clone <repository-url>
cd HAMSTER
conda env create -f hamster_env.yml
conda activate hamster_env
```

### Manual installation

Using the provided Conda environment is recommended, because HAMSTER depends on both Python packages and external bioinformatics command-line tools. If HAMSTER is installed without Conda, the following software must be available in the active environment and accessible through `$PATH`.

**Required Python packages**
- Python >= 3.8
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- pyrodigal
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- pyrodigal

**Required external tools for functionality** 
  - **DIAMOND** — for fast protein BLAST-like searches (https://github.com/bbuchfink/diamond)
  - **mmseqs2** — for fast and sensitive protein sequence searching and clustering (https://github.com/soedinglab/MMseqs2)
  - **MAFFT** — for multiple sequence alignment (https://mafft.cbrc.jp/alignment/software/)
  - **trimAl** — for automated alignment trimming (http://trimal.cgenomics.org/downloads)
  - **HMMER** — for HMM generation and testing (http://hmmer.org/download.html)



## Inputs

- **FASTA files**:  (required) Directory containing genome assemblies with .fna suffix or .faa with corresponding .gff files.
- **Query file**: (required) FASTA file with protein sequences. These are assumed to be commonly encoded in a syntenic gene cluster or occur in the same genome. Header should not include whitespaces.  
- **Results directory**: (optional) Output directory. If an existing results folder is provided resume the analysis.

HAMSTER can scale to > 300 000 input genomes. It is recommended to include genomes that do encode homologs of the query
proteins. If too few sequences are found that assumed, it might be necessary to adjust the minimal collinear syntenic block size (--min-csb-size)
or minimal sequence identity for singletons (--singleton-identity-cutoff) to select hits with small or without a conserved synteny.
Adjustments of the DIAMOND blastp search parameters can also increase/decrease sensitivity and number of seed sequences. 

## Outputs
For each new attempt to generate hidden Markov Models a new results folder is created that organizes the output files.
An example output folder structure with the most important output files is given below. The HMMs are located in the Hidden markov models
subdirectory and together with the cutoffs in the _ini_cutoffs.txt file. The training data for each HMM are located in the Sequences
subdirectory and a detailed report on each sequence in each training dataset in the Reports folder. HMMs, training data fasta files
and report files that belong together share the same basename.

### Main results

The main output directory contains the final HMMs, selected training sequences, and detailed validation reports.


- **Hidden_markov_models/**  
  Contains the generated profile Hidden Markov Models (HMMs) for each protein family and selection rule, together with the corresponding cutoff information.

  - **all_cutoffs.txt**  
    Summary table of the optimized, trusted, and noise cutoffs for all generated HMMs.

  - **all_performance.txt**  
    Summary table of HMM performance during classification of the underlying training data.

  - **cv_cutoff_performance.txt**  
    HMM performance during cross-validation

- **Sequences/**  
  Seed alignments and unaligned sequences for the HMM generation

- **Reports/**  
  Contains detailed reports for each sequence set, including selected sequences, classification results, and information on their genomic neighborhood.
- **Collinear_syntenic_blocks/**  
  Gene cluster pattern that were detected 

Other subdirectories include temporary files from the HAMSTER run.

## Citation

If you use HAMSTER in your research, please cite:

> [Placeholder citation]

---
