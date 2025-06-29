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

- **FASTA files**: Directory containing genome assemblies (required)
- **Query file**: FASTA file with protein or gene sequences of interest (required)
- **GFF files**: (optional) for annotated features
- **Results directory**: Resume or continue an analysis by providing an existing results folder

## Outputs

All outputs are organized into a results directory, including:
- Gene clusters and syntenic blocks
- Protein sequences and alignments
- Statistical and cross-validation reports
- Plots and figures (PDF, requires R for some scripts)
- Logs and pipeline configuration files

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
- Python packages: `numpy`, `pandas`, `scikit-learn`, `sqlite3`, and others (see `requirements.txt`)
- **R** (for plotting scripts)
- External tools may be required for full functionality (e.g., [DIAMOND](https://github.com/bbuchfink/diamond), [Prodigal](https://github.com/hyattpd/Prodigal), [mmseqs2](https://github.com/soedinglab/MMseqs2)).

## Command-line Arguments

HAMSTER supports both essential and advanced configuration options.
- Use `--help` for a concise summary (recommended for most users)
- Use `--help-all` to see all available advanced options

```bash
python hamster.py --help
python hamster.py --help-all
```

## Documentation

Detailed documentation, including parameter explanations and advanced usage, is available at:
[https://github.com/TSTanabe/HAMSTER](https://github.com/TSTanabe/HAMSTER)

## Citation

If you use HAMSTER in your research, please cite:

> [Add your publication/preprint/citation here]

## License

[Specify your license, e.g., MIT, GPL-3.0, or custom license]

## Contact

For questions, bug reports, or feature requests, please use the [GitHub Issues](https://github.com/TSTanabe/HAMSTER/issues) page  
or contact: [your.email@domain.com] (replace with your actual email)

---

### Acknowledgments

- Developed by [YOUR NAME/LAB/INSTITUTION]
- Powered by open-source tools including Python, R, scikit-learn, pandas, and more.

---
