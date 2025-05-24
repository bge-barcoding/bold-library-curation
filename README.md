# BOLD Library Curation Pipeline

![Perl CI](https://github.com/FabianDeister/Library_curation_BOLD/actions/workflows/ci.yml/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11267905.svg)](https://doi.org/10.5281/zenodo.11267905)
[![DOI](https://zenodo.org/badge/DOI/10.48546/WORKFLOWHUB.WORKFLOW.833.1.svg)](https://doi.org/10.48546/WORKFLOWHUB.WORKFLOW.833.1)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<img src="doc/IBOL_LOGO_TRANSPARENT.png" style="width: 200px;" alt="IBOL Logo">

## Overview

An automated pipeline for curating [BOLD Systems](https://boldsystems.org) barcode reference sequence data, implementing classification criteria developed by the [Biodiversity Genomics Europe (BGE)](https://biodiversitygenomics.eu) consortium. This tool processes BOLD data dumps in BCDM TSV format to classify barcode sequences based on standardized quality criteria.

### Key Features

- **Automated Quality Assessment**: Evaluates barcode sequences against multiple criteria including collection metadata, voucher information, and sequence quality
- **FAIR Compliance**: Built with reproducibility and provenance tracking using Snakemake workflows
- **Standardized Classification**: Implements BGE consortium criteria for barcode reference sequence classification
- **Scalable Processing**: Supports both local execution and HPC cluster deployment

## Background

The classification criteria being implemented are actively developed by the BGE consortium and documented in this [living document](https://docs.google.com/document/d/18m-7UnoJTG49TbvTsq_VncKMYZbYVbau98LE_q4rQvA/edit). The pipeline evaluates sequences based on:

- Collection date and collector information
- Geographic coordinates and location data
- Specimen voucher and museum information
- Sequence quality metrics
- Taxonomic identification details
- Image availability

## Requirements

- **Operating System**: Linux, macOS, or Windows with WSL
- **Package Manager**: Mamba (recommended) or Conda
- **Memory**: Minimum 8GB RAM (24GB recommended for large datasets)
- **Storage**: Sufficient space for BOLD data dumps (typically several GB)

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/bge-barcoding/bold-library-curation.git
   cd bold-library-curation
   ```

2. **Set up the environment:**
   ```bash
   mamba env create -f environment.yml
   mamba activate bold-curation
   ```

## Configuration

Before running the pipeline, configure your analysis in `config/config.yml`:

```yaml
# Input data
BOLD_TSV: resources/BOLD_Public.05-Apr-2024.tsv

# Analysis parameters
TAXON_LEVEL: "species"
KINGDOM: "Animalia"
CRITERIA: "COLLECTION_DATE COLLECTORS COORD COUNTRY..."

# Project settings
PROJECT_NAME: "bold-curation"
LOG_LEVEL: "INFO"
```

### Key Configuration Options

- **BOLD_TSV**: Path to your BOLD data dump file
- **TARGET_LIST**: Path to species and synonymy data (CSV format)
- **CRITERIA**: Space-separated list of quality criteria to evaluate
- **TAXON_LEVEL**: Taxonomic level for analysis (species, genus, etc.)
- **KINGDOM**: Target taxonomic kingdom

## Usage

### Basic Execution

Run the complete pipeline with default settings:

```bash
snakemake -p -c 4
```

### HPC Cluster Execution (SLURM)

For high-performance computing environments:

```bash
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_curate_bold_%j.out
#SBATCH --error=job_curate_bold_%j.err
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4

source activate bold-curation

# Clean previous runs (backup results first if needed)
snakemake --unlock -p -c 4 clean

# Execute pipeline
snakemake -p -c 4

echo "Pipeline complete!"
```

### Incremental Execution

The pipeline supports incremental execution. To run specific steps or restart from a checkpoint:

```bash
# Run only database creation
snakemake -p -c 4 results/bold.db

# Clean previous results (with unlock for interrupted runs)
snakemake --unlock -p -c 4 clean
```

## Output

The pipeline generates several output files in the `results/` directory:

- **bold.db**: SQLite database containing processed BOLD data
- **taxonomy_check.tsv**: Taxonomic validation results  
- **criteria_indexed**: Quality criteria evaluation results
- **Log files**: Detailed execution logs for troubleshooting

## Project Structure

```
bold-library-curation/
├── config/                 # Configuration files
│   └── config.yml         # Main configuration
├── workflow/              # Snakemake workflow
│   ├── Snakefile         # Main workflow definition
│   ├── scripts/          # Analysis scripts
│   └── envs/             # Environment specifications
├── resources/            # Input data and references
├── results/              # Output files
├── lib/                  # Library files
└── doc/                  # Documentation
```

## Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes with appropriate tests
4. Submit a pull request

## License

This project is licensed under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.

## Citation

If you use this tool in your research, please cite:

```
[Citation information - to be added when published]
```

**DOI**: [10.5281/zenodo.11267905](https://doi.org/10.5281/zenodo.11267905)

## Support

- **Issues**: Report bugs and request features via [GitHub Issues](https://github.com/bge-barcoding/bold-library-curation/issues)
- **Documentation**: Additional documentation available in the `doc/` directory
- **BGE Consortium**: Visit [biodiversitygenomics.eu](https://biodiversitygenomics.eu) for project updates

## Acknowledgments

This work is supported by the Biodiversity Genomics Europe (BGE) consortium and contributes to the International Barcode of Life (IBOL) initiative. Biodiversity Genomics Europe (Grant no.101059492) is funded by Horizon Europe under the Biodiversity, Circular Economy and Environment call (REA.B.3); co-funded by the Swiss State Secretariat for Education, Research and Innovation (SERI) under contract numbers 22.00173 and 24.00054; and by the UK Research and Innovation (UKRI) under the Department for Business, Energy and Industrial Strategy’s Horizon Europe Guarantee Scheme
