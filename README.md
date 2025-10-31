# BOLD Library Curation Pipeline

![Perl CI](https://github.com/FabianDeister/Library_curation_BOLD/actions/workflows/ci.yml/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17495973.svg)](https://doi.org/10.5281/zenodo.17495973)
[![DOI](https://zenodo.org/badge/DOI/10.48546/WORKFLOWHUB.WORKFLOW.833.1.svg)](https://doi.org/10.48546/WORKFLOWHUB.WORKFLOW.833.1)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<p align="left">
  <img src="doc/IBOL_LOGO_TRANSPARENT.png" style="width: 200px;" alt="IBOL Logo">
  <img src="doc/boldetective_logo.png" style="width: 200px;" alt="BOLDetective Logo">
  <img src="doc/Logo BGE.png" style="width: 200px;" alt="BGE Logo">
</p>

## Overview

A comprehensive Snakemake pipeline for processing and curating [BOLD Systems](https://boldsystems.org) barcode reference sequence data. This tool implements standardized quality assessment criteria developed by the [Biodiversity Genomics Europe (BGE)](https://biodiversitygenomics.eu) consortium to evaluate and rank DNA barcode sequences for library curation.

### Key Features

- **Comprehensive Quality Assessment**: Evaluates specimens against 16 standardized criteria including metadata completeness, voucher information, sequence quality, and phylogenetic analyses
- **Advanced Phylogenetic Analysis**: Includes haplotype identification and OTU clustering for genetic diversity assessment
- **BAGS Species Assessment**: Automated species-level quality grading system with subspecies inheritance
- **Geographic Representation**: Country representative selection for balanced geographic sampling
- **Scalable Architecture**: Family-level database splitting for efficient analysis of large datasets
- **FAIR Compliance**: Built with reproducibility and provenance tracking using Snakemake workflows

## Background

The classification criteria are actively developed by the BGE consortium and documented in this [living document](https://docs.google.com/document/d/18m-7UnoJTG49TbvTsq_VncKMYZbYVbau98LE_q4rQvA/edit). The pipeline evaluates sequences based on multiple quality dimensions to support evidence-based curation decisions.

## Pipeline Workflow
[DETAILED WORKFLOW](https://github.com/bge-barcoding/bold-library-curation/blob/main/workflow/README.md)

[Example RUN OUTPUT](https://github.com/bge-barcoding/bold-library-curation/blob/main/workflow/example_run_report.md)

The pipeline processes BOLD data through six main phases:

1. **Data Preparation**: Optional pre-filtering by taxa, geography, or genetic markers
2. **Database Setup**: SQLite database creation with taxonomic enrichment
3. **Quality Assessment**: Evaluation against 16 standardized criteria plus phylogenetic analyses
4. **BAGS Assessment**: Species-level quality grading with database optimization
5. **Data Integration**: Ranking system combining all assessments with country representative selection
6. **Family Splitting**: Creation of family-level databases for scalable downstream analysis

### Assessment Criteria

**Specimen Metadata**: Collection date, collectors, identifier, identification method  
**Geographic Data**: Country, region, site, sector, coordinates  
**Repository Info**: Institution, museum ID, public voucher  
**Sequence Quality**: DNA quality metrics, species ID, type specimen, images  
**Phylogenetic**: Haplotype identification, OTU clustering

## Requirements

- **Operating System**: Linux, macOS, or Windows with WSL
- **Package Manager**: Mamba (recommended) or Conda
- **Memory**: Minimum 8GB RAM (16GB+ recommended for large datasets)
- **Storage**: Sufficient space for BOLD data and family databases
- **Dependencies**: VSEARCH (for OTU clustering), SQLite, Perl, Python
 Data**: Country, region, site, sector, coordinates  
**Repository Info**: Institution, museum ID, public voucher  
**Sequence Quality**: DNA quality metrics, species ID, type specimen, images  
**Phylogenetic**: Haplotype identification, OTU clustering

## Requirements

- **Operating System**: Linux, macOS, or Windows with WSL
- **Package Manager**: Mamba (recommended) or Conda
- **Memory**: Minimum 8GB RAM (16GB+ recommended for large datasets)
- **Storage**: Sufficient space for BOLD data and family databases
- **Dependencies**: VSEARCH (for OTU clustering), SQLite, Perl, Python

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

## Quick Start

1. **Configure your analysis** in `config/config.yml`:
   ```yaml
   # Input data
   BOLD_TSV: "resources/your_bold_data.tsv"
   
   # Optional filtering
   ENABLE_PRESCORING_FILTER: false
   USE_TARGET_LIST: false
   
   # Analysis parameters
   OTU_CLUSTERING_THRESHOLD: 0.99
   FAMILY_SIZE_THRESHOLD: 10000
   ```

2. **Run the pipeline:**
   ```bash
   snakemake --cores 4 --use-conda
   ```

3. **Check results** in `results/result_output.tsv` and `results/family_databases/`

## Configuration Options

### Core Settings
- **BOLD_TSV**: Path to BOLD data dump (BCDM TSV format)
- **RESULTS_DIR/LOG_DIR**: Customizable output directories
- **TAXONOMY_CHUNK_SIZE**: Memory optimization (default: 10,000)

### Optional Filtering
- **Pre-scoring Filter**: Early dataset reduction by taxa, countries, markers, or BIN sharing
- **Target Lists**: Focus on specific species of interest
- **Geographic Filtering**: Country-based specimen filtering

### Analysis Parameters
- **OTU_CLUSTERING_THRESHOLD**: Genetic similarity for clustering (default: 0.99)
- **OTU_CLUSTERING_THREADS**: Parallel processing threads (default: 8)
- **FAMILY_SIZE_THRESHOLD**: Minimum records for individual family databases (default: 10,000)

## Usage Examples

### Basic Analysis
```bash
# Standard pipeline execution
snakemake --cores 8 --use-conda
```

### Large Dataset with Pre-filtering
```bash
# Configure pre-filtering in config.yml first
# ENABLE_PRESCORING_FILTER: true
# FILTER_TAXA: true
# FILTER_TAXA_LIST: "resources/target_taxa.txt"
snakemake --cores 16 --use-conda --resources mem_mb=32000
```

### HPC Cluster Execution (SLURM)
```bash
#!/bin/bash
#SBATCH --partition=day
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

source activate bold-curation
snakemake --cores 16 --use-conda
```

### Target Species Focus
```bash
# Enable target list in config.yml
# USE_TARGET_LIST: true
# TARGET_LIST: "resources/target_species.csv"
snakemake --cores 8 --use-conda
```

## Output Files

### Primary Outputs
- **`results/result_output.tsv`**: Final scored and ranked specimens with all assessments
- **`results/family_databases/`**: Family-level SQLite databases organized by phylum
- **`results/pipeline_summary.txt`**: Comprehensive execution summary

### Assessment Results
- **Individual criteria files**: `assessed_*.tsv` for each quality criterion
- **BAGS assessment**: Species-level quality grades with BIN sharing analysis
- **Phylogenetic analyses**: Haplotype IDs and OTU clustering results
- **Database**: Complete SQLite database with specialized tables for complex queries

### Performance Monitoring
- **Comprehensive logging**: Step-by-step execution logs in `logs/` directory
- **Progress tracking**: Real-time monitoring for long-running operations
- **Error handling**: Detailed debugging information for troubleshooting

## Performance Guidelines

### Dataset Size Recommendations
- **Small** (< 10K records): 30-60 minutes, 8GB RAM
- **Medium** (10K-100K records): 2-6 hours, 16GB RAM  
- **Large** (100K+ records): 6-24 hours, 32GB+ RAM
- **Very large** (1M+ records): 1-3 days, consider pre-filtering

### Optimization Tips
1. **Use pre-scoring filter** for very large datasets
2. **Adjust memory settings** based on available resources
3. **Configure threading** for OTU clustering based on CPU cores
4. **Use SSD storage** for database operations when possible

## Project Structure

```
bold-library-curation/
├── config/                 # Configuration files
├── workflow/              # Snakemake workflow
│   ├── bold-ranker.smk   # Main workflow definition
│   ├── scripts/          # Analysis scripts
│   ├── envs/             # Conda environments
│   └── README.md         # Detailed workflow documentation
├── resources/            # Input data and references
├── results/              # Output files and databases
└── logs/                 # Execution logs
```

## Advanced Features

### Phylogenetic Integration
- **Haplotype Analysis**: Identifies genetic variants within species/BIN groups
- **OTU Clustering**: VSEARCH-based clustering with configurable similarity thresholds
- **Multi-threaded Processing**: Parallel execution for computationally intensive steps

### Geographic Representation
- **Country Representatives**: Systematic selection of optimal specimens per region
- **Balanced Sampling**: Maintains geographic diversity while optimizing quality

### Database Architecture
- **Hierarchical Organization**: Family-level databases organized by phylum
- **Specialized Tables**: Optimized schema for different data types
- **Efficient Querying**: Comprehensive indexing for complex analyses

## Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) and:

1. Fork the repository
2. Create a feature branch
3. Make changes with appropriate tests
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
- **Documentation**: Detailed workflow documentation in `workflow/README.md`
- **BGE Consortium**: Visit [biodiversitygenomics.eu](https://biodiversitygenomics.eu) for project updates

## Acknowledgments

This work is supported by the Biodiversity Genomics Europe (BGE) consortium and contributes to the International Barcode of Life (IBOL) initiative. Biodiversity Genomics Europe (Grant no.101059492) is funded by Horizon Europe under the Biodiversity, Circular Economy and Environment call (REA.B.3); co-funded by the Swiss State Secretariat for Education, Research and Innovation (SERI) under contract numbers 22.00173 and 24.00054; and by the UK Research and Innovation (UKRI) under the Department for Business, Energy and Industrial Strategy's Horizon Europe Guarantee Scheme.
