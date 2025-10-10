# BOLD Library Curation Pipeline

This Snakemake workflow processes BOLD sequence data through comprehensive quality assessment criteria, performs phylogenetic analyses, and splits results into family-level databases for efficient analysis and curation.

## Pipeline Overview

The pipeline follows a systematic approach to evaluate BOLD specimens across multiple quality criteria, perform advanced phylogenetic analyses using parallel SLURM job arrays, and produce ranked outputs optimized for library curation. The workflow includes optional pre-filtering capabilities, OTU clustering, species-level BAGS assessment, parallel phylogenetic tree generation with curation checklists, and family-level database splitting for scalable downstream analysis.

![BOLD Pipeline workflow](/doc/pipeline_summary.png)

## Workflow Phases

### Phase 1: Data Preparation and Filtering

#### Pre-scoring Filter (Optional)
**Rules:** `prescoring_filter` or `skip_prescoring_filter`
- **Optional step**: Filters the input BOLD dataset before detailed processing
- Can filter by taxa, countries, genetic markers, and/or BIN sharing criteria
- **When enabled** (`ENABLE_PRESCORING_FILTER: true`): Runs filtering script to create filtered dataset
- **When disabled** (`ENABLE_PRESCORING_FILTER: false`): Skips filtering entirely - downstream steps use original file directly
- **Optimization**: No file copying when filtering is disabled, improving efficiency
- **Configuration**: Controlled by `ENABLE_PRESCORING_FILTER` in config
- **Input:** Original BOLD TSV file (when filtering enabled)
- **Output:** Filtered BOLD TSV file (when filtering enabled) or marker file only (when disabled)
- **Script:** `workflow/scripts/prescoring_filter.py` (when enabled)

**Available Filtering Options:**
- **Taxa filtering**: Filter by specific taxonomic groups (`FILTER_TAXA: true`, `FILTER_TAXA_LIST`)
- **Country filtering**: Filter by countries of collection (`FILTER_COUNTRIES: true`, `FILTER_COUNTRY_LIST`)
- **Marker filtering**: Filter by specific genetic markers (`MARKER`)
- **BIN sharing**: Enable BIN sharing analysis (`FILTER_BINS: true`)
- **Species filtering**: Filter for species-level identifications only (`FILTER_SPECIES: true`)

### Phase 2: Database Creation and Setup

#### Database Creation and Data Loading
**Rule:** `create_load_db`
- Creates SQLite database from BOLD TSV data using optimized fast_simple loader
- Uses dynamic input selection: filtered file (when filtering enabled) or original file (when filtering disabled)
- Loads the Barcode Core Data Model (BCDM) schema
- **Input:** BOLD TSV file (selected automatically based on filtering configuration), database schema
- **Output:** SQLite database file
- **Script:** `workflow/scripts/load_bcdm_fast_simple.pl`

#### Criteria Configuration Loading
**Rule:** `load_criteria`
- Imports assessment criteria definitions into the database
- **Input:** `resources/criteria.tsv`
- **Output:** Database with criteria table loaded

#### Database Indexing
**Rule:** `apply_indexes`
- Applies database indexes for improved query performance
- **Input:** Index definitions SQL file
- **Output:** Indexed database file

#### Taxonomy Loading
**Rule:** `load_taxonomy`
- Loads NCBI taxonomic information into database using optimized chunked loading
- Enriches specimen records with taxonomic hierarchy
- **Script:** `workflow/scripts/load_taxonomy_faster.pl`
- **Configuration:** `TAXONOMY_CHUNK_SIZE` for memory optimization

#### Target List Import (Optional)
**Rules:** `import_target_list` or `skip_target_list`
- **Conditional step**: Only runs if `USE_TARGET_LIST: true` in config
- When enabled: Filters specimens to only target species from provided list
- When disabled: Creates pass-through dependency for all specimens
- **Input:** Target species CSV file (when enabled)
- **Script:** `workflow/scripts/load_targetlist.pl`
- **Purpose:** Allows focused curation on specific species of interest

### Phase 3: Quality Criteria Assessment

The pipeline evaluates specimens against multiple quality criteria. Each criterion is assessed independently:

#### Specimen Metadata Criteria
- **COLLECTION_DATE**: Validates collection date completeness and validity
- **COLLECTORS**: Assesses collector information completeness
- **IDENTIFIER**: Evaluates taxonomic identifier information
- **ID_METHOD**: Checks identification method documentation

#### Geographic Information Criteria
- **COUNTRY**: Validates country information completeness and standardization
- **REGION**: Validates geographic region information
- **SITE**: Validates collection site information
- **SECTOR**: Validates geographic sector information
- **COORD**: Validates geographic coordinate completeness and precision

#### Institutional and Repository Criteria
- **INSTITUTION**: Assesses institutional affiliation completeness
- **MUSEUM_ID**: Validates museum/institution specimen ID
- **PUBLIC_VOUCHER**: Checks public voucher specimen availability

#### Specimen Quality Criteria
- **SEQ_QUALITY**: Evaluates DNA sequence quality metrics
- **SPECIES_ID**: Validates species identification completeness and accuracy
- **TYPE_SPECIMEN**: Identifies type specimen designation and documentation
- **HAS_IMAGE**: Checks specimen image availability using CAOS API

#### Advanced Molecular Analyses
- **OTU_CLUSTERING**: Performs VSEARCH-based clustering to identify Operational Taxonomic Units
  - Configurable similarity threshold (`OTU_CLUSTERING_THRESHOLD`, default: 0.99)
  - Multi-threaded processing (`OTU_CLUSTERING_THREADS`)
  - Temporary file management for large datasets
  - Exports unassigned sequences for quality review

Each assessment rule:
- **Input:** Database and taxonomy/target list dependency (conditional)
- **Output:** TSV file with assessment results
- **Script:** `workflow/scripts/assess_criteria.pl` (standard criteria), specialized scripts for images and OTUs
- **Dependency:** Uses `get_taxonomy_dependency()` function to determine assessment dependencies

### Phase 4: BAGS Assessment and Optimization

#### BAGS Database Optimization
**Rule:** `optimize_bags_database`
- Applies BAGS-specific database optimizations for improved performance
- **Optimizations:** WAL mode, memory settings, BAGS-specific indexes
- **Script:** `workflow/scripts/bags_indexes.sql`

#### BAGS Assessment (Species-level Quality Grading)
**Rule:** `BAGS`
- Performs comprehensive species-level assessment using simplified output format
- **BAGS Grading System**: Evaluates species based on taxonomic agreement within BINs
  - **Grade A**: High confidence - strong taxonomic-genetic concordance
  - **Grade B**: Moderate confidence - some discordance but resolvable
  - **Grade C**: Low confidence - significant taxonomic-genetic conflict
  - **Grade D**: Very low confidence - major discordance or insufficient data
- **Output format:** 4 essential columns (taxonid, BAGS_grade, BIN_URL, sharers)
- **Performance:** Optimized with progress tracking and simplified column structure
- **Script:** `workflow/scripts/assess_taxa_simplified.pl`

#### BAGS Data Integration
**Rules:** `import_bags`, `inherit_subspecies_bags`
- Imports BAGS results into database
- Inherits BAGS grades for subspecies from parent species
- Enables complex queries combining BAGS with other criteria

### Phase 5: Parallel Phylogenetic Analysis (Optional)

#### Phylogenetic Batch Preparation
**Rule:** `prepare_phylo_batches`
- Analyzes database to identify families suitable for phylogenetic analysis
- Creates job batches for parallel processing via SLURM job arrays
- **Configuration:**
  - `PHYLO_ENABLED`: Enable/disable phylogenetic analysis (default: false)
  - `PHYLO_NUM_JOBS`: Number of parallel jobs (default: 50)
  - `PHYLO_FAMILIES_PER_JOB`: Families per job batch (optional)
  - `PHYLO_MIN_OTUS`: Minimum OTUs required for analysis (default: 3)
  - `PHYLO_KINGDOMS`: Target kingdoms for analysis (default: ["all"])
- **Script:** `workflow/scripts/prepare_phylo_batches.py`

#### Parallel Phylogenetic Analysis
**Rule:** `run_phylogenetic_analysis_parallel`
- Executes phylogenetic analysis across families using SLURM job arrays
- **Features:**
  - Species-BIN representative selection with enhanced outgroup strategy
  - Hierarchical outgroup selection (order → class fallback)
  - Multiple alignment methods: MAFFT (default), MUSCLE
  - Multiple tree methods: FastTree (default), IQ-TREE
  - Bootstrap support analysis (configurable iterations)
  - BAGS grade integration and Grade C monophyly checking
  - Automated PDF tree visualization with ETE3
  - Curation checklist PDF generation for Grade C species

- **Configuration Options:**
  - `PHYLO_ALIGNMENT_METHOD`: "mafft" or "muscle" (default: "mafft")
  - `PHYLO_TREE_METHOD`: "fasttree" or "iqtree" (default: "fasttree")
  - `PHYLO_BOOTSTRAP`: Bootstrap iterations (default: 1000)
  - `PHYLO_NUM_OUTGROUPS`: Target outgroup count (default: 3)
  - `PHYLO_GENERATE_PDFS`: Enable tree PDF generation (default: true)
  - `PHYLO_JOB_MEMORY`: Memory per job (default: "32G")
  - `PHYLO_JOB_TIME`: Job time limit (default: "24:00:00")
  - `PHYLO_MAX_CONCURRENT`: Maximum concurrent jobs (default: 10)
  - `PHYLO_PARTITION`: SLURM partition (default: "week")
  - `PHYLO_CUSTOM_PARAMETERS`: Tool-specific parameters (dict, e.g., {"mafft": ["--maxiterate", "1000"], "fasttree": ["-fastest"]})

- **Output Products:**
  - Phylogenetic trees in Newick format
  - High-quality PDF visualizations with bootstrap support
  - Curation checklist PDFs for problematic species (Grade C)
  - Family-level analysis summaries with statistics
  - Failed family reports for troubleshooting

- **Scripts:** 
  - `workflow/scripts/phylo_array_job.sh` (SLURM array job wrapper)
  - `workflow/scripts/phylo_pipeline.py` (main phylogenetic analysis)
  - `workflow/scripts/consolidate_phylo_results.py` (result consolidation)

#### Enhanced Outgroup Selection Strategy
The phylogenetic pipeline uses a sophisticated outgroup selection approach:
1. **Primary Strategy**: Select outgroups from other families within the same taxonomic order
2. **Fallback Strategy**: If insufficient outgroups found, expand to other orders within the same class
3. **Flexible Requirements**: Accept 1-2 outgroups instead of requiring 3 when necessary
4. **Quality Filtering**: Prioritize high-quality specimens based on BOLD assessment criteria

This enhancement significantly improves success rates for families that previously failed due to outgroup limitations.

### Phase 6: Data Integration and Output

#### Results Consolidation
**Rule:** `concatenate`
- Combines individual criteria assessment results (excluding OTUs)
- **Input:** All standard criterion TSV files
- **Output:** `results/CONCATENATED.tsv`
- **Script:** `workflow/scripts/concat_tsvs.pl`

#### Specialized Data Import
**Rules:** `import_concatenated`, `import_otus`
- Imports consolidated results and OTU assignments into database
- Creates specialized tables: `bold_criteria`, `bold_otus`
- **Scripts:** Specialized loaders for each data type

#### Ranking and Selection
**Rules:** `create_ranks_schema`, `apply_ranking_indexes`, `calculate_store_ranks`
- Creates comprehensive ranking system combining all assessments
- Applies optimized indexes for ranking queries
- **Script:** `workflow/scripts/calculate_store_ranks.sql`

#### Country Representative Selection
**Rule:** `select_country_representatives`
- Selects best representative record per species per OTU per country
- **Selection criteria:** ranking ASC, sumscore DESC, recordid ASC
- **Filter:** Species-level identification only
- **Script:** `workflow/scripts/select_country_representatives.sql`

#### Manual Curation Support
**Rule:** `populate_manual_curation`
- Generates BOLD portal URLs for all specimens to facilitate manual curation
- **Output:** Database table with direct links to BOLD specimen pages
- **Script:** `workflow/scripts/populate_manual_curation.sql`

#### Final Output Generation
**Rule:** `output_filtered_data`
- Generates final scored and ranked output with all assessments including OTUs
- **Output:** `results/result_output.tsv`
- **Script:** `workflow/scripts/ranking_with_stored_ranks_otu.sql`

### Phase 7: Family-Level Database Creation (Parallel)

#### Family Database Splitting
**Rule:** `split_families`
- Splits main database into family-level databases using parallel SLURM job arrays
- **Performance Features:**
  - Parallel processing via SLURM job arrays
  - Configurable batch sizes and worker threads
  - Kingdom-specific export capabilities
  - Intelligent family grouping for small families
- **Configuration:**
  - `FAMILY_SIZE_THRESHOLD`: Minimum records for individual family database (default: 10,000)
  - `FAMILY_ARRAY_SIZE`: Number of parallel jobs (default: 64)
  - `WORKERS_PER_JOB`: Worker threads per job (default: 4)
  - `EXPORT_KINGDOMS`: Target kingdoms (default: ["all"])
  - `JOB_MEMORY`: Memory per job (default: "8G")
  - `JOB_TIME`: Job time limit (default: "04:00:00")
- **Output:** Individual SQLite databases per family with hierarchical organization
- **Scripts:** 
  - `workflow/scripts/prepare_family_batches.py` (batch preparation)
  - `workflow/scripts/family_split_array.sh` (SLURM array job)
  - `workflow/scripts/consolidate_results.py` (result consolidation)

#### Phylogenetic Results Integration
**Rule:** `integrate_phylogenetic_results`
- Integrates phylogenetic tree PDFs and curation checklist PDFs into family directories
- **Integration Strategy:**
  - Recursively scans family database directories
  - Copies tree PDF files to family directory root
  - Copies curation checklist PDFs for easy access
  - Maintains directory structure for compressed archives
- **Output:** Family directories containing databases, tree PDFs, and curation checklists

#### Database Compression
**Rule:** `compress_family_databases`
- Compresses family directories (including databases, PDFs, and other files) for storage efficiency
- **Compression Modes:**
  - `family_directories`: Compress entire family directories (default)
  - `individual_files`: Compress only database files (legacy mode)
- **Configuration:**
  - `COMPRESSION_WORKERS`: Parallel compression threads (default: 16)
  - `COMPRESSION_MODE`: Compression strategy (default: "family_directories")
- **Performance:** Multi-threaded compression with progress tracking
- **Script:** `workflow/scripts/zip_databases.py`

### Phase 8: Statistical Reporting and Archiving

#### Comprehensive Statistics Report
**Rule:** `generate_stats_report`
- Generates detailed PDF report with database statistics and quality metrics
- **Report Features:**
  - Dataset overview and processing statistics
  - Quality criteria distribution analysis
  - Geographic and taxonomic coverage maps
  - BAGS grade distribution and analysis
  - OTU clustering statistics
  - Phylogenetic analysis summary (when enabled)
  - Performance benchmarks and recommendations
- **Configuration:**
  - `SKIP_DETAILED_TAXONOMY`: Skip detailed taxonomic analysis for speed (default: true)
  - `STATS_REPORT_FAST_MODE`: Enable fast mode for large datasets (default: false)
  - `STATS_REPORT_SAMPLE_SIZE`: Sample size for analysis (default: 100,000)
  - `STATS_REPORT_CPU`: CPU cores for analysis (default: 8)
  - `STATS_REPORT_MEMORY_MB`: Memory allocation (default: 32,000 MB)
  - `STATS_REPORT_RUNTIME`: Time limit in minutes (default: 120)
- **Output:** `results/bold_database_statistics.pdf`
- **Script:** `workflow/scripts/bold_stats_report_hpc.py`

#### Final Results Archiving
**Rule:** `archive_final_results`
- Creates timestamped archive of key pipeline outputs for final delivery
- **Archived Components:**
  - Compressed main database (`bold.db.gz`)
  - Compressed final results (`result_output.tsv.gz`)
  - Compressed family databases archive (`family_databases_compressed.tar.gz`)
  - Complete log archive (`logs.tar.gz`)
  - Statistics report PDF (`bold_database_statistics.pdf`)
- **Archive Structure:** `results/final_results_YYYYMMDD_HHMMSS/`
- **Benefits:** Streamlined delivery package with all essential outputs

## Configuration

The pipeline is configured through `config/config.yml` which defines:

### Core Configuration
- Input file paths (BOLD TSV, schema, indexes)
- Database file locations and directory paths
- Logging levels (`LOG_LEVEL`)
- Library paths (`LIBS`)
- Results and log directories (`RESULTS_DIR`, `LOG_DIR`)

### Pre-scoring Filter Configuration
- `ENABLE_PRESCORING_FILTER`: Enable/disable pre-filtering (default: false)
- `PRESCORING_FILTERED_OUTPUT`: Output filename for filtered data
- `FILTER_TAXA`: Enable taxonomic filtering
- `FILTER_TAXA_LIST`: Path to taxa list file
- `FILTER_COUNTRIES`: Enable country filtering
- `FILTER_COUNTRY_LIST`: Path to countries list file
- `MARKER`: Specific genetic marker filtering
- `FILTER_BINS`: Enable BIN sharing analysis
- `FILTER_SPECIES`: Filter for species-level identifications only

### Target List Configuration
- `USE_TARGET_LIST`: Enable/disable target species filtering
- `TARGET_LIST`: Path to target species CSV file
- `PROJECT_NAME`: Project identifier for target list
- `TAXON_LEVEL`: Taxonomic level for target matching
- `KINGDOM`: Kingdom scope for target filtering

### Advanced Analysis Configuration
- `TAXONOMY_CHUNK_SIZE`: Memory optimization for taxonomy loading (default: 10,000)
- `OTU_CLUSTERING_THRESHOLD`: Similarity threshold for OTU clustering (default: 0.99)
- `OTU_CLUSTERING_THREADS`: Thread count for OTU clustering (default: 8)

### Phylogenetic Analysis Configuration
- `PHYLO_ENABLED`: Enable phylogenetic analysis (default: false)
- `PHYLO_NUM_JOBS`: Number of parallel jobs (default: 50)
- `PHYLO_FAMILIES_PER_JOB`: Families per job batch (optional)
- `PHYLO_MIN_OTUS`: Minimum OTUs for analysis (default: 3)
- `PHYLO_KINGDOMS`: Target kingdoms (default: ["all"])
- `PHYLO_ALIGNMENT_METHOD`: Alignment tool ("mafft" or "muscle", default: "mafft")
- `PHYLO_TREE_METHOD`: Tree method ("fasttree" or "iqtree", default: "fasttree")
- `PHYLO_BOOTSTRAP`: Bootstrap iterations (default: 1000)
- `PHYLO_NUM_OUTGROUPS`: Target outgroup count (default: 3)
- `PHYLO_GENERATE_PDFS`: Enable tree PDF generation (default: true)
- `PHYLO_JOB_MEMORY`: Memory per job (default: "32G")
- `PHYLO_JOB_TIME`: Job time limit (default: "24:00:00")
- `PHYLO_MAX_CONCURRENT`: Maximum concurrent jobs (default: 10)
- `PHYLO_PARTITION`: SLURM partition (default: "week")
- `PHYLO_OUTPUT_DIR`: Output directory name (default: "phylogenies")
- `PHYLO_CUSTOM_PARAMETERS`: Tool-specific parameters (dict, e.g., {"mafft": ["--maxiterate", "1000"]})

### Family Database Configuration
- `FAMILY_SIZE_THRESHOLD`: Minimum records for individual family database (default: 10,000)
- `FAMILY_ARRAY_SIZE`: Number of parallel jobs (default: 64)
- `WORKERS_PER_JOB`: Worker threads per job (default: 4)
- `EXPORT_KINGDOMS`: Target kingdoms (default: ["all"])
- `JOB_MEMORY`: Memory per job (default: "8G")
- `JOB_TIME`: Job time limit (default: "04:00:00")

### Compression and Output Configuration
- `COMPRESSION_WORKERS`: Parallel compression threads (default: 16)
- `COMPRESSION_MODE`: Compression strategy (default: "family_directories")

### Statistics Report Configuration
- `SKIP_DETAILED_TAXONOMY`: Skip detailed taxonomic analysis (default: true)
- `STATS_REPORT_FAST_MODE`: Enable fast mode (default: false)
- `STATS_REPORT_SAMPLE_SIZE`: Sample size for analysis (default: 100,000)
- `STATS_REPORT_CPU`: CPU cores (default: 8)
- `STATS_REPORT_MEMORY_MB`: Memory allocation (default: 32,000)
- `STATS_REPORT_RUNTIME`: Time limit in minutes (default: 120)

## Environment Requirements

The pipeline uses conda environments for different steps:
- `create_load_db.yaml`: Database creation and BCDM loading
- `sqlite.yaml`: SQLite operations
- `load_taxonomy.yaml`: Taxonomy processing
- `assess_criteria.yaml`: Standard criteria assessment
- `assess_images.yaml`: Image assessment via CAOS API
- `otu_clustering.yaml`: VSEARCH-based OTU clustering
- `phylogenetic_analysis.yaml`: Phylogenetic analysis with ETE3, MUSCLE, MAFFT, FastTree, IQ-TREE
- `prescoring_filter.yaml`: Pre-filtering operations (when enabled)
- `compress_databases.yaml`: Database compression utilities
- `stats_report.yaml`: Statistical analysis and PDF report generation

## Usage

### Basic Usage (No Pre-filtering, No Phylogenetics)
Run the complete pipeline without pre-filtering, target lists, or phylogenetic analysis:
```bash
snakemake --cores [number_of_cores] --use-conda
```

### With Pre-scoring Filter
Enable pre-filtering for large datasets:
```bash
# First, enable prescoring filter in config/config.yml:
# ENABLE_PRESCORING_FILTER: true
# Configure desired filtering options (taxa, countries, markers, bins)
snakemake --cores [number_of_cores] --use-conda
```

### With Phylogenetic Analysis (Recommended for Complete Analysis)
```bash
# First, enable phylogenetic analysis in config/config.yml:
# PHYLO_ENABLED: true
# Configure phylogenetic parameters as needed
snakemake --cores [number_of_cores] --use-conda
```

### With Target List Filtering
```bash
# First, set USE_TARGET_LIST: true in config/config.yml
# Ensure TARGET_LIST path is specified
snakemake --cores [number_of_cores] --use-conda
```

### Complete Configuration Example (Recommended)
```yaml
# Enable comprehensive analysis with all features:
ENABLE_PRESCORING_FILTER: true
USE_TARGET_LIST: false
PHYLO_ENABLED: true

# Pre-filtering options
FILTER_TAXA: true
FILTER_TAXA_LIST: "resources/target_taxa.txt"
FILTER_SPECIES: true
MARKER: "COI-5P"

# OTU clustering
OTU_CLUSTERING_THRESHOLD: 0.97
OTU_CLUSTERING_THREADS: 16

# Phylogenetic analysis
PHYLO_NUM_JOBS: 100
PHYLO_MIN_OTUS: 4
PHYLO_ALIGNMENT_METHOD: "mafft"
PHYLO_TREE_METHOD: "fasttree"
PHYLO_BOOTSTRAP: 1000
PHYLO_GENERATE_PDFS: true
PHYLO_MAX_CONCURRENT: 20

# Family databases
FAMILY_SIZE_THRESHOLD: 5000
FAMILY_ARRAY_SIZE: 128
WORKERS_PER_JOB: 8
EXPORT_KINGDOMS: ["Animalia", "Plantae"]

# Performance optimization
COMPRESSION_WORKERS: 32
STATS_REPORT_CPU: 16
STATS_REPORT_MEMORY_MB: 64000
```

### High-Performance Computing Usage
```bash
# For large datasets with comprehensive analysis on HPC systems:
snakemake --cores 32 --use-conda --resources mem_mb=128000 \
  --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads} --time=24:00:00" \
  --jobs 50
```

### Clean Intermediate Files
```bash
snakemake clean
```

## Filtering Strategy

The pipeline offers multiple complementary filtering approaches:

### 1. Pre-scoring Filter (Early Stage)
- **Purpose:** Reduce dataset size early in the pipeline for efficiency
- **When to use:** Large datasets that need broad filtering before detailed processing
- **Filters by:** Taxa, geography, genetic markers, BIN characteristics, species-level IDs
- **Advantage:** Reduces computational load for all downstream steps

### 2. Target List Filter (Post-taxonomy)
- **Purpose:** Focus curation on specific species of interest
- **When to use:** Project-specific curation targeting known species lists
- **Filters by:** Species matches against provided target list
- **Advantage:** Precise species-level targeting after full taxonomic processing

### 3. Country Representative Selection (Final Stage)
- **Purpose:** Select optimal representatives per country per species per OTU
- **When to use:** Reducing redundancy while maintaining geographic representation
- **Filters by:** Quality ranking within geographic and phylogenetic groups
- **Advantage:** Balanced representation across geographic regions

## Final Target Rules

### Main Pipeline Target (`rule all`)
The pipeline's main target includes all essential outputs for complete analysis:

```python
rule all:
    input:
        f"{get_results_dir()}/result_output.tsv",                    # Final scored output
        f"{get_results_dir()}/families_split.ok",                   # Family databases
        f"{get_results_dir()}/phylogenetic_results_integrated.ok",  # Phylogenetic integration
        f"{get_results_dir()}/family_databases_compressed.ok",      # Compressed archives
        f"{get_results_dir()}/country_representatives_selected.ok", # Geographic representatives
        f"{get_results_dir()}/stats_report_generated.ok",          # Statistical report
        f"{get_results_dir()}/final_results_archived.ok"           # Final archive package
```

This target ensures:
- **Complete specimen assessment** with all quality criteria and BAGS grades
- **Family-level database organization** for efficient downstream analysis
- **Phylogenetic analysis integration** when enabled
- **Compressed deliverables** optimized for storage and transfer
- **Geographic representation** maintaining sampling diversity
- **Comprehensive documentation** via statistical reports
- **Timestamped final archive** for streamlined delivery

## Output

### Primary Outputs
- **`results/result_output.tsv`**: Final scored and ranked specimens with all assessments
- **`results/family_databases/`**: Family-level SQLite databases for efficient analysis
- **`results/family_databases_compressed/`**: Compressed family directories with all associated files
- **`results/bold_database_statistics.pdf`**: Comprehensive statistical analysis report
- **`results/final_results_[timestamp]/`**: Timestamped archive of key deliverables

### Phylogenetic Analysis Outputs (When Enabled)
- **`results/phylogenies/[family]/`**: Family-level phylogenetic analysis directories containing:
  - Phylogenetic trees in Newick format
  - High-quality PDF tree visualizations with bootstrap support
  - Curation checklist PDFs for Grade C species requiring attention
  - Family analysis summaries and statistics
  - Alignment files and intermediate results

### Assessment Outputs
- Individual criteria assessment files (`assessed_*.tsv`)
- BAGS assessment results with species-level grades
- OTU clustering results with sequence assignments
- Country representative selection results

### Database Structure
The final database includes specialized tables:
- `bold`: Core specimen data
- `bold_criteria`: Quality criteria assessments
- `bold_otus`: OTU clustering results
- `bold_ranks`: Comprehensive ranking scores
- `bags`: Species-level BAGS grades
- `country_representatives`: Selected representatives per country/species/OTU
- `manual_curation`: BOLD portal URLs for manual curation

## Performance Benchmarks

### Typical Processing Times (Without Phylogenetics)
- **Small datasets** (< 10K records): 30-60 minutes
- **Medium datasets** (10K-100K records): 2-6 hours
- **Large datasets** (100K+ records): 6-24 hours
- **Very large datasets** (1M+ records): 1-3 days

### With Phylogenetic Analysis
- **Additional time**: 50-200% increase depending on family sizes and configuration
- **Parallel efficiency**: Scales effectively with available SLURM nodes
- **Typical family processing**: 5-60 minutes per family depending on size and complexity

### Memory Requirements
- **Minimum**: 16GB RAM (32GB recommended with phylogenetics)
- **Recommended**: 32GB+ RAM for large datasets
- **High-performance**: 64GB+ RAM with SSD storage for optimal performance
- **Phylogenetic analysis**: 32GB per job (configurable)
- **Statistics report**: 32GB default (configurable)

### Storage Requirements
- **Raw database**: 0.5-2x input file size
- **Family databases**: 1-1.5x main database size
- **Phylogenetic outputs**: 100-500MB per family (with PDFs)
- **Compressed archives**: 60-80% size reduction
- **Final archive**: Complete deliverable package

### SLURM Configuration Requirements
For phylogenetic analysis and family splitting:
- **Job arrays**: Support for arrays up to 1000+ jobs
- **Memory allocation**: 8-32GB per job depending on phase
- **Time limits**: 4-24 hours per job
- **Concurrent jobs**: 10-50 depending on cluster capacity
- **Partitions**: Access to medium/long partitions for extended analyses

### Optimization Tips
1. **Use pre-scoring filter** for very large datasets to reduce processing time
2. **Adjust `TAXONOMY_CHUNK_SIZE`** based on available system memory
3. **Enable target lists** when focusing on specific species
4. **Configure OTU clustering threads** based on available CPU cores
5. **Use SSD storage** for database operations when possible
6. **Set appropriate `FAMILY_SIZE_THRESHOLD`** for downstream analysis needs
7. **Optimize SLURM parameters** based on cluster specifications and queue times
8. **Enable phylogenetic analysis** for comprehensive taxonomic validation
9. **Use kingdom filtering** to focus phylogenetic analysis on target groups
10. **Configure compression workers** based on available CPU cores for final archiving

## Citation

If you use this pipeline in your research, please cite:

[Citation information will be provided when the pipeline is published]

## Support and Documentation

- **Pipeline issues**: Submit GitHub issues for bug reports and feature requests
- **Configuration help**: Refer to example configuration files in `config/examples/`
- **Performance tuning**: See performance benchmarks and optimization tips above
- **SLURM setup**: Consult your HPC administrator for optimal cluster configuration
