# BOLD Library Curation Pipeline

This Snakemake workflow processes BOLD sequence data through comprehensive quality assessment criteria and splits results into family-level databases for efficient analysis and curation.

## Pipeline Overview

The pipeline follows a systematic approach to evaluate BOLD specimens across multiple quality criteria, perform advanced phylogenetic analyses, and produce ranked outputs optimized for library curation. The workflow includes optional pre-filtering capabilities, haplotype analysis, OTU clustering, BAGS assessment, and family-level database splitting for scalable downstream analysis.

![BOLD Pipeline DAG](/doc/bold_pipeline_dag.svg)

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
 BOLD TSV file (when filtering enabled)
- **Output:** Filtered BOLD TSV file (when filtering enabled) or marker file only (when disabled)
- **Script:** `workflow/scripts/prescoring_filter.py` (when enabled)

**Available Filtering Options:**
- **Taxa filtering**: Filter by specific taxonomic groups (`FILTER_TAXA: true`, `FILTER_TAXA_LIST`)
- **Country filtering**: Filter by countries of collection (`FILTER_COUNTRIES: true`, `FILTER_COUNTRY_LIST`)
- **Marker filtering**: Filter by specific genetic markers (`MARKER`)
- **BIN sharing**: Enable BIN sharing analysis (`FILTER_BINS: true`)

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

#### Advanced Phylogenetic Analyses
- **HAPLOTYPE_ID**: Identifies unique haplotypes within each BIN and species group
- **OTU_CLUSTERING**: Performs VSEARCH-based clustering to identify Operational Taxonomic Units
  - Configurable similarity threshold (`OTU_CLUSTERING_THRESHOLD`, default: 0.99)
  - Multi-threaded processing (`OTU_CLUSTERING_THREADS`)
  - Temporary file management for large datasets

Each assessment rule:
- **Input:** Database and taxonomy/target list dependency (conditional)
- **Output:** TSV file with assessment results
- **Script:** `workflow/scripts/assess_criteria.pl` (standard criteria), specialized scripts for images, haplotypes, and OTUs
- **Dependency:** Uses `get_taxonomy_dependency()` function to determine assessment dependencies

### Phase 4: BAGS Assessment and Optimization

#### BAGS Database Optimization
**Rule:** `optimize_bags_database`
- Applies BAGS-specific database optimizations for improved performance
- **Optimizations:** WAL mode, memory settings, BAGS-specific indexes
- **Script:** `workflow/scripts/bags_indexes.sql`

#### BAGS Assessment
**Rule:** `BAGS`
- Performs species-level assessment using simplified output format
- **Output format:** 4 essential columns (taxonid, BAGS_grade, BIN_URL, sharers)
- **Performance:** Optimized with progress tracking and simplified column structure
- **Script:** `workflow/scripts/assess_taxa_simplified.pl`

#### BAGS Data Integration
**Rules:** `import_bags`, `inherit_subspecies_bags`
- Imports BAGS results into database
- Inherits BAGS grades for subspecies from parent species
- Enables complex queries combining BAGS with other criteria

### Phase 5: Data Integration and Output

#### Results Consolidation
**Rule:** `concatenate`
- Combines individual criteria assessment results (excluding haplotypes and OTUs)
- **Input:** All standard criterion TSV files
- **Output:** `results/CONCATENATED.tsv`
- **Script:** `workflow/scripts/concat_tsvs.pl`

#### Specialized Data Import
**Rules:** `import_concatenated`, `import_haplotypes`, `import_otus`
- Imports consolidated results, haplotype data, and OTU assignments into database
- Creates specialized tables: `bold_criteria`, `bold_haplotypes`, `bold_otus`
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

#### Final Output Generation
**Rule:** `output_filtered_data`
- Generates final scored and ranked output with all assessments including OTUs
- **Output:** `results/result_output.tsv`
- **Script:** `workflow/scripts/ranking_with_stored_ranks_otu.sql`

### Phase 6: Family-Level Database Creation

#### Family Database Splitting
**Rule:** `split_families`
- Splits main database into family-level databases for efficient analysis
- **Configuration:** `FAMILY_SIZE_THRESHOLD` (default: 10,000 records)
- **Organization:** Hierarchical structure by phylum/family
- **Output:** Individual SQLite databases per family or combined small families
- **Script:** `workflow/scripts/bold_family_splitter.py`

#### Pipeline Summary
**Rule:** `create_final_summary`
- Generates comprehensive pipeline execution summary
- **Output:** `results/pipeline_summary.txt`
- Includes processing statistics, family database information, and file structure overview

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
- `FAMILY_SIZE_THRESHOLD`: Minimum records for individual family database (default: 10,000)

## Environment Requirements

The pipeline uses conda environments for different steps:
- `create_load_db.yaml`: Database creation and BCDM loading
- `sqlite.yaml`: SQLite operations
- `load_taxonomy.yaml`: Taxonomy processing
- `assess_criteria.yaml`: Standard criteria assessment
- `assess_images.yaml`: Image assessment via CAOS API
- `haplotype_analysis.yaml`: Haplotype identification
- `otu_clustering.yaml`: VSEARCH-based OTU clustering
- `prescoring_filter.yaml`: Pre-filtering operations (when enabled)

## Usage

### Basic Usage (No Pre-filtering)
Run the complete pipeline without pre-filtering or target lists:
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

### With Target List Filtering
```bash
# First, set USE_TARGET_LIST: true in config/config.yml
# Ensure TARGET_LIST path is specified
snakemake --cores [number_of_cores] --use-conda
```

### Combined Configuration Example
```yaml
# Enable both pre-filtering and target lists for maximum efficiency:
ENABLE_PRESCORING_FILTER: true
USE_TARGET_LIST: true
FILTER_TAXA: true
FILTER_TAXA_LIST: "resources/target_taxa.txt"
MARKER: "COI-5P"
OTU_CLUSTERING_THRESHOLD: 0.97
FAMILY_SIZE_THRESHOLD: 5000
```

### Performance Tuning
```bash
# For large datasets with performance optimization:
snakemake --cores 16 --use-conda --resources mem_mb=32000
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
- **Filters by:** Taxa, geography, genetic markers, BIN characteristics
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

## Output

### Primary Outputs
- **`results/result_output.tsv`**: Final scored and ranked specimens with all assessments
- **`results/family_databases/`**: Family-level SQLite databases for efficient analysis
- **`results/pipeline_summary.txt`**: Comprehensive execution summary

### Assessment Outputs
- Individual criteria assessment files (`assessed_*.tsv`)
- BAGS assessment results with species-level grades
- Haplotype identification and OTU clustering results
- Country representative selection results

### Database Structure
The final database includes specialized tables:
- `bold`: Core specimen data
- `bold_criteria`: Quality criteria assessments
- `bold_haplotypes`: Haplotype assignments
- `bold_otus`: OTU clustering results
- `bold_ranks`: Comprehensive ranking scores
- `bags`: Species-level BAGS grades
- `country_representatives`: Selected representatives per country/species/OTU

## Performance Benchmarks

### Typical Processing Times
- **Small datasets** (< 10K records): 30-60 minutes
- **Medium datasets** (10K-100K records): 2-6 hours
- **Large datasets** (100K+ records): 6-24 hours
- **Very large datasets** (1M+ records): 1-3 days

### Memory Requirements
- **Minimum**: 8GB RAM
- **Recommended**: 16GB+ RAM for large datasets
- **High-performance**: 32GB+ RAM with SSD storage
- **Taxonomy chunk size**: Configurable based on available memory

### Optimization Tips
1. **Use pre-scoring filter** for very large datasets to reduce processing time
2. **Adjust `TAXONOMY_CHUNK_SIZE`** based on available system memory
3. **Enable target lists** when focusing on specific species
4. **Configure OTU clustering threads** based on available CPU cores
5. **Use SSD storage** for database operations when possible
6. **Set appropriate `FAMILY_SIZE_THRESHOLD`** for downstream analysis needs

## Troubleshooting

### Common Issues
1. **Memory errors during taxonomy loading**: Reduce `TAXONOMY_CHUNK_SIZE` in config
2. **OTU clustering failures**: Check VSEARCH installation and reduce thread count
3. **BAGS assessment slowdown**: Ensure database optimization completed successfully
4. **Family splitting errors**: Verify adequate disk space and write permissions
5. **Missing dependencies**: Verify all conda environments are properly installed

### Debug Mode
Enable debug logging by setting `LOG_LEVEL: "DEBUG"` in config for detailed execution information.

### Performance Monitoring
- Check individual log files for step-specific performance metrics
- Monitor disk space during family database creation
- Use system monitoring tools during resource-intensive operations (BAGS, OTU clustering)

## Version History

### v4.0 - Comprehensive Analysis Pipeline (Current)
- Added OTU clustering with VSEARCH integration
- Enhanced haplotype analysis with specialized database tables
- Country representative selection for geographic balance
- Family-level database splitting for scalable analysis
- Comprehensive ranking system with stored ranks
- Advanced performance optimizations throughout pipeline

### v3.0 - Enhanced BAGS Processing
- Optimized BAGS assessment with simplified output format
- BAGS database integration for complex querying
- Subspecies grade inheritance
- Performance optimizations for large datasets

### v2.0 - Dynamic File Handling
- Intelligent file handling based on filtering configuration
- Dynamic input selection functions
- Resource optimization for HPC environments

### v1.0 - Initial Release
- Core assessment pipeline with 16 quality criteria
- Basic BAGS assessment
- Pre-scoring filter capabilities
- Target list filtering

## Citation

If you use this pipeline in your research, please cite:
[Citation information to be added based on publication]

## Support

For issues, questions, or contributions, please:
1. Check the troubleshooting section above
2. Review log files in the `logs/` directory for specific error messages
3. Verify configuration settings match your data and system requirements
4. Open an issue in the project repository with relevant log excerpts
5. Contact the development team for complex issues

## License

[License information to be added]
