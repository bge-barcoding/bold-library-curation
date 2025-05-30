# BOLD Library Curation Pipeline

This Snakemake workflow processes BOLD records to assess and filter specimens based on multiple quality criteria for library curation purposes.

## Pipeline Overview

The pipeline follows a systematic approach to evaluate BOLD specimens across 16 different quality criteria, producing a ranked output of specimens suitable for subsequent manual library curation. The workflow includes optional pre-filtering capabilities to reduce dataset size before detailed assessment.

![BOLD Pipeline DAG](/doc/bold_pipeline_dag.svg)


## Workflow Steps

### 1. Pre-scoring Filter (Optional)
**Rules:** `prescoring_filter` or `skip_prescoring_filter`
- **Optional step**: Filters the input BOLD dataset before detailed processing
- Can filter by taxa, countries, and/or BIN sharing criteria
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
- **BIN sharing**: Enable BIN sharing analysis (`FILTER_BINS: true`)

### 2. Database Creation and Initial Data Loading
**Rule:** `create_load_db`
- Creates SQLite database from BOLD TSV data
- Uses dynamic input selection: filtered file (when filtering enabled) or original file (when filtering disabled)
- Loads the Barcode Core Data Model (BCDM) schema
- **Input:** BOLD TSV file (selected automatically based on filtering configuration), database schema
- **Output:** SQLite database file
- **Script:** `workflow/scripts/load_bcdm.pl`
- **Dependency:** Requires prescoring filter step completion (marker file)

### 3. Criteria Configuration Loading
**Rule:** `load_criteria`
- Imports assessment criteria definitions into the database
- **Input:** `resources/criteria.tsv`
- **Output:** Database with criteria table loaded

### 4. Database Indexing
**Rule:** `apply_indexes`
- Applies database indexes for improved query performance
- **Input:** Index definitions SQL file
- **Output:** Indexed database file

### 5. Taxonomy Loading
**Rule:** `load_taxonomy`
- Loads taxonomic information into the database
- Enriches specimen records with taxonomic hierarchy
- **Script:** `workflow/scripts/load_taxonomy.pl`

### 6. Target List Import (Optional)
**Rules:** `import_target_list` or `skip_target_list`
- **Conditional step**: Only runs if `USE_TARGET_LIST: true` in config
- When enabled: Filters specimens to only target species from provided list
- When disabled: Creates pass-through dependency for all specimens
- **Input:** Target species CSV file (when enabled)
- **Script:** `workflow/scripts/load_targetlist.pl`
- **Purpose:** Allows focused curation on specific species of interest

### 7. Quality Criteria Assessment
The pipeline evaluates specimens against 16 different quality criteria. Each criterion is assessed independently:

#### Specimen Metadata Criteria
- **COLLECTION_DATE**: Validates collection date present
- **COLLECTORS**: Assesses collector information present
- **IDENTIFIER**: Evaluates taxonomic identifier (BIN match vs person)
- **ID_METHOD**: Checks identification method (BIN match vs morphology)

#### Geographic Information Criteria
- **COUNTRY**: Validates country present
- **REGION**: Validates region present
- **SITE**: Validates site present
- **SECTOR**: Validates sector present
- **COORD**: Validates coordinates present

#### Institutional and Repository Criteria
- **INSTITUTION**: Assesses institutional affiliation
- **MUSEUM_ID**: Validates museum/collection identifier present
- **PUBLIC_VOUCHER**: Checks voucher specimen in public institution

#### Specimen Quality Criteria
- **SEQ_QUALITY**: Evaluates DNA sequence quality metrics
- **SPECIES_ID**: Validates species-level identification provided
- **TYPE_SPECIMEN**: Identifies type specimen status across multiple BCDM fields
- **HAS_IMAGE**: Checks for associated specimen images using CAOS api

Each assessment rule:
- **Input:** Database and taxonomy/target list dependency (conditional)
- **Output:** TSV file with assessment results
- **Script:** `workflow/scripts/assess_criteria.pl` (except HAS_IMAGE uses `assess_images.pl`)
- **Dependency:** Uses `get_taxonomy_dependency()` function to determine if specimens should be assessed after taxonomy loading alone or after target list filtering

### 8. Results Consolidation
**Rule:** `concatenate`
- Combines all individual criterion assessment results
- **Input:** All 16 criterion TSV files
- **Output:** `results/CONCATENATED.tsv`
- **Script:** `workflow/scripts/concat_tsvs.pl`

### 9. Results Import
**Rule:** `import_concatenated`
- Imports consolidated results back into the database
- Creates `bold_criteria` table with all assessments

### 10. Final Output Generation
**Rule:** `output_filtered_data`
- Applies ranking algorithm to prioritize specimens
- Generates final filtered and ranked specimen list
- **Output:** `results/result_output.tsv`
- **Script:** `workflow/scripts/ranking.sql`

## Configuration

The pipeline is configured through `config/config.yml` which defines:

### Core Configuration
- Input file paths (BOLD TSV, schema, indexes)
- Database file locations
- Logging levels (`LOG_LEVEL`)
- Library paths (`LIBS`)

### Pre-scoring Filter Configuration
- `ENABLE_PRESCORING_FILTER`: Enable/disable pre-filtering (default: false)
  - **When true**: Runs filtering script to create filtered dataset
  - **When false**: Skips filtering entirely, uses original file directly (no file copying)
- `PRESCORING_FILTERED_OUTPUT`: Output path for filtered data (only used when filtering enabled)
- `FILTER_TAXA`: Enable taxonomic filtering (requires `FILTER_TAXA_LIST`)
- `FILTER_TAXA_LIST`: Path to taxa list file
- `FILTER_COUNTRIES`: Enable country filtering (requires `FILTER_COUNTRY_LIST`)
- `FILTER_COUNTRY_LIST`: Path to countries list file
- `FILTER_BINS`: Enable BIN sharing analysis

### Target List Configuration
- `USE_TARGET_LIST`: Enable/disable target species filtering
- `TARGET_LIST`: Path to target species CSV file
- `PROJECT_NAME`: Project identifier for target list
- `TAXON_LEVEL`: Taxonomic level for target matching
- `KINGDOM`: Kingdom scope for target filtering

## Environment Requirements

The pipeline uses conda environments for different steps:
- `create_load_db.yaml`: Database creation and BCDM loading
- `sqlite.yaml`: SQLite operations
- `load_taxonomy.yaml`: Taxonomy processing
- `assess_criteria.yaml`: Criteria assessment
- `assess_images.yaml`: Image assessment
- `prescoring_filter.yaml`: Pre-filtering operations (when enabled)

## Usage

### Basic Usage (No Pre-filtering)
Run the complete pipeline without pre-filtering or target lists (most efficient for small-medium datasets):
```bash
snakemake --cores [number_of_cores] --use-conda
```
This configuration uses the original BOLD TSV file directly without any file copying or filtering overhead.

### With Pre-scoring Filter
Enable pre-filtering for large datasets that benefit from early-stage reduction:
```bash
# First, enable prescoring filter in config/config.yml:
# ENABLE_PRESCORING_FILTER: true
# Configure desired filtering options (taxa, countries, bins)
snakemake --cores [number_of_cores] --use-conda
```

### With Target List Filtering
```bash
# First, set USE_TARGET_LIST: true in config/config.yml
# Ensure TARGET_LIST path is specified
snakemake --cores [number_of_cores] --use-conda
```

### Combined Pre-filtering and Target Lists
```bash
# Enable both in config/config.yml:
# ENABLE_PRESCORING_FILTER: true
# USE_TARGET_LIST: true
# Configure all relevant filtering parameters
snakemake --cores [number_of_cores] --use-conda
```

### Clean Intermediate Files
```bash
snakemake clean
```

## Filtering Strategy

The pipeline offers two complementary filtering approaches with intelligent resource optimization:

### 1. Pre-scoring Filter (Early Stage)
- **Purpose:** Reduce dataset size early in the pipeline for efficiency
- **When to use:** Large datasets that need broad filtering before detailed processing
- **Filters by:** Taxa, geography, BIN characteristics
- **Advantage:** Reduces computational load for all downstream steps
- **Optimization:** When disabled, no file operations occur - original file is used directly

### 2. Target List Filter (Post-taxonomy)
- **Purpose:** Focus curation on specific species of interest
- **When to use:** Project-specific curation targeting known species lists
- **Filters by:** Species matches against provided target list
- **Advantage:** Precise species-level targeting after full taxonomic processing

### Resource Optimization
The pipeline automatically optimizes resource usage based on configuration:
- **No filtering**: Direct use of original file, no copying overhead
- **Pre-filtering only**: Creates filtered dataset, downstream steps use filtered file
- **Target filtering only**: All specimens processed, then filtered to target species
- **Combined filtering**: Maximum efficiency for focused curation projects

### Combined Strategy
Both filters can be used together for maximum efficiency:
1. Pre-scoring filter reduces initial dataset size
2. Target list filter provides species-specific focus
3. Result: Highly focused dataset optimized for targeted curation projects

## Output

The final output (`results/result_output.tsv`) contains filtered and ranked BOLD specimens that meet the quality criteria for library curation, with headers and tab-separated format suitable for further analysis or import into curation systems.

## Logging

All major steps generate log files in the `logs/` directory for troubleshooting and monitoring pipeline execution. Key log files include:
- `prescoring_filter.log`: Pre-filtering operations
- `create_load_db.log`: Database creation and loading
- `load_taxonomy.log`: Taxonomy processing
- `load_target_list.log`: Target list processing (when enabled)
- Individual assessment logs for each criterion
- `output_filtered_data.log`: Final output generation

## Pipeline Dependencies

The workflow uses a smart dependency system that adapts based on configuration:
- **Without target list:** Assessment steps depend on `taxonomy_loaded.ok`
- **With target list:** Assessment steps depend on `target_loaded.ok`
- **Pre-scoring filter enabled:** Filtering rule runs before database creation
- **Pre-scoring filter disabled:** Skip rule creates marker only, database uses original file directly
- **Helper functions:**
  - `get_taxonomy_dependency()`: Automatically determines assessment dependencies
  - `get_input_file()`: Automatically selects appropriate input file (original or filtered)

## Recent Optimizations

### File Handling Efficiency (v2.0)
The pipeline now intelligently handles file operations based on filtering configuration:
- **When filtering disabled**: No file copying occurs, original file used directly
- **When filtering enabled**: Filtered file created and used for downstream processing  
- **Dynamic input selection**: Automatic selection of appropriate input file throughout pipeline
- **Resource savings**: Eliminates unnecessary file operations for improved performance on HPC systems
### Enhanced BAGS Processing (v3.0 - Even Faster BAGS)
- **Faster Species Assessment**: Optimized `assess_taxa.pl` script for improved BAGS evaluation performance
- **Database Integration**: BAGS results now imported into database for complex querying capabilities
- **Advanced Ranking**: Enhanced `ranking_with_sumscore.sql` provides comprehensive scoring combining specimen and species-level data
- **Improved Dependency Management**: BAGS import integrated into workflow dependencies for consistent data availability

### Performance Improvements
- **Enhanced Taxonomy Loading**: Faster `load_taxonomy_faster.pl` with configurable chunking for memory optimization
- **Optimized Database Operations**: Improved indexing and import processes throughout pipeline
- **Smart Dependency Resolution**: Dynamic dependency management based on configuration options reduces unnecessary processing

### File Handling Efficiency (v2.0+)
The pipeline intelligently handles file operations based on filtering configuration:
- **When filtering disabled**: No file copying occurs, original file used directly
- **When filtering enabled**: Filtered file created and used for downstream processing  
- **Dynamic input selection**: Automatic selection of appropriate input file throughout pipeline
- **Resource savings**: Eliminates unnecessary file operations for improved performance on HPC systems

## Performance Benchmarks

### Typical Processing Times (Approximate)
- **Small datasets** (< 10K records): 15-30 minutes
- **Medium datasets** (10K-100K records): 1-3 hours
- **Large datasets** (100K+ records): 3-8 hours

### Memory Requirements
- **Minimum**: 4GB RAM
- **Recommended**: 8GB+ RAM for large datasets
- **Taxonomy chunk size**: Configurable based on available memory

### Optimization Tips
1. **Use pre-scoring filter** for very large datasets to reduce processing time
2. **Adjust `TAXONOMY_CHUNK_SIZE`** based on available system memory
3. **Enable target lists** when focusing on specific species to reduce assessment overhead
4. **Use SSD storage** for database operations when possible

## Troubleshooting

### Common Issues
1. **Memory errors during taxonomy loading**: Reduce `TAXONOMY_CHUNK_SIZE` in config
2. **Slow database operations**: Ensure adequate disk space and consider SSD storage
3. **Missing dependencies**: Verify all conda environments are properly installed
4. **Filtering errors**: Check that taxa and country list files exist and are properly formatted

### Debug Mode
Enable debug logging by setting `LOG_LEVEL: "DEBUG"` in config for detailed execution information.

## Version History

### v3.0 - Even Faster BAGS
- Enhanced BAGS assessment with optimized `assess_taxa.pl`
- BAGS database integration for improved querying
- Advanced ranking with `ranking_with_sumscore.sql`
- Performance optimizations throughout pipeline

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
2. Review log files in the `logs/` directory
3. Open an issue in the project repository
4. Contact the development team

## License

[License information to be added]
