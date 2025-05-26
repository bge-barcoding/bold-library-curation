# BOLD Library Curation Pipeline

This Snakemake workflow processes BOLD records to assess and filter specimens based on multiple quality criteria for library curation purposes.

## Pipeline Overview

The pipeline follows a systematic approach to evaluate BOLD specimens across 16 different quality criteria, producing a ranked output of specimens suitable for subsequent manual library curation.

## Workflow Steps

### 1. Database Creation and Initial Data Loading
**Rule:** `create_load_db`
- Creates SQLite database from BOLD TSV data
- Loads the Barcode Core Data Model (BCDM) schema
- **Input:** BOLD TSV file, database schema
- **Output:** SQLite database file
- **Script:** `workflow/scripts/load_bcdm.pl`

### 2. Criteria Configuration Loading
**Rule:** `load_criteria`
- Imports assessment criteria definitions into the database
- **Input:** `resources/criteria.tsv`
- **Output:** Database with criteria table loaded

### 3. Database Indexing
**Rule:** `apply_indexes`
- Applies database indexes for improved query performance
- **Input:** Index definitions SQL file
- **Output:** Indexed database file

### 4. Taxonomy Loading
**Rule:** `load_taxonomy`
- Loads taxonomic information into the database
- Enriches specimen records with taxonomic hierarchy
- **Script:** `workflow/scripts/load_taxonomy.pl`

### 4.5. Target List Import (Optional)
**Rules:** `import_target_list` or `skip_target_list`
- **Conditional step**: Only runs if `USE_TARGET_LIST: true` in config
- When enabled: Filters specimens to only target species from provided list
- When disabled: Creates pass-through dependency for all specimens
- **Input:** Target species CSV file (when enabled)
- **Script:** `workflow/scripts/load_targetlist.pl`
- **Purpose:** Allows focused curation on specific species of interest

### 5. Quality Criteria Assessment
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

### 6. Results Consolidation
**Rule:** `concatenate`
- Combines all individual criterion assessment results
- **Input:** All 16 criterion TSV files
- **Output:** `results/CONCATENATED.tsv`
- **Script:** `workflow/scripts/concat_tsvs.pl`

### 7. Results Import
**Rule:** `import_concatenated`
- Imports consolidated results back into the database
- Creates `bold_criteria` table with all assessments

### 8. Final Output Generation
**Rule:** `output_filtered_data`
- Applies ranking algorithm to prioritize specimens
- Generates final filtered and ranked specimen list
- **Output:** `results/result_output.tsv`
- **Script:** `workflow/scripts/ranking.sql`

## Configuration

The pipeline is configured through `config/config.yml` which defines:
- Input file paths (BOLD TSV, schema, indexes)
- Database file locations
- Logging levels
- Library paths

## Environment Requirements

The pipeline uses conda environments for different steps:
- `create_load_db.yaml`: Database creation and BCDM loading
- `sqlite.yaml`: SQLite operations
- `load_taxonomy.yaml`: Taxonomy processing
- `assess_criteria.yaml`: Criteria assessment
- `assess_images.yaml`: Image assessment

## Usage

Run the complete pipeline (without target list):
```bash
snakemake --cores [number_of_cores] --use-conda
```

Run with target list filtering enabled:
```bash
# First, set USE_TARGET_LIST: true in config/config.yml
snakemake --cores [number_of_cores] --use-conda
```

Clean intermediate files:
```bash
snakemake clean
```

## Target List Functionality

The pipeline can optionally filter specimens to only those matching a target species list:

**To enable target list filtering:**
1. Set `USE_TARGET_LIST: true` in `config/config.yml`
2. Ensure your target species CSV file is specified in `TARGET_LIST`
3. Run the pipeline normally

**Behavior differences:**
- **Without target list** (`USE_TARGET_LIST: false`): Processes all specimens in BOLD dataset
- **With target list** (`USE_TARGET_LIST: true`): Only processes specimens matching target species, making the pipeline more efficient for focused curation projects

## Output

The final output (`results/result_output.tsv`) contains filtered and ranked BOLD specimens that meet the quality criteria for library curation, with headers and tab-separated format suitable for further analysis or import into curation systems.

## Logging

All major steps generate log files in the `logs/` directory for troubleshooting and monitoring pipeline execution.
