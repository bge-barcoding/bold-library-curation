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
- **Input:** Database and taxonomy data
- **Output:** TSV file with assessment results
- **Script:** `workflow/scripts/assess_criteria.pl` (except HAS_IMAGE uses `assess_images.pl`)

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

Run the complete pipeline:
```bash
snakemake --cores [number_of_cores] --use-conda
```

Clean intermediate files:
```bash
snakemake clean
```

## Output

The final output (`results/result_output.tsv`) contains filtered and ranked BOLD specimens that meet the quality criteria for library curation, with headers and tab-separated format suitable for further analysis or import into curation systems.

## Logging

All major steps generate log files in the `logs/` directory for troubleshooting and monitoring pipeline execution.
