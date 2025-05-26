# Configuration Guide

This document explains all configuration options in `config.yml` for the BOLD Library Curation Pipeline.

## File Paths and Database Configuration

### Input Data Files
- **`BOLD_TSV`**: Path to the main BOLD public data file in TSV format
  - Default: `resources/BOLD_Public.05-Apr-2024.tsv`
  - This is the primary dataset containing all BOLD specimen records

### Database Configuration
- **`DB_FILE`**: Path for the main SQLite database file
  - Default: `results/bold.db`
  - Contains all loaded BOLD data and assessment results

- **`SCHEMA`**: SQL schema file for database structure
  - Default: `workflow/scripts/schema.sql`
  - Defines the Barcode Commons Data Model (BCDM) table structure

- **`INDEXES`**: SQL file containing database index definitions
  - Default: `workflow/scripts/indexes.sql`
  - Improves query performance for large datasets

- **`DB_FILE_INDEXED`**: Flag file indicating database indexing completion
  - Default: `results/bold_indexed.ok`
  - Touch file created after successful index application

## Library and Path Configuration

- **`LIBS`**: Directory containing Perl library modules
  - Default: `lib/`
  - Contains custom Perl modules for BCDM operations

- **`PATH`**: Additional path for library access
  - Default: `lib/`
  - Used for Perl module inclusion

## Taxonomy Configuration

- **`TAXONOMY_LOADED`**: Output file for taxonomy loading verification
  - Default: `results/taxonomy_check.tsv`
  - Contains taxonomy validation results

- **`KINGDOM`**: Taxonomic kingdom filter
  - Default: `"Animalia"`
  - Restricts processing to specific taxonomic kingdom

- **`TAXON_LEVEL`**: Taxonomic resolution level
  - Default: `"species"`
  - Defines the taxonomic level for analysis

## Assessment Criteria Configuration

- **`CRITERIA`**: Space-separated list of assessment criteria
  - Default: All 16 quality criteria
  - Complete list:
    - `COLLECTION_DATE`: Collection date completeness
    - `COLLECTORS`: Collector information quality
    - `COORD`: Geographic coordinate validation
    - `COUNTRY`: Country information accuracy
    - `REGION`: Regional location data
    - `SECTOR`: Collection sector/zone information
    - `HAS_IMAGE`: Specimen image availability
    - `IDENTIFIER`: Taxonomic identifier quality
    - `ID_METHOD`: Identification method documentation
    - `INSTITUTION`: Institutional affiliation
    - `MUSEUM_ID`: Museum collection identifiers
    - `PUBLIC_VOUCHER`: Public voucher availability
    - `SEQ_QUALITY`: DNA sequence quality metrics
    - `SITE`: Collection site specificity
    - `SPECIES_ID`: Species identification confidence
    - `TYPE_SPECIMEN`: Type specimen designation

- **`DB_CRITERIA_ADDED`**: Flag file for criteria loading completion
  - Default: `results/criteria_indexed`
  - Indicates successful criteria integration

## Project Configuration

- **`PROJECT_NAME`**: Identifier for the curation project
  - Default: `"bold-curation"`
  - Used for project tracking and output labeling

## Target List Configuration

- **`TARGET_LIST`**: CSV file containing target species list
  - Default: `resources/all_specs_and_syn.csv`
  - Contains species and synonyms for targeted curation

- **`USE_TARGET_LIST`**: Enable/disable target list filtering
  - Default: `false`
  - Set to `true` to filter specimens to only target species
  - When `false`: processes all specimens in the dataset
  - When `true`: only processes specimens matching target species list

## System Configuration

- **`LOG_LEVEL`**: Logging verbosity level
  - Default: `"INFO"`
  - Options: `DEBUG`, `INFO`, `WARN`, `ERROR`
  - Controls the amount of diagnostic information generated

## Customization Guidelines

### Modifying Input Data
To use a different BOLD dataset:
```yaml
BOLD_TSV: path/to/your/bold_data.tsv
```

### Adjusting Assessment Criteria
To assess only specific criteria, modify the CRITERIA string:
```yaml
CRITERIA: "COLLECTION_DATE COLLECTORS COORD COUNTRY SEQ_QUALITY"
```

### Changing Taxonomic Scope
To focus on different taxonomic groups:
```yaml
KINGDOM: "Plantae"
TAXON_LEVEL: "genus"
```

### Database Location
For different output locations:
```yaml
DB_FILE: /path/to/custom/database.db
```

### Logging Configuration
For detailed debugging:
```yaml
LOG_LEVEL: "DEBUG"
```

For minimal output:
```yaml
LOG_LEVEL: "ERROR"
```

## File Dependencies

Ensure these files exist before running the pipeline:
- Input TSV file specified in `BOLD_TSV`
- Schema file specified in `SCHEMA`
- Index file specified in `INDEXES`
- Library directory specified in `LIBS`

## Output Locations

The pipeline will create these output directories and files:
- `results/` - All output files and databases
- `logs/` - Log files for each pipeline step
- Files specified in `DB_FILE`, `DB_FILE_INDEXED`, `TAXONOMY_LOADED`

## Notes

- All relative paths are resolved from the workflow root directory
- The configuration supports both absolute and relative file paths
- Some configuration options (like `TARGET_LIST`) are prepared for future features
- Modify `LOG_LEVEL` to `DEBUG` for troubleshooting pipeline issues
