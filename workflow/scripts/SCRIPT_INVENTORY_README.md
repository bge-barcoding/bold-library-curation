# BOLD Library Curation Scripts Inventory

This directory contains all scripts used in the BOLD library curation workflow. There are currently three pipeline variants with different performance characteristics:

## Pipeline Variants

### even-faster
- Uses: `load_taxonomy_faster.pl`, `load_bcdm_fast_simple.pl`
- Snakefile: `Snakefile-even-faster`

### faster  
- Uses: `load_taxonomy.pl`, `load_bcdm_fast_simple.pl`
- Snakefile: `Snakefile-faster`

### standard
- Uses: `load_taxonomy.pl`, `load_bcdm.pl`
- Snakefile: `Snakefile`

## All Scripts in Directory

### Core Pipeline Scripts (Integrated)
- **load_taxonomy.pl** - Loads taxonomic data into database (used in faster & standard pipelines)
- **load_taxonomy_faster.pl** - Optimized version for loading taxonomic data (used in even-faster pipeline)
- **load_bcdm.pl** - Loads BCDM (Barcode Collection Data Model) data (used in standard pipeline)
- **load_bcdm_fast_simple.pl** - Simplified fast version for loading BCDM data (used in faster & even-faster pipelines)

### Assessment and Quality Control Scripts
- **assess_criteria.pl** - Assesses records against specified criteria (e.g., HAS_IMAGE)
- **assess_images.pl** - Validates and assesses image availability and quality for records
- **assess_taxa.pl** - Assesses taxonomic information and BAGS (BIN Analysis and Grouping System) data

### Data Processing and Filtering Scripts
- **autocurate.pl** - Automated curation of sequences based on BAGS ratings and taxonomic criteria
- **filter_output.pl** - Filters and processes output data from draft results
- **partition_result.pl** - Filters and partitions BCDM TSV files based on criteria and columns
- **prescoring_filter.py** - Pre-scoring filter to reduce dataset size before database creation
- **read_snapshot_data.pl** - Extracts and processes data from BOLD snapshot datasets

### Utility and Support Scripts
- **clean.pl** - Cleanup script that removes temporary files (*.tsv, *.ok, *.db) from results directory
- **concat_tsvs.pl** - Concatenates multiple TSV files into a single output file
- **get_statistics.pl** - Generates statistics from filtered information files
- **load_targetlist.pl** - Loads target species lists into the database
- **make_eu_specieslist.pl** - Creates European species lists based on ISO country codes

### Database and SQL Files
- **make_db.sh** - Shell script for database creation
- **schema.sql** - Database schema definition
- **indexes.sql** - Database index definitions
- **pivot.sql** - SQL for data pivoting operations
- **ranking.sql** - SQL for ranking operations
- **sumscore.sql** - SQL for score summation operations

### Documentation
- **README.md** - Existing documentation (primarily for read_snapshot_data.pl)

## Scripts NOT Currently Integrated into Any Pipeline

The following scripts are not referenced in any of the three Snakefiles (even-faster, faster, standard):

1. **assess_criteria.pl**
2. **assess_images.pl**  
3. **assess_taxa.pl**
4. **autocurate.pl**
5. **clean.pl**
6. **concat_tsvs.pl**
7. **filter_output.pl**
8. **get_statistics.pl**
9. **load_targetlist.pl**
10. **make_eu_specieslist.pl**
11. **partition_result.pl**
12. **prescoring_filter.py**
13. **read_snapshot_data.pl**

### Database/SQL files (likely used but not directly called):
- **indexes.sql**
- **make_db.sh** 
- **pivot.sql**
- **ranking.sql**
- **schema.sql**
- **sumscore.sql**

These unintegrated scripts appear to be utility tools for data preprocessing, quality assessment, post-processing analysis, and specialized curation tasks that may be run independently or as part of manual workflows outside the main automated pipelines.
