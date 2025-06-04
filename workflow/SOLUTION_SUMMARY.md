# BOLD Parallel Family Splitting Solution - CORRECTED

## Overview
This solution creates individual SQLite database files for each family with a single table called "records" containing the complete processed BOLD data with all quality assessments, rankings, and metadata - exactly matching the structure of your example CSV.

## Key Correction Made
The original implementation was incorrectly creating databases with just raw data from the `bold` table. The corrected version uses the complete ranking query from `ranking_with_stored_ranks_otu.sql` to ensure each family database contains:

- All quality assessment scores (SPECIES_ID, TYPE_SPECIMEN, SEQ_QUALITY, etc.)
- BAGS grades and rankings
- Manual curation data
- Haplotype and OTU IDs
- Complete taxonomic and specimen metadata

## File Structure
Each family database will have:
- **Single table**: `records`
- **Same columns**: Exactly matching your `records.csv` example (100+ columns)
- **Complete data**: All processed results, not just raw BOLD data

## Components Created

### Core Processing Scripts
1. **`prepare_family_batches.py`** - Analyzes families and creates batches for parallel processing
2. **`process_family_batch.py`** - Multi-threaded family processor (CORRECTED)
3. **`family_split_array.sh`** - SLURM job array script
4. **`consolidate_results.py`** - Results collection and reporting

### Integration Files
5. **`new_split_families_rule.smk`** - Updated Snakemake rule (replaces existing one)
6. **`parallel_splitting_config.yml`** - Configuration options for config.yaml
7. **`setup_parallel_splitting.sh`** - Setup script

## Configuration Options (add to config.yaml)
```yaml
# Family splitting configuration
FAMILY_SIZE_THRESHOLD: 10000    # Size threshold for subfamily splitting
FAMILY_ARRAY_SIZE: 64          # Number of parallel SLURM jobs
WORKERS_PER_JOB: 4             # Threads per job
JOB_MEMORY: "8G"               # Memory per job
JOB_TIME: "04:00:00"           # Time limit per job
```

## Output Structure
```
results/family_databases/
├── Animalia/
│   ├── Arthropoda/
│   │   ├── Insecta/
│   │   │   ├── Coleoptera/
│   │   │   │   ├── Tenebrionidae/
│   │   │   │   │   ├── Tenebrionidae.db           # Small family
│   │   │   │   │   ├── Tenebrionidae_Lagriinae.db # Large family split
│   │   │   │   │   └── Tenebrionidae_Alleculinae.db
```

## Database Content
Each `.db` file contains:
- **Table name**: `records`
- **Columns**: All 100+ columns from your CSV example
- **Data**: Complete processed BOLD records with rankings and assessments

## Integration Steps
1. Run `setup_parallel_splitting.sh` to make scripts executable
2. Add configuration from `parallel_splitting_config.yml` to your `config/config.yml`
3. Replace the existing `split_families` rule in `bold-ranker.smk` with the new rule from `new_split_families_rule.smk`

## Performance
- **Scalability**: 32-128+ cores based on configuration
- **Memory**: Streams data, loads individual families only
- **Speed**: Should process 3M records in 2-4 hours on HPC
- **Error handling**: Automatic retry with comprehensive logging

## Testing
Start with a small test database to validate the approach before running on the full 15GB dataset.

The corrected solution now ensures each family database contains the complete processed data structure you need, with a single `records` table matching your CSV example exactly.
