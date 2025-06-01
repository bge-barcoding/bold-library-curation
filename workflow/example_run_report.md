# BOLD Pipeline Example Run - Output Files Report

## Pipeline Configuration
- **Input File**: `resources/test_data/test_odonata_large.tsv`
- **Results Directory**: `results/` (configurable via `RESULTS_DIR`)
- **Logs Directory**: `logs/` (configurable via `LOG_DIR`)
- **Pre-scoring Filter**: Enabled
- **Target List**: Disabled
- **Family Size Threshold**: 10,000 records

## Expected Output Files

### Phase 1: Data Preparation
```
results/
├── prescoring_filter.ok                    # Filter completion marker
└── test_odonata_large_filtered.tsv        # Filtered input data

logs/
└── prescoring_filter.log                  # Pre-filtering process log
```

### Phase 2: Database Setup
```
results/
├── bold.db                                 # Main SQLite database
├── criteria_loaded.ok                     # Criteria loading marker
├── bold_indexed.ok                        # Database indexing marker
├── taxonomy_loaded.ok                     # Taxonomy loading marker
└── target_loaded.ok                       # Target list marker (passthrough)

logs/
├── create_load_db.log                     # Database creation log
├── load_criteria.log                      # Criteria loading log
├── apply_indexes.log                      # Indexing log
└── load_taxonomy.log                      # Taxonomy loading log
```

### Phase 3: Quality Criteria Assessment
```
results/
├── assessed_COLLECTION_DATE.tsv           # Collection date assessment
├── assessed_COLLECTORS.tsv                # Collectors information assessment
├── assessed_COUNTRY.tsv                   # Country information assessment
├── assessed_ID_METHOD.tsv                 # Identification method assessment
├── assessed_IDENTIFIER.tsv                # Identifier information assessment
├── assessed_INSTITUTION.tsv               # Institution affiliation assessment
├── assessed_COORD.tsv                     # Geographic coordinates assessment
├── assessed_MUSEUM_ID.tsv                 # Museum ID assessment
├── assessed_PUBLIC_VOUCHER.tsv            # Public voucher assessment
├── assessed_SEQ_QUALITY.tsv               # Sequence quality assessment
├── assessed_SITE.tsv                      # Collection site assessment
├── assessed_REGION.tsv                    # Geographic region assessment
├── assessed_SECTOR.tsv                    # Geographic sector assessment
├── assessed_SPECIES_ID.tsv                # Species identification assessment
├── assessed_TYPE_SPECIMEN.tsv             # Type specimen assessment
├── assessed_HAS_IMAGE.tsv                 # Image availability assessment
└── assessed_HAPLOTYPE_ID.tsv              # Haplotype identification assessment

logs/
├── assess_COLLECTION_DATE.log
├── assess_COLLECTORS.log
├── assess_COUNTRY.log
├── assess_ID_METHOD.log
├── assess_IDENTIFIER.log
├── assess_INSTITUTION.log
├── assess_COORD.log
├── assess_MUSEUM_ID.log
├── assess_PUBLIC_VOUCHER.log
├── assess_SEQ_QUALITY.log
├── assess_SITE.log
├── assess_REGION.log
├── assess_SECTOR.log
├── assess_SPECIES_ID.log
├── assess_TYPE_SPECIMEN.log
├── assess_HAS_IMAGE.log
└── assess_HAPLOTYPE_ID.log
```

### Phase 4: BAGS Assessment
```
results/
├── bags_optimized.ok                      # BAGS database optimization marker
├── assessed_BAGS.tsv                      # BAGS assessment results
├── bags_imported.ok                       # BAGS import completion marker
├── haplotypes_assigned.ok                 # Haplotype assignment marker
└── subspecies_bags_inherited.ok           # Subspecies inheritance marker

logs/
├── optimize_bags_database.log             # BAGS optimization log
├── assess_BAGS.log                        # BAGS assessment log
├── import_bags.log                        # BAGS import log
└── inherit_subspecies_bags.log            # Subspecies inheritance log
```

### Phase 5: Data Integration
```
results/
├── CONCATENATED.tsv                       # Combined criteria assessments
├── concatenated_imported.ok               # Import completion marker
└── result_output.tsv                      # Final scored and ranked output

logs/
├── import_concatenated.log                # Concatenation import log
└── output_filtered_data.log               # Final output generation log
```

### Phase 6: Family-Level Databases
```
results/
├── families_split.ok                      # Family splitting completion marker
├── pipeline_summary.txt                   # Comprehensive pipeline summary
└── family_databases/                      # Family-level database directory
    ├── splitting_report.txt               # Detailed splitting statistics
    ├── Arthropoda/                        # Phylum-level organization
    │   ├── Odonata_Aeshnidae.db          # Family database (>threshold)
    │   ├── Odonata_Libellulidae.db       # Family database (>threshold)
    │   └── Odonata_small_families.db     # Combined small families (<threshold)
    └── [other_phyla]/                     # Additional phyla as needed

logs/
└── split_families.log                     # Family splitting process log
```

## File Descriptions

### Key Output Files

#### `result_output.tsv`
Final scored dataset with all quality criteria assessments, BAGS grades, and ranking scores. Contains:
- Original specimen data
- Quality criteria scores (PASS/FAIL for each criterion)
- BAGS grades (A, B, C, D, F)
- Composite quality scores
- Ranking within taxonomic groups

#### `assessed_BAGS.tsv`
Species-level BAGS assessment with simplified 4-column format:
- `taxonid`: Taxonomic identifier
- `bags_grade`: BAGS grade (A/B/C/D/F)
- `bin_uri`: BIN identifier URL
- `sharers`: Number of species sharing the BIN

#### Family Databases
Individual SQLite databases for each family (or combined small families), containing:
- Complete specimen records for the family
- All quality assessments
- BAGS grades
- Optimized for family-specific queries and analysis

### Log Files
Comprehensive logging for each processing step, including:
- Processing statistics
- Error messages and warnings
- Performance metrics
- Progress indicators for long-running operations

## Configuration Recommendations

### For Custom Results/Logs Directories
```yaml
config:
  RESULTS_DIR: "/custom/path/results"
  LOG_DIR: "/custom/path/logs"
  # ... other configuration options
```
