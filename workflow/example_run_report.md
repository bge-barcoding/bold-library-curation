# BOLD Pipeline Example Run - Output Files Report

## Pipeline Configuration
- **Input File**: `resources/test_data/test_odonata_large.tsv`
- **Results Directory**: `results/` (configurable via `RESULTS_DIR`)
- **Logs Directory**: `logs/` (configurable via `LOG_DIR`)
- **Pre-scoring Filter**: Enabled
- **Target List**: Disabled
- **OTU Clustering**: Enabled (threshold: 0.99, threads: 8)
- **Haplotype Analysis**: Enabled
- **Family Size Threshold**: 10,000 records

## Expected Output Files

### Phase 1: Data Preparation
```
results/
├── prescoring_filter.ok                    # Filter completion marker
└── prescoring_filtered.tsv                # Filtered input data (configurable filename)

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
└── target_loaded.ok                       # Target list marker (passthrough when disabled)

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
├── assessed_HAPLOTYPE_ID.tsv              # Haplotype identification assessment
└── assessed_OTU_CLUSTERING.tsv            # OTU clustering results
```
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
├── assessed_HAPLOTYPE_ID.tsv              # Haplotype identification assessment
└── assessed_OTU_CLUSTERING.tsv            # OTU clustering results

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
├── assess_HAPLOTYPE_ID.log
└── assess_OTU_CLUSTERING.log
```

### Phase 4: BAGS Assessment
```
results/
├── bags_optimized.ok                      # BAGS database optimization marker
├── assessed_BAGS.tsv                      # BAGS assessment results (simplified 4-column format)
├── bags_imported.ok                       # BAGS import completion marker
└── subspecies_bags_inherited.ok           # Subspecies inheritance marker

logs/
├── optimize_bags_database.log             # BAGS optimization log
├── assess_BAGS.log                        # BAGS assessment log with progress tracking
├── import_bags.log                        # BAGS import log
└── inherit_subspecies_bags.log            # Subspecies inheritance log
```

### Phase 5: Data Integration and Ranking
```
results/
├── CONCATENATED.tsv                       # Combined criteria assessments (excluding haplotypes/OTUs)
├── concatenated_imported.ok               # Import completion marker
├── haplotypes_imported.ok                 # Haplotype data import marker
├── otus_imported.ok                       # OTU data import marker
├── schema_with_ranks_applied.ok           # Ranking schema marker
├── ranking_indexes_applied.ok             # Ranking indexes marker
├── ranks_calculated.ok                    # Rank calculation marker
├── country_representatives_selected.ok    # Country representative selection marker
└── result_output.tsv                      # Final scored and ranked output with OTUs

logs/
├── import_concatenated.log                # Concatenation import log
├── import_haplotypes.log                  # Haplotype import log
├── import_otus.log                        # OTU import log
├── create_ranks_schema.log                # Schema creation log
├── apply_ranking_indexes.log              # Ranking indexes log
├── calculate_store_ranks.log              # Rank calculation log
├── select_country_representatives.log     # Country representative selection log
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
Final scored dataset with all quality criteria assessments, BAGS grades, haplotype IDs, OTU assignments, and comprehensive ranking scores. Contains:
- Original specimen data from BOLD
- Quality criteria scores (PASS/FAIL for each of 16 criteria)
- BAGS grades (A, B, C, D, F) at species level
- Haplotype identifiers within BIN/species groups
- OTU clustering assignments at configurable similarity threshold
- Composite quality scores and rankings
- Country representative selection indicators

#### `assessed_BAGS.tsv`
Species-level BAGS assessment with simplified 4-column format for optimal performance:
- `taxonid`: Taxonomic identifier
- `BAGS_grade`: BAGS grade (A/B/C/D/F)
- `bin_uri`: BIN identifier URL
- `sharers`: Number of species sharing the BIN

#### `assessed_HAPLOTYPE_ID.tsv`
Haplotype identification results showing genetic diversity within species/BIN groups:
- `recordid`: Specimen identifier
- `HAPLOTYPE_ID`: Unique haplotype identifier within group
- Contains analysis of sequence variation patterns

#### `assessed_OTU_CLUSTERING.tsv`
OTU clustering results using VSEARCH at configurable similarity threshold:
- `recordid`: Specimen identifier
- `OTU_ID`: Operational Taxonomic Unit identifier
- Groups specimens by genetic similarity for phylogenetic analysis

#### Family Databases
Individual SQLite databases for each family (or combined small families), containing:
- Complete specimen records for the taxonomic group
- All quality assessments and criteria scores
- BAGS grades and haplotype assignments
- OTU clustering results
- Country representative selections
- Optimized for family-specific queries and analysis
- Hierarchical organization by phylum for logical structure

### Database Schema Enhancement
The final database includes specialized tables for advanced analysis:
- `bold`: Core specimen data with taxonomic information
- `bold_criteria`: Quality criteria assessments for all 16 criteria
- `bold_haplotypes`: Haplotype assignments with genetic diversity metrics
- `bold_otus`: OTU clustering results for phylogenetic grouping
- `bold_ranks`: Comprehensive ranking scores combining all assessments
- `bags`: Species-level BAGS grades with BIN sharing information
- `country_representatives`: Selected optimal representatives per country/species/OTU
- Enhanced indexes for efficient querying across all data types

## Advanced Analysis Features

### Phylogenetic Integration
- **Haplotype Analysis**: Identifies unique genetic variants within species/BIN groups
- **OTU Clustering**: Groups specimens by configurable genetic similarity thresholds
- **Multi-threaded Processing**: Parallel execution for computationally intensive phylogenetic analyses

### Geographic Representation
- **Country Representatives**: Systematic selection of optimal specimens per geographic region
- **Balanced Sampling**: Maintains geographic diversity while optimizing specimen quality
- **Ranking Integration**: Combines quality scores with geographic considerations

### Performance Optimizations
- **BAGS Database Optimization**: Specialized indexes and query optimization for species-level analysis
- **Chunked Processing**: Memory-efficient handling of large datasets
- **Progress Tracking**: Real-time monitoring for long-running operations

## Configuration Recommendations

### For Custom Results/Logs Directories
```yaml
config:
  RESULTS_DIR: "/custom/path/results"
  LOG_DIR: "/custom/path/logs"
  # ... other configuration options
```

### For Performance Tuning
```yaml
config:
  TAXONOMY_CHUNK_SIZE: 15000              # Increase for more memory
  OTU_CLUSTERING_THREADS: 16              # Match available CPU cores
  OTU_CLUSTERING_THRESHOLD: 0.97          # Adjust clustering stringency
  FAMILY_SIZE_THRESHOLD: 5000             # Lower for more granular family DBs
```
