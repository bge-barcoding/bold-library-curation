# BOLD Pipeline Example Run - Output Files Report

## Pipeline Configuration
- **Input File**: `resources/test_data/test_odonata_large.tsv`
- **Results Directory**: `results/` (configurable via `RESULTS_DIR`)
- **Logs Directory**: `logs/` (configurable via `LOG_DIR`)
- **Pre-scoring Filter**: Enabled
- **Target List**: Disabled
- **OTU Clustering**: Enabled (threshold: 0.99, threads: 8)
- **BAGS Assessment**: Enabled with species-level grading
- **Phylogenetic Analysis**: Enabled with parallel SLURM processing
- **Family Size Threshold**: 10,000 records
- **Phylogenetic Configuration**:
  - Alignment Method: MAFFT
  - Tree Method: FastTree
  - Bootstrap: 1000 iterations
  - PDF Generation: Enabled
  - Minimum OTUs: 3

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
└── assess_OTU_CLUSTERING.log
```

### Phase 4: BAGS Assessment and Optimization
```
results/
├── bags_optimized.ok                      # BAGS database optimization marker
├── assessed_BAGS.tsv                      # BAGS assessment results (simplified 4-column format)
├── bags_imported.ok                       # BAGS import completion marker
└── subspecies_bags_inherited.ok           # Subspecies inheritance marker

logs/
├── optimize_bags_database.log             # BAGS optimization log
├── assess_BAGS.log                        # BAGS assessment log with progress tracking
├── bags_full_debug.log                    # Detailed BAGS analysis debug log
├── import_bags.log                        # BAGS import log
└── inherit_subspecies_bags.log            # Subspecies inheritance log
```

### Phase 5: Parallel Phylogenetic Analysis (When Enabled)
```
results/
├── phylo_batches/                         # Phylogenetic batch preparation
│   ├── batch_1.json                      # Family batch 1 configuration
│   ├── batch_2.json                      # Family batch 2 configuration
│   ├── batch_N.json                      # Additional batches as needed
│   └── phylo_batch_summary.json          # Batch summary statistics
├── phylogenies/                          # Phylogenetic analysis results
│   ├── Aeshnidae/                        # Family-level analysis directory
│   │   ├── Aeshnidae_alignment.fasta     # Multiple sequence alignment
│   │   ├── Aeshnidae_tree.newick         # Phylogenetic tree
│   │   ├── Aeshnidae_tree.pdf            # Tree visualization PDF
│   │   ├── Aeshnidae_tree.svg            # Tree visualization SVG
│   │   ├── Aeshnidae_curation_checklist.pdf # Grade C species curation checklist
│   │   ├── Aeshnidae_analysis_summary.json # Analysis statistics and parameters
│   │   └── Aeshnidae_failed.log          # Error log (if analysis failed)
│   ├── Libellulidae/                     # Additional family directories
│   │   └── [similar structure]
│   └── phylogenetic_analysis_summary.json # Overall phylogenetic analysis summary
├── phylogenetic_analysis_parallel_completed.ok # Phylogenetic completion marker
├── phylo_custom_parameters.json          # Custom parameters used
└── phylogenetic_results_integrated.ok    # Results integration marker

logs/
├── prepare_phylo_batches.log              # Batch preparation log
├── phylogenetic_analysis_parallel.log    # Main phylogenetic analysis log
├── phylo_checkpoints/                    # SLURM job monitoring logs
├── phylo_1_123.out                       # SLURM job array output logs
├── phylo_1_123.err                       # SLURM job array error logs
└── phylo_N_N.out/.err                    # Additional job logs
```

### Phase 6: Data Integration and Ranking
```
results/
├── CONCATENATED.tsv                       # Combined criteria assessments (excluding OTUs)
├── concatenated_imported.ok               # Import completion marker
├── otus_imported.ok                       # OTU data import marker
├── schema_with_ranks_applied.ok           # Ranking schema marker
├── ranking_indexes_applied.ok             # Ranking indexes marker
├── ranks_calculated.ok                    # Rank calculation marker
├── country_representatives_selected.ok    # Country representative selection marker
├── manual_curation_populated.ok          # Manual curation URLs populated
└── result_output.tsv                      # Final scored and ranked output with OTUs

logs/
├── import_concatenated.log                # Concatenation import log
├── import_otus.log                        # OTU import log
├── create_ranks_schema.log                # Schema creation log
├── apply_ranking_indexes.log              # Ranking indexes log
├── calculate_store_ranks.log              # Rank calculation log
├── select_country_representatives.log     # Country representative selection log
├── populate_manual_curation.log          # Manual curation URL generation log
└── output_filtered_data.log               # Final output generation log
```

### Phase 7: Family-Level Database Creation (Parallel)
```
results/
├── family_batches/                        # Family splitting batch preparation
│   ├── batch_0.json                      # Family batch 0 configuration
│   ├── batch_1.json                      # Family batch 1 configuration
│   ├── batch_N.json                      # Additional batches as needed
│   └── batch_summary.json                # Batch summary with family assignments
├── families_split.ok                      # Family splitting completion marker
├── phylogenetic_results_integrated.ok    # Phylogenetic integration completion
├── family_databases_compressed.ok        # Compression completion marker
└── family_databases/                      # Family-level database directory
    ├── splitting_report.txt               # Detailed splitting statistics
    ├── Arthropoda/                        # Phylum-level organization
    │   ├── Odonata/                       # Order-level organization
    │   │   ├── Aeshnidae/                 # Family directory
    │   │   │   ├── Aeshnidae.db          # Family database
    │   │   │   ├── Aeshnidae_tree.pdf    # Phylogenetic tree PDF
    │   │   │   ├── Aeshnidae_tree.svg    # Phylogenetic tree SVG
    │   │   │   └── Aeshnidae_curation_checklist.pdf # Curation checklist
    │   │   ├── Libellulidae/              # Another family directory
    │   │   │   └── [similar structure]
    │   │   └── Odonata_small_families/    # Combined small families
    │   │       └── [combined files]
    │   └── [other_orders]/                # Additional orders as needed
    └── [other_phyla]/                     # Additional phyla as needed

logs/
├── split_families.log                     # Family splitting process log
├── integrate_phylogenetic_results.log    # Phylogenetic integration log
├── compress_family_databases.log         # Database compression log
├── family_splitting/                     # SLURM job array logs for family splitting
│   ├── family_split_123_0.out           # Job array output logs
│   ├── family_split_123_0.err           # Job array error logs
│   └── family_split_123_N.out/.err      # Additional job logs
└── consolidate_results.log               # Result consolidation log
```

### Phase 8: Statistical Reporting and Final Archiving
```
results/
├── bold_database_statistics.pdf          # Comprehensive statistical analysis report
├── stats_report_generated.ok             # Statistics report completion marker
├── final_results_archived.ok             # Final archiving completion marker
└── final_results_20241217_143022/        # Timestamped final archive
    ├── result_output.tsv.gz              # Compressed final results
    ├── bold.db.gz                        # Compressed main database
    ├── family_databases_compressed.tar.gz # Compressed family databases
    ├── logs.tar.gz                       # Compressed complete log archive
    └── bold_database_statistics.pdf      # Statistics report

logs/
├── stats_report.log                      # Statistics report generation log
└── archive_final_results.log             # Final archiving log
```

### Compressed Family Databases
```
results/
└── family_databases_compressed/          # Compressed family directories
    ├── Arthropoda_Odonata_Aeshnidae.zip # Compressed family directory with all files
    ├── Arthropoda_Odonata_Libellulidae.zip # Another compressed family
    ├── Arthropoda_Odonata_small_families.zip # Compressed combined families
    └── [additional_compressed_families].zip # More compressed families
```

## File Descriptions

### Key Output Files

#### `result_output.tsv`
Final scored dataset with all quality criteria assessments, BAGS grades, OTU assignments, and comprehensive ranking scores. Contains:
- Original specimen data from BOLD
- Quality criteria scores (PASS/FAIL for each of 16 criteria)
- BAGS grades (A, B, C, D, F) at species level
- OTU clustering assignments at configurable similarity threshold
- Composite quality scores and rankings
- Country representative selection indicators
- Manual curation URLs for direct BOLD portal access

#### `assessed_BAGS.tsv`
Species-level BAGS assessment with simplified 4-column format for optimal performance:
- `taxonid`: Taxonomic identifier
- `BAGS_grade`: BAGS grade (A/B/C/D/F)
- `bin_uri`: BIN identifier URL
- `sharers`: Number of species sharing the BIN

#### `assessed_OTU_CLUSTERING.tsv`
OTU clustering results using VSEARCH at configurable similarity threshold:
- `recordid`: Specimen identifier
- `OTU_ID`: Operational Taxonomic Unit identifier
- Groups specimens by genetic similarity for phylogenetic analysis

#### Family Databases with Integrated Phylogenetic Results
Individual SQLite databases for each family (or combined small families), containing:
- Complete specimen records for the taxonomic group
- All quality assessments and criteria scores
- BAGS grades and OTU assignments
- Country representative selections
- **Integrated phylogenetic tree PDFs** for visualization
- **Integrated phylogenetic tree SVGs** for scalable graphics
- **Curation checklist PDFs** for Grade C species requiring attention
- Optimized for family-specific queries and analysis
- Hierarchical organization by phylum/order/family for logical structure

#### Phylogenetic Analysis Outputs (When Enabled)

##### Family-Level Phylogenetic Trees
- **Newick Format Trees**: Standard phylogenetic tree format for analysis software
- **PDF Visualizations**: High-quality tree diagrams with bootstrap support values
- **SVG Visualizations**: Scalable vector graphics for web and print applications
- **Bootstrap Analysis**: Statistical support for tree topology (default: 1000 iterations)
- **Outgroup Selection**: Enhanced strategy with hierarchical fallback (order → class)

##### Curation Checklists
- **Grade C Species Focus**: Automated identification of species requiring taxonomic attention
- **PDF Format**: Professional curation checklists for systematic review
- **Monophyly Analysis**: Assessment of taxonomic coherence within phylogenetic context
- **Priority Recommendations**: Ranked suggestions for curation activities

##### Analysis Summaries
- **JSON Format Statistics**: Machine-readable analysis metadata
- **Processing Parameters**: Complete record of analysis settings
- **Success/Failure Reports**: Detailed logs for troubleshooting
- **Performance Metrics**: Processing times and resource utilization

#### Statistical Report (`bold_database_statistics.pdf`)
Comprehensive PDF report including:
- **Dataset Overview**: Input statistics and processing summary
- **Quality Criteria Analysis**: Distribution of assessment results across all criteria
- **BAGS Grade Distribution**: Species-level quality analysis with detailed breakdowns
- **Geographic Coverage**: Country and region representation analysis
- **Taxonomic Coverage**: Phylum, order, family distribution with diversity metrics
- **OTU Clustering Statistics**: Genetic diversity patterns and clustering effectiveness
- **Phylogenetic Analysis Summary**: Success rates, family coverage, and tree quality metrics
- **Performance Benchmarks**: Processing times, resource utilization, and optimization recommendations

#### Final Archive (`final_results_[timestamp]/`)
Timestamped delivery package containing:
- **Compressed Results**: All key outputs in space-efficient formats
- **Complete Database**: Full SQLite database with all assessments
- **Family Databases**: Compressed family-level databases with phylogenetic integration
- **Documentation**: Statistics report and complete processing logs
- **Metadata**: Processing parameters and pipeline version information
- **Size Optimization**: Typically 60-80% smaller than uncompressed outputs
- **Delivery Ready**: Single directory containing all essential pipeline outputs

### Database Schema Enhancement
The final database includes specialized tables for advanced analysis:
- `bold`: Core specimen data with taxonomic information
- `bold_criteria`: Quality criteria assessments for all 16 criteria
- `bold_otus`: OTU clustering results for phylogenetic grouping
- `bold_ranks`: Comprehensive ranking scores combining all assessments
- `bags`: Species-level BAGS grades with BIN sharing information
- `country_representatives`: Selected optimal representatives per country/species/OTU
- `manual_curation`: BOLD portal URLs for direct specimen access
- Enhanced indexes for efficient querying across all data types

## Advanced Analysis Features

### Parallel Phylogenetic Analysis
- **SLURM Job Arrays**: Scalable processing across HPC clusters
- **Enhanced Outgroup Selection**: Hierarchical strategy improving success rates
- **Multiple Analysis Methods**: Support for MAFFT/MUSCLE alignment and FastTree/IQ-TREE reconstruction
- **Automated Quality Control**: Bootstrap analysis and tree validation
- **PDF Visualization**: Professional tree diagrams with ETE3
- **Curation Integration**: Direct connection between phylogenetic patterns and taxonomic issues

### Geographic Representation
- **Country Representatives**: Systematic selection of optimal specimens per geographic region
- **Balanced Sampling**: Maintains geographic diversity while optimizing specimen quality
- **Ranking Integration**: Combines quality scores with geographic considerations

### Performance Optimizations
- **BAGS Database Optimization**: Specialized indexes and query optimization for species-level analysis
- **Parallel Processing**: SLURM job arrays for computationally intensive operations
- **Memory Management**: Chunked processing for large datasets
- **Progress Tracking**: Real-time monitoring for long-running operations
- **Resource Scaling**: Configurable memory and CPU allocation

### Quality Assurance
- **Comprehensive Logging**: Detailed logs for all pipeline phases
- **Error Recovery**: Robust handling of failed analyses with detailed reporting
- **Validation Checks**: Automated verification of output completeness
- **Performance Monitoring**: Resource utilization tracking and optimization recommendations

## Configuration Examples

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

### For Phylogenetic Analysis
```yaml
config:
  PHYLO_ENABLED: true                     # Enable phylogenetic analysis
  PHYLO_NUM_JOBS: 100                     # Number of parallel jobs
  PHYLO_MIN_OTUS: 4                       # Minimum OTUs for analysis
  PHYLO_ALIGNMENT_METHOD: "mafft"         # Alignment method
  PHYLO_TREE_METHOD: "fasttree"           # Tree reconstruction method
  PHYLO_BOOTSTRAP: 1000                   # Bootstrap iterations
  PHYLO_GENERATE_PDFS: true               # Enable PDF generation
  PHYLO_MAX_CONCURRENT: 20                # Maximum concurrent jobs
  PHYLO_JOB_MEMORY: "32G"                 # Memory per job
  PHYLO_JOB_TIME: "24:00:00"              # Time limit per job
```

### For Large-Scale Processing
```yaml
config:
  # Parallel processing optimization
  FAMILY_ARRAY_SIZE: 128                  # More parallel jobs for family splitting
  WORKERS_PER_JOB: 8                      # More workers per job
  COMPRESSION_WORKERS: 32                 # Parallel compression threads
  
  # Memory optimization
  STATS_REPORT_CPU: 16                    # More cores for statistics
  STATS_REPORT_MEMORY_MB: 64000           # More memory for large datasets
  
  # Quality filtering
  ENABLE_PRESCORING_FILTER: true          # Enable pre-filtering
  FILTER_SPECIES: true                    # Focus on species-level data
  FILTER_TAXA: true                       # Filter by target taxa
```

## Expected Processing Times

### Small Dataset (< 10K records)
- **Total Runtime**: 2-4 hours
- **Phylogenetic Analysis**: +1-2 hours (if enabled)
- **Family Database Creation**: 15-30 minutes

### Medium Dataset (10K-100K records)
- **Total Runtime**: 6-12 hours
- **Phylogenetic Analysis**: +3-8 hours (if enabled)
- **Family Database Creation**: 1-3 hours

### Large Dataset (100K+ records)
- **Total Runtime**: 12-48 hours
- **Phylogenetic Analysis**: +6-24 hours (if enabled)
- **Family Database Creation**: 3-8 hours

*Processing times depend on dataset complexity, HPC cluster capacity, and configuration parameters.*