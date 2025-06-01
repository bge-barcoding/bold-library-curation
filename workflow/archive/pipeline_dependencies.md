# Pipeline Dependencies and Execution Order (Even-Faster-New-BAGS Version)

## Dependency Chain Analysis

### Phase 1: Initial Setup (Conditional Pre-filtering)
```
1. prescoring_filter (if ENABLE_PRESCORING_FILTER=true)
   ├── Input: BOLD_TSV file
   ├── Conditional filters: taxa, country, marker, bins
   └── Output: PRESCORING_FILTERED_OUTPUT, results/prescoring_filter.ok
   
   OR
   
1. skip_prescoring_filter (if ENABLE_PRESCORING_FILTER=false)
   └── Output: results/prescoring_filter.ok (dummy marker)
```

### Phase 2: Database Creation
```
2. create_load_db
   ├── Input: schema, get_input_file() (either BOLD_TSV or PRESCORING_FILTERED_OUTPUT)
   ├── Dependency: results/prescoring_filter.ok
   └── Output: DB_FILE
```

### Phase 3: Database Setup
```
3. load_criteria
   ├── Input: criteria.tsv, DB_FILE
   └── Output: results/criteria_loaded.ok

4. apply_indexes
   ├── Input: indexes, DB_FILE, results/criteria_loaded.ok
   └── Output: DB_FILE_INDEXED
```

### Phase 4: Taxonomy Loading (Enhanced Performance)
```
5. load_taxonomy (using load_taxonomy_faster.pl)
   ├── Input: DB_FILE, DB_FILE_INDEXED
   ├── Enhanced: Configurable chunk size (TAXONOMY_CHUNK_SIZE)
   └── Output: results/taxonomy_loaded.ok
```

6. import_target_list (conditional - if USE_TARGET_LIST=true)
   ├── Input: DB_FILE, results/taxonomy_loaded.ok, TARGET_LIST
   └── Output: results/target_loaded.ok
   
   OR
   
6. skip_target_list (if USE_TARGET_LIST=false)
   ├── Input: results/taxonomy_loaded.ok
   └── Output: results/target_loaded.ok (pass-through)
```

### Phase 5A: BAGS Analysis Path (Heavily Optimized)
```
7A. optimize_bags_database **NEW OPTIMIZATION PHASE**
    ├── Input: DB_FILE, get_taxonomy_dependency() → results/target_loaded.ok
    ├── Optimizations: WAL mode, cache tuning, BAGS-specific indexes
    └── Output: results/bags_optimized.ok

8A. BAGS (Enhanced with Progress Tracking)
    ├── Input: DB_FILE, results/bags_optimized.ok
    ├── Enhanced: Progress tracking, performance monitoring, grade distribution
    ├── Script: assess_taxa_progress.pl (likely faster than standard assess_criteria.pl)
    └── Output: results/assessed_BAGS.tsv

9A. import_bags **NEW STEP**
    ├── Input: results/assessed_BAGS.tsv, DB_FILE
    ├── Creates: bags table in database
    └── Output: results/bags_imported.ok

10A. inherit_subspecies_bags **NEW FEATURE**
     ├── Input: DB_FILE, results/bags_imported.ok
     ├── Function: Propagates BAGS grades from species to subspecies
     └── Output: results/subspecies_bags_inherited.ok
```

### Phase 5B: Criteria Assessment Path (Parallel with BAGS)
```
7B-22B. All 16 criteria rules (running in parallel after taxonomy):
    ├── COLLECTION_DATE, COLLECTORS, COUNTRY, ID_METHOD, IDENTIFIER, INSTITUTION
    ├── COORD, MUSEUM_ID, PUBLIC_VOUCHER, SEQ_QUALITY, SITE, REGION
    ├── SECTOR, SPECIES_ID, TYPE_SPECIMEN, HAS_IMAGE (using assess_images.pl)
    
    Each rule:
    ├── Input: DB_FILE, get_taxonomy_dependency() → results/target_loaded.ok
    └── Output: results/assessed_[CRITERION].tsv
```

### Phase 6: Data Integration and Final Output
```
23. concatenate
    ├── Input: All 16 assessed_[CRITERION].tsv files
    ├── Script: concat_tsvs.pl
    └── Output: results/CONCATENATED.tsv

24. import_concatenated **UPDATED DEPENDENCY**
    ├── Input: results/CONCATENATED.tsv, DB_FILE
    ├── **NEW Dependency:** results/subspecies_bags_inherited.ok (ensures BAGS completion)
    └── Output: results/concatenated_imported.ok

25. output_filtered_data
    ├── Input: DB_FILE, results/concatenated_imported.ok
    ├── Query: ranking_with_sumscore.sql
    └── Output: results/result_output.tsv (final target)
```

## Key Improvements in Even-Faster-New-BAGS Version

### **1. Pre-processing Optimizations**
- **Conditional pre-filtering:** Optional prescoring_filter reduces dataset size early
- **Dynamic input selection:** get_input_file() chooses filtered vs. original data
- **Flexible filtering:** Support for taxa, country, marker, and BIN filtering

### **2. Database Performance Enhancements**
- **BAGS-specific optimization phase:** Dedicated database tuning before BAGS analysis
- **WAL mode:** Write-Ahead Logging for better concurrent performance
- **Chunked taxonomy loading:** Configurable TAXONOMY_CHUNK_SIZE parameter
- **Enhanced indexing:** BAGS-specific indexes via bags_indexes.sql

### **3. BAGS Analysis Improvements**
- **Progress tracking:** assess_taxa_progress.pl provides real-time feedback
- **Database integration:** BAGS results imported into database tables
- **Subspecies inheritance:** Automatic grade propagation to subspecies
- **Performance monitoring:** Detailed logging and statistics

### **4. Parallel Execution Optimization**
```
Timeline Optimization:
Phase 4 (Taxonomy) → Phase 5A (BAGS) + Phase 5B (Criteria) → Phase 6 (Output)
                           ↓                    ↓
                    BAGS optimization      16 criteria rules
                    BAGS analysis         (parallel execution)
                    BAGS import
                    Subspecies inherit
```

## Critical Path Analysis

### **Optimized Execution Flow:**
```
Time →  Phase 1-4  →        Phase 5A + 5B        →  Phase 6
       (Setup)     →    (BAGS + Criteria Parallel)  →  (Final)
```

### **Enhanced Dependency Graph:**
```
prescoring_filter.ok
    ↓
create_load_db (using dynamic input)
    ↓
load_criteria → apply_indexes
    ↓
load_taxonomy → target_loaded.ok
    ↓
    ├── optimize_bags_database → BAGS → import_bags → inherit_subspecies_bags
    └── 16 criteria rules (parallel)
            ↓
        concatenate → import_concatenated (waits for subspecies_bags_inherited.ok)
            ↓
        output_filtered_data
```

### **Performance Improvements Summary:**
1. **Reduced data volume** through optional pre-filtering
2. **Database optimization** specifically tuned for BAGS analysis  
3. **Progress visibility** during long-running BAGS analysis
4. **Integrated BAGS data** available for final ranking queries
5. **Subspecies coverage** through automatic grade inheritance
6. **Parallel execution** maintained between BAGS and criteria assessment

### **New Bottleneck Mitigation:**
- BAGS analysis now has dedicated optimization phase reducing runtime
- Progress tracking provides visibility into BAGS execution
- Database integration allows for more sophisticated final queries
- Subspecies inheritance ensures complete taxonomic coverage

## Key Dependencies for Execution Planning

**Critical path items that must complete sequentially:**
1. prescoring_filter.ok → create_load_db → load_criteria → apply_indexes → load_taxonomy → target_loaded.ok
2. target_loaded.ok → optimize_bags_database → BAGS → import_bags → inherit_subspecies_bags
3. subspecies_bags_inherited.ok + CONCATENATED.tsv → import_concatenated → output_filtered_data

**Parallel execution opportunities:**
- All 16 criteria rules can run simultaneously after target_loaded.ok
- BAGS optimization and analysis runs independently of criteria assessment
- Final integration waits for both paths to complete

## Major Changes from Previous Version

**New Rules Added:**
- `optimize_bags_database` - Database performance tuning for BAGS
- `import_bags` - Import BAGS results into database tables
- `inherit_subspecies_bags` - Propagate grades to subspecies

**Enhanced Rules:**
- `prescoring_filter` - Now conditional with flexible filtering options
- `load_taxonomy` - Uses faster script with configurable chunking
- `BAGS` - Enhanced with progress tracking and monitoring
- `import_concatenated` - Now depends on subspecies inheritance completion

**New Dependencies:**
- `results/bags_optimized.ok` - Ensures database is tuned before BAGS
- `results/bags_imported.ok` - BAGS data available in database
- `results/subspecies_bags_inherited.ok` - Complete taxonomic coverage

**Performance Optimizations:**
- Pre-filtering reduces dataset size early in pipeline
- BAGS-specific database optimizations improve analysis speed
- Progress tracking provides visibility into long-running operations
- Subspecies inheritance ensures complete coverage without re-analysis
