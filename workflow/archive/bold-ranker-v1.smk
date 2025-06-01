# BOLD Ranker Workflow
# ===============================================================================
# This Snakefile processes BOLD sequence data through quality assessment 
# criteria and splits results into family-level databases for efficient analysis.

# Configuration and helper functions
# ----------------------------------
import os

configfile: "config/config_optimized.yml"

def get_results_dir():
    """Return the configured results directory with default fallback"""
    return config.get("RESULTS_DIR", "results")

def get_db_file():
    """Return the database file path using the configured results directory"""
    results_dir = get_results_dir()
    return f"{results_dir}/bold.db"

def get_db_file_indexed():
    """Return the indexed database marker file path using the configured results directory"""
    results_dir = get_results_dir()
    return f"{results_dir}/bold_indexed.ok"

def get_prescoring_filtered_output():
    """Return the prescoring filtered output path using the configured results directory"""
    if config.get("ENABLE_PRESCORING_FILTER", False):
        results_dir = get_results_dir()
        # Config now contains just the filename, construct full path in results directory
        filename = config.get("PRESCORING_FILTERED_OUTPUT", "prescoring_filtered.tsv")
        return f"{results_dir}/{filename}"
    return f"{get_results_dir()}/prescoring_filter_disabled.tsv"  # Return dummy path when disabled

def get_input_file():
    """Return the appropriate input file based on whether filtering is enabled"""
    if config.get("ENABLE_PRESCORING_FILTER", False):
        return get_prescoring_filtered_output()
    else:
        return config["BOLD_TSV"]

def get_log_dir():
    """Return the configured log directory"""
    return config.get("LOG_DIR", "logs")

def get_taxonomy_dependency():
    """Return appropriate dependency based on target list usage"""
    results_dir = get_results_dir()
    if config.get("USE_TARGET_LIST", False):
        return f"{results_dir}/target_loaded.ok"
    else:
        return f"{results_dir}/taxonomy_loaded.ok"

# Ensure configured directories exist
# ----------------------------------
os.makedirs(get_results_dir(), exist_ok=True)
os.makedirs(get_log_dir(), exist_ok=True)
os.makedirs(f"{get_results_dir()}/family_databases", exist_ok=True)

# PHASE 1: DATA PREPARATION AND FILTERING
# =======================================

# Optional pre-scoring filter to reduce dataset size before main processing
if config.get("ENABLE_PRESCORING_FILTER", False):
    rule prescoring_filter:
        """Filter BOLD data by taxa, countries, markers, and BIN sharing before processing"""
        input:
            bold_tsv=config["BOLD_TSV"]
        output:
            filtered_file=get_prescoring_filtered_output(),
            marker=f"{get_results_dir()}/prescoring_filter.ok"
        params:
            taxa_arg=lambda wildcards: f"--taxa-list {config['FILTER_TAXA_LIST']}" if config.get("FILTER_TAXA", False) and config.get("FILTER_TAXA_LIST") else "",
            country_arg=lambda wildcards: f"--country-list {config['FILTER_COUNTRY_LIST']}" if config.get("FILTER_COUNTRIES", False) and config.get("FILTER_COUNTRY_LIST") else "",
            marker_arg=lambda wildcards: f"--marker {config['MARKER']}" if config.get("MARKER") else "",
            bin_arg="--enable-bin-sharing" if config.get("FILTER_BINS", False) else "",
            log_level=config['LOG_LEVEL']
        log: f"{get_log_dir()}/prescoring_filter.log"
        conda: "envs/prescoring_filter.yaml"
        shell:
            """
            python workflow/scripts/prescoring_filter.py \
                --input {input.bold_tsv} \
                --output {output.filtered_file} \
                --log-level {params.log_level} \
                {params.taxa_arg} \
                {params.country_arg} \
                {params.marker_arg} \
                {params.bin_arg} \
                2> {log}
            echo "Prescoring filter completed" > {output.marker}
            """
else:
    rule skip_prescoring_filter:
        """Create marker when prescoring filter is disabled"""
        input:
            original_file=config["BOLD_TSV"]
        output:
            marker=f"{get_results_dir()}/prescoring_filter.ok"
        shell:
            """
            echo 'Prescoring filter disabled - using original file: {input.original_file}' > {output.marker}
            """

# PHASE 2: DATABASE CREATION AND SETUP
# ====================================

rule create_load_db:
    """Create SQLite database and load BOLD data using optimized fast_simple loader"""
    input:
        schema=config["SCHEMA"],
        tsv_file=get_input_file(),
        prescoring_ok=f"{get_results_dir()}/prescoring_filter.ok"
    output: 
        get_db_file()
    params: 
        log_level=config['LOG_LEVEL']
    log: f"{get_log_dir()}/create_load_db.log"
    conda: "envs/create_load_db.yaml"
    shell:
        """
        perl workflow/scripts/load_bcdm_fast_simple.pl \
            --tsv {input.tsv_file} \
            --db {output} \
            --sql {input.schema} \
            --log {params.log_level} \
            --force 2> {log}
        """

rule load_criteria:
    """Load assessment criteria definitions into database"""
    input:
        criteria="resources/criteria.tsv",
        db=get_db_file()
    output:
        f"{get_results_dir()}/criteria_loaded.ok"
    log: f"{get_log_dir()}/load_criteria.log"
    conda: "envs/sqlite.yaml"
    shell:
        """
        sqlite3 {input.db} <<CRITERIA
.mode tabs
.import {input.criteria} criteria
.quit
CRITERIA
2> {log} && touch {output}
        """

rule apply_indexes:
    """Apply database indexes to optimize query performance"""
    input:
        indexes=config["INDEXES"],
        db=get_db_file(),
        criteria_ok=f"{get_results_dir()}/criteria_loaded.ok"
    output: 
        get_db_file_indexed()
    log: f"{get_log_dir()}/apply_indexes.log"
    conda: "envs/sqlite.yaml"
    shell:
        """
        sqlite3 {input.db} < {input.indexes} 2> {log} && touch {output}
        """

rule load_taxonomy:
    """Load NCBI taxonomy data into database using optimized chunked loading"""
    input:
        db=get_db_file(),
        index_ok=get_db_file_indexed()
    output:
        f"{get_results_dir()}/taxonomy_loaded.ok"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        chunk_size=config.get("TAXONOMY_CHUNK_SIZE", 10000)
    log: f"{get_log_dir()}/load_taxonomy.log"
    conda: "envs/load_taxonomy.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/load_taxonomy_faster.pl \
            --db {input.db} \
            --log {params.log_level} \
            --chunk {params.chunk_size} \
            2> {log} && touch {output}
        """

# Optional target list processing for focused analysis
if config.get("USE_TARGET_LIST", False):
    rule import_target_list:
        """Import specific target taxa list for focused assessment"""
        input:
            db=get_db_file(),
            taxonomy_ok=f"{get_results_dir()}/taxonomy_loaded.ok",
            targetlist=config["TARGET_LIST"]
        output:
            f"{get_results_dir()}/target_loaded.ok"
        params:
            log_level=config['LOG_LEVEL'],
            libs=config["LIBS"],
            project=config["PROJECT_NAME"],
            taxon=config["TAXON_LEVEL"],
            kingdom=config["KINGDOM"]
        log: f"{get_log_dir()}/load_target_list.log"
        conda: "envs/load_taxonomy.yaml"
        shell:
            """
            perl -I{params.libs} workflow/scripts/load_targetlist.pl \
                --list {input.targetlist} \
                --db {input.db} \
                --log {params.log_level} \
                --project {params.project} \
                --taxon {params.taxon} \
                --kingdom {params.kingdom} \
                2> {log} && touch {output}
            """
else:
    rule skip_target_list:
        """Pass-through rule when target list is not used"""
        input:
            f"{get_results_dir()}/taxonomy_loaded.ok"
        output:
            f"{get_results_dir()}/target_loaded.ok"
        shell:
            "cp {input} {output}"

# PHASE 3: QUALITY CRITERIA ASSESSMENT
# ====================================
# Each rule assesses specimens against specific data quality criteria

rule COLLECTION_DATE:
    """Assess specimen collection date completeness and validity"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="COLLECTION_DATE"
    output:
        tsv=f"{get_results_dir()}/assessed_COLLECTION_DATE.tsv"
    log: f"{get_log_dir()}/assess_COLLECTION_DATE.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule COLLECTORS:
    """Assess collector information completeness"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="COLLECTORS"
    output:
        tsv=f"{get_results_dir()}/assessed_COLLECTORS.tsv"
    log: f"{get_log_dir()}/assess_COLLECTORS.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule COUNTRY:
    """Assess country information completeness and standardization"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="COUNTRY"
    output:
        tsv=f"{get_results_dir()}/assessed_COUNTRY.tsv"
    log: f"{get_log_dir()}/assess_COUNTRY.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule ID_METHOD:
    """Assess taxonomic identification method documentation"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="ID_METHOD"
    output:
        tsv=f"{get_results_dir()}/assessed_ID_METHOD.tsv"
    log: f"{get_log_dir()}/assess_ID_METHOD.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule IDENTIFIER:
    """Assess taxonomic identifier information completeness"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="IDENTIFIER"
    output:
        tsv=f"{get_results_dir()}/assessed_IDENTIFIER.tsv"
    log: f"{get_log_dir()}/assess_IDENTIFIER.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule INSTITUTION:
    """Assess institutional affiliation completeness"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="INSTITUTION"
    output:
        tsv=f"{get_results_dir()}/assessed_INSTITUTION.tsv"
    log: f"{get_log_dir()}/assess_INSTITUTION.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule COORD:
    """Assess geographic coordinate completeness and precision"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="COORD"
    output:
        tsv=f"{get_results_dir()}/assessed_COORD.tsv"
    log: f"{get_log_dir()}/assess_COORD.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule MUSEUM_ID:
    """Assess museum/institution specimen ID completeness"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="MUSEUM_ID"
    output:
        tsv=f"{get_results_dir()}/assessed_MUSEUM_ID.tsv"
    log: f"{get_log_dir()}/assess_MUSEUM_ID.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule PUBLIC_VOUCHER:
    """Assess public voucher specimen availability"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="PUBLIC_VOUCHER"
    output:
        tsv=f"{get_results_dir()}/assessed_PUBLIC_VOUCHER.tsv"
    log: f"{get_log_dir()}/assess_PUBLIC_VOUCHER.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule SEQ_QUALITY:
    """Assess DNA sequence quality metrics"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="SEQ_QUALITY"
    output:
        tsv=f"{get_results_dir()}/assessed_SEQ_QUALITY.tsv"
    log: f"{get_log_dir()}/assess_SEQ_QUALITY.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule SITE:
    """Assess collection site information completeness"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="SITE"
    output:
        tsv=f"{get_results_dir()}/assessed_SITE.tsv"
    log: f"{get_log_dir()}/assess_SITE.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule REGION:
    """Assess geographic region information completeness"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="REGION"
    output:
        tsv=f"{get_results_dir()}/assessed_REGION.tsv"
    log: f"{get_log_dir()}/assess_REGION.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule SECTOR:
    """Assess geographic sector information completeness"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="SECTOR"
    output:
        tsv=f"{get_results_dir()}/assessed_SECTOR.tsv"
    log: f"{get_log_dir()}/assess_SECTOR.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule SPECIES_ID:
    """Assess species identification completeness and accuracy"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="SPECIES_ID"
    output:
        tsv=f"{get_results_dir()}/assessed_SPECIES_ID.tsv"
    log: f"{get_log_dir()}/assess_SPECIES_ID.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule TYPE_SPECIMEN:
    """Assess type specimen designation and documentation"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="TYPE_SPECIMEN"
    output:
        tsv=f"{get_results_dir()}/assessed_TYPE_SPECIMEN.tsv"
    log: f"{get_log_dir()}/assess_TYPE_SPECIMEN.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_criteria.pl \
            --db {input.db} \
            --log {params.log_level} \
            --criteria {params.criterion} \
            2> {log} > {output.tsv}
        """

rule HAS_IMAGE:
    """Assess specimen image availability"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"]
    output:
        tsv=f"{get_results_dir()}/assessed_HAS_IMAGE.tsv"
    log: f"{get_log_dir()}/assess_HAS_IMAGE.log"
    conda: "envs/assess_images.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_images.pl \
            --db {input.db} \
            --log {params.log_level} \
            2> {log} > {output.tsv}
        """

rule HAPLOTYPE_ID:
    """Identify unique haplotypes within each BIN and species group"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    output:
        tsv=f"{get_results_dir()}/assessed_HAPLOTYPE_ID.tsv"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"]
    log: f"{get_log_dir()}/assess_HAPLOTYPE_ID.log"
    conda: "envs/haplotype_analysis.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_haplotypes.pl \
            --db {input.db} \
            --log {params.log_level} \
            2> {log} > {output.tsv}
        """

# PHASE 4: BAGS ASSESSMENT AND OPTIMIZATION
# =========================================

rule optimize_bags_database:
    """Apply BAGS-specific database optimizations for improved performance"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    output:
        f"{get_results_dir()}/bags_optimized.ok"
    log: f"{get_log_dir()}/optimize_bags_database.log"
    conda: "envs/sqlite.yaml"
    shell:
        """
        echo "Applying BAGS-specific database optimizations..." > {log}
        
        sqlite3 {input.db} 2>> {log} <<OPTIMIZE
-- BAGS Performance Optimizations (Phase 1)
PRAGMA journal_mode = WAL;
PRAGMA synchronous = NORMAL;
PRAGMA cache_size = -64000;
PRAGMA temp_store = MEMORY;

-- Apply BAGS-specific indexes
.read workflow/scripts/bags_indexes.sql

-- Update query planner statistics
ANALYZE bold;
ANALYZE taxa;

-- Final optimization
PRAGMA optimize;
.quit
OPTIMIZE

        echo "BAGS database optimization completed on $(date)" >> {log}
        touch {output}
        """

# BAGS (species-level assessment) - SIMPLIFIED OUTPUT WITH ESSENTIAL COLUMNS ONLY
rule BAGS:
    input:
        db=get_db_file(),
        bags_optimized=f"{get_results_dir()}/bags_optimized.ok"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"]
    output:
        tsv=f"{get_results_dir()}/assessed_BAGS.tsv"
    log: f"{get_log_dir()}/assess_BAGS.log"
    conda: "envs/assess_criteria.yaml"
    shell:
        """
        echo "=== BAGS Analysis Started ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Database: {input.db}" >> {log}
        echo "Using simplified output format with 4 columns: taxonid, BAGS, BIN, sharers" >> {log}
        echo "" >> {log}
        
        # Run BAGS analysis with simplified output (no taxonomic rank columns)
        echo "Starting simplified BAGS analysis..." >> {log}
        echo "This eliminates column alignment issues by removing rank columns" >> {log}
        echo "" >> {log}
        
        # Run analysis and filter for progress lines
        perl -I{params.libs} workflow/scripts/assess_taxa_simplified.pl \
            --db {input.db} \
            --progress 50 \
            > {output.tsv} \
            2> >(tee logs/bags_full_debug.log | grep "PROGRESS:" >> {log}) || \
        # Fallback for systems without process substitution
        (perl -I{params.libs} workflow/scripts/assess_taxa_simplified.pl \
            --db {input.db} \
            --progress 50 \
            > {output.tsv} 2> logs/bags_all_output.log && \
         echo "Extracting progress information..." >> {log} && \
         grep "PROGRESS:" logs/bags_all_output.log >> {log} || true)
        
        echo "" >> {log}
        echo "=== BAGS Analysis Completed ===" >> {log}
        echo "End time: $(date)" >> {log}
        
        # Generate summary statistics
        if [ -f {output.tsv} ]; then
            TOTAL_RECORDS=$(tail -n +2 {output.tsv} | wc -l)
            UNIQUE_TAXA=$(tail -n +2 {output.tsv} | cut -f1 | sort -u | wc -l)
            echo "Summary:" >> {log}
            echo "  Total records: $TOTAL_RECORDS" >> {log}
            echo "  Unique taxa: $UNIQUE_TAXA" >> {log}
            echo "" >> {log}
            echo "Grade distribution:" >> {log}
            tail -n +2 {output.tsv} | cut -f2 | sort | uniq -c | sort -rn | while read count grade; do
                echo "  Grade $grade: $count records" >> {log}
            done
        fi
        
        echo "" >> {log}
        echo "Output file: {output.tsv}" >> {log}
        echo "Format: taxonid, BAGS_grade, BIN_URL, sharers" >> {log}
        """

# Modified rule to import simplified BAGS data into database
rule import_bags:
    input:
        bags_tsv=f"{get_results_dir()}/assessed_BAGS.tsv",
        db=get_db_file()
    output:
        f"{get_results_dir()}/bags_imported.ok"
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/import_bags.log"
    shell:
        """
        echo "Importing simplified BAGS data with 4 columns..." > {log}
        echo "Expected format: taxonid, BAGS_grade, BIN_URL, sharers" >> {log}
        echo "" >> {log}
        
        sqlite3 {input.db} 2>> {log} <<BAGS
-- Clear any existing data
DELETE FROM bags;

-- Import the simplified data
.mode tabs
.import {input.bags_tsv} bags_temp

-- Insert data, skipping header row
INSERT INTO bags SELECT * FROM bags_temp WHERE taxonid != 'taxonid';

-- Drop temporary table
DROP TABLE bags_temp;

-- Show import statistics
SELECT 'Total BAGS records imported: ' || COUNT(*) as result FROM bags;
SELECT 'Unique taxa with BAGS: ' || COUNT(DISTINCT taxonid) as result FROM bags;
SELECT 'Grade distribution:' as result;
SELECT bags_grade, COUNT(*) as count FROM bags GROUP BY bags_grade ORDER BY bags_grade;

.quit
BAGS

        echo "BAGS import completed on $(date)" >> {log}
        touch {output}
        """

rule inherit_subspecies_bags:
    """Inherit BAGS grades for subspecies from their parent species"""
    input:
        db=get_db_file(),
        bags_ok=f"{get_results_dir()}/bags_imported.ok"
    output:
        f"{get_results_dir()}/subspecies_bags_inherited.ok"
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/inherit_subspecies_bags.log"
    shell:
        """
        echo "Inheriting BAGS grades for subspecies from parent species..." > {log}
        
        sqlite3 {input.db} 2>> {log} <<INHERIT
-- Insert subspecies records with inherited BAGS grades from parent species
INSERT INTO bags (taxonid, bags_grade, bin_uri, sharers)
SELECT 
    s.taxonid,
    b.bags_grade,
    b.bin_uri,
    b.sharers
FROM taxa s
JOIN taxa parent ON s.parent_taxonid = parent.taxonid
JOIN bags b ON parent.taxonid = b.taxonid
WHERE s.level = 'subspecies'
  AND parent.level = 'species'
  AND s.taxonid NOT IN (SELECT taxonid FROM bags WHERE taxonid = s.taxonid);

-- Log subspecies inheritance statistics
SELECT 'Subspecies BAGS inheritance completed: ' || COUNT(*) || ' subspecies inherited grades' as result
FROM taxa s
JOIN taxa parent ON s.parent_taxonid = parent.taxonid  
JOIN bags b ON parent.taxonid = b.taxonid
WHERE s.level = 'subspecies' AND parent.level = 'species';
.quit
INHERIT

        echo "Subspecies BAGS inheritance completed on $(date)" >> {log}
        touch {output}
        """

# PHASE 5: DATA INTEGRATION AND OUTPUT
# ===================================

rule concatenate:
    """Combine all individual criteria assessment results into single file (excluding haplotypes)"""
    input:
        collection_date = rules.COLLECTION_DATE.output.tsv,
        collectors = rules.COLLECTORS.output.tsv,
        coord = rules.COORD.output.tsv,
        country = rules.COUNTRY.output.tsv,
        id_method = rules.ID_METHOD.output.tsv,
        identifier = rules.IDENTIFIER.output.tsv,
        institution = rules.INSTITUTION.output.tsv,
        museum_id = rules.MUSEUM_ID.output.tsv,
        public_voucher = rules.PUBLIC_VOUCHER.output.tsv,
        seq_quality = rules.SEQ_QUALITY.output.tsv,
        site = rules.SITE.output.tsv,
        region = rules.REGION.output.tsv,
        sector = rules.SECTOR.output.tsv,
        species_id = rules.SPECIES_ID.output.tsv,
        type_specimen = rules.TYPE_SPECIMEN.output.tsv,
        has_image = rules.HAS_IMAGE.output.tsv
    output:
        concat=f"{get_results_dir()}/CONCATENATED.tsv"
    shell:
        """
        perl workflow/scripts/concat_tsvs.pl \
            {input.collection_date} \
            {input.collectors} \
            {input.coord} \
            {input.country} \
            {input.id_method} \
            {input.identifier} \
            {input.institution} \
            {input.museum_id} \
            {input.public_voucher} \
            {input.seq_quality} \
            {input.site} \
            {input.region} \
            {input.sector} \
            {input.species_id} \
            {input.type_specimen} \
            {input.has_image} \
            > {output.concat}
        """

rule import_concatenated:
    """Import combined criteria assessment results into database"""
    input:
        concat=f"{get_results_dir()}/CONCATENATED.tsv",
        db=get_db_file(),
        subspecies_ok=f"{get_results_dir()}/subspecies_bags_inherited.ok"
    output:
        f"{get_results_dir()}/concatenated_imported.ok"
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/import_concatenated.log"
    shell:
        """
sqlite3 {input.db} 2> {log} <<IMPORT
.mode tabs
.import {input.concat} bold_criteria
.quit
IMPORT
touch {output}
        """

rule import_haplotypes:
    """Import haplotype data into specialized haplotype tables"""
    input:
        haplotype_tsv=rules.HAPLOTYPE_ID.output.tsv,
        db=get_db_file(),
        concatenated_ok=f"{get_results_dir()}/concatenated_imported.ok"
    output:
        f"{get_results_dir()}/haplotypes_imported.ok"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"]
    conda: "envs/assess_criteria.yaml"
    log: f"{get_log_dir()}/import_haplotypes.log"
    shell:
        """
        perl -I{params.libs} workflow/scripts/load_haplotypes.pl \
            --db {input.db} \
            --tsv {input.haplotype_tsv} \
            --log {params.log_level} \
            2> {log} && touch {output}
        """

rule output_filtered_data:
    """Generate final scored and ranked output with all assessments combined"""
    input:
        db=get_db_file(),
        import_ok=f"{get_results_dir()}/concatenated_imported.ok",
        haplotypes_ok=f"{get_results_dir()}/haplotypes_imported.ok"
    output:
        f"{get_results_dir()}/result_output.tsv"
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/output_filtered_data.log"
    shell:
        """
        sqlite3 {input.db} 2> {log} <<EOF
.headers ON        
.mode tabs
.output {output}
.read workflow/scripts/ranking_with_sumscore.sql
.quit
EOF
        """

# PHASE 6: FAMILY-LEVEL DATABASE CREATION
# =======================================

rule split_families:
    """Split main database into family-level databases for efficient analysis"""
    input:
        result_tsv=f"{get_results_dir()}/result_output.tsv",
        db=get_db_file()
    output:
        marker=f"{get_results_dir()}/families_split.ok",
        report=f"{get_results_dir()}/family_databases/splitting_report.txt"
    params:
        threshold=config.get("FAMILY_SIZE_THRESHOLD", 10000),
        output_dir=f"{get_results_dir()}/family_databases",
        script_path="workflow/scripts/bold_family_splitter.py"
    log: f"{get_log_dir()}/split_families.log"
    shell:
        """
        echo "Starting family database splitting..." > {log}
        echo "Input TSV: {input.result_tsv}" >> {log}
        echo "Input DB: {input.db}" >> {log}
        echo "Threshold: {params.threshold}" >> {log}
        echo "Output directory: {params.output_dir}" >> {log}
        echo "" >> {log}
        
        # Run the Python family splitter (try TSV first, fallback to DB)
        if [ -f "{input.result_tsv}" ] && [ -s "{input.result_tsv}" ]; then
            echo "Using TSV input: {input.result_tsv}" >> {log}
            python {params.script_path} "{input.result_tsv}" --output "{params.output_dir}" --threshold {params.threshold} 2>> {log}
        elif [ -f "{input.db}" ]; then
            echo "Using database input: {input.db}" >> {log}
            python {params.script_path} "{input.db}" --output "{params.output_dir}" --threshold {params.threshold} 2>> {log}
        else
            echo "ERROR: No valid input found" >> {log}
            exit 1
        fi
        
        # Verify output
        if [ -f "{output.report}" ]; then
            echo "Family splitting completed successfully" >> {log}
            echo "Report generated: {output.report}" >> {log}
            
            # Log summary statistics
            echo "" >> {log}
            echo "=== Summary ===" >> {log}
            DB_COUNT=$(find {params.output_dir} -name "*.db" | wc -l)
            echo "Total family databases created: $DB_COUNT" >> {log}
            
            # Show directory structure
            echo "" >> {log}
            echo "Output structure:" >> {log}
            find {params.output_dir} -type d | head -20 | while read dir; do
                echo "  $dir" >> {log}
            done
            
            touch {output.marker}
        else
            echo "ERROR: Family splitting failed - no report generated" >> {log}
            exit 1
        fi
        """

rule create_final_summary:
    """Generate comprehensive pipeline execution summary"""
    input:
        result_tsv=f"{get_results_dir()}/result_output.tsv",
        families_split=f"{get_results_dir()}/families_split.ok",
        split_report=f"{get_results_dir()}/family_databases/splitting_report.txt"
    output:
        summary=f"{get_results_dir()}/pipeline_summary.txt"
    params:
        db_file=get_db_file(),
        results_dir=get_results_dir()
    shell:
        """
        echo "BOLD Library Curation Pipeline - Final Summary" > {output.summary}
        echo "=============================================" >> {output.summary}
        echo "Completed: $(date)" >> {output.summary}
        echo "" >> {output.summary}
        
        # Main results
        echo "Main Results:" >> {output.summary}
        echo "- Processed database: {params.db_file}" >> {output.summary}
        echo "- Final scored data: {input.result_tsv}" >> {output.summary}
        if [ -f "{input.result_tsv}" ]; then
            TOTAL_RECORDS=$(tail -n +2 {input.result_tsv} | wc -l)
            echo "- Total records: $TOTAL_RECORDS" >> {output.summary}
        fi
        echo "" >> {output.summary}
        
        # Family databases
        echo "Family Databases:" >> {output.summary}
        if [ -f "{input.split_report}" ]; then
            grep -A 20 "Family Statistics:" {input.split_report} | head -20 >> {output.summary}
        fi
        echo "" >> {output.summary}
        
        # File structure
        echo "Output Structure:" >> {output.summary}
        echo "{params.results_dir}/" >> {output.summary}
        echo "├── bold.db                    # Main database" >> {output.summary}
        echo "├── result_output.tsv          # Final scored data" >> {output.summary}
        echo "├── assessed_*.tsv             # Individual criteria assessments" >> {output.summary}
        echo "└── family_databases/          # Split by taxonomy" >> {output.summary}
        
        DB_COUNT=$(find {params.results_dir}/family_databases -name "*.db" 2>/dev/null | wc -l)
        echo "    ├── $DB_COUNT family/subfamily databases" >> {output.summary}
        
        PHYLA_COUNT=$(find {params.results_dir}/family_databases -maxdepth 1 -type d 2>/dev/null | wc -l)
        echo "    └── organized in $PHYLA_COUNT phyla" >> {output.summary}
        
        echo "" >> {output.summary}
        echo "For detailed family splitting report, see: {input.split_report}" >> {output.summary}
        """

# FINAL TARGET RULES
# ==================

rule all:
    """Main pipeline target - produces final scored output and family databases"""
    input:
        f"{get_results_dir()}/result_output.tsv",
        f"{get_results_dir()}/families_split.ok",
        f"{get_results_dir()}/pipeline_summary.txt"
    default_target: True

# UTILITY RULES
# =============

rule clean:
    """Remove intermediate output files to free disk space"""
    shell:
        """
        perl workflow/scripts/clean.pl
        """