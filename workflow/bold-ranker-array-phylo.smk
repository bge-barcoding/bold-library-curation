# BOLD Ranker Workflow
# ===============================================================================
# This Snakefile processes BOLD sequence data through quality assessment 
# criteria and splits results into family-level databases for efficient analysis.

# Configuration and helper functions
# ----------------------------------
import os

configfile: "config/config.yml"

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
        base_dep = f"{results_dir}/target_loaded.ok"
    else:
        base_dep = f"{results_dir}/taxonomy_loaded.ok"
    
    # All assessment rules wait for database optimization to complete
    return [base_dep, f"{results_dir}/bags_optimized.ok"]

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
            species_arg="--filter-species" if config.get("FILTER_SPECIES", False) else "",
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
                {params.species_arg} \
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

rule OTU_CLUSTERING:
    """Identify OTUs (Operational Taxonomic Units) using VSEARCH clustering"""
    input:
        db=get_db_file(),
        taxonomy_ok=get_taxonomy_dependency()
    output:
        tsv=f"{get_results_dir()}/assessed_OTU_CLUSTERING.tsv"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        threshold=config.get("OTU_CLUSTERING_THRESHOLD", 0.99),  # Default 99% similarity
        temp_dir=f"{get_results_dir()}/temp_otu"
    log: f"{get_log_dir()}/assess_OTU_CLUSTERING.log"
    conda: "envs/otu_clustering.yaml"
    threads: config.get("OTU_CLUSTERING_THREADS", 8)
    shell:
        """
        echo "Starting OTU clustering analysis..." > {log}
        echo "Threshold: {params.threshold}" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "Temp directory: {params.temp_dir}" >> {log}
        echo "Conda environment: $CONDA_PREFIX" >> {log}
        echo "" >> {log}
        
        # Check if vsearch is available
        if command -v vsearch >/dev/null 2>&1; then
            echo "VSEARCH found at: $(which vsearch)" >> {log}
        else
            echo "ERROR: VSEARCH not found in PATH" >> {log}
            echo "Available binaries in conda environment:" >> {log}
            ls -la $CONDA_PREFIX/bin/ | grep -i vsearch >> {log} || echo "No vsearch found in conda bin directory" >> {log}
            exit 1
        fi
        echo "" >> {log}
        
        # Create temporary directory
        mkdir -p {params.temp_dir}
        
        # Run the OTU clustering script
        perl -I{params.libs} workflow/scripts/assess_otu_clustering.pl \
            --db {input.db} \
            --log {params.log_level} \
            --threshold {params.threshold} \
            --threads {threads} \
            --temp-dir {params.temp_dir} \
            2>> {log} > {output.tsv}
        
        # Clean up temporary directory
        rm -rf {params.temp_dir}
        
        echo "OTU clustering completed on $(date)" >> {log}
        """

# PHASE 4: BAGS ASSESSMENT AND OPTIMIZATION
# =========================================

rule optimize_bags_database:
    """Apply BAGS-specific database optimizations for improved performance"""
    input:
        db=get_db_file(),
        taxonomy_ok=lambda wildcards: f"{get_results_dir()}/target_loaded.ok" if config.get("USE_TARGET_LIST", False) else f"{get_results_dir()}/taxonomy_loaded.ok"
    output:
        f"{get_results_dir()}/bags_optimized.ok"
    log: f"{get_log_dir()}/optimize_bags_database.log"
    conda: "envs/sqlite.yaml"
    shell:
        """
        echo "Applying BAGS-specific database optimizations..." > {log}
        
        sqlite3 {input.db} < workflow/scripts/bags_indexes.sql 2>> {log}

        echo "BAGS database optimization completed on $(date)" >> {log}
        touch {output}
        """

# BAGS (species-level assessment)
rule BAGS:
    input:
        db=get_db_file(),
        bags_optimized=f"{get_results_dir()}/bags_optimized.ok"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        log_dir=get_log_dir()
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
            2> >(tee {params.log_dir}/bags_full_debug.log | grep "PROGRESS:" >> {log}) || \
        # Fallback for systems without process substitution
        (perl -I{params.libs} workflow/scripts/assess_taxa_simplified.pl \
            --db {input.db} \
            --progress 50 \
            > {output.tsv} 2> {params.log_dir}/bags_all_output.log && \
         echo "Extracting progress information..." >> {log} && \
         grep "PROGRESS:" {params.log_dir}/bags_all_output.log >> {log} || true)
        
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
    """Combine all individual criteria assessment results into single file"""
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

rule import_otus:
    """Import OTU data into specialized OTU table"""
    input:
        otu_tsv=rules.OTU_CLUSTERING.output.tsv,
        db=get_db_file(),
        concatenated_ok=f"{get_results_dir()}/concatenated_imported.ok"  # Add dependency to prevent concurrent database access
    output:
        f"{get_results_dir()}/otus_imported.ok"
    params:
        log_level=config['LOG_LEVEL']
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/import_otus.log"
    shell:
        """
        echo "Importing OTU data..." > {log}
        echo "Waiting for database to be available..." >> {log}
        
        # Add retry logic with timeout to handle transient database locks
        for i in {{1..5}}; do
            echo "Import attempt $i..." >> {log}
            
            if sqlite3 {input.db} 2>> {log} <<OTU
-- Clear any existing OTU data
DELETE FROM bold_otus;

-- Import the OTU data
.mode tabs
.import {input.otu_tsv} otu_temp

-- Insert data, skipping header row
INSERT INTO bold_otus (recordid, otu_id) 
SELECT recordid, OTU_ID 
FROM otu_temp 
WHERE recordid != 'recordid' 
AND recordid IS NOT NULL 
AND OTU_ID IS NOT NULL
AND OTU_ID != 'UNASSIGNED';

-- Drop temporary table
DROP TABLE otu_temp;

-- Show import statistics
SELECT 'Total OTU assignments imported: ' || COUNT(*) as result FROM bold_otus;
SELECT 'Unique OTUs: ' || COUNT(DISTINCT otu_id) as result FROM bold_otus;

-- Show top 10 largest OTUs
SELECT 'Top 10 largest OTUs:' as result;
SELECT otu_id, COUNT(*) as member_count 
FROM bold_otus 
GROUP BY otu_id 
ORDER BY member_count DESC 
LIMIT 10;

.quit
OTU
            then
                echo "OTU import completed successfully on $(date)" >> {log}
                break
            else
                echo "Import attempt $i failed, waiting before retry..." >> {log}
                sleep $((i * 2))  # Exponential backoff: 2, 4, 6, 8, 10 seconds
            fi
        done
        
        touch {output}
        """

# PHYLOGENETIC ANALYSIS - PARALLEL VERSION
# ========================================

rule prepare_phylo_batches:
    """Prepare family batches for parallel phylogenetic analysis"""
    input:
        db=get_db_file(),
        otus_ok=f"{get_results_dir()}/otus_imported.ok"
    output:
        batch_dir=directory(f"{get_results_dir()}/phylo_batches"),
        summary=f"{get_results_dir()}/phylo_batches/phylo_batch_summary.json"
    params:
        phylo_enabled=config.get("PHYLO_ENABLED", False),
        num_jobs=config.get("PHYLO_NUM_JOBS", 50),
        families_per_job=config.get("PHYLO_FAMILIES_PER_JOB", None),
        min_otus=config.get("PHYLO_MIN_OTUS", 3),
        export_kingdoms=config.get("PHYLO_KINGDOMS", ["all"]),
        kingdoms_arg=lambda wildcards: f"--export-kingdoms {' '.join(config.get('PHYLO_KINGDOMS', ['all']))}"
    conda: "envs/phylogenetic_analysis.yaml"
    log: f"{get_log_dir()}/prepare_phylo_batches.log"
    shell:
        """
        if [ "{params.phylo_enabled}" = "True" ]; then
            echo "=== Preparing Phylogenetic Analysis Batches ===" > {log}
            echo "Start time: $(date)" >> {log}
            echo "Database: {input.db}" >> {log}
            echo "Number of jobs: {params.num_jobs}" >> {log}
            echo "Minimum OTUs: {params.min_otus}" >> {log}
            echo "Export kingdoms: {params.export_kingdoms}" >> {log}
            echo "" >> {log}
            
            # Create batch directory
            mkdir -p {output.batch_dir}
            
            # Prepare batches
            python workflow/scripts/prepare_phylo_batches.py \
                {input.db} \
                --output-dir {output.batch_dir} \
                --num-jobs {params.num_jobs} \
                --min-otus {params.min_otus} \
                {params.kingdoms_arg} \
                2>> {log}
            
            echo "Batch preparation completed at $(date)" >> {log}
        else
            echo "Phylogenetic analysis disabled - creating empty batch directory" > {log}
            mkdir -p {output.batch_dir}
            echo '{{"total_batches": 0, "total_families": 0, "phylo_disabled": true}}' > {output.summary}
        fi
        """

rule run_phylogenetic_analysis_parallel:
    """Run phylogenetic analysis in parallel using SLURM job arrays"""
    input:
        batch_dir=f"{get_results_dir()}/phylo_batches",
        summary=f"{get_results_dir()}/phylo_batches/phylo_batch_summary.json",
        db=get_db_file()
    output:
        marker=f"{get_results_dir()}/phylogenetic_analysis_parallel_completed.ok",
        custom_params=f"{get_results_dir()}/phylo_custom_parameters.json"
    params:
        phylo_enabled=config.get("PHYLO_ENABLED", False),
        output_dir=f"{get_results_dir()}/{config.get('PHYLO_OUTPUT_DIR', 'phylogenies')}",
        min_otus=config.get("PHYLO_MIN_OTUS", 3),
        alignment_method=config.get("PHYLO_ALIGNMENT_METHOD", "mafft"),
        tree_method=config.get("PHYLO_TREE_METHOD", "fasttree"),
        bootstrap=config.get("PHYLO_BOOTSTRAP", 1000),
        num_outgroups=config.get("PHYLO_NUM_OUTGROUPS", 3),
        generate_pdfs=config.get("PHYLO_GENERATE_PDFS", True),
        job_memory=config.get("PHYLO_JOB_MEMORY", "32G"),
        job_time=config.get("PHYLO_JOB_TIME", "24:00:00"),
        max_concurrent=config.get("PHYLO_MAX_CONCURRENT", 10),
        partition=config.get("PHYLO_PARTITION", "week"),
        custom_parameters=config.get("PHYLO_CUSTOM_PARAMETERS", {}),
        log_dir=get_log_dir()
    log: f"{get_log_dir()}/phylogenetic_analysis_parallel.log"
    shell:
        """
        if [ "{params.phylo_enabled}" = "True" ]; then
            echo "=== Starting Parallel Phylogenetic Analysis ===" > {log}
            echo "Start time: $(date)" >> {log}
            echo "Database: {input.db}" >> {log}
            echo "Batch directory: {input.batch_dir}" >> {log}
            echo "Output directory: {params.output_dir}" >> {log}
            echo "Alignment method: {params.alignment_method}" >> {log}
            echo "Tree method: {params.tree_method}" >> {log}
            echo "Generate PDFs: {params.generate_pdfs}" >> {log}
            echo "" >> {log}
            
            # Create custom parameters file
            echo "Creating custom parameters file..." >> {log}
            python3 -c "
import json
import sys

# Custom parameters from Snakemake config
custom_params = {params.custom_parameters}

# Write to JSON file
with open('{output.custom_params}', 'w') as f:
    json.dump(custom_params, f, indent=2)

print('Custom parameters saved to {output.custom_params}')
if custom_params:
    for tool, tool_params in custom_params.items():
        if tool_params:  # Only show non-empty parameter lists
            print(f'  {{tool}}: {{tool_params}}')
else:
    print('  Using optimized defaults (no custom parameters)')
" 2>> {log}
            
            # Create output and log directories
            mkdir -p {params.output_dir}
            mkdir -p {params.log_dir}/phylo_checkpoints
            
            # Check if batches exist
            BATCH_SUMMARY=$(cat {input.summary})
            TOTAL_BATCHES=$(echo "$BATCH_SUMMARY" | python3 -c "import sys, json; data=json.load(sys.stdin); print(data.get('total_batches', 0))")
            
            if [ "$TOTAL_BATCHES" -eq 0 ]; then
                echo "No batches to process - phylogenetic analysis skipped" >> {log}
                touch {output.marker}
                exit 0
            fi
            
            echo "Processing $TOTAL_BATCHES batches..." >> {log}
            
            # Determine array range
            ARRAY_RANGE="1-$TOTAL_BATCHES%{params.max_concurrent}"
            echo "SLURM array range: $ARRAY_RANGE" >> {log}
            
            # Submit SLURM job array
            echo "Submitting SLURM job array..." >> {log}
            
            JOB_ID=$(sbatch \
                --array=$ARRAY_RANGE \
                --job-name=phylo_parallel \
                --output={params.log_dir}/phylo_%A_%a.out \
                --error={params.log_dir}/phylo_%A_%a.err \
                --time={params.job_time} \
                --mem={params.job_memory} \
                --cpus-per-task=8 \
                --partition={params.partition} \
                workflow/scripts/phylo_array_job.sh \
                {input.db} \
                {input.batch_dir} \
                {params.output_dir} \
                {params.min_otus} \
                {params.alignment_method} \
                {params.tree_method} \
                {params.bootstrap} \
                {params.num_outgroups} \
                {params.generate_pdfs} \
                {output.custom_params} \
                {params.log_dir} | awk '{{print $4}}')
            
            if [ -z "$JOB_ID" ]; then
                echo "ERROR: Failed to submit SLURM job array" >> {log}
                exit 1
            fi
            
            echo "Submitted job array with ID: $JOB_ID" >> {log}
            
            # Wait for job completion
            echo "Waiting for job completion..." >> {log}
            
            while squeue -j $JOB_ID -h >/dev/null 2>&1; do
                RUNNING=$(squeue -j $JOB_ID -h | wc -l)
                echo "$(date): $RUNNING tasks still running..." >> {log}
                sleep 60  # Check every minute
            done
            
            echo "Job array completed at $(date)" >> {log}
            
            # Wait a bit for file system to catch up
            sleep 30
            
            # Consolidate results
            echo "Consolidating results..." >> {log}
            
            python workflow/scripts/consolidate_phylo_results.py \
                {input.batch_dir} \
                {params.output_dir} \
                --output-file {params.output_dir}/phylogenetic_analysis_summary.json \
                --failed-families-file {params.output_dir}/failed_families.txt \
                2>> {log}
            
            echo "Phylogenetic analysis completed at $(date)" >> {log}
            touch {output.marker}
        else
            echo "Phylogenetic analysis disabled" > {log}
            echo '{{}}' > {output.custom_params}  # Create empty JSON file
            touch {output.marker}
        fi
        """

# RANKING ETC
# ===========

rule create_ranks_schema:
    input:
        db=get_db_file(),
        otus_ok=f"{get_results_dir()}/otus_imported.ok",
        concatenated_ok=f"{get_results_dir()}/concatenated_imported.ok",
        phylo_ok=f"{get_results_dir()}/phylogenetic_analysis_parallel_completed.ok" if config.get("PHYLO_ENABLED", False) else []
    output:
        f"{get_results_dir()}/schema_with_ranks_applied.ok"
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/create_ranks_schema.log"
    shell:
        """
        sqlite3 {input.db} 2> {log} <<SCHEMA
CREATE TABLE IF NOT EXISTS "bold_ranks" (
    "recordid" INTEGER NOT NULL,
    "ranking" INTEGER NOT NULL,
    "sumscore" INTEGER NOT NULL,
    "calculated_at" TEXT DEFAULT (datetime('now')),
    PRIMARY KEY (recordid),
    FOREIGN KEY(recordid) REFERENCES bold(recordid)
);
.quit
SCHEMA
        touch {output}
        """

rule apply_ranking_indexes:
    input:
        db=get_db_file(),
        schema_ok=f"{get_results_dir()}/schema_with_ranks_applied.ok"
    output: 
        f"{get_results_dir()}/ranking_indexes_applied.ok"
    log: f"{get_log_dir()}/apply_ranking_indexes.log"
    conda: "envs/sqlite.yaml"
    shell:
        """
        sqlite3 {input.db} < workflow/scripts/ranking_indexes.sql 2> {log}
        touch {output}
        """

rule calculate_store_ranks:
    input:
        db=get_db_file(),
        indexes_ok=f"{get_results_dir()}/ranking_indexes_applied.ok"
    output:
        f"{get_results_dir()}/ranks_calculated.ok"
    log: f"{get_log_dir()}/calculate_store_ranks.log"
    conda: "envs/sqlite.yaml"
    shell:
        """
        sqlite3 {input.db} < workflow/scripts/calculate_store_ranks.sql 2> {log}
        touch {output}
        """

rule select_country_representatives:
    """Select best representative record per species per OTU per country"""
    input:
        db=get_db_file(),
        ranks_ok=f"{get_results_dir()}/ranks_calculated.ok"
    output:
        f"{get_results_dir()}/country_representatives_selected.ok"
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/select_country_representatives.log"
    shell:
        """
        echo "=== Starting Country Representative Selection ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Database: {input.db}" >> {log}
        echo "" >> {log}
        
        echo "Selection criteria:" >> {log}
        echo "- Group by: country_iso + species + otu_id" >> {log}
        echo "- Select by: ranking ASC, sumscore DESC, recordid ASC" >> {log}
        echo "- Filter: species-level identification only" >> {log}
        echo "" >> {log}
        
        # Execute country representative selection
        sqlite3 {input.db} < workflow/scripts/select_country_representatives.sql 2>> {log}
        
        echo "" >> {log}
        echo "=== Country Representative Selection Completed ===" >> {log}
        echo "End time: $(date)" >> {log}
        
        touch {output}
        """

rule populate_manual_curation:
    """Populate manual_curation table with URLs for all BOLD records"""
    input:
        db=get_db_file(),
        rep_ok=f"{get_results_dir()}/country_representatives_selected.ok"
    output:
        f"{get_results_dir()}/manual_curation_populated.ok"
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/populate_manual_curation.log"
    shell:
        """
        echo "=== Starting Manual Curation Table Population ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Database: {input.db}" >> {log}
        echo "" >> {log}
        
        # Populate manual_curation table with URLs
        sqlite3 {input.db} < workflow/scripts/populate_manual_curation.sql 2>> {log}
        
        # Log statistics
        echo "" >> {log}
        echo "Population statistics:" >> {log}
        sqlite3 {input.db} "SELECT COUNT(*) as total_curation_records FROM manual_curation;" 2>> {log} | while read line; do echo "Total curation records: $line" >> {log}; done
        sqlite3 {input.db} "SELECT COUNT(*) as records_with_urls FROM manual_curation WHERE url LIKE '%/record/None';" 2>> {log} | while read line; do echo "Records with NULL/empty processid: $line" >> {log}; done
        
        echo "" >> {log}
        echo "=== Manual Curation Table Population Completed ===" >> {log}
        echo "End time: $(date)" >> {log}
        
        touch {output}
        """

rule output_filtered_data:
    """Generate final scored and ranked output with all assessments combined including OTUs"""
    input:
        db=get_db_file(),
        curation_ok=f"{get_results_dir()}/manual_curation_populated.ok"
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
.read workflow/scripts/ranking_with_stored_ranks_otu.sql
.quit
EOF
        """

# PHASE 6: PARALLEL FAMILY-LEVEL DATABASE CREATION
# ================================================

rule split_families:
    """Split main database into family-level databases using parallel job array processing"""
    input:
        result_tsv=f"{get_results_dir()}/result_output.tsv",
        db=get_db_file()
    output:
        marker=f"{get_results_dir()}/families_split.ok",
        report=f"{get_results_dir()}/family_databases/splitting_report.txt"
    params:
        threshold=config.get("FAMILY_SIZE_THRESHOLD", 10000),
        array_size=config.get("FAMILY_ARRAY_SIZE", 64),
        workers_per_job=config.get("WORKERS_PER_JOB", 4),
        job_memory=config.get("JOB_MEMORY", "8G"),
        job_time=config.get("JOB_TIME", "04:00:00"),
        export_kingdoms=config.get("EXPORT_KINGDOMS", ["all"]),
        kingdoms_arg=lambda wildcards: f"--export-kingdoms {' '.join(config.get('EXPORT_KINGDOMS', ['all']))}",
        output_dir=f"{get_results_dir()}/family_databases",
        batch_dir=f"{get_results_dir()}/family_batches",
        log_dir=get_log_dir(),
        script_dir="workflow/scripts"
    log: f"{get_log_dir()}/split_families.log"
    shell:
        """
        echo "=== Starting Parallel Family Database Splitting ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Database: {input.db}" >> {log}
        echo "Family size threshold: {params.threshold}" >> {log}
        echo "Job array size: {params.array_size}" >> {log}
        echo "Workers per job: {params.workers_per_job}" >> {log}
        echo "Export kingdoms: {params.export_kingdoms}" >> {log}
        echo "Output directory: {params.output_dir}" >> {log}
        echo "Batch directory: {params.batch_dir}" >> {log}
        echo "" >> {log}
        
        # Create necessary directories
        mkdir -p {params.output_dir}
        mkdir -p {params.batch_dir}
        mkdir -p {params.log_dir}
        mkdir -p {params.log_dir}/family_splitting
        
        # Step 1: Prepare family batches
        echo "Step 1: Preparing family batches..." >> {log}
        python {params.script_dir}/prepare_family_batches.py \
            {input.db} \
            --output-dir {params.batch_dir} \
            --num-jobs {params.array_size} \
            --threshold {params.threshold} \
            {params.kingdoms_arg} \
            2>> {log}
        
        if [ ! -f "{params.batch_dir}/batch_summary.json" ]; then
            echo "ERROR: Failed to create family batches" >> {log}
            exit 1
        fi
        
        # Determine actual number of batches created (exclude summary file)
        ACTUAL_BATCHES=$(find {params.batch_dir} -name "batch_[0-9]*.json" | wc -l)
        echo "Created $ACTUAL_BATCHES batch files" >> {log}
        
        if [ $ACTUAL_BATCHES -eq 0 ]; then
            echo "No batches created - no families to process" >> {log}
            touch {output.marker}
            echo "No families found to process for specified kingdoms: {params.export_kingdoms}" > {output.report}
            exit 0
        fi
        
        # Calculate array range (0-based indexing) - Fix: use actual count
        MAX_ARRAY_INDEX=$((ACTUAL_BATCHES - 1))
        echo "Job array range: 0-$MAX_ARRAY_INDEX" >> {log}
        
        # Step 2: Submit SLURM job array
        echo "Step 2: Submitting SLURM job array..." >> {log}
        
        JOB_ID=$(sbatch \
            --array=0-$MAX_ARRAY_INDEX \
            --job-name=bold_family_split \
            --output={params.log_dir}/family_splitting/family_split_%A_%a.out \
            --error={params.log_dir}/family_splitting/family_split_%A_%a.err \
            --time={params.job_time} \
            --mem={params.job_memory} \
            --cpus-per-task={params.workers_per_job} \
            {params.script_dir}/family_split_array.sh \
            {input.db} \
            {params.batch_dir} \
            {params.output_dir} \
            {params.threshold} \
            {params.workers_per_job} | awk '{{print $4}}')
        
        if [ -z "$JOB_ID" ]; then
            echo "ERROR: Failed to submit SLURM job array" >> {log}
            exit 1
        fi
        
        echo "Submitted job array with ID: $JOB_ID" >> {log}
        
        # Step 3: Wait for job completion
        echo "Step 3: Waiting for job completion..." >> {log}
        echo "Monitoring job $JOB_ID..." >> {log}
        
        # Wait for job to complete
        while squeue -j $JOB_ID -h >/dev/null 2>&1; do
            RUNNING=$(squeue -j $JOB_ID -h | wc -l)
            if [ $RUNNING -eq 0 ]; then
                echo "$(date): All tasks completed" >> {log}
                break
            fi
            echo "$(date): $RUNNING tasks still running..." >> {log}
            sleep 30
        done
        
        echo "Job array completed at $(date)" >> {log}
        
        # Step 4: Consolidate results
        echo "Step 4: Consolidating results..." >> {log}
        echo "Working directory: $(pwd)" >> {log}
        echo "Batch directory: {params.batch_dir}" >> {log}
        echo "Output directory: {params.output_dir}" >> {log}
        
        # Check what files exist
        echo "Files in batch directory:" >> {log}
        ls -la {params.batch_dir}/ >> {log} 2>&1 || echo "Batch directory not found" >> {log}
        echo "Files in output directory:" >> {log}
        ls -la {params.output_dir}/ >> {log} 2>&1 || echo "Output directory not found" >> {log}
        
        # Give jobs a moment to write their result files
        sleep 10
        
        # Check for result files specifically
        echo "Looking for result files:" >> {log}
        find {params.output_dir} -name "*result.json" >> {log} 2>&1 || echo "No result files found" >> {log}
        
        # Run consolidation with timeout
        timeout 300 python {params.script_dir}/consolidate_results.py \
            {params.batch_dir} \
            {params.output_dir} \
            --wait-timeout 60 \
            2>> {log} || {{
            echo "WARNING: Consolidation timed out or failed, proceeding anyway..." >> {log}
            # Create a basic report if consolidation fails
            echo "Consolidation failed or timed out" > {output.report}
            DB_COUNT=$(find {params.output_dir} -name "*.db" 2>/dev/null | wc -l || echo "0")
            echo "Found $DB_COUNT database files" >> {output.report}
        }}
        
        # Verify completion
        if [ -f "{output.report}" ]; then
            echo "Family splitting completed successfully" >> {log}
            echo "Report generated: {output.report}" >> {log}
            
            # Log summary statistics
            echo "" >> {log}
            echo "=== Final Summary ===" >> {log}
            DB_COUNT=$(find {params.output_dir} -name "*.db" 2>/dev/null | wc -l)
            echo "Total family databases created: $DB_COUNT" >> {log}
            echo "Exported kingdoms: {params.export_kingdoms}" >> {log}
            
            # Show sample of directory structure
            echo "" >> {log}
            echo "Sample output structure:" >> {log}
            find {params.output_dir} -type d | head -10 | while read dir; do
                echo "  $dir" >> {log}
            done
            
            touch {output.marker}
            echo "Parallel family splitting completed successfully at $(date)" >> {log}
        else
            echo "ERROR: Family splitting failed - no report generated" >> {log}
            exit 1
        fi
        """
	
rule compress_family_databases:
    """Compress entire family directories (including databases, PDFs, and other files) in parallel for storage efficiency"""
    input:
        families_split=f"{get_results_dir()}/families_split.ok",
        phylo_integrated=f"{get_results_dir()}/phylogenetic_results_integrated.ok"
    output:
        marker=f"{get_results_dir()}/family_databases_compressed.ok"
    params:
        source_dir=f"{get_results_dir()}/family_databases",
        output_dir=f"{get_results_dir()}/family_databases_compressed",
        workers=config.get("COMPRESSION_WORKERS", 16),
        compression_mode=config.get("COMPRESSION_MODE", "family_directories"),  # New parameter
        log_dir=get_log_dir()
    log: f"{get_log_dir()}/compress_family_databases.log"
    conda: "envs/compress_databases.yaml"
    threads: config.get("COMPRESSION_WORKERS", 16)
    shell:
        """
        echo "=== Starting Family Database Compression ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Compression mode: {params.compression_mode}" >> {log}
        echo "Source directory: {params.source_dir}" >> {log}
        echo "Output directory: {params.output_dir}" >> {log}
        echo "Workers: {params.workers}" >> {log}
        echo "" >> {log}
        
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Run compression script with specified mode
        python workflow/scripts/zip_databases.py \
            --source {params.source_dir} \
            --output {params.output_dir} \
            --log-dir {params.log_dir} \
            --workers {params.workers} \
            --mode {params.compression_mode} \
            --extensions .db \
            2>> {log}
        
        # Check if compression was successful based on compression mode
        COMPRESSED_COUNT=$(find {params.output_dir} -name "*.zip" | wc -l)
        
        if [ "{params.compression_mode}" = "family_directories" ]; then
            # Count family directories containing .db files
            ORIGINAL_COUNT=$(find {params.source_dir} -name "*.db" -exec dirname {{}} \; | sort -u | wc -l)
            ITEM_TYPE="family directories"
        else
            # Count individual .db files (legacy mode)
            ORIGINAL_COUNT=$(find {params.source_dir} -name "*.db" | wc -l)
            ITEM_TYPE="database files"
        fi
        
        echo "" >> {log}
        echo "Compression Summary:" >> {log}
        echo "Original $ITEM_TYPE: $ORIGINAL_COUNT" >> {log}
        echo "Compressed .zip files: $COMPRESSED_COUNT" >> {log}
        
        if [ "$COMPRESSED_COUNT" -eq "$ORIGINAL_COUNT" ]; then
            echo "SUCCESS: All $ITEM_TYPE compressed successfully" >> {log}
            
            touch {output.marker}
        else
            echo "ERROR: Compression incomplete. Expected $ORIGINAL_COUNT, got $COMPRESSED_COUNT" >> {log}
            exit 1
        fi
        
        echo "=== Compression Completed ===" >> {log}
        echo "End time: $(date)" >> {log}
        """

rule integrate_phylogenetic_results:
    """Integrate phylogenetic tree PDFs and curation checklist PDFs into family directories"""
    input:
        families_split=f"{get_results_dir()}/families_split.ok",
        phylo_done=f"{get_results_dir()}/phylogenetic_analysis_parallel_completed.ok" if config.get("PHYLO_ENABLED", False) else []
    output:
        marker=f"{get_results_dir()}/phylogenetic_results_integrated.ok"
    params:
        phylo_source_dir=f"{get_results_dir()}/{config.get('PHYLO_OUTPUT_DIR', 'phylogenies')}",
        family_db_dir=f"{get_results_dir()}/family_databases",
        phylo_enabled=config.get("PHYLO_ENABLED", False)
    log: f"{get_log_dir()}/integrate_phylogenetic_results.log"
    shell:
        """
        echo "=== Starting Phylogenetic PDF Integration ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Phylogenies source: {params.phylo_source_dir}" >> {log}
        echo "Family databases directory: {params.family_db_dir}" >> {log}
        echo "Phylogenetic analysis enabled: {params.phylo_enabled}" >> {log}
        echo "Integration mode: PDF files only (tree PDFs and curation checklists)" >> {log}
        echo "" >> {log}
        
        # Initialize counters
        TOTAL_FAMILIES=0
        PHYLO_MERGED=0
        CURATION_MERGED=0
        
        # Function to extract family name from various sources
        get_family_name() {{
            local family_dir="$1"
            local family_name=""
            
            # Try to get family name from directory name first
            family_name=$(basename "$family_dir")
            
            # If directory has a .db file, extract family name from filename
            if [ -z "$family_name" ] || [ "$family_name" = "." ]; then
                db_file=$(find "$family_dir" -maxdepth 1 -name "*.db" | head -1)
                if [ -n "$db_file" ]; then
                    family_name=$(basename "$db_file" .db)
                fi
            fi
            
            echo "$family_name"
        }}
        
        # Recursively find all directories containing .db files (family directories)
        if [ -d "{params.family_db_dir}" ]; then
            echo "Scanning for family directories containing .db files..." >> {log}
            
            # Create a temporary file to store family directories
            TEMP_FAMILIES=$(mktemp)
            find "{params.family_db_dir}" -name "*.db" -exec dirname {{}} \; | sort -u > "$TEMP_FAMILIES"
            
            # Count total families found
            TOTAL_FAMILIES=$(wc -l < "$TEMP_FAMILIES")
            echo "Found $TOTAL_FAMILIES family directories" >> {log}
            
            # Process each family directory
            while IFS= read -r family_dir; do
                if [ -d "$family_dir" ]; then
                    # Extract family name
                    family_name=$(get_family_name "$family_dir")
                    
                    if [ -n "$family_name" ] && [ "$family_name" != "." ]; then
                        echo "Processing family: $family_name (in $family_dir)" >> {log}
                        
                        # Integration successful flag
                        INTEGRATION_SUCCESS=false
                        
                        # Copy only PDF files if they exist and phylo is enabled
                        if [ "{params.phylo_enabled}" = "True" ] && [ -d "{params.phylo_source_dir}/$family_name" ]; then
                            echo "  Integrating PDF files for $family_name" >> {log}
                            
                            # Look for tree PDF files and copy them directly to family directory root
                            TREE_PDFS=$(find "{params.phylo_source_dir}/$family_name" -name "*_tree.pdf" 2>/dev/null)
                            
                            if [ -n "$TREE_PDFS" ]; then
                                echo "    Copying tree PDF files..." >> {log}
                                echo "$TREE_PDFS" | while read tree_pdf; do
                                    if [ -f "$tree_pdf" ]; then
                                        cp "$tree_pdf" "$family_dir/" 2>/dev/null
                                        echo "    ✓ Copied tree PDF: $(basename "$tree_pdf")" >> {log}
                                    fi
                                done
                                PHYLO_MERGED=$((PHYLO_MERGED + 1))
                                INTEGRATION_SUCCESS=true
                            else
                                echo "    ⚠ No tree PDF files found for $family_name" >> {log}
                            fi
                            
                            # Look for curation checklist PDFs and copy them directly to family directory root
                            CURATION_PDFS=$(find "{params.phylo_source_dir}/$family_name" -name "*curation_checklist.pdf" 2>/dev/null)
                            
                            if [ -n "$CURATION_PDFS" ]; then
                                echo "    Copying curation checklist PDF files..." >> {log}
                                echo "$CURATION_PDFS" | while read curation_pdf; do
                                    if [ -f "$curation_pdf" ]; then
                                        cp "$curation_pdf" "$family_dir/" 2>/dev/null
                                        echo "    ✓ Copied checklist PDF: $(basename "$curation_pdf")" >> {log}
                                    fi
                                done
                                CURATION_MERGED=$((CURATION_MERGED + 1))
                                INTEGRATION_SUCCESS=true
                            else
                                echo "    ⚠ No curation checklist PDF found for $family_name" >> {log}
                            fi
                        else
                            if [ "{params.phylo_enabled}" != "True" ]; then
                                echo "  Phylogenetic analysis disabled, skipping PDF integration" >> {log}
                            else
                                echo "  No phylogeny directory found for $family_name" >> {log}
                            fi
                        fi
                        
                        # Log final directory contents for this family
                        if [ "$INTEGRATION_SUCCESS" = "true" ]; then
                            echo "    Final family directory contents:" >> {log}
                            ls -la "$family_dir" | while read line; do
                                echo "      $line" >> {log}
                            done
                        fi
                    else
                        echo "  Warning: Could not determine family name for directory: $family_dir" >> {log}
                    fi
                fi
            done < "$TEMP_FAMILIES"
            
            # Clean up temporary file
            rm -f "$TEMP_FAMILIES"
        else
            echo "ERROR: Family databases directory not found: {params.family_db_dir}" >> {log}
            exit 1
        fi
        
        echo "" >> {log}
        echo "=== Integration Summary ===" >> {log}
        echo "Total families processed: $TOTAL_FAMILIES" >> {log}
        echo "Families with tree PDFs integrated: $PHYLO_MERGED" >> {log}
        echo "Families with curation checklist PDFs integrated: $CURATION_MERGED" >> {log}
        
        # Verify results
        if [ "$TOTAL_FAMILIES" -gt 0 ]; then
            echo "Phylogenetic results integration completed successfully" >> {log}
            
            # Create integration summary
            echo "Integration Summary:" >> {log}
            echo "  - Searched recursively in: {params.family_db_dir}" >> {log}
            echo "  - Source phylogenies from: {params.phylo_source_dir}" >> {log}
            echo "  - Families found: $TOTAL_FAMILIES" >> {log}
            echo "  - Tree PDFs integrated: $PHYLO_MERGED" >> {log}
            echo "  - Curation checklist PDFs integrated: $CURATION_MERGED" >> {log}
            
            if [ "$PHYLO_MERGED" -gt 0 ] || [ "$CURATION_MERGED" -gt 0 ]; then
                echo "  - PDF integration rate: $(echo "scale=1; ($PHYLO_MERGED + $CURATION_MERGED) * 50 / $TOTAL_FAMILIES" | bc -l)%" >> {log}
            fi
            
            touch {output.marker}
        else
            echo "ERROR: No families found to process" >> {log}
            exit 1
        fi
        
        echo "=== Phylogenetic PDF Integration Completed ===" >> {log}
        echo "End time: $(date)" >> {log}
        """

# PHASE 7: RUN SUMMARY
# ======================
rule generate_stats_report:
    """Generate comprehensive database statistics report as final pipeline step"""
    input:
        db=get_db_file(),
        families_split=f"{get_results_dir()}/families_split.ok",
        compressed=f"{get_results_dir()}/family_databases_compressed.ok"
    output:
        pdf_report=f"{get_results_dir()}/bold_database_statistics.pdf",
        marker=f"{get_results_dir()}/stats_report_generated.ok"
    params:
        db_path=get_db_file(),
        log_path=f"{get_log_dir()}/prescoring_filter.log",
        skip_detailed_taxonomy="--skip-detailed-taxonomy" if config.get("SKIP_DETAILED_TAXONOMY", True) else "",
        fast_mode="--fast-mode" if config.get("STATS_REPORT_FAST_MODE", False) else "",
        sample_size=config.get("STATS_REPORT_SAMPLE_SIZE", 100000),
        max_workers=config.get("STATS_REPORT_CPU", 8)
    resources:
        mem_mb=config.get("STATS_REPORT_MEMORY_MB", 32000),
        runtime=config.get("STATS_REPORT_RUNTIME", 120)  # 2 hours
    threads: config.get("STATS_REPORT_CPU", 8)
    conda: "envs/stats_report.yaml"
    log: f"{get_log_dir()}/stats_report.log"
    shell:
        """
        echo "=== Starting BOLD Database Statistics Report Generation ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Database: {params.db_path}" >> {log}
        echo "Memory allocated: {resources.mem_mb} MB" >> {log}
        echo "CPUs: {threads}" >> {log}
        echo "Skip detailed taxonomy: {params.skip_detailed_taxonomy}" >> {log}
        echo "Fast mode: {params.fast_mode}" >> {log}
        echo "" >> {log}
        
        # Set memory-efficient Python options for large datasets
        export PYTHONHASHSEED=0
        export OMP_NUM_THREADS={threads}
        export SQLITE_TMPDIR=${{TMPDIR:-/tmp}}
        
        # Run the statistics report generation
        python workflow/scripts/bold_stats_report_hpc.py \
            --database {params.db_path} \
            --log {params.log_path} \
            --output {output.pdf_report} \
            {params.skip_detailed_taxonomy} \
            {params.fast_mode} \
            --sample-size {params.sample_size} \
            --max-workers {params.max_workers} \
            2>> {log}
        
        if [ $? -eq 0 ]; then
            echo "Statistics report generated successfully: {output.pdf_report}" >> {log}
            echo "Report size: $(du -h {output.pdf_report} | cut -f1)" >> {log}
            touch {output.marker}
        else
            echo "ERROR: Failed to generate statistics report" >> {log}
            exit 1
        fi
        
        echo "=== Statistics Report Generation Completed ===" >> {log}
        echo "End time: $(date)" >> {log}
        """

# PHASE 8: ARCHIVE FINAL RESULTS
# ==============================
rule archive_final_results:
    """Archive key pipeline outputs with timestamp for final delivery"""
    input:
        stats_report=f"{get_results_dir()}/bold_database_statistics.pdf",
        stats_marker=f"{get_results_dir()}/stats_report_generated.ok",
        result_tsv=f"{get_results_dir()}/result_output.tsv",
        database=get_db_file(),
        compressed_marker=f"{get_results_dir()}/family_databases_compressed.ok"
    output:
        archive_marker=f"{get_results_dir()}/final_results_archived.ok"
    params:
        results_dir=get_results_dir(),
        log_dir=get_log_dir()
    log: f"{get_log_dir()}/archive_final_results.log"
    shell:
        """
        echo "=== Starting Final Results Archiving ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Results directory: {params.results_dir}" >> {log}
        echo "Log directory: {params.log_dir}" >> {log}
        echo "" >> {log}
        
        # Generate timestamp for archive folder
        TIMESTAMP=$(date +%Y%m%d_%H%M%S)
        ARCHIVE_DIR="{params.results_dir}/final_results_$TIMESTAMP"
        
        echo "Creating archive directory: $ARCHIVE_DIR" >> {log}
        mkdir -p "$ARCHIVE_DIR"
        
        # Initialize success/failure counters
        SUCCESS_COUNT=0
        FAILURE_COUNT=0
        
        echo "" >> {log}
        echo "=== Archiving Files ===" >> {log}
        
        # Archive result_output.tsv
        echo "Archiving result_output.tsv..." >> {log}
        if [ -f "{input.result_tsv}" ]; then
            if gzip -c "{input.result_tsv}" > "$ARCHIVE_DIR/result_output.tsv.gz"; then
                echo "SUCCESS: result_output.tsv compressed to result_output.tsv.gz" >> {log}
                echo "  Original size: $(du -h {input.result_tsv} | cut -f1)" >> {log}
                echo "  Compressed size: $(du -h $ARCHIVE_DIR/result_output.tsv.gz | cut -f1)" >> {log}
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            else
                echo "FAILURE: Failed to compress result_output.tsv" >> {log}
                FAILURE_COUNT=$((FAILURE_COUNT + 1))
            fi
        else
            echo "FAILURE: result_output.tsv not found at {input.result_tsv}" >> {log}
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
        fi
        
        # Archive bold.db
        echo "Archiving bold.db..." >> {log}
        if [ -f "{input.database}" ]; then
            if gzip -c "{input.database}" > "$ARCHIVE_DIR/bold.db.gz"; then
                echo "SUCCESS: bold.db compressed to bold.db.gz" >> {log}
                echo "  Original size: $(du -h {input.database} | cut -f1)" >> {log}
                echo "  Compressed size: $(du -h $ARCHIVE_DIR/bold.db.gz | cut -f1)" >> {log}
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            else
                echo "FAILURE: Failed to compress bold.db" >> {log}
                FAILURE_COUNT=$((FAILURE_COUNT + 1))
            fi
        else
            echo "FAILURE: bold.db not found at {input.database}" >> {log}
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
        fi
        
        # Archive family_databases_compressed folder
        echo "Archiving family_databases_compressed directory..." >> {log}
        FAMILY_DB_DIR="{params.results_dir}/family_databases_compressed"
        if [ -d "$FAMILY_DB_DIR" ]; then
            if tar -czf "$ARCHIVE_DIR/family_databases_compressed.tar.gz" -C "{params.results_dir}" "family_databases_compressed"; then
                echo "SUCCESS: family_databases_compressed directory archived to family_databases_compressed.tar.gz" >> {log}
                echo "  Directory size: $(du -sh $FAMILY_DB_DIR | cut -f1)" >> {log}
                echo "  Archive size: $(du -h $ARCHIVE_DIR/family_databases_compressed.tar.gz | cut -f1)" >> {log}
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            else
                echo "FAILURE: Failed to archive family_databases_compressed directory" >> {log}
                FAILURE_COUNT=$((FAILURE_COUNT + 1))
            fi
        else
            echo "FAILURE: family_databases_compressed directory not found at $FAMILY_DB_DIR" >> {log}
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
        fi
        
        # Archive logs directory
        echo "Archiving logs directory..." >> {log}
        if [ -d "{params.log_dir}" ]; then
            # Use relative path from parent directory to avoid including full path in archive
            PARENT_DIR=$(dirname "{params.log_dir}")
            LOG_DIR_NAME=$(basename "{params.log_dir}")
            if tar -czf "$ARCHIVE_DIR/logs.tar.gz" -C "$PARENT_DIR" "$LOG_DIR_NAME"; then
                echo "SUCCESS: logs directory archived to logs.tar.gz" >> {log}
                echo "  Directory size: $(du -sh {params.log_dir} | cut -f1)" >> {log}
                echo "  Archive size: $(du -h $ARCHIVE_DIR/logs.tar.gz | cut -f1)" >> {log}
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            else
                echo "FAILURE: Failed to archive logs directory" >> {log}
                FAILURE_COUNT=$((FAILURE_COUNT + 1))
            fi
        else
            echo "FAILURE: logs directory not found at {params.log_dir}" >> {log}
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
        fi
        
        # Copy statistics PDF (no compression needed)
        echo "Copying statistics report PDF..." >> {log}
        if [ -f "{input.stats_report}" ]; then
            if cp "{input.stats_report}" "$ARCHIVE_DIR/bold_database_statistics.pdf"; then
                echo "SUCCESS: bold_database_statistics.pdf copied" >> {log}
                echo "  File size: $(du -h {input.stats_report} | cut -f1)" >> {log}
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            else
                echo "FAILURE: Failed to copy bold_database_statistics.pdf" >> {log}
                FAILURE_COUNT=$((FAILURE_COUNT + 1))
            fi
        else
            echo "FAILURE: bold_database_statistics.pdf not found at {input.stats_report}" >> {log}
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
        fi
        
        echo "" >> {log}
        echo "=== Archive Summary ===" >> {log}
        echo "Archive directory: $ARCHIVE_DIR" >> {log}
        echo "Successful operations: $SUCCESS_COUNT" >> {log}
        echo "Failed operations: $FAILURE_COUNT" >> {log}
        echo "Total archive size: $(du -sh $ARCHIVE_DIR | cut -f1)" >> {log}
        
        # List final archive contents
        echo "" >> {log}
        echo "Archive contents:" >> {log}
        ls -lh "$ARCHIVE_DIR/" >> {log}
        
        echo "" >> {log}
        echo "=== Final Results Archiving Completed ===" >> {log}
        echo "End time: $(date)" >> {log}
        
        # Always create the marker file, even if some operations failed
        echo "Archive created at: $ARCHIVE_DIR" > {output.archive_marker}
        echo "Successful operations: $SUCCESS_COUNT" >> {output.archive_marker}
        echo "Failed operations: $FAILURE_COUNT" >> {output.archive_marker}
        echo "Completed at: $(date)" >> {output.archive_marker}
        """

# FINAL TARGET RULES
# ==================

rule all:
    """Main pipeline target - produces final scored output and family databases with phylogenies"""
    input:
        f"{get_results_dir()}/result_output.tsv",
        f"{get_results_dir()}/families_split.ok",
        f"{get_results_dir()}/phylogenetic_results_integrated.ok",
        f"{get_results_dir()}/family_databases_compressed.ok",
        f"{get_results_dir()}/country_representatives_selected.ok",
        f"{get_results_dir()}/stats_report_generated.ok",
        f"{get_results_dir()}/final_results_archived.ok"
    default_target: True
