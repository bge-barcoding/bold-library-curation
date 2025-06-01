# BOLD Ranker Workflow
# ===============================================================================
# This Snakefile processes BOLD sequence data through quality assessment 
# criteria and splits results into family-level databases for efficient analysis.

# Configuration and helper functions
# ----------------------------------
configfile: "config/config_optimized.yml"

def get_input_file():
    """Return the appropriate input file based on whether filtering is enabled"""
    if config.get("ENABLE_PRESCORING_FILTER", False):
        return config["PRESCORING_FILTERED_OUTPUT"]
    else:
        return config["BOLD_TSV"]

def get_taxonomy_dependency():
    """Return appropriate dependency based on target list usage"""
    if config.get("USE_TARGET_LIST", False):
        return "results/target_loaded.ok"
    else:
        return "results/taxonomy_loaded.ok"

# PHASE 1: DATA PREPARATION AND FILTERING
# =======================================

# Optional pre-scoring filter to reduce dataset size before main processing
if config.get("ENABLE_PRESCORING_FILTER", False):
    rule prescoring_filter:
        """Filter BOLD data by taxa, countries, markers, and BIN sharing before processing"""
        input:
            bold_tsv=config["BOLD_TSV"]
        output:
            filtered_file=config["PRESCORING_FILTERED_OUTPUT"],
            marker="results/prescoring_filter.ok"
        params:
            taxa_arg=lambda wildcards: f"--taxa-list {config['FILTER_TAXA_LIST']}" if config.get("FILTER_TAXA", False) and config.get("FILTER_TAXA_LIST") else "",
            country_arg=lambda wildcards: f"--country-list {config['FILTER_COUNTRY_LIST']}" if config.get("FILTER_COUNTRIES", False) and config.get("FILTER_COUNTRY_LIST") else "",
            marker_arg=lambda wildcards: f"--marker {config['MARKER']}" if config.get("MARKER") else "",
            bin_arg="--enable-bin-sharing" if config.get("FILTER_BINS", False) else "",
            log_level=config['LOG_LEVEL']
        log: "logs/prescoring_filter.log"
        conda: "envs/prescoring_filter.yaml"
        shell:
            """
            mkdir -p $(dirname {output.filtered_file})
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
        output:
            marker="results/prescoring_filter.ok"
        shell:
            "echo 'Prescoring filter disabled - using original file' > {output.marker}"

# PHASE 2: DATABASE CREATION AND SETUP
# ====================================

rule create_load_db:
    """Create SQLite database and load BOLD data using optimized fast_simple loader"""
    input:
        schema=config["SCHEMA"],
        tsv_file=get_input_file(),
        prescoring_ok="results/prescoring_filter.ok"
    output: 
        config["DB_FILE"]
    params: 
        log_level=config['LOG_LEVEL']
    log: "logs/create_load_db.log"
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
        db=config["DB_FILE"]
    output:
        "results/criteria_loaded.ok"
    log: "logs/load_criteria.log"
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
        db=config["DB_FILE"],
        criteria_ok="results/criteria_loaded.ok"
    output: 
        config["DB_FILE_INDEXED"]
    log: "logs/apply_indexes.log"
    conda: "envs/sqlite.yaml"
    shell:
        """
        sqlite3 {input.db} < {input.indexes} 2> {log} && touch {output}
        """

rule load_taxonomy:
    """Load NCBI taxonomy data into database using optimized chunked loading"""
    input:
        db=config["DB_FILE"],
        index_ok=config["DB_FILE_INDEXED"]
    output:
        "results/taxonomy_loaded.ok"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        chunk_size=config.get("TAXONOMY_CHUNK_SIZE", 10000)
    log: "logs/load_taxonomy.log"
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
            db=config["DB_FILE"],
            taxonomy_ok="results/taxonomy_loaded.ok",
            targetlist=config["TARGET_LIST"]
        output:
            "results/target_loaded.ok"
        params:
            log_level=config['LOG_LEVEL'],
            libs=config["LIBS"],
            project=config["PROJECT_NAME"],
            taxon=config["TAXON_LEVEL"],
            kingdom=config["KINGDOM"]
        log: "logs/load_target_list.log"
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
            "results/taxonomy_loaded.ok"
        output:
            "results/target_loaded.ok"
        shell:
            "cp {input} {output}"

# PHASE 3: QUALITY CRITERIA ASSESSMENT
# ====================================
# Each rule assesses specimens against specific data quality criteria

rule COLLECTION_DATE:
    """Assess specimen collection date completeness and validity"""
    input:
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="COLLECTION_DATE"
    output:
        tsv="results/assessed_COLLECTION_DATE.tsv"
    log: "logs/assess_COLLECTION_DATE.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="COLLECTORS"
    output:
        tsv="results/assessed_COLLECTORS.tsv"
    log: "logs/assess_COLLECTORS.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="COUNTRY"
    output:
        tsv="results/assessed_COUNTRY.tsv"
    log: "logs/assess_COUNTRY.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="ID_METHOD"
    output:
        tsv="results/assessed_ID_METHOD.tsv"
    log: "logs/assess_ID_METHOD.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="IDENTIFIER"
    output:
        tsv="results/assessed_IDENTIFIER.tsv"
    log: "logs/assess_IDENTIFIER.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="INSTITUTION"
    output:
        tsv="results/assessed_INSTITUTION.tsv"
    log: "logs/assess_INSTITUTION.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="COORD"
    output:
        tsv="results/assessed_COORD.tsv"
    log: "logs/assess_COORD.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="MUSEUM_ID"
    output:
        tsv="results/assessed_MUSEUM_ID.tsv"
    log: "logs/assess_MUSEUM_ID.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="PUBLIC_VOUCHER"
    output:
        tsv="results/assessed_PUBLIC_VOUCHER.tsv"
    log: "logs/assess_PUBLIC_VOUCHER.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="SEQ_QUALITY"
    output:
        tsv="results/assessed_SEQ_QUALITY.tsv"
    log: "logs/assess_SEQ_QUALITY.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="SITE"
    output:
        tsv="results/assessed_SITE.tsv"
    log: "logs/assess_SITE.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="REGION"
    output:
        tsv="results/assessed_REGION.tsv"
    log: "logs/assess_REGION.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="SECTOR"
    output:
        tsv="results/assessed_SECTOR.tsv"
    log: "logs/assess_SECTOR.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="SPECIES_ID"
    output:
        tsv="results/assessed_SPECIES_ID.tsv"
    log: "logs/assess_SPECIES_ID.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"],
        criterion="TYPE_SPECIMEN"
    output:
        tsv="results/assessed_TYPE_SPECIMEN.tsv"
    log: "logs/assess_TYPE_SPECIMEN.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"]
    output:
        tsv="results/assessed_HAS_IMAGE.tsv"
    log: "logs/assess_HAS_IMAGE.log"
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
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    output:
        tsv="results/assessed_HAPLOTYPE_ID.tsv",
        haplotypes_ok="results/haplotypes_assigned.ok"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"]
    log: "logs/assess_HAPLOTYPE_ID.log"
    conda: "envs/haplotype_analysis.yaml"
    shell:
        """
        perl -I{params.libs} workflow/scripts/assess_haplotypes.pl \
            --db {input.db} \
            --log {params.log_level} \
            2> {log} > {output.tsv}
        
        echo "Haplotype analysis completed on $(date)" >> {log}
        touch {output.haplotypes_ok}
        """

# PHASE 4: BAGS ASSESSMENT AND OPTIMIZATION
# =========================================

rule optimize_bags_database:
    """Apply BAGS-specific database optimizations for improved performance"""
    input:
        db=config["DB_FILE"],
        taxonomy_ok=get_taxonomy_dependency()
    output:
        "results/bags_optimized.ok"
    log: "logs/optimize_bags_database.log"
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
        db=config["DB_FILE"],
        bags_optimized="results/bags_optimized.ok"
    params:
        log_level=config['LOG_LEVEL'],
        libs=config["LIBS"]
    output:
        tsv="results/assessed_BAGS.tsv"
    log: "logs/assess_BAGS.log"
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
        bags_tsv="results/assessed_BAGS.tsv",
        db=config["DB_FILE"]
    output:
        "results/bags_imported.ok"
    conda: "envs/sqlite.yaml"
    log: "logs/import_bags.log"
    shell:
        """
        echo "Importing simplified BAGS data with 4 columns..." > {log}
        echo "Expected format: taxonid, BAGS_grade, BIN_URL, sharers" >> {log}
        echo "" >> {log}
        
        sqlite3 {input.db} 2>> {log} <<BAGS
CREATE TABLE IF NOT EXISTS bags (
    taxonid INTEGER,
    bags_grade TEXT,
    bin_uri TEXT,
    sharers TEXT
);

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
        db=config["DB_FILE"],
        bags_ok="results/bags_imported.ok"
    output:
        "results/subspecies_bags_inherited.ok"
    conda: "envs/sqlite.yaml"
    log: "logs/inherit_subspecies_bags.log"
    shell:
        """
        echo "Inheriting BAGS grades for subspecies from parent species..." > {log}
        
        sqlite3 {input.db} 2>> {log} <<INHERIT
-- Insert subspecies records with inherited BAGS grades from parent species
INSERT INTO bags (taxonid, order_name, family_name, genus_name, species_name, bags_grade, bin_uri, sharers)
SELECT 
    s.taxonid,
    b.order_name,
    b.family_name, 
    b.genus_name,
    s.name as species_name,
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
        has_image = rules.HAS_IMAGE.output.tsv,
        haplotype_id = rules.HAPLOTYPE_ID.output.tsv
    output:
        concat="results/CONCATENATED.tsv"
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
            {input.haplotype_id} \
            > {output.concat}
        """

rule import_concatenated:
    """Import combined criteria assessment results into database"""
    input:
        concat="results/CONCATENATED.tsv",
        db=config["DB_FILE"],
        subspecies_ok="results/subspecies_bags_inherited.ok",
        haplotypes_ok="results/haplotypes_assigned.ok"
    output:
        "results/concatenated_imported.ok"
    conda: "envs/sqlite.yaml"
    log: "logs/import_concatenated.log"
    shell:
        """
sqlite3 {input.db} 2> {log} <<IMPORT
.mode tabs
.import {input.concat} bold_criteria
.quit
IMPORT
touch {output}
        """

rule output_filtered_data:
    """Generate final scored and ranked output with all assessments combined"""
    input:
        db=config["DB_FILE"],
        import_ok="results/concatenated_imported.ok"
    output:
        "results/result_output.tsv"
    conda: "envs/sqlite.yaml"
    log: "logs/output_filtered_data.log"
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
        result_tsv="results/result_output.tsv",
        db=config["DB_FILE"]
    output:
        marker="results/families_split.ok",
        report="results/family_databases/splitting_report.txt"
    params:
        threshold=config.get("FAMILY_SIZE_THRESHOLD", 10000),
        output_dir="results/family_databases",
        script_path="workflow/scripts/bold_family_splitter.js"
    log: "logs/split_families.log"
    conda: "envs/family_splitter.yaml"
    shell:
        """
        echo "Starting family database splitting..." > {log}
        echo "Input TSV: {input.result_tsv}" >> {log}
        echo "Input DB: {input.db}" >> {log}
        echo "Threshold: {params.threshold}" >> {log}
        echo "Output directory: {params.output_dir}" >> {log}
        echo "" >> {log}
        
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Run the family splitter (try TSV first, fallback to DB)
        if [ -f "{input.result_tsv}" ] && [ -s "{input.result_tsv}" ]; then
            echo "Using TSV input: {input.result_tsv}" >> {log}
            node {params.script_path} "{input.result_tsv}" "{params.output_dir}" {params.threshold} 2>> {log}
        elif [ -f "{input.db}" ]; then
            echo "Using database input: {input.db}" >> {log}
            node {params.script_path} "{input.db}" "{params.output_dir}" {params.threshold} 2>> {log}
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
        result_tsv="results/result_output.tsv",
        families_split="results/families_split.ok",
        split_report="results/family_databases/splitting_report.txt"
    output:
        summary="results/pipeline_summary.txt"
    shell:
        """
        echo "BOLD Library Curation Pipeline - Final Summary" > {output.summary}
        echo "=============================================" >> {output.summary}
        echo "Completed: $(date)" >> {output.summary}
        echo "" >> {output.summary}
        
        # Main results
        echo "Main Results:" >> {output.summary}
        echo "- Processed database: results/bold.db" >> {output.summary}
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
        echo "results/" >> {output.summary}
        echo "├── bold.db                    # Main database" >> {output.summary}
        echo "├── result_output.tsv          # Final scored data" >> {output.summary}
        echo "├── assessed_*.tsv             # Individual criteria assessments" >> {output.summary}
        echo "└── family_databases/          # Split by taxonomy" >> {output.summary}
        
        DB_COUNT=$(find results/family_databases -name "*.db" 2>/dev/null | wc -l)
        echo "    ├── $DB_COUNT family/subfamily databases" >> {output.summary}
        
        PHYLA_COUNT=$(find results/family_databases -maxdepth 1 -type d 2>/dev/null | wc -l)
        echo "    └── organized in $PHYLA_COUNT phyla" >> {output.summary}
        
        echo "" >> {output.summary}
        echo "For detailed family splitting report, see: {input.split_report}" >> {output.summary}
        """

# FINAL TARGET RULES
# ==================

rule all:
    """Main pipeline target - produces final scored output and family databases"""
    input:
        "results/result_output.tsv",
        "results/families_split.ok",
        "results/pipeline_summary.txt"
    default_target: True

# UTILITY RULES
# =============

rule clean:
    """Remove intermediate output files to free disk space"""
    shell:
        """
        perl workflow/scripts/clean.pl
        """