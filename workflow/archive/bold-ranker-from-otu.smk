# BOLD Ranker Workflow - Starting from import_otus
# ===============================================================================
# Minimal version that starts from import_otus and completes the pipeline

import os

configfile: "config/config.yml"

def get_results_dir():
    """Return the configured results directory with default fallback"""
    return config.get("RESULTS_DIR", "results")

def get_db_file():
    """Return the database file path using the configured results directory"""
    results_dir = get_results_dir()
    return f"{results_dir}/bold.db"

def get_log_dir():
    """Return the configured log directory"""
    return config.get("LOG_DIR", "logs")

# Ensure configured directories exist
os.makedirs(get_results_dir(), exist_ok=True)
os.makedirs(get_log_dir(), exist_ok=True)
os.makedirs(f"{get_results_dir()}/family_databases", exist_ok=True)

# START FROM HERE: Import OTUs
rule import_otus:
    """Import OTU data into specialized OTU table"""
    input:
        otu_tsv=f"{get_results_dir()}/assessed_OTU_CLUSTERING.tsv",
        db=get_db_file(),
        concatenated_ok=f"{get_results_dir()}/concatenated_imported.ok"
    output:
        f"{get_results_dir()}/otus_imported.ok"
    params:
        log_level=config['LOG_LEVEL']
    conda: "envs/sqlite.yaml"
    log: f"{get_log_dir()}/import_otus.log"
    shell:
        """
        echo "Importing OTU data..." > {log}
        
        sqlite3 {input.db} 2>> {log} <<OTU
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

        echo "OTU import completed on $(date)" >> {log}
        touch {output}
        """

rule create_ranks_schema:
    input:
        db=get_db_file(),
        otus_ok=f"{get_results_dir()}/otus_imported.ok"
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
        script_path="workflow/scripts/bold_family_splitter.py",
        results_dir=get_results_dir()
    log: f"{get_log_dir()}/split_families.log"
    shell:
        """
        echo "Starting family database splitting..." > {log}
        echo "Results directory: {params.results_dir}" >> {log}
        echo "Threshold: {params.threshold}" >> {log}
        echo "Output directory: {params.output_dir}" >> {log}
        echo "" >> {log}
        
        # Run the Python family splitter with directory path
        # The script will automatically find the best input file in this order:
        # 1. bold.db (highest priority)
        # 2. Any other .db file
        # 3. result_output.tsv  
        # 4. Other fallbacks
        echo "Script will auto-resolve input file with priority: bold.db → other .db → result_output.tsv" >> {log}
        python {params.script_path} "{params.results_dir}" --output "{params.output_dir}" --threshold {params.threshold} 2>> {log}
        
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

# Default rule - runs everything from import_otus to the end
rule all:
    """Main pipeline target - produces final scored output and family databases"""
    input:
        f"{get_results_dir()}/result_output.tsv",
        f"{get_results_dir()}/families_split.ok",
        f"{get_results_dir()}/pipeline_summary.txt",
        f"{get_results_dir()}/country_representatives_selected.ok"
    default_target: True