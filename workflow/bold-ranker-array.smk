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
    """Compress all family database files in parallel for storage efficiency"""
    input:
        families_split=f"{get_results_dir()}/families_split.ok"
    output:
        marker=f"{get_results_dir()}/family_databases_compressed.ok",
        compression_report=f"{get_results_dir()}/family_databases/compression_report.txt"
    params:
        source_dir=f"{get_results_dir()}/family_databases",
        output_dir=f"{get_results_dir()}/family_databases_compressed",
        workers=config.get("COMPRESSION_WORKERS", 16),
        log_dir=get_log_dir()
    log: f"{get_log_dir()}/compress_family_databases.log"
    conda: "envs/compress_databases.yaml"
    threads: config.get("COMPRESSION_WORKERS", 16)
    shell:
        """
        echo "=== Starting Family Database Compression ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Source directory: {params.source_dir}" >> {log}
        echo "Output directory: {params.output_dir}" >> {log}
        echo "Workers: {params.workers}" >> {log}
        echo "" >> {log}
        
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Run compression script
        python workflow/scripts/zip_databases.py \
            --source {params.source_dir} \
            --output {params.output_dir} \
            --workers {params.workers} \
            --extensions .db \
            2>> {log}
        
        # Check if compression was successful
        COMPRESSED_COUNT=$(find {params.output_dir} -name "*.zip" | wc -l)
        ORIGINAL_COUNT=$(find {params.source_dir} -name "*.db" | wc -l)
        
        echo "" >> {log}
        echo "Compression Summary:" >> {log}
        echo "Original .db files: $ORIGINAL_COUNT" >> {log}
        echo "Compressed .zip files: $COMPRESSED_COUNT" >> {log}
        
        if [ "$COMPRESSED_COUNT" -eq "$ORIGINAL_COUNT" ]; then
            echo "SUCCESS: All database files compressed successfully" >> {log}
            
            # Generate compression report
            echo "Family Database Compression Report" > {output.compression_report}
            echo "=================================" >> {output.compression_report}
            echo "Completed: $(date)" >> {output.compression_report}
            echo "" >> {output.compression_report}
            echo "Files processed: $ORIGINAL_COUNT" >> {output.compression_report}
            echo "Files compressed: $COMPRESSED_COUNT" >> {output.compression_report}
            echo "" >> {output.compression_report}
            
            # Calculate size savings
            ORIGINAL_SIZE=$(du -sb {params.source_dir} | cut -f1)
            COMPRESSED_SIZE=$(du -sb {params.output_dir} | cut -f1)
            if [ "$ORIGINAL_SIZE" -gt 0 ]; then
                SAVINGS=$(echo "scale=1; (1 - $COMPRESSED_SIZE / $ORIGINAL_SIZE) * 100" | bc -l)
                echo "Original size: $(numfmt --to=iec $ORIGINAL_SIZE)" >> {output.compression_report}
                echo "Compressed size: $(numfmt --to=iec $COMPRESSED_SIZE)" >> {output.compression_report}
                echo "Space savings: $SAVINGS%" >> {output.compression_report}
            fi
            
            touch {output.marker}
        else
            echo "ERROR: Compression incomplete. Expected $ORIGINAL_COUNT, got $COMPRESSED_COUNT" >> {log}
            exit 1
        fi
        
        echo "=== Compression Completed ===" >> {log}
        echo "End time: $(date)" >> {log}
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
        f"{get_results_dir()}/pipeline_summary.txt",
        f"{get_results_dir()}/family_databases_compressed.ok",
        f"{get_results_dir()}/country_representatives_selected.ok"
    default_target: True
