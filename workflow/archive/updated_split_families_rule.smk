# Updated split_families rule with kingdom filtering support
# Replace the existing split_families rule in bold-ranker.smk with this version

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
        kingdoms_arg=lambda wildcards: " ".join([f"--export-kingdoms {' '.join(config.get('EXPORT_KINGDOMS', ['all']))}"]),
        output_dir=f"{get_results_dir()}/family_databases",
        batch_dir=f"{get_results_dir()}/family_batches",
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
        mkdir -p logs
        
        # Step 1: Prepare family batches with kingdom filtering
        echo "Step 1: Preparing family batches with kingdom filtering..." >> {log}
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
        
        # Determine actual number of batches created
        ACTUAL_BATCHES=$(find {params.batch_dir} -name "batch_*.json" | wc -l)
        echo "Created $ACTUAL_BATCHES batch files" >> {log}
        
        if [ $ACTUAL_BATCHES -eq 0 ]; then
            echo "No batches created - no families to process" >> {log}
            touch {output.marker}
            echo "No families found to process for specified kingdoms" > {output.report}
            exit 0
        fi
        
        # Calculate array range (0-based indexing)
        MAX_ARRAY_INDEX=$((ACTUAL_BATCHES - 1))
        echo "Job array range: 0-$MAX_ARRAY_INDEX" >> {log}
        
        # Step 2: Submit SLURM job array
        echo "Step 2: Submitting SLURM job array..." >> {log}
        
        JOB_ID=$(sbatch \
            --array=0-$MAX_ARRAY_INDEX \
            --job-name=bold_family_split \
            --output=logs/family_split_%A_%a.out \
            --error=logs/family_split_%A_%a.err \
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
            sleep 30
            RUNNING=$(squeue -j $JOB_ID -h | wc -l)
            echo "$(date): $RUNNING tasks still running..." >> {log}
        done
        
        echo "Job array completed at $(date)" >> {log}
        
        # Step 4: Consolidate results
        echo "Step 4: Consolidating results..." >> {log}
        python {params.script_dir}/consolidate_results.py \
            {params.batch_dir} \
            {params.output_dir} \
            --wait-timeout 300 \
            2>> {log}
        
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
