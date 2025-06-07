#!/bin/bash
#SBATCH --job-name=haplotype-analysis
#SBATCH --partition=day
#SBATCH --output=job_haplotype_%j.out
#SBATCH --error=job_haplotype_%j.err
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=b.price@nhm.ac.uk
#SBATCH --mail-type=ALL

# Optimized SLURM script for independent haplotype analysis
# Designed for maximum performance on 3M+ records
# Adjust resources based on your HPC cluster capabilities

echo "=================================================="
echo "BOLD Haplotype Analysis - Optimized Version"
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: $SLURM_MEM_PER_NODE MB"
echo ""

# Load required modules (adjust for your HPC environment)
# module load perl/5.30.0
# module load sqlite/3.36.0

# Activate conda environment
source activate bold-curation

# Set environment variables for optimal performance
export PERL_HASH_SEED=0  # For reproducible results
export TMPDIR=$SLURM_TMPDIR  # Use local fast storage for temp files
export SQLITE_TMPDIR=$SLURM_TMPDIR

# Performance parameters (auto-configured based on allocated resources)
PROCESSES=$SLURM_CPUS_PER_TASK
CHUNK_SIZE=10000  # Large chunks for high-memory environment
CACHE_SIZE=500000  # Large cache for 120GB memory

# Input/output paths (MODIFY THESE FOR YOUR SETUP)
DB_FILE="results/bold.db"
OUTPUT_FILE="results/assessed_HAPLOTYPE_ID_optimized.tsv"
LOG_FILE="logs/haplotype_analysis_optimized_${SLURM_JOB_ID}.log"

# Create directories if they don't exist
mkdir -p logs
mkdir -p results

echo "Configuration:"
echo "  Database: $DB_FILE"
echo "  Output: $OUTPUT_FILE"
echo "  Processes: $PROCESSES"
echo "  Chunk size: $CHUNK_SIZE"
echo "  Cache size: $CACHE_SIZE"
echo ""

# Check if database exists
if [[ ! -f "$DB_FILE" ]]; then
    echo "ERROR: Database file not found: $DB_FILE"
    echo "Please ensure the database has been created by the main pipeline first."
    exit 1
fi

# Check database size and record count
echo "Database Information:"
DB_SIZE=$(du -h "$DB_FILE" | cut -f1)
echo "  Database size: $DB_SIZE"

# Get record count (this might take a moment for large databases)
echo "  Counting records..."
RECORD_COUNT=$(sqlite3 "$DB_FILE" "SELECT COUNT(*) FROM bold WHERE nuc IS NOT NULL AND nuc != '' AND species IS NOT NULL AND species != '';")
echo "  Records to process: $RECORD_COUNT"

if [[ $RECORD_COUNT -eq 0 ]]; then
    echo "ERROR: No valid records found in database"
    exit 1
fi

echo ""

# Estimate runtime based on record count
if [[ $RECORD_COUNT -gt 0 ]]; then
    # Rough estimate: optimized script processes ~1000 records/second with good parallelization
    ESTIMATED_SECONDS=$((RECORD_COUNT / 1000 / PROCESSES * 2))  # Conservative estimate
    ESTIMATED_MINUTES=$((ESTIMATED_SECONDS / 60))
    echo "Estimated runtime: ~$ESTIMATED_MINUTES minutes"
    echo ""
fi

# Pre-flight checks
echo "Pre-flight checks:"

# Check available disk space
DISK_AVAIL=$(df -BG . | tail -1 | awk '{print $4}' | sed 's/G//')
echo "  Available disk space: ${DISK_AVAIL}GB"

if [[ $DISK_AVAIL -lt 10 ]]; then
    echo "  WARNING: Low disk space. Consider cleaning up temporary files."
fi

# Check if Parallel::ForkManager is available
perl -e "use Parallel::ForkManager; print 'Parallel::ForkManager: OK\n';" 2>/dev/null || {
    echo "  ERROR: Parallel::ForkManager not available"
    echo "  Installing via CPAN..."
    cpan -T Parallel::ForkManager || {
        echo "  Failed to install Parallel::ForkManager"
        exit 1
    }
}

# Check Digest::MD5
perl -e "use Digest::MD5; print 'Digest::MD5: OK\n';" 2>/dev/null || {
    echo "  ERROR: Digest::MD5 not available"
    exit 1
}

echo "  All dependencies OK"
echo ""

# Set up resource monitoring
echo "Starting resource monitoring..."
{
    while true; do
        echo "$(date): CPU: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}') | MEM: $(free -h | grep Mem | awk '{print $3"/"$2}') | LOAD: $(uptime | awk -F'load average:' '{print $2}')"
        sleep 60
    done
} > "logs/resource_monitor_${SLURM_JOB_ID}.log" &
MONITOR_PID=$!

# Main execution
echo "Starting haplotype analysis..."
echo "Command: perl assess_haplotypes_optimized.pl --db \"$DB_FILE\" --log INFO --processes $PROCESSES --chunk-size $CHUNK_SIZE --cache-size $CACHE_SIZE"
echo ""

START_TIME=$(date +%s)

# Run the optimized haplotype analysis
perl assess_haplotypes_optimized.pl \
    --db "$DB_FILE" \
    --log INFO \
    --processes $PROCESSES \
    --chunk-size $CHUNK_SIZE \
    --cache-size $CACHE_SIZE \
    > "$OUTPUT_FILE" \
    2> "$LOG_FILE"

EXIT_CODE=$?
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Stop resource monitoring
kill $MONITOR_PID 2>/dev/null

echo ""
echo "=================================================="
echo "EXECUTION COMPLETED"
echo "=================================================="
echo "End time: $(date)"
echo "Duration: $((DURATION / 60)) minutes $((DURATION % 60)) seconds"
echo "Exit code: $EXIT_CODE"

if [[ $EXIT_CODE -eq 0 ]]; then
    echo "Status: SUCCESS ✓"
    
    # Analyze results
    if [[ -f "$OUTPUT_FILE" ]]; then
        OUTPUT_LINES=$(wc -l < "$OUTPUT_FILE")
        ASSIGNMENTS=$((OUTPUT_LINES - 1))  # Subtract header
        echo "Results:"
        echo "  Output file: $OUTPUT_FILE"
        echo "  Haplotype assignments: $ASSIGNMENTS"
        echo "  Processing rate: $((RECORD_COUNT / DURATION)) records/second"
        
        # Show sample of results
        echo ""
        echo "Sample output (first 5 lines):"
        head -6 "$OUTPUT_FILE" | column -t -s $'\t'
        
        # Count unique haplotypes
        UNIQUE_HAPLOTYPES=$(tail -n +2 "$OUTPUT_FILE" | cut -f2 | sort -u | wc -l)
        echo ""
        echo "Summary statistics:"
        echo "  Total assignments: $ASSIGNMENTS"
        echo "  Unique haplotypes: $UNIQUE_HAPLOTYPES"
        echo "  Average sequences per haplotype: $((ASSIGNMENTS / UNIQUE_HAPLOTYPES))"
        
        # Disk usage
        OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
        echo "  Output file size: $OUTPUT_SIZE"
    else
        echo "WARNING: Output file not created"
        EXIT_CODE=1
    fi
else
    echo "Status: FAILED ✗"
    echo ""
    echo "Check the log file for details: $LOG_FILE"
    
    if [[ -f "$LOG_FILE" ]]; then
        echo ""
        echo "Last 20 lines of log:"
        tail -20 "$LOG_FILE"
    fi
fi

echo ""
echo "Log files:"
echo "  Main log: $LOG_FILE"
echo "  Resource monitor: logs/resource_monitor_${SLURM_JOB_ID}.log"
echo "  SLURM output: job_haplotype_${SLURM_JOB_ID}.out"
echo "  SLURM error: job_haplotype_${SLURM_JOB_ID}.err"

# Clean up temporary files
if [[ -n "$TMPDIR" && -d "$TMPDIR" ]]; then
    echo ""
    echo "Cleaning up temporary files in $TMPDIR..."
    find "$TMPDIR" -name "haplotype_*" -delete 2>/dev/null || true
fi

# Final system resource check
echo ""
echo "Final resource usage:"
echo "  Peak memory: $(grep VmPeak /proc/$$/status 2>/dev/null | awk '{print $2 $3}' || echo 'N/A')"
echo "  Disk space after: $(df -BG . | tail -1 | awk '{print $4}')GB available"

echo ""
echo "Haplotype analysis completed!"

exit $EXIT_CODE