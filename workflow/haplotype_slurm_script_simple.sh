#!/bin/bash
#SBATCH --job-name=haplotype-analysis
#SBATCH --partition=day
#SBATCH --output=job_haplotype_%j.out
#SBATCH --error=job_haplotype_%j.err
#SBATCH --mem=60G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --mail-user=b.price@nhm.ac.uk
#SBATCH --mail-type=ALL

echo "BOLD Haplotype Analysis - Optimized Version"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo ""

# Activate conda environment
source activate bold-curation

# Configuration
PROCESSES=$SLURM_CPUS_PER_TASK
CHUNK_SIZE=5000
CACHE_SIZE=100000

# Input/output paths
DB_FILE="results/bold.db"
OUTPUT_FILE="results/assessed_HAPLOTYPE_ID_optimized.tsv"
LOG_FILE="logs/haplotype_${SLURM_JOB_ID}.log"

# Create output directories
mkdir -p logs results

echo "Configuration:"
echo "  Database: $DB_FILE"
echo "  Output: $OUTPUT_FILE"  
echo "  Processes: $PROCESSES"
echo "  Chunk size: $CHUNK_SIZE"
echo ""

# Basic checks
if [[ ! -f "$DB_FILE" ]]; then
    echo "ERROR: Database file not found: $DB_FILE"
    exit 1
fi

# Get record count
RECORD_COUNT=$(sqlite3 "$DB_FILE" "SELECT COUNT(*) FROM bold WHERE nuc IS NOT NULL AND nuc != '' AND species IS NOT NULL AND species != '';")
echo "Records to process: $RECORD_COUNT"

if [[ $RECORD_COUNT -eq 0 ]]; then
    echo "ERROR: No valid records found in database"
    exit 1
fi

# Check dependencies
perl -e "use Parallel::ForkManager; use Digest::MD5;" 2>/dev/null || {
    echo "ERROR: Missing Perl dependencies (Parallel::ForkManager, Digest::MD5)"
    exit 1
}

echo ""
echo "Starting haplotype analysis..."

START_TIME=$(date +%s)

# Run the analysis
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

echo ""
echo "Execution completed"
echo "Duration: $((DURATION / 60)) minutes $((DURATION % 60)) seconds"
echo "Exit code: $EXIT_CODE"

if [[ $EXIT_CODE -eq 0 ]]; then
    echo "Status: SUCCESS"
    
    if [[ -f "$OUTPUT_FILE" ]]; then
        ASSIGNMENTS=$(( $(wc -l < "$OUTPUT_FILE") - 1 ))
        UNIQUE_HAPLOTYPES=$(tail -n +2 "$OUTPUT_FILE" | cut -f2 | sort -u | wc -l)
        
        echo "Results:"
        echo "  Haplotype assignments: $ASSIGNMENTS"
        echo "  Unique haplotypes: $UNIQUE_HAPLOTYPES"
        echo "  Processing rate: $((RECORD_COUNT / DURATION)) records/second"
        echo "  Output file: $OUTPUT_FILE"
    else
        echo "WARNING: Output file not created"
        EXIT_CODE=1
    fi
else
    echo "Status: FAILED"
    echo "Check log file: $LOG_FILE"
fi

echo "Completed at: $(date)"
exit $EXIT_CODE