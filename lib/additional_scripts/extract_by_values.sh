#!/bin/bash
#SBATCH --job-name=tsv_extract_values
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=tsv_extract_values_%j.out
#SBATCH --error=tsv_extract_values_%j.err
#SBATCH --partition=day
#SBATCH --mail-user=b.price@nhm.ac.uk
#SBATCH --mail-type=ALL

# HPC-optimized script for extracting records by multiple values from large TSV files
# This script uses extract_by_values.py to filter data based on a list of values

# =============================================================================
# CONFIGURATION - MODIFY THESE VARIABLES
# =============================================================================

# Input file path
INPUT_FILE="/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/BenPrice/bold_ranking/bold-library-curation/resources/BOLD_Public.16-May-2025.tsv"

# Column name to filter on (e.g., "processid", "order", "species_name", etc.)
COLUMN_NAME="processid"

# Option 1: Provide values as comma-separated string
# VALUES="ABC123,DEF456,GHI789,JKL012"

# Option 2: Provide path to file containing values (one per line)
VALUES_FILE="/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/BenPrice/bold_ranking/bold-library-curation/resources/chironominae_norway_processids.tsv"

# Output file path
OUTPUT_FILE="/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/BenPrice/bold_ranking/bold-library-curation/resources/extracted_chiro_${COLUMN_NAME}_${SLURM_JOB_ID}.tsv"

# Script path
SCRIPT_PATH="/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/BenPrice/bold_ranking/bold-library-curation/extract_by_values.py"

# Optional: Log file path
LOG_FILE="/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/BenPrice/bold_ranking/bold-library-curation/extraction_${SLURM_JOB_ID}.log"

# Performance tuning parameters
CHUNK_SIZE=10000
PROGRESS_INTERVAL=100000

# =============================================================================
# SCRIPT EXECUTION - DO NOT MODIFY BELOW THIS LINE
# =============================================================================

echo "=================================================="
echo "TSV Value Extraction Job Started"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=================================================="

# Load required modules (adjust based on your cluster)
# module load python/3.9
# module load gcc/11.2.0

# Create output and log directories
mkdir -p "$(dirname "$OUTPUT_FILE")"
mkdir -p "$(dirname "$LOG_FILE")"

# Build command
CMD="python3 $SCRIPT_PATH --input $INPUT_FILE --column $COLUMN_NAME --output $OUTPUT_FILE"

# Add values parameter (choose one)
if [ ! -z "$VALUES" ]; then
    CMD="$CMD --values \"$VALUES\""
elif [ ! -z "$VALUES_FILE" ]; then
    CMD="$CMD --values-file $VALUES_FILE"
else
    echo "ERROR: Either VALUES or VALUES_FILE must be set"
    exit 1
fi

# Add optional parameters
CMD="$CMD --chunk-size $CHUNK_SIZE --progress-interval $PROGRESS_INTERVAL"

if [ ! -z "$LOG_FILE" ]; then
    CMD="$CMD --log $LOG_FILE"
fi

echo "Command to execute:"
echo "$CMD"
echo ""

# Execute the extraction
eval $CMD

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "Job completed successfully!"
    echo "End time: $(date)"
    echo "Output file: $OUTPUT_FILE"
    
    # Display file information
    if [ -f "$OUTPUT_FILE" ]; then
        echo "Output file size: $(ls -lh "$OUTPUT_FILE" | awk '{print $5}')"
        echo "Output line count: $(wc -l < "$OUTPUT_FILE")"
    fi
    
    echo "=================================================="
else
    echo ""
    echo "=================================================="
    echo "Job failed with exit code: $?"
    echo "End time: $(date)"
    echo "Check error log for details."
    echo "=================================================="
    exit 1
fi
