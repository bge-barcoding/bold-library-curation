#!/bin/bash
#SBATCH --job-name=bold_family_split
#SBATCH --output=logs/family_split_%A_%a.out
#SBATCH --error=logs/family_split_%A_%a.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=day

# SLURM Job Array Script for BOLD Family Database Splitting
# Usage: sbatch --array=0-N family_split_array.sh <source_db> <batch_dir> <output_dir> <threshold>

# Arguments
SOURCE_DB="$1"
BATCH_DIR="$2" 
OUTPUT_DIR="$3"
THRESHOLD="$4"
MAX_WORKERS="${5:-4}"  # Default to 4 workers per job

# Validate arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <source_db> <batch_dir> <output_dir> <threshold> [max_workers]"
    exit 1
fi

# Set up environment
echo "=== BOLD Family Splitting Job Array Task ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Start Time: $(date)"
echo "Source Database: ${SOURCE_DB}"
echo "Batch Directory: ${BATCH_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Threshold: ${THRESHOLD}"
echo "Max Workers: ${MAX_WORKERS}"
echo ""

# Find batch file for this array task
BATCH_FILE="${BATCH_DIR}/batch_${SLURM_ARRAY_TASK_ID}.json"

if [ ! -f "${BATCH_FILE}" ]; then
    echo "ERROR: Batch file not found: ${BATCH_FILE}"
    echo "Available batch files:"
    ls -la "${BATCH_DIR}"/batch_*.json 2>/dev/null || echo "No batch files found"
    exit 1
fi

echo "Processing batch file: ${BATCH_FILE}"

# Create output and log directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p logs

# Set log file for this task
LOG_FILE="${OUTPUT_DIR}/batch_${SLURM_ARRAY_TASK_ID}.log"

# Load conda environment if available
if [ -n "${CONDA_DEFAULT_ENV}" ]; then
    echo "Using conda environment: ${CONDA_DEFAULT_ENV}"
elif command -v conda &> /dev/null; then
    echo "Activating conda base environment"
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate base
fi

# Run the batch processor
echo "Starting family batch processing..."
python workflow/scripts/process_family_batch.py \
    "${BATCH_FILE}" \
    --source-db "${SOURCE_DB}" \
    --output-dir "${OUTPUT_DIR}" \
    --threshold "${THRESHOLD}" \
    --max-workers "${MAX_WORKERS}" \
    --log-file "${LOG_FILE}"

EXIT_CODE=$?

echo ""
echo "=== Job Completed ==="
echo "Exit Code: ${EXIT_CODE}"
echo "End Time: $(date)"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ Batch processing completed successfully"
else
    echo "✗ Batch processing failed"
fi

exit ${EXIT_CODE}
