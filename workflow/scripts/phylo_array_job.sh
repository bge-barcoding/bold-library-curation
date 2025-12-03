#!/bin/bash
#SBATCH --job-name=phylo_parallel
#SBATCH --output=logs/phylo_%A_%a.out
#SBATCH --error=logs/phylo_%A_%a.err
#SBATCH --array=1-50%10
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=day

# partitions day/medium or week/long (NHM/Crop)

# Phylogenetic Analysis Parallel Job Array
# =========================================
# Processes families in parallel using SLURM job arrays
# Each job processes a batch of families defined in JSON files

echo "==============================================="
echo "PHYLOGENETIC ANALYSIS - ARRAY JOB ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "==============================================="
echo "Started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "CPUs allocated: ${SLURM_CPUS_PER_TASK}"
echo "Memory allocated: ${SLURM_MEM_PER_NODE}MB"
echo "Working directory: $(pwd)"
echo ""

# Set up error handling
set -e
set -u

# Function to handle errors
error_exit() {
    echo "ERROR: $1" >&2
    echo "Job failed at: $(date)" >&2
    exit 1
}

# Function to handle cleanup on exit
cleanup() {
    echo ""
    echo "Cleaning up temporary files..."
    
    # Stop virtual display if we started it
    if [[ -n "${XVFB_PID:-}" ]]; then
        echo "Stopping Xvfb process $XVFB_PID..."
        kill $XVFB_PID >/dev/null 2>&1 || true
    fi
    
    echo "Job completed at: $(date)"
}
trap cleanup EXIT

# Validate required environment variables
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    error_exit "SLURM_ARRAY_TASK_ID not set"
fi

# Configuration from command line arguments or defaults
DATABASE_PATH="${1:-results/bold.db}"
BATCH_DIR="${2:-results/phylo_batches}"
OUTPUT_DIR="${3:-results/phylogenies}"
MIN_OTUS="${4:-3}"
ALIGNMENT_METHOD="${5:-mafft}"
TREE_METHOD="${6:-fasttree}"
BOOTSTRAP="${7:-1000}"
NUM_OUTGROUPS="${8:-3}"
GENERATE_PDFS="${9:-false}"
CUSTOM_PARAMS_FILE="${10:-}"
LOG_DIR="${11:-logs}"
CLEANUP_INTERMEDIATES="${12:-true}"

echo "Configuration:"
echo "  Database: ${DATABASE_PATH}"
echo "  Batch directory: ${BATCH_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  Min OTUs: ${MIN_OTUS}"
echo "  Alignment method: ${ALIGNMENT_METHOD}"
echo "  Tree method: ${TREE_METHOD}"
echo "  Bootstrap: ${BOOTSTRAP}"
echo "  Outgroups: ${NUM_OUTGROUPS}"
echo "  Generate PDFs: ${GENERATE_PDFS}"
echo "  Custom parameters file: ${CUSTOM_PARAMS_FILE:-'none (using optimized defaults)'}"
echo "  Log directory: ${LOG_DIR}"
echo "  Cleanup intermediates: ${CLEANUP_INTERMEDIATES}"
echo ""

# Validate inputs
if [[ ! -f "${DATABASE_PATH}" ]]; then
    error_exit "Database file not found: ${DATABASE_PATH}"
fi

if [[ ! -d "${BATCH_DIR}" ]]; then
    error_exit "Batch directory not found: ${BATCH_DIR}"
fi

# Construct batch file name (zero-padded to 3 digits)
BATCH_FILE="${BATCH_DIR}/phylo_batch_$(printf "%03d" ${SLURM_ARRAY_TASK_ID}).json"

if [[ ! -f "${BATCH_FILE}" ]]; then
    error_exit "Batch file not found: ${BATCH_FILE}"
fi

echo "Processing batch file: ${BATCH_FILE}"

# Load batch information
BATCH_INFO=$(cat "${BATCH_FILE}")
FAMILY_COUNT=$(echo "${BATCH_INFO}" | python3 -c "import sys, json; data=json.load(sys.stdin); print(data['family_count'])")

echo "Batch contains ${FAMILY_COUNT} families"
echo ""

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}/phylo_checkpoints"

# Set up checkpoint file
CHECKPOINT_FILE="${LOG_DIR}/phylo_checkpoints/checkpoint_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"

# Activate conda environment
echo "Activating conda environment..."
if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "  Already in conda environment: ${CONDA_DEFAULT_ENV}"
else
    if command -v conda >/dev/null 2>&1; then
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate bold-curation || error_exit "Failed to activate conda environment 'bold-curation'"
        echo "  Activated conda environment: bold-curation"
    else
        error_exit "Conda not found in PATH"
    fi
fi

# Set environment variables for headless operation (for PDF generation)
export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
export DISPLAY=:99

# Additional ETE3/Qt environment variables for headless operation
export QT_DEBUG_PLUGINS=0
export QTWEBENGINE_DISABLE_SANDBOX=1

# Start virtual display if needed for ETE3 (suppress errors if already running)
echo "Setting up virtual display for ETE3..."
if command -v Xvfb >/dev/null 2>&1; then
    # Check if virtual display is already running
    if ! xdpyinfo -display :99 >/dev/null 2>&1; then
        echo "Starting Xvfb virtual display..."
        Xvfb :99 -screen 0 1024x768x24 >/dev/null 2>&1 &
        XVFB_PID=$!
        sleep 2  # Give Xvfb time to start
        echo "Xvfb started with PID: $XVFB_PID"
    else
        echo "Virtual display :99 already running"
    fi
else
    echo "Warning: Xvfb not available, ETE3 PDF generation may fail"
fi

# Check required tools
echo "Checking required tools..."
TOOLS_OK=true

check_tool() {
    local tool=$1
    local package=$2
    if command -v "${tool}" >/dev/null 2>&1; then
        echo "  ✓ ${tool} found"
    else
        echo "  ✗ ${tool} not found (install with: conda install -c bioconda ${package})"
        TOOLS_OK=false
    fi
}

check_tool "python3" "python"

# Check alignment tool
if [[ "${ALIGNMENT_METHOD}" == "mafft" ]]; then
    check_tool "mafft" "mafft"
elif [[ "${ALIGNMENT_METHOD}" == "muscle" ]]; then
    check_tool "muscle" "muscle"
fi

# Check tree building tool
if [[ "${TREE_METHOD}" == "iqtree" ]]; then
    check_tool "iqtree" "iqtree"
elif [[ "${TREE_METHOD}" == "fasttree" ]]; then
    if command -v fasttree >/dev/null 2>&1; then
        echo "  ✓ fasttree found"
    elif command -v FastTree >/dev/null 2>&1; then
        echo "  ✓ FastTree found"
    else
        echo "  ✗ fasttree/FastTree not found (install with: conda install -c bioconda fasttree)"
        TOOLS_OK=false
    fi
fi

if [[ "${TOOLS_OK}" != "true" ]]; then
    error_exit "Required tools missing"
fi

echo ""

# Prepare command arguments
PHYLO_ARGS=(
    "--database" "${DATABASE_PATH}"
    "--output-dir" "${OUTPUT_DIR}"
    "--family-list" "${BATCH_FILE}"
    "--min-otus" "${MIN_OTUS}"
    "--alignment-method" "${ALIGNMENT_METHOD}"
    "--tree-method" "${TREE_METHOD}"
    "--bootstrap" "${BOOTSTRAP}"
    "--num-outgroups" "${NUM_OUTGROUPS}"
    "--batch-id" "${SLURM_ARRAY_TASK_ID}"
    "--checkpoint-file" "${CHECKPOINT_FILE}"
)

# Add custom parameters file if provided
if [[ -n "${CUSTOM_PARAMS_FILE}" && -f "${CUSTOM_PARAMS_FILE}" ]]; then
    PHYLO_ARGS+=("--custom-parameters" "${CUSTOM_PARAMS_FILE}")
    echo "Using custom parameters from: ${CUSTOM_PARAMS_FILE}"
fi

# Add PDF generation if requested
if [[ "${GENERATE_PDFS}" == "true" ]] || [[ "${GENERATE_PDFS}" == "True" ]]; then
    PHYLO_ARGS+=("--generate-pdfs")
fi

# Add cleanup intermediates if requested
if [[ "${CLEANUP_INTERMEDIATES}" == "true" ]] || [[ "${CLEANUP_INTERMEDIATES}" == "True" ]]; then
    PHYLO_ARGS+=("--cleanup-intermediates")
fi

# Add BIN conflict analysis (always enabled for now)
PHYLO_ARGS+=("--bin-conflict-analysis")

echo "Running phylogenetic analysis..."
echo "Command: python3 workflow/scripts/phylo_pipeline.py ${PHYLO_ARGS[*]}"
echo ""

# Run the phylogenetic analysis
START_TIME=$(date +%s)

python3 workflow/scripts/phylo_pipeline.py "${PHYLO_ARGS[@]}" || error_exit "Phylogenetic analysis failed"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo ""
echo "==============================================="
echo "BATCH PROCESSING COMPLETE"
echo "==============================================="
echo "Batch: ${SLURM_ARRAY_TASK_ID}"
echo "Families in batch: ${FAMILY_COUNT}"
echo "Runtime: ${RUNTIME} seconds ($(date -d@${RUNTIME} -u +%H:%M:%S))"
echo "Completed at: $(date)"

# Save batch completion info
COMPLETION_INFO="{
  \"job_id\": \"${SLURM_JOB_ID}\",
  \"array_job_id\": \"${SLURM_ARRAY_JOB_ID}\",
  \"array_task_id\": \"${SLURM_ARRAY_TASK_ID}\",
  \"batch_file\": \"${BATCH_FILE}\",
  \"family_count\": ${FAMILY_COUNT},
  \"runtime_seconds\": ${RUNTIME},
  \"completed_at\": \"$(date -Iseconds)\",
  \"node\": \"$(hostname)\",
  \"success\": true
}"

echo "${COMPLETION_INFO}" > "${LOG_DIR}/phylo_checkpoints/completion_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.json"

echo "Batch completion info saved to ${LOG_DIR}/phylo_checkpoints/"
echo "Job array task ${SLURM_ARRAY_TASK_ID} completed successfully!"
