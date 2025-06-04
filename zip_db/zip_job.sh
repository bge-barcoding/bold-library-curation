#!/bin/bash
#SBATCH --job-name=zip_databases
#SBATCH --output=zip_databases_%j.out
#SBATCH --error=zip_databases_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=day
#SBATCH --mail-user=b.price@nhm.ac.uk
#SBATCH --mail-type=ALL

# Load required modules (adjust for your HPC system)
# module load python/3.9

# Set working directory
# cd $SLURM_SUBMIT_DIR

# Set paths (update these for your actual data locations)
SOURCE_DIR="/path/to/your/family_databases"
OUTPUT_DIR="/path/to/your/family_zipped"

echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: ${SLURM_MEM_PER_NODE}MB"
echo

# Run the compression script
python3 zip_databases.py \
    --source "$SOURCE_DIR" \
    --output "$OUTPUT_DIR" \
    --workers $SLURM_CPUS_PER_TASK

echo "Job completed on $(date)"
