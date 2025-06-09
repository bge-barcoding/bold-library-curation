#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_phylo_parallel_%j.out
#SBATCH --error=job_phylo_parallel_%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00

# BOLD Phylogenetic Analysis - Parallel Processing
# ================================================
# This script runs the parallel phylogenetic analysis workflow
# It prepares batches and submits them as SLURM job arrays

echo "==============================================="
echo "BOLD PHYLOGENETIC ANALYSIS - PARALLEL MODE"
echo "==============================================="
echo "Started at: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo ""

source activate bold-curation

echo "Configuration:"
echo "  Conda environment: ${CONDA_DEFAULT_ENV}"
echo "  Working directory: $(pwd)"
echo "  Python version: $(python --version)"
echo ""

# Enable phylogenetic analysis in the workflow
export PHYLO_ENABLED=True

# Run Snakemake workflow with phylogenetic analysis rules
echo "Running Snakemake workflow with parallel phylogenetic analysis..."
echo ""

snakemake \
    -s workflow/bold-ranker-array-phylo.smk \
    -p \
    -c 4 \
    --target run_phylogenetic_analysis_parallel \
    --rerun-incomplete \
    --printshellcmds

SNAKEMAKE_EXIT_CODE=$?

echo ""
echo "==============================================="
if [ $SNAKEMAKE_EXIT_CODE -eq 0 ]; then
    echo "PHYLOGENETIC ANALYSIS WORKFLOW COMPLETED SUCCESSFULLY"
    echo ""
    echo "Results can be found in:"
    echo "  - results/phylogenies/ (phylogenetic trees and PDFs)"
    echo "  - results/phylogenies/phylogenetic_analysis_summary.json (summary report)"
    echo "  - logs/phylo_* (detailed logs)"
    echo ""
    echo "To view the summary, run:"
    echo "  python workflow/scripts/consolidate_phylo_results.py results/phylo_batches results/phylogenies"
else
    echo "WORKFLOW FAILED WITH EXIT CODE: $SNAKEMAKE_EXIT_CODE"
    echo "Check the logs for details:"
    echo "  - job_phylo_parallel_${SLURM_JOB_ID}.err"
    echo "  - logs/phylogenetic_analysis_parallel.log"
fi
echo "==============================================="
echo "Completed at: $(date)"

exit $SNAKEMAKE_EXIT_CODE
