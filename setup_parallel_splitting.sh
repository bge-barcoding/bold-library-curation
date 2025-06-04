#!/bin/bash
# Setup script for parallel family splitting
# Run this once to prepare the environment

echo "Setting up BOLD parallel family splitting..."

# Make scripts executable
chmod +x workflow/scripts/family_split_array.sh
chmod +x workflow/scripts/prepare_family_batches.py
chmod +x workflow/scripts/process_family_batch.py
chmod +x workflow/scripts/consolidate_results.py

echo "✓ Made scripts executable"

# Create necessary directories
mkdir -p logs
mkdir -p results/family_databases
mkdir -p results/family_batches

echo "✓ Created directories"

# Check if SLURM is available
if command -v sbatch &> /dev/null; then
    echo "✓ SLURM detected"
else
    echo "⚠ Warning: SLURM not found - parallel processing may not work"
fi

# Check Python availability
if command -v python &> /dev/null; then
    echo "✓ Python detected"
else
    echo "⚠ Warning: Python not found"
fi

echo ""
echo "Setup complete!"
echo ""
echo "To integrate with your existing pipeline:"
echo "1. Add the configuration from config/parallel_splitting_config.yml to your config/config.yml"
echo "2. Replace the existing split_families rule in bold-ranker.smk with the new rule from new_split_families_rule.smk"
echo "3. Test with a small database first"
echo ""
echo "Example usage:"
echo "  snakemake --cores 1 split_families"
echo ""
