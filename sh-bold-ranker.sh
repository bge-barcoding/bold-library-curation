#!/bin/bash
#SBATCH --partition=medium
#SBATCH --output=%j_out.out
#SBATCH --error=%j_err.err
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=example@email.com
#SBATCH --mail-type=ALL

# Activate conda env 

# OPTION 1 (remove #)
# source ~/.bashrc  
# conda activate bold-curation

# OPTION 2 (remove #)
# source activate bold-curation

# Cluster partitions for SLURM scheduler
# day / medium 
# week / long

# snakemake useful things
# --unlock
# --forcerun specify as "--forcerun rule" (run rule -> completion) or "rule --forcerun rule" (= run the rule then stop)

snakemake -s workflow/bold-ranker-array-phylo.smk -p -c 16

echo Complete!