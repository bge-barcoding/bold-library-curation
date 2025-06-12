#!/bin/bash
#SBATCH --partition=medium
#SBATCH --output=%j_curate_bold.out
#SBATCH --error=%j_curate_bold.err
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=example@email.com
#SBATCH --mail-type=ALL

# Activate conda env (Crop)
source ~/.bashrc  
conda activate bold-curation

# Activate conda env (NHM)
source activate bold-curation

# partitions
# day / medium and week / long

# snakemake useful things
# --unlock
# --forcerun specify as "--forcerun rule" (run rule -> completion) or "rule --forcerun rule" (= run the rule then stop)

snakemake -s workflow/bold-ranker-array-phylo.smk -p -c 16

echo Complete!