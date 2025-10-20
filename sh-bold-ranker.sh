#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=%j_out.out
#SBATCH --error=%j_err.err
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=example@email.com
#SBATCH --mail-type=ALL

# Activate conda env 

# OPTION 1 (remove #)
# source ~/.bashrc  
# conda activate boldetective

# OPTION 2 (remove #)
# source activate boldetective

# Set PERL5LIB to include cpanm installations
export PERL5LIB=$CONDA_PREFIX/lib/perl5:$HOME/perl5/lib/perl5:$PERL5LIB

# Cluster partitions for SLURM scheduler
# day / medium 
# week / long

# snakemake useful things
# --unlock
# --forcerun specify as "--forcerun rule" (run rule -> completion) or "rule --forcerun rule" (= run the rule then stop)

snakemake -s workflow/bold-ranker-array-phylo.smk -p -c 16

echo Complete!