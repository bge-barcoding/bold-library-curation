#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_curate_bold_%j.out
#SBATCH --error=job_curate_bold_%j.err
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=example@email.com
#SBATCH --mail-type=ALL

source activate bold-curation

# snakemake useful things
# --unlock
# --forcerun specify as "--forcerun rule" (run rule -> completion) or "rule --forcerun rule" (= run the rule then stop)

snakemake -s bold-ranker.smk -p -c 16

echo Complete!