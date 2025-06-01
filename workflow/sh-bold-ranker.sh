#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_curate_bold_%j.out
#SBATCH --error=job_curate_bold_%j.err
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=b.price@nhm.ac.uk
#SBATCH --mail-type=ALL

source activate bold-curation

# use "clean" when running again but backup the results folder first if you want to keep it
# I've added the unlock function to the clean because I'm tired of doing this manually when something goes wrong

# snakemake useful things
# --unlock
# --forcerun specify as --forcerun rule (rule -> completion) or rule --forcerun rule (= only the rule then stop)

snakemake -s workflow/bold-ranker.smk -p -c 16

echo Complete!