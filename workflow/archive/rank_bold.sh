#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_curate_bold_%j.out
#SBATCH --error=job_curate_bold_%j.err
#SBATCH --mem=80G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=example@email.com
#SBATCH --mail-type=ALL

source activate bold-curation

# use "clean" when running again but backup the results folder first if you want to keep it
# I've added the unlock function to the clean because I'm tired of doing this manually when something goes wrong

# unlock implementation - uncomment to run
# snakemake --unlock -s workflow/Snakefile -p -c 10 clean

# original implementation - uncomment to run
# snakemake -s workflow/Snakefile -p -c 10

# faster implementation with BAGS grading - uncomment to run
snakemake -s workflow/Snakefile-even-faster-bags -p -c 10

echo Complete!