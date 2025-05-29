#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_curate_bold_%j.out
#SBATCH --error=job_curate_bold_%j.err
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=example@email.com
#SBATCH --mail-type=ALL

source activate bold-curation

# use "clean" when running again but backup the results folder first if you want to keep it
# I've added the unlock function to the clean because I'm tired of doing this manually when something goes wrong

# snakemake --unlock -s workflow/Snakefile -p -c 4 clean

# snakemake -s workflow/Snakefile -p -c 4

snakemake -s workflow/Snakefile-faster -p -c 4

echo Complete!