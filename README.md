![Perl CI](https://github.com/FabianDeister/Library_curation_BOLD/actions/workflows/ci.yml/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11267905.svg)](https://doi.org/10.5281/zenodo.11267905)
[![DOI](https://zenodo.org/badge/DOI/10.48546/WORKFLOWHUB.WORKFLOW.833.1.svg)](https://doi.org/10.48546/WORKFLOWHUB.WORKFLOW.833.1)

# Library curation BOLD

![alt text](doc/IBOL_LOGO_TRANSPARENT.png)

This repository contains scripts and synonymy data for pipelining the 
automated curation of [BOLD](https://boldsystems.org) data dumps in 
BCDM TSV format. The goal is to implement the classification of barcode 
reference sequences as is being developed by the 
[BGE](https://biodiversitygenomics.eu) consortium. A living document
in which these criteria are being developed is located
[here](https://docs.google.com/document/d/18m-7UnoJTG49TbvTsq_VncKMYZbYVbau98LE_q4rQvA/edit).

A further goal of this project is to develop the code in this repository
according to the standards developed by the community in terms of automation,
reproducibility, and provenance. In practice, this means including the
scripts in a pipeline system such as [snakemake](https://snakemake.readthedocs.io/),
adopting an environment configuration system such as
[conda](https://docs.conda.io/), and organizing the folder structure
in compliance with the requirements of
[WorkFlowHub](https://workflowhub.eu/). The latter will provide it with 
a DOI and will help generate [RO-crate](https://www.researchobject.org/ro-crate/)
documents, which means the entire tool chain is FAIR compliant according
to the current state of the art.

## Install
Clone the repo:
```{shell}
git clone https://github.com/bge-barcoding/bold-library-curation.git
```
Change directory: 
```{shell}
cd bold-library-curation
```
The code in this repo depends on various tools. These are managed using
the `mamba` program (a drop-in replacement of `conda`). The following
sets up an environment in which all needed tools are installed:

```{shell}
mamba env create -f environment.yml
```

Once set up, this is activated like so:

```{shell}
mamba activate bold-curation
```

## How to run

### snakemake

Follow the installation instructions above.

Update config/config.yml to define your input data.

in *bold-library-curation* type:
```{shell}
snakemake -p -c {number of cores}
```

If running on an HPC cluster with a SLURM scheduler you could use a bash script like this one:
```{shell}
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_curate_bold_%j.out
#SBATCH --error=job_curate_bold_%j.err
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4

source activate bold-curation

snakemake --unlock -p -c 4 clean

snakemake -p -c 4

# use "clean" when running again but backup the results folder first if you want to keep it
# I've added the unlock function to the clean because I'm tired of doing this manually when something goes wrong

echo Complete!
```
