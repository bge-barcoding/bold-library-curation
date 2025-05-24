# BOLD Record Reducer

A Perl script that filters BOLD (Barcode of Life Data) public database records based on a predefined species list.

## Overview

This script reduces a large BOLD public dataset by filtering records to include only species that are present in a reference species list. It's designed to work with BOLD database exports and create a more focused dataset for analysis.

## Files

- `reduce.pl` - Main Perl script
- `all_specs_and_syn.csv` - Input species reference file (semicolon-separated)
- `BOLD_Public.04-Apr-2025.tsv` - Input BOLD dataset (tab-separated)
- `reduced_BOLD_Public.04-Apr-2025.tsv` - Output filtered dataset

## How It Works

1. **Load Reference Species**: Reads `all_specs_and_syn.csv` and builds a hash of available species names
   - File format: semicolon-separated values
   - Each value becomes a key in the lookup hash

2. **Filter BOLD Records**: Processes the main BOLD dataset line by line
   - Preserves the original header row
   - Checks column 22 (index 21) for species names
   - Only outputs records where the species name exists in the reference list

3. **Generate Output**: Creates a reduced TSV file with only matching records

## Usage

```bash
perl reduce.pl
```

### Prerequisites

- Perl with standard libraries
- Input files must be present in the same directory:
  - `all_specs_and_syn.csv`
  - `BOLD_Public.04-Apr-2025.tsv`
### Input File Formats

**all_specs_and_syn.csv**:
```
species1;species2;species3;synonym1;synonym2
```

**BOLD_Public.04-Apr-2025.tsv**:
- Tab-separated values
- Species name expected in column 22 (0-indexed as column 21)
- Standard BOLD public dataset format

## Output

The script generates `reduced_BOLD_Public.04-Apr-2025.tsv` containing:
- Original header row
- Only records where column 22 matches a species in the reference list
- Same tab-separated format as input

## Notes

- The output file is recreated each time (existing file is deleted)
- Species matching is exact (case-sensitive)
- Column 22 is assumed to contain the species name for filtering
- No error handling for missing input files

## Data Processing Context

This script is part of the BOLD library curation workflow, designed to create focused datasets from the comprehensive BOLD public database by filtering for specific species of interest.