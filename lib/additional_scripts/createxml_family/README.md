# createxml_family.pl

A Perl script for processing taxonomic data and generating organized XML files by taxonomic families from BOLD (Barcode of Life Data) system records.

## Overview

This script takes a tab-separated values (TSV) file containing species/specimen data and converts it into organized XML files grouped by taxonomic families. It provides special handling for specified subfamilies and tracks Barcode Index Numbers (BINs) that span multiple taxonomic groups.

## Features

- **Family-based Organization**: Creates separate XML files for each taxonomic family
- **Subfamily Support**: Special handling for specified subfamilies with genus-level organization
- **Data Validation**: Filters out incomplete or invalid records
- **BIN Tracking**: Identifies and reports BINs shared across multiple subfamilies
- **XML Formatting**: Converts tabular data to structured XML format with proper escaping

## Requirements

- Perl 5 with the following modules:
  - `strict`
  - `warnings` 
  - `Getopt::Long`
- Input files:
  - Main data file (TSV format)
  - `all_specs_and_syn.csv` (reference file for species validation)

## Usage

```bash
perl createxml_family.pl -in <input_file> -sub <subfamily_list> [-db <database_flag>]
```

### Parameters

- `-in` (required): Path to the input TSV file containing taxonomic data
- `-sub` (required): Comma-separated list of subfamily names that need special processing
- `-db` (optional): Flag to define if database files will be created (currently not fully implemented)

### Example

```bash
perl createxml_family.pl -in specimen_data.tsv -sub "Lepidoptera,Diptera,Coleoptera"
```

## Input File Format

The script expects a tab-separated file with the following requirements:

### Required Columns (by index)
- Column 9: BIN (Barcode Index Number) - must not be "None" or empty
- Column 20: Family name
- Column 21: Genus name  
- Column 24: Species name - must not be "None" or empty
- Column 94: Must be present (content validation)

### Data Quality Requirements
- Records with "None" values in key taxonomic fields are skipped
- Species must exist in the reference file `all_specs_and_syn.csv`
- Complete taxonomic hierarchy required for processing

## Output Structure

The script creates the following output structure:

```
family_output/
├── [Family1].xml
├── [Family2].xml
├── [SpecialSubfamily1]/
│   ├── [Genus1].xml
│   ├── [Genus2].xml
│   └── subfamily_shared_BINs.csv
└── [SpecialSubfamily2]/
    ├── [Genus3].xml
    └── subfamily_shared_BINs.csv
```

### File Types

1. **Family XML files**: Contains all records for families not in the subfamily list
2. **Genus XML files**: For special subfamilies, contains records organized by genus
3. **subfamily_shared_BINs.csv**: Lists BINs that appear in multiple genera within a subfamily

## XML Output Format

Each record is converted to XML format:

```xml
<records>
    <record>
        <Keep>0</Keep>
        <field_name>field_value</field_name>
        <!-- ... additional fields ... -->
    </record>
    <!-- ... more records ... -->
</records>
```

### Field Processing

- **Field Renaming**: 
  - `COLLECTORS` → `HAS_COLLECTOR`
  - `bold_recordset_code_arr` → `recordset_code_arr`
  - `country/ocean` → `country_ocean`
  - `province/state` → `province_state`

- **Character Escaping**:
  - `&` characters are converted to ` et ` when between names
  - `<` characters are converted to ` ;`

- **Excluded Fields**: The following fields are filtered out:
  - `nuc`
  - `elev_accuracy` 
  - `primers_forward`
  - `primers_reverse`

## BIN Analysis

The script tracks Barcode Index Numbers (BINs) and identifies cases where:
- A single BIN appears across multiple genera within a subfamily
- These shared BINs are recorded in `subfamily_shared_BINs.csv` files

CSV format: `BIN_ID;Genus1,Genus2;Subfamily`

## Error Handling

The script includes validation for:
- Missing required parameters
- File access permissions
- Data completeness checks
- Reference file availability

Common error messages:
- `"Input file is not defined. Use --in=FILENAME"`
- `"Subfamily is not defined. Use --sub=SUBFAMILY"`
- File access errors for input files

## Performance Notes

- The script includes a debugging limit (5000 records) that can be uncommented for testing
- Processing time depends on input file size and number of families
- Output directory is recreated on each run (existing `family_output` is removed)

## Dependencies

Ensure the reference file `all_specs_and_syn.csv` is present in the same directory as the script. This file contains valid species names used for filtering.

## Troubleshooting

1. **"define input" error**: Ensure both `-in` and `-sub` parameters are provided
2. **Missing output**: Check that input species exist in `all_specs_and_syn.csv`
3. **Empty XML files**: Verify input data quality and column positions
4. **Permission errors**: Ensure write permissions in the script directory

## License

Part of the bold-library-curation toolkit for processing BOLD database records.