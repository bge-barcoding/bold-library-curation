# Gap Analysis Script

## Overview

The `gap_analysis.py` script performs a comprehensive gap analysis for BOLD library curation by comparing:

1. **Input species list** - Your target taxa list (from `FILTER_TAXA_LIST` in config)
2. **Result output** - All BOLD records processed (`result_output.tsv`)
3. **BAGS assessment** - Quality grades for each taxonid (`assessed_BAGS.tsv`)

## What it Does

The script identifies:
- ✅ Coverage of species from your input list
- 🔍 Additional species found in BOLD results but not in your input list (missed species)
- 📊 BAGS grades (A-F quality assessment) for all species
- 🧬 BIN information and sharing species
- 📈 Record counts per taxonid
- 🌳 Full taxonomic hierarchy

## Output Format

The script produces a TSV file with the following columns:

| Column | Description |
|--------|-------------|
| `valid_species` | Species name (title case) |
| `synonyms` | Pipe-separated list of synonyms (if provided in input) |
| `missed_species` | "Yes" if found in results but not in input list; "No" otherwise |
| `BAGS_grade` | Quality grade (A, B, C, D, E, F) from BAGS assessment |
| `BIN_uri` | Barcode Index Numbers (pipe-separated if multiple) |
| `sharers` | Other species sharing the same BIN (pipe-separated) |
| `BIN_record_count` | Number of BOLD records for this taxonid |
| `kingdom` | Kingdom (from result_output.tsv) |
| `phylum` | Phylum (from result_output.tsv) |
| `class` | Class (from result_output.tsv) |
| `order` | Order (from result_output.tsv) |
| `family` | Family (from result_output.tsv) |
| `genus` | Genus (from result_output.tsv) |

## Input Species List Format

The script supports two formats:

### Simple Format (one species per line)
```
Acisoma ascalaphoides
Acisoma attenboroughi
Acisoma inflatum
```

### With Synonyms (semicolon-separated)
```
valid_species;synonym1;synonym2;synonym3
Chironomus acidophilus;Chironomus meigeni;Chironomus (Chironomus) acidophilus
Chironomus agilis
```

**Note:** Species name matching is **case-insensitive**.

## Usage Examples

### Using Config File

The script reads `FILTER_TAXA_LIST` from your config file:

```bash
python workflow/scripts/gap_analysis.py \
    --config config/config.yml \
    --result-output results/result_output.tsv \
    --assessed-bags results/assessed_BAGS.tsv \
    --output results/gap_analysis.tsv
```

### Override Species List via CLI

```bash
python workflow/scripts/gap_analysis.py \
    --config config/config.yml \
    --species-list custom_species.csv \
    --result-output results/result_output.tsv \
    --assessed-bags results/assessed_BAGS.tsv \
    --output results/gap_analysis.tsv
```

### Without Config File (specify species list directly)

```bash
python workflow/scripts/gap_analysis.py \
    --species-list resources/test_data/test_odonata_spec.csv \
    --result-output results/result_output.tsv \
    --assessed-bags results/assessed_BAGS.tsv \
    --output results/gap_analysis.tsv
```

### With Debug Logging

```bash
python workflow/scripts/gap_analysis.py \
    --config config/config.yml \
    --result-output results/result_output.tsv \
    --assessed-bags results/assessed_BAGS.tsv \
    --output results/gap_analysis.tsv \
    --log-level DEBUG
```

## Command-Line Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--config` | No* | Path to config.yml (to read FILTER_TAXA_LIST) |
| `--species-list` | No* | Path to species list CSV (overrides config) |
| `--result-output` | Yes | Path to result_output.tsv |
| `--assessed-bags` | Yes | Path to assessed_BAGS.tsv |
| `--output` | Yes | Path to output gap analysis TSV |
| `--log-level` | No | Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO) |

*Either `--config` or `--species-list` must be provided.

## How It Works

### 1. Load Input Species List
- Reads species names from CSV file
- Supports optional synonyms (semicolon-separated)
- Case-insensitive matching

### 2. Load Result Output
- Parses all BOLD records from `result_output.tsv`
- Groups by taxonid
- Counts records per taxonid
- Extracts taxonomic hierarchy

### 3. Load BAGS Assessment
- Links taxonid to BAGS grades
- Preserves BIN URIs and sharing species (pipe-separated format)

### 4. Merge & Analyze
- Links species names → taxonid → BAGS data
- Identifies species from input list that have BOLD records
- Flags "missed species" (in BOLD but not in input list)
- Handles species with multiple taxonids (creates separate rows)

### 5. Generate Report
- One row per unique species-taxonid combination
- Complete taxonomic hierarchy
- BAGS quality assessment
- BIN information and record counts

## Output Statistics


The script logs useful summary statistics:

```
Gap analysis complete: 215 total entries
  - Input list species: 200
  - Missed species: 15
  - Species with BAGS assessment: 180
```

## Example Scenarios

### Scenario 1: Species in Input List with BOLD Records
```
valid_species: Acisoma inflatum
synonyms: 
missed_species: No
BAGS_grade: E
BIN_uri: BOLD:ABU8638
sharers: Acisoma panorpoides,Acisoma variegatum,Acisoma sp.
BIN_record_count: 45
kingdom: Animalia
phylum: Arthropoda
class: Insecta
order: Odonata
family: Libellulidae
genus: Acisoma
```

### Scenario 2: Species with Synonyms
```
valid_species: Chironomus acidophilus
synonyms: Chironomus meigeni|Chironomus (Chironomus) acidophilus
missed_species: No
BAGS_grade: B
BIN_uri: BOLD:AAA1234
sharers: 
BIN_record_count: 12
...
```

### Scenario 3: Missed Species (in BOLD, not in input)
```
valid_species: Acisoma variegatum
synonyms: 
missed_species: Yes
BAGS_grade: E
BIN_uri: BOLD:ABU8638
sharers: Acisoma inflatum,Acisoma panorpoides,Acisoma sp.
BIN_record_count: 23
...
```

### Scenario 4: Input Species with No BOLD Records
```
valid_species: Rare species name
synonyms: 
missed_species: No
BAGS_grade: 
BIN_uri: 
sharers: 
BIN_record_count: 0
kingdom: 
phylum: 
class: 
order: 
family: 
genus: 
```

## Use Cases

1. **Coverage Assessment**: Identify which species from your target list have BOLD records
2. **Quality Review**: See BAGS grades for your target species
3. **Discovery**: Find additional related species in BOLD that you might want to include
4. **BIN Analysis**: Understand which species share Barcode Index Numbers
5. **Data Completeness**: Identify species needing more records or better quality data

## Technical Notes

- **Multiple taxonids per species**: If a species has multiple taxonids in result_output.tsv, the script creates separate rows for each
- **Case-insensitive matching**: "Acisoma inflatum" matches "ACISOMA INFLATUM" and "acisoma inflatum"
- **Pipe-separated format preserved**: BIN_uri and sharers maintain the "|" delimiter from assessed_BAGS.tsv
- **Empty values**: Missing data is represented as empty strings in the TSV output

## Integration with Workflow

This script can be run after the main BOLD curation workflow completes. It requires:
1. `result_output.tsv` (generated by the main workflow)
2. `assessed_BAGS.tsv` (generated by BAGS assessment step)
3. Your input species list (from `FILTER_TAXA_LIST` config)

## Troubleshooting

### File Not Found Errors
Ensure all input files exist and paths are correct:
```bash
ls -l config/config.yml
ls -l results/result_output.tsv
ls -l results/assessed_BAGS.tsv
```

### Empty Output
- Check that species names in input list match those in result_output.tsv
- Use `--log-level DEBUG` to see detailed processing information

### Missing BAGS Grades
- Some species may not have BAGS assessments if they don't meet certain criteria
- Check `assessed_BAGS.tsv` to verify which taxonids have assessments

## Dependencies

- Python 3.10+
- pandas (already in environment.yml)
- pyyaml (added to environment.yml for this script)

Install via conda:
```bash
conda env update -f environment.yml
```
