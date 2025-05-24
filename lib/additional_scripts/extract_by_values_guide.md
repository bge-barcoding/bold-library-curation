# TSV Value Extractor - Usage Guide

## Overview
This enhanced toolkit provides an optimized solution for extracting records from large TSV files based on multiple values in any specified column. Unlike the original script that only filtered by the "order" column, this version allows flexible filtering on any column with a list of target values.

## Key Features
- **Flexible column filtering**: Filter on any column by name
- **Multiple value support**: Provide values via command line or file
- **High performance**: Optimized for HPC clusters with 19M+ row datasets
- **Memory efficient**: Chunked processing to handle large files
- **Progress tracking**: Real-time progress reporting
- **SLURM integration**: Ready for cluster job submission

## Files Included
1. `extract_by_values.py` - Main Python extraction script
2. `extract_by_values.sh` - SLURM job submission script
3. This usage guide

## Quick Start Examples

### 1. Extract by Process IDs (Command Line Values)
```bash
python3 extract_by_values.py \
    --input /path/to/data.tsv \
    --column "processid" \
    --values "ABC123,DEF456,GHI789" \
    --output /path/to/filtered_data.tsv
```

### 2. Extract by Order Values (From File)
Create a file `orders.txt` with one value per line:
```
Odonata
Lepidoptera
Diptera
Coleoptera
```

Then run:
```bash
python3 extract_by_values.py \
    --input /path/to/data.tsv \
    --column "order" \
    --values-file orders.txt \
    --output /path/to/filtered_orders.tsv \
    --log extraction.log
```

### 3. Extract by Species Names
```bash
python3 extract_by_values.py \
    --input /path/to/biodiversity.tsv \
    --column "species_name" \
    --values "Homo sapiens,Mus musculus,Drosophila melanogaster" \
    --output /path/to/target_species.tsv
```

## Command Line Options

### Required Parameters
- `--input, -i`: Path to input TSV file
- `--column, -c`: Name of the column to filter on
- `--output, -o`: Path to output TSV file
- **Values (choose one)**:
  - `--values, -v`: Comma-separated list of values
  - `--values-file, -f`: Path to file with values (one per line)

### Optional Parameters
- `--chunk-size`: Rows processed in memory at once (default: 10,000)
- `--progress-interval`: Progress logging frequency (default: 100,000)
- `--log`: Log file path (optional, logs to stdout if not specified)

## SLURM Job Submission

### 1. Configure the Job Script
Edit `extract_by_values.sh` and modify these variables:
```bash
INPUT_FILE="/path/to/your/data.tsv"
COLUMN_NAME="processid"
VALUES_FILE="/path/to/values.txt"  # OR set VALUES="val1,val2,val3"
OUTPUT_FILE="/path/to/output/extracted_${COLUMN_NAME}_${SLURM_JOB_ID}.tsv"
SCRIPT_PATH="/path/to/extract_by_values.py"
```

### 2. Submit the Job
```bash
sbatch extract_by_values.sh
```

### 3. Monitor Progress
```bash
# Check job status
squeue -u $USER

# View real-time output
tail -f tsv_extract_values_JOBID.out

# Check for errors
tail -f tsv_extract_values_JOBID.err
```
## Creating Values Files

### Format
Values files should contain one value per line:
```
ABC123
DEF456
GHI789
JKL012
```

### From Excel/CSV
If you have values in Excel or CSV format:
```bash
# Extract a column from CSV (column 2 in this example)
cut -d',' -f2 your_values.csv > values.txt

# Remove header if present
tail -n +2 values.txt > values_clean.txt
```

### From Database Query
```bash
# Example: MySQL query output
mysql -u user -p database_name -e "SELECT processid FROM samples WHERE status='active';" > values.txt
```

## Performance Guidelines

### Memory and Resource Allocation

| File Size | Recommended SLURM Settings | Expected Time |
|-----------|---------------------------|---------------|
| < 1M rows | 2GB RAM, 1 CPU, 30 min | 1-5 minutes |
| 1-5M rows | 4GB RAM, 2 CPU, 1 hour | 5-15 minutes |
| 5-20M rows | 8GB RAM, 4 CPU, 2 hours | 15-45 minutes |
| > 20M rows | 16GB RAM, 4-8 CPU, 4 hours | 1-2 hours |

### Chunk Size Optimization
- **Small files (<1M rows)**: `--chunk-size 50000`
- **Medium files (1-10M rows)**: `--chunk-size 10000` (default)
- **Large files (>10M rows)**: `--chunk-size 5000`

## Advanced Usage Examples

### 1. Process Multiple Columns Separately
```bash
#!/bin/bash
COLUMNS=("processid" "order" "family" "genus")
INPUT_FILE="/path/to/data.tsv"

for col in "${COLUMNS[@]}"; do
    python3 extract_by_values.py \
        --input "$INPUT_FILE" \
        --column "$col" \
        --values-file "${col}_values.txt" \
        --output "results/${col}_filtered.tsv" \
        --log "logs/${col}_extraction.log"
done
```

### 2. Extract Large Lists Efficiently
For very large value lists (>10,000 items):
```bash
# Split large values file into smaller chunks
split -l 1000 large_values.txt chunk_

# Process each chunk separately
for chunk in chunk_*; do
    python3 extract_by_values.py \
        --input data.tsv \
        --column "processid" \
        --values-file "$chunk" \
        --output "results/$(basename $chunk).tsv"
done

# Combine results (keeping header from first file only)
head -1 results/chunk_aa.tsv > combined_results.tsv
tail -n +2 -q results/chunk_*.tsv >> combined_results.tsv
```

### 3. Validate Results
```bash
#!/bin/bash
INPUT_FILE="data.tsv"
OUTPUT_FILE="filtered_data.tsv"
COLUMN_NAME="processid"

# Count unique values in target column of output
echo "Unique values found in output:"
cut -f$(head -1 "$INPUT_FILE" | tr '\t' '\n' | grep -n "$COLUMN_NAME" | cut -d: -f1) "$OUTPUT_FILE" | tail -n +2 | sort | uniq -c

# Verify no unexpected values
echo "Sample of values in output:"
cut -f$(head -1 "$INPUT_FILE" | tr '\t' '\n' | grep -n "$COLUMN_NAME" | cut -d: -f1) "$OUTPUT_FILE" | tail -n +2 | sort | uniq | head -10
```
## Troubleshooting

### Common Issues

1. **"Column not found"**
   ```bash
   # Check available columns
   head -1 input.tsv | tr '\t' '\n' | nl
   ```

2. **Empty output file**
   ```bash
   # Check if values exist in the column
   cut -f<column_number> input.tsv | grep -F -f values.txt | head -5
   ```

3. **Memory errors**
   - Reduce `--chunk-size` parameter
   - Increase SLURM memory allocation
   - Split input file into smaller pieces

4. **Performance issues**
   ```bash
   # Check file system performance
   time dd if=input.tsv of=/dev/null bs=1M count=100
   
   # Monitor job resources
   sstat -j $SLURM_JOB_ID --format=AveCPU,AvePages,AveRSS,AveVMSize
   ```

### Debugging Tips

1. **Test with small sample**
   ```bash
   # Create test file with first 1000 lines
   head -1000 large_file.tsv > test_file.tsv
   
   # Test extraction
   python3 extract_by_values.py -i test_file.tsv -c "processid" -v "test_value" -o test_output.tsv
   ```

2. **Verify column content**
   ```bash
   # Show unique values in column (first 20)
   cut -f<column_number> input.tsv | tail -n +2 | sort | uniq | head -20
   ```

3. **Check value format**
   ```bash
   # Look for whitespace or special characters
   cut -f<column_number> input.tsv | tail -n +2 | head -5 | xxd
   ```

## Output Verification

After extraction, verify results:
```bash
# Basic file information
ls -lh output.tsv
wc -l output.tsv

# Verify header preservation
head -1 output.tsv

# Check sample records
head -10 output.tsv

# Verify all records match criteria
COLUMN_NUM=$(head -1 output.tsv | tr '\t' '\n' | grep -n "processid" | cut -d: -f1)
cut -f$COLUMN_NUM output.tsv | tail -n +2 | sort | uniq -c
```

## Integration with Other Tools

### Converting to Other Formats
```bash
# Convert to CSV
python3 -c "
import csv
with open('output.tsv', 'r') as f_in, open('output.csv', 'w') as f_out:
    reader = csv.reader(f_in, delimiter='\t')
    writer = csv.writer(f_out)
    for row in reader:
        writer.writerow(row)
"

# Convert to JSON
python3 -c "
import csv, json
with open('output.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    data = list(reader)
with open('output.json', 'w') as f:
    json.dump(data, f, indent=2)
"
```

### Database Import
```bash
# Import to SQLite
sqlite3 results.db <<EOF
.mode tabs
.import output.tsv extracted_data
EOF

# Import to PostgreSQL
psql -d database_name -c "\COPY extracted_data FROM 'output.tsv' WITH DELIMITER E'\t' CSV HEADER;"
```

## Support and Best Practices

### Best Practices
1. **Always test with small samples first**
2. **Use appropriate chunk sizes for your file size**
3. **Monitor memory usage during runs**
4. **Keep backup copies of original data**
5. **Use descriptive output filenames with timestamps**
6. **Log all operations for reproducibility**

### Getting Help
- Check cluster documentation for SLURM-specific configurations
- Contact your HPC support team for resource limits and best practices
- Log files contain detailed timing and progress information for optimization
- Use the troubleshooting section above for common issues

### Example Workflows

**Workflow 1: Extract specific biological samples**
```bash
# 1. Get list of target process IDs from database
mysql -u user -p -e "SELECT processid FROM samples WHERE project='arctic_study';" > target_processids.txt

# 2. Extract matching records
python3 extract_by_values.py \
    --input BOLD_data.tsv \
    --column "processid" \
    --values-file target_processids.txt \
    --output arctic_samples.tsv \
    --log arctic_extraction.log

# 3. Verify results
echo "Extracted $(tail -n +2 arctic_samples.tsv | wc -l) samples"
```

**Workflow 2: Multi-stage filtering**
```bash
# 1. First filter by order
python3 extract_by_values.py -i full_data.tsv -c "order" -v "Lepidoptera" -o step1_lepidoptera.tsv

# 2. Then filter by family within that subset
python3 extract_by_values.py -i step1_lepidoptera.tsv -c "family" -v "Noctuidae,Sphingidae" -o final_moths.tsv
```

This enhanced script provides the flexibility you requested while maintaining the performance optimizations of the original. You can now filter on any column and provide values either as a comma-separated list or from a file containing one value per line.
