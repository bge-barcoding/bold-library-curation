# Database File Compression Utility

This utility efficiently compresses database files in parallel while preserving directory structure. It's optimized for HPC environments and can handle hundreds of files across multiple subdirectories.

## Features

- **Parallel Processing**: Uses multiple CPU cores for fast compression
- **Structure Preservation**: Maintains exact directory structure in output
- **Progress Monitoring**: Real-time progress updates and detailed logging
- **Compression Statistics**: Reports file sizes and compression ratios
- **Error Handling**: Robust error handling with detailed logging
- **HPC Ready**: Includes SLURM job script for cluster environments
- **Flexible**: Configurable file extensions, worker count, and paths

## Files

- `zip_databases.py` - Main Python script
- `run_zip.bat` - Windows batch script
- `run_zip.sh` - Linux/Unix shell script  
- `submit_zip_job.slurm` - SLURM job submission script
- `zip_databases.log` - Log file (created during execution)

## Quick Start

### Windows
```cmd
python zip_databases.py
```
or
```cmd
run_zip.bat
```

### Linux/HPC
```bash
python3 zip_databases.py
```
or
```bash
./run_zip.sh
```

### SLURM Cluster
```bash
sbatch submit_zip_job.slurm
```

## Command Line Options

```
python zip_databases.py [options]

Options:
  --source PATH         Source directory (default: family_databases)
  --output PATH         Output directory (default: family_zipped)
  --workers N           Number of parallel workers (default: min(CPU count, 16))
  --extensions EXT...   File extensions to compress (default: .db)
  --dry-run            Show what would be done without compressing
  -h, --help           Show help message
```

## Examples

### Basic usage with custom paths:
```bash
python3 zip_databases.py \
    --source /data/databases \
    --output /data/compressed \
    --workers 8
```

### Multiple file extensions:
```bash
python3 zip_databases.py \
    --extensions .db .sqlite .sqlite3 \
    --workers 12
```

### Dry run to preview operations:
```bash
python3 zip_databases.py --dry-run
```

## Performance Results (Test Dataset)

- **Files processed**: 16 database files
- **Original size**: 16.11 MB
- **Compressed size**: 1.37 MB  
- **Compression ratio**: 91.5%
- **Processing time**: 0.25 seconds
- **Throughput**: 64 files/second

## HPC Configuration

For large datasets on HPC systems, adjust the SLURM script parameters:

```bash
#SBATCH --cpus-per-task=32    # More CPUs for large datasets
#SBATCH --mem=64G             # More memory if needed
#SBATCH --time=04:00:00       # Longer time limit
```

## Requirements

- Python 3.6+
- Standard library only (no external dependencies)
- Sufficient disk space for compressed files

## Error Handling

- Failed compressions are logged with detailed error messages
- Partial failures don't stop the overall process
- Exit codes: 0 = success, 1 = some failures occurred

## Logging

All operations are logged to both console and `zip_databases.log` file:
- File discovery progress
- Compression progress (every 10 files)
- Error messages for failed operations
- Final summary with statistics

## Notes

- Each database file is compressed individually (not grouped)
- Directory structure is exactly preserved
- Compressed files get `.zip` extension added to original name
- For very large datasets (1000+ files), consider using higher worker counts
- The script automatically creates output directories as needed
