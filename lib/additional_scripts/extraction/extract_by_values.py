#!/usr/bin/env python3
"""
HPC-optimized script for extracting records with specific column values from large TSV files.
This version accepts a list of values and their corresponding column name for flexible filtering.
Designed for SLURM-managed clusters with 19M+ row datasets.

Usage:
    python extract_by_values.py --input /path/to/input.tsv --column "processid" --values "ABC123,DEF456,GHI789" --output /path/to/output.tsv
    python extract_by_values.py --input /path/to/input.tsv --column "order" --values-file values.txt --output /path/to/output.tsv
    
SLURM job submission example:
    sbatch --job-name=tsv_extract --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=02:00:00 extract_by_values.py
"""

import argparse
import csv
import sys
import os
import time
import logging
from pathlib import Path
from typing import Optional, List, Set
import mmap

# Configure logging for HPC environment
def setup_logging(log_file: Optional[str] = None) -> None:
    """Setup logging configuration optimized for HPC."""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    
    if log_file:
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

def get_file_line_count(file_path: str) -> int:
    """Efficiently count lines in large file using memory mapping."""
    try:
        with open(file_path, 'rb') as f:
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                return mm.count(b'\n')
    except Exception as e:
        logging.warning(f"Could not count lines efficiently: {e}. Using fallback method.")
        # Fallback method
        with open(file_path, 'r', encoding='utf-8') as f:
            return sum(1 for _ in f)

def find_column_index(header_row: List[str], column_name: str) -> int:
    """Find the index of the specified column in the header."""
    try:
        return header_row.index(column_name)
    except ValueError:
        # Case-insensitive search
        for i, col in enumerate(header_row):
            if col.lower() == column_name.lower():
                return i
        raise ValueError(f"No '{column_name}' column found in the header")

def load_values_from_file(file_path: str) -> Set[str]:
    """Load values from a text file (one value per line)."""
    values = set()
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                value = line.strip()
                if value:  # Skip empty lines
                    values.add(value)
        logging.info(f"Loaded {len(values)} values from {file_path}")
        return values
    except Exception as e:
        raise ValueError(f"Error reading values file {file_path}: {e}")

def parse_values_input(values_str: Optional[str], values_file: Optional[str]) -> Set[str]:
    """Parse values from either command line string or file."""
    if values_file:
        return load_values_from_file(values_file)
    elif values_str:
        # Split by comma and strip whitespace
        values = {v.strip() for v in values_str.split(',') if v.strip()}
        logging.info(f"Loaded {len(values)} values from command line")
        return values
    else:
        raise ValueError("Either --values or --values-file must be provided")

def extract_records_by_values(
    input_file: str,
    column_name: str,
    target_values: Set[str],
    output_file: str,
    chunk_size: int = 10000,
    progress_interval: int = 100000
) -> tuple[int, int]:
    """
    Extract records matching any of the specified values in the given column.
    
    Args:
        input_file: Path to input TSV file
        column_name: Name of the column to filter on
        target_values: Set of values to filter for
        output_file: Path to output TSV file
        chunk_size: Number of rows to process in memory at once
        progress_interval: How often to log progress
        
    Returns:
        Tuple of (total_rows_processed, matching_rows_found)
    """
    
    # Set CSV field size limit to handle large biological sequences
    # Use platform-safe maximum value to avoid Windows C long conversion issues
    try:
        csv.field_size_limit(sys.maxsize)
    except OverflowError:
        # On Windows, sys.maxsize might be too large for C long
        # Use a large but safe value (1GB)
        csv.field_size_limit(1024 * 1024 * 1024)
        logging.info("Using platform-safe CSV field size limit (1GB)")
    
    logging.info(f"Starting extraction for column: {column_name}")
    logging.info(f"Filtering for {len(target_values)} unique values")
    logging.info(f"Input file: {input_file}")
    logging.info(f"Output file: {output_file}")
    
    # Log sample of target values (first 10)
    sample_values = list(target_values)[:10]
    logging.info(f"Sample target values: {sample_values}")
    if len(target_values) > 10:
        logging.info(f"... and {len(target_values) - 10} more")
    
    # Get total line count for progress tracking
    total_lines = get_file_line_count(input_file)
    logging.info(f"Total lines in input file: {total_lines:,}")
    
    start_time = time.time()
    total_processed = 0
    matching_found = 0
    column_idx = None
    
    # Create output directory if it doesn't exist
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    try:
        with open(input_file, 'r', encoding='utf-8', newline='') as infile, \
             open(output_file, 'w', encoding='utf-8', newline='') as outfile:
            
            # Use csv.reader for proper TSV handling
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')
            
            # Process header
            try:
                header = next(reader)
                column_idx = find_column_index(header, column_name)
                writer.writerow(header)
                logging.info(f"Column '{column_name}' found at index: {column_idx}")
                total_processed += 1
            except StopIteration:
                raise ValueError("Input file is empty or has no header")
            
            # Process data rows in chunks
            chunk = []
            for row_num, row in enumerate(reader, start=2):  # Start at 2 since header is row 1
                
                # Skip empty rows
                if not row or len(row) <= column_idx:
                    continue
                
                chunk.append(row)
                
                # Process chunk when full
                if len(chunk) >= chunk_size:
                    matches = process_chunk(chunk, column_idx, target_values, writer)
                    matching_found += matches
                    total_processed += len(chunk)
                    chunk = []
                    
                    # Log progress
                    if total_processed % progress_interval == 0:
                        elapsed = time.time() - start_time
                        rate = total_processed / elapsed
                        progress_pct = (total_processed / total_lines) * 100
                        logging.info(
                            f"Progress: {total_processed:,}/{total_lines:,} ({progress_pct:.1f}%) "
                            f"- {matching_found:,} matches - {rate:.0f} rows/sec"
                        )
            
            # Process remaining chunk
            if chunk:
                matches = process_chunk(chunk, column_idx, target_values, writer)
                matching_found += matches
                total_processed += len(chunk)
    
    except Exception as e:
        logging.error(f"Error during extraction: {e}")
        raise
    
    elapsed = time.time() - start_time
    logging.info(f"Extraction completed in {elapsed:.2f} seconds")
    logging.info(f"Total rows processed: {total_processed:,}")
    logging.info(f"Matching rows found: {matching_found:,}")
    logging.info(f"Processing rate: {total_processed/elapsed:.0f} rows/sec")
    
    return total_processed, matching_found

def process_chunk(
    chunk: List[List[str]], 
    column_idx: int, 
    target_values: Set[str], 
    writer: csv.writer
) -> int:
    """Process a chunk of rows and write matching ones to output."""
    matches = 0
    for row in chunk:
        if len(row) > column_idx and row[column_idx].strip() in target_values:
            writer.writerow(row)
            matches += 1
    return matches

def validate_inputs(input_file: str, output_file: str, values_file: Optional[str] = None) -> None:
    """Validate input parameters and file accessibility."""
    
    # Check input file exists and is readable
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    if not os.access(input_file, os.R_OK):
        raise PermissionError(f"Cannot read input file: {input_file}")
    
    # Check values file if provided
    if values_file:
        if not os.path.exists(values_file):
            raise FileNotFoundError(f"Values file not found: {values_file}")
        if not os.access(values_file, os.R_OK):
            raise PermissionError(f"Cannot read values file: {values_file}")
    
    # Check output directory is writable
    output_dir = Path(output_file).parent
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise PermissionError(f"Cannot create output directory {output_dir}: {e}")
    
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_dir}")

def get_slurm_info() -> dict:
    """Get SLURM job information if running on SLURM cluster."""
    slurm_info = {}
    slurm_vars = [
        'SLURM_JOB_ID', 'SLURM_JOB_NAME', 'SLURM_NTASKS', 
        'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_NODE', 'SLURM_NODELIST'
    ]
    
    for var in slurm_vars:
        value = os.environ.get(var)
        if value:
            slurm_info[var] = value
    
    return slurm_info

def main():
    parser = argparse.ArgumentParser(
        description="Extract records with specified column values from large TSV files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract by process IDs (comma-separated)
  python extract_by_values.py -i data.tsv -c "processid" -v "ABC123,DEF456,GHI789" -o filtered_data.tsv
  
  # Extract by order values from a file
  python extract_by_values.py --input large_file.tsv --column "order" --values-file orders.txt --output results/filtered.tsv
  
  # Extract by any column with logging
  python extract_by_values.py -i data.tsv -c "species_name" -v "Homo sapiens,Mus musculus" -o results.tsv --log extraction.log
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to input TSV file'
    )
    
    parser.add_argument(
        '-c', '--column',
        required=True,
        help='Name of the column to filter on'
    )
    
    # Values input - either from command line or file
    values_group = parser.add_mutually_exclusive_group(required=True)
    values_group.add_argument(
        '-v', '--values',
        help='Comma-separated list of values to filter for'
    )
    values_group.add_argument(
        '-f', '--values-file',
        help='Path to file containing values (one per line)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to output TSV file'
    )
    
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=10000,
        help='Number of rows to process in memory at once (default: 10000)'
    )
    
    parser.add_argument(
        '--progress-interval',
        type=int,
        default=100000,
        help='How often to log progress (default: every 100000 rows)'
    )
    
    parser.add_argument(
        '--log',
        help='Path to log file (optional, logs to stdout if not specified)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log)
    
    # Log SLURM information if available
    slurm_info = get_slurm_info()
    if slurm_info:
        logging.info("SLURM Job Information:")
        for key, value in slurm_info.items():
            logging.info(f"  {key}: {value}")
    
    # Log system information
    logging.info(f"Python version: {sys.version}")
    logging.info(f"Working directory: {os.getcwd()}")
    
    try:
        # Validate inputs
        validate_inputs(args.input, args.output, args.values_file)
        
        # Parse target values
        target_values = parse_values_input(args.values, args.values_file)
        
        if not target_values:
            raise ValueError("No target values provided")
        
        # Perform extraction
        total_processed, matching_found = extract_records_by_values(
            input_file=args.input,
            column_name=args.column,
            target_values=target_values,
            output_file=args.output,
            chunk_size=args.chunk_size,
            progress_interval=args.progress_interval
        )
        
        # Summary
        logging.info("="*50)
        logging.info("EXTRACTION SUMMARY")
        logging.info("="*50)
        logging.info(f"Input file: {args.input}")
        logging.info(f"Target column: {args.column}")
        logging.info(f"Number of target values: {len(target_values)}")
        logging.info(f"Output file: {args.output}")
        logging.info(f"Total rows processed: {total_processed:,}")
        logging.info(f"Matching rows found: {matching_found:,}")
        
        if total_processed > 0:
            match_percentage = (matching_found / total_processed) * 100
            logging.info(f"Match percentage: {match_percentage:.2f}%")
        
        # Check output file size
        if os.path.exists(args.output):
            output_size = os.path.getsize(args.output)
            logging.info(f"Output file size: {output_size:,} bytes")
        
        logging.info("Extraction completed successfully!")
        
    except Exception as e:
        logging.error(f"Extraction failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
