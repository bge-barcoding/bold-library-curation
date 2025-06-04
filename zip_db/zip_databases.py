#!/usr/bin/env python3
"""
Efficient parallel database file compression script for HPC environments.
Zips all .db files while preserving directory structure.
"""

import os
import zipfile
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import time
from typing import List, Tuple

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('zip_databases.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def find_database_files(source_dir: str, extensions: List[str] = ['.db']) -> List[Tuple[str, str]]:
    """
    Find all database files in the source directory and calculate their target paths.
    
    Args:
        source_dir: Source directory to search
        extensions: List of file extensions to include
    
    Returns:
        List of tuples (source_file_path, target_zip_path)
    """
    database_files = []
    source_path = Path(source_dir)
    
    for file_path in source_path.rglob('*'):
        if file_path.is_file() and file_path.suffix.lower() in extensions:
            # Calculate relative path from source directory
            rel_path = file_path.relative_to(source_path)
            # Create target zip path with .zip extension
            target_zip = rel_path.with_suffix(rel_path.suffix + '.zip')
            database_files.append((str(file_path), str(target_zip)))
    
    return database_files


def zip_single_file(args: Tuple[str, str, str]) -> Tuple[str, bool, str]:
    """
    Zip a single database file.
    
    Args:
        args: Tuple of (source_file, target_zip_relative_path, output_base_dir)
    
    Returns:
        Tuple of (filename, success, error_message)
    """
    source_file, target_zip_rel, output_base_dir = args
    
    try:
        # Create full target path
        target_zip_path = Path(output_base_dir) / target_zip_rel
        
        # Create target directory if it doesn't exist
        target_zip_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create zip file
        with zipfile.ZipFile(target_zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as zip_file:
            # Add file to zip with just the filename (not full path)
            zip_file.write(source_file, Path(source_file).name)
        
        # Verify the zip file was created successfully
        if target_zip_path.exists() and target_zip_path.stat().st_size > 0:
            return (str(target_zip_path), True, "")
        else:
            return (str(target_zip_path), False, "Zip file not created or empty")
            
    except Exception as e:
        return (str(target_zip_path), False, str(e))


def main():
    parser = argparse.ArgumentParser(description='Efficiently zip database files in parallel')
    parser.add_argument('--source', 
                        default=r'C:\GitHub\bold-library-curation\results\results_otu_filter4\family_databases',
                        help='Source directory containing database files')
    parser.add_argument('--output', 
                        default=r'C:\GitHub\bold-library-curation\results\results_otu_filter4\family_zipped',
                        help='Output directory for zipped files')
    parser.add_argument('--workers', type=int, 
                        default=min(cpu_count(), 16),
                        help='Number of parallel workers (default: min(CPU count, 16))')
    parser.add_argument('--extensions', nargs='+', 
                        default=['.db'],
                        help='File extensions to zip (default: .db)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without actually zipping')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.exists(args.source):
        logger.error(f"Source directory does not exist: {args.source}")
        return 1
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    logger.info(f"Source directory: {args.source}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"File extensions: {args.extensions}")
    logger.info(f"Workers: {args.workers}")
    
    # Find all database files
    logger.info("Scanning for database files...")
    start_time = time.time()
    
    database_files = find_database_files(args.source, args.extensions)
    
    if not database_files:
        logger.warning("No database files found!")
        return 0
    
    logger.info(f"Found {len(database_files)} database files in {time.time() - start_time:.2f} seconds")
    
    if args.dry_run:
        logger.info("DRY RUN - Files that would be processed:")
        for source, target in database_files:
            logger.info(f"  {source} -> {os.path.join(args.output, target)}")
        return 0
    
    # Prepare arguments for parallel processing
    zip_args = [(source, target, args.output) for source, target in database_files]
    
    # Process files in parallel
    logger.info(f"Starting parallel compression with {args.workers} workers...")
    start_time = time.time()
    
    successful = 0
    failed = 0
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        # Submit all tasks
        future_to_file = {executor.submit(zip_single_file, arg): arg[0] for arg in zip_args}
        
        # Process completed tasks
        for future in as_completed(future_to_file):
            target_zip, success, error_msg = future.result()
            
            if success:
                successful += 1
                if successful % 10 == 0:  # Log progress every 10 files
                    logger.info(f"Compressed {successful}/{len(database_files)} files...")
            else:
                failed += 1
                logger.error(f"Failed to compress {future_to_file[future]}: {error_msg}")
    
    elapsed_time = time.time() - start_time
    
    # Summary
    logger.info("=" * 50)
    logger.info("COMPRESSION SUMMARY")
    logger.info("=" * 50)
    logger.info(f"Total files processed: {len(database_files)}")
    logger.info(f"Successfully compressed: {successful}")
    logger.info(f"Failed: {failed}")
    logger.info(f"Total time: {elapsed_time:.2f} seconds")
    logger.info(f"Average time per file: {elapsed_time/len(database_files):.3f} seconds")
    logger.info(f"Files per second: {len(database_files)/elapsed_time:.2f}")
    
    # Calculate compression statistics
    try:
        total_original_size = 0
        total_compressed_size = 0
        
        for source, target in database_files:
            if os.path.exists(source):
                total_original_size += os.path.getsize(source)
            
            compressed_path = os.path.join(args.output, target)
            if os.path.exists(compressed_path):
                total_compressed_size += os.path.getsize(compressed_path)
        
        if total_original_size > 0:
            compression_ratio = (1 - total_compressed_size / total_original_size) * 100
            logger.info(f"Original size: {total_original_size / (1024**2):.2f} MB")
            logger.info(f"Compressed size: {total_compressed_size / (1024**2):.2f} MB")
            logger.info(f"Compression ratio: {compression_ratio:.1f}%")
    
    except Exception as e:
        logger.warning(f"Could not calculate compression statistics: {e}")
    
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    exit(main())
