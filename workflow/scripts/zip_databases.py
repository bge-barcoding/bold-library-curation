#!/usr/bin/env python3
"""
Efficient parallel database file compression script for HPC environments.
Supports two modes:
1. Individual file compression: Zips each .db file separately
2. Family directory compression: Zips entire family directories containing .db files and related materials
"""

import os
import zipfile
import argparse
from pathlib import Path
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import time
from typing import List, Tuple, Dict, Union

# Logger will be configured in main() after parsing arguments
logger = logging.getLogger(__name__)


def find_database_files(source_dir: str, extensions: List[str] = ['.db']) -> List[Tuple[str, str]]:
    """
    Find all database files in the source directory and calculate their target paths.
    Used for individual file compression mode.
    
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


def find_family_directories(source_dir: str, extensions: List[str] = ['.db']) -> List[Tuple[str, str, str]]:
    """
    Find all directories containing database files (family directories).
    Used for family directory compression mode.
    
    Args:
        source_dir: Source directory to search
        extensions: List of file extensions that identify family directories
    
    Returns:
        List of tuples (family_directory_path, family_name, target_zip_path)
    """
    family_dirs = []
    source_path = Path(source_dir)
    processed_dirs = set()
    
    # Find all files with the specified extensions
    for file_path in source_path.rglob('*'):
        if file_path.is_file() and file_path.suffix.lower() in extensions:
            family_dir = file_path.parent
            
            # Skip if we've already processed this directory
            if family_dir in processed_dirs:
                continue
            
            processed_dirs.add(family_dir)
            
            # Use the directory name as the family name
            family_name = family_dir.name
            
            # Calculate relative path from source directory to preserve hierarchy
            rel_path = family_dir.relative_to(source_path)
            
            # Create target zip path preserving the directory structure
            # Replace the family directory name with family_name.zip
            target_zip_dir = rel_path.parent
            target_zip = target_zip_dir / f"{family_name}.zip"
            
            family_dirs.append((str(family_dir), family_name, str(target_zip)))
    
    return family_dirs


def zip_single_file(args: Tuple[str, str, str]) -> Tuple[str, bool, str, Dict[str, Union[int, float]]]:
    """
    Zip a single database file.
    Used for individual file compression mode.
    
    Args:
        args: Tuple of (source_file, target_zip_relative_path, output_base_dir)
    
    Returns:
        Tuple of (filename, success, error_message, stats_dict)
    """
    source_file, target_zip_rel, output_base_dir = args
    
    try:
        # Create full target path
        target_zip_path = Path(output_base_dir) / target_zip_rel
        
        # Create target directory if it doesn't exist
        target_zip_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Get original file size
        original_size = Path(source_file).stat().st_size
        
        # Create zip file
        with zipfile.ZipFile(target_zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as zip_file:
            # Add file to zip with just the filename (not full path)
            zip_file.write(source_file, Path(source_file).name)
        
        # Verify the zip file was created successfully
        if target_zip_path.exists() and target_zip_path.stat().st_size > 0:
            compressed_size = target_zip_path.stat().st_size
            stats = {
                'original_size': original_size,
                'compressed_size': compressed_size,
                'files_count': 1
            }
            return (str(target_zip_path), True, "", stats)
        else:
            return (str(target_zip_path), False, "Zip file not created or empty", {})
            
    except Exception as e:
        return (str(target_zip_path), False, str(e), {})


def zip_family_directory(args: Tuple[str, str, str, str]) -> Tuple[str, bool, str, Dict[str, Union[int, float]]]:
    """
    Zip an entire family directory including all files.
    Used for family directory compression mode.
    
    Args:
        args: Tuple of (family_directory_path, family_name, target_zip_relative_path, output_base_dir)
    
    Returns:
        Tuple of (zip_filename, success, error_message, stats_dict)
    """
    family_dir_path, family_name, target_zip_rel, output_base_dir = args
    
    try:
        family_dir = Path(family_dir_path)
        target_zip_path = Path(output_base_dir) / target_zip_rel
        
        # Create output directory if it doesn't exist (preserving hierarchy)
        target_zip_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Collect statistics
        total_original_size = 0
        files_count = 0
        file_types = {}
        
        # Get all files in the family directory
        all_files = []
        for file_path in family_dir.rglob("*"):
            if file_path.is_file():
                all_files.append(file_path)
                file_size = file_path.stat().st_size
                total_original_size += file_size
                files_count += 1
                
                # Track file types
                ext = file_path.suffix.lower()
                file_types[ext] = file_types.get(ext, 0) + 1
        
        if files_count == 0:
            return (str(target_zip_path), False, "No files found in family directory", {})
        
        # Create zip file
        with zipfile.ZipFile(target_zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as zip_file:
            for file_path in all_files:
                # Store files with relative path from family directory
                arcname = file_path.relative_to(family_dir)
                zip_file.write(file_path, arcname)
        
        # Verify the zip file was created successfully
        if target_zip_path.exists() and target_zip_path.stat().st_size > 0:
            compressed_size = target_zip_path.stat().st_size
            stats = {
                'original_size': total_original_size,
                'compressed_size': compressed_size,
                'files_count': files_count,
                'file_types': file_types,
                'family_name': family_name
            }
            return (str(target_zip_path), True, "", stats)
        else:
            return (str(target_zip_path), False, "Zip file not created or empty", {})
            
    except Exception as e:
        return (str(target_zip_path), False, str(e), {})


def main():
    parser = argparse.ArgumentParser(description='Efficiently zip database files in parallel')
    parser.add_argument('--source', 
                        default=r'C:\GitHub\bold-library-curation\results\results_otu_filter4\family_databases',
                        help='Source directory containing database files')
    parser.add_argument('--output', 
                        default=r'C:\GitHub\bold-library-curation\results\results_otu_filter4\family_zipped',
                        help='Output directory for zipped files')
    parser.add_argument('--log-dir',
                        default='.',
                        help='Directory for log files (default: current directory)')
    parser.add_argument('--workers', type=int, 
                        default=min(cpu_count(), 16),
                        help='Number of parallel workers (default: min(CPU count, 16))')
    parser.add_argument('--extensions', nargs='+', 
                        default=['.db'],
                        help='File extensions to zip (default: .db)')
    parser.add_argument('--mode', 
                        choices=['individual_files', 'family_directories'],
                        default='individual_files',
                        help='Compression mode: individual_files (zip each file separately) or family_directories (zip entire family folders)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without actually zipping')
    
    args = parser.parse_args()
    
    # Setup logging with configurable log directory
    log_file = os.path.join(args.log_dir, 'zip_databases.log')
    os.makedirs(args.log_dir, exist_ok=True)  # Ensure log directory exists
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    # Validate input directory
    if not os.path.exists(args.source):
        logger.error(f"Source directory does not exist: {args.source}")
        return 1
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    logger.info(f"Compression mode: {args.mode}")
    logger.info(f"Source directory: {args.source}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"File extensions: {args.extensions}")
    logger.info(f"Workers: {args.workers}")
    
    # Find files or directories to process based on mode
    logger.info("Scanning for files/directories to compress...")
    start_time = time.time()
    
    if args.mode == 'individual_files':
        items_to_process = find_database_files(args.source, args.extensions)
        process_function = zip_single_file
        item_type = "database files"
    else:  # family_directories
        items_to_process = find_family_directories(args.source, args.extensions)
        process_function = zip_family_directory
        item_type = "family directories"
    
    if not items_to_process:
        logger.warning(f"No {item_type} found!")
        return 0
    
    logger.info(f"Found {len(items_to_process)} {item_type} in {time.time() - start_time:.2f} seconds")
    
    if args.dry_run:
        logger.info(f"DRY RUN - {item_type.title()} that would be processed:")
        for item in items_to_process:
            if args.mode == 'individual_files':
                source, target = item
                logger.info(f"  {source} -> {os.path.join(args.output, target)}")
            else:  # family_directories
                family_dir, family_name, target_zip = item
                logger.info(f"  {family_dir} ({family_name}) -> {os.path.join(args.output, target_zip)}")
        return 0
    
    # Prepare arguments for parallel processing
    if args.mode == 'individual_files':
        process_args = [(source, target, args.output) for source, target in items_to_process]
    else:  # family_directories
        process_args = [(family_dir, family_name, target_zip, args.output) 
                       for family_dir, family_name, target_zip in items_to_process]
    
    # Process items in parallel
    logger.info(f"Starting parallel compression with {args.workers} workers...")
    start_time = time.time()
    
    successful = 0
    failed = 0
    all_stats = []
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        # Submit all tasks
        future_to_item = {executor.submit(process_function, arg): arg for arg in process_args}
        
        # Process completed tasks
        for future in as_completed(future_to_item):
            result_path, success, error_msg, stats = future.result()
            
            if success:
                successful += 1
                all_stats.append(stats)
                if successful % 50 == 0:  # Log progress every 50 items
                    logger.info(f"Compressed {successful}/{len(items_to_process)} {item_type}...")
            else:
                failed += 1
                logger.error(f"Failed to compress {future_to_item[future]}: {error_msg}")
    
    elapsed_time = time.time() - start_time
    
    # Summary
    logger.info("=" * 50)
    logger.info("COMPRESSION SUMMARY")
    logger.info("=" * 50)
    logger.info(f"Total {item_type} processed: {len(items_to_process)}")
    logger.info(f"Successfully compressed: {successful}")
    logger.info(f"Failed: {failed}")
    logger.info(f"Total time: {elapsed_time:.2f} seconds")
    logger.info(f"Average time per item: {elapsed_time/len(items_to_process):.3f} seconds")
    logger.info(f"Items per second: {len(items_to_process)/elapsed_time:.2f}")
    
    # Calculate compression statistics
    if all_stats:
        total_original_size = sum(stats.get('original_size', 0) for stats in all_stats)
        total_compressed_size = sum(stats.get('compressed_size', 0) for stats in all_stats)
        total_files = sum(stats.get('files_count', 0) for stats in all_stats)
        
        if total_original_size > 0:
            compression_ratio = (1 - total_compressed_size / total_original_size) * 100
            logger.info(f"Total files compressed: {total_files}")
            logger.info(f"Original size: {total_original_size / (1024**3):.2f} GB")
            logger.info(f"Compressed size: {total_compressed_size / (1024**3):.2f} GB")
            logger.info(f"Space savings: {compression_ratio:.1f}%")
    
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    exit(main())
