#!/usr/bin/env python3
"""
Consolidate Results from Family Batch Processing
Checks completion status and generates final report
"""

import json
import argparse
import os
import sys
from pathlib import Path
import glob
import time

def collect_batch_results(batch_dir, output_dir):
    """Collect results from all batch processing jobs"""
    results = []
    failed_batches = []
    
    # Find all batch result files
    result_files = glob.glob(str(Path(output_dir) / "batch_*_result.json"))
    
    for result_file in result_files:
        try:
            with open(result_file, 'r') as f:
                result = json.load(f)
            results.append(result)
        except Exception as e:
            print(f"Warning: Could not read result file {result_file}: {e}")
            batch_name = Path(result_file).stem.replace('_result', '')
            failed_batches.append(batch_name)
    
    return results, failed_batches

def check_missing_batches(batch_dir, results):
    """Check for batches that didn't produce results"""
    # Find all batch files
    batch_files = glob.glob(str(Path(batch_dir) / "batch_*.json"))
    expected_batches = set()
    
    for batch_file in batch_files:
        batch_name = Path(batch_file).stem
        expected_batches.add(batch_name)
    
    # Find batches that completed
    completed_batches = set()
    for result in results:
        batch_file = result.get('batch_file', '')
        if batch_file:
            batch_name = Path(batch_file).stem
            completed_batches.add(batch_name)
    
    missing_batches = expected_batches - completed_batches
    return list(missing_batches)

def count_database_files(output_dir):
    """Count created database files"""
    db_files = glob.glob(str(Path(output_dir) / "**/*.db"), recursive=True)
    return len(db_files), db_files

def generate_report(results, failed_batches, missing_batches, output_dir, db_count, db_files):
    """Generate comprehensive processing report"""
    report = {
        'summary': {
            'total_batches_expected': len(results) + len(failed_batches) + len(missing_batches),
            'batches_completed': len(results),
            'batches_failed': len(failed_batches),
            'batches_missing': len(missing_batches),
            'total_families_processed': sum(r.get('processed', 0) for r in results),
            'total_families_failed': sum(r.get('failed', 0) for r in results),
            'total_databases_created': db_count,
            'total_processing_time': sum(r.get('processing_time', 0) for r in results)
        },
        'batch_details': results,
        'failed_batches': failed_batches,
        'missing_batches': missing_batches,
        'database_files': db_files[:100],  # Limit to first 100 for readability
        'generated_at': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    
    # Write detailed report
    report_file = Path(output_dir) / "splitting_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    # Write human-readable summary
    summary_file = Path(output_dir) / "splitting_report.txt"
    with open(summary_file, 'w') as f:
        f.write("BOLD Family Database Splitting Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("Summary:\n")
        f.write(f"  Total batches expected: {report['summary']['total_batches_expected']}\n")
        f.write(f"  Batches completed: {report['summary']['batches_completed']}\n")
        f.write(f"  Batches failed: {report['summary']['batches_failed']}\n")
        f.write(f"  Batches missing: {report['summary']['batches_missing']}\n")
        f.write(f"  Families processed: {report['summary']['total_families_processed']}\n")
        f.write(f"  Families failed: {report['summary']['total_families_failed']}\n")
        f.write(f"  Database files created: {report['summary']['total_databases_created']}\n")
        f.write(f"  Total processing time: {report['summary']['total_processing_time']:.2f} seconds\n\n")
        
        if failed_batches:
            f.write("Failed Batches:\n")
            for batch in failed_batches:
                f.write(f"  - {batch}\n")
            f.write("\n")
        
        if missing_batches:
            f.write("Missing Batches:\n")
            for batch in missing_batches:
                f.write(f"  - {batch}\n")
            f.write("\n")
        
        f.write("Family Statistics:\n")
        family_counts = {}
        for result in results:
            for file_path in result.get('created_files', []):
                # Extract family name from path
                path_parts = Path(file_path).parts
                if len(path_parts) >= 6:  # kingdom/phylum/class/order/family/file.db
                    family = path_parts[-2]
                    family_counts[family] = family_counts.get(family, 0) + 1
        
        for family, count in sorted(family_counts.items())[:20]:  # Top 20
            f.write(f"  {family}: {count} database(s)\n")
        
        if len(family_counts) > 20:
            f.write(f"  ... and {len(family_counts) - 20} more families\n")
    
    return report_file, summary_file

def main():
    parser = argparse.ArgumentParser(description='Consolidate family batch processing results')
    parser.add_argument('batch_dir', help='Directory containing batch files')
    parser.add_argument('output_dir', help='Output directory with results')
    parser.add_argument('--wait-timeout', type=int, default=0, help='Wait timeout in seconds for missing results')
    
    args = parser.parse_args()
    
    print(f"Consolidating results from: {args.output_dir}")
    print(f"Batch directory: {args.batch_dir}")
    
    # Collect results
    results, failed_batches = collect_batch_results(args.batch_dir, args.output_dir)
    missing_batches = check_missing_batches(args.batch_dir, results)
    
    # Wait for missing results if timeout specified
    if args.wait_timeout > 0 and missing_batches:
        print(f"Waiting up to {args.wait_timeout} seconds for missing batches...")
        start_time = time.time()
        
        while missing_batches and (time.time() - start_time) < args.wait_timeout:
            time.sleep(10)  # Check every 10 seconds
            results, failed_batches = collect_batch_results(args.batch_dir, args.output_dir)
            missing_batches = check_missing_batches(args.batch_dir, results)
            
            if missing_batches:
                print(f"Still waiting for {len(missing_batches)} batches...")
    
    # Count database files
    db_count, db_files = count_database_files(args.output_dir)
    
    # Generate report
    report_file, summary_file = generate_report(results, failed_batches, missing_batches, 
                                               args.output_dir, db_count, db_files)
    
    print(f"\nConsolidation complete:")
    print(f"  Processed: {len(results)} batches")
    print(f"  Failed: {len(failed_batches)} batches")
    print(f"  Missing: {len(missing_batches)} batches")
    print(f"  Database files: {db_count}")
    print(f"  Report: {summary_file}")
    
    # Exit with error if there were failures
    if failed_batches or missing_batches:
        print(f"\nWarning: {len(failed_batches + missing_batches)} batches incomplete")
        sys.exit(1)
    else:
        print("\n✓ All batches completed successfully")

if __name__ == "__main__":
    main()
