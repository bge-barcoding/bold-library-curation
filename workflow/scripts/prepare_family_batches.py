#!/usr/bin/env python3
"""
BOLD Family Batch Preparation Script
Analyzes the database and creates family batches for parallel processing
"""

import sqlite3
import argparse
import json
import os
import sys
from pathlib import Path

def get_family_counts(db_path, threshold=10000):
    """Get family counts and classify them by size"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Query to get family counts with complete taxonomy
    query = """
    SELECT 
        COALESCE(kingdom, 'Unknown') as kingdom,
        COALESCE(phylum, 'Unknown') as phylum,
        COALESCE("class", 'Unknown') as class,
        COALESCE("order", 'Unknown') as "order",
        COALESCE(family, 'Unknown') as family,
        COUNT(*) as record_count,
        COUNT(DISTINCT subfamily) as subfamily_count
    FROM (
        SELECT DISTINCT recordid, kingdom, phylum, "class", "order", family, subfamily
        FROM bold
        WHERE family IS NOT NULL AND family != ''
    )
    GROUP BY kingdom, phylum, "class", "order", family
    HAVING record_count > 0
    ORDER BY record_count DESC
    """
    
    cursor.execute(query)
    families = []
    
    for row in cursor.fetchall():
        kingdom, phylum, class_name, order, family, count, subfamily_count = row
        
        family_info = {
            'kingdom': kingdom,
            'phylum': phylum,
            'class': class_name,
            'order': order,
            'family': family,
            'record_count': count,
            'subfamily_count': subfamily_count,
            'needs_splitting': count > threshold and subfamily_count > 1
        }
        families.append(family_info)
    
    conn.close()
    return families

def create_batches(families, num_jobs):
    """Create balanced batches for job array processing"""
    if not families:
        return []
    
    # Don't create more batches than families
    actual_num_jobs = min(num_jobs, len(families))
    
    # Simple round-robin distribution
    batches = [[] for _ in range(actual_num_jobs)]
    
    # Sort families by record count (largest first) for better load balancing
    families.sort(key=lambda x: x['record_count'], reverse=True)
    
    # Distribute families round-robin style
    for i, family in enumerate(families):
        batch_idx = i % actual_num_jobs
        batches[batch_idx].append(family)
    
    # All batches should be non-empty now, but filter just in case
    return [batch for batch in batches if batch]

def main():
    parser = argparse.ArgumentParser(description='Prepare family batches for parallel processing')
    parser.add_argument('database', help='Path to BOLD database')
    parser.add_argument('--output-dir', required=True, help='Output directory for batch files')
    parser.add_argument('--num-jobs', type=int, default=64, help='Number of job array tasks')
    parser.add_argument('--threshold', type=int, default=10000, help='Family size threshold for splitting')
    parser.add_argument('--result-tsv', help='Alternative: use result TSV instead of database')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Analyzing families from: {args.database}")
    print(f"Family size threshold: {args.threshold}")
    print(f"Target job array size: {args.num_jobs}")
    
    # Get family information
    try:
        families = get_family_counts(args.database, args.threshold)
        print(f"Found {len(families)} families to process")
        
        if not families:
            print("No families found - creating empty batch file")
            empty_batch_file = Path(args.output_dir) / "batch_0.json"
            with open(empty_batch_file, 'w') as f:
                json.dump([], f, indent=2)
            return
        
        # Create batches
        batches = create_batches(families, args.num_jobs)
        actual_jobs = len(batches)
        
        print(f"Created {actual_jobs} batches")
        
        # Write batch files
        for i, batch in enumerate(batches):
            batch_file = Path(args.output_dir) / f"batch_{i}.json"
            with open(batch_file, 'w') as f:
                json.dump(batch, f, indent=2)
            
            total_records = sum(fam['record_count'] for fam in batch)
            print(f"Batch {i}: {len(batch)} families, {total_records} total records")
        
        # Write summary file
        summary = {
            'total_families': len(families),
            'total_batches': actual_jobs,
            'threshold': args.threshold,
            'large_families': len([f for f in families if f['needs_splitting']]),
            'total_records': sum(f['record_count'] for f in families)
        }
        
        summary_file = Path(args.output_dir) / "batch_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"Summary written to: {summary_file}")
        print(f"Large families requiring splitting: {summary['large_families']}")
        print(f"Total records to process: {summary['total_records']}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
