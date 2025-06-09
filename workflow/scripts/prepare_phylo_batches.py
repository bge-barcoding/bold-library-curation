#!/usr/bin/env python3
"""
Prepare Family Batches for Parallel Phylogenetic Analysis
=========================================================

This script queries the BOLD database to identify families suitable for phylogenetic
analysis and splits them into batches for parallel processing via SLURM job arrays.

Usage:
    python prepare_phylo_batches.py DATABASE_PATH --output-dir BATCH_DIR [options]
"""

import sqlite3
import pandas as pd
import json
import argparse
import sys
from pathlib import Path
from datetime import datetime

def get_families_for_phylo(db_path, min_otus=3, export_kingdoms=None):
    """
    Get families suitable for phylogenetic analysis.
    
    Args:
        db_path (str): Path to BOLD SQLite database
        min_otus (int): Minimum number of species-BIN combinations required
        export_kingdoms (list): List of kingdoms to include, or ['all'] for all kingdoms
        
    Returns:
        pandas.DataFrame: Families with their statistics
    """
    conn = sqlite3.connect(db_path)
    
    # Build kingdom filter and parameters in the correct order
    kingdom_filter = ""
    params = []
    
    if export_kingdoms and 'all' not in export_kingdoms:
        placeholders = ','.join(['?' for _ in export_kingdoms])
        kingdom_filter = f"AND b.kingdom IN ({placeholders})"
        params.extend(export_kingdoms)
    
    # Add min_otus parameter at the end for HAVING clause
    params.append(min_otus)
    
    query = f"""
    SELECT 
        b.family,
        b.`order`,
        b.kingdom,
        COUNT(DISTINCT b.species || '|' || b.bin_uri) as species_bin_count,
        COUNT(DISTINCT b.species) as species_count,
        COUNT(*) as total_records,
        MIN(LENGTH(b.nuc)) as min_length,
        MAX(LENGTH(b.nuc)) as max_length,
        AVG(LENGTH(b.nuc)) as avg_length
    FROM bold b
    WHERE b.family IS NOT NULL 
        AND b.`order` IS NOT NULL
        AND b.nuc IS NOT NULL 
        AND LENGTH(b.nuc) > 200
        AND b.genus IS NOT NULL
        AND b.species IS NOT NULL
        AND b.bin_uri IS NOT NULL
        AND b.bin_uri != ''
        AND b.bin_uri != 'None'
        AND b.bin_uri LIKE 'BOLD:%'
        {kingdom_filter}
    GROUP BY b.family, b.`order`, b.kingdom
    HAVING species_bin_count >= ?
    ORDER BY species_bin_count DESC
    """
    
    families_df = pd.read_sql_query(query, conn, params=params)
    conn.close()
    
    return families_df

def create_family_batches(families_df, num_jobs, families_per_job=None):
    """
    Split families into batches for parallel processing.
    
    Args:
        families_df (pandas.DataFrame): Families to process
        num_jobs (int): Number of job array slots
        families_per_job (int): Fixed number of families per job (optional)
        
    Returns:
        list: List of batches, each containing family information
    """
    total_families = len(families_df)
    
    if families_per_job:
        # Use fixed families per job
        actual_jobs = min(num_jobs, (total_families + families_per_job - 1) // families_per_job)
        families_per_batch = families_per_job
    else:
        # Distribute evenly across available jobs
        actual_jobs = min(num_jobs, total_families)
        families_per_batch = (total_families + actual_jobs - 1) // actual_jobs
    
    batches = []
    for i in range(actual_jobs):
        start_idx = i * families_per_batch
        end_idx = min((i + 1) * families_per_batch, total_families)
        
        if start_idx < total_families:
            batch_families = families_df.iloc[start_idx:end_idx].copy()
            
            batch_info = {
                'batch_id': i + 1,  # 1-based for SLURM array
                'families': batch_families[['family', 'order', 'kingdom']].astype(str).to_dict('records'),
                'family_count': len(batch_families),
                'total_species_bins': int(batch_families['species_bin_count'].sum()),
                'complexity_score': float(batch_families['species_bin_count'].mean())
            }
            batches.append(batch_info)
    
    return batches

def save_batches(batches, output_dir):
    """
    Save batch information to JSON files.
    
    Args:
        batches (list): List of batch information
        output_dir (Path): Output directory for batch files
        
    Returns:
        dict: Summary information
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save individual batch files
    for batch in batches:
        batch_file = output_dir / f"phylo_batch_{batch['batch_id']:03d}.json"
        with open(batch_file, 'w') as f:
            json.dump(batch, f, indent=2)
    
    # Create summary
    summary = {
        'created_at': datetime.now().isoformat(),
        'total_batches': len(batches),
        'total_families': sum(b['family_count'] for b in batches),
        'families_per_batch_range': [
            int(min(b['family_count'] for b in batches)),
            int(max(b['family_count'] for b in batches))
        ],
        'total_species_bins': int(sum(b['total_species_bins'] for b in batches)),
        'batch_files': [f"phylo_batch_{b['batch_id']:03d}.json" for b in batches],
        'slurm_array_range': f"1-{len(batches)}"
    }
    
    # Save summary
    summary_file = output_dir / "phylo_batch_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    return summary

def print_analysis_summary(families_df, batches):
    """Print summary of families and batching strategy."""
    print("="*70)
    print("PHYLOGENETIC ANALYSIS BATCH PREPARATION")
    print("="*70)
    
    print(f"Total families identified: {len(families_df)}")
    print(f"Total species-BIN combinations: {families_df['species_bin_count'].sum():,}")
    print(f"Average species-BIN combinations per family: {families_df['species_bin_count'].mean():.1f}")
    print()
    
    # Kingdom distribution
    if 'kingdom' in families_df.columns:
        print("Kingdom distribution:")
        kingdom_counts = families_df.groupby('kingdom').agg({
            'family': 'count',
            'species_bin_count': 'sum'
        }).round(1)
        for kingdom, row in kingdom_counts.iterrows():
            print(f"  {kingdom}: {row['family']} families, {row['species_bin_count']:,} species-BIN combinations")
        print()
    
    # Order distribution (top 10)
    print("Top 10 orders by family count:")
    order_counts = families_df.groupby('order').size().sort_values(ascending=False).head(10)
    for order, count in order_counts.items():
        print(f"  {order}: {count} families")
    print()
    
    # Complexity distribution
    print("Family complexity distribution (species-BIN combinations):")
    bins = [0, 5, 10, 25, 50, 100, 500, float('inf')]
    labels = ['3-4', '5-9', '10-24', '25-49', '50-99', '100-499', '500+']
    complexity_dist = pd.cut(families_df['species_bin_count'], bins=bins, labels=labels, right=False)
    for label, count in complexity_dist.value_counts().sort_index().items():
        print(f"  {label} species-BIN combinations: {count} families")
    print()
    
    print("BATCH CONFIGURATION:")
    print(f"Number of batches: {len(batches)}")
    print(f"Families per batch: {min(b['family_count'] for b in batches)} - {max(b['family_count'] for b in batches)}")
    print(f"SLURM array range: 1-{len(batches)}")
    print()
    
    # Estimated runtime
    avg_time_per_family = 5  # minutes (conservative estimate)
    max_batch_time = max(b['family_count'] for b in batches) * avg_time_per_family
    total_sequential_time = sum(b['family_count'] for b in batches) * avg_time_per_family
    
    print("ESTIMATED RUNTIME:")
    print(f"Sequential processing: {total_sequential_time/60:.1f} hours")
    print(f"Parallel processing: {max_batch_time/60:.1f} hours (max batch)")
    print(f"Speed improvement: ~{total_sequential_time/max_batch_time:.1f}x")
    print("="*70)

def main():
    parser = argparse.ArgumentParser(description="Prepare family batches for parallel phylogenetic analysis")
    parser.add_argument("database", help="Path to BOLD SQLite database")
    parser.add_argument("--output-dir", required=True, help="Output directory for batch files")
    parser.add_argument("--num-jobs", type=int, default=50, help="Number of SLURM array jobs (default: 50)")
    parser.add_argument("--families-per-job", type=int, help="Fixed number of families per job (overrides automatic distribution)")
    parser.add_argument("--min-otus", type=int, default=3, help="Minimum species-BIN combinations required (default: 3)")
    parser.add_argument("--export-kingdoms", nargs='+', default=['all'], 
                       help="Kingdoms to include (default: all)")
    parser.add_argument("--dry-run", action="store_true", help="Show analysis without creating files")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.database).exists():
        print(f"Error: Database file not found: {args.database}")
        return 1
    
    if args.families_per_job and args.families_per_job < 1:
        print("Error: families-per-job must be at least 1")
        return 1
    
    print(f"Analyzing families in database: {args.database}")
    print(f"Minimum species-BIN combinations: {args.min_otus}")
    print(f"Export kingdoms: {args.export_kingdoms}")
    print()
    
    # Get families suitable for analysis
    try:
        families_df = get_families_for_phylo(
            args.database, 
            min_otus=args.min_otus,
            export_kingdoms=args.export_kingdoms if args.export_kingdoms != ['all'] else None
        )
    except Exception as e:
        print(f"Error querying database: {e}")
        return 1
    
    if len(families_df) == 0:
        print("No families found matching criteria")
        return 1
    
    # Create batches
    batches = create_family_batches(families_df, args.num_jobs, args.families_per_job)
    
    # Print analysis
    print_analysis_summary(families_df, batches)
    
    if args.dry_run:
        print("DRY RUN: No files created")
        return 0
    
    # Save batches
    try:
        summary = save_batches(batches, args.output_dir)
        print(f"Batch files created in: {args.output_dir}")
        print(f"Summary file: {Path(args.output_dir) / 'phylo_batch_summary.json'}")
        print(f"Ready for SLURM array: --array={summary['slurm_array_range']}")
        return 0
    except Exception as e:
        print(f"Error creating batch files: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
