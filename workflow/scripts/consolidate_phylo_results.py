#!/usr/bin/env python3
"""
Consolidate Phylogenetic Analysis Results
=========================================

This script consolidates results from parallel phylogenetic analysis jobs,
collecting completion information, generating summary reports, and identifying
any failed or incomplete families.

Usage:
    python consolidate_phylo_results.py BATCH_DIR OUTPUT_DIR [options]
"""

import json
import argparse
import sys
from pathlib import Path
from datetime import datetime
import pandas as pd

def load_batch_summary(batch_dir):
    """Load the original batch summary"""
    summary_file = Path(batch_dir) / "phylo_batch_summary.json"
    if not summary_file.exists():
        raise FileNotFoundError(f"Batch summary not found: {summary_file}")
    
    with open(summary_file, 'r') as f:
        return json.load(f)

def collect_completion_info(batch_dir, log_dir="logs"):
    """Collect completion information from all batch jobs"""
    completion_files = list(Path(log_dir).glob("phylo_checkpoints/completion_*.json"))
    
    completions = []
    for file in completion_files:
        try:
            with open(file, 'r') as f:
                completion = json.load(f)
            completions.append(completion)
        except Exception as e:
            print(f"Warning: Could not read completion file {file}: {e}")
    
    return completions

def collect_checkpoint_info(batch_dir, log_dir="logs"):
    """Collect checkpoint information from all batch jobs"""
    checkpoint_files = list(Path(log_dir).glob("phylo_checkpoints/checkpoint_*.json"))
    
    all_completed_families = set()
    checkpoint_info = {}
    
    for file in checkpoint_files:
        try:
            with open(file, 'r') as f:
                checkpoint = json.load(f)
            
            batch_id = checkpoint.get('batch_id')
            completed_families = checkpoint.get('completed_families', [])
            
            all_completed_families.update(completed_families)
            checkpoint_info[batch_id] = {
                'completed_families': completed_families,
                'last_updated': checkpoint.get('last_updated'),
                'count': len(completed_families)
            }
            
        except Exception as e:
            print(f"Warning: Could not read checkpoint file {file}: {e}")
    
    return all_completed_families, checkpoint_info

def scan_output_directory(output_dir):
    """Scan output directory for generated files"""
    output_path = Path(output_dir)
    
    results = {
        'families_with_trees': set(),
        'families_with_pdfs': set(),
        'families_with_curation_pdfs': set(),
        'tree_files': [],
        'pdf_files': [],
        'curation_pdf_files': [],
        'metadata_files': []
    }
    
    if not output_path.exists():
        return results
    
    # Scan for family directories
    for family_dir in output_path.iterdir():
        if family_dir.is_dir():
            family_name = family_dir.name
            
            # Check for tree files
            tree_file = family_dir / f"{family_name}.treefile"
            if tree_file.exists():
                results['families_with_trees'].add(family_name)
                results['tree_files'].append(str(tree_file))
            
            # Check for PDF files
            pdf_file = family_dir / f"{family_name}_tree.pdf"
            if pdf_file.exists():
                results['families_with_pdfs'].add(family_name)
                results['pdf_files'].append(str(pdf_file))
            
            # Check for curation PDF files
            curation_pdf = family_dir / f"{family_name}_curation_checklist.pdf"
            if curation_pdf.exists():
                results['families_with_curation_pdfs'].add(family_name)
                results['curation_pdf_files'].append(str(curation_pdf))
            
            # Check for metadata files
            metadata_file = family_dir / "tree_metadata.json"
            if metadata_file.exists():
                results['metadata_files'].append(str(metadata_file))
    
    return results

def analyze_bags_grades_distribution(output_dir):
    """Analyze BAGS grades distribution across all families"""
    metadata_files = list(Path(output_dir).glob("*/tree_metadata.json"))
    
    grade_distribution = {}
    monophyly_stats = {
        'total_grade_c_species': 0,
        'monophyletic_count': 0,
        'non_monophyletic_count': 0,
        'untestable_count': 0
    }
    
    family_summaries = []
    
    for metadata_file in metadata_files:
        try:
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            
            family_name = metadata.get('family', 'Unknown')
            bags_analysis = metadata.get('bags_analysis', {})
            monophyly_analysis = metadata.get('monophyly_analysis', {})
            
            # Collect BAGS grade distribution
            unique_species_per_grade = bags_analysis.get('unique_species_per_grade', {})
            for grade, count in unique_species_per_grade.items():
                if grade not in grade_distribution:
                    grade_distribution[grade] = 0
                grade_distribution[grade] += count
            
            # Collect monophyly statistics
            if monophyly_analysis:
                monophyly_stats['total_grade_c_species'] += len(monophyly_analysis.get('grade_c_monophyly', {}))
                monophyly_stats['monophyletic_count'] += monophyly_analysis.get('monophyletic_count', 0)
                monophyly_stats['non_monophyletic_count'] += monophyly_analysis.get('non_monophyletic_count', 0)
                monophyly_stats['untestable_count'] += monophyly_analysis.get('untestable_count', 0)
            
            # Family summary
            family_summaries.append({
                'family': family_name,
                'total_species': bags_analysis.get('total_species', 0),
                'grade_distribution': unique_species_per_grade,
                'monophyletic_species': monophyly_analysis.get('monophyletic_count', 0),
                'non_monophyletic_species': monophyly_analysis.get('non_monophyletic_count', 0)
            })
            
        except Exception as e:
            print(f"Warning: Could not read metadata file {metadata_file}: {e}")
    
    return grade_distribution, monophyly_stats, family_summaries

def generate_consolidation_report(batch_summary, completions, checkpoint_info, output_results, 
                                grade_distribution, monophyly_stats, family_summaries, output_file):
    """Generate comprehensive consolidation report"""
    
    # Calculate statistics
    total_batches_expected = batch_summary['total_batches']
    total_families_expected = batch_summary['total_families']
    batches_completed = len(completions)
    
    families_completed = len(output_results['families_with_trees'])
    families_with_pdfs = len(output_results['families_with_pdfs'])
    families_with_curation = len(output_results['families_with_curation_pdfs'])
    
    completion_rate = (families_completed / total_families_expected * 100) if total_families_expected > 0 else 0
    
    # Calculate total runtime
    total_runtime = sum(c.get('runtime_seconds', 0) for c in completions)
    avg_runtime_per_batch = total_runtime / batches_completed if batches_completed > 0 else 0
    
    # Generate report
    report = {
        'consolidation_info': {
            'generated_at': datetime.now().isoformat(),
            'batch_directory': str(Path(batch_summary.get('batch_files', [''])[0]).parent) if batch_summary.get('batch_files') else '',
            'output_directory': str(output_file.parent)
        },
        'batch_summary': {
            'total_batches_expected': total_batches_expected,
            'batches_completed': batches_completed,
            'batch_completion_rate': (batches_completed / total_batches_expected * 100) if total_batches_expected > 0 else 0,
            'total_families_expected': total_families_expected,
            'families_completed': families_completed,
            'family_completion_rate': completion_rate
        },
        'output_statistics': {
            'families_with_trees': families_completed,
            'families_with_pdfs': families_with_pdfs,
            'families_with_curation_checklists': families_with_curation,
            'total_tree_files': len(output_results['tree_files']),
            'total_pdf_files': len(output_results['pdf_files']),
            'total_curation_pdfs': len(output_results['curation_pdf_files']),
            'total_metadata_files': len(output_results['metadata_files'])
        },
        'runtime_statistics': {
            'total_runtime_seconds': total_runtime,
            'total_runtime_hours': total_runtime / 3600,
            'average_runtime_per_batch_seconds': avg_runtime_per_batch,
            'average_runtime_per_family_seconds': total_runtime / families_completed if families_completed > 0 else 0
        },
        'bags_grade_distribution': grade_distribution,
        'monophyly_analysis': monophyly_stats,
        'top_families_by_species_count': sorted(family_summaries, key=lambda x: x['total_species'], reverse=True)[:20]
    }
    
    # Save report
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    return report

def print_consolidation_summary(report):
    """Print human-readable summary of consolidation results"""
    print("="*70)
    print("PHYLOGENETIC ANALYSIS CONSOLIDATION SUMMARY")
    print("="*70)
    
    batch_summary = report['batch_summary']
    output_stats = report['output_statistics']
    runtime_stats = report['runtime_statistics']
    bags_grades = report['bags_grade_distribution']
    monophyly = report['monophyly_analysis']
    
    print(f"Generated at: {report['consolidation_info']['generated_at']}")
    print()
    
    print("BATCH COMPLETION:")
    print(f"  Batches completed: {batch_summary['batches_completed']}/{batch_summary['total_batches_expected']} "
          f"({batch_summary['batch_completion_rate']:.1f}%)")
    print(f"  Families completed: {batch_summary['families_completed']}/{batch_summary['total_families_expected']} "
          f"({batch_summary['family_completion_rate']:.1f}%)")
    print()
    
    print("OUTPUT FILES GENERATED:")
    print(f"  Phylogenetic trees: {output_stats['families_with_trees']:,}")
    print(f"  Tree visualization PDFs: {output_stats['families_with_pdfs']:,}")
    print(f"  Curation checklist PDFs: {output_stats['families_with_curation_checklists']:,}")
    print(f"  Metadata files: {output_stats['total_metadata_files']:,}")
    print()
    
    print("RUNTIME STATISTICS:")
    print(f"  Total runtime: {runtime_stats['total_runtime_hours']:.1f} hours")
    print(f"  Average per batch: {runtime_stats['average_runtime_per_batch_seconds']:.0f} seconds")
    print(f"  Average per family: {runtime_stats['average_runtime_per_family_seconds']:.1f} seconds")
    print()
    
    print("BAGS GRADE DISTRIBUTION:")
    total_species = sum(bags_grades.values())
    for grade in sorted(bags_grades.keys()):
        count = bags_grades[grade]
        percentage = (count / total_species * 100) if total_species > 0 else 0
        print(f"  Grade {grade}: {count:,} species ({percentage:.1f}%)")
    print(f"  Total species analyzed: {total_species:,}")
    print()
    
    if monophyly['total_grade_c_species'] > 0:
        print("GRADE C MONOPHYLY ANALYSIS:")
        print(f"  Total Grade C species tested: {monophyly['total_grade_c_species']:,}")
        print(f"  Monophyletic: {monophyly['monophyletic_count']:,}")
        print(f"  Non-monophyletic: {monophyly['non_monophyletic_count']:,}")
        print(f"  Insufficient data: {monophyly['untestable_count']:,}")
        
        if monophyly['monophyletic_count'] + monophyly['non_monophyletic_count'] > 0:
            testable = monophyly['monophyletic_count'] + monophyly['non_monophyletic_count']
            mono_rate = (monophyly['monophyletic_count'] / testable * 100)
            print(f"  Monophyly rate: {mono_rate:.1f}%")
        print()
    
    print("TOP 10 FAMILIES BY SPECIES COUNT:")
    top_families = report['top_families_by_species_count'][:10]
    for i, family in enumerate(top_families, 1):
        print(f"  {i:2d}. {family['family']}: {family['total_species']} species")
        if family['non_monophyletic_species'] > 0:
            print(f"      (⚠ {family['non_monophyletic_species']} non-monophyletic Grade C species)")
    
    print("="*70)

def main():
    parser = argparse.ArgumentParser(description="Consolidate parallel phylogenetic analysis results")
    parser.add_argument("batch_dir", help="Directory containing batch files")
    parser.add_argument("output_dir", help="Directory containing phylogenetic analysis outputs")
    parser.add_argument("--output-file", help="Output file for consolidation report", 
                       default="phylogenetic_analysis_summary.json")
    parser.add_argument("--failed-families-file", help="File to save list of failed families",
                       default="failed_families.txt")
    parser.add_argument("--log-dir", help="Log directory (defaults to 'logs')", default="logs")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    batch_dir = Path(args.batch_dir)
    output_dir = Path(args.output_dir)
    
    if not batch_dir.exists():
        print(f"Error: Batch directory not found: {batch_dir}")
        return 1
    
    if not output_dir.exists():
        print(f"Warning: Output directory not found: {output_dir}")
        print("Creating output directory...")
        output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Consolidating phylogenetic analysis results...")
    print(f"Batch directory: {batch_dir}")
    print(f"Output directory: {output_dir}")
    print()
    
    try:
        # Load batch summary
        print("Loading batch summary...")
        batch_summary = load_batch_summary(batch_dir)
        
        # Collect completion information
        print("Collecting completion information...")
        completions = collect_completion_info(batch_dir, args.log_dir)
        
        # Collect checkpoint information
        print("Collecting checkpoint information...")
        completed_families, checkpoint_info = collect_checkpoint_info(batch_dir, args.log_dir)
        
        # Scan output directory
        print("Scanning output directory...")
        output_results = scan_output_directory(output_dir)
        
        # Analyze BAGS grades distribution
        print("Analyzing BAGS grades distribution...")
        grade_distribution, monophyly_stats, family_summaries = analyze_bags_grades_distribution(output_dir)
        
        # Generate consolidation report
        print("Generating consolidation report...")
        output_file = Path(args.output_file)
        report = generate_consolidation_report(
            batch_summary, completions, checkpoint_info, output_results,
            grade_distribution, monophyly_stats, family_summaries, output_file
        )
        
        # Print summary
        print_consolidation_summary(report)
        
        # Identify failed families
        expected_families = set()
        for batch_file in batch_summary.get('batch_files', []):
            batch_path = batch_dir / batch_file
            if batch_path.exists():
                with open(batch_path, 'r') as f:
                    batch_data = json.load(f)
                for family_info in batch_data.get('families', []):
                    expected_families.add(family_info['family'])
        
        completed_families_from_output = output_results['families_with_trees']
        failed_families = expected_families - completed_families_from_output
        
        if failed_families:
            print(f"\nFAILED/INCOMPLETE FAMILIES ({len(failed_families)}):")
            failed_list = sorted(list(failed_families))
            for family in failed_list[:20]:  # Show first 20
                print(f"  - {family}")
            if len(failed_families) > 20:
                print(f"  ... and {len(failed_families) - 20} more")
            
            # Save failed families list
            with open(args.failed_families_file, 'w') as f:
                for family in failed_list:
                    f.write(f"{family}\n")
            print(f"\nComplete list saved to: {args.failed_families_file}")
        else:
            print("\n✓ All families completed successfully!")
        
        print(f"\nConsolidation report saved to: {output_file}")
        return 0
        
    except Exception as e:
        print(f"Error during consolidation: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
