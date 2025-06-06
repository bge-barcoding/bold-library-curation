#!/usr/bin/env python3
"""
Phylogenetic Analysis Pipeline for BOLD Database
================================================

This script generates phylogenies for each family using OTU representatives,
with outgroups selected from other families within the same order.

Requirements:
- biopython
- pandas
- muscle (alignment)
- iqtree (phylogeny reconstruction)
"""

import sqlite3
import pandas as pd
import os
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import random

class PhylogeneticPipeline:
    def __init__(self, db_path, output_dir="phylogenies"):
        """
        Initialize the phylogenetic pipeline.
        
        Args:
            db_path (str): Path to the BOLD SQLite database
            output_dir (str): Output directory for results
        """
        self.db_path = db_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        self.fasta_dir = self.output_dir / "fasta_files"
        self.alignment_dir = self.output_dir / "alignments"
        self.tree_dir = self.output_dir / "trees"
        
        for dir_path in [self.fasta_dir, self.alignment_dir, self.tree_dir]:
            dir_path.mkdir(exist_ok=True)
    
    def get_otu_representatives(self):
        """
        Get representative sequences for each OTU from the database.
        
        Returns:
            dict: Family -> list of (otu_id, sequence, record_info) tuples
        """
        conn = sqlite3.connect(self.db_path)
        
        # Get OTU representatives with sequences
        query = """
        SELECT DISTINCT
            b.family,
            b.`order`,
            o.otu_id,
            b.nuc,
            b.genus,
            b.species,
            b.processid,
            b.recordid
        FROM bold b
        JOIN bold_otus o ON b.recordid = o.recordid
        WHERE b.family IS NOT NULL 
            AND b.`order` IS NOT NULL
            AND b.nuc IS NOT NULL 
            AND b.nuc != ''
            AND LENGTH(b.nuc) > 200  -- Filter very short sequences
        ORDER BY b.family, o.otu_id
        """
        
        df = pd.read_sql_query(query, conn)
        conn.close()
        
        # Group by family and select one representative per OTU
        family_otus = defaultdict(list)
        
        for family, group in df.groupby('family'):
            # For each OTU in the family, select the longest sequence as representative
            for otu_id, otu_group in group.groupby('otu_id'):
                # Select longest sequence as representative
                rep = otu_group.loc[otu_group['nuc'].str.len().idxmax()]
                
                record_info = {
                    'otu_id': rep['otu_id'],
                    'genus': rep['genus'],
                    'species': rep['species'],
                    'processid': rep['processid'],
                    'recordid': rep['recordid'],
                    'order': rep['order']
                }
                
                family_otus[family].append((rep['otu_id'], rep['nuc'], record_info))
        
        return family_otus
    
    def select_outgroups(self, family_otus, target_family, num_outgroups=3):
        """
        Select outgroup sequences from other families in the same order.
        
        Args:
            family_otus (dict): All family OTU data
            target_family (str): Family to find outgroups for
            num_outgroups (int): Number of outgroup sequences to select
            
        Returns:
            list: List of (otu_id, sequence, record_info) tuples for outgroups
        """
        if target_family not in family_otus:
            return []
        
        # Get the order of the target family
        target_order = family_otus[target_family][0][2]['order']
        
        # Find other families in the same order
        other_families = [
            family for family in family_otus.keys() 
            if family != target_family and 
            family_otus[family][0][2]['order'] == target_order
        ]
        
        outgroups = []
        for other_family in other_families:
            if len(outgroups) >= num_outgroups:
                break
            
            # Randomly select one OTU from this family as outgroup
            if family_otus[other_family]:
                outgroup = random.choice(family_otus[other_family])
                outgroups.append(outgroup)
        
        return outgroups
    
    def create_fasta_file(self, family, ingroup_sequences, outgroup_sequences):
        """
        Create a FASTA file for a family with ingroup and outgroup sequences.
        
        Args:
            family (str): Family name
            ingroup_sequences (list): List of (otu_id, sequence, record_info) tuples
            outgroup_sequences (list): List of (otu_id, sequence, record_info) tuples
            
        Returns:
            Path: Path to the created FASTA file
        """
        fasta_file = self.fasta_dir / f"{family}.fasta"
        
        records = []
        
        # Add ingroup sequences
        for otu_id, sequence, info in ingroup_sequences:
            seq_id = f"{family}_{otu_id}_{info['genus']}_{info['species']}_{info['processid']}"
            record = SeqRecord(
                Seq(sequence),
                id=seq_id,
                description=f"Ingroup | {info['genus']} {info['species']}"
            )
            records.append(record)
        
        # Add outgroup sequences
        for otu_id, sequence, info in outgroup_sequences:
            seq_id = f"OUTGROUP_{otu_id}_{info['genus']}_{info['species']}_{info['processid']}"
            record = SeqRecord(
                Seq(sequence),
                id=seq_id,
                description=f"Outgroup | {info['genus']} {info['species']}"
            )
            records.append(record)
        
        SeqIO.write(records, fasta_file, "fasta")
        return fasta_file
    
    def align_sequences(self, fasta_file):
        """
        Align sequences using MUSCLE.
        
        Args:
            fasta_file (Path): Path to input FASTA file
            
        Returns:
            Path: Path to alignment file
        """
        family = fasta_file.stem
        alignment_file = self.alignment_dir / f"{family}_aligned.fasta"
        
        cmd = [
            "muscle",
            "-in", str(fasta_file),
            "-out", str(alignment_file)
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"✓ Aligned sequences for {family}")
            return alignment_file
        except subprocess.CalledProcessError as e:
            print(f"✗ Error aligning {family}: {e}")
            return None
    
    def build_tree(self, alignment_file, bootstrap=1000):
        """
        Build phylogenetic tree using IQ-TREE.
        
        Args:
            alignment_file (Path): Path to alignment file
            bootstrap (int): Number of bootstrap replicates
            
        Returns:
            Path: Path to tree file
        """
        family = alignment_file.stem.replace("_aligned", "")
        tree_prefix = self.tree_dir / family
        
        cmd = [
            "iqtree",
            "-s", str(alignment_file),
            "-pre", str(tree_prefix),
            "-bb", str(bootstrap),
            "-nt", "AUTO"
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"✓ Built tree for {family}")
            return tree_prefix.with_suffix(".treefile")
        except subprocess.CalledProcessError as e:
            print(f"✗ Error building tree for {family}: {e}")
            return None
    
    def run_pipeline(self, min_otus=3, max_families=None):
        """
        Run the complete phylogenetic analysis pipeline.
        
        Args:
            min_otus (int): Minimum number of OTUs required for analysis
            max_families (int): Maximum number of families to process (for testing)
        """
        print("Starting phylogenetic analysis pipeline...")
        print(f"Database: {self.db_path}")
        print(f"Output directory: {self.output_dir}")
        print("="*60)
        
        # Get OTU representatives
        print("1. Extracting OTU representatives...")
        family_otus = self.get_otu_representatives()
        
        print(f"Found {len(family_otus)} families:")
        for family, otus in family_otus.items():
            print(f"  - {family}: {len(otus)} OTUs")
        
        # Filter families with sufficient OTUs
        families_to_analyze = {
            family: otus for family, otus in family_otus.items() 
            if len(otus) >= min_otus
        }
        
        if max_families:
            families_to_analyze = dict(list(families_to_analyze.items())[:max_families])
        
        print(f"\nAnalyzing {len(families_to_analyze)} families with ≥{min_otus} OTUs")
        print("="*60)
        
        results = {}
        
        for i, (family, otus) in enumerate(families_to_analyze.items(), 1):
            print(f"\n{i}. Processing {family} ({len(otus)} OTUs)")
            
            # Select outgroups
            outgroups = self.select_outgroups(family_otus, family)
            print(f"   Selected {len(outgroups)} outgroups")
            
            if len(outgroups) == 0:
                print(f"   ⚠ No suitable outgroups found for {family}")
                continue
            
            # Create FASTA file
            fasta_file = self.create_fasta_file(family, otus, outgroups)
            print(f"   Created FASTA: {fasta_file.name}")
            
            # Align sequences
            alignment_file = self.align_sequences(fasta_file)
            if alignment_file is None:
                continue
            
            # Build tree
            tree_file = self.build_tree(alignment_file)
            if tree_file is None:
                continue
            
            results[family] = {
                'fasta': fasta_file,
                'alignment': alignment_file,
                'tree': tree_file,
                'num_ingroups': len(otus),
                'num_outgroups': len(outgroups)
            }
        
        print("\n" + "="*60)
        print("PIPELINE COMPLETE")
        print("="*60)
        print(f"Successfully processed {len(results)} families:")
        
        for family, result in results.items():
            print(f"  - {family}: {result['num_ingroups']} ingroups, "
                  f"{result['num_outgroups']} outgroups")
        
        return results
    
    def generate_summary_report(self, results):
        """
        Generate a summary report of the analysis.
        
        Args:
            results (dict): Results from run_pipeline()
        """
        report_file = self.output_dir / "analysis_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write("Phylogenetic Analysis Summary\n")
            f.write("="*50 + "\n\n")
            f.write(f"Database: {self.db_path}\n")
            f.write(f"Output directory: {self.output_dir}\n")
            f.write(f"Total families analyzed: {len(results)}\n\n")
            
            for family, result in results.items():
                f.write(f"Family: {family}\n")
                f.write(f"  Ingroup OTUs: {result['num_ingroups']}\n")
                f.write(f"  Outgroups: {result['num_outgroups']}\n")
                f.write(f"  Tree file: {result['tree'].name}\n\n")
        
        print(f"Summary report saved to: {report_file}")


def main():
    """Main function to run the pipeline."""
    
    # Configuration
    DB_PATH = r"C:\GitHub\bold-library-curation\results\__test_6\results\bold.db"
    OUTPUT_DIR = "phylogenies"
    MIN_OTUS = 3  # Minimum OTUs per family for analysis
    MAX_FAMILIES = None  # Set to a number for testing, None for all families
    
    # Check if required external programs are available
    required_programs = ['muscle', 'iqtree']
    missing_programs = []
    
    for program in required_programs:
        try:
            subprocess.run([program, '--help'], capture_output=True)
        except FileNotFoundError:
            missing_programs.append(program)
    
    if missing_programs:
        print(f"Missing required programs: {', '.join(missing_programs)}")
        print("\nInstallation instructions:")
        print("- MUSCLE: https://www.drive5.com/muscle/downloads.htm")
        print("- IQ-TREE: http://www.iqtree.org/")
        return
    
    # Run pipeline
    pipeline = PhylogeneticPipeline(DB_PATH, OUTPUT_DIR)
    results = pipeline.run_pipeline(min_otus=MIN_OTUS, max_families=MAX_FAMILIES)
    
    if results:
        pipeline.generate_summary_report(results)
    else:
        print("No families were successfully processed.")


if __name__ == "__main__":
    main()
