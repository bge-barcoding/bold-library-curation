#!/usr/bin/env python3
"""
Alternative Phylogenetic Pipeline using MAFFT + FastTree
=======================================================

This version uses more commonly available tools:
- MAFFT for sequence alignment
- FastTree for phylogeny reconstruction
- Optional: RAxML-NG as alternative to FastTree

Supports both local installation and conda environments.
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
import json

class SimplePhylogeneticPipeline:
    def __init__(self, db_path, output_dir="phylogenies_simple"):
        self.db_path = db_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        self.fasta_dir = self.output_dir / "sequences"
        self.alignment_dir = self.output_dir / "alignments"
        self.tree_dir = self.output_dir / "trees"
        self.log_dir = self.output_dir / "logs"
        
        for dir_path in [self.fasta_dir, self.alignment_dir, self.tree_dir, self.log_dir]:
            dir_path.mkdir(exist_ok=True)
        
        # Initialize metadata storage
        self.metadata = {}
    
    def check_dependencies(self):
        """Check for available phylogenetic software."""
        tools = {}
        
        # Check for alignment tools
        for aligner in ['mafft', 'muscle', 'clustalw']:
            try:
                result = subprocess.run([aligner, '--help'], 
                                      capture_output=True, text=True)
                tools[aligner] = True
                print(f"✓ Found {aligner}")
            except FileNotFoundError:
                tools[aligner] = False
        
        # Check for tree building tools
        for tree_builder in ['fasttree', 'FastTree', 'iqtree', 'raxml-ng']:
            try:
                if tree_builder in ['fasttree', 'FastTree']:
                    result = subprocess.run([tree_builder], 
                                          capture_output=True, text=True)
                else:
                    result = subprocess.run([tree_builder, '--help'], 
                                          capture_output=True, text=True)
                tools[tree_builder] = True
                print(f"✓ Found {tree_builder}")
            except FileNotFoundError:
                tools[tree_builder] = False
        
        return tools
    
    def get_family_statistics(self):
        """Get comprehensive statistics about families and OTUs."""
        conn = sqlite3.connect(self.db_path)
        
        query = """
        SELECT 
            b.family,
            b.`order`,
            COUNT(DISTINCT o.otu_id) as num_otus,
            COUNT(*) as num_records,
            MIN(LENGTH(b.nuc)) as min_seq_length,
            MAX(LENGTH(b.nuc)) as max_seq_length,
            AVG(LENGTH(b.nuc)) as avg_seq_length,
            COUNT(DISTINCT b.genus) as num_genera,
            COUNT(DISTINCT b.species) as num_species
        FROM bold b
        JOIN bold_otus o ON b.recordid = o.recordid
        WHERE b.family IS NOT NULL 
            AND b.`order` IS NOT NULL
            AND b.nuc IS NOT NULL 
            AND b.nuc != ''
            AND LENGTH(b.nuc) > 200
        GROUP BY b.family, b.`order`
        ORDER BY num_otus DESC
        """
        
        df = pd.read_sql_query(query, conn)
        conn.close()
        
        return df
    
    def extract_otu_representatives(self, strategy='longest'):
        """
        Extract representative sequences for each OTU.
        
        Args:
            strategy (str): Strategy for selecting representatives
                          'longest' - longest sequence per OTU
                          'random' - random sequence per OTU
                          'best_quality' - highest quality sequence per OTU
        """
        conn = sqlite3.connect(self.db_path)
        
        query = """
        SELECT 
            b.family,
            b.`order`,
            o.otu_id,
            b.nuc,
            b.genus,
            b.species,
            b.processid,
            b.recordid,
            LENGTH(b.nuc) as seq_length,
            CASE 
                WHEN b.nuc LIKE '%N%' THEN 1 
                ELSE 0 
            END as has_ambiguous
        FROM bold b
        JOIN bold_otus o ON b.recordid = o.recordid
        WHERE b.family IS NOT NULL 
            AND b.`order` IS NOT NULL
            AND b.nuc IS NOT NULL 
            AND b.nuc != ''
            AND LENGTH(b.nuc) > 200
        ORDER BY b.family, o.otu_id, seq_length DESC
        """
        
        df = pd.read_sql_query(query, conn)
        conn.close()
        
        family_data = defaultdict(lambda: defaultdict(list))
        
        # Group sequences by family and OTU
        for _, row in df.iterrows():
            family_data[row['family']][row['otu_id']].append(row.to_dict())
        
        # Select representatives based on strategy
        representatives = defaultdict(list)
        
        for family, otus in family_data.items():
            for otu_id, sequences in otus.items():
                if strategy == 'longest':
                    rep = max(sequences, key=lambda x: x['seq_length'])
                elif strategy == 'random':
                    rep = random.choice(sequences)
                elif strategy == 'best_quality':
                    # Prefer sequences without ambiguous bases, then longest
                    no_ambiguous = [s for s in sequences if s['has_ambiguous'] == 0]
                    if no_ambiguous:
                        rep = max(no_ambiguous, key=lambda x: x['seq_length'])
                    else:
                        rep = max(sequences, key=lambda x: x['seq_length'])
                
                representatives[family].append(rep)
        
        return representatives
    
    def smart_outgroup_selection(self, representatives, target_family, 
                               num_outgroups=3, selection_strategy='diverse'):
        """
        Intelligent outgroup selection with multiple strategies.
        
        Args:
            selection_strategy (str): 'diverse' - maximize taxonomic diversity
                                    'sister' - select closest sister families
                                    'random' - random selection
        """
        if target_family not in representatives:
            return []
        
        # Get target order
        target_order = representatives[target_family][0]['order']
        
        # Find candidate outgroup families in same order
        candidate_families = [
            family for family in representatives.keys()
            if family != target_family and 
            representatives[family][0]['order'] == target_order
        ]
        
        if not candidate_families:
            return []
        
        outgroups = []
        
        if selection_strategy == 'diverse':
            # Select from different families, prioritizing larger families
            family_sizes = {
                family: len(representatives[family]) 
                for family in candidate_families
            }
            sorted_families = sorted(family_sizes.items(), 
                                   key=lambda x: x[1], reverse=True)
            
            for family, _ in sorted_families[:num_outgroups]:
                # Select one representative from this family
                outgroup = random.choice(representatives[family])
                outgroups.append(outgroup)
        
        elif selection_strategy == 'random':
            # Randomly select from all available outgroup sequences
            all_outgroup_seqs = []
            for family in candidate_families:
                all_outgroup_seqs.extend(representatives[family])
            
            outgroups = random.sample(all_outgroup_seqs, 
                                    min(num_outgroups, len(all_outgroup_seqs)))
        
        return outgroups
    
    def create_alignment_file(self, family, ingroup_seqs, outgroup_seqs):
        """Create FASTA file and perform alignment."""
        # Create FASTA file
        fasta_file = self.fasta_dir / f"{family}.fasta"
        
        records = []
        
        # Add ingroup sequences
        for seq_data in ingroup_seqs:
            seq_id = f"{family}_{seq_data['otu_id']}_{seq_data['genus']}_{seq_data['species']}"
            # Clean sequence ID to avoid alignment software issues
            seq_id = seq_id.replace(' ', '_').replace('(', '').replace(')', '')
            
            record = SeqRecord(
                Seq(seq_data['nuc']),
                id=seq_id,
                description=f"Ingroup|{seq_data['genus']}_{seq_data['species']}"
            )
            records.append(record)
        
        # Add outgroup sequences
        for seq_data in outgroup_seqs:
            seq_id = f"OUT_{seq_data['family']}_{seq_data['otu_id']}_{seq_data['genus']}_{seq_data['species']}"
            seq_id = seq_id.replace(' ', '_').replace('(', '').replace(')', '')
            
            record = SeqRecord(
                Seq(seq_data['nuc']),
                id=seq_id,
                description=f"Outgroup|{seq_data['genus']}_{seq_data['species']}"
            )
            records.append(record)
        
        SeqIO.write(records, fasta_file, "fasta")
        
        # Store metadata
        self.metadata[family] = {
            'num_ingroups': len(ingroup_seqs),
            'num_outgroups': len(outgroup_seqs),
            'ingroup_otus': [s['otu_id'] for s in ingroup_seqs],
            'outgroup_families': list(set(s['family'] for s in outgroup_seqs))
        }
        
        return fasta_file
    
    def align_with_mafft(self, fasta_file, algorithm='auto'):
        """
        Perform sequence alignment using MAFFT.
        
        Args:
            algorithm (str): MAFFT algorithm ('auto', 'linsi', 'ginsi', 'einsi')
        """
        family = fasta_file.stem
        alignment_file = self.alignment_dir / f"{family}_aligned.fasta"
        log_file = self.log_dir / f"{family}_mafft.log"
        
        if algorithm == 'auto':
            cmd = ['mafft', '--auto', str(fasta_file)]
        elif algorithm == 'linsi':
            cmd = ['mafft', '--localpair', '--maxiterate', '1000', str(fasta_file)]
        elif algorithm == 'ginsi':
            cmd = ['mafft', '--globalpair', '--maxiterate', '1000', str(fasta_file)]
        elif algorithm == 'einsi':
            cmd = ['mafft', '--genafpair', '--maxiterate', '1000', str(fasta_file)]
        
        try:
            with open(alignment_file, 'w') as out_f, open(log_file, 'w') as log_f:
                result = subprocess.run(cmd, stdout=out_f, stderr=log_f, 
                                      check=True, text=True)
            
            print(f"✓ Aligned {family} using MAFFT ({algorithm})")
            return alignment_file
        except subprocess.CalledProcessError as e:
            print(f"✗ MAFFT alignment failed for {family}: {e}")
            return None
    
    def build_tree_fasttree(self, alignment_file, model='gtr'):
        """
        Build phylogenetic tree using FastTree.
        
        Args:
            model (str): Substitution model ('gtr', 'jc', 'nt')
        """
        family = alignment_file.stem.replace('_aligned', '')
        tree_file = self.tree_dir / f"{family}_fasttree.newick"
        log_file = self.log_dir / f"{family}_fasttree.log"
        
        # Determine FastTree command
        fasttree_cmd = None
        for cmd in ['fasttree', 'FastTree']:
            try:
                subprocess.run([cmd], capture_output=True)
                fasttree_cmd = cmd
                break
            except FileNotFoundError:
                continue
        
        if not fasttree_cmd:
            print(f"✗ FastTree not found for {family}")
            return None
        
        if model == 'gtr':
            cmd = [fasttree_cmd, '-nt', '-gtr', str(alignment_file)]
        else:
            cmd = [fasttree_cmd, '-nt', str(alignment_file)]
        
        try:
            with open(tree_file, 'w') as out_f, open(log_file, 'w') as log_f:
                result = subprocess.run(cmd, stdout=out_f, stderr=log_f, 
                                      check=True, text=True)
            
            print(f"✓ Built tree for {family} using FastTree")
            return tree_file
        except subprocess.CalledProcessError as e:
            print(f"✗ FastTree failed for {family}: {e}")
            return None
    
    def build_tree_iqtree(self, alignment_file, bootstrap=1000):
        """Build tree using IQ-TREE with model selection."""
        family = alignment_file.stem.replace('_aligned', '')
        tree_prefix = self.tree_dir / family
        
        cmd = [
            'iqtree',
            '-s', str(alignment_file),
            '-pre', str(tree_prefix),
            '-bb', str(bootstrap),
            '-nt', 'AUTO',
            '-m', 'MFP'  # Model finder plus
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"✓ Built tree for {family} using IQ-TREE")
            return tree_prefix.with_suffix('.treefile')
        except subprocess.CalledProcessError as e:
            print(f"✗ IQ-TREE failed for {family}: {e}")
            return None
    
    def run_analysis(self, min_otus=3, max_families=None, 
                    alignment_method='mafft', tree_method='fasttree',
                    num_outgroups=3, selection_strategy='longest'):
        """
        Run the complete phylogenetic analysis.
        
        Args:
            min_otus (int): Minimum OTUs per family
            max_families (int): Limit number of families (for testing)
            alignment_method (str): 'mafft' or 'muscle'
            tree_method (str): 'fasttree' or 'iqtree'
            num_outgroups (int): Number of outgroup sequences
            selection_strategy (str): OTU selection strategy
        """
        print("Simple Phylogenetic Analysis Pipeline")
        print("="*50)
        
        # Check dependencies
        print("Checking software dependencies...")
        tools = self.check_dependencies()
        
        if not any(tools.get(tool, False) for tool in ['mafft', 'muscle']):
            print("⚠ No alignment software found. Install MAFFT or MUSCLE.")
            return {}
        
        if not any(tools.get(tool, False) for tool in ['fasttree', 'FastTree', 'iqtree']):
            print("⚠ No tree building software found. Install FastTree or IQ-TREE.")
            return {}
        
        # Get family statistics
        print("\nAnalyzing database...")
        stats = self.get_family_statistics()
        print(f"Found {len(stats)} families with sequences")
        
        # Extract representatives
        print(f"Extracting OTU representatives using '{selection_strategy}' strategy...")
        representatives = self.extract_otu_representatives(selection_strategy)
        
        # Filter families by minimum OTU count
        families_to_analyze = {
            family: seqs for family, seqs in representatives.items()
            if len(seqs) >= min_otus
        }
        
        if max_families:
            families_to_analyze = dict(list(families_to_analyze.items())[:max_families])
        
        print(f"Analyzing {len(families_to_analyze)} families with ≥{min_otus} OTUs")
        print("="*50)
        
        results = {}
        
        for i, (family, ingroup_seqs) in enumerate(families_to_analyze.items(), 1):
            print(f"\n{i}. Processing {family} ({len(ingroup_seqs)} OTUs)")
            
            # Select outgroups
            outgroup_seqs = self.smart_outgroup_selection(
                representatives, family, num_outgroups, 'diverse'
            )
            
            if not outgroup_seqs:
                print(f"   ⚠ No outgroups available for {family}")
                continue
            
            print(f"   Selected {len(outgroup_seqs)} outgroups from families: "
                  f"{set(s['family'] for s in outgroup_seqs)}")
            
            # Create FASTA file
            fasta_file = self.create_alignment_file(family, ingroup_seqs, outgroup_seqs)
            
            # Perform alignment
            if alignment_method == 'mafft' and tools.get('mafft', False):
                alignment_file = self.align_with_mafft(fasta_file)
            else:
                print(f"   ⚠ {alignment_method} not available, skipping {family}")
                continue
            
            if not alignment_file:
                continue
            
            # Build tree
            if tree_method == 'fasttree' and any(tools.get(t, False) for t in ['fasttree', 'FastTree']):
                tree_file = self.build_tree_fasttree(alignment_file)
            elif tree_method == 'iqtree' and tools.get('iqtree', False):
                tree_file = self.build_tree_iqtree(alignment_file)
            else:
                print(f"   ⚠ {tree_method} not available, skipping {family}")
                continue
            
            if tree_file:
                results[family] = {
                    'fasta': fasta_file,
                    'alignment': alignment_file,
                    'tree': tree_file,
                    'metadata': self.metadata[family]
                }
        
        # Save metadata
        metadata_file = self.output_dir / "analysis_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)
        
        self.generate_final_report(results, stats)
        return results
    
    def generate_final_report(self, results, stats):
        """Generate comprehensive analysis report."""
        report_file = self.output_dir / "phylogenetic_analysis_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("BOLD Database Phylogenetic Analysis Report\n")
            f.write("="*60 + "\n\n")
            
            f.write(f"Database: {self.db_path}\n")
            f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Output directory: {self.output_dir}\n\n")
            
            f.write("SUMMARY STATISTICS\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total families in database: {len(stats)}\n")
            f.write(f"Families successfully analyzed: {len(results)}\n")
            f.write(f"Total OTUs analyzed: {sum(r['metadata']['num_ingroups'] for r in results.values())}\n")
            f.write(f"Total outgroups used: {sum(r['metadata']['num_outgroups'] for r in results.values())}\n\n")
            
            f.write("FAMILY DETAILS\n")
            f.write("-" * 30 + "\n")
            for family, result in results.items():
                meta = result['metadata']
                f.write(f"\n{family}:\n")
                f.write(f"  Ingroup OTUs: {meta['num_ingroups']}\n")
                f.write(f"  Outgroups: {meta['num_outgroups']}\n")
                f.write(f"  Outgroup families: {', '.join(meta['outgroup_families'])}\n")
                f.write(f"  Tree file: {result['tree'].name}\n")
            
            f.write(f"\nALL FILES SAVED TO: {self.output_dir}\n")
            f.write("- sequences/: FASTA files\n")
            f.write("- alignments/: Aligned sequences\n")
            f.write("- trees/: Phylogenetic trees\n")
            f.write("- logs/: Analysis logs\n")
        
        print(f"\n{'='*60}")
        print("ANALYSIS COMPLETE!")
        print(f"{'='*60}")
        print(f"Successfully analyzed {len(results)} families")
        print(f"Report saved to: {report_file}")
        print(f"All outputs in: {self.output_dir}")


def main():
    """Main execution function with configuration options."""
    
    # CONFIGURATION - Modify these parameters as needed
    CONFIG = {
        'db_path': r"C:\GitHub\bold-library-curation\results\__test_6\results\bold.db",
        'output_dir': "phylogenies_simple",
        'min_otus': 3,              # Minimum OTUs per family
        'max_families': 5,          # Set to None for all families
        'alignment_method': 'mafft', # 'mafft' or 'muscle'
        'tree_method': 'fasttree',  # 'fasttree' or 'iqtree'
        'num_outgroups': 3,         # Number of outgroup sequences
        'selection_strategy': 'best_quality'  # 'longest', 'random', 'best_quality'
    }
    
    print("BOLD Database Phylogenetic Analysis")
    print("Configuration:")
    for key, value in CONFIG.items():
        print(f"  {key}: {value}")
    print()
    
    # Initialize and run pipeline
    pipeline = SimplePhylogeneticPipeline(CONFIG['db_path'], CONFIG['output_dir'])
    
    results = pipeline.run_analysis(
        min_otus=CONFIG['min_otus'],
        max_families=CONFIG['max_families'],
        alignment_method=CONFIG['alignment_method'],
        tree_method=CONFIG['tree_method'],
        num_outgroups=CONFIG['num_outgroups'],
        selection_strategy=CONFIG['selection_strategy']
    )
    
    if not results:
        print("No phylogenies were successfully generated.")
        print("Check software dependencies and try again.")
    else:
        print(f"\nGenerated {len(results)} phylogenetic trees!")


if __name__ == "__main__":
    main()

        