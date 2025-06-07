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
    
    def get_otu_sequences(self):
        """
        Get sequences for each OTU from the database, including multiple species per OTU.
        
        Returns:
            dict: Family -> list of (otu_id, sequence, record_info) tuples
        """
        conn = sqlite3.connect(self.db_path)
        
        # Get all sequences with OTU information, including bin_uri
        query = """
        SELECT DISTINCT
            b.family,
            b.`order`,
            o.otu_id,
            b.nuc,
            b.genus,
            b.species,
            b.processid,
            b.recordid,
            b.bin_uri
        FROM bold b
        JOIN bold_otus o ON b.recordid = o.recordid
        WHERE b.family IS NOT NULL 
            AND b.`order` IS NOT NULL
            AND b.nuc IS NOT NULL 
            AND b.nuc != ''
            AND LENGTH(b.nuc) > 200  -- Filter very short sequences
            AND b.genus IS NOT NULL
            AND b.species IS NOT NULL
        ORDER BY b.family, o.otu_id, b.genus, b.species
        """
        
        df = pd.read_sql_query(query, conn)
        conn.close()
        
        # Group by family and include multiple species per OTU
        family_sequences = defaultdict(list)
        
        for family, group in df.groupby('family'):
            # For each OTU, include all unique species
            for otu_id, otu_group in group.groupby('otu_id'):
                # Group by genus-species combination to get unique species
                for (genus, species), species_group in otu_group.groupby(['genus', 'species']):
                    # Select the longest sequence for each genus-species combination within the OTU
                    rep = species_group.loc[species_group['nuc'].str.len().idxmax()]
                    
                    record_info = {
                        'otu_id': rep['otu_id'],
                        'genus': rep['genus'],
                        'species': rep['species'],
                        'processid': rep['processid'],
                        'recordid': rep['recordid'],
                        'order': rep['order'],
                        'family': family,  # Add family information for outgroup labeling
                        'bin_uri': rep['bin_uri'] if rep['bin_uri'] else 'NA'
                    }
                    
                    family_sequences[family].append((rep['otu_id'], rep['nuc'], record_info))
        
        return family_sequences
    
    def select_outgroups(self, family_sequences, target_family, num_outgroups=3):
        """
        Select outgroup sequences from other families in the same order.
        
        Args:
            family_sequences (dict): All family sequence data
            target_family (str): Family to find outgroups for
            num_outgroups (int): Number of outgroup sequences to select
            
        Returns:
            list: List of (otu_id, sequence, record_info) tuples for outgroups
        """
        if target_family not in family_sequences:
            return []
        
        # Get the order of the target family
        target_order = family_sequences[target_family][0][2]['order']
        
        # Find other families in the same order
        other_families = [
            family for family in family_sequences.keys() 
            if family != target_family and 
            family_sequences[family][0][2]['order'] == target_order
        ]
        
        outgroups = []
        for other_family in other_families:
            if len(outgroups) >= num_outgroups:
                break
            
            # Randomly select one sequence from this family as outgroup
            if family_sequences[other_family]:
                outgroup = random.choice(family_sequences[other_family])
                outgroups.append(outgroup)
        
        return outgroups
    
    def create_fasta_file(self, family, ingroup_sequences, outgroup_sequences):
        """
        Create a FASTA file for a family with ingroup and outgroup sequences.
        Uses the new branch labeling format: Genus-species-otu-bin_uri-processid
        
        Args:
            family (str): Family name
            ingroup_sequences (list): List of (otu_id, sequence, record_info) tuples
            outgroup_sequences (list): List of (otu_id, sequence, record_info) tuples
            
        Returns:
            Path: Path to the created FASTA file
        """
        fasta_file = self.fasta_dir / f"{family}.fasta"
        
        records = []
        
        # Add ingroup sequences with new labeling format
        for otu_id, sequence, info in ingroup_sequences:
            # Clean values to handle None or empty values
            genus = str(info.get('genus', 'Unknown')).replace(' ', '_')
            species = str(info.get('species', 'sp')).replace(' ', '_')
            otu = str(info.get('otu_id', 'NA')).replace(' ', '_')
            bin_uri = str(info.get('bin_uri', 'NA')).replace(' ', '_')
            processid = str(info.get('processid', 'NA')).replace(' ', '_')
            
            # Create sequence ID in format: Species-otu-bin_uri-processid
            seq_id = f"{species}-{otu}-{bin_uri}-{processid}"
            
            # Clean any problematic characters for phylogenetic software
            seq_id = seq_id.replace('(', '').replace(')', '').replace('[', '').replace(']', '')
            seq_id = seq_id.replace(',', '').replace(';', '').replace(':', '').replace('|', '-')
            
            record = SeqRecord(
                Seq(sequence),
                id=seq_id,
                description=f"Ingroup|{genus}_{species}|{family}"
            )
            records.append(record)
        
        # Add outgroup sequences with new labeling format
        for otu_id, sequence, info in outgroup_sequences:
            # Clean values to handle None or empty values
            genus = str(info.get('genus', 'Unknown')).replace(' ', '_')
            species = str(info.get('species', 'sp')).replace(' ', '_')
            otu = str(info.get('otu_id', 'NA')).replace(' ', '_')
            bin_uri = str(info.get('bin_uri', 'NA')).replace(' ', '_')
            processid = str(info.get('processid', 'NA')).replace(' ', '_')
            
            # Create sequence ID with OUTGROUP prefix
            seq_id = f"OUTGROUP-{species}-{otu}-{bin_uri}-{processid}"
            
            # Clean any problematic characters for phylogenetic software
            seq_id = seq_id.replace('(', '').replace(')', '').replace('[', '').replace(']', '')
            seq_id = seq_id.replace(',', '').replace(';', '').replace(':', '').replace('|', '-')
            
            record = SeqRecord(
                Seq(sequence),
                id=seq_id,
                description=f"Outgroup|{genus}_{species}|{info.get('family', 'Unknown')}"
            )
            records.append(record)
        
        SeqIO.write(records, fasta_file, "fasta")
        return fasta_file
    
    def align_sequences(self, fasta_file, method="mafft"):
        """
        Align sequences using MAFFT or MUSCLE.
        MAFFT is optimized for protein-coding genes with:
        - Automatic reverse complement detection (--adjustdirection)
        - Lower gap penalties suitable for coding sequences
        - L-INS-i algorithm for high accuracy alignment
        
        Args:
            fasta_file (Path): Path to input FASTA file
            method (str): Alignment method ("mafft" or "muscle")
            
        Returns:
            Path: Path to alignment file
        """
        family = fasta_file.stem
        alignment_file = self.alignment_dir / f"{family}_aligned.fasta"
        
        if method.lower() == "mafft":
            cmd = [
                "mafft",
                "--adjustdirection",  # Detect and correct reverse complement sequences
                "--op", "10.0",       # Much higher gap opening penalty to discourage gaps
                "--ep", "1.0",        # Higher gap extension penalty
                "--maxiterate", "1000",  # Allow sufficient refinement iterations
                "--localpair",        # Use L-INS-i for higher accuracy
                "--quiet",            # Reduce output
                str(fasta_file)
            ]
            
            try:
                with open(alignment_file, 'w') as output_file:
                    result = subprocess.run(cmd, check=True, stdout=output_file, 
                                          stderr=subprocess.PIPE, text=True)
                print(f"✓ Aligned sequences for {family} using MAFFT")
                return alignment_file
            except subprocess.CalledProcessError as e:
                print(f"✗ Error aligning {family} with MAFFT: {e}")
                if e.stderr:
                    print(f"   stderr: {e.stderr}")
                return None
                
        elif method.lower() == "muscle":
            # Check MUSCLE version and use appropriate syntax
            try:
                # Test MUSCLE version
                version_result = subprocess.run(["muscle", "-version"], 
                                              capture_output=True, text=True, timeout=10)
                
                if version_result.returncode == 0 and "5." in version_result.stdout:
                    # MUSCLE 5.x syntax
                    cmd = [
                        "muscle",
                        "-align", str(fasta_file),
                        "-output", str(alignment_file)
                    ]
                else:
                    # MUSCLE 3.x/4.x syntax (fallback)
                    cmd = [
                        "muscle",
                        "-in", str(fasta_file),
                        "-out", str(alignment_file)
                    ]
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError):
                # If version check fails, try MUSCLE 5.x syntax first (most likely on modern systems)
                cmd = [
                    "muscle",
                    "-align", str(fasta_file),
                    "-output", str(alignment_file)
                ]
            
            try:
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                print(f"✓ Aligned sequences for {family} using MUSCLE")
                return alignment_file
            except subprocess.CalledProcessError as e:
                # If MUSCLE 5.x syntax failed, try the old syntax
                if "-align" in cmd:
                    print(f"   MUSCLE 5.x syntax failed, trying legacy syntax...")
                    cmd_legacy = [
                        "muscle",
                        "-in", str(fasta_file),
                        "-out", str(alignment_file)
                    ]
                    try:
                        result = subprocess.run(cmd_legacy, check=True, capture_output=True, text=True)
                        print(f"✓ Aligned sequences for {family} using MUSCLE (legacy syntax)")
                        return alignment_file
                    except subprocess.CalledProcessError as e2:
                        print(f"✗ Error aligning {family} with both MUSCLE syntaxes:")
                        print(f"   MUSCLE 5.x error: {e}")
                        print(f"   Legacy error: {e2}")
                        if e.stderr:
                            print(f"   MUSCLE 5.x stderr: {e.stderr}")
                        if e2.stderr:
                            print(f"   Legacy stderr: {e2.stderr}")
                        return None
                else:
                    print(f"✗ Error aligning {family} with MUSCLE: {e}")
                    if e.stderr:
                        print(f"   stderr: {e.stderr}")
                    return None
        else:
            print(f"✗ Unknown alignment method: {method}")
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
    
    def run_pipeline(self, min_otus=3, max_families=None, alignment_method="mafft"):
        """
        Run the complete phylogenetic analysis pipeline.
        
        Args:
            min_otus (int): Minimum number of OTUs required for analysis
            max_families (int): Maximum number of families to process (for testing)
            alignment_method (str): Alignment method ("mafft" or "muscle")
        """
        print("Starting phylogenetic analysis pipeline...")
        print(f"Database: {self.db_path}")
        print(f"Output directory: {self.output_dir}")
        print("="*60)
        
        # Get sequences for each family (including multiple species per OTU)
        print("1. Extracting sequences from each OTU (including multiple species)...")
        family_sequences = self.get_otu_sequences()
        
        print(f"Found {len(family_sequences)} families:")
        for family, sequences in family_sequences.items():
            # Count unique OTUs and unique species
            unique_otus = len(set(seq[2]['otu_id'] for seq in sequences))
            unique_species = len(set(f"{seq[2]['genus']} {seq[2]['species']}" for seq in sequences))
            print(f"  - {family}: {len(sequences)} sequences, {unique_otus} OTUs, {unique_species} species")
        
        # Filter families with sufficient sequences
        families_to_analyze = {
            family: sequences for family, sequences in family_sequences.items() 
            if len(set(seq[2]['otu_id'] for seq in sequences)) >= min_otus
        }
        
        if max_families:
            families_to_analyze = dict(list(families_to_analyze.items())[:max_families])
        
        print(f"\nAnalyzing {len(families_to_analyze)} families with ≥{min_otus} OTUs")
        print("="*60)
        
        results = {}
        
        for i, (family, sequences) in enumerate(families_to_analyze.items(), 1):
            unique_otus = len(set(seq[2]['otu_id'] for seq in sequences))
            unique_species = len(set(f"{seq[2]['genus']} {seq[2]['species']}" for seq in sequences))
            
            print(f"\n{i}. Processing {family} ({len(sequences)} sequences, {unique_otus} OTUs, {unique_species} species)")
            
            # Select outgroups
            outgroups = self.select_outgroups(family_sequences, family)
            print(f"   Selected {len(outgroups)} outgroups")
            
            if len(outgroups) == 0:
                print(f"   ⚠ No suitable outgroups found for {family}")
                continue
            
            # Create FASTA file
            fasta_file = self.create_fasta_file(family, sequences, outgroups)
            print(f"   Created FASTA: {fasta_file.name}")
            
            # Align sequences
            alignment_file = self.align_sequences(fasta_file, alignment_method)
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
                'num_sequences': len(sequences),
                'num_otus': unique_otus,
                'num_species': unique_species,
                'num_outgroups': len(outgroups)
            }
        
        print("\n" + "="*60)
        print("PIPELINE COMPLETE")
        print("="*60)
        print(f"Successfully processed {len(results)} families:")
        
        for family, result in results.items():
            print(f"  - {family}: {result['num_sequences']} sequences, "
                  f"{result['num_otus']} OTUs, {result['num_species']} species, "
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
                f.write(f"  Total sequences: {result['num_sequences']}\n")
                f.write(f"  Unique OTUs: {result['num_otus']}\n")
                f.write(f"  Unique species: {result['num_species']}\n")
                f.write(f"  Outgroups: {result['num_outgroups']}\n")
                f.write(f"  Tree file: {result['tree'].name}\n\n")
        
        print(f"Summary report saved to: {report_file}")


def main():
    """Main function to run the pipeline with command-line arguments."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Phylogenetic Analysis Pipeline for BOLD Database")
    parser.add_argument("--database", required=True, help="Path to SQLite database file")
    parser.add_argument("--output-dir", default="phylogenies", help="Output directory for results")
    parser.add_argument("--min-otus", type=int, default=3, help="Minimum OTUs per family for analysis")
    parser.add_argument("--max-families", type=int, help="Maximum families to process (for testing)")
    parser.add_argument("--selection-strategy", default="best_quality", 
                       choices=["longest", "random", "best_quality"],
                       help="Strategy for selecting OTU representatives")
    parser.add_argument("--num-outgroups", type=int, default=3, help="Number of outgroup sequences")
    parser.add_argument("--alignment-method", default="mafft", 
                       choices=["mafft", "muscle"], help="Alignment method")
    parser.add_argument("--tree-method", default="fasttree",
                       choices=["fasttree", "iqtree"], help="Tree building method")
    parser.add_argument("--bootstrap", type=int, default=1000, help="Bootstrap replicates (IQ-TREE)")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPU threads")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    
    args = parser.parse_args()
    
    # Print system and tool information for debugging
    print("System Information:")
    print(f"Python version: {sys.version}")
    print(f"Working directory: {os.getcwd()}")
    print(f"PATH: {os.environ.get('PATH', 'Not set')}")
    
    # Check if required external programs are available
    required_programs = []
    if args.alignment_method == "mafft":
        required_programs.append("mafft")
    elif args.alignment_method == "muscle":
        required_programs.append("muscle")
    
    if args.tree_method == "iqtree":
        required_programs.append("iqtree")
    elif args.tree_method == "fasttree":
        required_programs.extend(["fasttree", "FastTree"])  # Try both names
    
    missing_programs = []
    program_versions = {}
    
    for program in required_programs:
        if program in ["fasttree", "FastTree"]:
            # Special handling for FastTree - try both names
            found = False
            for ft_name in ["fasttree", "FastTree"]:
                try:
                    result = subprocess.run([ft_name], capture_output=True, timeout=5, text=True)
                    found = True
                    program_versions[program] = f"{ft_name} found"
                    break
                except (FileNotFoundError, subprocess.TimeoutExpired):
                    continue
            if not found:
                missing_programs.append("fasttree/FastTree")
        else:
            try:
                if program == "muscle":
                    # Check MUSCLE version specifically
                    version_result = subprocess.run([program, "-version"], capture_output=True, timeout=5, text=True)
                    if version_result.returncode == 0:
                        program_versions[program] = version_result.stdout.strip()
                    else:
                        # Try alternative version commands
                        help_result = subprocess.run([program, "-h"], capture_output=True, timeout=5, text=True)
                        program_versions[program] = f"Available (help accessible)"
                else:
                    subprocess.run([program, '--help'], capture_output=True, timeout=5)
                    program_versions[program] = "Available"
            except (FileNotFoundError, subprocess.TimeoutExpired):
                missing_programs.append(program)
    
    print("\nTool Availability:")
    for program, version in program_versions.items():
        print(f"  {program}: {version}")
    
    if missing_programs:
        print(f"\nMissing required programs: {', '.join(missing_programs)}")
        print("\nInstallation instructions:")
        if "mafft" in missing_programs:
            print("- MAFFT: conda install -c bioconda mafft")
        if "muscle" in missing_programs:
            print("- MUSCLE: conda install -c bioconda muscle")
        if "iqtree" in missing_programs:
            print("- IQ-TREE: conda install -c bioconda iqtree")
        if "fasttree/FastTree" in missing_programs:
            print("- FastTree: conda install -c bioconda fasttree")
        return 1
    
    # Run pipeline
    pipeline = PhylogeneticPipeline(args.database, args.output_dir)
    results = pipeline.run_pipeline(
        min_otus=args.min_otus, 
        max_families=args.max_families,
        alignment_method=args.alignment_method
    )
    
    if results:
        pipeline.generate_summary_report(results)
        print(f"Phylogenetic analysis completed successfully!")
        print(f"Generated {len(results)} family phylogenies")
        print(f"Output directory: {args.output_dir}")
        return 0
    else:
        print("No families were successfully processed.")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
