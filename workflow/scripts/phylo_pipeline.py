#!/usr/bin/env python3
"""
Phylogenetic Analysis Pipeline for BOLD Database
================================================

This script generates phylogenies for each family using OTU representatives,
with outgroups selected from other families within the same order.
Includes BIN-species conflict analysis and PDF visualization.

Requirements:
- biopython
- pandas
- ete3 (for tree manipulation and PDF generation)
- muscle/mafft (alignment)
- iqtree/fasttree (phylogeny reconstruction)
"""

import sqlite3
import pandas as pd
import os
import subprocess
import sys
import json
import gc
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import random
from datetime import datetime

class PhylogeneticPipeline:
    def __init__(self, db_path, output_dir="phylogenies", generate_pdfs=False, bin_conflict_analysis=False):
        """
        Initialize the phylogenetic pipeline.
        
        Args:
            db_path (str): Path to the BOLD SQLite database
            output_dir (str): Output directory for results
            generate_pdfs (bool): Whether to generate PDF visualizations
            bin_conflict_analysis (bool): Whether to analyze BIN-species conflicts
        """
        self.db_path = db_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.generate_pdfs = generate_pdfs
        self.bin_conflict_analysis = bin_conflict_analysis
        
        # Import ETE3 only if PDF generation is requested
        if self.generate_pdfs:
            try:
                # Set multiple environment variables for headless operation
                import os
                os.environ['QT_QPA_PLATFORM'] = 'offscreen'
                os.environ['DISPLAY'] = ':99'  # Virtual display
                os.environ['MPLBACKEND'] = 'Agg'  # Force matplotlib non-GUI backend
                
                # Import matplotlib first with headless backend
                import matplotlib
                matplotlib.use('Agg')
                
                from ete3 import Tree, TreeStyle, TextFace, CircleFace
                self.Tree = Tree
                self.TreeStyle = TreeStyle
                self.TextFace = TextFace
                self.CircleFace = CircleFace
                print("✓ ETE3 loaded with headless backend")
            except ImportError:
                print("Warning: ETE3 not available. PDF generation disabled.")
                self.generate_pdfs = False
            except Exception as e:
                print(f"Warning: ETE3 backend issue: {e}. PDF generation disabled.")
                self.generate_pdfs = False
        
        # Performance settings
        self.cleanup_intermediates = True
        self.batch_size = 50  # Process families in batches for memory efficiency
    
    def get_families_to_process(self, min_otus=3):
        """
        Get families meeting criteria without loading all data.
        HPC-optimized for large datasets.
        
        Args:
            min_otus (int): Minimum number of OTUs required for analysis
            
        Returns:
            pandas.DataFrame: Families with their statistics
        """
        conn = sqlite3.connect(self.db_path)
        
        query = """
        SELECT 
            b.family,
            b.`order`,
            COUNT(DISTINCT o.otu_id) as otu_count,
            COUNT(DISTINCT b.species) as species_count,
            COUNT(*) as total_records
        FROM bold b
        JOIN bold_otus o ON b.recordid = o.recordid
        WHERE b.family IS NOT NULL 
            AND b.`order` IS NOT NULL
            AND b.nuc IS NOT NULL 
            AND LENGTH(b.nuc) > 200
            AND b.genus IS NOT NULL
            AND b.species IS NOT NULL
        GROUP BY b.family, b.`order`
        HAVING otu_count >= ?
        ORDER BY otu_count DESC
        """
        
        families_df = pd.read_sql_query(query, conn, params=[min_otus])
        conn.close()
        
        print(f"Found {len(families_df)} families with ≥{min_otus} OTUs")
        return families_df
    
    def get_family_sequences_optimized(self, family_name):
        """
        Get sequences for a single family with optimized query.
        
        Args:
            family_name (str): Name of the family
            
        Returns:
            pandas.DataFrame: Family sequence data
        """
        conn = sqlite3.connect(self.db_path)
        
        # Single family query with proper indexing
        query = """
        SELECT DISTINCT
            o.otu_id,
            b.nuc,
            b.genus,
            b.species,
            b.processid,
            b.recordid,
            b.bin_uri,
            b.`order`,
            b.family
        FROM bold b
        JOIN bold_otus o ON b.recordid = o.recordid
        WHERE b.family = ?
            AND b.nuc IS NOT NULL 
            AND LENGTH(b.nuc) > 200
            AND b.genus IS NOT NULL
            AND b.species IS NOT NULL
        ORDER BY o.otu_id, LENGTH(b.nuc) DESC
        """
        
        df = pd.read_sql_query(query, conn, params=[family_name])
        conn.close()
        return df
    
    def select_representatives_optimized(self, family_df):
        """
        Select representative sequences (one per OTU-species combination).
        
        Args:
            family_df (pandas.DataFrame): Family sequence data
            
        Returns:
            pandas.DataFrame: Representative sequences
        """
        representatives = []
        
        # Group by OTU and species combination
        for (otu_id, species), group in family_df.groupby(['otu_id', 'species']):
            # Select the longest sequence for each OTU-species combination
            rep = group.loc[group['nuc'].str.len().idxmax()]
            representatives.append(rep)
        
        return pd.DataFrame(representatives)
    
    def get_outgroups_optimized(self, family_order, family_name, num_outgroups=3):
        """
        Select outgroup sequences from other families in the same order.
        
        Args:
            family_order (str): Order of the target family
            family_name (str): Family to find outgroups for
            num_outgroups (int): Number of outgroup sequences to select
            
        Returns:
            pandas.DataFrame: Outgroup sequences
        """
        conn = sqlite3.connect(self.db_path)
        
        query = """
        SELECT DISTINCT
            o.otu_id,
            b.nuc,
            b.genus,
            b.species,
            b.processid,
            b.recordid,
            b.bin_uri,
            b.family
        FROM bold b
        JOIN bold_otus o ON b.recordid = o.recordid
        WHERE b.`order` = ?
            AND b.family != ?
            AND b.family IS NOT NULL
            AND b.nuc IS NOT NULL 
            AND LENGTH(b.nuc) > 200
            AND b.genus IS NOT NULL
            AND b.species IS NOT NULL
        ORDER BY RANDOM()
        LIMIT ?
        """
        
        outgroups_df = pd.read_sql_query(query, conn, params=[family_order, family_name, num_outgroups * 5])
        conn.close()
        
        # Select diverse outgroups (one per family if possible)
        outgroups = []
        used_families = set()
        
        for _, row in outgroups_df.iterrows():
            if len(outgroups) >= num_outgroups:
                break
            
            if row['family'] not in used_families:
                outgroups.append(row)
                used_families.add(row['family'])
        
        # If we don't have enough diverse families, add more from same families
        if len(outgroups) < num_outgroups:
            for _, row in outgroups_df.iterrows():
                if len(outgroups) >= num_outgroups:
                    break
                if row.name not in [og.name for og in outgroups]:
                    outgroups.append(row)
        
        return pd.DataFrame(outgroups[:num_outgroups]) if outgroups else pd.DataFrame()
    
    def parse_sequence_metadata(self, tree):
        """
        Extract metadata from sequence IDs for BIN-species analysis.
        
        Args:
            tree: ETE3 Tree object
            
        Returns:
            dict: Metadata for each sequence
        """
        metadata = {}
        
        for leaf in tree:
            if leaf.name.startswith('OUTGROUP-'):
                # Parse: OUTGROUP-species-otu-bin_uri-processid
                parts = leaf.name.split('-', 4)  # Split on first 4 hyphens
                is_outgroup = True
                species = parts[1] if len(parts) > 1 else 'Unknown'
                otu = parts[2] if len(parts) > 2 else 'NA'
                bin_uri = parts[3] if len(parts) > 3 else 'NA'
                processid = parts[4] if len(parts) > 4 else 'NA'
            else:
                # Parse: species-otu-bin_uri-processid
                parts = leaf.name.split('-', 3)  # Split on first 3 hyphens
                is_outgroup = False
                species = parts[0] if len(parts) > 0 else 'Unknown'
                otu = parts[1] if len(parts) > 1 else 'NA'
                bin_uri = parts[2] if len(parts) > 2 else 'NA'
                processid = parts[3] if len(parts) > 3 else 'NA'
            
            metadata[leaf.name] = {
                'species': species,
                'otu': otu,
                'bin_uri': bin_uri,
                'processid': processid,
                'is_outgroup': is_outgroup
            }
        
        return metadata
    
    def analyze_bin_species_discordance(self, metadata):
        """
        Identify BIN-species conflicts.
        
        Args:
            metadata (dict): Sequence metadata
            
        Returns:
            dict: Analysis of conflicts
        """
        # Group by BIN to find multi-species BINs
        bins_to_species = defaultdict(set)
        # Group by species to find multi-BIN species
        species_to_bins = defaultdict(set)
        
        for seq_id, data in metadata.items():
            if not data['is_outgroup'] and data['bin_uri'] != 'NA':
                bins_to_species[data['bin_uri']].add(data['species'])
                species_to_bins[data['species']].add(data['bin_uri'])
        
        # Identify conflicts
        multi_species_bins = {bin_uri: species_set 
                             for bin_uri, species_set in bins_to_species.items() 
                             if len(species_set) > 1}
        
        multi_bin_species = {species: bin_set 
                            for species, bin_set in species_to_bins.items() 
                            if len(bin_set) > 1}
        
        return {
            'multi_species_bins': multi_species_bins,
            'multi_bin_species': multi_bin_species,
            'bins_to_species': dict(bins_to_species),
            'species_to_bins': dict(species_to_bins)
        }
    
    def create_family_output_dir(self, family_name):
        """
        Create family-specific output directory structure.
        
        Args:
            family_name (str): Name of the family
            
        Returns:
            dict: Paths for family outputs
        """
        family_dir = self.output_dir / family_name
        family_dir.mkdir(exist_ok=True)
        
        return {
            'family_dir': family_dir,
            'fasta_file': family_dir / f"{family_name}.fasta",
            'alignment_file': family_dir / f"{family_name}_aligned.fasta",
            'tree_file': family_dir / f"{family_name}.treefile",
            'pdf_file': family_dir / f"{family_name}_tree.pdf",
            'metadata_file': family_dir / f"tree_metadata.json"
        }
    def create_fasta_file(self, representatives_df, outgroups_df, fasta_path):
        """
        Create a FASTA file for a family with ingroup and outgroup sequences.
        Uses the branch labeling format: species-otu-bin_uri-processid
        
        Args:
            representatives_df (pandas.DataFrame): Ingroup representative sequences
            outgroups_df (pandas.DataFrame): Outgroup sequences
            fasta_path (Path): Path to output FASTA file
            
        Returns:
            bool: Success status
        """
        try:
            records = []
            
            # Add ingroup sequences with labeling format
            for _, row in representatives_df.iterrows():
                # Clean values to handle None or empty values
                species = str(row.get('species', 'Unknown')).replace(' ', '_')
                otu = str(row.get('otu_id', 'NA')).replace(' ', '_')
                bin_uri = str(row.get('bin_uri', 'NA')).replace(' ', '_')
                processid = str(row.get('processid', 'NA')).replace(' ', '_')
                
                # Create sequence ID in format: species-otu-bin_uri-processid
                seq_id = f"{species}-{otu}-{bin_uri}-{processid}"
                
                # Clean any problematic characters for phylogenetic software
                seq_id = seq_id.replace('(', '').replace(')', '').replace('[', '').replace(']', '')
                seq_id = seq_id.replace(',', '').replace(';', '').replace(':', '').replace('|', '-')
                
                record = SeqRecord(
                    Seq(row['nuc']),
                    id=seq_id,
                    description=f"Ingroup|{row.get('genus', '')}_{species}|{row.get('family', '')}"
                )
                records.append(record)
            
            # Add outgroup sequences with OUTGROUP prefix
            for _, row in outgroups_df.iterrows():
                # Clean values to handle None or empty values
                species = str(row.get('species', 'Unknown')).replace(' ', '_')
                otu = str(row.get('otu_id', 'NA')).replace(' ', '_')
                bin_uri = str(row.get('bin_uri', 'NA')).replace(' ', '_')
                processid = str(row.get('processid', 'NA')).replace(' ', '_')
                
                # Create sequence ID with OUTGROUP prefix
                seq_id = f"OUTGROUP-{species}-{otu}-{bin_uri}-{processid}"
                
                # Clean any problematic characters for phylogenetic software
                seq_id = seq_id.replace('(', '').replace(')', '').replace('[', '').replace(']', '')
                seq_id = seq_id.replace(',', '').replace(';', '').replace(':', '').replace('|', '-')
                
                record = SeqRecord(
                    Seq(row['nuc']),
                    id=seq_id,
                    description=f"Outgroup|{row.get('genus', '')}_{species}|{row.get('family', 'Unknown')}"
                )
                records.append(record)
            
            SeqIO.write(records, fasta_path, "fasta")
            return True
            
        except Exception as e:
            print(f"Error creating FASTA file: {e}")
            return False
    
    def generate_bin_conflict_pdf(self, tree_file, pdf_file, family_name, metadata_file):
        """
        Generate PDF with BIN-species conflict visualization.
        
        Args:
            tree_file (Path): Path to tree file
            pdf_file (Path): Path to output PDF
            family_name (str): Name of the family
            metadata_file (Path): Path to metadata file
            
        Returns:
            bool: Success status
        """
        if not self.generate_pdfs:
            return True  # Skip if PDF generation disabled
            
        try:
            # Load and root tree
            tree = self.Tree(str(tree_file))
            
            # Find outgroup and root
            outgroup_nodes = [leaf for leaf in tree if leaf.name.startswith('OUTGROUP-')]
            if outgroup_nodes:
                tree.set_outgroup(outgroup_nodes[0])
            
            # Parse metadata and analyze conflicts
            metadata = self.parse_sequence_metadata(tree)
            conflicts = self.analyze_bin_species_discordance(metadata)
            
            # Style tree for conflicts
            self._style_tree_for_bin_conflicts(tree, metadata, conflicts)
            
            # Generate PDF
            success = self._generate_annotated_pdf(tree, family_name, conflicts, pdf_file)
            
            # Save metadata
            if success:
                self._save_tree_metadata(metadata, conflicts, family_name, metadata_file)
            
            return success
            
        except Exception as e:
            print(f"Error generating PDF for {family_name}: {e}")
            return False
    
    def _style_tree_for_bin_conflicts(self, tree, metadata, conflicts):
        """Apply visual styling to highlight BIN-species conflicts"""
        
        # Define color schemes
        colors = {
            'multi_species_bin': '#D32F2F',    # Dark red for BIN sharing
            'multi_bin_species': '#1976D2',    # Blue for BIN splitting  
            'both_conflicts': '#F57C00',       # Orange for both issues
            'normal': '#388E3C',               # Dark green for normal
            'outgroup': '#7B1FA2'              # Purple for outgroups
        }
        
        for leaf in tree:
            data = metadata[leaf.name]
            
            if data['is_outgroup']:
                leaf_color = colors['outgroup']
                conflict_type = "Outgroup"
            else:
                species = data['species']
                bin_uri = data['bin_uri']
                
                is_multi_species_bin = bin_uri in conflicts['multi_species_bins']
                is_multi_bin_species = species in conflicts['multi_bin_species']
                
                if is_multi_species_bin and is_multi_bin_species:
                    leaf_color = colors['both_conflicts']
                    conflict_type = "BIN sharing + BIN splitting"
                elif is_multi_species_bin:
                    leaf_color = colors['multi_species_bin']
                    conflict_type = "BIN sharing"
                elif is_multi_bin_species:
                    leaf_color = colors['multi_bin_species']
                    conflict_type = "BIN splitting"
                else:
                    leaf_color = colors['normal']
                    conflict_type = "Normal"
            
            # Style the leaf
            leaf.add_features(
                bgcolor=leaf_color,
                conflict_type=conflict_type,
                species_clean=data['species'],
                bin_clean=data['bin_uri']
            )
    
    def _generate_annotated_pdf(self, tree, family_name, conflicts, pdf_file):
        """Generate PDF with detailed annotations"""
        try:
            # Create tree style
            ts = self.TreeStyle()
            ts.show_leaf_name = False  # We'll add custom labels
            ts.show_branch_support = True
            ts.title.add_face(self.TextFace(f"Family: {family_name}", fsize=16), column=0)
            
            # Add conflict summary
            summary_text = f"""Conflicts Summary:
• Multi-species BINs: {len(conflicts['multi_species_bins'])}
• Multi-BIN species: {len(conflicts['multi_bin_species'])}"""
            ts.title.add_face(self.TextFace(summary_text, fsize=12), column=0)
            
            # Custom node styling function
            def layout(node):
                if node.is_leaf():
                    # Add colored circle
                    circle = self.CircleFace(radius=8, color=node.bgcolor, style="sphere")
                    node.add_face(circle, column=0, position="branch-right")
                    
                    # Add species name
                    species_face = self.TextFace(f" {node.species_clean}", fsize=10)
                    node.add_face(species_face, column=1, position="branch-right")
                    
                    # Add BIN info
                    bin_face = self.TextFace(f" [{node.bin_clean}]", fsize=8, fgcolor="gray")
                    node.add_face(bin_face, column=2, position="branch-right")
                    
                    # Add conflict indicator
                    if node.conflict_type not in ["Normal", "Outgroup"]:
                        conflict_face = self.TextFace(f" ⚠ {node.conflict_type}", fsize=8, fgcolor="red")
                        node.add_face(conflict_face, column=3, position="branch-right")
            
            ts.layout_fn = layout
            
            # Generate PDF
            tree.render(str(pdf_file), tree_style=ts, dpi=300)
            return True
            
        except Exception as e:
            print(f"Error rendering PDF: {e}")
            return False
    
    def _save_tree_metadata(self, metadata, conflicts, family_name, metadata_file):
        """Save tree metadata to JSON file"""
        try:
            tree_metadata = {
                'family': family_name,
                'generated_at': datetime.now().isoformat(),
                'total_sequences': len(metadata),
                'ingroup_sequences': len([m for m in metadata.values() if not m['is_outgroup']]),
                'outgroup_sequences': len([m for m in metadata.values() if m['is_outgroup']]),
                'conflicts': {
                    'multi_species_bins': len(conflicts['multi_species_bins']),
                    'multi_bin_species': len(conflicts['multi_bin_species']),
                    'multi_species_bin_details': {k: list(v) for k, v in conflicts['multi_species_bins'].items()},
                    'multi_bin_species_details': {k: list(v) for k, v in conflicts['multi_bin_species'].items()}
                }
            }
            
            with open(metadata_file, 'w') as f:
                json.dump(tree_metadata, f, indent=2)
                
        except Exception as e:
            print(f"Error saving metadata: {e}")
    
    def process_single_family(self, family_name, family_order, min_otus, num_outgroups=3, 
                             alignment_method="mafft", tree_method="iqtree", bootstrap=1000):
        """
        Process a single family efficiently.
        
        Args:
            family_name (str): Name of the family
            family_order (str): Order of the family
            min_otus (int): Minimum OTUs required
            num_outgroups (int): Number of outgroups
            alignment_method (str): Alignment method
            tree_method (str): Tree building method
            bootstrap (int): Bootstrap replicates
            
        Returns:
            dict: Processing results or None if failed
        """
        # Get output paths
        paths = self.create_family_output_dir(family_name)
        
        # Skip if already processed (optional optimization)
        if (paths['pdf_file'].exists() if self.generate_pdfs else paths['tree_file'].exists()):
            print(f"   Already processed, skipping {family_name}")
            return None
        
        try:
            # 1. Get sequences for this family
            family_df = self.get_family_sequences_optimized(family_name)
            
            if len(family_df) == 0:
                print(f"   No sequences found for {family_name}")
                return None
            
            # 2. Select representatives (one per OTU-species combination)
            representatives = self.select_representatives_optimized(family_df)
            
            # 3. Check if we have enough OTUs
            unique_otus = representatives['otu_id'].nunique()
            if unique_otus < min_otus:
                print(f"   Insufficient OTUs: {unique_otus} < {min_otus}")
                return None
            
            # 4. Get outgroups from same order
            outgroups = self.get_outgroups_optimized(family_order, family_name, num_outgroups)
            
            if len(outgroups) == 0:
                print(f"   No outgroups found for {family_name}")
                return None
            
            # 5. Create FASTA
            if not self.create_fasta_file(representatives, outgroups, paths['fasta_file']):
                return None
            
            # 6. Align sequences
            if not self.align_sequences_optimized(paths['fasta_file'], paths['alignment_file'], alignment_method):
                return None
            
            # 7. Build tree
            if not self.build_tree_optimized(paths['alignment_file'], paths['tree_file'], tree_method, bootstrap):
                return None
            
            # 8. Generate PDF if requested
            if self.generate_pdfs:
                if not self.generate_bin_conflict_pdf(paths['tree_file'], paths['pdf_file'], 
                                                    family_name, paths['metadata_file']):
                    print(f"   Warning: PDF generation failed for {family_name}")
            
            # 9. Clean up intermediate files to save space
            if self.cleanup_intermediates:
                paths['fasta_file'].unlink(missing_ok=True)
                paths['alignment_file'].unlink(missing_ok=True)
            
            return {
                'family': family_name,
                'order': family_order,
                'num_otus': unique_otus,
                'num_species': representatives['species'].nunique(),
                'num_sequences': len(representatives),
                'num_outgroups': len(outgroups),
                'tree_file': paths['tree_file'],
                'pdf_file': paths['pdf_file'] if self.generate_pdfs else None
            }
            
        except Exception as e:
            print(f"   Error processing {family_name}: {e}")
            return None
            
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
    
    def run_pipeline_hpc_optimized(self, min_otus=3, max_families=None, alignment_method="mafft", 
                                   tree_method="iqtree", bootstrap=1000, num_outgroups=3):
        """
        Run the HPC-optimized phylogenetic analysis pipeline.
        
        Args:
            min_otus (int): Minimum number of OTUs required for analysis
            max_families (int): Maximum number of families to process (for testing)
            alignment_method (str): Alignment method ("mafft" or "muscle")
            tree_method (str): Tree building method ("iqtree" or "fasttree")
            bootstrap (int): Bootstrap replicates
            num_outgroups (int): Number of outgroup sequences
        """
        print("Starting HPC-optimized phylogenetic analysis pipeline...")
        print(f"Database: {self.db_path}")
        print(f"Output directory: {self.output_dir}")
        print(f"Generate PDFs: {self.generate_pdfs}")
        print(f"BIN conflict analysis: {self.bin_conflict_analysis}")
        print("="*60)
        
        # Get families to process (single optimized query)
        print("1. Discovering families to process...")
        families_df = self.get_families_to_process(min_otus)
        
        if max_families and max_families > 0:
            families_df = families_df.head(max_families)
            print(f"   Limited to first {max_families} families for testing")
        
        print(f"Processing {len(families_df)} families:")
        for _, row in families_df.iterrows():
            print(f"  - {row['family']}: {row['otu_count']} OTUs, {row['species_count']} species")
        
        print("="*60)
        
        results = []
        failed_families = []
        
        # Process families in batches for memory efficiency
        for i, (_, family_row) in enumerate(families_df.iterrows(), 1):
            family_name = family_row['family']
            family_order = family_row['order']
            
            print(f"\n[{i}/{len(families_df)}] Processing {family_name}")
            
            try:
                result = self.process_single_family(
                    family_name, family_order, min_otus, num_outgroups,
                    alignment_method, tree_method, bootstrap
                )
                
                if result:
                    results.append(result)
                    print(f"✓ Success: {family_name}")
                else:
                    failed_families.append(family_name)
                    print(f"⚠ Skipped: {family_name}")
                    
            except Exception as e:
                failed_families.append(family_name)
                print(f"✗ Error: {family_name} - {e}")
                
            # Memory cleanup every 10 families
            if i % 10 == 0:
                gc.collect()
                print(f"   Memory cleanup after {i} families")
        
        # Generate summary
        self.generate_pipeline_summary(results, failed_families)
        return results
    
    def align_sequences_optimized(self, fasta_file, alignment_file, method="mafft"):
        """Optimized sequence alignment"""
        try:
            if method.lower() == "mafft":
                cmd = [
                    "mafft",
                    "--adjustdirection",  # Detect reverse complement
                    "--op", "10.0",       # Gap opening penalty
                    "--ep", "1.0",        # Gap extension penalty
                    "--maxiterate", "1000",
                    "--localpair",
                    "--quiet",
                    str(fasta_file)
                ]
                
                with open(alignment_file, 'w') as output_file:
                    result = subprocess.run(cmd, check=True, stdout=output_file, 
                                          stderr=subprocess.PIPE, text=True)
                return True
                
            elif method.lower() == "muscle":
                # Try MUSCLE 5.x syntax first
                cmd = ["muscle", "-align", str(fasta_file), "-output", str(alignment_file)]
                try:
                    subprocess.run(cmd, check=True, capture_output=True, text=True)
                    return True
                except subprocess.CalledProcessError:
                    # Fallback to legacy syntax
                    cmd = ["muscle", "-in", str(fasta_file), "-out", str(alignment_file)]
                    subprocess.run(cmd, check=True, capture_output=True, text=True)
                    return True
            
        except Exception as e:
            print(f"   Alignment error: {e}")
            return False
    
    def build_tree_optimized(self, alignment_file, tree_file, method="iqtree", bootstrap=1000):
        """Optimized tree building with enhanced error handling"""
        try:
            # Pre-flight checks for alignment file
            if not alignment_file.exists():
                print(f"   Error: Alignment file not found: {alignment_file}")
                return False
            
            # Check alignment file content
            try:
                with open(alignment_file, 'r') as f:
                    content = f.read().strip()
                    if not content:
                        print(f"   Error: Empty alignment file: {alignment_file}")
                        return False
                    
                    # Count sequences
                    seq_count = content.count('>')
                    if seq_count < 4:
                        print(f"   Error: Insufficient sequences for tree building ({seq_count} < 4): {alignment_file}")
                        return False
                    
                    print(f"   Alignment check: {seq_count} sequences found")
                    
            except Exception as e:
                print(f"   Error reading alignment file: {e}")
                return False
            
            if method.lower() == "iqtree":
                tree_prefix = tree_file.parent / tree_file.stem
                cmd = [
                    "iqtree",
                    "-s", str(alignment_file),
                    "-pre", str(tree_prefix),
                    "-nt", "AUTO"
                ]
                
                # Add bootstrap only if >= 1000 (IQ-TREE requirement)
                if bootstrap >= 1000:
                    cmd.extend(["-bb", str(bootstrap)])
                    print(f"   Using {bootstrap} bootstrap replicates")
                elif bootstrap > 0:
                    print(f"   Warning: Bootstrap {bootstrap} < 1000, using standard bootstrap instead")
                    cmd.extend(["-b", str(bootstrap)])
                else:
                    print(f"   Bootstrap disabled")
                
                # Run with error capture
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    print(f"   IQ-TREE error (exit code {result.returncode}):")
                    if result.stderr:
                        print(f"   STDERR: {result.stderr}")
                    if result.stdout:
                        print(f"   STDOUT: {result.stdout}")
                    return False
                
                # IQ-TREE creates .treefile, move to expected location
                iqtree_output = tree_prefix.with_suffix(".treefile")
                if iqtree_output.exists() and iqtree_output != tree_file:
                    iqtree_output.rename(tree_file)
                
                return True
                
            elif method.lower() == "fasttree":
                cmd = ["fasttree", "-nt", str(alignment_file)]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    print(f"   FastTree error (exit code {result.returncode}):")
                    if result.stderr:
                        print(f"   STDERR: {result.stderr}")
                    return False
                
                with open(tree_file, 'w') as output_file:
                    output_file.write(result.stdout)
                    
                return True
                
        except Exception as e:
            print(f"   Tree building error: {e}")
            return False
    
    def generate_pipeline_summary(self, results, failed_families):
        """Generate a summary report of the analysis"""
        summary_file = self.output_dir / "phylogenetic_analysis_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("HPC-Optimized Phylogenetic Analysis Summary\n")
            f.write("="*60 + "\n\n")
            f.write(f"Database: {self.db_path}\n")
            f.write(f"Output directory: {self.output_dir}\n")
            f.write(f"Generated at: {datetime.now().isoformat()}\n")
            f.write(f"PDF generation: {self.generate_pdfs}\n")
            f.write(f"BIN conflict analysis: {self.bin_conflict_analysis}\n\n")
            
            f.write(f"RESULTS:\n")
            f.write(f"Successfully processed: {len(results)} families\n")
            f.write(f"Failed/skipped: {len(failed_families)} families\n\n")
            
            if results:
                f.write("Successful families:\n")
                for result in results:
                    f.write(f"  - {result['family']}: {result['num_otus']} OTUs, "
                           f"{result['num_species']} species, {result['num_outgroups']} outgroups\n")
            
            if failed_families:
                f.write(f"\nFailed/skipped families:\n")
                for family in failed_families:
                    f.write(f"  - {family}\n")
        
        print(f"\nSummary report saved to: {summary_file}")
        print(f"Successfully processed {len(results)} families")
        if failed_families:
            print(f"Failed/skipped {len(failed_families)} families")
    
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
    parser.add_argument("--generate-pdfs", action="store_true", help="Generate PDF visualizations")
    parser.add_argument("--bin-conflict-analysis", action="store_true", help="Analyze BIN-species conflicts")
    parser.add_argument("--cleanup-intermediates", action="store_true", default=True, 
                       help="Clean up intermediate files to save space")
    
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
    
    # Run HPC-optimized pipeline
    pipeline = PhylogeneticPipeline(
        args.database, 
        args.output_dir,
        generate_pdfs=args.generate_pdfs,
        bin_conflict_analysis=args.bin_conflict_analysis
    )
    
    # Set cleanup option
    pipeline.cleanup_intermediates = args.cleanup_intermediates
    
    results = pipeline.run_pipeline_hpc_optimized(
        min_otus=args.min_otus, 
        max_families=args.max_families,
        alignment_method=args.alignment_method,
        tree_method=args.tree_method,
        bootstrap=args.bootstrap,
        num_outgroups=args.num_outgroups
    )
    
    if results:
        print(f"Phylogenetic analysis completed successfully!")
        print(f"Generated {len(results)} family phylogenies")
        if pipeline.generate_pdfs:
            pdf_count = len([r for r in results if r.get('pdf_file')])
            print(f"Generated {pdf_count} PDF visualizations")
        print(f"Output directory: {args.output_dir}")
        return 0
    else:
        print("No families were successfully processed.")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
