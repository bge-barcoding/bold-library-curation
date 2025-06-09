#!/usr/bin/env python3
"""
Phylogenetic Analysis Pipeline for BOLD Database
================================================

This script generates phylogenies for each family using species-BIN representatives,
with outgroups selected from other families within the same order.
Includes BAGS grade analysis, Grade C monophyly checking, and curation checklists.

Requirements:
- biopython
- pandas
- ete3 (for tree manipulation and PDF generation)
- muscle/mafft (alignment)
- iqtree/fasttree (phylogeny reconstruction)
- reportlab (for curation checklist PDFs)
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

# Import ReportLab for curation checklist PDFs
try:
    from reportlab.lib.pagesizes import letter, A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT, TA_CENTER
    REPORTLAB_AVAILABLE = True
except ImportError:
    print("Warning: ReportLab not available. Curation checklist PDFs will be disabled.")
    REPORTLAB_AVAILABLE = False

class PhylogeneticPipeline:
    def __init__(self, db_path, output_dir="phylogenies", generate_pdfs=False, bin_conflict_analysis=False, custom_parameters=None):
        """
        Initialize the phylogenetic pipeline.
        
        Args:
            db_path (str): Path to the BOLD SQLite database
            output_dir (str): Output directory for results
            generate_pdfs (bool): Whether to generate PDF visualizations
            bin_conflict_analysis (bool): Whether to analyze BIN-species conflicts (deprecated - now uses BAGS grades)
            custom_parameters (dict): Custom parameters for phylogenetic software
        """
        self.db_path = db_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.generate_pdfs = generate_pdfs
        self.bin_conflict_analysis = bin_conflict_analysis  # Kept for compatibility but not used
        self.custom_parameters = custom_parameters or {}
        
        # Print custom parameter info
        if self.custom_parameters:
            print(f"Custom parameters loaded:")
            for tool, params in self.custom_parameters.items():
                if params:  # Only show non-empty parameter lists
                    print(f"  {tool}: {params}")
        else:
            print("Using optimized defaults for phylogenetic software")
        
        # Curation checklist generation is always enabled if ReportLab is available
        self.generate_curation_checklists = REPORTLAB_AVAILABLE
        if not REPORTLAB_AVAILABLE:
            print("Warning: ReportLab not available. Curation checklists will be skipped.")
        
        # Import ETE3 only if PDF generation is requested
        print(f"DEBUG: PDF generation requested: {self.generate_pdfs}")
        if self.generate_pdfs:
            print("DEBUG: Starting ETE3 initialization...")
            try:
                # Set multiple environment variables for headless operation
                import os
                print("DEBUG: Setting environment variables...")
                os.environ['QT_QPA_PLATFORM'] = 'offscreen'
                os.environ['DISPLAY'] = ':99'  # Virtual display
                os.environ['MPLBACKEND'] = 'Agg'  # Force matplotlib non-GUI backend
                print("DEBUG: Environment variables set")
                
                # Import matplotlib first with headless backend
                print("DEBUG: Importing matplotlib...")
                import matplotlib
                matplotlib.use('Agg')
                print("DEBUG: Matplotlib imported successfully")
                
                print("DEBUG: Importing ETE3 components...")
                from ete3 import Tree, TreeStyle, TextFace, CircleFace
                print("DEBUG: ETE3 imports successful")
                
                self.Tree = Tree
                self.TreeStyle = TreeStyle
                self.TextFace = TextFace
                self.CircleFace = CircleFace
                print("✓ ETE3 loaded with headless backend")
            except ImportError as import_error:
                print(f"WARNING: ETE3 not available ({import_error}). PDF generation disabled.")
                self.generate_pdfs = False
            except Exception as e:
                print(f"WARNING: ETE3 backend issue ({type(e).__name__}: {e}). Continuing with PDF generation...")
                try:
                    # Try a simpler import without full initialization
                    print("DEBUG: Attempting fallback ETE3 import...")
                    from ete3 import Tree, TreeStyle, TextFace, CircleFace
                    self.Tree = Tree
                    self.TreeStyle = TreeStyle
                    self.TextFace = TextFace
                    self.CircleFace = CircleFace
                    print("✓ ETE3 loaded with basic backend (warnings ignored)")
                except Exception as fallback_error:
                    print(f"WARNING: ETE3 fallback failed ({fallback_error}). PDF generation disabled.")
                    self.generate_pdfs = False
        
        print(f"DEBUG: Final PDF generation status: {self.generate_pdfs}")
        
        # Performance settings
        self.cleanup_intermediates = True
        self.batch_size = 50  # Process families in batches for memory efficiency
    
    def get_families_to_process(self, min_otus=3):
        """
        Get families meeting criteria without loading all data.
        HPC-optimized for large datasets. Now counts distinct species-BIN combinations.
        
        Args:
            min_otus (int): Minimum number of species-BIN combinations required for analysis
            
        Returns:
            pandas.DataFrame: Families with their statistics
        """
        conn = sqlite3.connect(self.db_path)
        
        query = """
        SELECT 
            b.family,
            b.`order`,
            COUNT(DISTINCT b.species || '|' || b.bin_uri) as species_bin_count,
            COUNT(DISTINCT b.species) as species_count,
            COUNT(*) as total_records
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
        GROUP BY b.family, b.`order`
        HAVING species_bin_count >= ?
        ORDER BY species_bin_count DESC
        """
        
        families_df = pd.read_sql_query(query, conn, params=[min_otus])
        conn.close()
        
        print(f"Found {len(families_df)} families with ≥{min_otus} species-BIN combinations")
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
        
        # Single family query with proper indexing - now focusing on BIN_URI instead of OTU
        query = """
        SELECT DISTINCT
            b.nuc,
            b.genus,
            b.species,
            b.processid,
            b.recordid,
            b.bin_uri,
            b.`order`,
            b.family
        FROM bold b
        WHERE b.family = ?
            AND b.nuc IS NOT NULL 
            AND LENGTH(b.nuc) > 200
            AND b.genus IS NOT NULL
            AND b.species IS NOT NULL
            AND b.bin_uri IS NOT NULL
            AND b.bin_uri != ''
            AND b.bin_uri != 'None'
            AND b.bin_uri LIKE 'BOLD:%'
        ORDER BY b.species, b.bin_uri, LENGTH(b.nuc) DESC
        """
        
        df = pd.read_sql_query(query, conn, params=[family_name])
        conn.close()
        return df
    
    def select_representatives_optimized(self, family_df):
        """
        Select representative sequences (one per species-BIN combination).
        
        Args:
            family_df (pandas.DataFrame): Family sequence data
            
        Returns:
            pandas.DataFrame: Representative sequences
        """
        representatives = []
        
        # Group by species and BIN combination
        for (species, bin_uri), group in family_df.groupby(['species', 'bin_uri']):
            # Select the longest sequence for each species-BIN combination
            rep = group.loc[group['nuc'].str.len().idxmax()]
            representatives.append(rep)
        
        return pd.DataFrame(representatives)
    
    def collect_family_curation_data(self, family_name):
        """
        Collect all data needed for curation checklist organized by genus and BAGS grade.
        
        Args:
            family_name (str): Name of the family
            
        Returns:
            dict: Organized curation data by genus and grade
        """
        conn = sqlite3.connect(self.db_path)
        
        # Get all species in family with their BAGS grades, BINs, and genus information
        query = """
        SELECT DISTINCT
            b.genus,
            b.species,
            b.bin_uri,
            bags.bags_grade
        FROM bold b
        JOIN taxa t ON b.taxonid = t.taxonid
        LEFT JOIN bags ON t.taxonid = bags.taxonid
        WHERE b.family = ?
            AND b.genus IS NOT NULL
            AND b.species IS NOT NULL
            AND b.bin_uri IS NOT NULL
            AND b.bin_uri != ''
            AND b.bin_uri != 'None'
            AND b.bin_uri LIKE 'BOLD:%'
            AND t.level = 'species'
        ORDER BY b.genus, b.species, b.bin_uri
        """
        
        df = pd.read_sql_query(query, conn, params=[family_name])
        conn.close()
        
        if df.empty:
            return {'family': family_name, 'genera': {}}
        
        # Organize data by genus
        curation_data = {
            'family': family_name,
            'generated_at': datetime.now().isoformat(),
            'genera': {}
        }
        
        for genus in df['genus'].unique():
            genus_df = df[df['genus'] == genus]
            curation_data['genera'][genus] = self._organize_genus_by_grades(genus_df)
        
        return curation_data
    
    def _organize_genus_by_grades(self, genus_df):
        """
        Organize genus data by BAGS grades according to curation checklist format.
        
        Args:
            genus_df (pandas.DataFrame): Data for single genus
            
        Returns:
            dict: Organized data by grade categories
        """
        genus_data = {
            'grade_A_B_D': [],  # Grades A, B, D together
            'grade_C': [],      # Grade C separate
            'grade_E': [],      # Grade E separate 
            'grade_F': []       # Grade F separate
        }
        
        # Group by species to aggregate BINs
        for species in genus_df['species'].unique():
            species_df = genus_df[genus_df['species'] == species]
            
            # Get primary BAGS grade (should be same for all records of this species)
            bags_grade = species_df['bags_grade'].iloc[0] if not species_df['bags_grade'].isna().all() else 'Unknown'
            
            # Get all BINs for this species
            bins = sorted(species_df['bin_uri'].unique())
            
            species_data = {
                'species': species,
                'bins': bins,
                'grade': bags_grade
            }
            
            # Categorize by grade
            if bags_grade in ['A', 'B', 'D']:
                genus_data['grade_A_B_D'].append(species_data)
            elif bags_grade == 'C':
                genus_data['grade_C'].append(species_data)
            elif bags_grade == 'E':
                genus_data['grade_E'].append(species_data)
            elif bags_grade == 'F':
                genus_data['grade_F'].append(species_data)
        
        # For Grade E, we need to collect BIN sharing information
        if genus_data['grade_E']:
            genus_data['grade_E'] = self._add_bin_sharing_info(genus_data['grade_E'])
        
        return genus_data
    
    def _add_bin_sharing_info(self, grade_e_species):
        """
        Add BIN sharing information for Grade E species.
        
        Args:
            grade_e_species (list): List of Grade E species data
            
        Returns:
            list: Enhanced species data with sharing information
        """
        conn = sqlite3.connect(self.db_path)
        
        # For each species, check what other species share their BINs
        for species_data in grade_e_species:
            species_name = species_data['species']
            species_bins = species_data['bins']
            
            sharing_info = []
            
            for bin_uri in species_bins:
                # Find all other species that share this BIN
                query = """
                SELECT DISTINCT b.species
                FROM bold b
                JOIN taxa t ON b.taxonid = t.taxonid
                WHERE b.bin_uri = ?
                    AND b.species != ?
                    AND b.species IS NOT NULL
                    AND t.level = 'species'
                ORDER BY b.species
                """
                
                sharing_df = pd.read_sql_query(query, conn, params=[bin_uri, species_name])
                sharing_species = sharing_df['species'].tolist() if not sharing_df.empty else []
                
                sharing_info.append({
                    'bin_uri': bin_uri,
                    'sharing_species': sharing_species if sharing_species else None
                })
            
            species_data['sharing_info'] = sharing_info
        
        conn.close()
        return grade_e_species

    def get_bags_grades(self, family_df):
        """
        Get BAGS grades for sequences in family from the database.
        
        Args:
            family_df (pandas.DataFrame): Family sequence data
            
        Returns:
            dict: Mapping of species to BAGS grades
        """
        if family_df.empty:
            return {}
            
        conn = sqlite3.connect(self.db_path)
        
        # Get unique species in the family
        species_list = family_df['species'].unique().tolist()
        placeholders = ','.join(['?' for _ in species_list])
        
        query = f"""
        SELECT t.name as species, bags.bags_grade
        FROM taxa t
        JOIN bags ON t.taxonid = bags.taxonid
        WHERE t.level = 'species' 
        AND t.name IN ({placeholders})
        """
        
        grades_df = pd.read_sql_query(query, conn, params=species_list)
        conn.close()
        
        # Convert to dictionary for easy lookup
        grades_dict = dict(zip(grades_df['species'], grades_df['bags_grade']))
        return grades_dict

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
            b.nuc,
            b.genus,
            b.species,
            b.processid,
            b.recordid,
            b.bin_uri,
            b.family
        FROM bold b
        WHERE b.`order` = ?
            AND b.family != ?
            AND b.family IS NOT NULL
            AND b.nuc IS NOT NULL 
            AND LENGTH(b.nuc) > 200
            AND b.genus IS NOT NULL
            AND b.species IS NOT NULL
            AND b.bin_uri IS NOT NULL
            AND b.bin_uri != ''
            AND b.bin_uri != 'None'
            AND b.bin_uri LIKE 'BOLD:%'
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
    
    def parse_sequence_metadata(self, tree, bags_grades=None):
        """
        Extract metadata from sequence IDs including BAGS grades.
        
        Args:
            tree: ETE3 Tree object
            bags_grades (dict): Mapping of species to BAGS grades
            
        Returns:
            dict: Metadata for each sequence
        """
        metadata = {}
        bags_grades = bags_grades or {}
        
        for leaf in tree:
            if leaf.name.startswith('OUTGROUP-'):
                # Parse: OUTGROUP-species-bin_uri-processid
                parts = leaf.name.split('-', 3)  # Split on first 3 hyphens
                is_outgroup = True
                species = parts[1] if len(parts) > 1 else 'Unknown'
                bin_uri = parts[2] if len(parts) > 2 else 'NA'
                processid = parts[3] if len(parts) > 3 else 'NA'
            else:
                # Parse: species-bin_uri-processid
                parts = leaf.name.split('-', 2)  # Split on first 2 hyphens
                is_outgroup = False
                species = parts[0] if len(parts) > 0 else 'Unknown'
                bin_uri = parts[1] if len(parts) > 1 else 'NA'
                processid = parts[2] if len(parts) > 2 else 'NA'
            
            # Get BAGS grade for the species
            bags_grade = bags_grades.get(species.replace('_', ' '), 'Unknown') if not is_outgroup else 'Outgroup'
            
            metadata[leaf.name] = {
                'species': species,
                'bin_uri': bin_uri,
                'processid': processid,
                'is_outgroup': is_outgroup,
                'bags_grade': bags_grade
            }
        
        return metadata
    
    def generate_curation_checklist_pdf(self, curation_data, pdf_file):
        """
        Generate curation checklist PDF with proper formatting and checkboxes.
        
        Args:
            curation_data (dict): Organized curation data
            pdf_file (Path): Output PDF file path
            
        Returns:
            bool: Success status
        """
        if not REPORTLAB_AVAILABLE:
            print("ReportLab not available - skipping curation checklist PDF")
            return False
            
        try:
            # Create PDF document
            doc = SimpleDocTemplate(
                str(pdf_file),
                pagesize=A4,
                rightMargin=0.5*inch,
                leftMargin=0.5*inch,
                topMargin=0.75*inch,
                bottomMargin=0.75*inch
            )
            
            # Define styles
            styles = getSampleStyleSheet()
            title_style = ParagraphStyle(
                'CustomTitle',
                parent=styles['Heading1'],
                fontSize=16,
                spaceAfter=30,
                alignment=TA_CENTER
            )
            
            genus_style = ParagraphStyle(
                'GenusTitle',
                parent=styles['Heading2'],
                fontSize=14,
                spaceAfter=12,
                spaceBefore=20
            )
            
            grade_style = ParagraphStyle(
                'GradeTitle',
                parent=styles['Heading3'],
                fontSize=12,
                spaceAfter=8,
                spaceBefore=16
            )
            
            # Build document content
            story = []
            
            # Title
            family_name = curation_data['family']
            title = f"PDF CHECKLIST FOR CURATION - FAMILY {family_name.upper()}"
            story.append(Paragraph(title, title_style))
            
            # Process each genus
            genera = curation_data['genera']
            for genus_idx, (genus, genus_data) in enumerate(genera.items()):
                # Add page break if not first genus
                if genus_idx > 0:
                    story.append(PageBreak())
                
                # Genus heading
                story.append(Paragraph(f"<b>{genus}:</b>", genus_style))
                
                # Process each grade category
                self._add_grade_a_b_d_section(story, genus_data['grade_A_B_D'], grade_style)
                self._add_grade_c_section(story, genus_data['grade_C'], grade_style)
                self._add_grade_e_section(story, genus_data['grade_E'], grade_style)
                self._add_grade_f_section(story, genus_data['grade_F'], grade_style)
            
            # Build PDF
            doc.build(story)
            return True
            
        except Exception as e:
            print(f"Error generating curation checklist PDF: {e}")
            return False
    
    def _add_grade_a_b_d_section(self, story, species_list, grade_style):
        """Add BAGS Grade A/B/D section to PDF"""
        if not species_list:
            return
            
        # Section header
        header_text = """<b>BAGS Grade A / B / D:</b><br/>
        These species don't have any taxonomic conflict (1 species in 1 BIN), but can be checked to assess validity."""
        story.append(Paragraph(header_text, grade_style))
        story.append(Spacer(1, 6))
        
        # Create table data
        table_data = [['Species', 'BINs', 'Curated (tick)']]
        
        for species_data in species_list:
            bins_text = '; '.join(species_data['bins'])
            table_data.append([
                species_data['species'],
                bins_text,
                '☐'  # Empty checkbox
            ])
        
        # Create and style table
        table = Table(table_data, colWidths=[3*inch, 2.5*inch, 1*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
        ]))
        
        story.append(table)
        story.append(Spacer(1, 12))
    
    def _add_grade_c_section(self, story, species_list, grade_style):
        """Add BAGS Grade C section to PDF"""
        if not species_list:
            return
            
        # Section header
        header_text = """<b>BAGS Grade C:</b><br/>
        These species don't have any taxonomic conflict, but have multiple BINs so should be checked to assess validity of each BIN."""
        story.append(Paragraph(header_text, grade_style))
        story.append(Spacer(1, 6))
        
        # Create table data with each BIN on its own row
        table_data = [['Species', 'BIN', 'Curated (tick)']]
        
        for species_data in species_list:
            species_name = species_data['species']
            bins = species_data['bins']
            
            for i, bin_uri in enumerate(bins):
                if i == 0:
                    # First row includes species name and checkbox
                    species_cell = species_name
                    checkbox_cell = '☐'
                else:
                    # Subsequent rows have empty species name and checkbox
                    species_cell = ''
                    checkbox_cell = ''
                
                table_data.append([
                    species_cell,
                    bin_uri,
                    checkbox_cell
                ])
        
        # Create and style table
        table = Table(table_data, colWidths=[2.5*inch, 2.5*inch, 1*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
        ]))
        
        story.append(table)
        story.append(Spacer(1, 12))
    
    def _add_grade_e_section(self, story, species_list, grade_style):
        """Add BAGS Grade E section to PDF"""
        if not species_list:
            return
            
        # Section header
        header_text = """<b>BAGS Grade E:</b><br/>
        These species show taxonomic conflict, they must be checked to assess validity of these conflicts vs mis-identifications. 
        Note some species do share BINs so the aim is to resolve misidentifications rather than always make 1 species = 1 BIN."""
        story.append(Paragraph(header_text, grade_style))
        story.append(Spacer(1, 6))
        
        # Create table data with species-centric view
        table_data = [['Species', 'BIN URI', 'Sharing Species', 'Curated (tick)']]
        
        for species_data in species_list:
            species_name = species_data['species']
            sharing_info = species_data.get('sharing_info', [])
            
            for i, bin_info in enumerate(sharing_info):
                if i == 0:
                    # First row includes species name
                    species_cell = species_name
                    checkbox_cell = '☐'
                else:
                    # Subsequent rows have empty species name
                    species_cell = ''
                    checkbox_cell = ''
                
                bin_uri = bin_info['bin_uri']
                sharing_species = bin_info['sharing_species']
                
                if sharing_species:
                    sharing_text = '<br/>'.join(sharing_species)
                else:
                    sharing_text = 'N/A'
                
                table_data.append([
                    species_cell,
                    bin_uri,
                    Paragraph(sharing_text, getSampleStyleSheet()['Normal']),
                    checkbox_cell
                ])
        
        # Create and style table
        table = Table(table_data, colWidths=[2*inch, 1.5*inch, 2.5*inch, 0.75*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
        ]))
        
        story.append(table)
        story.append(Spacer(1, 12))
    
    def _add_grade_f_section(self, story, species_list, grade_style):
        """Add BAGS Grade F section to PDF"""
        if not species_list:
            return
            
        # Section header
        header_text = """<b>BAGS Grade F:</b><br/>
        These species have no valid BIN assignments and require sequence quality assessment."""
        story.append(Paragraph(header_text, grade_style))
        story.append(Spacer(1, 6))
        
        # Create table data
        table_data = [['Species', 'Status', 'Curated (tick)']]
        
        for species_data in species_list:
            table_data.append([
                species_data['species'],
                'No valid BINs',
                '☐'  # Empty checkbox
            ])
        
        # Create and style table
        table = Table(table_data, colWidths=[3*inch, 2.5*inch, 1*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 10),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
        ]))
        
        story.append(table)
        story.append(Spacer(1, 12))

    def analyze_grade_c_monophyly(self, tree, metadata):
        """
        Check monophyly of Grade C species using strict criteria.
        
        Args:
            tree: ETE3 Tree object
            metadata (dict): Sequence metadata including BAGS grades
            
        Returns:
            dict: Boolean monophyly results for each Grade C species
        """
        results = {}
        
        # Group sequences by species, focusing on Grade C
        grade_c_species = {}
        for seq_id, data in metadata.items():
            if not data['is_outgroup'] and data['bags_grade'] == 'C':
                species = data['species']
                if species not in grade_c_species:
                    grade_c_species[species] = []
                grade_c_species[species].append(seq_id)
        
        # Check monophyly for each Grade C species with ≥2 sequences
        for species, seq_ids in grade_c_species.items():
            if len(seq_ids) >= 2:  # Minimum requirement for meaningful test
                try:
                    # Find corresponding tree leaves
                    leaves = []
                    for seq_id in seq_ids:
                        leaf_nodes = tree.search_nodes(name=seq_id)
                        if leaf_nodes:
                            leaves.append(leaf_nodes[0])
                    
                    if len(leaves) >= 2:
                        # Check strict monophyly - no exceptions allowed
                        is_monophyletic = tree.check_monophyly(
                            values=[leaf.name for leaf in leaves], 
                            target_attr="name"
                        )[0]
                        
                        results[species] = {
                            'is_monophyletic': is_monophyletic,
                            'num_sequences': len(leaves),
                            'testable': True
                        }
                    else:
                        # Couldn't find enough leaves in tree
                        results[species] = {
                            'is_monophyletic': None,
                            'num_sequences': len(leaves),
                            'testable': False
                        }
                        
                except Exception as e:
                    print(f"Error checking monophyly for {species}: {e}")
                    results[species] = {
                        'is_monophyletic': None,
                        'num_sequences': len(seq_ids),
                        'testable': False
                    }
            else:
                # Insufficient sequences for monophyly test
                results[species] = {
                    'is_monophyletic': None,
                    'num_sequences': len(seq_ids),
                    'testable': False
                }
        
        return results

    def analyze_bags_grades(self, metadata):
        """
        Analyze BAGS grades distribution instead of BIN-species conflicts.
        
        Args:
            metadata (dict): Sequence metadata including BAGS grades
            
        Returns:
            dict: Analysis of BAGS grades
        """
        # Count BAGS grades (excluding outgroups)
        grade_counts = {}
        species_grades = {}
        
        for seq_id, data in metadata.items():
            if not data['is_outgroup']:
                grade = data['bags_grade']
                species = data['species']
                
                # Count grades
                grade_counts[grade] = grade_counts.get(grade, 0) + 1
                
                # Track species grades (one per species)
                species_grades[species] = grade
        
        # Count unique species per grade
        unique_species_per_grade = {}
        for grade in grade_counts.keys():
            unique_species_per_grade[grade] = len([s for s, g in species_grades.items() if g == grade])
        
        return {
            'grade_counts': grade_counts,
            'species_grades': species_grades,
            'unique_species_per_grade': unique_species_per_grade,
            'total_species': len(species_grades)
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
            'curation_pdf_file': family_dir / f"{family_name}_curation_checklist.pdf",
            'metadata_file': family_dir / f"tree_metadata.json",
            'curation_data_file': family_dir / f"curation_data.json"
        }
    def create_fasta_file(self, representatives_df, outgroups_df, fasta_path):
        """
        Create a FASTA file for a family with ingroup and outgroup sequences.
        Uses the branch labeling format: species-bin_uri-processid
        
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
                bin_uri = str(row.get('bin_uri', 'NA')).replace(' ', '_')
                processid = str(row.get('processid', 'NA')).replace(' ', '_')
                
                # Create sequence ID in format: species-bin_uri-processid
                seq_id = f"{species}-{bin_uri}-{processid}"
                
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
                bin_uri = str(row.get('bin_uri', 'NA')).replace(' ', '_')
                processid = str(row.get('processid', 'NA')).replace(' ', '_')
                
                # Create sequence ID with OUTGROUP prefix
                seq_id = f"OUTGROUP-{species}-{bin_uri}-{processid}"
                
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
        Generate PDF with BAGS grade visualization and Grade C monophyly analysis.
        
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
            
            # Get BAGS grades for this family
            # We need to get the family data to look up BAGS grades
            family_df = self.get_family_sequences_optimized(family_name)
            bags_grades = self.get_bags_grades(family_df)
            
            # Parse metadata and analyze BAGS grades
            metadata = self.parse_sequence_metadata(tree, bags_grades)
            bags_analysis = self.analyze_bags_grades(metadata)
            
            # Perform Grade C monophyly analysis
            monophyly_results = self.analyze_grade_c_monophyly(tree, metadata)
            
            # Style tree for BAGS grades including monophyly results
            self._style_tree_for_bags_grades(tree, metadata, bags_analysis, monophyly_results)
            
            # Generate PDF with monophyly information
            success = self._generate_annotated_pdf(tree, family_name, bags_analysis, pdf_file, monophyly_results)
            
            # Save metadata including monophyly results
            if success:
                self._save_tree_metadata(metadata, bags_analysis, family_name, metadata_file, monophyly_results)
                
                # Print summary of monophyly results for this family
                if monophyly_results:
                    mono_count = sum(1 for r in monophyly_results.values() 
                                   if r['testable'] and r['is_monophyletic'])
                    non_mono_count = sum(1 for r in monophyly_results.values() 
                                       if r['testable'] and not r['is_monophyletic'])
                    if non_mono_count > 0:
                        print(f"   ⚠ Found {non_mono_count} non-monophyletic Grade C species in {family_name}")
                    if mono_count > 0:
                        print(f"   ✓ Found {mono_count} monophyletic Grade C species in {family_name}")
            
            return success
            
        except Exception as e:
            print(f"Error generating PDF for {family_name}: {e}")
            return False
    
    def _style_tree_for_bags_grades(self, tree, metadata, bags_analysis, monophyly_results=None):
        """Apply visual styling based on BAGS grades and Grade C monophyly"""
        
        # Define color schemes for BAGS grades with Grade C monophyly distinction
        colors = {
            'A': '#2E7D32',           # Dark green - high quality (≥10 specimens, 1 unshared BIN)
            'B': '#388E3C',           # Green - good quality (3-9 specimens, 1 unshared BIN)
            'C_monophyletic': '#FFA000',    # Orange - BIN splitting but monophyletic
            'C_non_monophyletic': '#D32F2F', # Red - BIN splitting and non-monophyletic
            'C_untestable': '#FF8F00',      # Darker orange - Grade C but insufficient data for test
            'D': '#1976D2',           # Blue - low specimen count (<3 specimens, 1 unshared BIN)
            'E': '#7B1FA2',           # Purple - BIN sharing (shared with other species)
            'F': '#424242',           # Grey - no valid BINs
            'Unknown': '#9E9E9E',     # Light grey - unknown grade
            'Outgroup': '#9E9E9E'     # Light grey - outgroups
        }
        
        grade_descriptions = {
            'A': 'Grade A: ≥10 specimens, 1 unshared BIN',
            'B': 'Grade B: 3-9 specimens, 1 unshared BIN', 
            'C_monophyletic': 'Grade C: >1 unshared BIN, monophyletic',
            'C_non_monophyletic': 'Grade C: >1 unshared BIN, NON-MONOPHYLETIC',
            'C_untestable': 'Grade C: >1 unshared BIN, insufficient data for monophyly test',
            'D': 'Grade D: <3 specimens, 1 unshared BIN',
            'E': 'Grade E: BIN sharing with other species',
            'F': 'Grade F: No valid BINs',
            'Unknown': 'Unknown BAGS grade',
            'Outgroup': 'Outgroup sequence'
        }
        
        monophyly_results = monophyly_results or {}
        
        for leaf in tree:
            data = metadata[leaf.name]
            base_grade = data['bags_grade']
            
            # Special handling for Grade C with monophyly information
            if base_grade == 'C' and not data['is_outgroup']:
                species = data['species'].replace('_', ' ')  # Convert back to proper species name
                
                if species in monophyly_results:
                    mono_data = monophyly_results[species]
                    if mono_data['testable']:
                        if mono_data['is_monophyletic']:
                            grade_key = 'C_monophyletic'
                        else:
                            grade_key = 'C_non_monophyletic'
                    else:
                        grade_key = 'C_untestable'
                else:
                    grade_key = 'C_untestable'
            else:
                grade_key = base_grade
            
            leaf_color = colors.get(grade_key, colors['Unknown'])
            description = grade_descriptions.get(grade_key, 'Unknown grade')
            
            # Style the leaf
            leaf.add_features(
                bgcolor=leaf_color,
                bags_grade=base_grade,
                grade_key=grade_key,
                grade_description=description,
                species_clean=data['species'],
                bin_clean=data['bin_uri']
            )
    
    def _generate_annotated_pdf(self, tree, family_name, bags_analysis, pdf_file, monophyly_results=None):
        """Generate PDF with BAGS grade annotations and Grade C monophyly results"""
        try:
            # Create tree style
            ts = self.TreeStyle()
            ts.show_leaf_name = False  # We'll add custom labels
            ts.show_branch_support = True
            ts.title.add_face(self.TextFace(f"Family: {family_name}", fsize=16), column=0)
            
            # Add BAGS grade summary
            grade_summary = []
            for grade, count in sorted(bags_analysis['unique_species_per_grade'].items()):
                if grade == 'C' and monophyly_results:
                    # Break down Grade C by monophyly
                    mono_count = sum(1 for r in monophyly_results.values() 
                                   if r['testable'] and r['is_monophyletic'])
                    non_mono_count = sum(1 for r in monophyly_results.values() 
                                       if r['testable'] and not r['is_monophyletic'])
                    untestable_count = sum(1 for r in monophyly_results.values() 
                                         if not r['testable'])
                    
                    if mono_count > 0:
                        grade_summary.append(f"Grade C (monophyletic): {mono_count}")
                    if non_mono_count > 0:
                        grade_summary.append(f"Grade C (NON-monophyletic): {non_mono_count}")
                    if untestable_count > 0:
                        grade_summary.append(f"Grade C (untestable): {untestable_count}")
                else:
                    grade_summary.append(f"Grade {grade}: {count}")
            
            summary_text = f"""BAGS Grades Summary:
• {' • '.join(grade_summary)}
• Total species: {bags_analysis['total_species']}"""
            
            if monophyly_results:
                testable_grade_c = sum(1 for r in monophyly_results.values() if r['testable'])
                if testable_grade_c > 0:
                    summary_text += f"\n• Grade C species tested for monophyly: {testable_grade_c}"
            
            ts.title.add_face(self.TextFace(summary_text, fsize=12), column=0)
            
            # Custom node styling function
            def layout(node):
                if node.is_leaf():
                    # Add colored circle based on BAGS grade (including monophyly)
                    circle = self.CircleFace(radius=8, color=node.bgcolor, style="sphere")
                    node.add_face(circle, column=0, position="branch-right")
                    
                    # Add species name
                    species_face = self.TextFace(f" {node.species_clean}", fsize=10)
                    node.add_face(species_face, column=1, position="branch-right")
                    
                    # Add BIN info
                    bin_face = self.TextFace(f" [{node.bin_clean}]", fsize=8, fgcolor="gray")
                    node.add_face(bin_face, column=2, position="branch-right")
                    
                    # Add BAGS grade with monophyly info
                    if hasattr(node, 'bags_grade') and node.bags_grade not in ['Unknown', 'Outgroup']:
                        grade_text = node.bags_grade
                        if hasattr(node, 'grade_key') and 'C_non_monophyletic' in node.grade_key:
                            grade_text += " ⚠"  # Warning symbol for non-monophyletic
                        
                        grade_face = self.TextFace(f" {grade_text}", fsize=8, fgcolor="black")
                        node.add_face(grade_face, column=3, position="branch-right")
            
            ts.layout_fn = layout
            
            # Generate PDF
            tree.render(str(pdf_file), tree_style=ts, dpi=300)
            return True
            
        except Exception as e:
            print(f"Error rendering PDF: {e}")
            return False
    
    def _save_tree_metadata(self, metadata, bags_analysis, family_name, metadata_file, monophyly_results=None):
        """Save tree metadata with BAGS analysis and Grade C monophyly results to JSON file"""
        try:
            tree_metadata = {
                'family': family_name,
                'generated_at': datetime.now().isoformat(),
                'total_sequences': len(metadata),
                'ingroup_sequences': len([m for m in metadata.values() if not m['is_outgroup']]),
                'outgroup_sequences': len([m for m in metadata.values() if m['is_outgroup']]),
                'bags_analysis': {
                    'grade_counts': bags_analysis['grade_counts'],
                    'unique_species_per_grade': bags_analysis['unique_species_per_grade'],
                    'total_species': bags_analysis['total_species'],
                    'species_grades': bags_analysis['species_grades']
                }
            }
            
            # Add Grade C monophyly results if available
            if monophyly_results:
                monophyly_summary = {
                    'grade_c_monophyly': {},
                    'monophyletic_count': 0,
                    'non_monophyletic_count': 0,
                    'untestable_count': 0
                }
                
                for species, result in monophyly_results.items():
                    # Store boolean results only as requested
                    monophyly_summary['grade_c_monophyly'][species] = result['is_monophyletic']
                    
                    if result['testable']:
                        if result['is_monophyletic']:
                            monophyly_summary['monophyletic_count'] += 1
                        else:
                            monophyly_summary['non_monophyletic_count'] += 1
                    else:
                        monophyly_summary['untestable_count'] += 1
                
                tree_metadata['monophyly_analysis'] = monophyly_summary
            
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
            
            # 3. Check if we have enough species-BIN combinations
            unique_species_bins = representatives.groupby(['species', 'bin_uri']).ngroups
            if unique_species_bins < min_otus:
                print(f"   Insufficient species-BIN combinations: {unique_species_bins} < {min_otus}")
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
            
            # 9. Generate curation checklist PDF
            print(f"   Generating curation checklist...")
            curation_data = self.collect_family_curation_data(family_name)
            
            # Save curation data as JSON
            try:
                with open(paths['curation_data_file'], 'w') as f:
                    json.dump(curation_data, f, indent=2)
            except Exception as e:
                print(f"   Warning: Failed to save curation data: {e}")
            
            # Generate curation checklist PDF
            if not self.generate_curation_checklist_pdf(curation_data, paths['curation_pdf_file']):
                print(f"   Warning: Curation checklist PDF generation failed for {family_name}")
            else:
                print(f"   ✓ Curation checklist PDF generated: {paths['curation_pdf_file'].name}")
            
            # 10. Clean up intermediate files to save space
            if self.cleanup_intermediates:
                paths['fasta_file'].unlink(missing_ok=True)
                paths['alignment_file'].unlink(missing_ok=True)
            
            return {
                'family': family_name,
                'order': family_order,
                'num_species_bins': unique_species_bins,
                'num_species': representatives['species'].nunique(),
                'num_sequences': len(representatives),
                'num_outgroups': len(outgroups),
                'tree_file': paths['tree_file'],
                'pdf_file': paths['pdf_file'] if self.generate_pdfs else None,
                'curation_pdf_file': paths['curation_pdf_file'],
                'curation_data_file': paths['curation_data_file'],
                'total_family_species': len(curation_data['genera']) if curation_data else 0
            }
            
        except Exception as e:
            print(f"   Error processing {family_name}: {e}")
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
                                   tree_method="iqtree", bootstrap=1000, num_outgroups=3, families_to_process=None):
        """
        Run the HPC-optimized phylogenetic analysis pipeline.
        
        Args:
            min_otus (int): Minimum number of OTUs required for analysis
            max_families (int): Maximum number of families to process (for testing)
            alignment_method (str): Alignment method ("mafft" or "muscle")
            tree_method (str): Tree building method ("iqtree" or "fasttree")
            bootstrap (int): Bootstrap replicates
            num_outgroups (int): Number of outgroup sequences
            families_to_process (list): Pre-defined list of families to process (for batch mode)
        """
        batch_prefix = f"[Batch {getattr(self, 'batch_id', '?')}] " if hasattr(self, 'batch_id') and self.batch_id else ""
        
        print(f"{batch_prefix}Starting HPC-optimized phylogenetic analysis pipeline...")
        print(f"Database: {self.db_path}")
        print(f"Output directory: {self.output_dir}")
        print(f"Generate PDFs: {self.generate_pdfs}")
        print(f"BIN conflict analysis: {self.bin_conflict_analysis}")
        print("="*60)
        
        if families_to_process:
            # Batch processing mode - use provided family list
            print(f"{batch_prefix}Batch processing mode: {len(families_to_process)} families provided")
            
            # Verify families exist and get their statistics
            verified_families = []
            for family_info in families_to_process:
                family_name = family_info['family']
                
                # Skip if already completed (checkpoint)
                if hasattr(self, 'completed_families') and family_name in self.completed_families:
                    print(f"   Skipping {family_name} (already completed)")
                    continue
                
                # Verify family exists in database and get stats
                family_df = self.get_family_sequences_optimized(family_name)
                if len(family_df) == 0:
                    print(f"   Warning: No sequences found for {family_name}")
                    continue
                
                # Check species-BIN count
                representatives = self.select_representatives_optimized(family_df)
                species_bin_count = representatives.groupby(['species', 'bin_uri']).ngroups
                
                if species_bin_count < min_otus:
                    print(f"   Skipping {family_name}: insufficient species-BIN combinations ({species_bin_count} < {min_otus})")
                    continue
                
                # Get order information
                family_order = family_df['order'].iloc[0] if not family_df.empty else family_info.get('order')
                
                verified_families.append({
                    'family': family_name,
                    'order': family_order,
                    'species_bin_count': species_bin_count,
                    'species_count': representatives['species'].nunique()
                })
            
            families_df = pd.DataFrame(verified_families)
            print(f"{batch_prefix}Verified {len(families_df)} families for processing")
            
        else:
            # Discovery mode - find all suitable families
            print(f"{batch_prefix}Discovery mode: scanning database for suitable families...")
            families_df = self.get_families_to_process(min_otus)
            
            if max_families and max_families > 0:
                families_df = families_df.head(max_families)
                print(f"   Limited to first {max_families} families for testing")
        
        print(f"Processing {len(families_df)} families:")
        for _, row in families_df.iterrows():
            print(f"  - {row['family']}: {row['species_bin_count']} species-BIN combinations, {row['species_count']} species")
        
        print("="*60)
        
        results = []
        failed_families = []
        
        # Process families with checkpointing
        for i, (_, family_row) in enumerate(families_df.iterrows(), 1):
            family_name = family_row['family']
            family_order = family_row['order']
            
            print(f"\n{batch_prefix}[{i}/{len(families_df)}] Processing {family_name}")
            
            try:
                result = self.process_single_family(
                    family_name, family_order, min_otus, num_outgroups,
                    alignment_method, tree_method, bootstrap
                )
                
                if result:
                    results.append(result)
                    print(f"   ✓ Success: {family_name}")
                    
                    # Save checkpoint
                    self._save_checkpoint(family_name)
                else:
                    failed_families.append(family_name)
                    print(f"   ⚠ Skipped: {family_name}")
                    
            except Exception as e:
                failed_families.append(family_name)
                print(f"   ✗ Error: {family_name} - {e}")
                
            # Memory cleanup every 5 families
            if i % 5 == 0:
                gc.collect()
                print(f"   Memory cleanup after {i} families")
        
        # Generate summary
        self.generate_pipeline_summary(results, failed_families, batch_prefix)
        self.generate_summary_file(results, failed_families, batch_prefix)
        return results
    
    def _save_checkpoint(self, completed_family):
        """Save checkpoint after completing a family"""
        if not hasattr(self, 'checkpoint_file') or not self.checkpoint_file:
            return
            
        try:
            # Load existing checkpoint or create new one
            checkpoint_data = {}
            if Path(self.checkpoint_file).exists():
                with open(self.checkpoint_file, 'r') as f:
                    checkpoint_data = json.load(f)
            
            # Add completed family
            completed_families = set(checkpoint_data.get('completed_families', []))
            completed_families.add(completed_family)
            
            # Update checkpoint
            checkpoint_data.update({
                'completed_families': list(completed_families),
                'last_updated': datetime.now().isoformat(),
                'batch_id': getattr(self, 'batch_id', None)
            })
            
            # Save checkpoint
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
                
        except Exception as e:
            print(f"   Warning: Could not save checkpoint: {e}")
    
    def generate_pipeline_summary(self, results, failed_families, batch_prefix=""):
        """Generate comprehensive summary of pipeline results"""
        print("\n" + "="*60)
        print(f"{batch_prefix}PHYLOGENETIC ANALYSIS COMPLETE")
        print("="*60)
        
        total_attempted = len(results) + len(failed_families)
        success_rate = (len(results) / total_attempted * 100) if total_attempted > 0 else 0
        
        print(f"Families attempted: {total_attempted}")
        print(f"Successfully processed: {len(results)} ({success_rate:.1f}%)")
        print(f"Failed/skipped: {len(failed_families)}")
        
        if results:
            print(f"\nSuccessfully processed families:")
            total_species_bins = sum(r['num_species_bins'] for r in results)
            total_species = sum(r['num_species'] for r in results)
            total_sequences = sum(r['num_sequences'] for r in results)
            
            print(f"  Total species-BIN combinations: {total_species_bins:,}")
            print(f"  Total unique species: {total_species:,}")
            print(f"  Total sequences processed: {total_sequences:,}")
            
            for result in results[:10]:  # Show first 10
                print(f"  - {result['family']}: {result['num_species_bins']} species-BINs, "
                      f"{result['num_species']} species, {result['num_sequences']} sequences")
            
            if len(results) > 10:
                print(f"  ... and {len(results) - 10} more families")
        
        if failed_families:
            print(f"\nFailed/skipped families:")
            for family in failed_families[:20]:  # Show first 20
                print(f"  - {family}")
            if len(failed_families) > 20:
                print(f"  ... and {len(failed_families) - 20} more")
        
        print("="*60)
    
    def align_sequences_optimized(self, fasta_file, alignment_file, method="mafft"):
        """Optimized sequence alignment with custom parameter support"""
        try:
            if method.lower() == "mafft":
                # Optimized defaults for COI sequences - high gap penalties to prevent spurious indels
                default_params = ["--quiet", "--adjustdirection", "--op", "10.0", "--ep", "1.0", "--maxiterate", "1000", "--localpair"]
                
                # Check for custom parameters from config
                custom_params = getattr(self, 'custom_parameters', {}).get('mafft', None)
                if custom_params:
                    # Use custom parameters but ensure essential ones are present
                    cmd = ["mafft"] + custom_params + [str(fasta_file)]
                    print(f"   Using custom MAFFT parameters: {custom_params}")
                else:
                    # Use optimized defaults
                    cmd = ["mafft"] + default_params + [str(fasta_file)]
                    print(f"   Using optimized MAFFT defaults for COI sequences")
                
                with open(alignment_file, 'w') as output_file:
                    result = subprocess.run(cmd, check=True, stdout=output_file, 
                                          stderr=subprocess.PIPE, text=True)
                return True
                
            elif method.lower() == "muscle":
                # MUSCLE with optimized defaults (no custom params for now)
                # Try MUSCLE 5.x syntax first
                cmd = ["muscle", "-align", str(fasta_file), "-output", str(alignment_file)]
                try:
                    subprocess.run(cmd, check=True, capture_output=True, text=True)
                    print(f"   Using MUSCLE 5.x syntax")
                    return True
                except subprocess.CalledProcessError:
                    # Fallback to legacy syntax
                    cmd = ["muscle", "-in", str(fasta_file), "-out", str(alignment_file)]
                    subprocess.run(cmd, check=True, capture_output=True, text=True)
                    print(f"   Using MUSCLE legacy syntax")
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
                
                # Optimized defaults for COI sequences - fast mode with appropriate model
                default_params = ["--fast", "-m", "GTR+I+G"]
                
                # Check for custom parameters from config
                custom_params = getattr(self, 'custom_parameters', {}).get('iqtree', None)
                if custom_params:
                    # Use custom parameters
                    cmd = ["iqtree", "-s", str(alignment_file), "-pre", str(tree_prefix), "-nt", "AUTO"] + custom_params
                    print(f"   Using custom IQ-TREE parameters: {custom_params}")
                else:
                    # Use optimized defaults
                    cmd = ["iqtree", "-s", str(alignment_file), "-pre", str(tree_prefix), "-nt", "AUTO"] + default_params
                    print(f"   Using optimized IQ-TREE defaults (GTR+I+G model, fast mode)")
                
                # Add bootstrap only if >= 1000 (IQ-TREE requirement) and not overridden by custom params
                if bootstrap >= 1000 and not any('-b' in str(p) for p in (custom_params or [])):
                    cmd.extend(["-bb", str(bootstrap)])
                    print(f"   Using {bootstrap} ultrafast bootstrap replicates")
                elif bootstrap > 0 and not any('-b' in str(p) for p in (custom_params or [])):
                    print(f"   Warning: Bootstrap {bootstrap} < 1000, using standard bootstrap instead")
                    cmd.extend(["-b", str(bootstrap)])
                else:
                    print(f"   Bootstrap handled by custom parameters or disabled")
                
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
                # Optimized defaults for FastTree
                default_params = ["-gtr", "-gamma"]
                
                # Check for custom parameters from config
                custom_params = getattr(self, 'custom_parameters', {}).get('fasttree', None)
                if custom_params:
                    # Use custom parameters
                    cmd = ["fasttree"] + custom_params + [str(alignment_file)]
                    print(f"   Using custom FastTree parameters: {custom_params}")
                else:
                    # Use optimized defaults
                    cmd = ["fasttree"] + default_params + [str(alignment_file)]
                    print(f"   Using optimized FastTree defaults (GTR+Gamma model)")
                
                # Try both FastTree and fasttree command names
                for ft_cmd in ["fasttree", "FastTree"]:
                    try:
                        cmd[0] = ft_cmd
                        result = subprocess.run(cmd, capture_output=True, text=True)
                        
                        if result.returncode != 0:
                            print(f"   {ft_cmd} error (exit code {result.returncode}):")
                            if result.stderr:
                                print(f"   STDERR: {result.stderr}")
                            continue
                        
                        with open(tree_file, 'w') as output_file:
                            output_file.write(result.stdout)
                        
                        print(f"   Tree built successfully using {ft_cmd}")
                        return True
                    except (subprocess.CalledProcessError, FileNotFoundError):
                        continue
                
                print(f"   Error: Neither 'fasttree' nor 'FastTree' command found")
                return False
                
        except Exception as e:
            print(f"   Tree building error: {e}")
            return False
    
    def generate_summary_file(self, results, failed_families, batch_prefix=""):
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
                    f.write(f"  - {result['family']}: {result['num_species_bins']} species-BINs, "
                           f"{result['num_species']} species, {result['num_sequences']} sequences, "
                           f"{result['num_outgroups']} outgroups\n")
            
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
    parser.add_argument("--family-list", help="JSON file containing specific families to process (for batch processing)")
    parser.add_argument("--family-names", nargs='+', help="Specific family names to process")
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
    parser.add_argument("--batch-id", help="Batch identifier for logging purposes")
    parser.add_argument("--checkpoint-file", help="File to save/load processing checkpoint")
    parser.add_argument("--custom-parameters", help="JSON file with custom parameters for phylogenetic software")
    
    args = parser.parse_args()
    
    # Handle batch processing vs. discovery mode
    families_to_process = None
    
    if args.family_list:
        # Batch processing mode - load families from JSON file
        try:
            with open(args.family_list, 'r') as f:
                batch_data = json.load(f)
            families_to_process = batch_data.get('families', [])
            batch_info = f"Batch {batch_data.get('batch_id', '?')} ({len(families_to_process)} families)"
            print(f"Batch processing mode: {batch_info}")
        except Exception as e:
            print(f"Error loading family list from {args.family_list}: {e}")
            return 1
    elif args.family_names:
        # Specific families mode
        families_to_process = [{'family': name, 'order': None, 'kingdom': None} for name in args.family_names]
        print(f"Processing specific families: {args.family_names}")
    else:
        # Discovery mode - find all suitable families
        print("Discovery mode: scanning database for suitable families")
    
    # Load checkpoint if provided
    completed_families = set()
    checkpoint_data = {}
    if args.checkpoint_file and Path(args.checkpoint_file).exists():
        try:
            with open(args.checkpoint_file, 'r') as f:
                checkpoint_data = json.load(f)
            completed_families = set(checkpoint_data.get('completed_families', []))
            print(f"Loaded checkpoint: {len(completed_families)} families already completed")
        except Exception as e:
            print(f"Warning: Could not load checkpoint file: {e}")
    
    # Print system and tool information for debugging
    batch_prefix = f"[Batch {args.batch_id}] " if args.batch_id else ""
    print(f"{batch_prefix}System Information:")
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
    
    print(f"\n{batch_prefix}Tool Availability:")
    for program, version in program_versions.items():
        print(f"  {program}: {version}")
    
    if missing_programs:
        print(f"\n{batch_prefix}Missing required programs: {', '.join(missing_programs)}")
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
    
    # Load custom parameters if provided
    custom_parameters = {}
    if args.custom_parameters:
        try:
            with open(args.custom_parameters, 'r') as f:
                custom_parameters = json.load(f)
            print(f"Loaded custom parameters from: {args.custom_parameters}")
        except Exception as e:
            print(f"Warning: Could not load custom parameters from {args.custom_parameters}: {e}")
            print("Using default parameters instead")
    
    # Run HPC-optimized pipeline
    pipeline = PhylogeneticPipeline(
        args.database, 
        args.output_dir,
        generate_pdfs=args.generate_pdfs,
        bin_conflict_analysis=args.bin_conflict_analysis,
        custom_parameters=custom_parameters
    )
    
    # Set cleanup option
    pipeline.cleanup_intermediates = args.cleanup_intermediates
    
    # Add batch information to pipeline
    pipeline.batch_id = args.batch_id
    pipeline.checkpoint_file = args.checkpoint_file
    pipeline.completed_families = completed_families
    
    results = pipeline.run_pipeline_hpc_optimized(
        min_otus=args.min_otus, 
        max_families=args.max_families,
        alignment_method=args.alignment_method,
        tree_method=args.tree_method,
        bootstrap=args.bootstrap,
        num_outgroups=args.num_outgroups,
        families_to_process=families_to_process
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
