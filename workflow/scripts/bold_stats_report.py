#!/usr/bin/env python3
"""
BOLD Database Statistics Report Generator

This script generates comprehensive statistics from a BOLD database
and creates a PDF report with visualizations.

Usage:
    python bold_stats_report.py -d /path/to/bold.db -l /path/to/prescoring_filter.log -o report.pdf
"""

import argparse
import sqlite3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter
import re
from pathlib import Path

# Set matplotlib style
plt.style.use('default')
sns.set_palette("husl")

class BOLDStatsGenerator:
    def __init__(self, db_path, log_path=None):
        self.db_path = db_path
        self.log_path = log_path
        self.conn = sqlite3.connect(db_path)
        self.stats = {}
        
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn:
            self.conn.close()
    
    def get_table_info(self):
        """Get information about available tables and columns"""
        cursor = self.conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()
        
        print("Available tables:")
        for table in tables:
            table_name = table[0]
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = cursor.fetchall()
            print(f"\n{table_name}:")
            for col in columns:
                print(f"  - {col[1]} ({col[2]})")
    
    def extract_records_processed(self):
        """Extract number of records processed from pipeline summary or log file"""
        # First try to read from pipeline_summary.txt
        summary_path = Path(self.db_path).parent / "pipeline_summary.txt"
        if summary_path.exists():
            try:
                with open(summary_path, 'r') as f:
                    content = f.read()
                    match = re.search(r'Total records:\s*(\d+)', content, re.IGNORECASE)
                    if match:
                        return int(match.group(1))
            except Exception as e:
                print(f"Error reading pipeline summary: {e}")
        
        # Fallback to log file if provided
        if not self.log_path or not Path(self.log_path).exists():
            print(f"Warning: Neither pipeline summary nor log file found")
            return 0
            
        try:
            with open(self.log_path, 'r') as f:
                content = f.read()
                # Look for patterns like "processed X records" or similar
                patterns = [
                    r'processed\s+(\d+)\s+records',
                    r'(\d+)\s+records\s+processed',
                    r'total\s+records:\s*(\d+)',
                    r'input\s+records:\s*(\d+)'
                ]
                
                for pattern in patterns:
                    match = re.search(pattern, content, re.IGNORECASE)
                    if match:
                        return int(match.group(1))
                        
            print("Warning: Could not find records processed count in log file")
            return 0
        except Exception as e:
            print(f"Error reading log file: {e}")
            return 0    
    def get_basic_counts(self):
        """Get basic counts from the database using proper table structure"""
        queries = {
            'total_records': "SELECT COUNT(*) FROM bold",
            'species_count': "SELECT COUNT(DISTINCT species) FROM bold WHERE species IS NOT NULL",
            'bins_count': "SELECT COUNT(DISTINCT bin_uri) FROM bold WHERE bin_uri IS NOT NULL AND bin_uri != ''",
            'otus_count': "SELECT COUNT(DISTINCT otu_id) FROM bold_otus",
            'haplotypes_count': "SELECT COUNT(DISTINCT haplotype_id) FROM bold_haplotypes",
            'records_without_bin': "SELECT COUNT(*) FROM bold WHERE bin_uri IS NULL OR bin_uri = ''"
        }
        
        for key, query in queries.items():
            try:
                result = pd.read_sql_query(query, self.conn)
                self.stats[key] = result.iloc[0, 0]
            except Exception as e:
                print(f"Error executing query for {key}: {e}")
                self.stats[key] = 0
    
    def get_bags_distribution(self):
        """Get BAGS grade distribution using proper joins"""
        query = """
        SELECT b.bags_grade, COUNT(DISTINCT bold.species) as count 
        FROM bold 
        JOIN bags b ON bold.taxonid = b.taxonid 
        WHERE b.bags_grade IS NOT NULL AND bold.species IS NOT NULL
        GROUP BY b.bags_grade
        ORDER BY b.bags_grade
        """
        try:
            df = pd.read_sql_query(query, self.conn)
            self.stats['bags_distribution'] = df
        except Exception as e:
            print(f"Error getting BAGS distribution: {e}")
            self.stats['bags_distribution'] = pd.DataFrame()
    
    def get_rank_distribution(self):
        """Get rank distribution using proper joins"""
        query = """
        SELECT br.ranking, COUNT(*) as count 
        FROM bold 
        JOIN bold_ranks br ON bold.recordid = br.recordid 
        WHERE br.ranking IS NOT NULL 
        GROUP BY br.ranking 
        ORDER BY br.ranking
        """
        try:
            df = pd.read_sql_query(query, self.conn)
            self.stats['rank_distribution'] = df
        except Exception as e:
            print(f"Error getting rank distribution: {e}")
            self.stats['rank_distribution'] = pd.DataFrame()
    
    def get_sumscore_distribution(self):
        """Get sumscore distribution using proper joins"""
        query = """
        SELECT br.sumscore, COUNT(*) as count 
        FROM bold 
        JOIN bold_ranks br ON bold.recordid = br.recordid 
        WHERE br.sumscore IS NOT NULL 
        GROUP BY br.sumscore 
        ORDER BY br.sumscore
        """
        try:
            df = pd.read_sql_query(query, self.conn)
            self.stats['sumscore_distribution'] = df
        except Exception as e:
            print(f"Error getting sumscore distribution: {e}")
            self.stats['sumscore_distribution'] = pd.DataFrame()
    
    def get_criteria_distribution(self):
        """Get distribution of records against each criteria using normalized structure"""
        # Get all criteria from the criteria table
        criteria_query = "SELECT criterionid, name FROM criteria ORDER BY criterionid"
        
        try:
            criteria_df = pd.read_sql_query(criteria_query, self.conn)
            criteria_stats = {}
            
            for _, criterion in criteria_df.iterrows():
                criterion_id = criterion['criterionid']
                criterion_name = criterion['name']
                
                # Count records that pass (status=1) vs fail (status=0) this criterion
                query = f"""
                SELECT bc.status, COUNT(*) as count
                FROM bold_criteria bc
                WHERE bc.criterionid = {criterion_id}
                GROUP BY bc.status
                """
                
                try:
                    df = pd.read_sql_query(query, self.conn)
                    # Add human-readable labels
                    df['status_label'] = df['status'].map({0: 'Fail', 1: 'Pass'})
                    criteria_stats[criterion_name] = df
                except Exception as e:
                    print(f"Error getting criterion {criterion_name}: {e}")
                    continue
            
            self.stats['criteria_distribution'] = criteria_stats
            
        except Exception as e:
            print(f"Error getting criteria distribution: {e}")
            self.stats['criteria_distribution'] = {}
    
    def get_bags_rank_matrix(self):
        """Get matrix of BAGS grades vs ranks with species counts"""
        query = """
        WITH species_best_rank AS (
            SELECT 
                bold.species, 
                bags.bags_grade, 
                MIN(bold_ranks.ranking) as best_rank
            FROM bold 
            JOIN bags ON bold.taxonid = bags.taxonid
            JOIN bold_ranks ON bold.recordid = bold_ranks.recordid
            WHERE bold.species IS NOT NULL 
                AND bags.bags_grade IS NOT NULL 
                AND bold_ranks.ranking IS NOT NULL
            GROUP BY bold.species, bags.bags_grade
        )
        SELECT bags_grade, best_rank, COUNT(*) as species_count
        FROM species_best_rank
        GROUP BY bags_grade, best_rank
        ORDER BY bags_grade, best_rank
        """
        try:
            df = pd.read_sql_query(query, self.conn)
            # Pivot to create matrix
            matrix = df.pivot(index='bags_grade', columns='best_rank', values='species_count').fillna(0)
            self.stats['bags_rank_matrix'] = matrix
        except Exception as e:
            print(f"Error creating BAGS-rank matrix: {e}")
            self.stats['bags_rank_matrix'] = pd.DataFrame()
    
    def get_curation_categories(self):
        """Get species counts for different curation categories"""
        queries = {
            'auto_curable': """
                SELECT COUNT(DISTINCT bold.species) 
                FROM bold 
                JOIN bags ON bold.taxonid = bags.taxonid
                JOIN bold_ranks ON bold.recordid = bold_ranks.recordid
                WHERE bags.bags_grade IN ('A', 'B', 'D') 
                    AND bold_ranks.ranking IN (1, 2, 3)
                    AND bold.species IS NOT NULL
            """,
            'needs_attention': """
                SELECT COUNT(DISTINCT bold.species) 
                FROM bold 
                JOIN bags ON bold.taxonid = bags.taxonid
                JOIN bold_ranks ON bold.recordid = bold_ranks.recordid
                WHERE bags.bags_grade IN ('A', 'B', 'D') 
                    AND bold_ranks.ranking >= 4
                    AND bold.species IS NOT NULL
            """,
            'manual_intervention': """
                SELECT COUNT(DISTINCT bold.species) 
                FROM bold 
                JOIN bags ON bold.taxonid = bags.taxonid
                WHERE bags.bags_grade IN ('C', 'E')
                    AND bold.species IS NOT NULL
            """
        }
        
        for key, query in queries.items():
            try:
                result = pd.read_sql_query(query, self.conn)
                self.stats[key] = result.iloc[0, 0]
            except Exception as e:
                print(f"Error getting {key}: {e}")
                self.stats[key] = 0
    
    def get_taxonomic_breakdown(self):
        """Get breakdown across taxonomic hierarchy"""
        taxonomic_levels = ['phylum', 'class', 'order', 'family']
        breakdown = {}
        
        for level in taxonomic_levels:
            # Records count by taxonomic level
            records_query = f"""
                SELECT `{level}`, COUNT(*) as record_count 
                FROM bold 
                WHERE `{level}` IS NOT NULL 
                GROUP BY `{level}` 
                ORDER BY record_count DESC
            """
            
            # Species count by taxonomic level
            species_query = f"""
                SELECT `{level}`, COUNT(DISTINCT species) as species_count 
                FROM bold 
                WHERE `{level}` IS NOT NULL AND species IS NOT NULL
                GROUP BY `{level}` 
                ORDER BY species_count DESC
            """
            
            try:
                records_df = pd.read_sql_query(records_query, self.conn)
                species_df = pd.read_sql_query(species_query, self.conn)
                
                # Merge the dataframes
                combined = pd.merge(records_df, species_df, on=level, how='outer').fillna(0)
                breakdown[level] = combined
                
            except Exception as e:
                print(f"Error getting {level} breakdown: {e}")
                breakdown[level] = pd.DataFrame()
        
        self.stats['taxonomic_breakdown'] = breakdown
    
    def get_country_representatives_stats(self):
        """Get statistics about country representatives"""
        query = """
        SELECT 
            COUNT(*) as total_representatives,
            COUNT(DISTINCT country_iso) as countries_represented,
            COUNT(DISTINCT species) as species_represented,
            AVG(sumscore) as avg_sumscore,
            AVG(ranking) as avg_ranking
        FROM country_representatives
        """
        try:
            result = pd.read_sql_query(query, self.conn)
            self.stats['country_representatives'] = result.iloc[0].to_dict()
        except Exception as e:
            print(f"Error getting country representatives stats: {e}")
            self.stats['country_representatives'] = {}
    
    def generate_all_stats(self):
        """Generate all statistics"""
        print("Extracting records processed...")
        self.stats['records_processed'] = self.extract_records_processed()
        
        print("Getting basic counts...")
        self.get_basic_counts()
        
        print("Getting BAGS distribution...")
        self.get_bags_distribution()
        
        print("Getting rank distribution...")
        self.get_rank_distribution()
        
        print("Getting sumscore distribution...")
        self.get_sumscore_distribution()
        
        print("Getting criteria distribution...")
        self.get_criteria_distribution()
        
        print("Getting BAGS-rank matrix...")
        self.get_bags_rank_matrix()
        
        print("Getting curation categories...")
        self.get_curation_categories()
        
        print("Getting taxonomic breakdown...")
        self.get_taxonomic_breakdown()
        
        print("Getting country representatives stats...")
        self.get_country_representatives_stats()
        
        print("Statistics extraction complete!")
    
    def create_pdf_report(self, output_path):
        """Create comprehensive PDF report"""
        with PdfPages(output_path) as pdf:
            # Title page
            self._create_title_page(pdf)
            
            # Summary statistics
            self._create_summary_page(pdf)
            
            # BAGS distribution
            self._create_bags_page(pdf)
            
            # Rank distribution
            self._create_rank_page(pdf)
            
            # Sumscore distribution
            self._create_sumscore_page(pdf)
            
            # Criteria distribution
            self._create_criteria_pages(pdf)
            
            # BAGS-Rank matrix
            self._create_bags_rank_matrix_page(pdf)
            
            # Curation categories
            self._create_curation_page(pdf)
            
            # Taxonomic breakdown
            self._create_taxonomic_pages(pdf)
            
            # Country representatives
            self._create_country_representatives_page(pdf)
            
            # Additional statistics
            self._create_additional_stats_page(pdf)
    
    def _create_title_page(self, pdf):
        """Create title page"""
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.7, 'BOLD Database Statistics Report', 
                ha='center', va='center', fontsize=24, fontweight='bold')
        fig.text(0.5, 0.6, f'Database: {Path(self.db_path).name}', 
                ha='center', va='center', fontsize=14)
        fig.text(0.5, 0.5, f'Generated: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}', 
                ha='center', va='center', fontsize=12)
        plt.axis('off')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_summary_page(self, pdf):
        """Create summary statistics page"""
        fig, ax = plt.subplots(figsize=(8.5, 11))
        
        country_reps = self.stats.get('country_representatives', {})
        
        summary_text = f"""
SUMMARY STATISTICS

Records processed: {self.stats.get('records_processed', 'N/A'):,}
Total records: {self.stats.get('total_records', 0):,}
Unique species: {self.stats.get('species_count', 0):,}
Unique BINs: {self.stats.get('bins_count', 0):,}
Unique OTUs: {self.stats.get('otus_count', 0):,}
Unique Haplotypes: {self.stats.get('haplotypes_count', 0):,}
Records without BIN: {self.stats.get('records_without_bin', 0):,}

CURATION CATEGORIES

Auto-curable species: {self.stats.get('auto_curable', 0):,}
Needs attention: {self.stats.get('needs_attention', 0):,}
Manual intervention: {self.stats.get('manual_intervention', 0):,}

COUNTRY REPRESENTATIVES

Total representatives: {country_reps.get('total_representatives', 'N/A'):,}
Countries represented: {country_reps.get('countries_represented', 'N/A'):,}
Species represented: {country_reps.get('species_represented', 'N/A'):,}
        """
        
        ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=11,
                verticalalignment='top', fontfamily='monospace')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.title('Summary Statistics', fontsize=16, fontweight='bold', pad=20)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_bags_page(self, pdf):
        """Create BAGS distribution page"""
        if self.stats['bags_distribution'].empty:
            return
            
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 8.5))
        
        df = self.stats['bags_distribution']
        bars = ax1.bar(df['bags_grade'], df['count'])
        ax1.set_xlabel('BAGS Grade')
        ax1.set_ylabel('Number of Species')
        ax1.set_title('BAGS Grade Distribution')
        
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom')
        
        ax2.axis('tight')
        ax2.axis('off')
        table_data = [[grade, count] for grade, count in zip(df['bags_grade'], df['count'])]
        table = ax2.table(cellText=table_data, colLabels=['BAGS Grade', 'Species Count'],
                         cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1.2, 1.5)
        
        plt.suptitle('BAGS Grade Distribution', fontsize=16, fontweight='bold')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_rank_page(self, pdf):
        """Create rank distribution page"""
        if self.stats['rank_distribution'].empty:
            return
            
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 8.5))
        
        df = self.stats['rank_distribution']
        bars = ax1.bar(df['ranking'], df['count'])
        ax1.set_xlabel('Rank')
        ax1.set_ylabel('Number of Records')
        ax1.set_title('Rank Distribution')
        
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom')
        
        ax2.axis('tight')
        ax2.axis('off')
        table_data = [[rank, count] for rank, count in zip(df['ranking'], df['count'])]
        table = ax2.table(cellText=table_data, colLabels=['Rank', 'Record Count'],
                         cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1.2, 1.5)
        
        plt.suptitle('Rank Distribution', fontsize=16, fontweight='bold')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_sumscore_page(self, pdf):
        """Create sumscore distribution page"""
        if self.stats['sumscore_distribution'].empty:
            return
            
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 8.5))
        
        df = self.stats['sumscore_distribution']
        wedges, texts, autotexts = ax1.pie(df['count'], labels=df['sumscore'], autopct='%1.1f%%')
        ax1.set_title('Sumscore Distribution')
        
        ax2.axis('tight')
        ax2.axis('off')
        table_data = [[score, count, f"{count/df['count'].sum()*100:.1f}%"] 
                     for score, count in zip(df['sumscore'], df['count'])]
        table = ax2.table(cellText=table_data, colLabels=['Sumscore', 'Count', 'Percentage'],
                         cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        
        plt.suptitle('Sumscore Distribution (0-16)', fontsize=16, fontweight='bold')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_criteria_pages(self, pdf):
        """Create criteria distribution pages"""
        if not self.stats['criteria_distribution']:
            return
        
        criteria_items = list(self.stats['criteria_distribution'].items())
        
        for page_start in range(0, len(criteria_items), 4):
            page_criteria = criteria_items[page_start:page_start+4]
            
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
            axes = axes.flatten()
            
            for i, (criterion_name, df) in enumerate(page_criteria):
                if i >= 4:
                    break
                    
                ax = axes[i]
                
                if not df.empty:
                    colors = ['red', 'green']
                    wedges, texts, autotexts = ax.pie(df['count'], 
                                                     labels=df['status_label'], 
                                                     autopct='%1.1f%%',
                                                     colors=colors)
                ax.set_title(criterion_name, fontsize=10)
            
            for i in range(len(page_criteria), 4):
                axes[i].axis('off')
            
            plt.suptitle('Criteria Distribution (Pass/Fail)', fontsize=16, fontweight='bold')
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
    
    def _create_bags_rank_matrix_page(self, pdf):
        """Create BAGS-Rank matrix heatmap"""
        if self.stats['bags_rank_matrix'].empty:
            return
            
        fig, ax = plt.subplots(figsize=(10, 6))
        
        matrix = self.stats['bags_rank_matrix']
        proportions = matrix.div(matrix.sum().sum()) * 100
        
        sns.heatmap(proportions, annot=matrix.astype(int), fmt='d', 
                   cmap='Blues', ax=ax, cbar_kws={'label': 'Proportion (%)'})
        
        ax.set_xlabel('Best Rank')
        ax.set_ylabel('BAGS Grade')
        ax.set_title('Species Distribution: BAGS Grade vs Best Rank')
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_curation_page(self, pdf):
        """Create curation categories visualization"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 8.5))
        
        categories = ['Auto-curable', 'Needs attention', 'Manual intervention']
        counts = [self.stats['auto_curable'], self.stats['needs_attention'], 
                 self.stats['manual_intervention']]
        
        colors = ['green', 'orange', 'red']
        wedges, texts, autotexts = ax1.pie(counts, labels=categories, autopct='%1.1f%%', colors=colors)
        ax1.set_title('Curation Categories')
        
        bars = ax2.bar(categories, counts, color=colors)
        ax2.set_ylabel('Number of Species')
        ax2.set_title('Species by Curation Category')
        ax2.tick_params(axis='x', rotation=45)
        
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom')
        
        plt.suptitle('Curation Analysis', fontsize=16, fontweight='bold')
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_taxonomic_pages(self, pdf):
        """Create taxonomic breakdown pages"""
        for level, df in self.stats['taxonomic_breakdown'].items():
            if df.empty:
                continue
                
            df_top = df.head(15)
            
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 8.5))
            
            bars1 = ax1.barh(range(len(df_top)), df_top['record_count'])
            ax1.set_yticks(range(len(df_top)))
            ax1.set_yticklabels(df_top[level], fontsize=8)
            ax1.set_xlabel('Number of Records')
            ax1.set_title(f'Records by {level.title()}')
            ax1.invert_yaxis()
            
            bars2 = ax2.barh(range(len(df_top)), df_top['species_count'])
            ax2.set_yticks(range(len(df_top)))
            ax2.set_yticklabels(df_top[level], fontsize=8)
            ax2.set_xlabel('Number of Species')
            ax2.set_title(f'Species by {level.title()}')
            ax2.invert_yaxis()
            
            ax3.axis('tight')
            ax3.axis('off')
            table_data = [[taxon, records, species] for taxon, records, species in 
                         zip(df_top[level], df_top['record_count'], df_top['species_count'])]
            table = ax3.table(cellText=table_data, 
                             colLabels=[level.title(), 'Records', 'Species'],
                             cellLoc='left', loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
            
            total_records = df['record_count'].sum()
            total_species = df['species_count'].sum()
            unique_taxa = len(df)
            
            summary_text = f"""
{level.title()} Summary:
• Unique {level}s: {unique_taxa}
• Total records: {total_records:,}
• Total species: {total_species:,}
• Avg records per {level}: {total_records/unique_taxa:.1f}
• Avg species per {level}: {total_species/unique_taxa:.1f}
            """
            
            ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=10,
                    verticalalignment='top', fontfamily='monospace')
            ax4.set_xlim(0, 1)
            ax4.set_ylim(0, 1)
            ax4.axis('off')
            
            plt.suptitle(f'Taxonomic Breakdown: {level.title()}', fontsize=16, fontweight='bold')
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
    
    def _create_country_representatives_page(self, pdf):
        """Create country representatives analysis page"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 8.5))
        
        query = """
        SELECT country_iso, COUNT(*) as rep_count
        FROM country_representatives
        GROUP BY country_iso
        ORDER BY rep_count DESC
        LIMIT 15
        """
        try:
            country_counts = pd.read_sql_query(query, self.conn)
            
            bars1 = ax1.barh(range(len(country_counts)), country_counts['rep_count'])
            ax1.set_yticks(range(len(country_counts)))
            ax1.set_yticklabels(country_counts['country_iso'], fontsize=8)
            ax1.set_xlabel('Number of Representatives')
            ax1.set_title('Top Countries by Representatives')
            ax1.invert_yaxis()
        except:
            ax1.text(0.5, 0.5, 'No country data available', ha='center', va='center')
            ax1.set_title('Top Countries by Representatives')
        
        country_reps = self.stats.get('country_representatives', {})
        summary_text = f"""
Country Representatives Summary:

• Total representatives: {country_reps.get('total_representatives', 'N/A'):,}
• Countries represented: {country_reps.get('countries_represented', 'N/A'):,}
• Species with representatives: {country_reps.get('species_represented', 'N/A'):,}
• Average sumscore: {country_reps.get('avg_sumscore', 0):.1f}
• Average ranking: {country_reps.get('avg_ranking', 0):.1f}
        """
        
        ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace')
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        ax4.axis('off')
        
        ax2.axis('off')
        ax3.axis('off')
        
        plt.suptitle('Country Representatives Analysis', fontsize=16, fontweight='bold')
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_additional_stats_page(self, pdf):
        """Create additional statistics page"""
        fig, ax = plt.subplots(figsize=(8.5, 11))
        
        total_records = self.stats.get('total_records', 0)
        total_species = self.stats.get('species_count', 0)
        
        completion_rate = "N/A"
        if self.stats.get('records_processed', 0) > 0:
            completion_rate = f"{(total_records / self.stats['records_processed']) * 100:.1f}%"
        
        bin_coverage = "N/A"
        if total_records > 0:
            records_with_bin = total_records - self.stats.get('records_without_bin', 0)
            bin_coverage = f"{(records_with_bin / total_records) * 100:.1f}%"
        
        avg_records_per_species = total_records/total_species if total_species > 0 else 0
        avg_records_text = f"{avg_records_per_species:.1f}" if total_species > 0 else "N/A"
        
        additional_text = f"""
ADDITIONAL STATISTICS

Data Completeness:
• Processing completion rate: {completion_rate}
• BIN assignment coverage: {bin_coverage}
• Average records per species: {avg_records_text}

Quality Metrics:
• Most common BAGS grade: {self._get_most_common_bags_grade()}
• Most common rank: {self._get_most_common_rank()}
• Criteria pass rate: {self._get_overall_criteria_pass_rate()}

Database Information:
• Database file: {Path(self.db_path).name}
• File size: {self._get_file_size()}
• Available tables: {self._get_table_count()}
        """
        
        ax.text(0.1, 0.9, additional_text, transform=ax.transAxes, fontsize=11,
                verticalalignment='top', fontfamily='monospace')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.title('Additional Statistics', fontsize=16, fontweight='bold', pad=20)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _get_most_common_bags_grade(self):
        """Get most common BAGS grade"""
        if not self.stats['bags_distribution'].empty:
            df = self.stats['bags_distribution']
            most_common = df.loc[df['count'].idxmax(), 'bags_grade']
            return f"{most_common}"
        return "N/A"
    
    def _get_most_common_rank(self):
        """Get most common rank"""
        if not self.stats['rank_distribution'].empty:
            df = self.stats['rank_distribution']
            most_common = df.loc[df['count'].idxmax(), 'ranking']
            return f"{most_common}"
        return "N/A"
    
    def _get_file_size(self):
        """Get database file size"""
        try:
            size_bytes = Path(self.db_path).stat().st_size
            size_mb = size_bytes / (1024 * 1024)
            return f"{size_mb:.1f} MB"
        except:
            return "N/A"
    
    def _get_table_count(self):
        """Get number of tables in database"""
        try:
            cursor = self.conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table';")
            return cursor.fetchone()[0]
        except:
            return "N/A"
    
    def _get_overall_criteria_pass_rate(self):
        """Get overall criteria pass rate"""
        try:
            query = """
            SELECT 
                SUM(CASE WHEN status = 1 THEN 1 ELSE 0 END) * 100.0 / COUNT(*) as pass_rate
            FROM bold_criteria
            """
            result = pd.read_sql_query(query, self.conn)
            return f"{result.iloc[0, 0]:.1f}%"
        except:
            return "N/A"


def main():
    parser = argparse.ArgumentParser(description='Generate BOLD database statistics report')
    parser.add_argument('-d', '--database', required=True, help='Path to BOLD database file')
    parser.add_argument('-l', '--log', help='Path to prescoring_filter.log file (optional)')
    parser.add_argument('-o', '--output', default='bold_stats_report.pdf', help='Output PDF file path')
    parser.add_argument('--info', action='store_true', help='Show database table information and exit')
    
    args = parser.parse_args()
    
    if not Path(args.database).exists():
        print(f"Error: Database file not found: {args.database}")
        return 1
    
    try:
        with BOLDStatsGenerator(args.database, args.log) as generator:
            if args.info:
                generator.get_table_info()
                return 0
            
            print("Generating BOLD database statistics report...")
            generator.generate_all_stats()
            
            print(f"Creating PDF report: {args.output}")
            generator.create_pdf_report(args.output)
            
            print(f"Report generated successfully: {args.output}")
            return 0
            
    except Exception as e:
        print(f"Error generating report: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
