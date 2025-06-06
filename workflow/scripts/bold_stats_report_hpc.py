#!/usr/bin/env python3
"""
BOLD Database Statistics Report Generator (HPC Optimized)

This script generates comprehensive statistics from a BOLD database
and creates a PDF report with visualizations, optimized for large datasets
and HPC environments.

Usage:
    python bold_stats_report_hpc.py -d /path/to/bold.db -l /path/to/prescoring_filter.log -o report.pdf
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
import gc
import time
import psutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import warnings
from PIL import Image

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set matplotlib style
plt.style.use('default')
sns.set_palette("husl")

class MemoryMonitor:
    """Monitor memory usage throughout the process"""
    def __init__(self):
        self.peak_memory = 0
        self.process = psutil.Process()
    
    def get_memory_usage(self):
        """Get current memory usage in MB"""
        return self.process.memory_info().rss / (1024 * 1024)
    
    def log_memory(self, stage):
        """Log memory usage at a specific stage"""
        current_memory = self.get_memory_usage()
        self.peak_memory = max(self.peak_memory, current_memory)
        print(f"Memory usage at {stage}: {current_memory:.1f} MB (Peak: {self.peak_memory:.1f} MB)")


class HPCBOLDStatsGenerator:
    def __init__(self, db_path, log_path=None, skip_detailed_taxonomy=False, 
                 fast_mode=False, sample_size=100000, max_workers=4):
        self.db_path = db_path
        self.log_path = log_path
        self.skip_detailed_taxonomy = skip_detailed_taxonomy
        self.fast_mode = fast_mode
        self.sample_size = sample_size
        self.max_workers = max_workers
        self.conn = sqlite3.connect(db_path)
        self.stats = {}
        self.start_time = time.time()
        self.memory_monitor = MemoryMonitor()
        
        # Optimize SQLite for read-heavy workloads
        self.conn.execute("PRAGMA journal_mode = WAL")
        self.conn.execute("PRAGMA synchronous = NORMAL")
        self.conn.execute("PRAGMA cache_size = 10000")
        self.conn.execute("PRAGMA temp_store = MEMORY")
        
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn:
            self.conn.close()
    
    def create_performance_indexes(self):
        """Create temporary indexes for better performance"""
        print("Creating performance indexes...")
        try:
            # Only create if they don't exist
            indexes = [
                "CREATE INDEX IF NOT EXISTS idx_bold_species ON bold(species)",
                "CREATE INDEX IF NOT EXISTS idx_bold_bin_uri ON bold(bin_uri)",
                "CREATE INDEX IF NOT EXISTS idx_bold_recordid ON bold(recordid)",
                "CREATE INDEX IF NOT EXISTS idx_bags_taxonid ON bags(taxonid)",
                "CREATE INDEX IF NOT EXISTS idx_ranks_recordid ON bold_ranks(recordid)"
            ]
            
            for idx in indexes:
                self.conn.execute(idx)
            self.conn.commit()
            print("✓ Performance indexes created")
        except Exception as e:
            print(f"Warning: Could not create some indexes: {e}")
    
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
        """Get basic counts from the database with optimizations"""
        print("Getting basic counts...")
        
        queries = {
            'total_records': "SELECT COUNT(*) FROM bold",
            'species_count': "SELECT COUNT(DISTINCT species) FROM bold WHERE species IS NOT NULL",
            'bins_count': "SELECT COUNT(DISTINCT bin_uri) FROM bold WHERE bin_uri IS NOT NULL AND bin_uri != ''",
            'records_without_bin': "SELECT COUNT(*) FROM bold WHERE bin_uri IS NULL OR bin_uri = ''"
        }
        
        # Add OTU and haplotype counts if tables exist
        try:
            otu_count = pd.read_sql_query("SELECT COUNT(DISTINCT otu_id) FROM bold_otus", self.conn).iloc[0, 0]
            self.stats['otus_count'] = otu_count
        except:
            self.stats['otus_count'] = 0
            
        try:
            haplo_count = pd.read_sql_query("SELECT COUNT(DISTINCT haplotype_id) FROM bold_haplotypes", self.conn).iloc[0, 0]
            self.stats['haplotypes_count'] = haplo_count
        except:
            self.stats['haplotypes_count'] = 0
        
        for key, query in queries.items():
            try:
                if self.fast_mode and key in ['species_count', 'bins_count']:
                    # Use sampling for expensive queries in fast mode
                    sample_query = query.replace("FROM bold", f"FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {self.sample_size}) bold")
                    result = pd.read_sql_query(sample_query, self.conn)
                else:
                    result = pd.read_sql_query(query, self.conn)
                self.stats[key] = result.iloc[0, 0]
            except Exception as e:
                print(f"Error executing query for {key}: {e}")
                self.stats[key] = 0
        
        self.memory_monitor.log_memory("basic_counts")
    
    def get_bags_distribution(self):
        """Get BAGS grade distribution with optimizations"""
        print("Getting BAGS distribution...")
        query = """
        SELECT b.bags_grade, COUNT(DISTINCT bold.species) as count 
        FROM bold 
        JOIN bags b ON bold.taxonid = b.taxonid 
        WHERE b.bags_grade IS NOT NULL AND bold.species IS NOT NULL
        GROUP BY b.bags_grade
        ORDER BY b.bags_grade
        """
        try:
            if self.fast_mode:
                sample_query = """
                SELECT b.bags_grade, COUNT(DISTINCT bold.species) as count 
                FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {}) bold
                JOIN bags b ON bold.taxonid = b.taxonid 
                WHERE b.bags_grade IS NOT NULL AND bold.species IS NOT NULL
                GROUP BY b.bags_grade
                ORDER BY b.bags_grade
                """.format(self.sample_size)
                df = pd.read_sql_query(sample_query, self.conn)
            else:
                df = pd.read_sql_query(query, self.conn)
            self.stats['bags_distribution'] = df
        except Exception as e:
            print(f"Error getting BAGS distribution: {e}")
            self.stats['bags_distribution'] = pd.DataFrame()
        
        self.memory_monitor.log_memory("bags_distribution")
    
    def get_rank_distribution(self):
        """Get rank distribution with optimizations"""
        print("Getting rank distribution...")
        query = """
        SELECT br.ranking, COUNT(*) as count 
        FROM bold JOIN bold_ranks br ON bold.recordid = br.recordid 
        WHERE br.ranking IS NOT NULL 
        GROUP BY br.ranking ORDER BY br.ranking
        """
        try:
            if self.fast_mode:
                sample_query = query.replace("FROM bold", f"FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {self.sample_size}) bold")
                df = pd.read_sql_query(sample_query, self.conn)
            else:
                df = pd.read_sql_query(query, self.conn)
            self.stats['rank_distribution'] = df
        except Exception as e:
            print(f"Error getting rank distribution: {e}")
            self.stats['rank_distribution'] = pd.DataFrame()
        self.memory_monitor.log_memory("rank_distribution")
    
    def get_sumscore_distribution(self):
        """Get sumscore distribution using proper joins"""
        print("Getting sumscore distribution...")
        query = """
        SELECT br.sumscore, COUNT(*) as count 
        FROM bold JOIN bold_ranks br ON bold.recordid = br.recordid 
        WHERE br.sumscore IS NOT NULL 
        GROUP BY br.sumscore ORDER BY br.sumscore
        """
        try:
            if self.fast_mode:
                sample_query = query.replace("FROM bold", f"FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {self.sample_size}) bold")
                df = pd.read_sql_query(sample_query, self.conn)
            else:
                df = pd.read_sql_query(query, self.conn)
            self.stats['sumscore_distribution'] = df
        except Exception as e:
            print(f"Error getting sumscore distribution: {e}")
            self.stats['sumscore_distribution'] = pd.DataFrame()
        self.memory_monitor.log_memory("sumscore_distribution")
    
    def get_criteria_distribution(self):
        """Get distribution of records against each criteria"""
        print("Getting criteria distribution...")
        try:
            criteria_df = pd.read_sql_query("SELECT criterionid, name FROM criteria ORDER BY criterionid", self.conn)
            criteria_stats = {}
            
            for _, criterion in criteria_df.iterrows():
                criterion_id = criterion['criterionid']
                criterion_name = criterion['name']
                
                if self.fast_mode:
                    query = f"""
                    SELECT bc.status, COUNT(*) as count
                    FROM (SELECT * FROM bold_criteria ORDER BY RANDOM() LIMIT {self.sample_size}) bc
                    WHERE bc.criterionid = {criterion_id} GROUP BY bc.status
                    """
                else:
                    query = f"""
                    SELECT bc.status, COUNT(*) as count
                    FROM bold_criteria bc WHERE bc.criterionid = {criterion_id} GROUP BY bc.status
                    """
                
                try:
                    df = pd.read_sql_query(query, self.conn)
                    df['status_label'] = df['status'].map({0: 'Fail', 1: 'Pass'})
                    criteria_stats[criterion_name] = df
                except Exception as e:
                    print(f"Error getting criterion {criterion_name}: {e}")
                    continue
            
            self.stats['criteria_distribution'] = criteria_stats
        except Exception as e:
            print(f"Error getting criteria distribution: {e}")
            self.stats['criteria_distribution'] = {}
        self.memory_monitor.log_memory("criteria_distribution")
    
    def get_bags_rank_matrix(self):
        """Get matrix of BAGS grades vs ranks with species counts"""
        print("Getting BAGS-rank matrix...")
        if self.fast_mode:
            query = """
            SELECT bags.bags_grade, bold_ranks.ranking as best_rank, COUNT(DISTINCT bold.species) as species_count
            FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {}) bold
            JOIN bags ON bold.taxonid = bags.taxonid
            JOIN bold_ranks ON bold.recordid = bold_ranks.recordid
            WHERE bold.species IS NOT NULL AND bags.bags_grade IS NOT NULL AND bold_ranks.ranking IS NOT NULL
            GROUP BY bags.bags_grade, bold_ranks.ranking
            ORDER BY bags.bags_grade, bold_ranks.ranking
            """.format(self.sample_size)
        else:
            query = """
            WITH species_best_rank AS (
                SELECT bold.species, bags.bags_grade, MIN(bold_ranks.ranking) as best_rank
                FROM bold 
                JOIN bags ON bold.taxonid = bags.taxonid
                JOIN bold_ranks ON bold.recordid = bold_ranks.recordid
                WHERE bold.species IS NOT NULL AND bags.bags_grade IS NOT NULL AND bold_ranks.ranking IS NOT NULL
                GROUP BY bold.species, bags.bags_grade
            )
            SELECT bags_grade, best_rank, COUNT(*) as species_count
            FROM species_best_rank GROUP BY bags_grade, best_rank
            ORDER BY bags_grade, best_rank
            """
        try:
            df = pd.read_sql_query(query, self.conn)
            matrix = df.pivot(index='bags_grade', columns='best_rank', values='species_count').fillna(0)
            self.stats['bags_rank_matrix'] = matrix
        except Exception as e:
            print(f"Error creating BAGS-rank matrix: {e}")
            self.stats['bags_rank_matrix'] = pd.DataFrame()
        self.memory_monitor.log_memory("bags_rank_matrix")
    
    def get_curation_categories(self):
        """Get species counts for different curation categories"""
        print("Getting curation categories...")
        queries = {
            'auto_curatable': """
                SELECT COUNT(DISTINCT bold.species) 
                FROM bold JOIN bags ON bold.taxonid = bags.taxonid
                JOIN bold_ranks ON bold.recordid = bold_ranks.recordid
                WHERE bags.bags_grade IN ('A', 'B', 'D') AND bold_ranks.ranking IN (1, 2, 3)
                    AND bold.species IS NOT NULL
            """,
            'needs_attention': """
                SELECT COUNT(DISTINCT bold.species) 
                FROM bold JOIN bags ON bold.taxonid = bags.taxonid
                JOIN bold_ranks ON bold.recordid = bold_ranks.recordid
                WHERE bags.bags_grade IN ('A', 'B', 'D') AND bold_ranks.ranking >= 4
                    AND bold.species IS NOT NULL
            """,
            'manual_intervention': """
                SELECT COUNT(DISTINCT bold.species) 
                FROM bold JOIN bags ON bold.taxonid = bags.taxonid
                WHERE bags.bags_grade IN ('C', 'E') AND bold.species IS NOT NULL
            """
        }
        
        for key, query in queries.items():
            try:
                if self.fast_mode:
                    sample_query = query.replace("FROM bold", f"FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {self.sample_size}) bold")
                    result = pd.read_sql_query(sample_query, self.conn)
                else:
                    result = pd.read_sql_query(query, self.conn)
                self.stats[key] = result.iloc[0, 0]
            except Exception as e:
                print(f"Error getting {key}: {e}")
                self.stats[key] = 0
        self.memory_monitor.log_memory("curation_categories")
    
    def get_taxonomic_breakdown(self):
        """Get breakdown across taxonomic hierarchy"""
        print("Getting taxonomic breakdown...")
        taxonomic_levels = ['phylum', 'class', 'order', 'family']
        breakdown = {}
        
        for level in taxonomic_levels:
            print(f"  Processing {level}...")
            records_query = f"""
                SELECT `{level}`, COUNT(*) as record_count 
                FROM bold WHERE `{level}` IS NOT NULL 
                GROUP BY `{level}` ORDER BY record_count DESC
            """
            species_query = f"""
                SELECT `{level}`, COUNT(DISTINCT species) as species_count 
                FROM bold WHERE `{level}` IS NOT NULL AND species IS NOT NULL
                GROUP BY `{level}` ORDER BY species_count DESC
            """
            
            try:
                if self.fast_mode:
                    records_sample = records_query.replace("FROM bold", f"FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {self.sample_size}) bold")
                    species_sample = species_query.replace("FROM bold", f"FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {self.sample_size}) bold")
                    records_df = pd.read_sql_query(records_sample, self.conn)
                    species_df = pd.read_sql_query(species_sample, self.conn)
                else:
                    records_df = pd.read_sql_query(records_query, self.conn)
                    species_df = pd.read_sql_query(species_query, self.conn)
                
                combined = pd.merge(records_df, species_df, on=level, how='outer').fillna(0)
                breakdown[level] = combined
            except Exception as e:
                print(f"Error getting {level} breakdown: {e}")
                breakdown[level] = pd.DataFrame()
        
        self.stats['taxonomic_breakdown'] = breakdown
        self.memory_monitor.log_memory("taxonomic_breakdown")
    
    def get_country_representatives_stats(self):
        """Get statistics about country representatives"""
        print("Getting country representatives stats...")
        query = """
        SELECT COUNT(*) as total_representatives, COUNT(DISTINCT cr.country_iso) as countries_represented,
               COUNT(DISTINCT cr.species) as species_represented, AVG(cr.sumscore) as avg_sumscore,
               AVG(cr.ranking) as avg_ranking
        FROM country_representatives cr
        """
        try:
            if self.fast_mode:
                sample_query = f"""
                SELECT COUNT(*) as total_representatives, COUNT(DISTINCT cr.country_iso) as countries_represented,
                       COUNT(DISTINCT cr.species) as species_represented, AVG(cr.sumscore) as avg_sumscore,
                       AVG(cr.ranking) as avg_ranking
                FROM (SELECT * FROM country_representatives ORDER BY RANDOM() LIMIT {self.sample_size}) cr
                """
                result = pd.read_sql_query(sample_query, self.conn)
            else:
                result = pd.read_sql_query(query, self.conn)
            self.stats['country_representatives'] = result.iloc[0].to_dict()
        except Exception as e:
            print(f"Error getting country representatives stats: {e}")
            self.stats['country_representatives'] = {}
        self.memory_monitor.log_memory("country_representatives")
    
    def get_geographic_distribution(self):
        """Get geographic distribution of records by country"""
        print("Getting geographic distribution...")
        country_columns = ['country/ocean', 'country', 'country_ocean']
        
        for col in country_columns:
            try:
                test_query = f"SELECT `{col}` FROM bold LIMIT 1"
                pd.read_sql_query(test_query, self.conn)
                
                query = f"""
                SELECT COALESCE(`{col}`, 'Unknown') as country_name,
                       COUNT(*) as record_count, COUNT(DISTINCT species) as species_count
                FROM bold WHERE `{col}` IS NOT NULL AND `{col}` != ''
                GROUP BY `{col}` ORDER BY record_count DESC
                """
                
                if self.fast_mode:
                    sample_query = query.replace("FROM bold", f"FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {self.sample_size}) bold")
                    df = pd.read_sql_query(sample_query, self.conn)
                else:
                    df = pd.read_sql_query(query, self.conn)
                
                self.stats['geographic_distribution'] = df
                self.stats['country_column_used'] = col
                self.memory_monitor.log_memory("geographic_distribution")
                return
            except Exception:
                continue
        
        # Fallback
        try:
            query = """
            SELECT COALESCE(country, 'Unknown') as country_name,
                   COUNT(*) as record_count, COUNT(DISTINCT species) as species_count
            FROM bold WHERE country IS NOT NULL AND country != ''
            GROUP BY country ORDER BY record_count DESC
            """
            if self.fast_mode:
                sample_query = query.replace("FROM bold", f"FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {self.sample_size}) bold")
                df = pd.read_sql_query(sample_query, self.conn)
            else:
                df = pd.read_sql_query(query, self.conn)
            self.stats['geographic_distribution'] = df
            self.stats['country_column_used'] = 'country'
        except Exception as e:
            print(f"Error getting geographic distribution: {e}")
            self.stats['geographic_distribution'] = pd.DataFrame()
            self.stats['country_column_used'] = None
        
        self.memory_monitor.log_memory("geographic_distribution")
    
    def generate_all_stats_parallel(self):
        """Generate all statistics with parallel processing where possible"""
        print("=== Starting HPC-Optimized Statistics Generation ===")
        print(f"Fast mode: {self.fast_mode}")
        print(f"Sample size: {self.sample_size:,}" if self.fast_mode else "Full dataset")
        print(f"Max workers: {self.max_workers}")
        print(f"Skip detailed taxonomy: {self.skip_detailed_taxonomy}")
        
        # Create performance indexes first
        if not self.fast_mode:
            self.create_performance_indexes()
        
        # Sequential tasks that must run first
        print("\n=== Phase 1: Basic Statistics ===")
        self.memory_monitor.log_memory("start")
        
        self.stats['records_processed'] = self.extract_records_processed()
        self.get_basic_counts()
        
        # Distribution analysis tasks
        print("\n=== Phase 2: Distribution Analysis ===")
        
        analysis_tasks = [
            ("BAGS distribution", self.get_bags_distribution),
            ("Rank distribution", self.get_rank_distribution),
            ("Sumscore distribution", self.get_sumscore_distribution),
            ("Criteria distribution", self.get_criteria_distribution),
            ("BAGS-rank matrix", self.get_bags_rank_matrix),
            ("Curation categories", self.get_curation_categories),
            ("Taxonomic breakdown", self.get_taxonomic_breakdown),
            ("Country representatives stats", self.get_country_representatives_stats),
            ("Geographic distribution", self.get_geographic_distribution)
        ]
        
        # Run tasks sequentially for now (can be parallelized later)
        for task_name, task_func in analysis_tasks:
            try:
                print(f"Processing: {task_name}")
                task_func()
                print(f"✓ Completed: {task_name}")
                gc.collect()  # Force garbage collection after each task
            except Exception as e:
                print(f"✗ Error in {task_name}: {e}")
        
        # Final memory report
        elapsed_time = time.time() - self.start_time
        print(f"\n=== Statistics Generation Complete ===")
        print(f"Total time: {elapsed_time/60:.1f} minutes")
        self.memory_monitor.log_memory("final")
        print(f"Peak memory usage: {self.memory_monitor.peak_memory:.1f} MB")
    
    def create_pdf_report_optimized(self, output_path):
        """Create comprehensive PDF report"""
        print(f"\n=== Creating PDF Report: {output_path} ===")
        
        try:
            with PdfPages(output_path) as pdf:
                # Title page
                self._create_title_page(pdf)
                
                # Summary statistics
                self._create_summary_page(pdf)
                
                # BAGS distribution
                if not self.stats.get('bags_distribution', pd.DataFrame()).empty:
                    self._create_bags_page(pdf)
                
                # Rank distribution
                if not self.stats.get('rank_distribution', pd.DataFrame()).empty:
                    self._create_rank_page(pdf)
                
                # Sumscore distribution
                if not self.stats.get('sumscore_distribution', pd.DataFrame()).empty:
                    self._create_sumscore_page(pdf)
                
                # Criteria distribution
                if self.stats.get('criteria_distribution'):
                    self._create_criteria_pages(pdf)
                
                # BAGS-Rank matrix
                if not self.stats.get('bags_rank_matrix', pd.DataFrame()).empty:
                    self._create_bags_rank_matrix_page(pdf)
                
                # Curation categories
                self._create_curation_page(pdf)
                
                # Country representatives
                if self.stats.get('country_representatives'):
                    self._create_country_representatives_page(pdf)
                
                # Geographic distribution
                if not self.stats.get('geographic_distribution', pd.DataFrame()).empty:
                    self._create_geographic_distribution_page(pdf)
                
                # Taxonomic breakdown (conditional)
                if not self.skip_detailed_taxonomy:
                    print("Generating detailed taxonomic breakdown pages...")
                    self._create_taxonomic_pages(pdf)
                else:
                    print("Skipping detailed taxonomic breakdown (--skip-detailed-taxonomy flag used)")
                    self._create_basic_taxonomic_summary_page(pdf)
                
                # Additional statistics
                self._create_additional_stats_page(pdf)
                
                # Performance metrics
                self._create_performance_page(pdf)
                
            print(f"✓ PDF report created successfully: {output_path}")
            file_size = Path(output_path).stat().st_size / (1024 * 1024)  # MB
            print(f"Report size: {file_size:.1f} MB")
            
        except Exception as e:
            print(f"Error creating PDF report: {e}")
            raise
    
    def _create_title_page(self, pdf):
        """Create title page with system information"""
        fig = plt.figure(figsize=(8.5, 11))
        
        mode_text = "Fast Mode (Sampled)" if self.fast_mode else "Full Analysis"
        sample_text = f"Sample size: {self.sample_size:,}" if self.fast_mode else ""
        
        fig.text(0.5, 0.7, 'BOLD Database Statistics Report', 
                ha='center', va='center', fontsize=24, fontweight='bold')
        fig.text(0.5, 0.6, f'Database: {Path(self.db_path).name}', 
                ha='center', va='center', fontsize=14)
        fig.text(0.5, 0.55, f'Mode: {mode_text}', 
                ha='center', va='center', fontsize=12)
        if sample_text:
            fig.text(0.5, 0.5, sample_text, 
                    ha='center', va='center', fontsize=12)
        fig.text(0.5, 0.45, f'Generated: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}', 
                ha='center', va='center', fontsize=12)
        fig.text(0.5, 0.4, f'Processing time: {(time.time() - self.start_time)/60:.1f} minutes', 
                ha='center', va='center', fontsize=10)
        fig.text(0.5, 0.35, f'Peak memory: {self.memory_monitor.peak_memory:.1f} MB', 
                ha='center', va='center', fontsize=10)
        
        plt.axis('off')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_summary_page(self, pdf):
        """Create enhanced summary statistics page"""
        fig, ax = plt.subplots(figsize=(8.5, 11))
        
        mode_note = "(Estimated from sample)" if self.fast_mode else ""
        country_reps = self.stats.get('country_representatives', {})
        
        summary_text = f"""
SUMMARY STATISTICS {mode_note}

Records processed: {self.stats.get('records_processed', 'N/A')}
Total records: {self.stats.get('total_records', 0):,}
Unique species: {self.stats.get('species_count', 0):,}
Unique BINs: {self.stats.get('bins_count', 0):,}
Unique OTUs: {self.stats.get('otus_count', 0):,}
Unique Haplotypes: {self.stats.get('haplotypes_count', 0):,}
Records without BIN: {self.stats.get('records_without_bin', 0):,}

CURATION CATEGORIES

Auto-curatable species: {self.stats.get('auto_curatable', 0):,}
Needs attention: {self.stats.get('needs_attention', 0):,}
Manual intervention: {self.stats.get('manual_intervention', 0):,}

COUNTRY REPRESENTATIVES

Total representatives: {country_reps.get('total_representatives', 'N/A')}
Countries represented: {country_reps.get('countries_represented', 'N/A')}
Species represented: {country_reps.get('species_represented', 'N/A')}

PROCESSING PERFORMANCE

Total processing time: {(time.time() - self.start_time)/60:.1f} minutes
Peak memory usage: {self.memory_monitor.peak_memory:.1f} MB
Processing mode: {'Fast (sampled)' if self.fast_mode else 'Full dataset'}
Database size: {self._get_file_size()}
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
            ax1.text(bar.get_x() + bar.get_width()/2., height, f'{int(height)}', ha='center', va='bottom')
        
        # Replace table with BAGS scoring summary text
        ax2.axis('off')
        bags_summary_text = """BAGS scoring summary:

Grade A = Species has > 10 specimens in 1 BIN

Grade B = Species has 3 - 10 specimens in 1 BIN

Grade C = Species has > 1 BIN (BIN splitting)

Grade D = Species has < 3 specimens in 1 BIN

Grade E = BIN with > 1 species (BIN sharing, includes Genus sp. or similar)

Grade F = other"""
        
        ax2.text(0.1, 0.9, bags_summary_text, transform=ax2.transAxes, fontsize=11,
                verticalalignment='top', fontfamily='monospace')
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
        
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
            ax1.text(bar.get_x() + bar.get_width()/2., height, f'{int(height)}', ha='center', va='bottom')
        
        # Replace table with image
        ax2.axis('off')
        try:
            # Load and display the rank image
            from PIL import Image
            img_path = r"C:\GitHub\bold-library-curation\doc\bold-curation-ranks.png"
            img = Image.open(img_path)
            ax2.imshow(img)
            ax2.set_title('Rank Categories')
        except Exception as e:
            # Fallback if image can't be loaded
            ax2.text(0.5, 0.5, f'Image not available:\n{img_path}', 
                    ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('Rank Categories')
        
        plt.suptitle('Rank Distribution', fontsize=16, fontweight='bold')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_sumscore_page(self, pdf):
        """Create sumscore distribution page"""
        if self.stats['sumscore_distribution'].empty:
            return
        fig, ax = plt.subplots(figsize=(11, 8.5))
        df = self.stats['sumscore_distribution']
        bars = ax.bar(df['sumscore'], df['count'])
        ax.set_xlabel('Sumscore')
        ax.set_ylabel('Number of Records')
        ax.set_title('Sumscore Distribution (0-16)')
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{int(height)}', ha='center', va='bottom')
        
        # Add figure caption at the bottom
        caption_text = "The sum of all passed criteria, out of a total of 16 criteria assessed for each record. A higher score indicates a record with better specimen data."
        fig.text(0.5, 0.02, caption_text, ha='center', va='bottom', fontsize=10, 
                style='italic', wrap=True)
        
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)  # Make room for caption
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    # Essential placeholder methods for PDF generation
    def _create_criteria_pages(self, pdf):
        """Create criteria distribution pages"""
        if not self.stats['criteria_distribution']:
            return
        
        # Define the specific order for criteria across two pages
        page_1_criteria = [
            'SPECIES_ID', 'TYPE_SPECIMEN', 'SEQ_QUALITY', 'HAS_IMAGE',
            'ID_METHOD', 'PUBLIC_VOUCHER', 'INSTITUTION', 'MUSEUM_ID'
        ]
        
        page_2_criteria = [
            'COUNTRY', 'REGION', 'SECTOR', 'SITE',
            'COORD', 'COLLECTORS', 'COLLECTION_DATE', 'IDENTIFIER'
        ]
        
        # Create pages with 8 criteria each
        for page_num, criteria_list in enumerate([page_1_criteria, page_2_criteria], 1):
            fig, axes = plt.subplots(2, 4, figsize=(16, 8.5))
            axes = axes.flatten()
            
            for i, criterion_name in enumerate(criteria_list):
                ax = axes[i]
                
                # Check if we have data for this criterion
                if criterion_name in self.stats['criteria_distribution']:
                    df = self.stats['criteria_distribution'][criterion_name]
                    if not df.empty:
                        colors = ['red', 'green']
                        wedges, texts, autotexts = ax.pie(df['count'], 
                                                         labels=df['status_label'], 
                                                         autopct='%1.1f%%',
                                                         colors=colors)
                    else:
                        ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                else:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                
                ax.set_title(criterion_name, fontsize=10)
            
            plt.suptitle(f'Criteria Distribution (Pass/Fail) - Page {page_num}', fontsize=16, fontweight='bold')
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
        
        # Add figure caption at the bottom
        caption_text = "The rank of the best specimen record within each species within each BAGS grade. Top left box indicates the number of species with a BAGS grade A which have a specimen of rank of 1 (= type specimens)."
        fig.text(0.5, 0.02, caption_text, ha='center', va='bottom', fontsize=10, 
                style='italic', wrap=True)
        
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)  # Make room for caption
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    def _create_curation_page(self, pdf):
        """Create curation categories visualization"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 8.5))
        
        categories = ['Auto-curatable', 'Needs attention', 'Manual intervention']
        counts = [self.stats['auto_curatable'], self.stats['needs_attention'], 
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
        
        # Add figure caption at the bottom
        caption_text = """Distribution of species which are possible to auto-curate, probably need attention and require manual intervention.
Auto-curatable = BAGS grade A, B or D and specimen ranks 1, 2, or 3
Need attention = BAGS grade A, B or D and specimen ranks >=4
Manual intervention = BAGS grade C or E"""
        fig.text(0.5, 0.02, caption_text, ha='center', va='bottom', fontsize=9, 
                style='italic', wrap=True)
        
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.2)  # Make room for caption
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    def _create_country_representatives_page(self, pdf):
        """Create country representatives analysis page"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 8.5))
        
        query = """
        SELECT COALESCE(b.`country/ocean`, cr.country_iso) as country_name, COUNT(*) as rep_count
        FROM country_representatives cr
        LEFT JOIN bold b ON cr.recordid = b.recordid
        GROUP BY COALESCE(b.`country/ocean`, cr.country_iso)
        ORDER BY rep_count DESC
        LIMIT 15
        """
        try:
            if self.fast_mode:
                sample_query = f"""
                SELECT COALESCE(b.`country/ocean`, cr.country_iso) as country_name, COUNT(*) as rep_count
                FROM (SELECT * FROM country_representatives ORDER BY RANDOM() LIMIT {self.sample_size}) cr
                LEFT JOIN bold b ON cr.recordid = b.recordid
                GROUP BY COALESCE(b.`country/ocean`, cr.country_iso)
                ORDER BY rep_count DESC
                LIMIT 15
                """
                country_counts = pd.read_sql_query(sample_query, self.conn)
            else:
                country_counts = pd.read_sql_query(query, self.conn)
            
            bars1 = ax1.barh(range(len(country_counts)), country_counts['rep_count'])
            ax1.set_yticks(range(len(country_counts)))
            ax1.set_yticklabels(country_counts['country_name'], fontsize=8)
            ax1.set_xlabel('Number of Representatives')
            ax1.set_title('Top Countries by Representatives')
            ax1.invert_yaxis()
        except:
            ax1.text(0.5, 0.5, 'No country data available', ha='center', va='center')
            ax1.set_title('Top Countries by Representatives')
        
        country_reps = self.stats.get('country_representatives', {})
        summary_text = f"""
Country Representatives Summary:

• Total representatives: {country_reps.get('total_representatives', 'N/A')}
• Countries represented: {country_reps.get('countries_represented', 'N/A')}
• Species with representatives: {country_reps.get('species_represented', 'N/A')}
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
    def _create_geographic_distribution_page(self, pdf):
        """Create geographic distribution analysis page"""
        fig = plt.figure(figsize=(16, 11))
        
        # Create a grid layout for multiple plots
        gs = fig.add_gridspec(3, 3, height_ratios=[2, 2, 1], width_ratios=[2, 2, 1])
        
        geo_df = self.stats.get('geographic_distribution', pd.DataFrame())
        
        if geo_df.empty:
            # Create a simple message if no geographic data
            ax = fig.add_subplot(gs[:, :])
            ax.text(0.5, 0.5, 'No geographic data available', 
                   ha='center', va='center', fontsize=16)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
        else:
            # Top 20 countries by record count
            ax1 = fig.add_subplot(gs[0, :2])
            top_countries = geo_df.head(20)
            
            bars = ax1.barh(range(len(top_countries)), top_countries['record_count'])
            ax1.set_yticks(range(len(top_countries)))
            ax1.set_yticklabels(top_countries['country_name'], fontsize=8)
            ax1.set_xlabel('Number of Records')
            ax1.set_title('Top 20 Countries by Record Count')
            ax1.invert_yaxis()
            
            # Add value labels on bars for top 10
            for i, bar in enumerate(bars[:10]):
                width = bar.get_width()
                ax1.text(width, bar.get_y() + bar.get_height()/2., 
                        f'{int(width):,}', ha='left', va='center', fontsize=8)
            
            # Top 15 countries by species diversity
            ax2 = fig.add_subplot(gs[1, :2])
            top_species = geo_df.nlargest(15, 'species_count')
            
            bars2 = ax2.barh(range(len(top_species)), top_species['species_count'])
            ax2.set_yticks(range(len(top_species)))
            ax2.set_yticklabels(top_species['country_name'], fontsize=8)
            ax2.set_xlabel('Number of Species')
            ax2.set_title('Top 15 Countries by Species Diversity')
            ax2.invert_yaxis()
            
            # Add value labels on bars for top 10
            for i, bar in enumerate(bars2[:10]):
                width = bar.get_width()
                ax2.text(width, bar.get_y() + bar.get_height()/2., 
                        f'{int(width):,}', ha='left', va='center', fontsize=8)
            
            # Summary statistics
            ax3 = fig.add_subplot(gs[:2, 2])
            ax3.axis('off')
            
            total_countries = len(geo_df)
            total_records = geo_df['record_count'].sum()
            total_species = geo_df['species_count'].sum()
            avg_records_per_country = total_records / total_countries if total_countries > 0 else 0
            avg_species_per_country = total_species / total_countries if total_countries > 0 else 0
            
            # Get top countries info safely
            if not top_countries.empty:
                top_record_country = top_countries.iloc[0]['country_name']
                top_record_count = f"{top_countries.iloc[0]['record_count']:,}"
            else:
                top_record_country = 'N/A'
                top_record_count = 'N/A'
                
            if not top_species.empty:
                top_species_country = top_species.iloc[0]['country_name']
                top_species_count = f"{top_species.iloc[0]['species_count']:,}"
            else:
                top_species_country = 'N/A'
                top_species_count = 'N/A'
            
            summary_text = f"""
Geographic Summary:

• Total countries: {total_countries:,}
• Total records: {total_records:,}
• Total species: {total_species:,}

• Avg records/country: {avg_records_per_country:.0f}
• Avg species/country: {avg_species_per_country:.0f}

Top Contributors:
• Most records: 
  {top_record_country}
  ({top_record_count} records)

• Most species: 
  {top_species_country}
  ({top_species_count} species)

Data source: {self.stats.get('country_column_used', 'Unknown')} column
            """
            
            ax3.text(0.05, 0.95, summary_text, transform=ax3.transAxes, fontsize=10,
                    verticalalignment='top', fontfamily='monospace')
            
            # Distribution chart (records vs species)
            ax4 = fig.add_subplot(gs[2, :2])
            
            # Scatter plot of records vs species for top 30 countries
            scatter_data = geo_df.head(30)
            scatter = ax4.scatter(scatter_data['record_count'], scatter_data['species_count'], 
                                alpha=0.6, s=50)
            
            ax4.set_xlabel('Number of Records')
            ax4.set_ylabel('Number of Species')
            ax4.set_title('Records vs Species Diversity (Top 30 Countries)')
            
            # Add trend line
            if len(scatter_data) > 1:
                z = np.polyfit(scatter_data['record_count'], scatter_data['species_count'], 1)
                p = np.poly1d(z)
                ax4.plot(scatter_data['record_count'], p(scatter_data['record_count']), 
                        "r--", alpha=0.8, linewidth=1)
            
            # Annotate a few top countries
            for i, row in scatter_data.head(5).iterrows():
                ax4.annotate(row['country_name'], 
                           (row['record_count'], row['species_count']),
                           xytext=(5, 5), textcoords='offset points', 
                           fontsize=8, alpha=0.7)
        
        plt.suptitle('Geographic Distribution of Records', fontsize=16, fontweight='bold')
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    def _create_taxonomic_pages(self, pdf):
        """Create taxonomic breakdown pages"""
        # First create pages for phylum, class, and order
        for level in ['phylum', 'class', 'order']:
            df = self.stats['taxonomic_breakdown'].get(level, pd.DataFrame())
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
        
        # Only create family breakdown pages if detailed taxonomy is enabled
        if not self.skip_detailed_taxonomy:
            self._create_family_by_order_pages(pdf)
    
    def _create_family_by_order_pages(self, pdf):
        """Create family breakdown pages for each order"""
        # Get all orders from the database
        orders_query = """
        SELECT DISTINCT `order`
        FROM bold 
        WHERE `order` IS NOT NULL AND family IS NOT NULL
        ORDER BY `order`
        """
        
        try:
            if self.fast_mode:
                # Limit the number of orders processed in fast mode
                orders_query += " LIMIT 20"
            
            orders_df = pd.read_sql_query(orders_query, self.conn)
            
            for _, order_row in orders_df.iterrows():
                order_name = order_row['order']
                
                # Get family breakdown for this order
                family_query = f"""
                SELECT family, COUNT(*) as record_count 
                FROM bold 
                WHERE `order` = '{order_name}' AND family IS NOT NULL 
                GROUP BY family 
                ORDER BY record_count DESC
                """
                
                species_query = f"""
                SELECT family, COUNT(DISTINCT species) as species_count 
                FROM bold 
                WHERE `order` = '{order_name}' AND family IS NOT NULL AND species IS NOT NULL
                GROUP BY family 
                ORDER BY species_count DESC
                """
                
                try:
                    if self.fast_mode:
                        # Limit results in fast mode
                        family_query += " LIMIT 20"
                        species_query += " LIMIT 20"
                    
                    records_df = pd.read_sql_query(family_query, self.conn)
                    species_df = pd.read_sql_query(species_query, self.conn)
                    
                    # Merge the dataframes
                    combined = pd.merge(records_df, species_df, on='family', how='outer').fillna(0)
                    
                    if combined.empty:
                        continue
                    
                    # Limit to top 15 families for readability
                    df_top = combined.head(15)
                    
                    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 8.5))
                    
                    bars1 = ax1.barh(range(len(df_top)), df_top['record_count'])
                    ax1.set_yticks(range(len(df_top)))
                    ax1.set_yticklabels(df_top['family'], fontsize=8)
                    ax1.set_xlabel('Number of Records')
                    ax1.set_title(f'Records by Family')
                    ax1.invert_yaxis()
                    
                    bars2 = ax2.barh(range(len(df_top)), df_top['species_count'])
                    ax2.set_yticks(range(len(df_top)))
                    ax2.set_yticklabels(df_top['family'], fontsize=8)
                    ax2.set_xlabel('Number of Species')
                    ax2.set_title(f'Species by Family')
                    ax2.invert_yaxis()
                    
                    ax3.axis('tight')
                    ax3.axis('off')
                    table_data = [[family, records, species] for family, records, species in 
                                 zip(df_top['family'], df_top['record_count'], df_top['species_count'])]
                    table = ax3.table(cellText=table_data, 
                                     colLabels=['Family', 'Records', 'Species'],
                                     cellLoc='left', loc='center')
                    table.auto_set_font_size(False)
                    table.set_fontsize(8)
                    table.scale(1, 1.5)
                    
                    total_records = combined['record_count'].sum()
                    total_species = combined['species_count'].sum()
                    unique_families = len(combined)
                    
                    summary_text = f"""
Family Summary ({order_name}):
• Unique families: {unique_families}
• Total records: {total_records:,}
• Total species: {total_species:,}
• Avg records per family: {total_records/unique_families:.1f}
• Avg species per family: {total_species/unique_families:.1f}
                    """
                    
                    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=10,
                            verticalalignment='top', fontfamily='monospace')
                    ax4.set_xlim(0, 1)
                    ax4.set_ylim(0, 1)
                    ax4.axis('off')
                    
                    plt.suptitle(f'Taxonomic Breakdown: Family ({order_name})', fontsize=16, fontweight='bold')
                    plt.tight_layout()
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close()
                    
                except Exception as e:
                    print(f"Error creating family breakdown for order {order_name}: {e}")
                    continue
                    
        except Exception as e:
            print(f"Error getting orders for family breakdown: {e}")
    def _create_basic_taxonomic_summary_page(self, pdf):
        """Create a single page with basic taxonomic summaries (faster alternative)"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 8.5))
        
        # Summary counts by taxonomic level
        levels = ['phylum', 'class', 'order', 'family']
        level_counts = []
        
        for level in levels:
            df = self.stats['taxonomic_breakdown'].get(level, pd.DataFrame())
            if not df.empty:
                level_counts.append({
                    'level': level.title(),
                    'unique_count': len(df),
                    'total_records': df['record_count'].sum(),
                    'total_species': df['species_count'].sum()
                })
        
        if level_counts:
            counts_df = pd.DataFrame(level_counts)
            
            # Bar chart of unique taxa by level
            bars1 = ax1.bar(counts_df['level'], counts_df['unique_count'])
            ax1.set_ylabel('Number of Unique Taxa')
            ax1.set_title('Taxonomic Diversity by Level')
            ax1.tick_params(axis='x', rotation=45)
            
            for bar in bars1:
                height = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2., height,
                        f'{int(height)}', ha='center', va='bottom')
            
            # Replace table with text summary
            ax2.axis('off')
            summary_text = "Taxonomic Summary Statistics:\n\n"
            for _, row in counts_df.iterrows():
                summary_text += f"{row['level']}:\n"
                summary_text += f"  • Unique taxa: {row['unique_count']}\n"
                summary_text += f"  • Total records: {row['total_records']:,}\n"
                summary_text += f"  • Total species: {row['total_species']:,}\n\n"
            
            ax2.text(0.1, 0.9, summary_text, transform=ax2.transAxes, fontsize=10,
                    verticalalignment='top', fontfamily='monospace')
            ax2.set_title('Taxonomic Summary Statistics')
        
        # Top 10 most diverse orders (by species count)
        order_df = self.stats['taxonomic_breakdown'].get('order', pd.DataFrame())
        if not order_df.empty:
            top_orders = order_df.nlargest(10, 'species_count')
            
            bars3 = ax3.barh(range(len(top_orders)), top_orders['species_count'])
            ax3.set_yticks(range(len(top_orders)))
            ax3.set_yticklabels(top_orders['order'], fontsize=8)
            ax3.set_xlabel('Number of Species')
            ax3.set_title('Top 10 Orders by Species Diversity')
            ax3.invert_yaxis()
        else:
            ax3.text(0.5, 0.5, 'No order data available', ha='center', va='center')
            ax3.set_title('Top 10 Orders by Species Diversity')
        
        # General taxonomic statistics
        total_orders = len(order_df) if not order_df.empty else 0
        family_df = self.stats['taxonomic_breakdown'].get('family', pd.DataFrame())
        total_families = len(family_df) if not family_df.empty else 0
        
        # Get top order info safely
        if not order_df.empty and 'top_orders' in locals():
            most_diverse_order = top_orders.iloc[0]['order']
            top_order_records = top_orders.iloc[0]['record_count']
        else:
            most_diverse_order = 'N/A'
            top_order_records = 'N/A'
        
        # Format the record count properly
        if isinstance(top_order_records, int):
            records_display = f"{top_order_records:,}"
        else:
            records_display = str(top_order_records)
        
        stats_text = f"""
Taxonomic Overview:

• Total Orders: {total_orders:,}
• Total Families: {total_families:,}
• Most diverse order: {most_diverse_order}
• Records in top order: {records_display}

Note: Detailed family-by-order 
breakdown was skipped for faster 
processing. Use without 
--skip-detailed-taxonomy for 
complete analysis.
        """
        
        ax4.text(0.1, 0.9, stats_text, transform=ax4.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace')
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        ax4.axis('off')
        
        plt.suptitle('Taxonomic Overview (Summary)', fontsize=16, fontweight='bold')
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
    def _create_performance_page(self, pdf):
        """Create performance metrics page"""
        fig, ax = plt.subplots(figsize=(8.5, 11))
        
        elapsed_time = time.time() - self.start_time
        
        performance_text = f"""
PERFORMANCE METRICS

Processing Configuration:
• HPC optimized mode: Enabled
• Fast mode: {self.fast_mode}
• Sample size: {self.sample_size:,} records (fast mode only)
• Max workers: {self.max_workers}
• Skip detailed taxonomy: {self.skip_detailed_taxonomy}

Resource Usage:
• Total processing time: {elapsed_time/60:.1f} minutes
• Peak memory usage: {self.memory_monitor.peak_memory:.1f} MB
• Current memory usage: {self.memory_monitor.get_memory_usage():.1f} MB

Database Information:
• Database file: {Path(self.db_path).name}
• File size: {self._get_file_size()}
• SQLite optimizations: Enabled
• Performance indexes: {'Created' if not self.fast_mode else 'Skipped (fast mode)'}

System Information:
• Available CPU cores: {psutil.cpu_count()}
• Total system memory: {psutil.virtual_memory().total / (1024**3):.1f} GB
• Available memory: {psutil.virtual_memory().available / (1024**3):.1f} GB

Processing Rates:
• Records per minute: {(self.stats.get('total_records', 0) / elapsed_time * 60):,.0f}
• Memory efficiency: {(self.stats.get('total_records', 0) / self.memory_monitor.peak_memory):,.0f} records/MB
        """
        
        ax.text(0.1, 0.9, performance_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.title('Performance Metrics', fontsize=16, fontweight='bold', pad=20)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    # Helper methods
    def _get_most_common_bags_grade(self):
        """Get most common BAGS grade"""
        if not self.stats.get('bags_distribution', pd.DataFrame()).empty:
            df = self.stats['bags_distribution']
            most_common = df.loc[df['count'].idxmax(), 'bags_grade']
            return f"{most_common}"
        return "N/A"
    
    def _get_most_common_rank(self):
        """Get most common rank"""
        if not self.stats.get('rank_distribution', pd.DataFrame()).empty:
            df = self.stats['rank_distribution']
            most_common = df.loc[df['count'].idxmax(), 'ranking']
            return f"{most_common}"
        return "N/A"
    
    def _get_file_size(self):
        """Get database file size"""
        try:
            size_bytes = Path(self.db_path).stat().st_size
            size_gb = size_bytes / (1024 * 1024 * 1024)
            return f"{size_gb:.2f} GB"
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
            if self.fast_mode:
                query = f"""
                SELECT SUM(CASE WHEN status = 1 THEN 1 ELSE 0 END) * 100.0 / COUNT(*) as pass_rate
                FROM (SELECT * FROM bold_criteria ORDER BY RANDOM() LIMIT {self.sample_size}) bc
                """
            else:
                query = """
                SELECT SUM(CASE WHEN status = 1 THEN 1 ELSE 0 END) * 100.0 / COUNT(*) as pass_rate
                FROM bold_criteria
                """
            result = pd.read_sql_query(query, self.conn)
            return f"{result.iloc[0, 0]:.1f}%"
        except:
            return "N/A"


def main():
    parser = argparse.ArgumentParser(description='Generate HPC-optimized BOLD database statistics report')
    parser.add_argument('-d', '--database', required=True, help='Path to BOLD database file')
    parser.add_argument('-l', '--log', help='Path to prescoring_filter.log file (optional)')
    parser.add_argument('-o', '--output', default='bold_stats_report.pdf', help='Output PDF file path')
    parser.add_argument('--info', action='store_true', help='Show database table information and exit')
    parser.add_argument('--skip-detailed-taxonomy', action='store_true', 
                       help='Skip detailed family-by-order taxonomy pages (faster execution)')
    parser.add_argument('--fast-mode', action='store_true',
                       help='Use fast mode with sampling for very large datasets')
    parser.add_argument('--sample-size', type=int, default=100000,
                       help='Sample size for fast mode (default: 100000)')
    parser.add_argument('--max-workers', type=int, default=4,
                       help='Maximum number of worker threads (default: 4)')
    
    args = parser.parse_args()
    
    if not Path(args.database).exists():
        print(f"Error: Database file not found: {args.database}")
        return 1
    
    try:
        with HPCBOLDStatsGenerator(
            args.database, 
            args.log, 
            args.skip_detailed_taxonomy,
            args.fast_mode,
            args.sample_size,
            args.max_workers
        ) as generator:
            
            if args.info:
                generator.get_table_info()
                return 0
            
            print("=== HPC-Optimized BOLD Database Statistics Report ===")
            print(f"Database: {args.database}")
            print(f"Output: {args.output}")
            
            generator.generate_all_stats_parallel()
            generator.create_pdf_report_optimized(args.output)
            
            return 0
            
    except Exception as e:
        print(f"Error generating report: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
