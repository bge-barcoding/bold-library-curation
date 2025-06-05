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
        
        if self.fast_mode:
            # Use approximate counts for very large datasets
            queries = {
                'total_records': "SELECT COUNT(*) FROM bold",
                'species_count': "SELECT COUNT(DISTINCT species) FROM bold WHERE species IS NOT NULL",
                'bins_count': "SELECT COUNT(DISTINCT bin_uri) FROM bold WHERE bin_uri IS NOT NULL AND bin_uri != ''",
            }
            
            # Get approximate counts for OTUs and haplotypes
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
        else:
            queries = {
                'total_records': "SELECT COUNT(*) FROM bold",
                'species_count': "SELECT COUNT(DISTINCT species) FROM bold WHERE species IS NOT NULL",
                'bins_count': "SELECT COUNT(DISTINCT bin_uri) FROM bold WHERE bin_uri IS NOT NULL AND bin_uri != ''",
                'otus_count': "SELECT COUNT(DISTINCT otu_id) FROM bold_otus",
                'haplotypes_count': "SELECT COUNT(DISTINCT haplotype_id) FROM bold_haplotypes",
            }
        
        for key, query in queries.items():
            try:
                result = pd.read_sql_query(query, self.conn)
                self.stats[key] = result.iloc[0, 0]
            except Exception as e:
                print(f"Error executing query for {key}: {e}")
                self.stats[key] = 0
        
        # Calculate records without BIN
        self.stats['records_without_bin'] = self.stats['total_records'] - self.stats['bins_count']
        
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
                # Use sample for very large datasets
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
        FROM bold 
        JOIN bold_ranks br ON bold.recordid = br.recordid 
        WHERE br.ranking IS NOT NULL 
        GROUP BY br.ranking 
        ORDER BY br.ranking
        """
        try:
            if self.fast_mode:
                # Use sample for very large datasets
                sample_query = """
                SELECT br.ranking, COUNT(*) as count 
                FROM (SELECT * FROM bold ORDER BY RANDOM() LIMIT {}) bold
                JOIN bold_ranks br ON bold.recordid = br.recordid 
                WHERE br.ranking IS NOT NULL 
                GROUP BY br.ranking 
                ORDER BY br.ranking
                """.format(self.sample_size)
                df = pd.read_sql_query(sample_query, self.conn)
            else:
                df = pd.read_sql_query(query, self.conn)
            self.stats['rank_distribution'] = df
        except Exception as e:
            print(f"Error getting rank distribution: {e}")
            self.stats['rank_distribution'] = pd.DataFrame()
        
        self.memory_monitor.log_memory("rank_distribution")
    
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
        
        # Simple additional tasks
        print("\n=== Phase 2: Distribution Analysis ===")
        
        simple_tasks = [
            ("BAGS distribution", self.get_bags_distribution),
            ("Rank distribution", self.get_rank_distribution)
        ]
        
        # Run tasks sequentially for now (can be parallelized later)
        for task_name, task_func in simple_tasks:
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
        """Create optimized PDF report for large datasets"""
        print(f"\n=== Creating PDF Report: {output_path} ===")
        
        try:
            with PdfPages(output_path) as pdf:
                # Title page
                self._create_title_page(pdf)
                
                # Summary statistics
                self._create_summary_page(pdf)
                
                # Core distribution pages
                if not self.stats['bags_distribution'].empty:
                    self._create_bags_page(pdf)
                
                if not self.stats['rank_distribution'].empty:
                    self._create_rank_page(pdf)
                
                # Performance and system info
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
        
        summary_text = f"""
SUMMARY STATISTICS {mode_note}

Dataset Information:
Records processed: {self.stats.get('records_processed', 'N/A'):,}
Total records: {self.stats.get('total_records', 0):,}
Unique species: {self.stats.get('species_count', 0):,}
Unique BINs: {self.stats.get('bins_count', 0):,}
Unique OTUs: {self.stats.get('otus_count', 0):,}
Unique Haplotypes: {self.stats.get('haplotypes_count', 0):,}
Records without BIN: {self.stats.get('records_without_bin', 0):,}

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
    
    def _get_file_size(self):
        """Get database file size"""
        try:
            size_bytes = Path(self.db_path).stat().st_size
            size_gb = size_bytes / (1024 * 1024 * 1024)
            return f"{size_gb:.2f} GB"
        except:
            return "N/A"


def main():
    parser = argparse.ArgumentParser(description='Generate HPC-optimized BOLD database statistics report')
    parser.add_argument('-d', '--database', required=True, help='Path to BOLD database file')
    parser.add_argument('-l', '--log', help='Path to prescoring_filter.log file (optional)')
    parser.add_argument('-o', '--output', default='bold_stats_report.pdf', help='Output PDF file path')
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
