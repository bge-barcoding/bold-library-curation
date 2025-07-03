#!/usr/bin/env python3

"""
BOLD Family Database Splitter (Python Version)

Splits BOLD database or TSV export into family-level databases
organized by taxonomic hierarchy with configurable size thresholds
"""

import os
import sys
import sqlite3
import csv
import argparse
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime
import re

# Configuration
CONFIG = {
    # Thresholds
    'FAMILY_SIZE_THRESHOLD': 10000,  # Split families larger than this by subfamily
    
    # Fallback values for missing taxonomy
    'UNKNOWN_PHYLUM': 'Unknown_Phylum',
    'UNKNOWN_CLASS': 'Unknown_Class', 
    'UNKNOWN_ORDER': 'Unknown_Order',
    'UNKNOWN_FAMILY': 'Unknown_Family',
    'UNKNOWN_SUBFAMILY': 'Unknown_Subfamily'
}

class BoldFamilySplitter:
    def __init__(self, input_path, output_base='taxonomic_output', threshold=10000):
        self.input_path = Path(input_path)
        self.output_base = Path(output_base)
        self.threshold = threshold
        
        # Handle directory input - look for files in the directory
        if self.input_path.is_dir():
            print(f"Directory provided: {self.input_path}")
            self.input_path, self.is_database = self._resolve_input_file(self.input_path / "dummy")
        else:
            # Resolve specific file input with new priority logic
            self.input_path, self.is_database = self._resolve_input_file(input_path)
        
        self.family_stats = {}
        self.taxonomy_cache = {}
        
        # Ensure output directory exists
        self.output_base.mkdir(parents=True, exist_ok=True)
        
        # Update threshold in config
        CONFIG['FAMILY_SIZE_THRESHOLD'] = threshold

    def _resolve_input_file(self, input_path):
        """
        Resolve input file with priority: bold.db → any other .db → result_output.tsv → original input
        Returns tuple of (resolved_path, is_database)
        """
        input_path = Path(input_path)
        input_dir = input_path.parent
        
        print(f"Resolving input file from: {input_path}")
        
        # Priority 1: Look for bold.db in the same directory
        bold_db_path = input_dir / "bold.db"
        if bold_db_path.exists():
            print(f"Found bold.db, using: {bold_db_path}")
            return bold_db_path, True
        
        # Priority 2: Look for any other .db file in the same directory
        db_files = list(input_dir.glob("*.db"))
        if db_files:
            # Use the first .db file found (excluding bold.db which we already checked)
            other_db = db_files[0]
            print(f"Found other database file, using: {other_db}")
            return other_db, True
        
        # Priority 3: Look for result_output.tsv in the same directory
        result_output_path = input_dir / "result_output.tsv"
        if result_output_path.exists():
            print(f"Found result_output.tsv, using: {result_output_path}")
            return result_output_path, False
        
        # Priority 4: If input is explicitly a database and exists, use it
        if str(input_path).endswith('.db') and input_path.exists():
            print(f"Using specified database file: {input_path}")
            return input_path, True
        
        # Priority 5: If input is a TSV and exists, use it
        if (str(input_path).endswith('.tsv') or str(input_path).endswith('.txt')) and input_path.exists():
            print(f"Using specified TSV file: {input_path}")
            return input_path, False
        
        # Priority 6: Try to detect file type if it exists without clear extension
        if input_path.exists():
            # Check if it's actually a database by trying to open it
            try:
                with sqlite3.connect(input_path) as conn:
                    cursor = conn.cursor()
                    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1")
                    cursor.fetchone()
                print(f"Detected as database file: {input_path}")
                return input_path, True
            except:
                print(f"Treating as TSV file: {input_path}")
                return input_path, False
        
        # Priority 7: Look for alternatives based on input base name
        base_path = input_path.with_suffix('')
        
        # Try .db variants
        for db_name in ['bold.db', f'{base_path.name}.db']:
            db_path = input_dir / db_name
            if db_path.exists():
                print(f"Found alternative database file: {db_path}")
                return db_path, True
        
        # Fall back to TSV variants
        for ext in ['.tsv', '.txt', '.csv']:
            tsv_path = base_path.with_suffix(ext)
            if tsv_path.exists():
                print(f"Found alternative TSV file: {tsv_path}")
                return tsv_path, False
        
        # If nothing found, return original path (will cause error later)
        print(f"WARNING: No suitable file found in directory {input_dir}")
        print(f"Falling back to original input: {input_path}")
        return input_path, str(input_path).endswith('.db')

    def run(self):
        """Main execution method"""
        print(f"Processing {'database' if self.is_database else 'TSV'}: {self.input_path}")
        print(f"Output directory: {self.output_base}")
        print(f"Family size threshold: {self.threshold} records")
        
        try:
            # Step 1: Analyze family sizes and get taxonomy info
            self.analyze_families()
            
            # Step 2: Create family databases with proper folder structure
            self.create_family_databases()
            
            # Step 3: Generate summary report
            self.generate_report()
            
            print('\n*** Processing completed successfully!')
            
        except Exception as error:
            print(f'ERROR: {error}')
            raise error

    def analyze_families(self):
        """Analyze family sizes and collect taxonomy information"""
        print('\nAnalyzing family sizes and taxonomy...')
        
        if self.is_database:
            self._analyze_families_from_db()
        else:
            self._analyze_families_from_tsv()
        
        print('\nFamily analysis results:')
        for family, stats in self.family_stats.items():
            strategy = 'subfamily' if stats['count'] > self.threshold else 'family'
            print(f"  {family}: {stats['count']} records (split by {strategy})")

    def _analyze_families_from_db(self):
        """Analyze families from SQLite database"""
        with sqlite3.connect(self.input_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()
            
            # Get family counts
            cursor.execute(f"""
                SELECT 
                    COALESCE(family, '{CONFIG['UNKNOWN_FAMILY']}') as family,
                    COUNT(*) as count
                FROM bold 
                GROUP BY family
                ORDER BY count DESC
            """)
            
            family_counts = cursor.fetchall()
            
            # Get taxonomy information for each family
            for row in family_counts:
                family = row['family']
                
                # Get sample taxonomy for this family
                cursor.execute(f"""
                    SELECT 
                        COALESCE(phylum, '{CONFIG['UNKNOWN_PHYLUM']}') as phylum,
                        COALESCE(class, '{CONFIG['UNKNOWN_CLASS']}') as class,
                        COALESCE("order", '{CONFIG['UNKNOWN_ORDER']}') as order_name,
                        COALESCE(family, '{CONFIG['UNKNOWN_FAMILY']}') as family,
                        subfamily
                    FROM bold 
                    WHERE COALESCE(family, '{CONFIG['UNKNOWN_FAMILY']}') = ?
                    LIMIT 1
                """, (family,))
                
                taxonomy_row = cursor.fetchone()
                taxonomy = dict(taxonomy_row) if taxonomy_row else {}
                
                # Get subfamily counts if family is large
                subfamilies = []
                if row['count'] > self.threshold:
                    cursor.execute(f"""
                        SELECT 
                            COALESCE(subfamily, '{CONFIG['UNKNOWN_SUBFAMILY']}') as subfamily,
                            COUNT(*) as count
                        FROM bold 
                        WHERE COALESCE(family, '{CONFIG['UNKNOWN_FAMILY']}') = ?
                        GROUP BY subfamily
                        ORDER BY count DESC
                    """, (family,))
                    
                    subfamilies = [dict(row) for row in cursor.fetchall()]
                
                self.family_stats[family] = {
                    'count': row['count'],
                    'taxonomy': {
                        'phylum': taxonomy.get('phylum', CONFIG['UNKNOWN_PHYLUM']),
                        'class': taxonomy.get('class', CONFIG['UNKNOWN_CLASS']),
                        'order': taxonomy.get('order_name', CONFIG['UNKNOWN_ORDER']),
                        'family': family
                    },
                    'subfamilies': subfamilies
                }

    def _analyze_families_from_tsv(self):
        """Analyze families from TSV file"""
        family_counts = Counter()
        taxonomy_examples = {}
        
        # First pass: count families and collect taxonomy examples
        with open(self.input_path, 'r', encoding='utf-8', errors='replace') as file:
            reader = csv.DictReader(file, delimiter='\t')
            
            # Find relevant columns (case-insensitive)
            headers = [h.lower() for h in reader.fieldnames]
            family_col = self._find_column(reader.fieldnames, ['family'])
            phylum_col = self._find_column(reader.fieldnames, ['phylum'])
            class_col = self._find_column(reader.fieldnames, ['class'])
            order_col = self._find_column(reader.fieldnames, ['order'])
            subfamily_col = self._find_column(reader.fieldnames, ['subfamily'])
            
            if family_col is None:
                raise ValueError('Family column not found in TSV')
            
            line_count = 0
            for row in reader:
                family = row.get(family_col, CONFIG['UNKNOWN_FAMILY']) or CONFIG['UNKNOWN_FAMILY']
                
                # Count families
                family_counts[family] += 1
                
                # Store taxonomy example for each family
                if family not in taxonomy_examples:
                    taxonomy_examples[family] = {
                        'phylum': (row.get(phylum_col) if phylum_col else None) or CONFIG['UNKNOWN_PHYLUM'],
                        'class': (row.get(class_col) if class_col else None) or CONFIG['UNKNOWN_CLASS'],
                        'order': (row.get(order_col) if order_col else None) or CONFIG['UNKNOWN_ORDER'],
                        'subfamily': row.get(subfamily_col) if subfamily_col else None
                    }
                
                line_count += 1
                if line_count % 100000 == 0:
                    print(f"\rProcessed {line_count} records...", end='', flush=True)
        
        print(f"\nProcessed {line_count} total records")
        
        # Analyze subfamily distribution for large families
        for family, count in family_counts.items():
            taxonomy = taxonomy_examples[family]
            subfamilies = []
            
            if count > self.threshold:
                # Re-scan to get subfamily counts for this family
                subfamilies = self._get_subfamily_counts_from_tsv(family, family_col, subfamily_col)
            
            self.family_stats[family] = {
                'count': count,
                'taxonomy': taxonomy,
                'subfamilies': subfamilies
            }

    def _get_subfamily_counts_from_tsv(self, target_family, family_col, subfamily_col):
        """Get subfamily counts for a specific family from TSV"""
        if subfamily_col is None:
            return []
        
        subfamily_counts = Counter()
        
        with open(self.input_path, 'r', encoding='utf-8', errors='replace') as file:
            reader = csv.DictReader(file, delimiter='\t')
            
            for row in reader:
                family = row.get(family_col, CONFIG['UNKNOWN_FAMILY']) or CONFIG['UNKNOWN_FAMILY']
                
                if family == target_family:
                    subfamily = row.get(subfamily_col, CONFIG['UNKNOWN_SUBFAMILY']) or CONFIG['UNKNOWN_SUBFAMILY']
                    subfamily_counts[subfamily] += 1
        
        return [{'subfamily': subfamily, 'count': count} 
                for subfamily, count in subfamily_counts.items()]

    def _find_column(self, fieldnames, possible_names):
        """Find column by possible names (case-insensitive)"""
        for name in possible_names:
            for field in fieldnames:
                if field.lower() == name.lower():
                    return field
        return None

    def create_family_databases(self):
        """Create family-level databases with taxonomic organization"""
        print('\nCreating family databases...')
        
        for family, stats in self.family_stats.items():
            taxonomy = stats['taxonomy']
            subfamilies = stats['subfamilies']
            
            # Create taxonomic folder structure
            family_dir = self.output_base / self._sanitize_filename(taxonomy['phylum']) / \
                        self._sanitize_filename(taxonomy['class']) / \
                        self._sanitize_filename(taxonomy['order'])
            
            family_dir.mkdir(parents=True, exist_ok=True)
            
            if subfamilies:
                # Split by subfamily
                print(f"  Creating subfamily databases for {family}...")
                self._create_subfamily_databases(family, family_dir, subfamilies)
            else:
                # Create single family database
                db_path = family_dir / f"{self._sanitize_filename(family)}.db"
                self._create_single_family_database(family, db_path)
                print(f"  *** Created {db_path}")

    def _create_subfamily_databases(self, family, family_dir, subfamilies):
        """Create separate databases for each subfamily"""
        subfamily_dir = family_dir / self._sanitize_filename(family)
        subfamily_dir.mkdir(parents=True, exist_ok=True)
        
        for subfamily_info in subfamilies:
            subfamily = subfamily_info['subfamily']
            db_path = subfamily_dir / f"{self._sanitize_filename(subfamily)}.db"
            
            if self.is_database:
                self._create_subfamily_database_from_db(family, subfamily, db_path)
            else:
                self._create_subfamily_database_from_tsv(family, subfamily, db_path)
            
            print(f"    *** Created {db_path} ({subfamily_info['count']} records)")
        
        # Create shared BINs analysis for this family
        self._analyze_shared_bins(family, subfamily_dir)

    def _create_single_family_database(self, family, db_path):
        """Create a single database for a family"""
        if self.is_database:
            self._create_family_database_from_db(family, db_path)
        else:
            self._create_family_database_from_tsv(family, db_path)

    def _create_family_database_from_db(self, family, db_path):
        """Create family database from source database"""
        # Connect to source and target databases
        source_conn = sqlite3.connect(self.input_path)
        target_conn = sqlite3.connect(db_path)
        
        try:
            # Copy schema
            source_cursor = source_conn.cursor()
            source_cursor.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='bold'")
            schema_row = source_cursor.fetchone()
            
            if schema_row:
                target_conn.execute(schema_row[0])
            
            # Copy data
            source_cursor.execute(f"""
                SELECT * FROM bold 
                WHERE COALESCE(family, '{CONFIG['UNKNOWN_FAMILY']}') = ?
            """, (family,))
            
            # Insert data into target database
            target_cursor = target_conn.cursor()
            rows = source_cursor.fetchall()
            
            if rows:
                # Get column info to build insert statement
                source_cursor.execute("PRAGMA table_info(bold)")
                columns = [row[1] for row in source_cursor.fetchall()]
                placeholders = ','.join(['?' for _ in columns])
                
                target_cursor.executemany(f"INSERT INTO bold VALUES ({placeholders})", rows)
            
            # Create indexes
            self._create_indexes(target_conn)
            target_conn.commit()
            
        finally:
            source_conn.close()
            target_conn.close()

    def _create_subfamily_database_from_db(self, family, subfamily, db_path):
        """Create subfamily database from source database"""
        source_conn = sqlite3.connect(self.input_path)
        target_conn = sqlite3.connect(db_path)
        
        try:
            # Copy schema
            source_cursor = source_conn.cursor()
            source_cursor.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='bold'")
            schema_row = source_cursor.fetchone()
            
            if schema_row:
                target_conn.execute(schema_row[0])
            
            # Copy data
            source_cursor.execute(f"""
                SELECT * FROM bold 
                WHERE COALESCE(family, '{CONFIG['UNKNOWN_FAMILY']}') = ? 
                AND COALESCE(subfamily, '{CONFIG['UNKNOWN_SUBFAMILY']}') = ?
            """, (family, subfamily))
            
            # Insert data into target database
            target_cursor = target_conn.cursor()
            rows = source_cursor.fetchall()
            
            if rows:
                # Get column info to build insert statement
                source_cursor.execute("PRAGMA table_info(bold)")
                columns = [row[1] for row in source_cursor.fetchall()]
                placeholders = ','.join(['?' for _ in columns])
                
                target_cursor.executemany(f"INSERT INTO bold VALUES ({placeholders})", rows)
            
            # Create indexes
            self._create_indexes(target_conn)
            target_conn.commit()
            
        finally:
            source_conn.close()
            target_conn.close()

    def _create_family_database_from_tsv(self, family, db_path):
        """Create family database from TSV file"""
        with sqlite3.connect(db_path) as conn:
            # Create table structure
            self._create_table_from_tsv(conn)
            
            # Insert data
            def filter_func(row, family_col, subfamily_col=None):
                record_family = row.get(family_col, CONFIG['UNKNOWN_FAMILY']) or CONFIG['UNKNOWN_FAMILY']
                return record_family == family
            
            self._insert_tsv_data(conn, filter_func)
            
            # Create indexes
            self._create_indexes(conn)

    def _create_subfamily_database_from_tsv(self, family, subfamily, db_path):
        """Create subfamily database from TSV file"""
        with sqlite3.connect(db_path) as conn:
            # Create table structure
            self._create_table_from_tsv(conn)
            
            # Insert data
            def filter_func(row, family_col, subfamily_col=None):
                record_family = row.get(family_col, CONFIG['UNKNOWN_FAMILY']) or CONFIG['UNKNOWN_FAMILY']
                record_subfamily = row.get(subfamily_col, CONFIG['UNKNOWN_SUBFAMILY']) or CONFIG['UNKNOWN_SUBFAMILY'] if subfamily_col else CONFIG['UNKNOWN_SUBFAMILY']
                return record_family == family and record_subfamily == subfamily
            
            self._insert_tsv_data(conn, filter_func)
            
            # Create indexes
            self._create_indexes(conn)

    def _create_table_from_tsv(self, conn):
        """Create table structure based on TSV headers"""
        with open(self.input_path, 'r', encoding='utf-8', errors='replace') as file:
            reader = csv.DictReader(file, delimiter='\t')
            headers = reader.fieldnames
        
        # Handle duplicate column names by adding suffixes
        seen_columns = {}
        unique_headers = []
        
        for header in headers:
            clean_header = header.replace('"', '""')
            if clean_header.lower() in seen_columns:
                # Add suffix to make it unique
                suffix = seen_columns[clean_header.lower()] + 1
                seen_columns[clean_header.lower()] = suffix
                unique_header = f"{clean_header}_v{suffix}"
            else:
                seen_columns[clean_header.lower()] = 0
                unique_header = clean_header
            unique_headers.append(unique_header)
        
        # Create table with all columns as TEXT
        columns_sql = ',\n'.join([f'"{header}" TEXT' for header in unique_headers])
        create_sql = f"CREATE TABLE bold (\n{columns_sql}\n)"
        
        conn.execute(create_sql)
        
        # Store the header mapping for later use
        self._header_mapping = dict(zip(headers, unique_headers))

    def _insert_tsv_data(self, conn, filter_function):
        """Insert filtered TSV data into database"""
        with open(self.input_path, 'r', encoding='utf-8', errors='replace') as file:
            reader = csv.DictReader(file, delimiter='\t')
            headers = reader.fieldnames
            
            # Find column indices for filtering
            family_col = self._find_column(headers, ['family'])
            subfamily_col = self._find_column(headers, ['subfamily'])
            
            # Use mapped headers if available, otherwise original headers
            if hasattr(self, '_header_mapping'):
                mapped_headers = [self._header_mapping[h] for h in headers]
            else:
                mapped_headers = headers
            
            # Prepare insert statement
            placeholders = ','.join(['?' for _ in mapped_headers])
            insert_sql = f"INSERT INTO bold VALUES ({placeholders})"
            
            # Insert data in batches
            batch_size = 10000
            batch = []
            inserted_count = 0
            
            conn.execute('BEGIN TRANSACTION')
            
            for row in reader:
                if filter_function(row, family_col, subfamily_col):
                    values = [row.get(header, '') for header in headers]
                    batch.append(values)
                    inserted_count += 1
                    
                    if len(batch) >= batch_size:
                        conn.executemany(insert_sql, batch)
                        batch = []
            
            # Insert remaining records
            if batch:
                conn.executemany(insert_sql, batch)
            
            conn.execute('COMMIT')
            return inserted_count

    def _analyze_shared_bins(self, family, output_dir):
        """Analyze shared BINs across subfamilies"""
        print(f"    > Analyzing shared BINs for {family}...")
        
        if not self.is_database:
            print(f"    > Skipping BIN analysis for TSV input")
            return
        
        with sqlite3.connect(self.input_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()
            
            cursor.execute(f"""
                SELECT 
                    bin_uri,
                    GROUP_CONCAT(DISTINCT COALESCE(subfamily, '{CONFIG['UNKNOWN_SUBFAMILY']}')) as subfamilies,
                    COUNT(DISTINCT COALESCE(subfamily, '{CONFIG['UNKNOWN_SUBFAMILY']}')) as subfamily_count
                FROM bold 
                WHERE COALESCE(family, '{CONFIG['UNKNOWN_FAMILY']}') = ? 
                    AND bin_uri IS NOT NULL 
                    AND bin_uri != 'None' 
                    AND bin_uri != ''
                GROUP BY bin_uri 
                HAVING subfamily_count > 1
            """, (family,))
            
            shared_bins = cursor.fetchall()
        
        if shared_bins:
            csv_path = output_dir / 'subfamily_shared_BINs.csv'
            with open(csv_path, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=';')
                for row in shared_bins:
                    writer.writerow([row['bin_uri'], row['subfamilies'], family])
            
            print(f"    > Saved {len(shared_bins)} shared BINs to {csv_path}")

    def _create_indexes(self, conn):
        """Create database indexes for performance"""
        indexes = [
            'CREATE INDEX IF NOT EXISTS idx_family ON bold(family)',
            'CREATE INDEX IF NOT EXISTS idx_genus ON bold(genus)',
            'CREATE INDEX IF NOT EXISTS idx_species ON bold(species)',
            'CREATE INDEX IF NOT EXISTS idx_bin ON bold(bin_uri)',
            'CREATE INDEX IF NOT EXISTS idx_processid ON bold(processid)'
        ]
        
        for index_sql in indexes:
            try:
                conn.execute(index_sql)
            except sqlite3.OperationalError:
                # Ignore errors for missing columns
                pass

    def _sanitize_filename(self, name):
        """Sanitize filename for filesystem compatibility"""
        if not name or name in ('None', ''):
            return 'Unknown'
        # Replace problematic characters with underscores
        sanitized = re.sub(r'[<>:"/\\|?*]', '_', str(name))
        sanitized = re.sub(r'\s+', '_', sanitized)
        return sanitized

    def generate_report(self):
        """Generate comprehensive splitting report"""
        report_path = self.output_base / 'splitting_report.txt'
        
        lines = []
        lines.append('BOLD Database Splitting Report')
        lines.append('=' * 50)
        lines.append(f'Input: {self.input_path}')
        lines.append(f'Output: {self.output_base}')
        lines.append(f'Threshold: {self.threshold} records')
        lines.append(f'Generated: {datetime.now().isoformat()}')
        lines.append('')
        
        lines.append('Family Statistics:')
        lines.append('-' * 50)
        
        for family, stats in self.family_stats.items():
            strategy = 'subfamily' if stats['subfamilies'] else 'family'
            lines.append(f'{family}: {stats["count"]} records (split by {strategy})')
            
            if stats['subfamilies']:
                for subfam in stats['subfamilies']:
                    lines.append(f'  +-- {subfam["subfamily"]}: {subfam["count"]} records')
        
        lines.append('')
        lines.append('Output Structure:')
        lines.append('-' * 50)
        
        total_families = len(self.family_stats)
        split_families = sum(1 for stats in self.family_stats.values() if stats['subfamilies'])
        
        lines.append(f'Total families processed: {total_families}')
        lines.append(f'Families split by subfamily: {split_families}')
        lines.append(f'Regular family databases: {total_families - split_families}')
        
        with open(report_path, 'w') as f:
            f.write('\n'.join(lines))
        
        print(f'\nReport saved to: {report_path}')


def main():
    """Main command-line interface"""
    parser = argparse.ArgumentParser(
        description='BOLD Family Database Splitter (Python Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python bold_family_splitter.py bold.db
  python bold_family_splitter.py result_output.tsv --output my_output --threshold 5000
  python bold_family_splitter.py "C:\\GitHub\\bold-library-curation\\results\\results_test3"
  python bold_family_splitter.py "C:\\GitHub\\bold-library-curation\\results\\_bags_sub-3\\bold.db"
  python bold_family_splitter.py data_file  # Will auto-detect .db or .tsv

File Selection Priority:
  1. bold.db in the same directory (highest priority)
  2. Any other .db file in the same directory  
  3. result_output.tsv in the same directory
  4. The explicitly specified file (if it exists)
  5. Alternative files with same base name
  - Database files are preferred for better performance with large datasets

Output Structure:
  output_dir/
  +-- Phylum1/
      +-- Class1/
          +-- Order1/
              +-- SmallFamily.db
              +-- LargeFamily/
                  +-- Subfamily1.db
                  +-- Subfamily2.db
                  +-- subfamily_shared_BINs.csv
        """)
    
    parser.add_argument('input', help='Path to input database (.db), TSV file (.tsv), or directory containing data files. When a directory is provided, files are auto-resolved by priority.')
    parser.add_argument('--output', '-o', default='taxonomic_output', 
                       help='Output directory (default: taxonomic_output)')
    parser.add_argument('--threshold', '-t', type=int, default=10000,
                       help='Family size threshold for subfamily splitting (default: 10000)')
    
    args = parser.parse_args()
    
    # Validate input
    if not os.path.exists(args.input):
        print(f'ERROR: Input file not found: {args.input}')
        sys.exit(1)
    
    try:
        splitter = BoldFamilySplitter(args.input, args.output, args.threshold)
        splitter.run()
    except Exception as error:
        print(f'ERROR: Processing failed: {error}')
        sys.exit(1)


if __name__ == '__main__':
    main()
