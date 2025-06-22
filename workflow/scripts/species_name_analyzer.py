#!/usr/bin/env python3
"""
Species Name Analyzer for BOLD BAGS Dataset
Handles | and , separators for BIN groups and species names
"""

import pandas as pd
import re
import logging
from datetime import datetime
from pathlib import Path

# Set up logging
log_filename = f"species_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class SpeciesNameAnalyzer:
    def __init__(self):
        """Initialize the analyzer with regex patterns for non-binomial names"""
        
        # Patterns that indicate NON-proper species names
        self.non_binomial_patterns = [
            # Basic "sp." patterns
            r'\bsp\.\s*$',                          # "Genus sp."
            r'\bsp\.\s+\d+',                        # "Genus sp. 1", "Genus sp. 22"
            r'\bsp\.\s+[A-Z]+\d*',                  # "Genus sp. ABC", "Genus sp. MS5619"
            r'\bsp\.\s+[a-z]+\d*',                  # "Genus sp. morph1", "sp. code"
            r'\bsp\.\s+nov\.',                      # "sp. nov."
            
            # Geographic/collection codes  
            r'\bsp\.\s+[A-Z]{2,}_\d+',             # "sp. USA_72", "sp. VIE4"
            r'\bsp\.\s+[A-Z]{3,}\d+',              # "sp. USCM5903"
            r'\bsp\.\s+[A-Z]+\d+[A-Z]*',           # "sp. SAA15", "sp. CHR116"
            
            # Morphotype designations
            r'\bsp\.\s+[A-Z]*[Mm]orph\d*',         # "sp. Morph0046", "sp. SAEVG Morph0118"
            
            # cf./aff. patterns
            r'\bcf\.\s+',                           # "Genus cf. species"
            r'\baff\.\s+',                          # "Genus aff. species"
            r'\bnr\.\s+',                           # "Genus nr. species" (near)
            
            # Complex/group designations
            r'\bcomplex\s*$',                       # "Genus species complex"
            r'\bgroup\s*$',                         # "Genus species group"
            r'\bss\.\s*$',                          # "species ss." (sensu stricto)
            r'-group$',                             # "species-group"
            
            # Hybrid designations
            r'\sx\s',                               # "Species x Species"
            r'\shybrid\b',                          # "Species hybrid"
            r'\striploid\s+hybrid',                 # "Species triploid hybrid"
            r'\stetraploid\s+hybrid',               # "Species tetraploid hybrid"
            
            # Uncertain designations
            r'\?\s*$',                              # "Genus species?"
            
            # Code patterns at end (fixed - only match actual codes)
            r'\s+[A-Z]{2,}\d+[A-Z]*$',             # " ABC123", " CHR116" (must have digits)
            r'\s+\d+[A-Z]+\d*$',                    # " 1A", " 22Ko4"
            r'\s+[A-Z]+_\d+$',                      # " USA_72"
            
            # Numbered forms
            r'\s+\d+\s*$',                          # " 1", " 22" at end
            r'_\d+$',                               # "_1", "_22" at end
            
            # Special collection codes
            r'BOLD:',                               # BOLD identifiers shouldn't be in species names
            r'USNM\s+IZ',                          # Museum collection codes
            r'ZMBN\s+\d+',                         # Museum codes
            r'NTNUVM\s+\d+',                       # Collection codes
            
            # Other non-binomial indicators
            r'\bssp\.\s+\d+',                      # "ssp. 1" - subspecies numbers
            r'\[.*\]',                             # Things in brackets
        ]
        
        # Compile patterns for efficiency
        self.compiled_patterns = [re.compile(pattern, re.IGNORECASE) for pattern in self.non_binomial_patterns]
        
        # Statistics tracking
        self.stats = {
            'total_entries': 0,
            'entries_with_sharers': 0,
            'total_species_names': 0,
            'proper_species_names': 0,
            'non_binomial_names': 0,
            'tautonyms_found': 0,
            'pattern_matches': {},
            'proper_examples': [],
            'non_binomial_examples': [],
            'tautonym_examples': [],
            'empty_bin_groups': 0
        }
    
    def is_proper_species_name(self, name):
        """
        Check if a species name is a proper binomial name
        Returns (is_proper, reason_if_not_proper)
        """
        if not name or pd.isna(name):
            return False, "Empty/NaN"
        
        name = name.strip()
        
        if not name:  # Empty after stripping
            return False, "Empty after stripping"
        
        # Check for non-binomial patterns first
        for i, pattern in enumerate(self.compiled_patterns):
            if pattern.search(name):
                pattern_name = self.non_binomial_patterns[i]
                self.stats['pattern_matches'][pattern_name] = \
                    self.stats['pattern_matches'].get(pattern_name, 0) + 1
                return False, f"Matches pattern: {pattern_name}"
        
        # Basic binomial check: should have exactly 2 words (Genus species)
        words = name.split()
        if len(words) != 2:
            return False, f"Not binomial ({len(words)} words)"
        
        genus, species = words
        
        # Genus should start with capital letter
        if not genus[0].isupper():
            return False, "Genus not capitalized"
        
        # Species should start with lowercase letter  
        if not species[0].islower():
            return False, "Species not lowercase"
        
        # Both should be mostly alphabetic (allow for some special chars like hyphens)
        if not re.match(r'^[A-Za-z][A-Za-z-]*$', genus):
            return False, "Genus contains invalid characters"
        
        if not re.match(r'^[a-z][a-z-]*$', species):
            return False, "Species contains invalid characters"
        
        # Check if it's a tautonym (genus == species) - THESE ARE VALID!
        if genus.lower() == species.lower():
            self.stats['tautonyms_found'] += 1
            if len(self.stats['tautonym_examples']) < 50:
                self.stats['tautonym_examples'].append(name)
            return True, "Proper binomial (tautonym)"
        
        return True, "Proper binomial"
    
    def parse_sharers_text(self, sharers_text):
        """
        Parse the sharers text accounting for | (BIN separators) and , (species separators)
        Returns list of all individual species names
        """
        if not sharers_text or pd.isna(sharers_text):
            return []
        
        # First split by | to get BIN groups
        bin_groups = sharers_text.split('|')
        
        all_species = []
        for bin_group in bin_groups:
            bin_group = bin_group.strip()
            if bin_group:  # Skip empty BIN groups
                # Split by comma to get individual species in this BIN
                species_in_bin = [species.strip() for species in bin_group.split(',') if species.strip()]
                all_species.extend(species_in_bin)
            else:
                self.stats['empty_bin_groups'] += 1
        
        return all_species
    
    def analyze_sharers_cell(self, sharers_text):
        """
        Analyze a single sharers cell that may contain multiple BIN groups and species names
        Returns dict with analysis results
        """
        if not sharers_text or pd.isna(sharers_text):
            return {
                'has_sharers': False,
                'species_count': 0,
                'proper_count': 0,
                'non_binomial_count': 0,
                'has_any_proper': False,
                'all_proper': True,  # vacuously true for empty
                'species_list': [],
                'proper_list': [],
                'non_binomial_list': []
            }
        
        # Parse the sharers text correctly
        species_names = self.parse_sharers_text(sharers_text)
        
        if not species_names:  # No species found after parsing
            return {
                'has_sharers': False,
                'species_count': 0,
                'proper_count': 0,
                'non_binomial_count': 0,
                'has_any_proper': False,
                'all_proper': True,
                'species_list': [],
                'proper_list': [],
                'non_binomial_list': []
            }
        
        proper_names = []
        non_binomial_names = []
        
        for name in species_names:
            is_proper, reason = self.is_proper_species_name(name)
            if is_proper:
                proper_names.append(name)
                # Store examples for logging
                if len(self.stats['proper_examples']) < 100:
                    self.stats['proper_examples'].append(name)
            else:
                non_binomial_names.append((name, reason))
                # Store examples for logging
                if len(self.stats['non_binomial_examples']) < 100:
                    self.stats['non_binomial_examples'].append((name, reason))
        
        return {
            'has_sharers': True,
            'species_count': len(species_names),
            'proper_count': len(proper_names),
            'non_binomial_count': len(non_binomial_names),
            'has_any_proper': len(proper_names) > 0,
            'all_proper': len(non_binomial_names) == 0 and len(proper_names) > 0,
            'species_list': species_names,
            'proper_list': proper_names,
            'non_binomial_list': non_binomial_names
        }
    
    def analyze_dataset(self, filepath):
        """
        Analyze the entire dataset
        """
        logger.info(f"Loading dataset from {filepath}")
        
        try:
            df = pd.read_csv(filepath, sep='\t', dtype=str)
            logger.info(f"Loaded {len(df)} rows with columns: {list(df.columns)}")
        except Exception as e:
            logger.error(f"Error loading file: {e}")
            return None
        
        # Check for required columns
        if 'sharers' not in df.columns:
            logger.error("Column 'sharers' not found in dataset")
            return None
        
        self.stats['total_entries'] = len(df)
        
        # Analyze each row
        results = []
        for idx, row in df.iterrows():
            if idx % 5000 == 0:
                logger.info(f"Processing row {idx:,}/{len(df):,}")
                
            analysis = self.analyze_sharers_cell(row.get('sharers', ''))
            
            # Update statistics
            if analysis['has_sharers']:
                self.stats['entries_with_sharers'] += 1
                self.stats['total_species_names'] += analysis['species_count']
                self.stats['proper_species_names'] += analysis['proper_count']
                self.stats['non_binomial_names'] += analysis['non_binomial_count']
            
            # Add proper species column value
            if not analysis['has_sharers']:
                proper_species_value = 'no_data'
            elif analysis['all_proper']:
                proper_species_value = 'yes'
            elif analysis['has_any_proper']:
                proper_species_value = 'mixed'
            else:
                proper_species_value = 'no'
            
            results.append({
                'taxonid': row.get('taxonid', ''),
                'BAGS': row.get('BAGS', ''),
                'BIN': row.get('BIN', ''),
                'sharers': row.get('sharers', ''),
                'proper_species': proper_species_value,
                'species_count': analysis['species_count'],
                'proper_count': analysis['proper_count'],
                'non_binomial_count': analysis['non_binomial_count']
            })
        
        result_df = pd.DataFrame(results)
        
        # Log statistics
        self.log_statistics()
        
        return result_df
    
    def log_statistics(self):
        """Log analysis statistics"""
        logger.info("=" * 70)
        logger.info("FINAL ANALYSIS STATISTICS (with tautonyms fixed)")
        logger.info("=" * 70)
        logger.info(f"Total entries: {self.stats['total_entries']:,}")
        logger.info(f"Entries with sharers: {self.stats['entries_with_sharers']:,}")
        logger.info(f"Empty BIN groups encountered: {self.stats['empty_bin_groups']:,}")
        logger.info(f"Total species names: {self.stats['total_species_names']:,}")
        logger.info(f"Proper species names: {self.stats['proper_species_names']:,}")
        logger.info(f"  - Including tautonyms: {self.stats['tautonyms_found']:,}")
        logger.info(f"Non-binomial names: {self.stats['non_binomial_names']:,}")
        
        if self.stats['total_species_names'] > 0:
            proper_pct = (self.stats['proper_species_names'] / self.stats['total_species_names']) * 100
            logger.info(f"Percentage proper: {proper_pct:.1f}%")
        
        logger.info("\n" + "=" * 70)
        logger.info("TOP NON-BINOMIAL PATTERNS")
        logger.info("=" * 70)
        sorted_patterns = sorted(self.stats['pattern_matches'].items(), 
                               key=lambda x: x[1], reverse=True)
        for pattern, count in sorted_patterns[:15]:
            logger.info(f"{pattern}: {count:,} matches")
        
        logger.info("\n" + "=" * 70)
        logger.info("SAMPLE PROPER SPECIES NAMES")
        logger.info("=" * 70)
        for i, example in enumerate(self.stats['proper_examples'][:30]):
            logger.info(f"{i+1:2d}. {example}")
        
        logger.info("\n" + "=" * 70)
        logger.info("SAMPLE TAUTONYMS (valid species names)")
        logger.info("=" * 70)
        for i, example in enumerate(self.stats['tautonym_examples'][:20]):
            logger.info(f"{i+1:2d}. {example}")
        
        logger.info("\n" + "=" * 70)
        logger.info("SAMPLE NON-BINOMIAL NAMES")
        logger.info("=" * 70)
        for i, (example, reason) in enumerate(self.stats['non_binomial_examples'][:30]):
            logger.info(f"{i+1:2d}. '{example}' -> {reason}")

def main():
    """Main execution function"""
    # File paths
    input_file = r"C:\_claude_files\projects\bold-ranker-testing\BAGS\assessed_BAGS.tsv"
    output_file = r"C:\_claude_files\projects\bold-ranker-testing\BAGS\assessed_BAGS_with_proper_species_FINAL.tsv"
    
    # Initialize analyzer
    analyzer = SpeciesNameAnalyzer()
    
    # Analyze dataset
    result_df = analyzer.analyze_dataset(input_file)
    
    if result_df is not None:
        # Save results
        logger.info(f"Saving results to {output_file}")
        result_df.to_csv(output_file, sep='\t', index=False)
        logger.info("Analysis complete!")
        
        # Show sample results
        logger.info("\n" + "=" * 60)
        logger.info("SAMPLE RESULTS")
        logger.info("=" * 60)
        print(result_df[['taxonid', 'sharers', 'proper_species', 'species_count', 'proper_count']].head(15))
        
        # Summary by proper_species category
        logger.info("\n" + "=" * 60)
        logger.info("SUMMARY BY CATEGORY")
        logger.info("=" * 60)
        summary = result_df['proper_species'].value_counts()
        print(summary)
        
        # Show some examples of different categories
        logger.info("\n" + "=" * 60)
        logger.info("EXAMPLES BY CATEGORY")
        logger.info("=" * 60)
        
        for category in ['yes', 'mixed', 'no']:
            if category in result_df['proper_species'].values:
                logger.info(f"\nCategory '{category}' examples:")
                examples = result_df[result_df['proper_species'] == category]['sharers'].dropna().head(5)
                for i, example in enumerate(examples):
                    logger.info(f"  {i+1}. {example}")
    
    logger.info(f"\nLog file saved as: {log_filename}")

if __name__ == "__main__":
    main()
