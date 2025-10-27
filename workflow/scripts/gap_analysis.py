#!/usr/bin/env python3
"""
Gap analysis for BOLD library curation workflow.

This script performs a gap analysis by comparing:
1. Input species list (from config FILTER_TAXA_LIST)
2. Result output data (result_output.tsv)
3. BAGS assessment data (assessed_BAGS.tsv)

The script identifies:
- Species from the input list and their coverage
- Additional species found in results but not in input list (missed species)
- BAGS grades and BIN information for all species
- Record counts per taxonid

Output: A comprehensive TSV file with gap analysis results.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict
import yaml

# Increase CSV field size limit to handle large TSV files
# Use a large but safe value for Windows compatibility
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    # Windows workaround
    csv.field_size_limit(2147483647)


def setup_logging(log_level: str = "INFO") -> None:
    """Set up logging configuration."""
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')
    
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def load_config(config_path: Path) -> Dict:
    """Load configuration from YAML file."""
    logging.info(f"Loading configuration from {config_path}")
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        logging.error(f"Failed to load configuration: {e}")
        sys.exit(1)


def format_species_name(species_lower: str) -> str:
    """
    Format species name to proper binomial nomenclature: Genus species
    (capitalize only the first word/genus, keep species lowercase)
    
    Examples:
        'anax parthenope' -> 'Anax parthenope'
        'acisoma inflatum' -> 'Acisoma inflatum'
    """
    parts = species_lower.split()
    if len(parts) >= 2:
        # Capitalize first word (genus), keep rest lowercase
        return parts[0].capitalize() + ' ' + ' '.join(parts[1:])
    elif len(parts) == 1:
        # Just genus name
        return parts[0].capitalize()
    return species_lower


def parse_species_list(species_file: Path) -> Tuple[Dict[str, List[str]], Set[str], Dict[str, str]]:
    """
    Parse species list with optional synonyms.
    
    Format: valid_species OR valid_species;synonym1;synonym2;etc
    
    Returns:
        Tuple of:
        - Dictionary mapping valid species (lowercase) to list of synonyms
        - Set of all valid species names (original case)
        - Dictionary mapping synonym (lowercase) to valid species (lowercase)
    """
    logging.info(f"Loading species list from {species_file}")
    species_synonyms = {}  # lowercase valid species -> list of synonyms
    valid_species_set = set()  # original case valid species
    synonym_to_valid = {}  # lowercase synonym -> lowercase valid species
    
    try:
        with open(species_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split(';')
                valid_species = parts[0].strip()
                
                if not valid_species:
                    logging.warning(f"Line {line_num}: Empty species name, skipping")
                    continue
                
                # Store with lowercase key for case-insensitive matching
                key = valid_species.lower()
                
                # Get synonyms (if any)
                synonyms = [s.strip() for s in parts[1:] if s.strip()]
                
                if key in species_synonyms:
                    logging.warning(f"Duplicate species found: {valid_species}")
                
                species_synonyms[key] = synonyms
                valid_species_set.add(valid_species)
                
                # Build synonym to valid species mapping
                for syn in synonyms:
                    synonym_to_valid[syn.lower()] = key
        
        logging.info(f"Loaded {len(species_synonyms)} species from input list")
        total_synonyms = sum(len(syns) for syns in species_synonyms.values())
        logging.info(f"Total synonyms: {total_synonyms}")
        
        return species_synonyms, valid_species_set, synonym_to_valid
        
    except Exception as e:
        logging.error(f"Failed to load species list: {e}")
        sys.exit(1)


def load_assessed_bags(bags_file: Path) -> Dict[str, Dict]:
    """
    Load BAGS assessment data.
    
    Handles both UTF-8 and Latin-1 encoded files automatically.
    
    Returns dictionary: taxonid -> {BAGS: grade, BIN: uri, sharers: list}
    """
    logging.info(f"Loading BAGS assessments from {bags_file}")
    bags_data = {}
    
    # Try UTF-8 first, fall back to Latin-1 if decoding fails
    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(bags_file, 'r', encoding=encoding) as f:
                reader = csv.DictReader(f, delimiter='\t')
                
                for row in reader:
                    taxonid = row.get('taxonid', '').strip()
                    if not taxonid:
                        continue
                    
                    bags_data[taxonid] = {
                        'BAGS': row.get('BAGS', '').strip(),
                        'BIN': row.get('BIN', '').strip(),
                        'sharers': row.get('sharers', '').strip()
                    }
            
            # If we successfully read the file, log and return
            if encoding == 'latin-1':
                logging.info(f"Successfully read file using {encoding} encoding")
            logging.info(f"Loaded BAGS data for {len(bags_data)} taxonids")
            return bags_data
            
        except UnicodeDecodeError as e:
            if encoding == 'utf-8':
                logging.warning(f"UTF-8 decoding failed for BAGS file, trying Latin-1 encoding")
                bags_data = {}  # Reset for retry
                continue
            else:
                logging.error(f"Failed to load BAGS data with Latin-1 encoding: {e}")
                sys.exit(1)
        except Exception as e:
            logging.error(f"Failed to load BAGS data: {e}")
            sys.exit(1)
    
    # Should never reach here
    logging.error("Failed to load BAGS data with any supported encoding")
    sys.exit(1)


def load_result_output(result_file: Path) -> Tuple[Dict, Dict]:
    """
    Load result_output.tsv and extract species information.
    
    Includes subspecies in the full species name (e.g., "Genus species subspecies")
    
    Handles both UTF-8 and Latin-1 encoded files automatically.
    
    Returns:
        Tuple of:
        - species_taxonid_map: dict mapping species (lowercase) -> list of (taxonid, taxonomy_dict)
        - taxonid_record_count: dict mapping taxonid -> count of records
    """
    logging.info(f"Loading result output from {result_file}")
    
    species_taxonid_map = defaultdict(list)
    taxonid_record_count = defaultdict(int)
    subspecies_count = 0
    
    # Try UTF-8 first, fall back to Latin-1 if decoding fails
    for encoding in ['utf-8', 'latin-1']:
        try:
            with open(result_file, 'r', encoding=encoding) as f:
                reader = csv.DictReader(f, delimiter='\t')
                
                for row in reader:
                    species_full = row.get('species', '').strip()
                    subspecies = row.get('subspecies', '').strip()
                    taxonid = row.get('taxonid', '').strip()
                    
                    if not species_full or not taxonid:
                        continue
                    
                    # If subspecies exists, include it in the full species name
                    if subspecies and subspecies.lower() not in ['none', 'null', '']:
                        # Subspecies field may contain full trinomial (e.g., "Genus species subspecies")
                        # Extract just the subspecies epithet (last word)
                        subspecies_parts = subspecies.split()
                        if len(subspecies_parts) >= 3:
                            # Full trinomial provided, take the last word
                            subspecies_epithet = subspecies_parts[-1]
                        else:
                            # Just subspecies epithet provided
                            subspecies_epithet = subspecies
                        
                        # Create full trinomial: "Genus species subspecies"
                        species = f"{species_full} {subspecies_epithet}"
                        subspecies_count += 1
                    else:
                        species = species_full
                    
                    # Count records per taxonid
                    taxonid_record_count[taxonid] += 1
                    
                    # Extract taxonomy (only add once per species-taxonid combo)
                    species_lower = species.lower()
                    
                    # Check if we already have this species-taxonid combination
                    existing_taxonids = [t for t, _ in species_taxonid_map[species_lower]]
                    if taxonid not in existing_taxonids:
                        taxonomy = {
                            'kingdom': row.get('kingdom', '').strip(),
                            'phylum': row.get('phylum', '').strip(),
                            'class': row.get('class', '').strip(),
                            'order': row.get('order', '').strip(),
                            'family': row.get('family', '').strip(),
                            'genus': row.get('genus', '').strip()
                        }
                        species_taxonid_map[species_lower].append((taxonid, taxonomy))
            
            # If we successfully read the file, log success and return
            if encoding == 'latin-1':
                logging.info(f"Successfully read file using {encoding} encoding")
            logging.info(f"Loaded {len(species_taxonid_map)} unique species from result output")
            logging.info(f"Found {subspecies_count} subspecies records")
            logging.info(f"Processed {sum(taxonid_record_count.values())} total records")
            logging.info(f"Found {len(taxonid_record_count)} unique taxonids")
            
            return dict(species_taxonid_map), dict(taxonid_record_count)
            
        except UnicodeDecodeError as e:
            if encoding == 'utf-8':
                logging.warning(f"UTF-8 decoding failed at position {e.start}, trying Latin-1 encoding")
                # Reset counters for retry
                species_taxonid_map = defaultdict(list)
                taxonid_record_count = defaultdict(int)
                subspecies_count = 0
                continue
            else:
                # Latin-1 should handle any byte sequence, so this shouldn't happen
                logging.error(f"Failed to load result output with Latin-1 encoding: {e}")
                sys.exit(1)
        except Exception as e:
            logging.error(f"Failed to load result output: {e}")
            sys.exit(1)
    
    # Should never reach here
    logging.error("Failed to load result output with any supported encoding")
    sys.exit(1)


def perform_gap_analysis(
    species_synonyms: Dict[str, List[str]],
    valid_species_set: Set[str],
    synonym_to_valid: Dict[str, str],
    species_taxonid_map: Dict,
    taxonid_record_count: Dict,
    bags_data: Dict
) -> List[Dict]:
    """
    Perform gap analysis by merging all data sources.
    
    Returns list of dictionaries, one per unique species-taxonid combination.
    """
    logging.info("Performing gap analysis...")
    
    results = []
    processed_species = set()
    
    # Process species from input list
    for species_lower, synonyms in species_synonyms.items():
        processed_species.add(species_lower)
        
        # Find matching species in result output
        if species_lower in species_taxonid_map:
            # Species found in results
            for taxonid, taxonomy in species_taxonid_map[species_lower]:
                # Get BAGS data for this taxonid
                bags_info = bags_data.get(taxonid, {})
                
                result = {
                    'species': format_species_name(species_lower),  # Proper format: Genus species (or Genus species subspecies)
                    'synonyms': '|'.join(synonyms) if synonyms else '',
                    'gaplist_species': 'Valid',  # Matched valid species in input gap list
                    'BAGS_grade': bags_info.get('BAGS', ''),
                    'BIN_uri': bags_info.get('BIN', ''),
                    'sharers': bags_info.get('sharers', ''),
                    'total_record_count': taxonid_record_count.get(taxonid, 0),
                    'kingdom': taxonomy.get('kingdom', ''),
                    'phylum': taxonomy.get('phylum', ''),
                    'class': taxonomy.get('class', ''),
                    'order': taxonomy.get('order', ''),
                    'family': taxonomy.get('family', ''),
                    'genus': taxonomy.get('genus', '')
                }
                results.append(result)
        else:
            # Species in input list but not in results
            result = {
                'species': format_species_name(species_lower),  # Proper format: Genus species
                'synonyms': '|'.join(synonyms) if synonyms else '',
                'gaplist_species': 'Valid',  # Matched valid species in input gap list
                'BAGS_grade': '',
                'BIN_uri': '',
                'sharers': '',
                'total_record_count': 0,
                'kingdom': '',
                'phylum': '',
                'class': '',
                'order': '',
                'family': '',
                'genus': ''
            }
            results.append(result)
    
    # Process species found in results but NOT in input list (check for synonyms or mark as No)
    for species_lower, taxonid_list in species_taxonid_map.items():
        if species_lower not in processed_species:
            # Check if this species is a synonym
            if species_lower in synonym_to_valid:
                # This species matches a synonym
                gaplist_value = 'Synonym'
            else:
                # This species is not in the input gap list at all
                gaplist_value = 'No'
            
            for taxonid, taxonomy in taxonid_list:
                bags_info = bags_data.get(taxonid, {})
                
                result = {
                    'species': format_species_name(species_lower),  # Proper format: Genus species (or Genus species subspecies)
                    'synonyms': '',
                    'gaplist_species': gaplist_value,  # Either 'Synonym' or 'No'
                    'BAGS_grade': bags_info.get('BAGS', ''),
                    'BIN_uri': bags_info.get('BIN', ''),
                    'sharers': bags_info.get('sharers', ''),
                    'total_record_count': taxonid_record_count.get(taxonid, 0),
                    'kingdom': taxonomy.get('kingdom', ''),
                    'phylum': taxonomy.get('phylum', ''),
                    'class': taxonomy.get('class', ''),
                    'order': taxonomy.get('order', ''),
                    'family': taxonomy.get('family', ''),
                    'genus': taxonomy.get('genus', '')
                }
                results.append(result)
    
    logging.info(f"Gap analysis complete: {len(results)} total entries")
    
    # Summary statistics
    valid_species_count = len([r for r in results if r['gaplist_species'] == 'Valid'])
    synonym_species_count = len([r for r in results if r['gaplist_species'] == 'Synonym'])
    non_gaplist_species_count = len([r for r in results if r['gaplist_species'] == 'No'])
    species_with_bags = len([r for r in results if r['BAGS_grade']])
    
    logging.info(f"  - Valid species (from gap list): {valid_species_count}")
    logging.info(f"  - Synonym species: {synonym_species_count}")
    logging.info(f"  - Non-gap list species: {non_gaplist_species_count}")
    logging.info(f"  - Species with BAGS assessment: {species_with_bags}")
    
    return results


def write_gap_analysis(results: List[Dict], output_file: Path) -> None:
    """Write gap analysis results to TSV file."""
    logging.info(f"Writing gap analysis to {output_file}")
    
    fieldnames = [
        'species',
        'synonyms',
        'gaplist_species',
        'BAGS_grade',
        'BIN_uri',
        'sharers',
        'total_record_count',
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus'
    ]
    
    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
        
        logging.info(f"Successfully wrote {len(results)} entries to {output_file}")
        
    except Exception as e:
        logging.error(f"Failed to write output file: {e}")
        sys.exit(1)


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Perform gap analysis for BOLD library curation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using config file (reads FILTER_TAXA_LIST from config)
  python gap_analysis.py --config config/config.yml \\
      --result-output results/result_output.tsv \\
      --assessed-bags results/assessed_BAGS.tsv \\
      --output results/gap_analysis.tsv
  
  # Override species list via CLI
  python gap_analysis.py --config config/config.yml \\
      --species-list custom_species.csv \\
      --result-output results/result_output.tsv \\
      --assessed-bags results/assessed_BAGS.tsv \\
      --output results/gap_analysis.tsv
        """
    )
    
    parser.add_argument(
        '--config',
        type=Path,
        help='Path to config.yml file (to read FILTER_TAXA_LIST)'
    )
    
    parser.add_argument(
        '--species-list',
        type=Path,
        help='Path to species list CSV (overrides config FILTER_TAXA_LIST)'
    )
    
    parser.add_argument(
        '--result-output',
        type=Path,
        required=True,
        help='Path to result_output.tsv file'
    )
    
    parser.add_argument(
        '--assessed-bags',
        type=Path,
        required=True,
        help='Path to assessed_BAGS.tsv file'
    )
    
    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Path to output gap analysis TSV file'
    )
    
    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level)
    
    logging.info("=" * 80)
    logging.info("BOLD Library Curation - Gap Analysis")
    logging.info("=" * 80)
    
    # Determine species list source
    species_list_path = None
    
    if args.species_list:
        species_list_path = args.species_list
        logging.info(f"Using species list from CLI argument: {species_list_path}")
    elif args.config:
        config = load_config(args.config)
        filter_taxa_list = config.get('FILTER_TAXA_LIST')
        if filter_taxa_list:
            species_list_path = Path(filter_taxa_list)
            logging.info(f"Using species list from config: {species_list_path}")
        else:
            logging.error("FILTER_TAXA_LIST not found in config file")
            sys.exit(1)
    else:
        logging.error("Must provide either --config or --species-list")
        parser.print_help()
        sys.exit(1)
    
    # Validate input files exist
    if not species_list_path.exists():
        logging.error(f"Species list file not found: {species_list_path}")
        sys.exit(1)
    
    if not args.result_output.exists():
        logging.error(f"Result output file not found: {args.result_output}")
        sys.exit(1)
    
    if not args.assessed_bags.exists():
        logging.error(f"Assessed BAGS file not found: {args.assessed_bags}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    # Load all data
    species_synonyms, valid_species_set, synonym_to_valid = parse_species_list(species_list_path)
    bags_data = load_assessed_bags(args.assessed_bags)
    species_taxonid_map, taxonid_record_count = load_result_output(args.result_output)
    
    # Perform gap analysis
    results = perform_gap_analysis(
        species_synonyms,
        valid_species_set,
        synonym_to_valid,
        species_taxonid_map,
        taxonid_record_count,
        bags_data
    )
    
    # Write output
    write_gap_analysis(results, args.output)
    
    logging.info("=" * 80)
    logging.info("Gap analysis complete!")
    logging.info("=" * 80)


if __name__ == '__main__':
    main()
