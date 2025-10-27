#!/usr/bin/env python3
"""
Pre-scoring filter for BOLD library curation workflow.

This script filters BOLD TSV data before database creation to reduce dataset size
and improve performance for small target lists.

The filtering process:
1. Finds records matching taxa, country, and/or marker criteria
2. Optionally expands the dataset by including all records sharing BIN_URIs with the matched records
3. Always preserves records that match the initial criteria, even if they lack BIN_URIs

Records without BIN_URIs that match the taxa/country/marker criteria are always kept - 
the purpose of BIN expansion is to add additional related records, not to exclude existing ones.

When multiple filters are used (taxa, countries, marker), the marker filter acts as an additional 
constraint, keeping only records that match ALL specified criteria.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Set, Optional, Dict, Any
import time

# Increase CSV field size limit to handle large TSV files
# Use platform-specific max size (Windows has different limits)
import platform
if platform.system() == 'Windows':
    csv.field_size_limit(2147483647)  # Max int on Windows
else:
    csv.field_size_limit(sys.maxsize)


def sanitize_field(field: str) -> str:
    """
    Sanitize a field value to prevent parsing issues caused by unescaped quotes.
    
    Args:
        field: Raw field value from TSV
        
    Returns:
        Sanitized field value with quotes removed and whitespace cleaned
    """
    if not field:
        return ''
    
    # Remove embedded newlines and carriage returns that cause row merging
    field = field.replace('\r', ' ').replace('\n', ' ')
    
    # Remove all double quotes to prevent CSV parsing issues
    # The BOLD data contains unescaped quotes like "Syn. that break parsers
    field = field.replace('"', '')
    
    return field.strip()


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


def has_valid_species(species_value: str) -> bool:
    """
    Check if species value is valid (not None, empty, or "None").
    
    Args:
        species_value: Raw species value from TSV
        
    Returns:
        True if species has a valid value, False otherwise
    """
    if not species_value:
        return False
    
    species_clean = species_value.strip().lower()
    return species_clean and species_clean != 'none'


def load_taxa_list(file_path: str) -> Set[str]:
    """
    Load species names from taxa list file (semicolon separated format).
    
    Args:
        file_path: Path to taxa list file
        
    Returns:
        Set of species names (lowercased for case-insensitive matching)
    """
    taxa_set = set()
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                    
                # Split by semicolon and add all variants
                species_variants = [s.strip() for s in line.split(';')]
                for species in species_variants:
                    if species:  # Skip empty strings
                        taxa_set.add(species.lower())
        
        logging.info(f"Loaded {len(taxa_set)} unique species names from {file_path}")
        return taxa_set
        
    except FileNotFoundError:
        logging.error(f"Taxa list file not found: {file_path}")
        raise
    except Exception as e:
        logging.error(f"Error loading taxa list from {file_path}: {e}")
        raise


def load_country_list(file_path: str) -> Set[str]:
    """
    Load country codes from country list file (semicolon separated format).
    
    Args:
        file_path: Path to country list file (semicolon-separated format like taxa)
        
    Returns:
        Set of country codes (uppercased)
    """
    country_set = set()
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                    
                # Split by semicolon and add all variants (ISO code + country names)
                country_variants = [c.strip() for c in line.split(';')]
                for country in country_variants:
                    if country:  # Skip empty strings
                        country_set.add(country.upper())
        
        logging.info(f"Loaded {len(country_set)} unique country codes/names from {file_path}")
        return country_set
        
    except FileNotFoundError:
        logging.error(f"Country list file not found: {file_path}")
        raise
    except Exception as e:
        logging.error(f"Error loading country list from {file_path}: {e}")
        raise


def filter_by_taxa(input_file: str, taxa_set: Set[str], species_column: str) -> Set[str]:
    """
    Return set of processids matching taxa list.
    
    Args:
        input_file: Path to input TSV file
        taxa_set: Set of species names to match (lowercased)
        species_column: Name of species column
        
    Returns:
        Set of processids that match the taxa criteria
    """
    matching_processids = set()
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        
        for row_num, row in enumerate(reader, 1):
            if row_num % 10000 == 0:
                logging.info(f"Taxa filtering: processed {row_num} rows")
            
            # Sanitize all fields to prevent parsing issues
            row = {k: sanitize_field(v) if v else '' for k, v in row.items()}
                
            species_raw = row.get(species_column, '') or ''
            species = species_raw.strip().lower() if species_raw else ''
            if species and species in taxa_set:
                processid_raw = row.get('processid', '') or ''
                processid = processid_raw.strip() if processid_raw else ''
                if processid:
                    matching_processids.add(processid)
    
    logging.info(f"Taxa filter matched {len(matching_processids)} processids")
    return matching_processids


def filter_by_countries(input_file: str, country_set: Set[str], country_column: str) -> Set[str]:
    """
    Return set of processids matching country list.
    
    Args:
        input_file: Path to input TSV file
        country_set: Set of country codes to match (uppercased)
        country_column: Name of country column
        
    Returns:
        Set of processids that match the country criteria
    """
    matching_processids = set()
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        
        for row_num, row in enumerate(reader, 1):
            if row_num % 10000 == 0:
                logging.info(f"Country filtering: processed {row_num} rows")
            
            # Sanitize all fields to prevent parsing issues
            row = {k: sanitize_field(v) if v else '' for k, v in row.items()}
                
            country_raw = row.get(country_column, '') or ''
            country = country_raw.strip().upper() if country_raw else ''
            if country and country in country_set:
                processid_raw = row.get('processid', '') or ''
                processid = processid_raw.strip() if processid_raw else ''
                if processid:
                    matching_processids.add(processid)
    
    logging.info(f"Country filter matched {len(matching_processids)} processids")
    return matching_processids


def filter_by_marker(input_file: str, marker_code: str, marker_column: str) -> Set[str]:
    """
    Return set of processids matching the specified marker code.
    
    Args:
        input_file: Path to input TSV file
        marker_code: Marker code to match (e.g., "COI-5P")
        marker_column: Name of marker column
        
    Returns:
        Set of processids that match the marker criteria
    """
    matching_processids = set()
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        
        for row_num, row in enumerate(reader, 1):
            if row_num % 10000 == 0:
                logging.info(f"Marker filtering: processed {row_num} rows")
            
            # Sanitize all fields to prevent parsing issues
            row = {k: sanitize_field(v) if v else '' for k, v in row.items()}
                
            marker_raw = row.get(marker_column, '') or ''
            marker = marker_raw.strip() if marker_raw else ''
            if marker and marker == marker_code:
                processid_raw = row.get('processid', '') or ''
                processid = processid_raw.strip() if processid_raw else ''
                if processid:
                    matching_processids.add(processid)
    
    logging.info(f"Marker filter matched {len(matching_processids)} processids")
    return matching_processids


def get_bins_for_processids(processids: Set[str], input_file: str) -> Set[str]:
    """
    Extract BIN_URIs for given processids.
    
    Args:
        processids: Set of processids to look up
        input_file: Path to input TSV file
        
    Returns:
        Set of BIN_URIs associated with the processids
    """
    bins = set()
    processids_with_bins = 0
    processids_without_bins = 0
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        
        for row in reader:
            # Sanitize all fields to prevent parsing issues
            row = {k: sanitize_field(v) if v else '' for k, v in row.items()}
            
            processid_raw = row.get('processid', '') or ''
            processid = processid_raw.strip() if processid_raw else ''
            if processid in processids:
                bin_uri_raw = row.get('bin_uri', '') or ''
                bin_uri = bin_uri_raw.strip() if bin_uri_raw else ''
                if bin_uri and bin_uri.lower() != 'none':
                    bins.add(bin_uri)
                    processids_with_bins += 1
                else:
                    processids_without_bins += 1
    
    logging.info(f"Found {len(bins)} unique BIN_URIs for {len(processids)} processids")
    logging.info(f"  - {processids_with_bins} processids have BIN_URIs")
    logging.info(f"  - {processids_without_bins} processids have no BIN_URI (will be preserved)")
    return bins


def get_processids_for_bins(bins: Set[str], input_file: str) -> Set[str]:
    """
    Find all processids sharing the given BIN_URIs.
    
    Args:
        bins: Set of BIN_URIs to match
        input_file: Path to input TSV file
        
    Returns:
        Set of processids that share the BIN_URIs
    """
    processids = set()
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        
        for row in reader:
            # Sanitize all fields to prevent parsing issues
            row = {k: sanitize_field(v) if v else '' for k, v in row.items()}
            
            bin_uri_raw = row.get('bin_uri', '') or ''
            bin_uri = bin_uri_raw.strip() if bin_uri_raw else ''
            if bin_uri and bin_uri in bins:
                processid_raw = row.get('processid', '') or ''
                processid = processid_raw.strip() if processid_raw else ''
                if processid:
                    processids.add(processid)
    
    logging.info(f"Found {len(processids)} processids sharing {len(bins)} BIN_URIs")
    return processids


def expand_by_bin_sharing(initial_processids: Set[str], input_file: str) -> Set[str]:
    """
    Iteratively expand processid set based on BIN_URI sharing.
    
    Args:
        initial_processids: Starting set of processids
        input_file: Path to input TSV file
        
    Returns:
        Expanded set of processids after BIN sharing
    """
    current_processids = initial_processids.copy()
    previous_size = 0
    iteration = 0
    
    logging.info(f"Starting BIN expansion with {len(current_processids)} processids")
    
    while len(current_processids) > previous_size:
        iteration += 1
        previous_size = len(current_processids)
        
        # Get all BIN_URIs for current processids
        current_bins = get_bins_for_processids(current_processids, input_file)
        
        if not current_bins:
            logging.info("No BIN_URIs found, stopping expansion")
            break
            
        # Get all processids sharing those BIN_URIs
        current_processids = get_processids_for_bins(current_bins, input_file)
        
        logging.info(f"BIN expansion iteration {iteration}: {len(current_processids)} processids")
        
        # Safety check to prevent infinite loops
        if iteration > 50:
            logging.warning("BIN expansion reached maximum iterations (50), stopping")
            break
    
    expansion_factor = len(current_processids) / len(initial_processids) if initial_processids else 0
    logging.info(f"BIN expansion completed: {len(initial_processids)} -> {len(current_processids)} "
                f"processids (factor: {expansion_factor:.2f})")
    
    return current_processids


def filter_by_kingdom(input_file: str, kingdom_set: Set[str], kingdom_column: str) -> Set[str]:
    """
    Return set of processids matching the specified kingdoms.
    
    Args:
        input_file: Path to input TSV file
        kingdom_set: Set of kingdom names to match (case-insensitive)
        kingdom_column: Name of kingdom column
        
    Returns:
        Set of processids that match the kingdom criteria
    """
    matching_processids = set()
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        
        for row_num, row in enumerate(reader, 1):
            if row_num % 10000 == 0:
                logging.info(f"Kingdom filtering: processed {row_num} rows")
            
            # Sanitize all fields to prevent parsing issues
            row = {k: sanitize_field(v) if v else '' for k, v in row.items()}
                
            kingdom_raw = row.get(kingdom_column, '') or ''
            kingdom = kingdom_raw.strip() if kingdom_raw else ''
            if kingdom and kingdom.lower() in kingdom_set:
                processid_raw = row.get('processid', '') or ''
                processid = processid_raw.strip() if processid_raw else ''
                if processid:
                    matching_processids.add(processid)
    
    logging.info(f"Kingdom filter matched {len(matching_processids)} processids")
    return matching_processids


def detect_column_names(input_file: str) -> Dict[str, str]:
    """
    Auto-detect column names in the TSV file.
    
    Args:
        input_file: Path to input TSV file
        
    Returns:
        Dictionary mapping standard names to actual column names
    """
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        headers = reader.fieldnames or []
    
    column_map = {}
    
    # Species column detection
    for col in headers:
        if col.lower() in ['species', 'species_name', 'taxon']:
            column_map['species'] = col
            break
    
    # Country column detection  
    for col in headers:
        if col.lower() in ['country_iso', 'country', 'country_code']:
            column_map['country'] = col
            break
    
    # BIN column detection
    for col in headers:
        if col.lower() in ['bin_uri', 'bin', 'bin_code']:
            column_map['bin'] = col
            break
    
    # Marker column detection
    for col in headers:
        if col.lower() in ['marker_code', 'marker', 'gene']:
            column_map['marker'] = col
            break
    
    # Kingdom column detection
    for col in headers:
        if col.lower() in ['kingdom']:
            column_map['kingdom'] = col
            break
    
    logging.info(f"Detected columns: {column_map}")
    return column_map


def prescoring_filter(
    input_tsv: str,
    output_tsv: str,
    taxa_list: Optional[str] = None,
    country_list: Optional[str] = None,
    marker_code: Optional[str] = None,
    kingdom_list: Optional[list] = None,
    enable_bin_sharing: bool = False,
    filter_species: bool = False
) -> Dict[str, Any]:
    """
    Main filtering function with multi-criteria support.
    
    Args:
        input_tsv: Path to input BOLD TSV file
        output_tsv: Path to output filtered TSV file
        taxa_list: Path to taxa list file (optional)
        country_list: Path to country list file (optional)
        marker_code: Marker code to filter by (e.g., "COI-5P") (optional)
        kingdom_list: List of kingdom names to filter by (optional)
        enable_bin_sharing: Whether to include BIN_URI sharing expansion
        
    Returns:
        Dictionary with filtering statistics
    """
    start_time = time.time()
    
    # Validate inputs
    if not Path(input_tsv).exists():
        raise FileNotFoundError(f"Input file not found: {input_tsv}")
    
    if not taxa_list and not country_list and not marker_code and not kingdom_list:
        raise ValueError("At least one filter (taxa_list, country_list, marker_code, or kingdom_list) must be provided")
    
    # Auto-detect column names
    column_map = detect_column_names(input_tsv)
    
    # Count total rows for progress reporting
    with open(input_tsv, 'r', encoding='utf-8') as f:
        total_rows = sum(1 for _ in f) - 1  # Subtract header
    logging.info(f"Input file has {total_rows} data rows")
    
    # Load filter lists
    taxa_set = set()
    country_set = set()
    
    if taxa_list:
        taxa_set = load_taxa_list(taxa_list)
    
    if country_list:
        country_set = load_country_list(country_list)
    
    # Apply filters
    matching_processids = set()
    
    if taxa_list and 'species' in column_map:
        taxa_matches = filter_by_taxa(input_tsv, taxa_set, column_map['species'])
        matching_processids.update(taxa_matches)
    elif taxa_list:
        logging.warning("Taxa list provided but species column not found")
    
    if country_list and 'country' in column_map:
        country_matches = filter_by_countries(input_tsv, country_set, column_map['country'])
        matching_processids.update(country_matches)
    elif country_list:
        logging.warning("Country list provided but country column not found")
    
    initial_matches = len(matching_processids)
    logging.info(f"Initial filtering matched {initial_matches} processids")
    
    # Store original matches (including those without bin_uri) to preserve them
    original_matches = matching_processids.copy()
    
    # BIN sharing expansion
    if enable_bin_sharing and matching_processids and 'bin' in column_map:
        expanded_matches = expand_by_bin_sharing(matching_processids, input_tsv)
        # Always include original matches even if they don't have bin_uri
        matching_processids = original_matches.union(expanded_matches)
        logging.info(f"Combined original matches with BIN expansion: {len(matching_processids)} total processids")
    elif enable_bin_sharing and not matching_processids:
        logging.warning("BIN sharing requested but no initial matches found")
    elif enable_bin_sharing:
        logging.warning("BIN sharing requested but bin_uri column not found")
    
    # Apply marker filter at the end - this removes records that don't match the specified marker
    if marker_code and 'marker' in column_map:
        marker_matches = filter_by_marker(input_tsv, marker_code, column_map['marker'])
        if matching_processids:
            # Intersect final results with marker filter
            before_marker_count = len(matching_processids)
            matching_processids = matching_processids.intersection(marker_matches)
            logging.info(f"Marker filter: {before_marker_count} -> {len(matching_processids)} processids")
        else:
            # If no other filters, use marker filter as the only filter
            matching_processids = marker_matches
            logging.info(f"Using marker filter only: {len(matching_processids)} processids")
    elif marker_code:
        logging.warning("Marker code provided but marker column not found")
    
    # Apply kingdom filter at the end - this removes records that don't match the specified kingdoms
    if kingdom_list and 'kingdom' in column_map:
        kingdom_set = {k.lower() for k in kingdom_list}  # Convert to lowercase set
        kingdom_matches = filter_by_kingdom(input_tsv, kingdom_set, column_map['kingdom'])
        if matching_processids:
            # Intersect final results with kingdom filter
            before_kingdom_count = len(matching_processids)
            matching_processids = matching_processids.intersection(kingdom_matches)
            logging.info(f"Kingdom filter: {before_kingdom_count} -> {len(matching_processids)} processids")
        else:
            # If no other filters, use kingdom filter as the only filter
            matching_processids = kingdom_matches
            logging.info(f"Using kingdom filter only: {len(matching_processids)} processids")
    elif kingdom_list:
        logging.warning("Kingdom list provided but kingdom column not found")
    
    final_matches = len(matching_processids)
    
    # Write filtered output
    rows_written = 0
    rows_skipped_no_processid = 0
    rows_skipped_not_matching = 0
    rows_skipped_wrong_marker = 0
    rows_skipped_no_species = 0
    unique_processids_written = set()
    
    with open(input_tsv, 'r', encoding='utf-8') as infile, \
         open(output_tsv, 'w', encoding='utf-8', newline='') as outfile:
        
        reader = csv.DictReader(infile, delimiter='\t', quoting=csv.QUOTE_NONE)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter='\t', 
                               quoting=csv.QUOTE_NONE)
        writer.writeheader()
        
        for row_num, row in enumerate(reader, 1):
            if row_num % 10000 == 0:
                logging.info(f"Output writing: processed {row_num}/{total_rows} rows")
            
            # Sanitize all fields FIRST to prevent parsing issues
            row = {k: sanitize_field(v) if v else '' for k, v in row.items()}
                
            processid_raw = row.get('processid', '') or ''
            processid = processid_raw.strip() if processid_raw else ''
            if not processid:
                rows_skipped_no_processid += 1
                continue
                
            # Check if processid matches our filtering criteria
            if processid not in matching_processids:
                rows_skipped_not_matching += 1
                continue
            
            # Additional marker check: if marker filtering is enabled, 
            # check this specific record's marker, not just the processid
            if marker_code and 'marker' in column_map:
                record_marker_raw = row.get(column_map['marker'], '') or ''
                record_marker = record_marker_raw.strip() if record_marker_raw else ''
                if record_marker != marker_code:
                    rows_skipped_wrong_marker += 1
                    continue
            
            # Species filter check (if enabled)
            if filter_species and 'species' in column_map:
                species_raw = row.get(column_map['species'], '') or ''
                if not has_valid_species(species_raw):
                    rows_skipped_no_species += 1
                    continue
            
            # If we get here, the record passes all filters
            # Filter out any fields not in fieldnames to prevent DictWriter errors
            valid_keys = set(writer.fieldnames)
            row = {k: v for k, v in row.items() if k in valid_keys}
            writer.writerow(row)
            rows_written += 1
            unique_processids_written.add(processid)
    
    logging.info(f"Output summary: {rows_written} rows written, {rows_skipped_no_processid} rows skipped (no processid), {rows_skipped_not_matching} rows skipped (processid not matching), {rows_skipped_wrong_marker} rows skipped (wrong marker), {rows_skipped_no_species} rows skipped (no valid species)")
    logging.info(f"Unique processids written: {len(unique_processids_written)}")
    
    # Note about the counts
    if marker_code:
        logging.info(f"Note: When marker filtering is enabled, final processid count ({final_matches}) represents processids that have at least one record with the specified marker, but final row count ({rows_written}) represents only the actual records with that marker.")
    
    # Calculate statistics
    end_time = time.time()
    processing_time = end_time - start_time
    reduction_factor = (total_rows - rows_written) / total_rows if total_rows > 0 else 0
    
    stats = {
        'input_rows': total_rows,
        'output_rows': rows_written,
        'reduction_factor': reduction_factor,
        'reduction_percentage': reduction_factor * 100,
        'initial_matches': initial_matches,
        'final_matches': final_matches,
        'bin_expansion_factor': final_matches / initial_matches if initial_matches > 0 else 0,
        'processing_time_seconds': processing_time,
        'taxa_filter_used': bool(taxa_list),
        'country_filter_used': bool(country_list),
        'marker_filter_used': bool(marker_code),
        'kingdom_filter_used': bool(kingdom_list),
        'bin_sharing_used': enable_bin_sharing
    }
    
    logging.info(f"Filtering completed in {processing_time:.2f} seconds")
    logging.info(f"Reduced dataset from {total_rows:,} to {rows_written:,} rows "
                f"({reduction_factor*100:.1f}% reduction)")
    
    return stats


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Pre-scoring filter for BOLD library curation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter by taxa list only
  python prescoring_filter.py --input data.tsv --output filtered.tsv --taxa-list species.csv
  
  # Filter by country list only  
  python prescoring_filter.py --input data.tsv --output filtered.tsv --country-list countries.txt
  
  # Filter by marker only
  python prescoring_filter.py --input data.tsv --output filtered.tsv --marker COI-5P
  
  # Combined filtering with BIN sharing
  python prescoring_filter.py --input data.tsv --output filtered.tsv \\
    --taxa-list species.csv --country-list countries.txt --marker COI-5P --enable-bin-sharing
        """
    )
    
    parser.add_argument('--input', required=True, 
                        help='Input BOLD TSV file')
    parser.add_argument('--output', required=True,
                        help='Output filtered TSV file')
    parser.add_argument('--taxa-list',
                        help='Taxa list file path (semicolon-separated format)')
    parser.add_argument('--country-list',
                        help='Country list file path (semicolon-separated format like taxa)')
    parser.add_argument('--marker',
                        help='Marker code to filter by (e.g., COI-5P)')
    parser.add_argument('--kingdom-list', nargs='+',
                        help='Kingdom names to filter by (e.g., --kingdom-list Animalia Plantae)')
    parser.add_argument('--enable-bin-sharing', action='store_true',
                        help='Include BIN_URI sharing expansion')
    parser.add_argument('--filter-species', action='store_true',
                        help='Only include records with valid species names (not null, empty, or "None")')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='Logging level (default: INFO)')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log_level)
    
    try:
        # Run the filtering
        stats = prescoring_filter(
            input_tsv=args.input,
            output_tsv=args.output,
            taxa_list=args.taxa_list,
            country_list=args.country_list,
            marker_code=args.marker,
            kingdom_list=args.kingdom_list,
            enable_bin_sharing=args.enable_bin_sharing,
            filter_species=args.filter_species
        )
        
        # Print summary statistics
        print(f"\nFiltering Summary:")
        print(f"  Input rows: {stats['input_rows']:,}")
        print(f"  Output rows: {stats['output_rows']:,}")
        print(f"  Reduction: {stats['reduction_percentage']:.1f}%")
        print(f"  Processing time: {stats['processing_time_seconds']:.2f}s")
        
        if stats['bin_sharing_used']:
            print(f"  BIN expansion: {stats['initial_matches']:,} -> {stats['final_matches']:,} "
                  f"(factor: {stats['bin_expansion_factor']:.2f})")
        
        logging.info("Pre-scoring filter completed successfully")
        
    except Exception as e:
        logging.error(f"Pre-scoring filter failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
