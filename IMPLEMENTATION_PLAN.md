# Simplified Pre-Scoring Filter Implementation Plan

## Overview

Create a new pre-scoring filter that runs before database creation, keeping the existing target list implementation completely untouched. This filter will reduce dataset size dramatically before expensive database operations.

## Configuration Changes

### Add to `config/config.yml`:
```yaml
# Pre-scoring filter options
ENABLE_PRESCORING_FILTER: false  # Master switch
FILTER_TAXA: false
FILTER_TAXA_LIST: ""
FILTER_COUNTRIES: false  
FILTER_COUNTRY_LIST: ""
FILTER_BINS: false  # Include BIN_URI sharing
PRESCORING_FILTERED_OUTPUT: "results/prescoring_filtered.tsv"
```

## Implementation Phases

### Phase 1: Core Filter Script

**File**: `workflow/scripts/prescoring_filter.py`

**Key Functions**:
```python
def load_taxa_list(file_path: str) -> Set[str]:
    """Load species names from taxa list file (semicolon or comma separated)"""
    
def load_country_list(file_path: str) -> Set[str]: 
    """Load country codes from country list file"""
    
def filter_by_taxa(input_file: str, taxa_set: Set[str], species_column: str) -> Set[str]:
    """Return set of processids matching taxa list"""
    
def filter_by_countries(input_file: str, country_set: Set[str], country_column: str) -> Set[str]:
    """Return set of processids matching country list"""
    
def expand_by_bin_sharing(processids: Set[str], input_file: str, bin_column: str) -> Set[str]:
    """Iteratively expand processid set based on BIN_URI sharing"""
    
def prescoring_filter(
    input_tsv: str,
    output_tsv: str,
    taxa_list: Optional[str] = None,
    country_list: Optional[str] = None, 
    enable_bin_sharing: bool = False
) -> dict:
    """Main filtering function with multi-criteria support"""
```

**Command Line Interface**:
```python
parser.add_argument('--input', required=True, help='Input BOLD TSV file')
parser.add_argument('--output', required=True, help='Output filtered TSV file')
parser.add_argument('--taxa-list', help='Taxa list file path')
parser.add_argument('--country-list', help='Country list file path')
parser.add_argument('--enable-bin-sharing', action='store_true', help='Include BIN_URI sharing')
parser.add_argument('--log-level', default='INFO', help='Logging level')
```

**Checklist**:
- [ ] Basic script structure with argument parsing
- [ ] Taxa list loading (support semicolon and comma separation)
- [ ] Country list loading (CSV format)
- [ ] Single-pass filtering for taxa and countries
- [ ] Output TSV with proper header preservation
- [ ] Progress reporting and logging
- [ ] Error handling for missing files/columns

### Phase 2: BIN Sharing Logic

**BIN Sharing Algorithm**:
```python
def expand_by_bin_sharing(initial_processids: Set[str], input_file: str) -> Set[str]:
    """
    Iteratively expand processid set based on BIN_URI sharing:
    1. Get BIN_URIs from current processids
    2. Find all processids sharing those BIN_URIs
    3. Repeat until no new processids found
    """
    current_processids = initial_processids.copy()
    previous_size = 0
    iteration = 0
    
    while len(current_processids) > previous_size:
        iteration += 1
        previous_size = len(current_processids)
        
        # Get all BIN_URIs for current processids
        current_bins = get_bins_for_processids(current_processids, input_file)
        
        # Get all processids sharing those BIN_URIs
        current_processids = get_processids_for_bins(current_bins, input_file)
        
        logging.info(f"BIN expansion iteration {iteration}: {len(current_processids)} processids")
    
    return current_processids

def get_bins_for_processids(processids: Set[str], input_file: str) -> Set[str]:
    """Extract BIN_URIs for given processids"""
    
def get_processids_for_bins(bins: Set[str], input_file: str) -> Set[str]:
    """Find all processids sharing the given BIN_URIs"""
```

**Checklist**:
- [ ] Iterative BIN_URI expansion algorithm
- [ ] Memory-efficient multi-pass processing
- [ ] BIN_URI extraction and matching functions
- [ ] Progress tracking for expansion iterations
- [ ] Validation of BIN_URI format and existence

### Phase 3: Snakemake Integration

**Add to `workflow/Snakefile`**:
```python
# Conditional rule for pre-scoring filter
if config.get("ENABLE_PRESCORING_FILTER", False):
    rule prescoring_filter:
        input:
            bold_tsv=config["BOLD_TSV"],
            taxa_list=config.get("FILTER_TAXA_LIST") if config.get("FILTER_TAXA", False) else [],
            country_list=config.get("FILTER_COUNTRY_LIST") if config.get("FILTER_COUNTRIES", False) else []
        output:
            config["PRESCORING_FILTERED_OUTPUT"]
        params:
            taxa_arg=lambda wildcards, input: f"--taxa-list {input.taxa_list}" if config.get("FILTER_TAXA", False) and input.taxa_list else "",
            country_arg=lambda wildcards, input: f"--country-list {input.country_list}" if config.get("FILTER_COUNTRIES", False) and input.country_list else "",
            bin_arg="--enable-bin-sharing" if config.get("FILTER_BINS", False) else "",
            log_level=config['LOG_LEVEL']
        log: "logs/prescoring_filter.log"
        conda: "envs/prescoring_filter.yaml"
        shell:
            """
            python workflow/scripts/prescoring_filter.py \
                --input {input.bold_tsv} \
                --output {output} \
                {params.taxa_arg} \
                {params.country_arg} \
                {params.bin_arg} \
                --log-level {params.log_level} \
                2> {log}
            """
```

**Update database creation rule**:
```python
rule create_load_db:
    input:
        bold_tsv=config["PRESCORING_FILTERED_OUTPUT"] if config.get("ENABLE_PRESCORING_FILTER", False) else config["BOLD_TSV"],
        schema=config["SCHEMA"],
        filter_done=config["PRESCORING_FILTERED_OUTPUT"] if config.get("ENABLE_PRESCORING_FILTER", False) else []
    output: config["DB_FILE"]
    params: log_level=config['LOG_LEVEL']
    log: "logs/create_load_db.log"
    conda: "envs/create_load_db.yaml"
    shell:
        "perl workflow/scripts/load_bcdm.pl \
            --tsv {input.bold_tsv} \
            --db {output} \
            --sql {input.schema} \
            --log {params.log_level} \
            --force 2> {log}"
```

**Environment file** `workflow/envs/prescoring_filter.yaml`:
```yaml
name: prescoring_filter
channels:
  - conda-forge
dependencies:
  - python=3.9
  - pandas>=1.3
  - numpy>=1.20
```

**Checklist**:
- [ ] Conditional Snakemake rule with proper input dependencies
- [ ] Parameter handling for optional arguments
- [ ] Database creation rule updated to use filtered data
- [ ] Conda environment file created
- [ ] Rule dependencies and file flow validated

### Phase 4: Testing and Validation

**Test Configuration Examples**:
```yaml
# Example 1: Taxa filtering only
ENABLE_PRESCORING_FILTER: true
FILTER_TAXA: true
FILTER_TAXA_LIST: "resources/test_chironominae_norway_spec_and_syn.csv"
FILTER_COUNTRIES: false
FILTER_BINS: true

# Example 2: Country filtering only  
ENABLE_PRESCORING_FILTER: true
FILTER_TAXA: false
FILTER_COUNTRIES: true
FILTER_COUNTRY_LIST: "resources/european_countries.csv"
FILTER_BINS: false

# Example 3: Combined filtering
ENABLE_PRESCORING_FILTER: true
FILTER_TAXA: true
FILTER_TAXA_LIST: "resources/target_species.csv"
FILTER_COUNTRIES: true
FILTER_COUNTRY_LIST: "resources/target_countries.csv" 
FILTER_BINS: true
```

**Test Scenarios**:
```python
def test_taxa_filtering():
    """Test filtering with species list"""
    
def test_country_filtering():
    """Test filtering with country list"""
    
def test_bin_sharing():
    """Test BIN_URI sharing expansion"""
    
def test_combined_filtering():
    """Test all filters together"""
    
def test_edge_cases():
    """Test empty lists, missing columns, malformed data"""
```

**Checklist**:
- [ ] Test with existing test datasets
- [ ] Validate taxa filtering accuracy
- [ ] Validate country filtering accuracy  
- [ ] Validate BIN sharing expansion
- [ ] Test combined filtering scenarios
- [ ] Performance comparison (before/after filtering)
- [ ] Edge case testing (empty files, missing columns)
- [ ] Integration testing with full Snakemake pipeline

## Processing Logic Flow

```
Input: BOLD TSV (19M rows)
  ↓
Step 1: Load filter lists (taxa, countries)
  ↓
Step 2: Single pass - identify matching processids
  ├── Taxa matches → Set A
  └── Country matches → Set B  
  ↓
Step 3: Combine sets → Set C = A ∪ B
  ↓
Step 4 (if FILTER_BINS=true): BIN sharing expansion
  ├── Get BIN_URIs for Set C
  ├── Find all processids sharing those BIN_URIs
  ├── Repeat until no new processids found
  └── Final Set D
  ↓
Step 5: Second pass - output all records with processids in final set
  ↓
Output: Filtered TSV (reduced size)
```

## Column Assumptions

Based on existing code analysis:
- **Species column**: `species_name` or similar (will auto-detect)
- **Country column**: `country_iso` or `country` (will auto-detect)  
- **BIN column**: `bin_uri` (will auto-detect)
- **ProcessID column**: `processid` (primary key for filtering)

## File Format Support

- **Taxa lists**: CSV with semicolon separation (existing format) or comma separation
- **Country lists**: CSV with country codes
- **Input/Output**: Standard BOLD TSV format (tab-separated, UTF-8)

## Success Criteria

- [ ] Non-invasive: Existing target list workflow unchanged
- [ ] Optional: Can be completely disabled via configuration
- [ ] Performance: >10x speedup for small target lists
- [ ] Accuracy: No data loss, correct filtering logic
- [ ] Maintainable: Clean code structure, good error handling
- [ ] Compatible: Works with existing Snakemake pipeline