#!/bin/bash

# Setup script for independent haplotype analysis
# Run this before submitting the SLURM job to ensure everything is ready

echo "=============================================="
echo "BOLD Haplotype Analysis - Setup & Validation"
echo "=============================================="
echo ""

# Configuration (modify these paths for your setup)
SCRIPT_DIR="workflow/scripts"
DB_FILE="results/bold.db"
CONDA_ENV="bold-curation"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Helper functions
check_ok() {
    echo -e "${GREEN}✓${NC} $1"
}

check_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
}

check_error() {
    echo -e "${RED}✗${NC} $1"
}

# 1. Check if we're in the right directory
echo "1. Checking working directory..."
if [[ -f "workflow/bold-ranker.smk" ]]; then
    check_ok "Found BOLD ranker Snakemake file"
else
    check_error "Not in BOLD library curation directory"
    echo "   Please run this script from the bold-library-curation root directory"
    exit 1
fi

# 2. Check if optimized script exists
echo ""
echo "2. Checking optimized haplotype script..."
if [[ -f "assess_haplotypes_optimized.pl" ]]; then
    check_ok "Found optimized script: assess_haplotypes_optimized.pl"
else
    check_error "Optimized script not found: assess_haplotypes_optimized.pl"
    echo "   Please ensure you've downloaded the optimized script to the current directory"
    exit 1
fi

# 3. Check database
echo ""
echo "3. Checking database..."
if [[ -f "$DB_FILE" ]]; then
    DB_SIZE=$(du -h "$DB_FILE" | cut -f1)
    check_ok "Database found: $DB_FILE ($DB_SIZE)"
    
    # Check if database has required tables and data
    echo "   Validating database structure..."
    
    # Check if bold table exists
    BOLD_TABLE=$(sqlite3 "$DB_FILE" "SELECT name FROM sqlite_master WHERE type='table' AND name='bold';" 2>/dev/null)
    if [[ "$BOLD_TABLE" == "bold" ]]; then
        check_ok "Table 'bold' exists"
    else
        check_error "Table 'bold' not found"
        exit 1
    fi
    
    # Check record count
    RECORD_COUNT=$(sqlite3 "$DB_FILE" "SELECT COUNT(*) FROM bold WHERE nuc IS NOT NULL AND nuc != '' AND species IS NOT NULL AND species != '';" 2>/dev/null)
    if [[ $RECORD_COUNT -gt 0 ]]; then
        check_ok "Found $RECORD_COUNT records with sequences and species"
    else
        check_error "No valid records found in database"
        exit 1
    fi
    
    # Check for required columns
    COLUMNS=$(sqlite3 "$DB_FILE" "PRAGMA table_info(bold);" | cut -d'|' -f2)
    REQUIRED_COLS=("recordid" "species" "nuc" "bin_uri")
    for col in "${REQUIRED_COLS[@]}"; do
        if echo "$COLUMNS" | grep -q "$col"; then
            check_ok "Column '$col' present"
        else
            check_error "Required column '$col' missing"
            exit 1
        fi
    done
    
else
    check_error "Database not found: $DB_FILE"
    echo "   Please run the main BOLD pipeline first to create the database"
    exit 1
fi

# 4. Check conda environment
echo ""
echo "4. Checking conda environment..."
if command -v conda >/dev/null 2>&1; then
    check_ok "Conda available"
    
    # Check if environment exists
    if conda env list | grep -q "$CONDA_ENV"; then
        check_ok "Environment '$CONDA_ENV' exists"
    else
        check_error "Environment '$CONDA_ENV' not found"
        echo "   Available environments:"
        conda env list | grep -v "^#" | sed 's/^/     /'
        exit 1
    fi
else
    check_error "Conda not available"
    echo "   Please ensure conda/miniconda is installed and in PATH"
    exit 1
fi

# 5. Check Perl dependencies
echo ""
echo "5. Checking Perl dependencies..."

# Activate conda environment for testing
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$CONDA_ENV" 2>/dev/null

REQUIRED_MODULES=("DBI" "Getopt::Long" "Time::HiRes" "Parallel::ForkManager" "Digest::MD5")

for module in "${REQUIRED_MODULES[@]}"; do
    if perl -e "use $module; print 'OK';" >/dev/null 2>&1; then
        check_ok "Perl module '$module' available"
    else
        check_error "Perl module '$module' missing"
        
        if [[ "$module" == "Parallel::ForkManager" ]]; then
            echo "   To install: cpan Parallel::ForkManager"
            echo "   Or: conda install perl-parallel-forkmanager"
        elif [[ "$module" == "Digest::MD5" ]]; then
            echo "   To install: cpan Digest::MD5"
            echo "   Or: conda install perl-digest-md5"
        fi
    fi
done

# 6. Check system resources
echo ""
echo "6. Checking system resources..."

# Check available memory
if command -v free >/dev/null 2>&1; then
    TOTAL_MEM=$(free -g | grep Mem | awk '{print $2}')
    AVAIL_MEM=$(free -g | grep Mem | awk '{print $7}')
    check_ok "System memory: ${AVAIL_MEM}GB available / ${TOTAL_MEM}GB total"
    
    if [[ $AVAIL_MEM -lt 8 ]]; then
        check_warn "Low available memory. Consider reducing --processes parameter"
    fi
else
    check_warn "Cannot check memory usage (free command not available)"
fi

# Check CPU cores
if [[ -f /proc/cpuinfo ]]; then
    CPU_CORES=$(grep -c ^processor /proc/cpuinfo)
    check_ok "CPU cores: $CPU_CORES"
else
    check_warn "Cannot determine CPU core count"
fi

# Check disk space
DISK_AVAIL=$(df -BG . | tail -1 | awk '{print $4}' | sed 's/G//')
if [[ $DISK_AVAIL -gt 20 ]]; then
    check_ok "Disk space: ${DISK_AVAIL}GB available"
else
    check_warn "Low disk space: ${DISK_AVAIL}GB available"
fi

# 7. Test script execution (dry run)
echo ""
echo "7. Testing script execution..."

# Quick syntax check
if perl -c assess_haplotypes_optimized.pl >/dev/null 2>&1; then
    check_ok "Script syntax valid"
else
    check_error "Script has syntax errors"
    perl -c assess_haplotypes_optimized.pl
    exit 1
fi

# Test with --help flag
if perl assess_haplotypes_optimized.pl --help >/dev/null 2>&1; then
    check_ok "Script help function works"
else
    check_warn "Script help function not available"
fi

# 8. Create necessary directories
echo ""
echo "8. Preparing output directories..."

DIRS=("logs" "results")
for dir in "${DIRS[@]}"; do
    if [[ ! -d "$dir" ]]; then
        mkdir -p "$dir"
        check_ok "Created directory: $dir"
    else
        check_ok "Directory exists: $dir"
    fi
done

# 9. Generate recommended SLURM parameters
echo ""
echo "9. Recommended SLURM parameters for your system:"
echo ""

# Calculate optimal settings based on record count and system
if [[ -n "$RECORD_COUNT" && -n "$CPU_CORES" ]]; then
    # Suggest memory based on record count (rough estimate: 40MB per 1000 records)
    SUGGESTED_MEM=$(( (RECORD_COUNT / 1000 * 40) + 20 ))  # Add 20GB base
    if [[ $SUGGESTED_MEM -lt 60 ]]; then SUGGESTED_MEM=60; fi
    if [[ $SUGGESTED_MEM -gt 200 ]]; then SUGGESTED_MEM=200; fi
    
    # Suggest CPU cores (use all available, but cap at 32 for efficiency)
    SUGGESTED_CPU=$CPU_CORES
    if [[ $SUGGESTED_CPU -gt 32 ]]; then SUGGESTED_CPU=32; fi
    
    # Estimate time based on record count
    EST_MINUTES=$(( RECORD_COUNT / 1000 / SUGGESTED_CPU * 2 ))  # Conservative estimate
    if [[ $EST_MINUTES -lt 60 ]]; then EST_MINUTES=60; fi  # Minimum 1 hour
    EST_HOURS=$(( (EST_MINUTES + 59) / 60 ))  # Round up to nearest hour
    
    echo "   For $RECORD_COUNT records:"
    echo "   --mem=${SUGGESTED_MEM}G"
    echo "   --cpus-per-task=$SUGGESTED_CPU"
    echo "   --time=${EST_HOURS}:00:00"
    echo ""
    echo "   Edit these values in sh-haplotype-optimized.sh before submitting"
else
    echo "   Unable to calculate optimal parameters"
fi

# 10. Final summary
echo ""
echo "=============================================="
echo "SETUP VALIDATION COMPLETE"
echo "=============================================="

# Count any errors or warnings
ERROR_COUNT=$(grep -c "✗" <<< "$OUTPUT" 2>/dev/null || echo 0)
WARN_COUNT=$(grep -c "⚠" <<< "$OUTPUT" 2>/dev/null || echo 0)

if [[ $ERROR_COUNT -eq 0 ]]; then
    echo -e "${GREEN}✓ All checks passed!${NC}"
    echo ""
    echo "Ready to submit SLURM job:"
    echo "   sbatch sh-haplotype-optimized.sh"
    echo ""
    echo "Monitor progress with:"
    echo "   squeue -u \$USER"
    echo "   tail -f job_haplotype_<jobid>.out"
else
    echo -e "${RED}✗ $ERROR_COUNT errors found${NC}"
    echo "Please fix the errors above before proceeding"
    exit 1
fi

if [[ $WARN_COUNT -gt 0 ]]; then
    echo -e "${YELLOW}⚠ $WARN_COUNT warnings - review before proceeding${NC}"
fi

conda deactivate 2>/dev/null