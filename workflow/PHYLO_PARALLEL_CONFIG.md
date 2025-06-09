# Parallel Phylogenetic Analysis Configuration
# ============================================

# This file documents the configuration options for the parallel phylogenetic analysis system.
# Add these parameters to your config/config.yml file to customize the analysis.

# Basic Settings
# --------------
PHYLO_ENABLED: True                    # Enable/disable phylogenetic analysis
PHYLO_OUTPUT_DIR: "phylogenies"       # Output directory name (within results/)
PHYLO_MIN_OTUS: 3                     # Minimum species-BIN combinations per family

# Parallelization Settings
# ------------------------
PHYLO_NUM_JOBS: 50                    # Number of SLURM array jobs
PHYLO_FAMILIES_PER_JOB: null          # Fixed families per job (null = auto-distribute)
PHYLO_MAX_CONCURRENT: 10              # Maximum concurrent jobs (%10 in array)
PHYLO_KINGDOMS: ["all"]               # Kingdoms to process (["all"] or ["Animalia", "Plantae"])

# Job Resource Settings
# ---------------------
PHYLO_JOB_MEMORY: "32G"               # Memory per job
PHYLO_JOB_TIME: "24:00:00"            # Time limit per job
PHYLO_PARTITION: "week"               # SLURM partition

# Analysis Methods
# ----------------
PHYLO_ALIGNMENT_METHOD: "mafft"       # Alignment method: "mafft" or "muscle"
PHYLO_TREE_METHOD: "fasttree"         # Tree method: "fasttree" or "iqtree"
PHYLO_BOOTSTRAP: 1000                 # Bootstrap replicates (for iqtree)
PHYLO_NUM_OUTGROUPS: 3               # Number of outgroup sequences

# Output Options
# --------------
PHYLO_GENERATE_PDFS: True             # Generate tree visualization PDFs
PHYLO_CLEANUP_INTERMEDIATES: True     # Clean up intermediate files

# Performance Recommendations
# ===========================

# For 3000 families:
# ------------------
# PHYLO_NUM_JOBS: 50-100               # Balance between parallelism and queue limits
# PHYLO_MAX_CONCURRENT: 10-20          # Prevent overwhelming the cluster
# PHYLO_JOB_MEMORY: "32G"              # Sufficient for most families
# PHYLO_JOB_TIME: "24:00:00"           # Allow for complex families

# For testing (small subset):
# ---------------------------
# PHYLO_NUM_JOBS: 5
# PHYLO_MAX_CONCURRENT: 5
# PHYLO_FAMILIES_PER_JOB: 2

# For high-performance clusters:
# ------------------------------
# PHYLO_NUM_JOBS: 100
# PHYLO_MAX_CONCURRENT: 30
# PHYLO_JOB_MEMORY: "64G"
# PHYLO_PARTITION: "compute"

# Method Selection Guidelines
# ===========================

# Alignment Methods:
# -----------------
# MAFFT: Faster, good for most cases, handles reverse complements
# MUSCLE: Slower but sometimes more accurate for divergent sequences

# Tree Methods:
# ------------
# FastTree: Much faster (~10x), good for large datasets, approximate ML
# IQ-TREE: Slower but more accurate, model selection, better bootstrap

# Kingdom Selection:
# -----------------
# ["all"]: Process all kingdoms (recommended for comprehensive analysis)
# ["Animalia"]: Process only animal families
# ["Plantae"]: Process only plant families
# ["Animalia", "Fungi"]: Process multiple specific kingdoms

# Expected Performance
# ===================

# With recommended settings (50 jobs, 32G memory, FastTree):
# - Total runtime: 4-8 hours for 3000 families
# - Individual job runtime: 2-6 hours depending on family complexity
# - Speed improvement: ~20-30x faster than sequential processing
# - Resource usage: 50 concurrent jobs × 32G = 1.6TB peak memory across cluster

# Output Files Generated
# =====================

# Per family:
# -----------
# - {family}.treefile                           # Newick format phylogenetic tree
# - {family}_tree.pdf                          # Tree visualization with BAGS grades
# - {family}_curation_checklist.pdf            # Curation checklist organized by BAGS grades
# - tree_metadata.json                         # Analysis metadata and results
# - curation_data.json                         # Raw curation data

# Global summaries:
# ----------------
# - phylogenetic_analysis_summary.json         # Complete analysis summary
# - failed_families.txt                        # List of families that failed processing

# Usage Examples
# ==============

# 1. Full analysis with default settings:
# ---------------------------------------
# sbatch workflow/sh-phylo-parallel.sh

# 2. Test run with subset:
# -----------------------
# Edit config.yml:
#   PHYLO_NUM_JOBS: 5
#   PHYLO_FAMILIES_PER_JOB: 2
# sbatch workflow/sh-phylo-parallel.sh

# 3. High-performance run:
# -----------------------
# Edit config.yml:
#   PHYLO_NUM_JOBS: 100
#   PHYLO_MAX_CONCURRENT: 30
#   PHYLO_JOB_MEMORY: "64G"
#   PHYLO_TREE_METHOD: "iqtree"
# sbatch workflow/sh-phylo-parallel.sh

# 4. Animals only:
# ---------------
# Edit config.yml:
#   PHYLO_KINGDOMS: ["Animalia"]
# sbatch workflow/sh-phylo-parallel.sh

# Monitoring and Troubleshooting
# ==============================

# Monitor job progress:
# --------------------
# squeue -u $USER | grep phylo                 # Check running jobs
# ls logs/phylo_*.out | wc -l                  # Count completed jobs
# tail -f logs/phylogenetic_analysis_parallel.log  # Monitor main workflow

# Check results:
# -------------
# ls results/phylogenies/ | wc -l              # Count processed families
# python workflow/scripts/consolidate_phylo_results.py results/phylo_batches results/phylogenies

# Resume failed jobs:
# ------------------
# The system automatically uses checkpoints to resume failed families
# Simply rerun: sbatch workflow/sh-phylo-parallel.sh

# Troubleshooting common issues:
# -----------------------------
# 1. "No batches to process": Check PHYLO_ENABLED=True and database exists
# 2. Jobs fail immediately: Check conda environment and tool availability
# 3. Memory errors: Increase PHYLO_JOB_MEMORY or reduce families per job
# 4. Time limit exceeded: Increase PHYLO_JOB_TIME or use FastTree instead of IQ-TREE
# 5. Queue limits: Reduce PHYLO_MAX_CONCURRENT
