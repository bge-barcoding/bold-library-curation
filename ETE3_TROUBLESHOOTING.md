# ETE3 PDF Generation Troubleshooting Guide

## Problem
The phylogenetic analysis pipeline completes successfully but PDF tree visualizations are not generated. Log shows "Generate PDFs: False" even when configured as True.

## Root Cause
ETE3 library fails to import or initialize in headless HPC environments due to missing Qt/X11 dependencies or improper display setup.

## Quick Diagnosis

Run the validation script to test your environment:
```bash
cd /path/to/bold-library-curation
conda activate bold-curation
python workflow/scripts/validate_ete3.py
```

For detailed diagnosis:
```bash
python workflow/scripts/test_ete3_import.py
```

## Common Issues and Solutions

### 1. ETE3 Not Installed
**Symptoms:** ImportError: No module named 'ete3'
**Solution:**
```bash
conda install -c etetoolkit ete3
```

### 2. Missing Qt/X11 Dependencies
**Symptoms:** Qt/PyQt import errors, display connection failures
**Solution:**
```bash
# Install missing packages
conda install qt-main xorg-libx11 xorg-libxext xorg-libxrender
conda install xvfb fontconfig freetype
```

### 3. Virtual Display Issues
**Symptoms:** "cannot connect to X server" errors
**Solution:**
Start virtual display before running:
```bash
# Start Xvfb in background
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
```

### 4. Environment Variable Setup
**Symptoms:** Qt platform plugin errors
**Solution:**
Ensure these are set:
```bash
export QT_QPA_PLATFORM=offscreen
export DISPLAY=:99
export MPLBACKEND=Agg
export QT_DEBUG_PLUGINS=0
export QTWEBENGINE_DISABLE_SANDBOX=1
```

## Updated Environment Files

The following files have been updated to include all necessary dependencies:

### environment.yml
Added packages:
- `qt-main` - Qt platform plugins
- `xorg-libx11` - Core X11 library  
- `xvfb` - Virtual framebuffer
- `fontconfig` - Font configuration
- `freetype` - Font rendering

### workflow/envs/phylogenetic_analysis.yaml
Same additions as environment.yml for consistency.

## Testing Your Fix

After updating the environment:

1. **Recreate the conda environment:**
   ```bash
   conda env remove -n bold-curation
   conda env create -f environment.yml
   conda activate bold-curation
   ```

2. **Test ETE3 setup:**
   ```bash
   python workflow/scripts/validate_ete3.py
   ```

3. **Run a small phylogenetic test:**
   ```bash
   # Test on a single family
   python workflow/scripts/phylo_pipeline.py \
     --database results/bold.db \
     --family-names Coenagrionidae \
     --generate-pdfs \
     --output-dir test_phylo
   ```

## Verification

If PDF generation is working, you should see:
- Log shows "✓ ETE3 loaded with headless backend successfully"
- Tree visualization PDFs are created in family directories
- No "Generate PDFs: False" in the logs

## Additional Notes

- The pipeline generates TWO types of PDFs:
  - **Tree visualization PDFs** (ETE3) - these were failing
  - **Curation checklist PDFs** (ReportLab) - these work fine
- SLURM jobs use the main `bold-curation` environment, not the phylogenetic_analysis.yaml
- Virtual display setup is now automated in phylo_array_job.sh

## Still Having Issues?

Check the detailed diagnostic output from:
```bash
python workflow/scripts/test_ete3_import.py
```

This will show exactly which component is failing and provide specific fix suggestions.
