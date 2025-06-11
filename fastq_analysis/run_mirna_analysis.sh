#!/bin/bash

# =============================================================================
# Complete miRNA-seq Analysis Runner
# =============================================================================
# This script runs the complete analysis pipeline step by step
# Usage: /path/to/scripts/run_mirna_analysis.sh
# =============================================================================

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Get the directory where the script was called from (data directory)
DATA_DIR="$(pwd)"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to print step headers
print_step() {
    echo -e "\n${BLUE}════════════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}STEP $1: $2${NC}"
    echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}\n"
}

# Function to check if command succeeded
check_status() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ $1 completed successfully${NC}\n"
    else
        echo -e "${RED}✗ $1 failed${NC}\n"
        exit 1
    fi
}

# Start
clear
echo -e "${GREEN}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║           miRNA-seq Analysis Pipeline Runner                  ║${NC}"
echo -e "${GREEN}║                                                               ║${NC}"
echo -e "${GREEN}║  Samples:                                                     ║${NC}"
echo -e "${GREEN}║  - IL21602-001: Electroporated (Let-7-5p knockdown)         ║${NC}"
echo -e "${GREEN}║  - IL21603-001: Not electroporated (Control)                 ║${NC}"
echo -e "${GREEN}╚═══════════════════════════════════════════════════════════════╝${NC}"
echo
echo -e "${BLUE}Script directory: ${SCRIPT_DIR}${NC}"
echo -e "${BLUE}Data directory: ${DATA_DIR}${NC}"
echo

# =============================================================================
# STEP 1: Check Prerequisites
# =============================================================================
print_step "1" "Checking Prerequisites"

echo "Checking required tools..."

# Check for R
if command -v R >/dev/null 2>&1; then
    echo -e "✓ R is installed ($(R --version | head -n1))"
else
    echo -e "${RED}✗ R is not installed. Please install R first.${NC}"
    exit 1
fi

# Check for Python tools
if command -v fastqc >/dev/null 2>&1; then
    echo "✓ FastQC is installed"
else
    echo -e "${YELLOW}⚠ FastQC not found. Install with: brew install fastqc${NC}"
fi

if command -v multiqc >/dev/null 2>&1 || python3 -m multiqc --version >/dev/null 2>&1; then
    echo "✓ MultiQC is installed"
else
    echo -e "${YELLOW}⚠ MultiQC not found. Install with: pip3 install multiqc${NC}"
fi

if command -v cutadapt >/dev/null 2>&1 || python3 -m cutadapt --version >/dev/null 2>&1; then
    echo "✓ Cutadapt is installed"
else
    echo -e "${YELLOW}⚠ Cutadapt not found. Install with: pip3 install cutadapt${NC}"
fi

# =============================================================================
# STEP 2: Install R Packages
# =============================================================================
print_step "2" "Installing/Checking R Packages"

echo "Running R package installation check..."
Rscript "${SCRIPT_DIR}/01_install_packages.R"
check_status "R package installation"

# =============================================================================
# STEP 3: Extract Data
# =============================================================================
print_step "3" "Extracting and Organizing Data"

# Check for TAR file in data directory
if [ -f "${DATA_DIR}/XTRIA.20250601_163957.ILLUMINA_DATA.1-of-1.tar" ]; then
    echo "Extracting TAR file..."
    cd "${DATA_DIR}"
    tar -xf XTRIA.20250601_163957.ILLUMINA_DATA.1-of-1.tar
    check_status "Data extraction"
    
    echo "Files extracted. Checking structure..."
    
    # Check what was extracted
    if [ -d "Electroporated" ] && [ -d "Notelectroporated" ]; then
        echo "✓ Found sample directories"
        
        # List files in each directory
        echo "Electroporated samples:"
        ls -la Electroporated/ILLUMINA_DATA/*.fastq.gz 2>/dev/null | head -5
        echo "Notelectroporated samples:"
        ls -la Notelectroporated/ILLUMINA_DATA/*.fastq.gz 2>/dev/null | head -5
    else
        echo "Checking for FASTQ files in main directory..."
        ls -la *.fastq.gz 2>/dev/null | head -10
    fi
else
    echo -e "${YELLOW}TAR archive not found. Assuming data is already extracted.${NC}"
fi

# Copy mature.fa if it exists in data directory
if [ -f "${DATA_DIR}/mature.fa" ] && [ ! -f "${DATA_DIR}/complete_analysis/4_mapping/mirbase/mature.fa" ]; then
    echo "Found mature.fa, will copy during analysis..."
fi

# =============================================================================
# STEP 4: Run Preprocessing (FastQC/MultiQC)
# =============================================================================
print_step "4" "Running Quality Control (FastQC/MultiQC)"

if command -v fastqc >/dev/null 2>&1; then
    Rscript "${SCRIPT_DIR}/02_run_preprocessing.R" "${DATA_DIR}"
    check_status "Quality control analysis"
else
    echo -e "${YELLOW}Skipping FastQC/MultiQC (tools not installed)${NC}"
fi

# =============================================================================
# STEP 5: Run Main Analysis Pipeline
# =============================================================================
print_step "5" "Running miRNA Analysis Pipeline"

echo "Starting main analysis..."
Rscript "${SCRIPT_DIR}/03_run_main_analysis.R" "${DATA_DIR}" "${SCRIPT_DIR}"
check_status "miRNA analysis"

# =============================================================================
# STEP 6: Generate Visualizations
# =============================================================================
print_step "6" "Generating Visualization Plots"

echo "Creating visualization plots..."
Rscript "${SCRIPT_DIR}/04_create_visualizations.R" "${DATA_DIR}"
check_status "Plot generation"

# =============================================================================
# STEP 7: Generate Summary Report
# =============================================================================
print_step "7" "Generating Summary Report"

echo "Creating HTML summary report..."
Rscript "${SCRIPT_DIR}/05_generate_report.R" "${DATA_DIR}"
check_status "Report generation"

# =============================================================================
# FINAL SUMMARY
# =============================================================================
echo
echo -e "${GREEN}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║                    ANALYSIS COMPLETE!                         ║${NC}"
echo -e "${GREEN}╚═══════════════════════════════════════════════════════════════╝${NC}"
echo
echo -e "${BLUE}Results saved in: ${DATA_DIR}/complete_analysis/${NC}"
echo
echo "Key files:"
echo "  • Excel results: complete_analysis/complete_analysis_results.xlsx"
echo "  • HTML report: complete_analysis/analysis_report.html"
echo "  • Plots: complete_analysis/plots/"
echo "  • QC reports: complete_analysis/*/multiqc_*.html"
echo
echo -e "${GREEN}Open the Excel file to explore all results!${NC}"
echo

# Open results folder (macOS)
if [[ "$OSTYPE" == "darwin"* ]]; then
    open "${DATA_DIR}/complete_analysis/"
fi