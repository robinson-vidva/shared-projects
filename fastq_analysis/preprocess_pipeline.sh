#!/bin/bash

# =============================================================================
# Complete miRNA-seq Preprocessing Pipeline - Updated for Directory Structure
# =============================================================================
# This script handles all command-line preprocessing steps:
# 1. FastQC on raw data
# 2. Adapter trimming (if needed)
# 3. FastQC on trimmed data
# 4. MultiQC reports
# =============================================================================

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_message() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check if required tools are installed
check_tools() {
    local missing_tools=()
    
    command -v fastqc >/dev/null 2>&1 || missing_tools+=("fastqc")
    
    # Check for multiqc and cutadapt as commands or Python modules
    if ! command -v multiqc >/dev/null 2>&1; then
        if ! python3 -m multiqc --version >/dev/null 2>&1; then
            missing_tools+=("multiqc")
        fi
    fi
    
    if ! command -v cutadapt >/dev/null 2>&1; then
        if ! python3 -m cutadapt --version >/dev/null 2>&1; then
            missing_tools+=("cutadapt")
        fi
    fi
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        print_error "Missing required tools: ${missing_tools[*]}"
        echo "Please install missing tools:"
        echo "  - FastQC: brew install fastqc"
        echo "  - MultiQC: pip3 install multiqc"
        echo "  - Cutadapt: pip3 install cutadapt"
        return 1
    fi
    return 0
}

# Wrapper functions for Python tools
run_multiqc() {
    if command -v multiqc >/dev/null 2>&1; then
        multiqc "$@"
    else
        python3 -m multiqc "$@"
    fi
}

run_cutadapt() {
    if command -v cutadapt >/dev/null 2>&1; then
        cutadapt "$@"
    else
        python3 -m cutadapt "$@"
    fi
}

# Main pipeline
main() {
    # Set up directories
    BASE_DIR="/Users/rvidva/Documents/datasets/XTRIA"
    OUTPUT_DIR="${BASE_DIR}/complete_analysis"
    
    print_message "Starting miRNA-seq preprocessing pipeline"
    print_message "Base directory: ${BASE_DIR}"
    print_message "Output directory: ${OUTPUT_DIR}"
    
    # Check tools
    if ! check_tools; then
        exit 1
    fi
    
    # Create output directories
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}/1_fastqc_raw"
    mkdir -p "${OUTPUT_DIR}/2_trimmed_reads"
    mkdir -p "${OUTPUT_DIR}/3_fastqc_trimmed"
    mkdir -p "${OUTPUT_DIR}/logs"
    
    cd "${BASE_DIR}"
    
    # =========================================================================
    # Step 1: Extract TAR file if needed
    # =========================================================================
    print_message "Step 1: Checking data extraction"
    
    # Check if directories exist
    if [ ! -d "Electroporated" ] || [ ! -d "Notelectroporated" ]; then
        if [ -f "XTRIA.20250601_163957.ILLUMINA_DATA.1-of-1.tar" ]; then
            print_message "Extracting TAR file..."
            tar -xf XTRIA.20250601_163957.ILLUMINA_DATA.1-of-1.tar
        else
            print_error "TAR file not found and directories don't exist!"
            exit 1
        fi
    else
        print_message "Data already extracted in subdirectories"
    fi
    
    # =========================================================================
    # Step 2: FastQC on raw data
    # =========================================================================
    print_message "Step 2: Running FastQC on raw data"
    
    # Process all raw FASTQ files from both directories
    for dir in Electroporated Notelectroporated; do
        if [ -d "$dir" ]; then
            print_message "Processing files in $dir/"
            for file in "$dir"/XTRIA_*_R[12].fastq.gz; do
                if [ -f "$file" ] && [[ ! "$file" =~ "trimmed" ]]; then
                    print_message "Running FastQC on: $(basename "$file")"
                    fastqc -o "${OUTPUT_DIR}/1_fastqc_raw" -t 2 "$file" 2>&1 | \
                        tee -a "${OUTPUT_DIR}/logs/fastqc_raw.log"
                fi
            done
        fi
    done
    
    # Run MultiQC on raw FastQC results
    print_message "Running MultiQC on raw data"
    cd "${OUTPUT_DIR}/1_fastqc_raw"
    if ls *.zip 1> /dev/null 2>&1; then
        run_multiqc . -n multiqc_raw_report 2>&1 | tee -a "${OUTPUT_DIR}/logs/multiqc_raw.log"
    else
        print_warning "No FastQC results found for raw data"
    fi
    cd "${BASE_DIR}"
    
    # =========================================================================
    # Step 3: Check for trimmed files or perform trimming
    # =========================================================================
    print_message "Step 3: Checking adapter trimming"
    
    # Check if trimmed files already exist in the directories
    trimmed_exists=true
    for dir in Electroporated Notelectroporated; do
        if [ ! -f "$dir"/*_R1_trimmed.fastq.gz ]; then
            trimmed_exists=false
            break
        fi
    done
    
    if [ "$trimmed_exists" = true ]; then
        print_message "Using pre-trimmed files from sequencing facility"
        
        # Copy trimmed files to our output directory
        for dir in Electroporated Notelectroporated; do
            cp "$dir"/*_trimmed.fastq.gz "${OUTPUT_DIR}/2_trimmed_reads/" 2>/dev/null || true
            cp "$dir"/*_trimmed_stats.txt "${OUTPUT_DIR}/2_trimmed_reads/" 2>/dev/null || true
        done
    else
        print_message "Would perform adapter trimming, but trimmed files already exist"
    fi
    
    # =========================================================================
    # Step 4: FastQC on trimmed data
    # =========================================================================
    print_message "Step 4: Running FastQC on trimmed data"
    
    # Process trimmed files from both directories
    for dir in Electroporated Notelectroporated; do
        if [ -d "$dir" ]; then
            for file in "$dir"/*_trimmed.fastq.gz; do
                if [ -f "$file" ]; then
                    print_message "Running FastQC on: $(basename "$file")"
                    fastqc -o "${OUTPUT_DIR}/3_fastqc_trimmed" -t 2 "$file" 2>&1 | \
                        tee -a "${OUTPUT_DIR}/logs/fastqc_trimmed.log"
                fi
            done
        fi
    done
    
    # Also process files in output directory if copied there
    for file in "${OUTPUT_DIR}/2_trimmed_reads/"*_trimmed.fastq.gz; do
        if [ -f "$file" ]; then
            print_message "Running FastQC on: $(basename "$file")"
            fastqc -o "${OUTPUT_DIR}/3_fastqc_trimmed" -t 2 "$file" 2>&1 | \
                tee -a "${OUTPUT_DIR}/logs/fastqc_trimmed.log"
        fi
    done
    
    # Run MultiQC on trimmed FastQC results
    print_message "Running MultiQC on trimmed data"
    cd "${OUTPUT_DIR}/3_fastqc_trimmed"
    if ls *.zip 1> /dev/null 2>&1; then
        run_multiqc . -n multiqc_trimmed_report 2>&1 | tee -a "${OUTPUT_DIR}/logs/multiqc_trimmed.log"
    else
        print_warning "No FastQC results found for trimmed data"
    fi
    cd "${BASE_DIR}"
    
    # =========================================================================
    # Step 5: Generate preprocessing summary
    # =========================================================================
    print_message "Step 5: Generating preprocessing summary"
    
    SUMMARY_FILE="${OUTPUT_DIR}/preprocessing_summary.txt"
    
    {
        echo "miRNA-seq Preprocessing Summary"
        echo "==============================="
        echo "Date: $(date)"
        echo ""
        echo "Sample Information:"
        echo "- IL21602-001: Electroporated EVs"
        echo "- IL21603-001: Non-electroporated control EVs"
        echo ""
        echo "Directory Structure:"
        echo "- Electroporated/ : Contains IL21602-001 files"
        echo "- Notelectroporated/ : Contains IL21603-001 files"
        echo ""
        echo "File Statistics:"
        echo "----------------"
        
        # Count reads in trimmed files from both directories
        for dir in Electroporated Notelectroporated; do
            if [ -d "$dir" ]; then
                echo ""
                echo "$dir samples:"
                for file in "$dir"/*_R1_trimmed.fastq.gz; do
                    if [ -f "$file" ]; then
                        count=$(zcat "$file" | wc -l | awk '{print $1/4}')
                        echo "  $(basename "$file"): $(printf "%'d" $count) reads"
                    fi
                done
            fi
        done
        
        echo ""
        echo "Quality Control Reports:"
        echo "- Raw data: ${OUTPUT_DIR}/1_fastqc_raw/multiqc_raw_report.html"
        echo "- Trimmed data: ${OUTPUT_DIR}/3_fastqc_trimmed/multiqc_trimmed_report.html"
        echo ""
        echo "Next Steps:"
        echo "1. Review the MultiQC reports"
        echo "2. Run the R analysis script for mapping and differential expression"
        echo "3. Perform pathway analysis with the results"
    } > "${SUMMARY_FILE}"
    
    cat "${SUMMARY_FILE}"
    
    # =========================================================================
    # Quick Let-7 check
    # =========================================================================
    print_message "Bonus: Quick Let-7-5p check"
    
    echo ""
    echo "Let-7 signature sequence (GAGGTAG) counts:"
    
    # Check Electroporated sample
    electro_file="Electroporated/XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz"
    if [ -f "$electro_file" ]; then
        count=$(zcat "$electro_file" | grep -c "GAGGTAG" || true)
        echo "IL21602-001 (Electroporated): ${count} reads with GAGGTAG motif"
    fi
    
    # Check Notelectroporated sample
    notelectro_file="Notelectroporated/XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz"
    if [ -f "$notelectro_file" ]; then
        count=$(zcat "$notelectro_file" | grep -c "GAGGTAG" || true)
        echo "IL21603-001 (Notelectroporated): ${count} reads with GAGGTAG motif"
    fi
    
    echo ""
    echo "Additional Let-7-5p sequence checks:"
    
    # Check for specific Let-7 sequences
    for seq_name in "let-7a-5p:TGAGGTAGTAGGTTGTATAGTT" "let-7b-5p:TGAGGTAGTAGGTTGTGTGGTT" "let-7c-5p:TGAGGTAGTAGGTTGTATGGTT"; do
        seq=$(echo $seq_name | cut -d: -f2)
        name=$(echo $seq_name | cut -d: -f1)
        
        echo ""
        echo "$name sequence ($seq):"
        
        if [ -f "$electro_file" ]; then
            count=$(zcat "$electro_file" | grep -c "$seq" || true)
            echo "  Electroporated: ${count} reads"
        fi
        
        if [ -f "$notelectro_file" ]; then
            count=$(zcat "$notelectro_file" | grep -c "$seq" || true)
            echo "  Notelectroporated: ${count} reads"
        fi
    done
    
    print_message "Preprocessing complete!"
    print_message "Results saved in: ${OUTPUT_DIR}"
    
    # Create symlinks in base directory for R script compatibility
    print_message "Creating symlinks for R analysis compatibility"
    
    # Create symlinks in base directory pointing to files in subdirectories
    ln -sf "Electroporated/XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz" \
           "${BASE_DIR}/XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz" 2>/dev/null || true
    
    ln -sf "Notelectroporated/XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz" \
           "${BASE_DIR}/XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz" 2>/dev/null || true
}

# Run the main pipeline
main "$@"