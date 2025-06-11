#!/bin/bash

# =============================================================================
# Proper adapter trimming for NextFlex Small RNA v4
# =============================================================================

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}=== NextFlex v4 Small RNA Adapter Trimming ===${NC}"
echo
echo "This script properly trims NextFlex v4 adapters for miRNA analysis"
echo

# Check if cutadapt is available
if ! command -v cutadapt &> /dev/null; then
    echo -e "${RED}Error: cutadapt not found. Install with: pip install cutadapt${NC}"
    exit 1
fi

# Set directories
DATA_DIR="${1:-.}"
cd "$DATA_DIR"

# Create output directory
OUTPUT_DIR="properly_trimmed"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

# NextFlex v4 adapter sequences
ADAPTER_3P="TGGAATTCTCGGGTGCCAAGG"  # 3' adapter
ADAPTER_5P="GTTCAGAGTTCTACAGTCCGACGATC"  # 5' adapter (if needed)

echo -e "${BLUE}Trimming parameters:${NC}"
echo "  3' adapter: $ADAPTER_3P"
echo "  Minimum length: 18 nt (for miRNAs)"
echo "  Maximum length: 35 nt (for miRNAs)"
echo "  Minimum quality: 20"
echo

# Find original (non-trimmed) R1 files
echo -e "${BLUE}Finding original R1 files...${NC}"

# Look in multiple possible locations
R1_FILES=""
for dir in "original_files" "." "Electroporated/ILLUMINA_DATA" "Notelectroporated/ILLUMINA_DATA"; do
    if [ -d "$dir" ]; then
        FILES=$(find "$dir" -name "*_R1.fastq.gz" ! -name "*trimmed*" -type f 2>/dev/null)
        if [ -n "$FILES" ]; then
            R1_FILES="$R1_FILES$FILES"$'\n'
        fi
    fi
done

# Remove empty lines and duplicates based on filename
R1_FILES=$(echo "$R1_FILES" | grep -v "^$" | while read -r file; do
    basename "$file"
done | sort | uniq | while read -r filename; do
    # Find the first occurrence of this filename
    echo "$R1_FILES" | grep -v "^$" | grep "/$filename$" | head -1
done)

if [ -z "$R1_FILES" ]; then
    echo -e "${RED}No original R1 files found!${NC}"
    echo "Looked in: original_files/, current directory, and sample subdirectories"
    echo
    echo "If files are still in TAR format, run:"
    echo "  ./extract_original_files.sh"
    exit 1
fi

# Show what we found
echo "Found files:"
echo "$R1_FILES" | while read -r file; do
    if [ -n "$file" ]; then
        echo "  - $file"
    fi
done
echo

# Process each file
echo "$R1_FILES" | while read -r R1; do
    if [ -z "$R1" ]; then
        continue
    fi
    
    # Get base name and sample info
    BASENAME=$(basename "$R1" .fastq.gz)
    DIRNAME=$(dirname "$R1")
    SAMPLE_NAME=$(echo "$BASENAME" | grep -oE "IL[0-9]+-[0-9]+")
    
    echo
    echo -e "${GREEN}Processing: $BASENAME${NC}"
    echo "Sample: $SAMPLE_NAME"
    
    # Define output file
    OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_R1_properly_trimmed.fastq.gz"
    LOG_FILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}_cutadapt.log"
    
    # Skip if already processed
    if [ -f "$OUTPUT_FILE" ]; then
        echo "  Already processed, skipping..."
        continue
    fi
    
    # Run cutadapt with proper parameters for NextFlex v4
    echo "Running cutadapt..."
    
    cutadapt \
        -a "$ADAPTER_3P" \
        --quality-cutoff 20 \
        --minimum-length 18 \
        --maximum-length 35 \
        --trim-n \
        --max-n 0 \
        -o "$OUTPUT_FILE" \
        "$R1" > "$LOG_FILE" 2>&1
    
    # Check if successful
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Successfully trimmed${NC}"
        
        # Extract key stats from log
        TOTAL_READS=$(grep "Total reads processed:" "$LOG_FILE" | awk '{print $NF}')
        READS_WITH_ADAPTERS=$(grep "Reads with adapters:" "$LOG_FILE" | awk '{print $4}')
        READS_TOO_SHORT=$(grep "Reads that were too short:" "$LOG_FILE" | awk '{print $6}')
        READS_WRITTEN=$(grep "Reads written (passing filters):" "$LOG_FILE" | awk '{print $5}')
        
        echo "  Total reads: $TOTAL_READS"
        echo "  Reads with adapters: $READS_WITH_ADAPTERS"
        echo "  Reads too short: $READS_TOO_SHORT"
        echo "  Reads kept: $READS_WRITTEN"
    else
        echo -e "${RED}✗ Trimming failed! Check log: $LOG_FILE${NC}"
    fi
done

# Quick length check on trimmed files
echo
echo -e "${BLUE}=== Checking trimmed read lengths ===${NC}"

# List what files were actually created
echo "Files created in $OUTPUT_DIR:"
ls -la "$OUTPUT_DIR"/*.fastq.gz 2>/dev/null || echo "No .fastq.gz files found"

# Check each file that exists
for TRIMMED in "$OUTPUT_DIR"/*_properly_trimmed.fastq.gz; do
    if [ -f "$TRIMMED" ]; then
        SAMPLE=$(basename "$TRIMMED" | cut -d'_' -f1)
        echo
        echo "Sample: $SAMPLE"
        echo "File: $(basename "$TRIMMED")"
        echo "Length distribution (top 10):"
        
        # Use gunzip -c instead of zcat for macOS compatibility
        gunzip -c "$TRIMMED" | \
            awk 'NR%4==2 {print length($0)}' | \
            sort | uniq -c | sort -nr | head -10
        
        # Calculate stats
        TOTAL=$(gunzip -c "$TRIMMED" | awk 'NR%4==2' | wc -l)
        if [ "$TOTAL" -gt 0 ]; then
            IN_RANGE=$(gunzip -c "$TRIMMED" | awk 'NR%4==2 {if(length($0)>=18 && length($0)<=30) count++} END {print count+0}')
            PCT=$(awk "BEGIN {printf \"%.1f\", $IN_RANGE * 100 / $TOTAL}")
            
            echo "Total reads: $TOTAL"
            echo "Reads 18-30nt: $IN_RANGE ($PCT%)"
        fi
    fi
done

echo
echo -e "${GREEN}=== Trimming Complete ===${NC}"
echo
echo "Properly trimmed files are in: $OUTPUT_DIR/"
echo "Cutadapt logs are in: $OUTPUT_DIR/logs/"
echo
echo -e "${YELLOW}Next steps:${NC}"
echo "1. Check the length distributions above"
echo "2. If >80% reads are 18-30nt, proceed with alignment"
echo "3. Use these files instead of the pre-trimmed ones"
echo

# Also process R2 if doing paired-end
echo -e "${BLUE}Note: This script only processed R1 files.${NC}"
echo "For paired-end analysis, R2 files would need different handling."