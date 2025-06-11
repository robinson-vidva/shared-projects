#!/bin/bash

# =============================================================================
# Extract only original (non-trimmed) files from TAR archives
# =============================================================================

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}=== Extracting Original Files from TAR Archives ===${NC}"
echo

# Find TAR files
TAR_FILES=$(find . -name "*.tar" -type f)

if [ -z "$TAR_FILES" ]; then
    echo -e "${RED}No TAR files found!${NC}"
    exit 1
fi

echo -e "${BLUE}Found TAR files:${NC}"
echo "$TAR_FILES"
echo

# Create directory for original files
ORIG_DIR="original_files"
mkdir -p "$ORIG_DIR"

# Extract only non-trimmed R1 and R2 files
for TAR in $TAR_FILES; do
    echo -e "${BLUE}Processing: $(basename "$TAR")${NC}"
    
    # List contents to see what's inside
    echo "Contents:"
    tar -tf "$TAR" | grep -E "_R[12]\.fastq\.gz$" | grep -v "trimmed" | head -10
    
    # Extract only original R1 and R2 files (not trimmed)
    echo "Extracting original files..."
    
    # For macOS compatibility, extract specific files without wildcards
    tar -tf "$TAR" | grep -E "_R[12]\.fastq\.gz$" | grep -v "trimmed" | while read -r file; do
        echo "  Extracting: $file"
        tar -xf "$TAR" "$file"
    done
    
    echo -e "${GREEN}✓ Done${NC}"
    echo
done

# Move files to a clean directory structure
echo -e "${BLUE}Organizing files...${NC}"

# Find all extracted original files
find . -name "*_R[12].fastq.gz" ! -name "*trimmed*" ! -path "./$ORIG_DIR/*" -type f | while read -r file; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        echo "Moving: $filename to $ORIG_DIR/"
        cp "$file" "$ORIG_DIR/"
    fi
done

echo
echo -e "${GREEN}=== Extraction Complete ===${NC}"
echo
echo "Original files are in: $ORIG_DIR/"
echo
echo -e "${BLUE}Files extracted:${NC}"
if [ -d "$ORIG_DIR" ] && [ "$(ls -A "$ORIG_DIR"/*.fastq.gz 2>/dev/null | wc -l)" -gt 0 ]; then
    ls -lh "$ORIG_DIR"/*.fastq.gz | awk '{print $9, "(" $5 ")"}'
else
    echo "No files found in $ORIG_DIR/"
    echo
    echo "Checking if files exist in their original locations:"
    find . -name "*_R[12].fastq.gz" ! -name "*trimmed*" -type f | while read -r file; do
        size=$(ls -lh "$file" | awk '{print $5}')
        echo "  $file ($size)"
    done
fi

echo
echo -e "${YELLOW}Next step:${NC}"
echo "Run the trimming script again:"
echo "  ~/GitHub/shared-projects/fastq_analysis/trim_nextflex_v4.sh"