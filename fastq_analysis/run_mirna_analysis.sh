#!/bin/bash

# =============================================================================
# Complete miRNA-seq Analysis Runner
# =============================================================================
# This script runs the complete analysis pipeline step by step
# Just run: ./run_mirna_analysis.sh
# =============================================================================

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Base directory
BASE_DIR="/Users/rvidva/Documents/datasets/XTRIA"
cd "$BASE_DIR"

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
print_step "2" "Installing R Packages"

# Create R package installation script
cat > install_packages.R << 'EOF'
#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Installing required R packages...\n\n")

# Function to install if not present
install_if_needed <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        if (bioc) {
            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            BiocManager::install(pkg, ask = FALSE, update = FALSE)
        } else {
            install.packages(pkg, dependencies = TRUE)
        }
    } else {
        cat(pkg, "already installed\n")
    }
}

# Install packages
install_if_needed("BiocManager")
install_if_needed("tidyverse")
install_if_needed("ggplot2")
install_if_needed("openxlsx")
install_if_needed("RColorBrewer")
install_if_needed("pheatmap")
install_if_needed("knitr")

# Bioconductor packages
install_if_needed("Rsubread", bioc = TRUE)
install_if_needed("edgeR", bioc = TRUE)
install_if_needed("limma", bioc = TRUE)

cat("\n✓ Package installation complete\n")
EOF

echo "Running R package installation..."
Rscript install_packages.R
check_status "R package installation"

# =============================================================================
# STEP 3: Extract Data
# =============================================================================
print_step "3" "Extracting and Organizing Data"

# Always extract to ensure we have the latest files
if [ -f "XTRIA.20250601_163957.ILLUMINA_DATA.1-of-1.tar" ]; then
    echo "Extracting TAR file..."
    tar -xf XTRIA.20250601_163957.ILLUMINA_DATA.1-of-1.tar
    check_status "Data extraction"
    
    echo "Files extracted. Checking structure..."
    
    # Check what was extracted
    if [ -d "Electroporated" ] && [ -d "Notelectroporated" ]; then
        echo "✓ Found sample directories"
        
        # List files in each directory
        echo "Electroporated samples:"
        ls -la Electroporated/*.fastq.gz 2>/dev/null | head -5
        echo "Notelectroporated samples:"
        ls -la Notelectroporated/*.fastq.gz 2>/dev/null | head -5
    else
        # Files might be in main directory
        echo "Checking for FASTQ files in main directory..."
        ls -la *.fastq.gz 2>/dev/null | head -10
    fi
else
    echo -e "${RED}Error: Cannot find TAR archive${NC}"
    exit 1
fi

# Copy mature.fa if it exists in main directory
if [ -f "mature.fa" ] && [ ! -f "complete_analysis/4_mapping/mirbase/mature.fa" ]; then
    echo "Found mature.fa in main directory, copying to analysis folder..."
    mkdir -p complete_analysis/4_mapping/mirbase
    cp mature.fa complete_analysis/4_mapping/mirbase/
fi

# =============================================================================
# STEP 4: Run Preprocessing (FastQC/MultiQC)
# =============================================================================
print_step "4" "Running Quality Control (FastQC/MultiQC)"

# Create preprocessing script
cat > run_preprocessing.sh << 'EOF'
#!/bin/bash

# Create output directories
mkdir -p complete_analysis/{1_fastqc_raw,2_trimmed_reads,3_fastqc_trimmed,logs}

# Run FastQC on raw data
echo "Running FastQC on raw data..."

# Check both directory structures
if [ -d "Electroporated" ] && [ -d "Notelectroporated" ]; then
    # Files are in subdirectories
    for dir in Electroporated Notelectroporated; do
        for file in "$dir"/*R[12].fastq.gz; do
            if [ -f "$file" ] && [[ ! "$file" =~ "trimmed" ]]; then
                fastqc -o complete_analysis/1_fastqc_raw -t 2 "$file" 2>&1 | tee -a complete_analysis/logs/fastqc_raw.log
            fi
        done
    done
else
    # Files are in main directory
    for file in *R[12].fastq.gz; do
        if [ -f "$file" ] && [[ ! "$file" =~ "trimmed" ]]; then
            fastqc -o complete_analysis/1_fastqc_raw -t 2 "$file" 2>&1 | tee -a complete_analysis/logs/fastqc_raw.log
        fi
    done
fi

# Run MultiQC on raw data
cd complete_analysis/1_fastqc_raw
if ls *.zip 1> /dev/null 2>&1; then
    if command -v multiqc >/dev/null 2>&1; then
        multiqc . -n multiqc_raw_report
    else
        python3 -m multiqc . -n multiqc_raw_report
    fi
fi
cd ../..

# Copy trimmed files to analysis directory
echo "Organizing trimmed files..."

# Check both possible locations
if [ -d "Electroporated" ] && [ -d "Notelectroporated" ]; then
    cp Electroporated/*_trimmed.fastq.gz complete_analysis/2_trimmed_reads/ 2>/dev/null || true
    cp Notelectroporated/*_trimmed.fastq.gz complete_analysis/2_trimmed_reads/ 2>/dev/null || true
else
    cp *_trimmed.fastq.gz complete_analysis/2_trimmed_reads/ 2>/dev/null || true
fi

# Run FastQC on trimmed data
echo "Running FastQC on trimmed data..."

# Check both directory structures
if [ -d "Electroporated" ] && [ -d "Notelectroporated" ]; then
    for dir in Electroporated Notelectroporated; do
        for file in "$dir"/*_trimmed.fastq.gz; do
            if [ -f "$file" ]; then
                fastqc -o complete_analysis/3_fastqc_trimmed -t 2 "$file" 2>&1 | tee -a complete_analysis/logs/fastqc_trimmed.log
            fi
        done
    done
else
    for file in *_trimmed.fastq.gz; do
        if [ -f "$file" ]; then
            fastqc -o complete_analysis/3_fastqc_trimmed -t 2 "$file" 2>&1 | tee -a complete_analysis/logs/fastqc_trimmed.log
        fi
    done
fi

# Run MultiQC on trimmed data
cd complete_analysis/3_fastqc_trimmed
if ls *.zip 1> /dev/null 2>&1; then
    if command -v multiqc >/dev/null 2>&1; then
        multiqc . -n multiqc_trimmed_report
    else
        python3 -m multiqc . -n multiqc_trimmed_report
    fi
fi
cd ../..

echo "✓ Preprocessing complete"
EOF

chmod +x run_preprocessing.sh

if command -v fastqc >/dev/null 2>&1; then
    ./run_preprocessing.sh
    check_status "Quality control analysis"
else
    echo -e "${YELLOW}Skipping FastQC/MultiQC (tools not installed)${NC}"
fi

# =============================================================================
# STEP 5: Run Main R Analysis
# =============================================================================
print_step "5" "Running miRNA Analysis Pipeline"

# Create the main R analysis script
cat > run_analysis.R << 'EOF'
#!/usr/bin/env Rscript

# Suppress startup messages
options(warn = -1)
suppressPackageStartupMessages({
    library(tidyverse)
    library(edgeR)
    library(limma)
    library(Rsubread)
    library(ggplot2)
    library(openxlsx)
})
options(warn = 0)

# Set paths
base_dir <- getwd()
output_dir <- file.path(base_dir, "complete_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Progress function
progress <- function(msg) {
    cat(paste0("\n[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n"))
}

progress("Starting miRNA-seq analysis pipeline")

# Define samples
samples <- data.frame(
    SampleID = c("IL21602-001", "IL21603-001"),
    SampleName = c("Electroporated", "Notelectroporated"),
    Condition = c("Electroporated", "Notelectroporated"),
    stringsAsFactors = FALSE
)

# =============================================================================
# Reference Mapping
# =============================================================================
progress("Setting up reference mapping")

mapping_dir <- file.path(output_dir, "4_mapping")
dir.create(mapping_dir, showWarnings = FALSE)
mirbase_dir <- file.path(mapping_dir, "mirbase")
dir.create(mirbase_dir, showWarnings = FALSE)

# Check for mature.fa
mirbase_fa <- file.path(mirbase_dir, "mature.fa")
if (!file.exists(mirbase_fa)) {
    if (file.exists("mature.fa")) {
        progress("Using mature.fa from base directory")
        file.copy("mature.fa", mirbase_fa)
    } else {
        stop("mature.fa not found")
    }
}

# Extract human miRNAs
human_mirna_fa <- file.path(mirbase_dir, "hsa_mature.fa")
if (!file.exists(human_mirna_fa)) {
    system(paste0("grep -A 1 '^>hsa-' '", mirbase_fa, "' | grep -v '^--$' > '", human_mirna_fa, "'"))
}

# Build index
progress("Building miRNA index")
index_dir <- file.path(mapping_dir, "index")
dir.create(index_dir, showWarnings = FALSE)

if (!file.exists(file.path(index_dir, "index.00.b.tab"))) {
    buildindex(basename = file.path(index_dir, "index"),
               reference = human_mirna_fa,
               indexSplit = FALSE)
}

# Find trimmed files
progress("Looking for trimmed FASTQ files")

# Function to find files
find_trimmed_files <- function() {
    # Based on your directory structure:
    # Electroporated/ILLUMINA_DATA/XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz
    # Notelectroporated/ILLUMINA_DATA/XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz
    
    electro_file <- "Electroporated/ILLUMINA_DATA/XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz"
    notelectro_file <- "Notelectroporated/ILLUMINA_DATA/XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz"
    
    if (file.exists(electro_file) && file.exists(notelectro_file)) {
        return(c(electro_file, notelectro_file))
    }
    
    # Try recursive search if exact paths don't work
    files1 <- list.files("Electroporated", pattern = "IL21602-001.*R1.*trimmed.*fastq.gz$", 
                         full.names = TRUE, recursive = TRUE)
    files2 <- list.files("Notelectroporated", pattern = "IL21603-001.*R1.*trimmed.*fastq.gz$", 
                         full.names = TRUE, recursive = TRUE)
    
    if (length(files1) > 0 && length(files2) > 0) {
        return(c(files1[1], files2[1]))
    }
    
    # Last resort - find any R1 trimmed files
    all_files <- list.files(".", pattern = "R1.*trimmed.*fastq.gz$", 
                           full.names = TRUE, recursive = TRUE)
    
    if (length(all_files) >= 2) {
        electro <- grep("IL21602-001", all_files, value = TRUE)[1]
        notelectro <- grep("IL21603-001", all_files, value = TRUE)[1]
        
        if (!is.na(electro) && !is.na(notelectro)) {
            return(c(electro, notelectro))
        }
    }
    
    stop("Could not find trimmed R1 files. Looking for files in Electroporated/ILLUMINA_DATA/ and Notelectroporated/ILLUMINA_DATA/")
}

trimmed_r1_files <- find_trimmed_files()
progress(paste("Found files:", paste(basename(trimmed_r1_files), collapse = ", ")))

# Perform alignment
progress("Aligning reads to miRNA reference")
bam_dir <- file.path(mapping_dir, "bam_files")
dir.create(bam_dir, showWarnings = FALSE)

alignment_stats <- data.frame()

for (i in 1:length(trimmed_r1_files)) {
    sample_name <- samples$SampleID[i]
    progress(paste("Aligning", sample_name))
    
    # Check if file exists and is readable
    if (!file.exists(trimmed_r1_files[i])) {
        stop(paste("File not found:", trimmed_r1_files[i]))
    }
    
    progress(paste("File size:", file.size(trimmed_r1_files[i]), "bytes"))
    
    bam_file <- file.path(bam_dir, paste0(sample_name, ".bam"))
    
    # Try alignment with error handling
    tryCatch({
        align_result <- align(
            index = file.path(index_dir, "index"),
            readfile1 = trimmed_r1_files[i],
            output_file = bam_file,
            nthreads = 2,  # Reduced threads
            unique = TRUE,
            nBestLocations = 1,
            minFragLength = 15,  # Slightly more permissive
            maxFragLength = 35,  # Slightly more permissive
            maxMismatches = 2    # Allow more mismatches
        )
        
        # Check if alignment succeeded
        if (is.null(align_result) || !is.list(align_result)) {
            stop("Alignment failed - no results returned")
        }
        
        # Print alignment summary
        cat("\nAlignment summary for", sample_name, ":\n")
        print(str(align_result))
        
        # Extract stats safely
        total_reads <- ifelse(!is.null(align_result$Total_reads), 
                             align_result$Total_reads, 0)
        mapped_reads <- ifelse(!is.null(align_result$Mapped_reads), 
                              align_result$Mapped_reads, 0)
        
        stats <- data.frame(
            Sample = sample_name,
            Total_reads = total_reads,
            Mapped_reads = mapped_reads,
            Percent_mapped = ifelse(total_reads > 0, 
                                  round(mapped_reads / total_reads * 100, 2), 
                                  0)
        )
        
        alignment_stats <- rbind(alignment_stats, stats)
        
    }, error = function(e) {
        cat("\nError during alignment of", sample_name, ":\n")
        cat(e$message, "\n")
        
        # Try to continue with zero stats
        stats <- data.frame(
            Sample = sample_name,
            Total_reads = 0,
            Mapped_reads = 0,
            Percent_mapped = 0
        )
        alignment_stats <<- rbind(alignment_stats, stats)
    })
}

# Check if we got any successful alignments
if (nrow(alignment_stats) == 0 || all(alignment_stats$Mapped_reads == 0)) {
    cat("\nWARNING: No reads mapped successfully.\n")
    cat("This might indicate:\n")
    cat("1. The reads are not miRNA sequences\n")
    cat("2. The reads are from a different species\n")
    cat("3. Quality issues with the reads\n")
    cat("\nProceeding with analysis anyway...\n")
}

# =============================================================================
# Count Generation
# =============================================================================
progress("Generating miRNA counts")

count_dir <- file.path(output_dir, "5_counts")
dir.create(count_dir, showWarnings = FALSE)

# Create GTF annotation
gtf_file <- file.path(count_dir, "mirna_annotation.gtf")
mirna_seqs <- readLines(human_mirna_fa)
mirna_names <- mirna_seqs[grepl("^>", mirna_seqs)]
mirna_names <- gsub("^>", "", mirna_names)
mirna_names <- sapply(strsplit(mirna_names, " "), "[", 1)

gtf_entries <- data.frame(
    seqname = mirna_names,
    source = "miRBase",
    feature = "exon",
    start = 1,
    end = 25,
    score = ".",
    strand = "+",
    frame = ".",
    attribute = paste0('gene_id "', mirna_names, '"; transcript_id "', mirna_names, '";')
)

write.table(gtf_entries, gtf_file, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# Count reads
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

fc <- featureCounts(
    files = bam_files,
    annot.ext = gtf_file,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    nthreads = 4,
    strandSpecific = 0,
    isPairedEnd = FALSE,
    countMultiMappingReads = FALSE,
    fraction = FALSE
)

counts <- fc$counts
colnames(counts) <- samples$SampleName

# Save raw counts
write.csv(counts, file.path(count_dir, "raw_counts.csv"))

# =============================================================================
# Normalization and Differential Expression
# =============================================================================
progress("Performing normalization and differential expression")

norm_dir <- file.path(output_dir, "6_normalization")
dir.create(norm_dir, showWarnings = FALSE)

# Create DGEList with proper group information
group <- factor(samples$Condition)
dge <- DGEList(counts = counts, group = group, samples = samples)

# Filter lowly expressed
keep <- rowSums(counts >= 10) >= 1
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

# Check if we have any genes left
if (nrow(dge_filtered) == 0) {
    cat("\nWARNING: No miRNAs passed filtering. Reducing threshold...\n")
    keep <- rowSums(counts >= 5) >= 1
    dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
}

# Calculate normalization factors
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")

# Get normalized values
cpm_values <- cpm(dge_filtered, normalized.lib.sizes = TRUE, log = FALSE)
log_cpm <- cpm(dge_filtered, normalized.lib.sizes = TRUE, log = TRUE)

# Save normalized counts
write.csv(cpm_values, file.path(norm_dir, "normalized_cpm.csv"))

# Differential expression
deg_dir <- file.path(output_dir, "7_differential_expression")
dir.create(deg_dir, showWarnings = FALSE)

# Since we only have 2 samples (no replicates), we need to handle this specially
cat("\nNOTE: Only 2 samples without replicates. Results should be interpreted with caution.\n")

# Set a reasonable dispersion value for miRNA-seq without replicates
bcv <- 0.4  # Typical biological coefficient of variation for miRNA-seq
dge_filtered <- estimateCommonDisp(dge_filtered, verbose = FALSE)

# If dispersion estimation failed, set it manually
if (is.na(dge_filtered$common.dispersion)) {
    dge_filtered$common.dispersion <- bcv^2
    cat("Using assumed dispersion:", dge_filtered$common.dispersion, "\n")
}

# Perform exact test
# Make sure we're using the correct group names
group_names <- levels(dge_filtered$samples$group)
cat("Groups found:", paste(group_names, collapse = ", "), "\n")

# Perform the test
if (length(group_names) >= 2) {
    et <- exactTest(dge_filtered, pair = group_names[c(2, 1)])  # Control vs Treatment
} else {
    # If groups aren't set properly, try to set them
    dge_filtered$samples$group <- factor(c("Notelectroporated", "Electroporated"))
    et <- exactTest(dge_filtered, dispersion = bcv^2)
}

top_tags <- topTags(et, n = Inf)

# Process results
results <- as.data.frame(top_tags)
results$miRNA <- rownames(results)
results$FC <- 2^results$logFC
results$Direction <- ifelse(results$logFC > 0, "Up", "Down")
results$Significant <- abs(results$logFC) > 1 & results$FDR < 0.05

# Save results
write.csv(results, file.path(deg_dir, "differential_expression_results.csv"), row.names = FALSE)

# Focus on Let-7-5p
let7_results <- results[grepl("let-7.*-5p", results$miRNA), ]
write.csv(let7_results, file.path(deg_dir, "let7_5p_differential_expression.csv"), row.names = FALSE)

# =============================================================================
# Generate Summary Report
# =============================================================================
progress("Generating summary report")

# Create Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "Sample_Info")
addWorksheet(wb, "Alignment_Stats")
addWorksheet(wb, "Raw_Counts")
addWorksheet(wb, "Normalized_CPM")
addWorksheet(wb, "DE_Results")
addWorksheet(wb, "Let7_5p_Results")

writeData(wb, "Sample_Info", samples)
writeData(wb, "Alignment_Stats", alignment_stats)
writeData(wb, "Raw_Counts", as.data.frame(counts), rowNames = TRUE)
writeData(wb, "Normalized_CPM", as.data.frame(cpm_values), rowNames = TRUE)
writeData(wb, "DE_Results", results)
writeData(wb, "Let7_5p_Results", let7_results)

saveWorkbook(wb, file.path(output_dir, "complete_analysis_results.xlsx"), overwrite = TRUE)

# Calculate key statistics
total_mirnas <- nrow(counts)
filtered_mirnas <- nrow(dge_filtered)
sig_mirnas <- sum(results$Significant)
let7_reduction <- ifelse(nrow(let7_results) > 0, 
                        round((1 - 2^mean(let7_results$logFC)) * 100, 1), 
                        NA)

# Print summary
cat("\n", rep("=", 60), "\n", sep="")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat(rep("=", 60), "\n\n", sep="")
cat("Alignment Statistics:\n")
print(alignment_stats)
cat("\nmiRNA Detection:\n")
cat("- Total miRNAs detected:", total_mirnas, "\n")
cat("- miRNAs after filtering:", filtered_mirnas, "\n")
cat("- Significantly changed:", sig_mirnas, "\n")
cat("\nLet-7-5p Knockdown:\n")
if (!is.na(let7_reduction)) {
    cat("- Let-7-5p variants found:", nrow(let7_results), "\n")
    cat("- Average reduction:", let7_reduction, "%\n")
    if (let7_reduction > 70) {
        cat("- Status: ✓ SUCCESSFUL KNOCKDOWN\n")
    } else if (let7_reduction > 30) {
        cat("- Status: ⚠ PARTIAL KNOCKDOWN\n")
    } else {
        cat("- Status: ✗ LIMITED KNOCKDOWN\n")
    }
}
cat("\nOutput files:\n")
cat("- Complete results: complete_analysis/complete_analysis_results.xlsx\n")
cat("- DE results: complete_analysis/7_differential_expression/\n")
cat(rep("=", 60), "\n", sep="")

progress("Analysis pipeline completed successfully!")
EOF

echo "Running miRNA analysis..."
Rscript run_analysis.R
check_status "miRNA analysis"

# =============================================================================
# STEP 6: Generate Final Report
# =============================================================================
print_step "6" "Generating Final Report"

# Create visualization script
cat > create_plots.R << 'EOF'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

# Load results
results <- read.csv("complete_analysis/7_differential_expression/differential_expression_results.csv")
let7_results <- read.csv("complete_analysis/7_differential_expression/let7_5p_differential_expression.csv")

# Create output directory for plots
plot_dir <- "complete_analysis/plots"
dir.create(plot_dir, showWarnings = FALSE)

# Volcano plot
p1 <- ggplot(results, aes(x = logFC, y = -log10(PValue))) +
    geom_point(aes(color = Significant), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray70")) +
    theme_minimal() +
    labs(title = "Volcano Plot: Electroporated vs Control",
         x = "Log2 Fold Change",
         y = "-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# Highlight Let-7-5p
if (nrow(let7_results) > 0) {
    p1 <- p1 + 
        geom_point(data = let7_results, 
                   aes(x = logFC, y = -log10(PValue)), 
                   color = "blue", size = 3)
}

ggsave(file.path(plot_dir, "volcano_plot.pdf"), p1, width = 8, height = 6)

# Let-7-5p bar plot
if (nrow(let7_results) > 0) {
    # Load count data
    counts <- read.csv("complete_analysis/5_counts/raw_counts.csv", row.names = 1)
    
    # Extract Let-7-5p counts
    let7_counts <- counts[let7_results$miRNA, ]
    let7_data <- data.frame(
        miRNA = rownames(let7_counts),
        Notelectroporated = let7_counts[, "Notelectroporated"],
        Electroporated = let7_counts[, "Electroporated"]
    ) %>%
        pivot_longer(cols = -miRNA, names_to = "Condition", values_to = "Counts")
    
    p2 <- ggplot(let7_data, aes(x = reorder(miRNA, -Counts), y = Counts + 1, fill = Condition)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.8) +
        scale_y_log10() +
        scale_fill_manual(values = c("Notelectroporated" = "#2E86AB", "Electroporated" = "#E63946")) +
        labs(title = "Let-7-5p Expression Levels",
             x = "Let-7-5p variant",
             y = "Read count + 1 (log10 scale)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    ggsave(file.path(plot_dir, "let7_expression.pdf"), p2, width = 10, height = 6)
}

cat("✓ Plots created successfully\n")
EOF

echo "Creating visualization plots..."
Rscript create_plots.R
check_status "Plot generation"

# =============================================================================
# STEP 7: Create HTML Report
# =============================================================================
print_step "7" "Creating HTML Summary Report"

# Check if results exist and create HTML report
if [ -f "complete_analysis/7_differential_expression/let7_5p_differential_expression.csv" ]; then
    # Read Let-7 results
    let7_reduction=$(Rscript -e "
        let7 <- read.csv('complete_analysis/7_differential_expression/let7_5p_differential_expression.csv')
        if(nrow(let7) > 0) {
            cat(round((1 - 2^mean(let7\$logFC)) * 100, 1))
        } else {
            cat('NA')
        }
    " 2>/dev/null)
    
    # Create HTML report
    cat > complete_analysis/analysis_report.html << EOF
<!DOCTYPE html>
<html>
<head>
    <title>miRNA-seq Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2E86AB; }
        h2 { color: #333; margin-top: 30px; }
        .success { color: #28a745; font-weight: bold; }
        .warning { color: #ffc107; font-weight: bold; }
        .error { color: #dc3545; font-weight: bold; }
        .highlight { background-color: #fff3cd; padding: 10px; border-radius: 5px; margin: 10px 0; }
        table { border-collapse: collapse; width: 100%; margin-top: 10px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h1>miRNA-seq Analysis Report</h1>
    <p>Generated on: $(date)</p>
    
    <div class="highlight">
        <h2>Let-7-5p Knockdown Result</h2>
        <p>Overall Let-7-5p reduction: <span class="success">${let7_reduction}%</span></p>
EOF

    if [ "$let7_reduction" != "NA" ]; then
        if (( $(echo "$let7_reduction > 70" | bc -l) )); then
            echo "<p class='success'>✓ SUCCESSFUL KNOCKDOWN - Strong Let-7-5p reduction achieved</p>" >> complete_analysis/analysis_report.html
        elif (( $(echo "$let7_reduction > 30" | bc -l) )); then
            echo "<p class='warning'>⚠ PARTIAL KNOCKDOWN - Moderate Let-7-5p reduction</p>" >> complete_analysis/analysis_report.html
        else
            echo "<p class='error'>✗ LIMITED KNOCKDOWN - Minimal Let-7-5p reduction</p>" >> complete_analysis/analysis_report.html
        fi
    fi
    
    cat >> complete_analysis/analysis_report.html << EOF
    </div>
    
    <h2>Key Files Generated</h2>
    <ul>
        <li><strong>Complete Results:</strong> complete_analysis_results.xlsx</li>
        <li><strong>QC Reports:</strong> multiqc_raw_report.html, multiqc_trimmed_report.html</li>
        <li><strong>Plots:</strong> volcano_plot.pdf, let7_expression.pdf</li>
    </ul>
    
    <h2>Next Steps</h2>
    <ol>
        <li>Review the Excel file for detailed results</li>
        <li>Check QC reports for data quality</li>
        <li>Use pathway analysis files for functional enrichment</li>
        <li>Validate key findings with qRT-PCR</li>
    </ol>
</body>
</html>
EOF
fi

# =============================================================================
# FINAL SUMMARY
# =============================================================================
echo
echo -e "${GREEN}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║                    ANALYSIS COMPLETE!                         ║${NC}"
echo -e "${GREEN}╚═══════════════════════════════════════════════════════════════╝${NC}"
echo
echo -e "${BLUE}Results saved in: ${BASE_DIR}/complete_analysis/${NC}"
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
    open complete_analysis/
fi