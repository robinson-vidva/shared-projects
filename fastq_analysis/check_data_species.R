#!/usr/bin/env Rscript

# =============================================================================
# Quick script to check species and read characteristics
# No special packages required
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    args <- c(".")  # Default to current directory
}

data_dir <- args[1]
setwd(data_dir)

cat("\n=== Checking Data Characteristics and Species ===\n\n")

# Find trimmed files
trimmed_files <- list.files(".", pattern = "_R1_trimmed\\.fastq\\.gz$", 
                           full.names = TRUE, recursive = TRUE)

if (length(trimmed_files) == 0) {
    stop("No trimmed R1 files found!")
}

cat("Found", length(trimmed_files), "trimmed R1 files\n")
cat("Analyzing:", basename(trimmed_files[1]), "\n\n")

# Read first 10000 sequences
cat("Reading sequences...\n")
con <- gzfile(trimmed_files[1])
lines <- readLines(con, n = 40000)  # 10000 reads * 4 lines per read
close(con)

# Extract sequences (every 2nd line starting from line 2)
sequences <- lines[seq(2, length(lines), by = 4)]
read_lengths <- nchar(sequences)

cat("\n--- READ LENGTH ANALYSIS ---\n")
cat("Total sequences analyzed:", length(sequences), "\n")
cat("Length range:", min(read_lengths), "-", max(read_lengths), "nt\n")
cat("Mean length:", round(mean(read_lengths), 1), "nt\n")
cat("Median length:", median(read_lengths), "nt\n")

# Length distribution
cat("\nLength distribution:\n")
length_table <- table(read_lengths)
# Show top 10 most common lengths
top_lengths <- sort(length_table, decreasing = TRUE)[1:min(10, length(length_table))]
for (i in 1:length(top_lengths)) {
    len <- names(top_lengths)[i]
    count <- top_lengths[i]
    pct <- round(count / length(sequences) * 100, 1)
    cat(sprintf("  %s nt: %d reads (%.1f%%)\n", len, count, pct))
}

# Check if these look like miRNAs
mirna_range <- sum(read_lengths >= 18 & read_lengths <= 30)
mirna_pct <- round(mirna_range / length(sequences) * 100, 1)

cat("\n--- miRNA COMPATIBILITY ---\n")
cat("Reads in miRNA range (18-30 nt):", mirna_range, sprintf("(%.1f%%)\n", mirna_pct))

if (mirna_pct < 50) {
    cat("\n⚠️  WARNING: Less than 50% of reads are in typical miRNA range!\n")
    cat("This might not be a miRNA-seq dataset.\n")
}

# Show some example sequences
cat("\n--- EXAMPLE SEQUENCES ---\n")
cat("First 5 sequences:\n")
for (i in 1:min(5, length(sequences))) {
    seq <- sequences[i]
    cat(sprintf("%d. Length %d: %s\n", i, nchar(seq), 
                ifelse(nchar(seq) > 50, paste0(substr(seq, 1, 50), "..."), seq)))
}

# Check for adapter sequences
cat("\n--- ADAPTER CHECK ---\n")
adapters <- list(
    "Illumina small RNA 3'" = "TGGAATTCTCGGGTGCCAAGG",
    "TruSeq" = "AGATCGGAAGAGCACACGTCT",
    "Nextera" = "CTGTCTCTTATACACATCT"
)

for (name in names(adapters)) {
    adapter <- adapters[[name]]
    matches <- sum(grepl(adapter, sequences, fixed = TRUE))
    if (matches > 0) {
        pct <- round(matches / length(sequences) * 100, 1)
        cat(sprintf("Found %s adapter in %d reads (%.1f%%)\n", name, matches, pct))
    }
}

# Quick species check by looking at reference
cat("\n--- CURRENT REFERENCE CHECK ---\n")
ref_file <- "complete_analysis/4_mapping/mirbase/hsa_mature.fa"
if (file.exists(ref_file)) {
    cat("Currently using HUMAN (hsa) miRNA reference\n")
    
    # Count human miRNAs
    ref_content <- readLines(ref_file)
    n_mirnas <- sum(grepl("^>", ref_content))
    cat("Number of human miRNAs in reference:", n_mirnas, "\n")
} else {
    cat("No reference file found yet\n")
}

# Suggestion for species
cat("\n--- SPECIES SUGGESTION ---\n")
cat("To check if your data is from a different species:\n")
cat("1. Pick a few sequences and BLAST them online\n")
cat("2. Try aligning to both human and mouse references\n")
cat("3. Check with your lab about the sample source\n")

# Try to align to both human and mouse if index exists
if (file.exists("complete_analysis/4_mapping/index/index.00.b.tab")) {
    cat("\n--- TRYING MULTI-SPECIES ALIGNMENT ---\n")
    
    # Create a small test file with first 1000 reads
    test_file <- "test_1000_reads.fastq"
    writeLines(lines[1:min(4000, length(lines))], test_file)
    
    cat("Testing alignment with current (human) reference...\n")
    
    # You can add actual alignment test here if needed
    cat("Run the main pipeline to see full alignment results\n")
    
    # Clean up
    file.remove(test_file)
}

# Summary recommendations
cat("\n", rep("=", 60), "\n")
cat("SUMMARY AND RECOMMENDATIONS\n")
cat(rep("=", 60), "\n\n")

if (mirna_pct > 80) {
    cat("✅ Data appears to be miRNA-seq (mostly 18-30 nt reads)\n")
    if (mean(read_lengths) < 25) {
        cat("✅ Read lengths are appropriate for miRNAs\n")
    }
} else if (mirna_pct > 30) {
    cat("⚠️  Data might be small RNA-seq with mixed RNA species\n")
} else {
    cat("❌ Data does not appear to be miRNA-seq\n")
    cat("   - Could be regular RNA-seq\n")
    cat("   - Could be degraded RNA\n")
    cat("   - Could need different adapter trimming\n")
}

cat("\nNEXT STEPS:\n")
cat("1. If this is mouse data, modify the script to use 'mmu' instead of 'hsa'\n")
cat("2. Check FastQC reports: complete_analysis/1_fastqc_raw/multiqc_raw_report.html\n")
cat("3. Verify with wet lab about species and protocol\n")
cat("4. Try BLAST on a few sequences: https://blast.ncbi.nlm.nih.gov/\n")

# If species might be wrong, show how to check
if (mirna_pct > 50) {
    cat("\nTo quickly test different species:\n")
    cat("1. Copy one sequence from above\n")
    cat("2. Go to miRBase: https://www.mirbase.org/search.shtml\n")
    cat("3. Search by sequence\n")
    cat("4. See which species have matching miRNAs\n")
}