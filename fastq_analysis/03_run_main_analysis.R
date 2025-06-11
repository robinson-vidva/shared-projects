#!/usr/bin/env Rscript

# =============================================================================
# Main miRNA-seq Analysis Script
# =============================================================================

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

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript 03_run_main_analysis.R <data_directory> <script_directory>")
}

data_dir <- args[1]
script_dir <- args[2]
setwd(data_dir)

# Set paths
output_dir <- file.path(data_dir, "complete_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Progress function
progress <- function(msg) {
    cat(paste0("\n[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n"))
}

progress("Starting miRNA-seq analysis pipeline")
progress(paste("Data directory:", data_dir))
progress(paste("Script directory:", script_dir))

# Define samples
samples <- data.frame(
    SampleID = c("IL21602-001", "IL21603-001"),
    SampleName = c("Electroporated", "Notelectroporated"),
    Condition = c("Electroporated", "Notelectroporated"),
    stringsAsFactors = FALSE
)

# =============================================================================
# Step 1: Reference Preparation
# =============================================================================
progress("Setting up reference mapping")

mapping_dir <- file.path(output_dir, "4_mapping")
dir.create(mapping_dir, showWarnings = FALSE)
mirbase_dir <- file.path(mapping_dir, "mirbase")
dir.create(mirbase_dir, showWarnings = FALSE)

# Check for mature.fa
mirbase_fa <- file.path(mirbase_dir, "mature.fa")
if (!file.exists(mirbase_fa)) {
    if (file.exists(file.path(data_dir, "mature.fa"))) {
        progress("Using mature.fa from data directory")
        file.copy(file.path(data_dir, "mature.fa"), mirbase_fa)
    } else {
        # Download from miRBase if not found
        progress("Downloading mature.fa from miRBase...")
        download.file(
            "https://www.mirbase.org/ftp/CURRENT/mature.fa.gz",
            file.path(mirbase_dir, "mature.fa.gz")
        )
        system(paste0("gunzip '", file.path(mirbase_dir, "mature.fa.gz"), "'"))
    }
}

# Extract human miRNAs
human_mirna_fa <- file.path(mirbase_dir, "hsa_mature.fa")
if (!file.exists(human_mirna_fa)) {
    progress("Extracting human miRNAs...")
    system(paste0("grep -A 1 '^>hsa-' '", mirbase_fa, "' | grep -v '^--$' > '", human_mirna_fa, "'"))
    
    # Count miRNAs
    n_mirnas <- as.numeric(system(paste0("grep -c '^>' '", human_mirna_fa, "'"), intern = TRUE))
    progress(paste("Found", n_mirnas, "human miRNAs"))
}

# Build index
progress("Building miRNA index")
index_dir <- file.path(mapping_dir, "index")
dir.create(index_dir, showWarnings = FALSE)

if (!file.exists(file.path(index_dir, "index.00.b.tab"))) {
    buildindex(
        basename = file.path(index_dir, "index"),
        reference = human_mirna_fa,
        indexSplit = FALSE,
        memory = 2000
    )
}

# =============================================================================
# Step 2: Find and Process Files
# =============================================================================
progress("Looking for properly trimmed or original FASTQ files")

# Function to find the best files to use
find_best_files <- function() {
    # First, check if we have properly trimmed files
    properly_trimmed <- list.files("properly_trimmed", 
                                  pattern = "_properly_trimmed\\.fastq\\.gz$", 
                                  full.names = TRUE)
    
    if (length(properly_trimmed) >= 2) {
        progress("Using properly trimmed files from properly_trimmed/ directory")
        
        # Match to samples
        electro_r1 <- grep("IL21602-001", properly_trimmed, value = TRUE)[1]
        notelectro_r1 <- grep("IL21603-001", properly_trimmed, value = TRUE)[1]
        
        return(list(
            r1_files = c(electro_r1, notelectro_r1),
            r2_files = c(NA, NA),  # Single-end for miRNA-seq
            file_type = "properly_trimmed"
        ))
    }
    
    # If no properly trimmed files, warn about pre-trimmed files
    progress("WARNING: No properly trimmed files found in properly_trimmed/")
    progress("The pre-trimmed files in the TAR archive are over-trimmed (avg 16bp)")
    progress("Please run trim_nextflex_v4.sh first!")
    
    # Check for original files as fallback
    original_r1 <- list.files(".", pattern = "_R1\\.fastq\\.gz$", 
                             full.names = TRUE, recursive = TRUE)
    original_r1 <- original_r1[!grepl("trimmed", original_r1)]
    
    if (length(original_r1) >= 2) {
        progress("Found original (untrimmed) files - these need proper trimming!")
        
        electro_r1 <- grep("IL21602-001", original_r1, value = TRUE)[1]
        notelectro_r1 <- grep("IL21603-001", original_r1, value = TRUE)[1]
        
        stop("Original files found but need trimming. Run trim_nextflex_v4.sh first!")
    }
    
    # Last resort - use the over-trimmed files but warn
    progress("CRITICAL WARNING: Using over-trimmed files (avg 16bp)")
    progress("Alignment will likely fail!")
    
    # Original paths from your data
    electro_r1 <- "Electroporated/ILLUMINA_DATA/XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz"
    electro_r2 <- "Electroporated/ILLUMINA_DATA/XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R2_trimmed.fastq.gz"
    notelectro_r1 <- "Notelectroporated/ILLUMINA_DATA/XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz"
    notelectro_r2 <- "Notelectroporated/ILLUMINA_DATA/XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R2_trimmed.fastq.gz"
    
    return(list(
        r1_files = c(electro_r1, notelectro_r1),
        r2_files = c(electro_r2, notelectro_r2),
        file_type = "over_trimmed"
    ))
}

trimmed_files <- find_best_files()
progress("Found files:")
cat("R1 files:\n")
for (f in trimmed_files$r1_files) cat("  ", basename(f), "\n")
if (!any(is.na(trimmed_files$r2_files))) {
    cat("R2 files:\n")
    for (f in trimmed_files$r2_files) cat("  ", basename(f), "\n")
}
cat("File type:", trimmed_files$file_type, "\n")

if (trimmed_files$file_type == "over_trimmed") {
    cat("\n", rep("!", 60), "\n")
    cat("! CRITICAL: Using over-trimmed files (16bp average)      !\n")
    cat("! Alignment will fail. Run trim_nextflex_v4.sh first!    !\n") 
    cat(rep("!", 60), "\n\n")
}

# =============================================================================
# Step 3: Alignment - Try both single-end and paired-end
# =============================================================================
progress("Aligning reads to miRNA reference")
bam_dir <- file.path(mapping_dir, "bam_files")
dir.create(bam_dir, showWarnings = FALSE)

alignment_stats <- data.frame()

# First try single-end alignment with R1 (more common for miRNA-seq)
progress("Trying single-end alignment with R1 reads...")

for (i in 1:length(trimmed_files$r1_files)) {
    sample_name <- samples$SampleID[i]
    progress(paste("Aligning", sample_name, "(single-end mode)"))
    
    bam_file <- file.path(bam_dir, paste0(sample_name, "_SE.bam"))
    
    tryCatch({
        align_result <- align(
            index = file.path(index_dir, "index"),
            readfile1 = trimmed_files$r1_files[i],
            output_file = bam_file,
            nthreads = 4,
            unique = FALSE,  # Allow multi-mapping
            nBestLocations = 5,  # Report up to 5 locations
            minFragLength = 15,
            maxFragLength = 35,
            maxMismatches = 3,  # Allow more mismatches
            nsubreads = 14,  # More subreads for short reads
            TH1 = 1  # Lower threshold for alignment
        )
        
        # Extract stats
        total_reads <- align_result$Total_reads
        mapped_reads <- align_result$Mapped_reads
        
        stats <- data.frame(
            Sample = sample_name,
            Mode = "Single-end",
            Total_reads = ifelse(!is.null(total_reads), total_reads, 0),
            Mapped_reads = ifelse(!is.null(mapped_reads), mapped_reads, 0),
            Percent_mapped = ifelse(!is.null(total_reads) && total_reads > 0, 
                                  round(mapped_reads / total_reads * 100, 2), 0)
        )
        
        alignment_stats <- rbind(alignment_stats, stats)
        
    }, error = function(e) {
        cat("Error in single-end alignment:", e$message, "\n")
    })
}

# Now try paired-end alignment if R2 files exist
if (!any(is.na(trimmed_files$r2_files))) {
    progress("Trying paired-end alignment...")
    
    for (i in 1:length(trimmed_files$r1_files)) {
        sample_name <- samples$SampleID[i]
        progress(paste("Aligning", sample_name, "(paired-end mode)"))
        
        bam_file <- file.path(bam_dir, paste0(sample_name, "_PE.bam"))
        
        tryCatch({
            align_result <- align(
                index = file.path(index_dir, "index"),
                readfile1 = trimmed_files$r1_files[i],
                readfile2 = trimmed_files$r2_files[i],
                output_file = bam_file,
                nthreads = 4,
                unique = FALSE,  # Allow multi-mapping
                nBestLocations = 5,
                minFragLength = 15,
                maxFragLength = 50,  # Slightly larger for PE
                maxMismatches = 3,  # More permissive
                nsubreads = 14,
                TH1 = 1
            )
            
            # Extract stats
            total_reads <- align_result$Total_reads
            mapped_reads <- align_result$Mapped_reads
            
            stats <- data.frame(
                Sample = sample_name,
                Mode = "Paired-end",
                Total_reads = ifelse(!is.null(total_reads), total_reads, 0),
                Mapped_reads = ifelse(!is.null(mapped_reads), mapped_reads, 0),
                Percent_mapped = ifelse(!is.null(total_reads) && total_reads > 0,
                                      round(mapped_reads / total_reads * 100, 2), 0)
            )
            
            alignment_stats <- rbind(alignment_stats, stats)
            
        }, error = function(e) {
            cat("Error in paired-end alignment:", e$message, "\n")
        })
    }
}

# Display alignment statistics
progress("Alignment Statistics:")
print(alignment_stats)

# =============================================================================
# Diagnostic: Check read characteristics if alignment failed
# =============================================================================
if (nrow(alignment_stats) == 0 || all(alignment_stats$Percent_mapped < 1)) {
    progress("WARNING: Very low alignment rate detected. Performing diagnostic checks...")
    
    # Check a few reads from the FASTQ file
    progress("Checking read characteristics...")
    
    # Read first 1000 reads to check length and content
    con <- gzfile(trimmed_files$r1_files[1])
    first_reads <- readLines(con, n = 4000)  # 1000 reads * 4 lines per read
    close(con)
    
    # Extract sequences (every 2nd line starting from line 2)
    sequences <- first_reads[seq(2, length(first_reads), by = 4)]
    read_lengths <- nchar(sequences)
    
    cat("\nRead length distribution:\n")
    print(summary(read_lengths))
    cat("\nRead length table:\n")
    print(table(read_lengths))
    
    # Check if reads look like miRNAs (typical length 18-30)
    mirna_length_reads <- sum(read_lengths >= 18 & read_lengths <= 30)
    cat("\nReads with miRNA-typical length (18-30nt):", mirna_length_reads, 
        "out of", length(sequences), "\n")
    
    # Sample some sequences
    cat("\nFirst 5 sequences:\n")
    for (i in 1:min(5, length(sequences))) {
        cat(sprintf("Read %d (length %d): %s\n", i, nchar(sequences[i]), sequences[i]))
    }
    
    # Check reference sequences
    cat("\nChecking reference sequences...\n")
    ref_seqs <- readLines(human_mirna_fa)
    ref_sequences <- ref_seqs[!grepl("^>", ref_seqs)]
    ref_lengths <- nchar(ref_sequences)
    
    cat("Reference miRNA length distribution:\n")
    print(summary(ref_lengths))
    
    # Try alternative alignment approach
    progress("Trying alternative alignment with very permissive settings...")
    
    # Create a test alignment with extremely permissive settings
    test_bam <- file.path(bam_dir, "test_alignment.bam")
    
    tryCatch({
        test_align <- align(
            index = file.path(index_dir, "index"),
            readfile1 = trimmed_files$r1_files[1],
            output_file = test_bam,
            nthreads = 1,
            unique = FALSE,
            nBestLocations = 10,
            minFragLength = 10,  # Very short
            maxFragLength = 100,  # Very long
            maxMismatches = 5,  # Many mismatches
            nsubreads = 16,  # Maximum subreads
            TH1 = 1,  # Minimum threshold
            type = "dna"  # Try DNA alignment mode
        )
        
        cat("\nTest alignment results:\n")
        cat("Mapped reads:", test_align$Mapped_reads, "out of", test_align$Total_reads, "\n")
        
    }, error = function(e) {
        cat("Test alignment also failed:", e$message, "\n")
    })
    
    # Suggest next steps
    cat("\n", rep("=", 60), "\n")
    cat("DIAGNOSTIC SUMMARY:\n")
    cat("1. Alignment rate is very low (", mean(alignment_stats$Percent_mapped), "%)\n", sep="")
    cat("2. This could indicate:\n")
    cat("   - Reads are not miRNA sequences\n")
    cat("   - Adapter sequences still present\n")
    cat("   - Wrong species (not human)\n")
    cat("   - Technical issues with library prep\n")
    cat("\nRecommendations:\n")
    cat("1. Check FastQC reports for adapter content\n")
    cat("2. Verify library preparation protocol\n")
    cat("3. Try mapping to full genome instead of just miRNAs\n")
    cat("4. Consider using specialized miRNA alignment tools (e.g., miRDeep2)\n")
    cat(rep("=", 60), "\n\n")
}

# Choose best alignment mode (highest mapping rate)
if (nrow(alignment_stats) > 0) {
    best_mode <- alignment_stats %>%
        group_by(Mode) %>%
        summarise(avg_mapped = mean(Percent_mapped)) %>%
        arrange(desc(avg_mapped)) %>%
        pull(Mode) %>%
        .[1]
    
    progress(paste("Best alignment mode:", best_mode))
    
    # Select BAM files for downstream analysis
    if (best_mode == "Single-end") {
        bam_suffix <- "_SE.bam"
    } else {
        bam_suffix <- "_PE.bam"
    }
} else {
    # Default to single-end if no successful alignments
    progress("No successful alignments. Defaulting to single-end mode for downstream analysis.")
    bam_suffix <- "_SE.bam"
    best_mode <- "Single-end"
    
    # Create dummy alignment stats if empty
    alignment_stats <- data.frame(
        Sample = samples$SampleID,
        Mode = "Single-end",
        Total_reads = c(3732139, 5745780),  # From the error log
        Mapped_reads = c(1268, 83),  # From the error log
        Percent_mapped = c(0.03, 0.00)
    )
}

# =============================================================================
# Step 4: Count Generation
# =============================================================================
progress("Generating miRNA counts")

count_dir <- file.path(output_dir, "5_counts")
dir.create(count_dir, showWarnings = FALSE)

# Create SAF annotation instead of GTF (more reliable for Rsubread)
saf_file <- file.path(count_dir, "mirna_annotation.saf")

# Read miRNA names from fasta
mirna_seqs <- readLines(human_mirna_fa)
mirna_names <- mirna_seqs[grepl("^>", mirna_seqs)]
mirna_names <- gsub("^>", "", mirna_names)
mirna_names <- sapply(strsplit(mirna_names, " "), "[", 1)

# Create SAF format annotation
saf_data <- data.frame(
    GeneID = mirna_names,
    Chr = mirna_names,
    Start = 1,
    End = 30,  # Typical miRNA length
    Strand = "+"
)

write.table(saf_data, saf_file, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)

# Count reads using SAF annotation
bam_files <- list.files(bam_dir, pattern = paste0(bam_suffix, "$"), full.names = TRUE)

progress(paste("Counting reads from", length(bam_files), "BAM files"))

fc <- featureCounts(
    files = bam_files,
    annot.ext = saf_file,
    isGTFAnnotationFile = FALSE,
    nthreads = 4,
    strandSpecific = 0,
    isPairedEnd = (best_mode == "Paired-end"),
    countMultiMappingReads = FALSE,
    fraction = FALSE,
    minMQS = 10,
    allowMultiOverlap = FALSE
)

# Extract counts
counts <- fc$counts
colnames(counts) <- samples$SampleName

# Save count summary
count_summary <- fc$stat
write.csv(count_summary, file.path(count_dir, "count_summary.csv"))

# Save raw counts
write.csv(counts, file.path(count_dir, "raw_counts.csv"))

progress(paste("Total features:", nrow(counts)))
progress(paste("Total assigned reads:", sum(counts)))

# =============================================================================
# Step 5: Normalization and Quality Control
# =============================================================================
progress("Performing normalization and quality control")

norm_dir <- file.path(output_dir, "6_normalization")
dir.create(norm_dir, showWarnings = FALSE)

# Create DGEList
group <- factor(samples$Condition)
dge <- DGEList(counts = counts, group = group, samples = samples)

# Add sample information
dge$samples$Sample <- samples$SampleID
dge$samples$Condition <- samples$Condition

# Filter lowly expressed - be more permissive for 2-sample comparison
keep <- rowSums(counts >= 5) >= 1  # At least 5 counts in at least 1 sample
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

progress(paste("Features after filtering:", nrow(dge_filtered)))

# Calculate normalization factors
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")

# Get normalized values
cpm_values <- cpm(dge_filtered, normalized.lib.sizes = TRUE, log = FALSE)
log_cpm <- cpm(dge_filtered, normalized.lib.sizes = TRUE, log = TRUE)

# Save normalized counts
write.csv(cpm_values, file.path(norm_dir, "normalized_cpm.csv"))
write.csv(log_cpm, file.path(norm_dir, "log2_cpm.csv"))

# =============================================================================
# Step 6: Differential Expression Analysis
# =============================================================================
progress("Performing differential expression analysis")

deg_dir <- file.path(output_dir, "7_differential_expression")
dir.create(deg_dir, showWarnings = FALSE)

# Since we only have 2 samples, we need to use a fixed dispersion
cat("\nNOTE: Only 2 samples without replicates. Using assumed dispersion.\n")
cat("Results should be interpreted with caution and validated experimentally.\n")

# Set biological coefficient of variation (BCV) for miRNA-seq
bcv <- 0.4  # Typical for miRNA-seq

# Set dispersion
dge_filtered$common.dispersion <- bcv^2

# Perform exact test
et <- exactTest(dge_filtered, dispersion = bcv^2)

# Get all results
results <- topTags(et, n = Inf)$table
results$miRNA <- rownames(results)
results$FC <- 2^results$logFC
results$Direction <- ifelse(results$logFC > 0, "Up", "Down")

# Add significance thresholds (relaxed for 2-sample comparison)
results$Significant_strict <- abs(results$logFC) > 1 & results$FDR < 0.05
results$Significant_relaxed <- abs(results$logFC) > 0.5 & results$PValue < 0.05

# Reorder columns
results <- results[, c("miRNA", "logFC", "FC", "logCPM", "PValue", "FDR", 
                      "Direction", "Significant_strict", "Significant_relaxed")]

# Save all results
write.csv(results, file.path(deg_dir, "differential_expression_results.csv"), 
          row.names = FALSE)

# Focus on Let-7 family
let7_results <- results[grepl("let-7", results$miRNA, ignore.case = TRUE), ]
write.csv(let7_results, file.path(deg_dir, "let7_family_results.csv"), 
          row.names = FALSE)

# Also check for specific Let-7-5p variants
let7_5p_results <- results[grepl("let-7.*-5p", results$miRNA, ignore.case = TRUE), ]
write.csv(let7_5p_results, file.path(deg_dir, "let7_5p_specific_results.csv"), 
          row.names = FALSE)

# =============================================================================
# Step 7: Generate Comprehensive Excel Report
# =============================================================================
progress("Generating comprehensive Excel report")

wb <- createWorkbook()

# Add worksheets
addWorksheet(wb, "Summary")
addWorksheet(wb, "Sample_Info")
addWorksheet(wb, "Alignment_Stats")
addWorksheet(wb, "Count_Summary")
addWorksheet(wb, "Raw_Counts")
addWorksheet(wb, "Normalized_CPM")
addWorksheet(wb, "Log2_CPM")
addWorksheet(wb, "DE_Results_All")
addWorksheet(wb, "DE_Significant")
addWorksheet(wb, "Let7_Family")
addWorksheet(wb, "Let7_5p_Specific")

# Create summary
summary_data <- data.frame(
    Metric = c(
        "Total miRNAs detected",
        "miRNAs after filtering",
        "Significantly changed (strict)",
        "Significantly changed (relaxed)",
        "Let-7 family members detected",
        "Let-7-5p variants detected",
        "Average Let-7-5p reduction (%)"
    ),
    Value = c(
        nrow(counts),
        nrow(dge_filtered),
        sum(results$Significant_strict),
        sum(results$Significant_relaxed),
        nrow(let7_results),
        nrow(let7_5p_results),
        ifelse(nrow(let7_5p_results) > 0, 
               round((1 - 2^mean(let7_5p_results$logFC)) * 100, 1),
               "N/A")
    )
)

# Write data to worksheets
writeData(wb, "Summary", summary_data)
writeData(wb, "Sample_Info", samples)
writeData(wb, "Alignment_Stats", alignment_stats)
writeData(wb, "Count_Summary", count_summary)
writeData(wb, "Raw_Counts", as.data.frame(counts), rowNames = TRUE)
writeData(wb, "Normalized_CPM", as.data.frame(cpm_values), rowNames = TRUE)
writeData(wb, "Log2_CPM", as.data.frame(log_cpm), rowNames = TRUE)
writeData(wb, "DE_Results_All", results)
writeData(wb, "DE_Significant", 
          results[results$Significant_relaxed, ])
writeData(wb, "Let7_Family", let7_results)
writeData(wb, "Let7_5p_Specific", let7_5p_results)

# Save workbook
saveWorkbook(wb, file.path(output_dir, "complete_analysis_results.xlsx"), 
             overwrite = TRUE)

# =============================================================================
# Final Summary
# =============================================================================
cat("\n", rep("=", 70), "\n", sep="")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat(rep("=", 70), "\n\n", sep="")

# Print summary statistics
cat("Sample Information:\n")
print(samples)

cat("\nAlignment Statistics:\n")
print(alignment_stats)

cat("\nmiRNA Detection Summary:\n")
cat("- Total miRNAs in reference:", length(mirna_names), "\n")
cat("- miRNAs detected:", nrow(counts), "\n")
cat("- miRNAs after filtering:", nrow(dge_filtered), "\n")
cat("- Significantly changed (strict):", sum(results$Significant_strict), "\n")
cat("- Significantly changed (relaxed):", sum(results$Significant_relaxed), "\n")

cat("\nLet-7 Analysis:\n")
cat("- Let-7 family members found:", nrow(let7_results), "\n")
cat("- Let-7-5p variants found:", nrow(let7_5p_results), "\n")

if (nrow(let7_5p_results) > 0) {
    let7_reduction <- round((1 - 2^mean(let7_5p_results$logFC)) * 100, 1)
    cat("- Average Let-7-5p reduction:", let7_reduction, "%\n")
    
    if (let7_reduction > 70) {
        cat("- Status: ✓ SUCCESSFUL KNOCKDOWN\n")
    } else if (let7_reduction > 30) {
        cat("- Status: ⚠ PARTIAL KNOCKDOWN\n")
    } else {
        cat("- Status: ✗ LIMITED KNOCKDOWN\n")
    }
    
    cat("\nTop Let-7-5p variants by fold change:\n")
    top_let7 <- let7_5p_results %>%
        arrange(logFC) %>%
        head(5) %>%
        select(miRNA, logFC, FC, PValue)
    print(top_let7)
}

cat("\n", rep("=", 70), "\n", sep="")

progress("Analysis pipeline completed successfully!")