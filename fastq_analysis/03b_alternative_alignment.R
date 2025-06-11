#!/usr/bin/env Rscript

# =============================================================================
# Alternative Alignment Script for miRNA-seq
# Uses different alignment strategies when standard approach fails
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(Rsubread)
    library(ShortRead)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 03b_alternative_alignment.R <data_directory>")
}

data_dir <- args[1]
setwd(data_dir)

cat("\n=== Alternative miRNA Alignment Approach ===\n")

# =============================================================================
# Step 1: Analyze Read Content
# =============================================================================
cat("\nStep 1: Analyzing read content...\n")

# Find trimmed files
trimmed_files <- list.files("complete_analysis/2_trimmed_reads", 
                           pattern = "_R1_trimmed\\.fastq\\.gz$", 
                           full.names = TRUE)

if (length(trimmed_files) == 0) {
    trimmed_files <- list.files(".", pattern = "_R1_trimmed\\.fastq\\.gz$", 
                               full.names = TRUE, recursive = TRUE)
}

# Read a sample of sequences
cat("Reading sample sequences from:", basename(trimmed_files[1]), "\n")

# Use ShortRead to properly read FASTQ
fq_sample <- readFastq(trimmed_files[1], n = 10000)
sequences <- as.character(sread(fq_sample))
qualities <- as.character(quality(quality(fq_sample)))

# Analyze sequences
seq_lengths <- width(fq_sample)
cat("\nSequence length statistics:\n")
print(summary(seq_lengths))

# Check for adapter sequences
common_adapters <- c(
    "TGGAATTCTCGGGTGCCAAGG",  # Illumina small RNA 3' adapter
    "AGATCGGAAGAGCACACGTCT",  # TruSeq adapter
    "CTGTCTCTTATACACATCT"     # Nextera adapter
)

cat("\nChecking for adapter contamination:\n")
for (adapter in common_adapters) {
    matches <- sum(grepl(adapter, sequences, fixed = TRUE))
    if (matches > 0) {
        cat(sprintf("Found %d reads containing adapter: %s\n", matches, adapter))
    }
}

# =============================================================================
# Step 2: Try Alignment with Different Parameters
# =============================================================================
cat("\nStep 2: Trying different alignment strategies...\n")

mapping_dir <- "complete_analysis/4_mapping"
index_dir <- file.path(mapping_dir, "index")
bam_dir <- file.path(mapping_dir, "bam_files_alt")
dir.create(bam_dir, showWarnings = FALSE)

# Strategy 1: Very short seed length
cat("\nStrategy 1: Short seed alignment (for degraded or modified RNA)...\n")

align_strategies <- list(
    list(
        name = "short_seed",
        params = list(
            nsubreads = 20,  # Many short seeds
            TH1 = 1,  # Minimum voting threshold
            TH2 = 1,
            maxMismatches = 5,
            unique = FALSE,
            nBestLocations = 10
        )
    ),
    list(
        name = "local_align",
        params = list(
            type = "dna",  # Use DNA mode
            nsubreads = 14,
            maxMismatches = 3,
            indels = 5,  # Allow indels
            unique = FALSE
        )
    ),
    list(
        name = "seed_and_vote",
        params = list(
            seed = 10,  # Very short seed
            nsubreads = 30,
            TH1 = 2,
            maxMismatches = 4,
            unique = FALSE
        )
    )
)

results <- data.frame()

for (i in 1:length(trimmed_files)) {
    sample_name <- gsub(".*/(.*?)_.*", "\\1", trimmed_files[i])
    cat("\nProcessing sample:", sample_name, "\n")
    
    for (strategy in align_strategies) {
        cat("  Trying strategy:", strategy$name, "\n")
        
        bam_file <- file.path(bam_dir, paste0(sample_name, "_", strategy$name, ".bam"))
        
        # Prepare parameters
        align_params <- c(
            list(
                index = file.path(index_dir, "index"),
                readfile1 = trimmed_files[i],
                output_file = bam_file,
                nthreads = 2
            ),
            strategy$params
        )
        
        # Try alignment
        tryCatch({
            align_result <- do.call(align, align_params)
            
            result_row <- data.frame(
                Sample = sample_name,
                Strategy = strategy$name,
                Total_reads = align_result$Total_reads,
                Mapped_reads = align_result$Mapped_reads,
                Percent_mapped = round(align_result$Mapped_reads / align_result$Total_reads * 100, 2)
            )
            
            results <- rbind(results, result_row)
            
            cat(sprintf("    Mapped: %d / %d (%.2f%%)\n", 
                       align_result$Mapped_reads, 
                       align_result$Total_reads,
                       result_row$Percent_mapped))
            
        }, error = function(e) {
            cat("    Failed:", e$message, "\n")
        })
    }
}

# =============================================================================
# Step 3: Try Mapping to Genome Instead
# =============================================================================
cat("\nStep 3: Checking if reads map to human genome (not just miRNAs)...\n")
cat("This would help determine if reads are human but not miRNAs\n")
cat("(Requires human genome index - skipping for now)\n")

# =============================================================================
# Step 4: Check for Contamination
# =============================================================================
cat("\nStep 4: Checking sequence composition...\n")

# Get GC content
gc_content <- alphabetFrequency(fq_sample)[, c("G", "C")]
gc_percent <- rowSums(gc_content) / width(fq_sample) * 100

cat("\nGC content statistics:\n")
print(summary(gc_percent))

# Check for poly-A tails (RNA degradation)
polyA_reads <- sum(grepl("AAAAAAA", sequences))
cat("\nReads with poly-A sequences:", polyA_reads, "\n")

# Check for rRNA contamination patterns
rRNA_patterns <- c("GGCTGGTCCTC", "CCTGGCTCTA", "ACCGCGGCTGCTGG")
rRNA_matches <- 0
for (pattern in rRNA_patterns) {
    rRNA_matches <- rRNA_matches + sum(grepl(pattern, sequences))
}
cat("Potential rRNA contamination:", rRNA_matches, "reads\n")

# =============================================================================
# Summary and Recommendations
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("ALTERNATIVE ALIGNMENT SUMMARY\n")
cat(rep("=", 70), "\n\n")

if (nrow(results) > 0) {
    cat("Alignment results by strategy:\n")
    print(results)
    
    best_strategy <- results[which.max(results$Percent_mapped), ]
    cat("\nBest performing strategy:", best_strategy$Strategy, 
        "with", best_strategy$Percent_mapped, "% mapped\n")
}

cat("\nRECOMMENDATIONS:\n")

if (all(results$Percent_mapped < 1)) {
    cat("1. Very low mapping suggests reads are not standard miRNAs\n")
    cat("2. Consider these possibilities:\n")
    cat("   - IsomiRs with non-templated additions\n")
    cat("   - Cross-species contamination\n")
    cat("   - Library prep artifacts\n")
    cat("   - Wrong adapter trimming\n")
    cat("\n3. Next steps:\n")
    cat("   - Try miRNA-specific tools: miRDeep2, sRNAbench\n")
    cat("   - Map to full genome to check read origin\n")
    cat("   - BLAST a few sequences to identify source\n")
    cat("   - Re-examine library prep protocol\n")
} else if (any(results$Percent_mapped > 10)) {
    cat("1. Some alignment success with alternative parameters\n")
    cat("2. Use the best performing strategy for full analysis\n")
    cat("3. Results suggest non-canonical miRNA sequences\n")
}

# Save results
write.csv(results, "complete_analysis/4_mapping/alternative_alignment_results.csv", 
          row.names = FALSE)

cat("\nResults saved to: complete_analysis/4_mapping/alternative_alignment_results.csv\n")
cat(rep("=", 70), "\n")