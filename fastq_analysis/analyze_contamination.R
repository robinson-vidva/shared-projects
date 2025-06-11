#!/usr/bin/env Rscript

# =============================================================================
# Contamination Analysis for Small RNA-seq Data
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(Biostrings)
    library(Rsubread)
})

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    args <- c(".")
}

data_dir <- args[1]
setwd(data_dir)

cat("\n=== Small RNA Contamination Analysis ===\n")
cat("This will check for common contaminants in small RNA-seq data\n\n")

# =============================================================================
# Step 1: Analyze unmapped reads
# =============================================================================
cat("Step 1: Extracting and analyzing unmapped reads...\n")

# Read a sample of properly trimmed sequences
fastq_file <- "properly_trimmed/IL21602-001_R1_properly_trimmed.fastq.gz"
if (!file.exists(fastq_file)) {
    stop("Properly trimmed files not found. Run trimming first.")
}

# Read sequences
cat("Reading sequences...\n")
con <- gzfile(fastq_file)
lines <- readLines(con, n = 40000)  # 10000 reads
close(con)

# Extract sequences
headers <- lines[seq(1, length(lines), by = 4)]
sequences <- lines[seq(2, length(lines), by = 4)]
qualities <- lines[seq(4, length(lines), by = 4)]

# Create a data frame
read_data <- data.frame(
    header = headers,
    sequence = sequences,
    length = nchar(sequences),
    stringsAsFactors = FALSE
)

# =============================================================================
# Step 2: Check for common contaminants
# =============================================================================
cat("\nStep 2: Checking for known contaminants...\n")

# Common contaminant sequences
contaminants <- list(
    # Adapter sequences
    "TruSeq_Adapter" = "AGATCGGAAGAGC",
    "Illumina_SmallRNA_3p" = "TGGAATTCTCGGGTGCCAAGG",
    "NextFlex_3p" = "TGGAATTCTCGGGTGCCAAGG",
    
    # rRNA fragments (human)
    "5S_rRNA_1" = "GTCTACGGCCATACCACCCTGAA",
    "5S_rRNA_2" = "CGCCTGGGAACACGC",
    "5.8S_rRNA" = "CGACTCTTAGCGGTGGATCACT",
    "18S_rRNA_1" = "AAACGGCTACCACATCC",
    "18S_rRNA_2" = "CGGGTCATAAGCTTGC",
    "28S_rRNA_1" = "CGATGAAGAACGCAGC",
    "28S_rRNA_2" = "TCGGAGGGAACCAGCTACTA",
    
    # tRNA fragments
    "tRNA_Gly" = "GCATTGGTGGTTCAGTGG",
    "tRNA_Glu" = "TCCCTGGTGGTCTAGTGG",
    "tRNA_Val" = "GTTTCCGTAGTGTAGTGG",
    "tRNA_Lys" = "GCCCGGCTAGCTCAGTCGG",
    
    # snRNA/snoRNA
    "U1_snRNA" = "ATACTTACCTGGCAGGGGAGA",
    "U2_snRNA" = "ATCGCTTCTCGGCCTTTTGGC",
    "U6_snRNA" = "GTGCTCGCTTCGGCAGCACA",
    
    # PolyA/PolyT
    "PolyA" = "AAAAAAAAAA",
    "PolyT" = "TTTTTTTTTT",
    
    # Common primers
    "RT_primer" = "CAAGCAGAAGACGGCATACGAGAT",
    "PCR_primer" = "AATGATACGGCGACCACCGAGAT"
)

# Check for exact and partial matches
contamination_results <- data.frame()

for (cont_name in names(contaminants)) {
    cont_seq <- contaminants[[cont_name]]
    
    # Exact matches
    exact_matches <- sum(grepl(cont_seq, read_data$sequence, fixed = TRUE))
    
    # Partial matches (at least 10bp)
    if (nchar(cont_seq) > 10) {
        partial_seq <- substr(cont_seq, 1, 10)
        partial_matches <- sum(grepl(partial_seq, read_data$sequence, fixed = TRUE))
    } else {
        partial_matches <- exact_matches
    }
    
    if (exact_matches > 0 || partial_matches > 0) {
        contamination_results <- rbind(contamination_results, data.frame(
            Contaminant = cont_name,
            Exact_Matches = exact_matches,
            Partial_Matches = partial_matches,
            Percent_Exact = round(exact_matches / nrow(read_data) * 100, 2),
            Percent_Partial = round(partial_matches / nrow(read_data) * 100, 2)
        ))
    }
}

if (nrow(contamination_results) > 0) {
    cat("\nContaminants found:\n")
    print(contamination_results)
} else {
    cat("\nNo known contaminants detected in this sample.\n")
}

# =============================================================================
# Step 3: Sequence composition analysis
# =============================================================================
cat("\nStep 3: Analyzing sequence composition...\n")

# GC content
gc_content <- sapply(read_data$sequence, function(seq) {
    bases <- strsplit(seq, "")[[1]]
    gc_count <- sum(bases %in% c("G", "C"))
    gc_count / length(bases) * 100
})

read_data$gc_percent <- gc_content

# Categorize by GC content
cat("\nGC content distribution:\n")
cat("Very low GC (<20%):", sum(gc_content < 20), "reads\n")
cat("Low GC (20-40%):", sum(gc_content >= 20 & gc_content < 40), "reads\n")
cat("Normal GC (40-60%):", sum(gc_content >= 40 & gc_content < 60), "reads\n")
cat("High GC (60-80%):", sum(gc_content >= 60 & gc_content < 80), "reads\n")
cat("Very high GC (>80%):", sum(gc_content >= 80), "reads\n")

# Low complexity sequences
low_complexity <- sapply(read_data$sequence, function(seq) {
    bases <- table(strsplit(seq, "")[[1]])
    max(bases) / nchar(seq) > 0.5  # More than 50% same base
})

cat("\nLow complexity sequences:", sum(low_complexity), 
    sprintf("(%.1f%%)\n", sum(low_complexity) / nrow(read_data) * 100))

# =============================================================================
# Step 4: Identify over-represented sequences
# =============================================================================
cat("\nStep 4: Finding over-represented sequences...\n")

# Count unique sequences
seq_counts <- table(read_data$sequence)
seq_counts <- sort(seq_counts, decreasing = TRUE)

# Get top 20 most common sequences
top_seqs <- head(seq_counts, 20)

cat("\nTop 20 most abundant sequences:\n")
top_seq_df <- data.frame(
    Sequence = names(top_seqs),
    Count = as.numeric(top_seqs),
    Percent = round(as.numeric(top_seqs) / nrow(read_data) * 100, 2),
    Length = nchar(names(top_seqs))
)

print(top_seq_df)

# Save for BLAST analysis
write.table(
    data.frame(
        header = paste0(">seq", 1:20, "_count", top_seq_df$Count),
        sequence = top_seq_df$Sequence
    ),
    "top_sequences_for_blast.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\n"
)

cat("\nTop sequences saved to: top_sequences_for_blast.txt\n")
cat("BLAST these sequences at: https://blast.ncbi.nlm.nih.gov/\n")

# =============================================================================
# Step 5: Create filtered FASTQ
# =============================================================================
cat("\nStep 5: Creating contamination-filtered FASTQ files...\n")

# Define sequences to remove
remove_patterns <- c()

# Add highly repetitive sequences
if (nrow(contamination_results) > 0) {
    high_contam <- contamination_results[contamination_results$Percent_Partial > 5, ]
    if (nrow(high_contam) > 0) {
        for (i in 1:nrow(high_contam)) {
            cont_name <- high_contam$Contaminant[i]
            remove_patterns <- c(remove_patterns, contaminants[[cont_name]])
        }
    }
}

# Add low complexity sequences
remove_low_complexity <- TRUE

# Function to filter reads
filter_reads <- function(input_file, output_file, patterns, remove_low_comp = TRUE) {
    con_in <- gzfile(input_file, "r")
    con_out <- gzfile(output_file, "w")
    
    removed_count <- 0
    kept_count <- 0
    
    while (length(lines <- readLines(con_in, n = 4)) == 4) {
        header <- lines[1]
        sequence <- lines[2]
        plus <- lines[3]
        quality <- lines[4]
        
        # Check if should remove
        remove <- FALSE
        
        # Check contamination patterns
        for (pattern in patterns) {
            if (grepl(pattern, sequence, fixed = TRUE)) {
                remove <- TRUE
                break
            }
        }
        
        # Check low complexity
        if (!remove && remove_low_comp) {
            bases <- table(strsplit(sequence, "")[[1]])
            if (max(bases) / nchar(sequence) > 0.5) {
                remove <- TRUE
            }
        }
        
        # Check polyA/polyT
        if (!remove && (grepl("AAAAAAA", sequence) || grepl("TTTTTTT", sequence))) {
            remove <- TRUE
        }
        
        # Write or skip
        if (!remove) {
            writeLines(c(header, sequence, plus, quality), con_out)
            kept_count <- kept_count + 1
        } else {
            removed_count <- removed_count + 1
        }
    }
    
    close(con_in)
    close(con_out)
    
    return(list(kept = kept_count, removed = removed_count))
}

# Create filtered directory
dir.create("contamination_filtered", showWarnings = FALSE)

# Process each file
for (sample in c("IL21602-001", "IL21603-001")) {
    input_file <- paste0("properly_trimmed/", sample, "_R1_properly_trimmed.fastq.gz")
    output_file <- paste0("contamination_filtered/", sample, "_R1_filtered.fastq.gz")
    
    if (file.exists(input_file)) {
        cat("\nFiltering", sample, "...\n")
        result <- filter_reads(input_file, output_file, remove_patterns, remove_low_complexity)
        
        cat("Kept:", result$kept, "reads\n")
        cat("Removed:", result$removed, "reads\n")
        cat("Removal rate:", round(result$removed / (result$kept + result$removed) * 100, 1), "%\n")
    }
}

# =============================================================================
# Step 6: Summary and Recommendations
# =============================================================================
cat("\n", rep("=", 60), "\n")
cat("CONTAMINATION ANALYSIS SUMMARY\n")
cat(rep("=", 60), "\n\n")

cat("Key findings:\n")

if (nrow(contamination_results) > 0) {
    cat("1. Found contamination from:\n")
    print(contamination_results[, c("Contaminant", "Percent_Partial")])
}

cat("\n2. Sequence characteristics:\n")
cat("   - Low complexity sequences:", 
    sprintf("%.1f%%", sum(low_complexity) / nrow(read_data) * 100), "\n")
cat("   - Mean GC content:", round(mean(gc_content), 1), "%\n")

cat("\n3. Next steps:\n")
cat("   a) BLAST the sequences in 'top_sequences_for_blast.txt'\n")
cat("   b) Try aligning filtered reads from 'contamination_filtered/'\n")
cat("   c) Consider using a different RNA database (tRNA, rRNA, etc.)\n")
cat("   d) Try miRNA discovery tools (miRDeep2) for novel miRNAs\n")

cat("\nTo re-run the pipeline with filtered data:\n")
cat("1. Copy filtered files to properly_trimmed/\n")
cat("2. Run the main pipeline again\n")

cat(rep("=", 60), "\n")