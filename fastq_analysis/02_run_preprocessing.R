#!/usr/bin/env Rscript

# =============================================================================
# Preprocessing Script: Quality Control with FastQC and MultiQC
# =============================================================================

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 02_run_preprocessing.R <data_directory>")
}

data_dir <- args[1]
setwd(data_dir)

cat("\n=== Running Preprocessing and Quality Control ===\n")
cat("Working directory:", getwd(), "\n\n")

# Create output directories
output_dirs <- c(
    "complete_analysis/1_fastqc_raw",
    "complete_analysis/2_trimmed_reads", 
    "complete_analysis/3_fastqc_trimmed",
    "complete_analysis/logs"
)

for (dir in output_dirs) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# Function to run FastQC
run_fastqc <- function(files, output_dir, log_file) {
    if (length(files) == 0) {
        cat("No files found for FastQC\n")
        return()
    }
    
    cat("Running FastQC on", length(files), "files...\n")
    
    for (file in files) {
        cat("  Processing:", basename(file), "\n")
        cmd <- sprintf("fastqc -o '%s' -t 2 '%s' >> '%s' 2>&1", 
                      output_dir, file, log_file)
        system(cmd)
    }
}

# Function to run MultiQC
run_multiqc <- function(input_dir, report_name) {
    cat("Running MultiQC in", input_dir, "\n")
    
    # Check if there are any FastQC results
    fastqc_files <- list.files(input_dir, pattern = "_fastqc\\.zip$", full.names = TRUE)
    
    if (length(fastqc_files) == 0) {
        cat("No FastQC results found for MultiQC\n")
        return()
    }
    
    # Try different ways to run MultiQC
    if (system("command -v multiqc", ignore.stdout = TRUE) == 0) {
        cmd <- sprintf("cd '%s' && multiqc . -n %s --quiet", input_dir, report_name)
    } else {
        cmd <- sprintf("cd '%s' && python3 -m multiqc . -n %s --quiet", input_dir, report_name)
    }
    
    system(cmd)
}

# =============================================================================
# Step 1: FastQC on Raw Data
# =============================================================================
cat("\n--- Step 1: FastQC on Raw Data ---\n")

# Find raw FASTQ files (both R1 and R2)
raw_files <- c()

# Check both directory structures
if (dir.exists("Electroporated") && dir.exists("Notelectroporated")) {
    # Files are in subdirectories
    raw_files <- c(
        list.files("Electroporated/ILLUMINA_DATA", pattern = "\\.fastq\\.gz$", 
                   full.names = TRUE, recursive = TRUE),
        list.files("Notelectroporated/ILLUMINA_DATA", pattern = "\\.fastq\\.gz$", 
                   full.names = TRUE, recursive = TRUE)
    )
} else {
    # Files in main directory
    raw_files <- list.files(".", pattern = "\\.fastq\\.gz$", 
                           full.names = TRUE, recursive = FALSE)
}

# Filter out trimmed files
raw_files <- raw_files[!grepl("trimmed", raw_files)]

cat("Found", length(raw_files), "raw FASTQ files\n")

if (length(raw_files) > 0) {
    # Show first few files
    cat("Files:\n")
    for (i in 1:min(6, length(raw_files))) {
        cat("  ", basename(raw_files[i]), "\n")
    }
    if (length(raw_files) > 6) cat("  ... and", length(raw_files) - 6, "more\n")
    
    run_fastqc(raw_files, "complete_analysis/1_fastqc_raw", 
               "complete_analysis/logs/fastqc_raw.log")
    run_multiqc("complete_analysis/1_fastqc_raw", "multiqc_raw_report")
}

# =============================================================================
# Step 2: Organize Trimmed Files
# =============================================================================
cat("\n--- Step 2: Organizing Trimmed Files ---\n")

# Find trimmed files
trimmed_files <- c()

if (dir.exists("Electroporated") && dir.exists("Notelectroporated")) {
    trimmed_files <- c(
        list.files("Electroporated", pattern = "_trimmed\\.fastq\\.gz$", 
                   full.names = TRUE, recursive = TRUE),
        list.files("Notelectroporated", pattern = "_trimmed\\.fastq\\.gz$", 
                   full.names = TRUE, recursive = TRUE)
    )
} else {
    trimmed_files <- list.files(".", pattern = "_trimmed\\.fastq\\.gz$", 
                               full.names = TRUE, recursive = TRUE)
}

cat("Found", length(trimmed_files), "trimmed files\n")

# Copy trimmed files to analysis directory
if (length(trimmed_files) > 0) {
    for (file in trimmed_files) {
        dest <- file.path("complete_analysis/2_trimmed_reads", basename(file))
        if (!file.exists(dest)) {
            file.copy(file, dest)
            cat("  Copied:", basename(file), "\n")
        }
    }
}

# =============================================================================
# Step 3: FastQC on Trimmed Data
# =============================================================================
cat("\n--- Step 3: FastQC on Trimmed Data ---\n")

if (length(trimmed_files) > 0) {
    run_fastqc(trimmed_files, "complete_analysis/3_fastqc_trimmed", 
               "complete_analysis/logs/fastqc_trimmed.log")
    run_multiqc("complete_analysis/3_fastqc_trimmed", "multiqc_trimmed_report")
}

# =============================================================================
# Step 4: Adapter Trimming (if needed)
# =============================================================================
cat("\n--- Step 4: Checking for Adapter Trimming ---\n")

# If no trimmed files exist, create them
if (length(trimmed_files) == 0 && length(raw_files) > 0) {
    cat("No trimmed files found. Running adapter trimming...\n")
    
    # Common Illumina adapters for small RNA-seq
    adapter_3prime <- "TGGAATTCTCGGGTGCCAAGG"
    
    # Process each R1 file
    r1_files <- raw_files[grepl("_R1", raw_files)]
    
    for (r1_file in r1_files) {
        output_file <- sub("\\.fastq\\.gz$", "_trimmed.fastq.gz", basename(r1_file))
        output_path <- file.path("complete_analysis/2_trimmed_reads", output_file)
        
        cat("  Trimming:", basename(r1_file), "\n")
        
        # Cutadapt command for miRNA-seq
        cmd <- sprintf(
            "cutadapt -a %s -m 18 -M 30 --trim-n -o '%s' '%s' > '%s' 2>&1",
            adapter_3prime,
            output_path,
            r1_file,
            file.path("complete_analysis/logs", paste0(basename(r1_file), "_cutadapt.log"))
        )
        
        if (system("command -v cutadapt", ignore.stdout = TRUE) == 0) {
            system(cmd)
        } else {
            # Try with python3 -m cutadapt
            cmd <- sub("cutadapt", "python3 -m cutadapt", cmd)
            system(cmd)
        }
    }
    
    # Re-run FastQC on newly trimmed files
    new_trimmed <- list.files("complete_analysis/2_trimmed_reads", 
                             pattern = "_trimmed\\.fastq\\.gz$", 
                             full.names = TRUE)
    
    if (length(new_trimmed) > 0) {
        run_fastqc(new_trimmed, "complete_analysis/3_fastqc_trimmed", 
                   "complete_analysis/logs/fastqc_trimmed_new.log")
        run_multiqc("complete_analysis/3_fastqc_trimmed", "multiqc_trimmed_report")
    }
}

cat("\n✓ Preprocessing complete\n")

# Summary
cat("\n=== Preprocessing Summary ===\n")
cat("Raw FASTQ files processed:", length(raw_files), "\n")
cat("Trimmed files available:", length(list.files("complete_analysis/2_trimmed_reads", 
                                                pattern = "_trimmed\\.fastq\\.gz$")), "\n")

if (file.exists("complete_analysis/1_fastqc_raw/multiqc_raw_report.html")) {
    cat("✓ Raw data MultiQC report generated\n")
}

if (file.exists("complete_analysis/3_fastqc_trimmed/multiqc_trimmed_report.html")) {
    cat("✓ Trimmed data MultiQC report generated\n")
}

cat("\nCheck the MultiQC reports for detailed quality metrics.\n")