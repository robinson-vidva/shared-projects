#!/usr/bin/env Rscript

# =============================================================================
# COMPLETE miRNA-seq WORKFLOW PIPELINE
# =============================================================================
# Following the exact workflow:
# 1. FastQC/MultiQC on raw data
# 2. Adapter trimming
# 3. FastQC/MultiQC on trimmed data
# 4. Reference mapping and counting
# 5. Normalization and DEG analysis
# 6. Pathway analysis
# =============================================================================

# Set up environment
options(repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(Rsubread)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplot2)
  library(knitr)
})

# Set paths
base_dir <- "/Users/rvidva/Documents/datasets/XTRIA"
output_dir <- file.path(base_dir, "complete_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create log file
log_file <- file.path(output_dir, paste0("analysis_log_", Sys.Date(), ".txt"))
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
}

log_message("Starting complete miRNA-seq analysis pipeline")

# =============================================================================
# STEP 1: Data Extraction and Initial Setup
# =============================================================================
log_message("STEP 1: Extracting data")

setwd(base_dir)
tar_file <- "XTRIA.20250601_163957.ILLUMINA_DATA.1-of-1.tar"

if (!file.exists("XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1.fastq.gz")) {
  log_message("Extracting TAR file...")
  system(paste0("tar -xf ", tar_file))
}

# Define samples
samples <- data.frame(
  SampleID = c("IL21602-001", "IL21603-001"),
  SampleName = c("Electroporated", "Notelectroporated"),
  Condition = c("Electroporated", "Notelectroporated"),
  FileName_R1 = c("XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1.fastq.gz",
                  "XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1.fastq.gz"),
  FileName_R2 = c("XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R2.fastq.gz",
                  "XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R2.fastq.gz"),
  stringsAsFactors = FALSE
)

write.csv(samples, file.path(output_dir, "sample_sheet.csv"), row.names = FALSE)

# =============================================================================
# STEP 2: FastQC on Raw Data
# =============================================================================
log_message("STEP 2: Running FastQC on raw data")

fastqc_raw_dir <- file.path(output_dir, "1_fastqc_raw")
dir.create(fastqc_raw_dir, showWarnings = FALSE)

# Run FastQC on raw files
raw_files <- c(samples$FileName_R1, samples$FileName_R2)
for (file in raw_files) {
  if (file.exists(file)) {
    log_message(paste("Running FastQC on", basename(file)))
    cmd <- paste0("fastqc -o '", fastqc_raw_dir, "' -t 2 '", file, "' 2>&1")
    system(cmd)
  }
}

# Run MultiQC
log_message("Running MultiQC on raw FastQC results")
setwd(fastqc_raw_dir)
system("multiqc . -n multiqc_raw_report 2>&1")
setwd(base_dir)

# =============================================================================
# STEP 3: Adapter Trimming
# =============================================================================
log_message("STEP 3: Adapter trimming")

# Note: Your data already includes trimmed files, but here's how to trim if needed
trimmed_dir <- file.path(output_dir, "2_trimmed_reads")
dir.create(trimmed_dir, showWarnings = FALSE)

# Check if trimmed files already exist
trimmed_exists <- all(file.exists(c(
  "XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz",
  "XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R2_trimmed.fastq.gz"
)))

if (trimmed_exists) {
  log_message("Using pre-trimmed files from sequencing facility")
  # Copy or link trimmed files to our directory
  trimmed_files <- list.files(".", pattern = "_trimmed\\.fastq\\.gz$", full.names = TRUE)
  for (file in trimmed_files) {
    file.copy(file, trimmed_dir)
  }
} else {
  log_message("Trimming adapters with cutadapt")
  # Example cutadapt command for NextFlex Small RNA adapters
  # Adjust adapter sequences based on your library prep kit
  
  for (i in 1:nrow(samples)) {
    # Trim R1
    cmd_r1 <- paste0(
      "cutadapt ",
      "-a TGGAATTCTCGGGTGCCAAGG ",  # NextFlex 3' adapter
      "-m 18 -M 30 ",  # min and max length for miRNAs
      "-o '", file.path(trimmed_dir, gsub("\\.fastq\\.gz", "_trimmed.fastq.gz", samples$FileName_R1[i])), "' ",
      "'", samples$FileName_R1[i], "' ",
      "> '", file.path(trimmed_dir, paste0(samples$SampleID[i], "_R1_trim.log")), "' 2>&1"
    )
    system(cmd_r1)
    
    # Trim R2 if doing paired-end analysis
    cmd_r2 <- paste0(
      "cutadapt ",
      "-a GATCGTCGGACTGTAGAACTCTGAAC ",  # NextFlex 5' adapter
      "-m 18 -M 30 ",
      "-o '", file.path(trimmed_dir, gsub("\\.fastq\\.gz", "_trimmed.fastq.gz", samples$FileName_R2[i])), "' ",
      "'", samples$FileName_R2[i], "' ",
      "> '", file.path(trimmed_dir, paste0(samples$SampleID[i], "_R2_trim.log")), "' 2>&1"
    )
    system(cmd_r2)
  }
}

# =============================================================================
# STEP 4: FastQC on Trimmed Data
# =============================================================================
log_message("STEP 4: Running FastQC on trimmed data")

fastqc_trimmed_dir <- file.path(output_dir, "3_fastqc_trimmed")
dir.create(fastqc_trimmed_dir, showWarnings = FALSE)

# Run FastQC on trimmed files
trimmed_files <- list.files(pattern = "_trimmed\\.fastq\\.gz$", full.names = TRUE)
for (file in trimmed_files) {
  if (file.exists(file)) {
    log_message(paste("Running FastQC on", basename(file)))
    cmd <- paste0("fastqc -o '", fastqc_trimmed_dir, "' -t 2 '", file, "' 2>&1")
    system(cmd)
  }
}

# Run MultiQC
log_message("Running MultiQC on trimmed FastQC results")
setwd(fastqc_trimmed_dir)
system("multiqc . -n multiqc_trimmed_report 2>&1")
setwd(base_dir)

# =============================================================================
# STEP 5: Reference Mapping
# =============================================================================
log_message("STEP 5: Reference mapping")

mapping_dir <- file.path(output_dir, "4_mapping")
dir.create(mapping_dir, showWarnings = FALSE)

# Download and prepare miRNA reference
log_message("Downloading miRBase reference")
mirbase_dir <- file.path(mapping_dir, "mirbase")
dir.create(mirbase_dir, showWarnings = FALSE)

mirbase_fa <- file.path(mirbase_dir, "mature.fa")
if (!file.exists(mirbase_fa)) {
  download.file("https://www.mirbase.org/ftp/CURRENT/mature.fa.gz",
                paste0(mirbase_fa, ".gz"))
  system(paste0("gunzip '", mirbase_fa, ".gz'"))
}

# Extract human miRNAs
human_mirna_fa <- file.path(mirbase_dir, "hsa_mature.fa")
system(paste0("grep -A 1 '^>hsa-' '", mirbase_fa, "' | grep -v '^--$' > '", human_mirna_fa, "'"))

# Build Rsubread index
log_message("Building miRNA index")
index_dir <- file.path(mapping_dir, "index")
dir.create(index_dir, showWarnings = FALSE)

if (!file.exists(file.path(index_dir, "index.00.b.tab"))) {
  buildindex(basename = file.path(index_dir, "index"),
             reference = human_mirna_fa,
             indexSplit = FALSE)
}

# Perform alignment
log_message("Aligning reads to miRNA reference")
bam_dir <- file.path(mapping_dir, "bam_files")
dir.create(bam_dir, showWarnings = FALSE)

# Get trimmed R1 files for miRNA analysis
trimmed_r1_files <- c(
  "XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz",
  "XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz"
)

alignment_stats <- data.frame()
for (i in 1:length(trimmed_r1_files)) {
  sample_name <- samples$SampleID[i]
  log_message(paste("Aligning", sample_name))
  
  bam_file <- file.path(bam_dir, paste0(sample_name, ".bam"))
  
  align_result <- align(
    index = file.path(index_dir, "index"),
    readfile1 = trimmed_r1_files[i],
    output_file = bam_file,
    nthreads = 4,
    unique = TRUE,
    nBestLocations = 1,
    minFragLength = 18,
    maxFragLength = 30,
    nMismatch = 1
  )
  
  # Collect alignment statistics
  stats <- data.frame(
    Sample = sample_name,
    Total_reads = align_result$Total_reads,
    Mapped_reads = align_result$Mapped_reads,
    Uniquely_mapped = align_result$Uniquely_mapped_reads,
    Percent_mapped = round(align_result$Mapped_reads / align_result$Total_reads * 100, 2)
  )
  alignment_stats <- rbind(alignment_stats, stats)
}

write.csv(alignment_stats, file.path(mapping_dir, "alignment_statistics.csv"), row.names = FALSE)
log_message("Alignment statistics saved")

# =============================================================================
# STEP 6: Count Generation
# =============================================================================
log_message("STEP 6: Generating miRNA counts")

count_dir <- file.path(output_dir, "5_counts")
dir.create(count_dir, showWarnings = FALSE)

# Create GTF annotation
log_message("Creating miRNA annotation file")
gtf_file <- file.path(count_dir, "mirna_annotation.gtf")

mirna_seqs <- readLines(human_mirna_fa)
mirna_names <- mirna_seqs[grepl("^>", mirna_seqs)]
mirna_names <- gsub("^>", "", mirna_names)
mirna_names <- sapply(strsplit(mirna_names, " "), "[", 1)

# Create GTF entries
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

# Count reads using featureCounts
log_message("Counting reads per miRNA")
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

# Extract count matrix
counts <- fc$counts
colnames(counts) <- samples$SampleName

# Save raw counts
write.csv(counts, file.path(count_dir, "raw_counts.csv"))
write.csv(fc$stat, file.path(count_dir, "count_statistics.csv"))

# =============================================================================
# STEP 7: Normalization and Filtering
# =============================================================================
log_message("STEP 7: Normalization and filtering")

norm_dir <- file.path(output_dir, "6_normalization")
dir.create(norm_dir, showWarnings = FALSE)

# Create DGEList object
dge <- DGEList(counts = counts, samples = samples)

# Filter lowly expressed miRNAs
# Keep miRNAs with at least 10 counts in at least one sample
keep <- rowSums(counts >= 10) >= 1
dge_filtered <- dge[keep, ]

log_message(paste("Filtered from", nrow(dge), "to", nrow(dge_filtered), "miRNAs"))

# Calculate normalization factors (TMM normalization)
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")

# Get normalized expression values
cpm_values <- cpm(dge_filtered, normalized.lib.sizes = TRUE, log = FALSE)
log_cpm <- cpm(dge_filtered, normalized.lib.sizes = TRUE, log = TRUE)

# Save normalized counts
write.csv(cpm_values, file.path(norm_dir, "normalized_cpm.csv"))
write.csv(log_cpm, file.path(norm_dir, "normalized_log2_cpm.csv"))

# Create QC plots
pdf(file.path(norm_dir, "normalization_qc_plots.pdf"), width = 10, height = 8)

# Library size plot
par(mfrow = c(2, 2))
barplot(dge_filtered$samples$lib.size / 1e6, 
        names = dge_filtered$samples$SampleName,
        las = 2, main = "Library Sizes", ylab = "Millions of reads")

# Box plots before and after normalization
boxplot(log2(counts + 1), las = 2, main = "Raw counts (log2)")
boxplot(log_cpm, las = 2, main = "Normalized log2 CPM")

# MDS plot
plotMDS(dge_filtered, main = "MDS Plot")

dev.off()

# =============================================================================
# STEP 8: Differential Expression Analysis
# =============================================================================
log_message("STEP 8: Differential expression analysis")

deg_dir <- file.path(output_dir, "7_differential_expression")
dir.create(deg_dir, showWarnings = FALSE)

# Since we only have 2 samples (no replicates), we'll use exact test
# Note: This is not ideal - biological replicates are recommended

# Estimate dispersion (will use default value due to no replicates)
dge_filtered <- estimateCommonDisp(dge_filtered, verbose = TRUE)
dge_filtered$common.dispersion <- 0.4  # Typical value for miRNA-seq

# Perform exact test
et <- exactTest(dge_filtered, pair = c("Notelectroporated", "Electroporated"))
top_tags <- topTags(et, n = Inf)

# Add additional statistics
results <- as.data.frame(top_tags)
results$miRNA <- rownames(results)
results$FC <- 2^results$logFC
results$Direction <- ifelse(results$logFC > 0, "Up", "Down")

# Identify significant miRNAs (adjust these criteria as needed)
results$Significant <- abs(results$logFC) > 1 & results$FDR < 0.05

# Save results
write.csv(results, file.path(deg_dir, "differential_expression_results.csv"), row.names = FALSE)

# Focus on Let-7-5p
let7_results <- results[grepl("let-7.*-5p", results$miRNA), ]
write.csv(let7_results, file.path(deg_dir, "let7_5p_differential_expression.csv"), row.names = FALSE)

# Create volcano plot
pdf(file.path(deg_dir, "volcano_plot.pdf"), width = 8, height = 6)
plot(results$logFC, -log10(results$PValue),
     xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano Plot: Electroporated vs Control",
     pch = 20, col = ifelse(results$Significant, "red", "gray"))
abline(v = c(-1, 1), lty = 2)
abline(h = -log10(0.05), lty = 2)

# Highlight Let-7-5p
let7_points <- results[grepl("let-7.*-5p", results$miRNA), ]
points(let7_points$logFC, -log10(let7_points$PValue), 
       pch = 20, col = "blue", cex = 1.5)
dev.off()

# MA plot
pdf(file.path(deg_dir, "ma_plot.pdf"), width = 8, height = 6)
plotSmear(et, de.tags = rownames(results)[results$Significant])
# Highlight Let-7-5p
points(let7_results$logCPM, let7_results$logFC, 
       pch = 20, col = "blue", cex = 1.5)
legend("topright", legend = c("Significant", "Let-7-5p"), 
       col = c("red", "blue"), pch = 20)
dev.off()

# =============================================================================
# STEP 9: Pathway Analysis
# =============================================================================
log_message("STEP 9: Pathway analysis preparation")

pathway_dir <- file.path(output_dir, "8_pathway_analysis")
dir.create(pathway_dir, showWarnings = FALSE)

# Prepare data for pathway analysis tools
# Get significantly changed miRNAs
sig_mirnas <- results[results$Significant, ]
sig_mirnas <- sig_mirnas[order(abs(sig_mirnas$logFC), decreasing = TRUE), ]

# Save lists for pathway analysis tools
# For DIANA-miRPath
mirpath_input <- data.frame(
  miRNA = gsub("hsa-", "", sig_mirnas$miRNA),
  logFC = sig_mirnas$logFC,
  FDR = sig_mirnas$FDR
)
write.csv(mirpath_input, file.path(pathway_dir, "mirpath_input.csv"), row.names = FALSE)

# For miRSystem
write.table(gsub("hsa-", "", sig_mirnas$miRNA), 
            file.path(pathway_dir, "mirsystem_input.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Create summary report
log_message("Creating summary report")

report <- paste0(
  "# miRNA-seq Analysis Report\n",
  "Date: ", Sys.Date(), "\n\n",
  "## Summary Statistics\n",
  "- Total miRNAs detected: ", nrow(counts), "\n",
  "- miRNAs after filtering: ", nrow(dge_filtered), "\n",
  "- Significantly changed miRNAs: ", sum(results$Significant), "\n",
  "- Upregulated in Electroporated: ", sum(results$Significant & results$Direction == "Up"), "\n",
  "- Downregulated in Electroporated: ", sum(results$Significant & results$Direction == "Down"), "\n\n",
  "## Let-7-5p Analysis\n",
  "- Total Let-7-5p variants: ", nrow(let7_results), "\n"
)

if (nrow(let7_results) > 0) {
  avg_let7_fc <- mean(let7_results$logFC)
  report <- paste0(report,
    "- Average Let-7-5p log2FC: ", round(avg_let7_fc, 2), "\n",
    "- Let-7-5p knockdown efficiency: ", 
    round((1 - 2^avg_let7_fc) * 100, 1), "%\n\n"
  )
}

report <- paste0(report,
  "## Next Steps for Pathway Analysis\n",
  "1. Upload 'mirpath_input.csv' to DIANA-miRPath (http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=mirpath/index)\n",
  "2. Use 'mirsystem_input.txt' for miRSystem (http://mirsystem.cgm.ntu.edu.tw/)\n",
  "3. Consider using miRNet for network analysis (https://www.mirnet.ca/)\n",
  "4. Validate key findings with qRT-PCR\n"
)

writeLines(report, file.path(output_dir, "analysis_summary_report.txt"))

# =============================================================================
# STEP 10: Generate Final Summary
# =============================================================================
log_message("STEP 10: Generating final summary")

# Create a comprehensive Excel-like summary
library(openxlsx)
wb <- createWorkbook()

# Add sheets
addWorksheet(wb, "Sample_Info")
addWorksheet(wb, "Alignment_Stats")
addWorksheet(wb, "Raw_Counts")
addWorksheet(wb, "Normalized_CPM")
addWorksheet(wb, "DE_Results")
addWorksheet(wb, "Let7_5p_Results")

# Write data to sheets
writeData(wb, "Sample_Info", samples)
writeData(wb, "Alignment_Stats", alignment_stats)
writeData(wb, "Raw_Counts", as.data.frame(counts), rowNames = TRUE)
writeData(wb, "Normalized_CPM", as.data.frame(cpm_values), rowNames = TRUE)
writeData(wb, "DE_Results", results)
writeData(wb, "Let7_5p_Results", let7_results)

# Save workbook
saveWorkbook(wb, file.path(output_dir, "complete_analysis_results.xlsx"), overwrite = TRUE)

log_message("Analysis complete!")
log_message(paste("All results saved in:", output_dir))

# Print final summary
cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Results directory:", output_dir, "\n")
cat("\nKey findings:\n")
if (nrow(let7_results) > 0) {
  cat("- Let-7-5p knockdown efficiency:", 
      round((1 - 2^mean(let7_results$logFC)) * 100, 1), "%\n")
}
cat("- Total significantly changed miRNAs:", sum(results$Significant), "\n")
cat("\nCheck the following files:\n")
cat("- Complete results: complete_analysis_results.xlsx\n")
cat("- Summary report: analysis_summary_report.txt\n")
cat("- QC reports: */multiqc_*_report.html\n")
cat("=======================================\n")