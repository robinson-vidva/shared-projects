#!/usr/bin/env Rscript

# =============================================================================
# Quick Let-7-5p Analysis with Correct Sample Names
# =============================================================================
# IL21602-001: Electroporated
# IL21603-001: Notelectroporated
# =============================================================================

cat("=== Quick Let-7-5p Knockdown Analysis ===\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# Set paths
base_dir <- "/Users/rvidva/Documents/datasets/XTRIA"
output_dir <- file.path(base_dir, "let7_quick_analysis")
dir.create(output_dir, showWarnings = FALSE)

setwd(base_dir)

# Define samples with correct names
samples <- data.frame(
  SampleID = c("IL21602-001", "IL21603-001"),
  SampleName = c("Electroporated", "Notelectroporated"),
  Description = c("Electroporated EVs (expected reduced Let-7-5p)", 
                  "Not electroporated EVs (control)"),
  stringsAsFactors = FALSE
)

cat("Sample Information:\n")
print(samples)
cat("\n")

# Quick check if processed data exists
if (file.exists("complete_analysis/5_counts/raw_counts.csv")) {
  cat("Loading existing count data...\n")
  counts <- read.csv("complete_analysis/5_counts/raw_counts.csv", row.names = 1)
  
  # Rename columns to match our sample names
  colnames(counts) <- samples$SampleName
  
} else if (file.exists("analysis_output/raw_mirna_counts.csv")) {
  cat("Loading count data from previous analysis...\n")
  counts <- read.csv("analysis_output/raw_mirna_counts.csv", row.names = 1)
  
  # Rename columns
  colnames(counts) <- samples$SampleName
  
} else {
  cat("No processed data found. Running quick FASTQ search...\n\n")
  
  # Quick Let-7 motif search
  r1_files <- c(
    "XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz",
    "XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz"
  )
  
  let7_counts <- data.frame()
  
  for (i in 1:length(r1_files)) {
    if (file.exists(r1_files[i])) {
      cat("Searching", samples$SampleName[i], "...\n")
      
      # Count common Let-7 sequences
      cmd_let7a <- paste0("zcat '", r1_files[i], "' | grep -c 'TGAGGTAGTAGGTTGTATAGTT' || echo 0")
      cmd_let7b <- paste0("zcat '", r1_files[i], "' | grep -c 'TGAGGTAGTAGGTTGTGTGGTT' || echo 0")
      cmd_let7c <- paste0("zcat '", r1_files[i], "' | grep -c 'TGAGGTAGTAGGTTGTATGGTT' || echo 0")
      cmd_gaggtag <- paste0("zcat '", r1_files[i], "' | grep -c 'GAGGTAG' || echo 0")
      
      count_let7a <- as.numeric(system(cmd_let7a, intern = TRUE))
      count_let7b <- as.numeric(system(cmd_let7b, intern = TRUE))
      count_let7c <- as.numeric(system(cmd_let7c, intern = TRUE))
      count_gaggtag <- as.numeric(system(cmd_gaggtag, intern = TRUE))
      
      result <- data.frame(
        Sample = samples$SampleName[i],
        let7a_5p = count_let7a,
        let7b_5p = count_let7b,
        let7c_5p = count_let7c,
        GAGGTAG_motif = count_gaggtag,
        Total_let7 = count_let7a + count_let7b + count_let7c
      )
      
      let7_counts <- rbind(let7_counts, result)
    }
  }
  
  if (nrow(let7_counts) > 0) {
    cat("\nLet-7 sequence counts:\n")
    print(let7_counts)
    
    # Calculate reduction
    if (let7_counts$Total_let7[2] > 0) {
      reduction <- (1 - let7_counts$Total_let7[1] / let7_counts$Total_let7[2]) * 100
      cat(sprintf("\nEstimated Let-7 reduction: %.1f%%\n", reduction))
    }
  }
  
  stop("Please run the complete pipeline for accurate quantification")
}

# If we have count data, analyze Let-7-5p
cat("\n=== Analyzing Let-7-5p Expression ===\n")

# Find Let-7-5p miRNAs
let7_5p <- rownames(counts)[grepl("let-7.*-5p", rownames(counts), ignore.case = TRUE)]
cat("Found", length(let7_5p), "Let-7-5p variants\n\n")

if (length(let7_5p) > 0) {
  # Extract Let-7-5p data
  let7_data <- data.frame(
    miRNA = let7_5p,
    Notelectroporated = counts[let7_5p, "Notelectroporated"],
    Electroporated = counts[let7_5p, "Electroporated"],
    stringsAsFactors = FALSE
  )
  
  # Calculate fold change and reduction
  let7_data$Log2_FC <- log2((let7_data$Electroporated + 1) / (let7_data$Notelectroporated + 1))
  let7_data$Percent_Reduction <- pmax(0, (1 - (let7_data$Electroporated + 1) / (let7_data$Notelectroporated + 1)) * 100)
  
  # Sort by expression in control
  let7_data <- let7_data[order(let7_data$Notelectroporated, decreasing = TRUE), ]
  
  # Calculate overall statistics
  total_notelectro <- sum(let7_data$Notelectroporated)
  total_electro <- sum(let7_data$Electroporated)
  overall_reduction <- (1 - total_electro / total_notelectro) * 100
  
  # Print results
  cat("Let-7-5p Expression Summary:\n")
  cat("---------------------------\n")
  cat(sprintf("Total Let-7-5p in Notelectroporated: %d reads\n", total_notelectro))
  cat(sprintf("Total Let-7-5p in Electroporated: %d reads\n", total_electro))
  cat(sprintf("Overall Let-7-5p reduction: %.1f%%\n\n", overall_reduction))
  
  cat("Top Let-7-5p variants:\n")
  print(head(let7_data[, c("miRNA", "Notelectroporated", "Electroporated", "Percent_Reduction")], 10))
  
  # Save results
  write.csv(let7_data, file.path(output_dir, "let7_5p_expression.csv"), row.names = FALSE)
  
  # Create visualization
  cat("\nCreating visualization...\n")
  
  # Prepare data for plotting
  plot_data <- let7_data %>%
    select(miRNA, Notelectroporated, Electroporated) %>%
    pivot_longer(cols = -miRNA, names_to = "Condition", values_to = "Counts") %>%
    mutate(Condition = factor(Condition, 
                             levels = c("Notelectroporated", "Electroporated")))
  
  # Create bar plot
  p1 <- ggplot(plot_data, aes(x = reorder(miRNA, -Counts), y = Counts + 1, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    scale_y_log10() +
    scale_fill_manual(values = c("Notelectroporated" = "#2E86AB", "Electroporated" = "#E63946")) +
    labs(title = "Let-7-5p Expression: Electroporated vs Not Electroporated",
         subtitle = sprintf("Overall Let-7-5p reduction: %.1f%%", overall_reduction),
         x = "Let-7-5p variant",
         y = "Read count + 1 (log10 scale)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "top")
  
  ggsave(file.path(output_dir, "let7_5p_barplot.pdf"), p1, width = 12, height = 8)
  
  # Create scatter plot
  p2 <- ggplot(let7_data, aes(x = Notelectroporated + 1, y = Electroporated + 1)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", size = 1) +
    geom_abline(intercept = 0, slope = 0.1, linetype = "dotted", color = "red", size = 1) +
    geom_point(size = 4, alpha = 0.8, color = "#2E86AB") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Let-7-5p Expression Scatter Plot",
         subtitle = "Each point represents a Let-7-5p variant",
         x = "Not electroporated (counts + 1)",
         y = "Electroporated (counts + 1)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold")) +
    annotate("text", x = 100, y = 10, label = "90% reduction line", 
             angle = 12, color = "red", size = 3)
  
  ggsave(file.path(output_dir, "let7_5p_scatter.pdf"), p2, width = 8, height = 8)
  
  # Summary report
  report <- paste0(
    "# Let-7-5p Knockdown Analysis Report\n\n",
    "## Summary\n",
    "- Date: ", Sys.Date(), "\n",
    "- Electroporated sample: IL21602-001\n",
    "- Not electroporated sample: IL21603-001\n\n",
    "## Results\n",
    "- Total Let-7-5p variants detected: ", length(let7_5p), "\n",
    "- Total Let-7-5p reads in not electroporated: ", format(total_notelectro, big.mark = ","), "\n",
    "- Total Let-7-5p reads in electroporated: ", format(total_electro, big.mark = ","), "\n",
    "- **Overall Let-7-5p reduction: ", sprintf("%.1f%%", overall_reduction), "**\n\n",
    "## Interpretation\n"
  )
  
  if (overall_reduction > 70) {
    report <- paste0(report, 
      "✅ **SUCCESS**: Strong Let-7-5p knockdown achieved. ",
      "The electroporation effectively reduced Let-7-5p expression.\n")
  } else if (overall_reduction > 30) {
    report <- paste0(report,
      "⚠️ **PARTIAL**: Moderate Let-7-5p knockdown achieved. ",
      "Consider optimization for stronger effect.\n")
  } else {
    report <- paste0(report,
      "❌ **LIMITED**: Minimal Let-7-5p knockdown observed. ",
      "Review electroporation protocol.\n")
  }
  
  writeLines(report, file.path(output_dir, "let7_analysis_report.md"))
  
  cat("\n=== Analysis Complete ===\n")
  cat("Results saved in:", output_dir, "\n")
  cat("\nConclusion: ")
  if (overall_reduction > 70) {
    cat("✅ Successful Let-7-5p knockdown!\n")
  } else if (overall_reduction > 30) {
    cat("⚠️ Partial knockdown achieved.\n")
  } else {
    cat("❌ Limited knockdown effect.\n")
  }
  
} else {
  cat("No Let-7-5p variants found in the data.\n")
  cat("Please check the alignment and counting steps.\n")
}