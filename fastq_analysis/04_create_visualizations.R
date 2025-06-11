#!/usr/bin/env Rscript

# =============================================================================
# Visualization Script for miRNA-seq Analysis
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(RColorBrewer)
    library(pheatmap)
    library(scales)
    library(gridExtra)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 04_create_visualizations.R <data_directory>")
}

data_dir <- args[1]
setwd(data_dir)

cat("\n=== Creating Visualizations ===\n")

# Create output directory
plot_dir <- file.path("complete_analysis", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading analysis results...\n")

# Load DE results
de_results <- read.csv("complete_analysis/7_differential_expression/differential_expression_results.csv")
let7_results <- read.csv("complete_analysis/7_differential_expression/let7_family_results.csv")
let7_5p_results <- read.csv("complete_analysis/7_differential_expression/let7_5p_specific_results.csv")

# Load count data
raw_counts <- read.csv("complete_analysis/5_counts/raw_counts.csv", row.names = 1)
cpm_data <- read.csv("complete_analysis/6_normalization/normalized_cpm.csv", row.names = 1)
log_cpm <- read.csv("complete_analysis/6_normalization/log2_cpm.csv", row.names = 1)

# Set theme
theme_set(theme_minimal(base_size = 12))

# Color palette
colors <- c("Notelectroporated" = "#2E86AB", "Electroporated" = "#E63946")

# =============================================================================
# 1. Volcano Plot
# =============================================================================
cat("Creating volcano plot...\n")

# Add labels for top genes
de_results$label <- ifelse(
    de_results$Significant_strict | grepl("let-7", de_results$miRNA, ignore.case = TRUE),
    de_results$miRNA,
    ""
)

p_volcano <- ggplot(de_results, aes(x = logFC, y = -log10(PValue))) +
    geom_point(aes(color = Significant_relaxed), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray70"),
                      name = "Significant") +
    labs(title = "Volcano Plot: Electroporated vs Control",
         subtitle = "Let-7-5p knockdown experiment",
         x = "Log2 Fold Change",
         y = "-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

# Highlight Let-7 family
if (nrow(let7_results) > 0) {
    p_volcano <- p_volcano + 
        geom_point(data = let7_results, 
                   aes(x = logFC, y = -log10(PValue)), 
                   color = "blue", size = 3, shape = 17) +
        geom_text(aes(label = label), size = 3, hjust = -0.1, vjust = 0.5,
                  check_overlap = TRUE)
}

ggsave(file.path(plot_dir, "volcano_plot.pdf"), p_volcano, width = 10, height = 8)
ggsave(file.path(plot_dir, "volcano_plot.png"), p_volcano, width = 10, height = 8, dpi = 300)

# =============================================================================
# 2. MA Plot
# =============================================================================
cat("Creating MA plot...\n")

p_ma <- ggplot(de_results, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = Significant_relaxed), alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray70"),
                      name = "Significant") +
    labs(title = "MA Plot: Electroporated vs Control",
         x = "Average log CPM",
         y = "Log2 Fold Change") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# Highlight Let-7 family
if (nrow(let7_results) > 0) {
    p_ma <- p_ma + 
        geom_point(data = let7_results, 
                   aes(x = logCPM, y = logFC), 
                   color = "blue", size = 3, shape = 17)
}

ggsave(file.path(plot_dir, "ma_plot.pdf"), p_ma, width = 10, height = 8)
ggsave(file.path(plot_dir, "ma_plot.png"), p_ma, width = 10, height = 8, dpi = 300)

# =============================================================================
# 3. Let-7 Family Expression Plot
# =============================================================================
if (nrow(let7_results) > 0) {
    cat("Creating Let-7 family expression plot...\n")
    
    # Get Let-7 counts
    let7_counts <- raw_counts[let7_results$miRNA, , drop = FALSE]
    let7_cpm <- cpm_data[let7_results$miRNA, , drop = FALSE]
    
    # Prepare data for plotting
    let7_plot_data <- let7_cpm %>%
        rownames_to_column("miRNA") %>%
        pivot_longer(cols = -miRNA, names_to = "Condition", values_to = "CPM") %>%
        left_join(let7_results[, c("miRNA", "logFC", "PValue")], by = "miRNA")
    
    # Order by fold change
    mirna_order <- let7_results %>%
        arrange(logFC) %>%
        pull(miRNA)
    
    let7_plot_data$miRNA <- factor(let7_plot_data$miRNA, levels = mirna_order)
    
    p_let7 <- ggplot(let7_plot_data, aes(x = miRNA, y = CPM + 1, fill = Condition)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.8) +
        scale_y_log10(labels = comma) +
        scale_fill_manual(values = colors) +
        labs(title = "Let-7 Family Expression Levels",
             subtitle = "Normalized counts (CPM + 1)",
             x = "Let-7 variant",
             y = "CPM + 1 (log10 scale)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))
    
    ggsave(file.path(plot_dir, "let7_expression.pdf"), p_let7, 
           width = max(10, 0.5 * nrow(let7_results)), height = 8)
    ggsave(file.path(plot_dir, "let7_expression.png"), p_let7, 
           width = max(10, 0.5 * nrow(let7_results)), height = 8, dpi = 300)
    
    # Create focused Let-7-5p plot if available
    if (nrow(let7_5p_results) > 0) {
        cat("Creating Let-7-5p specific plot...\n")
        
        let7_5p_cpm <- cpm_data[let7_5p_results$miRNA, , drop = FALSE]
        
        let7_5p_plot_data <- let7_5p_cpm %>%
            rownames_to_column("miRNA") %>%
            pivot_longer(cols = -miRNA, names_to = "Condition", values_to = "CPM") %>%
            left_join(let7_5p_results[, c("miRNA", "logFC", "PValue")], by = "miRNA")
        
        mirna_order_5p <- let7_5p_results %>%
            arrange(logFC) %>%
            pull(miRNA)
        
        let7_5p_plot_data$miRNA <- factor(let7_5p_plot_data$miRNA, levels = mirna_order_5p)
        
        p_let7_5p <- ggplot(let7_5p_plot_data, aes(x = miRNA, y = CPM + 1, fill = Condition)) +
            geom_bar(stat = "identity", position = "dodge", width = 0.8) +
            scale_y_log10(labels = comma) +
            scale_fill_manual(values = colors) +
            labs(title = "Let-7-5p Variants: Knockdown Efficiency",
                 subtitle = "Target miRNAs for knockdown",
                 x = "Let-7-5p variant",
                 y = "CPM + 1 (log10 scale)") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5))
        
        ggsave(file.path(plot_dir, "let7_5p_knockdown.pdf"), p_let7_5p, 
               width = max(8, 0.5 * nrow(let7_5p_results)), height = 8)
        ggsave(file.path(plot_dir, "let7_5p_knockdown.png"), p_let7_5p, 
               width = max(8, 0.5 * nrow(let7_5p_results)), height = 8, dpi = 300)
    }
}

# =============================================================================
# 4. Top Differentially Expressed miRNAs
# =============================================================================
cat("Creating top DE miRNAs plot...\n")

# Get top 20 up and down regulated
top_up <- de_results %>%
    filter(logFC > 0) %>%
    arrange(desc(abs(logFC))) %>%
    head(20)

top_down <- de_results %>%
    filter(logFC < 0) %>%
    arrange(desc(abs(logFC))) %>%
    head(20)

top_de <- rbind(top_up, top_down)

# Get expression data for top DE miRNAs
top_de_cpm <- cpm_data[top_de$miRNA, , drop = FALSE]

top_de_plot_data <- top_de_cpm %>%
    rownames_to_column("miRNA") %>%
    pivot_longer(cols = -miRNA, names_to = "Condition", values_to = "CPM") %>%
    left_join(top_de[, c("miRNA", "logFC")], by = "miRNA")

# Order by fold change
top_de_plot_data$miRNA <- factor(top_de_plot_data$miRNA, 
                                 levels = top_de %>% arrange(logFC) %>% pull(miRNA))

p_top_de <- ggplot(top_de_plot_data, aes(x = miRNA, y = CPM + 1, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    scale_y_log10(labels = comma) +
    scale_fill_manual(values = colors) +
    labs(title = "Top Differentially Expressed miRNAs",
         subtitle = "Top 20 up- and down-regulated miRNAs",
         x = "miRNA",
         y = "CPM + 1 (log10 scale)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)) +
    coord_flip()

ggsave(file.path(plot_dir, "top_de_mirnas.pdf"), p_top_de, width = 10, height = 12)
ggsave(file.path(plot_dir, "top_de_mirnas.png"), p_top_de, width = 10, height = 12, dpi = 300)

# =============================================================================
# 5. Heatmap of Top Variable miRNAs
# =============================================================================
cat("Creating heatmap...\n")

# Calculate variance across samples
log_cpm_matrix <- as.matrix(log_cpm)
row_vars <- apply(log_cpm_matrix, 1, var)

# Get top 50 most variable miRNAs
top_var_mirnas <- names(sort(row_vars, decreasing = TRUE))[1:min(50, length(row_vars))]
heatmap_data <- log_cpm_matrix[top_var_mirnas, ]

# Scale rows for heatmap
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Create annotation for Let-7 family
row_annotation <- data.frame(
    Let7_Family = ifelse(grepl("let-7", rownames(heatmap_data), ignore.case = TRUE), 
                        "Let-7", "Other"),
    row.names = rownames(heatmap_data)
)

# Define colors
ann_colors <- list(
    Let7_Family = c("Let-7" = "blue", "Other" = "gray90")
)

# Create heatmap
pdf(file.path(plot_dir, "heatmap_top_variable.pdf"), width = 8, height = 12)
pheatmap(heatmap_data_scaled,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_row = row_annotation,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 50 Most Variable miRNAs\n(scaled log2 CPM)",
         fontsize_row = 8)
dev.off()

# Also save as PNG
png(file.path(plot_dir, "heatmap_top_variable.png"), width = 800, height = 1200, res = 100)
pheatmap(heatmap_data_scaled,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_row = row_annotation,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 50 Most Variable miRNAs\n(scaled log2 CPM)",
         fontsize_row = 8)
dev.off()

# =============================================================================
# 6. Sample Correlation Plot
# =============================================================================
cat("Creating sample correlation plot...\n")

# Calculate correlation
sample_cor <- cor(log_cpm_matrix, method = "pearson")

# Create correlation plot
pdf(file.path(plot_dir, "sample_correlation.pdf"), width = 6, height = 5)
pheatmap(sample_cor,
         display_numbers = TRUE,
         number_format = "%.3f",
         color = colorRampPalette(c("white", "steelblue"))(100),
         main = "Sample Correlation Matrix\n(Pearson correlation of log2 CPM)",
         fontsize = 12)
dev.off()

# =============================================================================
# 7. Summary Statistics Plot
# =============================================================================
cat("Creating summary statistics plot...\n")

# Prepare summary data
summary_stats <- data.frame(
    Category = c("Total\nDetected", "After\nFiltering", "Significant\n(strict)", 
                 "Significant\n(relaxed)", "Let-7\nFamily", "Let-7-5p\nVariants"),
    Count = c(
        nrow(raw_counts),
        nrow(de_results),
        sum(de_results$Significant_strict),
        sum(de_results$Significant_relaxed),
        nrow(let7_results),
        nrow(let7_5p_results)
    )
)

p_summary <- ggplot(summary_stats, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_text(aes(label = Count), vjust = -0.5, size = 5) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "miRNA Detection and Analysis Summary",
         x = "",
         y = "Number of miRNAs") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text.x = element_text(size = 10))

ggsave(file.path(plot_dir, "summary_statistics.pdf"), p_summary, width = 10, height = 8)
ggsave(file.path(plot_dir, "summary_statistics.png"), p_summary, width = 10, height = 8, dpi = 300)

# =============================================================================
# Summary
# =============================================================================
cat("\n✓ Visualizations created successfully!\n")
cat("\nPlots saved in:", file.path(data_dir, plot_dir), "\n")
cat("\nGenerated plots:\n")
cat("- volcano_plot.pdf/png: Volcano plot showing differential expression\n")
cat("- ma_plot.pdf/png: MA plot of fold changes vs expression levels\n")
if (nrow(let7_results) > 0) {
    cat("- let7_expression.pdf/png: Let-7 family expression levels\n")
}
if (nrow(let7_5p_results) > 0) {
    cat("- let7_5p_knockdown.pdf/png: Let-7-5p knockdown efficiency\n")
}
cat("- top_de_mirnas.pdf/png: Top differentially expressed miRNAs\n")
cat("- heatmap_top_variable.pdf/png: Heatmap of most variable miRNAs\n")
cat("- sample_correlation.pdf: Sample correlation matrix\n")
cat("- summary_statistics.pdf/png: Analysis summary statistics\n")