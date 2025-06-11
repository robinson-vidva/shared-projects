#!/usr/bin/env Rscript

# =============================================================================
# R Package Installation Script for miRNA-seq Analysis
# =============================================================================

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Checking and installing required R packages...\n\n")

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
        cat("✓", pkg, "already installed\n")
    }
}

# CRAN packages
install_if_needed("BiocManager")
install_if_needed("tidyverse")
install_if_needed("ggplot2")
install_if_needed("openxlsx")
install_if_needed("RColorBrewer")
install_if_needed("pheatmap")
install_if_needed("knitr")
install_if_needed("scales")
install_if_needed("gridExtra")

# Bioconductor packages
install_if_needed("Rsubread", bioc = TRUE)
install_if_needed("edgeR", bioc = TRUE)
install_if_needed("limma", bioc = TRUE)
install_if_needed("DESeq2", bioc = TRUE)
install_if_needed("biomaRt", bioc = TRUE)
install_if_needed("ShortRead", bioc = TRUE)

cat("\n✓ Package installation check complete\n")