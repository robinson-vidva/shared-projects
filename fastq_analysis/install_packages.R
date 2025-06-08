#!/usr/bin/env Rscript

# =============================================================================
# Package Installation Script for Let-7-5p Analysis
# =============================================================================
# Run this script first to install all required packages
# =============================================================================

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("=== Installing Required Packages for miRNA Analysis ===\n\n")

# Function to check if a package is installed
is_installed <- function(pkg) {
  pkg %in% installed.packages()[, "Package"]
}

# Install BiocManager if needed
if (!is_installed("BiocManager")) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}

# Bioconductor packages
bioc_packages <- c("Rsubread", "edgeR", "limma")
cat("\nChecking Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!is_installed(pkg)) {
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(pkg, "is already installed.\n")
  }
}

# CRAN packages
cran_packages <- c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer", 
                   "knitr", "ggrepel", "scales")
cat("\nChecking CRAN packages...\n")
for (pkg in cran_packages) {
  if (!is_installed(pkg)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(pkg, "is already installed.\n")
  }
}

# Verify installation
cat("\n=== Verifying Installation ===\n")
all_packages <- c(bioc_packages, cran_packages)
installed <- sapply(all_packages, is_installed)

if (all(installed)) {
  cat("\n✅ All packages successfully installed!\n")
  cat("You can now run the analysis pipeline.\n")
} else {
  cat("\n❌ Some packages failed to install:\n")
  failed <- all_packages[!installed]
  cat(paste("  -", failed), sep = "\n")
  cat("\nTry installing these manually:\n")
  
  failed_bioc <- failed[failed %in% bioc_packages]
  if (length(failed_bioc) > 0) {
    cat("\nFor Bioconductor packages:\n")
    cat("BiocManager::install(c(", paste0('"', failed_bioc, '"', collapse = ", "), "))\n")
  }
  
  failed_cran <- failed[failed %in% cran_packages]
  if (length(failed_cran) > 0) {
    cat("\nFor CRAN packages:\n")
    cat("install.packages(c(", paste0('"', failed_cran, '"', collapse = ", "), "))\n")
  }
}

# Test loading packages
cat("\n=== Testing Package Loading ===\n")
load_success <- TRUE
for (pkg in all_packages) {
  tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    cat("✓", pkg, "loads successfully\n")
  }, error = function(e) {
    cat("✗", pkg, "failed to load:", e$message, "\n")
    load_success <<- FALSE
  })
}

if (load_success) {
  cat("\n🎉 All packages are ready! You can now run the analysis.\n")
} else {
  cat("\n⚠️  Some packages failed to load. Please check the errors above.\n")
}

# Session info
cat("\n=== Session Information ===\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("BiocManager version:", as.character(packageVersion("BiocManager")), "\n")