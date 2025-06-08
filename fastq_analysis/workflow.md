# Complete miRNA-seq Workflow Execution Guide

## Overview
This guide follows your exact workflow requirements:
1. FastQC/MultiQC on raw data
2. Adapter trimming
3. FastQC/MultiQC on trimmed data
4. Reference mapping and counting
5. Normalization and differential expression
6. Pathway analysis

## Prerequisites Installation

### 1. Install Command-Line Tools

```bash
# On macOS with Homebrew:
brew install fastqc
pip3 install multiqc cutadapt

# On Linux:
sudo apt-get install fastqc
pip install multiqc cutadapt
```

### 2. Install R and Required Packages

In R or RStudio:
```r
# Set CRAN mirror
options(repos = "https://cloud.r-project.org")

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
BiocManager::install(c("Rsubread", "edgeR", "limma"))
install.packages(c("tidyverse", "ggplot2", "pheatmap", 
                   "RColorBrewer", "knitr", "openxlsx"))
```

## Step-by-Step Execution

### Step 1: Navigate to Your Data Directory
```bash
cd /Users/rvidva/Documents/datasets/XTRIA
```

### Step 2: Run Preprocessing (FastQC, Trimming, MultiQC)

Save the bash script as `preprocess_pipeline.sh` and run:
```bash
# Make executable
chmod +x preprocess_pipeline.sh

# Run preprocessing
./preprocess_pipeline.sh
```

This will:
- Extract your TAR file if needed
- Run FastQC on raw reads
- Handle adapter trimming (or use pre-trimmed files)
- Run FastQC on trimmed reads
- Generate MultiQC reports
- Create a preprocessing summary

### Step 3: Review Quality Control Reports

Open the MultiQC HTML reports in your browser:
- Raw data: `complete_analysis/1_fastqc_raw/multiqc_raw_report.html`
- Trimmed data: `complete_analysis/3_fastqc_trimmed/multiqc_trimmed_report.html`

Check for:
- Read quality scores (should be >30)
- Adapter content (should be minimal after trimming)
- Read length distribution (18-30nt for miRNAs)

### Step 4: Run Complete R Analysis

Save the R script as `complete_workflow.R` and run:
```bash
Rscript complete_workflow.R
```

Or in RStudio:
```r
setwd("/Users/rvidva/Documents/datasets/XTRIA")
source("complete_workflow.R")
```

This will perform:
- Reference mapping to miRBase
- Read counting
- TMM normalization
- Differential expression analysis
- Pathway analysis preparation
- Generate comprehensive results

### Step 5: Review Results

The analysis creates a `complete_analysis` folder with:

```
complete_analysis/
├── 1_fastqc_raw/               # Raw data QC
├── 2_trimmed_reads/             # Trimmed sequences
├── 3_fastqc_trimmed/            # Trimmed data QC
├── 4_mapping/                   # Alignment files
│   ├── alignment_statistics.csv # Mapping rates
│   └── bam_files/              # Aligned reads
├── 5_counts/                    # Raw count data
├── 6_normalization/             # Normalized data
│   ├── normalized_cpm.csv      # CPM values
│   └── normalization_qc_plots.pdf
├── 7_differential_expression/   # DE results
│   ├── differential_expression_results.csv
│   ├── let7_5p_differential_expression.csv
│   ├── volcano_plot.pdf
│   └── ma_plot.pdf
├── 8_pathway_analysis/          # Files for pathway tools
│   ├── mirpath_input.csv       # For DIANA-miRPath
│   └── mirsystem_input.txt     # For miRSystem
├── complete_analysis_results.xlsx  # All results in Excel
└── analysis_summary_report.txt     # Summary report
```

### Step 6: Interpret Key Results

1. **Check Let-7-5p knockdown efficiency**:
   - Open `let7_5p_differential_expression.csv`
   - Look for negative logFC values (indicating reduction)
   - Calculate knockdown: (1 - 2^logFC) × 100%

2. **Review differential expression**:
   - Open `complete_analysis_results.xlsx`
   - Check "DE_Results" sheet for all miRNAs
   - Significant: |logFC| > 1 and FDR < 0.05

3. **Examine QC plots**:
   - `normalization_qc_plots.pdf` shows data quality
   - `volcano_plot.pdf` shows overall changes
   - `ma_plot.pdf` highlights significant miRNAs

### Step 7: Pathway Analysis

Use the prepared files for web-based tools:

1. **DIANA-miRPath**:
   - Go to: http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=mirpath/index
   - Upload `mirpath_input.csv`
   - Select pathways of interest

2. **miRSystem**:
   - Go to: http://mirsystem.cgm.ntu.edu.tw/
   - Upload `mirsystem_input.txt`
   - Analyze target genes and pathways

3. **miRNet**:
   - Go to: https://www.mirnet.ca/
   - Input your significant miRNAs
   - Create network visualizations

## Troubleshooting

### If FastQC fails:
- Check file permissions: `ls -la *.fastq.gz`
- Ensure enough disk space: `df -h`

### If R analysis fails:
- Check package installation
- Reduce threads if memory limited
- Check file paths are correct

### If no reads map:
- Verify you're using human samples
- Check read length distribution
- Ensure adapter trimming worked

## Expected Outcomes

For successful Let-7-5p knockdown:
- **Mapping rate**: >60% to miRNAs
- **Let-7-5p reduction**: >70% in electroporated
- **Volcano plot**: Let-7-5p variants below y=0 line
- **Pathway analysis**: Affected pathways related to Let-7 targets

## Quick Validation

For immediate Let-7-5p check without full pipeline:
```bash
# Electroporated sample
echo "Electroporated (IL21602-001):"
zcat XTRIA_20250520_A00904_IL21602-001_PE2V4-B02_L003_R1_trimmed.fastq.gz | \
  grep -c "TGAGGTAGTAGGTT"

# Control sample  
echo "Control (IL21603-001):"
zcat XTRIA_20250520_A00904_IL21603-001_PE2V4-C02_L003_R1_trimmed.fastq.gz | \
  grep -c "TGAGGTAGTAGGTT"
```

The control should have significantly more hits.