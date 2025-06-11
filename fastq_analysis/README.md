# miRNA-seq Analysis Pipeline

A professional, modular pipeline for analyzing miRNA-seq data with a focus on Let-7-5p knockdown experiments.

## Overview

This pipeline provides a complete workflow for miRNA-seq analysis including:
- Quality control (FastQC/MultiQC)
- Read alignment to miRBase
- Expression quantification
- Differential expression analysis
- Comprehensive visualization
- HTML and Excel reporting

## Directory Structure

```
scripts/                    # Place all scripts here
├── run_mirna_analysis.sh   # Main runner script
├── 01_install_packages.R   # R package installation
├── 02_run_preprocessing.R  # QC and preprocessing
├── 03_run_main_analysis.R  # Main analysis pipeline
├── 04_create_visualizations.R  # Generate plots
└── 05_generate_report.R    # Create HTML report

data_directory/             # Your data location
├── XTRIA.*.tar            # Raw data archive
├── mature.fa              # miRBase reference (optional)
├── Electroporated/        # Sample directory
├── Notelectroporated/     # Sample directory
└── complete_analysis/     # Output directory (created)
    ├── 1_fastqc_raw/
    ├── 2_trimmed_reads/
    ├── 3_fastqc_trimmed/
    ├── 4_mapping/
    ├── 5_counts/
    ├── 6_normalization/
    ├── 7_differential_expression/
    ├── plots/
    ├── complete_analysis_results.xlsx
    └── analysis_report.html
```

## Prerequisites

### Required Software
- R (>= 4.0)
- Python 3
- FastQC (optional but recommended)
- MultiQC (optional but recommended)
- Cutadapt (optional, for adapter trimming)

### Installation
```bash
# macOS with Homebrew
brew install r
brew install fastqc
pip3 install multiqc cutadapt

# Linux
sudo apt-get install r-base
sudo apt-get install fastqc
pip3 install multiqc cutadapt
```

### R Packages
The pipeline will automatically install required R packages:
- tidyverse, ggplot2, openxlsx
- Rsubread, edgeR, limma (Bioconductor)
- pheatmap, RColorBrewer, scales

## Usage

### Basic Usage
1. Place all scripts in a common directory (e.g., `/path/to/scripts/`)
2. Navigate to your data directory
3. Run the main script:

```bash
cd /path/to/your/data
/path/to/scripts/run_mirna_analysis.sh
```

### Script Details

#### Main Runner Script
`run_mirna_analysis.sh` - Orchestrates the entire pipeline
- Checks prerequisites
- Installs R packages
- Extracts data
- Runs all analysis steps
- Opens results on completion

#### Individual R Scripts

1. **01_install_packages.R**
   - Checks and installs required R packages
   - Sets up Bioconductor if needed

2. **02_run_preprocessing.R**
   - Runs FastQC on raw and trimmed reads
   - Generates MultiQC reports
   - Performs adapter trimming if needed
   - Handles both R1 and R2 files

3. **03_run_main_analysis.R**
   - Main analysis workflow:
     - Downloads/prepares miRBase reference
     - Builds alignment index
     - Aligns reads (tries both SE and PE modes)
     - Quantifies miRNA expression
     - Normalizes counts (TMM)
     - Performs differential expression
     - Generates Excel report

4. **04_create_visualizations.R**
   - Creates publication-ready plots:
     - Volcano plot
     - MA plot
     - Let-7 expression plots
     - Top DE miRNAs
     - Expression heatmap
     - Sample correlation

5. **05_generate_report.R**
   - Generates comprehensive HTML report
   - Summarizes all findings
   - Includes embedded visualizations

### Advanced Usage

Run individual steps:
```bash
# Just preprocessing
Rscript /path/to/scripts/02_run_preprocessing.R /path/to/data

# Main analysis only
Rscript /path/to/scripts/03_run_main_analysis.R /path/to/data /path/to/scripts

# Generate plots
Rscript /path/to/scripts/04_create_visualizations.R /path/to/data

# Create HTML report
Rscript /path/to/scripts/05_generate_report.R /path/to/data
```

## Input Data Requirements

### Expected File Structure
- **TAR archive**: `XTRIA.*.tar` containing:
  - `Electroporated/ILLUMINA_DATA/*.fastq.gz`
  - `Notelectroporated/ILLUMINA_DATA/*.fastq.gz`

### File Naming Convention
- R1 files: `*_R1.fastq.gz` or `*_R1_trimmed.fastq.gz`
- R2 files: `*_R2.fastq.gz` or `*_R2_trimmed.fastq.gz`

### Reference Files
- `mature.fa` from miRBase (downloaded automatically if not present)

## Output Files

### Main Results
- `complete_analysis_results.xlsx` - Multi-sheet Excel workbook with all results
- `analysis_report.html` - Interactive HTML summary report

### Quality Control
- `1_fastqc_raw/multiqc_raw_report.html` - Raw data QC
- `3_fastqc_trimmed/multiqc_trimmed_report.html` - Trimmed data QC

### Analysis Results
- `5_counts/raw_counts.csv` - Raw count matrix
- `6_normalization/normalized_cpm.csv` - Normalized expression values
- `7_differential_expression/differential_expression_results.csv` - All DE results
- `7_differential_expression/let7_5p_specific_results.csv` - Let-7-5p focused results

### Visualizations
All plots are saved in both PDF and PNG formats:
- Volcano plot
- MA plot
- Let-7 family expression
- Top differentially expressed miRNAs
- Expression heatmap
- Sample correlation matrix

## Interpreting Results

### Knockdown Efficiency
- **Successful**: >70% reduction in Let-7-5p
- **Partial**: 30-70% reduction
- **Limited**: <30% reduction

### Statistical Considerations
With only 2 samples (no replicates):
- Fixed dispersion (BCV = 0.4) is assumed
- Results should be validated with additional experiments
- Consider "relaxed" significance criteria for exploratory analysis

### Key Metrics
1. **Alignment rate**: Should be >50% for good quality miRNA-seq
2. **Total miRNAs detected**: Typically 200-800 for human samples
3. **Let-7 family**: Should detect multiple Let-7 variants

## Troubleshooting

### Low Alignment Rate
- Check adapter trimming parameters
- Verify species (human miRNAs expected)
- Consider contamination or degradation

### No Let-7-5p Detection
- Check if Let-7 is expressed in your cell type
- Verify knockdown method efficiency
- Consider sequencing depth

### Script Errors
- Ensure all prerequisites are installed
- Check file paths and permissions
- Verify data file naming conventions

## Citation

If you use this pipeline, please cite:
- Rsubread: Liao et al. (2019) Bioinformatics
- edgeR: Robinson et al. (2010) Bioinformatics
- miRBase: Kozomara et al. (2019) Nucleic Acids Research

## Support

For issues or questions:
1. Check the HTML report for warnings
2. Review QC reports for data quality
3. Examine log files in `complete_analysis/logs/`
4. Validate key findings experimentally

## Version History

- v1.0: Initial release with modular structure
- Supports both single-end and paired-end miRNA-seq
- Automated report generation
- Comprehensive visualization suite