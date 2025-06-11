#!/usr/bin/env Rscript

# =============================================================================
# HTML Report Generation Script for miRNA-seq Analysis
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(knitr)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 05_generate_report.R <data_directory>")
}

data_dir <- args[1]
setwd(data_dir)

cat("\n=== Generating HTML Summary Report ===\n")

# Function to safely read CSV files
safe_read_csv <- function(file_path, default = data.frame()) {
    if (file.exists(file_path)) {
        tryCatch({
            read.csv(file_path, check.names = FALSE)
        }, error = function(e) {
            cat("Warning: Could not read", file_path, "\n")
            default
        })
    } else {
        cat("Warning: File not found:", file_path, "\n")
        default
    }
}

# Load results with error handling
de_results <- safe_read_csv("complete_analysis/7_differential_expression/differential_expression_results.csv")
let7_results <- safe_read_csv("complete_analysis/7_differential_expression/let7_family_results.csv")
let7_5p_results <- safe_read_csv("complete_analysis/7_differential_expression/let7_5p_specific_results.csv")

# Try to load or create alignment statistics
alignment_stats <- safe_read_csv("complete_analysis/alignment_statistics.csv")
if (nrow(alignment_stats) == 0) {
    # Create placeholder with properly trimmed read counts
    cat("Creating alignment statistics from available data...\n")
    alignment_stats <- data.frame(
        Sample = c("IL21602-001", "IL21603-001"),
        Mode = c("Single-end", "Single-end"),
        Total_reads = c(1028879, 1399238),  # From properly trimmed files
        Mapped_reads = c(NA, NA),
        Percent_mapped = c(NA, NA)
    )
    
    # Try to get actual alignment stats from count summary if available
    count_summary <- safe_read_csv("complete_analysis/5_counts/count_summary.csv")
    if (nrow(count_summary) > 0 && "Assigned" %in% rownames(count_summary)) {
        assigned_reads <- as.numeric(count_summary["Assigned", ])
        if (length(assigned_reads) >= 2) {
            alignment_stats$Mapped_reads <- assigned_reads[1:2]
            alignment_stats$Percent_mapped <- round(assigned_reads[1:2] / alignment_stats$Total_reads * 100, 2)
        }
    }
}

# Calculate summary statistics
total_mirnas <- nrow(de_results)
sig_strict <- if ("Significant_strict" %in% colnames(de_results)) sum(de_results$Significant_strict) else 0
sig_relaxed <- if ("Significant_relaxed" %in% colnames(de_results)) sum(de_results$Significant_relaxed) else 0
let7_family_count <- nrow(let7_results)
let7_5p_count <- nrow(let7_5p_results)

# Calculate Let-7-5p reduction
let7_reduction <- if (nrow(let7_5p_results) > 0 && "logFC" %in% colnames(let7_5p_results)) {
    round((1 - 2^mean(let7_5p_results$logFC)) * 100, 1)
} else {
    NA
}

# Determine knockdown status
knockdown_status <- if (!is.na(let7_reduction)) {
    if (let7_reduction > 70) {
        list(status = "SUCCESSFUL", class = "success", icon = "✓")
    } else if (let7_reduction > 30) {
        list(status = "PARTIAL", class = "warning", icon = "⚠")
    } else {
        list(status = "LIMITED", class = "error", icon = "✗")
    }
} else {
    list(status = "NOT DETECTED", class = "error", icon = "✗")
}

# Get current date and time
report_date <- format(Sys.time(), "%B %d, %Y at %I:%M %p")

# Create HTML report
html_content <- sprintf('<!DOCTYPE html>
<html>
<head>
    <title>miRNA-seq Analysis Report</title>
    <meta charset="UTF-8">
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            color: #2E86AB;
            text-align: center;
            margin-bottom: 10px;
        }
        h2 {
            color: #333;
            margin-top: 30px;
            border-bottom: 2px solid #2E86AB;
            padding-bottom: 10px;
        }
        h3 {
            color: #555;
            margin-top: 20px;
        }
        .date {
            text-align: center;
            color: #666;
            margin-bottom: 30px;
        }
        .summary-box {
            background-color: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
        }
        .highlight-box {
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
        }
        .warning-box {
            background-color: #f8d7da;
            border: 1px solid #f5c6cb;
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
        }
        .success {
            color: #28a745;
            font-weight: bold;
        }
        .warning {
            color: #ffc107;
            font-weight: bold;
        }
        .error {
            color: #dc3545;
            font-weight: bold;
        }
        .metric {
            display: flex;
            justify-content: space-between;
            padding: 10px 0;
            border-bottom: 1px solid #eee;
        }
        .metric:last-child {
            border-bottom: none;
        }
        .metric-label {
            font-weight: 600;
            color: #555;
        }
        .metric-value {
            font-weight: bold;
            color: #2E86AB;
        }
        table {
            border-collapse: collapse;
            width: 100%%;
            margin-top: 20px;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }
        th {
            background-color: #2E86AB;
            color: white;
            font-weight: bold;
        }
        tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        .plot-gallery {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        .plot-item {
            text-align: center;
        }
        .plot-item img {
            max-width: 100%%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 5px;
        }
        .file-list {
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin-top: 10px;
        }
        .file-list ul {
            margin: 10px 0;
        }
        .file-list li {
            margin: 5px 0;
            color: #555;
        }
        .file-list code {
            background-color: #e9ecef;
            padding: 2px 5px;
            border-radius: 3px;
            font-family: monospace;
        }
        .next-steps {
            background-color: #e7f3ff;
            border: 1px solid #b8daff;
            border-radius: 8px;
            padding: 20px;
            margin-top: 30px;
        }
        .next-steps ol {
            margin: 10px 0;
        }
        .next-steps li {
            margin: 10px 0;
        }
        .footer {
            text-align: center;
            color: #666;
            margin-top: 50px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>miRNA-seq Analysis Report</h1>
        <p class="date">Generated on %s</p>',
    report_date
)

# Add warning if analysis is incomplete
if (total_mirnas == 0) {
    html_content <- paste0(html_content, '
        <div class="warning-box">
            <h2>⚠️ Analysis May Be Incomplete</h2>
            <p>Some analysis results are missing. This could be due to:</p>
            <ul>
                <li>Low alignment rate</li>
                <li>Processing errors</li>
                <li>Incomplete pipeline execution</li>
            </ul>
            <p>Check the log files for more information.</p>
        </div>')
}

# Add main results section
html_content <- paste0(html_content, sprintf('
        <div class="highlight-box">
            <h2>🎯 Let-7-5p Knockdown Result</h2>
            <div class="metric">
                <span class="metric-label">Overall Let-7-5p reduction:</span>
                <span class="metric-value">%s%%</span>
            </div>
            <div class="metric">
                <span class="metric-label">Knockdown status:</span>
                <span class="%s">%s %s KNOCKDOWN</span>
            </div>
            <div class="metric">
                <span class="metric-label">Let-7-5p variants detected:</span>
                <span class="metric-value">%d</span>
            </div>
        </div>
        
        <h2>📊 Analysis Overview</h2>
        <div class="summary-box">
            <h3>Sample Information</h3>
            <table>
                <tr>
                    <th>Sample ID</th>
                    <th>Condition</th>
                    <th>Description</th>
                </tr>
                <tr>
                    <td>IL21602-001</td>
                    <td>Electroporated</td>
                    <td>Let-7-5p knockdown</td>
                </tr>
                <tr>
                    <td>IL21603-001</td>
                    <td>Not electroporated</td>
                    <td>Control</td>
                </tr>
            </table>
        </div>
        
        <div class="summary-box">
            <h3>Detection Summary</h3>
            <div class="metric">
                <span class="metric-label">Total miRNAs detected:</span>
                <span class="metric-value">%d</span>
            </div>
            <div class="metric">
                <span class="metric-label">Significantly changed (strict criteria):</span>
                <span class="metric-value">%d</span>
            </div>
            <div class="metric">
                <span class="metric-label">Significantly changed (relaxed criteria):</span>
                <span class="metric-value">%d</span>
            </div>
            <div class="metric">
                <span class="metric-label">Let-7 family members detected:</span>
                <span class="metric-value">%d</span>
            </div>
        </div>',
    ifelse(is.na(let7_reduction), "N/A", paste0(let7_reduction)),
    knockdown_status$class,
    knockdown_status$icon,
    knockdown_status$status,
    let7_5p_count,
    total_mirnas,
    sig_strict,
    sig_relaxed,
    let7_family_count
))

# Add alignment statistics if available
if (nrow(alignment_stats) > 0) {
    html_content <- paste0(html_content, '
        <div class="summary-box">
            <h3>Alignment Statistics</h3>
            <table>
                <tr>
                    <th>Sample</th>
                    <th>Mode</th>
                    <th>Total Reads</th>
                    <th>Mapped Reads</th>
                    <th>Percent Mapped</th>
                </tr>')
    
    for (i in 1:nrow(alignment_stats)) {
        mapped_str <- ifelse(is.na(alignment_stats$Mapped_reads[i]), 
                            "N/A", 
                            format(alignment_stats$Mapped_reads[i], big.mark = ","))
        percent_str <- ifelse(is.na(alignment_stats$Percent_mapped[i]), 
                             "N/A", 
                             paste0(alignment_stats$Percent_mapped[i], "%"))
        
        html_content <- paste0(html_content, sprintf('
                <tr>
                    <td>%s</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>%s</td>
                </tr>',
            alignment_stats$Sample[i],
            alignment_stats$Mode[i],
            format(alignment_stats$Total_reads[i], big.mark = ","),
            mapped_str,
            percent_str
        ))
    }
    
    html_content <- paste0(html_content, '
            </table>
        </div>')
}

# Add top Let-7-5p results if available
if (nrow(let7_5p_results) > 0 && all(c("miRNA", "logFC", "FC", "PValue") %in% colnames(let7_5p_results))) {
    top_let7 <- let7_5p_results %>%
        arrange(logFC) %>%
        head(5)
    
    html_content <- paste0(html_content, '
        <div class="summary-box">
            <h3>Top Let-7-5p Variants (Most Downregulated)</h3>
            <table>
                <tr>
                    <th>miRNA</th>
                    <th>Log2 FC</th>
                    <th>Fold Change</th>
                    <th>P-value</th>
                    <th>Reduction (%)</th>
                </tr>')
    
    for (i in 1:nrow(top_let7)) {
        reduction_pct <- round((1 - top_let7$FC[i]) * 100, 1)
        html_content <- paste0(html_content, sprintf('
                <tr>
                    <td>%s</td>
                    <td>%.2f</td>
                    <td>%.3f</td>
                    <td>%.2e</td>
                    <td>%s%%</td>
                </tr>',
            top_let7$miRNA[i],
            top_let7$logFC[i],
            top_let7$FC[i],
            top_let7$PValue[i],
            reduction_pct
        ))
    }
    
    html_content <- paste0(html_content, '
            </table>
        </div>')
}

# Add visualization gallery
html_content <- paste0(html_content, '
        <h2>📈 Visualizations</h2>
        <div class="plot-gallery">')

# Check which plots exist and add them
plots <- list(
    c("volcano_plot.png", "Volcano Plot"),
    c("ma_plot.png", "MA Plot"),
    c("let7_expression.png", "Let-7 Family Expression"),
    c("let7_5p_knockdown.png", "Let-7-5p Knockdown"),
    c("top_de_mirnas.png", "Top DE miRNAs"),
    c("summary_statistics.png", "Summary Statistics")
)

plots_found <- FALSE
for (plot_info in plots) {
    plot_file <- file.path("complete_analysis/plots", plot_info[1])
    if (file.exists(plot_file)) {
        plots_found <- TRUE
        html_content <- paste0(html_content, sprintf('
            <div class="plot-item">
                <h4>%s</h4>
                <img src="plots/%s" alt="%s">
            </div>',
            plot_info[2],
            plot_info[1],
            plot_info[2]
        ))
    }
}

if (!plots_found) {
    html_content <- paste0(html_content, '
            <p>No plots generated yet. Check if the visualization script completed successfully.</p>')
}

html_content <- paste0(html_content, '
        </div>')

# Add file list
html_content <- paste0(html_content, '
        <h2>📁 Generated Files</h2>
        <div class="file-list">
            <h3>Main Results</h3>
            <ul>')

# Check which files actually exist
if (file.exists("complete_analysis/complete_analysis_results.xlsx")) {
    html_content <- paste0(html_content, '
                <li><strong>Excel Report:</strong> <code>complete_analysis_results.xlsx</code> - Comprehensive results with multiple sheets</li>')
}

html_content <- paste0(html_content, '
                <li><strong>HTML Report:</strong> <code>analysis_report.html</code> - This summary report</li>
            </ul>
            
            <h3>Quality Control</h3>
            <ul>')

if (file.exists("complete_analysis/1_fastqc_raw/multiqc_raw_report.html")) {
    html_content <- paste0(html_content, '
                <li><strong>Raw data QC:</strong> <code>1_fastqc_raw/multiqc_raw_report.html</code></li>')
}

if (file.exists("complete_analysis/3_fastqc_trimmed/multiqc_trimmed_report.html")) {
    html_content <- paste0(html_content, '
                <li><strong>Trimmed data QC:</strong> <code>3_fastqc_trimmed/multiqc_trimmed_report.html</code></li>')
}

html_content <- paste0(html_content, '
            </ul>
            
            <h3>Analysis Results</h3>
            <ul>')

# List actual result files
result_files <- list(
    c("5_counts/raw_counts.csv", "Count data"),
    c("6_normalization/normalized_cpm.csv", "Normalized data"),
    c("7_differential_expression/differential_expression_results.csv", "DE results"),
    c("7_differential_expression/let7_5p_specific_results.csv", "Let-7 results")
)

for (file_info in result_files) {
    if (file.exists(file.path("complete_analysis", file_info[1]))) {
        html_content <- paste0(html_content, sprintf('
                <li><strong>%s:</strong> <code>%s</code></li>',
            file_info[2], file_info[1]
        ))
    }
}

html_content <- paste0(html_content, '
            </ul>
            
            <h3>Visualizations</h3>
            <ul>
                <li><strong>All plots:</strong> <code>plots/</code> directory (PDF and PNG formats)</li>
            </ul>
        </div>')

# Add next steps
html_content <- paste0(html_content, '
        <div class="next-steps">
            <h2>🚀 Next Steps</h2>
            <ol>')

if (file.exists("complete_analysis/complete_analysis_results.xlsx")) {
    html_content <- paste0(html_content, '
                <li><strong>Review the Excel file</strong> for detailed results across all miRNAs</li>')
} else {
    html_content <- paste0(html_content, '
                <li><strong>Check analysis logs</strong> to troubleshoot any errors</li>')
}

html_content <- paste0(html_content, '
                <li><strong>Check QC reports</strong> to ensure data quality meets standards</li>
                <li><strong>Validate key findings</strong> with qRT-PCR, especially for Let-7-5p variants</li>
                <li><strong>Perform pathway analysis</strong> on significantly changed miRNAs</li>
                <li><strong>Consider biological replicates</strong> for more robust statistical analysis</li>
            </ol>
        </div>')

# Add footer
html_content <- paste0(html_content, '
        <div class="footer">
            <p>Generated by miRNA-seq Analysis Pipeline<br>
            For questions or issues, please refer to the documentation</p>
        </div>
    </div>
</body>
</html>')

# Write HTML file
output_file <- file.path("complete_analysis", "analysis_report.html")
writeLines(html_content, output_file)

# Save alignment statistics if we created it
if (!file.exists("complete_analysis/alignment_statistics.csv")) {
    write.csv(alignment_stats, 
              file.path("complete_analysis", "alignment_statistics.csv"), 
              row.names = FALSE)
}

cat("\n✓ HTML report generated successfully!\n")
cat("Report saved as:", output_file, "\n\n")

# Print summary to console
cat("=== Analysis Summary ===\n")

if (total_mirnas > 0) {
    cat("Let-7-5p knockdown:", ifelse(is.na(let7_reduction), "Not detected", paste0(let7_reduction, "%")), "\n")
    cat("Status:", knockdown_status$icon, knockdown_status$status, "KNOCKDOWN\n")
    cat("\nKey findings:\n")
    cat("- Total miRNAs detected:", total_mirnas, "\n")
    cat("- Significantly changed:", sig_relaxed, "\n")
    cat("- Let-7-5p variants:", let7_5p_count, "\n")
} else {
    cat("WARNING: No analysis results found.\n")
    cat("Check if the main analysis completed successfully.\n")
    cat("\nPossible issues:\n")
    cat("- Very low alignment rate\n")
    cat("- Pipeline stopped early\n")
    cat("- Wrong species reference\n")
}

cat("\nGenerated files:\n")
cat("- HTML report: complete_analysis/analysis_report.html\n")

if (file.exists("complete_analysis/complete_analysis_results.xlsx")) {
    cat("- Excel report: complete_analysis/complete_analysis_results.xlsx\n")
}

cat("\nOpen the HTML report in a web browser for detailed results.\n")