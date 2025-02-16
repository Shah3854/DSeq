# Differential Expression Analysis using DESeq2

## Overview
This repository contains an R script for performing differential gene expression analysis using the DESeq2 package. The workflow includes data preprocessing, normalization, statistical testing, and visualization of results through various plots.

## Requirements
### R Packages
Ensure you have the following R packages installed before running the script:

```r
install.packages(c("ggplot2", "gplots", "pheatmap", "RColorBrewer", "geneplotter", "Rsubread", "dplyr", "genefilter", "stringr", "LSD", "apeglm", "tidyverse", "matrixStats"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "DESeq2", "EnhancedVolcano"))
```

## Files and Structure
```
├── PCA_raw.txt             # Raw count matrix (tab-separated)
├── condition.csv           # Sample metadata
├── result.csv              # DESeq2 full results
├── filter.csv              # Filtered results (log2FC > 1.5 & p-value <= 0.05)
├── norm_counts.csv         # Normalized count data
├── maplot.pdf              # MA plot
├── boxplot.pdf             # Boxplot of normalized counts
├── scatterplot.pdf         # Scatter plot of PC vs Treatment
├── volcano.pdf             # Basic volcano plot
├── heatmap.pdf             # Heatmap of top 50 differentially expressed genes
├── enhanced_volcano.pdf    # Enhanced volcano plot
├── script.R                # Main R script
└── README.md               # Documentation
```

## Usage
1. Prepare your input files:
   - `PCA_raw.txt`: A tab-separated file with gene counts, rows as genes, columns as samples.
   - `condition.csv`: A CSV file with sample names and their respective conditions.

2. Run the script:
   ```sh
   Rscript script.R
   ```

3. Output files, including DE results and plots, will be generated in the specified `Path/` directory.

## Plots and Visualization
- **MA Plot**: Shows the relationship between mean expression and fold-change.
- **Boxplot**: Visualizes the distribution of normalized counts across samples.
- **Scatter Plot**: Compares average counts of the positive control and treatment.
- **Volcano Plot**: Highlights significantly up/downregulated genes.
- **Heatmap**: Displays expression patterns of the top 50 differentially expressed genes.
- **Enhanced Volcano Plot**: A more detailed volcano plot using `EnhancedVolcano` package.
