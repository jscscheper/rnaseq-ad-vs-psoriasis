# RNA-seq Differential Expression Analysis

![R](https://img.shields.io/badge/R-4.x-276DC3?logo=r)
![Python](https://img.shields.io/badge/Python-3.x-blue?logo=python)
![DESeq2](https://img.shields.io/badge/DESeq2-1.x-green)
![License](https://img.shields.io/badge/License-Educational-lightgrey)

A comprehensive RNA-seq differential expression analysis project comparing gene expression in atopic dermatitis and psoriasis samples.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Technology Stack](#technology-stack)
- [Data](#data)
- [Usage](#usage)
- [Author](#author)

---

## Overview

This project analyzes RNA-seq data from the GEO dataset GSE121212, comparing gene expression profiles between atopic dermatitis (AD) and psoriasis (PSO) conditions. The analysis includes:

- **Data preprocessing** and quality control
- **Differential expression analysis** using DESeq2
- **Pathway enrichment** and functional annotation
- **Visualization** including heatmaps and volcano plots

---

## Project Structure

```
rna_diff/
├── analysis/
│   └── ad_pso_analysis.Rmd             # Main R Markdown analysis
├── data/
│   ├── GSE121212_raw_counts*.tsv       # Raw count matrices
│   ├── GSE121212_series_matrix.txt     # Sample metadata
│   ├── drug_targets.csv                # Drug target annotations
│   ├── celltype_h3k27ac_gene_sets.rds  # Gene set data
│   └── roadmap_H3K27ac/                # Epigenetic annotations
├── results/
│   └── ad_pso_analysis.pdf             # Generated report
├── src/
│   └── helpers.R                       # Helper functions
└── README.md                           # This file
```

---

## Technology Stack

| Component | Technology | Purpose |
|-----------|------------|---------|
| **Language** | R 4.x | Statistical analysis |
| **DE Analysis** | DESeq2 | Differential expression |
| **Visualization** | ggplot2 | Plots and figures |
| **Data Source** | GEO (GSE121212) | RNA-seq counts |

---

## Data

The analysis uses publicly available data from NCBI GEO:

- **Dataset**: [GSE121212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121212)
- **Samples**: Skin biopsies from AD and PSO patients
- **Format**: Raw read counts (GRCh38.p13)

---

## Usage

1. **Open the analysis in RStudio**
   ```r
   # Navigate to rna_diff/analysis/
   # Open ad_pso_analysis.Rmd
   ```

2. **Install required packages**
   ```r
   install.packages("BiocManager")
   BiocManager::install(c("DESeq2", "ggplot2", "pheatmap"))
   ```

3. **Knit the document** to generate the report

---

## Author

**Jamie Scheper**

- GitHub: [@jscscheper](https://github.com/jscscheper)
