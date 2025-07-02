# Transcriptomic Analysis of Young and Aged Human MSC

This repository contains the full R code used in the manuscript:

## 📝 Associated Manuscript

**Title:** Exploring the Involvement of Ferroptosis-Associated Genes and Pathways in Mesenchymal Stem Cell Aging Through Bioinformatics Analysis  
**Authors:** Laleh Mavaddatiyan¹, Yasaman Khamineh¹, Mahmood Talkhabi¹*  
**Affiliation:** Department of Animal Sciences and Biotechnology, Faculty of Life Sciences and Biotechnology, Shahid Beheshti University, Tehran, Iran  
**Corresponding Author:** m_talkhabi@sbu.ac.ir

---

## 🧬 Overview

Mesenchymal stem cells (MSCs) derived from fetal, adult, and aged tissues exhibit distinct transcriptomic profiles associated with aging. To investigate these differences, we analyzed gene expression data from three publicly available GEO datasets:

- **GSE68374**: Fetal vs Adult MSCs
- **GSE97311**: Fetal vs Aged MSCs
- **GSE119987**: Early vs Late passage MSCs 

The pipeline includes:
- Data acquisition from GEO
- Normalization and transformation
- Differential expression analysis using `limma`
- Gene annotation via `biomaRt`
- Gene Set Enrichment Analysis (GSEA) using `fgsea` and `msigdbr`
- Heatmap visualizations of selected gene sets

---

## 📁 Repository Structure

- `MSC_analysis.R`: Complete R script including all steps for:
  - DEG identification
  - GSEA (KEGG and Hallmark pathways)
  - Heatmap generation
- `README.md`: Documentation of the analysis and instructions
- Output files:
  - Differential gene expression tables (`.txt`)
  - Enrichment plots (`.png`, `.pdf`)

---

## 📦 Required R Packages

Install these packages before running the analysis:

```r
install.packages(c("pheatmap", "RColorBrewer", "ggplot2", "cowplot",
                   "data.table", "dplyr", "tidyr", "gridExtra", "magrittr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "fgsea", "biomaRt", "msigdbr", "enrichplot"))
```

---

## 🚀 Usage Instructions

1. Clone or download this repository.
2. Open `MSC_analysis.R` in R or RStudio.
3. Ensure all input files (e.g., expression matrices, platform annotations) are correctly referenced.
4. Run the script section by section for each dataset:
   - **GSE68374**: Fetal vs Adult MSCs
   - **GSE97311**: Fetal vs Aged MSCs
   - **GSE119987**: Early vs Late passage MSCs

Each section includes:
- Normalization & DEG analysis
- GSEA with KEGG & Hallmark gene sets
- Heatmap visualizations

---

## 🔬 Output Examples

- `heatmap_sample_groups.pdf`: Heatmap of selected genes
- `final_deg.txt`: Filtered differentially expressed genes
- `Enrichment plot1.png` ... `Hallmark_enrichment_50.png`: GSEA plots

---

## ⚠️ Notes

- Some objects (`annot`, `annot2`, `illumina`, `heatmap.raw`, etc.) must be preloaded or preprocessed from the GEO platform data or user-provided data.
- Be sure to adjust file paths in `setwd()` or use `here::here()` for reproducibility.

---

## 📄 Citation

If you use this code in your work, please cite our manuscript.
---

## 📬 Contact

For questions or feedback, please contact:

**Corresponding Author**  
[Mahmood Talkhabi]  
[Department of Animal Sciences and Biotechnology, Faculty of Life Sciences and Biotechnology, Shahid
Beheshti University, Tehran, Iran]  
Email: [m_talkhabi@sbu.ac.ir]
