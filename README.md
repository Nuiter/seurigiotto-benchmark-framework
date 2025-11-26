# seurigiotto-benchmark-framework

![seurigiotto-framework](https://img.shields.io/badge/seurigiotto-framework-9cf?colorA=7A009C&colorB=4CAF50&logo=r)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17628633.svg)](https://doi.org/10.5281/zenodo.17628633)
![R Version](https://img.shields.io/badge/R-4.4%2B-blue?logo=r)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/status-active-success.svg)
![GitHub last commit](https://img.shields.io/github/last-commit/Nuiter/seurigiotto-benchmark-framework)

![Seurat](https://img.shields.io/badge/Built%20with-Seurat-2777B6?logo=r)
![Giotto](https://img.shields.io/badge/Built%20with-Giotto-4CAF50?logo=r)


**Optimized pipelines for Spatial Transcriptomics (ST) data analysis using Seurat & Giotto, designed for reproducible benchmarking and biological insight.**

*MSc Thesis Project | Bioinformatics & Biostatistics, UB-UOC | 2025*

---

### Project Overview

This repository provides a comprehensive and reproducible suite of R scripts for analyzing Spatial Transcriptomics (ST) data. It compares the performance and results of **Seurat** and **Giotto**, and introduces a novel **Optimized Hybrid Pipeline** that leverages the strengths of both.

Developed as part of a Master's Thesis, this project is structured as a robust resource for the wider bioinformatics community, focusing on lung cancer research using 10x Genomics Visium HD data.

### Key Features

- **Standalone Pipelines:** Fully documented, independent scripts for Seurat and Giotto using harmonized parameters for fair comparison.
- **Hybrid Pipeline:** A novel workflow that combines Seurat’s efficiency for initial analysis with Giotto’s granularity for detailed clustering, governed by an automated decision step.
- **Rigorous Benchmarking:** Measures key performance metrics including memory usage (RAM), execution time, and clustering resolution.
- **In-depth Downstream Analysis:** Includes modules for marker gene detection, cell-type signature scoring, simulated deconvolution, and advanced spatial visualization.
- **100% Reproducible:** All code is open, well-documented, and designed for reproducibility.

### Note on Resolution & Scalability:

- The provided scripts and benchmarks are configured for 16 µm resolution binning. This resolution was chosen to balance spatial granularity with accessibility, allowing the full workflow (including Giotto) to run on standard workstations (32GB - 64GB RAM).
- However, the framework is resolution-agnostic. To analyze finer resolutions (8 µm or 2 µm), simply provide the corresponding input matrix. Caution: High-resolution bins (2 µm) scale memory requirements exponentially.

### Repository Structure
```
📁 scripts/
    📜 01_seurat_pipeline.R
    📜 02_giotto_pipeline.R
    📜 03_compare_pipelines.R
    📜 4_Optimized_pipeline.R
    📜 5_Integrated_analysis.R
    📜 6_native_visualizations.R
    📜 7_final_report_generation.R
    📜 8_deconvolution_analysis.R
    📜 9_balanced_deconvolution.R
    
📁 data/
    📄 Tumor.csv
    📄 Fibroblasts.csv
    📄 ... (and other signature files))
    📁 visium_dataset/
        *(Raw spatial dataset - to be downloaded separately)*
        
📁 results/ # (Git-ignored) Where all outputs are saved.

📄 .gitignore
📄 README.md

```
### Getting Started
```
#### 1. Clone the Repository
git clone https://github.com/Nuiter/seurigiotto-benchmark-framework.git
cd seurigiotto-benchmark-framework
```
#### 2. Install Dependencies
## Reproducibility snapshot
This repository reflects the core code and environment used for the MSc thesis defence (June 2025).  
Key versions:

| Component   | Version |
|:---|:---|
| R           | 4.4.2   |
| Seurat      | 5.3.0   |
| Giotto      | 4.2.1   |
| OS          | Linux   |


- **R version:** 4.4.2 or higher recommended.  
- **Key R Packages:**  
  Seurat, Giotto, ggplot2, dplyr, patchwork, pryr, arrow, data.table, scales, RColorBrewer, viridis, pheatmap, tidyr.

You can install them in your R session with:

```
install.packages(c("Seurat", "Giotto", "ggplot2", "dplyr", "patchwork", "pryr", "arrow", "data.table", "scales", "RColorBrewer", "viridis", "pheatmap", "tidyr"))
```
(Note: Seurat and Giotto may require installation from Bioconductor or specific sources. Please refer to their official documentation.)


#### 3. Data Setup

The analysis scripts are configured to use public data from 10x Genomics (Visium HD Human Lung Cancer).

- The cell signature CSV files are included in the /data directory.  
- The main spatial dataset must be downloaded separately. Place the dataset folder (e.g., visium_dataset) inside the /data directory.
- Link: [Visium HD Human Lung Cancer Dataset](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-lung-cancer-fixed-frozen)

#### 4. Run the Analysis Workflow

The scripts are designed to be run sequentially from the command line.  
The core benchmarking workflow consists of scripts 1, 2, and 3.  
The others perform more advanced, integrated analyses.

```bash
# === Core Benchmarking Workflow ===
# 1. Run the standalone Seurat pipeline
Rscript scripts/01_seurat_pipeline.R

# 2. Run the standalone Giotto pipeline
Rscript scripts/02_giotto_pipeline.R

# 3. Run the comparative analysis
Rscript scripts/03_compare_pipelines.R

# === Advanced & Integrated Analysis ===
# These scripts build upon the initial results for deeper insights.
Rscript scripts/4_Optimized_pipeline.R
Rscript scripts/5_Integrated_analysis.R
Rscript scripts/6_native_visualizations.R
Rscript scripts/7_final_report_generation.R
Rscript scripts/8_deconvolution_analysis.R
Rscript scripts/9_balanced_deconvolution.R

```

All outputs (figures, tables, logs) will be saved in the /results directory in organized subfolders.

### Citation

If you use this workflow or find the code helpful in your research, please cite this repository.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17628633.svg)](https://doi.org/10.5281/zenodo.17628633)

Author: Ángel I. Pérez Santiago  
Project: MSc in Bioinformatics & Biostatistics (University of Barcelona / UOC, 2025)  
GitHub: https://github.com/nuiter

### License

This project is licensed under the MIT License. See the LICENSE file for full details.

### Acknowledgements

- Dataset: Visium HD CytAssist Gene Expression of Human Lung Cancer (Fixed Frozen) by 10x Genomics.  
- Supervision: Dr. Alfonso Saera Vila (MSc Thesis Advisor).  
- Core Frameworks: The Satija Lab for [Seurat](https://satijalab.org/seurat/) and the Dries Lab for [Giotto](https://giottosuite.com/).

---

For questions, collaborations, or feedback, please open an issue in this repository.
