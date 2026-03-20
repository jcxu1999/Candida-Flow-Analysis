# Candida albicans Multi-color Flow Cytometry Analysis

This repository contains an automated R-based pipeline for de-convolving complex cell mixtures in ***Candida albicans*** using a **Combinatorial Labeling (Fluorescence Barcoding)** strategy.

## Analysis Pipeline
The core script (`analyze_mixed_populations.R`) performs high-dimensional data processing:

1.  **Preprocessing**: Arcsinh Transformation (Cofactor = 150).
2.  **Dimensionality Reduction**: **UMAP** for 2D visualization of 21 clusters.
3.  **Classification**: **k-Nearest Neighbors (kNN, k=15)** to predict cell identities.
