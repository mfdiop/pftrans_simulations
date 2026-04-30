# Benchmarking Framework: Simulation and IBD/IBS Validation Pipeline

## 1. Overview

This document provides a complete description of the simulation and benchmarking framework used to evaluate the performance of different genetic relatedness inference methods (IBD, IBS, and phylogenetic inference) in accurately recovering transmission links. The framework integrates **SLiM**, **pyslim**, **msprime**, **tskit**, and **tskibd**, with validation and visualization components in **Python**.

The central goal is to determine under which biological and sampling conditions these methods most accurately reconstruct transmission dynamics, using *Plasmodium falciparum* as the main recombinant pathogen model. The framework remains extensible to non-recombinant pathogens.

---

## 2. Simulation Framework and Core Workflow

The pipeline consists of two main scripts:

1. **SLiM Simulation Script** – defines population structure, recombination, selection, mutation, and sampling logic. The output is a tree sequence (`.trees`) containing complete genealogical information.
2. **Python Post-Processing Script** – processes the `.trees` file using **pyslim**, **tskit**, and **tskibd** to extract genotypes (VCF), compute IBD segments, introduce experimental noise, and produce benchmarking datasets.

### Simulation Parameters
- **Population size**: adjustable via `N` (e.g., 1000).
- **Genome length**: typically 1e6–1e7 bp.
- **Mutation rate (μ)**: varied between 1e−8 to 1e−6 per base per generation.
- **Recombination rate (r)**: varied between 1e−9 to 1e−6 per base per generation.
- **Selection coefficient (s)**: optional, defining beneficial or deleterious mutations.
- **Sampling scheme**: defined by generation, coverage, and proportion.

---

## 3. Validation Strategy

The validation strategy introduces systematic variation in simulation parameters and evaluates the ability of different inference methods to recover true genetic relationships.

### Variables Tested
1. **Mutation rate** – affects SNP density and information content for relatedness inference.
2. **Recombination rate** – determines haplotype shuffling and IBD fragment resolution.
3. **Sample size and coverage** – impacts signal-to-noise ratio in relatedness estimation.
4. **Selection regime** – explores distortion of relatedness under linked selection.

### Optional Extensions
- **Noise injection** – simulate genotyping error (allele flip probability) and sequencing depth variation.
- **Temporal sampling** – label samples by generation to assess degradation of relatedness inference across time.
- **Transmission mapping** – extract true adjacency (parent–offspring) matrix from pedigree for direct link recovery evaluation.

---

## 4. Benchmarking and Evaluation Metrics

### Overview
Benchmarking focuses on the agreement between true relatedness (pedigree or simulated IBD) and inferred relatedness (IBD, IBS, or tree-based distances). Performance is summarized across parameter combinations using the following core metrics.

### 4.1 Correlation Metrics
- **Pearson correlation (r)** and **Spearman correlation (ρ)** between true and inferred IBD proportions.
- Computed across all pairs of individuals.
- Mean and variance of correlations are tracked per simulation condition.

### 4.2 Recovery Accuracy
- Define the *top-k related pairs* based on true IBD (e.g., top 5% highest relatedness).
- Compute the fraction of those pairs recovered among the top-k inferred pairs.
- Accuracy = |intersection(true, inferred)| / |true|.

### 4.3 Segment Length Error Distribution
- For all true IBD segments, compute the absolute difference between true and inferred segment lengths.
- Summarize via mean absolute error (MAE), variance, and quantiles (e.g., 5th–95th percentile range).

### 4.4 Transmission Recovery Metrics
- Represent the true parent–offspring links as an adjacency matrix.
- Infer pairwise connectivity using IBD/IBS thresholds.
- Compute **Precision**, **Recall**, and **F1-score** for recovered transmission pairs:
  - Precision = TP / (TP + FP)
  - Recall = TP / (TP + FN)
  - F1 = 2 × (Precision × Recall) / (Precision + Recall)

### 4.5 Summary Statistics
For each simulation condition, summarize:
- Mean correlation (IBD true vs inferred)
- Variance of correlation
- Top-related pair recovery accuracy
- Segment length MAE
- Transmission precision, recall, and F1-score

All metrics are stored in a structured results file (`results_metrics.csv`) for reproducibility.

---

## 5. Visualization and Reporting

### 5.1 Correlation Plots
- Scatterplots of inferred vs. true IBD proportions with regression lines.
- Heatmaps showing mean correlation by mutation and recombination rate combinations.

### 5.2 Segment Length Distributions
- Histograms comparing true and inferred IBD segment length distributions.
- Density plots highlighting error magnitude.

### 5.3 Transmission Accuracy Visualization
- ROC and Precision–Recall curves for recovered transmission pairs.
- Network visualizations overlaying true and inferred links (colored by accuracy category).

### 5.4 Summary Heatmaps
- Correlation heatmap: mutation rate × recombination rate.
- Accuracy heatmap: sample size × sequencing coverage.
- Temporal degradation plot: correlation vs. sampling interval.

---

## 6. Reproducibility

- All simulations use deterministic random seeds.
- Parameter grids are logged in `params.json`.
- Output structure:
  ```
  /simulation_results/
  ├── params.json
  ├── VCF_files/
  ├── IBD_truth/
  ├── IBD_inferred/
  ├── IBS_matrices/
  ├── transmission_links/
  └── results_metrics.csv
  ```
- Metadata include timestamp, seed, mutation/recombination rates, selection mode, and sampling regime.

---

## 7. Next Steps

1. Integrate the Python benchmarking driver to automate metric computation.
2. Extend the framework with cross-method comparisons (IBD vs. IBS vs. phylogenetic distance).
3. Develop an R Shiny or Python dashboard for interactive visualization of benchmarking outcomes.
4. Incorporate real genomic data for empirical validation of the simulation-informed conclusions.

---

## 8. Summary

This framework creates a reproducible, parameterized simulation environment to evaluate the accuracy and robustness of transmission inference methods. It directly connects ground-truth genealogical data from tree sequences to observed genomic patterns, providing a biologically realistic, statistically rigorous, and computationally efficient foundation for transmission benchmarking research.

