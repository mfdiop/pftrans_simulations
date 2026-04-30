# Transmission Discrimination Analysis Under Recombination

## Overview

This analysis evaluates whether **high recombination rates improve IBD-based methods' ability to discriminate recent transmission from distant ancestry**, even as overall detection sensitivity decreases due to shorter IBD segments.

## Research Question

**Does recombination help IBD methods separate the signal of recent transmission from ancestral relatedness?**

### Key Hypotheses

1. **H1 (Detection)**: IBD detection sensitivity decreases with recombination
   - Mechanism: Shorter segments harder to detect
   - Prediction: AUPR for overall relatedness ↓ as recombination ↑

2. **H2 (Discrimination)**: IBD transmission discrimination improves with recombination
   - Mechanism: Distant ancestry loses detectable segments faster than recent transmission
   - Prediction: AUPR for transmission vs. ancestry ↑ or stable as recombination ↑

3. **H3 (Enrichment)**: Among detectable pairs, recent transmission is enriched at high recombination
   - Mechanism: Recombination filters out distant ancestry preferentially
   - Prediction: Enrichment factor in top K pairs ↑ as recombination ↑

## Ground Truth Definitions

| Category | Generation Distance | Symbol | Interpretation |
|----------|-------------------|---------|----------------|
| **Transmission** | ≤5 generations | T | Recent, direct transmission chains |
| **Recent Ancestry** | 6-15 generations | A₁ | Detectable relatedness, not transmission |
| **Distant Ancestry** | 16-25 generations | A₂ | May become undetectable at high recombination |
| **Unrelated** | >25 generations | U | Background population |

## Evaluation Tasks

The analysis evaluates three distinct classification tasks:

### Task 1: Overall Relatedness Detection
- **Question**: Can we detect any relatedness?
- **Positive class**: T ∪ A₁ ∪ A₂ (all ≤25 gen)
- **Negative class**: U (>25 gen)
- **Expected**: IBD performs worse at high recombination (shorter segments)

### Task 2: Transmission Detection
- **Question**: Can we identify recent transmission?
- **Positive class**: T (≤5 gen)
- **Negative class**: A₁ ∪ A₂ ∪ U (>5 gen)
- **Expected**: IBD may be stable or improve (better separation)

### Task 3: Transmission Discrimination
- **Question**: Among related pairs, can we distinguish transmission from ancestry?
- **Positive class**: T (≤5 gen)
- **Negative class**: A₁ ∪ A₂ (6-25 gen)
- **Expected**: IBD improves at high recombination (ancestry filtered out)

## Methods Compared

| Method | Metric Type | Segment-Dependent | Expected Pattern |
|--------|------------|-------------------|------------------|
| **IBD (HMM)** | Similarity (total IBD bp) | Yes | Decreasing overall, stable/improving for transmission |
| **IBS** | Similarity (proportion shared) | No | Stable across recombination |
| **Phylogenetic** | Similarity (max_dist - patristic) | No | Stable across recombination |

### Critical: Metric Orientation

All methods are converted to **similarity metrics** where:
- **Higher score** = **more related**
- **Lower score** = **less related**

This ensures AUPR/AUROC metrics are comparable across methods.

#### Phylogenetic Distance Conversion

Patristic distances are **inverted** to similarity:

```r
# Patristic distance: smaller = more related
max_dist <- max(patristic_distances)
similarity <- max_dist - patristic_distance

# Now: larger similarity = more related ✓
```

## Files

### Main Analysis Scripts

1. **`transmission_analysis_v2.R`** - Main analysis pipeline
   - Loads data (IBD, IBS, phylogenetic)
   - Evaluates all three tasks
   - Computes enrichment metrics
   - Generates all visualizations

2. **`diagnostic_validation.R`** - Pre-analysis validation
   - Checks file existence
   - Validates metric directions
   - Tests on one replicate before full run
   - **Run this first!**

### Input Data Structure

```
simulations/multiple_runs/
├── inferred/
│   ├── rep1/
│   │   ├── run1_rec1e09_chr1/
│   │   │   ├── true_ibd_summary.tsv    # Ground truth
│   │   │   ├── inferred_ibd_hmm.rds    # HMM-IBD predictions
│   │   │   └── ibs_matrix.rds          # IBS matrix
│   │   ├── run1_rec1e08_chr1/
│   │   └── ...
│   ├── rep2/
│   └── ...
└── phylo_results/
    ├── rep1/
    │   ├── *_rec1e09_*modelfinder.treefile
    │   └── ...
    └── ...
```

### Output Files

```
simulations/multiple_runs/transmission_evaluation/
├── tables/
│   ├── results_all_replicates.csv       # Raw results per replicate
│   ├── results_aggregated.csv           # Mean ± SE across replicates
│   ├── detectability_by_relationship.csv
│   ├── detectability_summary.csv
│   ├── enrichment_all_replicates.csv
│   └── enrichment_summary.csv
└── figures/
    ├── 00_composite_main_figure.png     # Combined overview
    ├── 01_overall_vs_transmission_aupr.png
    ├── 02_transmission_discrimination_aupr.png
    ├── 03_detectability_by_relationship.png
    ├── 04_transmission_enrichment_topk.png
    └── 05_performance_heatmap.png
```

## Usage

### Step 1: Run Diagnostic Validation

```r
source("diagnostic_validation.R")
```

**Expected output:**
```
✅ Files: PASS
✅ IBD: PASS
✅ IBS: PASS
✅ Phylo: PASS

🎉 ALL CHECKS PASSED! 🎉
```

If any checks fail, fix the issues before proceeding.

### Step 2: Configure Paths

Edit `transmission_analysis_v2.R` configuration:

```r
CONFIG <- list(
  ROOT_DIR = "simulations/multiple_runs",  # Your base directory
  INFERRED_SUBPATH = "inferred",
  PHYLO_SUBPATH = "phylo_results",
  OUTDIR = "simulations/multiple_runs/transmission_evaluation",
  
  REC_RATES = c("1e-09", "1e-08", "1e-07", "1e-06"),
  GENOME_BP = 640851,  # Your genome size
  
  GEN_TRANSMISSION = 5,      # Transmission threshold
  GEN_RECENT_ANCESTRY = 15,
  GEN_DISTANT_ANCESTRY = 25,
  
  ALPHA_DETECT = 0.01,  # False positive rate for detectability
  SAVE_PLOTS = TRUE,
  SAVE_TABLES = TRUE
)
```

### Step 3: Run Main Analysis

```r
source("transmission_analysis_v2.R")

# Or run interactively:
analysis_results <- run_transmission_analysis()

# Access results:
results_dt <- analysis_results$results
agg_results <- analysis_results$agg_results
enrichment <- analysis_results$enrichment_summary
```

## Interpreting Results

### Key Metrics

#### AUPR (Area Under Precision-Recall Curve)
- Measures classification performance
- Range: [0, 1], higher is better
- >0.5 indicates better than random
- Compare across recombination rates for each task

#### Enrichment Factor
- (Observed % transmission in top K) / (Baseline % transmission)
- \>1 indicates enrichment
- \>1 and increasing with recombination → supports H3

#### Detectability
- Proportion of pairs with IBD ≥ threshold
- Threshold set at baseline rate (1% FPR)
- Track separately for each relationship type

### Expected Patterns

If the hypothesis is correct, you should see:

| Metric | Low Recombination | High Recombination |
|--------|------------------|-------------------|
| **AUPR (Overall)** | High (~0.85) | Low (~0.50) |
| **AUPR (Transmission)** | Medium (~0.75) | Stable or Higher (~0.80) |
| **AUPR (Discrimination)** | Low (~0.60) | Higher (~0.80) |
| **Enrichment (top 50)** | Low (~1.5x) | High (~3-5x) |
| **Detectability (transmission)** | High (~95%) | Medium (~70%) |
| **Detectability (distant)** | High (~90%) | Low (~20%) |

### Biological Interpretation

**If IBD transmission detection improves with recombination:**

1. **Practical implication**: In high-transmission settings (where recombination is high due to many infection cycles), IBD methods may be *more* useful for identifying direct transmission chains, not less.

2. **Methodological insight**: The "failure" of IBD at high recombination is actually a *feature* for transmission detection - it preferentially loses signal from distant ancestry while preserving recent transmission.

3. **Study design**: Transmission studies should consider using higher recombination as a natural filter to enrich for recent transmission events.

## Troubleshooting

### Common Issues

#### "AUPR < 0.5 for phylogenetic method"
- **Cause**: Distance not converted to similarity
- **Fix**: Check `load_phylo()` function has: `score := max_dist - distance`

#### "Negative Spearman correlation"
- **Cause**: Metric inverted (distance vs. similarity mismatch)
- **Fix**: Run diagnostic script to identify which metric

#### "No phylogenetic data found"
- **Cause**: Tree files not matching pattern or in wrong directory
- **Fix**: Check `PHYLO_SUBPATH` and file naming pattern

#### "Very low enrichment factors"
- **Cause**: Not enough transmission pairs, or threshold too low
- **Fix**: Adjust `GEN_TRANSMISSION` or check ground truth distribution

### Diagnostic Checks

The diagnostic script validates:

1. **File existence**: All required files are present
2. **Metric direction**: All scores higher for related pairs
3. **AUPR > 0.5**: All methods better than random
4. **Positive correlation**: All methods correlate positively with IBD

Run diagnostic after any changes to data loading functions.

## Differences from Original Script

### What Changed

1. **Phylogenetic handling**: Proper distance → similarity conversion
2. **Multiple ground truths**: Transmission, recent ancestry, distant ancestry
3. **Multiple tasks**: Overall, transmission detection, discrimination
4. **Enrichment analysis**: Transmission enrichment in top K pairs
5. **Stratified detectability**: Separate for each relationship type
6. **Clearer organization**: Modular functions, extensive comments

### What Stayed the Same

- Data loading logic (file paths, formats)
- AUPR/AUROC computation (PRROC, pROC packages)
- Replicate aggregation (mean ± SE)
- Overall structure (load → evaluate → visualize)

## References

### Key Concepts

- **Patristic distance**: Sum of branch lengths between tips in phylogenetic tree
- **Cophenetic distance**: Alternative name for patristic distance from `ape::cophenetic.phylo()`
- **IBD segment**: Continuous genomic region inherited from common ancestor without recombination
- **HMM-IBD**: Hidden Markov Model for detecting IBD segments

### Related Methods

This analysis framework can be adapted for:
- Other IBD detection methods (iLASH, GERMLINE, RaPID)
- Other genetic distances (FST, kinship, pi)
- Other recombination rate proxies (genetic map distance)
- Other infectious disease systems (malaria, TB, etc.)

## Citation

If you use this analysis framework, please cite:

```
[Your paper citation here]
```

## Contact

For questions or issues:
- [Your contact information]
- GitHub issues: [repository link]

---

**Last updated**: 2025-02-08  
**Version**: 2.0
