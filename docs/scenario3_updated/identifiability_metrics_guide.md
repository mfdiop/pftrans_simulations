# QUICK REFERENCE: Identifiability Metrics for Transmission Inference

## RECOMMENDED DEFINITION

**"Identifiable" = ROC-AUC ≥ 0.80 AND Sensitivity ≥ 0.60 at 90% specificity**

---

## 1. PRIMARY METRIC: ROC-AUC

### Calculation
```python
from sklearn.metrics import roc_auc_score, roc_curve

# For each scenario + method combination:
y_true = [1 if pair represents direct transmission, else 0]
y_score = [genetic distance metric for each pair]

# Calculate AUC
auc = roc_auc_score(y_true, y_score)

# Find sensitivity at 90% specificity
fpr, tpr, thresholds = roc_curve(y_true, y_score)
idx_90spec = np.argmin(np.abs(fpr - 0.10))  # FPR = 1 - specificity
sensitivity_at_90spec = tpr[idx_90spec]
```

### Identifiability Thresholds

| AUC Range | Sens @ 90% Spec | Verdict | Confidence |
|-----------|-----------------|---------|------------|
| **≥0.90** | ≥0.70 | **IDENTIFIABLE** | High |
| **0.80-0.89** | ≥0.60 | **IDENTIFIABLE** | Acceptable |
| **0.70-0.79** | 0.40-0.59 | **MARGINAL** | Low - context dependent |
| **0.60-0.69** | <0.40 | **NOT IDENTIFIABLE** | Poor discrimination |
| **<0.60** | <0.30 | **FAILURE** | Worse than informed guess |

### Interpretation
- **AUC = 0.90**: Perfect separation of transmission pairs from unrelated pairs
- **AUC = 0.80**: Good discrimination; transmission pairs have systematically higher/lower distance
- **AUC = 0.70**: Fair discrimination; signal exists but noisy
- **AUC = 0.50**: Random performance; no signal

---

## 2. SECONDARY METRICS (For Sensitivity Analysis)

### A. Traditional Classification Metrics

If you choose a fixed distance threshold:

```python
# Set threshold (e.g., IBD >0.5 = transmission)
threshold = 0.5
y_pred = [1 if distance > threshold else 0]

# Calculate metrics
from sklearn.metrics import confusion_matrix, classification_report

tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

sensitivity = tp / (tp + fn)  # True positive rate
specificity = tn / (tn + fp)  # True negative rate
PPV = tp / (tp + fp)          # Positive predictive value
F1 = 2 * (PPV * sensitivity) / (PPV + sensitivity)
```

**Identifiable when**:
- Sensitivity ≥ 0.70 AND Specificity ≥ 0.80 AND PPV ≥ 0.60

---

### B. Clustering Metrics

If you group samples into transmission clusters:

```python
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

true_clusters = [ground_truth_cluster_id for each sample]
inferred_clusters = [method_cluster_id for each sample]

ARI = adjusted_rand_score(true_clusters, inferred_clusters)
NMI = normalized_mutual_info_score(true_clusters, inferred_clusters)
```

**Identifiable when**:
- ARI ≥ 0.60 OR NMI ≥ 0.65

**Interpretation**:
- ARI = 1.0: Perfect clustering agreement
- ARI = 0.0: Random clustering
- ARI < 0.0: Worse than random

---

## 3. GROUND TRUTH DEFINITION

### Pairwise Relationship Labels

```python
# For each pair of sampled individuals (i, j):
relationships = {
    'direct_transmission': 1,      # j directly infected by i (1 step)
    'indirect_2steps': 0,          # 2 transmission steps between i and j
    'indirect_3+steps': 0,         # ≥3 steps (same outbreak cluster)
    'same_population': 0,          # Same deme, distant in transmission chain
    'different_population': 0,     # From different populations
    'unrelated': 0                 # No shared recent ancestor
}

# For binary classification (recommended):
y_true = {
    'direct_transmission': 1,      # Positive class
    'all_others': 0               # Negative class
}
```

**Alternative** (if testing cluster detection):
```python
# Multi-level relationships:
y_true = {
    'direct_transmission': 1,
    '2-3_steps': 2,
    'distant_cluster': 3,
    'unrelated': 4
}
# Then use ARI/NMI instead of AUC
```

---

## 4. METHOD COMPARISON FRAMEWORK (RQ4)

### Detecting Methodological vs. Fundamental Failure

```python
# For each scenario:
auc_ibd = calculate_auc(ground_truth, ibd_distances)
auc_ibs = calculate_auc(ground_truth, ibs_distances)
auc_phylo = calculate_auc(ground_truth, phylo_distances)

# Method performance range
auc_range = max(auc_ibd, auc_ibs, auc_phylo) - min(auc_ibd, auc_ibs, auc_phylo)

# Classification:
if max(auc_ibd, auc_ibs, auc_phylo) < 0.70:
    if auc_range < 0.10:
        failure_type = "FUNDAMENTAL (all methods fail equally)"
    else:
        failure_type = "MIXED (some methods work, identifiability possible)"
        best_method = identify_best(auc_ibd, auc_ibs, auc_phylo)
        
elif min(auc_ibd, auc_ibs, auc_phylo) < 0.70:
    if auc_range >= 0.15:
        failure_type = "METHODOLOGICAL (method choice matters)"
        failed_methods = [m for m in methods if auc[m] < 0.70]
        working_methods = [m for m in methods if auc[m] >= 0.80]
```

### Interpretation Matrix

| Max AUC | Min AUC | AUC Range | Interpretation |
|---------|---------|-----------|----------------|
| >0.80 | >0.80 | <0.10 | All methods work; identifiable |
| >0.80 | <0.70 | >0.15 | **Methodological failure**; method choice critical |
| <0.70 | <0.70 | <0.10 | **Fundamental failure**; insufficient signal |
| 0.70-0.80 | 0.70-0.80 | <0.10 | Marginal; context-dependent |

---

## 5. REGRESSION ANALYSIS (RQ2)

### Quantifying Parameter Effects

```python
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

# Compile results across all scenarios
results = pd.DataFrame({
    'rec_rate': [...],           # Log-transform: log10(rec_rate)
    'bottleneck': [...],
    'est': [...],
    'migration': [...],          # Log-transform: log10(migration_rate)
    'sampling': [...],
    'auc_ibd': [...],
    'auc_ibs': [...],
    'auc_phylo': [...]
})

# Standardize predictors for effect size comparison
X = results[['rec_rate', 'bottleneck', 'est', 'migration', 'sampling']]
X_scaled = StandardScaler().fit_transform(np.log10(X + 1e-10))

# Fit for each method
for method in ['ibd', 'ibs', 'phylo']:
    y = results[f'auc_{method}']
    model = LinearRegression().fit(X_scaled, y)
    
    # Standardized coefficients = effect sizes
    print(f"\n{method.upper()} Effect Sizes:")
    for param, coef in zip(X.columns, model.coef_):
        print(f"  {param}: {coef:.3f}")
```

**Interpret standardized coefficients**:
- |coef| > 0.3: Strong effect
- |coef| 0.1-0.3: Moderate effect  
- |coef| < 0.1: Weak effect

**Expected pattern for RQ2**:
- IBD: Large negative coefficient for `rec_rate` (recombination degrades IBD)
- IBS: Large positive coefficient for `rec_rate` (works better clonally)
- Phylo: Large negative coefficient for `rec_rate` + `bottleneck`

---

## 6. MIGRATION THRESHOLD ANALYSIS (RQ3)

### Detecting Importation Events

```python
# For pairs from different populations:
imported_pairs = pairs[pairs['true_pop_i'] != pairs['true_pop_j']]
local_pairs = pairs[pairs['true_pop_i'] == pairs['true_pop_j']]

# Calculate discrimination ability
for migration_rate in [0.001, 0.01, 0.05]:
    scenario_data = results[results['migration'] == migration_rate]
    
    # Can we distinguish imported from local?
    y_true = scenario_data['is_imported']
    y_score = scenario_data['genetic_distance']
    
    auc_import = roc_auc_score(y_true, y_score)
    
    # Threshold: accuracy <70% = can't distinguish
    if auc_import < 0.70:
        print(f"Migration {migration_rate}: CANNOT distinguish importation")
    else:
        print(f"Migration {migration_rate}: Can detect imports (AUC={auc_import:.2f})")
```

---

## 7. RECOMMENDED FIGURES FOR THESIS

### Figure 1: Method Performance Heatmap
```
          Recombination Rate
          1e-9  1e-8  1e-7  1e-6
IBD       0.68  0.75  0.85  0.90
IBS       0.88  0.82  0.72  0.65
Phylo     0.85  0.78  0.68  0.60

Color scale: Red (<0.70) → Yellow (0.70-0.80) → Green (>0.80)
```

### Figure 2: ROC Curves by Scenario
- Panel per recombination rate
- Overlay IBD, IBS, phylo curves
- Shows where methods diverge

### Figure 3: Identifiability Phase Diagram
- X-axis: Recombination rate (log scale)
- Y-axis: Sampling proportion
- Color: AUC (best method)
- Contour line at AUC = 0.80 (identifiability threshold)

### Figure 4: Parameter Effect Sizes
- Forest plot of standardized regression coefficients
- Separate panel for IBD, IBS, phylo
- Shows which parameters matter most for each method

---

## 8. VALIDATION CHECKLIST

Before analyzing full dataset:

- [ ] All three methods (IBD, IBS, phylo) use **identical pairs** from ground truth
- [ ] Ground truth labels verified on 2-3 test scenarios
- [ ] Sanity check: AUC >0.85 for "easy" scenario (low recomb, high sampling)
- [ ] Sanity check: AUC <0.60 for "hard" scenario (high recomb, sparse sampling)
- [ ] Hypothesis test: IBD < IBS for clonal scenarios (recomb = 1e-9)
- [ ] Hypothesis test: IBD > IBS for recombining scenarios (recomb = 1e-6)

If any sanity check fails → debug before scaling up

---

## 9. SUMMARY: ANSWERING "IS IT IDENTIFIABLE?"

**For each scenario**, calculate:

1. **AUC** for each method (primary metric)
2. **Sensitivity at 90% specificity** (practical utility)
3. **Best method** (IBD, IBS, or phylo)

**Then classify**:

```python
def classify_identifiability(auc, sens_at_90spec):
    if auc >= 0.90 and sens_at_90spec >= 0.70:
        return "IDENTIFIABLE - High confidence"
    elif auc >= 0.80 and sens_at_90spec >= 0.60:
        return "IDENTIFIABLE - Acceptable"
    elif auc >= 0.70:
        return "MARGINAL - Context dependent"
    else:
        return "NOT IDENTIFIABLE"
```

**Report**:
- % of scenarios identifiable overall
- % identifiable per recombination rate
- % identifiable per sampling proportion
- Identify parameter combinations where NO method works (fundamental limit)
- Identify where method choice matters (methodological)

---

## FINAL RECOMMENDATION

**Primary metric**: ROC-AUC with threshold ≥0.80  
**Validation metric**: Sensitivity ≥0.60 at 90% specificity  
**Comparison framework**: AUC range across methods to detect methodological vs. fundamental failure

This gives you:
- **Interpretable** (AUC is standard in epidemiology)
- **Comparable** (all three methods produce distance/similarity scores)
- **Defensible** (published thresholds for "good" discrimination)
- **Thesis-ready** (answers all four RQs with quantitative metrics)
