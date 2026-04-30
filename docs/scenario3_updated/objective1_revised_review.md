# REVISED Critical Review: PhD Objective 1 with Inference Methods

## EXECUTIVE SUMMARY

**Status**: Design is now **strong and complete** with IBD, IBS, and phylogenetic distance implemented separately.

**Core thesis**: "IBD is not universally superior for malaria transmission inference" - this is publishable if you can demonstrate the conditions under which it fails.

**Main tasks remaining**: 
1. Define quantitative identifiability thresholds
2. Justify recombination-as-MOI proxy
3. Ensure inference methods use identical ground truth for comparison

---

## 1. UPDATED ASSESSMENT: INFERENCE METHODS

### 1.1 Your Three-Method Approach

**IBD (Identity-by-Descent)**
- Detects genomic segments inherited from recent common ancestor
- Gold standard for malaria due to high recombination breaking up background LD
- **Hypothesis**: Will fail when recombination is too low (clonal) or too high (signal fragmented)

**IBS (Identity-by-State)** 
- Measures raw sequence similarity
- Simpler, doesn't account for recombination structure
- **Hypothesis**: May outperform IBD in low-recombination settings (clonal transmission)

**Phylogenetic Distance**
- Genetic distance on maximum likelihood tree
- Standard for bacterial outbreaks
- **Hypothesis**: Will fail in high-recombination scenarios (conflicting phylogenies)

### 1.2 Why This Design is Strong

**This is exactly right for your RQ4**: "methodological vs. fundamental failure"

You're comparing methods with different assumptions:
- **IBD**: Assumes recombination breaks up haplotypes (true for high-transmission malaria)
- **IBS**: Assumes recent transmission = high similarity (true for clonal pathogens)  
- **Phylogenetic**: Assumes vertical inheritance (violated by recombination)

**Predicted outcome** (testable hypothesis):
- Low recombination (1e-9): IBS and phylogenetic distance outperform IBD
- High recombination (1e-6): IBD outperforms others
- Intermediate: Methods converge in performance

**This would be a strong first-author paper**: "Rethinking IBD for malaria transmission inference: performance depends on population recombination rate"

---

## 2. DEFINING "IDENTIFIABLE"

### 2.1 Core Concept

**Identifiable** = the genetic data contains sufficient information to distinguish true transmission links from background relatedness with acceptable accuracy.

### 2.2 Operational Definitions (Choose Based on Your Inference Output)

#### **Option A: Pairwise Transmission Link Classification**

If your methods classify pairs as "direct transmission" vs. "unrelated":

**Identifiable when**:
- **Sensitivity (Recall) ≥ 0.70**: Detect 70%+ of true transmission pairs
- **Specificity ≥ 0.80**: Avoid 80%+ of false positive pairs
- **Positive Predictive Value (PPV) ≥ 0.60**: 60%+ of predicted links are true

**Failure threshold**:
- Sensitivity <0.50 OR Specificity <0.70 = "not identifiable"
- This means the method performs worse than a moderately informed guess

**Metrics to calculate**:
```
For each scenario:
  True Positives (TP) = direct transmission pairs correctly identified
  False Positives (FP) = unrelated pairs misclassified as transmission
  True Negatives (TN) = unrelated pairs correctly classified
  False Negatives (FN) = transmission pairs missed
  
  Sensitivity = TP / (TP + FN)
  Specificity = TN / (TN + FP)  
  PPV = TP / (TP + FP)
  F1-score = 2 × (PPV × Sensitivity) / (PPV + Sensitivity)
  
  IDENTIFIABLE if: Sensitivity ≥ 0.70 AND Specificity ≥ 0.80
```

---

#### **Option B: Transmission Cluster Reconstruction**

If your methods group samples into transmission clusters:

**Identifiable when**:
- **Adjusted Rand Index (ARI) ≥ 0.60**: Clustering matches true transmission chains
- **Normalized Mutual Information (NMI) ≥ 0.65**: Shared information between inferred and true clusters
- **Cluster Purity ≥ 0.75**: 75%+ of samples in each inferred cluster share true transmission source

**Failure threshold**:
- ARI <0.40 OR NMI <0.50 = "not identifiable"  

**Metrics to calculate**:
```python
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# For each scenario:
true_clusters = [cluster_id for each sample based on ground truth]
inferred_clusters = [cluster_id from your IBD/IBS/phylo method]

ARI = adjusted_rand_score(true_clusters, inferred_clusters)
NMI = normalized_mutual_info_score(true_clusters, inferred_clusters)

# Cluster purity
for each inferred_cluster:
    purity = (count of most common true_cluster) / cluster_size
    
IDENTIFIABLE if: ARI ≥ 0.60 AND NMI ≥ 0.65
```

---

#### **Option C: Distance-Based Threshold Performance**

If your methods use genetic distance thresholds to define transmission:

**Identifiable when**:
- **ROC-AUC ≥ 0.85**: Strong discriminative ability between transmission pairs vs. unrelated
- **Optimal threshold is stable**: Best cutoff doesn't vary wildly across replicates

**Failure threshold**:
- ROC-AUC <0.70 = "barely better than random"
- Optimal threshold CV >30% = "unstable, not identifiable"

**Metrics to calculate**:
```python
from sklearn.metrics import roc_curve, roc_auc_score

# For each scenario:
y_true = [1 if pair is direct transmission, 0 otherwise]
y_score = [genetic distance for each pair]  # From IBD, IBS, or phylo

auc = roc_auc_score(y_true, y_score)
fpr, tpr, thresholds = roc_curve(y_true, y_score)

# Find optimal threshold (max Youden's J)
J = tpr - fpr
optimal_threshold = thresholds[np.argmax(J)]

IDENTIFIABLE if: AUC ≥ 0.85
```

---

### 2.3 RECOMMENDED APPROACH FOR YOUR STUDY

**Use Option C (Distance-Based + ROC)** because:
1. Works for all three methods (IBD, IBS, phylogenetic distance)
2. Allows direct comparison of method performance
3. ROC-AUC is interpretable and standard in epidemiology

**Define identifiability thresholds**:

| AUC Range | Interpretation | Identifiable? |
|-----------|----------------|---------------|
| **0.90-1.0** | Excellent discrimination | **YES - High confidence** |
| **0.80-0.89** | Good discrimination | **YES - Acceptable** |
| **0.70-0.79** | Fair discrimination | **MARGINAL - Context-dependent** |
| **0.50-0.69** | Poor discrimination | **NO - Not identifiable** |
| **<0.50** | Worse than random | **NO - Complete failure** |

**Additional criterion**: Require **sensitivity ≥ 0.60 at 90% specificity**
- This ensures practical utility (detect most true cases without false alarms)
- Calculate from ROC curve: find sensitivity where specificity = 0.90

---

### 2.4 Ground Truth Definition

**Critical**: Ensure your simulations export pairwise relationships:

```python
# For each pair of sampled individuals (i, j):
ground_truth = {
    'sample_i': i,
    'sample_j': j,
    'relationship': 'direct_transmission',  # or 'same_cluster', 'unrelated'
    'transmission_distance': 1,  # Steps in transmission chain (1 = direct)
    'generation_gap': 3,  # Time between infections
    'shared_population': True  # Same deme or different
}
```

**Relationship categories**:
- `direct_transmission`: j infected by i (1 step)
- `2nd_degree`: j infected by someone i infected (2 steps)  
- `same_cluster`: j and i in same transmission chain (≤5 steps)
- `same_population`: j and i from same deme but distant transmission
- `imported`: j and i from different populations
- `unrelated`: No recent common ancestor in simulation

---

## 3. RECOMBINATION AS MOI PROXY: JUSTIFICATION

### 3.1 Your Approach

**What you did**: Vary recombination rate (1e-9 to 1e-6) as a proxy for MOI (multiplicity of infection)

**Why this is reasonable**:
1. **Biological correlation**: High-transmission areas have both:
   - High MOI (chronic polyclonal infections)
   - More opportunities for recombination (multiple genotypes per mosquito)

2. **Genetic consequence is similar**:
   - High MOI → diverse haplotypes within hosts
   - High recombination → mosaic genomes in offspring
   - **Both erode IBD signals** by shuffling haplotype blocks

### 3.2 How to Justify in Your Thesis

**Section in Methods**:

> "Due to technical constraints in implementing superinfection dynamics in SLiM, we used recombination rate as a proxy for transmission intensity and multiplicity of infection (MOI). This approach is biologically justified because:
>
> 1. **Epidemiological correlation**: High-transmission settings exhibit both elevated MOI (chronic polyclonal infections) and increased recombination opportunities (mosquitoes ingesting multiple genotypes) (References: Nkhoma et al. 2020, Chang et al. 2019).
>
> 2. **Equivalent genetic signatures**: From an inference perspective, high MOI and high recombination produce similar genetic patterns—fragmented haplotype blocks and elevated within-host diversity—both of which challenge IBD-based transmission reconstruction (References: Taylor et al. 2019, Schaffner et al. 2018).
>
> 3. **Parameter calibration**: Our recombination rates (1e-9 to 1e-6) span values estimated from field isolates in low-transmission (predominantly clonal, r~1e-9) to high-transmission (frequent outcrossing, r~1e-6) settings (References: Miles et al. 2016 PMID:27428910, MalariaGEN).
>
> While this simplification does not capture all dynamics of superinfection (e.g., within-host competition, immune selection), it isolates the key challenge for genomic inference: the erosion of linkage disequilibrium through recombination."

### 3.3 Limitations to Acknowledge

**In Discussion**:
- Your model doesn't capture **temporal accumulation** of MOI (reinfection over weeks)
- Doesn't model **within-host selection** (some clones outcompete others)
- Bottleneck size (1-20) partially captures MOI but not chronic polyclonal infections (MOI 5-10 stable for months)

**Mitigation**: Phrase findings as "recombination rate-dependent" rather than "MOI-dependent" when presenting results. This is more precise.

---

## 4. REVISED RESEARCH QUESTIONS (WITH METRICS)

### **RQ1 (quantified)**:
"Under what combinations of recombination rate, transmission bottleneck, mutation accumulation, and sampling density can IBD, IBS, and phylogenetic distance methods achieve identifiable transmission inference (ROC-AUC ≥0.80, sensitivity ≥0.60 at 90% specificity)?"

**Analysis plan**:
- Calculate AUC for each method × scenario × replicate
- Plot AUC as function of: recombination rate, sampling proportion, migration rate
- Identify "phase transitions" where AUC drops below 0.80

---

### **RQ2 (quantified)**:
"Which biological processes (recombination, bottleneck size, mutation accumulation) most strongly degrade inference accuracy, and do different methods (IBD vs. IBS vs. phylogenetic) show differential sensitivity to these processes?"

**Analysis plan**:
- Fit regression model: `AUC ~ recombination + bottleneck + est + migration + sampling`
- Compare standardized coefficients to rank effect sizes
- Test for **method × parameter interactions** (e.g., does IBD degrade faster than IBS as recombination increases?)

**Hypothesis**: IBD will show strongest negative coefficient for recombination rate; phylogenetic distance for bottleneck size.

---

### **RQ3 (quantified)**:
"At what migration rate does spatial structure become insufficient to distinguish local transmission from importation (classification accuracy <70% for imported cases)?"

**Analysis plan**:
- For each scenario, classify pairs as "same population" vs. "imported"
- Calculate sensitivity for detecting imported cases using genetic distance
- Plot accuracy vs. migration rate
- Find critical threshold where accuracy <0.70

**Hypothesis**: At migration >0.01 (1% per generation), populations are too homogeneous to infer importation.

---

### **RQ4 (quantified)**:
"When transmission inference fails (AUC <0.70), is failure uniform across methods (fundamental information loss) or method-specific (methodological limitation)?"

**Analysis plan**:
```
For each failed scenario (AUC <0.70):
  Calculate: range_AUC = max(AUC_IBD, AUC_IBS, AUC_phylo) - min(...)
  
  If range_AUC < 0.10:
    → "Fundamental failure" (all methods fail equally)
  
  If range_AUC ≥ 0.15:
    → "Methodological failure" (one method outperforms others)
    → Identify which method works best under those conditions
```

**Hypothesis**: 
- Low recombination failures are **methodological** (IBS/phylo work, IBD fails)
- High recombination failures are **fundamental** (all methods fail)
- Sparse sampling failures are **fundamental** (insufficient data)

---

## 5. CRITICAL NEXT STEPS

### 5.1 Validate Inference Scripts Against Known Truth

Before running 900 scenarios, test on 2-3 replicates:

**Sanity check**:
```python
# For a simple scenario (low recombination, high sampling):
assert AUC_all_methods > 0.85, "Methods should work on easy case"

# For impossible scenario (high recombination, sparse sampling):  
assert AUC_all_methods < 0.60, "Methods should fail on hard case"

# For your hypothesis:
# At recombination = 1e-9 (clonal):
assert AUC_IBS > AUC_IBD, "IBS should beat IBD in clonal settings"

# At recombination = 1e-6 (high outcrossing):
assert AUC_IBD > AUC_IBS, "IBD should beat IBS with recombination"
```

If these don't hold, **debug before scaling up**.

---

### 5.2 Ensure Consistent Ground Truth

**Critical question**: Are your IBD, IBS, and phylogenetic scripts using the **identical pairwise labels** from simulation truth?

**Required for fair comparison**:
```python
# All three methods must evaluate against same pairs:
ground_truth.csv:
  sample_i, sample_j, is_direct_transmission
  
ibd_results.csv:
  sample_i, sample_j, ibd_proportion, ibd_distance
  
ibs_results.csv:
  sample_i, sample_j, ibs_similarity
  
phylo_results.csv:
  sample_i, sample_j, patristic_distance
  
# Then merge and calculate AUC for each method
```

**Common mistake**: Each method uses different sample subsets (e.g., phylogenetic drops low-quality samples). Ensure same n for all methods.

---

### 5.3 Parameter Range Additions (Still Recommended)

**A. Add high mutation scenarios**:
- `est = 1.0, 2.0` to test if rapid evolution masks recent transmission

**B. Add sparse sampling**:
- `sampling_proportion = 0.01, 0.02` to probe real-world surveillance limits

**Justification**: These are the regimes where you'd expect fundamental failure (not enough data) vs. methodological differences disappearing.

---

## 6. EXPECTED RESULTS (TESTABLE PREDICTIONS)

Based on your design, here's what I predict you'll find:

### 6.1 Main Effects (RQ1-2)

**Recombination rate** (strongest predictor):
- 1e-9 (clonal): IBS/phylo AUC ~0.90, IBD AUC ~0.70
- 1e-6 (high recomb): IBD AUC ~0.85, IBS/phylo AUC ~0.65
- **Crossover point**: ~1e-8 where all methods converge

**Sampling proportion**:
- 5%: All methods AUC ~0.75 (marginal)
- 20%: All methods AUC >0.85 (identifiable)
- **Critical threshold**: ~10% for reliable inference

**Migration rate** (RQ3):
- <0.01: Can distinguish populations (importation accuracy >80%)
- >0.03: Populations homogeneous (importation accuracy ~50%, random guess)

**Bottleneck size**:
- Tight (1): Phylogenetic distance works well (clonal descent)
- Loose (20): IBD works better (captures recombinant haplotypes)

### 6.2 Interactions (RQ4)

**Key interaction**: `recombination × method`
- Low recombination: Method choice matters (IBS best)
- High recombination: Method choice matters (IBD best)
- **Sparse sampling**: Method choice doesn't matter (all fail)

This pattern would support: **"Failure is methodological at extremes, fundamental at low sampling"**

---

## 7. PUBLISHABILITY ASSESSMENT

### 7.1 Potential Outcomes → Papers

**If you find IBD fails in low-transmission settings**:
- **Title**: "Rethinking identity-by-descent for malaria transmission inference in elimination settings"
- **Journal**: PLOS Computational Biology, Genetics, or MBE
- **Impact**: Changes how people use IBD for low-transmission malaria

**If you find no method works below X% sampling**:
- **Title**: "Genomic surveillance thresholds for infectious disease transmission inference"  
- **Journal**: PLOS Pathogens, eLife
- **Impact**: Informs surveillance design for elimination programs

**If you find recombination is the dominant factor**:
- **Title**: "The identifiability crisis in transmission genomics: recombination erodes inference across methods"
- **Journal**: Nature Communications, Molecular Biology & Evolution
- **Impact**: Broad implications beyond malaria

### 7.2 Thesis Contribution

**This objective alone could yield**:
- **1 first-author paper** (method comparison + conditions for identifiability)
- **1 methods paper** (if you develop a new composite metric combining IBD+IBS)
- **Framework** for Objectives 2-3 (you've established the "ground truth" testing paradigm)

---

## 8. FINAL VERDICT (REVISED)

### Objective Statement: A- 
With the revised understanding of your methods, the objective is strong. Minor suggestion:

**Add method specificity**:
"...quantified by comparing identity-by-descent, identity-by-state, and phylogenetic distance methods against simulated transmission networks."

### Research Questions: A
With the quantitative metrics defined above (ROC-AUC thresholds, sensitivity/specificity), these are now testable and defensible.

### Simulation Design: A
- 45 scenarios × 20 replicates = robust
- Recombination-as-MOI is a reasonable proxy (with caveats acknowledged)
- Would benefit from expanded est and sampling ranges

### Inference Framework: A
- Three complementary methods with different assumptions
- Tests a falsifiable hypothesis (IBD superiority)
- Missing: just need to validate they use identical ground truth

---

## 9. IMMEDIATE ACTION ITEMS

**Before running full design**:

1. ✅ **Define identifiability**: Use ROC-AUC ≥0.80 + sensitivity ≥0.60 @ 90% specificity

2. ✅ **Test hypothesis on 3 pilot scenarios**:
   - Low recomb (1e-9): Expect IBS > IBD
   - High recomb (1e-6): Expect IBD > IBS  
   - Sparse sampling (5%): Expect all methods fail

3. ✅ **Verify ground truth consistency**: Same pairwise labels for IBD, IBS, phylo

4. ⚠️ **Add parameter ranges**: est=1.0, 2.0 and sampling=1%, 2%

5. ✅ **Document recombination-MOI justification**: Add to methods section

6. ✅ **Pre-register analysis plan**: Commit to RQ1-4 metrics before seeing results

**Timeline suggestion**:
- Week 1: Pilot validation (items 1-3)
- Week 2: Expand design if needed (item 4)
- Weeks 3-8: Run full simulation + inference pipeline
- Weeks 9-12: Analysis and write-up

---

## CONCLUSION

**Your design is thesis-ready with implementation of the identifiability metrics above.**

The core hypothesis—that IBD is not universally superior for malaria—is scientifically sound and publishable. Your three-method comparison allows you to distinguish methodological from fundamental failures (RQ4), which is the hardest but most valuable contribution.

**Key strengths**:
1. Falsifiable hypothesis about IBD performance
2. Well-justified recombination rate proxy for MOI
3. Comprehensive parameter space (45 scenarios)
4. Three complementary inference methods

**Main risk**: If all three methods always agree (always work or always fail together), RQ4 becomes less interesting. But I predict you'll find method-specific failures at the extremes (low/high recombination), which is the publishable finding.

**You're in excellent shape. Focus on validating the pilot scenarios before scaling up.**
