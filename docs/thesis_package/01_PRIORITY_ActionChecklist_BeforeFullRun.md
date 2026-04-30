# IMMEDIATE ACTION CHECKLIST
## Before Running Full 900 Simulations

---

## ✅ WHAT YOU ALREADY HAVE (EXCELLENT)

1. ✅ **AUROC calculation** - Line 478 of evaluate_recombination_effects_v1_1.r
2. ✅ **Three inference methods** - IBD, IBS, Phylogenetic distance
3. ✅ **Ground truth extraction** - From SLiM pedigrees
4. ✅ **Comprehensive evaluation** - AUPR, Spearman, Brier, Precision@K
5. ✅ **Systematic replication** - Multiple replicates per scenario
6. ✅ **Proper canonicalization** - Consistent pair ordering (id1 <= id2)

**You're in great shape!** Just need a few clarifications and additions.

---

## ⚠️ CRITICAL: VERIFY BEFORE FULL RUN (2-3 days)

### 1. Ground Truth Definition (HIGHEST PRIORITY)

**Current setting**: `gen_threshold = 25` (line 42, evaluate_recombination_effects_v1_1.r)

**Questions**:
- What does `gen_distance` measure? (Parent-child steps in pedigree?)
- Why 25 generations?
- Is this "direct transmission" (1-2 steps) or "transmission cluster" (5-10 steps)?

**Action needed**:
```r
# OPTION A: If studying direct transmission
gen_threshold = 1  # Change from 25 to 1

# OPTION B: If studying transmission clusters
gen_threshold = 5  # Or 3, 10 depending on definition

# OPTION C: Test sensitivity
gen_thresholds = c(1, 2, 3, 5, 10, 25)  # Run analysis with multiple thresholds
```

**Why this matters**:
- Your RQ asks about "transmission links"
- 25 generations ≈ 8-25 months (very distant relationships)
- For malaria elimination, direct transmission (1-2 steps) is more relevant

**Deliverable**: Document in methods section why you chose this threshold.

---

### 2. Add Sensitivity at 90% Specificity (30 minutes)

**Current**: Only AUROC calculated
**Needed**: Sensitivity at 90% specificity (validation metric from my recommendations)

**Where to add**: In `evaluate_method_on_table()` function, after line 478

```r
# AFTER line 481 (res$auroc <- ...), ADD:

if (!is.null(roc_obj)) {
  # Get all sensitivity/specificity coordinates
  coords <- pROC::coords(roc_obj, x = "all", 
                        ret = c("sensitivity", "specificity"))
  
  # Find sensitivity where specificity >= 0.90
  idx_90spec <- which(coords$specificity >= 0.90)
  
  if (length(idx_90spec) > 0) {
    res$sensitivity_at_90spec <- max(coords$sensitivity[idx_90spec])
  } else {
    res$sensitivity_at_90spec <- 0  # Can't achieve 90% specificity
  }
} else {
  res$sensitivity_at_90spec <- NA_real_
}
```

**Then store it** (line 762-774):
```r
row <- data.table(
  replicate = rep,
  method = m,
  rate = rate_numeric,
  rate_label = rate,
  aupr = eval_res$aupr,
  auroc = eval_res$auroc,
  sensitivity_at_90spec = eval_res$sensitivity_at_90spec,  # ADD THIS LINE
  spearman = eval_res$spearman,
  ...
)
```

**Why**: This is your identifiability validation metric:
- AUROC ≥ 0.80 AND Sensitivity ≥ 0.60 at 90% specificity = "Identifiable"

---

### 3. Verify IBS Direction (15 minutes)

**Check**: Line 198 of recombination_evaluation_pipeline_single_run.R

```r
long[, score := distance]  # Is "distance" actually similarity?
```

**Test**:
```r
# For one replicate, check IBS values:
ibs_mat <- readRDS("path/to/ibs_matrix.rds")

# Look at values:
summary(as.vector(ibs_mat))

# Expected for IBS SIMILARITY:
# - Range: 0 to 1
# - Related pairs: High values (>0.8)
# - Unrelated pairs: Low values (<0.5)

# If opposite (related = low, unrelated = high):
# Your matrix is DISTANCE, not similarity
# CHANGE line 198 to: long[, score := -distance]
```

**Why this matters**: If direction is wrong, IBS results will be completely inverted.

---

### 4. Check Pair Coverage (15 minutes)

**Add diagnostics** after line 724 in evaluate_recombination_effects_v1_1.r:

```r
# AFTER merging all methods (line 720), ADD:
message("  → Method coverage:")
for (m in names(method_tables)) {
  score_col <- paste0("score_", m)
  n_scored <- sum(!is.na(merged[[score_col]]))
  n_total <- nrow(merged)
  pct <- round(100 * n_scored / n_total, 1)
  message("    • ", m, ": ", n_scored, " / ", n_total, 
          " (", pct, "%)")
}

# Check if all methods score the same pairs
score_cols <- grep("^score_", names(merged), value = TRUE)
complete_pairs <- rowSums(!is.na(merged[, ..score_cols])) == length(score_cols)
n_complete <- sum(complete_pairs)
message("  → Pairs scored by ALL methods: ", n_complete, " / ", nrow(merged),
        " (", round(100 * n_complete / nrow(merged), 1), "%)")
```

**Why**: If IBD has 95% coverage but Phylo has 80%, you're comparing on different data.

**If coverage differs >10%**: Document and report as limitation, or subset to complete pairs.

---

## 📋 MEDIUM PRIORITY: ENHANCE ANALYSIS (1-2 days)

### 5. Test Hypothesis on Pilot Data

**Before running 900 simulations**, test on 3 scenarios:

**Scenario 1: Easy (Low recombination)**
- rec_rate = 1e-9
- sampling = 20%
- **Expected**: IBS > IBD (AUC_IBS should be higher)

**Scenario 2: Hard (High recombination)**
- rec_rate = 1e-6  
- sampling = 5%
- **Expected**: IBD > IBS (AUC_IBD should be higher)

**Scenario 3: Realistic**
- rec_rate = 1e-7
- sampling = 10%
- **Expected**: Methods differ (AUC range >0.15)

**Verification checklist**:
```r
# For each scenario, check:
- All methods AUROC >0.85 for Easy? ✓ / ✗
- All methods AUROC <0.60 for Hard? ✓ / ✗
- IBS > IBD for Low recombination? ✓ / ✗
- IBD > IBS for High recombination? ✓ / ✗
```

**If any fail** → Debug before running full design

---

### 6. Add Method Comparison Metrics (For RQ4)

**Current**: You calculate AUROC per method separately
**Needed**: Direct comparison to detect methodological vs. fundamental failure

**Add to final analysis**:
```r
# After all methods evaluated, compute:
method_comparison <- results[, .(
  auc_ibd = mean(auroc[method == "IBD"], na.rm = TRUE),
  auc_ibs = mean(auroc[method == "IBS"], na.rm = TRUE),
  auc_phylo = mean(auroc[method == "Phylo"], na.rm = TRUE)
), by = .(rate, replicate)]

# Calculate range
method_comparison[, auc_range := pmax(auc_ibd, auc_ibs, auc_phylo, na.rm = TRUE) - 
                                  pmin(auc_ibd, auc_ibs, auc_phylo, na.rm = TRUE)]

# Classify failure type (RQ4)
method_comparison[, failure_type := fcase(
  pmax(auc_ibd, auc_ibs, auc_phylo, na.rm = TRUE) < 0.70 & auc_range < 0.10,
    "FUNDAMENTAL (all methods fail)",
  pmax(auc_ibd, auc_ibs, auc_phylo, na.rm = TRUE) >= 0.80 & auc_range >= 0.15,
    "METHODOLOGICAL (method choice matters)",
  default = "MIXED"
)]

# Report
method_comparison[, .N, by = failure_type]
```

---

### 7. Document Ground Truth Extraction

**Create a methods note** explaining:

```markdown
# Ground Truth Extraction

## Pedigree Information
- Source: SLiM `keepPedigrees=T` output
- Metric: `gen_distance` = number of parent-child steps between samples
- Range: 0 (self) to MAX_GENERATIONS (unrelated)

## Labeling Scheme
- Direct transmission: gen_distance = 1 (parent directly infected child)
- 2-step transmission: gen_distance = 2 (grandparent-grandchild)
- Transmission cluster: gen_distance ≤ [THRESHOLD]
- Unrelated: gen_distance > [THRESHOLD]

## IBD Calculation
- Source: Tree sequence true IBD (not inferred)
- Metric: Sum of shared genomic segments (bp)
- Use: Detectability analysis only (not used for truth labels)

## Rationale for Threshold
[YOUR JUSTIFICATION HERE - e.g., "We use gen_distance ≤ 1 to define direct 
transmission because this represents parent-child infection within a single 
mosquito-human cycle, which is the biological definition of a transmission link 
in malaria epidemiology."]
```

---

## 🔧 OPTIONAL: NICE-TO-HAVE ENHANCEMENTS

### 8. Save ROC Curves for Plotting
```r
# In evaluate_method_on_table(), after calculating roc_obj:
if (!is.null(roc_obj)) {
  roc_coords <- pROC::coords(roc_obj, x = "all", 
                             ret = c("threshold", "sensitivity", "specificity"))
  # Save for later visualization
  res$roc_curve <- as.data.table(roc_coords)
}
```

### 9. Add Regression Analysis for RQ2
```r
library(lme4)

# Effect of recombination on each method
model_ibd <- lmer(auroc ~ log10(rate) + (1|replicate), 
                  data = results[method == "IBD"])
model_ibs <- lmer(auroc ~ log10(rate) + (1|replicate), 
                  data = results[method == "IBS"])

# Extract coefficients (effect sizes)
summary(model_ibd)  # Expect positive (higher rec → better IBD)
summary(model_ibs)  # Expect negative (higher rec → worse IBS)
```

---

## 📅 SUGGESTED TIMELINE

### Week 1: Verification (5 days)

**Day 1**: 
- [ ] Review ground truth definition in simulation code
- [ ] Decide on appropriate `gen_threshold` (1, 3, 5, 10, or 25?)
- [ ] Document rationale

**Day 2**:
- [ ] Add sensitivity at 90% specificity calculation
- [ ] Test on one replicate
- [ ] Verify output includes new metric

**Day 3**:
- [ ] Check IBS direction (similarity vs. distance)
- [ ] Add pair coverage diagnostics
- [ ] Run diagnostics on one replicate

**Day 4-5**:
- [ ] Run 3 pilot scenarios (easy, hard, realistic)
- [ ] Verify hypothesis:
  - Easy: IBS > IBD ✓
  - Hard: IBD > IBS ✓
  - Realistic: Methods differ ✓
- [ ] Check all AUROC values are reasonable (0.5-1.0)

### Week 2: Enhancement (3-5 days)

**Day 6-7**:
- [ ] Add method comparison framework
- [ ] Test on pilot data
- [ ] Verify failure type classification works

**Day 8-9**:
- [ ] Document ground truth extraction
- [ ] Write methods section text
- [ ] Add any optional enhancements

**Day 10**:
- [ ] Final review of all code
- [ ] Run one complete replicate end-to-end
- [ ] Verify all outputs as expected

### Week 3+: Full Run
- [ ] Launch full 900 simulations
- [ ] Monitor progress
- [ ] Begin analysis as results come in

---

## ✅ FINAL CHECKLIST BEFORE FULL RUN

Before submitting 900 simulation jobs:

- [ ] **Ground truth threshold justified or changed** (gen_threshold = ?)
- [ ] **Sensitivity at 90% specificity added** to evaluation
- [ ] **IBS direction verified** (high score = more related)
- [ ] **Pair coverage checked** (all methods ~same coverage)
- [ ] **Pilot validated** (3 scenarios show expected patterns)
- [ ] **Method comparison framework ready** for RQ4
- [ ] **All code documented** (comments, methods notes)
- [ ] **Output directory structure tested**

**If all checked** → **Proceed with confidence!**

---

## 🎯 SUCCESS CRITERIA

After full run, you should be able to answer:

**RQ1**: "Under [threshold] generations and [parameters], IBD/IBS/Phylo achieve 
identifiable inference (AUROC ≥0.80) in X% of scenarios."

**RQ2**: "Recombination rate has the strongest effect (β = X, p < 0.001), 
followed by sampling (β = Y), with IBD showing 3× larger effect than IBS."

**RQ3**: "At migration rate >0.01, importation detection accuracy drops below 
70%, indicating loss of spatial structure signal."

**RQ4**: "In Z% of failed scenarios (AUROC <0.70), failure is fundamental 
(all methods fail), but in W% failure is methodological (some methods work)."

**Hypothesis confirmed**: "IBS outperforms IBD in low-recombination scenarios 
(1e-9: IBS=0.88, IBD=0.68, p<0.001), while IBD outperforms in high-recombination 
scenarios (1e-6: IBD=0.90, IBS=0.65, p<0.001)."

---

## 🚀 YOU'RE ALMOST THERE!

Your code is **high quality**. With these minor verifications and additions, 
you'll have a robust, defensible analysis pipeline.

**Main strengths**:
- Already calculating AUROC ✅
- Comprehensive metrics ✅
- Proper replication ✅
- Good code structure ✅

**Main needs**:
- Verify/justify ground truth threshold ⚠️
- Add one missing metric (sensitivity at 90% spec) ⚠️
- Test on pilots before scaling ⚠️

**Estimated time to fully ready**: 5-10 days (including pilot validation)

Good luck!
