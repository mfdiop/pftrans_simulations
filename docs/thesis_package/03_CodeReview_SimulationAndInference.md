# CODE REVIEW: Simulation & Inference Pipeline
## Critical Assessment of Implementation for PhD Objective 1

---

## EXECUTIVE SUMMARY

**Overall Assessment: A- (Excellent with minor refinements needed)**

Your implementation is **strong and well-structured**. You have:
- ✅ Comprehensive simulation framework (SLiM + msprime)
- ✅ Three inference methods (IBD, IBS, phylogenetic)
- ✅ **AUROC already calculated** alongside AUPR
- ✅ Proper ground truth extraction
- ✅ Systematic evaluation across replicates

**Main gaps identified**:
1. ⚠️ Ground truth definition needs clarification (gen_threshold = 25?)
2. ⚠️ Methods may not be using identical pair sets (need verification)
3. ⚠️ Missing: Sensitivity at 90% specificity calculation
4. ✓ AUROC is already implemented (line 478 of evaluate script)

---

## 1. SIMULATION PIPELINE REVIEW

### 1.1 Python Script (01_run_simulation_replicates.py)

**STRENGTHS:**
- ✅ **Proper pedigree tracking** (`keepPedigrees=T` in SLiM script, line 108)
- ✅ **Ground truth extraction** with genealogical relationships
- ✅ **VCF output** for downstream analysis
- ✅ **Realistic Pf parameters** (chromosome lengths, mutation rates)
- ✅ **Multiple replicates** with recombination sweeps

**OBSERVATIONS:**

**Line 694:** Ground truth statistics computation
```python
stats, genetic_dist = self.compute_ground_truth_stats(ts_mutate, chr_idx)
```
**Question**: What genealogical metric is used? Is it generation distance (parent-child steps)?

**Line 707:** Genetic distance matrix saved
```python
np.savetxt(dist_file, genetic_dist, fmt='%.6f')
```
**Good**: You're saving pairwise genetic distances for validation.

**Line 732:** Ground truth saved per chromosome
```python
stats_df.to_csv(stats_output, sep='\t', index=False)
```
**Recommendation**: Ensure this includes:
- `id1, id2` (sample pair)
- `gen_distance` (generation steps in transmission tree)
- `is_direct_transmission` (boolean: 1 step = direct)
- `total_ibd_bp` (true IBD from tree sequence)

### 1.2 Missing from Simulation Code (View Limitation)

I couldn't see the full `compute_ground_truth_stats()` function (line 694). **Critical to verify**:

**Required ground truth fields** (for your analysis):
```python
# For each sample pair (i, j):
ground_truth = {
    'id1': sample_i,
    'id2': sample_j,
    'gen_distance': steps_in_transmission_tree,  # e.g., 1 = direct parent-child
    'is_direct_transmission': (gen_distance == 1),
    'total_ibd_bp': sum_of_ibd_segments,
    'max_segment_bp': longest_ibd_segment,
    'true_population': population_id  # For migration scenarios
}
```

**Why this matters**: Your evaluation script (line 651) uses:
```r
truth_dt <- define_ground_truth(truth_dt, gen_threshold = GEN_THRESHOLD)
```
If `GEN_THRESHOLD = 25`, you're labeling pairs ≤25 generations apart as "positive".

**CRITICAL QUESTION**: Is 25 generations biologically meaningful for malaria?
- 25 generations ≈ 25-50 mosquito cycles (assuming 10-15 days per cycle)
- This is ~250-750 days (8-25 months)
- For **direct transmission inference**, you likely want `gen_threshold = 1` or `2`
- For **cluster detection**, 5-10 generations might be appropriate

**Recommendation**: Test multiple thresholds:
- `gen_threshold = 1`: Direct transmission only
- `gen_threshold = 2-3`: Within 2-3 transmission steps
- `gen_threshold = 5-10`: Transmission cluster
- `gen_threshold = 25`: Your current setting (rationale needed)

---

## 2. EVALUATION PIPELINE REVIEW

### 2.1 Main Evaluation Script (evaluate_recombination_effects_v1_1.r)

**STRENGTHS:**
- ✅ **AUROC calculated** (line 478, 768)
- ✅ **AUPR calculated** (line 466)
- ✅ **Multiple metrics**: Spearman, Brier, Precision@K
- ✅ **Proper baseline calibration** (lines 732-738)
- ✅ **Detectability framework** (lines 656-682)
- ✅ **ROC threshold optimization** (lines 172-173)

**CRITICAL OBSERVATIONS:**

### Line 42: Ground Truth Definition
```r
define_ground_truth <- function(dt, gen_col = "gen_distance", gen_threshold = 25)
```

**Issue**: `gen_threshold = 25` is **not direct transmission**.

**Why this is important for your RQs**:
- **RQ1** asks "when can methods recover transmission links"
  - If "link" = direct transmission → use `gen_threshold = 1`
  - If "link" = transmission cluster → use `gen_threshold = 5-10`
  - 25 generations includes very distant relationships

**Recommendation**: 
```r
# Test multiple definitions:
direct_transmission <- define_ground_truth(dt, gen_threshold = 1)    # Your core analysis
cluster_2_3_steps <- define_ground_truth(dt, gen_threshold = 3)      # Secondary analysis  
transmission_chain <- define_ground_truth(dt, gen_threshold = 10)    # Cluster detection
```

### Line 473-481: AUROC Calculation (EXCELLENT)
```r
roc_obj <- tryCatch(
  pROC::roc(dt[[truth_col]], dt$score_internal, quiet = TRUE),
  error = function(e) NULL
)
res$auroc <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
```

**Good**: You're already calculating AUROC. This is the **primary metric I recommended**.

**Missing**: Sensitivity at 90% specificity.

**How to add**:
```r
# After line 478, add:
if (!is.null(roc_obj)) {
  coords <- pROC::coords(roc_obj, x = "all", ret = c("sensitivity", "specificity"))
  
  # Find sensitivity where specificity ≥ 0.90
  idx_90spec <- which(coords$specificity >= 0.90)
  if (length(idx_90spec) > 0) {
    res$sensitivity_at_90spec <- max(coords$sensitivity[idx_90spec])
  } else {
    res$sensitivity_at_90spec <- 0
  }
} else {
  res$sensitivity_at_90spec <- NA_real_
}
```

Then store in results (line 762-774):
```r
row <- data.table(
  replicate = rep,
  method = m,
  rate = rate_numeric,
  rate_label = rate,
  aupr = eval_res$aupr,
  auroc = eval_res$auroc,
  sensitivity_at_90spec = eval_res$sensitivity_at_90spec,  # ADD THIS
  spearman = eval_res$spearman,
  ...
)
```

### Line 656-682: Detectability Framework

**Observation**: You're computing "detectability" based on IBD proportion:
```r
ibd_eps_global <- define_ibd_eps(truth_dt, alpha = ALPHA_DETECT)
truth_dt[, detectable := as.integer(ibd_prop >= ibd_eps_global)]
```

**This is interesting but orthogonal to identifiability**.

**What you're measuring**: "What fraction of true transmission pairs have enough IBD signal to be detected?"

**What your RQ asks**: "Under what conditions can methods achieve accurate inference?"

**These are related but different**:
- **Detectability** (your current metric): Do true pairs have measurable IBD?
- **Identifiability** (your RQ): Can methods distinguish true from false pairs?

**Example where they diverge**:
- High recombination → low detectability (little IBD)
- But if what little IBD exists is highly specific → identifiability could still be good
- Method: IBS might work where IBD fails

**Recommendation**: Keep detectability as a secondary metric, but your primary outcome should be AUROC.

---

## 3. INFERENCE METHODS REVIEW

### 3.1 IBD Loading (recombination_evaluation_pipeline_single_run.R)

**Line 119-146: Load IBD**
```r
load_inferred_ibd <- function(rep_dir) {
  fp <- file.path(rep_dir, INFERRED_IBD_FILE)
  dt <- fread(fp)
  # ... canonicalize pair_key ...
  out <- unique(dt[, .(method = "IBD", id1, id2, score, pair_key)])
}
```

**Good**: Proper canonicalization of pair keys.

**Question**: What IBD inference tool are you using? (hmmIBD, isoRelate, Refinetti?)

### 3.2 IBS Loading

**Line 154-202: Load IBS Matrix**
```r
load_ibs_matrix <- function(rep_dir) {
  obj <- readRDS(fp)
  mat <- as.matrix(obj)
  # ... convert to long format ...
  long[, score := distance]  # Line 198
}
```

**CRITICAL ISSUE**: Line 198 says `score := distance`.

**If this is IBS (Identity-by-State)**, it should be a **similarity** (higher = more related), not distance.

**Check**:
- If your IBS matrix contains **proportion of shared alleles** (0-1, higher = similar) → score should be `distance` (keep as is)
- If your IBS matrix contains **genetic distance** (higher = less related) → score should be `-distance` or `1 - distance`

**How IBS is typically defined**:
```r
IBS_similarity = (number_of_matching_genotypes) / (total_genotypes)
# Range: 0 (completely different) to 1 (identical)
```

**Verify your IBS calculation**. If computed correctly, high IBS should = high relatedness.

### 3.3 Phylogenetic Distance Loading

**Line 205-255: Load Tree**
```r
load_tree <- function(rep_dir) {
  dt <- as.data.table(phylo_long(fp, method = "patristic"))
  # Converts distance to similarity using multiple approaches
}
```

**Good**: You're converting phylogenetic distance to similarity (lines 218-230).

**Line 750 in main script**:
```r
is_distance <- (m == "Phylo")
```

**Good**: You're flagging phylogenetic as distance and inverting it (line 449-453).

---

## 4. GROUND TRUTH CONSISTENCY CHECK

### **CRITICAL: Are all three methods using the same pairs?**

**Current implementation**:

**Line 712-720** (evaluate_recombination_effects_v1_1.r):
```r
merged <- unique(truth_dt[, .(pair_key, id1, id2, ibd_prop, 
                              true_transmission, detectable, 
                              gen_distance)])

for (m in names(method_tables)) {
  mt <- unique(method_tables[[m]][, .(pair_key, score)])
  setnames(mt, "score", paste0("score_", m))
  merged <- merge(merged, mt, by = "pair_key", all.x = TRUE)  # LEFT JOIN
}
```

**Analysis**:
- ✅ **Good**: You merge methods onto `truth_dt` with `all.x = TRUE`
- ✅ **This means**: All true pairs are retained even if method didn't predict them
- ✅ **Missing predictions** will have `NA` scores
- ✅ **AUROC calculation** (line 473-481) will drop NAs: `dt <- dt[!is.na(get(score_col))]`

**Potential issue**: If different methods drop different samples (e.g., phylogenetic tree drops low-quality samples), your AUROC for each method is calculated on **different subsets**.

**How to verify consistency**:
```r
# After line 724, add diagnostic:
for (m in names(method_tables)) {
  score_col <- paste0("score_", m)
  n_scored <- sum(!is.na(merged[[score_col]]))
  n_total <- nrow(merged)
  message("    • ", m, ": ", n_scored, " / ", n_total, 
          " pairs (", round(100 * n_scored / n_total, 1), "%)")
}
```

**If coverage differs substantially** (e.g., IBD has 95% coverage, phylo has 80%):
- This could bias comparisons
- Recommendation: Report coverage as a metric
- For fair comparison: subset to pairs scored by all methods

**How to fix** (if needed):
```r
# Find pairs scored by all methods
score_cols <- grep("^score_", names(merged), value = TRUE)
complete_pairs <- rowSums(!is.na(merged[, ..score_cols])) == length(score_cols)

# Report
message("  → Complete scoring: ", sum(complete_pairs), " / ", nrow(merged), " pairs")

# Option 1: Analyze complete set only (strict)
merged_complete <- merged[complete_pairs]

# Option 2: Report both (recommended)
# - Main analysis: all pairs (current approach)
# - Sensitivity analysis: complete pairs only
```

---

## 5. KEY RECOMMENDATIONS (PRIORITY ORDER)

### **HIGH PRIORITY (Do before final analysis)**

#### 1. Clarify ground truth definition
**Current**: `gen_threshold = 25` (very liberal)
**Recommended**: 
```r
# Primary analysis:
gen_threshold = 1  # Direct transmission (parent-child)

# Sensitivity analysis:
gen_threshold = c(1, 2, 3, 5, 10)  # Test multiple thresholds
```

**Justification needed**: Why is 25 generations biologically meaningful for Pf transmission inference?

#### 2. Add sensitivity at 90% specificity
**Implementation** (see Section 2.1 above):
```r
# In evaluate_method_on_table(), after line 478:
coords <- pROC::coords(roc_obj, x = "all", ret = c("sensitivity", "specificity"))
idx_90spec <- which(coords$specificity >= 0.90)
res$sensitivity_at_90spec <- if (length(idx_90spec) > 0) {
  max(coords$sensitivity[idx_90spec])
} else {
  0
}
```

**Why**: This is your identifiability validation metric (from my recommendations).

#### 3. Verify IBS direction
**Check**: Line 198 in recombination_evaluation_pipeline_single_run.R
```r
long[, score := distance]  # Is this correct?
```

**If IBS is similarity** (high = related) → keep as is
**If IBS is distance** (low = related) → change to `score := -distance`

**How to verify**:
```r
# Look at IBS values for known related pairs vs. unrelated
# Related pairs should have HIGH scores
```

#### 4. Check pair coverage consistency
**Add diagnostics** (see Section 4):
```r
# Report how many pairs each method scored
# Verify all methods evaluate the same pairs (or document differences)
```

---

### **MEDIUM PRIORITY (Enhance analysis quality)**

#### 5. Document ground truth extraction
**Recommendation**: Add comments or separate document explaining:
- How `gen_distance` is calculated (MRCA? parent-child steps?)
- What `total_ibd_bp` represents (from true tree sequence)
- How superinfection is handled (if at all)

#### 6. Add method comparison metrics (RQ4)
**Current**: You calculate AUROC per method separately
**Needed for RQ4**: Direct comparison

**Add after all methods evaluated**:
```r
# For each scenario:
method_aucs <- merged_summary[, .(
  auc_ibd = mean(auroc[method == "IBD"]),
  auc_ibs = mean(auroc[method == "IBS"]),
  auc_phylo = mean(auroc[method == "Phylo"])
), by = .(rate)]

method_aucs[, auc_range := pmax(auc_ibd, auc_ibs, auc_phylo) - 
                          pmin(auc_ibd, auc_ibs, auc_phylo)]

method_aucs[, failure_type := fcase(
  pmax(auc_ibd, auc_ibs, auc_phylo) < 0.70 & auc_range < 0.10, 
    "FUNDAMENTAL (all fail)",
  pmax(auc_ibd, auc_ibs, auc_phylo) >= 0.80 & pmin(auc_ibd, auc_ibs, auc_phylo) < 0.70, 
    "METHODOLOGICAL (method choice matters)",
  default = "MIXED"
)]
```

#### 7. Test hypothesis explicitly
**From your thesis**: IBD should fail at low recombination, work at high recombination

**Add statistical test**:
```r
# Hypothesis: IBD performance increases with recombination rate
# IBS performance decreases with recombination rate

library(lme4)
model_ibd <- lmer(auroc ~ log10(rate) + (1|replicate), 
                  data = results[method == "IBD"])
model_ibs <- lmer(auroc ~ log10(rate) + (1|replicate), 
                  data = results[method == "IBS"])

# Expect:
# - Positive coefficient for IBD (higher rec → higher AUC)
# - Negative coefficient for IBS (higher rec → lower AUC)
```

---

### **LOW PRIORITY (Nice to have)**

#### 8. Expand detectability analysis
**Current**: You report detectability separately (line 674)
**Enhancement**: Plot detectability vs. identifiability

```r
# For each scenario:
ggplot(summary_data, aes(x = detectability, y = auroc, color = method)) +
  geom_point() +
  facet_wrap(~rate) +
  labs(title = "Detectability vs. Identifiability",
       x = "Fraction of true pairs with detectable IBD",
       y = "AUROC") +
  theme_minimal()
```

**Expected pattern**:
- Low detectability + low AUROC → Fundamental failure
- Low detectability + high AUROC for IBS → Methodological (IBS works when IBD fails)

#### 9. Add ROC curves to output
**Current**: You calculate AUROC but don't save full ROC curves
**Enhancement**:
```r
# Save ROC coordinates for plotting
if (!is.null(roc_obj)) {
  roc_coords <- coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"))
  # Save for later plotting
}
```

---

## 6. VALIDATION CHECKLIST

Before running full 900 simulations, verify:

- [ ] **Ground truth contains**: `id1, id2, gen_distance, is_direct_transmission, total_ibd_bp`
- [ ] **Ground truth threshold justified**: Why gen_threshold = 25? (or change to 1?)
- [ ] **IBS direction correct**: High IBS score = more related (verify on known pairs)
- [ ] **Phylogenetic inversion working**: Line 449-453 correctly inverts distance
- [ ] **All methods use same pairs**: Check coverage diagnostics
- [ ] **AUROC calculated correctly**: Values between 0.5-1.0 (not 0-0.5)
- [ ] **Sensitivity at 90% specificity added**: Per my recommendations
- [ ] **Method comparison framework**: Ready for RQ4 analysis

---

## 7. EXAMPLE ANALYSIS WORKFLOW

Here's how your analysis should flow:

```r
# 1. Load results
results <- fread("all_metrics_combined.csv")

# 2. Filter to primary analysis (adjust gen_threshold if needed)
results_primary <- results  # If gen_threshold already correct

# 3. Calculate identifiability
results_primary[, identifiable := (auroc >= 0.80 & sensitivity_at_90spec >= 0.60)]

# 4. RQ1: Percentage identifiable by parameter
rq1_summary <- results_primary[, .(
  pct_identifiable = 100 * mean(identifiable),
  mean_auroc = mean(auroc),
  mean_sens = mean(sensitivity_at_90spec)
), by = .(method, rate)]

# 5. RQ2: Parameter effects
library(lme4)
model <- lmer(auroc ~ log10(rate) + bottleneck + sampling + 
              (1|replicate), data = results_primary)
summary(model)  # Standardized coefficients = effect sizes

# 6. RQ3: Migration analysis
# (Your migration scenario, not in these scripts)

# 7. RQ4: Methodological vs. fundamental
method_comparison <- results_primary[, .(
  auc_range = max(auroc) - min(auroc),
  max_auc = max(auroc),
  min_auc = min(auroc),
  best_method = method[which.max(auroc)]
), by = .(rate, replicate)]

method_comparison[, failure_type := fcase(
  max_auc < 0.70 & auc_range < 0.10, "FUNDAMENTAL",
  max_auc >= 0.80 & auc_range >= 0.15, "METHODOLOGICAL",
  default = "MIXED"
)]

# 8. Hypothesis test: IBD vs. IBS across recombination
t.test(auroc ~ method, data = results_primary[rate == 1e-9 & method %in% c("IBD", "IBS")])
t.test(auroc ~ method, data = results_primary[rate == 1e-6 & method %in% c("IBD", "IBS")])
```

---

## 8. FINAL ASSESSMENT

### Code Quality: A

**Strengths**:
- Professional structure and documentation
- Proper error handling (`tryCatch`)
- Comprehensive metrics (AUROC, AUPR, Spearman, etc.)
- Good logging and progress tracking

**Areas for improvement**:
- Ground truth definition needs justification or revision
- Missing sensitivity at 90% specificity (easy to add)
- Could benefit from method comparison framework

### Implementation Completeness: A-

**What you have**:
- ✅ Simulation framework
- ✅ Three inference methods
- ✅ AUROC calculation (my primary recommendation)
- ✅ Ground truth extraction
- ✅ Replicate evaluation

**What's missing/unclear**:
- ⚠️ Ground truth threshold (gen = 25?) needs justification
- ⚠️ Sensitivity at 90% specificity not calculated
- ⚠️ Direct method comparison for RQ4 (can add in analysis)

### Readiness for Full Analysis: B+

**Ready to proceed IF**:
1. You verify/justify `gen_threshold = 25` OR change to 1-3
2. You add sensitivity at 90% specificity calculation
3. You verify IBS direction and pair coverage

**Once these are addressed** → A (fully ready)

---

## 9. PRIORITY ACTIONS BEFORE FULL RUN

### Week 1 Tasks:

**Day 1-2**: Verify ground truth
- [ ] Check what `gen_distance` actually represents
- [ ] Decide on appropriate threshold (1 for direct transmission?)
- [ ] Document rationale

**Day 3**: Add missing metric
- [ ] Implement sensitivity at 90% specificity (code provided above)
- [ ] Test on one replicate

**Day 4**: Verify consistency
- [ ] Check IBS direction (should high = more related)
- [ ] Add pair coverage diagnostics
- [ ] Verify all methods use same pairs

**Day 5**: Pilot validation
- [ ] Run 3 test scenarios (easy, hard, realistic)
- [ ] Verify:
  - Easy scenario: All AUROC >0.85
  - Hard scenario: All AUROC <0.60
  - Realistic: Method differences (AUROC range >0.15)
  - Hypothesis: IBS > IBD at low recomb, IBD > IBS at high recomb

**If all pass** → Proceed with full 900 simulations
**If any fail** → Debug before scaling

---

## CONCLUSION

Your code is **high quality and nearly ready**. The main concerns are:

1. **Ground truth definition** (gen_threshold = 25 seems too liberal for direct transmission)
2. **Missing metric** (sensitivity at 90% specificity - easy to add)
3. **Verification needs** (IBS direction, pair coverage)

**Once these are addressed, you're fully ready to run and analyze your full design.**

**Your existing AUROC calculation is excellent** - this was my primary recommendation, and you already have it implemented.

Good work! With minor refinements, this will be a strong thesis chapter.
