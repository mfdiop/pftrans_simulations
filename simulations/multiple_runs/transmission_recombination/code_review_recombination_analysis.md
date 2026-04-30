# Code Review: Recombination Effects Analysis
## Critical Issues with Distance Metrics and Method Comparison

---

## 🚨 **CRITICAL ISSUE #1: Phylogenetic Distance-to-Similarity Conversion Not Implemented**

### Location
- **Lines 368-399**: `load_phylo_patristic()` function
- **Line 390**: Comment claims "Distances are converted to similarity"
- **Lines 395-396**: Score assignment

### The Problem
```r
# Line 390: Comment says distances → similarity
dt <- as.data.table(phylo_long(tf[1], method = "patristic")) # Distances are converted to similarity

# Lines 395-396: But no conversion actually happens!
if ("phylo" %in% names(dt)) dt[, score := phylo] 
else if ("distance" %in% names(dt)) dt[, score := distance] 
else dt[, score := 0]
```

**The code directly assigns patristic distances as scores WITHOUT any transformation.** Patristic distances are DISTANCES (larger values = less related), not similarities.

### Why This Breaks Your Analysis
1. **IBD scores**: Higher values = more IBD sharing = MORE related ✅
2. **IBS scores**: Likely higher values = more sharing = MORE related ✅
3. **Phylogenetic distances**: Higher values = MORE evolutionary distance = LESS related ❌

You're comparing apples and oranges. When you compute AUPR/AUROC, the metrics assume higher scores indicate positive class membership. For phylogenetic distances, **the relationship is inverted**.

### Expected Impact on Results
- **Phylogenetic method will appear to perform poorly** because:
  - Related pairs (true positives) have SMALL distances
  - Unrelated pairs (true negatives) have LARGE distances
  - Your ROC/PR curve calculations expect the opposite
  - This will produce AUPR < 0.5 (worse than random!)

---

## 🚨 **CRITICAL ISSUE #2: Missing `invert_dist` Parameter in Function Calls**

### Location
- **Lines 770-775**: Call to `evaluate_method_on_table()`
- **Lines 405-573**: Function definition with `invert_dist` parameter

### The Problem
```r
# Line 770-775: Function call WITHOUT invert_dist parameter
eval_res <- evaluate_method_on_table(
  merged, 
  score_col = score_col, 
  truth_col = "true_transmission", 
  ks = TOP_K
  # invert_dist parameter is MISSING!
)

# Line 412: Function has this parameter with default = FALSE
evaluate_method_on_table <- function(
    ...
    invert_dist = FALSE  # Default assumes scores are similarities
)
```

### Why This Matters
The function HAS logic to handle distance metrics (lines 438-443):
```r
if (invert_dist) {
  # higher value should mean *more related*
  s_use <- -s
} else {
  s_use <- s
}
```

**But you never set `invert_dist = TRUE` for phylogenetic distances!**

This means:
- Phylogenetic distances are treated as similarities
- Related pairs get LOW scores
- Unrelated pairs get HIGH scores
- All downstream metrics (AUPR, AUROC, Precision@K) are inverted/incorrect

---

## 🚨 **CRITICAL ISSUE #3: IBS Metric Direction Unclear**

### Location
- **Lines 312-366**: `load_ibs_matrix()` function
- **Lines 363-365**: Comment acknowledges the ambiguity

### The Problem
```r
# Line 363-365
# If IBS is a distance (smaller=more related), 
# we may invert later; here keep raw score and mark method
unique(dt[, .(pair_key, id1, id2, score)])
```

The comment says "we may invert later" but:
1. You never actually invert it
2. You don't track whether IBS is a distance or similarity
3. Different IBS implementations use different conventions!

### Common IBS Metrics
- **Proportion of shared alleles** (0-1): Similarity (higher = more related) ✅
- **Hamming distance**: Distance (higher = less related) ❌
- **1 - IBS proportion**: Distance (higher = less related) ❌

**You need to know which one you're using!**

---

## 🚨 **ISSUE #4: Inconsistent Distance Handling in Precision@K**

### Location
- **Lines 502-516**: Precision@K calculation

### The Problem
```r
# Lines 504-510
if (invert_dist) {
  # similarity-like; sort descending
  ord <- dt[order(-score_internal)]
} else {
  # distance-like; sort ascending
  ord <- dt[order(score_internal)]
}
```

This logic is BACKWARDS for distance metrics when `invert_dist = FALSE`:
- For distances: Smaller values = more related → should rank ASCENDING ✅
- For similarities: Larger values = more related → should rank DESCENDING ✅

**BUT** since phylogenetic scores are distances and `invert_dist = FALSE` (default), the code sorts ascending, which is correct by accident. However, the logic at lines 438-443 doesn't negate the score, so AUPR/AUROC calculations are still wrong!

---

## 🚨 **ISSUE #5: Cophenetic vs Patristic Distance Confusion**

### Location
- **Title/Comments**: Mentions "patristic/cophenetic"
- **Line 390**: Uses `method = "patristic"`

### Conceptual Issue
These are different metrics:

| Metric | Definition | Use Case |
|--------|-----------|----------|
| **Patristic** | Sum of branch lengths along path | Evolutionary distance (time/mutations) |
| **Cophenetic** | Height of node where taxa join | Clustering/UPGMA distance |

For phylogenies:
- **Patristic distances** reflect total evolutionary change
- **Cophenetic distances** reflect tree topology/hierarchy

**Which one should you use?**
- If your tree branch lengths are in substitutions/time: **Patristic** ✅
- If analyzing ultrametric trees (molecular clock): **Either** ✅
- If tree is non-ultrametric: **Patristic preferred** ✅

The code uses patristic (line 390) but the title mentions both. Clarify which you intended!

---

## 🔧 **RECOMMENDED FIXES**

### Fix #1: Properly Handle Phylogenetic Distances
```r
load_phylo_patristic <- function(rep_phylo_dir, rate_str) {
  # ... existing code ...
  
  dt <- as.data.table(phylo_long(tf[1], method = "patristic"))
  
  # ... existing code ...
  
  # CONVERT DISTANCE TO SIMILARITY
  if ("phylo" %in% names(dt)) {
    max_dist <- max(dt$phylo, na.rm = TRUE)
    dt[, score := max_dist - phylo]  # Invert: smaller distance → higher similarity
  } else if ("distance" %in% names(dt)) {
    max_dist <- max(dt$distance, na.rm = TRUE)
    dt[, score := max_dist - distance]  # Invert: smaller distance → higher similarity
  } else {
    dt[, score := 0]
  }
  
  # OR: Add metadata to track metric type
  dt[, metric_type := "similarity"]  # Now it's a similarity after conversion
  
  unique(dt[, .(pair_key, id1, id2, score)])
}
```

### Fix #2: Pass Metric Type Through Pipeline
```r
# Modify method_tables to include metadata
method_tables <- list()
if (!is.null(ibd_hmm_dt)) {
  method_tables$IBD <- list(data = ibd_hmm_dt, type = "similarity")
}
if (!is.null(ibs_dt)) {
  method_tables$IBS <- list(data = ibs_dt, type = "similarity")  # OR "distance" - verify!
}
if (!is.null(phylo_dt)) {
  method_tables$Phylo <- list(data = phylo_dt, type = "distance")
}

# Update evaluation call
for (m in names(method_tables)) {
  score_col <- paste0("score_", m)
  is_distance <- (method_tables[[m]]$type == "distance")
  
  eval_res <- evaluate_method_on_table(
    merged, 
    score_col = score_col, 
    truth_col = "true_transmission", 
    ks = TOP_K,
    invert_dist = is_distance  # ✅ Now correctly set!
  )
  # ...
}
```

### Fix #3: Verify IBS Metric Type
Check your IBS matrix generation code to determine if it's:
```r
# Option A: Similarity (proportion shared)
ibs_matrix[i,j] <- sum(geno[i,] == geno[j,]) / ncol(geno)

# Option B: Distance (Hamming)
ibs_matrix[i,j] <- sum(geno[i,] != geno[j,]) / ncol(geno)

# Option C: Distance (1 - proportion)
ibs_matrix[i,j] <- 1 - (sum(geno[i,] == geno[j,]) / ncol(geno))
```

Then set the appropriate type in `method_tables`.

### Fix #4: Alternative - Normalize All Metrics to [0,1] Similarity
```r
# In evaluate_method_on_table, BEFORE computing metrics:
if (invert_dist) {
  # Convert distance to similarity [0,1]
  s_use <- 1 - norm_to_01(s)  # Min dist → max similarity
} else {
  # Already similarity, just normalize
  s_use <- norm_to_01(s)
}
```

This ensures all metrics are on the same scale and direction.

---

## ⚠️ **VALIDATION CHECKS TO RUN**

### Check #1: Score Distributions
```r
# After merging, before evaluation
summary(merged$score_IBD)     # Should be high for related pairs
summary(merged$score_IBS)     # Should be high for related pairs
summary(merged$score_Phylo)   # Currently: HIGH for UNRELATED pairs (WRONG!)

# Split by truth
merged[true_transmission == 1, summary(score_Phylo)]  # Should be HIGH (currently LOW)
merged[true_transmission == 0, summary(score_Phylo)]  # Should be LOW (currently HIGH)
```

### Check #2: AUPR Values
If phylogenetic AUPR < 0.5, you have the distance/similarity inversion problem!

```r
# After running analysis
agg_dt[method == "Phylo", .(rate, aupr_mean)]
# If aupr_mean < 0.5 consistently → inverted metric confirmed
```

### Check #3: Correlation Signs
```r
# All methods should have POSITIVE Spearman correlation with IBD proportion
metrics_dt[, .(method, median_spearman = median(spearman, na.rm = TRUE)), by = method]

# If Phylo has NEGATIVE correlation → inverted metric confirmed
```

---

## 📊 **WHY YOU EXPECT VARIATION ACROSS RECOMBINATION RATES**

Your expectation is correct:
1. **Higher recombination** → IBD segments broken into smaller pieces
2. **Smaller IBD segments** → harder for HMM-IBD to detect
3. **Detection threshold issues** → affects Sensitivity/FDR

However, **if phylogenetic distances are inverted**, you won't see the expected patterns because:
- The metric is measuring the OPPOSITE of what you think
- AUPR will decrease (or stay artificially low) across all rates
- Method rankings will be incorrect

---

## ✅ **SUMMARY: What Needs to Change**

| Issue | Current Behavior | Required Fix |
|-------|------------------|--------------|
| Phylo distances | Used as-is (distance) | Convert to similarity OR set `invert_dist=TRUE` |
| `invert_dist` param | Never passed | Pass based on metric type |
| IBS type | Unknown | Verify and document if distance/similarity |
| Metric comparison | Apples vs oranges | Ensure all metrics have same direction |
| AUPR interpretation | May be inverted for Phylo | Validate with correlation checks |

---

## 🎯 **IMMEDIATE ACTION ITEMS**

1. **Add diagnostic checks** before running full analysis:
   ```r
   # Check score distributions by true class
   for (m in c("IBD", "IBS", "Phylo")) {
     cat("\n", m, ":\n")
     print(merged[, .(
       mean_true_pos = mean(get(paste0("score_", m))[true_transmission == 1], na.rm=TRUE),
       mean_true_neg = mean(get(paste0("score_", m))[true_transmission == 0], na.rm=TRUE)
     )])
   }
   # For IBD and IBS: mean_true_pos should be HIGHER
   # For Phylo (if distance): mean_true_pos should be LOWER
   ```

2. **Implement proper metric type handling** (see Fix #2)

3. **Re-run analysis** and check if:
   - All AUPR values > 0.5
   - Phylogenetic method shows expected degradation with recombination
   - Correlation signs are all positive

4. **Document metric definitions** in your methods section

---

## 📚 **CONCEPTUAL BACKGROUND**

### Why Metric Direction Matters for AUPR/AUROC

These metrics assume scores are **discriminant functions** where:
- Higher score → more likely to be positive class
- Lower score → more likely to be negative class

**For distance metrics**, this assumption is violated:
- Related pairs (positive class) have LOW distances
- Unrelated pairs (negative class) have HIGH distances

**Without inversion**, ROC curves will be below diagonal (AUC < 0.5)!

### Why This Affects Your Biological Question

Recombination breaks up IBD segments but doesn't change:
- IBS at the SNP level (still similar genotypes)
- Phylogenetic distance (still evolutionarily related)

**Expected pattern**:
- IBD-based methods: Sensitivity ↓ as recombination ↑ (segments too small)
- IBS-based methods: Relatively stable (uses all SNPs)
- Phylogenetic methods: Relatively stable (uses all SNPs for tree)

**But you can't observe these patterns if metrics are inverted!**

---

## 📝 **BOTTOM LINE**

Your analysis has a systematic error in how phylogenetic distances are handled. The metric is being used as-is (distance) when the evaluation framework expects similarities. This will:

1. Make phylogenetic methods appear to perform poorly (artificially low AUPR)
2. Produce incorrect method rankings
3. Obscure the true effects of recombination on different relatedness metrics

**Fix the distance/similarity handling BEFORE interpreting results.**
