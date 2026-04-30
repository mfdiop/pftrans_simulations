# Corrected Code Snippets for Recombination Analysis
## Ready-to-use fixes for distance metric handling

---

## 📝 **CORRECTED FUNCTION: load_phylo_patristic()**

Replace lines 368-399 with this corrected version:

```r
# Load patristic distances and convert to similarity scores
load_phylo_patristic <- function(rep_phylo_dir, rate_str) {
  
  # Find treefile matching rec rate pattern
  patt <- paste0("rec", rate_str)
  files <- list.files(rep_phylo_dir, full.names = TRUE, pattern = patt)
  
  # Pick modelfinder one if available
  tf <- files[grepl("_modelfinder.treefile$", files)]
  
  if (length(tf) == 0 && length(files) > 0) tf <- files[1]
  if (length(tf) == 0) return(NULL)
  
  # Expect function phylo_long() available
  if (!exists("phylo_long")) {
    helper <- file.path("simulations", "R", "patristic_distances.R")
    if (file.exists(helper)) source(helper) else {
      warning("phylo_long not found; skipping patristic for ", rep_phylo_dir)
      return(NULL)
    }
  }
  
  # Load patristic distances
  dt <- as.data.table(phylo_long(tf[1], method = "patristic"))
  id_cols <- intersect(c("id1", "id2", "Id1", "Id2", "sample1", "sample2"), names(dt))
  
  if (length(id_cols) >= 2) {
    setnames(dt, old = id_cols[1:2], new = c("id1", "id2"), skip_absent = TRUE)
  }
  
  # ✅ CRITICAL FIX: Convert distance to similarity
  # Patristic distances: smaller = more related
  # We need: higher score = more related (for AUPR/AUROC)
  
  if ("phylo" %in% names(dt)) {
    max_dist <- max(dt$phylo, na.rm = TRUE)
    # Option 1: Subtract from max (preserves scale)
    dt[, score := max_dist - phylo]
    
    # Option 2 (alternative): Inverse transform
    # dt[, score := 1 / (1 + phylo)]
    
    # Option 3 (alternative): Negative (for use with invert_dist=TRUE)
    # dt[, score := -phylo]
    
  } else if ("distance" %in% names(dt)) {
    max_dist <- max(dt$distance, na.rm = TRUE)
    dt[, score := max_dist - distance]
    
  } else {
    warning("No distance column found in phylogenetic data")
    dt[, score := 0]
  }
  
  if (!("pair_key" %in% names(dt))) {
    dt[, pair_key := mapply(canonical_pair, id1, id2)]
  }
  
  # Add metadata for tracking (optional but recommended)
  dt[, metric_original_type := "distance"]
  dt[, metric_converted_to := "similarity"]
  
  message("    • Phylo distances converted to similarity (range: ",
          round(min(dt$score, na.rm = TRUE), 4), " - ",
          round(max(dt$score, na.rm = TRUE), 4), ")")
  
  unique(dt[, .(pair_key, id1, id2, score)])
}
```

**Key Changes:**
1. Added actual distance-to-similarity conversion (line with `max_dist - phylo`)
2. Added informative message showing score range
3. Removed misleading comment that claimed conversion was happening
4. Added optional metadata tracking

---

## 📝 **ALTERNATIVE APPROACH: Track Metric Type & Use invert_dist**

If you prefer to keep distances as distances and handle them in evaluation:

### Modified load_phylo_patristic():
```r
load_phylo_patristic <- function(rep_phylo_dir, rate_str) {
  # ... [same file finding code] ...
  
  dt <- as.data.table(phylo_long(tf[1], method = "patristic"))
  
  # ... [same column renaming code] ...
  
  # Keep as distance (don't convert)
  if ("phylo" %in% names(dt)) {
    dt[, score := phylo]
  } else if ("distance" %in% names(dt)) {
    dt[, score := distance]
  } else {
    dt[, score := 0]
  }
  
  if (!("pair_key" %in% names(dt))) {
    dt[, pair_key := mapply(canonical_pair, id1, id2)]
  }
  
  # ✅ Add attribute to track this is a distance metric
  attr(dt, "metric_type") <- "distance"
  
  unique(dt[, .(pair_key, id1, id2, score)])
}
```

### Modified data loading section (lines 715-760):
```r
# Load methods
message("  Loading method predictions...")

ibd_hmm_dt <- load_inferred_ibd(rep_dir, rate)
if (!is.null(ibd_hmm_dt)) {
  attr(ibd_hmm_dt, "metric_type") <- "similarity"  # IBD is similarity
  message(" ✓ HMM-IBD: ", nrow(ibd_hmm_dt), " pairs")
}

ibs_dt <- load_ibs_matrix(rep_dir, rate)
if (!is.null(ibs_dt)) {
  # ⚠️ VERIFY YOUR IBS CALCULATION!
  # If it's proportion shared: similarity
  # If it's Hamming distance or 1-proportion: distance
  attr(ibs_dt, "metric_type") <- "similarity"  # ← CHANGE THIS IF NEEDED
  message(" ✓ IBS: ", nrow(ibs_dt), " pairs")
}

phylo_dt <- NULL
phylo_rep_dir <- file.path(PHYLO_ROOT, rep)
if (dir.exists(phylo_rep_dir)) {
  phylo_dt <- load_phylo_patristic(phylo_rep_dir, rate)
  if (!is.null(phylo_dt)) {
    # Already has metric_type = "distance" from function
    message(" ✓ Phylo: ", nrow(phylo_dt), " pairs")
  }
}

# Store method tables with metadata
method_tables <- list()
if (!is.null(ibd_hmm_dt)) {
  method_tables$IBD <- list(
    data = ibd_hmm_dt,
    type = attr(ibd_hmm_dt, "metric_type")
  )
}
if (!is.null(ibs_dt)) {
  method_tables$IBS <- list(
    data = ibs_dt,
    type = attr(ibs_dt, "metric_type")
  )
}
if (!is.null(phylo_dt)) {
  method_tables$Phylo <- list(
    data = phylo_dt,
    type = attr(phylo_dt, "metric_type")
  )
}
```

### Modified merging section (lines 740-760):
```r
# Merge methods with truth
message("  Merging data...")

# Base merged table (truth)
merged <- unique(truth_dt[, .(pair_key, id1, id2, ibd_prop, 
                              true_transmission, detectable, gen_distance)])

# Track metric types for later use
metric_types <- list()

for (m in names(method_tables)) {
  mt <- unique(method_tables[[m]]$data[, .(pair_key, score)])
  setnames(mt, "score", paste0("score_", m))
  merged <- merge(merged, mt, by = "pair_key", all.x = TRUE)
  
  # Store metric type
  metric_types[[m]] <- method_tables[[m]]$type
}

# Add metric types as attribute
attr(merged, "metric_types") <- metric_types

message("    ✓ Merged table: ", nrow(merged), " pairs × ", ncol(merged), " columns")
message("    • Metric types: ", paste(names(metric_types), "=", unlist(metric_types), collapse = ", "))
```

### Modified evaluation loop (lines 763-830):
```r
# Evaluate each method
message("  Evaluating methods...")

metric_types <- attr(merged, "metric_types")

for (m in names(method_tables)) {
  score_col <- paste0("score_", m)
  
  # ✅ CRITICAL FIX: Set invert_dist based on metric type
  is_distance <- (metric_types[[m]] == "distance")
  
  message("    • Evaluating ", m, " (type: ", metric_types[[m]], ")")
  
  eval_res <- evaluate_method_on_table(
    merged, 
    score_col = score_col, 
    truth_col = "true_transmission", 
    ks = TOP_K,
    invert_dist = is_distance  # ✅ NOW CORRECTLY SET!
  )
  
  # [rest of the code remains the same]
  row <- data.table(
    replicate = rep,
    method = m,
    rate = rate_numeric,
    rate_label = rate,
    aupr = eval_res$aupr,
    auroc = eval_res$auroc,
    spearman = eval_res$spearman,
    brier = eval_res$brier,
    detectability = detectability_rate,
    n_detectable = n_detectable,
    n_true_positives = n_true_pos,
    metric_type = metric_types[[m]]  # Add this for tracking
  )
  
  for (k in TOP_K) {
    row[[paste0("prec_at_", k)]] <- eval_res[[paste0("prec_at_", k)]]
  }
  
  all_metrics_rows[[length(all_metrics_rows) + 1]] <- row
  
  # ... [rest of storage code] ...
  
  message("    ✓ ", m, " (", metric_types[[m]], "): ",
          "AUPR = ", round(eval_res$aupr, 3),
          ", AUROC = ", round(eval_res$auroc, 3),
          ", Prec@10 = ", round(eval_res$prec_at_10, 3))
}
```

---

## 📝 **DIAGNOSTIC CODE TO ADD BEFORE EVALUATION**

Insert this after merging (around line 758):

```r
# ============================================================================
# DIAGNOSTIC: Verify metric directions before evaluation
# ============================================================================
message("\n  [DIAGNOSTIC] Checking metric directions...")

for (m in names(method_tables)) {
  score_col <- paste0("score_", m)
  
  if (!score_col %in% names(merged)) next
  
  # Calculate mean scores for true positives vs negatives
  mean_pos <- merged[true_transmission == 1, mean(get(score_col), na.rm = TRUE)]
  mean_neg <- merged[true_transmission == 0, mean(get(score_col), na.rm = TRUE)]
  median_pos <- merged[true_transmission == 1, median(get(score_col), na.rm = TRUE)]
  median_neg <- merged[true_transmission == 0, median(get(score_col), na.rm = TRUE)]
  
  # Calculate correlation with IBD proportion (should be positive for all)
  cor_ibd <- cor(merged[[score_col]], merged$ibd_prop, 
                 method = "spearman", use = "complete.obs")
  
  # Determine expected direction
  metric_type <- if (exists("metric_types")) metric_types[[m]] else "unknown"
  
  message("\n    ", m, " (", metric_type, "):")
  message("      • Mean: positives = ", round(mean_pos, 4), 
          ", negatives = ", round(mean_neg, 4))
  message("      • Median: positives = ", round(median_pos, 4), 
          ", negatives = ", round(median_neg, 4))
  message("      • Correlation with IBD proportion: ", round(cor_ibd, 4))
  
  # Validation checks
  if (metric_type == "similarity") {
    if (mean_pos < mean_neg) {
      warning("      ⚠️ WARNING: Similarity metric has LOWER values for positives!")
      warning("         This suggests metric is actually a distance, not similarity")
    } else {
      message("      ✓ Direction correct (positives > negatives)")
    }
  } else if (metric_type == "distance") {
    if (mean_pos > mean_neg) {
      warning("      ⚠️ WARNING: Distance metric has HIGHER values for positives!")
      warning("         This suggests metric is actually a similarity, not distance")
      warning("         OR invert_dist parameter is not being applied correctly")
    } else {
      message("      ✓ Direction correct (positives < negatives for distance)")
    }
  }
  
  if (abs(cor_ibd) < 0.3) {
    warning("      ⚠️ Low correlation with IBD (|r| = ", round(abs(cor_ibd), 3), ")")
  } else if (cor_ibd < 0 && metric_type == "similarity") {
    warning("      ⚠️ NEGATIVE correlation for similarity metric - likely inverted!")
  }
}

message("\n  [DIAGNOSTIC] Complete. Check warnings above!\n")
```

---

## 📝 **ENHANCED evaluate_method_on_table() FUNCTION**

Replace the function definition (lines 405-573) with this enhanced version:

```r
#' Evaluate a single method on detectability-conditioned data
#' 
#' @param dt_merge data.table with scores and truth
#' @param score_col Name of column containing scores
#' @param truth_col Name of truth column (default: "true_transmission")
#' @param threshold Optional threshold (NULL = use Youden's index)
#' @param ks Vector of K values for Precision@K
#' @param invert_dist Logical. If TRUE, scores are distances (negate for evaluation)
#' @return List of evaluation metrics
evaluate_method_on_table <- function(
    dt_merge,
    score_col = "score",
    truth_col = "true_transmission",
    threshold = NULL,
    ks = c(5, 10, 25, 50),
    invert_dist = FALSE
) 
{
  res <- list()
  res$invert_dist <- invert_dist  # Track what was used
  
  # If empty input
  if (is.null(dt_merge) || nrow(dt_merge) == 0 || !score_col %in% names(dt_merge)) {
    null_metrics <- c("aupr", "auroc", "spearman", "brier")
    for (m in null_metrics) res[[m]] <- NA_real_
    for (k in ks) res[[paste0("prec_at_", k)]] <- NA_real_
    res$TP <- res$TN <- res$FP <- res$FN <- NA_integer_
    res$precision <- res$recall <- res$f1 <- NA_real_
    res$confusion_matrix <- NULL
    return(res)
  }
  
  dt <- copy(dt_merge)
  dt <- dt[!is.na(get(score_col))]
  
  if (nrow(dt) == 0) return(res)
  
  # ----------------------------
  # HANDLE DISTANCE OR SIMILARITY
  # ----------------------------
  s <- dt[[score_col]]
  
  if (invert_dist) {
    # Score is a DISTANCE (smaller = more related)
    # Convert to similarity-like (larger = more related)
    # Option 1: Simple negation
    s_use <- -s
    
    # Option 2: Max - distance (better for interpretability)
    # max_s <- max(s, na.rm = TRUE)
    # s_use <- max_s - s
  } else {
    # Score is already a SIMILARITY (larger = more related)
    s_use <- s
  }
  
  dt[, score_internal := s_use]
  
  # Store raw and transformed scores for diagnostics
  res$score_raw_range <- range(s, na.rm = TRUE)
  res$score_internal_range <- range(s_use, na.rm = TRUE)
  
  # ----------------------------
  # AUPR (Precision-Recall)
  # ----------------------------
  pos_scores <- dt[get(truth_col) == 1, score_internal]
  neg_scores <- dt[get(truth_col) == 0, score_internal]
  
  res$n_pos <- length(pos_scores)
  res$n_neg <- length(neg_scores)
  
  if (length(pos_scores) > 0 && length(neg_scores) > 0) {
    pr <- tryCatch(
      PRROC::pr.curve(scores.class0 = pos_scores,
                      scores.class1 = neg_scores,
                      curve = TRUE),
      error = function(e) {
        warning("PR curve calculation failed: ", e$message)
        NULL
      }
    )
    
    if (!is.null(pr)) {
      res$pr <- pr
      res$aupr <- pr$auc.integral
      
      # Diagnostic: AUPR should be > 0.5 for useful classifier
      if (pr$auc.integral < 0.5) {
        warning("⚠️ AUPR < 0.5 (", round(pr$auc.integral, 3), 
                ") suggests inverted metric!")
      }
    } else {
      res$pr <- NA
      res$aupr <- NA_real_
    }
  } else {
    res$pr <- NA
    res$aupr <- NA_real_
  }
  
  # ----------------------------
  # AUROC (Receiver Operating Characteristic)
  # ----------------------------
  if (length(unique(dt[[truth_col]])) > 1) {
    roc_obj <- tryCatch(
      pROC::roc(dt[[truth_col]], dt$score_internal, quiet = TRUE),
      error = function(e) {
        warning("ROC calculation failed: ", e$message)
        NULL
      }
    )
    
    if (!is.null(roc_obj)) {
      res$roc <- roc_obj
      res$auroc <- pROC::auc(roc_obj)
      
      # Diagnostic: AUROC should be > 0.5 for useful classifier
      if (pROC::auc(roc_obj) < 0.5) {
        warning("⚠️ AUROC < 0.5 (", round(pROC::auc(roc_obj), 3), 
                ") suggests inverted metric!")
      }
    } else {
      res$roc <- NA
      res$auroc <- NA_real_
      roc_obj <- NULL
    }
  } else {
    res$roc <- NA
    res$auroc <- NA_real_
    roc_obj <- NULL
  }
  
  # ----------------------------
  # Spearman correlation vs continuous truth (IBD proportion)
  # ----------------------------
  if ("ibd_prop" %in% names(dt)) {
    res$spearman <- suppressWarnings(
      cor(dt$score_internal, dt$ibd_prop, method = "spearman", use = "complete.obs")
    )
    
    # Diagnostic: Should be positive for all methods
    if (!is.na(res$spearman) && res$spearman < 0) {
      warning("⚠️ Negative Spearman correlation (", round(res$spearman, 3), 
              ") with IBD proportion!")
      if (!invert_dist) {
        warning("   Metric appears to be a distance - set invert_dist=TRUE")
      }
    }
  } else {
    res$spearman <- NA_real_
  }
  
  # ----------------------------
  # Brier score (calibration metric)
  # ----------------------------
  rng <- range(dt$score_internal, na.rm = TRUE)
  if (diff(rng) == 0) {
    pred_prob <- rep(0.5, nrow(dt))
  } else {
    pred_prob <- (dt$score_internal - rng[1]) / diff(rng)
  }
  res$brier <- mean((pred_prob - dt[[truth_col]])^2, na.rm = TRUE)
  
  # ----------------------------
  # Precision@K (ranking metric)
  # ----------------------------
  # For similarity: sort descending (highest scores first)
  # For distance (after inversion): also sort descending
  ord <- dt[order(-score_internal)]  # Always descending after handling invert_dist
  
  for (k in ks) {
    topk <- head(ord, k)
    res[[paste0("prec_at_", k)]] <-
      if (nrow(topk) == 0) NA_real_ else mean(topk[[truth_col]] == 1)
  }
  
  # ----------------------------
  # Confusion matrix at optimal threshold
  # ----------------------------
  if (!is.null(roc_obj)) {
    # Use Youden's index for optimal threshold
    coords <- pROC::coords(roc_obj, "best", best.method = "youden")
    thr <- coords$threshold
    res$threshold_used <- thr
    res$threshold_method <- "youden"
    res$threshold_sensitivity <- coords$sensitivity
    res$threshold_specificity <- coords$specificity
  } else {
    # Fallback to median if ROC failed
    thr <- median(dt$score_internal, na.rm = TRUE)
    res$threshold_used <- thr
    res$threshold_method <- "median"
    res$threshold_sensitivity <- NA
    res$threshold_specificity <- NA
  }
  
  pred <- as.integer(dt$score_internal >= thr)
  
  TP <- sum(pred == 1 & dt[[truth_col]] == 1)
  FP <- sum(pred == 1 & dt[[truth_col]] == 0)
  TN <- sum(pred == 0 & dt[[truth_col]] == 0)
  FN <- sum(pred == 0 & dt[[truth_col]] == 1)
  
  res$TP <- TP
  res$FP <- FP
  res$TN <- TN
  res$FN <- FN
  
  precision <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
  recall <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0)
    NA_real_ else 2 * precision * recall / (precision + recall)
  
  res$precision <- precision
  res$recall <- recall
  res$f1 <- f1
  
  # Diagnostic summary
  res$diagnostic_summary <- paste0(
    "Metric type: ", ifelse(invert_dist, "distance (inverted)", "similarity"),
    " | AUPR: ", round(res$aupr, 3),
    " | AUROC: ", round(res$auroc, 3),
    " | Spearman: ", round(res$spearman, 3)
  )
  
  return(res)
}
```

**Key enhancements:**
1. Added diagnostic checks for AUPR < 0.5
2. Added warnings for negative correlations
3. Stores raw and transformed score ranges
4. Simplified Precision@K logic (always descending after inversion)
5. Added diagnostic summary string
6. Better error handling with tryCatch

---

## 📝 **VALIDATION SCRIPT TO RUN BEFORE FULL ANALYSIS**

Save this as `validate_metrics.R` and run before your main analysis:

```r
# validate_metrics.R
# Quick validation of metric directions before running full analysis

library(data.table)
library(PRROC)
library(pROC)

# Test data: create pairs with known relatedness
set.seed(123)
n_related <- 100
n_unrelated <- 400

test_dt <- data.table(
  pair_id = 1:(n_related + n_unrelated),
  true_transmission = c(rep(1, n_related), rep(0, n_unrelated))
)

# Simulate scores
# IBD: similarity (0-1), higher for related
test_dt[, score_IBD := ifelse(true_transmission == 1,
                               runif(.N, 0.3, 0.9),  # Related: high
                               runif(.N, 0.0, 0.2))] # Unrelated: low

# IBS: similarity (0-1), higher for related (with noise)
test_dt[, score_IBS := ifelse(true_transmission == 1,
                               runif(.N, 0.6, 0.95),  # Related: high
                               runif(.N, 0.4, 0.75))] # Unrelated: medium-low

# Phylo: DISTANCE, smaller for related
test_dt[, score_Phylo_dist := ifelse(true_transmission == 1,
                                      runif(.N, 0.0, 0.3),  # Related: small
                                      runif(.N, 0.5, 1.0))] # Unrelated: large

# Phylo converted to similarity
test_dt[, score_Phylo_sim := max(score_Phylo_dist) - score_Phylo_dist]

# Test evaluation
test_similarity <- function(dt, score_col, label) {
  pos <- dt[true_transmission == 1, get(score_col)]
  neg <- dt[true_transmission == 0, get(score_col)]
  
  pr <- pr.curve(scores.class0 = pos, scores.class1 = neg)
  roc <- roc(dt$true_transmission, dt[[score_col]], quiet = TRUE)
  
  cat("\n", label, ":\n")
  cat("  Mean: positives =", round(mean(pos), 3), 
      ", negatives =", round(mean(neg), 3), "\n")
  cat("  AUPR =", round(pr$auc.integral, 3), "\n")
  cat("  AUROC =", round(auc(roc), 3), "\n")
  cat("  Expected: AUPR > 0.5, AUROC > 0.5, mean(pos) > mean(neg)\n")
  
  if (pr$auc.integral < 0.5 || auc(roc) < 0.5) {
    cat("  ❌ FAIL: Metric appears inverted!\n")
  } else if (mean(pos) > mean(neg)) {
    cat("  ✅ PASS: Correctly oriented similarity metric\n")
  } else {
    cat("  ⚠️  WARNING: Unexpected pattern\n")
  }
}

cat("\n========================================")
cat("\n   METRIC DIRECTION VALIDATION TEST")
cat("\n========================================\n")

test_similarity(test_dt, "score_IBD", "IBD (similarity)")
test_similarity(test_dt, "score_IBS", "IBS (similarity)")
test_similarity(test_dt, "score_Phylo_dist", "Phylo WRONG (raw distance)")
test_similarity(test_dt, "score_Phylo_sim", "Phylo CORRECT (converted)")

cat("\n========================================")
cat("\n   INTERPRETATION:")
cat("\n========================================\n")
cat("If 'Phylo WRONG' shows AUPR < 0.5 → confirms inversion problem\n")
cat("If 'Phylo CORRECT' shows AUPR > 0.5 → confirms fix works\n")
cat("\nApply the same conversion in your load_phylo_patristic() function!\n\n")
```

Run this and verify:
- IBD ✅
- IBS ✅  
- Phylo WRONG ❌ (should have AUPR < 0.5)
- Phylo CORRECT ✅ (should have AUPR > 0.5)

---

## 📋 **CHECKLIST BEFORE RUNNING FULL ANALYSIS**

- [ ] ✅ Phylogenetic distances converted to similarity in `load_phylo_patristic()`
- [ ] ✅ IBS metric type verified (distance or similarity)
- [ ] ✅ Metric types tracked through pipeline
- [ ] ✅ `invert_dist` parameter passed correctly
- [ ] ✅ Diagnostic code added before evaluation
- [ ] ✅ Validation script run and passed
- [ ] ✅ All methods show AUPR > 0.5 on test data
- [ ] ✅ All methods show positive correlation with IBD proportion
- [ ] ✅ Score distributions checked (positives should score higher)

Once all items checked, your analysis should correctly compare methods! 🎉

---

## 🔍 **HOW TO VERIFY THE FIX WORKED**

After rerunning with fixes, check these key indicators:

### 1. AUPR Values
```r
# All should be > 0.5 (better than random)
agg_dt[, .(method, rate, aupr_mean)]

# If Phylo still < 0.5, conversion didn't work
```

### 2. Spearman Correlations
```r
# All should be positive
metrics_dt[, .(method, median_spearman = median(spearman, na.rm = TRUE)), 
           by = method]

# Negative correlation = still inverted
```

### 3. Recombination Effect Pattern
```r
# Plot AUPR across rates - should see expected patterns
ggplot(agg_dt, aes(x = rate, y = aupr_mean, color = method)) +
  geom_line() +
  geom_point() +
  scale_x_log10()

# Expected:
# - IBD: decreases with higher recombination (segments broken)
# - IBS/Phylo: more stable (use all SNPs)
```

### 4. Method Ranking
```r
# At baseline (lowest recombination)
agg_dt[rate == min(rate), .(method, aupr_mean)][order(-aupr_mean)]

# Phylo should be competitive with IBS, not far below
```

If all 4 checks pass → fix successful! 🎯
