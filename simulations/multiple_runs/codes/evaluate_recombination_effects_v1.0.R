# evaluate_recombination_full_pipeline_CORRECTED.R
# ============================================================================
# CORRECTED Pipeline for IBD Method Evaluation Under Recombination
# ============================================================================
# 
# Key Corrections:
# 1. Fixed rate conversion bug (1/as.numeric(rate) -> proper lookup)
# 2. Fixed rank_norm() direction (removed minus sign)
# 3. Separated biological truth from detectability
# 4. Calibrated thresholds at baseline only
# 5. Added comprehensive logging and progress tracking
# 6. Reports detectability alongside AUPR
#
# Authors: [Your Name]
# Date: 2025-01-XX
# ============================================================================

library(data.table)
library(tidyverse)
library(PRROC)
library(pROC)
library(gridExtra)
library(cowplot)
library(viridis)
library(ggridges)
library(GGally)
library(patchwork)
library(mgcv)

set.seed(42)

# ============================================================================
# SECTION 1: HELPER FUNCTIONS - DETECTABILITY FRAMEWORK
# ============================================================================

#' Define biological ground truth based on genealogical distance
#' 
#' @param dt data.table with genealogical information
#' @param gen_col Name of column containing generation distance
#' @param gen_threshold Maximum generations to consider "related"
#' @return data.table with new column 'true_transmission'
define_ground_truth <- function(dt, gen_col = "gen_distance", 
                                gen_threshold = 25) {
  message("  → Defining biological ground truth (gen_threshold = ", gen_threshold, ")")
  
  if (!gen_col %in% names(dt)) {
    stop("Column '", gen_col, "' not found in data")
  }
  
  dt[, true_transmission := as.integer(get(gen_col) <= gen_threshold)]
  
  n_pos <- sum(dt$true_transmission == 1, na.rm = TRUE)
  n_total <- nrow(dt)
  
  message("    ✓ True positives: ", n_pos, " / ", n_total, 
          " (", round(100 * n_pos / n_total, 2), "%)")
  
  return(dt)
}

#' Define detectability threshold based on false positive rate
#' (For descriptive statistics only - NOT used for filtering)
#' 
#' @param dt data.table with IBD proportions
#' @param ibd_col Name of column containing IBD proportion
#' @param truth_col Name of column containing ground truth
#' @param alpha False positive rate threshold (default: 1%)
#' @return Numeric threshold value
define_ibd_eps <- function(dt, ibd_col = "ibd_prop", 
                           truth_col = "true_transmission", 
                           alpha = 0.01) {
  message("  → Computing detectability threshold (alpha = ", alpha, ")\n")
  
  negatives <- dt[get(truth_col) == 0 & !is.na(get(ibd_col)), get(ibd_col)]
  positives <- dt[get(truth_col) == 1 & !is.na(get(ibd_col)), get(ibd_col)]
  
  if (length(negatives) == 0) {
    warning("    ⚠ No negatives found, using 0 as threshold")
    return(0)
  }
  
  # Diagnostic output
  message("    • Negative IBD: min = ", round(min(negatives, na.rm = TRUE), 4),
          ", median=", round(median(negatives, na.rm = TRUE), 4),
          ", 99th% = ", round(quantile(negatives, 0.99, na.rm = TRUE), 4),
          ", max = ", round(max(negatives, na.rm = TRUE), 4))
  
  if (length(positives) > 0) {
    message("    • Positive IBD: min = ", round(min(positives, na.rm = TRUE), 4),
            ", median = ", round(median(positives, na.rm = TRUE), 4),
            ", max = ", round(max(positives, na.rm = TRUE), 4))
  }
  
  # Compute threshold from negatives
  threshold_neg <- quantile(negatives, probs = 1 - alpha, na.rm = TRUE)
  
  # Also check 5th percentile of positives as a sanity check
  if (length(positives) > 0) {
    threshold_pos_5th <- quantile(positives, probs = 0.05, na.rm = TRUE)
    
    # If negatives' 99th percentile exceeds positives' 5th percentile,
    # there's no separation → use the minimum
    if (threshold_neg > threshold_pos_5th) {
      warning("    ⚠ Poor separation: 99th% of negatives (", round(threshold_neg, 4),
              ") > 5th% of positives (", round(threshold_pos_5th, 4), ")")
      message("    → Using 5th percentile of positives as threshold")
      ibd_eps <- threshold_pos_5th
    } else {
      ibd_eps <- threshold_neg
    }
  } else {
    ibd_eps <- threshold_neg
  }
  
  # Cap at reasonable maximum (90% IBD)
  if (ibd_eps > 0.9) {
    warning("    ⚠ Computed threshold (", round(ibd_eps, 4), 
            ") > 0.9, capping at 0.1")
    ibd_eps <- 0.2
  }
  message("  → Estimated detectability threshold: ", threshold_neg, "\n")
  message("  → Chosen detectability threshold: ", ibd_eps, "\n")
  message("    ✓ Detectability threshold: ", round(ibd_eps, 6), 
          " (", round(100 * (1 - alpha), 1), "th percentile)")
  
  return(as.numeric(ibd_eps))
}

#' Calibrate classification thresholds using ROC analysis
#' 
#' @param dt data.table with scores and truth
#' @param score_cols Vector of score column names
#' @param truth_col Name of truth column
#' @param method Threshold selection method ("youden" or "f1")
#' @return Named list of thresholds
calibrate_thresholds <- function(dt, score_cols, 
                                 truth_col = "identifiable_truth", 
                                 method = "youden") 
{
  message("  → Calibrating thresholds (method = ", method, ")")
  
  thresholds <- list()
  
  for (sc in score_cols) {
    if (!sc %in% names(dt)) {
      warning(" ⚠ Score column '", sc, "' not found, skipping")
      next
    }
    
    scores <- dt[[sc]]
    truth <- dt[[truth_col]]
    
    # Remove NAs
    valid <- !is.na(scores) & !is.na(truth)
    scores <- scores[valid]
    truth <- truth[valid]
    
    if (length(unique(truth)) < 2) {
      warning(" ⚠ Only one class present for '", sc, "', skipping")
      next
    }
    
    roc <- tryCatch(
      pROC::roc(truth, scores, quiet = TRUE),
      error = function(e) {
        warning(" ⚠ ROC failed for '", sc, "': ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(roc)) {
      thresholds[[sc]] <- median(scores, na.rm = TRUE)
      next
    }
    
    coords <- pROC::coords(roc, x = "best", 
                           best.method = ifelse(method == "youden", "youden", "closest.topleft"))
    
    thresholds[[sc]] <- coords$threshold
    
    message(" ✓ ", sc, ": threshold = ", round(coords$threshold, 4),
            " (sens = ", round(coords$sensitivity, 3), 
            ", spec = ", round(coords$specificity, 3), ")")
  }
  
  return(thresholds)
}

# ============================================================================
# SECTION 2: UTILITY FUNCTIONS
# ============================================================================

#' Canonical pair key (ensures id1 <= id2)
canonical_pair <- function(a,b) {
  if (is.na(a) || is.na(b)) return(NA_character_)
  if (a <= b) paste(a,b,sep="--") else paste(b,a,sep="--")
}

#' Safe file reading with error handling
safe_fread <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(fread(fp), error = function(e) { 
    warning("fread failed: ", fp, " - ", e$message)
    NULL 
  })
}

safe_readRDS <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(readRDS(fp), error = function(e) { 
    warning("readRDS failed: ", fp, " - ", e$message)
    NULL 
  })
}

#' Min-max normalization to [0, 1]
norm_to_01 <- function(x) {
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

#' Rank normalization (CORRECTED - removed minus sign)
rank_norm <- function(x) {
  if (all(is.na(x))) return(x)
  # Higher similarity -> higher rank -> higher normalized score
  r <- rank(x, ties.method = "average", na.last = "keep")
  (r - 1) / (max(r, na.rm = TRUE) - 1)
}


# ============================================================================
# SECTION 3: DATA LOADING FUNCTIONS
# ============================================================================

#' Load true IBD summary
load_true_ibd_summary <- function(rep_dir, rate_str) {
  
  repnum <- gsub("rep", "", basename(rep_dir))
  candidate <- file.path(rep_dir, paste0("run", repnum, "_rec", rate_str, "_chr1"))
  
  if (!dir.exists(candidate)) {
    cand <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
    cand <- cand[grepl(paste0("rec", rate_str), cand)]
    if (length(cand) >= 1) candidate <- cand[1] else candidate <- rep_dir
  }
  
  fp <- file.path(candidate, "true_ibd_summary.tsv")
  dt <- safe_fread(fp)
  
  if (is.null(dt)) return(NULL)
  
  # Modify sample IDs
  dt <- dt %>% 
    mutate(Id1 = paste0("tsk_", Id1), Id2 = paste0("tsk_", Id2)) %>% 
    as.data.table()
  
  # Standardize columns
  id_cols <- intersect(c("id1", "id2", "Id1", "Id2", "sample1", "sample2"), names(dt))
  if (length(id_cols) >= 2) setnames(dt, old = id_cols[1:2], new = c("id1", "id2"), skip_absent = TRUE)
  
  if (!("pair_key" %in% names(dt))) dt[, pair_key := mapply(canonical_pair, id1, id2)]
  if (!("total_ibd_bp" %in% names(dt))) dt[, total_ibd_bp := NA_real_]
  
  if (!("true_ibd_prop" %in% names(dt))) dt[, true_ibd_prop := total_ibd_bp / GENOME_BP]
  
  if (!("ibd_prop" %in% names(dt))) dt[, ibd_prop := true_ibd_prop]
  
  # Generation distance
  if (!("gen_distance" %in% names(dt))) {
    possible_gens <- intersect(c("generations", "t_generations", "min_tmrca"), names(dt))
    if (length(possible_gens) >= 1) setnames(dt, old = possible_gens[1], new = "gen_distance", skip_absent = TRUE)
  }
  
  dt[, .(pair_key, id1, id2, total_ibd_bp, ibd_prop, 
         max_seg_bp = ifelse("max_segment_bp" %in% names(dt), max_segment_bp, NA_real_), 
         gen_distance)]
}

# load inferred HMM-IBD predictions
load_inferred_ibd <- function(rep_dir, rate_str) {
  
  repnum <- gsub("rep","", basename(rep_dir))
  candidate <- file.path(rep_dir, paste0("run", repnum, "_rec", rate_str, "_chr1"))
  
  if (!dir.exists(candidate)) {
    cand <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
    cand <- cand[grepl(paste0("rec", rate_str), cand)]
    if (length(cand) >= 1) candidate <- cand[1] else candidate <- rep_dir
  }
  
  fp <- file.path(candidate, "inferred_ibd_hmm.rds")
  dt <- as.data.table(safe_readRDS(fp))
  
  if (is.null(dt)) return(NULL)
  
  id_cols <- intersect(c("id1","id2","Id1","Id2","sample1","sample2"), names(dt))
  
  if (length(id_cols) >= 2) setnames(dt, old = id_cols[1:2], new = c("id1","id2"), skip_absent = TRUE)
  # pick score column
  score_candidates <- intersect(c("total_ibd_bp","total_ibd","score","ibd","ibd_prop","hmm","posterior"), names(dt))
  
  if (length(score_candidates) > 0) dt[, score := as.numeric(get(score_candidates[1]))] else {
    if ("n_segments" %in% names(dt)) dt[, score := as.numeric(n_segments)] else dt[, score := 1]
  }
  
  if (!("pair_key" %in% names(dt))) dt[, pair_key := mapply(canonical_pair, id1, id2)]
  unique(dt[, .(pair_key, id1, id2, score)])
}

# load IBS matrix in RDS to long form
load_ibs_matrix <- function(rep_dir, rate_str) {
  
  repnum <- gsub("rep","", basename(rep_dir))
  candidate <- file.path(rep_dir, paste0("run", repnum, "_rec", rate_str, "_chr1"))
  
  if (!dir.exists(candidate)) {
    cand <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
    cand <- cand[grepl(paste0("rec", rate_str), cand)]
    if (length(cand) >= 1) candidate <- cand[1] else candidate <- rep_dir
  }
  
  fp <- file.path(candidate, "ibs_matrix.rds")
  
  if (!file.exists(fp)) return(NULL)
  
  obj <- tryCatch(readRDS(fp), error = function(e) { 
    warning("readRDS failed: ", fp) 
    NULL 
  })
  
  if (is.null(obj)) return(NULL)
  mat <- as.matrix(obj)
  ids <- rownames(mat)
  
  if (is.null(ids)) {
    # maybe first column contains ids
    if (is.data.frame(obj)) {
      ids <- as.character(obj[[1]])
      mat <- as.matrix(obj[, -1])
      rownames(mat) <- ids; colnames(mat) <- ids
    } else return(NULL)
  }
  
  n <- length(ids)
  rows <- vector("list", n*(n-1)/2)
  idx <- 1
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      rows[[idx]] <- list(
        pair_key = canonical_pair(ids[i], ids[j]), 
        id1 = ids[i], 
        id2 = ids[j], 
        score = as.numeric(mat[i,j])
      )
      idx <- idx + 1
    }
  }
  dt <- rbindlist(rows)

  # If IBS is a distance (smaller=more related), 
  # we may invert later; here keep raw score and mark method
  unique(dt[, .(pair_key, id1, id2, score)])
}

# load patristic distances using phylo_long helper (must exist)
load_phylo_patristic <- function(rep_phylo_dir, rate_str) {
  
  # find treefile matching rec rate pattern
  patt <- paste0("rec", rate_str)
  files <- list.files(rep_phylo_dir, full.names = TRUE, pattern = patt)
  
  # pick modelfinder one if available
  tf <- files[grepl("_chr1.treefile$", files)]
  
  if (length(tf) == 0 && length(files)>0) tf <- files[1]
  if (length(tf)==0) return(NULL)
  
  # expect function phylo_long() available
  if (!exists("phylo_long")) {
    helper <- file.path("simulations", "R", "patristic_distances.R")
    if (file.exists(helper)) source(helper) else {
      warning("phylo_long not found; skipping patristic for ", rep_phylo_dir)
      return(NULL)
    }
  }
  
  dt <- as.data.table(phylo_long(tf[1], method = "patristic")) # Distances are converted to similarity
  id_cols <- intersect(c("id1","id2","Id1","Id2","sample1","sample2"), names(dt))
  
  if (length(id_cols) >= 2) setnames(dt, old = id_cols[1:2], new = c("id1","id2"), skip_absent = TRUE)
  
  # choose distance col and convert to similarity (higher = more related)
  if ("phylo" %in% names(dt)) dt[, score := phylo] else if ("distance" %in% names(dt)) dt[, score := distance] else dt[, score := 0]
  if (!("pair_key" %in% names(dt))) dt[, pair_key := mapply(canonical_pair, id1, id2)]
  unique(dt[, .(pair_key, id1, id2, score)])
}

# ============================================================================
# SECTION 4: EVALUATION FUNCTION
# ============================================================================

#' Evaluate a single method on detectability-conditioned data
evaluate_method_on_table <- function(
    dt_merge,
    score_col = "score",
    truth_col = "true_transmission",  # Changed default to genealogical truth
    threshold = NULL,
    ks = c(5, 10, 25, 50),
    invert_dist  = FALSE    # assumes higher == more related; invert if score is a distance
) 
{
  res <- list()
  
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
  
  # nothing left after filtering
  if (nrow(dt) == 0) return(res)
  
  # ----------------------------
  # HANDLE DISTANCE OR SIMILARITY
  # ----------------------------
  s <- dt[[score_col]]
  if (invert_dist) {
    # higher value should mean *more related*
    s_use <- -s
  } else {
    s_use <- s
  }
  dt[, score_internal := s_use]
  
  # ----------------------------
  # AUPR
  # ----------------------------
  pos_scores <- dt[get(truth_col) == 1, score_internal]
  neg_scores <- dt[get(truth_col) == 0, score_internal]
  
  if (length(pos_scores) > 0 && length(neg_scores) > 0) {
    pr <- PRROC::pr.curve(scores.class0 = pos_scores,
                          scores.class1 = neg_scores,
                          curve = TRUE)
    
    res$pr <- pr
    res$aupr <- pr$auc.integral
  } else {
    res$pr <- NA_real_
    res$aupr <- NA_real_
  }
  
  # ----------------------------
  # AUROC
  # ----------------------------
  if (length(unique(dt[[truth_col]])) > 1) {
    roc_obj <- tryCatch(
      pROC::roc(dt[[truth_col]], dt$score_internal, quiet = TRUE),
      error = function(e) NULL
    )
    res$roc <- if (!is.null(roc_obj)) roc_obj else NA_real_
    res$auroc <- if (!is.null(roc_obj)) pROC::auc(roc_obj) else NA_real_
  } else {
    res$roc <- NA_real_
    res$auroc <- NA_real_
  }
  
  # ----------------------------
  # Spearman correlation vs continuous truth
  # ----------------------------
  if ("ibd_prop" %in% names(dt)) {
    res$spearman <- suppressWarnings(
      cor(dt$score_internal, dt$ibd_prop, method = "spearman", use = "complete.obs")
    )
  } else {
    res$spearman <- NA_real_
  }
  
  # ----------------------------
  # Brier score (scaled to 0–1)
  # ----------------------------
  rng <- range(dt$score_internal, na.rm = TRUE)
  if (diff(rng) == 0) {
    pred_prob <- rep(0.5, nrow(dt))
  } else {
    pred_prob <- (dt$score_internal - rng[1]) / diff(rng)
  }
  res$brier <- mean((pred_prob - dt[[truth_col]])^2, na.rm = TRUE)
  
  # ----------------------------
  # Precision@K
  # ----------------------------
  if (invert_dist) {
    # similarity-like; sort descending
    ord <- dt[order(-score_internal)]
  } else {
    # distance-like; sort ascending
    ord <- dt[order(score_internal)]
  }
  
  for (k in ks) {
    topk <- head(ord, k)
    res[[paste0("prec_at_", k)]] <-
      if (nrow(topk) == 0) NA_real_ else mean(topk[[truth_col]] == 1)
  }
  
  # ----------------------------
  # CONFUSION MATRIX at 0 threshold
  # ----------------------------
  # predicted positive = score_internal > 0 (or use median, but 0 is fine)
  # pred <- ifelse(dt$score_internal > 0, 1L, 0L)
  
  # Arbitrary threshold: The median is data-dependent and has no biological justification
  # Class imbalance ignored: If 90% of pairs are negatives, median threshold produces ~50% positives (wrong)
  # Metric interpretation: Precision/recall at median threshold don't reflect real-world performance
  # Biological Impact: F1-scores are meaningless when the threshold doesn't respect the true positive rate in the population.
  
  # Use Youden's index (optimal threshold from ROC curve)
  if (!is.null(roc_obj)) {
    coords <- pROC::coords(roc_obj, "best", best.method = "youden")
    thr <- coords$threshold
  } else {
    thr <- median(dt$score_internal, na.rm = TRUE)
  }
  
  message("    ✓ ", score_col, ": threshold = ", 
          round(coords$threshold, 4),
          " (sens = ", round(coords$sensitivity, 3), 
          ", spec = ", round(coords$specificity, 3), ")")

  # thr <- median(dt$score_internal, na.rm = TRUE) # Median
  pred <- as.integer(dt$score_internal >= thr)
  
  # # Confusion matrix at calibrated threshold
  # if (!is.null(threshold)) {
  #   if (invert_dist) {
  #     pred <- as.integer(dt$score_internal <= -threshold)  # Note: negated
  #   } else {
  #     pred <- as.integer(dt$score_internal >= threshold)
  #   }
  
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
  
    return(res)
}

# ============================================================================
# SECTION 5: MAIN ANALYSIS PIPELINE
# ============================================================================

run_recombination_analysis <- function(
    ROOT_DIR = "simulations/multiple_runs/metrics",
    INFERRED_SUBPATH = "inferred",
    PHYLO_ROOT = file.path(ROOT_DIR, "phylo_results"),
    REC_RATES = c("1e-09","1e-08","1e-07","1e-06"),
    REPLICATES = NULL,            # NULL => auto-detect rep* under INFERRED_SUBPATH
    GEN_THRESHOLD = 25,           # binary truth: gen_distance <= GEN_THRESHOLD => positive
    TOP_K = c(5,10,25,50),
    GENOME_BP = 640851,      # for ibd_prop computation fallback
    OUTDIR = file.path(ROOT_DIR, "recombination_evaluation"),
    # OUTDIR = file.path(ROOT_DIR, "test"),
    SAVE_PLOTS = TRUE,
    ALPHA_DETECT = 0.01  # False positive rate for detectability threshold
) 
{
  
  message("\n" , rep("=", 60), "\n",
          "\tIBD METHOD EVALUATION PIPELINE \n",
          rep("=", 60))
  
  # Setup directories
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(OUTDIR,"tables"), showWarnings = FALSE)
  dir.create(file.path(OUTDIR,"figures"), showWarnings = FALSE)
  
  INFERRED_ROOT <- file.path(ROOT_DIR, INFERRED_SUBPATH)
  
  # auto-detect replicates as directories under INFERRED_ROOT
  if (is.null(REPLICATES)) {
    REPLICATES <- list.dirs(INFERRED_ROOT, recursive = FALSE, full.names = FALSE)
    REPLICATES <- REPLICATES[grepl("^rep", REPLICATES)]
  }
  
  message("\n[CONFIG]")
  message("  • Replicates: ", paste(REPLICATES, collapse = ", "))
  message("  • Recombination rates: ", paste(REC_RATES, collapse = ", "))
  message("  • Generation threshold: ", GEN_THRESHOLD)
  message("  • Detectability alpha: ", ALPHA_DETECT)
  message("  • Output directory: ", OUTDIR)
  
  # Create rate lookup table (CORRECTED)
  REC_RATES_LABEL <- gsub("-", "", REC_RATES)
  rate_lookup <- setNames(as.numeric(REC_RATES), REC_RATES_LABEL)
  
  message("\n[RATE LOOKUP TABLE]")
  for (i in seq_along(REC_RATES_LABEL)) {
    message("  ", REC_RATES_LABEL[i], " -> ", rate_lookup[REC_RATES_LABEL[i]])
  }
  
  
  # Storage
  all_metrics_rows <- list()
  merged_scores_store <- list()
  all_performance <- list()
  all_prcurve <- list()
  detectability_rows <- list()
  
  # Global variables for calibration
  global_thresholds <- NULL
  ibd_eps_global <- NULL
  
  # ============================================================================
  # MAIN LOOP: Per-replicate, per-rate evaluation
  # ============================================================================
  
  message("\n", rep("=", 60))
  message("\tSTEP 1: DATA LOADING AND EVALUATION")
  message(rep("=", 60))
  
  # --------- Step 1: per-replicate, per-rate evaluation -------------
  for (rep_idx in seq_along(REPLICATES)) {
    
    pr_curve <- NULL
    rep <- REPLICATES[rep_idx]
    rep_dir <- file.path(INFERRED_ROOT, rep)
    
    message("\n[", rep_idx, "/", length(REPLICATES), "] Processing replicate: ", rep)
    
    if (!dir.exists(rep_dir)) {
      warning("  ⚠ Replicate dir missing: ", rep_dir)
      next
    }
    
    for (rate_idx in seq_along(REC_RATES_LABEL)) {
      rate <- REC_RATES_LABEL[rate_idx]
      rate_numeric <- rate_lookup[rate]
      is_baseline <- (rate_idx == 1)
      
      message("\n  [Rate ", rate_idx, "/", 
              length(REC_RATES_LABEL), "] ", 
              rate, " (", rate_numeric, " bp⁻¹)", 
              ifelse(is_baseline, " [BASELINE]", ""))
      
      # -------------------------------------------------------------------
      # Load ground truth
      # -------------------------------------------------------------------
      message("  Loading ground truth...")
      truth_dt <- load_true_ibd_summary(rep_dir, rate)
      
      if (is.null(truth_dt) || nrow(truth_dt) == 0) {
        warning("    ⚠ Missing truth data, skipping")
        next
      }
      
      message("    ✓ Loaded ", nrow(truth_dt), " pairs")
      
      # -------------------------------------------------------------------
      # Define biological truth (NO detectability conditioning)
      # -------------------------------------------------------------------
      truth_dt <- define_ground_truth(truth_dt, gen_threshold = GEN_THRESHOLD)
      
      if (is_baseline) {
        ibd_eps_global <- define_ibd_eps(truth_dt, alpha = ALPHA_DETECT)
      }
      
      # Mark detectable pairs (for descriptive statistics only)
      truth_dt[, detectable := as.integer(ibd_prop >= ibd_eps_global)]
      
      # Store detectability stats
      n_true_pos <- sum(truth_dt$true_transmission == 1, na.rm = TRUE)
      n_detectable <- sum(truth_dt$true_transmission == 1 & truth_dt$detectable == 1, na.rm = TRUE)
      detectability_rate <- if (n_true_pos > 0) n_detectable / n_true_pos else 0
      
      message("  → Detectability statistics (for reporting only):")
      message("    • True positives: ", n_true_pos)
      message("    • With detectable IBD (≥", round(ibd_eps_global, 4), "): ", 
              n_detectable, " (", round(100 * detectability_rate, 1), "%)")
      message("    ⓘ Evaluation uses ALL positives (genealogical truth)")
      
      detectability_rows[[length(detectability_rows) + 1]] <- data.table(
        replicate = rep,
        rate = rate_numeric,
        n_pairs = nrow(truth_dt),
        n_true_positives = n_true_pos,
        n_detectable = n_detectable,
        detectability = detectability_rate,
        ibd_eps = ibd_eps_global
      )
      
      # -------------------------------------------------------------------
      # Load methods
      # -------------------------------------------------------------------
      message("  Loading method predictions...")
      
      ibd_hmm_dt <- load_inferred_ibd(rep_dir, rate)
      if (!is.null(ibd_hmm_dt)) 
        message(" ✓ HMM-IBD: ", nrow(ibd_hmm_dt), " pairs")
      
      ibs_dt <- load_ibs_matrix(rep_dir, rate)
      if (!is.null(ibs_dt)) message(" ✓ IBS: ", nrow(ibs_dt), " pairs")
      
      phylo_dt <- NULL
      phylo_rep_dir <- file.path(PHYLO_ROOT, rep)
      if (dir.exists(phylo_rep_dir)) {
        phylo_dt <- load_phylo_patristic(phylo_rep_dir, rate)
        if (!is.null(phylo_dt)) 
          message(" ✓ Phylo: ", nrow(phylo_dt), " pairs")
      }
      
      # -------------------------------------------------------------------
      # Merge methods with truth
      # -------------------------------------------------------------------
      message("  Merging data...")

      method_tables <- list()
      if (!is.null(ibd_hmm_dt)) method_tables$IBD <- ibd_hmm_dt
      if (!is.null(ibs_dt)) method_tables$IBS <- ibs_dt
      if (!is.null(phylo_dt)) method_tables$Phylo <- phylo_dt
      
      # base merged table (truth)
      merged <- unique(truth_dt[, .(pair_key, id1, id2, ibd_prop, 
                                    true_transmission, detectable, gen_distance)])
      
      for (m in names(method_tables)) {
        mt <- unique(method_tables[[m]][, .(pair_key, score)])
        setnames(mt, "score", paste0("score_", m))
        merged <- merge(merged, mt, by = "pair_key", all.x = TRUE)
      }
      
      message("    ✓ Merged table: ", 
              nrow(merged), " pairs × ", 
              ncol(merged), " columns")
      
      # Store merged table
      merged_scores_store[[paste(rep, rate, sep = "__")]] <- merged

      # -------------------------------------------------------------------
      # Evaluate each method
      # -------------------------------------------------------------------
      message("  Evaluating methods...")
      for (m in names(method_tables)) {
        
        score_col <- paste0("score_", m)
        eval_res <- evaluate_method_on_table(
          merged, 
          score_col = score_col, 
          truth_col = "true_transmission", 
          ks = TOP_K
        )
        
        # Store metrics
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
          n_true_positives = n_true_pos
        )
        
        for (k in TOP_K) {
          row[[paste0("prec_at_", k)]] <- eval_res[[paste0("prec_at_", k)]]
        }
        
        all_metrics_rows[[length(all_metrics_rows) + 1]] <- row
        
        # Store performance metrics
        performance <- data.table(
          replicate = rep,
          method = m,
          rate = rate_numeric,
          TP = eval_res$TP,
          FP = eval_res$FP,
          TN = eval_res$TN,
          FN = eval_res$FN,
          Precision = eval_res$precision,
          Recall = eval_res$recall,
          f1 = eval_res$f1
        )
        
        all_performance[[length(all_performance) + 1]] <- performance
        
        # Store PR curve
        if (!is.null(eval_res$pr)) { #  && inherits(eval_res$pr, "precrec")
          pr_curve <- rbind.data.frame(pr_curve, data.table(
            replicate = rep,
            method = m,
            rate = rate_numeric,
            eval_res$pr$curve
          ))
          # names(pr_curve)[4:6] <- c("recall", "precision", "threshold")
        }
        
        message("    ✓ ", m, ": AUPR = ", round(eval_res$aupr, 3),
                ", AUROC = ", round(eval_res$auroc, 3),
                ", Prec@10 = ", round(eval_res$prec_at_10, 3))
      }
    }
    all_prcurve[[rep]] <- pr_curve
    message("\n  ✓ Completed replicate ", rep)
  } 
  
  # end replicates loop
  
  # ============================================================================
  # STEP 2: AGGREGATE RESULTS
  # ============================================================================
  
  message("\n", rep("=", 60))
  message("\tSTEP 2: AGGREGATING RESULTS ACROSS REPLICATES")
  message(rep("=", 60))
  
  metrics_dt <- rbindlist(all_metrics_rows, fill = TRUE)
  performance_dt <- rbindlist(all_performance, fill = TRUE)
  prcurve_dt <- rbindlist(all_prcurve, fill = TRUE)
  detectability_dt <- rbindlist(detectability_rows, fill = TRUE)
  
  message("\n  • Total evaluations: ", nrow(metrics_dt))
  message("  • Methods: ", paste(unique(metrics_dt$method), collapse = ", "))
  message("  • Rates: ", paste(unique(metrics_dt$rate), collapse = ", "))
  
  # Compute summary statistics
  agg_dt <- metrics_dt[, .(
    aupr_median = median(aupr, na.rm = TRUE),
    aupr_mean = mean(aupr, na.rm = TRUE),
    aupr_sd = sd(aupr, na.rm = TRUE),
    aupr_se = sd(aupr, na.rm = TRUE) / sqrt(.N),
    auroc_median = median(auroc, na.rm = TRUE),
    auroc_mean = mean(auroc, na.rm = TRUE),
    auroc_sd = sd(auroc, na.rm = TRUE),
    auroc_se = sd(auroc, na.rm = TRUE) / sqrt(.N),
    spearman_median = median(spearman, na.rm = TRUE),
    spearman_mean = mean(spearman, na.rm = TRUE),
    spearman_sd = sd(spearman, na.rm = TRUE),
    detectability_mean = mean(detectability, na.rm = TRUE),
    n_replicates = .N
  ), by = .(rate, method)]
  
  # Compute detectability summary
  detect_agg <- detectability_dt[, .(
    detectability_mean = mean(detectability, na.rm = TRUE),
    detectability_sd = sd(detectability, na.rm = TRUE),
    detectability_se = sd(detectability, na.rm = TRUE) / sqrt(.N),
    ibd_eps_mean = mean(ibd_eps, na.rm = TRUE)
  ), by = rate]
  
  # Save tables
  message("\n  Saving tables...")
  write_tsv(metrics_dt, file.path(OUTDIR, "tables", "metrics_by_rate_replicate_method.tsv"))
  writexl::write_xlsx(metrics_dt, file.path(OUTDIR, "tables", "metrics_by_rate_replicate_method.xlsx"))
  
  writexl::write_xlsx(agg_dt, file.path(OUTDIR, "tables", "agg_metrics_by_rate_method.xlsx"))
  writexl::write_xlsx(performance_dt, file.path(OUTDIR, "tables", "evaluation_metrics.xlsx"))
  writexl::write_xlsx(detectability_dt, file.path(OUTDIR, "tables", "detectability_by_rate.xlsx"))
  
  write_tsv(prcurve_dt, file.path(OUTDIR, "tables", "pr_curve.tsv"))
  
  message("    ✓ All tables saved")
  
  # ============================================================================
  # STEP 3: VISUALIZATION
  # ============================================================================
  message("\n", rep("=", 60))
  message("\t STEP 3: GENERATING VISUALIZATIONS")
  message(rep("=", 60))
  
  theme_set(theme_bw(base_size = 14)) # , base_family = "Arial"
  
  custom_theme <- theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black", face = "bold"),
    # plot.title = element_text(color = "black", face = "bold", size = 15, hjust = 0.5),
    # plot.subtitle = element_text(color = "gray30", size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm"),
    strip.text = element_text(size = 14, color = "black", face = "bold"),
    panel.grid.minor = element_blank()
  )
  
  # -------------------------------------------------------------------
  # Plot 1: AUPR vs Recombination Rate
  # -------------------------------------------------------------------
  message("\n  [1/7] AUPR vs recombination rate...")
  
  # p_aupr <- ggplot(agg_dt, aes(x = rate, y = aupr_mean, color = method, group = method)) +
  #   geom_line(linewidth = 1.2) +
  #   geom_point(size = 3) +
  #   geom_errorbar(aes(ymin = aupr_mean - aupr_se, ymax = aupr_mean + aupr_se), 
  #                 width = 0.1, linewidth = 0.8) +
  #   scale_x_log10(
  #     breaks = as.numeric(REC_RATES),
  #     labels = REC_RATES
  #   ) +
  #   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  #   labs(
  #     title = "Method Performance vs Recombination Rate",
  #     subtitle = "AUPR (Area Under Precision-Recall Curve)",
  #     x = "Recombination Rate (bp⁻¹, log scale)",
  #     y = "AUPR (mean ± SE)"
  #   ) +
  #   custom_theme
  # 
  # if (SAVE_PLOTS) {
  #   ggsave(file.path(OUTDIR, "figures", "AUPR_vs_recrate.png"), p_aupr, 
  #          width = 10, height = 7, dpi = 600)
  #   message("    ✓ Saved")
  # }  
  # p_auroc <- ggplot(agg_dt, aes(x = rate_num, y = auroc_mean, color = method, group = method)) +
  #   geom_line(linewidth = 1.2) + 
  #   geom_point() +
  #   geom_errorbar(aes(ymin = auroc_mean - auroc_sd, ymax = auroc_mean + auroc_sd), width = 0.08) +
  #   scale_x_log10() +
  #   labs(x = "", y = "AUROC (mean ± sd)", title = "AUROC vs recombination rate") +
  #   custom
  # 
  # # if (SAVE_PLOTS) ggsave(file.path(OUTDIR,"figures","AUROC_vs_recrate.pdf"), p_auroc, width = 8, height = 5)
  # if (SAVE_PLOTS) ggsave(file.path(OUTDIR,"figures","AUROC_vs_recrate.png"), p_auroc, width = 12, height = 8)
  
  # -------------------------------------------------------------------
  # Plot 2: Detectability vs Recombination Rate
  # -------------------------------------------------------------------
  message("  [2/7] Detectability vs recombination rate...")
  
  # p_detect <- ggplot(detect_agg, aes(x = rate, y = detectability_mean)) +
  #   geom_line(linewidth = 1.2, color = "darkred") +
  #   geom_point(size = 3, color = "darkred") +
  #   geom_errorbar(aes(ymin = detectability_mean - detectability_se,
  #                     ymax = detectability_mean + detectability_se),
  #                 width = 0.1, linewidth = 0.8, color = "darkred") +
  #   scale_x_log10(
  #     breaks = as.numeric(REC_RATES),
  #     labels = REC_RATES
  #   ) +
  #   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
  #                      labels = scales::percent_format()) +
  #   labs(
  #     title = "IBD Detectability vs Recombination Rate",
  #     subtitle = paste0("Proportion of related pairs (≤", GEN_THRESHOLD, 
  #                       " generations) with detectable IBD"),
  #     x = "Recombination Rate (bp⁻¹, log scale)",
  #     y = "Detectability (% of true positives)"
  #   ) +
  #   custom_theme +
  #   theme(legend.position = "none")
  # 
  # if (SAVE_PLOTS) {
  #   ggsave(file.path(OUTDIR, "figures", "Detectability_vs_recrate.png"), p_detect,
  #          width = 10, height = 7, dpi = 600)
  #   message("    ✓ Saved")
  # }
  
  # -------------------------------------------------------------------
  # Plot 3: Combined AUPR + Detectability
  # -------------------------------------------------------------------
  message("  [3/7] Combined AUPR and detectability...")
  
  # Normalize detectability to same scale as AUPR for dual-axis
  p_combined <- ggplot() +
    geom_line(data = agg_dt, aes(x = rate, y = aupr_mean, color = method, group = method),
              linewidth = 1.2) +
    geom_point(data = agg_dt, aes(x = rate, y = aupr_mean, color = method),
               size = 3) +
    geom_errorbar(data = agg_dt,
      aes(x = rate,
        ymin = aupr_mean - aupr_se,
        ymax = aupr_mean + aupr_se,
        color = method),
      width = 0.1, linewidth = 0.8) +
  
    geom_line(data = detect_agg, aes(x = rate, y = detectability_mean),
              linetype = "dashed", linewidth = 1.5, color = "black") +
    geom_point(data = detect_agg, aes(x = rate, y = detectability_mean),
               size = 4, shape = 18, color = "black") +
    scale_x_log10(
      breaks = as.numeric(REC_RATES),
      labels = REC_RATES) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Method Performance and IBD Detectability",
      # subtitle = "Solid lines = AUPR by method | Dashed line = % detectable IBD",
      x = "",
      y = "AUPR / Detectability") +
    custom_theme + 
    theme(plot.title = element_text(color = "gray50", face = "bold", size = 15),
          plot.title.position = "plot") +
    annotate("text", x = max(detect_agg$rate) * 0.8, 
             y = min(detect_agg$detectability_mean) + 0.03,
             label = "Detectability -->", hjust = 1, size = 4, fontface = "bold")
  
  if (SAVE_PLOTS) {
    ggsave(file.path(OUTDIR, "figures", "AUPR_Detectability_combined.png"), p_combined,
           width = 12, height = 7, dpi = 600)
    message("    ✓ Saved")
  }
  
  # -------------------------------------------------------------------
  # Plot 4: Spearman Correlation
  # -------------------------------------------------------------------
  message("  [4/7] Spearman correlation...")
  
  p_spear <- ggplot(agg_dt, aes(x = rate, y = spearman_mean, color = method, group = method)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = spearman_mean - spearman_sd, 
                      ymax = spearman_mean + spearman_sd),
                  width = 0.1, linewidth = 0.8) +
    scale_x_log10(breaks = as.numeric(REC_RATES), labels = REC_RATES) +
    labs(
      title = "Spearman Correlation vs Recombination Rate",
      # subtitle = "Correlation between predicted score and true IBD proportion",
      x = "",
      y = "Spearman ρ (mean ± SD)") + 
    custom_theme + 
    theme(plot.title = element_text(color = "gray50", face = "bold", size = 15),
          plot.title.position = "plot") +
  
  if (SAVE_PLOTS) {
    ggsave(file.path(OUTDIR, "figures", "Spearman_vs_recrate.png"), 
           p_spear, width = 10, height = 7, dpi = 600)
    message("    ✓ Saved")
  }
  
  # -------------------------------------------------------------------
  # Plot 5: Precision@K
  # -------------------------------------------------------------------
  message("  [5/7] Precision@K...")
  
  preck_cols <- grep("^prec_at_", names(metrics_dt), value = TRUE)
  prec_dt <- reshape2::melt(metrics_dt, 
                  id.vars = c("replicate", "rate", "method"),
                  measure.vars = preck_cols,
                  variable.name = "metric", 
                  value.name = "value")
  
  prec_dt <- as.data.table(prec_dt)
  prec_dt[, K := as.integer(sub("^prec_at_", "", metric))]
  
  agg_prec <- prec_dt[, .(
    mean_prec = mean(value, na.rm = TRUE),
    sd_prec = sd(value, na.rm = TRUE)
  ), by = .(rate, method, K)]
  
  p_precK <- ggplot(agg_prec, aes(x = rate, y = mean_prec, color = method, group = method)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    facet_wrap(~K, scales = "free_y", labeller = label_both) +
    scale_x_log10(
      breaks = as.numeric(REC_RATES),
      labels = REC_RATES) +
    labs(
      title = "Precision@K vs Recombination Rate",
      # subtitle = "Proportion of true positives in top K predictions",
      x = "Recombination Rate",
      y = "Precision@K (mean)") + custom_theme
  
  if (SAVE_PLOTS) {
    ggsave(file.path(OUTDIR, "figures", "Precision_at_K_vs_recrate.png"), 
           p_precK, width = 12, height = 8, dpi = 600)
    message("    ✓ Saved")
  }
  
  # -------------------------------------------------------------------
  # Plot 6: Delta AUPR Heatmap
  # -------------------------------------------------------------------
  message("  [6/7] Delta AUPR heatmap...")
  
  # Delta-AUPR heatmap (method x rate)
  baseline_rate <- as.numeric(REC_RATES[1])
  base <- agg_dt[rate == baseline_rate, .(method, base_aupr = aupr_mean)]
  delta <- merge(agg_dt, base, by = "method", all.x = TRUE)
  delta <- as.data.table(delta)
  delta[, delta_aupr := aupr_mean - base_aupr]
  
  heat_dt <- delta[, .(method, rate, delta_aupr)]
  heat_dt[, rate_factor := factor(rate, levels = as.numeric(REC_RATES), labels = REC_RATES)]
  
  p_heat_diverg <- ggplot(heat_dt, aes(x = rate_factor, y = method, fill = delta_aupr)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", delta_aupr)),
              color = ifelse(abs(heat_dt$delta_aupr) > 0.25, "white", "black"),
              size = 5, fontface = "bold") +
    scale_fill_gradient2(low = "darkred", mid = "gold", high = "#1a9850", # "#d73027"
      midpoint = 0, name = "Δ AUPR",
      limits = c(min(heat_dt$delta_aupr, na.rm = TRUE), 
                 max(heat_dt$delta_aupr, na.rm = TRUE))) +
    labs(
      title = "Change in AUPR Relative to Baseline",
      # subtitle = paste0("Baseline = ", REC_RATES[1], " bp⁻¹ | Green = improved, Red = degraded"),
      x = "Recombination Rate", y = "") +
    theme_minimal(base_size = 14) + 
    theme(
      axis.text = element_text(size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 14, color = "black", face = "bold", margin = margin(t = 10)),
      legend.text = element_text(size = 12, color = "black", face = "bold"),
      legend.title = element_text(size = 13, color = "black", face = "bold"),
      plot.title = element_text(size = 16, color = "gray50", face = "bold"),
      plot.title.position = "plot",
      panel.grid = element_blank())
  
  if (SAVE_PLOTS) {
    ggsave(file.path(OUTDIR, "figures", "DeltaAUPR_heatmap.png"), p_heat_diverg,
           width = 10, height = 6, dpi = 600)
    message("    ✓ Saved")
  }
  
  # -------------------------------------------------------------------
  # Composite Figure
  # -------------------------------------------------------------------
  message("\n  Creating composite figure...")
  
  if (SAVE_PLOTS) {
    composite <- (p_combined | p_spear ) / p_heat_diverg +  #(p_aupr | p_detect)
      plot_layout(heights = c(2, 1)) +
      plot_annotation(
        title = "Method Performance Under Increasing Recombination",
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))
    
    ggsave(file.path(OUTDIR, "figures", "Composite_main_figure.png"), composite,
           width = 16, height = 12, dpi = 600)
    message("    ✓ Saved composite figure")
  }


  # -------------------------------------------------------------------
  # Plot 7: IBD Decay by Generation
  # -------------------------------------------------------------------
  # message("  [7/7] IBD decay by generation...")
  # 
  # decay_rows <- list()
  # for (key in names(merged_scores_store)) {
  #   merged <- merged_scores_store[[key]]
  #   parts <- strsplit(key, "__")[[1]]
  #   repname <- parts[1]
  #   rrate <- parts[2]
  #   
  #   if (!("gen_distance" %in% names(merged))) next
  #   
  #   tmp <- unique(merged[, .(pair_key, gen_distance, ibd_prop)])
  #   tmp[, replicate := repname]
  #   tmp[, rate := rate_lookup[rrate]]
  #   tmp[, rate_label := rrate]
  #   decay_rows[[length(decay_rows) + 1]] <- tmp
  # }
  # 
  # if (length(decay_rows) > 0) {
  #   decay_dt <- rbindlist(decay_rows)
  #   
  #   # Downsample for plotting
  #   set.seed(123)
  #   decay_sub <- decay_dt[, .SD[sample(.N, min(20000, .N))], by = rate]
  #   
  #   # Bin and summarize
  #   decay_dt[, gen_bin := floor(gen_distance)]
  #   sum_dt <- decay_dt[, .(
  #     ibd_med = median(ibd_prop, na.rm = TRUE),
  #     ibd_lo = quantile(ibd_prop, 0.25, na.rm = TRUE),
  #     ibd_hi = quantile(ibd_prop, 0.75, na.rm = TRUE),
  #     n = .N
  #   ), by = .(rate, gen_bin)]
  #   
  #   sum_dt[, rate_label := factor(rate, 
  #                                 levels = as.numeric(REC_RATES),
  #                                 labels = REC_RATES)]
  #   
  #   p_decay <- ggplot(sum_dt, aes(x = gen_bin, y = ibd_med, color = rate_label)) +
  #     geom_ribbon(aes(ymin = ibd_lo, ymax = ibd_hi, fill = rate_label), 
  #                 alpha = 0.2, color = NA) +
  #     geom_line(linewidth = 1.2) +
  #     scale_color_viridis_d(name = "Recombination\nRate (bp⁻¹)") +
  #     scale_fill_viridis_d(name = "Recombination\nRate (bp⁻¹)") +
  #     labs(
  #       title = "IBD Decay with Genealogical Distance",
  #       # subtitle = "Median ± IQR across all replicates",
  #       x = "Generations Apart",
  #       y = "IBD Proportion"
  #     ) +
  #     custom_theme
  #   
  #   if (SAVE_PLOTS) {
  #     ggsave(file.path(OUTDIR, "figures", "IBD_decay_by_generation.png"), p_decay,
  #            width = 10, height = 7, dpi = 600)
  #     message("    ✓ Saved")
  #   }
  # }
  # 
  # # =====================================================
  # # Replace LOESS with GAM (mgcv::gam) — massively faster
  # library(mgcv)
  # p_decay_sub <- ggplot(decay_sub, aes(gen_distance, ibd_prop)) +
  #   geom_point(alpha = 0.3, size = 0.6) +
  #   geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp"), se = TRUE) +
  #   facet_wrap(~ rate) + # , scales = "free_y"
  #   labs(title = "Decay of IBD proportion with generation by recombination rate", 
  #        x = "Generations apart", y = "IBD proportion") +
  #   theme(legend.position = "bottom",
  #         legend.title = element_blank(),
  #         legend.text = element_text(size = 14, color = 'black', face = 'bold'),
  #         plot.title = element_text(color = 'black', face = 'bold'),
  #         # plot.title.position = "plot",
  #         axis.title = element_text(size = 16, color = 'black', face = 'bold'),
  #         axis.text = element_text(size = 14, color = 'black', face = 'bold'),
  #         strip.text = element_text(size = 14, color = 'black', face = 'bold'),
  #         axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
  #         axis.ticks = element_line(color = 'black', linewidth = .7),
  #         axis.ticks.length = unit(.22, "cm"))
  # 
  # if (SAVE_PLOTS){
  #   ggsave(file.path(OUTDIR, "figures", "IBD_decay_by_recrate.png"), 
  #          p_decay_sub, width = 12, height = 8)
  # }
  
  # setDT(decay_dt)
  # 
  # # bin generation distances (e.g. 1-gen bins, or wider if needed)
  # decay_dt[, gen_bin := floor(gen_distance)]
  # 
  # # summarize: median + IQR per bin
  # sum_dt <- decay_dt[
  #   , .(
  #     ibd_med = median(ibd_prop, na.rm = TRUE),
  #     ibd_lo  = quantile(ibd_prop, 0.25, na.rm = TRUE),
  #     ibd_hi  = quantile(ibd_prop, 0.75, na.rm = TRUE),
  #     n = .N
  #   ),
  #   by = .(rate, gen_bin)]
  # 
  # p_decay_med <- ggplot(sum_dt, aes(x = gen_bin, y = ibd_med)) +
  #   geom_ribbon(aes(ymin = ibd_lo, ymax = ibd_hi), alpha = 0.2) +
  #   geom_line(color = "steelblue", linewidth = 1.3) +
  #   facet_wrap(~ rate, scales = "free_y") +
  #   labs(title = "IBD decay (median and IQR) by recombination rate",
  #        x = "Generations apart", y = "IBD proportion") +
  #   theme(legend.position = "bottom",
  #         legend.title = element_blank(),
  #         legend.text = element_text(size = 14, color = 'black', face = 'bold'),
  #         plot.title = element_text(color = 'black', face = 'bold'),
  #         axis.title = element_text(size = 16, color = 'black', face = 'bold'),
  #         axis.text = element_text(size = 14, color = 'black', face = 'bold'))
  # 
  # if (SAVE_PLOTS){
  #   ggsave(file.path(OUTDIR, "figures", "IBD_decay_by_recrate_v2.png"), 
  #          p_decay_med, width = 12, height = 10)
  # }


  # -------------------------------------------------------------------
  # Plot 8: AUPR across recombination rate
  # -------------------------------------------------------------------
  message("  [8/7] AUPR across recombination rate ...")

  pr_rows <- list()
  for (key in names(merged_scores_store)) {
    merged <- merged_scores_store[[key]]
    # parse key to rep and rate
    parts <- strsplit(key, "__")[[1]]
    repname <- parts[1]; rrate <- parts[2]
    method_cols <- grep("^score_", names(merged), value = TRUE)
    
    for (sc in method_cols) {
      method_name <- sub("^score_", "", sc)
      # need both positives and negatives
      if (sum(merged$true_transmission == 1, na.rm = TRUE) > 0 && sum(merged$true_transmission == 0, na.rm = TRUE) > 0) {
        pos <- merged[true_transmission == 1 & !is.na(get(sc)), get(sc)]
        neg <- merged[true_transmission == 0 & !is.na(get(sc)), get(sc)]
        if (length(pos) > 0 && length(neg) > 0) {
          # Try/catch to avoid PRROC crashing
          probj <- tryCatch(
            pr.curve(scores.class0 = pos, scores.class1 = neg, curve = TRUE),
            error = function(e) NULL
          )
          
          if (!is.null(probj)) {
            prmat <- probj$curve
            prdt <- as.data.table(prmat)
            setnames(prdt, c("recall","precision","threshold"))
            prdt[, replicate := repname]
            prdt[, rate := rate_lookup[rrate]]
            prdt[, method := method_name]
            pr_rows[[length(pr_rows)+1]] <- prdt
          }
        }
      }
    }
  }
  
  pr_all <- rbindlist(pr_rows, fill = TRUE)
  
  if (nrow(pr_all) > 0) {
    # interpolate to common recall grid per rate-method across replicates and compute mean±95%CI
    recall_grid <- seq(0,1,length.out = 200)
    pr_summary_rows <- list(); idx <- 0
    for (rr in unique(pr_all$rate)) {
      for (mm in unique(pr_all$method)) {
        sub <- pr_all[rate == rr & method == mm]
        if (nrow(sub) == 0) next
        # per replicate interpolation
        recs_by_rep <- split(sub, by = "replicate")
        interp_mat <- sapply(recs_by_rep, function(dt) {
          # remove duplicate recall values
          ok <- !duplicated(dt$recall)
          approx(x = dt$recall[ok], y = dt$precision[ok], xout = recall_grid, rule = 2)$y
        })
        if (is.null(interp_mat)) next
        mean_prec <- rowMeans(interp_mat, na.rm = TRUE)
        lo <- apply(interp_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
        hi <- apply(interp_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
        tmp <- data.table(rate = rr, 
                          method = mm, 
                          recall = recall_grid, 
                          precision_mean = mean_prec, 
                          precision_lo = lo, 
                          precision_hi = hi)
        
        pr_summary_rows[[length(pr_summary_rows)+1]] <- tmp
      }
    }
    pr_summary <- rbindlist(pr_summary_rows)
    
    p_pr_ribbon <- ggplot(pr_summary, aes(x = recall, y = precision_mean, color = method, fill = method)) +
      geom_ribbon(aes(ymin = precision_lo, ymax = precision_hi), alpha = 0.15, color = NA) +
      geom_line(linewidth = 1.5) +
      facet_wrap(~ rate, ncol = 2) +
      labs(title = paste0("Area Under Precision–Recall Curve across recombination rates (G=",GEN_THRESHOLD,")"), x = "Recall", y = "Precision") +
      theme(axis.text = element_text(size = 14, color = 'black', face = 'bold'), # angle = 45, hjust = 1,
            legend.text = element_text(size = 14, color = 'black', face = 'bold'),
            legend.title = element_blank(),
            plot.title = element_text(color = 'black', face = 'bold'),
            axis.title = element_text(size = 16, color = 'black', face = 'bold'),
            strip.text = element_text(size = 16, color = 'black', face = 'bold'))
    
    if (SAVE_PLOTS) {

      ggsave(file.path(OUTDIR, "figures", "PR_ribbons_by_rate.png"),
             p_pr_ribbon, width = 10, height = 8, dpi = 600)
      message("    ✓ Saved")
    }
  }

  # -------------------------------------------------------------------
  # Plot 8: Method ranking panel (AUPR and AUROC)
  # -------------------------------------------------------------------
  message("  [8/7] Method ranking panel (AUPR and AUROC) ...")
  
  rank_dt <- melt(agg_dt, 
                  id.vars = c("rate", "method"), 
                  measure.vars = c("aupr_mean","auroc_mean"), 
                  variable.name = "metric", 
                  value.name = "value")
  
  p_rank <- rank_dt %>% 
    filter(metric != "auroc_mean") %>% 
    ggplot( aes(x = rate, 
                y = value, 
                color = method, 
                group = method)) +
    geom_line(linewidth = 1.3) + 
    geom_point(size = 3) +
    # facet_wrap(~ metric, scales = "free_y") +
    scale_x_log10() +
    labs(x = "", y = "AUPR") +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14, color = 'black', face = 'bold'),
          plot.title = element_text(color = 'black', face = 'bold'),
          axis.title = element_text(size = 16, color = 'black', face = 'bold'),
          axis.text = element_text(size = 13, color = 'black', face = 'bold'))
  
  if (SAVE_PLOTS) {
    ggsave(file.path(OUTDIR, "figures", "Method_ranking_AUPR.png"), 
           p_rank, width = 12, height = 10, dpi = 600)
  }
  
  # -------------------------------------------------------------------
  # Save workspace
  # -------------------------------------------------------------------
  message("\n  Saving workspace...")
  saveRDS(
    list(
      metrics_dt = metrics_dt,
      agg_dt = agg_dt,
      detectability_dt = detectability_dt,
      merged_scores_store = merged_scores_store,
      # global_thresholds = global_thresholds,
      ibd_eps_global = ibd_eps_global,
      config = list(
        REC_RATES = REC_RATES,
        GEN_THRESHOLD = GEN_THRESHOLD,
        ALPHA_DETECT = ALPHA_DETECT)
    ),
    file.path(OUTDIR, "tables", "recombination_workspace.rds")
  )
  
  message("    ✓ Workspace saved")
  
  # ============================================================================
  # FINAL SUMMARY
  # ============================================================================
  
  message("\n", rep("=", 40))
  message("\t ANALYSIS COMPLETE")
  message(rep("=", 40))
  message("\n[SUMMARY STATISTICS]")
  message("  • Total evaluations: ", nrow(metrics_dt))
  message("  • Baseline detectability: ", 
          round(100 * detect_agg[rate == baseline_rate, detectability_mean], 1), "%")
  message("  • Highest rate detectability: ",
          round(100 * detect_agg[rate == max(rate), detectability_mean], 1), "%")
  message("\n[BEST PERFORMING METHOD (by AUPR at baseline)]")
  best_method <- agg_dt[rate == baseline_rate][which.max(aupr_mean)]
  message("  • Method: ", best_method$method)
  message("  • AUPR: ", round(best_method$aupr_mean, 3), " ± ", round(best_method$aupr_se, 3))
  message("\n[OUTPUT LOCATION]")
  message("  ", OUTDIR)
  message("\n", rep("=", 60), "\n")
  
  return(list(
    metrics = metrics_dt,
    agg = agg_dt,
    detectability = detectability_dt,
    merged = merged_scores_store,
    # thresholds = global_thresholds,
    ibd_eps = ibd_eps_global
  ))
}

# If you source this script, you can run:
res <- run_recombination_analysis()
