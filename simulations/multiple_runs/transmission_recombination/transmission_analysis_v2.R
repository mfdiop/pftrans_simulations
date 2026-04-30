# ============================================================================
# IBD Method Evaluation: Transmission Discrimination Under Recombination
# ============================================================================
# 
# RESEARCH QUESTION:
# Does high recombination improve IBD's ability to discriminate recent 
# transmission from distant ancestry, even as overall detection sensitivity 
# decreases?
#
# KEY HYPOTHESES:
# H1: IBD detection sensitivity decreases with recombination (shorter segments)
# H2: IBD transmission discrimination improves with recombination (better separation)
# H3: Among detectable pairs, recent transmission is enriched at high recombination
#
# GROUND TRUTH DEFINITIONS:
# - Transmission: gen_distance ≤ 5 (recent, direct transmission chains)
# - Recent Ancestry: 6 ≤ gen_distance ≤ 15 (detectable but not transmission)
# - Distant Ancestry: 16 ≤ gen_distance ≤ 25 (may be undetectable)
# - Unrelated: gen_distance > 25 (background)
#
# METRICS EVALUATED:
# - IBD (HMM-based): Segment length dependent
# - IBS (genetic similarity): All SNPs, segment independent
# - Phylogenetic (patristic distance): All SNPs, evolutionary distance
#
# Author: [Your Name]
# Date: 2025-02-08
# ============================================================================

library(data.table)
library(tidyverse)
library(PRROC)
library(pROC)
library(ape)        # For phylogenetic tree handling
library(phangorn)   # For patristic distances
library(gridExtra)
library(patchwork)
library(viridis)

set.seed(42)

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIG <- list(
  # Paths
  ROOT_DIR = "simulations/multiple_runs",
  INFERRED_SUBPATH = "inferred",
  PHYLO_SUBPATH = "phylo_results",
  OUTDIR = "simulations/multiple_runs/transmission_evaluation",
  
  # Analysis parameters
  REC_RATES = c("1e-09", "1e-08", "1e-07", "1e-06"),
  GENOME_BP = 640851,
  CHR = "chr1",
  
  # Ground truth thresholds (generations)
  GEN_TRANSMISSION = 5,      # Recent transmission
  GEN_RECENT_ANCESTRY = 15,  # Recent ancestry
  GEN_DISTANT_ANCESTRY = 25, # Distant ancestry
  
  # Detectability
  ALPHA_DETECT = 0.01,  # False positive rate for detectability threshold
  
  # Evaluation
  TOP_K = c(5, 10, 25, 50, 100),
  
  # Output
  SAVE_PLOTS = TRUE,
  SAVE_TABLES = TRUE
)

# Create output directories
dir.create(CONFIG$OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CONFIG$OUTDIR, "tables"), showWarnings = FALSE)
dir.create(file.path(CONFIG$OUTDIR, "figures"), showWarnings = FALSE)

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Create canonical pair key (ensures id1 <= id2)
canonical_pair <- function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_character_)
  if (a <= b) paste(a, b, sep = "--") else paste(b, a, sep = "--")
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

# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

#' Load ground truth IBD data
#' @return data.table with pair_key, id1, id2, total_ibd_bp, ibd_prop, gen_distance
load_ground_truth <- function(rep_dir, rate_str, genome_bp = CONFIG$GENOME_BP) {
  
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
  
  # Standardize sample IDs
  dt <- dt %>% 
    mutate(Id1 = paste0("tsk_", Id1), Id2 = paste0("tsk_", Id2)) %>% 
    as.data.table()
  
  # Standardize column names
  id_cols <- intersect(c("id1", "id2", "Id1", "Id2"), names(dt))
  if (length(id_cols) >= 2) {
    setnames(dt, old = id_cols[1:2], new = c("id1", "id2"), skip_absent = TRUE)
  }
  
  # Create pair key
  if (!("pair_key" %in% names(dt))) {
    dt[, pair_key := mapply(canonical_pair, id1, id2)]
  }
  
  # Ensure IBD proportion
  if (!("total_ibd_bp" %in% names(dt))) dt[, total_ibd_bp := NA_real_]
  if (!("true_ibd_prop" %in% names(dt))) {
    dt[, true_ibd_prop := total_ibd_bp / genome_bp]
  }
  if (!("ibd_prop" %in% names(dt))) dt[, ibd_prop := true_ibd_prop]
  
  # Ensure generation distance
  if (!("gen_distance" %in% names(dt))) {
    possible_gens <- intersect(c("generations", "t_generations", "min_tmrca"), names(dt))
    if (length(possible_gens) >= 1) {
      setnames(dt, old = possible_gens[1], new = "gen_distance", skip_absent = TRUE)
    }
  }
  
  # Get max segment
  if ("max_segment_bp" %in% names(dt)) {
    dt[, max_seg_bp := max_segment_bp]
  } else {
    dt[, max_seg_bp := NA_real_]
  }
  
  return(dt[, .(pair_key, id1, id2, total_ibd_bp, ibd_prop, max_seg_bp, gen_distance)])
}

#' Load inferred IBD (HMM-IBD predictions)
#' @return data.table with pair_key, id1, id2, score (similarity)
load_ibd_hmm <- function(rep_dir, rate_str) {
  
  repnum <- gsub("rep", "", basename(rep_dir))
  candidate <- file.path(rep_dir, paste0("run", repnum, "_rec", rate_str, "_chr1"))
  
  if (!dir.exists(candidate)) {
    cand <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
    cand <- cand[grepl(paste0("rec", rate_str), cand)]
    if (length(cand) >= 1) candidate <- cand[1] else candidate <- rep_dir
  }
  
  fp <- file.path(candidate, "inferred_ibd_hmm.rds")
  dt <- as.data.table(safe_readRDS(fp))
  
  if (is.null(dt)) return(NULL)
  
  # Standardize column names
  id_cols <- intersect(c("id1", "id2", "Id1", "Id2"), names(dt))
  if (length(id_cols) >= 2) {
    setnames(dt, old = id_cols[1:2], new = c("id1", "id2"), skip_absent = TRUE)
  }
  
  # Extract score (total IBD or other metric - higher = more related)
  score_candidates <- c("total_ibd_bp", "total_ibd", "ibd_prop", "score", "posterior")
  score_col <- intersect(score_candidates, names(dt))
  
  if (length(score_col) > 0) {
    dt[, score := as.numeric(get(score_col[1]))]
  } else if ("n_segments" %in% names(dt)) {
    dt[, score := as.numeric(n_segments)]
  } else {
    dt[, score := 1]
  }
  
  # Create pair key
  if (!("pair_key" %in% names(dt))) {
    dt[, pair_key := mapply(canonical_pair, id1, id2)]
  }
  
  return(unique(dt[, .(pair_key, id1, id2, score)]))
}

#' Load IBS matrix
#' @return data.table with pair_key, id1, id2, score (similarity)
load_ibs <- function(rep_dir, rate_str) {
  
  repnum <- gsub("rep", "", basename(rep_dir))
  candidate <- file.path(rep_dir, paste0("run", repnum, "_rec", rate_str, "_chr1"))
  
  if (!dir.exists(candidate)) {
    cand <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
    cand <- cand[grepl(paste0("rec", rate_str), cand)]
    if (length(cand) >= 1) candidate <- cand[1] else candidate <- rep_dir
  }
  
  fp <- file.path(candidate, "ibs_matrix.rds")
  if (!file.exists(fp)) return(NULL)
  
  obj <- safe_readRDS(fp)
  if (is.null(obj)) return(NULL)
  
  mat <- as.matrix(obj)
  ids <- rownames(mat)
  
  if (is.null(ids)) {
    if (is.data.frame(obj)) {
      ids <- as.character(obj[[1]])
      mat <- as.matrix(obj[, -1])
      rownames(mat) <- ids
      colnames(mat) <- ids
    } else {
      return(NULL)
    }
  }
  
  # Convert to long format
  n <- length(ids)
  rows <- vector("list", n * (n - 1) / 2)
  idx <- 1
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      rows[[idx]] <- list(
        pair_key = canonical_pair(ids[i], ids[j]),
        id1 = ids[i],
        id2 = ids[j],
        score = as.numeric(mat[i, j])
      )
      idx <- idx + 1
    }
  }
  
  dt <- rbindlist(rows)
  
  # Assume IBS is similarity (proportion shared alleles)
  # If your IBS is a distance, you need to invert it here:
  # dt[, score := max(score) - score] or dt[, score := 1 - score]
  
  return(unique(dt[, .(pair_key, id1, id2, score)]))
}

#' Load phylogenetic tree and compute patristic distances
#' Converts distances to similarity scores
#' @return data.table with pair_key, id1, id2, score (similarity)
load_phylo <- function(phylo_dir, rate_str, chr) {
  
  # Find tree file matching recombination rate
  patt <- paste0("rec", rate_str, "_", chr)
  files <- list.files(phylo_dir, full.names = TRUE, pattern = patt)
  
  # Prefer modelfinder treefile
  tf <- files[grepl(".treefile$", files)]
  
  if (length(tf) == 0 && length(files) > 0) tf <- files[1]
  if (length(tf) == 0) return(NULL)
  
  message("      Reading tree: ", basename(tf[1]))
  
  # Read tree
  tree <- tryCatch(
    read.tree(tf[1]),
    error = function(e) {
      warning("Failed to read tree: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(tree)) return(NULL)
  
  # Compute patristic distances
  message("      Computing patristic distances...")
  dist_mat <- tryCatch(
    cophenetic.phylo(tree),  # This gives patristic distances
    error = function(e) {
      warning("Failed to compute patristic distances: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(dist_mat)) return(NULL)
  
  # Convert distance matrix to long format
  ids <- rownames(dist_mat)
  n <- length(ids)
  rows <- vector("list", n * (n - 1) / 2)
  idx <- 1
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      rows[[idx]] <- list(
        pair_key = canonical_pair(ids[i], ids[j]),
        id1 = ids[i],
        id2 = ids[j],
        distance = as.numeric(dist_mat[i, j])
      )
      idx <- idx + 1
    }
  }
  
  dt <- rbindlist(rows)
  
  # ✅ CRITICAL: Convert distance to similarity
  # Patristic distance: smaller = more related
  # We need: higher score = more related (for AUPR/AUROC)
  max_dist <- max(dt$distance, na.rm = TRUE)
  # dt[, score := max_dist - distance]  # Invert: max distance becomes 0, min becomes max
  
  # Better alternatives (ordered from minimal to ideal)
  # ✅ Option 1: Normalize by tree diameter (minimal fix)
  dt[, score := 1 - (distance / max_dist)]
  
  # # ✅ Option 2: Exponential decay kernel (recommended)
  # # Convert patristic distance to similarity using exponential kernel
  # alpha <- 1 / median(dt$distance, na.rm = TRUE)
  # dt[, score := exp(-alpha * distance)]
  
  message("      Patristic distance range: [", 
          round(min(dt$distance, na.rm = TRUE), 4), ", ",
          round(max(dt$distance, na.rm = TRUE), 4), "]")
  message("      Converted similarity range: [",
          round(min(dt$score, na.rm = TRUE), 4), ", ",
          round(max(dt$score, na.rm = TRUE), 4), "]")
  
  return(dt[, .(pair_key, id1, id2, score)])
}

# ============================================================================
#                       GROUND TRUTH ANNOTATION
# ============================================================================

#' Annotate pairs with relationship categories and ground truths
annotate_relationships <- function(dt, 
                                   gen_transmission = CONFIG$GEN_TRANSMISSION,
                                   gen_recent = CONFIG$GEN_RECENT_ANCESTRY,
                                   gen_distant = CONFIG$GEN_DISTANT_ANCESTRY) {
  
  # Relationship categories
  dt[, relationship := fcase(
    gen_distance <= gen_transmission, "transmission",
    gen_distance <= gen_recent, "recent_ancestry",
    gen_distance <= gen_distant, "distant_ancestry",
    default = "unrelated"
  )]
  
  dt[, relationship := factor(relationship, 
                              levels = c("transmission", "recent_ancestry", 
                                       "distant_ancestry", "unrelated"))]
  
  # Binary ground truths for different tasks
  dt[, truth_transmission := as.integer(gen_distance <= gen_transmission)]
  dt[, truth_recent_ancestry := as.integer(gen_distance <= gen_recent)]
  dt[, truth_any_relatedness := as.integer(gen_distance <= gen_distant)]
  
  # For transmission discrimination (among related pairs only)
  dt[, is_transmission_vs_ancestry := fcase(
    gen_distance <= gen_transmission, 1L,
    gen_distance > gen_transmission & gen_distance <= gen_distant, 0L,
    default = NA_integer_
  )]
  
  message("  Relationship breakdown:")
  message("    Transmission (≤", gen_transmission, " gen): ", 
          sum(dt$relationship == "transmission"))
  message("    Recent ancestry (", gen_transmission + 1, "-", gen_recent, " gen): ",
          sum(dt$relationship == "recent_ancestry"))
  message("    Distant ancestry (", gen_recent + 1, "-", gen_distant, " gen): ",
          sum(dt$relationship == "distant_ancestry"))
  message("    Unrelated (>", gen_distant, " gen): ",
          sum(dt$relationship == "unrelated"))
  
  return(dt)
}

#' Define detectability threshold based on false positive rate
define_detectability_threshold <- function(dt, truth_col = "truth_any_relatedness",
                                          alpha = CONFIG$ALPHA_DETECT) {
  
  negatives <- dt[get(truth_col) == 0 & !is.na(ibd_prop), ibd_prop]
  positives <- dt[get(truth_col) == 1 & !is.na(ibd_prop), ibd_prop]
  
  if (length(negatives) == 0) {
    warning("  No negatives found, using 0 as threshold")
    return(0)
  }
  
  # Use (1 - alpha) quantile of negatives
  threshold <- quantile(negatives, probs = 1 - alpha, na.rm = TRUE)
  
  message("  Detectability threshold (", (1-alpha)*100, "th percentile of negatives):")
  message("    Threshold: ", round(threshold, 6))
  message("    Negatives > threshold: ", 
          round(100 * mean(negatives > threshold, na.rm = TRUE), 2), "%")
  
  if (length(positives) > 0) {
    message("    Positives > threshold: ",
            round(100 * mean(positives > threshold, na.rm = TRUE), 2), "%")
  }
  
  return(as.numeric(threshold))
}

# ============================================================================
#                         EVALUATION FUNCTIONS
# ============================================================================

#' Compute evaluation metrics for binary classification
#' All scores assumed to be similarities (higher = more related)
evaluate_binary_classification <- function(dt, score_col, truth_col, 
                                          ks = CONFIG$TOP_K) {
  
  res <- list(score_col = score_col, truth_col = truth_col)
  
  # Remove NAs
  valid <- !is.na(dt[[score_col]]) & !is.na(dt[[truth_col]])
  if (sum(valid) == 0) {
    res$aupr <- res$auroc <- NA_real_
    res$TP <- res$TN <- res$FP <- res$FN <- NA_integer_
    res$precision <- res$recall <- res$f1 <- NA_real_
    res$confusion_matrix <- NULL
    return(res)
  }
  
  
  dt_eval <- dt[valid]
  scores <- dt_eval[[score_col]]
  truth <- dt_eval[[truth_col]]
  
  # Check for class imbalance
  n_pos <- sum(truth == 1)
  n_neg <- sum(truth == 0)
  
  res$n_pos <- n_pos
  res$n_neg <- n_neg
  res$prevalence <- n_pos / (n_pos + n_neg)
  
  if (n_pos == 0 || n_neg == 0) {
    res$aupr <- res$auroc <- NA_real_
    res$TP <- res$TN <- res$FP <- res$FN <- NA_integer_
    res$precision <- res$recall <- res$f1 <- NA_real_
    res$confusion_matrix <- NULL
    return(res)
  }
  
  # AUPR
  pos_scores <- scores[truth == 1]
  neg_scores <- scores[truth == 0]
  
  pr_obj <- tryCatch(
    PRROC::pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = TRUE),
    error = function(e) {
      warning("PR curve failed: ", e$message)
      NULL
    }
  )
  
  if (!is.null(pr_obj)) {
    res$aupr <- pr_obj$auc.integral
    res$pr_curve <- as.data.table(pr_obj$curve)
    setnames(res$pr_curve, c("recall", "precision", "threshold"))
  } else {
    res$aupr <- NA_real_
  }
  
  # AUROC
  roc_obj <- tryCatch(
    pROC::roc(truth, scores, quiet = TRUE),
    error = function(e) {
      warning("ROC failed: ", e$message)
      NULL
    }
  )
  
  if (!is.null(roc_obj)) {
    res$auroc <- as.numeric(pROC::auc(roc_obj))
    
    # Get optimal threshold (Youden's index)
    coords <- pROC::coords(roc_obj, "best", best.method = "youden")
    res$optimal_threshold <- coords$threshold
    res$optimal_sensitivity <- coords$sensitivity
    res$optimal_specificity <- coords$specificity
  } else {
    res$auroc <- NA_real_
  }
  
  # Spearman correlation with continuous IBD
  if ("ibd_prop" %in% names(dt_eval)) {
    res$spearman <- cor(scores, dt_eval$ibd_prop, method = "spearman", 
                       use = "complete.obs")
  }
  
  # Precision@K
  dt_ranked <- dt_eval[order(-get(score_col))]  # Descending order
  
  for (k in ks) {
    if (nrow(dt_ranked) >= k) {
      topk <- head(dt_ranked, k)
      prec <- mean(topk[[truth_col]] == 1)
      res[[paste0("prec_at_", k)]] <- prec
    } else {
      res[[paste0("prec_at_", k)]] <- NA_real_
    }
  }
  
  # ---------------------------
  # Confusion matrix at top-K
  # ---------------------------
  
  for (k in ks) {
    if (nrow(dt_ranked) >= k) {
      
      dt_ranked[, pred := 0L]
      dt_ranked[1:k, pred := 1L]
      
      TP <- sum(dt_ranked$pred == 1 & dt_ranked[[truth_col]] == 1)
      FP <- sum(dt_ranked$pred == 1 & dt_ranked[[truth_col]] == 0)
      FN <- sum(dt_ranked$pred == 0 & dt_ranked[[truth_col]] == 1)
      TN <- sum(dt_ranked$pred == 0 & dt_ranked[[truth_col]] == 0)
      
      precision <- ifelse(TP + FP > 0, TP / (TP + FP), NA_real_)
      recall    <- ifelse(TP + FN > 0, TP / (TP + FN), NA_real_)
      f1        <- ifelse(!is.na(precision + recall) && (precision + recall) > 0,
                          2 * precision * recall / (precision + recall),
                          NA_real_)
      
      res[[paste0("TP_at_", k)]] <- TP
      res[[paste0("FP_at_", k)]] <- FP
      res[[paste0("FN_at_", k)]] <- FN
      res[[paste0("TN_at_", k)]] <- TN
      res[[paste0("precision_at_", k)]] <- precision
      res[[paste0("recall_at_", k)]] <- recall
      res[[paste0("f1_at_", k)]] <- f1
      
    } else {
      res[[paste0("TP_at_", k)]] <- NA_integer_
      res[[paste0("FP_at_", k)]] <- NA_integer_
      res[[paste0("FN_at_", k)]] <- NA_integer_
      res[[paste0("TN_at_", k)]] <- NA_integer_
      res[[paste0("precision_at_", k)]] <- NA_real_
      res[[paste0("recall_at_", k)]] <- NA_real_
      res[[paste0("f1_at_", k)]] <- NA_real_
    }
  }
  
  return(res)
}

#' Evaluate enrichment of transmission pairs in top K
evaluate_transmission_enrichment <- function(dt, score_col, 
                                            ks = CONFIG$TOP_K,
                                            gen_transmission = CONFIG$GEN_TRANSMISSION) {
  
  valid <- !is.na(dt[[score_col]])
  dt_eval <- dt[valid]
  
  if (nrow(dt_eval) == 0) return(NULL)
  
  # Rank by score
  dt_ranked <- dt_eval[order(-get(score_col))]
  
  results <- list()
  
  for (k in ks) {
    if (nrow(dt_ranked) >= k) {
      topk <- head(dt_ranked, k)
      
      # Count relationship types in top K
      n_transmission <- sum(topk$gen_distance <= gen_transmission, na.rm = TRUE)
      n_ancestry <- sum(topk$gen_distance > gen_transmission & 
                       topk$gen_distance <= CONFIG$GEN_DISTANT_ANCESTRY, na.rm = TRUE)
      n_unrelated <- sum(topk$gen_distance > CONFIG$GEN_DISTANT_ANCESTRY, na.rm = TRUE)
      
      # Baseline (expected if random sampling)
      baseline_transmission <- mean(dt_eval$gen_distance <= gen_transmission, na.rm = TRUE)
      
      # Enrichment
      observed_transmission <- n_transmission / k
      enrichment <- observed_transmission / baseline_transmission
      
      results[[paste0("k", k)]] <- data.table(
        k = k,
        n_transmission = n_transmission,
        n_ancestry = n_ancestry,
        n_unrelated = n_unrelated,
        prop_transmission = observed_transmission,
        baseline_transmission = baseline_transmission,
        enrichment = enrichment
      )
    }
  }
  
  return(rbindlist(results))
}

# ============================================================================
#                         MAIN ANALYSIS PIPELINE
# ============================================================================

run_transmission_analysis <- function() {
  
  message("\n", rep("=", 70))
  message("  TRANSMISSION DISCRIMINATION UNDER RECOMBINATION")
  message(rep("=", 70), "\n")
  
  message("[CONFIG]")
  message("  Root directory: ", CONFIG$ROOT_DIR)
  message("  Recombination rates: ", paste(CONFIG$REC_RATES, collapse = ", "))
  message("  Transmission threshold: ≤", CONFIG$GEN_TRANSMISSION, " generations")
  message("  Recent ancestry: ≤", CONFIG$GEN_RECENT_ANCESTRY, " generations")
  message("  Distant ancestry: ≤", CONFIG$GEN_DISTANT_ANCESTRY, " generations")
  message("  Output: ", CONFIG$OUTDIR)
  
  # Detect replicates
  inferred_root <- file.path(CONFIG$ROOT_DIR, CONFIG$INFERRED_SUBPATH)
  replicates <- list.dirs(inferred_root, recursive = FALSE, full.names = FALSE)
  replicates <- replicates[grepl("^rep", replicates)]
  
  message("\n  Detected replicates: ", paste(replicates, collapse = ", "))
  
  # Create rate lookup
  rec_rates_label <- gsub("-", "", CONFIG$REC_RATES)
  rate_lookup <- setNames(as.numeric(CONFIG$REC_RATES), rec_rates_label)
  
  # Storage
  all_results <- list()
  all_enrichment <- list()
  all_detectability <- list()
  merged_data_store <- list()
  all_eval_overall <- list()
  all_eval_transmission <- list()
  all_eval_discrimination <- list()
  
  baseline_rate <- min(as.numeric(CONFIG$REC_RATES))
  ibd_eps_global <- NULL
  chr <- CONFIG$CHR
  
  # ============================================================================
  #               STEP 1: DATA LOADING AND EVALUATION
  # ============================================================================
  
  message("\n", rep("=", 70))
  message("  STEP 1: LOADING DATA AND COMPUTING METRICS")
  message(rep("=", 70))
  
  for (rep_idx in seq_along(replicates)) {
    rep <- replicates[rep_idx]
    rep_dir <- file.path(inferred_root, rep)
    phylo_dir <- file.path(CONFIG$ROOT_DIR, CONFIG$PHYLO_SUBPATH, rep)
    
    message("\n[", rep_idx, "/", length(replicates), "] Replicate: ", rep)
    
    if (!dir.exists(rep_dir)) {
      warning("  Directory not found: ", rep_dir)
      next
    }
    
    for (rate_idx in seq_along(rec_rates_label)) {
      rate_label <- rec_rates_label[rate_idx]
      rate_numeric <- rate_lookup[rate_label]
      is_baseline <- (rate_numeric == baseline_rate)
      
      message("\n  Rate ", rate_idx, "/", length(rec_rates_label), ": ",
              rate_label, " (", rate_numeric, " bp⁻¹)",
              ifelse(is_baseline, " [BASELINE]", ""))
      
      # ----------------------------------------------------------------------
      # Load data
      # ----------------------------------------------------------------------
      message("  Loading data...")
      
      # Ground truth
      truth <- load_ground_truth(rep_dir, rate_label)
      if (is.null(truth) || nrow(truth) == 0) {
        warning("    No ground truth data")
        next
      }
      message("    ✓ Ground truth: ", nrow(truth), " pairs")
      
      # Annotate relationships
      truth <- annotate_relationships(truth)
      
      # Define detectability threshold (only at baseline)
      if (is_baseline) {
        ibd_eps_global <- define_detectability_threshold(truth)
      }
      
      # Mark detectable pairs
      truth[, detectable := as.integer(ibd_prop >= ibd_eps_global)]
      
      # IBD method
      ibd <- load_ibd_hmm(rep_dir, rate_label)
      if (!is.null(ibd)) {
        message("    ✓ IBD: ", nrow(ibd), " pairs")
      }
      
      # IBS
      ibs <- load_ibs(rep_dir, rate_label)
      if (!is.null(ibs)) {
        message("    ✓ IBS: ", nrow(ibs), " pairs")
      }
      
      # Phylogenetic
      phylo <- NULL
      if (dir.exists(phylo_dir)) {
        phylo <- load_phylo(phylo_dir, rate_label, chr)
        if (!is.null(phylo)) {
          message("    ✓ Phylo: ", nrow(phylo), " pairs")
        }
      }
      
      # ----------------------------------------------------------------------
      # Merge data
      # ----------------------------------------------------------------------
      message("  Merging methods...")
      
      merged <- truth[, .(pair_key, id1, id2, ibd_prop, gen_distance, 
                         relationship, truth_transmission, truth_recent_ancestry,
                         truth_any_relatedness, is_transmission_vs_ancestry,
                         detectable, max_seg_bp)]
      
      if (!is.null(ibd)) {
        ibd_scores <- ibd[, .(pair_key, score_IBD = score)]
        merged <- merge(merged, ibd_scores, by = "pair_key", all.x = TRUE)
      }
      
      if (!is.null(ibs)) {
        ibs_scores <- ibs[, .(pair_key, score_IBS = score)]
        merged <- merge(merged, ibs_scores, by = "pair_key", all.x = TRUE)
      }
      
      if (!is.null(phylo)) {
        phylo_scores <- phylo[, .(pair_key, score_Phylo = score)]
        merged <- merge(merged, phylo_scores, by = "pair_key", all.x = TRUE)
      }
      
      message("    ✓ Merged: ", nrow(merged), " pairs × ", ncol(merged), " columns")
      
      # Store merged data
      merged_data_store[[paste(rep, rate_label, sep = "__")]] <- merged
      
      # ----------------------------------------------------------------------
      # Detectability analysis by relationship type
      # ----------------------------------------------------------------------
      message("  Analyzing detectability...")
      
      detect_stats <- merged[, .(
        n_pairs = .N,
        n_detectable = sum(detectable == 1, na.rm = TRUE),
        detectability = mean(detectable, na.rm = TRUE),
        mean_ibd = mean(ibd_prop, na.rm = TRUE),
        median_ibd = median(ibd_prop, na.rm = TRUE)
      ), by = relationship]
      
      detect_stats[, `:=`(
        replicate = rep,
        rate = rate_numeric,
        rate_label = rate_label
      )]
      
      all_detectability[[length(all_detectability) + 1]] <- detect_stats
      
      message("  Detectability by relationship:")
      for (i in 1:nrow(detect_stats)) {
        message("    ", detect_stats$relationship[i], ": ",
                round(100 * detect_stats$detectability[i], 1), "% (",
                detect_stats$n_detectable[i], "/", detect_stats$n_pairs[i], ")")
      }
      
      # ----------------------------------------------------------------------
      # Evaluate methods on multiple tasks
      # ----------------------------------------------------------------------
      message("  Evaluating methods...")
      
      methods <- grep("^score_", names(merged), value = TRUE)
      methods <- gsub("^score_", "", methods)
      
      for (method in methods) {
        score_col <- paste0("score_", method)
        
        message("\n    Method: ", method)
        
        # Task 1: Overall relatedness detection (≤25 gen vs >25 gen)
        eval_overall <- evaluate_binary_classification(
          merged, score_col, "truth_any_relatedness"
        )
        
        # Task 2: Transmission detection (≤5 gen vs >5 gen)
        eval_transmission <- evaluate_binary_classification(
          merged, score_col, "truth_transmission"
        )
        
        # Task 3: Transmission discrimination (≤5 vs 6-25, among related only)
        merged_related <- merged[!is.na(is_transmission_vs_ancestry)]
        eval_discrimination <- if (nrow(merged_related) > 0) {
          evaluate_binary_classification(
            merged_related, score_col, "is_transmission_vs_ancestry"
          )
        } else {
          list(aupr = NA_real_, auroc = NA_real_)
        }
        
        # Task 4: Transmission enrichment in top K
        enrichment <- evaluate_transmission_enrichment(merged, score_col)
        if (!is.null(enrichment)) {
          enrichment[, `:=`(
            replicate = rep,
            rate = rate_numeric,
            rate_label = rate_label,
            method = method
          )]
          all_enrichment[[length(all_enrichment) + 1]] <- enrichment
        }
        
        # Store results
        result_row <- data.table(
          replicate = rep,
          rate = rate_numeric,
          rate_label = rate_label,
          method = method,
          
          # Overall relatedness
          aupr_overall = eval_overall$aupr,
          auroc_overall = eval_overall$auroc,
          n_pos_overall = eval_overall$n_pos,
          n_neg_overall = eval_overall$n_neg,
          
          # Transmission detection
          aupr_transmission = eval_transmission$aupr,
          auroc_transmission = eval_transmission$auroc,
          n_pos_transmission = eval_transmission$n_pos,
          n_neg_transmission = eval_transmission$n_neg,
          
          # Transmission discrimination (among related)
          aupr_discrimination = eval_discrimination$aupr,
          auroc_discrimination = eval_discrimination$auroc,
          
          # Spearman correlation with IBD
          spearman = eval_overall$spearman,
          
          # Precision@K for transmission
          prec_at_10_trans = eval_transmission$prec_at_10,
          prec_at_25_trans = eval_transmission$prec_at_25,
          prec_at_50_trans = eval_transmission$prec_at_50
        )
        
        all_results[[length(all_results) + 1]] <- result_row
        all_eval_overall[[length(all_eval_overall) + 1]] <- eval_overall
        all_eval_transmission[[length(all_eval_transmission) + 1]] <- eval_transmission
        all_eval_discrimination[[length(all_eval_discrimination) + 1]] <- eval_discrimination
        
        eval_overall$replicate  <- rep
        eval_overall$rate       <- as.numeric(rate_numeric)
        eval_overall$rate_label <-  rate_label
        eval_overall$method     <- method
        
        
        message("      Overall relatedness AUPR: ", round(eval_overall$aupr, 3))
        message("      Transmission detection AUPR: ", round(eval_transmission$aupr, 3))
        message("      Transmission discrimination AUPR: ", round(eval_discrimination$aupr, 3))
      }
    }
  }
  
  # ============================================================================
  #                 STEP 2: AGGREGATE AND SUMMARIZE
  # ============================================================================
  
  message("\n", rep("=", 70))
  message("  STEP 2: AGGREGATING RESULTS")
  message(rep("=", 70))
  
  results_dt <- rbindlist(all_results, fill = TRUE)
  detectability_dt <- rbindlist(all_detectability, fill = TRUE)
  enrichment_dt <- rbindlist(all_enrichment, fill = TRUE)
  
  # Compute summary statistics
  agg_results <- results_dt[, .(
    # Overall relatedness
    aupr_overall_mean = mean(aupr_overall, na.rm = TRUE),
    aupr_overall_se = sd(aupr_overall, na.rm = TRUE) / sqrt(.N),
    auroc_overall_mean = mean(auroc_overall, na.rm = TRUE),
    
    # Transmission detection
    aupr_transmission_mean = mean(aupr_transmission, na.rm = TRUE),
    aupr_transmission_se = sd(aupr_transmission, na.rm = TRUE) / sqrt(.N),
    auroc_transmission_mean = mean(auroc_transmission, na.rm = TRUE),
    
    # Transmission discrimination
    aupr_discrimination_mean = mean(aupr_discrimination, na.rm = TRUE),
    aupr_discrimination_se = sd(aupr_discrimination, na.rm = TRUE) / sqrt(.N),
    auroc_discrimination_mean = mean(auroc_discrimination, na.rm = TRUE),
    
    # Correlation
    spearman_mean = mean(spearman, na.rm = TRUE),
    
    # Sample size
    n_replicates = .N
  ), by = .(rate, method)]
  
  # Detectability summary
  detect_summary <- detectability_dt[, .(
    detectability_mean = mean(detectability, na.rm = TRUE),
    detectability_se = sd(detectability, na.rm = TRUE) / sqrt(.N),
    mean_ibd_mean = mean(mean_ibd, na.rm = TRUE),
    n_pairs_mean = mean(n_pairs, na.rm = TRUE)
  ), by = .(rate, relationship)]
  
  # Enrichment summary
  if (nrow(enrichment_dt) > 0) {
    enrichment_summary <- enrichment_dt[, .(
      enrichment_mean = mean(enrichment, na.rm = TRUE),
      enrichment_se = sd(enrichment, na.rm = TRUE) / sqrt(.N),
      prop_transmission_mean = mean(prop_transmission, na.rm = TRUE),
      baseline_transmission_mean = mean(baseline_transmission, na.rm = TRUE)
    ), by = .(rate, method, k)]
  } else {
    enrichment_summary <- data.table()
  }
  
  # ============================================================================
  #                           STEP 3: SAVE RESULTS
  # ============================================================================
  
  message("\n", rep("=", 70))
  message("  STEP 3: SAVING RESULTS")
  message(rep("=", 70))
  
  if (CONFIG$SAVE_TABLES) {
    message("\n  Saving tables...")
    
    fwrite(results_dt, 
           file.path(CONFIG$OUTDIR, "tables", "results_all_replicates.csv"))
    fwrite(agg_results,
           file.path(CONFIG$OUTDIR, "tables", "results_aggregated.csv"))
    fwrite(detectability_dt,
           file.path(CONFIG$OUTDIR, "tables", "detectability_by_relationship.csv"))
    fwrite(detect_summary,
           file.path(CONFIG$OUTDIR, "tables", "detectability_summary.csv"))
    
    readRDS(all_eval_overall, file.path(CONFIG$OUTDIR, "tables", "detectability_summary.csv"))
    
    if (nrow(enrichment_dt) > 0) {
      fwrite(enrichment_dt,
             file.path(CONFIG$OUTDIR, "tables", "enrichment_all_replicates.csv"))
      fwrite(enrichment_summary,
             file.path(CONFIG$OUTDIR, "tables", "enrichment_summary.csv"))
    }
    
    message("    ✓ Tables saved")
  }
  
  # ============================================================================
  #                       STEP 4: VISUALIZATIONS
  # ============================================================================
  
  message("\n", rep("=", 70))
  message("  STEP 4: CREATING VISUALIZATIONS")
  message(rep("=", 70))
  
  if (CONFIG$SAVE_PLOTS) {
    source_code_for_plots(agg_results, detect_summary, enrichment_summary)
  }
  
  # ============================================================================
  # FINAL SUMMARY
  # ============================================================================
  
  message("\n", rep("=", 70))
  message("  ANALYSIS COMPLETE")
  message(rep("=", 70))
  
  message("\n[SUMMARY]")
  message("  Total evaluations: ", nrow(results_dt))
  message("  Methods: ", paste(unique(results_dt$method), collapse = ", "))
  message("  Rates: ", paste(unique(results_dt$rate), collapse = ", "))
  message("  Replicates: ", length(replicates))
  
  # Test hypothesis: Does transmission enrichment increase with recombination?
  if (nrow(enrichment_summary) > 0) {
    message("\n[HYPOTHESIS TEST: Transmission Enrichment]")
    
    for (m in unique(enrichment_summary$method)) {
      enrich_data <- enrichment_summary[method == m & k == 50]
      
      if (nrow(enrich_data) > 1) {
        # Correlation between rate and enrichment
        cor_test <- cor.test(log10(enrich_data$rate), 
                            enrich_data$enrichment_mean,
                            method = "spearman")
        
        message("  ", m, " (top 50):")
        message("    Baseline enrichment: ", 
                round(enrich_data[rate == min(rate), enrichment_mean], 2), "x")
        message("    Highest rate enrichment: ",
                round(enrich_data[rate == max(rate), enrichment_mean], 2), "x")
        message("    Correlation (rate vs enrichment): rho = ",
                round(cor_test$estimate, 3), 
                ", p = ", format.pval(cor_test$p.value, digits = 3))
      }
    }
  }
  
  message("\n  Output directory: ", CONFIG$OUTDIR)
  message("\n", rep("=", 70), "\n")
  
  return(list(
    results = results_dt,
    agg_results = agg_results,
    detectability = detectability_dt,
    detect_summary = detect_summary,
    enrichment = enrichment_dt,
    enrichment_summary = enrichment_summary,
    merged_data = merged_data_store
  ))
}

# ============================================================================
#                     VISUALIZATION FUNCTIONS
# ============================================================================

source_code_for_plots <- function(agg_results, detect_summary, enrichment_summary) {
  
  message("\n  Generating plots...")
  
  theme_set(theme_bw(base_size = 14))
  
  custom_theme <- theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black", face = "bold"),
    plot.title = element_text(size = 15, color = "black", face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm"),
    strip.text = element_text(size = 13, color = "black", face = "bold"),
    panel.grid.minor = element_blank()
  )
  
  # ------------------------------------------------------------------------
  # Plot 1: Overall vs Transmission AUPR
  # ------------------------------------------------------------------------
  message("    [1/6] Overall vs Transmission AUPR comparison...")
  
  plot_data_overall <- agg_results[, .(rate, method, 
                                       AUPR = aupr_overall_mean,
                                       SE = aupr_overall_se,
                                       Task = "Overall Relatedness")]
  
  plot_data_trans <- agg_results[, .(rate, method,
                                     AUPR = aupr_transmission_mean,
                                     SE = aupr_transmission_se,
                                     Task = "Transmission Detection")]
  
  plot_data <- rbind(plot_data_overall, plot_data_trans)
  
  p1 <- ggplot(plot_data, aes(x = rate, y = AUPR, color = method, 
                              linetype = Task, group = interaction(method, Task))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = AUPR - SE, ymax = AUPR + SE), 
                  width = 0.05, linewidth = 0.7) +
    scale_x_log10(breaks = unique(plot_data$rate),
                  labels = function(x) format(x, scientific = TRUE)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_linetype_manual(values = c("Overall Relatedness" = "solid",
                                    "Transmission Detection" = "dashed")) +
    labs(title = "Method Performance: Overall vs Transmission Detection",
         subtitle = "Overall = detect ≤25 gen; Transmission = detect ≤5 gen",
         x = "Recombination Rate (bp⁻¹, log scale)",
         y = "AUPR (mean ± SE)") +
    custom_theme +
    guides(color = guide_legend(order = 1),
           linetype = guide_legend(order = 2))
  
  ggsave(file.path(CONFIG$OUTDIR, "figures", "01_overall_vs_transmission_aupr.png"),
         p1, width = 10, height = 6, dpi = 300)
  
  # ------------------------------------------------------------------------
  # Plot 2: Transmission Discrimination (among related pairs)
  # ------------------------------------------------------------------------
  message("    [2/6] Transmission discrimination AUPR...")
  
  p2 <- ggplot(agg_results, aes(x = rate, y = aupr_discrimination_mean, 
                                color = method, group = method)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = aupr_discrimination_mean - aupr_discrimination_se,
                      ymax = aupr_discrimination_mean + aupr_discrimination_se),
                  width = 0.05, linewidth = 0.8) +
    scale_x_log10(breaks = unique(agg_results$rate),
                  labels = function(x) format(x, scientific = TRUE)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Transmission Discrimination Among Related Pairs",
         subtitle = "Task: Distinguish transmission (≤5 gen) from ancestry (6-25 gen)",
         x = "Recombination Rate (bp⁻¹, log scale)",
         y = "AUPR (mean ± SE)") +
    custom_theme
  
  ggsave(file.path(CONFIG$OUTDIR, "figures", "02_transmission_discrimination_aupr.png"),
         p2, width = 10, height = 6, dpi = 300)
  
  # ------------------------------------------------------------------------
  # Plot 3: Detectability by Relationship Type
  # ------------------------------------------------------------------------
  message("    [3/6] Detectability by relationship...")
  
  p3 <- ggplot(detect_summary[relationship != "unrelated"], 
               aes(x = rate, y = detectability_mean, 
                   color = relationship, group = relationship)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = detectability_mean - detectability_se,
                      ymax = detectability_mean + detectability_se),
                  width = 0.05, linewidth = 0.8) +
    scale_x_log10(breaks = unique(detect_summary$rate),
                  labels = function(x) format(x, scientific = TRUE)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       labels = scales::percent_format()) +
    scale_color_viridis_d(option = "D", end = 0.8) +
    labs(title = "IBD Detectability by Relationship Type",
         subtitle = paste0("Threshold: ", CONFIG$ALPHA_DETECT * 100, 
                          "% FPR among unrelated pairs"),
         x = "Recombination Rate (bp⁻¹, log scale)",
         y = "Detectability (% with IBD ≥ threshold)") +
    custom_theme +
    theme(legend.position = "right")
  
  ggsave(file.path(CONFIG$OUTDIR, "figures", "03_detectability_by_relationship.png"),
         p3, width = 10, height = 6, dpi = 300)
  
  # ------------------------------------------------------------------------
  # Plot 4: Transmission Enrichment in Top K
  # ------------------------------------------------------------------------
  if (nrow(enrichment_summary) > 0) {
    message("    [4/6] Transmission enrichment...")
    
    p4 <- ggplot(enrichment_summary[k %in% c(10, 25, 50)],
                 aes(x = rate, y = enrichment_mean, 
                     color = method, linetype = factor(k), 
                     group = interaction(method, k))) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = enrichment_mean - enrichment_se,
                        ymax = enrichment_mean + enrichment_se),
                    width = 0.05, linewidth = 0.7) +
      scale_x_log10(breaks = unique(enrichment_summary$rate),
                    labels = function(x) format(x, scientific = TRUE)) +
      labs(title = "Transmission Enrichment in Top-Ranked Pairs",
           subtitle = "Enrichment = (observed % transmission) / (baseline % transmission)",
           x = "Recombination Rate (bp⁻¹, log scale)",
           y = "Enrichment Factor") +
      custom_theme +
      guides(color = guide_legend(title = "Method", order = 1),
             linetype = guide_legend(title = "Top K", order = 2))
    
    ggsave(file.path(CONFIG$OUTDIR, "figures", "04_transmission_enrichment_topk.png"),
           p4, width = 10, height = 6, dpi = 300)
  }
  
  # ------------------------------------------------------------------------
  # Plot 5: Heatmap of Performance Metrics
  # ------------------------------------------------------------------------
  message("    [5/6] Performance heatmap...")
  
  heatmap_data <- melt(agg_results,
                       id.vars = c("rate", "method"),
                       measure.vars = c("aupr_overall_mean", 
                                       "aupr_transmission_mean",
                                       "aupr_discrimination_mean"),
                       variable.name = "metric",
                       value.name = "aupr")
  
  heatmap_data[, metric := fcase(
    metric == "aupr_overall_mean", "Overall Relatedness",
    metric == "aupr_transmission_mean", "Transmission Detection",
    metric == "aupr_discrimination_mean", "Transmission Discrimination"
  )]
  
  p5 <- ggplot(heatmap_data, aes(x = factor(rate), y = method, fill = aupr)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = round(aupr, 2)), size = 3.5, fontface = "bold") +
    facet_wrap(~ metric, ncol = 1) +
    scale_fill_viridis_c(option = "plasma", limits = c(0, 1),
                        na.value = "gray90") +
    labs(title = "Performance Heatmap: AUPR Across Tasks",
         x = "Recombination Rate (bp⁻¹)",
         y = "Method",
         fill = "AUPR") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  ggsave(file.path(CONFIG$OUTDIR, "figures", "05_performance_heatmap.png"),
         p5, width = 10, height = 8, dpi = 300)
  
  # ------------------------------------------------------------------------
  # Plot 6: Composite Figure
  # ------------------------------------------------------------------------
  message("    [6/6] Creating composite figure...")
  
  composite <- (p1 / p2) | p3
  composite <- composite + plot_annotation(
    title = "Transmission Discrimination Under Recombination",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
  
  ggsave(file.path(CONFIG$OUTDIR, "figures", "00_composite_main_figure.png"),
         composite, width = 16, height = 10, dpi = 300)
  
  message("    ✓ All plots saved")
}

# ============================================================================
# RUN ANALYSIS
# ============================================================================

if (interactive() || !exists("analysis_results")) {
  analysis_results <- run_transmission_analysis()
  
  # Save output list
  message(" Save output list...")
  saveRDS(analysis_results, file.path(CONFIG$OUTDIR, "tables", "recombination_sweeps.rds"))
}
