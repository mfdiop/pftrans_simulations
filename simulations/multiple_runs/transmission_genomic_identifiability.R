# ============================================================================
# transmission_genomic_identifiability.R
# ============================================================================
#
# TITLE:   Quantifying the Limits of Genomic Identifiability Under
#          Realistic Biology and Surveillance
#
# GOAL:    Evaluate how recombination rate affects the ability of genomic
#          methods (IBD, IBS, Phylogenetic) to identify transmission events,
#          separating biological detectability from method discriminability.
#
# KEY HYPOTHESES:
#   H1: IBD detection sensitivity decreases with recombination rate
#       (shorter segments → fewer detectable pairs)
#   H2: Transmission discrimination may improve with recombination
#       (old ancestry IBD erodes faster → better signal/noise for recent events)
#   H3: Among detectable pairs, recent transmission is enriched at
#       high recombination rates relative to background ancestry
#
# GROUND TRUTH STRATA (generations apart):
#   - Transmission       : gen_distance ≤ GEN_TRANSMISSION  (e.g. ≤ 5)
#   - Recent ancestry    : GEN_TRANSMISSION < gen ≤ GEN_RECENT_ANCESTRY (e.g. 6–15)
#   - Distant ancestry   : GEN_RECENT_ANCESTRY < gen ≤ GEN_DISTANT_ANCESTRY (e.g. 16–25)
#   - Unrelated          : gen_distance > GEN_DISTANT_ANCESTRY (e.g. > 25)
#
# TASKS EVALUATED (per method, per recombination rate):
#   Task 1 – Overall relatedness     : ≤ 25 gen vs > 25 gen
#   Task 2 – Transmission detection  : ≤  5 gen vs everything else
#   Task 3 – Transmission discrimination (among related only):
#             ≤  5 gen vs 6–25 gen
#   Task 4 – Transmission enrichment in top-K ranked pairs
#
# DESIGN PRINCIPLES (from v1 corrections):
#   - Detectability threshold calibrated ONCE at the baseline rate only,
#     then held fixed for all higher rates → fair cross-rate comparison
#   - Classification thresholds (for confusion matrix) calibrated ONCE
#     at baseline via Youden's index, held fixed across rates
#   - Detectability reported alongside AUPR but NOT used to filter pairs;
#     evaluation always uses full genealogical ground truth
#   - Delta-AUPR reported relative to baseline rate
#
# METHODS EVALUATED:
#   IBD  – HMM-based IBD inference  (score = total IBD bp or proportion)
#   IBS  – Identity-by-state matrix (score = proportion shared alleles)
#   Phylo– Patristic distance from ML tree (converted to similarity via
#          exponential decay kernel)
#
# Authors : [Your Name]
# Date    : 2025-02-08 (unified version)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(PRROC)
  library(pROC)
  library(ape)
  library(phangorn)
  library(gridExtra)
  library(patchwork)
  library(viridis)
  library(writexl)
})

set.seed(42)

# ============================================================================
# SECTION 1: CONFIGURATION
# ============================================================================

CONFIG <- list(

  # --- Paths ------------------------------------------------------------------
  ROOT_DIR        = "simulations/multiple_runs/metrics",
  INFERRED_SUBPATH= "inferred",
  PHYLO_SUBPATH   = "phylo_results",
  OUTDIR          = "simulations/multiple_runs/genomic_evaluation",  # genomic_identifiability_evaluation

  # --- Simulation parameters --------------------------------------------------
  REC_RATES  = c("1e-09", "1e-08", "1e-07", "1e-06"),  # bp⁻¹, ascending
  GENOME_BP  = 640851,
  CHR        = "chr1",

  # --- Ground truth generation thresholds ------------------------------------
  GEN_TRANSMISSION    = 5,   # ≤ this  → "transmission" (positive for Task 2 & 3)
  GEN_RECENT_ANCESTRY = 15,  # ≤ this  → "recent ancestry"
  GEN_DISTANT_ANCESTRY= 25,  # ≤ this  → "distant ancestry" (positive for Task 1)

  # --- Detectability ----------------------------------------------------------
  ALPHA_DETECT = 0.01,  # 1% FPR among unrelated pairs → detectability threshold

  # --- Evaluation -------------------------------------------------------------
  TOP_K = c(5, 10, 25, 50, 100),

  # --- Output -----------------------------------------------------------------
  SAVE_PLOTS  = TRUE,
  SAVE_TABLES = TRUE
)

# ============================================================================
# SECTION 2: UTILITY FUNCTIONS
# ============================================================================

#' Create canonical pair key ensuring id_a ≤ id_b (avoids duplicate pairs)
canonical_pair <- function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_character_)
  if (a <= b) paste(a, b, sep = "--") else paste(b, a, sep = "--")
}

#' Safe data.table file read with informative warning on failure
safe_fread <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(fread(fp), error = function(e) {
    warning("fread failed: ", fp, " — ", e$message); NULL
  })
}

#' Safe RDS read with informative warning on failure
safe_readRDS <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(readRDS(fp), error = function(e) {
    warning("readRDS failed: ", fp, " — ", e$message); NULL
  })
}

#' Min-max normalization to [0, 1]
norm_to_01 <- function(x) {
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

#' Rank normalization: higher value → higher normalized score
#' (corrected: no minus sign — higher similarity gets higher rank)
rank_norm <- function(x) {
  if (all(is.na(x))) return(x)
  r <- rank(x, ties.method = "average", na.last = "keep")
  (r - 1) / (max(r, na.rm = TRUE) - 1)
}

# ============================================================================
# SECTION 3: DATA LOADING FUNCTIONS
# ============================================================================

#' Resolve the per-replicate per-rate subdirectory
#' Tries the canonical naming convention first, then falls back to grep.
resolve_run_dir <- function(rep_dir, rate_str) {
  repnum    <- gsub("rep", "", basename(rep_dir))
  candidate <- file.path(rep_dir, paste0("run", repnum, "_rec", rate_str, "_chr1"))
  if (dir.exists(candidate)) return(candidate)
  alts <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
  alts <- alts[grepl(paste0("rec", rate_str), alts)]
  if (length(alts) >= 1) return(alts[1])
  return(rep_dir)
}

# ----------------------------------------------------------------------------

#' Load ground truth IBD summary for one replicate × rate
#'
#' Returns: data.table(pair_key, id1, id2, total_ibd_bp, ibd_prop,
#'                     max_seg_bp, gen_distance)
load_ground_truth <- function(rep_dir, rate_str,
                              genome_bp = CONFIG$GENOME_BP) {
  run_dir <- resolve_run_dir(rep_dir, rate_str)
  fp      <- file.path(run_dir, "true_ibd_summary.tsv")
  dt      <- safe_fread(fp)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  # Standardize sample IDs
  dt <- dt %>%
    mutate(Id1 = paste0("tsk_", Id1),
           Id2 = paste0("tsk_", Id2)) %>%
    as.data.table()

  # Standardize ID column names
  id_cols <- intersect(c("id1", "id2", "Id1", "Id2"), names(dt))
  if (length(id_cols) >= 2)
    setnames(dt, old = id_cols[1:2], new = c("id1", "id2"), skip_absent = TRUE)

  # Canonical pair key
  if (!("pair_key" %in% names(dt)))
    dt[, pair_key := mapply(canonical_pair, id1, id2)]

  # IBD proportion: use true_ibd_prop if present, else compute from bp
  if (!("total_ibd_bp" %in% names(dt))) dt[, total_ibd_bp := NA_real_]
  if ("true_ibd_prop" %in% names(dt)) {
    dt[, ibd_prop := true_ibd_prop]
  } else {
    dt[, ibd_prop := total_ibd_bp / genome_bp]
  }

  # Sanity-check: ibd_prop must be in [0, 1]
  max_ibd <- max(dt$ibd_prop, na.rm = TRUE)
  if (!is.na(max_ibd) && max_ibd > 1.0) {
    warning("  ⚠ ibd_prop > 1 detected (max = ", round(max_ibd, 3),
            "), normalizing")
    dt[, ibd_prop := pmin(ibd_prop / max_ibd, 1.0)]
  }

  # Generation distance
  if (!("gen_distance" %in% names(dt))) {
    possible <- intersect(c("generations", "t_generations", "min_tmrca"),
                          names(dt))
    if (length(possible) >= 1)
      setnames(dt, old = possible[1], new = "gen_distance", skip_absent = TRUE)
  }

  # Max segment
  if ("max_segment_bp" %in% names(dt)) {
    dt[, max_seg_bp := max_segment_bp]
  } else {
    dt[, max_seg_bp := NA_real_]
  }

  dt[, .(pair_key, id1, id2, total_ibd_bp, ibd_prop, max_seg_bp, gen_distance)]
}

# ----------------------------------------------------------------------------

#' Load HMM-IBD inferred output for one replicate × rate
#'
#' Returns: data.table(pair_key, id1, id2, score)  score = total IBD (higher = more related)
load_ibd_hmm <- function(rep_dir, rate_str) {
  run_dir <- resolve_run_dir(rep_dir, rate_str)
  fp      <- file.path(run_dir, "inferred_ibd_hmm.rds")
  dt      <- as.data.table(safe_readRDS(fp))
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  id_cols <- intersect(c("id1", "id2", "Id1", "Id2"), names(dt))
  if (length(id_cols) >= 2)
    setnames(dt, old = id_cols[1:2], new = c("id1", "id2"), skip_absent = TRUE)

  score_candidates <- c("total_ibd_bp", "total_ibd", "ibd_prop", "score", "posterior")
  sc <- intersect(score_candidates, names(dt))
  if (length(sc) > 0) {
    dt[, score := as.numeric(get(sc[1]))]
  } else if ("n_segments" %in% names(dt)) {
    dt[, score := as.numeric(n_segments)]
  } else {
    dt[, score := 1]
  }

  if (!("pair_key" %in% names(dt)))
    dt[, pair_key := mapply(canonical_pair, id1, id2)]

  unique(dt[, .(pair_key, id1, id2, score)])
}

# ----------------------------------------------------------------------------

#' Load IBS matrix for one replicate × rate
#'
#' Returns: data.table(pair_key, id1, id2, score)  score = IBS proportion (higher = more similar)
load_ibs <- function(rep_dir, rate_str) {
  run_dir <- resolve_run_dir(rep_dir, rate_str)
  fp      <- file.path(run_dir, "ibs_matrix.rds")
  if (!file.exists(fp)) return(NULL)

  obj <- safe_readRDS(fp)
  if (is.null(obj)) return(NULL)

  mat <- as.matrix(obj)
  ids <- rownames(mat)

  if (is.null(ids)) {
    if (is.data.frame(obj)) {
      ids <- as.character(obj[[1]])
      mat <- as.matrix(obj[, -1])
      rownames(mat) <- ids; colnames(mat) <- ids
    } else return(NULL)
  }

  n    <- length(ids)
  rows <- vector("list", n * (n - 1) / 2)
  idx  <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      rows[[idx]] <- list(
        pair_key = canonical_pair(ids[i], ids[j]),
        id1      = ids[i],
        id2      = ids[j],
        score    = as.numeric(mat[i, j])
      )
      idx <- idx + 1
    }
  }

  unique(rbindlist(rows)[, .(pair_key, id1, id2, score)])
}

# ----------------------------------------------------------------------------

#' Load ML phylogenetic tree and convert patristic distances to similarity scores.
#'
#' Uses an exponential decay kernel: score = exp(-alpha * distance)
#' where alpha = 1 / median(distance).  This is biologically motivated
#' (similarity decays exponentially with evolutionary distance) and handles
#' skewed branch-length distributions better than linear rescaling.
#'
#' Returns: data.table(pair_key, id1, id2, score)  score ∈ (0, 1], higher = more related
load_phylo <- function(phylo_dir, rate_str, chr) {
  patt  <- paste0("rec", rate_str, "_", chr, ".treefile")
  files <- list.files(phylo_dir, full.names = TRUE, pattern = patt)

  tf <- files[grepl("\\.treefile$", files)]
  if (length(tf) == 0 && length(files) > 0) tf <- files[1]
  if (length(tf) == 0) return(NULL)

  message("      Reading tree: ", basename(tf[1]))

  tree <- tryCatch(ape::read.tree(tf[1]),
                   error = function(e) { warning("read.tree failed: ", e$message); NULL })
  if (is.null(tree)) return(NULL)

  message("      Computing patristic distances...")
  dist_mat <- tryCatch(ape::cophenetic.phylo(tree),
                       error = function(e) { warning("cophenetic failed: ", e$message); NULL })
  if (is.null(dist_mat)) return(NULL)

  ids <- rownames(dist_mat)
  n   <- length(ids)
  rows <- vector("list", n * (n - 1) / 2)
  idx  <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      rows[[idx]] <- list(
        pair_key = canonical_pair(ids[i], ids[j]),
        id1      = ids[i],
        id2      = ids[j],
        distance = as.numeric(dist_mat[i, j])
      )
      idx <- idx + 1
    }
  }
  dt <- rbindlist(rows)

  # Exponential decay kernel: exp(-d / median_d)
  # Handles skewed branch-length distributions; score → 1 as d → 0
  med_d <- median(dt$distance, na.rm = TRUE)
  if (is.na(med_d) || med_d == 0) med_d <- 1  # safety fallback
  alpha <- 1 / med_d
  dt[, score := exp(-alpha * distance)]

  message("      Patristic distance range : [",
          round(min(dt$distance, na.rm = TRUE), 6), ", ",
          round(max(dt$distance, na.rm = TRUE), 6), "]")
  message("      Similarity score range   : [",
          round(min(dt$score, na.rm = TRUE), 6), ", ",
          round(max(dt$score, na.rm = TRUE), 6), "]")

  unique(dt[, .(pair_key, id1, id2, score)])
}

# ============================================================================
# SECTION 4: GROUND TRUTH ANNOTATION
# ============================================================================

#' Annotate all pairs with epidemiological relationship strata and
#' binary ground truth labels for each evaluation task.
annotate_relationships <- function(dt,
                                   gen_transmission = CONFIG$GEN_TRANSMISSION,
                                   gen_recent       = CONFIG$GEN_RECENT_ANCESTRY,
                                   gen_distant      = CONFIG$GEN_DISTANT_ANCESTRY) {
  # Categorical strata
  dt[, relationship := fcase(
    gen_distance <= gen_transmission,                             "transmission",
    gen_distance >  gen_transmission & gen_distance <= gen_recent,"recent_ancestry",
    gen_distance >  gen_recent       & gen_distance <= gen_distant,"distant_ancestry",
    default = "unrelated"
  )]
  dt[, relationship := factor(relationship,
                               levels = c("transmission", "recent_ancestry",
                                          "distant_ancestry", "unrelated"))]

  # Binary ground truths
  dt[, truth_transmission    := as.integer(gen_distance <= gen_transmission)]  # Task 2
  dt[, truth_any_relatedness := as.integer(gen_distance <= gen_distant)]       # Task 1

  # Task 3 label: only defined for pairs that are related (≤ gen_distant)
  # 1 = transmission, 0 = ancestry background; NA = unrelated (excluded)
  dt[, is_transmission_vs_ancestry := fcase(
    gen_distance <= gen_transmission,                              1L,
    gen_distance >  gen_transmission & gen_distance <= gen_distant, 0L,
    default = NA_integer_
  )]

  message("  Relationship breakdown:")
  message("    Transmission      (≤ ", gen_transmission, " gen) : ",
          sum(dt$relationship == "transmission"))
  message("    Recent ancestry   (", gen_transmission + 1, "–", gen_recent, " gen) : ",
          sum(dt$relationship == "recent_ancestry"))
  message("    Distant ancestry  (", gen_recent + 1, "–", gen_distant, " gen) : ",
          sum(dt$relationship == "distant_ancestry"))
  message("    Unrelated         (> ", gen_distant, " gen) : ",
          sum(dt$relationship == "unrelated"))

  return(dt)
}

# ============================================================================
# SECTION 5: DETECTABILITY THRESHOLD
# ============================================================================

#' Compute the IBD detectability threshold at alpha FPR.
#'
#' Calibrated ONCE on the baseline (lowest) recombination rate and then
#' held fixed for all higher rates — this is the correct design for
#' measuring how detectability degrades as recombination increases.
#'
#' Logic:
#'   threshold = (1 - alpha) quantile of unrelated pairs' IBD distribution.
#'   If the threshold exceeds the 5th percentile of positive IBD
#'   (poor separation), we warn and use the 5th percentile of positives
#'   instead so the threshold remains meaningful.
#'
#' @param dt        data.table with ibd_prop and truth_any_relatedness columns
#' @param alpha     false positive rate (default 0.01)
#' @return          scalar threshold value
define_detectability_threshold <- function(dt,
                                           truth_col = "truth_any_relatedness",
                                           alpha     = CONFIG$ALPHA_DETECT) {
  message("  → Computing detectability threshold (alpha = ", alpha, ")")

  negatives <- dt[get(truth_col) == 0 & !is.na(ibd_prop), ibd_prop]
  positives <- dt[get(truth_col) == 1 & !is.na(ibd_prop), ibd_prop]

  if (length(negatives) == 0) {
    warning("  ⚠ No negatives found, returning threshold = 0")
    return(0)
  }

  message("    Negative IBD : min = ", round(min(negatives), 6),
          "  median = ", round(median(negatives), 6),
          "  99th% = ", round(quantile(negatives, 0.99), 6),
          "  max = ", round(max(negatives), 6))

  if (length(positives) > 0) {
    message("    Positive IBD : min = ", round(min(positives), 6),
            "  median = ", round(median(positives), 6),
            "  max = ", round(max(positives), 6))
  }

  threshold_neg <- quantile(negatives, probs = 1 - alpha, na.rm = TRUE)

  # Check for poor separation
  if (length(positives) > 0) {
    threshold_pos5 <- quantile(positives, probs = 0.05, na.rm = TRUE)
    if (threshold_neg > threshold_pos5) {
      warning("  ⚠ Poor separation: 99th% negatives (", round(threshold_neg, 6),
              ") > 5th% positives (", round(threshold_pos5, 6), ")")
      message("    → Falling back to 5th percentile of positives")
      threshold_neg <- threshold_pos5
    }
  }

  # Hard cap at 0.9 with fallback to 0.1 — avoids degenerate thresholds
  if (!is.na(threshold_neg) && threshold_neg > 0.9) {
    warning("  ⚠ Threshold > 0.9, capping at 0.1")
    threshold_neg <- 0.1
  }

  message("    ✓ Detectability threshold : ", round(as.numeric(threshold_neg), 6))
  return(as.numeric(threshold_neg))
}

# ============================================================================
# SECTION 6: CLASSIFICATION THRESHOLD CALIBRATION
# ============================================================================

#' Calibrate per-method classification thresholds using Youden's index on
#' the baseline recombination rate data.
#'
#' These thresholds are fixed and re-applied at all higher rates, allowing
#' fair comparisons of how confusion-matrix metrics (precision, recall, F1)
#' degrade as recombination increases.
#'
#' @param dt         merged data.table at baseline
#' @param score_cols character vector of score column names (e.g. "score_IBD")
#' @param truth_col  binary truth column to calibrate against
#' @return           named list of scalar thresholds
calibrate_classification_thresholds <- function(dt,
                                                 score_cols,
                                                 truth_col = "truth_any_relatedness") {
  message("  → Calibrating classification thresholds (Youden's index on baseline)")

  thresholds <- list()

  for (sc in score_cols) {
    if (!sc %in% names(dt)) {
      warning("  ⚠ Score column '", sc, "' not found, skipping"); next
    }

    valid  <- !is.na(dt[[sc]]) & !is.na(dt[[truth_col]])
    scores <- dt[[sc]][valid]
    truth  <- dt[[truth_col]][valid]

    if (length(unique(truth)) < 2) {
      warning("  ⚠ Only one class for '", sc, "', skipping"); next
    }

    roc_obj <- tryCatch(
      pROC::roc(truth, scores, quiet = TRUE),
      error = function(e) { warning("ROC failed for '", sc, "': ", e$message); NULL }
    )

    if (is.null(roc_obj)) {
      thresholds[[sc]] <- median(scores, na.rm = TRUE)
      next
    }

    coords <- pROC::coords(roc_obj, x = "best", best.method = "youden", transpose = FALSE)
    # thresholds[[sc]] <- as.numeric(coords$threshold)
    
    # coords may return multiple rows when ties exist at the Youden optimum,
    # and pROC appends sentinel rows with threshold = -Inf / Inf.
    # Keep only finite rows; if multiple remain, choose the one with the
    # highest Youden index (sensitivity + specificity - 1).
    if (is.data.frame(coords) && nrow(coords) > 1) {
      finite_rows <- coords[is.finite(coords$threshold), , drop = FALSE]
      if (nrow(finite_rows) == 0) {
        # Truly no finite optimum — fall back to median score
        warning("  ⚠ No finite Youden threshold for '", sc,
                "', using score median")
        thresholds[[sc]] <- median(scores, na.rm = TRUE)
        next
      }
      finite_rows$youden <- finite_rows$sensitivity + finite_rows$specificity - 1
      coords <- finite_rows[which.max(finite_rows$youden), , drop = FALSE]
    }
    
    # Enforce scalar — take first value as a safety net
    thresh_val <- as.numeric(coords$threshold[1])
    
    thresholds[[sc]] <- thresh_val

    message("    ✓ ", sc, " : threshold = ", round(coords$threshold, 6),
            "  (sensitivity = ", round(coords$sensitivity, 3),
            ", specificity = ", round(coords$specificity, 3), ")")
  }

  return(thresholds)
}

# ============================================================================
# SECTION 7: EVALUATION FUNCTIONS
# ============================================================================

#' Compute the full suite of metrics for one method × one binary task.
#'
#' Assumes all scores are similarities (higher = more related).
#' IBS and IBD scores already satisfy this; Phylo distances are converted
#' to similarities in load_phylo().
#'
#' @param dt          merged data.table
#' @param score_col   name of the score column
#' @param truth_col   name of the binary truth column
#' @param threshold   pre-calibrated classification threshold (or NULL)
#' @param ks          integer vector of K values for Precision@K
#' @return            named list of metric values
evaluate_binary_task <- function(dt,
                                 score_col,
                                 truth_col,
                                 threshold = NULL,
                                 ks        = CONFIG$TOP_K) {
  res <- list(score_col = score_col, truth_col = truth_col)

  valid <- !is.na(dt[[score_col]]) & !is.na(dt[[truth_col]])
  if (sum(valid) == 0) {
    res[c("aupr","auroc","spearman","brier","n_pos","n_neg","prevalence")] <-
      list(NA_real_, NA_real_, NA_real_, NA_real_, NA_integer_, NA_integer_, NA_real_)
    for (k in ks) {
      res[[paste0("prec_at_",  k)]] <- NA_real_
      res[[paste0("recall_at_",k)]] <- NA_real_
      res[[paste0("f1_at_",    k)]] <- NA_real_
    }
    return(res)
  }

  dt_eval <- dt[valid]
  scores  <- dt_eval[[score_col]]
  truth   <- dt_eval[[truth_col]]

  n_pos <- sum(truth == 1); n_neg <- sum(truth == 0)
  res$n_pos      <- n_pos
  res$n_neg      <- n_neg
  res$prevalence <- n_pos / (n_pos + n_neg)

  if (n_pos == 0 || n_neg == 0) {
    res[c("aupr","auroc","spearman","brier")] <-
      list(NA_real_, NA_real_, NA_real_, NA_real_)
    for (k in ks) {
      res[[paste0("prec_at_",  k)]] <- NA_real_
      res[[paste0("recall_at_",k)]] <- NA_real_
      res[[paste0("f1_at_",    k)]] <- NA_real_
    }
    return(res)
  }

  # ------ AUCPR ---------------------------------------------------------------
  pos_scores <- scores[truth == 1]
  neg_scores <- scores[truth == 0]
  pr_obj <- tryCatch(
    PRROC::pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores,
                    curve = TRUE),
    error = function(e) { warning("PR curve failed: ", e$message); NULL }
  )
  if (!is.null(pr_obj)) {
    res$aupr     <- pr_obj$auc.integral
    res$pr_curve <- as.data.table(pr_obj$curve) %>%
      setnames(c("recall", "precision", "threshold"))
  } else {
    res$aupr     <- NA_real_
    res$pr_curve <- NULL
  }

  # ------ AUROC --------------------------------------------------------------
  roc_obj <- tryCatch(
    pROC::roc(truth, scores, quiet = TRUE),
    error = function(e) { warning("ROC failed: ", e$message); NULL }
  )
  if (!is.null(roc_obj)) {
    res$auroc <- as.numeric(pROC::auc(roc_obj))
    # Also store Youden optimal point for reference
    coords <- tryCatch(pROC::coords(roc_obj, "best", best.method = "youden"),
                       error = function(e) NULL)
    if (!is.null(coords)) {
      res$youden_threshold   <- as.numeric(coords$threshold)
      res$youden_sensitivity <- as.numeric(coords$sensitivity)
      res$youden_specificity <- as.numeric(coords$specificity)
    }
  } else {
    res$auroc <- NA_real_
  }

  # ------ Spearman correlation with continuous IBD ---------------------------
  if ("ibd_prop" %in% names(dt_eval)) {
    res$spearman <- suppressWarnings(
      cor(scores, dt_eval$ibd_prop, method = "spearman", use = "complete.obs")
    )
  } else {
    res$spearman <- NA_real_
  }

  # ------ Brier score --------------------------------------------------------
  rng <- range(scores, na.rm = TRUE)
  pred_prob <- if (diff(rng) == 0) rep(0.5, length(scores)) else
    (scores - rng[1]) / diff(rng)
  res$brier <- mean((pred_prob - truth)^2, na.rm = TRUE)

  # ------ Precision@K, Recall@K, F1@K ----------------------------------------
  # Sort descending: higher score = more related
  dt_ranked <- dt_eval[order(-get(score_col))]

  for (k in ks) {
    if (nrow(dt_ranked) < k) {
      res[[paste0("prec_at_",  k)]] <- NA_real_
      res[[paste0("recall_at_",k)]] <- NA_real_
      res[[paste0("f1_at_",    k)]] <- NA_real_
      next
    }
    # Assign top-K as predicted positive, rest as negative
    dt_ranked[, pred := 0L]
    dt_ranked[seq_len(k), pred := 1L]

    TP <- sum(dt_ranked$pred == 1 & dt_ranked[[truth_col]] == 1)
    FP <- sum(dt_ranked$pred == 1 & dt_ranked[[truth_col]] == 0)
    FN <- sum(dt_ranked$pred == 0 & dt_ranked[[truth_col]] == 1)

    prec   <- ifelse(TP + FP > 0, TP / (TP + FP), NA_real_)
    recall <- ifelse(TP + FN > 0, TP / (TP + FN), NA_real_)
    f1     <- ifelse(!is.na(prec + recall) && (prec + recall) > 0,
                     2 * prec * recall / (prec + recall), NA_real_)

    res[[paste0("prec_at_",   k)]] <- prec
    res[[paste0("recall_at_", k)]] <- recall
    res[[paste0("f1_at_",     k)]] <- f1

    # Reset pred column for next k
    dt_ranked[, pred := NULL]
  }

  # ------ Confusion matrix at calibrated fixed threshold ---------------------
  if (!is.null(threshold) && !is.na(threshold)) {
    pred_thresh <- as.integer(scores >= threshold)
    TP  <- sum(pred_thresh == 1 & truth == 1)
    FP  <- sum(pred_thresh == 1 & truth == 0)
    TN  <- sum(pred_thresh == 0 & truth == 0)
    FN  <- sum(pred_thresh == 0 & truth == 1)
    pr  <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
    rc  <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
    f1  <- if (!is.na(pr + rc) && (pr + rc) > 0)
      2 * pr * rc / (pr + rc) else NA_real_
    res$thresh_TP        <- TP
    res$thresh_FP        <- FP
    res$thresh_TN        <- TN
    res$thresh_FN        <- FN
    res$thresh_precision <- pr
    res$thresh_recall    <- rc
    res$thresh_f1        <- f1
  }

  return(res)
}

# ----------------------------------------------------------------------------

#' Evaluate how well top-ranked pairs are enriched for true transmission events.
#'
#' Enrichment = (observed proportion of transmission pairs in top-K) /
#'              (baseline proportion of transmission pairs in the full dataset)
#' Enrichment > 1 means the method is recovering transmission pairs at better
#' than chance; enrichment = 1 is random; < 1 is anti-enriched.
evaluate_transmission_enrichment <- function(dt, score_col,
                                              ks = CONFIG$TOP_K,
                                              gen_trans = CONFIG$GEN_TRANSMISSION,
                                              gen_dist  = CONFIG$GEN_DISTANT_ANCESTRY) {
  valid   <- !is.na(dt[[score_col]])
  dt_eval <- dt[valid]
  if (nrow(dt_eval) == 0) return(NULL)

  dt_ranked  <- dt_eval[order(-get(score_col))]
  baseline_p <- mean(dt_eval$gen_distance <= gen_trans, na.rm = TRUE)

  results <- lapply(ks, function(k) {
    if (nrow(dt_ranked) < k) return(NULL)
    topk <- head(dt_ranked, k)
    n_trans    <- sum(topk$gen_distance <= gen_trans,   na.rm = TRUE)
    n_ancestry <- sum(topk$gen_distance >  gen_trans &
                      topk$gen_distance <= gen_dist,    na.rm = TRUE)
    n_unrel    <- sum(topk$gen_distance >  gen_dist,    na.rm = TRUE)
    obs_p      <- n_trans / k
    data.table(
      k                   = k,
      n_transmission      = n_trans,
      n_ancestry          = n_ancestry,
      n_unrelated         = n_unrel,
      prop_transmission   = obs_p,
      baseline_proportion = baseline_p,
      enrichment          = if (baseline_p > 0) obs_p / baseline_p else NA_real_
    )
  })

  rbindlist(Filter(Negate(is.null), results))
}

# ----------------------------------------------------------------------------

#' Evaluate IBD as a sparse DETECTOR (not a ranker).
#'
#' Strategy rationale
#' ------------------
#' Two valid ways to evaluate the IBD HMM exist and answer different questions:
#'
#' 1. RANKER evaluation  (used in evaluate_binary_task after zero-fill):
#'    Universe = all N*(N-1)/2 pairs.  Undetected pairs get score = 0.
#'    Metric: AUPR/AUROC over the full population.
#'    Question: "How well does IBD rank ALL pairs by relatedness?"
#'
#' 2. DETECTOR evaluation  (this function):
#'    Universe = only the pairs the HMM reported (score > 0).
#'    Metric: PPV (precision), recall, F1 of the detection call itself.
#'    Question: "Among the pairs IBD chose to flag, how accurate is it?"
#'
#' These are complementary.  At low recombination, IBD flags almost everything
#' (low specificity as a detector but high recall).  At high recombination, IBD
#' flags fewer pairs but those it does flag are more precise.
#'
#' @param merged     merged data.table with score_IBD (zero-filled) and truth columns
#' @param truth_col  binary ground-truth column to evaluate against
#' @return           data.table with detection PPV, recall, F1, n_detected
evaluate_ibd_detection_performance <- function(merged,
                                               truth_col = "truth_any_relatedness") {
  if (!"score_IBD" %in% names(merged)) return(NULL)
  if (!truth_col  %in% names(merged)) return(NULL)
  
  valid    <- !is.na(merged[[truth_col]])
  dt       <- merged[valid]
  
  # Detected = HMM reported at least one segment (score > 0 after zero-fill)
  detected <- dt[score_IBD > 0]
  n_total  <- nrow(dt)
  n_det    <- nrow(detected)
  
  if (n_det == 0) {
    return(data.table(
      n_total = n_total, n_detected = 0L,
      detection_rate = 0, ppv = NA_real_,
      detection_recall = NA_real_, detection_f1 = NA_real_,
      n_true_pos_detected = NA_integer_,
      n_true_pos_total    = sum(dt[[truth_col]] == 1, na.rm = TRUE)
    ))
  }
  
  n_true_pos_total    <- sum(dt[[truth_col]] == 1, na.rm = TRUE)
  n_true_pos_detected <- sum(detected[[truth_col]] == 1, na.rm = TRUE)
  n_false_pos_detected<- n_det - n_true_pos_detected
  
  ppv    <- n_true_pos_detected / n_det                                     # precision
  recall <- if (n_true_pos_total > 0)
    n_true_pos_detected / n_true_pos_total else NA_real_         # recall
  f1     <- if (!is.na(ppv) && !is.na(recall) && (ppv + recall) > 0)
    2 * ppv * recall / (ppv + recall) else NA_real_
  
  data.table(
    n_total              = n_total,
    n_detected           = n_det,
    detection_rate       = n_det / n_total,
    n_true_pos_detected  = n_true_pos_detected,
    n_true_pos_total     = n_true_pos_total,
    ppv                  = ppv,
    detection_recall     = recall,
    detection_f1         = f1
  )
}

# ============================================================================
# SECTION 8: MAIN ANALYSIS PIPELINE
# ============================================================================

run_genomic_identifiability_analysis <- function() {

  message("\n", strrep("=", 75))
  message("  GENOMIC IDENTIFIABILITY UNDER RECOMBINATION — UNIFIED PIPELINE")
  message(strrep("=", 75))

  # --- Setup output directories -----------------------------------------------
  dir.create(CONFIG$OUTDIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(CONFIG$OUTDIR, "tables"),  showWarnings = FALSE)
  dir.create(file.path(CONFIG$OUTDIR, "figures"), showWarnings = FALSE)

  # --- Rate lookup: label ("1e09") → numeric (1e-9) ---------------------------
  rec_rates_label <- gsub("-", "", CONFIG$REC_RATES)   # "1e-09" → "1e09"
  rate_lookup     <- setNames(as.numeric(CONFIG$REC_RATES), rec_rates_label)
  baseline_rate   <- min(as.numeric(CONFIG$REC_RATES))

  message("\n[RATE LOOKUP TABLE]")
  for (i in seq_along(rec_rates_label))
    message("  ", rec_rates_label[i], " → ", rate_lookup[rec_rates_label[i]], " bp⁻¹")

  # --- Auto-detect replicates -------------------------------------------------
  inferred_root <- file.path(CONFIG$ROOT_DIR, CONFIG$INFERRED_SUBPATH)
  replicates    <- list.dirs(inferred_root, recursive = FALSE, full.names = FALSE)
  replicates    <- replicates[grepl("^rep", replicates)]

  message("\n[CONFIG]")
  message("  Replicates            : ", paste(replicates, collapse = ", "))
  message("  Recombination rates   : ", paste(CONFIG$REC_RATES, collapse = ", "))
  message("  Transmission threshold: ≤ ", CONFIG$GEN_TRANSMISSION, " generations")
  message("  Distant ancestry      : ≤ ", CONFIG$GEN_DISTANT_ANCESTRY, " generations")
  message("  Detectability alpha   : ", CONFIG$ALPHA_DETECT)
  message("  Output directory      : ", CONFIG$OUTDIR)

  # --- Storage ----------------------------------------------------------------
  all_results          <- list()
  all_enrichment       <- list()
  all_detectability    <- list()
  all_ibd_detection    <- list()   # IBD as sparse detector (PPV / recall of flagged pairs)
  merged_data_store    <- list()
  all_pr_curve          <- list()

  # Global calibration values — set at baseline, reused for all higher rates
  ibd_eps_global          <- NULL   # detectability threshold
  global_class_thresholds <- NULL   # Youden classification thresholds

  # ============================================================================
  # STEP 1: DATA LOADING AND EVALUATION
  # ============================================================================

  message("\n", strrep("=", 75))
  message("  STEP 1: LOADING DATA AND COMPUTING METRICS")
  message(strrep("=", 75))

  for (rep_idx in seq_along(replicates)) {
    
    prcurve <- NULL
    rep     <- replicates[rep_idx]
    rep_dir <- file.path(inferred_root, rep)

    message("\n[", rep_idx, "/", length(replicates), "] Replicate: ", rep)

    if (!dir.exists(rep_dir)) {
      warning("  ⚠ Directory not found: ", rep_dir); next
    }

    for (rate_idx in seq_along(rec_rates_label)) {
      rate_label   <- rec_rates_label[rate_idx]
      rate_numeric <- rate_lookup[rate_label]
      is_baseline  <- (rate_numeric == baseline_rate)

      message("\n  Rate ", rate_idx, "/", length(rec_rates_label), ": ",
              rate_label, " (", rate_numeric, " bp⁻¹)",
              if (is_baseline) "  [BASELINE — thresholds calibrated here]" else "")

      # ------------------ Load ground truth ------------------------------------
      message("  Loading ground truth...")
      truth <- load_ground_truth(rep_dir, rate_label)

      if (is.null(truth) || nrow(truth) == 0) {
        warning("    ⚠ No ground truth data, skipping"); next
      }
      message("    ✓ ", nrow(truth), " pairs loaded")

      # Annotate relationship strata and binary labels
      truth <- annotate_relationships(truth)

      # # ------------------ Calibrate thresholds ONCE at baseline ---------------
      # if (is_baseline) {
      #   message("  Calibrating detectability threshold (baseline only)...")
      #   ibd_eps_global <- define_detectability_threshold(truth)
      # }
      # 
      # # Mark detectable pairs (descriptive only — not used to filter evaluation)
      # truth[, detectable := as.integer(ibd_prop >= ibd_eps_global)]
      # 
      # # ------------------ Detectability by relationship type ------------------
      # detect_stats <- truth[, .(
      #   n_pairs      = .N,
      #   n_detectable = sum(detectable == 1, na.rm = TRUE),
      #   detectability= mean(detectable,     na.rm = TRUE),
      #   mean_ibd     = mean(ibd_prop,        na.rm = TRUE),
      #   median_ibd   = median(ibd_prop,      na.rm = TRUE)
      # ), by = relationship]
      # 
      # detect_stats[, `:=`(replicate   = rep,
      #                      rate        = rate_numeric,
      #                      rate_label  = rate_label)]
      # all_detectability[[length(all_detectability) + 1]] <- detect_stats
      # 
      # message("  Detectability by relationship:")
      # for (i in seq_len(nrow(detect_stats))) {
      #   message("    ", detect_stats$relationship[i], " : ",
      #           round(100 * detect_stats$detectability[i], 1), "% (",
      #           detect_stats$n_detectable[i], "/", detect_stats$n_pairs[i], ")")
      # }

      # ------------------ Load inference methods ------------------------------
      message("  Loading method predictions...")
      ibd   <- load_ibd_hmm(rep_dir, rate_label)
      if (!is.null(ibd))   message("    ✓ IBD  : ", nrow(ibd),  " pairs")

      ibs   <- load_ibs(rep_dir, rate_label)
      if (!is.null(ibs))   message("    ✓ IBS  : ", nrow(ibs),  " pairs")

      phylo <- NULL
      phylo_dir <- file.path(CONFIG$ROOT_DIR, CONFIG$PHYLO_SUBPATH, rep)
      if (dir.exists(phylo_dir)) {
        phylo <- load_phylo(phylo_dir, rate_label, CONFIG$CHR)
        if (!is.null(phylo)) message("    ✓ Phylo: ", nrow(phylo), " pairs")
      }

      # ------------------ Merge all methods with ground truth -----------------
      message("  Merging...")
      merged <- truth[, .(pair_key, id1, id2, ibd_prop, gen_distance,
                           relationship, truth_transmission, truth_any_relatedness,
                           is_transmission_vs_ancestry, max_seg_bp)] # detectable, 

      if (!is.null(ibd))   merged <- merge(merged, ibd[,  .(pair_key, score_IBD   = score)], by = "pair_key", all.x = TRUE)
      if (!is.null(ibs))   merged <- merge(merged, ibs[,  .(pair_key, score_IBS   = score)], by = "pair_key", all.x = TRUE)
      if (!is.null(phylo)) merged <- merge(merged, phylo[, .(pair_key, score_Phylo = score)], by = "pair_key", all.x = TRUE)
      
      # -----------------------------------------------------------------------
      # CRITICAL: IBD HMM only outputs DETECTED pairs (those with ≥1 segment).
      # Pairs absent from the HMM file were assigned zero segments by the HMM
      # — that is a real score of 0, not missing data.
      # Leaving them as NA drops all true negatives → AUROC = 0.5 exactly.
      # IBS and Phylo are dense matrices (all pairs), so their NAs are genuine.
      # -----------------------------------------------------------------------
      if ("score_IBD" %in% names(merged)) {
        n_ibd_na <- sum(is.na(merged$score_IBD))
        merged[is.na(score_IBD), score_IBD := 0]
        if (n_ibd_na > 0)
          message("    ✓ IBD: ", n_ibd_na, " undetected pairs assigned score = 0")
      }

      message("    ✓ Merged table: ", nrow(merged), " rows × ", ncol(merged), " columns")
      # merged_data_store[[paste(rep, rate_label, sep = "__")]] <- merged
      # 
      # # ------------------ Calibrate classification thresholds at baseline -----
      # score_cols <- grep("^score_", names(merged), value = TRUE)
      # 
      # if (is_baseline) {
      #   message("  Calibrating classification thresholds (baseline only)...")
      #   global_class_thresholds <- calibrate_classification_thresholds(
      #     merged, score_cols, truth_col = "truth_any_relatedness"
      #   )
      # }
      
      # ------------------------------------------------------------------------
      # ------------------ Calibrate ALL thresholds ONCE at baseline -----------
      # Both detectability and classification thresholds use the merged table
      # (post zero-fill) so that IBD score distribution is correct.
      score_cols <- grep("^score_", names(merged), value = TRUE)
      
      if (is_baseline) {
        # Detectability: uses inferred IBD score distribution, not tree-sequence ibd_prop
        message("  Calibrating detectability threshold (baseline only)...")
        ibd_eps_global <- define_detectability_threshold(merged)
        
        # Classification thresholds: Youden's index on each method
        message("  Calibrating classification thresholds (baseline only)...")
        global_class_thresholds <- calibrate_classification_thresholds(
          merged, score_cols, truth_col = "truth_any_relatedness"
        )
      }
      
      # ------------------ Mark detectable pairs (descriptive only) -----------
      # Uses the inferred IBD score threshold (not tree-sequence ibd_prop)
      if ("score_IBD" %in% names(merged)) {
        merged[, detectable := as.integer(score_IBD > ibd_eps_global)]
      } else {
        merged[, detectable := 0L]
      }
      
      # ------------------ Detectability by relationship type ------------------
      detect_stats <- merged[, .(
        n_pairs       = .N,
        n_detectable  = sum(detectable == 1, na.rm = TRUE),
        detectability = mean(detectable,     na.rm = TRUE),
        mean_ibd_score  = mean(score_IBD,    na.rm = TRUE),
        median_ibd_score= median(score_IBD,  na.rm = TRUE)
      ), by = relationship]
      
      detect_stats[, `:=`(replicate  = rep,
                          rate       = rate_numeric,
                          rate_label = rate_label)]
      all_detectability[[length(all_detectability) + 1]] <- detect_stats
      
      message("  Detectability by relationship:")
      for (i in seq_len(nrow(detect_stats))) {
        message("    ", detect_stats$relationship[i], " : ",
                round(100 * detect_stats$detectability[i], 1), "% (",
                detect_stats$n_detectable[i], "/", detect_stats$n_pairs[i], ")")
      }
      
      # Store merged table (includes detectable column)
      merged_data_store[[paste(rep, rate_label, sep = "__")]] <- merged

      # ------------------ Evaluate each method on all tasks -------------------
      message("  Evaluating methods...")
      methods <- gsub("^score_", "", score_cols)

      for (method in methods) {
        sc <- paste0("score_", method)
        thresh <- global_class_thresholds[[sc]]

        message("    Method: ", method)

        # Task 1 — Overall relatedness: ≤25 gen vs >25 gen
        ev1 <- evaluate_binary_task(merged, sc, "truth_any_relatedness",
                                     threshold = thresh)

        # Task 2 — Transmission detection: ≤5 gen vs everything else
        ev2 <- evaluate_binary_task(merged, sc, "truth_transmission",
                                     threshold = thresh)

        # Task 3 — Transmission discrimination: ≤5 gen vs 6–25 gen
        #          (evaluated only on related pairs)
        merged_related <- merged[!is.na(is_transmission_vs_ancestry)]
        ev3 <- if (nrow(merged_related) > 0)
          evaluate_binary_task(merged_related, sc, "is_transmission_vs_ancestry",
                                threshold = thresh)
        else
          list(aupr = NA_real_, auroc = NA_real_, spearman = NA_real_,
               brier = NA_real_, n_pos = NA_integer_, n_neg = NA_integer_)

        # Task 4 — Transmission enrichment in top-K
        enrichment <- evaluate_transmission_enrichment(merged, sc)
        if (!is.null(enrichment) && nrow(enrichment) > 0) {
          enrichment[, `:=`(replicate  = rep,
                             rate       = rate_numeric,
                             rate_label = rate_label,
                             method     = method)]
          all_enrichment[[length(all_enrichment) + 1]] <- enrichment
        }
        
        # Store PR curve
        if (!is.null(ev1$pr_curve)) {
          prcurve <- rbind.data.frame(prcurve,
                                      data.table(replicate = rep, method = method, rate = rate_numeric, ev1$pr_curve)
                                      )
            names(prcurve)[4:6] <- c("recall", "precision", "threshold")
        }

        # --- Collate result row -----------------------------------------------
        row <- data.table(
          replicate  = rep,
          rate       = rate_numeric,
          rate_label = rate_label,
          method     = method,

          # Task 1
          aupr_overall  = ev1$aupr,
          auroc_overall = ev1$auroc,
          n_pos_overall = ev1$n_pos,

          # Task 2
          aupr_transmission  = ev2$aupr,
          auroc_transmission = ev2$auroc,
          n_pos_transmission = ev2$n_pos,

          # Task 3
          aupr_discrimination  = ev3$aupr,
          auroc_discrimination = ev3$auroc,

          # Correlation
          spearman_overall      = ev1$spearman,
          spearman_transmission = ev2$spearman,

          # Brier
          brier_overall      = ev1$brier,
          brier_transmission = ev2$brier,

          # Detectability (reference only)
          detectability_all = detect_stats[relationship != "unrelated",
                                           mean(detectability, na.rm = TRUE)],

          # Precision@K (transmission task)
          prec_at_10_trans = ev2$prec_at_10,
          prec_at_25_trans = ev2$prec_at_25,
          prec_at_50_trans = ev2$prec_at_50,

          # F1@K (transmission task)
          f1_at_10_trans   = ev2$f1_at_10,
          f1_at_25_trans   = ev2$f1_at_25,
          f1_at_50_trans   = ev2$f1_at_50,

          # Confusion matrix at calibrated threshold (task 1)
          thresh_precision = ev1$thresh_precision,
          thresh_recall    = ev1$thresh_recall,
          thresh_f1        = ev1$thresh_f1
        )

        all_results[[length(all_results) + 1]] <- row

        message("      Task 1 (overall)         AUPR = ", round(ev1$aupr, 3),
                "  AUROC = ", round(ev1$auroc, 3))
        message("      Task 2 (transmission)    AUPR = ", round(ev2$aupr, 3),
                "  AUROC = ", round(ev2$auroc, 3))
        message("      Task 3 (discrimination)  AUPR = ", round(ev3$aupr, 3),
                "  AUROC = ", round(ev3$auroc, 3))
      }
    }
    
    all_pr_curve[[rep]] <- prcurve
    message("\n  ✓ Completed replicate: ", rep)
  }

  # ============================================================================
  # STEP 2: AGGREGATE RESULTS
  # ============================================================================

  message("\n", strrep("=", 75))
  message("  STEP 2: AGGREGATING RESULTS ACROSS REPLICATES")
  message(strrep("=", 75))

  results_dt       <- rbindlist(all_results,       fill = TRUE)
  detectability_dt <- rbindlist(all_detectability, fill = TRUE)
  enrichment_dt    <- rbindlist(all_enrichment,    fill = TRUE)
  prcurve_dt       <- rbindlist(all_pr_curve,    fill = TRUE)

  message("  Total evaluations : ", nrow(results_dt))
  message("  Methods           : ", paste(unique(results_dt$method), collapse = ", "))
  message("  Rates             : ", paste(unique(results_dt$rate),   collapse = ", "))

  # Aggregate across replicates
  agg_results <- results_dt[, .(
    aupr_overall_mean     = mean(aupr_overall,         na.rm = TRUE),
    aupr_overall_se       = sd(aupr_overall,           na.rm = TRUE) / sqrt(.N),
    auroc_overall_mean    = mean(auroc_overall,         na.rm = TRUE),

    aupr_transmission_mean= mean(aupr_transmission,    na.rm = TRUE),
    aupr_transmission_se  = sd(aupr_transmission,      na.rm = TRUE) / sqrt(.N),
    auroc_transmission_mean = mean(auroc_transmission, na.rm = TRUE),

    aupr_discrimination_mean = mean(aupr_discrimination, na.rm = TRUE),
    aupr_discrimination_se   = sd(aupr_discrimination,   na.rm = TRUE) / sqrt(.N),
    auroc_discrimination_mean= mean(auroc_discrimination,na.rm = TRUE),

    spearman_overall_mean = mean(spearman_overall,     na.rm = TRUE),
    spearman_transm_mean  = mean(spearman_transmission,na.rm = TRUE),

    brier_overall_mean    = mean(brier_overall,        na.rm = TRUE),
    brier_transm_mean     = mean(brier_transmission,   na.rm = TRUE),

    n_replicates          = .N
  ), by = .(rate, method)]

  # Delta AUPR relative to baseline (from Script 1 design)
  for (task_col in c("aupr_overall_mean", "aupr_transmission_mean",
                      "aupr_discrimination_mean")) {
    base_vals <- agg_results[rate == baseline_rate, .(method, base = get(task_col))]
    agg_results <- merge(agg_results, base_vals, by = "method", all.x = TRUE)
    delta_col   <- gsub("_mean$", "_delta", task_col)
    agg_results[, (delta_col) := get(task_col) - base]
    agg_results[, base := NULL]
  }

  # Detectability summary by relationship and rate
  detect_summary <- detectability_dt[, .(
    detectability_mean = mean(detectability, na.rm = TRUE),
    detectability_se   = sd(detectability,   na.rm = TRUE) / sqrt(.N),
    mean_ibd_mean      = mean(mean_ibd_score,       na.rm = TRUE),
    n_pairs_mean       = mean(n_pairs,         na.rm = TRUE)
  ), by = .(rate, relationship)]

  # Enrichment summary
  enrichment_summary <- if (nrow(enrichment_dt) > 0)
    enrichment_dt[, .(
      enrichment_mean         = mean(enrichment,         na.rm = TRUE),
      enrichment_se           = sd(enrichment,           na.rm = TRUE) / sqrt(.N),
      prop_transmission_mean  = mean(prop_transmission,  na.rm = TRUE),
      baseline_prop_mean      = mean(baseline_proportion,na.rm = TRUE)
    ), by = .(rate, method, k)]
  else
    data.table()

  # IBD decay data (for Plot 7, from Script 1 design)
  decay_rows <- lapply(names(merged_data_store), function(key) {
    merged <- merged_data_store[[key]]
    parts  <- strsplit(key, "__")[[1]]
    tmp    <- unique(merged[!is.na(gen_distance), .(pair_key, gen_distance, ibd_prop)])
    tmp[, `:=`(replicate   = parts[1],
               rate        = rate_lookup[parts[2]],
               rate_label  = parts[2])]
    tmp
  })
  
  decay_dt <- rbindlist(decay_rows, fill = TRUE)

  # ============================================================================
  # STEP 3: SAVE TABLES
  # ============================================================================

  message("\n", strrep("=", 75))
  message("  STEP 3: SAVING RESULTS")
  message(strrep("=", 75))

  if (CONFIG$SAVE_TABLES) {
    tdir <- file.path(CONFIG$OUTDIR, "tables")

    fwrite(results_dt,       file.path(tdir, "results_all_replicates.csv"))
    fwrite(agg_results,      file.path(tdir, "results_aggregated.csv"))
    fwrite(detectability_dt, file.path(tdir, "detectability_by_relationship.csv"))
    fwrite(detect_summary,   file.path(tdir, "detectability_summary.csv"))

    write_xlsx(list(
      results_all       = as.data.frame(results_dt),
      results_agg       = as.data.frame(agg_results),
      detectability     = as.data.frame(detectability_dt),
      detect_summary    = as.data.frame(detect_summary),
      enrichment        = as.data.frame(enrichment_dt),
      enrichment_summary= as.data.frame(enrichment_summary)
    ), file.path(tdir, "all_results.xlsx"))

    saveRDS(
      list(
        results_dt         = results_dt,
        agg_results        = agg_results,
        detectability_dt   = detectability_dt,
        detect_summary     = detect_summary,
        enrichment_dt      = enrichment_dt,
        enrichment_summary = enrichment_summary,
        decay_dt           = decay_dt,
        merged_data_store  = merged_data_store,
        ibd_eps_global     = ibd_eps_global,
        global_class_thresholds = global_class_thresholds,
        config             = CONFIG
      ),
      file.path(tdir, "workspace.rds")
    )

    message("  ✓ All tables and workspace saved to: ", tdir)
  }

  # ============================================================================
  # STEP 4: VISUALIZATIONS
  # ============================================================================

  if (CONFIG$SAVE_PLOTS) {
    message("\n", strrep("=", 75))
    message("  STEP 4: GENERATING VISUALIZATIONS")
    message(strrep("=", 75))

    generate_plots(agg_results, detect_summary, enrichment_summary,
                   decay_dt, baseline_rate)
  }

  # ============================================================================
  # STEP 5: HYPOTHESIS TESTS
  # ============================================================================

  message("\n", strrep("=", 75))
  message("  STEP 5: HYPOTHESIS TESTS")
  message(strrep("=", 75))

  # H1: Detectability decreases with recombination
  message("\n[H1] IBD detection sensitivity vs recombination rate")
  det_trend <- detect_summary[relationship == "transmission"]
  for (m in unique(agg_results$method)) {
    cor_res <- cor.test(log10(det_trend$rate),
                        det_trend$detectability_mean,
                        method = "spearman")
    message("  Transmission detectability (all methods) rho = ",
            round(cor_res$estimate, 3), "  p = ",
            format.pval(cor_res$p.value, digits = 3))
    break   # detectability is method-independent; report once
  }

  # H2 & H3: Transmission discrimination and enrichment vs recombination
  if (nrow(enrichment_summary) > 0) {
    message("\n[H2/H3] Transmission discrimination and enrichment vs recombination rate")
    for (m in unique(enrichment_summary$method)) {
      ed <- enrichment_summary[method == m & k == 50]
      if (nrow(ed) > 1) {
        cor_disc <- tryCatch(
          cor.test(log10(agg_results[method == m, rate]),
                   agg_results[method == m, aupr_discrimination_mean],
                   method = "spearman"),
          error = function(e) NULL
        )
        cor_enr <- cor.test(log10(ed$rate), ed$enrichment_mean,
                            method = "spearman")

        message("\n  Method: ", m)
        if (!is.null(cor_disc))
          message("    H2 — Discrimination AUPR vs rate : rho = ",
                  round(cor_disc$estimate, 3), "  p = ",
                  format.pval(cor_disc$p.value, digits = 3))
        message("    H3 — Enrichment@50 vs rate        : rho = ",
                round(cor_enr$estimate, 3), "  p = ",
                format.pval(cor_enr$p.value, digits = 3))
        message("         Baseline enrichment  = ",
                round(ed[rate == min(rate), enrichment_mean], 2), "×")
        message("         Max-rate enrichment  = ",
                round(ed[rate == max(rate), enrichment_mean], 2), "×")
      }
    }
  }

  # ============================================================================
  # FINAL SUMMARY
  # ============================================================================

  message("\n", strrep("=", 75))
  message("  ANALYSIS COMPLETE")
  message(strrep("=", 75))
  message("  Total pair evaluations : ", nrow(results_dt))
  message("  Methods evaluated      : ", paste(unique(results_dt$method), collapse = ", "))
  message("  Recombination rates    : ", paste(sort(unique(results_dt$rate)), collapse = ", "))
  message("  Replicates             : ", length(replicates))
  message("  Output                 : ", CONFIG$OUTDIR)
  message(strrep("=", 75), "\n")

  invisible(list(
    results            = results_dt,
    agg_results        = agg_results,
    detectability      = detectability_dt,
    detect_summary     = detect_summary,
    enrichment         = enrichment_dt,
    enrichment_summary = enrichment_summary,
    decay_dt           = decay_dt,
    merged_data        = merged_data_store,
    ibd_eps            = ibd_eps_global,
    class_thresholds   = global_class_thresholds
  ))
}

# ============================================================================
# SECTION 9: VISUALIZATION
# ============================================================================

generate_plots <- function(agg_results, detect_summary, enrichment_summary,
                            decay_dt, baseline_rate) {

  fdir <- file.path(CONFIG$OUTDIR, "figures")

  custom_theme <- theme_bw(base_size = 14) +
    theme(
      legend.position      = "bottom",
      legend.title         = element_blank(),
      legend.text          = element_text(size = 12, color = "black", face = "bold"),
      plot.title           = element_text(size = 15, color = "black", face = "bold",   hjust = 0.5),
      plot.subtitle        = element_text(size = 10, color = "gray30", hjust = 0.5),
      axis.title           = element_text(size = 14, color = "black", face = "bold"),
      axis.text            = element_text(size = 12, color = "black", face = "bold"),
      axis.line            = element_line(linewidth = 1, colour = "black", lineend = "square"),
      axis.ticks           = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length    = unit(0.22, "cm"),
      strip.text           = element_text(size = 13, color = "black", face = "bold"),
      panel.grid.minor     = element_blank()
    )

  x_scale_log <- scale_x_log10(
    breaks = sort(unique(agg_results$rate)),
    labels = function(x) format(x, scientific = TRUE)
  )

  # --------------------------------------------------------------------------
  # Plot 1: AUPR — Overall vs Transmission vs Discrimination
  # --------------------------------------------------------------------------
  message("  [1/8] Overall vs Transmission vs Discrimination AUPR...")

  plot_data <- rbind(
    agg_results[, .(rate, method, AUPR = aupr_overall_mean,      SE = aupr_overall_se,      Task = "Overall Relatedness")],
    agg_results[, .(rate, method, AUPR = aupr_transmission_mean, SE = aupr_transmission_se, Task = "Transmission Detection")],
    agg_results[, .(rate, method, AUPR = aupr_discrimination_mean,SE= aupr_discrimination_se,Task= "Transmission Discrimination")]
  )

  p1 <- ggplot(plot_data, aes(x = rate, y = AUPR, color = method,
                               linetype = Task,
                               group = interaction(method, Task))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = AUPR - SE, ymax = AUPR + SE),
                  width = 0.05, linewidth = 0.7) +
    x_scale_log +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_linetype_manual(values = c("Overall Relatedness"       = "solid",
                                     "Transmission Detection"    = "dashed",
                                     "Transmission Discrimination"= "dotted")) +
    labs(title    = "Method Performance Across Tasks",
         subtitle = "Overall ≤25 gen | Transmission ≤5 gen | Discrimination: ≤5 vs 6–25 gen",
         x = "Recombination Rate (bp⁻¹, log scale)",
         y = "AUPR (mean ± SE)") +
    custom_theme +
    guides(color    = guide_legend(order = 1),
           linetype = guide_legend(order = 2))

  ggsave(file.path(fdir, "01_AUPR_all_tasks.png"), p1,
         width = 11, height = 7, dpi = 300)

  # --------------------------------------------------------------------------
  # Plot 2: Transmission Discrimination AUPR (Task 3 alone — key result)
  # --------------------------------------------------------------------------
  message("  [2/8] Transmission discrimination AUPR...")

  p2 <- ggplot(agg_results, aes(x = rate, y = aupr_discrimination_mean,
                                 color = method, group = method)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = aupr_discrimination_mean - aupr_discrimination_se,
                      ymax = aupr_discrimination_mean + aupr_discrimination_se),
                  width = 0.05, linewidth = 0.8) +
    x_scale_log +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title    = "Transmission Discrimination Among Related Pairs",
         subtitle = "Task: distinguish transmission (≤5 gen) from ancestry (6–25 gen)",
         x = "Recombination Rate (bp⁻¹, log scale)",
         y = "AUPR (mean ± SE)") +
    custom_theme

  ggsave(file.path(fdir, "02_transmission_discrimination_AUPR.png"), p2,
         width = 10, height = 6, dpi = 300)

  # --------------------------------------------------------------------------
  # Plot 3: Detectability by Relationship Type
  # --------------------------------------------------------------------------
  message("  [3/8] Detectability by relationship type...")

  p3 <- ggplot(detect_summary[relationship != "unrelated"],
               aes(x = rate, y = detectability_mean,
                   color = relationship, group = relationship)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = detectability_mean - detectability_se,
                      ymax = detectability_mean + detectability_se),
                  width = 0.05, linewidth = 0.8) +
    x_scale_log +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       labels = scales::percent_format()) +
    scale_color_viridis_d(option = "D", end = 0.85) +
    labs(title    = "IBD Detectability by Relationship Type",
         subtitle = paste0("Threshold calibrated at baseline: ",
                           CONFIG$ALPHA_DETECT * 100, "% FPR among unrelated pairs"),
         x = "Recombination Rate (bp⁻¹, log scale)",
         y = "Detectability (% with IBD ≥ threshold)") +
    custom_theme +
    theme(legend.position = "right")

  ggsave(file.path(fdir, "03_detectability_by_relationship.png"), p3,
         width = 10, height = 6, dpi = 300)

  # --------------------------------------------------------------------------
  # Plot 4: Delta AUPR Heatmap (vs baseline) — from Script 1
  # --------------------------------------------------------------------------
  message("  [4/8] Delta AUPR heatmap...")

  heat_dt <- melt(
    agg_results[, .(rate, method,
                    `Overall Relatedness`        = aupr_overall_delta,
                    `Transmission Detection`     = aupr_transmission_delta,
                    `Transmission Discrimination`= aupr_discrimination_delta)],
    id.vars      = c("rate", "method"),
    variable.name= "Task",
    value.name   = "delta_aupr"
  )
  heat_dt[, rate_label := format(rate, scientific = TRUE)]
  heat_dt$rate_label <- factor(heat_dt$rate_label, levels = c("1e-09", "1e-08", "1e-07", "1e-06"))

  p4 <- ggplot(heat_dt, aes(x = rate_label, y = method, fill = delta_aupr)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", delta_aupr)),
              color = ifelse(abs(heat_dt$delta_aupr) > 0.20, "white", "black"),
              size = 4, fontface = "bold") +
    facet_wrap(~ Task, ncol = 1) +
    scale_fill_gradient2(
      low      = "darkred",
      mid      = "gold",
      high     = "#1a9850",
      midpoint = 0,
      name     = "Δ AUPR"
    ) +
    labs(title    = "Change in AUPR Relative to Baseline Recombination Rate",
         subtitle = paste0("Baseline = ", format(baseline_rate, scientific = TRUE),
                           " bp⁻¹  |  Green = improved, Red = degraded"),
         x = "Recombination Rate (bp⁻¹)",
         y = "") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text        = element_text(size = 12, color = "black", face = "bold"),
      axis.title.x     = element_text(size = 13, color = "black", face = "bold",
                                      margin = margin(t = 8)),
      strip.text       = element_text(size = 12, color = "black", face = "bold"),
      legend.text      = element_text(size = 11),
      legend.title     = element_text(size = 12, face = "bold"),
      plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 10, color = "gray30", hjust = 0.5),
      panel.grid       = element_blank()
    )

  ggsave(file.path(fdir, "04_delta_AUPR_heatmap.png"), p4,
         width = 10, height = 10, dpi = 300)

  # --------------------------------------------------------------------------
  # Plot 5: Transmission Enrichment in Top-K
  # --------------------------------------------------------------------------
  if (nrow(enrichment_summary) > 0) {
    message("  [5/8] Transmission enrichment in top-K...")

    p5 <- ggplot(enrichment_summary[k %in% c(10, 25, 50)],
                 aes(x = rate, y = enrichment_mean,
                     color = method, linetype = factor(k),
                     group = interaction(method, k))) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = enrichment_mean - enrichment_se,
                        ymax = enrichment_mean + enrichment_se),
                    width = 0.05, linewidth = 0.7) +
      x_scale_log +
      labs(title    = "Transmission Enrichment in Top-Ranked Pairs",
           subtitle = "Enrichment = (observed % transmission) / (background %); dashed line = random",
           x = "Recombination Rate (bp⁻¹, log scale)",
           y = "Enrichment Factor") +
      custom_theme +
      guides(color    = guide_legend(title = "Method", order = 1),
             linetype = guide_legend(title = "Top K",  order = 2))

    ggsave(file.path(fdir, "05_transmission_enrichment_topK.png"), p5,
           width = 10, height = 6, dpi = 300)
  }

  # --------------------------------------------------------------------------
  # Plot 6: IBD Decay by Generation Distance — from Script 1
  # --------------------------------------------------------------------------
  message("  [6/8] IBD decay with genealogical distance...")

  if (!is.null(decay_dt) && nrow(decay_dt) > 0) {
    decay_dt[, gen_bin := floor(gen_distance)]

    sum_decay <- decay_dt[, .(
      ibd_med = median(ibd_prop, na.rm = TRUE),
      ibd_lo  = quantile(ibd_prop, 0.25, na.rm = TRUE),
      ibd_hi  = quantile(ibd_prop, 0.75, na.rm = TRUE)
    ), by = .(rate, gen_bin)]

    sum_decay[, rate_label := factor(format(rate, scientific = TRUE))]

    p6 <- ggplot(sum_decay, aes(x = gen_bin, y = ibd_med, color = rate_label)) +
      geom_ribbon(aes(ymin = ibd_lo, ymax = ibd_hi, fill = rate_label),
                  alpha = 0.2, color = NA) +
      geom_line(linewidth = 1.2) +
      scale_color_viridis_d(name = "Recombination\nRate (bp⁻¹)") +
      scale_fill_viridis_d(name  = "Recombination\nRate (bp⁻¹)") +
      labs(title    = "IBD Decay with Genealogical Distance",
           subtitle = "Median ± IQR across all replicates",
           x = "Generations Apart",
           y = "IBD Proportion") +
      custom_theme

    ggsave(file.path(fdir, "06_IBD_decay_by_generation.png"), p6,
           width = 10, height = 7, dpi = 300)
  }

  # --------------------------------------------------------------------------
  # Plot 7: Performance Heatmap (AUPR across methods × rates)
  # --------------------------------------------------------------------------
  message("  [7/8] Performance summary heatmap...")

  heatmap_data <- melt(
    agg_results[, .(rate, method,
                    `Overall Relatedness`        = aupr_overall_mean,
                    `Transmission Detection`     = aupr_transmission_mean,
                    `Transmission Discrimination`= aupr_discrimination_mean)],
    id.vars      = c("rate", "method"),
    variable.name= "Task",
    value.name   = "aupr"
  )
  heatmap_data[, rate_label := format(rate, scientific = TRUE)]
  heatmap_data$rate_label <- factor(heatmap_data$rate_label, levels = c("1e-09", "1e-08", "1e-07", "1e-06"))

  p7 <- ggplot(heatmap_data, aes(x = rate_label, y = method, fill = aupr)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = round(aupr, 2)), size = 4, fontface = "bold") +
    facet_wrap(~ Task, ncol = 1) +
    scale_fill_viridis_c(option = "plasma", limits = c(0, 1), na.value = "gray90") +
    labs(title = "AUPR Summary Across Tasks and Recombination Rates",
         x     = "Recombination Rate (bp⁻¹)",
         y     = "Method",
         fill  = "AUPR") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  ggsave(file.path(fdir, "07_AUPR_summary_heatmap.png"), p7,
         width = 10, height = 9, dpi = 300)
  
  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  pr_rows <- list()
  for (key in names(merged_data_store)) {
    merged <- merged_data_store[[key]]
    # parse key to rep and rate
    parts <- strsplit(key, "__")[[1]]
    repname <- parts[1]; rrate <- parts[2]
    method_cols <- grep("^score_", names(merged), value = TRUE)
    
    for (sc in method_cols) {
      method_name <- sub("^score_", "", sc)
      # need both positives and negatives
      if (sum(merged$truth_transmission == 1, na.rm = TRUE) > 0 && sum(merged$truth_transmission == 0, na.rm = TRUE) > 0) {
        pos <- merged[truth_transmission == 1 & !is.na(get(sc)), get(sc)]
        neg <- merged[truth_transmission == 0 & !is.na(get(sc)), get(sc)]
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
      labs(title = "Area Under Precision–Recall Curve across recombination rates ", x = "Recall", y = "Precision") +
      theme(axis.text = element_text(size = 14, color = 'black', face = 'bold'), # angle = 45, hjust = 1,
            legend.text = element_text(size = 14, color = 'black', face = 'bold'),
            legend.title = element_blank(),
            plot.title = element_text(color = 'black', face = 'bold'),
            axis.title = element_text(size = 16, color = 'black', face = 'bold'),
            strip.text = element_text(size = 16, color = 'black', face = 'bold'))
    
    ggsave(file.path(fdir, "08_AUPR_across_recombination.png"),
           p_pr_ribbon, width = 12, height = 10, dpi = 300)
  }

  # --------------------------------------------------------------------------
  # Plot 8: Composite Main Figure
  # --------------------------------------------------------------------------
  message("  [8/8] Composite main figure...")

  composite <- (p1 / p2) | p3
  composite <- composite +
    plot_annotation(
      title  = "Limits of Genomic Identifiability Under Recombination",
      theme  = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))
    )

  ggsave(file.path(fdir, "00_composite_main_figure.png"), composite,
         width = 18, height = 12, dpi = 300)

  message("  ✓ All 8 plots saved to: ", fdir)
}

# ============================================================================
# RUN
# ============================================================================

if (!interactive() || !exists("analysis_results")) {
  analysis_results <- run_genomic_identifiability_analysis()
}
