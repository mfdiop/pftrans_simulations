# ============================================================================
# Diagnostic Script: Validate Metric Directions and Data Quality
# ============================================================================
# 
# Run this BEFORE the main analysis to verify:
# 1. All metrics are properly oriented (higher = more related)
# 2. Phylogenetic distances are correctly converted to similarities
# 3. Data files are readable and contain expected columns
# 4. Ground truth relationships are distributed as expected
#
# Usage: source("diagnostic_validation.R")
# ============================================================================

library(tidyverse)
library(data.table)
library(conflicted)
library(ape)
library(PRROC)
library(pROC)

# ============================================================================
# CONFIGURATION (match your main script)
# ============================================================================

DIAGNOSTIC_CONFIG <- list(
  ROOT_DIR = "simulations/multiple_runs",
  INFERRED_SUBPATH = "inferred",
  PHYLO_SUBPATH = "phylo_results",
  TEST_REPLICATE = "rep1",  # Test on first replicate
  TEST_RATE = "1e09",        # Test on baseline rate
  GENOME_BP = 640851,
  GEN_TRANSMISSION = 5,
  CHR = "chr1"
)

# ============================================================================
# LOAD HELPER FUNCTIONS (simplified versions from main script)
# ============================================================================

canonical_pair <- function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_character_)
  if (a <= b) paste(a, b, sep = "--") else paste(b, a, sep = "--")
}

# ============================================================================
# DIAGNOSTIC FUNCTIONS
# ============================================================================

check_file_exists <- function(filepath, description) {
  if (file.exists(filepath)) {
    cat("  ✓", description, "exists\n")
    return(TRUE)
  } else {
    cat("  ✗", description, "NOT FOUND:", filepath, "\n")
    return(FALSE)
  }
}

# test_metric_direction <- function(scores, truth, metric_name) {
#   cat("\n[", metric_name, "]\n")
#   
#   # Split by true class
#   pos_scores <- scores[truth == 1]
#   neg_scores <- scores[truth == 0]
#   
#   # Summary statistics
#   cat("  Positive class (related):\n")
#   cat("    n =", length(pos_scores), "\n")
#   cat("    mean =", round(mean(pos_scores, na.rm = TRUE), 4), "\n")
#   cat("    median =", round(median(pos_scores, na.rm = TRUE), 4), "\n")
#   cat("    range = [", round(min(pos_scores, na.rm = TRUE), 4), ",",
#       round(max(pos_scores, na.rm = TRUE), 4), "]\n")
#   
#   cat("  Negative class (unrelated):\n")
#   cat("    n =", length(neg_scores), "\n")
#   cat("    mean =", round(mean(neg_scores, na.rm = TRUE), 4), "\n")
#   cat("    median =", round(median(neg_scores, na.rm = TRUE), 4), "\n")
#   cat("    range = [", round(min(neg_scores, na.rm = TRUE), 4), ",",
#       round(max(neg_scores, na.rm = TRUE), 4), "]\n")
#   
#   # Test direction
#   mean_diff <- mean(pos_scores, na.rm = TRUE) - mean(neg_scores, na.rm = TRUE)
#   cat("  Mean difference (pos - neg):", round(mean_diff, 4), "\n")
#   
#   # AUPR test
#   pr <- tryCatch(
#     pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE),
#     error = function(e) NULL
#   )
#   
#   roc <- tryCatch(
#     roc(truth, scores, quiet = TRUE),
#     error = function(e) NULL
#   )
#   
#   if (!is.null(pr)) {
#     cat("  AUCPR:", round(pr$auc.integral, 4), "\n")
#   }
#   
#   if (!is.null(roc)) {
#     cat("  AUROC:", round(auc(roc), 4), "\n")
#   }
#   
#   # Validation
#   issues <- c()
#   
#   if (mean_diff < 0) {
#     issues <- c(issues, "Negative mean difference (positives score LOWER than negatives)")
#   }
#   
#   if (!is.null(pr) && pr$auc.integral < 0.5) {
#     issues <- c(issues, "AUCPR < 0.5 (worse than random)")
#   }
#   
#   if (!is.null(roc) && auc(roc) < 0.5) {
#     issues <- c(issues, "AUROC < 0.5 (worse than random)")
#   }
#   
#   if (length(issues) > 0) {
#     cat("\n  ❌ VALIDATION FAILED:\n")
#     for (issue in issues) {
#       cat("     •", issue, "\n")
#     }
#     cat("  → This metric appears to be INVERTED (distance instead of similarity)\n")
#     return(FALSE)
#   } else {
#     cat("\n  ✅ VALIDATION PASSED: Metric correctly oriented\n")
#     return(TRUE)
#   }
# }

test_metric_direction <- function(scores, truth, metric_name) {
  cat("\n[", metric_name, "]\n")
  
  # Split by true class
  pos_scores <- scores[truth == 1]
  neg_scores <- scores[truth == 0]
  
  # Class balance
  n_pos <- length(pos_scores)
  n_neg <- length(neg_scores)
  n_total <- n_pos + n_neg
  prevalence <- n_pos / n_total
  
  cat("  Class distribution:\n")
  cat("    Positive (related): n =", n_pos, 
      "(", round(100 * prevalence, 2), "%)\n")
  cat("    Negative (unrelated): n =", n_neg,
      "(", round(100 * (1 - prevalence), 2), "%)\n")
  
  # Summary statistics
  cat("\n  Positive class (related):\n")
  cat("    mean =", round(mean(pos_scores, na.rm = TRUE), 4), "\n")
  cat("    median =", round(median(pos_scores, na.rm = TRUE), 4), "\n")
  cat("    range = [", round(min(pos_scores, na.rm = TRUE), 4), ",",
      round(max(pos_scores, na.rm = TRUE), 4), "]\n")
  
  cat("  Negative class (unrelated):\n")
  cat("    mean =", round(mean(neg_scores, na.rm = TRUE), 4), "\n")
  cat("    median =", round(median(neg_scores, na.rm = TRUE), 4), "\n")
  cat("    range = [", round(min(neg_scores, na.rm = TRUE), 4), ",",
      round(max(neg_scores, na.rm = TRUE), 4), "]\n")
  
  # Test direction (key metric for similarity vs distance)
  mean_diff <- mean(pos_scores, na.rm = TRUE) - mean(neg_scores, na.rm = TRUE)
  median_diff <- median(pos_scores, na.rm = TRUE) - median(neg_scores, na.rm = TRUE)
  
  cat("\n  Separation metrics:\n")
  cat("    Mean difference (pos - neg):", round(mean_diff, 4), "\n")
  cat("    Median difference (pos - neg):", round(median_diff, 4), "\n")
  
  # AUPR test (affected by class imbalance)
  pr <- tryCatch(
    pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE),
    error = function(e) NULL
  )
  
  # AUROC test (NOT affected by class imbalance)
  roc <- tryCatch(
    roc(truth, scores, quiet = TRUE),
    error = function(e) NULL
  )
  
  cat("\n  Performance metrics:\n")
  
  if (!is.null(pr)) {
    aupr <- pr$auc.integral
    # AUPR baseline is the prevalence (proportion of positives)
    # A random classifier would achieve AUPR ≈ prevalence
    aupr_baseline <- prevalence
    aupr_enrichment <- aupr / aupr_baseline
    
    cat("    AUPR:", round(aupr, 4), "\n")
    cat("    AUPR baseline (prevalence):", round(aupr_baseline, 4), "\n")
    cat("    AUPR enrichment (AUPR/baseline):", round(aupr_enrichment, 2), "x\n")
  }
  
  if (!is.null(roc)) {
    auroc <- as.numeric(auc(roc))
    cat("    AUROC:", round(auroc, 4), "(threshold: 0.5)\n")
  }
  
  # Validation checks
  issues <- c()
  warnings <- c()
  
  # PRIMARY CHECK: Mean difference
  # For similarity metrics: positives should score HIGHER than negatives
  if (mean_diff < 0) {
    issues <- c(issues, 
                paste0("Negative mean difference (", round(mean_diff, 4), 
                       ") - positives score LOWER than negatives"))
  }
  
  # SECONDARY CHECK: AUROC (NOT affected by imbalance)
  # Should always be > 0.5 for a useful classifier with correct orientation
  if (!is.null(roc)) {
    auroc <- as.numeric(auc(roc))
    if (auroc < 0.5) {
      issues <- c(issues,
                  paste0("AUROC < 0.5 (", round(auroc, 4), 
                         ") - worse than random"))
    } else if (auroc < 0.6) {
      warnings <- c(warnings,
                    paste0("Low AUROC (", round(auroc, 4), 
                           ") - weak discrimination"))
    }
  }
  
  # TERTIARY CHECK: AUPR enrichment
  # Should be > 1 (better than random), but this is less critical with imbalance
  if (!is.null(pr)) {
    aupr <- pr$auc.integral
    aupr_baseline <- prevalence
    aupr_enrichment <- aupr / aupr_baseline
    
    if (aupr_enrichment < 1) {
      issues <- c(issues,
                  paste0("AUPR enrichment < 1 (", round(aupr_enrichment, 2),
                         "x) - worse than random"))
    } else if (aupr_enrichment < 2) {
      warnings <- c(warnings,
                    paste0("Low AUPR enrichment (", round(aupr_enrichment, 2),
                           "x) - weak performance"))
    }
  }
  
  # Report results
  cat("\n")
  
  if (length(issues) > 0) {
    cat("  ❌ VALIDATION FAILED:\n")
    for (issue in issues) {
      cat("     •", issue, "\n")
    }
    cat("\n  → This metric appears to be INVERTED\n")
    cat("     (likely a distance metric that needs conversion to similarity)\n")
    return(FALSE)
  } else if (length(warnings) > 0) {
    cat("  ⚠️  VALIDATION PASSED WITH WARNINGS:\n")
    for (warning in warnings) {
      cat("     •", warning, "\n")
    }
    cat("\n  → Metric is correctly oriented but may have weak discrimination\n")
    cat("     (this could be expected for some methods)\n")
    return(TRUE)
  } else {
    cat("  ✅ VALIDATION PASSED: Metric correctly oriented\n")
    cat("     • Positives score higher than negatives ✓\n")
    cat("     • AUROC > 0.5 (good discrimination) ✓\n")
    if (!is.null(pr)) {
      aupr_enrichment <- pr$auc.integral / prevalence
      cat("     • AUPR enrichment:", round(aupr_enrichment, 2), "x ✓\n")
    }
    return(TRUE)
  }
}

# ============================================================================
# MAIN DIAGNOSTIC
# ============================================================================

run_diagnostic <- function() {
  
  cat("\n")
  cat(rep("=", 70), "\n")
  cat("  DIAGNOSTIC VALIDATION\n")
  cat(rep("=", 70), "\n\n")
  
  cat("[CONFIG]\n")
  cat("  Testing replicate:", DIAGNOSTIC_CONFIG$TEST_REPLICATE, "\n")
  cat("  Testing rate:", DIAGNOSTIC_CONFIG$TEST_RATE, "\n\n")
  
  # ============================================================================
  # STEP 1: Check file existence
  # ============================================================================
  
  cat(rep("-", 70), "\n")
  cat("STEP 1: File Existence Check\n")
  cat(rep("-", 70), "\n")
  
  inferred_root <- file.path(DIAGNOSTIC_CONFIG$ROOT_DIR, 
                            DIAGNOSTIC_CONFIG$INFERRED_SUBPATH)
  
  rep_dir <- file.path(inferred_root, DIAGNOSTIC_CONFIG$TEST_REPLICATE)
  
  repnum <- gsub("rep", "", DIAGNOSTIC_CONFIG$TEST_REPLICATE)
  
  data_dir <- file.path(rep_dir, paste0("run", repnum, "_rec", DIAGNOSTIC_CONFIG$TEST_RATE, "_chr1"))
  
  cat("\nData directory:", data_dir, "\n\n")
  
  # Check files
  files_ok <- TRUE
  files_ok <- check_file_exists(file.path(data_dir, "true_ibd_summary.tsv"), 
                               "Ground truth IBD") && files_ok
  files_ok <- check_file_exists(file.path(data_dir, "inferred_ibd_hmm.rds"),
                               "Inferred IBD (HMM)") && files_ok
  files_ok <- check_file_exists(file.path(data_dir, "ibs_matrix.rds"),
                               "IBS matrix") && files_ok
  
  # Check phylo
  phylo_dir <- file.path(DIAGNOSTIC_CONFIG$ROOT_DIR, 
                        DIAGNOSTIC_CONFIG$PHYLO_SUBPATH,
                        DIAGNOSTIC_CONFIG$TEST_REPLICATE)
  
  cat("\nPhylo directory:", phylo_dir, "\n\n")
  
  if (dir.exists(phylo_dir)) {
    phylo_files <- list.files(phylo_dir,
                              pattern = paste0("rec", DIAGNOSTIC_CONFIG$TEST_RATE, "_", DIAGNOSTIC_CONFIG$CHR),
                              full.names = TRUE)
    if (length(phylo_files) > 0) {
      cat("  ✓ Phylogenetic tree files found:", length(phylo_files), " for chromosome1", "\n")
      cat("    Files:", paste(head(basename(phylo_files), 3), collapse = ", "), "\n")
    } else {
      cat("  ✗ No phylogenetic tree files found\n")
      files_ok <- FALSE
    }
  } else {
    cat("  ✗ Phylo directory not found\n")
    files_ok <- FALSE
  }
  
  if (!files_ok) {
    cat("\n❌ File check FAILED. Fix paths before proceeding.\n")
    return(invisible(NULL))
  }
  
  # ============================================================================
  # STEP 2: Load and examine ground truth
  # ============================================================================
  
  cat("\n", rep("-", 70), "\n")
  cat("STEP 2: Ground Truth Examination\n")
  cat(rep("-", 70), "\n\n")
  
  truth_file <- file.path(data_dir, "true_ibd_summary.tsv")
  truth <- fread(truth_file)
  
  cat("Ground truth dimensions:", nrow(truth), "rows ×", ncol(truth), "columns\n")
  cat("Columns:", paste(names(truth), collapse = ", "), "\n\n")
  
  # Standardize
  truth <- truth %>% 
    mutate(Id1 = paste0("tsk_", Id1), Id2 = paste0("tsk_", Id2)) %>%
    as.data.table()
  
  if ("generations" %in% names(truth)) {
    setnames(truth, "generations", "gen_distance")
  } else if ("t_generations" %in% names(truth)) {
    setnames(truth, "t_generations", "gen_distance")
  } else if ("min_tmrca" %in% names(truth)) {
    setnames(truth, "min_tmrca", "gen_distance")
  }
  
  if ("total_ibd_bp" %in% names(truth)) {
    truth[, ibd_prop := total_ibd_bp / DIAGNOSTIC_CONFIG$GENOME_BP]
  }
  
  # Relationship distribution
  cat("Relationship distribution (generations):\n")
  breaks <- c(0, DIAGNOSTIC_CONFIG$GEN_TRANSMISSION, 15, 25, Inf)
  labels <- c("Transmission (≤5)", "Recent (6-15)", "Distant (16-25)", "Unrelated (>25)")
  truth[, rel_cat := cut(gen_distance, breaks = breaks, labels = labels, right = TRUE)]
  
  rel_summary <- truth[, .N, by = rel_cat]
  print(rel_summary)
  
  cat("\nIBD proportion summary:\n")
  print(summary(truth$ibd_prop))
  
  # Define binary truth
  truth[, is_related := as.integer(gen_distance <= 25)]
  truth[, is_transmission := as.integer(gen_distance <= DIAGNOSTIC_CONFIG$GEN_TRANSMISSION)]
  
  cat("\nBinary classifications:\n")
  cat("  Related (≤25 gen):", sum(truth$is_related), "/", nrow(truth), "\n")
  cat("  Transmission (≤", DIAGNOSTIC_CONFIG$GEN_TRANSMISSION, " gen):", 
      sum(truth$is_transmission), "/", nrow(truth), "\n")
  
  # ============================================================================
  # STEP 3: Test IBD method
  # ============================================================================
  
  cat("\n", rep("-", 70), "\n")
  cat("STEP 3: IBD Method Validation\n")
  cat(rep("-", 70), "\n")
  
  ibd_file <- file.path(data_dir, "inferred_ibd_hmm.rds")
  ibd <- readRDS(ibd_file)
  ibd <- as.data.table(ibd)
  
  cat("\nIBD data dimensions:", nrow(ibd), "rows ×", ncol(ibd), "columns\n")
  cat("Columns:", paste(names(ibd), collapse = ", "), "\n")
  
  # Get score
  if ("total_ibd_bp" %in% names(ibd)) {
    ibd[, score := total_ibd_bp]
  } else if ("ibd_prop" %in% names(ibd)) {
    ibd[, score := ibd_prop]
  }  else if ("hmm" %in% names(ibd)) {
    ibd[, score := hmm]
  } else if ("score" %in% names(ibd)) {
    # already there
  } else {
    cat("\n⚠️  Warning: No obvious score column, using first numeric column\n")
    numeric_cols <- names(ibd)[sapply(ibd, is.numeric)]
    ibd[, score := get(numeric_cols[1])]
  }
  
  # Merge with truth
  if ("id1" %in% names(ibd)) {
    # already good
  } else if ("Id1" %in% names(ibd)) {
    setnames(ibd, c("Id1", "Id2"), c("id1", "id2"))
  }
  
  ibd[, pair_key := mapply(canonical_pair, id1, id2)]
  truth[, pair_key := mapply(canonical_pair, Id1, Id2)]
  
  merged <- merge(truth[, .(pair_key, is_related, is_transmission, ibd_prop, gen_distance)],
                 ibd[, .(pair_key, score_ibd = score)],
                 by = "pair_key")
  
  cat("\nMerged data:", nrow(merged), "pairs\n")
  
  # Test
  ibd_ok <- test_metric_direction(merged$score_ibd, merged$is_related, "IBD (HMM-based)")
  
  # ============================================================================
  # STEP 4: Test IBS method
  # ============================================================================
  
  cat("\n", rep("-", 70), "\n")
  cat("STEP 4: IBS Method Validation\n")
  cat(rep("-", 70), "\n")
  
  ibs_file <- file.path(data_dir, "ibs_matrix.rds")
  ibs_mat <- readRDS(ibs_file)
  ibs_mat <- as.matrix(ibs_mat)
  
  cat("\nIBS matrix dimensions:", nrow(ibs_mat), "×", ncol(ibs_mat), "\n")
  
  # Convert to long
  ids <- rownames(ibs_mat)
  if (is.null(ids)) {
    cat("⚠️  Warning: No row names, using indices\n")
    ids <- paste0("tsk_", 1:nrow(ibs_mat))
    rownames(ibs_mat) <- ids
    colnames(ibs_mat) <- ids
  }
  
  n <- length(ids)
  ibs_rows <- list()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      ibs_rows[[length(ibs_rows) + 1]] <- list(
        pair_key = canonical_pair(ids[i], ids[j]),
        score_ibs = as.numeric(ibs_mat[i, j])
      )
    }
  }
  ibs_long <- rbindlist(ibs_rows)
  
  merged_ibs <- merge(truth[, .(pair_key, is_related)],
                     ibs_long,
                     by = "pair_key")
  
  cat("\nMerged IBS data:", nrow(merged_ibs), "pairs\n")
  cat("IBS score range: [", round(min(merged_ibs$score_ibs, na.rm = TRUE), 4), ",",
      round(max(merged_ibs$score_ibs, na.rm = TRUE), 4), "]\n")
  
  # Test
  ibs_ok <- test_metric_direction(merged_ibs$score_ibs, merged_ibs$is_related, 
                                  "IBS (Genetic Similarity)")
  
  # ============================================================================
  # STEP 5: Test Phylogenetic method
  # ============================================================================
  
  cat("\n", rep("-", 70), "\n")
  cat("STEP 5: Phylogenetic Method Validation\n")
  cat(rep("-", 70), "\n")
  
  # # Find tree file
  # phylo_files <- list.files(phylo_dir, 
  #                          pattern = paste0("rec", DIAGNOSTIC_CONFIG$TEST_RATE),
  #                          full.names = TRUE)
  
  if (length(phylo_files) == 0) {
    cat("\n⚠️  No phylogenetic tree files found, skipping\n")
    phylo_ok <- NA
  } else {
    # Prefer modelfinder
    tree_file <- phylo_files[grepl(".treefile$", phylo_files)]
    if (length(tree_file) == 0) tree_file <- phylo_files[1]
    
    cat("\nReading tree:", basename(tree_file), "\n")
    
    tree <- read.tree(tree_file)
    cat("Tree has", length(tree$tip.label), "tips\n")
    
    # Compute patristic distances
    cat("Computing patristic distances...\n")
    dist_mat <- cophenetic.phylo(tree)
    
    cat("Distance matrix:", nrow(dist_mat), "×", ncol(dist_mat), "\n")
    cat("Distance range: [", round(min(dist_mat), 6), ",", 
        round(max(dist_mat), 6), "]\n")
    
    # Convert to long format (distances)
    ids <- rownames(dist_mat)
    n <- length(ids)
    phylo_rows <- list()
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        phylo_rows[[length(phylo_rows) + 1]] <- list(
          pair_key = canonical_pair(ids[i], ids[j]),
          distance = as.numeric(dist_mat[i, j])
        )
      }
    }
    
    phylo_long <- rbindlist(phylo_rows)
    
    # ✅ Convert distance to similarity
    max_dist <- max(phylo_long$distance, na.rm = TRUE)
    # phylo_long[, score_phylo := max_dist - distance]
    
    # Better alternatives (ordered from minimal to ideal)
    # ✅ Option 1: Normalize by tree diameter (minimal fix)
    
    phylo_long[, score_phylo := 1 - (distance / max_dist)]
    
    # # ✅ Option 2: Exponential decay kernel (recommended)
    # # Convert patristic distance to similarity using exponential kernel
    # alpha <- 1 / median(phylo_long$distance, na.rm = TRUE)
    # phylo_long[, score_phylo := exp(-alpha * distance)]
    
    
    cat("\n✅ Converting distance to similarity:\n")
    cat("  Original distance range: [", round(min(phylo_long$distance), 6), ",",
        round(max(phylo_long$distance), 6), "]\n")
    cat("  Converted similarity range: [", round(min(phylo_long$score_phylo), 6), ",",
        round(max(phylo_long$score_phylo), 6), "]\n")
    
    # Merge
    merged_phylo <- merge(truth[, .(pair_key, is_related)],
                         phylo_long[, .(pair_key, score_phylo)],
                         by = "pair_key")
    
    cat("\nMerged phylo data:", nrow(merged_phylo), "pairs\n")
    
    # Test CONVERTED metric
    phylo_ok <- test_metric_direction(merged_phylo$score_phylo, merged_phylo$is_related,
                                      "Phylo (Patristic → Similarity)")
    
    # Also show what would happen WITHOUT conversion
    cat("\n⚠️  COMPARISON: If we did NOT convert (using raw distance):\n")
    test_metric_direction(merged_phylo$score_phylo * (-1) + max_dist,  # Flip back to distance
                         merged_phylo$is_related,
                         "Phylo (RAW Distance - WRONG)")
  }
  
  # ============================================================================
  # FINAL SUMMARY
  # ============================================================================
  
  cat("\n", rep("=", 70), "\n")
  cat("DIAGNOSTIC SUMMARY\n")
  cat(rep("=", 70), "\n\n")
  
  results <- c(
    "Files" = ifelse(files_ok, "PASS", "FAIL"),
    "IBD" = ifelse(ibd_ok, "PASS", "FAIL"),
    "IBS" = ifelse(ibs_ok, "PASS", "FAIL"),
    "Phylo" = ifelse(is.na(phylo_ok), "SKIPPED", ifelse(phylo_ok, "PASS", "FAIL"))
  )
  
  for (i in seq_along(results)) {
    status_symbol <- ifelse(results[i] == "PASS", "✅", 
                          ifelse(results[i] == "FAIL", "❌", "⚠️ "))
    cat(status_symbol, names(results)[i], ":", results[i], "\n")
  }
  
  cat("\n")
  
  if (all(results[results != "SKIPPED"] == "PASS")) {
    cat("🎉 ALL CHECKS PASSED! 🎉\n")
    cat("\nYou can proceed with the main analysis.\n")
    cat("All metrics are correctly oriented (higher = more related).\n")
  } else {
    cat("⚠️  SOME CHECKS FAILED ⚠️\n")
    cat("\nPlease review the failures above before running the main analysis.\n")
    cat("Common issues:\n")
    cat("  • IBS metric is a distance (not similarity) - needs inversion\n")
    cat("  • Phylogenetic distances not converted to similarity\n")
    cat("  • File paths incorrect\n")
  }
  
  cat("\n", rep("=", 70), "\n\n")
  
  return(invisible(results))
}

# ============================================================================
# RUN
# ============================================================================

diagnostic_results <- run_diagnostic()
