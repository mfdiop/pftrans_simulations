#!/usr/bin/env Rscript

## Load command-line options from options.R
source(file.path("wrapper_codes", "options.R"))
## - Supports: IBD (segment and pairwise summaries), SFS, and Ne comparisons
## Usage (example):
## Rscript 01_slim_msprime_evaluation.r --true_ibd analysis/true_ibd.all.tsv --inf_ibd inferred.ibd --outdir analysis/eval

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage:\n  Rscript 01_slim_msprime_evaluation.r --true_ibd <true_ibd.tsv> --inf_ibd <inferred.ibd> --true_ne <true.ne> --est_ne <est.ne> --true_sfs <true.sfs.tsv> --est_sfs <est.sfs.tsv> --outdir <outdir>\n")
  quit(status = 1)
}

parse_kv <- function(argv) {
  kv <- list()
  i <- 1
  while (i <= length(argv)) {
    if (startsWith(argv[i], "--")) {
      key <- sub("^--", "", argv[i])
      val <- if (i + 1 <= length(argv) && !startsWith(argv[i+1], "--")) argv[i+1] else TRUE
      kv[[key]] <- val
      i <- i + if (identical(val, TRUE)) 1 else 2
    } else i <- i + 1
  }
  kv
}

opts <- parse_kv(args)
outdir <- if (!is.null(opts$outdir)) opts$outdir else "analysis/eval"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Load helpers (if present)
helper_candidates <- c(
  file.path("wrapper_codes", "01_useful_functions.R"), "01_useful_functions.R",
  file.path(dirname(getwd()), "wrapper_codes", "01_useful_functions.R")
)

found_helper <- NULL

for (hc in helper_candidates) if (file.exists(hc)) { found_helper <- hc; break }
if (!is.null(found_helper)) source(found_helper)

library(tidyverse)

## Utilities
safe_read_tsv <- function(fn) {
  if (is.null(fn) || !file.exists(fn)) return(NULL)
  readr::read_tsv(fn, show_col_types = FALSE)
}

## Convert interval list to pairwise summaries
intervals_to_pairwise <- function(df) {
  # df should contain Id1, Id2, Start, End (bp) and optionally Chrom
  req <- c("Id1","Id2","Start","End")
  if (!all(req %in% names(df))) stop("intervals_to_pairwise: missing columns: ", paste(setdiff(req, names(df)), collapse=","))
  df <- df %>% mutate(len = End - Start) %>% filter(len > 0)
  sum_by_pair <- df %>% group_by(Id1 = as.character(Id1), Id2 = as.character(Id2)) %>%
    summarise(total_ibd_bp = sum(len, na.rm=TRUE), n_segments = n(), mean_len = mean(len), .groups = "drop")
  # ensure symmetric pairs aggregated: create ordered pair key
  sum_by_pair <- sum_by_pair %>% mutate(key = if_else(Id1 <= Id2, paste(Id1,Id2,sep='|'), paste(Id2,Id1,sep='|')))
  sum_by_pair <- sum_by_pair %>% group_by(key) %>%
    summarise(Id1 = first(strsplit(key, "\\|")[[1]][1]), Id2 = first(strsplit(key, "\\|")[[1]][2]),
              total_ibd_bp = sum(total_ibd_bp), n_segments = sum(n_segments), mean_len = mean(mean_len), .groups = 'drop')
  sum_by_pair %>% select(Id1, Id2, total_ibd_bp, n_segments, mean_len)
}

## Pairwise comparison summary
compare_pairwise <- function(true_int, inf_int, out_prefix) {
  tpair <- intervals_to_pairwise(true_int)
  ipair <- intervals_to_pairwise(inf_int)

  merged <- full_join(tpair, ipair, by = c("Id1","Id2"), suffix = c("_true","_inf")) %>%
    replace_na(list(total_ibd_bp_true = 0, n_segments_true = 0, mean_len_true = 0,
                    total_ibd_bp_inf = 0, n_segments_inf = 0, mean_len_inf = 0))

  # metrics
  merged <- merged %>% mutate(
    diff_bp = total_ibd_bp_inf - total_ibd_bp_true,
    abs_err = abs(diff_bp),
    rel_err = if_else(total_ibd_bp_true == 0, NA_real_, diff_bp / total_ibd_bp_true)
  )

  # global scores
  rmse <- sqrt(mean((merged$total_ibd_bp_inf - merged$total_ibd_bp_true)^2))
  cor_val <- cor(merged$total_ibd_bp_true, merged$total_ibd_bp_inf, method = "pearson")

  readr::write_tsv(merged, paste0(out_prefix, ".pairwise_summary.tsv"))

  list(summary = merged, rmse = rmse, cor = cor_val)
}

## Segment-level overlap: compute TP/FP/FN using overlap fraction
segment_confusion <- function(true_int, inf_int, min_overlap_frac = 0.5) {
  # For each inferred segment, check if it overlaps any true segment for same pair (unordered)
  require(GenomicRanges)
  if (!all(c('Id1','Id2','Start','End') %in% names(true_int))) stop('true_int missing cols')
  if (!all(c('Id1','Id2','Start','End') %in% names(inf_int))) stop('inf_int missing cols')

  # normalize pair ordering
  norm_pair <- function(df) {
    df <- df %>% mutate(a = pmin(as.character(Id1), as.character(Id2)),
                        b = pmax(as.character(Id1), as.character(Id2)))
    df
  }
  tdf <- norm_pair(true_int)
  idf <- norm_pair(inf_int)

  # build GRanges per pair and test overlaps
  pairs <- union(paste(tdf$a,tdf$b,sep='|'), paste(idf$a,idf$b,sep='|'))
  results <- tibble(pair = pairs, TP = 0L, FP = 0L, FN = 0L)

  for (p in pairs) {
    parts <- strsplit(p, '\\|')[[1]]
    a <- parts[1]; b <- parts[2]
    tsub <- tdf %>% filter(a == !!a & b == !!b)
    isub <- idf %>% filter(a == !!a & b == !!b)

    if (nrow(isub) == 0 && nrow(tsub) == 0) next

    if (nrow(isub) == 0) {
      results <- results %>% mutate(FN = if_else(pair == p, as.integer(nrow(tsub)), FN))
      next
    }
    if (nrow(tsub) == 0) {
      results <- results %>% mutate(FP = if_else(pair == p, as.integer(nrow(isub)), FP))
      next
    }

  chrom_t <- if ('Chrom' %in% names(tsub)) as.character(tsub$Chrom) else rep('1', nrow(tsub))
  chrom_i <- if ('Chrom' %in% names(isub)) as.character(isub$Chrom) else rep('1', nrow(isub))
  gr_t <- GRanges(seqnames = chrom_t, ranges = IRanges(start = tsub$Start, end = tsub$End))
  gr_i <- GRanges(seqnames = chrom_i, ranges = IRanges(start = isub$Start, end = isub$End))
    ov <- findOverlaps(gr_i, gr_t)

    # for each inferred segment, decide if it's TP (any overlap fraction >= min_overlap_frac)
    tp_count <- 0L
    fp_count <- 0L
    for (ii in seq_along(gr_i)) {
      hits <- subjectHits(ov[queryHits(ov) == ii])
      is_tp <- FALSE
      if (length(hits) > 0) {
        for (h in hits) {
          ol <- intersect(gr_i[ii], gr_t[h])
          ol_len <- width(ol)
          frac <- ol_len / width(gr_i[ii])
          if (!is.na(frac) && frac >= min_overlap_frac) { is_tp <- TRUE; break }
        }
      }
      if (is_tp) tp_count <- tp_count + 1L else fp_count <- fp_count + 1L
    }
    fn_count <- nrow(tsub) - tp_count
    results <- results %>% mutate(TP = if_else(pair == p, tp_count, TP), FP = if_else(pair == p, fp_count, FP), FN = if_else(pair == p, fn_count, FN))
  }
  results <- results %>% mutate(precision = TP / (TP + FP), recall = TP / (TP + FN))
  results
}

## SFS comparison: compute chi-square or KS between counts (aligned by maf bin)
compare_sfs <- function(true_sfs, est_sfs, out_prefix) {
  if (is.null(true_sfs) || is.null(est_sfs)) return(NULL)
  t <- readr::read_tsv(true_sfs, show_col_types = FALSE)
  e <- readr::read_tsv(est_sfs, show_col_types = FALSE)
  # align by maf or k
  if ('maf' %in% names(t) && 'maf' %in% names(e)) {
    comb <- full_join(t, e, by = 'maf', suffix = c('_true','_est')) %>% replace_na(list(count_true = 0, count_est = 0))
  } else {
    comb <- full_join(t, e, by = 'k', suffix = c('_true','_est')) %>% replace_na(list(count_true = 0, count_est = 0))
  }
  chi <- chisq.test(comb$count_true + 1, comb$count_est + 1)
  readr::write_tsv(comb, paste0(out_prefix, ".sfs_comp.tsv"))
  list(table = comb, chi_sq = chi)
}

## Ne comparison: RMSE and Pearson correlation
compare_ne <- function(true_ne_fn, est_ne_fn, out_prefix) {
  t <- safe_read_tsv(true_ne_fn)
  e <- safe_read_tsv(est_ne_fn)
  if (is.null(t) || is.null(e)) return(NULL)
  if (!('GEN' %in% names(t))) t <- rename(t, GEN = 1)
  if (!('GEN' %in% names(e))) e <- rename(e, GEN = 1)
  comb <- full_join(t, e, by = 'GEN', suffix = c('_true','_est')) %>% replace_na(list(NE_true = NA, NE_est = NA))
  comb <- comb %>% mutate(NE_true = coalesce(NE_true, NE), NE_est = coalesce(NE_est, NE))
  comb <- comb %>% filter(!is.na(NE_true) & !is.na(NE_est))
  rmse <- sqrt(mean((comb$NE_est - comb$NE_true)^2))
  corv <- cor(comb$NE_true, comb$NE_est)
  readr::write_tsv(comb, paste0(out_prefix, ".ne_comp.tsv"))
  list(table = comb, rmse = rmse, cor = corv)
}

## ----- Transmission Chain and Outbreak Analysis Functions -----

#' Compute confusion matrix metrics for transmission pairs
#' @param true_pairs data.frame with columns source,target representing true transmission pairs
#' @param inf_pairs data.frame with columns source,target representing inferred pairs
#' @param directed logical; if TRUE treat pairs as directed (source->target matters)
#' @return List with TP,FP,TN,FN counts and derived metrics (sensitivity, specificity, etc)
compute_confusion_metrics <- function(true_pairs, inf_pairs, directed = FALSE) {
  # Standardize inputs
  tp <- as.data.frame(true_pairs)
  ip <- as.data.frame(inf_pairs)
  if (!all(c('source','target') %in% names(tp))) stop('true_pairs missing source/target columns')
  if (!all(c('source','target') %in% names(ip))) stop('inf_pairs missing source/target columns')
  
  # Create pair keys (handle directed/undirected)
  if (directed) {
    make_key <- function(s,t) paste(s,t,sep='->')
  } else {
    make_key <- function(s,t) {
      k1 <- pmin(as.character(s), as.character(t))
      k2 <- pmax(as.character(s), as.character(t))
      paste(k1,k2,sep='|')
    }
  }
  
  tp$key <- make_key(tp$source, tp$target)
  ip$key <- make_key(ip$source, ip$target)
  
  # Basic counts
  TP <- sum(ip$key %in% tp$key)
  FP <- sum(!ip$key %in% tp$key)
  FN <- sum(!tp$key %in% ip$key)
  
  # Derived metrics
  sensitivity <- TP / (TP + FN)
  precision <- TP / (TP + FP)
  recall <- sensitivity  # same as sensitivity
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  list(
    counts = list(TP=TP, FP=FP, FN=FN),
    metrics = list(
      sensitivity=sensitivity, 
      precision=precision,
      recall=recall,
      f1=f1
    )
  )
}

#' Assess accuracy of transmission direction prediction
#' @param true_pairs data.frame with source->target columns for true transmission chains
#' @param inf_pairs data.frame with source->target columns for inferred chains
#' @return Fraction of correctly oriented edges among true positive predictions
assess_direction_accuracy <- function(true_pairs, inf_pairs) {
  # First find matching pairs (ignoring direction)
  make_undir_key <- function(s,t) {
    k1 <- pmin(as.character(s), as.character(t))
    k2 <- pmax(as.character(s), as.character(t))
    paste(k1,k2,sep='|')
  }
  
  tp <- as.data.frame(true_pairs)
  ip <- as.data.frame(inf_pairs)
  
  tp$ukey <- make_undir_key(tp$source, tp$target)
  ip$ukey <- make_undir_key(ip$source, ip$target)
  tp$dkey <- paste(tp$source, tp$target, sep='->')
  ip$dkey <- paste(ip$source, ip$target, sep='->')
  
  # Find common undirected edges
  common_pairs <- intersect(unique(tp$ukey), unique(ip$ukey))
  if (length(common_pairs) == 0) return(NA_real_)
  
  # For each common edge, check if direction matches
  n_correct <- 0
  n_total <- 0
  for (p in common_pairs) {
    true_dirs <- tp$dkey[tp$ukey == p]
    inf_dirs <- ip$dkey[ip$ukey == p]
    n_correct <- n_correct + sum(inf_dirs %in% true_dirs)
    n_total <- n_total + length(inf_dirs)
  }
  
  n_correct / n_total
}

#' Evaluate outbreak cluster prediction accuracy
#' @param true_clusters data.frame with id,cluster columns for true outbreak grouping
#' @param inf_clusters data.frame with id,cluster columns for inferred outbreaks
#' @return Adjusted Rand Index (ARI) between clusterings
evaluate_clustering <- function(true_clusters, inf_clusters) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    warning("mclust package required for ARI calculation")
    return(NA_real_)
  }
  
  # Align cluster assignments
  ids <- intersect(true_clusters$id, inf_clusters$id)
  if (length(ids) == 0) return(NA_real_)
  
  true_vec <- true_clusters$cluster[match(ids, true_clusters$id)]
  inf_vec <- inf_clusters$cluster[match(ids, inf_clusters$id)]
  
  mclust::adjustedRandIndex(true_vec, inf_vec)
}

#' Calibrate genetic distance thresholds for link prediction
#' @param inf_links data.frame with source,target,distance columns
#' @param true_pairs data.frame with source,target columns
#' @param thresholds numeric vector of distance thresholds to evaluate
#' @return List with ROC curve data and optimal threshold
calibrate_genetic_distances <- function(inf_links, true_pairs, thresholds = NULL) {
  if (is.null(thresholds)) {
    thresholds <- sort(unique(inf_links$distance))
  }
  
  # Evaluate each threshold
  results <- data.frame(
    threshold = thresholds,
    TPR = NA_real_,
    FPR = NA_real_,
    precision = NA_real_,
    F1 = NA_real_
  )
  
  for (i in seq_along(thresholds)) {
    t <- thresholds[i]
    pred_pairs <- inf_links[inf_links$distance <= t,]
    metrics <- compute_confusion_metrics(true_pairs, pred_pairs)
    
    results$TPR[i] <- metrics$metrics$sensitivity
    results$precision[i] <- metrics$metrics$precision
    results$F1[i] <- metrics$metrics$f1
  }
  
  # Find optimal threshold (maximum F1)
  best_idx <- which.max(results$F1)
  
  list(
    calibration = results,
    optimal = list(
      threshold = thresholds[best_idx],
      f1 = results$F1[best_idx],
      precision = results$precision[best_idx],
      sensitivity = results$TPR[best_idx]
    )
  )
}

## MAIN: wire options
if (!is.null(opts$true_ibd) && !is.null(opts$inf_ibd)) {
  true_ibd <- readr::read_tsv(opts$true_ibd, show_col_types = FALSE)
  inf_ibd <- readr::read_tsv(opts$inf_ibd, show_col_types = FALSE)
  res <- compare_pairwise(true_ibd, inf_ibd, file.path(outdir, "ibd"))
  cat(sprintf("IBD pairwise: RMSE=%.3f, COR=%.3f\n", res$rmse, res$cor))
  # segment-level confusion if GenomicRanges available
  if (requireNamespace('GenomicRanges', quietly = TRUE)) {
    conf <- segment_confusion(true_ibd, inf_ibd, min_overlap_frac = as.numeric(opts$min_overlap_frac %||% 0.5))
    readr::write_tsv(conf, file.path(outdir, "ibd_segment_confusion.tsv"))
  }
}

if (!is.null(opts$true_sfs) && !is.null(opts$est_sfs)) {
  sfs_res <- compare_sfs(opts$true_sfs, opts$est_sfs, file.path(outdir, "sfs"))
  if (!is.null(sfs_res)) cat(sprintf("SFS chi-sq p=%.3g\n", sfs_res$chi_sq$p.value))
}

if (!is.null(opts$true_ne) && !is.null(opts$est_ne)) {
  ne_res <- compare_ne(opts$true_ne, opts$est_ne, file.path(outdir, "ne"))
  if (!is.null(ne_res)) cat(sprintf("NE RMSE=%.3f COR=%.3f\n", ne_res$rmse, ne_res$cor))
}

# Transmission chain analysis if input files provided
if (!is.null(opts$true_pairs) && !is.null(opts$inf_pairs)) {
  true_pairs <- readr::read_tsv(opts$true_pairs, show_col_types = FALSE)
  inf_pairs <- readr::read_tsv(opts$inf_pairs, show_col_types = FALSE)
  
  # Basic confusion metrics
  conf <- compute_confusion_metrics(true_pairs, inf_pairs, 
                                  directed = as.logical(opts$directed %||% FALSE))
  readr::write_tsv(
    tibble::tibble(
      metric = c("TP","FP","FN","sensitivity","precision","recall","f1"),
      value = c(unlist(conf$counts), unlist(conf$metrics))
    ),
    file.path(outdir, "transmission_metrics.tsv")
  )
  cat(sprintf("Transmission metrics:\n  Sensitivity=%.3f\n  Precision=%.3f\n  F1=%.3f\n",
              conf$metrics$sensitivity, conf$metrics$precision, conf$metrics$f1))
  
  # Direction accuracy if analyzing as directed
  if (as.logical(opts$directed %||% FALSE)) {
    dir_acc <- assess_direction_accuracy(true_pairs, inf_pairs)
    cat(sprintf("Direction accuracy: %.3f\n", dir_acc))
  }
}

# Outbreak cluster comparison if provided
if (!is.null(opts$true_clusters) && !is.null(opts$inf_clusters)) {
  true_cl <- readr::read_tsv(opts$true_clusters, show_col_types = FALSE)
  inf_cl <- readr::read_tsv(opts$inf_clusters, show_col_types = FALSE)
  ari <- evaluate_clustering(true_cl, inf_cl)
  cat(sprintf("Cluster ARI: %.3f\n", ari))
}

# Genetic distance calibration if provided
if (!is.null(opts$inf_links) && !is.null(opts$true_pairs)) {
  inf_links <- readr::read_tsv(opts$inf_links, show_col_types = FALSE)
  true_pairs <- readr::read_tsv(opts$true_pairs, show_col_types = FALSE)
  
  if ("distance" %in% names(inf_links)) {
    cal <- calibrate_genetic_distances(inf_links, true_pairs)
    readr::write_tsv(as.data.frame(cal$calibration), 
                     file.path(outdir, "distance_calibration.tsv"))
    cat(sprintf("Optimal distance threshold: %.3f (F1=%.3f)\n",
                cal$optimal$threshold, cal$optimal$f1))
  }
}

cat("Done. Results written to:", normalizePath(outdir), "\n")
