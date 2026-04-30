# recombination_evaluation_pipeline_single_run.R
# Pipeline adapted to: ROOT_DIR = "simulations/single_run"
# Structure:
#   simulations/single_run/
#     inferred/
#       replicate_1/
#         inferred_ibd_hmm.tsv
#         ibs_matrix.rds
#         true_ibd_summary.tsv
#       replicate_2/ ...
#     phylo_output/
#       replicate_1/ ...    # not used by default in this script

# ----------------------
# Packages
# ----------------------
required_pkgs <- c("data.table", "PRROC", "tidyverse", "cowplot")
new_pkgs <- setdiff(required_pkgs, installed.packages()[, "Package"])

if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

library(data.table)
library(PRROC)
library(tidyverse)
library(cowplot)

# ----------------------
# User settings (edit if needed)
# ----------------------
ROOT_DIR <- "simulations/single_run"
INFERRED_ROOT <- file.path(ROOT_DIR, "inferred")   # where replicate_*/ live
PHYLO_ROOT    <- file.path(ROOT_DIR, "phylo_output") # kept for future
REPLICA_PATTERN <- "replicate"     # match folder names like replicate_1
INFERRED_IBD_FILE <- "inferred_ibd_hmm.tsv"
IBS_MATRIX_FILE   <- "ibs_matrix.rds"
TRUE_IBD_SUMMARY  <- "true_ibd_summary.tsv"
TREE_FILE <- "run1_chr1_modelfinder.treefile"

# recombination rate (Morgans per bp), genome length (bp)
r <- 0.01 / 15000              # as provided
genome_bp <- 100 * 15000       # as provided

# hybrid labeling parameters
S_MULT <- 1                    # total_ibd threshold multiplier (S_G = S_MULT * L_G_bp)
# Option C parameters
K_top <- 10                    # top-K per sample to include in eligible universe
# G grid (emphasis on G = 25)
G_values <- c(3, 5, 10, 15, 20, 25)
# minimal segment length to consider (if using segment-level files)
MIN_SEG_LEN_BP <- 30000 # 2 cM = 2 × 15,000 bp = 30,000 bp | 0.01/r
# output prefix
OUT_PREFIX <- "eval_single_run"

# ----------------------
# Helpers
# ----------------------
canonical_pair <- function(a, b) {
  # returns "idA--idB" where idA <= idB (string lexicographic)
  if (is.na(a) || is.na(b)) return(NA_character_)
  if (a <= b) paste(a, b, sep = "--") else paste(b, a, sep = "--")
}

split_pair_key <- function(pair_key) {
  cs <- tstrsplit(pair_key, "--")
  list(id1 = cs[[1]], id2 = cs[[2]])
}

# ----------------------
# Load truth summary per replicate
# Expect TRUE_IBD_SUMMARY to be per-pair summary with columns:
#   id1, id2, total_ibd_bp, max_seg_bp, n_segments   (but function is flexible)
# ----------------------
load_true_ibd_summary <- function(rep_dir) {
  fp <- file.path(rep_dir, TRUE_IBD_SUMMARY)
  # fp <- file.path(rep_inferred_dir, paste0(REPLICA_PATTERN, 1), TRUE_IBD_SUMMARY)
  if (!file.exists(fp)) {
    warning("Missing true_ibd_summary in ", rep_dir)
    return(NULL)
  }
  dt <- fread(fp)
  # canonicalize names (accept different conventions)
  # find id columns
  id_cols <- intersect(c("id1","id2","Id1","Id2","sample1","sample2"), names(dt))
  if (length(id_cols) < 2) stop("true_ibd_summary must have id1 and id2 columns")
  
  # rename first two matches to id1,id2
  setnames(dt, old = id_cols[1:2], new = c("id1","id2"))
  
  dt <- dt %>% 
    mutate(id1 = paste("tsk", id1, sep = '_'), id2 = paste("tsk", id2, sep = '_')) %>% 
    as.data.table()
  
  # ensure numeric fields exist; if not, create defaults
  if (!("total_ibd_bp" %in% names(dt))) dt[, total_ibd_bp := NA_real_]
  if (!("max_segment_bp" %in% names(dt))) dt[, max_seg_bp := NA_real_]
  if (!("n_segments" %in% names(dt))) dt[, n_segments := NA_integer_]
  
  # canonicalize pair_key
  dt[, pair_key := mapply(canonical_pair, id1, id2)]
  
  # deduplicate and aggregate if necessary (in case summary contains multiple rows)
  dt_pair <- dt[, .(
    total_ibd_bp = sum(as.numeric(total_ibd_bp), na.rm = TRUE),
    max_seg_bp   = max(as.numeric(max_segment_bp), na.rm = TRUE),
    n_segments   = sum(as.integer(n_segments), na.rm = TRUE)
  ), by = .(pair_key)]
  
  # split ids back
  ids <- split_pair_key(dt_pair$pair_key)
  dt_pair[, id1 := ids$id1]; dt_pair[, id2 := ids$id2]
  setcolorder(dt_pair, c("pair_key","id1","id2","total_ibd_bp","max_seg_bp","n_segments"))
  return(dt_pair)
}

# ----------------------
# Load inferred_ibd_hmm.tsv (long format) and convert to canonical pair, ensure a 'score' column
# Expect columns: id1,id2 and some score (if no score, create score=total_ibd_pred or 1)
# ----------------------
load_inferred_ibd <- function(rep_dir) {
  fp <- file.path(rep_dir, INFERRED_IBD_FILE)
  if (!file.exists(fp)) return(NULL)
  dt <- fread(fp)
  
  # allow different column names
  id_cols <- intersect(c("id1","id2","Id1","Id2","sample1","sample2"), names(dt))
  if (length(id_cols) < 2) stop("inferred_ibd file must have id columns")
  setnames(dt, old = id_cols[1:2], new = c("id1","id2"))
  
  # choose score: prefer total_ibd_bp, total_ibd, score or create constant
  score_candidates <- intersect(c("total_ibd_bp","total_ibd","score","ibd", "hmm", "ibd_score"), names(dt))
  if (length(score_candidates) > 0) {
    score_col <- score_candidates[1]
    dt[, score := as.numeric(get(score_col))]
  } else {
    # If inferred is a count of segments, use n_segments if available
    if ("n_segments" %in% names(dt)) dt[, score := as.numeric(n_segments)] else dt[, score := 1]
  }
  
  dt[, pair_key := mapply(canonical_pair, id1, id2)]
  
  # canonicalize id1,id2 so id1 < id2
  ids <- split_pair_key(dt$pair_key)
  dt[, id1 := ids$id1]; dt[, id2 := ids$id2]
  out <- unique(dt[, .(method = "IBD", id1, id2, score, pair_key)])
  return(out)
}

# ----------------------
# Load ibs_matrix.rds -> long format
# Assume RDS contains a matrix or data.frame where first column/rownames are sample IDs
# The numeric matrix is symmetric; convert to long (id1<id2) and transform distance -> score
# (we will convert distance to score: higher == stronger evidence of relatedness; so score = -distance)
# ----------------------
load_ibs_matrix <- function(rep_dir) {
  fp <- file.path(rep_dir, IBS_MATRIX_FILE)
  if (!file.exists(fp)) return(NULL)
  obj <- readRDS(fp)
  
  # obj might be matrix or data.frame
  if (is.data.frame(obj)) mat <- as.matrix(obj) else mat <- obj
  
  # If matrix has rownames, use them; else if first column is id column, handle that
  if (is.null(rownames(mat))) {
    # maybe it was saved as data.frame with first column being ids
    if (is.matrix(mat)) stop("IBS matrix has no rownames; cannot infer ids")
  }
  
  ids <- rownames(mat)
  
  if (is.null(ids)) {
    # try if obj was data.frame with first column as id
    if (is.data.frame(obj) && ncol(obj) >= 2) {
      ids <- as.character(obj[[1]])
      mat <- as.matrix(obj[, -1])
      rownames(mat) <- ids; colnames(mat) <- ids
    } else {
      stop("Cannot determine ids for IBS matrix in ", fp)
    }
  }
  
  # Melt upper triangle
  long <- data.table()
  n <- length(ids)
  rows <- vector("list", n*(n-1)/2)
  idx <- 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      rows[[idx]] <- list(id1 = ids[i], id2 = ids[j], distance = as.numeric(mat[i,j]))
      idx <- idx + 1
    }
  }
  
  long <- rbindlist(rows)
  # convert to score: smaller distance -> larger score
  # We'll create score = -distance (so higher score == more related). If your IBS is similarity (higher==more related),
  # you can set score = distance directly.
  # long[, score := -distance]
  long[, score := distance]
  long[, method := "IBS"]
  setcolorder(long, c("method","id1","id2","score"))
  return(long[, .(method,id1,id2,score)])
}


load_tree <- function(rep_dir) {
  source("simulations/R/patristic_distances.R")
  
  fp <- file.path(rep_dir, TREE_FILE)
  if (!file.exists(fp)) return(NULL) # stop("❌ Tree file not found: ", treefile) 
  
  dt <- as.data.table(phylo_long(fp, method = "patristic"))
  cat("   ✅ Phylogenetic distances calculated:", nrow(dt), "pairs\n")

  # Convert distance to similarity
  # Assuming your dataframe is a data.table called dt
  # Columns: id1, id2, phylo, score, pair_key
  
  dt[, `:=`(
    
    # 1. Inverse similarity: 1 / (1 + d)
    sim_inverse = 1 / (1 + phylo),
    
    # 2. Min–max normalized similarity
    #    (max distance -> 0, min distance -> 1)
    sim_minmax = 1 - (phylo - min(phylo)) / (max(phylo) - min(phylo)),
    
    # 3. Exponential similarity: exp(-d)
    sim_exp = exp(-phylo)
    
  )]
  
  # allow different column names
  id_cols <- intersect(c("id1","id2","Id1","Id2","sample1","sample2"), names(dt))
  if (length(id_cols) < 2) stop("inferred_ibd file must have id columns")
  setnames(dt, old = id_cols[1:2], new = c("id1","id2"))
  
  # choose score: prefer total_ibd_bp, total_ibd, score or create constant
  score_candidates <- intersect(c("score","sim_exp", "distance", "dist"), names(dt))
  
  if (length(score_candidates) > 0) {
    score_col <- score_candidates[1]
    dt[, score := as.numeric(get(score_col))]
  } else {
    # If inferred is a count of segments, use n_segments if available
    if ("n_segments" %in% names(dt)) dt[, score := as.numeric(n_segments)] else dt[, score := 1]
  }
  
  dt[, pair_key := mapply(canonical_pair, id1, id2)]
  
  # canonicalize id1,id2 so id1 < id2
  ids <- split_pair_key(dt$pair_key)
  dt[, id1 := ids$id1]; dt[, id2 := ids$id2]
  out <- unique(dt[, .(method = "patristic", id1, id2, score, pair_key)])
  return(out)
}

# ----------------------
# Build eligible universe Option C:
# union of all true-IBD pairs (pair_key) and top-K predicted pairs per sample across methods
# inferred_methods_list is a list of data.tables (method,id1,id2,score)
# ----------------------
build_eligible_universe_optionC <- function(true_pairs_dt, inferred_methods_list, K = K_top) {
  # true_pairs_dt: data.table with pair_key,id1,id2
  eligible_keys <- character()
  
  if (!is.null(true_pairs_dt) && nrow(true_pairs_dt) > 0) eligible_keys <- unique(true_pairs_dt$pair_key)
  
  # collect top-K keys per method
  for (dt in inferred_methods_list) {
    if (is.null(dt) || nrow(dt) == 0) next
    
    # ensure pair_key
    dt[, pair_key := mapply(canonical_pair, id1, id2)]
    
    # top K per focal sample (id1 and id2)
    setorder(dt, id1, -score)
    top1 <- dt[, head(.SD, K), by = id1]
    setorder(dt, id2, -score)
    top2 <- dt[, head(.SD, K), by = id2]
    cand <- unique(c(top1$pair_key, top2$pair_key))
    eligible_keys <- unique(c(eligible_keys, cand))
  }
  
  # assemble table
  if (length(eligible_keys) == 0) return(data.table())
  ids <- split_pair_key(eligible_keys)
  out <- data.table(pair_key = eligible_keys, id1 = ids$id1, id2 = ids$id2)
  return(out)
}

# ----------------------
# Label positives by hybrid rule using pair_summ (true pair summaries)
# pair_summ must have pair_key, max_seg_bp, total_ibd_bp
# maps G -> Lbp using Lbp = (1/G)/r
# ----------------------
# label_pairs_by_hybrid <- function(pair_summ, G, r, S_mult = S_MULT) {
#   if (is.null(pair_summ) || nrow(pair_summ) == 0) return(data.table())
#   Lbp <- (1 / G) / r
#   Sbp <- S_mult * Lbp
#   
#   # if max_seg_bp/total_ibd_bp are NA, fall back to total_ibd_bp > 0
#   pair_summ[, is_positive := 0L]
#   has_max <- !is.na(pair_summ$max_seg_bp)
#   has_total <- !is.na(pair_summ$total_ibd_bp)
#   
#   # positive if either condition satisfied; if no numeric data, use total_ibd_bp>0 fallback
#   if (any(has_max, na.rm = TRUE) || any(has_total, na.rm = TRUE)) {
#     pair_summ[has_max & (max_seg_bp >= Lbp), is_positive := 1L]
#     pair_summ[has_total & (total_ibd_bp >= Sbp), is_positive := 1L]
#   } else {
#     # fallback: if no numeric fields, treat pairs with total_ibd_bp > 0 as positive (if present)
#     pair_summ[!is.na(total_ibd_bp) & (total_ibd_bp > 0), is_positive := 1L]
#   }
#   pair_summ[, Lbp := Lbp]; pair_summ[, Sbp := Sbp]
#   return(pair_summ)
# }

label_pairs_by_ibd <- function(pair_summ,
                               genome_bp,
                               max_ibd_threshold = 5e5,   # 0.5 Mb
                               total_ibd_fraction = 0.02) # 2% genome
{
  if (is.null(pair_summ) || nrow(pair_summ) == 0) 
    return(data.table())

  pair_summ[, is_positive := 0L]
  
  # compute fraction of genome that is IBD
  if ("total_ibd_bp" %in% names(pair_summ)) {
    pair_summ[, f_ibd := total_ibd_bp / genome_size]
  } else {
    pair_summ[, f_ibd := NA_real_]
  }
  
  # biologically defensible rules:
  # 1. max IBD > threshold → strong recent relatedness
  pair_summ[!is.na(max_seg_bp) & (max_seg_bp >= max_ibd_threshold),
            is_positive := 1L]
  
  # 2. fraction IBD > threshold → high relatedness
  pair_summ[!is.na(f_ibd) & (f_ibd >= total_ibd_fraction),
            is_positive := 1L]
  
  # optional fallback: no data means "unknown", not positive
  return(pair_summ[])
}


# ----------------------
# Compute metric for a single method given eligible universe and truth labels
# inferred_dt: data.table(method,id1,id2,score)
# eligible_dt: data.table(pair_key,id1,id2)
# truth_dt: data.table(pair_key,is_positive)
# ----------------------
compute_metrics_for_method <- function(inferred_dt, eligible_dt, truth_dt) {
  
  if (is.null(inferred_dt) || nrow(inferred_dt) == 0) {
    return(list(
      TP = 0,
      FP = 0,
      FN = sum(truth_dt$is_positive == 1),
      TN = sum(truth_dt$is_positive == 0),
      precision = NA,
      recall = 0,
      f1 = NA,
      aupr = NA,
      pr_curve = NULL
    ))
  }
  
  # ensure canonical pair_key
  inferred_dt[, pair_key := mapply(canonical_pair, id1, id2)]
  
  # restrict to eligible universe only
  inferred_elig <- inferred_dt[pair_key %in% eligible_dt$pair_key]
  
  # truth lookup table
  label_map <- setNames(truth_dt$is_positive, truth_dt$pair_key)
  
  # assign truth label to predicted links
  inferred_elig[, true_label := as.integer(label_map[pair_key])]
  
  # TRUE POSITIVE: predicted positive & truly positive
  TP <- sum(inferred_elig$true_label == 1, na.rm = TRUE)
  
  # FALSE POSITIVE: predicted positive & truly negative
  FP <- sum(inferred_elig$true_label == 0, na.rm = TRUE)
  
  # lists used for FN / TN
  true_pos_keys <- truth_dt[is_positive == 1, pair_key]
  predicted_pos_keys <- unique(inferred_elig$pair_key)
  
  # FALSE NEGATIVE: true positive but not predicted positive
  FN <- length(setdiff(true_pos_keys, predicted_pos_keys))
  
  # TRUE NEGATIVE: all negatives - FP
  total_negatives <- nrow(eligible_dt) - length(true_pos_keys)
  TN <- total_negatives - FP
  
  # metrics
  precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
  recall    <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  f1        <- if (is.na(precision) || is.na(recall) || 
                   (precision + recall) == 0) NA
  else 2 * precision * recall / (precision + recall)
  
  aupr <- NA
  pr_curve <- NULL
  
  # PR curve only makes sense if scores exist
  if (nrow(inferred_elig) > 0 && ("score" %in% names(inferred_elig))) {
    pos_scores <- inferred_elig[true_label == 1, score]
    neg_scores <- inferred_elig[true_label == 0, score]
    
    if (length(pos_scores) > 0 && length(neg_scores) > 0) {
      pr <- pr.curve(scores.class0 = pos_scores,
                     scores.class1 = neg_scores,
                     curve = TRUE)
      aupr <- pr$auc.integral
      pr_curve <- as.data.table(pr$curve)
    }
  }
  
  return(list(
    TP = TP, FP = FP, FN = FN, TN = TN,
    precision = precision, recall = recall, f1 = f1,
    aupr = aupr, pr_curve = pr_curve
  ))
}


# ----------------------
# Evaluate one replicate
# ----------------------
evaluate_one_replicate <- function(rep_inferred_dir, tree_dir, G_values = G_values, K = K_top, r = r, S_mult = S_MULT) {
  # rep_inferred_dir: path to inferred/replicate_i
  # Load true summary
  true_dt <- load_true_ibd_summary(rep_inferred_dir)
  # xx <- load_true_ibd_summary(file.path(rep_inferred_dir, "replicate1"))
  if (is.null(true_dt) || nrow(true_dt) == 0) {
    warning("No true IBD summary for ", rep_inferred_dir); return(NULL)
  }
  
  # Load inferred methods
  cat("   Load inferred IBD dataframe...\n")
  ibd_pred <- tryCatch(load_inferred_ibd(rep_inferred_dir), error = function(e) { warning(e); NULL })
  
  cat("   Load inferred IBS dataframe...\n")
  ibs_pred <- tryCatch(load_ibs_matrix(rep_inferred_dir), error = function(e) { warning(e); NULL })
  
  # Calculate phylogenetic distances
  cat("   Calculating phylogenetic distances...\n")
  dist_pred <- tryCatch(load_tree(tree_dir), error = function(e) { warning(e); NULL })
  
  
  # collect all inferred method tables into a list (if present)
  inferred_list <- list()
  if (!is.null(ibd_pred)) inferred_list[[ibd_pred$method[1]]] <- ibd_pred
  if (!is.null(ibs_pred)) inferred_list[[ibs_pred$method[1]]] <- ibs_pred
  if (!is.null(dist_pred)) inferred_list[[dist_pred$method[1]]] <- dist_pred
  
  # build eligible universe (Option C)
  eligible_dt <- build_eligible_universe_optionC(true_dt, inferred_list, K = K)
  
  # if eligible is empty, fallback to all true pairs
  if (is.null(eligible_dt) || nrow(eligible_dt) == 0) {
    eligible_dt <- true_dt[, .(pair_key, id1, id2)]
  }
  
  # prepare pair_summ (true_dt already has pair_key,id1,id2,total_ibd_bp,max_seg_bp,n_segments)
  pair_summ <- copy(true_dt)
  
  # ensure numeric columns present
  if (!("total_ibd_bp" %in% names(pair_summ))) pair_summ[, total_ibd_bp := NA_real_]
  if (!("max_seg_bp" %in% names(pair_summ))) pair_summ[, max_seg_bp := NA_real_]
  if (!("n_segments" %in% names(pair_summ))) pair_summ[, n_segments := NA_integer_]
  
  # For each G compute labels and metrics per method
  res_per_G <- list()
  for (G in G_values) {
    # labeled_pairs <- label_pairs_by_hybrid(pair_summ, G = G, r = r, S_mult = S_mult)
    labeled_pairs <- label_pairs_by_hybrid(pair_summ, genome_bp)
    
    # truth labels for eligible universe
    eligible_labels <- merge(eligible_dt, labeled_pairs[, .(pair_key, is_positive)], by = "pair_key", all.x = TRUE)
    eligible_labels[is.na(is_positive), is_positive := 0]
    
    # For each inferred method compute metrics
    method_metrics <- list()
    for (mname in names(inferred_list)) {
      inferred_dt <- inferred_list[[mname]]
      
      # ensure the inferred dt has columns method,id1,id2,score
      mm <- compute_metrics_for_method(inferred_dt, eligible_labels, eligible_labels)
      method_metrics[[mname]] <- mm
    }
    # if there are inferred files saved as raw CSV (other methods), also try to read them:
    # attempt to read any other CSV/TSV in the replicate folder (exclude the known files)
    other_files <- list.files(rep_inferred_dir, pattern = "\\.(csv|tsv)$", full.names = TRUE)
    other_files <- setdiff(other_files, c(file.path(rep_inferred_dir, INFERRED_IBD_FILE), file.path(rep_inferred_dir, TRUE_IBD_SUMMARY)))
    
    for (fp in other_files) {
      dt <- tryCatch(fread(fp), error = function(e) NULL)
      if (is.null(dt)) next
      # try to standardize: need id1,id2,score
      if (!all(c("id1","id2") %in% names(dt))) {
        # try other names
        id_cols <- intersect(c("sample1","sample2","Id1","Id2"), names(dt))
        if (length(id_cols) >= 2) setnames(dt, old = id_cols[1:2], new = c("id1","id2"))
      }
      if (!("score" %in% names(dt))) {
        # attempt to identify a plausible score column
        sc <- intersect(c("total_ibd_bp","total_ibd","max_seg_bp","distance","dist","ibs","ibd"), names(dt))
        if (length(sc) > 0) dt[, score := as.numeric(get(sc[1]))] else dt[, score := 1]
      }
      # canonicalize pair ordering
      dt[, pair_key := mapply(canonical_pair, id1, id2)]
      ids <- split_pair_key(dt$pair_key)
      dt[, id1 := ids$id1]; dt[, id2 := ids$id2]
      dt_out <- unique(dt[, .(method = tools::file_path_sans_ext(basename(fp)), id1, id2, score)])
      method_metrics[[dt_out$method[1]]] <- compute_metrics_for_method(dt_out, eligible_labels, eligible_labels)
    }
    res_per_G[[as.character(G)]] <- list(Lbp = (1/G)/r, Sbp = S_mult * ((1/G)/r), metrics = method_metrics, eligible = eligible_labels)
  }
  return(res_per_G)
}

# ----------------------
# Evaluate all replicates under INFERRED_ROOT
# ----------------------
evaluate_all_replicates <- function(inferred_root = INFERRED_ROOT, tree_dir = tree_dir,
                                    G_values = G_values, K = K_top, r = r, S_mult = S_MULT) 
  {
  # find replicate folders under inferred_root
  if (!dir.exists(inferred_root)) stop("Inferred root does not exist: ", inferred_root)
  rep_dirs <- list.dirs(inferred_root, recursive = FALSE, full.names = TRUE)
  # tree_dirs <- list.dirs(tree_dir, recursive = FALSE, full.names = TRUE)
  
  # extract replicate numbers
  rep_nums <- as.numeric(sub(".*replicate", "", basename(rep_dirs)))
  
  # order by numeric value
  rep_dirs <- rep_dirs[order(rep_nums)]

  # filter by pattern
  #rep_dirs <- rep_dirs[grepl(REPLICA_PATTERN, basename(rep_dirs))]
  all_results <- list()
  
  for (rep in rep_dirs) {
    cat("Processing:", rep, "\n")
    res <- tryCatch(evaluate_one_replicate(rep, file.path(tree_dir, basename(rep)), G_values = G_values, K = K, r = r, S_mult = S_mult),
                    error = function(e) { message("Error in replicate ", basename(rep), ": ", e$message); NULL })
    all_results[[basename(rep)]] <- res
  }
  
  return(all_results)
}

# ----------------------
# Aggregation utilities for plotting
# ----------------------
aggregate_results <- function(all_results) {
  rows <- list(); k <- 0
  for (rep in names(all_results)) {
    rep_res <- all_results[[rep]]
    if (is.null(rep_res)) next
    for (G in names(rep_res)) {
      mm <- rep_res[[G]]$metrics
      for (m in names(mm)) {
        s <- mm[[m]]
        k <- k + 1
        rows[[k]] <- data.table(replicate = rep, G = as.integer(G), method = m,
                                TP = s$TP, FP = s$FP, FN = s$FN, TN = s$TN,
                                precision = s$precision, recall = s$recall, f1 = s$f1, aupr = s$aupr)
      }
    }
  }
  if (length(rows) == 0) return(data.table())
  return(rbindlist(rows, fill = TRUE))
}

# PR ribbon from saved pr_curves in results
build_pr_ribbon <- function(all_results, recall_grid = seq(0,1,length.out = 200)) {
  interp_rows <- list(); idx <- 0
  for (rep in names(all_results)) {
    rep_res <- all_results[[rep]]
    if (is.null(rep_res)) next
    for (G in names(rep_res)) {
      for (m in names(rep_res[[G]]$metrics)) {
        prc <- rep_res[[G]]$metrics[[m]]$pr_curve
        if (is.null(prc) || nrow(prc) == 0) next
        dt <- as.data.table(prc)
        setnames(dt, old = colnames(dt), new = c("recall","precision","threshold"))
        dt[, replicate := rep]; dt[, method := m]; dt[, G := as.integer(G)]
        ord <- order(dt$recall); dt <- dt[ord]
        # interpolate to recall_grid
        uniq <- !duplicated(dt$recall)
        rec <- dt$recall[uniq]; prec <- dt$precision[uniq]
        interp_prec <- approx(rec, prec, xout = recall_grid, rule = 2)$y
        for (i in seq_along(recall_grid)) {
          idx <- idx + 1
          interp_rows[[idx]] <- data.table(method = m, G = as.integer(G), replicate = rep, recall = recall_grid[i], precision = interp_prec[i])
        }
      }
    }
  }
  if (length(interp_rows) == 0) return(NULL)
  interp_dt <- rbindlist(interp_rows)
  summary_dt <- interp_dt[, .(precision_mean = mean(precision, na.rm = TRUE),
                              precision_lo = quantile(precision, 0.025, na.rm = TRUE),
                              precision_hi = quantile(precision, 0.975, na.rm = TRUE)), by = .(method, G, recall)]
  return(summary_dt)
}

# ----------------------
# Plotting functions
# ----------------------
plot_metrics_boxplots <- function(agg_dt, out_file) {
  png(out_file, width = 2200, height = 1200, res = 150)
  metrics <- c("precision","recall","f1")
  plots <- lapply(metrics, function(m) {
    ggplot(agg_dt, aes(x = "method", y = m, fill = "method")) +
      geom_boxplot(outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
      facet_wrap(~ G, scales = "free_x") +
      theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = NULL, y = toupper(m))
  })
  print(plot_grid(plotlist = plots, ncol = length(plots), labels = "AUTO"))
  dev.off()
}

plot_pr_ribbon <- function(pr_summary_dt, out_file) {
  if (is.null(pr_summary_dt)) return(NULL)
  png(out_file, width = 1200, height = 900, res = 150)
  p <- ggplot(pr_summary_dt, aes(x = recall, y = precision_mean, color = method, fill = method)) +
    geom_ribbon(aes(ymin = precision_lo, ymax = precision_hi), alpha = 0.2, color = NA) +
    geom_line(size = 1) +
    facet_wrap(~ G, ncol = 2) +
    theme_minimal() + labs(x = "Recall", y = "Precision", title = "PR curve (mean ±95% CI across replicates)")
  print(p)
  dev.off()
}

plot_performance_vs_G <- function(agg_dt, out_file) {
  perf_summary <- agg_dt[, .(precision_mean = mean(precision, na.rm = TRUE), precision_lo = quantile(precision, 0.025, na.rm = TRUE), precision_hi = quantile(precision, 0.975, na.rm = TRUE),
                             recall_mean = mean(recall, na.rm = TRUE), recall_lo = quantile(recall, 0.025, na.rm = TRUE), recall_hi = quantile(recall, 0.975, na.rm = TRUE),
                             f1_mean = mean(f1, na.rm = TRUE), f1_lo = quantile(f1, 0.025, na.rm = TRUE), f1_hi = quantile(f1, 0.975, na.rm = TRUE)), by = .(method, G)]
  png(out_file, width = 1400, height = 900, res = 150)
  p1 <- ggplot(perf_summary, aes(x = G, y = precision_mean, color = method, group = method)) + geom_line() + geom_ribbon(aes(ymin = precision_lo, ymax = precision_hi, fill = method), alpha = 0.15) + theme_minimal() + labs(y = "Precision", x = "G generations")
  p2 <- ggplot(perf_summary, aes(x = G, y = recall_mean, color = method, group = method)) + geom_line() + geom_ribbon(aes(ymin = recall_lo, ymax = recall_hi, fill = method), alpha = 0.15) + theme_minimal() + labs(y = "Recall", x = "G generations")
  p3 <- ggplot(perf_summary, aes(x = G, y = f1_mean, color = method, group = method)) + geom_line() + geom_ribbon(aes(ymin = f1_lo, ymax = f1_hi, fill = method), alpha = 0.15) + theme_minimal() + labs(y = "F1", x = "G generations")
  print(plot_grid(p1, p2, p3, ncol = 1, labels = "AUTO"))
  dev.off()
}

# ----------------------
# Top-level run function
# ----------------------
run_pipeline <- function(root_inferred = INFERRED_ROOT, tree_dir = PHYLO_ROOT, out_prefix = OUT_PREFIX, 
                         G_values = G_values, K = K_top, r = r, S_mult = S_MULT) {
  
  all_res <- evaluate_all_replicates(inferred_root = root_inferred, tree_dir = tree_dir, 
                                     G_values = G_values, K = K, r = r, S_mult = S_mult)
  
  saveRDS(all_res, "simulations/single_run/evaluation_results.rds")
  
  # Aggregate results between replicates
  agg <- aggregate_results(all_res)
  
  # Save aggregated results
  fwrite(agg, file = file.path("simulations/single_run", paste0(out_prefix, "_aggregated_metrics.csv")))
  
  # plots
  plot_metrics_boxplots(agg, file.path("simulations/single_run", paste0(out_prefix, "_metrics_boxplots.png")))
  
  pr_summary <- build_pr_ribbon(all_res)
  
  if (!is.null(pr_summary)) 
    plot_pr_ribbon(pr_summary, paste0(out_prefix, "_pr_ribbon.png"))
  
  plot_performance_vs_G(agg, paste0(out_prefix, "_performance_vs_G.png"))
  return(list(aggregated = agg, pr_summary = pr_summary, raw = all_res))
}

# ----------------------
# Run example
# ----------------------
# results <- run_pipeline(root_inferred = INFERRED_ROOT, out_prefix = "eval_single_run", G_values = c(3, 5, 10, 15, 20, 25), K = 10, r = 0.01/15000, S_mult = 1)

# End of script
# 
# Overview (short)

# Build a labelled training set (empirical linked pairs or simulated data).
# 
# Engineer features (max_seg_bp, total_ibd_bp, f_ibd, SNP distance, time difference, same_village, etc.).
# 
# Train a probabilistic classifier (start with logistic regression).
# 
# Calibrate and validate (cross-validation, AUPR, calibration curve).
# 
# Use the model to predict probabilities for inferred pairs.
# 
# Score/evaluate across thresholds (PR curve & AUPR) and produce final binary calls at a chosen threshold (or keep probabilities).
# 
# Plug probabilistic output into your compute_metrics_for_method by treating predictions with p >= t as predicted positives (or use whole PR sweep to compute AUPR).
# 
# Key modeling considerations (biological sanity)
# 
# Train on data resembling your study system (simulate malaria life cycle where possible).
# 
# Include covariates that affect IBD: sampling time difference, geography, multiplicity of infection (MOI), within-host diversity.
# 
# Address class imbalance (far fewer true transmission links than negatives): use stratified CV, class weights, or oversampling (SMOTE).
# 
# Quantify uncertainty: bootstrap or use Bayesian models (brms/rstanarm) if you want posterior probabilities.
# 
# Evaluate with PR curves (AUPR is a better metric than AUC when classes are imbalanced).
# 
# Keep probabilities (not only binary labels) for downstream epidemiological inference.



# packages
library(data.table)
library(caret)
library(PRROC)     # pr.curve
library(pROC)      # optional for ROC
set.seed(42)

# 1. Prepare training data ----------------------------------------------------
# Assumption: truth_dt contains pair_key and is_positive (1/0)
# and pair_summ contains pair_key and features: max_seg_bp, total_ibd_bp, snp_dist, time_diff, same_village, ...

prepare_training_table <- function(pair_summ, truth_dt, genome_size = NULL) {
  dt <- copy(pair_summ)
  dt <- merge(dt, truth_dt[, .(pair_key, is_positive)], by = "pair_key", all.x = FALSE, all.y = TRUE)
  # compute derived features
  if (!is.null(genome_size) && "total_ibd_bp" %in% names(dt)) {
    dt[, f_ibd := total_ibd_bp / genome_size]
  } else {
    dt[, f_ibd := NA_real_]
  }
  # add other derived features if available: log-transform, indicator of NA, etc.
  dt[, log_max_seg := ifelse(is.na(max_seg_bp), NA, log1p(max_seg_bp))]
  dt[, log_total_ibd := ifelse(is.na(total_ibd_bp), NA, log1p(total_ibd_bp))]
  # remove pairs with no truth label
  dt <- dt[!is.na(is_positive)]
  return(dt)
}

# 2. Train logistic regression with CV ---------------------------------------
train_prob_model <- function(train_dt, formula = is_positive ~ log_max_seg + log_total_ibd + f_ibd + time_diff + same_village,
                             folds = 5, repeat_cv = 2) {
  # caret preprocessing pipeline: center/scale, impute median for NAs
  ctrl <- trainControl(method = "repeatedcv", number = folds, repeats = repeat_cv,
                       classProbs = TRUE, summaryFunction = twoClassSummary,
                       savePredictions = "final", verboseIter = FALSE)
  # caret wants factor outcome "pos"/"neg"
  train_dt[, outcome := ifelse(is_positive == 1, "pos", "neg")]
  train_dt[, outcome := factor(outcome, levels = c("pos", "neg"))]
  
  preProc <- c("medianImpute", "center", "scale")
  # fit logistic regression (glmnet could be substituted)
  model <- train(as.formula(formula), data = train_dt,
                 method = "glm", family = "binomial",
                 trControl = ctrl, metric = "ROC", preProcess = preProc)
  return(model)
}

# 3. Predict probabilities on new inferred pairs ------------------------------
predict_probabilities <- function(model, pair_summ, genome_size = NULL) {
  dt <- copy(pair_summ)
  if (!is.null(genome_size) && "total_ibd_bp" %in% names(dt)) dt[, f_ibd := total_ibd_bp / genome_size]
  dt[, log_max_seg := ifelse(is.na(max_seg_bp), NA, log1p(max_seg_bp))]
  dt[, log_total_ibd := ifelse(is.na(total_ibd_bp), NA, log1p(total_ibd_bp))]
  probs <- predict(model, newdata = dt, type = "prob")
  dt[, prob := probs[,"pos"]]
  return(dt)
}

# 4. Evaluate model (PR curve & AUPR) ----------------------------------------
evaluate_prob_predictions <- function(pred_dt, truth_dt) {
  dt <- merge(pred_dt[, .(pair_key, prob)], truth_dt[, .(pair_key, is_positive)], by = "pair_key", all.x = FALSE)
  # PRROC wants scores for positives and negatives:
  pos_scores <- dt[is_positive == 1, prob]
  neg_scores <- dt[is_positive == 0, prob]
  if (length(pos_scores) > 0 && length(neg_scores) > 0) {
    pr <- pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = TRUE)
    return(list(pr = pr, AUPR = pr$auc.integral))
  } else {
    return(list(pr = NULL, AUPR = NA_real_))
  }
}

# 5. Convert probs to binary predictions at threshold t and compute confusion matrix
predict_binary_and_metrics <- function(pred_dt, truth_dt, threshold = 0.5, eligible_dt) {
  dt <- merge(pred_dt[, .(pair_key, prob)], eligible_dt[, .(pair_key)], by = "pair_key", all.y = TRUE)
  # pairs not present in pred_dt -> prob = 0 (or NA); here we set prob = 0 for "not predicted"
  dt[is.na(prob), prob := 0]
  dt[, pred_positive := as.integer(prob >= threshold)]
  # merge truth
  dt <- merge(dt, truth_dt[, .(pair_key, is_positive)], by = "pair_key", all.x = TRUE)
  dt[is.na(is_positive), is_positive := 0]  # if truth not present, treat as negative (careful)
  
  TP <- nrow(dt[pred_positive == 1 & is_positive == 1])
  FP <- nrow(dt[pred_positive == 1 & is_positive == 0])
  FN <- nrow(dt[pred_positive == 0 & is_positive == 1])
  TN <- nrow(dt[pred_positive == 0 & is_positive == 0])
  precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
  recall <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA else 2 * precision * recall / (precision + recall)
  return(list(TP = TP, FP = FP, FN = FN, TN = TN, precision = precision, recall = recall, f1 = f1))
}
