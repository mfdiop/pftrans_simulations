########################################################################
# evaluate_recombination_effects.R
# Purpose: evaluate how varying recombination rates affect metric performance
# Input layout (example):
# simulations/multiple_runs/inferred/rep1/run1rec1e-9_chr1/
#   - ibs_matrix.rds
#   - inferred_ibd_hmm.tsv
#   - true_ibd_summary.tsv
# simulations/multiple_runs/phylo_results/rep1/run1rec1e-9_chr1_modelfinder.treefile
#
# Adjust ROOT_DIR, REC_RATES and replicates as necessary.
########################################################################

library(data.table)
library(tidyverse)
library(PRROC)
library(pROC)
library(gridExtra)
library(cowplot)

# ---------------------------
# USER SETTINGS
# ---------------------------
ROOT_DIR <- "simulations/multiple_runs"        # change if needed
INFERRED_ROOT <- file.path(ROOT_DIR, "inferred")
PHYLO_ROOT    <- file.path(ROOT_DIR, "phylo_results") # contains treefiles per rate per rep
OUTDIR <- file.path(ROOT_DIR, "recombination_evaluation")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "tables"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "figures"), showWarnings = FALSE)

# recombination rates (character strings used in folder/file names)
REC_RATES <- c("1e-09", "1e-08", "1e-07", "1e-06")

# numeric versions for plotting (converted from strings)
REC_RATES_NUM <- as.numeric(sapply(REC_RATES, function(x) as.numeric(x)))

# Format recombination rate for filename matching
REC_RATES_LABEL = gsub("-", "", REC_RATES)  # Fix negative exponents


# replicates (use actual folder names under INFERRED_ROOT, e.g. rep1...rep5)
REPLICATES <- list.dirs(INFERRED_ROOT, recursive = FALSE, full.names = FALSE)

# Optionally filter by pattern like "rep" if other folders exist
REPLICATES <- REPLICATES[grepl("^rep", REPLICATES)]

GEN_THRESHOLD <- c(3, 5, 10, 15, 20, 25)   # your generational cutoff for positive label (24th cousins)
TOP_K <- c(1, 5, 10, 25, 50) # precision@K

# genome length for IBD proportion if needed (should be in truth summary ideally)
GENOME_BP <- 640851  # length of P.f chromosome1

# ---------------------------
# HELPER FUNCTIONS
# ---------------------------

canonical_pair <- function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_character_)
  if (a <= b) paste(a, b, sep = "--") else paste(b, a, sep = "--")
}

safe_fread <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  res <- tryCatch(fread(fp), error = function(e) {
    message("Failed to read ", fp, ": ", e$message)
    NULL
  })
  return(res)
}

# Load summarized tskit truth per replicate-rate folder
load_true_ibd_summary <- function(rep_dir, rate_str) {
  # expected file path: <rep_dir>/run<repnum>rec<RATE>_chr1/true_ibd_summary.tsv
  # or maybe directly under the folder for that run; we test both
  base_name <- basename(rep_dir)
  
  # try pattern search
  candidate_dir <- file.path(rep_dir, paste0("run", gsub("rep","", base_name), "_rec", rate_str, "_chr1"))
  
  # fallback: find any folder in rep_dir that contains rec{rate}
  if (!dir.exists(candidate_dir)) {
    cand <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
    cand_rate <- cand[grepl(paste0("rec", rate_str), cand)]
    if (length(cand_rate) > 0) candidate_dir <- cand_rate[1] else candidate_dir <- rep_dir
  }
  
  fp <- file.path(candidate_dir, "true_ibd_summary.tsv")
  dt <- safe_fread(fp)
  if (is.null(dt)) return(NULL)
  
  # Modify the sample IDs columns
  dt <- dt %>% 
    mutate(Id1 = paste0("tsk_", Id1),
           Id2 = paste0("tsk_", Id2)) %>% 
    as.data.table()
  
  # ensure pair_key present
  id_cols <- intersect(c("id1","id2","Id1","Id2","sample1","sample2"), names(dt))
  if (length(id_cols) >= 2) {
    setnames(dt, old = id_cols[1:2], new = c("id1","id2"))
  }
  
  # compute pair_key and ibd_prop if missing
  dt[, pair_key := mapply(canonical_pair, id1, id2)]
  if (!("total_ibd_bp" %in% names(dt))) dt[, total_ibd_bp := NA_real_]
  if (!("ibd_prop" %in% names(dt))) dt[, ibd_prop := total_ibd_bp / GENOME_BP]
  
  # if generations info present keep it
  if (!("gen_distance" %in% names(dt)) && ("min_tmrca" %in% names(dt))) dt[, gen_distance := min_tmrca]
  if (!("gen_distance" %in% names(dt))) dt[, gen_distance := NA_integer_]
  
  return(dt[, .(pair_key, id1, id2, total_ibd_bp, ibd_prop, max_seg_bp = ifelse("max_segment_bp" %in% names(dt), max_segment_bp, NA_real_), gen_distance)])
}

# Load inferred HMM-IBD file
load_inferred_ibd <- function(rep_dir, rate_str) {
  candidate_dir <- file.path(rep_dir, paste0("run", gsub("rep","", basename(rep_dir)), "_rec", rate_str, "_chr1"))
  
  if (!dir.exists(candidate_dir)) {
    cand <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
    cand_rate <- cand[grepl(paste0("rec", rate_str), cand)]
    if (length(cand_rate) > 0) candidate_dir <- cand_rate[1] else candidate_dir <- rep_dir
  }
  
  fp <- file.path(candidate_dir, "inferred_ibd_hmm.tsv")
  dt <- safe_fread(fp)
  
  if (is.null(dt)) return(NULL)
  id_cols <- intersect(c("id1","id2","Id1","Id2","sample1","sample2"), names(dt))
  
  if (length(id_cols) >= 2) setnames(dt, old = id_cols[1:2], new = c("id1","id2"))
  
  # choose score col
  score_candidates <- intersect(c("total_ibd_bp","total_ibd","score","ibd","ibd_prop","hmm"), names(dt))
  if (length(score_candidates) > 0) dt[, score := as.numeric(get(score_candidates[1]))] else {
    if ("n_segments" %in% names(dt)) dt[, score := as.numeric(n_segments)] else dt[, score := 1]
  }
  
  # Convert IBD proportion to distance
  dt[, score := 1-score]
  dt[, pair_key := mapply(canonical_pair, id1, id2)]
  return(dt[, .(pair_key, id1, id2, score)])
}

# Load IBS matrix rds -> long
load_ibs_matrix <- function(rep_dir, rate_str) {
  candidate_dir <- file.path(rep_dir, paste0("run", gsub("rep","", basename(rep_dir)), "_rec", rate_str, "_chr1"))
  
  if (!dir.exists(candidate_dir)) {
    cand <- list.dirs(rep_dir, recursive = FALSE, full.names = TRUE)
    cand_rate <- cand[grepl(paste0("rec", rate_str), cand)]
    if (length(cand_rate) > 0) candidate_dir <- cand_rate[1] else candidate_dir <- rep_dir
  }
  
  fp <- file.path(candidate_dir, "ibs_matrix.rds")
  
  if (!file.exists(fp)) return(NULL)
  obj <- tryCatch(readRDS(fp), error = function(e) { message("Failed readRDS: ", fp); NULL })
  
  if (is.null(obj)) return(NULL)
  
  # convert to pairwise long: assume matrix with rownames
  if (is.matrix(obj) || is.data.frame(obj)) {
    mat <- as.matrix(obj)
    ids <- rownames(mat)
    if (is.null(ids)) {
      # maybe data.frame with first col ids
      if (is.data.frame(obj)) {
        ids <- as.character(obj[[1]])
        mat <- as.matrix(obj[, -1])
        rownames(mat) <- ids; colnames(mat) <- ids
      } else return(NULL)
    }
    
    # melt upper triangle
    rows <- list(); k <- 0
    n <- length(ids)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        k <- k + 1
        rows[[k]] <- list(pair_key = canonical_pair(ids[i], ids[j]), id1 = ids[i], id2 = ids[j], distance = as.numeric(mat[i,j]))
      }
    }
    
    long <- rbindlist(rows)
    # convert to score: if distance is similarity bigger->more related, otherwise invert
    # attach column name 'score_ibs'
    long[, score := 1- distance] # Convert IBS to genetic distance
    return(long[, .(pair_key, id1, id2, score)])
  } else {
    return(NULL)
  }
}

# Load phylogenetic patristic distances from treefile (assumes you have a function phylo_long that returns pairwise distances)
load_phylo_patristic <- function(phylo_rep_dir, rate_str) {
  # rep_phylo_dir: directory for replicate under PHYLO_ROOT (e.g., PHYLO_ROOT/rep1)
  # treefile pattern: run<repnum>rec<RATE>_chr1_modelfinder.treefile
  # try to find file
  patt <- paste0("rec", rate_str, "_chr1")
  files <- list.files(phylo_rep_dir, pattern = patt, full.names = TRUE)
  
  # prefer the treefile as described
  treefile <- files[grepl("_modelfinder.treefile$", files)]
  if (length(treefile) == 0) {
    treefile <- files[1]
  }
  
  if (length(treefile) == 0) return(NULL)
  
  fp <- treefile[1]
  # user previously used phylo_long(fp, method = "patristic") in earlier pipeline
  # try to source that helper if exists
  
  if (!exists("phylo_long")) {
    # attempt to source a helper script in repository if exists
    helper <- file.path("simulations", "R", "patristic_distances.R")
    if (file.exists(helper)) source(helper) else {
      message("phylo_long helper not found; cannot load patristic distances for ", fp)
      return(NULL)
    }
  }
  
  dt <- as.data.table(phylo_long(fp, method = "patristic"))
  
  # ensure id1,id2 columns
  id_cols <- intersect(c("id1","id2","Id1","Id2","sample1","sample2"), names(dt))
  
  if (length(id_cols) >= 2) setnames(dt, old = id_cols[1:2], new = c("id1","id2"))
  # choose score column (smaller distance -> less related). Convert to similarity if desired
  
  # if ("phylo" %in% names(dt)) dt[, score := -phylo] else if ("distance" %in% names(dt)) dt[, score := -distance] else dt[, score := 1]
  # Let's keep the distance and convert IBD and IBS to distance using 1 - IBD/IBS
  if ("phylo" %in% names(dt)) dt[, score := phylo] else if ("distance" %in% names(dt)) dt[, score := distance] else dt[, score := 1]
  
  dt[, pair_key := mapply(canonical_pair, id1, id2)]
  return(dt[, .(pair_key, id1, id2, score)])
}


# Compute PR/AUROC/Spearman/Brier/precision@K for a given method on merged table
evaluate_method_on_table <- function(dt_merge, score_col = "score", truth_col = "is_positive", ks = TOP_K) {
  res <- list()
  
  # prepare scores and truth: remove NA scores
  dt <- copy(dt_merge)
  dt <- dt[!is.na(get(score_col))]
  
  # scale scores to [0,1] for calibration/Brier
  s <- dt[[score_col]]
  
  if (diff(range(s, na.rm = TRUE)) == 0) s_scaled <- rep(0.5, length(s)) else s_scaled <- (s - min(s, na.rm = TRUE)) / (max(s, na.rm = TRUE) - min(s, na.rm = TRUE))
  dt[, score_scaled := s_scaled]
  
  # AUPR
  pos_scores <- dt[get(truth_col) == 1, get(score_col)]
  neg_scores <- dt[get(truth_col) == 0, get(score_col)]
  
  if (length(pos_scores) > 0 && length(neg_scores) > 0) {
    pr <- pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE)
    res$aupr <- pr$auc.integral
  } else res$aupr <- NA_real_
  
  # AUROC
  if (length(unique(dt[[truth_col]]))>1) {
    roc_obj <- tryCatch(roc(dt[[truth_col]], dt[[score_col]]), error = function(e) NULL)
    res$auroc <- if (!is.null(roc_obj)) auc(roc_obj) else NA_real_
  } else res$auroc <- NA_real_
  
  # Spearman correlation with continuous truth (ibd_prop exists)
  if ("ibd_prop" %in% names(dt)) {
    res$spearman <- suppressWarnings(cor(dt[[score_col]], dt$ibd_prop, method = "spearman", use = "complete.obs"))
  } else res$spearman <- NA_real_
  
  # Brier score on scaled score vs truth
  if (!any(is.na(dt$score_scaled))) res$brier <- mean((dt$score_scaled - dt[[truth_col]])^2, na.rm = TRUE) else res$brier <- NA_real_
  
  # precision@K: rank by score desc, compute precision at top K
  ord <- dt[order(-get(score_col))]
  
  for (k in ks) {
    topk <- head(ord, k)
    if (nrow(topk) == 0) res[[paste0("prec_at_", k)]] <- NA_real_ else res[[paste0("prec_at_", k)]] <- mean(topk[[truth_col]] == 1, na.rm = TRUE)
  }
  return(res)
}

# -----------------------------------------------------
# MAIN: loop over replicates and recombination rates
# -----------------------------------------------------

all_rows <- list()
metrics_store <- list()  # nested list: metrics_store[[rate]][[replicate]][[method]] -> list of metrics
per_rate_merged_scores <- list()

for (rate_str in REC_RATES_LABEL) {
  metrics_store[[rate_str]] <- list()
  per_rate_merged_scores[[rate_str]] <- list()
  
  for (rep in REPLICATES) {
    cat("[rate=", rate_str, "] processing replicate ", rep, "...\n")
    
    rep_dir <- file.path(INFERRED_ROOT, rep)
    
    # load truth
    truth_dt <- load_true_ibd_summary(rep_dir, rate_str)
    
    if (is.null(truth_dt)) {
      message("Missing truth for rep ", rep, " rate ", rate_str); next
    }
    
    # construct binary truth using gen threshold if available, else use ibd_prop>0 threshold fallback
    if (!is.na(truth_dt$gen_distance[1])) {
      truth_dt[, is_positive := as.integer(gen_distance <= GEN)] # GEN_THRESHOLD
    } else {
      # fallback: ibd_prop > 0 and optionally > tiny epsilon (1e-8)
      truth_dt[, is_positive := as.integer(ibd_prop > 0)]
    }
    
    # load inferred method outputs
    ibd_dt  <- load_inferred_ibd(rep_dir, rate_str)   # HMM-IBD
    ibs_dt  <- load_ibs_matrix(rep_dir, rate_str)     # IBS
    
    phylo_dt <- NULL
    # phylo files are under PHYLO_ROOT/<rep> and contain runXrec<RATE>...treefile
    phylo_rep_dir <- file.path(PHYLO_ROOT, rep)
    if (dir.exists(phylo_rep_dir)) {
      phylo_dt <- load_phylo_patristic(phylo_rep_dir, rate_str)
    } else {
      message("No phylo dir for rep ", rep)
    }
    
    # standardize and merge: we want per pair scores from each method merged with truth
    # prepare method tables with column 'score' and pair_key
    method_tables <- list()
    if (!is.null(ibd_dt)) method_tables$IBD <- ibd_dt
    if (!is.null(ibs_dt)) method_tables$IBS <- ibs_dt
    if (!is.null(phylo_dt)) method_tables$Phylo <- phylo_dt
    
    # ensure all pair_keys in truth exist
    # merge each method with truth to evaluate
    metrics_store[[rate_str]][[rep]] <- list()
    merged_scores_rep <- data.table(pair_key = truth_dt$pair_key)
    merged_scores_rep <- merge(merged_scores_rep, truth_dt[, .(pair_key, is_positive, ibd_prop, gen_distance)], by = "pair_key", all.x = TRUE)
    
    for (m in names(method_tables)) {
      mt <- method_tables[[m]]
      # ensure unique pair_key & score
      mt <- unique(mt[, .(pair_key, score)])
      setnames(mt, "score", paste0("score_", m))
      merged_scores_rep <- merge(merged_scores_rep, mt, by = "pair_key", all.x = TRUE)
    }
    
    # store merged table for diagnostics / correlation
    per_rate_merged_scores[[rate_str]][[rep]] <- copy(merged_scores_rep)
    
    # evaluate each method
    for (m in names(method_tables)) {
      score_col <- paste0("score_", m)
      metrics_store[[rate_str]][[rep]][[m]] <- evaluate_method_on_table(merged_scores_rep, score_col = score_col, truth_col = "is_positive", ks = TOP_K)
      # also store AUPR/AUROC row in all_rows
      row <- data.table(rate = rate_str, rate_num = as.numeric(rate_str), replicate = rep, method = m,
                        aupr = metrics_store[[rate_str]][[rep]][[m]]$aupr,
                        auroc = metrics_store[[rate_str]][[rep]][[m]]$auroc,
                        spearman = metrics_store[[rate_str]][[rep]][[m]]$spearman,
                        brier = metrics_store[[rate_str]][[rep]][[m]]$brier)
      
      # add precision@k columns
      for (k in TOP_K) row[[paste0("prec_at_",k)]] <- metrics_store[[rate_str]][[rep]][[m]][[paste0("prec_at_",k)]]
      all_rows[[length(all_rows)+1]] <- row
    }
  } # end replicates loop
  
  # save per-rate merged scores for later inspection
  saveRDS(per_rate_merged_scores[[rate_str]], file.path(OUTDIR, "tables", paste0("merged_scores_rate_", rate_str, ".rds")))

} # end rates loop

# collate metrics table
metrics_dt <- rbindlist(all_rows, fill = TRUE)
fwrite(metrics_dt, file.path(OUTDIR, "tables", "metrics_by_rate_and_replicate.csv"))

# ---------------------------
# AGGREGATE ACROSS REPLICATES (summary per rate x method)
# ---------------------------
agg_by_rate_method <- metrics_dt[, .(
  aupr_mean = mean(aupr, na.rm = TRUE),
  aupr_sd   = sd(aupr, na.rm = TRUE),
  auroc_mean = mean(auroc, na.rm = TRUE),
  auroc_sd   = sd(auroc, na.rm = TRUE),
  spearman_mean = mean(spearman, na.rm = TRUE),
  spearman_sd = sd(spearman, na.rm = TRUE)
), by = .(rate, rate_num, method)]

fwrite(agg_by_rate_method, file.path(OUTDIR, "tables", "agg_by_rate_method.csv"))

# ---------------------------
# PLOTTING: metrics vs recombination rate
# ---------------------------
theme_set(theme_bw(base_size = 14))

# AUPR vs rate (log scale)
p_aupr <- ggplot(agg_by_rate_method, aes(x = rate_num, y = aupr_mean, color = method, group = method)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = aupr_mean - aupr_sd, ymax = aupr_mean + aupr_sd), width = 0.1) +
  scale_x_log10() + labs(x = "Recombination rate (bp^-1, log10)", y = "AUPR (mean ± sd)", title = "AUPR vs recombination rate") +
  theme(legend.position = "bottom")
ggsave(file.path(OUTDIR, "figures", "AUPR_vs_recrate.png"), p_aupr, width = 8, height = 5)

# AUROC vs rate
p_auroc <- ggplot(agg_by_rate_method, aes(x = rate_num, y = auroc_mean, color = method, group = method)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = auroc_mean - auroc_sd, ymax = auroc_mean + auroc_sd), width = 0.1) +
  scale_x_log10() + labs(x = "Recombination rate (bp^-1, log10)", y = "AUROC (mean ± sd)", title = "AUROC vs recombination rate") +
  theme(legend.position = "bottom")
ggsave(file.path(OUTDIR, "figures", "AUROC_vs_recrate.png"), p_auroc, width = 8, height = 5)

# Spearman vs rate
p_spear <- ggplot(agg_by_rate_method, aes(x = rate_num, y = spearman_mean, color = method, group = method)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = spearman_mean - spearman_sd, ymax = spearman_mean + spearman_sd), width = 0.1) +
  scale_x_log10() + labs(x = "Recombination rate (bp^-1, log10)", y = "Spearman correlation (mean ± sd)", title = "Spearman vs recombination rate") +
  theme(legend.position = "bottom")
ggsave(file.path(OUTDIR, "figures", "Spearman_vs_recrate.png"), p_spear, width = 8, height = 5)

# Precision@K vs rate (example K=10)
p_precK <- ggplot(metrics_dt, aes(x = as.numeric(rate), y = prec_at_10, color = method)) +
  stat_summary(fun = mean, geom = "line", aes(group = method)) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", aes(ymin = ..y.. - ..ymin.., ymax = ..y.. + ..ymax..), width = 0.05) +
  scale_x_log10() + labs(x = "Recombination rate", y = "Precision@10 (mean across replicates)", title = "Precision@10 vs recombination rate") +
  theme(legend.position = "bottom")
# save simple version
ggsave(file.path(OUTDIR, "figures", "Precision_at_10_vs_recrate.png"), p_precK, width = 8, height = 5)

# Heatmap of delta AUPR relative to smallest rate
# compute baseline (lowest recomb)
baseline <- agg_by_rate_method[rate == REC_RATES[1], .(method, base_aupr = aupr_mean)]
delta_dt <- merge(agg_by_rate_method, baseline, by = "method", all.x = TRUE)
delta_dt[, delta_aupr := aupr_mean - base_aupr]
# wide format
heat_dt <- dcast(delta_dt, method ~ rate_num, value.var = "delta_aupr")
mheat <- as.matrix(heat_dt[, -1, with = FALSE])
rownames(mheat) <- heat_dt$method
png(file.path(OUTDIR, "figures", "deltaAUPR_heatmap.png"), width = 800, height = 800)
par(mar = c(6,8,4,2))
heatmap(mheat, Rowv = NA, Colv = NA, scale = "none", col = viridis::viridis(100), margins = c(5,8), main = "Delta AUPR vs baseline recomb")
dev.off()

# ---------------------------
# PR ribbons for selected rates (overlay)
# ---------------------------
# We'll build PR curves for each replicate and method for a subset of rates (e.g., min, median, max)
selected_rates <- REC_RATES
# aggregate PR curves from saved merged_scores to compute pr.curve per replicate
pr_overlay_list <- list()
for (rate_str in selected_rates) {
  merged_list <- readRDS(file.path(OUTDIR, "tables", paste0("merged_scores_rate_", rate_str, ".rds")))
  # merged_list is a list of per-replicate dt
  for (rep in names(merged_list)) {
    dtm <- merged_list[[rep]]
    # for each method column score_<method>
    score_cols <- grep("^score_", names(dtm), value = TRUE)
    for (sc in score_cols) {
      method_name <- sub("^score_", "", sc)
      # if both positives and negatives exist
      if (sum(dtm$is_positive==1, na.rm = TRUE)>0 && sum(dtm$is_positive==0, na.rm = TRUE)>0) {
        pr_obj <- pr.curve(scores.class0 = dtm[get("is_positive")==1, get(sc)], scores.class1 = dtm[get("is_positive")==0, get(sc)], curve = TRUE)
        pr_dt <- as.data.table(pr_obj$curve)
        setnames(pr_dt, c("recall","precision","threshold"))
        pr_dt[, method := method_name]; pr_dt[, replicate := rep]; pr_dt[, rate := rate_str]
        pr_overlay_list[[length(pr_overlay_list)+1]] <- pr_dt
      }
    }
  }
}
pr_overlay_dt <- rbindlist(pr_overlay_list, fill = TRUE)
# compute mean precision per recall per method per rate
if (nrow(pr_overlay_dt)>0) {
  pr_summary <- pr_overlay_dt[, .(precision_mean = mean(precision, na.rm = TRUE), precision_lo = quantile(precision,0.025,na.rm=TRUE), precision_hi = quantile(precision,0.975,na.rm=TRUE)), by = .(rate, method, recall)]
  p_pr_overlay <- ggplot(pr_summary, aes(x = recall, y = precision_mean, color = method)) +
    geom_line(size = 1) +
    facet_wrap(~ rate, scales = "free_x") +
    labs(title = "PR curves across recombination rates", x = "Recall", y = "Precision") +
    theme_minimal()
  ggsave(file.path(OUTDIR, "figures", "PR_by_recrate_facet.png"), p_pr_overlay, width = 12, height = 8)
}

# ---------------------------
# Decay of IBDprop vs generation by recomb rate (if gen_distance present in truth)
# ---------------------------
decay_rows <- list()
for (rate_str in REC_RATES) {
  merged_list <- readRDS(file.path(OUTDIR, "tables", paste0("merged_scores_rate_", rate_str, ".rds")))
  for (rep in names(merged_list)) {
    dtm <- merged_list[[rep]]
    if (!("gen_distance" %in% names(dtm)) || all(is.na(dtm$gen_distance))) next
    tmp <- unique(dtm[, .(pair_key, gen_distance, ibd_prop)])
    tmp[, replicate := rep]; tmp[, rate := rate_str]
    decay_rows[[length(decay_rows)+1]] <- tmp
  }
}
if (length(decay_rows)>0) {
  decay_dt <- rbindlist(decay_rows)
  p_decay <- ggplot(decay_dt, aes(x = gen_distance, y = ibd_prop)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_smooth(method = "loess", se = TRUE, color = "steelblue") +
    facet_wrap(~ rate, scales = "free_y") +
    scale_x_continuous() +
    labs(title = "Decay of IBD proportion vs generation by recomb rate", x = "Generations apart", y = "IBD proportion") +
    theme_minimal()
  ggsave(file.path(OUTDIR, "figures", "IBD_decay_by_recrate.png"), p_decay, width = 12, height = 8)
}

# ---------------------------
# Save aggregated plots into a composite figure
# ---------------------------
png(file.path(OUTDIR, "figures", "composite_recrate_summary.png"), width = 1600, height = 1200, res = 150)
grid.arrange(p_aupr + ggtitle("AUPR vs recomb rate"), p_auroc + ggtitle("AUROC vs recomb rate"), p_spear + ggtitle("Spearman vs recomb rate"), ncol = 1)
dev.off()

# ---------------------------
# Save workspace for later inspection
# ---------------------------
saveRDS(list(metrics_dt = metrics_dt, agg_by_rate_method = agg_by_rate_method, per_rate_merged_scores = per_rate_merged_scores), file.path(OUTDIR, "tables", "recombination_evaluation_workspace.rds"))

message("Done. Outputs in: ", OUTDIR)
