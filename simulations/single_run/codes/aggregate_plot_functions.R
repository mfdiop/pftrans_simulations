
# ------------------------
# Utilities
# ------------------------

split_pair_key <- function(pair_key) {
  cs <- tstrsplit(pair_key, "--")
  list(id1 = cs[[1]], id2 = cs[[2]])
}

# safe get pr-curve table: expects PRROC::pr.curve$curve with columns (recall, precision, threshold)
pr_curve_to_dt <- function(pr) {
  if (is.null(pr)) return(NULL)
  # PRROC returns matrix with 3 cols: recall, precision, threshold
  dt <- as.data.table(pr$curve)
  setnames(dt, c("recall","precision","threshold"))
  return(dt)
}

# safe get roc table: if pROC::roc object -> produce fpr,tpr,threshold
roc_obj_to_dt <- function(roc_obj) {
  if (is.null(roc_obj)) return(NULL)
  # pROC::roc stores sensitivities (tpr) & specificities
  tpr <- roc_obj$sensitivities
  fpr <- 1 - roc_obj$specificities
  thr <- roc_obj$thresholds
  dt <- data.table(fpr = fpr, tpr = tpr, threshold = thr)
  # order by fpr ascending
  setorder(dt, fpr)
  return(dt)
}

# Interpolate curve (x,y) to a common grid x_grid using approx; rule=2 for extrapolation
interp_curve <- function(x, y, x_grid) {
  # remove duplicates in x
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 2) return(rep(NA_real_, length(x_grid)))
  ux <- x[ok]; uy <- y[ok]
  # ensure monotonic x for approx; approx expects increasing x
  ord <- order(ux)
  ux <- ux[ord]; uy <- uy[ord]
  approx(ux, uy, xout = x_grid, rule = 2)$y
}

# ------------------------
# Aggregate PR curves across replicates for a given G and method list
# Input: all_results (list of replicates), G_cut (integer), method_names (optional)
# ------------------------
aggregate_pr_ribbons <- function(all_results, G_cut = 25, recall_grid = seq(0,1,length.out = 200)) {
  rows <- list(); idx <- 0
  for (rep in names(all_results)) {
    rep_res <- all_results[[rep]]
    if (is.null(rep_res)) next
    strG <- as.character(G_cut)
    for (m in names(rep_res$curve_data)) {
      prc <- rep_res$curve_data[[m]]$pr$curve
      if (is.null(prc) || nrow(prc) == 0) next
      
      # prc may be stored as data.table or matrix with three columns (recall, precision, thr)
      dt <- as.data.table(prc)
      
      # name columns properly if needed
      if (ncol(dt) >= 3) setnames(dt, old = colnames(dt)[1:3], new = c("recall","precision","threshold"), skip_absent = TRUE)
      
      # ensure numeric
      dt[, recall := as.numeric(recall)]; dt[, precision := as.numeric(precision)]
      
      # unique recall
      uniq_idx <- !duplicated(dt$recall)
      rec <- dt$recall[uniq_idx]; prec <- dt$precision[uniq_idx]
      interp_prec <- interp_curve(rec, prec, recall_grid)
      
      for (i in seq_along(recall_grid)) {
        idx <- idx + 1
        rows[[idx]] <- data.table(replicate = rep, method = m, G = G_cut, recall = recall_grid[i], precision = interp_prec[i])
      }
    }
  }
  
  if (length(rows) == 0) return(NULL)
  interp_dt <- rbindlist(rows)
  summary_dt <- interp_dt[, .(
    precision_mean = mean(precision, na.rm = TRUE),
    precision_lo   = quantile(precision, 0.025, na.rm = TRUE),
    precision_hi   = quantile(precision, 0.975, na.rm = TRUE)
  ), by = .(method, G, recall)]
  return(summary_dt[])
}

# ------------------------
# Aggregate ROC curves across replicates
# ------------------------
aggregate_roc_ribbons <- function(all_results, G_cut = 25, fpr_grid = seq(0,1,length.out = 200)) {
  rows <- list(); idx <- 0
  for (rep in names(all_results)) {
    rep_res <- all_results[[rep]]
    if (is.null(rep_res)) next
    strG <- as.character(G_cut)
    for (m in names(rep_res$curve_data)) {
      # attempt to find stored roc, else compute from pr_curve if available by reconstructing scores not possible
      metr <- rep_res$curve_data[[m]]$roc
      roc_dt <- NULL
      if (!is.null(metr)) {
        # roc_dt <- as.data.table(metr$roc_curve)
        roc_dt <- as.data.table(roc_obj_to_dt(metr))
        # ensure columns
        if (ncol(roc_dt) >= 3) setnames(roc_dt, old = colnames(roc_dt)[1:3], new = c("fpr","tpr","threshold"), skip_absent = TRUE)
      } else if (!is.null(metr$pr_curve) && nrow(as.data.table(metr$pr_curve))>0) {
        # cannot reliably convert PR->ROC; skip
        roc_dt <- NULL
      } else {
        # fallback: skip if no stored ROC and no scores
        roc_dt <- NULL
      }
      if (is.null(roc_dt) || nrow(roc_dt)==0) next
      roc_dt[, fpr := as.numeric(fpr)]; roc_dt[, tpr := as.numeric(tpr)]
      uniq_idx <- !duplicated(roc_dt$fpr)
      fpr <- roc_dt$fpr[uniq_idx]; tpr <- roc_dt$tpr[uniq_idx]
      interp_tpr <- interp_curve(fpr, tpr, fpr_grid)
      for (i in seq_along(fpr_grid)) {
        idx <- idx + 1
        rows[[idx]] <- data.table(method = m, replicate = rep, G = G_cut, fpr = fpr_grid[i], tpr = interp_tpr[i])
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  interp_dt <- rbindlist(rows)
  summary_dt <- interp_dt[, .(
    tpr_mean = mean(tpr, na.rm = TRUE),
    tpr_lo   = quantile(tpr, 0.025, na.rm = TRUE),
    tpr_hi   = quantile(tpr, 0.975, na.rm = TRUE)
  ), by = .(method, G, fpr)]
  return(summary_dt[])
}

# ------------------------
# Collect per-replicate AUPR and AUROC for ranking (G_cut)
# ------------------------
collect_ranking_metrics <- function(all_results, G_cut = 25) {
  rows <- list(); k <- 0
  for (rep in names(all_results)) {
    rep_res <- all_results[[rep]]
    if (is.null(rep_res)) next
    for (m in names(rep_res$curve_data)) {
      k <- k + 1
      s <- rep_res$curve_data[[m]]
      # ensure aupr, auroc fields exist
      aupr <- if (!is.null(s$auc_pr)) s$auc_pr else NA_real_
      auroc <- if (!is.null(s$auc_roc)) s$auc_roc else NA_real_
      rows[[k]] <- data.table(replicate = rep, method = m, G = G_cut, aupr = aupr, auroc = auroc)
    }
  }
  if (length(rows) == 0) return(NULL)
  return(rbindlist(rows, fill = TRUE))
}

# ------------------------
# Correlation matrix: across replicates, compute Spearman correlation between method score and truth ibd_prop
# Expects each replicate at G has stored eligible table with columns pair_key, score per method? 
# We'll attempt to extract merged eligible table from results: rep_res[[G]]$eligible (as in your pipeline)

# compute_correlations_across_replicates <- function(all_results, G_cut = 25, methods = NULL) {
#         # build a list of data.tables with columns pair_key, ibd_prop, <method_scores...>
#         merged_list <- list()
#         for (rep in names(all_results)) {
#                 rep_res <- all_results[[rep]]
#                 if (is.null(rep_res)) next
#                 eligible_dt <- copy(rep_res$merged_data$conservative) # expected to contain pair_key, is_positive maybe scores
#                 if (is.null(eligible_dt) || nrow(eligible_dt)==0) next
#                 # Attempt: attach scores for each method from rep_res[[G]]$metrics? Usually metrics don't store raw scores.
#                 # So we expect your eligible_dt to have columns score_<method> or we cannot compute correlations.
#                 # To be flexible, look for columns that look like scores e.g., names containing "score" or method names
#                 dt <- as.data.table(copy(eligible_dt))
#                 
#                 # ensure ibd_prop present
#                 if (!("ibd" %in% names(dt)) && ("total_ibd_bp" %in% names(dt) && !is.na(dt$total_ibd_bp[1]))) {
#                         # need genome length to compute ibd_prop; if not available skip
#                         warning("eligible_dt has total_ibd_bp but not ibd_prop; skipping correlation for replicate ", rep)
#                         next
#                 }
#                 
#                 dt <- dt[, .(pair_key, ibd)]
#                 
#                 # now attempt to attach method scores: search rep_res[[strG]]$metrics entries for a stored 'scores_dt' if present
#                 for (m in names(rep_res[[strG]]$metrics)) {
#                         mobj <- rep_res[[strG]]$metrics[[m]]
#                         # when compute_metrics_for_method was called earlier, it didn't store raw scores. If you stored the 'inferred_dt' in another place, adapt here.
#                         # For now, try to find a stored column in eligible_dt named after method (score_<method> or <method>_score)
#                         cand_names <- c(paste0("score_", m), paste0(m, "_score"), m)
#                         found <- intersect(cand_names, names(rep_res[[strG]]$eligible))
#                         if (length(found)==1) {
#                                 dt[, (m) := rep_res[[strG]]$eligible[[found[1]]]]
#                         } else {
#                                 # try to extract scores from metrics object (if you attached 'scores_dt' earlier)
#                                 if (!is.null(mobj$scores_dt)) {
#                                         tmp <- copy(mobj$scores_dt[, .(pair_key, score)])
#                                         setnames(tmp, "score", m)
#                                         dt <- merge(dt, tmp, by = "pair_key", all.x = TRUE)
#                                 } else {
#                                         # cannot find scores; fill NA
#                                         dt[, (m) := NA_real_]
#                                 }
#                         }
#                 }
#                 merged_list[[rep]] <- dt
#         }
#         if (length(merged_list)==0) return(NULL)
#         # stack and compute Spearman correlations per method
#         combined <- rbindlist(merged_list, idcol = "replicate", fill = TRUE)
#         # pick methods columns (all except pair_key, ibd_prop, replicate)
#         method_cols <- setdiff(names(combined), c("pair_key","ibd_prop","replicate"))
#         # compute Spearman correlations
#         cors <- sapply(method_cols, function(mc) suppressWarnings(cor(combined[[mc]], combined$ibd_prop, method = "spearman", use = "complete.obs")))
#         cors_dt <- data.table(method = method_cols, spearman = cors)
#         return(list(combined = combined, cors_dt = cors_dt))
# }

#' Compute correlations between ground truth IBD and inferred metrics across replicates
#'
#' @param all_results A list of replicate results, each containing:
#'        - $merged_data$conservative: data.table with ground truth and inferred metrics
#'        - $merged_data$strict: alternative merging strategy
#'        - $merged_data$union: alternative merging strategy
#' @param G_cut Maximum generations cutoff (for filtering if needed)
#' @param methods Vector of method names to analyze (if NULL, uses all available)
#' @param merge_strategy Which merging strategy to use: "conservative", "strict", "union"
#' @param compute_ci Whether to compute confidence intervals via bootstrapping
#' @param n_boot Number of bootstrap samples for CI calculation
#' @param return_plot_data Whether to return data for plotting
#' @return List containing correlation results and optional plot data
compute_correlations_across_replicates <- function(all_results, 
                                                   G_cut = 25, 
                                                   methods = NULL,
                                                   merge_strategy = c("conservative", "strict", "union"),
                                                   compute_ci = FALSE,
                                                   n_boot = 100,
                                                   return_plot_data = TRUE) {
  
  # Validate inputs
  merge_strategy <- match.arg(merge_strategy)
  
  # Initialize storage
  combined_list <- list()
  method_correlations <- list()
  
  # Process each replicate
  for (rep_name in names(all_results)) {
    rep_data <- all_results[[rep_name]]
    
    # Skip if replicate has no data
    if (is.null(rep_data) || is.null(rep_data$merged_data)) next
    
    # Get the appropriate merged dataset
    merged_dt <- rep_data$merged_data[[merge_strategy]]
    
    if (is.null(merged_dt) || nrow(merged_dt) == 0) next
    
    # Ensure we have the required columns
    required_cols <- c("pair_key", "ibd")
    if (!all(required_cols %in% names(merged_dt))) {
      warning(sprintf("Replicate %s missing required columns. Skipping.", rep_name))
      next
    }
    
    # Create a clean data.table for this replicate
    dt <- data.table(
      replicate = rep_name,
      pair_key = merged_dt$pair_key,
      ibd_truth = merged_dt$ibd  # Ground truth IBD proportion
    )
    
    # Add inferred metrics (ensure they exist)
    inferred_metrics <- c("IBD", "IBS", "phylo")  # Your inferred metrics
    available_metrics <- intersect(inferred_metrics, names(merged_dt))
    
    if (length(available_metrics) == 0) {
      warning(sprintf("No inferred metrics found in replicate %s. Skipping.", rep_name))
      next
    }
    
    # Add each available metric
    for (metric in available_metrics) {
      dt[[metric]] <- merged_dt[[metric]]
    }
    
    # # Optional: Filter by generations if needed
    # if ("generations" %in% names(merged_dt) && !is.null(G_cut)) {
    #   dt <- dt[merged_dt$generations <= G_cut]
    # }
    
    # Optional: Add true_link indicator if needed
    if ("true_link" %in% names(merged_dt)) {
      dt$true_link <- merged_dt$true_link
    }
    
    # Store this replicate's data
    combined_list[[rep_name]] <- dt
  }
  
  # Check if we have any data
  if (length(combined_list) == 0) {
    message("No valid replicate data found.")
    return(NULL)
  }
  
  # Combine all replicates
  combined_dt <- rbindlist(combined_list, fill = TRUE, idcol = NULL)
  
  # Determine which methods to analyze
  if (is.null(methods)) {
    # Auto-detect methods (exclude metadata columns)
    metadata_cols <- c("replicate", "pair_key", "ibd_truth", "true_link", "generations")
    methods <- setdiff(names(combined_dt), metadata_cols)
  }
  
  # Initialize results storage
  results <- list(
    combined_data = combined_dt,
    merge_strategy = merge_strategy,
    n_replicates = length(combined_list),
    n_pairs = nrow(combined_dt),
    correlations = list()
  )
  
  # Compute correlations for each method
  for (method in methods) {
    if (!method %in% names(combined_dt)) {
      warning(sprintf("Method %s not found in combined data.", method))
      next
    }
    
    # Remove NA values for this method
    valid_data <- combined_dt[!is.na(get(method)) & !is.na(ibd_truth)]
    
    if (nrow(valid_data) < 2) {
      warning(sprintf("Insufficient non-NA data for method %s.", method))
      next
    }
    
    # Compute Spearman correlation
    cor_test <- suppressWarnings(
      cor.test(
        x = valid_data[[method]],
        y = valid_data$ibd_truth,
        method = "spearman",
        exact = FALSE,  # For large datasets
        use = "complete.obs"
      )
    )
    
    # Store basic correlation results
    method_results <- list(
      rho = cor_test$estimate,
      p_value = cor_test$p.value,
      n = nrow(valid_data),
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
    
    # Compute bootstrap confidence intervals if requested
    if (compute_ci && nrow(valid_data) >= 10) {
      set.seed(123)  # For reproducibility
      boot_results <- bootstrapped_correlation(
        x = valid_data[[method]],
        y = valid_data$ibd_truth,
        n_boot = n_boot
      )
      
      method_results$ci_lower <- boot_results$ci[1]
      method_results$ci_upper <- boot_results$ci[2]
      method_results$boot_samples <- boot_results$samples
    }
    
    # Compute correlation per replicate
    replicate_cors <- sapply(unique(combined_dt$replicate), function(rep) {
      rep_data <- combined_dt[replicate == rep]
      rep_data <- rep_data[!is.na(get(method)) & !is.na(ibd_truth)]
      if (nrow(rep_data) >= 2) {
        suppressWarnings(
          cor(rep_data[[method]], rep_data$ibd_truth, 
              method = "spearman", use = "complete.obs")
        )
      } else {
        NA_real_
      }
    })
    
    method_results$replicate_correlations <- replicate_cors
    method_results$mean_replicate_rho <- mean(replicate_cors, na.rm = TRUE)
    method_results$sd_replicate_rho <- sd(replicate_cors, na.rm = TRUE)
    
    # Store method results
    results$correlations[[method]] <- method_results
  }
  
  # Prepare summary table
  summary_dt <- data.table(
    method = character(),
    rho = numeric(),
    p_value = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    n_pairs = integer(),
    mean_replicate_rho = numeric(),
    sd_replicate_rho = numeric()
  )
  
  for (method in names(results$correlations)) {
    res <- results$correlations[[method]]
    summary_dt <- rbind(summary_dt, data.table(
      method = method,
      rho = res$rho,
      p_value = res$p_value,
      ci_lower = ifelse(is.null(res$ci_lower), NA_real_, res$ci_lower),
      ci_upper = ifelse(is.null(res$ci_upper), NA_real_, res$ci_upper),
      n_pairs = res$n,
      mean_replicate_rho = res$mean_replicate_rho,
      sd_replicate_rho = res$sd_replicate_rho
    ), fill = TRUE)
  }
  
  results$summary_table <- summary_dt[order(-abs(rho))]  # Sort by absolute correlation
  
  # Prepare plot data if requested
  if (return_plot_data) {
    results$plot_data <- prepare_correlation_plot_data(results, methods)
  }
  
  return(results)
}

#' Helper function for bootstrapped confidence intervals
bootstrapped_correlation <- function(x, y, n_boot = 1000, conf_level = 0.95) {
  n <- length(x)
  boot_samples <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    indices <- sample(1:n, n, replace = TRUE)
    boot_samples[i] <- suppressWarnings(
      cor(x[indices], y[indices], method = "spearman", use = "complete.obs")
    )
  }
  
  # Remove NA bootstrap samples
  boot_samples <- boot_samples[!is.na(boot_samples)]
  
  if (length(boot_samples) == 0) {
    return(list(samples = numeric(), ci = c(NA_real_, NA_real_)))
  }
  
  ci <- quantile(boot_samples, probs = c((1 - conf_level)/2, 1 - (1 - conf_level)/2))
  
  list(samples = boot_samples, ci = ci)
}

#' Prepare data for correlation visualization
prepare_correlation_plot_data <- function(results_list, methods) {
  plot_data <- list()
  
  cat("Prepare data for correlation for ", methods, "... \n")
  combined_dt <- results_list$combined_data
  
  # Prepare data for scatter plots
  scatter_data <- list()
  for (method in methods) {
    if (method %in% names(combined_dt)) {
      # Sample data for plotting (avoid too many points)
      plot_dt <- combined_dt[!is.na(get(method)) & !is.na(ibd_truth)]
      if (nrow(plot_dt) > 10000) {
        set.seed(123)
        plot_dt <- plot_dt[sample(.N, 10000)]
      }
      
      scatter_data[[method]] <- plot_dt[, .(ibd_truth, inferred = get(method), method)]
    }
  }
  
  cat("Prepare data for correlation distribution across replicates\n")
  plot_data$scatter <- rbindlist(scatter_data, use.names = TRUE)
  
  # Prepare data for correlation distribution across replicates
  rep_cor_data <- list()
  for (method in methods) {
    # This will be filled by the main function
    rep_cor_data[[method]] <- data.table(
      method = method,
      replicate = names(results_list$correlations[[method]]$replicate_correlations),
      rho = results_list$correlations[[method]]$replicate_correlations
    )
  }
  cat("Prepare data for correlation distribution across replicates\n")
  plot_data$replicate_cors <- rbindlist(rep_cor_data, use.names = TRUE)
  
  return(plot_data)
}

#' Create correlation visualization plots
plot_correlation_results <- function(correlation_results, 
                                     plot_type = c("scatter", "violin", "forest"),
                                     methods_to_plot = NULL,
                                     alpha = 0.3,
                                     point_size = 1) {
  
  if (is.null(correlation_results$plot_data)) {
    stop("Correlation results do not contain plot data. Run with return_plot_data = TRUE")
  }
  
  plot_type <- match.arg(plot_type)
  
  if (is.null(methods_to_plot)) {
    methods_to_plot <- names(correlation_results$correlations)
  }
  
  library(ggplot2)
  
  if (plot_type == "scatter") {
    # Scatter plots of IBD truth vs inferred metrics
    scatter_dt <- correlation_results$plot_data$scatter
    scatter_dt <- scatter_dt[method %in% methods_to_plot]
    
    p <- ggplot(scatter_dt, aes(x = ibd_truth, y = inferred)) +
      geom_point(alpha = alpha, size = point_size) +
      geom_smooth(method = "lm", formula = y ~ x, color = "red", se = TRUE) +
      facet_wrap(~ method, scales = "free_y") +
      labs(
        x = "Ground Truth IBD Proportion",
        y = "Inferred Metric Value",
        title = "Correlation between Ground Truth and Inferred Metrics",
        subtitle = sprintf("Merge Strategy: %s | %d Replicates", 
                           correlation_results$merge_strategy, 
                           correlation_results$n_replicates)
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", color = "black", size = 15), # hjust = 0, 
        strip.background = element_rect(fill = "gray90", color = "black", linetype = "solid", linewidth = 1.5),
        plot.title = element_text(size = 16, color = 'black', face = 'bold'),
        plot.subtitle = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black', face = 'bold'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square")
      )
    
  } else if (plot_type == "violin") {
    # Distribution of correlations across replicates
    rep_cor_dt <- correlation_results$plot_data$replicate_cors
    rep_cor_dt <- rep_cor_dt[method %in% methods_to_plot & !is.na(rho)]
    
    p <- ggplot(rep_cor_dt, aes(x = method, y = rho, fill = method)) +
      geom_violin(alpha = 0.6) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(
        x = "",
        y = "Spearman's ρ (per replicate)",
        title = "Distribution of Correlations Across Replicates",
        fill = "Method") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 16, color = 'black', face = 'bold'),
            axis.text.x = element_text(size = 12, color = 'black'), # angle = 45 , face = 'bold', hjust = 1
            axis.title = element_text(size = 14, color = 'black', face = 'bold'),
            axis.text.y = element_text(size = 12, color = 'black'),
            axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square")
      )
    
  } else if (plot_type == "forest") {
    # Forest plot of correlations with confidence intervals
    summary_dt <- correlation_results$summary_table
    summary_dt <- summary_dt[method %in% methods_to_plot]
    
    p <- ggplot(summary_dt, aes(x = rho, y = reorder(method, rho))) +
      geom_point(size = 3) +
      geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      labs(
        x = "Spearman's ρ (95% CI)",
        y = "Methods",
        title = "Correlation with Ground Truth IBD",
        subtitle = "Points show correlation, bars show 95% confidence intervals"
      ) +
      theme_minimal() +
      theme(panel.grid.major.y = element_blank(),
            plot.title = element_text(size = 16, color = 'black', face = 'bold'),
            axis.title = element_text(size = 14, color = 'black', face = 'bold'),
            axis.text = element_text(size = 12, color = 'black'),
            axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"))
  }
  
  return(p)
}

#' Quick summary function
summarize_correlation_results <- function(correlation_results) {
  cat("Correlation Analysis Summary\n")
  cat("============================\n")
  cat(sprintf("Merge Strategy: %s\n", correlation_results$merge_strategy))
  cat(sprintf("Number of Replicates: %d\n", correlation_results$n_replicates))
  cat(sprintf("Total Pairs: %d\n", correlation_results$n_pairs))
  cat("\nMethod Performance (sorted by |ρ|):\n")
  
  print(correlation_results$summary_table[, .(
    Method = method,
    ρ = round(rho, 3),
    `p-value` = ifelse(p_value < 0.001, "<0.001", round(p_value, 4)),
    `Mean Rep ρ` = round(mean_replicate_rho, 3),
    `SD Rep ρ` = round(sd_replicate_rho, 3),
    `N Pairs` = n_pairs
  )])
  
  cat("\nInterpretation:\n")
  cat("• Positive ρ: Higher inferred metric → Higher IBD (expected for IBD/IBS)\n")
  cat("• Negative ρ: Higher inferred metric → Lower IBD (expected for phylogenetic distance)\n")
  cat("• |ρ| > 0.7: Strong correlation\n")
  cat("• |ρ| 0.5-0.7: Moderate correlation\n")
  cat("• |ρ| 0.3-0.5: Weak correlation\n")
  cat("• |ρ| < 0.3: Negligible correlation\n")
}

# ------------------------
# Distribution plots of scores vs ground truth
# We'll create violin + jitter per method, separating positive vs negative under generational cutoff
# ------------------------
plot_score_distributions <- function(all_results, G_cut = 25, 
                                     merge_strategy = c("conservative", "strict", "union"),
                                     out_file = "score_distributions.png") {
  
  rows <- list(); k <- 0
  for (rep in names(all_results)) {
    rep_res <- all_results[[rep]]
    
    if (is.null(rep_res)) next
    
    eligible <- as.data.frame(copy(rep_res$merged_data[[merge_strategy]]))
    
    if (is.null(eligible) || nrow(eligible)==0) next
    
    # attempt to attach any score columns present for methods; we assume eligible has per-method score columns or single 'score' for a primary method
    # We'll reshape eligible so any column that is numeric and not pair_key/id becomes candidate score
    num_cols <- names(eligible)[sapply(eligible, is.numeric)]
    
    # we want columns besides total_ibd_bp, max_seg_bp, n_segments, is_positive, etc.
    skip_cols <- c("total_ibd_bp", "max_seg_bp", "n_segments","is_positive", "generation")
    cand_cols <- setdiff(num_cols, skip_cols)
    if (length(cand_cols)==0) next
    for (c in cand_cols) {
      k <- k + 1
      tmp <- as.data.table(eligible[, c("pair_key", "true_link", "ibd")])
      tmp[, method := c]
      tmp[, score := eligible[[c]]]
      tmp[, replicate := rep]
      rows[[k]] <- tmp
    }
  }
  
  if (length(rows)==0) return(NULL)
  df <- rbindlist(rows)
  df[, is_positive := factor(true_link, levels = c(0,1), labels = c("neg","pos"))]
  
  p <- ggplot(df, aes(x = method, y = score, fill = is_positive)) +
    geom_violin(position = position_dodge(width = 0.9), alpha = 0.6) +
    geom_jitter(aes(color = is_positive), size = 0.4, width = 0.15, alpha = 0.5) +
    facet_wrap(~ replicate, scales = "free_x") +
    theme_minimal() + 
    labs(y = "Score", x = "Method (score column name)", 
         title = paste0("Score distributions by method (G=", G_cut, ")")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(out_file, p, width = 14, height = 8)
  return(p)
}

# ------------------------
# Plot PR ribbons
# ------------------------
plot_pr_ribbon <- function(pr_summary_dt, out_file = "pr_ribbon.png") {
  if (is.null(pr_summary_dt) || nrow(pr_summary_dt)==0) return(NULL)
  p <- ggplot(pr_summary_dt, aes(x = recall, y = precision_mean, color = method, fill = method)) +
    geom_ribbon(aes(ymin = precision_lo, ymax = precision_hi), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.5) +
    theme_minimal() +
    labs(x = "Recall", y = "Precision", 
         title = "Precision–Recall", 
         subtitle = "Mean ±95% CI across replicates") +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 20, colour = 'black', face = 'bold'),
          plot.subtitle = element_text(size = 12, colour = 'black'),
          axis.title = element_text(size = 18, colour = 'black', face = 'bold'),
          axis.text = element_text(size = 15, colour = 'black'),
          axis.line = element_line(linewidth = 1.5, lineend = "square"),
          legend.text = element_text(size = 14, colour = 'black'),
          legend.title = element_blank()
    )
  ggsave(out_file, p, width = 12, height = 8)
  return(p)
}

# ------------------------
# Plot ROC ribbons
# ------------------------
plot_roc_ribbon <- function(roc_summary_dt, out_file = "roc_ribbon.png") {
  if (is.null(roc_summary_dt) || nrow(roc_summary_dt)==0) return(NULL)
  p <- ggplot(roc_summary_dt, aes(x = fpr, y = tpr_mean, color = method, fill = method)) +
    geom_ribbon(aes(ymin = tpr_lo, ymax = tpr_hi), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() +
    labs(x = "False positive rate", y = "True positive rate", title = "ROC (mean ±95% CI across replicates)") +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 20, colour = 'black', face = 'bold'),
          plot.subtitle = element_text(size = 12, colour = 'black'),
          axis.title = element_text(size = 18, colour = 'black', face = 'bold'),
          axis.text = element_text(size = 15, colour = 'black'),
          axis.line = element_line(linewidth = 1.5, lineend = "square"),
          legend.text = element_text(size = 14, colour = 'black'),
          legend.title = element_blank()
    )
  ggsave(out_file, p, width = 12, height = 8)
  return(p)
}

# ------------------------
# Correlation heatmap
# ------------------------
plot_correlation_heatmap <- function(cors_dt, out_file = "correlation_heatmap.png") {
  
  if (is.null(cors_dt) || nrow(cors_dt)==0) return(NULL)
  cors_dt[, method := factor(method, levels = method[order(-spearman)])]
  
  p <- ggplot(cors_dt, aes(x = 1, y = method, fill = spearman)) +
    geom_tile() + scale_fill_viridis_c(option = "D", limits = c(-1,1)) +
    geom_text(aes(label = round(spearman, 2))) +
    theme_minimal() + labs(x = "", y = "Method", title = "Spearman correlation vs IBD proportion (combined replicates)") +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 20, colour = 'black', face = 'bold'),
          plot.subtitle = element_text(size = 12, colour = 'black'),
          axis.title = element_text(size = 16, colour = 'black', face = 'bold'),
          axis.text.y = element_text(size = 14, colour = 'black'),
          axis.line = element_line(linewidth = 1.5, lineend = "square"),
          legend.text = element_text(size = 14, colour = 'black'),
          legend.title = element_blank()
    )
  
  ggsave(out_file, p, width = 6, height = 4)
  return(p)
}

# ------------------------
# Method ranking panel
# ------------------------
plot_method_ranking <- function(ranking_dt, out_file = "method_ranking.png") {
  if (is.null(ranking_dt) || nrow(ranking_dt)==0) return(NULL)
  agg <- ranking_dt[, .(
    aupr_mean = mean(aupr, na.rm = TRUE),
    aupr_lo = quantile(aupr, 0.025, na.rm = TRUE),
    aupr_hi = quantile(aupr, 0.975, na.rm = TRUE),
    auroc_mean = mean(auroc, na.rm = TRUE),
    auroc_lo = quantile(auroc, 0.025, na.rm = TRUE),
    auroc_hi = quantile(auroc, 0.975, na.rm = TRUE)
  ), by = method]
  
  # order methods by aupr_mean
  agg[, method := factor(method, levels = agg[order(-aupr_mean)]$method)]
  
  p1 <- ggplot(agg, aes(x = method, y = aupr_mean)) +
    geom_point(size = 3) + geom_errorbar(aes(ymin = aupr_lo, ymax = aupr_hi), width = 0.2) +
    theme_minimal() + labs(y = "AUPR (mean ±95% CI)", x = "", title = "Method ranking by AUPR") +
    theme(plot.title = element_text(size = 20, colour = 'black', face = 'bold'),
          plot.subtitle = element_text(size = 12, colour = 'black'),
          axis.title = element_text(size = 16, colour = 'black', face = 'bold'),
          axis.text.y = element_text(size = 14, colour = 'black'),
          axis.line = element_line(linewidth = 1.5, lineend = "square"),
          legend.text = element_text(size = 14, colour = 'black'),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 14, colour = 'black', angle = 45, hjust = 1))
  
  p2 <- ggplot(agg, aes(x = method, y = auroc_mean)) +
    geom_point(size = 3) + geom_errorbar(aes(ymin = auroc_lo, ymax = auroc_hi), width = 0.2) +
    theme_minimal() + labs(y = "AUROC (mean ±95% CI)", x = "", title = "Method ranking by AUROC") +
    theme(plot.title = element_text(size = 20, colour = 'black', face = 'bold'),
          plot.subtitle = element_text(size = 12, colour = 'black'),
          axis.title = element_text(size = 16, colour = 'black', face = 'bold'),
          axis.text.y = element_text(size = 14, colour = 'black'),
          axis.line = element_line(linewidth = 1.5, lineend = "square"),
          legend.text = element_text(size = 14, colour = 'black'),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 14, colour = 'black', angle = 45, hjust = 1))
  
  comb <- plot_grid(p1, p2, ncol = 1, labels = "AUTO")
  ggsave(out_file, comb, width = 10, height = 10)
  return(comb)
}

# ------------------------
# Composite figure: PR ribbon + ROC ribbon + correlation heatmap + ranking
# ------------------------
# assemble_composite_figure <- function(pr_plot, roc_plot, corr_plot, rank_plot, out_file = "composite_figure.png") {
#         # arrange
#         top <- plot_grid(pr_plot, roc_plot, ncol = 2, rel_widths = c(1,1))
#         bottom <- plot_grid(corr_plot, rank_plot, ncol = 2, rel_widths = c(0.5,1))
#         full <- plot_grid(top, bottom, ncol = 1, rel_heights = c(2,1))
#         ggsave(out_file, full, width = 16, height = 12)
#         return(full)
# }

assemble_composite_figure <- function(pr_plot, roc_plot, rank_plot, out_file = "composite_figure.png") {
  # arrange
  top <- plot_grid(pr_plot, rank_plot, ncol = 2, rel_widths = c(1,1))
  bottom <- plot_grid( roc_plot, ncol = 2, rel_widths = c(1,1))
  full <- plot_grid(top, bottom, ncol = 1, rel_heights = c(2,1))
  ggsave(out_file, full, width = 16, height = 12)
  return(full)
}


#' Plot correlation heatmap for method comparison
#' 
#' @param correlation_results Output from compute_correlations_across_replicates
#' @param out_file Output file path (optional)
#' @param show_ci Whether to show confidence intervals on the heatmap
#' @param sort_by Sorting method: "rho" (default), "mean_replicate", "method"
#' @param color_palette Color palette: "viridis" (default), "diverging", "sequential"
#' @param plot_type Type of heatmap: "single" (just ρ), "range" (with CI), "detailed"
#' @param font_size Base font size for plot elements
#' @return ggplot object
plot_correlation_heatmap <- function(correlation_results, 
                                     out_file = NULL, 
                                     show_ci = FALSE,
                                     sort_by = c("rho", "mean_replicate", "method"),
                                     color_palette = c("viridis", "diverging", "sequential"),
                                     plot_type = c("single", "range", "detailed"),
                                     font_size = 14) {
  
  # Validate input
  if (is.null(correlation_results)) {
    warning("Correlation results are NULL.")
    return(NULL)
  }
  
  if (is.null(correlation_results$summary_table) || 
      nrow(correlation_results$summary_table) == 0) {
    warning("No correlation data available.")
    return(NULL)
  }
  
  # Get parameters
  sort_by <- match.arg(sort_by)
  color_palette <- match.arg(color_palette)
  plot_type <- match.arg(plot_type)
  
  # Prepare data
  dt <- copy(correlation_results$summary_table)
  
  # Sort methods
  if (sort_by == "rho") {
    dt[, method := factor(method, levels = method[order(-abs(rho))])]
  } else if (sort_by == "mean_replicate") {
    dt[, method := factor(method, levels = method[order(-abs(mean_replicate_rho))])]
  } else {
    dt[, method := factor(method, levels = sort(method))]
  }
  
  # Determine color limits
  rho_range <- range(dt$rho, na.rm = TRUE)
  max_abs <- max(abs(rho_range))
  color_limits <- c(-max_abs, max_abs)
  
  # Choose color scale
  if (color_palette == "viridis") {
    fill_scale <- scale_fill_viridis_c(
      option = "D",
      limits = color_limits,
      na.value = "gray80"
    )
  } else if (color_palette == "diverging") {
    fill_scale <- scale_fill_gradient2(
      low = "#2166ac",  # Blue
      mid = "#f7f7f7",  # White
      high = "#b2182b", # Red
      midpoint = 0,
      limits = color_limits,
      na.value = "gray80"
    )
  } else {
    fill_scale <- scale_fill_gradientn(
      colors = c("#f7fbff", "#6baed6", "#08519c", "#08306b"),
      limits = c(0, max_abs),
      na.value = "gray80"
    )
  }
  
  # Create base plot based on plot_type
  if (plot_type == "single") {
    # Simple heatmap with correlation values
    p <- ggplot(dt, aes(x = 1, y = method, fill = rho)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = sprintf("%.3f\n(p=%s)", 
                                    rho, 
                                    ifelse(p_value < 0.001, "<0.001", 
                                           sprintf("%.3f", p_value)))), 
                size = font_size * 0.3, 
                fontface = "bold") +
      fill_scale +
      labs(
        title = "Correlation with Ground Truth IBD",
        subtitle = sprintf("Merge Strategy: %s | %d Replicates", 
                           correlation_results$merge_strategy,
                           correlation_results$n_replicates),
        x = "", 
        y = "",
        fill = "Spearman's ρ"
      )
    
  } else if (plot_type == "range") {
    # Heatmap with confidence intervals
    dt[, ci_range := ifelse(!is.na(ci_upper), 
                            sprintf("[%.3f, %.3f]", ci_lower, ci_upper), 
                            "NA")]
    
    p <- ggplot(dt, aes(x = 1, y = method, fill = rho)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = sprintf("%.3f\n%s", rho, ci_range)), 
                size = font_size * 0.25, 
                fontface = "bold") +
      fill_scale +
      labs(
        title = "Correlation with Ground Truth IBD (95% CI)",
        subtitle = sprintf("Merge Strategy: %s | %d Replicates", 
                           correlation_results$merge_strategy,
                           correlation_results$n_replicates),
        x = "", 
        y = "",
        fill = "Spearman's ρ"
      )
    
  } else if (plot_type == "detailed") {
    # Detailed heatmap with multiple metrics
    dt_long <- melt(dt, 
                    id.vars = "method",
                    measure.vars = c("rho", "mean_replicate_rho"),
                    variable.name = "metric",
                    value.name = "value")
    
    dt_long[, metric := factor(metric, 
                               levels = c("rho", "mean_replicate_rho"),
                               labels = c("Overall ρ", "Mean Replicate ρ"))]
    
    # Calculate max absolute value for scaling
    max_val <- max(abs(dt_long$value), na.rm = TRUE)
    
    p <- ggplot(dt_long, aes(x = metric, y = method, fill = value)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = sprintf("%.3f", value)), 
                size = font_size * 0.3, 
                fontface = "bold") +
      scale_fill_gradient2(
        low = "#2166ac",
        mid = "#f7f7f7",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-max_val, max_val),
        na.value = "gray80"
      ) +
      labs(
        title = "Detailed Correlation Analysis",
        subtitle = sprintf("Merge Strategy: %s | %d Replicates | %d Pairs", 
                           correlation_results$merge_strategy,
                           correlation_results$n_replicates,
                           correlation_results$n_pairs),
        x = "Metric", 
        y = "",
        fill = "Correlation"
      )
  }
  
  # Add common theme elements
  p <- p + theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(
        size = font_size + 4,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = font_size,
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      axis.title.y = element_text(
        size = font_size + 2,
        face = "bold",
        margin = margin(r = 10)
      ),
      axis.text.y = element_text(
        size = font_size,
        face = "bold"
      ),
      legend.title = element_text(
        size = font_size,
        face = "bold"
      ),
      legend.text = element_text(size = font_size - 2),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Save if output file specified
  if (!is.null(out_file)) {
    # Determine dimensions based on number of methods
    n_methods <- nrow(dt)
    height <- max(4, n_methods * 0.5 + 2)
    
    if (plot_type == "detailed") {
      width <- 8
    } else {
      width <- 6
    }
    
    ggsave(
      filename = out_file,
      plot = p,
      width = width,
      height = height,
      dpi = 300,
      bg = "white"
    )
    
    message(sprintf("Heatmap saved to: %s", out_file))
  }
  
  return(p)
}

#' Create a comparison heatmap across multiple merge strategies
#' 
#' @param results_list Named list of correlation results for different merge strategies
#' @param out_file Output file path (optional)
#' @param method_subset Subset of methods to include (NULL = all)
#' @param font_size Base font size
#' @return ggplot object
plot_strategy_comparison_heatmap <- function(results_list, 
                                             out_file = NULL,
                                             method_subset = NULL,
                                             font_size = 14) {
  
  # Validate input
  if (is.null(results_list) || length(results_list) == 0) {
    warning("No results provided.")
    return(NULL)
  }
  
  # Prepare combined data
  combined_data <- list()
  
  for (strategy_name in names(results_list)) {
    res <- results_list[[strategy_name]]
    if (is.null(res$summary_table)) next
    
    dt <- copy(res$summary_table)
    dt$strategy <- strategy_name
    combined_data[[strategy_name]] <- dt
  }
  
  if (length(combined_data) == 0) {
    warning("No valid data across strategies.")
    return(NULL)
  }
  
  # Combine and prepare
  dt_combined <- rbindlist(combined_data, fill = TRUE)
  
  # Filter methods if specified
  if (!is.null(method_subset)) {
    dt_combined <- dt_combined[method %in% method_subset]
  }
  
  # Order strategies and methods
  strategy_order <- names(results_list)
  dt_combined[, strategy := factor(strategy, levels = strategy_order)]
  
  # Order methods by average correlation across strategies
  method_avg <- dt_combined[, .(avg_rho = mean(abs(rho), na.rm = TRUE)), by = method]
  method_order <- method_avg[order(-avg_rho), method]
  dt_combined[, method := factor(method, levels = method_order)]
  
  # Create heatmap
  p <- ggplot(dt_combined, aes(x = strategy, y = method, fill = rho)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = sprintf("%.3f", rho)), 
              size = font_size * 0.3, 
              fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#f7f7f7",
      high = "#b2182b",
      midpoint = 0,
      na.value = "gray80",
      name = "Spearman's ρ"
    ) +
    labs(
      title = "Correlation Comparison Across Merge Strategies",
      subtitle = sprintf("%d Methods, %d Strategies", 
                         length(unique(dt_combined$method)),
                         length(unique(dt_combined$strategy))),
      x = "Merge Strategy",
      y = "Method"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        size = font_size + 4,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = font_size,
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      axis.title = element_text(
        size = font_size + 2,
        face = "bold"
      ),
      axis.text = element_text(size = font_size),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = font_size, face = "bold"),
      legend.text = element_text(size = font_size - 2),
      panel.grid = element_blank()
    )
  
  # Save if output file specified
  if (!is.null(out_file)) {
    n_methods <- length(unique(dt_combined$method))
    n_strategies <- length(unique(dt_combined$strategy))
    
    height <- max(4, n_methods * 0.5 + 2)
    width <- max(6, n_strategies * 1.5 + 2)
    
    ggsave(
      filename = out_file,
      plot = p,
      width = width,
      height = height,
      dpi = 300,
      bg = "white"
    )
    
    message(sprintf("Comparison heatmap saved to: %s", out_file))
  }
  
  return(p)
}




