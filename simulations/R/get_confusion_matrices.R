#' Calculate Confusion Matrices for Multiple Inference Metrics
#' 
#' @param data A dataframe containing ground truth labels and inference metrics
#' @param truth_col Name of the column containing ground truth labels (TRUE/FALSE or 1/0)
#' @param metric_cols Character vector of column names for the inference metrics to evaluate
#' @param threshold_method Method for determining classification threshold: 
#'        "fixed" (use specific cutoff), "optimal_youden" (maximize Youden's J statistic),
#'        "optimal_f1" (maximize F1 score), or "median" (use median of metric)
#' @param fixed_thresholds Named vector of fixed thresholds for each metric (if threshold_method = "fixed")
#' @param positive_class Value indicating positive class in truth_col (default: TRUE or 1)
#' 
#' @return A list containing:
#'   - confusion_matrices: List of confusion matrices for each metric
#'   - metrics_summary: Dataframe with performance metrics (accuracy, precision, recall, F1, etc.)
#'   - thresholds_used: Vector of thresholds actually used for each metric
#'   - detailed_results: List with full classification results for each metric
#'
#' @examples
#' results <- get_confusion_matrices(
#'   data = merge_metrics,
#'   truth_col = "ground_truth",
#'   metric_cols = c("IBD", "IBS", "phylo_dist"),
#'   threshold_method = "optimal_youden"
#' )


get_confusion_matrices <- function(data, 
                                   truth_col = "ground_truth", 
                                   metric_cols = c("IBD", "IBS", "phylo_dist"),
                                   threshold_method = "optimal_youden",
                                   fixed_thresholds = NULL,
                                   positive_class = TRUE) {
  
  # Input validation
  if (!truth_col %in% names(data)) {
    stop("Ground truth column '", truth_col, "' not found in data")
  }
  
  missing_metrics <- setdiff(metric_cols, names(data))
  if (length(missing_metrics) > 0) {
    stop("The following metric columns are missing: ", paste(missing_metrics, collapse = ", "))
  }
  
  if (!threshold_method %in% c("fixed", "optimal_youden", "optimal_f1", "median")) {
    stop("threshold_method must be one of: 'fixed', 'optimal_youden', 'optimal_f1', 'median'")
  }
  
  # Initialize results storage
  confusion_matrices <- list()
  performance_metrics <- data.frame()
  thresholds_used <- numeric()
  detailed_results <- list()
  
  # Extract ground truth (convert to logical if needed)
  ground_truth <- as.logical(data[[truth_col]])
  
  for (metric in metric_cols) {
    # Extract metric values
    metric_values <- data[[metric]]
    
    # Determine threshold for this metric
    threshold <- determine_threshold(
      metric_values = metric_values,
      ground_truth = ground_truth,
      metric_name = metric,
      method = threshold_method,
      fixed_thresholds = fixed_thresholds
    )
    
    thresholds_used[metric] <- threshold
    
    # Classify predictions based on threshold
    # Note: Assumes higher values indicate stronger relationship for all metrics
    # For phylogenetic distance (lower = closer relationship), we'll handle inversion
    predictions <- classify_relationship(
      metric_values = metric_values,
      threshold = threshold,
      metric_name = metric
    )
    
    # Create confusion matrix
    cm <- table(
      Actual = factor(ground_truth, levels = c(TRUE, FALSE)),
      Predicted = factor(predictions, levels = c(TRUE, FALSE))
    )
    
    confusion_matrices[[metric]] <- cm
    
    # Calculate performance metrics
    perf_metrics <- calculate_performance_metrics(cm)
    
    # Store detailed results
    detailed_results[[metric]] <- list(
      predictions = predictions,
      metric_values = metric_values,
      threshold = threshold,
      confusion_matrix = cm,
      performance_metrics = perf_metrics
    )
    
    # Add to summary dataframe
    perf_df <- data.frame(
      Metric = metric,
      Threshold = round(threshold, 4),
      Accuracy = round(perf_metrics$accuracy, 4),
      Precision = round(perf_metrics$precision, 4),
      Recall = round(perf_metrics$recall, 4),
      F1_Score = round(perf_metrics$f1, 4),
      Specificity = round(perf_metrics$specificity, 4),
      NPV = round(perf_metrics$npv, 4),  # Negative Predictive Value
      TP = cm["TRUE", "TRUE"],
      FP = cm["FALSE", "TRUE"],
      TN = cm["FALSE", "FALSE"],
      FN = cm["TRUE", "FALSE"]
    )
    
    performance_metrics <- rbind(performance_metrics, perf_df)
  }
  
  # Return comprehensive results
  return(list(
    confusion_matrices = confusion_matrices,
    metrics_summary = performance_metrics,
    thresholds_used = thresholds_used,
    detailed_results = detailed_results
  ))
}

# Helper function to determine classification threshold
determine_threshold <- function(metric_values, ground_truth, metric_name, 
                                method = "optimal_youden", fixed_thresholds = NULL) {
  
  if (method == "fixed") {
    if (is.null(fixed_thresholds) || !metric_name %in% names(fixed_thresholds)) {
      stop("Fixed threshold method requires fixed_thresholds named vector with entry for '", metric_name, "'")
    }
    return(fixed_thresholds[[metric_name]])
  }
  
  if (method == "median") {
    return(median(metric_values, na.rm = TRUE))
  }
  
  # For optimal threshold methods, use ROC analysis
  if (method %in% c("optimal_youden", "optimal_f1")) {
    
    # Handle metrics where lower values indicate stronger relationships (e.g., phylogenetic distance)
    if (grepl("phylo|distance", metric_name, ignore.case = TRUE)) {
      # Invert for ROC analysis (higher = better for ROC)
      roc_obj <- pROC::roc(ground_truth, -metric_values, quiet = TRUE)
      # roc_obj <- pROC::roc(ground_truth, metric_values, quiet = TRUE) # because the distance was transformed to similarity metric
    } else {
      roc_obj <- pROC::roc(ground_truth, metric_values, quiet = TRUE)
    }
    
    if (method == "optimal_youden") {
      optimal_point <- pROC::coords(roc_obj, "best", best.method = "youden")
    } else { # optimal_f1
      optimal_point <- pROC::coords(roc_obj, "best", best.method = "f1")
    }
    
    threshold <- optimal_point$threshold
    
    # Convert back to original scale if we inverted
    if (grepl("phylo|distance", metric_name, ignore.case = TRUE)) {
      threshold <- -threshold
      # threshold <- threshold # because the distance was transformed to similarity metric
    }
    
    return(threshold)
  }
}

# Helper function for classification
classify_relationship <- function(metric_values, threshold, metric_name) {
  # Handle direction of relationship
  # For most metrics (IBD, IBS), higher values = stronger relationship
  # For phylogenetic distance, lower values = stronger relationship
  
  if (grepl("phylo|distance", metric_name, ignore.case = TRUE)) {
    # Lower values indicate closer relationships
    return(metric_values <= threshold)
    # return(metric_values >= threshold) # distance was transformed to similarity metric
  } else {
    # Higher values indicate closer relationships
    return(metric_values >= threshold)
  }
}

# Helper function to calculate performance metrics
calculate_performance_metrics <- function(confusion_matrix) {
  TP <- confusion_matrix["TRUE", "TRUE"]
  FP <- confusion_matrix["FALSE", "TRUE"]
  TN <- confusion_matrix["FALSE", "FALSE"]
  FN <- confusion_matrix["TRUE", "FALSE"]
  
  accuracy <- (TP + TN) / sum(confusion_matrix)
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)  # Also called sensitivity
  specificity <- TN / (TN + FP)
  f1 <- 2 * (precision * recall) / (precision + recall)
  npv <- TN / (TN + FN)  # Negative Predictive Value
  
  return(list(
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    specificity = specificity,
    f1 = f1,
    npv = npv
  ))
}


# ROC and AUPRC for a single predictor column
compute_roc_aupr <- function(df, score_col = "IBD_prop", label_col = "is_true") {
  
  # Input validation
  if (!score_col %in% names(df)) {
    stop("Score column '", score_col, "' not found")
  }
  if (!label_col %in% names(df)) {
    stop("Label column '", label_col, "' not found")
  }
  
  scores <- df[[score_col]]
  labels <- df[[label_col]]
  
  # Remove NA values from both scores and labels
  complete_cases <- !is.na(scores) & !is.na(labels)
  scores_clean <- scores[complete_cases]
  labels_clean <- labels[complete_cases]
  
  cat("   Using", sum(complete_cases), "complete cases out of", length(scores), "total\n")
  
  # Check if we have enough data
  if (sum(complete_cases) < 10) {
    warning("Insufficient complete cases (n = ", sum(complete_cases), ") for reliable evaluation")
    return(list(roc = NULL, auc = NA, pr = NULL, aupr = NA, n_used = sum(complete_cases)))
  }
  
  # Check class balance
  class_table <- table(labels_clean)
  if (length(class_table) < 2) {
    warning("Only one class present in complete cases - cannot compute ROC/PR")
    return(list(roc = NULL, auc = NA, pr = NULL, aupr = NA, n_used = sum(complete_cases)))
  }
  
  if (min(class_table) < 5) {
    warning("Very small minority class (n = ", min(class_table), ") - results may be unstable")
  }
  
  tryCatch({
    # ROC / AUC
    roc_obj <- pROC::roc(labels_clean, scores_clean, quiet = TRUE, direction = "<")
    auc <- pROC::auc(roc_obj)
    
    # AUPRC with additional safety checks
    scores_class0 <- scores_clean[labels_clean == 1]
    scores_class1 <- scores_clean[labels_clean == 0]
    
    # Check if we have samples in both classes for PR curve
    if (length(scores_class0) == 0 || length(scores_class1) == 0) {
      warning("No samples in one of the classes for PR curve calculation")
      return(list(roc = roc_obj, auc = auc, pr = NULL, aupr = NA, n_used = sum(complete_cases)))
    }
    
    pr_obj <- PRROC::pr.curve(
      scores.class0 = scores_class0, 
      scores.class1 = scores_class1, 
      curve = TRUE
    )
    aupr <- pr_obj$auc.integral
    
    return(list(roc = roc_obj, auc = auc, pr = pr_obj, aupr = aupr, n_used = sum(complete_cases)))
    
  }, error = function(e) {
    warning("Error in ROC/PR calculation: ", e$message)
    return(list(roc = NULL, auc = NA, pr = NULL, aupr = NA, n_used = sum(complete_cases), error = e$message))
  })
}

# precision@k / recall@k (directional)
# We assume pairs are candidate edges id1 -> id2. Group by id2 (target infectee), rank candidates by score descending,
# compute for each id2 precision@k and recall@k across its true infectors (rare).
precision_recall_at_k <- function(pairs_df, score_col = "IBD_prop", label_col = "is_true", ks = c(1,3,5)) {
  df <- pairs_df %>%
    arrange(id2, desc(.data[[score_col]])) %>%
    group_by(id2) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  results <- tibble(k = ks, precision = NA_real_, recall = NA_real_)
  for (i in seq_along(ks)) {
    k <- ks[i]
    topk <- df %>% filter(rank <= k)
    # per-target precision: for each id2, precision = (# true among topk candidates)/k
    per_target <- topk %>% group_by(id2) %>%
      summarise(tp = sum(.data[[label_col]] == 1), n = n()) %>%
      mutate(prec = tp / k)
    # average precision across targets (macro-averaged)
    avg_prec <- mean(per_target$prec, na.rm = TRUE)
    # recall: need ground truth positives per target (#true infectors)
    true_per_target <- df %>% group_by(id2) %>% summarise(true_total = sum(.data[[label_col]] == 1))
    # compute recall per target as tp/top_true (avoid division by zero)
    per_target <- per_target %>% left_join(true_per_target, by = "id2") %>%
      mutate(rec = ifelse(true_total > 0, tp / true_total, NA))
    avg_rec <- mean(per_target$rec, na.rm = TRUE)
    results$precision[i] <- avg_prec
    results$recall[i] <- avg_rec
  }
  return(results)
}


evaluate_predictors <- function(df, predictors = c("IBD_prop","IBS_prop","phylo_sim_exp"), 
                                label_col = "is_true", 
                                ks = c(1,3,5),
                                ...) {
  out_list <- list()
  for (pred in predictors) {
    roc_pr <- compute_roc_aupr(df, score_col = pred, label_col = label_col)
    prk <- precision_recall_at_k(df, score_col = pred, label_col = label_col, ks = ks)
    out_list[[pred]] <- list(all = roc_pr,
                             auc = as.numeric(roc_pr$auc), 
                             aupr = as.numeric(roc_pr$aupr), 
                             prk = prk, 
                             roc = roc_pr$roc)
  }
  return(out_list)
}

# Evaluate per scenario and gather AUC/AUPR and precision@k
evaluate_across_scenarios <- function(pairs_df, predictors = c("IBD_prop","IBS_prop","phylo_sim_exp"), label_col = "is_true", ks = c(1,3,5)) {
  scenarios <- unique(pairs_df$scenario)
  results <- list()
  for (sc in scenarios) {
    df_sc <- filter(pairs_df, scenario == sc)
    evals <- evaluate_predictors(df_sc, predictors = predictors, label_col = label_col, ks = ks)
    # tidy summary for each predictor
    summary_rows <- map_dfr(names(evals), function(pred) {
      tibble(scenario = sc,
             predictor = pred,
             auc = evals[[pred]]$auc,
             aupr = evals[[pred]]$aupr) %>%
        bind_cols(bind_rows(evals[[pred]]$prk %>% mutate(.pred = pred)))
    })
    results[[sc]] <- summary_rows
  }
  final_df <- bind_rows(results)
  return(final_df)
}

compute_auc <- function(df) {
  out <- list()
  
  # Only retain pairs with known truth
  d <- df %>% filter(!is.na(direct_link))
  
  out$ROC_IBD   <- roc(d$direct_link, d$IBD)
  out$ROC_IBS   <- roc(d$direct_link, d$IBS)
  out$ROC_PHYLO <- roc(d$direct_link, d$phylo_sim)
  
  out$AUPRC_IBD   <- pr.curve(scores.class0 = d$IBD[d$direct_link==1],
                              scores.class1 = d$IBD[d$direct_link==0])
  out$AUPRC_IBS   <- pr.curve(scores.class0 = d$IBS[d$direct_link==1],
                              scores.class1 = d$IBS[d$direct_link==0])
  out$AUPRC_PHYLO <- pr.curve(scores.class0 = d$phylo_sim[d$direct_link==1],
                              scores.class1 = d$phylo_sim[d$direct_link==0])
  
  out
}

precision_at_k <- function(df, metric = "IBD", k = 3) {
  df2 <- df %>% filter(!is.na(direct_link))
  
  df_ranked <- df2 %>%
    arrange(desc(.data[[metric]])) %>%
    mutate(rank = row_number())
  
  topk <- df_ranked %>% slice(1:k)
  
  precision <- mean(topk$direct_link == 1)
  recall    <- sum(topk$direct_link == 1) / sum(df2$direct_link == 1)
  
  list(precision = precision, recall = recall)
}




compute_confusion <- function(df, score_col, truth_col = "true_link", threshold) {
  
  predicted <- df[[score_col]] >= threshold
  truth <- df[[truth_col]]
  
  TP <- sum(predicted & truth, na.rm = TRUE)
  FP <- sum(predicted & !truth, na.rm = TRUE)
  TN <- sum(!predicted & !truth, na.rm = TRUE)
  FN <- sum(!predicted & truth, na.rm = TRUE)
  
  list(
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    confusion_matrix = matrix(
      c(TP, FP, FN, TN),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("Pred_Pos", "Pred_Neg"),
        c("True_Pos", "True_Neg")
      )
    )
  )
}



compute_roc_pr <- function(df, score_col, truth_col = "is_true") {
  
  # ROC
  roc_obj <- pROC::roc(df[[truth_col]], df[[score_col]], quiet = TRUE)
  
  # PR Curve (needs scores for positives & negatives)
  fg <- df[df[[truth_col]] == TRUE,  ][[score_col]]
  bg <- df[df[[truth_col]] == FALSE, ][[score_col]]
  
  pr_obj <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
  
  list(
    roc = roc_obj,
    pr = pr_obj,
    auc_roc = pROC::auc(roc_obj),
    auc_pr = pr_obj$auc.integral
  )
}
