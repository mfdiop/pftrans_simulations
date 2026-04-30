
# -------------------------------------------------------
# Define a function to classify true transmission pairs
# -------------------------------------------------------
true_pairs_by_tmrca <- function(df, g = 5) {
  # Filter pairs with TMRCA <= g generations
  close <- df %>% filter(generations <= g)
  
  # Create set of sorted pairs
  pairs <- close %>%
    rowwise() %>%
    mutate(pair = list(sort(c(Id1, Id2)))) %>%
    pull(pair) %>%
    lapply(function(x) paste(x, collapse = "_"))
  
  return(unique(unlist(pairs)))
}

# -----------------------------
# Compute Precision-Recall AUC 
# -----------------------------
compute_pr <- function(true_metric, inferred_metric, index = "ibd", generation = 5, outdir){
  #' Compute precision-recall AUC for IBD detection
  #' 
  #' @param true_ibd DataFrame of true pairwise IBD segments
  #' @param obs_ibd DataFrame of observed pairwise IBD segments
  #' @return AUC score

  
  # -----------------------------------------
  # Define true transmission pairs
  # Use IBD as a proxy to true transmission
  # -----------------------------------------
  
  # True pairs for g = 25 generations (recent)
  GT_IBD_25 <- true_pairs_by_tmrca(true_metric, g = generation)
  
  # -----------------------------
  # 4. Build labeled dataframe
  # -----------------------------
  label_data <- inferred_metric %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Id1, Id2)), collapse = "_")) %>%
    ungroup() %>%
    mutate(is_true = as.integer(pair %in% GT_IBD_25))
  
  # ----------------------
  # 5. Compute PR curve
  # ----------------------
  pr_curve <- pr.curve(
    scores.class0 = label_data[label_data$is_true == 1, ][[index]],
    scores.class1 = label_data[label_data$is_true == 0, ][[index]],
    curve = TRUE
  )
  
  pr_auc <- pr_curve$auc.integral
  cat(sprintf("PR-AUC = %.3f\n", pr_auc))
  
  pr_df <- data.frame(
    recall = pr_curve$curve[, 1],
    precision = pr_curve$curve[, 2])
  
  # Find optimal threshold based on your priorities
  # Option 1: Maximize F1 score (balance precision and recall)
  pr_df <- pr_df %>%
    mutate(f1_score = 2 * (precision * recall) / (precision + recall))
  
  optimal_f1 <- pr_df[which.max(pr_df$f1_score), ]
  cat(sprintf("Optimal F1: Precision=%.3f, Recall=%.3f, Threshold=%.3f\n",
              optimal_f1$precision, optimal_f1$recall, optimal_f1$f1_score))
  
  # Option 2: Require minimum precision
  min_precision <- 0.90
  optimal_recall <- pr_df %>%
    filter(precision >= min_precision) %>%
    filter(recall == max(recall, na.rm = TRUE))
  
  # Option 3: Require minimum recall
  min_recall <- 0.95
  optimal_precision <- pr_df %>%
    filter(recall >= min_recall) %>%
    filter(precision == max(precision, na.rm = TRUE))
  
  # -------------------
  # 6. Plot PR curve
  # -------------------
  
  plot_curve <- ggplot(pr_df, aes(x = recall, y = precision)) +
    geom_line(color = "blue", linewidth = 2, lineend = 'round') +
    labs(
      x = "Recall", y = "Precision",
      title = sprintf("Transmission recoverability from observed %s metric", index),
      subtitle = sprintf("AUC = %.3f", pr_auc)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 19), # hjust = 0.5, 
      plot.subtitle = element_text(face = "bold", size = 15),
      axis.title = element_text(size = 15, color = 'black', face = "bold"),
      axis.text = element_text(size = 12, color = 'black'))
  
  ggsave(file.path(outdir, "pr_auc_", index, ".pdf"), plot = plot_curve,
         width = 12, height = 8) # , units = 'mm'
  
}

# -----------  perform ROC curve analysis --------------

# Function to compute and plot ROC curve
compute_roc <- function(obs_rel, metric_col = "metric") {
  library(pROC)
  library(PRROC)
  
  # Compute ROC curve using pROC package
  roc_obj <- roc(obs_rel$is_true, obs_rel[[metric_col]], 
                 levels = c(0, 1), 
                 direction = "<")  # Use "<" if higher values = more related
  # Use ">" if lower values = more related (e.g., distance)
  
  # Get AUC
  auc_value <- auc(roc_obj)
  
  # Print results
  cat(sprintf("ROC-AUC = %.3f\n", auc_value))
  
  # Plot ROC curve
  plot(roc_obj, 
       main = sprintf("ROC Curve (AUC = %.3f)", auc_value),
       col = "blue", 
       lwd = 2,
       print.auc = TRUE,
       print.auc.x = 0.6,
       print.auc.y = 0.4)
  
  # Add diagonal reference line (random classifier)
  abline(a = 0, b = 1, lty = 2, col = "gray")
  
  return(list(roc = roc_obj, auc = auc_value))
}

# --------------------------------------------------
# Compute precision-recall AUC for IBD detection
# --------------------------------------------------
compute_pr <- function(true_data, obs_data, all_pairs) {
  #' Compute precision-recall AUC for IBD detection
  #' 
  #' @param true_ibd DataFrame of true pairwise IBD segments
  #' @param obs_ibd DataFrame of observed pairwise IBD segments
  #' @param all_pairs List of all possible sample pairs (as list of vectors or data frame)
  #' @return AUC score
  
  library(PRROC)
  
  # Create set of true pairs
  true_pairs <- true_data %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Id1, Id2)), collapse = "_")) %>%
    pull(pair) %>%
    unique()
  
  # Create set of observed pairs
  obs_pairs <- obs_data %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Id1, Id2)), collapse = "_")) %>%
    pull(pair) %>%
    unique()
  
  # Convert all_pairs to the same format
  # Handle different input formats for all_pairs
  if (is.list(all_pairs) && !is.data.frame(all_pairs)) {
    # If all_pairs is a list of vectors
    all_pairs_str <- sapply(all_pairs, function(p) paste(sort(p), collapse = "_"))
  } else if (is.data.frame(all_pairs)) {
    # If all_pairs is a data frame with two columns
    all_pairs_str <- all_pairs %>%
      rowwise() %>%
      mutate(pair = paste(sort(c(.[1], .[2])), collapse = "_")) %>%
      pull(pair)
  } else if (is.matrix(all_pairs)) {
    # If all_pairs is a matrix
    all_pairs_str <- apply(all_pairs, 1, function(p) paste(sort(p), collapse = "_"))
  }
  
  # Create binary labels
  y_true <- as.integer(all_pairs_str %in% true_pairs)
  y_score <- as.integer(all_pairs_str %in% obs_pairs)
  
  # Compute PR curve
  pr_curve <- pr.curve(
    scores.class0 = y_score[y_true == 1],
    scores.class1 = y_score[y_true == 0],
    curve = TRUE
  )
  
  auc_score <- pr_curve$auc.integral
  
  return(auc_score)
}

# Alternative version using yardstick package (more similar to sklearn):

compute_pr_v1 <- function(true_ibd, obs_ibd, all_pairs) {
  library(yardstick)
  library(dplyr)
  
  # Create set of true pairs
  true_pairs <- true_ibd %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Id1, Id2)), collapse = "_")) %>%
    pull(pair) %>%
    unique()
  
  # Create set of observed pairs
  obs_pairs <- obs_ibd %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Id1, Id2)), collapse = "_")) %>%
    pull(pair) %>%
    unique()
  
  # Convert all_pairs to string format
  if (is.list(all_pairs) && !is.data.frame(all_pairs)) {
    all_pairs_str <- sapply(all_pairs, function(p) paste(sort(p), collapse = "_"))
  } else if (is.data.frame(all_pairs) || is.matrix(all_pairs)) {
    all_pairs_str <- apply(as.matrix(all_pairs), 1, function(p) paste(sort(p), collapse = "_"))
  }
  
  # Create data frame for evaluation
  eval_df <- data.frame(
    truth = factor(ifelse(all_pairs_str %in% true_pairs, "positive", "negative"),
                   levels = c("positive", "negative")),
    estimate = as.numeric(all_pairs_str %in% obs_pairs)
  )
  
  # Compute PR AUC
  auc_score <- pr_auc(eval_df, truth, estimate, event_level = "first")$.estimate
  
  return(auc_score)
}


topk_overlap <- function(true_matrix, inferred_matrix, k = 100) {
  # Get indices of top k values (flattened matrix)
  true_top <- order(as.vector(true_matrix), decreasing = TRUE)[1:k]
  inf_top <- order(as.vector(inferred_matrix), decreasing = TRUE)[1:k]
  
  # Calculate overlap
  overlap <- length(intersect(true_top, inf_top)) / k
  
  return(overlap)
}



# ROC Curve with ggplot2
compute_roc_ggplot <- function(obs_rel, metric_col = "metric", 
                               direction = "<") {
  
  # Compute ROC
  roc_obj <- roc(obs_rel$is_true, obs_rel[[metric_col]], 
                 levels = c(0, 1), 
                 direction = direction)
  
  auc_value <- auc(roc_obj)
  
  # Extract coordinates
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,  # False Positive Rate
    tpr = roc_obj$sensitivities,      # True Positive Rate (Recall)
    threshold = roc_obj$thresholds
  )
  
  # Create plot
  p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "blue", linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", color = "gray50") +
    annotate("text", x = 0.7, y = 0.3, 
             label = sprintf("AUC = %.3f", auc_value), 
             size = 5, fontface = "bold") +
    labs(
      title = "ROC Curve for IBD Detection through HMM",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity/Recall)",
      subtitle = "Diagonal line represents random classifier"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 15),
      axis.title = element_text(size = 15, , colour = 'black', face = 'bold'),
      axis.text = element_text(size = 12, colour = 'black')
    ) +
    coord_fixed()  # Equal aspect ratio
  
  print(p)
  
  ggsave(paste0("single_pop/out_true_ibd/roc_", metric_col, ".png"),
         plot = p, width = 12, height = 8)
  
  return(list(roc = roc_obj, auc = auc_value, data = roc_df))
}




# Compare Multiple Metrics (Multiple ROC Curves)
compare_roc_curves <- function(data_list, metric_names) {
  #' data_list: list of data frames, each with 'is_true' and a metric column
  #' metric_names: vector of names for each metric
  
  library(pROC)
  
  # Compute ROC for each metric
  roc_list <- list()
  auc_values <- c()
  
  for (i in seq_along(data_list)) {
    df <- data_list[[i]]
    roc_obj <- roc(df$is_true, df$metric, 
                   levels = c(0, 1), 
                   direction = "<")
    roc_list[[i]] <- roc_obj
    auc_values[i] <- auc(roc_obj)
  }
  
  # Combine data for plotting
  roc_combined <- map_df(seq_along(roc_list), function(i) {
    data.frame(
      fpr = 1 - roc_list[[i]]$specificities,
      tpr = roc_list[[i]]$sensitivities,
      method = metric_names[i],
      auc = auc_values[i]
    )
  })
  
  # Create plot
  p <- ggplot(roc_combined, aes(x = fpr, y = tpr, color = method)) +
    geom_line(size = 1.2) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", color = "gray50") +
    labs(
      title = "ROC Curve Comparison",
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = "Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    coord_fixed()
  
  # Add AUC values to legend
  auc_labels <- sprintf("%s (AUC=%.3f)", metric_names, auc_values)
  p <- p + scale_color_discrete(labels = auc_labels)
  
  print(p)
  
  return(list(roc_list = roc_list, auc_values = auc_values, plot_data = roc_combined))
}





compute_roc_prroc <- function(obs_rel, metric_col = "metric") {
  
  # Separate scores for positive and negative classes
  scores_pos <- obs_rel[[metric_col]][obs_rel$is_true == 1]
  scores_neg <- obs_rel[[metric_col]][obs_rel$is_true == 0]
  
  # Compute ROC curve
  roc_obj <- roc.curve(
    scores.class0 = scores_pos,
    scores.class1 = scores_neg,
    curve = TRUE
  )
  
  auc_value <- roc_obj$auc
  
  cat(sprintf("ROC-AUC = %.3f\n", auc_value))
  
  # Plot
  plot(roc_obj, 
       main = sprintf("ROC Curve (AUC = %.3f)", auc_value),
       col = "blue",
       lwd = 2)
  
  return(roc_obj)
}

# Usage
roc_result <- compute_roc_prroc(obs_rel, metric_col = "metric")

# Complete ROC Analysis Function
perform_roc_analysis <- function(true_ibd, obs_rel, g = 5, 
                                 metric_col = "metric",
                                 genome_length = 30140) {
  
  library(pROC)
  library(ggplot2)
  
  # 1. Create ground truth pairs
  GT_pairs <- true_ibd %>%
    filter(Tmrca <= g) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Id1, Id2)), collapse = "_")) %>%
    pull(pair) %>%
    unique()
  
  # 2. Label observed pairs
  obs_labeled <- obs_rel %>%
    mutate(
      pair = paste(pmin(Id1, Id2), pmax(Id1, Id2), sep = "_"),
      is_true = as.integer(pair %in% GT_pairs)
    )
  
  # 3. Check metric separation
  cat("=== Metric Separation ===\n")
  separation <- obs_labeled %>%
    group_by(is_true) %>%
    summarise(
      mean = mean(.data[[metric_col]]),
      median = median(.data[[metric_col]]),
      sd = sd(.data[[metric_col]]),
      .groups = 'drop'
    )
  print(separation)
  
  # 4. Compute ROC
  roc_obj <- roc(obs_labeled$is_true, obs_labeled[[metric_col]], 
                 levels = c(0, 1), 
                 direction = "<")  # Adjust based on your metric
  
  auc_value <- auc(roc_obj)
  
  # 5. Find optimal threshold (Youden's index)
  coords_all <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
  youden <- coords_all$sensitivity + coords_all$specificity - 1
  optimal_idx <- which.max(youden)
  optimal_threshold <- coords_all$threshold[optimal_idx]
  optimal_sens <- coords_all$sensitivity[optimal_idx]
  optimal_spec <- coords_all$specificity[optimal_idx]
  
  cat(sprintf("\n=== Optimal Threshold ===\n"))
  cat(sprintf("Threshold: %.4f\n", optimal_threshold))
  cat(sprintf("Sensitivity: %.3f\n", optimal_sens))
  cat(sprintf("Specificity: %.3f\n", optimal_spec))
  
  # 6. Plot ROC curve
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
  
  p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "blue", size = 1.2) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", color = "gray50") +
    geom_point(aes(x = 1 - optimal_spec, y = optimal_sens),
               color = "red", size = 3) +
    annotate("text", x = 0.7, y = 0.3, 
             label = sprintf("AUC = %.3f", auc_value), 
             size = 5) +
    annotate("text", x = 1 - optimal_spec + 0.1, y = optimal_sens + 0.05,
             label = "Optimal", color = "red", size = 3) +
    labs(
      title = sprintf("ROC Curve (TMRCA ≤ %d generations)", g),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal() +
    coord_fixed()
  
  print(p)
  
  # 7. Confusion matrix at optimal threshold
  obs_labeled <- obs_labeled %>%
    mutate(predicted = as.integer(.data[[metric_col]] >= optimal_threshold))
  
  conf_mat <- table(Predicted = obs_labeled$predicted, 
                    Actual = obs_labeled$is_true)
  
  cat("\n=== Confusion Matrix at Optimal Threshold ===\n")
  print(conf_mat)
  
  # Calculate performance metrics
  if (all(dim(conf_mat) == c(2, 2))) {
    tn <- conf_mat[1, 1]
    fp <- conf_mat[2, 1]
    fn <- conf_mat[1, 2]
    tp <- conf_mat[2, 2]
    
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    f1 <- 2 * (precision * recall) / (precision + recall)
    
    cat(sprintf("\nPrecision: %.3f\n", precision))
    cat(sprintf("Recall: %.3f\n", recall))
    cat(sprintf("F1-Score: %.3f\n", f1))
  }
  
  return(list(
    roc = roc_obj,
    auc = auc_value,
    optimal_threshold = optimal_threshold,
    data = obs_labeled,
    plot = p
  ))
}


compare_roc_pr <- function(obs_rel, metric_col = "metric") {
  
  library(pROC)
  library(PRROC)
  library(patchwork)
  
  # ROC curve
  roc_obj <- roc(obs_rel$is_true, obs_rel[[metric_col]], 
                 levels = c(0, 1), direction = "<")
  roc_auc <- auc(roc_obj)
  
  roc_df <- data.frame(
    x = 1 - roc_obj$specificities,
    y = roc_obj$sensitivities
  )
  
  p_roc <- ggplot(roc_df, aes(x = x, y = y)) +
    geom_line(color = "blue", size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    annotate("text", x = 0.7, y = 0.3, 
             label = sprintf("AUC = %.3f", roc_auc), size = 5) +
    labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal() +
    coord_fixed()
  
  # PR curve
  pr_obj <- pr.curve(
    scores.class0 = obs_rel[[metric_col]][obs_rel$is_true == 1],
    scores.class1 = obs_rel[[metric_col]][obs_rel$is_true == 0],
    curve = TRUE
  )
  pr_auc <- pr_obj$auc.integral
  
  pr_df <- data.frame(
    x = pr_obj$curve[, 1],  # Recall
    y = pr_obj$curve[, 2]   # Precision
  )
  
  baseline <- mean(obs_rel$is_true)
  
  p_pr <- ggplot(pr_df, aes(x = x, y = y)) +
    geom_line(color = "red", size = 1.2) +
    geom_hline(yintercept = baseline, linetype = "dashed", color = "gray") +
    annotate("text", x = 0.3, y = 0.7, 
             label = sprintf("AUC = %.3f", pr_auc), size = 5) +
    labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
    theme_minimal() +
    coord_fixed()
  
  # Combine plots
  combined_plot <- p_roc + p_pr + 
    plot_annotation(title = "ROC vs PR Curve Comparison")
  
  print(combined_plot)
  
  cat(sprintf("ROC-AUC: %.3f\n", roc_auc))
  cat(sprintf("PR-AUC: %.3f\n", pr_auc))
  cat(sprintf("Baseline (random): %.3f\n", baseline))
  
  return(list(roc_auc = roc_auc, pr_auc = pr_auc, 
              baseline = baseline, plot = combined_plot))
}




































