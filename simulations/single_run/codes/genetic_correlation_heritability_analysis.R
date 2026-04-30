
library(tidyverse)
library(data.table)

results_list <- readRDS("simulations/single_run/evaluation_results.rds")

split_pair_key <- function(df, col_name, keep_original = TRUE) {
  result_df <- df
  
  # Extract all tsk patterns (tsk followed by a number)
  # This regex captures "tsk" followed by optional separator and digits
  tsk_matches <- str_extract_all(
    as.character(result_df[[col_name]]), 
    "tsk[ _]*\\d+",
    simplify = FALSE
  )
  
  # Clean the matches (remove underscores and spaces)
  clean_tsk <- function(x) {
    # x <- gsub("[ _]", "", x)
    x <- gsub("[ ]", "", x)
    return(x)
  }
  
  # Apply cleaning and ensure we have at least 2 elements
  cleaned <- lapply(tsk_matches, function(x) {
    cleaned <- clean_tsk(x)
    length(cleaned) <- 2  # Ensure length 2
    return(cleaned)
  })
  
  # Add columns
  result_df[["id1"]] <- sapply(cleaned, `[`, 1)
  result_df[["id2"]] <- sapply(cleaned, `[`, 2)
  
  # Remove original if specified
  if (!keep_original) {
    result_df[[col_name]] <- NULL
  }
  
  return(result_df[, c("id1", "id2")])
}


build_sym_matrix <- function(pairs_df, value_col) {
  pairs_df <- as.data.table(pairs_df)
  
  pairs <- split_pair_key(pairs_df, "pair_key") 
  pairs_df[, id1 := pairs$id1]; pairs_df[, id2 := pairs$id2]
  
  ids <- unique(c(pairs$id1, pairs$id2))
  n <- length(ids)
  m <- matrix(NA, nrow = n, ncol = n, dimnames = list(ids, ids))
  
  for (r in seq_len(nrow(pairs_df))) {
    i <- as.character(pairs_df$id1[r]); j <- as.character(pairs_df$id2[r])
    if (i %in% ids & j %in% ids) {
      m[i,j] <- pairs_df[[value_col]][r]; m[j,i] <- m[i,j]
    }
  }
  diag(m) <- 1
  m[is.na(m)] <- 0
  return(m)
}

# 1. Genetic Correlation & Heritability Analysis
#' Analyze genetic parameters across replicates
#' 
#' @param results_list List of results from 20 replicate runs
#' @param trait_values Optional: Simulated trait values for heritability analysis

analyze_genetic_parameters <- function(results_list, trait_values = NULL) {
  
  cat("🧬 Analyzing Genetic Parameters Across Replicates...\n")
  
  # Extract correlation matrices and compute genetic parameters
  genetic_stats <- imap_dfr(results_list, function(replicate, replicate_name) {
    
    # Get the merged data (using conservative merge)
    merged_data <- replicate$merged_data$conservative
    
    # Print progress - use replicate_name (string) not replicate (list)
    cat("Merged data with ", nrow(merged_data), "pairs for replicate ", replicate_name, "\n")
    
    # Build relationship matrices
    true_ibd <- build_sym_matrix(merged_data, "ibd")
    K_ibd <- build_sym_matrix(merged_data, "IBD")
    K_ibs <- build_sym_matrix(merged_data, "IBS") 
    K_phylo <- build_sym_matrix(merged_data, "phylo")
    
    # Compute matrix correlations with ground truth
    # (Assuming you have ground truth relationship matrix)
    # NOTE: You might want to check these correlations - they seem to be correlating
    # different matrices with each other, not with ground truth
    ground_truth_cor <- c(
      IBD = cor(as.vector(K_ibd), as.vector(true_ibd), use = "complete.obs"),
      IBS = cor(as.vector(K_ibs), as.vector(true_ibd), use = "complete.obs"),
      Phylo = cor(as.vector(K_phylo), as.vector(true_ibd), use = "complete.obs")
    )
    
    # Compute matrix properties
    matrix_properties <- function(mat) {
      eigen_vals <- eigen(mat)$values
      list(
        effective_loci = sum(eigen_vals^2) / sum(eigen_vals)^2, # Effective number of loci
        condition_number = max(eigen_vals) / min(eigen_vals[eigen_vals > 1e-10]),
        mean_relatedness = mean(mat[upper.tri(mat)])
      )
    }
    
    props_ibd <- matrix_properties(K_ibd)
    props_ibs <- matrix_properties(K_ibs)
    props_phylo <- matrix_properties(K_phylo)
    
    # Return summary for this replicate
    tibble(
      replicate = replicate_name,  # Use the name from imap_dfr
      IBD_truth = ground_truth_cor["IBD"],
      IBS_truth = ground_truth_cor["IBS"],
      Phylo_truth = ground_truth_cor["Phylo"],
      IBD_effective_loci = props_ibd$effective_loci,
      IBS_effective_loci = props_ibs$effective_loci,
      Phylo_effective_loci = props_phylo$effective_loci,
      IBD_mean_relatedness = props_ibd$mean_relatedness,
      IBS_mean_relatedness = props_ibs$mean_relatedness,
      Phylo_mean_relatedness = props_phylo$mean_relatedness
    )
  })
  
  # Alternative: If you want numeric replicate indices
  # genetic_stats <- map_dfr(seq_along(results_list), function(i) {
  #   replicate <- results_list[[i]]
  #   # ... same code ...
  #   tibble(
  #     replicate = i,  # Use index number
  #     # ... rest of code ...
  #   )
  # })
  
  # Compute confidence intervals
  ci_summary <- genetic_stats %>%
    summarise(across(where(is.numeric), 
                     list(mean = ~mean(., na.rm = TRUE),
                          sd = ~sd(., na.rm = TRUE),
                          ci_lower = ~quantile(., 0.025, na.rm = TRUE),
                          ci_upper = ~quantile(., 0.975, na.rm = TRUE))))
  
  cat("✅ Genetic parameter analysis completed\n")
  
  return(list(
    genetic_stats = genetic_stats,
    confidence_intervals = ci_summary
  ))
}

genetic_parameters <- analyze_genetic_parameters(results_list)


#' 2. Heritability Estimation Simulation
#' Simulate traits and estimate heritability using different relationship matrices
#' 
#' @param results_list List of results from replicate runs
#' @param h2_simulated Vector of simulated heritability values to test
# simulate_heritability_analysis <- function(results_list, h2_simulated = c(0.3, 0.5, 0.7)) {
#   
#   cat("📊 Simulating Heritability Estimation...\n")
#   
#   heritability_results <- map_dfr(results_list, function(replicate, h2) {
#     
#     merged_data <- as.data.table(replicate$merged_data$conservative)
#     pairs <- split_pair_key(merged_data, "pair_key") 
#     merged_data[, id1 := pairs$id1]; merged_data[, id2 := pairs$id2]
#     
#     ids <- unique(c(pairs$id1, pairs$id2))
#     # ids <- unique(c(merged_data$id1, merged_data$id2))
#     n_individuals <- length(ids)
#     
#     # Build relationship matrices
#     true_ibd <- build_sym_matrix(merged_data, "ibd")
#     K_ibd <- build_sym_matrix(merged_data, "IBD")
#     K_ibs <- build_sym_matrix(merged_data, "IBS")
#     K_phylo <- build_sym_matrix(merged_data, "phylo")
#     
#     # Simulate genetic values and phenotypes
#     simulate_trait <- function(K, h2, n_ind = n_individuals) {
#       # Simulate genetic values
#       g <- MASS::mvrnorm(1, mu = rep(0, n_ind), Sigma = K * h2)
#       # Simulate environmental noise
#       e <- rnorm(n_ind, mean = 0, sd = sqrt(1 - h2))
#       # Phenotype
#       y <- g + e
#       return(y)
#     }
#     
#     # Test each heritability value
#     map_dfr(h2_simulated, function(true_h2) {
#       
#       # Simulate trait using IBD matrix as "true" genetic architecture
#       trait <- simulate_trait(true_ibd, true_h2)
#       
#       # Estimate heritability using each relationship matrix
#       estimate_h2 <- function(K, y) {
#         # Simple MLM for heritability estimation
#         n <- length(y)
#         K_scaled <- K / mean(diag(K))
#         
#         # Using EMMA/REML-like approach (simplified)
#         tryCatch({
#           model <- lm(y ~ 1)  # Intercept only
#           var_components <- lme4::VarCorr(lme4::lmer(y ~ 1 + (1|id), 
#                                                      data = data.frame(y = y, id = factor(1:n)),
#                                                      REML = TRUE))
#           vc <- as.data.frame(var_components)
#           h2_est <- vc$vcov[1] / sum(vc$vcov)
#           return(h2_est)
#         }, error = function(e) return(NA))
#       }
#       
#       h2_ibd <- estimate_h2(K_ibd, trait)
#       h2_ibs <- estimate_h2(K_ibs, trait) 
#       h2_phylo <- estimate_h2(K_phylo, trait)
#       
#       tibble(
#         true_heritability = true_h2,
#         estimated_ibd = h2_ibd,
#         estimated_ibs = h2_ibs,
#         estimated_phylo = h2_phylo,
#         bias_ibd = h2_ibd - true_h2,
#         bias_ibs = h2_ibs - true_h2,
#         bias_phylo = h2_phylo - true_h2
#       )
#     })
#   }, h2_simulated, .id = "replicate")
#   
#   # Summarize across replicates
#   h2_summary <- heritability_results %>%
#     group_by(true_heritability) %>%
#     summarise(
#       ibd_mean = mean(estimated_ibd, na.rm = TRUE),
#       ibd_bias = mean(bias_ibd, na.rm = TRUE),
#       ibd_rmse = sqrt(mean(bias_ibd^2, na.rm = TRUE)),
#       ibs_mean = mean(estimated_ibs, na.rm = TRUE),
#       ibs_bias = mean(bias_ibs, na.rm = TRUE),
#       ibs_rmse = sqrt(mean(bias_ibs^2, na.rm = TRUE)),
#       phylo_mean = mean(estimated_phylo, na.rm = TRUE),
#       phylo_bias = mean(bias_phylo, na.rm = TRUE),
#       phylo_rmse = sqrt(mean(bias_phylo^2, na.rm = TRUE))
#     )
#   
#   cat("✅ Heritability simulation completed\n")
#   
#   return(list(
#     detailed_results = heritability_results,
#     summary = h2_summary
#   ))
# }

simulate_heritability_analysis <- function(results_list, h2_simulated = c(0.3, 0.5, 0.7)) {
  
  cat("📊 Simulating Heritability Estimation...\n")
  
  heritability_results <- imap_dfr(results_list, function(replicate, replicate_name) {
    
    # Check if merged_data exists
    if (is.null(replicate$merged_data$conservative)) {
      cat("Warning: No conservative merged data for replicate", replicate_name, "\n")
      return(NULL)
    }
    
    merged_data <- as.data.table(replicate$merged_data$conservative)
    
    # Check if we have enough data
    if (nrow(merged_data) < 2) {
      cat("Warning: Insufficient data for replicate", replicate_name, "\n")
      return(NULL)
    }
    
    # Extract pairs
    if ("pair_key" %in% names(merged_data)) {
      pairs <- split_pair_key(merged_data, "pair_key")
      merged_data[, id1 := pairs$id1]
      merged_data[, id2 := pairs$id2]
    }
    
    ids <- unique(c(merged_data$id1, merged_data$id2))
    n_individuals <- length(ids)
    
    # Build relationship matrices
    true_ibd <- build_sym_matrix(merged_data, "ibd")
    K_ibd <- build_sym_matrix(merged_data, "IBD")
    K_ibs <- build_sym_matrix(merged_data, "IBS")
    K_phylo <- build_sym_matrix(merged_data, "phylo")
    
    # Function to make matrix positive definite
    make_positive_definite <- function(K, tol = 1e-6) {
      # Add small ridge to diagonal if not PD
      if (!is.matrix(K)) {
        # Convert to matrix if it's not
        if (is.null(dim(K))) {
          # If it's a vector, create a 1x1 matrix
          K <- as.matrix(K)
        } else {
          # Try to coerce to matrix
          K <- as.matrix(K)
        }
      }
      
      # Check if matrix is square
      if (nrow(K) != ncol(K)) {
        # Make it square (could be an issue with your build_sym_matrix function)
        n <- max(nrow(K), ncol(K))
        K_new <- diag(n)
        K_new[1:nrow(K), 1:ncol(K)] <- K
        K <- K_new
      }
      
      # Ensure symmetry
      K <- (K + t(K)) / 2
      
      # Check eigenvalues
      eig_vals <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
      
      if (any(eig_vals < tol)) {
        # Add ridge to make it positive definite
        ridge <- abs(min(eig_vals)) + tol
        K <- K + diag(ridge, nrow(K))
        cat("  Added ridge of", ridge, "to make matrix PD for replicate", replicate_name, "\n")
      }
      
      return(K)
    }
    
    # Ensure matrices are positive definite
    true_ibd <- make_positive_definite(true_ibd)
    K_ibd <- make_positive_definite(K_ibd)
    K_ibs <- make_positive_definite(K_ibs)
    K_phylo <- make_positive_definite(K_phylo)
    
    # Also ensure dimensions match
    if (nrow(true_ibd) != n_individuals) {
      cat("Warning: Matrix dimension mismatch for replicate", replicate_name, "\n")
      return(NULL)
    }
    
    # Simulate genetic values and phenotypes
    simulate_trait <- function(K, h2, n_ind = n_individuals) {
      # Ensure K is the right dimension
      if (nrow(K) != n_ind || ncol(K) != n_ind) {
        # Resize if needed
        K <- K[1:n_ind, 1:n_ind, drop = FALSE]
      }
      
      # Scale the matrix appropriately
      K_scaled <- K * h2
      
      # Add small ridge for numerical stability
      K_scaled <- K_scaled + diag(1e-6, n_ind)
      
      # Simulate genetic values using Cholesky decomposition (more stable)
      tryCatch({
        L <- chol(K_scaled)  # Cholesky decomposition
        g <- t(L) %*% rnorm(n_ind)
      }, error = function(e) {
        # If Cholesky fails, use eigen decomposition
        eig <- eigen(K_scaled, symmetric = TRUE)
        # Keep only positive eigenvalues
        pos_eig <- eig$values[eig$values > 1e-10]
        pos_vec <- eig$vectors[, eig$values > 1e-10]
        g <- pos_vec %*% (sqrt(pos_eig) * rnorm(length(pos_eig)))
      })
      
      # Simulate environmental noise
      e <- rnorm(n_ind, mean = 0, sd = sqrt(1 - h2))
      
      # Phenotype
      y <- as.vector(g) + e
      return(y)
    }
    
    # Test each heritability value
    map_dfr(h2_simulated, function(true_h2) {
      
      # Simulate trait using IBD matrix as "true" genetic architecture
      trait <- simulate_trait(true_ibd, true_h2)
      
      # Estimate heritability using each relationship matrix
      estimate_h2 <- function(K, y) {
        # Ensure K is positive definite for mixed model
        K <- make_positive_definite(K)
        
        n <- length(y)
        
        # Create data for mixed model
        df <- data.frame(
          y = y,
          id = factor(1:n),
          row.names = NULL
        )
        
        # Using lme4 for REML estimation
        tryCatch({
          # Scale the relationship matrix
          K_scaled <- K / mean(diag(K))
          
          # Create sparse matrix if large
          if (n > 1000) {
            K_scaled <- as(K_scaled, "sparseMatrix")
          }
          
          # Fit the model using the relationship matrix
          # Using a simplified approach - in practice you might use:
          # kinship2::gls() or other packages that accept covariance matrices
          
          # Alternative: Use EMMA/GEMMA approach or custom REML
          # For simplicity, using a correlation-based estimator
          
          # Method 1: Simple correlation-based estimator (for demonstration)
          K_cor <- cov2cor(K_scaled)
          # Estimate heritability using Haseman-Elston regression style approach
          # This is a simplified approach
          y_centered <- y - mean(y)
          vy <- var(y_centered)
          vg <- sum(K_cor * outer(y_centered, y_centered)) / (n^2)
          h2_est <- vg / vy
          
          # Ensure h2 is between 0 and 1
          h2_est <- max(0, min(1, h2_est))
          
          return(h2_est)
        }, error = function(e) {
          # cat("  Error in h2 estimation:", e$message, "\n")
          return(NA)
        })
      }
      
      h2_ibd <- estimate_h2(K_ibd, trait)
      h2_ibs <- estimate_h2(K_ibs, trait) 
      h2_phylo <- estimate_h2(K_phylo, trait)
      
      tibble(
        replicate = replicate_name,
        true_heritability = true_h2,
        estimated_ibd = h2_ibd,
        estimated_ibs = h2_ibs,
        estimated_phylo = h2_phylo,
        bias_ibd = h2_ibd - true_h2,
        bias_ibs = h2_ibs - true_h2,
        bias_phylo = h2_phylo - true_h2
      )
    })
  })
  
  # Check if we have results
  if (is.null(heritability_results) || nrow(heritability_results) == 0) {
    warning("No heritability results were generated!")
    return(NULL)
  }
  
  # Summarize across replicates
  h2_summary <- heritability_results %>%
    group_by(true_heritability) %>%
    summarise(
      n_replicates = n(),
      ibd_mean = mean(estimated_ibd, na.rm = TRUE),
      ibd_bias = mean(bias_ibd, na.rm = TRUE),
      ibd_rmse = sqrt(mean(bias_ibd^2, na.rm = TRUE)),
      ibs_mean = mean(estimated_ibs, na.rm = TRUE),
      ibs_bias = mean(bias_ibs, na.rm = TRUE),
      ibs_rmse = sqrt(mean(bias_ibs^2, na.rm = TRUE)),
      phylo_mean = mean(estimated_phylo, na.rm = TRUE),
      phylo_bias = mean(bias_phylo, na.rm = TRUE),
      phylo_rmse = sqrt(mean(bias_phylo^2, na.rm = TRUE))
    )
  
  cat("✅ Heritability simulation completed\n")
  cat("   Processed", n_distinct(heritability_results$replicate), "replicates\n")
  
  return(list(
    detailed_results = heritability_results,
    summary = h2_summary
  ))
}

heritability <- simulate_heritability_analysis(results_list = results_list)

# =======================================================================
visualize_genetic_parameters <- function(analysis_results) {
  
  library(tidyverse)
  library(patchwork)
  library(ggridges)
  library(gt)
  
  genetic_stats <- analysis_results$genetic_stats
  
  # DEBUG: Check column names
  cat("Genetic stats columns:\n")
  print(colnames(genetic_stats))
  
  cat("\nConfidence intervals columns:\n")
  print(colnames(analysis_results$confidence_intervals))
  
  # 1. Correlation Heatmap Across Replicates
  # Fix: Use correct column names - assuming they're like "IBD_truth", "IBS_truth", "Phylo_truth"
  p1 <- genetic_stats %>%
    select(ends_with("_truth")) %>%
    pivot_longer(everything(), names_to = "correlation_type", values_to = "correlation") %>%
    mutate(correlation_type = gsub("_truth", "", correlation_type)) %>%
    ggplot(aes(x = correlation_type, y = correlation, fill = correlation_type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    labs(title = "Matrix Correlation with Ground Truth",
         x = "Matrix Type", y = "Correlation") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 2. Effective Loci Distribution
  p2 <- genetic_stats %>%
    select(ends_with("_effective_loci")) %>%
    pivot_longer(everything(), names_to = "matrix_type", values_to = "effective_loci") %>%
    mutate(matrix_type = gsub("_effective_loci", "", matrix_type)) %>%
    ggplot(aes(x = effective_loci, y = matrix_type, fill = matrix_type)) +
    geom_density_ridges(alpha = 0.7) +
    labs(title = "Effective Number of Loci",
         x = "Effective Loci", y = "Matrix Type") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 3. Mean Relatedness Comparison
  p3 <- genetic_stats %>%
    select(ends_with("_mean_relatedness")) %>%
    pivot_longer(everything(), names_to = "matrix_type", values_to = "mean_relatedness") %>%
    mutate(matrix_type = gsub("_mean_relatedness", "", matrix_type)) %>%
    ggplot(aes(x = matrix_type, y = mean_relatedness, fill = matrix_type)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.9) +
    labs(title = "Mean Relatedness Across Matrices",
         x = "Matrix Type", y = "Mean Relatedness") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 4. Confidence Interval Plot - SIMPLIFIED
  # Based on your column pattern: IBD_truth_mean, IBD_truth_sd, etc.
  # Pattern is: METRIC_TYPE_stat
  
  # First, let's see what we're working with
  cat("\nProcessing confidence intervals...\n")
  
  ci_data <- analysis_results$confidence_intervals %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    # Extract the statistic from the end (mean, sd, ci_lower, ci_upper)
    mutate(
      stat = case_when(
        str_detect(parameter, "_mean$") ~ "mean",
        str_detect(parameter, "_sd$") ~ "sd",
        str_detect(parameter, "_ci_lower$") ~ "ci_lower",
        str_detect(parameter, "_ci_upper$") ~ "ci_upper",
        TRUE ~ NA_character_
      ),
      # Extract metric (everything before the last underscore)
      metric = str_remove(parameter, paste0("_", stat, "$"))
    ) %>%
    filter(!is.na(stat)) %>%  # Remove any rows where stat wasn't identified
    select(metric, stat, value) %>%
    # Now pivot to wide format
    pivot_wider(
      names_from = "stat",
      values_from = "value"
    )
  
  cat("Processed CI data structure:\n")
  print(str(ci_data))
  cat("\nCI data preview:\n")
  print(ci_data)
  
  # Check if we have the required columns
  required_cols <- c("mean", "ci_lower", "ci_upper")
  if (!all(required_cols %in% colnames(ci_data))) {
    warning("Missing required columns for CI plot. Available columns: ", 
            paste(colnames(ci_data), collapse = ", "))
    
    # Create a placeholder plot
    p4 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "CI data not in expected format", size = 6) +
      theme_void()
  } else {
    # Create the CI plot
    p4 <- ggplot(ci_data, aes(x = metric, y = mean, color = metric)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
      labs(title = "Parameter Estimates with 95% CI",
           x = "Parameter", y = "Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  }
  
  # 5. Additional diagnostic plot: Show all CI data
  ci_long <- analysis_results$confidence_intervals %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    separate(parameter, into = c("metric", "type", "stat"), sep = "_", fill = "right") %>%
    unite("metric_type", metric, type, sep = "_", remove = FALSE)
  
  p5 <- ggplot(ci_long, aes(x = stat, y = value, fill = metric)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    labs(title = "All Confidence Interval Statistics",
         x = "Statistic", y = "Value") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine plots - updated layout
  combined_plot <- (p1 | p2) / (p3 | p4)
  
  print(combined_plot)
  print(p5)
  
  # Create a summary table
  if (exists("ci_data") && nrow(ci_data) > 0) {
    summary_table <- ci_data %>%
      gt() %>%
      tab_header(
        title = "Genetic Parameter Confidence Intervals"
      ) %>%
      fmt_number(columns = where(is.numeric), decimals = 3) %>%
      cols_label(
        metric = "Parameter",
        mean = "Mean",
        sd = "SD",
        ci_lower = "CI Lower (2.5%)",
        ci_upper = "CI Upper (97.5%)"
      )
    
    print(summary_table)
  }
  
  # Return individual plots for saving
  return(list(
    correlation_plot = p1,
    effective_loci_plot = p2,
    relatedness_plot = p3,
    ci_plot = p4,
    diagnostic_plot = p5,
    combined_plot = combined_plot,
    ci_data = if(exists("ci_data")) ci_data else NULL,
    summary_table = if(exists("summary_table")) summary_table else NULL
  ))
}


visualize_heritability_results <- function(heritability_results) {
  
  library(tidyverse)
  library(patchwork)
  library(ggh4x)
  library(ggpubr)
  
  detailed <- heritability_results$detailed_results
  summary <- heritability_results$summary
  
  # 1. Bias Plot by True Heritability
  p1 <- detailed %>%
    pivot_longer(cols = starts_with("bias_"), 
                 names_to = "matrix_type", 
                 values_to = "bias") %>%
    mutate(matrix_type = gsub("bias_", "", matrix_type)) %>%
    ggplot(aes(x = factor(true_heritability), y = bias, fill = matrix_type)) +
    geom_boxplot(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Bias in Heritability Estimation",
         x = "True Heritability", y = "Bias (Estimated - True)",
         fill = "Matrix Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 2. Scatter Plot: True vs Estimated
  p2 <- detailed %>%
    pivot_longer(cols = starts_with("estimated_"), 
                 names_to = "matrix_type", 
                 values_to = "estimated") %>%
    mutate(matrix_type = gsub("estimated_", "", matrix_type)) %>%
    ggplot(aes(x = true_heritability, y = estimated, color = matrix_type)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ matrix_type) +
    labs(title = "True vs Estimated Heritability",
         x = "True Heritability", y = "Estimated Heritability",
         color = "Matrix Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 3. RMSE Comparison
  p3 <- summary %>%
    pivot_longer(cols = ends_with("_rmse"), 
                 names_to = "matrix_type", 
                 values_to = "rmse") %>%
    mutate(matrix_type = gsub("_rmse", "", matrix_type)) %>%
    ggplot(aes(x = factor(true_heritability), y = rmse, fill = matrix_type)) +
    geom_col(position = position_dodge()) +
    labs(title = "Root Mean Square Error (RMSE) by Method",
         x = "True Heritability", y = "RMSE",
         fill = "Matrix Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 4. Detailed Performance Summary Plot
  perf_data <- summary %>%
    pivot_longer(cols = c(ends_with("_bias"), ends_with("_rmse")),
                 names_to = "metric", values_to = "value") %>%
    separate(metric, into = c("matrix_type", "metric"), sep = "_") %>%
    pivot_wider(names_from = "metric", values_from = "value")
  
  p4 <- ggplot(perf_data, aes(x = factor(true_heritability), y = bias, 
                              size = rmse, color = matrix_type)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = "Performance Summary: Bias (color) and RMSE (size)",
         x = "True Heritability", y = "Bias",
         color = "Matrix Type", size = "RMSE") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 5. Distribution of Estimates (Ridge Plot)
  p5 <- detailed %>%
    pivot_longer(cols = starts_with("estimated_"), 
                 names_to = "matrix_type", 
                 values_to = "estimated") %>%
    mutate(matrix_type = gsub("estimated_", "", matrix_type)) %>%
    ggplot(aes(x = estimated, y = factor(true_heritability), fill = matrix_type)) +
    geom_density_ridges(alpha = 0.6, scale = 0.9) +
    facet_wrap(~ matrix_type, ncol = 1) +
    geom_vline(xintercept = unique(detailed$true_heritability), 
               linetype = "dashed", alpha = 0.5) +
    labs(title = "Distribution of Heritability Estimates",
         x = "Estimated Heritability", y = "True Heritability",
         fill = "Matrix Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Combine plots
  main_plots <- (p1 | p2) / (p3 | p4)
  
  print(main_plots)
  print(p5)
  
  # Create summary table
  summary_table <- summary %>%
    gt::gt() %>%
    gt::tab_header(
      title = "Heritability Simulation Summary"
    ) %>%
    gt::fmt_number(columns = where(is.numeric), decimals = 3)
  
  # Print summary table
  print(summary_table)
  
  return(list(
    bias_plot = p1,
    scatter_plot = p2,
    rmse_plot = p3,
    performance_plot = p4,
    distribution_plot = p5,
    summary_table = summary_table,
    main_layout = main_plots
  ))
}


# Combined Dashboard Function
create_analysis_dashboard <- function(genetic_analysis, heritability_analysis) {
  
  library(patchwork)
  library(tidyverse)
  library(gridExtra)
  library(cowplot)
  
  # Create individual visualizations
  genetic_plots <- visualize_genetic_parameters(genetic_analysis)
  heritability_plots <- visualize_heritability_results(heritability_analysis)
  
  # 1. Create a comprehensive dashboard
  dashboard <- plot_grid(
    genetic_plots$combined_plot,
    heritability_plots$main_layout,
    ncol = 1,
    labels = c("A. Genetic Parameter Analysis", "B. Heritability Simulation"),
    label_size = 14
  )
  
  print(dashboard)
  
  # 2. Create performance comparison across methods
  if (!is.null(heritability_analysis)) {
    # Extract key metrics for comparison
    performance_comparison <- heritability_analysis$detailed_results %>%
      pivot_longer(cols = starts_with("estimated_"), 
                   names_to = "method", 
                   values_to = "estimate") %>%
      mutate(method = gsub("estimated_", "", method)) %>%
      group_by(true_heritability, method) %>%
      summarise(
        mean_estimate = mean(estimate, na.rm = TRUE),
        sd_estimate = sd(estimate, na.rm = TRUE),
        bias = mean(estimate - true_heritability, na.rm = TRUE),
        rmse = sqrt(mean((estimate - true_heritability)^2, na.rm = TRUE)),
        .groups = "drop"
      )
    
    perf_plot <- ggplot(performance_comparison, 
                        aes(x = factor(true_heritability), y = bias, 
                            color = method, group = method)) +
      geom_line(size = 1.2, alpha = 0.7) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      facet_wrap(~ method, nrow = 1) +
      labs(title = "Method Performance Comparison",
           x = "True Heritability", y = "Bias",
           color = "Method") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    print(perf_plot)
  }
  
  # 3. Create summary report
  cat("\n" + strrep("=", 80) + "\n")
  cat("ANALYSIS SUMMARY REPORT\n")
  cat(strrep("=", 80) + "\n\n")
  
  cat("GENETIC PARAMETERS:\n")
  cat("- Number of replicates:", nrow(genetic_analysis$genetic_stats), "\n")
  cat("- Parameters analyzed:", 
      paste(names(genetic_analysis$genetic_stats)[-1], collapse = ", "), "\n\n")
  
  if (!is.null(heritability_analysis)) {
    cat("HERITABILITY SIMULATION:\n")
    cat("- Number of replicates:", 
        n_distinct(heritability_analysis$detailed_results$replicate), "\n")
    cat("- True h² values tested:", 
        paste(unique(heritability_analysis$detailed_results$true_heritability), 
              collapse = ", "), "\n")
    
    # Best performing method
    best_method <- heritability_analysis$summary %>%
      filter(true_heritability == median(true_heritability)) %>%
      pivot_longer(cols = ends_with("_rmse"), 
                   names_to = "method", 
                   values_to = "rmse") %>%
      slice_min(rmse) %>%
      pull(method)
    
    cat("- Best performing method (lowest RMSE at median h²):", 
        gsub("_rmse", "", best_method), "\n")
  }
  
  cat(strrep("=", 80) + "\n")
  
  # Return all visualization objects
  return(list(
    dashboard = dashboard,
    genetic_visualizations = genetic_plots,
    heritability_visualizations = heritability_plots,
    performance_comparison_plot = if(exists("perf_plot")) perf_plot else NULL
  ))
}

# 4. Usage Example
# Run analyses
genetic_results <- analyze_genetic_parameters(results_list)
heritability_results <- simulate_heritability_analysis(results_list)

# Create visualizations
vis1 <- visualize_genetic_parameters(genetic_parameters)
vis2 <- visualize_heritability_results(heritability)

# Save individual plots
ggsave("genetic_parameters.png", vis1$combined_plot, width = 12, height = 10)
ggsave("heritability_results.png", vis2$main_layout, width = 12, height = 10)

# Create comprehensive dashboard
# dashboard <- create_analysis_dashboard(genetic_results, heritability_results)
dashboard <- create_analysis_dashboard(genetic_parameters, heritability)
ggsave("analysis_dashboard.png", dashboard$dashboard, width = 16, height = 14)

# Create interactive HTML report (optional)
if (requireNamespace("plotly", quietly = TRUE)) {
  library(plotly)
  
  interactive_plot <- ggplotly(vis2$scatter_plot)
  htmlwidgets::saveWidget(interactive_plot, "heritability_interactive.html")
}

# 5. Quick Diagnostic Plot Function
quick_diagnostic_plots <- function(results_list, n_replicates = 3) {
  
  # Set up plotting area
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Calculate layout based on number of replicates
  n_reps <- min(n_replicates, length(results_list))
  n_cols <- min(3, n_reps)
  n_rows <- ceiling(n_reps / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 2))
  
  for (i in 1:n_reps) {
    replicate <- results_list[[i]]
    
    # Check if replicate has the expected structure
    if (is.null(replicate$merged_data$conservative)) {
      plot(0, 0, type = "n", xlab = "", ylab = "", main = paste("Replicate", i))
      text(0, 0, "No merged_data$conservative", col = "red")
      next
    }
    
    merged_data <- replicate$merged_data$conservative
    
    # Check if merged_data is a data frame with rows
    if (!is.data.frame(merged_data) && !is.data.table(merged_data)) {
      plot(0, 0, type = "n", xlab = "", ylab = "", main = paste("Replicate", i))
      text(0, 0, "Invalid data format", col = "red")
      next
    }
    
    if (nrow(merged_data) == 0) {
      plot(0, 0, type = "n", xlab = "", ylab = "", main = paste("Replicate", i))
      text(0, 0, "Empty data frame", col = "orange")
      next
    }
    
    # Plot 1: IBD distribution (with proper NA handling)
    if ("IBD" %in% names(merged_data)) {
      ibd_vals <- merged_data$IBD
      # Remove infinite and NA values for plotting
      ibd_vals_clean <- ibd_vals[is.finite(ibd_vals) & !is.na(ibd_vals)]
      
      if (length(ibd_vals_clean) > 0) {
        # Use tryCatch to handle any plotting errors
        tryCatch({
          hist(ibd_vals_clean, 
               main = paste("Rep", i, "- IBD (n =", length(ibd_vals_clean), ")"),
               xlab = "IBD", 
               col = "lightblue",
               breaks = 30,
               xlim = range(ibd_vals_clean, na.rm = TRUE))
          
          # Add summary stats as text
          stats_text <- paste(
            "Mean:", round(mean(ibd_vals_clean), 3),
            "\nSD:", round(sd(ibd_vals_clean), 3)
          )
          mtext(stats_text, side = 3, line = -2, cex = 0.7)
        }, error = function(e) {
          plot(0, 0, type = "n", xlab = "", ylab = "", 
               main = paste("Rep", i, "- IBD"))
          text(0, 0, paste("Plot error:", e$message), col = "red", cex = 0.8)
        })
      } else {
        plot(0, 0, type = "n", xlab = "", ylab = "", 
             main = paste("Rep", i, "- IBD"))
        text(0, 0, "No valid IBD values", col = "orange")
      }
    } else {
      plot(0, 0, type = "n", xlab = "", ylab = "", 
           main = paste("Rep", i, "- IBD"))
      text(0, 0, "No IBD column", col = "gray")
    }
  }
  
  par(mfrow = c(1, 1))
  
  # Print summary information
  cat("\n", strrep("=", 60), "\n")
  cat("DIAGNOSTIC SUMMARY\n")
  cat(strrep("=", 60), "\n")
  
  for (i in 1:n_reps) {
    replicate <- results_list[[i]]
    
    if (is.null(replicate$merged_data$conservative)) {
      cat(sprintf("Replicate %d: No conservative merged_data\n", i))
      next
    }
    
    merged_data <- replicate$merged_data$conservative
    
    if (nrow(merged_data) == 0) {
      cat(sprintf("Replicate %d: Empty data frame (0 rows)\n", i))
      next
    }
    
    # Summary statistics
    cat(sprintf("\nReplicate %d:\n", i))
    cat(sprintf("  Rows: %d, Columns: %d\n", nrow(merged_data), ncol(merged_data)))
    
    # Check key columns
    key_cols <- c("IBD", "IBS", "phylo", "pair_key")
    for (col in key_cols) {
      if (col %in% names(merged_data)) {
        n_missing <- sum(is.na(merged_data[[col]]))
        n_infinite <- sum(!is.finite(merged_data[[col]]))
        cat(sprintf("  %s: NAs = %d, Infinite = %d\n", 
                    col, n_missing, n_infinite))
      }
    }
  }
  cat(strrep("=", 60), "\n")
}


quick_diagnostic_plots(results_list = results_list, n_replicates = 10)
