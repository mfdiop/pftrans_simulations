

#' 3. Population Structure Analysis
#' Analyze population structure confounding
#' 
#' @param results_list List of results from replicate runs
#' @param genotype_data Optional genotype matrix for PCA
analyze_population_structure <- function(results_list, genotype_data = NULL) {
  
  cat("🏛️ Analyzing Population Structure Effects...\n")
  
  structure_results <- imap_dfr(results_list, function(replicate, replicate_name) {
    
    # Get the merged data (using conservative merge)
    merged_data <- replicate$merged_data$conservative
    
    # Print progress - use replicate_name (string) not replicate (list)
    cat("Merged data with ", nrow(merged_data), "pairs for replicate ", replicate_name, "\n")
    
    # Method 1: Compare relationship estimates with geographic distance
    # (Assuming you have some spatial coordinates or population labels)
    
    # Method 2: Perform Mantel tests between relationship matrices
    true_ibd <- build_sym_matrix(merged_data, "ibd")
    K_ibd <- build_sym_matrix(merged_data, "IBD")
    K_ibs <- build_sym_matrix(merged_data, "IBS")
    K_phylo <- build_sym_matrix(merged_data, "phylo")
    
    # Create a "population structure" matrix (example: based on first few PCs)
    if (!is.null(genotype_data)) {
      # Perform PCA on genotypes
      pca <- prcomp(genotype_data[ids, ], scale. = TRUE)
      # Use first PC as population structure proxy
      pc_dist <- dist(pca$x[, 1:2])  # Distance in PC space
      pc_matrix <- as.matrix(pc_dist)
      
      # Mantel tests between relationship matrices and population structure
      mantel_ibd <- vegan::mantel(as.dist(K_ibd), as.dist(pc_matrix), permutations = 999)
      mantel_ibs <- vegan::mantel(as.dist(K_ibs), as.dist(pc_matrix), permutations = 999)
      mantel_phylo <- vegan::mantel(as.dist(K_phylo), as.dist(pc_matrix), permutations = 999)
      
      mantel_results <- c(
        IBD_pop_cor = mantel_ibd$statistic,
        IBS_pop_cor = mantel_ibs$statistic, 
        Phylo_pop_cor = mantel_phylo$statistic
      )
    } else {
      mantel_results <- c(IBD_pop_cor = NA, IBS_pop_cor = NA, Phylo_pop_cor = NA)
    }
    
    # Method 3: Analyze clustering patterns
    # Perform hierarchical clustering on each relationship matrix
    cluster_analysis <- function(mat, method = "average") {
      hc <- hclust(as.dist(1 - mat), method = method)
      # Compute cophenetic correlation
      coph <- cor(as.dist(1 - mat), cophenetic(hc))
      return(coph)
    }
    
    coph_ibd <- cluster_analysis(K_ibd)
    coph_ibs <- cluster_analysis(K_ibs)
    coph_phylo <- cluster_analysis(K_phylo)
    
    tibble(
      replicate = replicate_name,
      IBD_population_correlation = mantel_results["IBD_pop_cor"],
      IBS_population_correlation = mantel_results["IBS_pop_cor"],
      Phylo_population_correlation = mantel_results["Phylo_pop_cor"],
      IBD_cophenetic_correlation = coph_ibd,
      IBS_cophenetic_correlation = coph_ibs,
      Phylo_cophenetic_correlation = coph_phylo
    )
  })
  
  cat("✅ Population structure analysis completed\n")
  
  return(structure_results)
}


analyze_population_structure(results_list)

#' 4. Comprehensive Replicate Analysis Framework
#' Master function to analyze all 20 replicates comprehensively
#' 
#' @param results_list List of results from 20 replicate runs
run_comprehensive_replicate_analysis <- function(results_list) {
  
  cat("🎯 Starting Comprehensive Analysis of", length(results_list), "Replicates...\n\n")
  
  # 1. Performance Metrics Consolidation
  cat("1. Consolidating Performance Metrics...\n")
  performance_metrics <- imap_dfr(results_list, function(repl, replicate_name) {
    
    # Print progress - use replicate_name (string) not replicate (list)
    cat(" Compute analysis for ", replicate_name, "\n")
    
    # Extract key metrics from each replicate
    perf_table <- repl$performance_summary
    
    # Add additional metrics from evaluations
    optimal_eval <- repl$evaluations$optimal_youden$metrics_summary
    
    tibble(
      replicate = replicate_name,
      IBD_AUC = perf_table$auc[perf_table$predictor == "IBD"],
      IBS_AUC = perf_table$auc[perf_table$predictor == "IBS"],
      Phylo_AUC = perf_table$auc[perf_table$predictor == "phylo"],
      IBD_AUPR = perf_table$aupr[perf_table$predictor == "IBD"],
      IBS_AUPR = perf_table$aupr[perf_table$predictor == "IBS"], 
      Phylo_AUPR = perf_table$aupr[perf_table$predictor == "phylo"]
    )
  })
  
  # 2. Statistical Testing Across Replicates
  cat("2. Performing Statistical Comparisons...\n")
  
  # Paired t-tests between methods
  performance_comparisons <- list(
    IBD_vs_IBS_AUC = t.test(performance_metrics$IBD_AUC, performance_metrics$IBS_AUC, paired = TRUE),
    IBD_vs_Phylo_AUC = t.test(performance_metrics$IBD_AUC, performance_metrics$Phylo_AUC, paired = TRUE),
    IBS_vs_Phylo_AUC = t.test(performance_metrics$IBS_AUC, performance_metrics$Phylo_AUC, paired = TRUE),
    
    IBD_vs_IBS_AUPR = t.test(performance_metrics$IBD_AUPR, performance_metrics$IBS_AUPR, paired = TRUE),
    IBD_vs_Phylo_AUPR = t.test(performance_metrics$IBD_AUPR, performance_metrics$Phylo_AUPR, paired = TRUE),
    IBS_vs_Phylo_AUPR = t.test(performance_metrics$IBS_AUPR, performance_metrics$Phylo_AUPR, paired = TRUE)
  )
  
  # 3. Effect Size Calculations
  cat("3. Calculating Effect Sizes...\n")
  effect_sizes <- performance_metrics %>%
    summarise(
      IBD_IBS_AUC_CohenD = (mean(IBD_AUC) - mean(IBS_AUC)) / sqrt((sd(IBD_AUC)^2 + sd(IBS_AUC)^2)/2),
      IBD_Phylo_AUC_CohenD = (mean(IBD_AUC) - mean(Phylo_AUC)) / sqrt((sd(IBD_AUC)^2 + sd(Phylo_AUC)^2)/2),
      IBS_Phylo_AUC_CohenD = (mean(IBS_AUC) - mean(Phylo_AUC)) / sqrt((sd(IBS_AUC)^2 + sd(Phylo_AUC)^2)/2)
    )
  
  # 4. Bootstrap Confidence Intervals
  cat("4. Computing Bootstrap Confidence Intervals...\n")
  bootstrap_cis <- map_dfr(c("IBD_AUC", "IBS_AUC", "Phylo_AUC", "IBD_AUPR", "IBS_AUPR", "Phylo_AUPR"), 
                           function(metric) {
                             values <- performance_metrics[[metric]]
                             boot_results <- boot::boot(values, function(x, i) mean(x[i]), R = 1000)
                             ci <- boot::boot.ci(boot_results, type = "bca")
                             
                             tibble(
                               metric = metric,
                               mean = mean(values),
                               sd = sd(values),
                               ci_lower = ci$bca[4],
                               ci_upper = ci$bca[5]
                             )
                           })
  
  # 5. Variance Component Analysis
  cat("5. Analyzing Variance Components...\n")
  # Fit linear model to understand sources of variation
  perf_long <- performance_metrics %>%
    pivot_longer(cols = c(IBD_AUC, IBS_AUC, Phylo_AUC, IBD_AUPR, IBS_AUPR, Phylo_AUPR),
                 names_to = "metric_method", values_to = "value") %>%
    separate(metric_method, into = c("method", "metric"), sep = "_")
  
  variance_components <- lme4::lmer(value ~ method + metric + (1|repl), data = perf_long)
  
  # 6. Rank Analysis
  cat("6. Performing Rank Analysis...\n")
  rank_analysis <- performance_metrics %>%
    mutate(
      AUC_rank_IBD = rank(-IBD_AUC),
      AUC_rank_IBS = rank(-IBS_AUC),
      AUC_rank_Phylo = rank(-Phylo_AUC),
      AUPR_rank_IBD = rank(-IBD_AUPR),
      AUPR_rank_IBS = rank(-IBS_AUPR),
      AUPR_rank_Phylo = rank(-Phylo_AUPR)
    ) %>%
    summarise(
      IBD_AUC_median_rank = median(AUC_rank_IBD),
      IBS_AUC_median_rank = median(AUC_rank_IBS),
      Phylo_AUC_median_rank = median(AUC_rank_Phylo),
      IBD_wins_AUC = sum(AUC_rank_IBD == 1) / n(),
      IBS_wins_AUC = sum(AUC_rank_IBS == 1) / n(),
      Phylo_wins_AUC = sum(AUC_rank_Phylo == 1) / n()
    )
  
  cat("✅ Comprehensive analysis completed!\n")
  
  return(list(
    performance_metrics = performance_metrics,
    statistical_tests = performance_comparisons,
    effect_sizes = effect_sizes,
    bootstrap_cis = bootstrap_cis,
    variance_components = summary(variance_components),
    rank_analysis = rank_analysis,
    performance_summary = list(
      AUC = performance_metrics %>% summarise(
        IBD = paste0(round(mean(IBD_AUC), 3), " (", round(sd(IBD_AUC), 3), ")"),
        IBS = paste0(round(mean(IBS_AUC), 3), " (", round(sd(IBS_AUC), 3), ")"),
        Phylo = paste0(round(mean(Phylo_AUC), 3), " (", round(sd(Phylo_AUC), 3), ")")
      ),
      AUPR = performance_metrics %>% summarise(
        IBD = paste0(round(mean(IBD_AUPR), 3), " (", round(sd(IBD_AUPR), 3), ")"),
        IBS = paste0(round(mean(IBS_AUPR), 3), " (", round(sd(IBS_AUPR), 3), ")"),
        Phylo = paste0(round(mean(Phylo_AUPR), 3), " (", round(sd(Phylo_AUPR), 3), ")")
      )
    )
  ))
}


# 5. Additional Robust Analysis Suggestions
# Sensitivity Analysis for Key Parameters
run_sensitivity_analysis <- function() {
  # Test different generational cutoffs
  cutoffs <- c(3, 5, 10, 15, 25)
  # Test different relationship thresholds
  thresholds <- seq(0.1, 0.9, 0.1)
  # Test different phylogenetic transformation parameters
  sigma_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)
}

# Power Analysis
calculate_statistical_power <- function(results_list, alpha = 0.05) {
  # Calculate power to detect differences between methods
  # Based on your 20 replicates, determine if you have sufficient power
  # to detect meaningful differences in performance metrics
}

# Meta-analysis across replicates
perform_meta_analysis <- function(results_list) {
  # Use meta-analysis techniques to combine results across replicates
  # Calculate overall effect sizes with random effects models
}


# Run comprehensive analysis on your 20 replicates
comprehensive_results <- run_comprehensive_replicate_analysis(results_list)

# Print key findings
cat("\n📊 KEY FINDINGS ACROSS 20 REPLICATES:\n")
cat("AUC (mean ± sd):\n")
print(comprehensive_results$performance_summary$AUC)
cat("\nAUPR (mean ± sd):\n") 
print(comprehensive_results$performance_summary$AUPR)
cat("\nEffect Sizes (Cohen's d):\n")
print(comprehensive_results$effect_sizes)
cat("\nRank Analysis:\n")
print(comprehensive_results$rank_analysis)

# Generate publication-ready tables and figures
generate_manuscript_tables(comprehensive_results)