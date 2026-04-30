

# =============================================================================
# IBD Analysis Pipeline - Relationship Inference Evaluation
# Description: This script analyzes and compares different relationship inference
#              methods (IBD, IBS, phylogenetic distance) against ground truth data.
# Author: Mouhamadou Fadel DIOP
# Date: [Date]
#' Run Full Relationship Inference Pipeline
#' 
#' This function executes a complete pipeline for comparing relationship inference methods
#' including IBD, IBS, and phylogenetic approaches against ground truth data.
#'
#' @param ibd_file   Path to IBD inference results file (TSV format)
#' @param ibs_file   Path to IBS matrix file (RDS format)  
#' @param treefile   Path to phylogenetic tree file
#' @param truth_file Path to ground truth relationships file (TSV format)
#' @param method     Phylogenetic distance method ("patristic": default, "linear" or "gaussian")
#' @param gen_cutoff Generational cutoff for defining close relatives (default: 5)
#' @param outdir     Output directory for results (default: "out")
#'
#' @return List containing analysis results and saves output files to disk
#' 
#' @examples
#' \dontrun{
#' results <- run_full_pipeline(
#'   ibd_file = "inferred/replicate1/inferred_ibd_hmm.tsv",
#'   ibs_file = "inferred/replicate1/ibs_matrix.rds", 
#'   treefile = "phylo_output/replicate1/run1_chr1_modelfinder.treefile",
#'   truth_file = "inferred/replicate1/true_ibd_summary.tsv"
#' )
#' }
#' =============================================================================
run_full_pipeline <- function(ibd_file, ibs_file, treefile, truth_file, 
                              method = "patristic", gen_cutoff = 5, 
                              outdir = "out") {
  
  # ============================================================================
  # INITIALIZATION AND SETUP
  # ============================================================================
  
  cat("🚀 Starting full pipeline execution...\n")
  cat("📁 Output directory:", outdir, "\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    cat("📂 Creating output directory:", outdir, "\n")
    dir.create(outdir, recursive = TRUE)
  }
  
  # Load required packages with error handling
  cat("📦 Loading required packages...\n")
  required_packages <- c("tidyverse", "hexbin", "gridExtra", "PRROC", "pROC", 
                         "vegan", "scales", "ape", "phangorn", "ggpubr", "data.table")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop("❌ Package ", pkg, " is not installed. Please install it first.")
    }
  }
  cat("✅ All packages loaded successfully\n")
  
  # ============================================================================
  # DATA LOADING AND PREPROCESSING
  # ============================================================================
  
  cat("\n📥 Loading input data...\n")
  
  # Load IBD data
  cat("   Reading IBD file:", ibd_file, "\n")
  if (!file.exists(ibd_file)) stop("❌ IBD file not found: ", ibd_file)
  ibd_df <- read_tsv(ibd_file, show_col_types = FALSE)
  cat("   ✅ IBD data loaded:", nrow(ibd_df), "rows\n")
  
  # Load and process IBS data
  cat("   Reading IBS file:", ibs_file, "\n")
  if (!file.exists(ibs_file)) stop("❌ IBS file not found: ", ibs_file)
  ibs_matrix <- readRDS(ibs_file)
  cat("   ✅ IBS matrix loaded:", nrow(ibs_matrix), "x", ncol(ibs_matrix), "matrix\n")
  
  # Convert IBS matrix to long format
  cat("   Converting IBS matrix to long format...\n")
  ibs_df <- matrix_to_long_efficient(ibs_matrix, value_name = "IBS") %>% 
    arrange(id1)
  cat("   ✅ IBS data converted:", nrow(ibs_df), "unique pairs\n")
  
  # Save IBS long format for future use
  ibs_output_file <- file.path(dirname(ibs_file), "inferred_ibs.tsv")
  cat("   Saving IBS long format to:", ibs_output_file, "\n")
  write_tsv(ibs_df, ibs_output_file)
  
  # Load and process ground truth data
  cat("   Reading ground truth file:", truth_file, "\n")
  if (!file.exists(truth_file)) stop("❌ Truth file not found: ", truth_file)
  true_links <- read_tsv(truth_file, show_col_types = FALSE) %>% 
    mutate(
      Id1 = paste0("tsk_", Id1), 
      Id2 = paste0("tsk_", Id2), 
      ibd = round(total_ibd_prop, 3)
    ) %>% 
    select(Id1, Id2, ibd, min_tmrca) %>% 
    rename(generations = min_tmrca)
  cat("   ✅ Ground truth data loaded:", nrow(true_links), "rows\n")
  
  # Calculate phylogenetic distances
  cat("   Calculating phylogenetic distances...\n")
  if (!file.exists(treefile)) stop("❌ Tree file not found: ", treefile)
  phylo_df <- phylo_long(treefile, method = "patristic")
  cat("   ✅ Phylogenetic distances calculated:", nrow(phylo_df), "pairs\n")
  
  # Define input file paths for a single simulation replicate
  # truth_file <- "simulations/single_run/inferred/replicate1/true_ibd_summary.tsv"    # Ground truth IBD data
  # ibd_file <- "simulations/single_run/inferred/replicate1/inferred_ibd_hmm.tsv"      # Inferred IBD segments
  # ibs_file <- "simulations/single_run/inferred/replicate1/ibs_matrix.rds"          # Inferred IBS segments
  # treefile <- "simulations/single_run/phylo_output/replicate1/run1_chr1_modelfinder.treefile"  # Phylogenetic tree
  
  # ============================================================================
  # DATA INTEGRATION
  # ============================================================================
  
  cat("\n🔗 Integrating datasets...\n")
  
  # Define ground truth relationships
  cat("   Defining ground truth with generational cutoff:", gen_cutoff, "\n")
  df_ground <- define_ground_truth(true_links, gen_cutoff = gen_cutoff) %>% 
    rename(id1 = Id1, id2 = Id2)
  cat("   ✅ Ground truth defined:", sum(df_ground$true_link), "true links out of", nrow(df_ground), "pairs\n")
  
  # Merge datasets using different strategies
  cat("   Merging datasets using different strategies...\n")
  
  # Conservative merge (only pairs with ground truth)
  cat("   - Conservative merge (ground truth pairs only)...\n")
  merged_conservative <- merge_metrics_with(
    df_ground = df_ground,
    ibd_df = ibd_df, 
    ibs_df = ibs_df, 
    phylo_df = phylo_df,
    merge_method = "ground_truth"
  )
  cat("     ✅ Conservative merge:", nrow(merged_conservative), "pairs\n")
  
  # Strict merge (intersection of all methods)
  cat("   - Strict merge (intersection of all methods)...\n")
  merged_strict <- merge_metrics_with(
    df_ground = df_ground,
    ibd_df = ibd_df,
    ibs_df = ibs_df,
    phylo_df = phylo_df, 
    merge_method = "intersection"
  )
  cat("     ✅ Strict merge:", nrow(merged_strict), "pairs\n")
  
  # Union merge (all pairs from any method)
  cat("   - Union merge (all pairs from any method)...\n")
  merged_union <- merge_metrics_with(
    df_ground = df_ground,
    ibd_df = ibd_df,
    ibs_df = ibs_df,
    phylo_df = phylo_df, 
    merge_method = "union"
  )
  cat("     ✅ Union merge:", nrow(merged_union), "pairs\n")
  
  # ============================================================================
  # VISUALIZATION
  # ============================================================================
  
  cat("\n📊 Generating visualizations...\n")
  
  # Generate base filename from input directory
  # base_name <- basename(dirname(ibs_file))
  # 
  # # Scatter plots
  # scatter_output <- file.path(outdir, paste0(base_name, "_metric_comparisons.png"))
  # cat("   Creating scatter plots:", scatter_output, "\n")
  # plot_scatter(merged_conservative, scatter_output)
  # cat("   ✅ Scatter plots saved\n")
  # 
  # # Hexbin plots
  # hex_output <- file.path(outdir, paste0(base_name, "_ibs_vs_ibd.png"))
  # cat("   Creating hexbin plots:", hex_output, "\n")
  # plot_hex(merged_conservative, hex_output)
  # cat("   ✅ Hexbin plots saved\n")
  
  # ============================================================================
  # STATISTICAL ANALYSIS
  # ============================================================================
  
  cat("\n📈 Performing statistical analysis...\n")
  
  # Correlation analysis
  cat("   Calculating correlations between IBD and IBS...\n")
  corr_res <- compute_correlations_v1(merged_conservative, "IBD", "IBS", transform_for_pearson = "logit")
  cat("   ✅ Correlation analysis completed:\n")
  cat("     - Pearson correlation:", round(corr_res$pearson$estimate, 4), "\n")
  cat("     - Spearman correlation:", round(corr_res$spearman$estimate, 4), "\n")
  
  # Save correlation results
  # corr_output <- file.path(outdir, paste0(base_name, "_corr_coef.rds"))
  # cat("   Saving correlation results:", corr_output, "\n")
  # saveRDS(corr_res, file = corr_output)
  # 
  # Prepare data for Mantel test
  cat("   Preparing matrices for Mantel test...\n")
  merged_metrics <- merged_conservative %>% 
    separate(pair_key, into = c("id1", "id2"), sep = "(?<=[0-9])_", extra = "merge")
  
  ids <- unique(c(merged_metrics$id1, merged_metrics$id2))
  cat("   Unique individuals for Mantel test:", length(ids), "\n")
  
  # Build symmetric matrices
  mat_ibd <- build_sym_matrix(merged_metrics, "IBD", ids)
  mat_ibs <- build_sym_matrix(merged_metrics, "IBS", ids)
  mat_phylo <- build_sym_matrix(merged_metrics, "phylo", ids)
  cat("   ✅ Symmetric matrices built for Mantel test\n")
  
  # # used mat_phylo instead of 1-mat_phylo because we have distance not similarity
  # vegan::mantel(as.dist(1-mat_ibd), as.dist(mat_phylo), permutations = 999, method="pearson") 
  # vegan::mantel(as.dist(1-mat_ibd), as.dist(1-mat_ibs), permutations = 999, method="pearson")
  # vegan::mantel(as.dist(1-mat_ibs), as.dist(mat_phylo), permutations = 999, method="pearson")
  
  # ============================================================================
  # MODEL EVALUATION
  # ============================================================================
  
  cat("\n🎯 Evaluating model performance...\n")
  
  # Evaluation with different threshold methods
  threshold_methods <- list(
    optimal_youden = list(method = "optimal_youden", fixed_thresholds = NULL),
    fixed = list(method = "fixed", fixed_thresholds = c(IBD = 0.125, IBS = 0.8, phylo = 0.5)),
    median = list(method = "median", fixed_thresholds = NULL)
  )
  
  all_results <- list()
  
  for (method_name in names(threshold_methods)) {
    cat("   Evaluating with", method_name, "thresholds...\n")
    
    results <- get_confusion_matrices(
      data = merged_metrics,
      truth_col = "true_link",
      metric_cols = c("IBD", "IBS", "phylo"),
      threshold_method = threshold_methods[[method_name]]$method,
      fixed_thresholds = threshold_methods[[method_name]]$fixed_thresholds
    )
    
    all_results[[method_name]] <- results
    cat("     ✅", method_name, "evaluation completed\n")
    
    # Print summary
    cat("     Thresholds used:\n")
    for (metric in names(results$thresholds_used)) {
      cat("       -", metric, ":", round(results$thresholds_used[metric], 4), "\n")
    }
  }
  
  # Save evaluation results
  # cat("   Saving evaluation results...\n")
  # base_output <- file.path(outdir, base_name)
  # 
  # # Save RDS files
  # saveRDS(all_results$optimal_youden, paste0(base_output, "_optimal_youden.rds"))
  # saveRDS(all_results$fixed, paste0(base_output, "_fixed.rds")) 
  # saveRDS(all_results$median, paste0(base_output, "_median.rds"))
  # 
  # # Save Excel summary
  # excel_output <- paste0(base_output, "_confusion_metrics.xlsx")
  # writexl::write_xlsx(
  #   list(
  #     optimal_youden = all_results$optimal_youden$metrics_summary,
  #     fixed = all_results$fixed$metrics_summary,
  #     median = all_results$median$metrics_summary
  #   ),
  #   excel_output
  # )
  # cat("   ✅ Evaluation results saved to:", excel_output, "\n")
  # 
  # # Generate confusion matrix plots
  # cat("   Generating confusion matrix plots...\n")
  # plot_confusion_matrices(all_results$optimal_youden, paste0(base_output, "_oy_confusion_metrics"))
  # cat("   ✅ Confusion matrix plots saved\n")
  
  # ============================================================================
  # CLASSIFICATION PERFORMANCE
  # ============================================================================
  
  cat("\n📊 Evaluating classification performance...\n")
  
  # Comprehensive predictor evaluation
  cat("   Running comprehensive predictor evaluation...\n")
  evals <- evaluate_predictors(merged_metrics, predictors = c("IBD", "IBS", "phylo"), label_col = "true_link")
  
  # # Save classification results
  # class_output <- file.path(outdir, paste0(base_name, "_classification.rds"))
  # saveRDS(evals, class_output)
  # cat("   ✅ Classification results saved:", class_output, "\n")
  
  # Extract and save AUC/AUPR metrics
  cat("   Extracting AUC/AUPR metrics...\n")
  performance_table <- map_dfr(names(evals), ~ {
    tibble(
      predictor = .x, 
      auc = evals[[.x]]$auc, 
      aupr = evals[[.x]]$aupr
    )
  })
  
  # table_output <- file.path(outdir, paste0(base_name, "_classification.xlsx"))
  # writexl::write_xlsx(performance_table, table_output)
  # cat("   ✅ Performance metrics saved:", table_output, "\n")
  
  # Generate PR curves for each predictor
  cat("   Generating precision-recall curves...\n")
  predictors <- c("IBD", "IBS", "phylo")
  precision_table <- list()
  curve_data <- list()
  
  for (variable in predictors) {
    cat("   - Processing", variable, "metric...\n")
    
    # Calculate precision@k
    precision_table[[variable]] <- precision_recall_at_k(merged_metrics, score_col = variable, label_col = "true_link", ks = c(1, 3, 5))
    
    # Compute ROC/PR curves
    curve_data[[variable]] <- compute_roc_pr(merged_metrics, variable, truth_col = "true_link")
    # pr_auc <- curve_data$pr$auc.integral
    
    # cat("     PR-AUC =", round(pr_auc, 3), "\n")
    # 
    # # Create PR curve plot
    # pr_df <- data.frame(
    #   recall = curve_data$pr$curve[, 1],
    #   precision = curve_data$pr$curve[, 2]
    # )
    # 
    # plot_curve <- ggplot(pr_df, aes(x = recall, y = precision)) +
    #   geom_line(color = "blue", linewidth = 2, lineend = 'round') +
    #   labs(
    #     x = "Recall", y = "Precision",
    #     title = sprintf("Precision-Recall Curve: %s", variable),
    #     subtitle = sprintf("AUC = %.3f", pr_auc)
    #   ) +
    #   theme_minimal() +
    #   theme(
    #     plot.title = element_text(face = "bold", size = 19), 
    #     plot.subtitle = element_text(face = "bold", size = 15),
    #     axis.title = element_text(size = 15, color = 'black', face = "bold"),
    #     axis.text = element_text(size = 12, color = 'black')
    #   )
    # 
    # # Save PR curve
    # curve_output <- paste0(base_output, "_prauc_", variable, ".pdf")
    # ggsave(curve_output, plot = plot_curve, width = 20, height = 15, units = "cm")
    # cat("     ✅ PR curve saved:", curve_output, "\n")
  }
  
  # ============================================================================
  # FINAL SUMMARY
  # ============================================================================
  
  cat("\n✅ Pipeline execution completed successfully!\n")
  cat("📁 All results saved to:", outdir, "\n")
  cat("📊 Performance summary:\n")
  
  # Print final performance summary
  for (i in 1:nrow(performance_table)) {
    row <- performance_table[i, ]
    cat("   -", row$predictor, ": AUC =", round(row$auc, 3), "AUPR =", round(row$aupr, 3), "\n")
  }
  
  # Return comprehensive results
  return(list(
    input_data = list(
      ground_truth = true_links,
      ibd = ibd_df,
      ibs = ibs_df,
      patristic = phylo_df
    ),
    merged_data = list(
      conservative = merged_conservative,
      strict = merged_strict,
      union = merged_union
    ),
    correlations = corr_res,
    evaluations = all_results,
    classification = evals,
    performance_summary = performance_table,
    precision_table = precision_table,
    curve_data = curve_data
  ))
}

# ========================================================================
#               MAIN (RUN BENCHMARKING ACROSS REPLICATES)
# ========================================================================              
# Define the base directory and replicate paths
library(tidyverse)

source("simulations/R/utils.R")
source_all_scripts()

base_dir <- "simulations/single_run"
replicates <- sprintf("replicate%d", 1:10)

results <- map(replicates, ~ {
  replicate_dir <- file.path(base_dir, "inferred", .x)
  phylo_dir <- file.path(base_dir, "phylo_output", .x)
  
  run_full_pipeline(
    ibd_file = file.path(replicate_dir, "inferred_ibd_hmm.tsv"),
    ibs_file = file.path(replicate_dir, "ibs_matrix.rds"),
    treefile = file.path(phylo_dir, paste0("run", str_extract(.x, "\\d+"), "_chr1_modelfinder.treefile")),
    truth_file = file.path(replicate_dir, "true_ibd_summary.tsv"),
    gen_cutoff = 25, 
    outdir = "simulations/single_run/evaluation_results"
  )
})

names(results) <- replicates
# Save output list across all replicates
saveRDS(results, "simulations/single_run/evaluation_results.rds")

# Now, perform the evaluation for different Generation cutoffs
gen <- c(3, 5, 10, 15, 20, 25)
all_results <- list()

for(G in gen){
  
  results <- map(replicates, ~ {
    replicate_dir <- file.path(base_dir, "inferred", .x)
    phylo_dir <- file.path(base_dir, "phylo_output", .x)
    
    run_full_pipeline(
      ibd_file = file.path(replicate_dir, "inferred_ibd_hmm.tsv"),
      ibs_file = file.path(replicate_dir, "ibs_matrix.rds"),
      treefile = file.path(phylo_dir, paste0("run", str_extract(.x, "\\d+"), "_chr1_modelfinder.treefile")),
      truth_file = file.path(replicate_dir, "true_ibd_summary.tsv"),
      gen_cutoff = G, 
      outdir = "simulations/single_run/evaluation_results"
    )
  })
  all_results[[paste0("gen", G)]] <- results
}


# =============================================================================
# NOTES FOR REPRODUCIBILITY:
# 1. Ensure all file paths are correct for your directory structure
# 2. Verify that custom functions (phylo_long, merge_metrics, etc.) are loaded
# 3. The gen_cutoff = 25 defines what constitutes a "close relative"
# 4. Logit transformation is used for Pearson correlation to improve normality
# 5. Sample IDs are standardized with "tsk_" prefix for consistent merging
# =============================================================================

# results <- map(1:20, ~ {
#   rep_num <- .x
#   replicate_dir <- file.path(base_dir, "inferred", paste0("replicate", rep_num))
#   phylo_dir <- file.path(base_dir, "phylo_output", paste0("replicate", rep_num))
#   
#   run_full_pipeline(
#     ibd_file = file.path(replicate_dir, "inferred_ibd_hmm.tsv"),
#     ibs_file = file.path(replicate_dir, "inferred_ibs.tsv"),
#     treefile = file.path(phylo_dir, paste0("run", rep_num, "_chr1_modelfinder.treefile")),
#     truth_file = file.path(replicate_dir, "true_ibd_summary.tsv")
#   )
# })




















