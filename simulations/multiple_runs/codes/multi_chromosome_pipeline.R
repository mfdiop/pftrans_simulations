
#' Helper function to convert matrix to long format (efficient)
matrix_to_long_efficient <- function(matrix_data, value_name = "value", id_cols = c("id1", "id2")) {
  
  if (!is.matrix(matrix_data)) {
    stop("Input must be a matrix")
  }
  
  n <- nrow(matrix_data)
  if (n != ncol(matrix_data)) {
    stop("Input matrix must be square")
  }
  
  # Get only upper triangle (excluding diagonal)
  upper_tri_indices <- which(upper.tri(matrix_data), arr.ind = TRUE)
  
  long_df <- data.frame(
    id1 = rownames(matrix_data)[upper_tri_indices[, 1]],
    id2 = colnames(matrix_data)[upper_tri_indices[, 2]], 
    value = matrix_data[upper_tri_indices],
    stringsAsFactors = FALSE
  )
  
  names(long_df) <- c(id_cols, value_name)
  
  return(long_df)
}

#' Helper function for phylogenetic distance calculation
phylo_long <- function(treefile, method = "patristic") {
  
  if (!file.exists(treefile)) {
    stop("Tree file not found: ", treefile)
  }
  
  tree <- read.tree(treefile)
  dist_matrix <- cophenetic(tree)
  
  # Convert to long format with unique pairs only
  phylo_long <- matrix_to_long_efficient(dist_matrix, value_name = "phylo")
  
  return(phylo_long)
}

#' Helper function to plot confusion matrices
plot_confusion_matrices <- function(confusion_results, output_base) {
  
  cm_list <- confusion_results$confusion_matrices
  plots <- list()
  
  for (metric_name in names(cm_list)) {
    cm <- as.data.frame(cm_list[[metric_name]])
    
    p <- ggplot(cm, aes(x = Predicted, y = Actual, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), vjust = 0.5, size = 4) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      labs(
        title = paste("Confusion Matrix:", metric_name),
        subtitle = paste("Threshold:", round(confusion_results$thresholds_used[metric_name], 4))
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    plots[[metric_name]] <- p
  }
  
  # Save combined plot
  combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = length(plots))
  ggsave(paste0(output_base, ".pdf"), combined_plot, width = 12, height = 4)
  
  # Save individual plots
  for (metric_name in names(plots)) {
    ggsave(
      paste0(output_base, "_", metric_name, ".pdf"),
      plots[[metric_name]], width = 6, height = 4
    )
  }
}

# Place these in a separate utils.R file or at the top of your script

#' Define ground truth relationships based on generational cutoff
define_ground_truth <- function(true_links, gen_cutoff = 5) {
  true_links %>%
    mutate(true_link = generations <= gen_cutoff)
}

#' Merge metrics with ground truth data
merge_metrics_with <- function(df_ground, ibd_df, ibs_df, phylo_df, merge_method = "ground_truth") {
  # Standardize all dataframes to use sorted pair_key
  standardize_df <- function(df, id_cols = c("id1", "id2")) {
    df %>%
      rowwise() %>%
      mutate(pair_key = paste(sort(c(!!sym(id_cols[1]), !!sym(id_cols[2]))), collapse = "_")) %>%
      ungroup() %>%
      distinct(pair_key, .keep_all = TRUE)
  }
  
  ground_std <- standardize_df(df_ground)
  ibd_std <- standardize_df(ibd_df)
  ibs_std <- standardize_df(ibs_df)
  phylo_std <- standardize_df(phylo_df)
  
  # Determine which pairs to keep based on merge method
  if (merge_method == "ground_truth") {
    keep_pairs <- ground_std$pair_key
  } else if (merge_method == "intersection") {
    keep_pairs <- Reduce(intersect, list(ground_std$pair_key, ibd_std$pair_key, 
                                         ibs_std$pair_key, phylo_std$pair_key))
  } else if (merge_method == "union") {
    keep_pairs <- unique(c(ground_std$pair_key, ibd_std$pair_key, 
                           ibs_std$pair_key, phylo_std$pair_key))
  }
  
  # Create base dataframe and merge
  base_df <- data.frame(pair_key = keep_pairs)
  
  safe_join <- function(left_df, right_df, right_name, by = "pair_key") {
    right_cols_to_join <- setdiff(names(right_df), c("id1", "id2", "pair_key"))
    left_df %>%
      left_join(select(right_df, all_of(c(by, right_cols_to_join))), by = by)
  }
  
  merged <- base_df %>%
    safe_join(ground_std, "ground_truth") %>%
    safe_join(ibd_std, "IBD") %>%
    safe_join(ibs_std, "IBS") %>%
    safe_join(phylo_std, "phylo")
  
  return(merged)
}

# Helper function for string concatenation
`%+%` <- function(a, b) paste0(a, b)

#' Run Single Chromosome Pipeline
#' 
#' Modified version of the original pipeline for integration with multi-chromosome framework
#' Returns results object instead of just saving files
#'
#' @param ibd_file Path to IBD inference results file (TSV format)
#' @param ibs_file Path to IBS matrix file (RDS format)  
#' @param treefile Path to phylogenetic tree file
#' @param truth_file Path to ground truth relationships file (TSV format)
#' @param gen_cutoff Generational cutoff for defining close relatives (default: 5)
#' @param outdir Output directory for results
#' @param combo_id Combination identifier for output files
#'
#' @return List containing analysis results
run_single_chromosome_pipeline <- function(ibd_file, ibs_file, treefile, truth_file, 
                                           gen_cutoff = 5, outdir = "out", combo_id = NULL) {
  
  source(normalizePath("simulations/R/utils.R", winslash = "\\"))
  source_all_scripts()
  
  # ============================================================================
  # INITIALIZATION
  # ============================================================================
  
  cat("   🧬 Starting single chromosome pipeline...\n")
  
  # Use combo_id for output filenames if provided
  if (is.null(combo_id)) {
    base_name <- tools::file_path_sans_ext(basename(ibd_file))
  } else {
    base_name <- combo_id
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # ============================================================================
  # DATA LOADING AND PREPROCESSING
  # ============================================================================
  
  cat("   📥 Loading input data...\n")
  
  # Load IBD data
  cat("     Reading IBD file...\n")
  ibd_df <- read_tsv(ibd_file, show_col_types = FALSE)
  cat("     ✅ IBD data loaded:", nrow(ibd_df), "rows\n")
  
  # Load and process IBS data
  cat("     Reading IBS file...\n")
  ibs_matrix <- readRDS(ibs_file)
  cat("     ✅ IBS matrix loaded:", nrow(ibs_matrix), "x", ncol(ibs_matrix), "matrix\n")
  
  # Convert IBS matrix to long format
  cat("     Converting IBS matrix to long format...\n")
  ibs_df <- matrix_to_long_efficient(ibs_matrix, value_name = "IBS") %>% 
    arrange(id1)
  cat("     ✅ IBS data converted:", nrow(ibs_df), "unique pairs\n")
  
  # Save IBS long format for future use
  ibs_output_file <- file.path(outdir, paste0(base_name, "_inferred_ibs.tsv"))
  write_tsv(ibs_df, ibs_output_file)
  
  # Load and process ground truth data
  cat("     Reading ground truth file...\n")
  true_links <- read_tsv(truth_file, show_col_types = FALSE) %>% 
    mutate(
      Id1 = paste0("tsk_", Id1), 
      Id2 = paste0("tsk_", Id2), 
      ibd = round(total_ibd_prop, 3)) %>% 
    select(Id1, Id2, ibd, min_tmrca) %>% 
    rename(generations = min_tmrca)
  cat("     ✅ Ground truth data loaded:", nrow(true_links), "rows\n")
  
  # Calculate phylogenetic distances
  cat("     Calculating phylogenetic distances...\n")
  phylo_df <- phylo_long(treefile, method = "patristic")
  cat("     ✅ Phylogenetic distances calculated:", nrow(phylo_df), "pairs\n")
  
  # ============================================================================
  # DATA INTEGRATION
  # ============================================================================
  
  cat("   🔗 Integrating datasets...\n")
  
  # Define ground truth relationships
  df_ground <- define_ground_truth(true_links, gen_cutoff = gen_cutoff) %>% 
    rename(id1 = Id1, id2 = Id2)
  cat("     ✅ Ground truth defined:", sum(df_ground$true_link), "true links\n")
  
  # Merge datasets (using conservative approach)
  merged_conservative <- merge_metrics_with(
    df_ground = df_ground,
    ibd_df = ibd_df, 
    ibs_df = ibs_df, 
    phylo_df = phylo_df,
    merge_method = "ground_truth"
  )
  cat("     ✅ Data merged:", nrow(merged_conservative), "pairs\n")
  
  # ============================================================================
  # STATISTICAL ANALYSIS
  # ============================================================================
  
  cat("   📈 Performing statistical analysis...\n")
  
  # Correlation analysis
  corr_res <- compute_correlations_v1(merged_conservative, "IBD", "IBS", transform_for_pearson = "logit")
  cat("     ✅ Correlation analysis completed\n")
  
  # Prepare for Mantel test
  merged_metrics <- merged_conservative %>% 
    separate(pair_key, into = c("id1", "id2"), sep = "(?<=[0-9])_", extra = "merge")
  
  # ============================================================================
  # MODEL EVALUATION
  # ============================================================================
  
  cat("   🎯 Evaluating model performance...\n")
  
  # Evaluation with different threshold methods
  threshold_methods <- list(
    optimal_youden = list(method = "optimal_youden", fixed_thresholds = NULL),
    fixed = list(method = "fixed", fixed_thresholds = c(IBD = 0.125, IBS = 0.8, phylo = 0.5)),
    median = list(method = "median", fixed_thresholds = NULL)
  )
  
  all_evaluations <- list()
  
  for (method_name in names(threshold_methods)) {
    cat("     Evaluating with", method_name, "thresholds...\n")
    
    results <- get_confusion_matrices(
      data = merged_metrics,
      truth_col = "true_link",
      metric_cols = c("IBD", "IBS", "phylo"),
      threshold_method = threshold_methods[[method_name]]$method,
      fixed_thresholds = threshold_methods[[method_name]]$fixed_thresholds
    )
    
    all_evaluations[[method_name]] <- results
  }
  
  # ============================================================================
  # CLASSIFICATION PERFORMANCE
  # ============================================================================
  
  cat("   📊 Evaluating classification performance...\n")
  
  # Comprehensive predictor evaluation
  evals <- evaluate_predictors(merged_metrics, predictors = c("IBD", "IBS", "phylo"), label_col = "true_link")
  
  # Extract performance metrics
  performance_table <- map_dfr(names(evals), ~ {
    tibble(
      predictor = .x, 
      auc = round(evals[[.x]]$auc, 3), 
      aupr = round(evals[[.x]]$aupr, 3)
    )
  })
  
  cat("     ✅ Classification evaluation completed\n")
  
  # ============================================================================
  # SAVE OUTPUT FILES (for individual inspection)
  # ============================================================================
  
  cat("   💾 Saving output files...\n")
  saveRDS(all_evaluations, file.path(outdir, paste0(base_name, "_evaluations.rds")))
  
  # Save evaluation results
  # saveRDS(all_evaluations$optimal_youden, file.path(outdir, paste0(base_name, "_optimal_youden.rds")))
  # saveRDS(all_evaluations$fixed, file.path(outdir, paste0(base_name, "_fixed.rds")))
  # saveRDS(all_evaluations$median, file.path(outdir, paste0(base_name, "_median.rds")))
  
  # Save Excel summary
  writexl::write_xlsx(
    list(
      optimal_youden = all_evaluations$optimal_youden$metrics_summary,
      fixed = all_evaluations$fixed$metrics_summary,
      median = all_evaluations$median$metrics_summary
    ),
    file.path(outdir, paste0(base_name, "_confusion_metrics.xlsx"))
  )
  
  # Save classification results
  saveRDS(evals, file.path(outdir, paste0(base_name, "_classification.rds")))
  writexl::write_xlsx(performance_table, file.path(outdir, paste0(base_name, "_classification.xlsx")))
  
  # ============================================================================
  # GENERATE VISUALIZATIONS
  # ============================================================================
  
  cat("   📊 Generating visualizations...\n")
  
  # Generate key plots
  tryCatch({
    # PR curves for each predictor
    predictors <- c("IBD", "IBS", "phylo")
    
    for (variable in predictors) {
      curve_data <- compute_roc_pr(merged_metrics, variable, truth_col = "true_link")
      pr_auc <- curve_data$pr$auc.integral
      
      pr_df <- data.frame(
        recall = curve_data$pr$curve[, 1],
        precision = curve_data$pr$curve[, 2]
      )
      
      plot_curve <- ggplot(pr_df, aes(x = recall, y = precision)) +
        geom_line(color = "blue", linewidth = 1, lineend = 'round') +
        labs(
          x = "Recall", y = "Precision",
          title = sprintf("PR Curve: %s", variable),
          subtitle = sprintf("AUC = %.3f", pr_auc)
        ) +
        theme_minimal()
      
      ggsave(
        file.path(outdir, paste0(base_name, "_prauc_", variable, ".pdf")),
        plot_curve, width = 8, height = 6
      )
    }
    
    # Confusion matrix plot
    plot_confusion_matrices(all_evaluations$optimal_youden, 
                            file.path(outdir, paste0(base_name, "_confusion_matrices")))
    
    cat("     ✅ Visualizations generated\n")
  }, error = function(e) {
    cat("     ⚠️  Visualization generation failed:", e$message, "\n")
  })
  
  # ============================================================================
  # RETURN RESULTS OBJECT
  # ============================================================================
  
  cat("   ✅ Single chromosome pipeline completed!\n")
  
  # Return comprehensive results object
  results <- list(
    metadata = list(
      base_name = base_name,
      n_pairs = nrow(merged_conservative),
      n_true_links = sum(df_ground$true_link),
      files = list(
        ibd = ibd_file,
        ibs = ibs_file,
        truth = truth_file,
        tree = treefile
      )
    ),
    performance_metrics = performance_table,
    correlations = corr_res,
    evaluations = all_evaluations,
    classification = evals,
    merged_data = list(
      conservative = merged_conservative,
      metrics = merged_metrics
    )
  )
  
  return(results)
}

# Summary function for multi-chromosome data
create_multi_chromosome_summary <- function(all_results, outdir) {
  
  cat("   Creating multi-chromosome summary...\n")
  
  summary_data <- map_dfr(all_results, function(result) {
    metadata <- result$metadata
    results <- result$results
    
    # Extract key metrics
    perf_metrics <- results$performance_metrics
    
    tibble(
      replicate = metadata$replicate,
      recombination_rate = metadata$recombination_rate,
      chromosome = metadata$chromosome,
      output_id = metadata$output_id,
      IBD_AUC = ifelse("IBD" %in% perf_metrics$predictor, 
                       perf_metrics$auc[perf_metrics$predictor == "IBD"], NA),
      IBS_AUC = ifelse("IBS" %in% perf_metrics$predictor,
                       perf_metrics$auc[perf_metrics$predictor == "IBS"], NA),
      Phylo_AUC = ifelse("phylo" %in% perf_metrics$predictor,
                         perf_metrics$auc[perf_metrics$predictor == "phylo"], NA),
      IBD_AUPR = ifelse("IBD" %in% perf_metrics$predictor,
                        perf_metrics$aupr[perf_metrics$predictor == "IBD"], NA),
      IBS_AUPR = ifelse("IBS" %in% perf_metrics$predictor,
                        perf_metrics$aupr[perf_metrics$predictor == "IBS"], NA),
      Phylo_AUPR = ifelse("phylo" %in% perf_metrics$predictor,
                          perf_metrics$aupr[perf_metrics$predictor == "phylo"], NA)
    )
  })
  
  # Save summary
  summary_file <- file.path(outdir, "multi_chromosome_summary.csv")
  write_csv(summary_data, summary_file)
  cat("   💾 Summary saved to:", summary_file, "\n")
  
  return(summary_data)
}

# # Comparative plotting function
# generate_comparative_plots <- function(all_results, outdir) {
#   
#   cat("   Generating comparative plots...\n")
#   
#   # This would create plots comparing performance across:
#   # - Different recombination rates
#   # - Different chromosomes
#   # - Different replicates
#   
#   # Implementation depends on your specific plotting needs
#   # You can adapt the plotting functions I provided earlier
#   
#   cat("   ✅ Comparative plots generated\n")
# }


#' Run Full Relationship Inference Pipeline for Multi-Chromosome, Multi-Rate Data
#' 
#' This function processes simulations across multiple chromosomes, recombination rates,
#' and replicates. Designed for SLURM batch processing with checkpointing.
#'
#' @param base_dir Base directory containing inferred and phylo_results folders
#' @param replicates Vector of replicate numbers (e.g., 1:5)
#' @param recombination_rates Vector of recombination rates (e.g., c(1e-9, 1e-8, 1e-7, 1e-6))
#' @param chromosomes Vector of chromosome numbers (e.g., 1:5)
#' @param outdir Base output directory for results
#' @param gen_cutoff Generational cutoff for defining close relatives (default: 5)
#' @param skip_processed Skip already processed combinations (default: TRUE)
#' @param slurm_job_id SLURM job ID for checkpointing (optional)
#'
#' @return List containing analysis results and saves output files to disk
#' 
#' @examples
#' \dontrun{
#' # For SLURM array job
#' results <- run_multi_chromosome_pipeline(
#'   base_dir = "simulations/multiple_runs",
#'   replicates = 1:5,
#'   recombination_rates = c(1e-9, 1e-8, 1e-7, 1e-6),
#'   chromosomes = 1:3,
#'   outdir = "results_multi_chrom",
#'   gen_cutoff = 25,
#'   slurm_job_id = Sys.getenv("SLURM_ARRAY_TASK_ID")
#' )
#' }


# Define input file paths for a single simulation replicate
truth_file <- "simulations/multiple_runs/inferred/rep1/run1_rec1e09_chr1/true_ibd_summary.tsv"    # Ground truth IBD data
ibd_file <- "simulations/multiple_runs/inferred/rep1/run1_rec1e09_chr1/inferred_ibd_hmm.tsv"      # Inferred IBD segments
ibs_file <- "simulations/multiple_runs/inferred/rep1/run1_rec1e09_chr1/ibs_matrix.rds"          # Inferred IBS segments
treefile <- "simulations/multiple_runs/phylo_results/rep1/run1_rec1e09_chr1_modelfinder.treefile"  # Phylogenetic tree


run_multi_chromosome_pipeline <- function(base_dir, replicates, recombination_rates, 
                                          chromosomes, outdir = "simulations/multiple_runs/results_multi_chrom",
                                          gen_cutoff = 5, skip_processed = TRUE,
                                          slurm_job_id = NULL) {
  # ============================================================================
  # INITIALIZATION AND SETUP
  # ============================================================================
  
  cat("🚀 Starting Multi-Chromosome/Multi-Rate Pipeline...\n")
  cat("📁 Base directory:", base_dir, "\n")
  cat("📁 Output directory:", outdir, "\n")
  cat("🔢 Replicates:", paste(replicates, collapse = ", "), "\n")
  cat("🔄 Recombination rates:", paste(recombination_rates, collapse = ", "), "\n")
  cat("🧬 Chromosomes:", paste(chromosomes, collapse = ", "), "\n")
  
  if (!is.null(slurm_job_id)) {
    cat("⚡ SLURM Job ID:", slurm_job_id, "\n")
  }
  
  # Create output directory structure
  if (!dir.exists(outdir)) {
    cat("📂 Creating output directory:", outdir, "\n")
    dir.create(outdir, recursive = TRUE)
  }
  
  # Load required packages
  cat("📦 Loading required packages...\n")
  required_packages <- c("tidyverse", "hexbin", "gridExtra", "PRROC", "pROC", 
                         "vegan", "scales", "ape", "phangorn", "ggpubr", "data.table")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop("❌ Package ", pkg, " is not installed. Please install it first.")
    }
  }
  cat("✅ All packages loaded successfully\n")
  
  # ============================================================================
  # GENERATE PROCESSING GRID
  # ============================================================================
  
  cat("\n📋 Generating processing grid...\n")
  
  # Create all combinations of parameters
  processing_grid <- expand.grid(
    replicate = replicates,
    recombination_rate = recombination_rates,
    chromosome = chromosomes,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      # Format recombination rate for filename matching
      rec_rate_label = sprintf("rec%.0e", recombination_rate),
      rec_rate_label = gsub("\\+", "", rec_rate_label),  # Remove + from scientific notation
      rec_rate_label = gsub("-0", "-", rec_rate_label),  # Fix negative exponents
      
      # Generate expected file patterns
      ibd_pattern = paste0("run", replicate, "_", rec_rate_label, "_chr", chromosome),
      tree_pattern = paste0("run", replicate, "_chr", chromosome, "_modelfinder.treefile"),
      
      # Generate output identifier
      output_id = paste0("rep", replicate, "_", rec_rate_label, "_chr", chromosome)
    )
  
  cat("📊 Total combinations to process:", nrow(processing_grid), "\n")
  
  # ============================================================================
  # CHECKPOINTING SETUP
  # ============================================================================
  
  if (skip_processed) {
    cat("🔍 Checking for already processed combinations...\n")
    
    processed_files <- list.files(outdir, pattern = "*_classification.xlsx", full.names = FALSE)
    processed_ids <- gsub("_classification\\.xlsx$", "", processed_files)
    
    processing_grid <- processing_grid %>%
      mutate(already_processed = output_id %in% processed_ids)
    
    n_processed <- sum(processing_grid$already_processed)
    n_remaining <- nrow(processing_grid) - n_processed
    
    cat("✅ Already processed:", n_processed, "combinations\n")
    cat("⏳ Remaining to process:", n_remaining, "combinations\n")
    
    if (n_remaining == 0) {
      cat("🎉 All combinations already processed!\n")
      return(NULL)
    }
  }
  
  # ============================================================================
  # SLURM ARRAY JOB SUPPORT
  # ============================================================================
  
  if (!is.null(slurm_job_id)) {
    slurm_job_id <- as.numeric(slurm_job_id)
    cat("⚡ Running as SLURM array job, processing combination:", slurm_job_id, "\n")
    
    if (slurm_job_id > nrow(processing_grid)) {
      cat("❌ SLURM job ID exceeds number of combinations\n")
      return(NULL)
    }
    
    # Process only the assigned combination
    processing_grid <- processing_grid[slurm_job_id, , drop = FALSE]
  }
  
  # ============================================================================
  # PROCESS EACH COMBINATION
  # ============================================================================
  
  all_results <- list()
  
  for (i in 1:nrow(processing_grid)) {
    current_combo <- processing_grid[i, ]
    
    cat("\n" %+% strrep("=", 60) %+% "\n")
    cat("  🔄 Processing Combination", i, "of", nrow(processing_grid), "\n")
    cat("  📊 Replicate:", current_combo$replicate, "\n")
    cat("  🔄 Recombination rate:", current_combo$recombination_rate, "\n")
    cat("  🧬 Chromosome:", current_combo$chromosome, "\n")
    cat("  🎯 Output ID:", current_combo$output_id, "\n")
    
    # Skip if already processed
    if (skip_processed && current_combo$already_processed) {
      cat("⏭️  Already processed, skipping...\n")
      next
    }
    
    tryCatch({
      # ========================================================================
      # FILE PATH CONSTRUCTION
      # ========================================================================
      
      # Construct file paths
      inferred_dir <- file.path(base_dir, "inferred", paste0("rep", current_combo$replicate))
      phylo_dir <- file.path(base_dir, "phylo_results", paste0("rep", current_combo$replicate))
      
      cat("📁 Inferred directory:", inferred_dir, "\n")
      cat("📁 Phylo directory:", phylo_dir, "\n")
      
      # Find matching files
      ibd_files <- list.files(inferred_dir, pattern = paste0(current_combo$ibd_pattern, ".*\\.tsv$"), full.names = TRUE)
      ibs_files <- list.files(inferred_dir, pattern = paste0(current_combo$ibd_pattern, ".*_ibs\\.rds$"), full.names = TRUE)
      truth_files <- list.files(inferred_dir, pattern = paste0(current_combo$ibd_pattern, ".*_true_ibd_summary\\.tsv$"), full.names = TRUE)
      tree_files <- list.files(phylo_dir, pattern = current_combo$tree_pattern, full.names = TRUE)
      
      # Check if all required files exist
      if (length(ibd_files) == 0) {
        stop("❌ No IBD file found for pattern: ", current_combo$ibd_pattern)
      }
      if (length(ibs_files) == 0) {
        stop("❌ No IBS file found for pattern: ", current_combo$ibd_pattern)
      }
      if (length(truth_files) == 0) {
        stop("❌ No truth file found for pattern: ", current_combo$ibd_pattern)
      }
      if (length(tree_files) == 0) {
        stop("❌ No tree file found for pattern: ", current_combo$tree_pattern)
      }
      
      # Use first matching file (should be only one per pattern)
      ibd_file <- ibd_files[1]
      ibs_file <- ibs_files[1]
      truth_file <- truth_files[1]
      tree_file <- tree_files[1]
      
      cat("📄 IBD file:", basename(ibd_file), "\n")
      cat("📄 IBS file:", basename(ibs_file), "\n")
      cat("📄 Truth file:", basename(truth_file), "\n")
      cat("📄 Tree file:", basename(tree_file), "\n")
      
      # ========================================================================
      # CREATE COMBINATION-SPECIFIC OUTPUT DIRECTORY
      # ========================================================================
      
      combo_outdir <- file.path(outdir, current_combo$output_id)
      if (!dir.exists(combo_outdir)) {
        dir.create(combo_outdir, recursive = TRUE)
      }
      cat("📂 Combination output directory:", combo_outdir, "\n")
      
      # ========================================================================
      # RUN SINGLE PIPELINE (Reusing your existing function)
      # ========================================================================
      
      cat("🔧 Running single pipeline for this combination...\n")
      
      single_result <- run_single_chromosome_pipeline(
        ibd_file = ibd_file,
        ibs_file = ibs_file, 
        treefile = tree_file,
        truth_file = truth_file,
        gen_cutoff = gen_cutoff,
        outdir = combo_outdir,
        combo_id = current_combo$output_id
      )
      
      # ========================================================================
      # STORE RESULTS WITH METADATA
      # ========================================================================
      
      result_with_metadata <- list(
        metadata = list(
          replicate = current_combo$replicate,
          recombination_rate = current_combo$recombination_rate,
          chromosome = current_combo$chromosome,
          output_id = current_combo$output_id,
          files = list(
            ibd = ibd_file,
            ibs = ibs_file,
            truth = truth_file,
            tree = tree_file
          )
        ),
        results = single_result
      )
      
      all_results[[current_combo$output_id]] <- result_with_metadata
      
      # Save individual combination results
      result_file <- file.path(combo_outdir, paste0(current_combo$output_id, "_full_results.rds"))
      saveRDS(result_with_metadata, result_file)
      cat("💾 Saved individual results to:", result_file, "\n")
      
      cat("✅ Successfully processed combination:", current_combo$output_id, "\n")
      
    }, error = function(e) {
      cat("❌ ERROR processing combination:", current_combo$output_id, "\n")
      cat("   Error message:", e$message, "\n")
      
      # Log error for debugging
      error_log <- file.path(outdir, "processing_errors.log")
      write_lines(
        paste(Sys.time(), current_combo$output_id, e$message, sep = "\t"),
        error_log,
        append = TRUE
      )
      cat("   📝 Error logged to:", error_log, "\n")
    })
  }
  
  # ============================================================================
  # AGGREGATE RESULTS ACROSS ALL COMBINATIONS
  # ============================================================================
  
  cat("\n" %+% strrep("=", 60) %+% "\n")
  cat("📊 Aggregating results across all combinations...\n")
  
  if (length(all_results) > 0) {
    # Create summary across all processed combinations
    summary_results <- create_multi_chromosome_summary(all_results, outdir)
    
    # Save aggregated results
    aggregated_file <- file.path(outdir, "aggregated_results_all_combinations.rds")
    saveRDS(list(individual_results = all_results, summary = summary_results), aggregated_file)
    cat("💾 Saved aggregated results to:", aggregated_file, "\n")
    
    # # Generate comparative plots across combinations
    # generate_comparative_plots(all_results, outdir)
    
    cat("✅ Pipeline completed successfully!\n")
    cat("📁 All results saved to:", outdir, "\n")
    cat("📊 Processed", length(all_results), "combinations\n")
    
    return(list(individual_results = all_results, summary = summary_results))
  } else {
    cat("⚠️  No combinations were successfully processed\n")
    return(NULL)
  }
}



# # Test one combination
# results <- run_multi_chromosome_pipeline(
#   base_dir = "simulations/single_run",
#   replicates = 1,
#   recombination_rates = 1e-8,
#   chromosomes = 1,
#   outdir = "test_results",
#   skip_processed = FALSE
# )






