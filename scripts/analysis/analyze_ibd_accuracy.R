#!/usr/bin/env Rscript

################################################################################
# Script: analyze_ibd_accuracy.R
# Purpose: Evaluate IBD detection accuracy across different recombination rates
# Author: [Your Name]
# Date: 2025-11-05
#
# Description:
#   Compares true IBD segments (from tskibd/tree sequences) with inferred IBD
#   segments (from hmmIBD, hapIBD, isoRelate, etc.) to calculate:
#   - False Negative Rate (FNR): Proportion of true IBD missed
#   - False Positive Rate (FPR): Proportion of called IBD that's false
#   - Precision, Recall, F1 score
#
# Input:
#   - True IBD segments from tskibd
#   - Inferred IBD segments from various callers
#
# Output:
#   - Summary statistics table
#   - Plots of FNR and FPR vs recombination rate
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(GenomicRanges)
  library(cowplot)
  library(scales)
})

################################################################################
# COMMAND LINE ARGUMENTS
################################################################################

option_list <- list(
  make_option("--input_dir", type = "character", default = "ibd_results",
              help = "Directory containing IBD results [default: %default]"),
  make_option("--output_dir", type = "character", default = "ibd_accuracy",
              help = "Output directory for results [default: %default]"),
  make_option("--min_length", type = "numeric", default = 0,
              help = "Minimum IBD segment length in bp [default: %default]"),
  make_option("--overlap_threshold", type = "numeric", default = 0.5,
              help = "Minimum overlap fraction for segment matching [default: %default]"),
  make_option("--methods", type = "character", default = "hmm,iso,ibs",
              help = "Comma-separated list of methods to evaluate [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Parse methods
methods_to_eval <- strsplit(opt$methods, ",")[[1]]

# Create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Convert IBD data frame to GRanges object
#' 
#' @param ibd_df Data frame with columns: Id1, Id2, Start, End, (Chr)
#' @param genome_length Total genome length for chromosome bounds
#' @return GRanges object with pair information in metadata
ibd_to_granges <- function(ibd_df, genome_length = NULL) {
  
  if (nrow(ibd_df) == 0) {
    return(GRanges())
  }
  
  # Ensure required columns exist
  required_cols <- c("Id1", "Id2", "Start", "End")
  missing_cols <- setdiff(required_cols, colnames(ibd_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Add chromosome if missing
  if (!"Chr" %in% colnames(ibd_df)) {
    ibd_df$Chr <- 1
  }
  
  # Create pair identifier
  ibd_df <- ibd_df %>%
    mutate(pair = paste(pmin(Id1, Id2), pmax(Id1, Id2), sep = "_"))
  
  # Create GRanges
  gr <- GRanges(
    seqnames = paste0("chr", ibd_df$Chr),
    ranges = IRanges(start = ibd_df$Start, end = ibd_df$End),
    pair = ibd_df$pair,
    Id1 = ibd_df$Id1,
    Id2 = ibd_df$Id2
  )
  
  # Set sequence lengths if provided
  if (!is.null(genome_length)) {
    seqlengths(gr) <- setNames(genome_length, unique(as.character(seqnames(gr))))
  }
  
  return(gr)
}

#' Calculate overlap between true and inferred IBD segments
#' 
#' @param true_gr GRanges of true IBD segments
#' @param inferred_gr GRanges of inferred IBD segments
#' @param overlap_threshold Minimum overlap fraction to consider a match
#' @return List with matched and unmatched segments
calculate_segment_overlap <- function(true_gr, inferred_gr, overlap_threshold = 0.5) {
  
  if (length(true_gr) == 0 || length(inferred_gr) == 0) {
    return(list(
      true_positive = 0,
      false_negative = length(true_gr),
      false_positive = length(inferred_gr),
      true_positive_bp = 0,
      false_negative_bp = sum(width(true_gr)),
      false_positive_bp = sum(width(inferred_gr))
    ))
  }
  
  # Find overlaps between true and inferred segments
  overlaps <- findOverlaps(true_gr, inferred_gr)
  
  # Calculate overlap fractions
  overlap_widths <- width(pintersect(true_gr[queryHits(overlaps)], 
                                     inferred_gr[subjectHits(overlaps)]))
  true_widths <- width(true_gr[queryHits(overlaps)])
  inferred_widths <- width(inferred_gr[subjectHits(overlaps)])
  
  # Overlap fraction relative to both segments
  overlap_frac_true <- overlap_widths / true_widths
  overlap_frac_inferred <- overlap_widths / inferred_widths
  
  # Consider match if overlap exceeds threshold for both segments
  is_match <- (overlap_frac_true >= overlap_threshold) & 
    (overlap_frac_inferred >= overlap_threshold)
  
  # Also check if pairs match
  true_pairs <- mcols(true_gr)$pair[queryHits(overlaps)]
  inferred_pairs <- mcols(inferred_gr)$pair[subjectHits(overlaps)]
  pairs_match <- true_pairs == inferred_pairs
  
  # Valid matches must have both overlap and matching pairs
  valid_matches <- is_match & pairs_match
  
  # Identify matched segments
  matched_true <- unique(queryHits(overlaps)[valid_matches])
  matched_inferred <- unique(subjectHits(overlaps)[valid_matches])
  
  # Calculate statistics
  true_positive <- length(matched_true)
  false_negative <- length(true_gr) - true_positive
  false_positive <- length(inferred_gr) - length(matched_inferred)
  
  true_positive_bp <- sum(width(true_gr[matched_true]))
  false_negative_bp <- sum(width(true_gr[-matched_true]))
  false_positive_bp <- sum(width(inferred_gr[-matched_inferred]))
  
  return(list(
    true_positive = true_positive,
    false_negative = false_negative,
    false_positive = false_positive,
    true_positive_bp = true_positive_bp,
    false_negative_bp = false_negative_bp,
    false_positive_bp = false_positive_bp,
    matched_true_idx = matched_true,
    matched_inferred_idx = matched_inferred
  ))
}

#' Calculate IBD detection metrics
#' 
#' @param overlap_results Results from calculate_segment_overlap
#' @return Data frame with sensitivity, precision, FNR, FPR, F1
calculate_metrics <- function(overlap_results) {
  
  tp <- overlap_results$true_positive
  fn <- overlap_results$false_negative
  fp <- overlap_results$false_positive
  
  tp_bp <- overlap_results$true_positive_bp
  fn_bp <- overlap_results$false_negative_bp
  fp_bp <- overlap_results$false_positive_bp
  
  # Segment-based metrics
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA
  f1_score <- if (!is.na(sensitivity) && !is.na(precision) && (sensitivity + precision) > 0) {
    2 * (precision * sensitivity) / (precision + sensitivity)
  } else NA
  
  fnr <- if ((tp + fn) > 0) fn / (tp + fn) else NA
  fpr <- if ((fp + tp) > 0) fp / (fp + tp) else NA
  
  # Base pair-based metrics
  total_true_bp <- tp_bp + fn_bp
  total_inferred_bp <- tp_bp + fp_bp
  
  sensitivity_bp <- if (total_true_bp > 0) tp_bp / total_true_bp else NA
  precision_bp <- if (total_inferred_bp > 0) tp_bp / total_inferred_bp else NA
  f1_score_bp <- if (!is.na(sensitivity_bp) && !is.na(precision_bp) && 
                     (sensitivity_bp + precision_bp) > 0) {
    2 * (precision_bp * sensitivity_bp) / (precision_bp + sensitivity_bp)
  } else NA
  
  fnr_bp <- if (total_true_bp > 0) fn_bp / total_true_bp else NA
  fpr_bp <- if (total_inferred_bp > 0) fp_bp / total_inferred_bp else NA
  
  return(data.frame(
    true_positive = tp,
    false_negative = fn,
    false_positive = fp,
    sensitivity = sensitivity,
    precision = precision,
    fnr = fnr,
    fpr = fpr,
    f1_score = f1_score,
    true_positive_bp = tp_bp,
    false_negative_bp = fn_bp,
    false_positive_bp = fp_bp,
    sensitivity_bp = sensitivity_bp,
    precision_bp = precision_bp,
    fnr_bp = fnr_bp,
    fpr_bp = fpr_bp,
    f1_score_bp = f1_score_bp
  ))
}

#' Load IBD data for a specific run/chromosome/method
#' 
#' @param input_dir Base input directory
#' @param run_id Run number
#' @param rec_rate Recombination rate string
#' @param chr Chromosome number
#' @param method Method name (hmm, iso, ibs) or "true"
#' @return Data frame with IBD segments
load_ibd_data <- function(input_dir, run_id, rec_rate, chr, method = "true") {
  
  # Construct directory path
  dir_path <- file.path(input_dir, sprintf("run%s_rec%s_chr%s", run_id, rec_rate, chr))
  
  if (!dir.exists(dir_path)) {
    warning(sprintf("Directory not found: %s", dir_path))
    return(NULL)
  }
  
  # Construct file path based on method
  if (method == "true") {
    file_path <- file.path(dir_path, "true_ibd_summary.tsv")
    
    # Read and convert to segment format
    if (file.exists(file_path)) {
      ibd_data <- read_tsv(file_path, show_col_types = FALSE)
      # True IBD summary doesn't have Start/End, only total_ibd_bp
      # We need the detailed segment file
      # Look for original tskibd output
      true_ibd_dir <- file.path("out_true_ibd", sprintf("run%s_rec%s_chr%s_ibd", run_id, rec_rate, chr))
      true_ibd_file <- list.files(true_ibd_dir, pattern = "\\.ibd$", full.names = TRUE)
      
      if (length(true_ibd_file) > 0) {
        ibd_data <- read_tsv(true_ibd_file[1], show_col_types = FALSE)
        # tskibd output format: Id1, Id2, Start, End, Tmrca, etc.
      } else {
        warning(sprintf("True IBD segment file not found in: %s", true_ibd_dir))
        return(NULL)
      }
    } else {
      return(NULL)
    }
  } else {
    file_path <- file.path(dir_path, sprintf("inferred_ibd_%s.tsv", method))
    
    if (file.exists(file_path)) {
      ibd_data <- read_tsv(file_path, show_col_types = FALSE)
      # Inferred methods might have different formats
      # Standardize column names if needed
    } else {
      return(NULL)
    }
  }
  
  return(ibd_data)
}

################################################################################
# MAIN ANALYSIS
################################################################################

message("Starting IBD accuracy analysis...")
message(sprintf("Input directory: %s", opt$input_dir))
message(sprintf("Output directory: %s", opt$output_dir))
message(sprintf("Methods to evaluate: %s", paste(methods_to_eval, collapse = ", ")))

# Define recombination rates to analyze
rec_rates <- data.frame(
  run_id = c(1, 2, 3, 4),
  rec_rate = c("1.0e-09", "1.0e-08", "1.0e-07", "1.0e-06"),
  rec_rate_numeric = c(1e-9, 1e-8, 1e-7, 1e-6)
)

# Chromosomes
chromosomes <- c(1, 2, 3)

# Initialize results storage
all_results <- list()
result_counter <- 1

# Loop through all combinations
for (i in 1:nrow(rec_rates)) {
  run_id <- rec_rates$run_id[i]
  rec_rate <- rec_rates$rec_rate[i]
  rec_rate_num <- rec_rates$rec_rate_numeric[i]
  
  message(sprintf("\nProcessing Run %d (rec_rate = %s)...", run_id, rec_rate))
  
  for (chr in chromosomes) {
    message(sprintf("  Chromosome %d...", chr))
    
    # Load true IBD segments
    true_ibd <- load_ibd_data(opt$input_dir, run_id, rec_rate, chr, method = "true")
    
    if (is.null(true_ibd) || nrow(true_ibd) == 0) {
      warning(sprintf("    No true IBD data for run%d_chr%d", run_id, chr))
      next
    }
    
    # Filter by minimum length if specified
    if (opt$min_length > 0) {
      true_ibd <- true_ibd %>%
        mutate(length = End - Start) %>%
        filter(length >= opt$min_length)
    }
    
    # Get genome length from true IBD data
    genome_length <- max(true_ibd$End, na.rm = TRUE)
    
    # Convert to GRanges
    true_gr <- ibd_to_granges(true_ibd, genome_length)
    
    message(sprintf("    True IBD segments: %d", length(true_gr)))
    
    # Evaluate each method
    for (method in methods_to_eval) {
      
      # Load inferred IBD segments
      inferred_ibd <- load_ibd_data(opt$input_dir, run_id, rec_rate, chr, method = method)
      
      if (is.null(inferred_ibd) || nrow(inferred_ibd) == 0) {
        warning(sprintf("    No %s data for run%d_chr%d", method, run_id, chr))
        next
      }
      
      # Filter by minimum length
      if (opt$min_length > 0 && "Start" %in% colnames(inferred_ibd) && 
          "End" %in% colnames(inferred_ibd)) {
        inferred_ibd <- inferred_ibd %>%
          mutate(length = End - Start) %>%
          filter(length >= opt$min_length)
      }
      
      # Convert to GRanges
      inferred_gr <- ibd_to_granges(inferred_ibd, genome_length)
      
      message(sprintf("    %s segments: %d", toupper(method), length(inferred_gr)))
      
      # Calculate overlap
      overlap_results <- calculate_segment_overlap(true_gr, inferred_gr, opt$overlap_threshold)
      
      # Calculate metrics
      metrics <- calculate_metrics(overlap_results)
      
      # Store results
      all_results[[result_counter]] <- data.frame(
        run_id = run_id,
        rec_rate = rec_rate,
        rec_rate_numeric = rec_rate_num,
        chr = chr,
        method = method,
        n_true_segments = length(true_gr),
        n_inferred_segments = length(inferred_gr),
        metrics,
        stringsAsFactors = FALSE
      )
      
      result_counter <- result_counter + 1
      
      message(sprintf("    %s: Sensitivity=%.3f, Precision=%.3f, FNR=%.3f, FPR=%.3f",
                      toupper(method), metrics$sensitivity, metrics$precision, 
                      metrics$fnr, metrics$fpr))
    }
  }
}

# Combine all results
if (length(all_results) == 0) {
  stop("No results generated. Check input files and paths.")
}

results_df <- bind_rows(all_results)

# Save detailed results
output_file <- file.path(opt$output_dir, "ibd_accuracy_detailed.tsv")
write_tsv(results_df, output_file)
message(sprintf("\nDetailed results saved to: %s", output_file))

# Calculate genome-wide averages (across chromosomes)
summary_df <- results_df %>%
  group_by(run_id, rec_rate, rec_rate_numeric, method) %>%
  summarise(
    n_chromosomes = n(),
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_fnr = mean(fnr, na.rm = TRUE),
    mean_fpr = mean(fpr, na.rm = TRUE),
    mean_f1 = mean(f1_score, na.rm = TRUE),
    mean_sensitivity_bp = mean(sensitivity_bp, na.rm = TRUE),
    mean_precision_bp = mean(precision_bp, na.rm = TRUE),
    mean_fnr_bp = mean(fnr_bp, na.rm = TRUE),
    mean_fpr_bp = mean(fpr_bp, na.rm = TRUE),
    mean_f1_bp = mean(f1_score_bp, na.rm = TRUE),
    .groups = 'drop'
  )

# Save summary results
summary_file <- file.path(opt$output_dir, "ibd_accuracy_summary.tsv")
write_tsv(summary_df, summary_file)
message(sprintf("Summary results saved to: %s", summary_file))

################################################################################
# VISUALIZATION
################################################################################

message("\nGenerating plots...")

# Color palette for methods
method_colors <- c("hmm" = "#E41A1C", "iso" = "#377EB8", "ibs" = "#4DAF4A")

# Plot 1: False Negative Rate vs Recombination Rate
p_fnr <- ggplot(summary_df, aes(x = rec_rate_numeric, y = mean_fnr, 
                                color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_color_manual(values = method_colors,
                     labels = c("hmm" = "hmmIBD", "iso" = "isoRelate", "ibs" = "IBS")) +
  labs(
    title = "False Negative Rate vs Recombination Rate",
    subtitle = "Proportion of true IBD segments missed by callers",
    x = "Recombination rate (per bp per generation)",
    y = "False Negative Rate (FNR)",
    color = "Method"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# Plot 2: False Positive Rate vs Recombination Rate
p_fpr <- ggplot(summary_df, aes(x = rec_rate_numeric, y = mean_fpr, 
                                color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_color_manual(values = method_colors,
                     labels = c("hmm" = "hmmIBD", "iso" = "isoRelate", "ibs" = "IBS")) +
  labs(
    title = "False Positive Rate vs Recombination Rate",
    subtitle = "Proportion of called IBD segments that are false",
    x = "Recombination rate (per bp per generation)",
    y = "False Positive Rate (FPR)",
    color = "Method"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# Combine plots
combined_plot <- plot_grid(p_fnr, p_fpr, ncol = 1, align = "v", labels = c("A", "B"))

# Save plots
ggsave(file.path(opt$output_dir, "ibd_accuracy_fnr_fpr.png"), 
       combined_plot, width = 10, height = 10, dpi = 300)
ggsave(file.path(opt$output_dir, "ibd_accuracy_fnr_fpr.pdf"), 
       combined_plot, width = 10, height = 10)

message(sprintf("Combined plot saved to: %s", 
                file.path(opt$output_dir, "ibd_accuracy_fnr_fpr.png")))

# Plot 3: Sensitivity and Precision
p_sens_prec <- ggplot(summary_df, aes(x = rec_rate_numeric, color = method, group = method)) +
  geom_line(aes(y = mean_sensitivity, linetype = "Sensitivity"), size = 1) +
  geom_point(aes(y = mean_sensitivity), size = 3) +
  geom_line(aes(y = mean_precision, linetype = "Precision"), size = 1) +
  geom_point(aes(y = mean_precision), size = 3, shape = 17) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = method_colors,
                     labels = c("hmm" = "hmmIBD", "iso" = "isoRelate", "ibs" = "IBS")) +
  scale_linetype_manual(values = c("Sensitivity" = "solid", "Precision" = "dashed")) +
  labs(
    title = "Sensitivity and Precision vs Recombination Rate",
    x = "Recombination rate (per bp per generation)",
    y = "Rate",
    color = "Method",
    linetype = "Metric"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave(file.path(opt$output_dir, "ibd_accuracy_sens_prec.png"), 
       p_sens_prec, width = 10, height = 6, dpi = 300)

message(sprintf("Sensitivity/Precision plot saved to: %s",
                file.path(opt$output_dir, "ibd_accuracy_sens_prec.png")))

# Plot 4: F1 Score
p_f1 <- ggplot(summary_df, aes(x = rec_rate_numeric, y = mean_f1, 
                               color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_color_manual(values = method_colors,
                     labels = c("hmm" = "hmmIBD", "iso" = "isoRelate", "ibs" = "IBS")) +
  labs(
    title = "F1 Score vs Recombination Rate",
    subtitle = "Harmonic mean of precision and sensitivity",
    x = "Recombination rate (per bp per generation)",
    y = "F1 Score",
    color = "Method"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

ggsave(file.path(opt$output_dir, "ibd_accuracy_f1.png"), 
       p_f1, width = 8, height = 6, dpi = 300)

message(sprintf("F1 score plot saved to: %s",
                file.path(opt$output_dir, "ibd_accuracy_f1.png")))

################################################################################
# PRINT SUMMARY
################################################################################

message("\n" , paste(rep("=", 80), collapse = ""))
message("IBD ACCURACY ANALYSIS SUMMARY")
message(paste(rep("=", 80), collapse = ""))

for (method in unique(summary_df$method)) {
  method_data <- summary_df %>% filter(method == !!method)
  
  message(sprintf("\n%s:", toupper(method)))
  message(sprintf("  Rec Rate          FNR       FPR    Sensitivity  Precision   F1 Score"))
  message(paste(rep("-", 80), collapse = ""))
  
  for (i in 1:nrow(method_data)) {
    message(sprintf("  %-15s  %6.2f%%  %6.2f%%     %6.2f%%     %6.2f%%    %6.2f%%",
                    method_data$rec_rate[i],
                    method_data$mean_fnr[i] * 100,
                    method_data$mean_fpr[i] * 100,
                    method_data$mean_sensitivity[i] * 100,
                    method_data$mean_precision[i] * 100,
                    method_data$mean_f1[i] * 100))
  }
}

message("\n" , paste(rep("=", 80), collapse = ""))
message("Analysis complete!")
message(paste(rep("=", 80), collapse = ""))
message(sprintf("\nAll results saved to: %s", opt$output_dir))
message("Files generated:")
message(sprintf("  - %s", "ibd_accuracy_detailed.tsv"))
message(sprintf("  - %s", "ibd_accuracy_summary.tsv"))
message(sprintf("  - %s", "ibd_accuracy_fnr_fpr.png"))
message(sprintf("  - %s", "ibd_accuracy_sens_prec.png"))
message(sprintf("  - %s", "ibd_accuracy_f1.png"))