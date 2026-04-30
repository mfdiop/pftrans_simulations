

library(tidyverse)
library(patchwork)
library(gt)
library(broom)

all_results <- readRDS("simulations/single_run/evaluation_results.rds")
out_prefix = "simulations/single_run/figures/"

# Helper function to calculate metrics from confusion matrix
calculate_confusion_metrics <- function(conf_matrix) {
  
  # Ensure it's a matrix
  if (!is.matrix(conf_matrix)) {
    conf_matrix <- as.matrix(conf_matrix)
  }
  
  n_classes <- nrow(conf_matrix)
  metrics <- list()
  
  # Overall metrics
  metrics$Accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  
  # Precision, Recall, F1 for each class
  class_wise <- matrix(NA, nrow = n_classes, ncol = 3,
                       dimnames = list(rownames(conf_matrix),
                                       c("Precision", "Recall", "F1")))
  
  for (i in 1:n_classes) {
    TP <- conf_matrix[i, i]
    FP <- sum(conf_matrix[, i]) - TP
    FN <- sum(conf_matrix[i, ]) - TP
    
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    Recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    F1 <- ifelse((Precision + Recall) > 0, 
                 2 * Precision * Recall / (Precision + Recall), 0)
    
    class_wise[i, ] <- c(Precision, Recall, F1)
  }
  
  # Macro-averaged metrics
  metrics$Macro_Precision <- mean(class_wise[, "Precision"], na.rm = TRUE)
  metrics$Macro_Recall <- mean(class_wise[, "Recall"], na.rm = TRUE)
  metrics$Macro_F1 <- mean(class_wise[, "F1"], na.rm = TRUE)
  
  # Weighted metrics (by class prevalence)
  class_weights <- rowSums(conf_matrix) / sum(conf_matrix)
  metrics$Weighted_Precision <- sum(class_wise[, "Precision"] * class_weights)
  metrics$Weighted_Recall <- sum(class_wise[, "Recall"] * class_weights)
  metrics$Weighted_F1 <- sum(class_wise[, "F1"] * class_weights)
  
  return(list(
    metrics = metrics,
    class_wise = class_wise
  ))
}

# 1. Basic Aggregation - Summing Across Replicates

plot_aggregated_confusion <- function(results_list, 
                                      metric_name = c("optimal_youden", "fixed", "median"),
                                      matrix_type = c("IBD", "IBS", "phylo")) {
  
  cat("Extracting confusion matrices for:", matrix_type, "-", metric_name, "\n")
  
  # Extract confusion matrices for the specific matrix type
  confusion_matrices <- map(results_list, function(replicate) {
    # Navigate through the nested structure
    if (!is.null(replicate$evaluation[[metric_name]])) {
      # Check structure - it might be a list with IBD, IBS, phylo
      eval_data <- replicate$evaluation[[metric_name]]
      
      # Option 1: eval_data is a list with matrix_type as element
      if (matrix_type %in% names(eval_data)) {
        return(eval_data[[matrix_type]]$confusion_matrices)
      }
      # Option 2: eval_data contains confusion_matrices as a list
      else if (!is.null(eval_data$confusion_matrices)) {
        return(eval_data$confusion_matrices[[matrix_type]])
      }
      # Option 3: eval_data is the confusion matrix itself
      else if (is.matrix(eval_data)) {
        return(eval_data)
      }
    }
    return(NULL)
  })
  
  # Filter out NULLs
  confusion_matrices <- compact(confusion_matrices)
  
  if (length(confusion_matrices) == 0) {
    # Try alternative structure
    cat("Trying alternative structure search...\n")
    confusion_matrices <- map(results_list, function(replicate) {
      # Deep search for confusion matrix
      if (!is.null(replicate$evaluation)) {
        # Check if there's a nested structure
        for (metric in names(replicate$evaluation)) {
          if (grepl(metric_name, metric, ignore.case = TRUE)) {
            metric_data <- replicate$evaluation[[metric]]
            if (!is.null(metric_data[[matrix_type]])) {
              if (is.list(metric_data[[matrix_type]]) && 
                  !is.null(metric_data[[matrix_type]]$confusion_matrix)) {
                return(metric_data[[matrix_type]]$confusion_matrix)
              } else if (is.matrix(metric_data[[matrix_type]])) {
                return(metric_data[[matrix_type]])
              }
            }
          }
        }
      }
      return(NULL)
    })
    confusion_matrices <- compact(confusion_matrices)
  }
  
  if (length(confusion_matrices) == 0) {
    warning("No confusion matrices found for ", matrix_type, " - ", metric_name)
    
    # Debug: Show available structure
    cat("\nDebug - Available structure in first replicate:\n")
    if (!is.null(results_list[[1]]$evaluation)) {
      cat("Evaluation metrics available:\n")
      print(names(results_list[[1]]$evaluation))
      
      if (!is.null(results_list[[1]]$evaluation[[metric_name]])) {
        cat("\nStructure of", metric_name, ":\n")
        print(str(results_list[[1]]$evaluation[[metric_name]]))
      }
    }
    return(NULL)
  }
  
  cat("Found", length(confusion_matrices), "confusion matrices for aggregation\n")
  
  # Check that all matrices have same dimensions
  dims <- map(confusion_matrices, dim)
  if (length(unique(dims)) > 1) {
    warning("Confusion matrices have different dimensions. Using first matrix dimensions.")
    # Use first matrix dimensions as reference
    ref_dims <- dims[[1]]
    # Resize matrices if needed (simplistic approach)
    confusion_matrices <- map(confusion_matrices, function(cm) {
      if (!all(dim(cm) == ref_dims)) {
        # Create empty matrix with correct dimensions
        new_cm <- matrix(0, nrow = ref_dims[1], ncol = ref_dims[2])
        rownames(new_cm) <- paste0("Class", 1:ref_dims[1])
        colnames(new_cm) <- paste0("Class", 1:ref_dims[2])
        
        # Fill with available data
        common_rows <- intersect(rownames(cm), rownames(new_cm))
        common_cols <- intersect(colnames(cm), colnames(new_cm))
        if (length(common_rows) > 0 && length(common_cols) > 0) {
          new_cm[common_rows, common_cols] <- cm[common_rows, common_cols]
        }
        return(new_cm)
      }
      return(cm)
    })
  }
  
  # Sum across replicates
  aggregated <- Reduce(`+`, confusion_matrices)
  
  # Convert to tidy format for plotting
  # Convert matrix to tidy format
  # If you have broom package

  agg_tidy <- tibble::as_tibble(aggregated) %>%
    rename(Count = n) %>%
    mutate(
      Actual = factor(Actual, levels = rownames(aggregated)),
      Predicted = factor(Predicted, levels = colnames(aggregated))
    )
  
  # Calculate percentages
  total <- sum(agg_tidy$Count)
  agg_tidy <- agg_tidy %>%
    mutate(Percentage = Count / total * 100)
  
  # Plot 1: Heatmap with counts
  p1 <- ggplot(agg_tidy, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count), color = "black", size = 6) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      title = paste("Aggregated Confusion Matrix -", matrix_type), # metric_name
      # subtitle = paste("Sum across", length(confusion_matrices), "replicates"),
      x = "Predicted Class", y = "Actual Class") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, colour = "black", face = "bold"),
      plot.title.position = "plot",
      axis.title = element_text(size = 14, colour = "black", face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      legend.position = "none") +
    coord_fixed()
  
  # Plot 2: Heatmap with percentages
  p2 <- ggplot(agg_tidy, aes(x = Predicted, y = Actual, fill = Percentage)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
              color = "black", size = 5) +
    scale_fill_gradient(low = "white", high = "darkred", 
                        limits = c(0, 100)) +
    labs(
      title = paste("Percentage Confusion Matrix -", matrix_type),
      # subtitle = "Percentage of total cases",
      x = "Predicted Class",
      y = "Actual Class"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, colour = "black", face = "bold"),
      plot.title.position = "plot",
      axis.title = element_text(size = 14, colour = "black", face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      legend.position = "none") +
    coord_fixed()
  
  # Calculate performance metrics
  performance <- calculate_confusion_metrics(aggregated)
  
  # Plot 3: Performance metrics bar plot
  perf_df <- data.frame(
    Metric = names(performance$metrics),
    Value = unlist(performance$metrics))
  
  # Remove any non-numeric metrics if present (like "Confusion Matrix")
  perf_df <- perf_df[!is.na(as.numeric(perf_df$Value)), ]
  
  # Convert Value to numeric
  perf_df$Value <- as.numeric(perf_df$Value)
  
  # Reorder metrics for better visualization
  desired_order <- c("Accuracy", "Macro_Precision", "Macro_Recall", "Macro_F1", 
                     "Weighted_Precision", "Weighted_Recall", "Weighted_F1")
  
  perf_df$Metric <- factor(perf_df$Metric, levels = desired_order)
  
  # Create the plot
  p3 <- ggplot(perf_df, aes(x = Metric, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", Value)), 
              vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Performance Metrics Summary",
      x = "", y = "Value") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()) +
    ylim(0, min(1.1, max(perf_df$Value) * 1.2))
  

  p3_1 <- ggplot(perf_df, aes(x = Metric, y = Value, fill = Metric)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.3f", Value)), 
              vjust = -0.5, size = 4, fontface = "bold") +
    labs(
      title = "Performance Metrics",
      x = "",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    ylim(0, 1.1)
  
  # Alternative: Horizontal bar plot (often better for many metrics)
  p5 <- ggplot(perf_df, aes(x = Value, y = reorder(Metric, Value), fill = Value)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", Value)), 
              hjust = -0.1, size = 4, fontface = "bold") +
    scale_fill_gradient(low = "#4B9CD3", high = "#13294B", limits = c(0, 1)) +
    labs(
      x = "", y = "") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 12),
      panel.grid.major.y = element_blank()
    ) +
    xlim(0, min(1.1, max(perf_df$Value) * 1.2))
  
  # # Alternative: Radar/spider plot for a different visualization
  # # (Requires fmsb or ggradar package)
  # if (requireNamespace("fmsb", quietly = TRUE)) {
  #   library(fmsb)
  #   
  #   # Prepare data for radar chart
  #   radar_data <- rbind(rep(1, nrow(perf_df)),  # Max value
  #                       rep(0, nrow(perf_df)),  # Min value
  #                       perf_df$Value)
  #   colnames(radar_data) <- perf_df$Metric
  #   
  #   # Create radar chart
  #   radarchart(
  #     radar_data,
  #     axistype = 1,
  #     pcol = "#13294B", 
  #     pfcol = rgb(0.2, 0.5, 0.8, 0.3),
  #     plwd = 2,
  #     cglcol = "grey",
  #     cglty = 1,
  #     axislabcol = "grey",
  #     caxislabels = seq(0, 1, 0.25),
  #     cglwd = 0.8,
  #     vlcex = 0.8,
  #     title = "Performance Metrics Radar Chart"
  #   )
  # }
  
  # Plot class-wise metrics (if available)
  if (!is.null(performance$class_wise)) {
    # Assuming class_wise is a matrix or data.frame
    class_metrics <- as.data.frame(performance$class_wise)
    class_metrics$Class <- rownames(class_metrics)
    
    # Convert to long format
    class_long <- reshape2::melt(class_metrics, id.vars = "Class", 
                                 variable.name = "Metric", value.name = "Value")
    
    # Or using tidyr if available:
    # class_long <- tidyr::pivot_longer(class_metrics, -Class, 
    #                                   names_to = "Metric", values_to = "Value")
    
    p4 <- ggplot(class_long, aes(x = Class, y = Value, fill = Metric)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = sprintf("%.2f", Value)), 
                position = position_dodge(width = 0.8), 
                vjust = -0.5, size = 3) +
      labs(
        title = "Class-wise Performance Metrics",
        x = "Class",
        y = "Value",
        fill = "Metric"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      ylim(0, min(1.1, max(class_long$Value, na.rm = TRUE) * 1.2))
    
    # print(p4)
  }
  
  # Combine plots
  # combined_plot <- (p1 | p2) / (p3 | p4)
  
  # Create summary table
  summary_table <- tibble(
    Metric = names(performance$metrics),
    Value = as.numeric(performance$metrics)
  ) %>%
    filter(!Metric %in% c("Confusion Matrix", "Class Metrics")) %>%
    gt() %>%
    tab_header(
      title = paste("Performance Summary -", matrix_type)
      # subtitle = paste("Aggregated across", length(confusion_matrices), "replicates")
    ) %>%
    fmt_number(columns = Value, decimals = 3) %>%
    cols_label(
      Metric = "Performance Metric",
      Value = "Value")
  
  print(summary_table)
  
  return(list(
    aggregated_matrix = aggregated,
    tidy_data = agg_tidy,
    performance = performance,
    plots = list(
      count_heatmap = p1,
      percent_heatmap = p2,
      metrics_bar = p3,
      metrics_bar_1 = p3_1,
      class_metrics = p4,
      horizontal_bar = p5), # ,combined = combined_plot
    summary_table = summary_table
  ))
}


# Aggregate All Matrix Types at Once
plot_all_confusion_matrices <- function(results_list, metric_name = c("optimal_youden", "fixed", "median")) {
  
  matrix_types <- c("IBD", "IBS", "phylo")
  all_results <- list()
  
  for (mt in matrix_types) {
    cat("\nProcessing:", mt, "...\n")
    result <- plot_aggregated_confusion(results_list, metric_name, mt)
    if (!is.null(result)) {
      all_results[[mt]] <- result
    }
  }
  
  # Create comparison plot
  if (length(all_results) > 0) {
    # Combine aggregated matrices
    agg_matrices <- map(all_results, "aggregated_matrix")
    
    # Create a comparison heatmap
    comparison_data <- map_dfr(names(agg_matrices), function(mt) {
      df <- as.data.frame(as.table(agg_matrices[[mt]]))
      colnames(df) <- c("Actual", "Predicted", "Count")
      total <- sum(df$Count)
      df$percentage <- round((df$Count / total) * 100, 1)
      df$Matrix_Type <- mt
      return(df)
    })
    
    # Plot comparison count
    p_compare <- ggplot(comparison_data, 
                        aes(x = Predicted, y = Actual, fill = Count)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.0f", Count)), 
                vjust = -0.5, size = 8, fontface = "bold") +
      facet_wrap(~ Matrix_Type, ncol = 3) +
      scale_fill_gradient(low = "white", high = "darkorange") +
      theme_minimal() +
      theme(legend.position = "none",
        axis.title = element_text(size = 16, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        strip.text = element_text(size = 18,, colour = "black", face = "bold")
        )
    
    # Plot comparison proportion
    p_prop <- ggplot(comparison_data, 
                        aes(x = Predicted, y = Actual, fill = percentage)) +
      geom_tile() +
      geom_text(aes(label = percentage), 
                vjust = -0.5, size = 8, fontface = "bold") +
      facet_wrap(~ Matrix_Type, ncol = 3) +
      scale_fill_gradient(low = "white", high = "darkred") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title = element_text(size = 16, colour = "black", face = "bold"),
            axis.text.x = element_text(size = 13, angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(size = 13, colour = "black"),
            strip.text = element_text(size = 18, colour = "black", face = "bold")
      )
    

    # Compare performance metrics
    perf_comparison <- map_dfr(names(all_results), function(mt) {
      perf <- all_results[[mt]]$performance$metrics
      tibble(
        Matrix_Type = mt,
        Accuracy = perf$Accuracy,
        Macro_F1 = perf$Macro_F1,
        Weighted_F1 = perf$Weighted_F1
      )
    })
    
    p_perf <- perf_comparison %>%
      pivot_longer(-Matrix_Type, names_to = "Metric", values_to = "Value") %>%
      ggplot(aes(x = Matrix_Type, y = Value, fill = Metric)) +
      geom_col(position = position_dodge()) +
      geom_text(aes(label = sprintf("%.3f", Value)), 
                position = position_dodge(width = 0.9),
                vjust = -0.5, fontface = "bold") +
      labs(x = "", y = "", title = "Performance Comparison") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 18, colour = "black", face = "bold"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text = element_text(size = 12, colour = "black")
      )
    
  }
  all_results$p_compare <- p_compare
  all_results$p_prop <- p_prop
  all_results$p_perf <- p_perf
  
  return(all_results)
}

# Then use the appropriate function
# 1. Basic aggregated confusion matrix
confusion_results <- plot_aggregated_confusion(all_results, "optimal_youden", "IBD")

# Option 1: Plot specific matrix type
confusion_ibd <- plot_aggregated_confusion(all_results, "optimal_youden", "IBD")
confusion_ibs <- plot_aggregated_confusion(all_results, "optimal_youden", "IBS")
confusion_phylo <- plot_aggregated_confusion(all_results, "optimal_youden", "phylo")

# Option 2: Plot all matrix types together
all_confusion <- plot_all_confusion_matrices(all_results, "optimal_youden")

# Save plot output confusion matrix
saveRDS(all_confusion, file.path(dirname(out_prefix), "aggregated_confusion_matrix.rds"))

ggsave(filename = paste0(out_prefix, "aggregated_confusion_matrix.png"),
      plot = all_confusion$p_compare,
      width = 12, height = 10, dpi = 600
       )

ggsave(filename = paste0(out_prefix, "aggregated_confusion_prop.png"),
       plot = all_confusion$p_prop,
       width = 12, height = 10, dpi = 600
)

ggsave(filename = paste0(out_prefix, "performance_matrix.png"),
       plot = all_confusion$p_perf,
       width = 12, height = 10, dpi = 600
)


ggsave(filename = paste0(out_prefix, "ibd_performance_class.png"),
       plot = all_confusion$IBD$plots$class_metrics,
       width = 12, height = 10, dpi = 600
)

ggsave(filename = paste0(out_prefix, "ibs_performance_class.png"),
       plot = all_confusion$IBS$plots$class_metrics,
       width = 12, height = 10, dpi = 600
)

ggsave(filename = paste0(out_prefix, "phylo_performance_class.png"),
       plot = all_confusion$phylo$plots$class_metrics,
       width = 12, height = 10, dpi = 600
)

# Save Summary Tables
write.csv(list(IBD = all_confusion$IBD$summary_table, 
               IBS = all_confusion$IBS$summary_table,
               Patristic_distance = all_confusion$phylo$summary_table),
          file = file.path(dirname(out_prefix), "tables/performance_summary.csv"),
          row.names = FALSE,
          na = "")

write_performance_summary <- function(results_list, out_file) {
  # Ensure the writexl package is available
  if (!requireNamespace("writexl", quietly = TRUE)) {
    install.packages("writexl")
  }
  
  # Create a named list for Excel sheets.
  # It filters for objects that are data frames and have a 'summary_table'.
  sheet_list <- list()
  
  # Loop through each method in your results (e.g., 'IBD', 'IBS', 'phylo')
  for (method_name in names(results_list)[1:3]) {
    method_obj <- results_list[[method_name]]
    
    # Check if this method object has a 'summary_table' that is a data frame
    if (!is.null(method_obj$summary_table)) { #  && is.data.frame(method_obj$summary_table)
      # Use the method name (e.g., "IBD") as the sheet name
      sheet_list[[method_name]] <- as.data.frame(method_obj$summary_table)
    }
  }
  
  # Check if we have anything to write
  if (length(sheet_list) == 0) {
    warning("No valid summary_table data frames found in the provided list.")
    return(invisible(NULL))
  }
  
  # Ensure the output directory exists
  out_dir <- dirname(out_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Write to Excel
  writexl::write_xlsx(sheet_list, path = out_file)
  message("File successfully written to: ", out_file)
}

# How to use the function with your 'all_confusion' object:
# First, specify your output file path (adjust 'out_prefix' as needed)
output_excel_path <- file.path(dirname(out_prefix), "tables", "performance_summary.xlsx")

# Then call the function
write_performance_summary(all_confusion, output_excel_path)

# # Option 3: Quick individual plots
# plot_confusion_by_type(all_results, "optimal_youden", "IBD")
# 
# # 2. Compare multiple metrics
# comparison <- compare_confusion_metrics(
#   results_list,
#   c("optimal_youden", "optimal_f1", "default_threshold")
# )
# 
# # 3. Create complete dashboard
# dashboard <- create_confusion_dashboard(results_list, "optimal_youden", 
#                                         save_plots = TRUE)
# 
# # 4. Interactive version
# interactive_plot <- plot_interactive_confusion(results_list, "optimal_youden")
# 
# # 5. Extract key statistics for reporting
# key_stats <- list(
#   accuracy = confusion_results$performance$metrics$Accuracy,
#   macro_f1 = confusion_results$performance$metrics$Macro_F1,
#   n_replicates = length(results_list),
#   best_class = names(which.max(rowSums(confusion_results$aggregated_matrix)))
# )


