

# Scatter: IBD_prop vs IBS_prop or vs phylo_sim
plot_scatter_pairs <- function(df,
                               x = "IBD_prop",
                               y = "IBS_prop",
                               color_by = c("is_true", "time_diff"),
                               title = NULL,
                               alpha = 0.6,
                               point_size = 1.2) {
  color_by <- match.arg(color_by)
  p <- ggplot(df, aes_string(x = x, y = y)) +
    geom_point(aes_string(color = color_by), alpha = alpha, size = point_size) +
    theme_minimal() +
    labs(title = title, x = x, y = y, color = color_by)
  
  # sensible coloring for binary is_true
  if (color_by == "is_true") {
    p <- p + scale_color_manual(values = c("0" = "grey70", "1" = "red"),
                                labels = c("0" = "not_true", "1" = "true"))
  } else {
    p <- p + scale_color_viridis_c()
  }
  return(p)
}

# Hexbin joint density (good for dense scatter)
plot_hex_density <- function(df, x = "IBD_prop", y = "phylo_sim_exp", bins = 40, title = NULL) {
  ggplot(df, aes_string(x = x, y = y)) +
    stat_binhex(bins = bins) +
    scale_fill_gradient(low = "lightyellow", high = "darkred", trans = "log10") +
    theme_minimal() +
    labs(title = title, x = x, y = y, fill = "count")
}


# Plot confusion matrices for comparison
plot_confusion_matrices <- function(confusion_results, output = "oy_confusion_metrics") {
  cm_list <- confusion_results$confusion_matrices
  plots <- list()
  
  for (metric_name in names(cm_list)) {
    cm <- as.data.frame(cm_list[[metric_name]])
    
    p <- ggplot(cm, aes(x = Predicted, y = Actual, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), vjust = 0.5, size = 7, color = 'black') +
      # scale_fill_gradient(low = "white", high = "steelblue") +
      scale_fill_gradient(low = "white", high = "darkred") + # low = "orange", high = "#009194"
      labs(
        title = paste("Confusion Matrix:", metric_name),
        subtitle = paste("Threshold:", round(confusion_results$thresholds_used[metric_name], 4))
      ) +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 18, face = 'bold', colour = 'black'),
            plot.subtitle = element_text(size = 14, face = 'bold'),
            axis.title = element_text(size = 16, face = 'bold', colour = 'black'),
            axis.text = element_text(size = 14, colour = 'black'))
    
    plots[[metric_name]] <- p
  }
  
  # Save as PDF
  pdf(paste0(output, ".pdf"), width = 22, height = 20)
  # Arrange all plots in a grid
  do.call(gridExtra::grid.arrange, c(plots, ncol = length(plots)))
  dev.off()
  
  # Save as PNG
  png(paste0(output, ".png"), width = 1200, height = 800, res = 300)
  do.call(gridExtra::grid.arrange, c(plots, ncol = length(plots)))
  dev.off()
  
  # Save with specific dimensions
  ggsave(paste0(output, "_v1.pdf"), 
         plot = do.call(gridExtra::grid.arrange, c(plots, ncol = length(plots))),
         width = 22, height = 20)
}


# Combined explorer for a scenario (IBD vs IBS and IBD vs phylo_sim)
plot_exploratory_for_scenario <- function(df, phylo_col = "phylo_sim_exp", scenario_name = NULL) {
  p1 <- plot_scatter_pairs(df, x = "IBD_prop", y = "IBS_prop", color_by = "is_true",
                           title = paste0("IBD_prop vs IBS_prop: ", scenario_name))
  p2 <- plot_scatter_pairs(df, x = "IBD_prop", y = phylo_col, color_by = "time_diff",
                           title = paste0("IBD_prop vs phylo_sim (", phylo_col, "): ", scenario_name))
  p3 <- plot_hex_density(df, x = "IBD_prop", y = "IBS_prop", title = paste0("Hex: IBD vs IBS: ", scenario_name))
  p4 <- plot_hex_density(df, x = "IBD_prop", y = phylo_col, title = paste0("Hex: IBD vs phylo_sim: ", scenario_name))
  
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}




# Simple plotting of performance vs scenario (scenario can encode recomb/mut/sampling)
plot_performance_vs_scenario <- function(summary_df, metric = c("auc","aupr"), predictor = NULL) {
  metric <- match.arg(metric)
  dfp <- summary_df %>% ungroup()
  p <- ggplot(dfp, aes(x = scenario, y = .data[[metric]], color = predictor, group = predictor)) +
    geom_point() + geom_line() + theme_minimal() +
    labs(title = paste0(toupper(metric), " across scenarios"), x = "scenario", y = metric) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

# If you have a scenario-level table linking scenario name -> parameter values, join it
# Then plot AUC/AUPR as function of numeric parameter using geom_smooth or geom_point
plot_metric_vs_param <- function(summary_df, scenario_meta, metric = "aupr", param = "sampling_frac") {
  df <- summary_df %>% left_join(scenario_meta, by = "scenario")
  ggplot(df, aes_string(x = param, y = metric, color = "predictor")) +
    geom_point() + geom_smooth(method = "loess", se = TRUE) +
    theme_minimal() +
    labs(title = paste0(metric, " vs ", param), x = param, y = metric)
}



### Scatter plots (IBD vs IBS, IBD vs phylogenetic similarity)

plot_scatter <- function(df, output) {
  p1 <- ggplot(df, aes(IBD, IBS, color = factor(true_link))) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(color = "Transmission") +
    theme(
          axis.title = element_text(size = 18, color = 'black', face = "bold"),
          axis.text = element_text(size = 14, color = 'black')
          )
  
  p2 <- ggplot(df, aes(IBD, phylo, color = factor(true_link))) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(color = "Transmission") +
    theme(
      axis.title = element_text(size = 18, color = 'black', face = "bold"),
      axis.text = element_text(size = 14, color = 'black')
    )
  
  p3 <- ggplot(df, aes(IBS, phylo, color = factor(true_link))) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(color = "Transmission") +
    theme(
      axis.title = element_text(size = 18, color = 'black', face = "bold"),
      axis.text = element_text(size = 14, color = 'black')
    )
  
  combined <- ggpubr::ggarrange(p1, p2, p3, ncol = 3, common.legend = TRUE)
  
  ggsave(output, plot = combined, width = 20, height = 12, dpi = 600)
}


### Hexbin / joint density

plot_hex <- function(df, output = "out.png") {
  hex <- ggplot(df, aes(IBD, IBS)) +
    geom_hex(bins = 40) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18, color = 'black', face = "bold"),
      axis.text = element_text(size = 14, color = 'black'))
  
  ggsave(output, plot = hex, width = 20, height = 12, dpi = 600)
}
