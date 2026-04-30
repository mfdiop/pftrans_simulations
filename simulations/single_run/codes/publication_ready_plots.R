
# These plots provide:
#   
#   Multi-dimensional comparisons (radar charts)
# 
#   Uncertainty visualization (confidence bands)
# 
#   Network perspectives (relationship structures)
# 
#   Trade-off analyses (practical considerations)
# 
#   Aggregate patterns (UpSet plots)

create_performance_dashboard <- function(comprehensive_results, output_file = "performance_dashboard.pdf") {
  
  # Create a comprehensive 2x2 dashboard
  p1 <- ggplot(comprehensive_results$performance_metrics) +
    geom_violin(aes(x = "IBD", y = IBD_AUC), fill = "#1f77b4", alpha = 0.7) +
    geom_violin(aes(x = "IBS", y = IBS_AUC), fill = "#ff7f0e", alpha = 0.7) +
    geom_violin(aes(x = "Phylo", y = Phylo_AUC), fill = "#2ca02c", alpha = 0.7) +
    geom_point(aes(x = "IBD", y = IBD_AUC), position = position_jitter(width = 0.2), alpha = 0.5) +
    geom_point(aes(x = "IBS", y = IBS_AUC), position = position_jitter(width = 0.2), alpha = 0.5) +
    geom_point(aes(x = "Phylo", y = Phylo_AUC), position = position_jitter(width = 0.2), alpha = 0.5) +
    labs(title = "AUC Distribution Across Replicates", x = "Method", y = "AUC") +
    theme_minimal()
  
  p2 <- ggplot(comprehensive_results$performance_metrics) +
    geom_density(aes(x = IBD_AUC, fill = "IBD"), alpha = 0.6) +
    geom_density(aes(x = IBS_AUC, fill = "IBS"), alpha = 0.6) +
    geom_density(aes(x = Phylo_AUC, fill = "Phylo"), alpha = 0.6) +
    labs(title = "Performance Density Distributions", x = "AUC", y = "Density") +
    scale_fill_manual(values = c("IBD" = "#1f77b4", "IBS" = "#ff7f0e", "Phylo" = "#2ca02c")) +
    theme_minimal()
  
  # Rank plot
  rank_data <- comprehensive_results$performance_metrics %>%
    mutate(replicate = row_number()) %>%
    pivot_longer(cols = c(IBD_AUC, IBS_AUC, Phylo_AUC), 
                 names_to = "method", values_to = "auc") %>%
    group_by(replicate) %>%
    mutate(rank = rank(-auc)) %>%
    ungroup()
  
  p3 <- ggplot(rank_data, aes(x = method, y = rank, fill = method)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.9) +
    scale_y_reverse() +  # So rank 1 is at top
    labs(title = "Method Ranking Distribution", x = "Method", y = "Rank (1 = Best)") +
    scale_fill_manual(values = c("IBD_AUC" = "#1f77b4", "IBS_AUC" = "#ff7f0e", "Phylo_AUC" = "#2ca02c")) +
    theme_minimal()
  
  # Performance trajectory across replicates
  p4 <- ggplot(comprehensive_results$performance_metrics) +
    geom_line(aes(x = as.numeric(replicate_id), y = IBD_AUC, color = "IBD"), size = 1) +
    geom_line(aes(x = as.numeric(replicate_id), y = IBS_AUC, color = "IBS"), size = 1) +
    geom_line(aes(x = as.numeric(replicate_id), y = Phylo_AUC, color = "Phylo"), size = 1) +
    geom_point(aes(x = as.numeric(replicate_id), y = IBD_AUC, color = "IBD"), size = 2) +
    geom_point(aes(x = as.numeric(replicate_id), y = IBS_AUC, color = "IBS"), size = 2) +
    geom_point(aes(x = as.numeric(replicate_id), y = Phylo_AUC, color = "Phylo"), size = 2) +
    labs(title = "Performance Trajectory Across Replicates", 
         x = "Replicate", y = "AUC", color = "Method") +
    scale_color_manual(values = c("IBD" = "#1f77b4", "IBS" = "#ff7f0e", "Phylo" = "#2ca02c")) +
    theme_minimal()
  
  # Combine all plots
  dashboard <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  ggsave(output_file, dashboard, width = 16, height = 12, dpi = 300)
  
  return(dashboard)
}


create_relationship_network <- function(merged_data, method = "IBD", threshold = 0.1, 
                                        output_file = "relationship_network.pdf") {
  
  library(igraph)
  library(ggraph)
  
  # Filter significant relationships
  network_data <- merged_data %>%
    filter(.data[[method]] > threshold) %>%
    select(id1, id2, weight = .data[[method]])
  
  # Create graph object
  g <- graph_from_data_frame(network_data, directed = FALSE)
  
  # Community detection
  communities <- cluster_louvain(g)
  
  # Create network plot
  network_plot <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(alpha = weight, width = weight), 
                   color = "gray70", show.legend = FALSE) +
    geom_node_point(aes(color = as.factor(communities$membership)), 
                    size = 3, show.legend = FALSE) +
    geom_node_text(aes(label = name), size = 2, repel = TRUE) +
    scale_edge_width_continuous(range = c(0.5, 2)) +
    scale_color_viridis_d() +
    labs(title = paste("Relationship Network:", method),
         subtitle = paste("Threshold =", threshold)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(output_file, network_plot, width = 10, height = 8, dpi = 300)
  return(network_plot)
}


create_method_agreement_heatmap <- function(results_list, output_file = "method_agreement_heatmap.pdf") {
  
  # Calculate agreement matrices between methods across replicates
  agreement_data <- map_dfr(results_list, function(repl) {
    correlations <- repl$correlations
    
    tibble(
      replicate = which(results_list == repl),
      IBD_IBS_cor = correlations$pearson$estimate,
      IBD_Phylo_cor = cor(repl$merged_data$conservative$IBD, 
                          repl$merged_data$conservative$phylo, use = "complete.obs"),
      IBS_Phylo_cor = cor(repl$merged_data$conservative$IBS, 
                          repl$merged_data$conservative$phylo, use = "complete.obs")
    )
  })
  
  # Create correlation heatmap
  cor_matrix <- agreement_data %>%
    select(-replicate) %>%
    cor() %>%
    as.data.frame() %>%
    rownames_to_column("method1") %>%
    pivot_longer(cols = -method1, names_to = "method2", values_to = "correlation")
  
  heatmap_plot <- ggplot(cor_matrix, aes(x = method1, y = method2, fill = correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(correlation, 3)), color = "white", size = 5) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(title = "Method Agreement Across Replicates",
         x = "Method", y = "Method", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_file, heatmap_plot, width = 8, height = 6, dpi = 300)
  return(heatmap_plot)
}


create_pr_landscape <- function(results_list, output_file = "pr_landscape.pdf") {
  
  # Extract all PR curves across replicates
  pr_landscape_data <- map_dfr(results_list, function(repl, rep_id) {
    
    map_dfr(c("IBD", "IBS", "phylo"), function(method) {
      curve_data <- compute_roc_pr(repl$merged_data$conservative, method, truth_col = "true_link")
      
      tibble(
        replicate = rep_id,
        method = method,
        recall = curve_data$pr$curve[, 1],
        precision = curve_data$pr$curve[, 2],
        aupr = curve_data$pr$auc.integral
      )
    })
  }, seq_along(results_list))
  
  # Create PR landscape with confidence bands
  pr_summary <- pr_landscape_data %>%
    group_by(method, recall) %>%
    summarise(
      mean_precision = mean(precision, na.rm = TRUE),
      sd_precision = sd(precision, na.rm = TRUE),
      ci_lower = quantile(precision, 0.025, na.rm = TRUE),
      ci_upper = quantile(precision, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  landscape_plot <- ggplot(pr_summary, aes(x = recall, color = method, fill = method)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
    geom_line(aes(y = mean_precision), size = 1.2) +
    scale_color_manual(values = c("IBD" = "#1f77b4", "IBS" = "#ff7f0e", "phylo" = "#2ca02c")) +
    scale_fill_manual(values = c("IBD" = "#1f77b4", "IBS" = "#ff7f0e", "phylo" = "#2ca02c")) +
    labs(title = "Precision-Recall Landscape with Confidence Bands",
         x = "Recall", y = "Precision", color = "Method", fill = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(output_file, landscape_plot, width = 10, height = 8, dpi = 300)
  return(landscape_plot)
}


create_performance_radar <- function(comprehensive_results, output_file = "performance_radar.pdf") {
  
  library(fmsb)
  
  # Prepare data for radar chart
  radar_data <- comprehensive_results$performance_metrics %>%
    summarise(
      IBD = c(
        AUC = mean(IBD_AUC),
        AUPR = mean(IBD_AUPR),
        Stability = 1 - sd(IBD_AUC)/mean(IBD_AUC),
        Consistency = mean(IBD_AUC > 0.8)  # Proportion of replicates with high performance
      ),
      IBS = c(
        AUC = mean(IBS_AUC),
        AUPR = mean(IBS_AUPR), 
        Stability = 1 - sd(IBS_AUC)/mean(IBS_AUC),
        Consistency = mean(IBS_AUC > 0.8)
      ),
      Phylo = c(
        AUC = mean(Phylo_AUC),
        AUPR = mean(Phylo_AUPR),
        Stability = 1 - sd(Phylo_AUC)/mean(Phylo_AUC),
        Consistency = mean(Phylo_AUC > 0.8)
      )
    ) %>%
    as.data.frame()
  
  # Add max and min for radar chart
  radar_data <- rbind(rep(1, 3), rep(0, 3), radar_data)
  rownames(radar_data) <- c("max", "min", "AUC", "AUPR", "Stability", "Consistency")
  
  # Create radar chart
  colors_border <- c("#1f77b4", "#ff7f0e", "#2ca02c")
  colors_in <- alpha(colors_border, 0.3)
  
  pdf(output_file, width = 8, height = 8)
  radarchart(radar_data, axistype = 1,
             # Customize
             pcol = colors_border, pfcol = colors_in, plwd = 2, plty = 1,
             # Customize the axis
             cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = seq(0, 1, 0.25),
             # Variable labels
             vlcex = 0.8,
             title = "Multi-dimensional Performance Comparison")
  legend(x = 1.2, y = 1.2, legend = colnames(radar_data), 
         bty = "n", pch = 20, col = colors_in, text.col = "grey", cex = 1, pt.cex = 3)
  dev.off()
}


create_tradeoff_analysis <- function(results_list, output_file = "tradeoff_analysis.pdf") {
  
  tradeoff_data <- map_dfr(results_list, function(repl, rep_id) {
    
    # Calculate computational efficiency (simulated - replace with actual metrics)
    efficiency_metrics <- tibble(
      replicate = rep_id,
      method = c("IBD", "IBS", "Phylo"),
      computation_time = c(1.0, 0.3, 0.7),  # Relative times
      memory_usage = c(0.8, 0.2, 0.5),     # Relative memory
      accuracy = c(mean(repl$performance_summary$IBD_AUC, na.rm = TRUE),
                   mean(repl$performance_summary$IBS_AUC, na.rm = TRUE), 
                   mean(repl$performance_summary$Phylo_AUC, na.rm = TRUE))
    )
    
    return(efficiency_metrics)
  }, seq_along(results_list))
  
  tradeoff_plot <- ggplot(tradeoff_data, aes(x = computation_time, y = accuracy, color = method)) +
    geom_point(aes(size = memory_usage), alpha = 0.7) +
    geom_text(aes(label = method), nudge_y = 0.02, size = 3) +
    scale_color_manual(values = c("IBD" = "#1f77b4", "IBS" = "#ff7f0e", "Phylo" = "#2ca02c")) +
    labs(title = "Accuracy vs Computational Efficiency Trade-off",
         x = "Relative Computation Time", 
         y = "Accuracy (AUC)",
         size = "Memory Usage",
         color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(output_file, tradeoff_plot, width = 10, height = 8, dpi = 300)
  return(tradeoff_plot)
}


create_upset_plot <- function(results_list, output_file = "method_agreement_upset.pdf") {
  
  library(UpSetR)
  
  # Create binary matrix of top-performing methods per replicate
  upset_data <- map_dfr(results_list, function(repl, rep_id) {
    metrics <- repl$performance_summary
    
    # Identify which methods achieve high performance
    tibble(
      replicate = rep_id,
      IBD = as.integer(metrics$IBD_AUC > 0.8),
      IBS = as.integer(metrics$IBS_AUC > 0.8),
      Phylo = as.integer(metrics$Phylo_AUC > 0.8)
    )
  }, seq_along(results_list)) %>%
    select(-replicate)
  
  # Create UpSet plot
  pdf(output_file, width = 10, height = 6)
  print(upset(upset_data, 
              sets = c("IBD", "IBS", "Phylo"),
              mb.ratio = c(0.6, 0.4),
              order.by = "freq"))
  dev.off()
}


# Usage for Manuscript:
# Generate all innovative plots for your manuscript
dashboard <- create_performance_dashboard(comprehensive_results)
network <- create_relationship_network(merged_data, "IBD")
heatmap <- create_method_agreement_heatmap(results_list)
pr_landscape <- create_pr_landscape(results_list)
radar <- create_performance_radar(comprehensive_results)
tradeoff <- create_tradeoff_analysis(results_list)
upset_plot <- create_upset_plot(results_list)

# Create a summary figure panel
summary_panel <- gridExtra::grid.arrange(
  dashboard, network, heatmap, pr_landscape,
  ncol = 2, nrow = 2)

ggsave("manuscript_summary_panel.pdf", summary_panel,
       width = 16, height = 12, dpi = 300)