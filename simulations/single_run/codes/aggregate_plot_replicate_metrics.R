
# Put this file in your project and source it after you have `all_results`.
# Requires: data.table, ggplot2, PRROC, pROC, cowplot, viridis (optional)
library(data.table)
library(tidyverse)
library(PRROC)
library(pROC)
library(cowplot)
library(viridis)

source("simulations/single_run/aggregate_plot_functions")

all_results <- readRDS("simulations/single_run/evaluation_results.rds")

# ------------------------
# Master function to generate all figures for a given G cutoff
# ------------------------
generate_all_figures <- function(all_results, G_cut = 25, out_prefix = "simulations/single_run/figures/fig") {
        dir.create(dirname(paste0(out_prefix,"_pr_ribbon.png")), showWarnings = FALSE, recursive = TRUE)
        pr_summary <- aggregate_pr_ribbons(all_results, G_cut)
        pr_plot <- plot_pr_ribbon(pr_summary, paste0(out_prefix,"_pr_ribbon.png"))
        
        roc_summary <- aggregate_roc_ribbons(all_results, G_cut)
        roc_plot <- plot_roc_ribbon(roc_summary, paste0(out_prefix,"_roc_ribbon.png"))
        
        ranking_dt <- collect_ranking_metrics(all_results, G_cut)
        rank_plot <- plot_method_ranking(ranking_dt, paste0(out_prefix,"_method_ranking.png"))
        
        dist_plot <- plot_score_distributions(all_results, G_cut, merge_strategy = "conservative",
                                              paste0(out_prefix,"_score_distributions.png"))
        
        # Compute correlations and compare across merge strategies
        results_conservative <- compute_correlations_across_replicates(all_results, 
                                                                       merge_strategy = "conservative",
                                                                       compute_ci = TRUE,
                                                                       n_boot = 1000,
                                                                       return_plot_data = TRUE)
        
        # results_strict <- compute_correlations_across_replicates(all_results, 
        #                                                          merge_strategy = "strict",
        #                                                          compute_ci = FALSE,
        #                                                          return_plot_data = TRUE)
        # 
        # results_union <- compute_correlations_across_replicates(all_results, 
        #                                                         merge_strategy = "union",
        #                                                         compute_ci = FALSE,
        #                                                         return_plot_data = TRUE)
        # View summary
        summarize_correlation_results(results_conservative)
        
        # Create plots
        plot_scatter <- plot_correlation_results(results_conservative, plot_type = "scatter")
        ggsave(paste0(out_prefix,"_correlation_scatter.png"), plot_scatter, width = 14, height = 10)
        
        plot_violin <- plot_correlation_results(results_conservative, plot_type = "violin")
        ggsave(paste0(out_prefix,"_correlation_violin.png"), plot_scatter, width = 14, height = 10)
        
        plot_forest <- plot_correlation_results(results_conservative, plot_type = "forest")
        ggsave(paste0(out_prefix,"_correlation_forest.png"), plot_scatter, width = 14, height = 10)
        
        # cor_res <- compute_correlations_across_replicates(all_results, G_cut)
        # corr_plot <- NULL
        # if (!is.null(cor_res)) 
        #       corr_plot <- plot_correlation_heatmap(cor_res$cors_dt, paste0(out_prefix,"_correlation_heatmap.png"))
        
        # Basic usage with single strategy
        plot_correlation_heatmap(results_conservative, plot_type = "single",
                                 out_file = paste0(out_prefix,"_correlation_single.png"))

        # With confidence intervals
        plot_correlation_heatmap(results_conservative, plot_type = "range", show_ci = TRUE,
                                 out_file = paste0(out_prefix,"_correlation_range.png"))
        
        # Save to file
        plot_correlation_heatmap(results_conservative, 
                                 out_file = paste0(out_prefix,"_correlation_detailed.png"),
                                 plot_type = "detailed",
                                 font_size = 12)
        
        
        
        strategy_list <- list(
                conservative = results_conservative,
                strict = results_strict,
                union = results_union
        )
        
        plot_strategy_comparison_heatmap(strategy_list, 
                                         out_file = "strategy_comparison.png",
                                         method_subset = c("IBD", "IBS", "phylo"))
        
        composite <- assemble_composite_figure(pr_plot, roc_plot, corr_plot, rank_plot, paste0(out_prefix,"_composite.png"))
        
        # return a list of ggplot objects for further inspection
        return(list(pr = pr_plot, roc = roc_plot, ranking = rank_plot, correlation = corr_plot, distribution = dist_plot, composite = composite))
}



figs <- generate_all_figures(all_results, G_cut = 25, out_prefix = "simulations/single_run/figures/fig")
