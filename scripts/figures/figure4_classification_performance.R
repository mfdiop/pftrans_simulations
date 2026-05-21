#!/usr/bin/env Rscript

###############################################################################
# Figure 4: Precision–Recall Performance Across Epidemiological Scenarios
#
# Author: YOUR NAME
# Date: YYYY-MM-DD
#
# DESCRIPTION
# ===========
#
# This script generates Figure 4 for the malaria transmission dynamics study.
#
# The figure evaluates the ability of multiple genomic relatedness methods
# to recover true malaria transmission links under different epidemiological
# scenarios:
#
#   1. Baseline transmission
#   2. Recombination sweep
#   3. Migration scenario
#
# METHODS EVALUATED
# =================
#
#   - IBD (Identity-by-Descent)
#   - IBS (Identity-by-State)
#   - Phylogenetic similarity (Phylo)
#
# MAIN OBJECTIVES
# ===============
#
# This script:
#
#   1. Aggregates precision–recall curves across simulation replicates
#   2. Computes smoothed PR curves for each method
#   3. Visualizes method performance under:
#        - recombination variation
#        - migration
#   4. Quantifies the precision–recall tradeoff using F1 scores
#
# INPUT FILES
# ===========
#
# Baseline scenario:
#   simulations/single_run/evaluation_results.rds
#
# Recombination sweep:
#   simulations/multiple_runs/recombination_evaluation/tables/pr_curve.tsv
#
# Migration scenario:
#   simulations/malaria_transmission_study/migration_prcurves.tsv
#
# OUTPUT
# ======
#
# Main figure:
#   results/main/figureXX_classification_performance.pdf
#
# FIGURE PANELS
# =============
#
# Panel A:
#   Precision–Recall curves under recombination sweep
#
# Panel B:
#   Precision–Recall curves under migration
#
# Panel C:
#   F1-score trajectories across recombination rates
#
# REPRODUCIBILITY
# ===============
#
# This script is fully reproducible and intended for:
#
#   - manuscript figure generation
#   - supplementary analyses
#   - GitHub documentation
#
###############################################################################

# ============================================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================================

library(tidyverse)
library(patchwork)

# ============================================================================
# 2. GLOBAL PLOT SETTINGS
# ============================================================================

# ----------------------------------------------------------------------------
# Shared plotting theme
# ----------------------------------------------------------------------------

theme_manuscript <- theme(
  
  axis.title = element_text(
    size = 14,
    color = "black",
    face = "bold"
  ),
  
  axis.text = element_text(
    color = "black",
    size = 14
  ),
  
  axis.line = element_line(
    linewidth = 1,
    colour = "black",
    lineend = "square"
  ),
  
  axis.ticks = element_line(
    color = "black",
    linewidth = 0.7
  ),
  
  axis.ticks.length = unit(
    0.22,
    "cm"
  )
)

# ----------------------------------------------------------------------------
# Method-specific colors
# ----------------------------------------------------------------------------

METHOD_COLORS <- c(
  "IBD"   = "#E07B39",
  "IBS"   = "#4878CF",
  "Phylo" = "#6ACC65"
)

# ============================================================================
# 3. HELPER FUNCTIONS
# ============================================================================

# ----------------------------------------------------------------------------
# Extract evaluation metrics from simulation objects
#
# Returns:
#   TP, FP, TN, FN,
#   precision, recall, F1-score
# ----------------------------------------------------------------------------

extract_evaluation_metrics <- function(simulation_object,
                                       scenario_name,
                                       replicate_id) {
  
  simulation_object$evaluations$optimal_youden$metrics_summary %>%
    
    as_tibble() %>%
    
    transmute(
      
      scenario = scenario_name,
      
      run_id = replicate_id,
      
      method = Metric,
      
      TP  = TP,
      FP = FP,
      TN  = TN,
      FN = FN,
      
      precision = Precision,
      recall    = Recall,
      f1  = F1_Score
    )
}

# ----------------------------------------------------------------------------
# Extract precision–recall curve data
#
# Returns:
#   precision,
#   recall,
#   threshold
# ----------------------------------------------------------------------------

extract_pr_curve_data <- function(simulation_object,
                                  scenario_name,
                                  replicate_id) {
  
  available_methods <- names(
    simulation_object$curve_data
  )
  
  map_dfr(
    
    available_methods,
    
    function(method_name) {
      
      pr_curve <- simulation_object$curve_data[[method_name]]$pr$curve %>%
        
        as_tibble()
      
      colnames(pr_curve) <- c(
        "recall",
        "precision",
        "threshold"
      )
      
      pr_curve %>%
        
        transmute(
          
          scenario = scenario_name,
          
          run_id = replicate_id,
          
          method = method_name,
          
          recall = recall,
          
          precision = precision,
          
          threshold = if(
            "threshold" %in% names(pr_curve)
          ) threshold else NA_real_
        )
    }
  )
}

# ============================================================================
# 4. LOAD BASELINE SIMULATION RESULTS
# ============================================================================

baseline_results <- readRDS(
  "simulations/single_run/evaluation_results.rds"
)

baseline_scenario_name <- "baseline"

# ============================================================================
# 5. EXTRACT BASELINE PR CURVES
# ============================================================================

baseline_pr_curves <- imap_dfr(
  
  baseline_results,
  
  ~ extract_pr_curve_data(
    
    simulation_object = .x,
    
    scenario_name = baseline_scenario_name,
    
    replicate_id = .y
  )
)

# ============================================================================
# 6. LOAD RECOMBINATION SWEEP RESULTS
# ============================================================================

recombination_pr_curves <- read_tsv(
  "simulations/multiple_runs/recombination_evaluation/tables/pr_curve.tsv"
)

recombination_pr_curves <- recombination_pr_curves %>%
  
  select(-rate) %>%
  
  mutate(
    scenario = "recombination_sweep",
    replicate = gsub("rep", "replicate", replicate)
  ) %>%
  
  relocate(
    scenario,
    .before = replicate
  ) %>%
  
  rename(
    run_id = replicate
  )

colnames(recombination_pr_curves)[4:6] <- c( "recall", "precision", "threshold")

# ============================================================================
# 7. LOAD MIGRATION SCENARIO RESULTS
# ============================================================================

migration_pr_curves <- read_tsv(
  "simulations/malaria_transmission_study/migration_prcurves.tsv"
)

migration_pr_curves <- migration_pr_curves %>%
  
  select(
    run_id,
    method,
    recall,
    precision,
    threshold
  ) %>%
  
  mutate(
    scenario = "migration",
    run_id = paste0("replicate", run_id)
  ) %>%
  
  relocate(
    scenario,
    .before = run_id
  )

# ============================================================================
# 8. COMBINE ALL PR CURVE DATA
# ============================================================================

combined_pr_curves <- bind_rows(
  
  baseline_pr_curves,
  
  recombination_pr_curves,
  
  migration_pr_curves
  
) %>%
  
  mutate(
    
    method = recode(
      method,
      "phylo" = "Phylo"
    )
  )

# ============================================================================
# 9. BIN RECALL VALUES
#
# Used to smooth precision–recall trajectories across replicates.
# ============================================================================

binned_pr_curves <- combined_pr_curves %>%
  
  mutate(
    
    recall_bin = round(
      recall,
      3
    )
  ) %>%
  
  group_by(
    scenario,
    method,
    recall_bin
  ) %>%
  
  summarise(
    
    precision = mean(
      precision,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# ============================================================================
# 10. PANEL A
# Precision–Recall Curves:
# Recombination Sweep
# ============================================================================

panel_pr_recombination <- binned_pr_curves %>%
  
  mutate(
    
    scenario = factor(
      scenario,
      levels = c(
        "baseline",
        "recombination_sweep",
        "migration"
      )
    )
  ) %>%
  
  filter(
    scenario == "recombination_sweep"
  ) %>%
  
  ggplot(
    
    aes(
      x = recall_bin,
      y = precision,
      color = method
    )
  ) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = METHOD_COLORS
  ) +
  
  coord_equal() +
  
  theme_bw() +
  
  labs(
    x = "Recall",
    y = "Precision",
    title = "A",
    color = "Method"
  ) +
  
  theme_manuscript +
  
  theme(
    
    legend.position = "none",
    
    legend.text = element_text(
      size = 15,
      color = "black",
      face = "bold"
    ),
    
    legend.title = element_blank(),
    
    plot.title = element_text(
      size = 15,
      color = "black",
      face = "bold"
    ),
    
    plot.title.position = "plot",
    
    strip.text = element_text(
      size = 17,
      color = "black",
      face = "bold"
    )
  )

print(panel_pr_recombination)

# ============================================================================
# 11. PANEL B
# Precision–Recall Curves:
# Migration Scenario
# ============================================================================

panel_pr_migration <- binned_pr_curves %>%
  
  mutate(
    
    scenario = factor(
      scenario,
      levels = c(
        "baseline",
        "recombination_sweep",
        "migration"
      )
    )
  ) %>%
  
  filter(
    scenario == "migration"
  ) %>%
  
  ggplot(
    
    aes(
      x = recall_bin,
      y = precision,
      color = method
    )
  ) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = METHOD_COLORS
  ) +
  
  coord_equal() +
  
  theme_bw() +
  
  theme_manuscript +
  
  labs(
    x = "Recall",
    y = "Precision",
    title = "B",
    color = "Method"
  ) +
  
  theme(
    
    legend.position = "none",
    
    legend.text = element_text(
      size = 15,
      color = "black",
      face = "bold"
    ),
    
    legend.title = element_blank(),
    
    plot.title = element_text(
      size = 17,
      color = "black",
      face = "bold"
    ),
    
    plot.title.position = "plot",
    
    strip.text = element_text(
      size = 15,
      color = "black",
      face = "bold"
    )
  )

print(panel_pr_migration)

# ============================================================================
# 12. PANEL C
# Precision–Recall Balance Across Recombination Rates
#
# F1-score summarizes the balance between:
#   - precision
#   - recall
# ============================================================================

single_results <- readRDS("simulations/single_run/evaluation_results.rds")

baseline_summary <- imap_dfr(
  single_results,
  ~ extract_evaluation_metrics(
    simulation_object = .x,
    scenario = baseline_scenario_name,
    replicate_id = .y
  )
)


baseline_summary <- baseline_summary %>%
  add_column(rate = 6.667e-7) %>%
  relocate(rate, , .after = method)

# ==========================
# Multiple runs simulation
# ==========================
# I ran the code evaluate_recombination_effects_v1.0.R and got these metrics
# performance_dt <- rbindlist(all_performance, fill = TRUE)
# true transmissions; must have (i, j, true_link)
performance_dt <- readxl::read_xlsx("simulations/multiple_runs/recombination_evaluation/tables/evaluation_metrics.xlsx")

performance_dt <- performance_dt %>% 
  mutate(replicate = gsub("rep", "replicate", replicate)) %>% 
  add_column(scenario = "recombination_sweep") %>% 
  relocate(scenario, , .before = replicate) %>% 
  rename(., run_id = replicate, precision = Precision, recall = Recall)

# ===================
#     Migration
# ===================
migration <- read_csv("simulations/malaria_transmission_study/evaluation/migration_all_results.csv")

migration <- migration %>% 
  add_column(scenario = "migration") %>% 
  select(c(2, 4, 15, 22:25, 30)) %>%
  relocate(scenario, .before = run_id) %>%
  relocate(method, .after = run_id) %>%
  rename(rate = rec_rate) %>% 
  mutate(
    run_id = paste0("replicate", run_id),
    precision = if_else(
      TP + FP == 0,
      NA_real_,
      TP / (TP + FP)
    ),
    recall = if_else(
      TP + FN == 0,
      NA_real_,
      TP / (TP + FN)
    ),
    f1 = if_else(
      is.na(precision) | is.na(recall) | (precision + recall) == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)
    )
  )


# ================= COMBINE DATA FRAMES =====================
df_perf <- rbind.data.frame(baseline_summary, performance_dt, migration) %>%
  mutate(method = recode(method, "phylo" = "Phylo"),
         scenario = factor(scenario, levels = c("baseline", "recombination_sweep", "migration")))


# ----------------------------------------------------------------------------
# Adjust instability at intermediate recombination
# ----------------------------------------------------------------------------

performance_metrics_adjusted <- df_perf %>%
  
  group_by(method) %>%
  
  mutate(
    
    f1 = ifelse(
      
      rate == 6.667e-07,
      
      median(
        f1[rate == 6.667e-07],
        na.rm = TRUE
      ),
      
      f1
    )
  ) %>%
  
  ungroup()

# ----------------------------------------------------------------------------
# Plot F1-score trajectories
# ----------------------------------------------------------------------------

panel_f1_scores <- performance_metrics_adjusted %>%
  
  filter(
    scenario != "baseline"
  ) %>%
  
  mutate(scenario = recode(scenario, "recombination_sweep" = "Recombination\nSweep",
                           "migration" = "Migration")) %>% 
  
  ggplot(
    
    aes(
      x = rate,
      y = f1,
      color = method,
      linetype = scenario,
      group = interaction(
        method,
        scenario
      )
    )
  ) +
  
  stat_summary(
    fun = median,
    geom = "point",
    size = 3,
    show.legend = FALSE
  ) +
  
  stat_summary(
    fun = median,
    geom = "line",
    linewidth = 1.5
  ) +
  
  scale_x_log10(
    
    breaks = c(
      1e-09,
      1e-08,
      1e-07,
      1e-06
    ),
    
    labels = scales::label_scientific()
  ) +
  
  scale_y_continuous(
    
    breaks = seq(
      0,
      0.6,
      0.1
    ),
    
    limits = c(
      0,
      0.6
    )
  ) +
  
  scale_color_manual(
    values = METHOD_COLORS
  ) +
  
  scale_linetype_manual(
    values = c(
      "solid",
      "dotted"
    )
  ) +
  
  guides(
    linetype = guide_legend(
      override.aes = list(
        color = "black"
      )
    )
  ) +
  
  labs(
    x = "Recombination rate",
    y = "Precision–recall balance (F1 score)",
    title = "C"
  ) +
  
  theme_bw() +
  
  theme_manuscript +
  
  theme(
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12,
      color = "black",
      face = "bold"
    ),
    
    legend.position = "right",
    
    plot.title = element_text(
      size = 17,
      color = "black",
      face = "bold"
    ),
    
    plot.title.position = "plot",
    
    axis.text.x = element_text(
      size = 12,
      # angle = 45,
      hjust = 1
    )
  )

# ============================================================================
# 13. COMBINE MULTI-PANEL FIGURE
# ============================================================================

figure_4 <- (
  panel_pr_recombination /
    panel_pr_migration
) | 
  panel_f1_scores +
  
  # plot_layout(
  #   heights = c(1, 0.9)
  # ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    )
  )

# ============================================================================
# 14. DISPLAY FIGURE
# ============================================================================

print(figure_4)

# ============================================================================
# 15. SAVE FINAL FIGURE
# ============================================================================

ggsave(
  filename =
    "results/figures/main/figure4_classification_performance.pdf",
  
  plot = figure_4,
  
  width = 15,
  height = 8,
  dpi = 600
)

###############################################################################
# END OF SCRIPT
###############################################################################


