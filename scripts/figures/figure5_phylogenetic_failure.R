#!/usr/bin/env Rscript

###############################################################################
# Figure 5: Failure of Phylogenetic Metrics to Resolve Direct Transmission
#
# Author: Mouhamadou Fadel DIOP
# Date: YYYY-MM-DD
#
# DESCRIPTION
# ===========
#
# This script evaluates whether phylogenetic distances can accurately recover
# direct malaria transmission events across varying recombination regimes.
#
# Specifically, this analysis investigates whether:
#
#   - cophenetic distances,
#   - phylogenetic similarity,
#   - and tree-derived clustering
#
# can distinguish:
#
#   - direct transmission,
#   - recent ancestry,
#   - and indirect shared ancestry.
#
# BIOLOGICAL RATIONALE
# ====================
#
# Phylogenetic approaches assume that recently transmitted parasites should
# cluster closely on phylogenetic trees.
#
# However, malaria recombination creates mosaic genomes that violate simple
# tree-like inheritance assumptions.
#
# Consequently:
#
#   - closely related parasites may not represent direct transmission,
#   - direct transmission pairs may appear phylogenetically distant,
#   - recombination progressively decouples ancestry from transmission.
#
# This script quantifies these effects across recombination rates.
#
# FIGURE PANELS
# =============
#
# Panel A
# -------
# Distribution of true TMRCA values across transmission classes.
#
# Panel B
# -------
# Distribution of phylogenetic (cophenetic) distances for:
#   - direct transmission,
#   - indirect transmission.
#
# Panel C
# -------
# Relationship between:
#   - true IBD sharing,
#   - phylogenetic distance.
#
# Panel D
# -------
# PRAUC of phylogenetic distance for predicting direct transmission
# across recombination rates.
#
# Supplementary Figure
# --------------------
# Misclassification heatmaps per replicate.
#
# INPUT FILES
# ===========
#
# Required directory structure:
#
# simulations/multiple_runs/
# ├── metrics/
# │   ├── phylo_results/
# │   └── inferred/
#
# OUTPUT FILES
# ============
#
# results/figures/main/figure5_phylogenetic_failure.pdf
# results/figures/supplementary/supplementary_figure5_heatmaps.pdf
#
# REPRODUCIBILITY
# ===============
#
# This script:
#
#   - automatically detects recombination scenarios,
#   - iterates across replicates,
#   - standardizes pairwise relationships,
#   - and produces publication-quality figures.
#
###############################################################################

# ============================================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================================

library(tidyverse)
library(ape)
library(phangorn)
library(pROC)
library(PRROC)
library(patchwork)
library(glue)

# ============================================================================
# 2. GLOBAL PLOT SETTINGS
# ============================================================================

theme_manuscript <- theme_bw() +
  
  theme(
    
    axis.title = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    
    axis.text = element_text(
      size = 13,
      color = "black"
    ),
    
    legend.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.text = element_text(
      size = 12, color = "black"
    ),
    
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    
    axis.line = element_line(
      linewidth = 0.9,
      color = "black",
      lineend = "square",
      linejoin = "round"
      
    ),
    
    axis.ticks = element_line(
      linewidth = 0.7,
      color = "black"
    ),
    
    axis.ticks.length = unit(
      0.22,
      "cm"
    )
  )

# ----------------------------------------------------------------------------
# Transmission class colors
# ----------------------------------------------------------------------------

TRANSMISSION_COLORS <- c(
  "Indirect" = "grey40",
  "Direct"   = "blue3"
)

# ============================================================================
# 3. IDENTIFY RECOMBINATION SCENARIOS
# ============================================================================

recombination_directories <- list.dirs(
  "simulations/multiple_runs/metrics/inferred",
  recursive = TRUE,
  full.names = TRUE
)

recombination_metadata <- tibble(
  directory = recombination_directories
) %>%
  
  mutate(
    
    folder_name = basename(directory),
    
    replicate_id = str_extract(
      folder_name,
      "(?<=run)\\d+"
    ),
    
    recombination_rate = str_extract(
      folder_name,
      "(?<=rec)[^_]+"
    ),
    
    chromosome = str_extract(
      folder_name,
      "(?<=chr)\\d+"
    )
  ) %>%
  
  filter(
    !is.na(recombination_rate),
    chromosome == 1
  )

# ============================================================================
# 4. HELPER FUNCTION:
# STANDARDIZE PAIR ORDER
# ============================================================================

standardize_pair_order <- function(dataframe,
                                   sample1 = "sample1",
                                   sample2 = "sample2") {
  
  dataframe %>%
    
    rowwise() %>%
    
    mutate(
      
      pair_id = paste(
        sort(c(
          get(sample1),
          get(sample2)
        )),
        collapse = "_"
      )
    ) %>%
    
    ungroup()
}

# ============================================================================
# 5. LOAD PHYLOGENETIC DATA
# ============================================================================

load_phylogenetic_replicate <- function(recombination_rate,
                                        replicate_id) {
  
  # --------------------------------------------------------------------------
  # File paths
  # --------------------------------------------------------------------------
  
  tree_file <- glue(
    "simulations/multiple_runs/metrics/phylo_results/",
    "rep{replicate_id}/",
    "run{replicate_id}_rec{recombination_rate}_chr1.treefile"
  )
  
  truth_file <- glue(
    "simulations/multiple_runs/metrics/inferred/",
    "rep{replicate_id}/",
    "run{replicate_id}_rec{recombination_rate}_chr1/",
    "true_ibd_summary.tsv"
  )
  
  inferred_file <- glue(
    "simulations/multiple_runs/metrics/inferred/",
    "rep{replicate_id}/",
    "run{replicate_id}_rec{recombination_rate}_chr1/",
    "inferred_ibd_hmm.tsv"
  )
  
  # --------------------------------------------------------------------------
  # Load phylogenetic tree
  # --------------------------------------------------------------------------
  
  phylogenetic_tree <- read.tree(tree_file)
  
  # --------------------------------------------------------------------------
  # Load true IBD relationships
  # --------------------------------------------------------------------------
  
  true_ibd <- read_tsv(
    truth_file,
    show_col_types = FALSE
  )
  
  # --------------------------------------------------------------------------
  # Load inferred IBD
  # --------------------------------------------------------------------------
  
  inferred_ibd <- read_tsv(
    inferred_file,
    show_col_types = FALSE
  ) %>%
    
    rename(
      sample_1 = Id1,
      sample_2 = Id2
    )
  
  # --------------------------------------------------------------------------
  # Compute pairwise cophenetic distances
  # --------------------------------------------------------------------------
  
  cophenetic_matrix <- cophenetic(
    phylogenetic_tree
  ) %>%
    
    as.data.frame() %>%
    
    rownames_to_column(
      "sample_1"
    )
  
  cophenetic_long <- cophenetic_matrix %>%
    
    pivot_longer(
      
      cols = -sample_1,
      
      names_to = "sample_2",
      
      values_to = "cophenetic_distance"
    )
  
  # --------------------------------------------------------------------------
  # Merge all pairwise relationships
  # --------------------------------------------------------------------------
  
  merged_data <- true_ibd %>%
    
    select(
      -Id1,
      -Id2
    ) %>%
    
    separate(
      pair,
      into = c(
        "sample_1",
        "sample_2"
      ),
      sep = "_"
    ) %>%
    
    mutate(
      
      sample_1 = paste0(
        "tsk_",
        sample_1
      ),
      
      sample_2 = paste0(
        "tsk_",
        sample_2
      )
    ) %>%
    
    standardize_pair_order() %>%
    
    left_join(
      
      inferred_ibd,
      
      by = c(
        "sample_1",
        "sample_2"
      )
    ) %>%
    
    left_join(
      
      cophenetic_long,
      
      by = c(
        "sample_1",
        "sample_2"
      )
    ) %>%
    
    mutate(
      
      recombination_rate = recombination_rate,
      
      recombination_rate = factor(
        recombination_rate,
        levels = c("1e09", "1e08", "1e07", "1e06"),
        labels = c("1e-09", "1e-08", "1e-07", "1e-06")
        ),
      
      replicate_id = replicate_id,
      
      is_direct = min_tmrca == 1,
      
      transmission_class = if_else(
        is_direct,
        "Direct",
        "Indirect"
      )
    )
  
  return(merged_data)
}

# ============================================================================
# 6. LOAD ALL RECOMBINATION SCENARIOS
# ============================================================================

phylogenetic_results <- pmap_dfr(
  
  list(
    recombination_metadata$recombination_rate,
    recombination_metadata$replicate_id
  ),
  
  load_phylogenetic_replicate
)

# ============================================================================
# 7. PANEL A
# TRUE TMRCA DISTRIBUTIONS
# ============================================================================

panel_tmrca_distribution <- phylogenetic_results %>%
  
  ggplot(
    
    aes(
      x = min_tmrca,
      fill = transmission_class
    )
  ) +
  
  geom_density(
    alpha = 0.55
  ) + scale_x_log10() +
  
  facet_wrap(
    ~ recombination_rate,
    scales = "free"
  ) +
  
  scale_fill_manual(
    values = TRANSMISSION_COLORS
  ) +
  
  labs(
    x = "True TMRCA (generations)",
    y = "Density",
    fill = "Transmission"
  ) +
  
  theme_manuscript

# ============================================================================
# 8. PANEL B
# COPHENETIC DISTANCE DISTRIBUTIONS
# ============================================================================

panel_phylogenetic_distance <- phylogenetic_results %>%
  
  ggplot(
    
    aes(
      x = cophenetic_distance,
      fill = transmission_class
    )
  ) +
  
  geom_density(
    alpha = 0.6
  ) +
  
  facet_wrap(
    ~ recombination_rate,
    scales = "free"
  ) +
  
  scale_fill_manual(
    values = TRANSMISSION_COLORS
  ) +
  
  labs(
    x = "Cophenetic distance",
    y = "Density",
    fill = "Transmission"
  ) +
  
  theme_manuscript

# ============================================================================
# 9. PANEL C
# IBD VS PHYLOGENETIC DISTANCE
# ============================================================================

panel_ibd_phylogenetic <- phylogenetic_results %>%
  
  ggplot(
    
    aes(
      x = true_ibd_prop,
      y = cophenetic_distance,
      color = transmission_class
    )
  ) +
  
  geom_point(
    alpha = 0.45,
    size = 1.2
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  
  facet_wrap(
    ~ recombination_rate
  ) +
  
  scale_color_manual(
    values = TRANSMISSION_COLORS
  ) +
  
  labs(
    x = "True IBD proportion",
    y = "Cophenetic distance",
    color = "Transmission"
  ) +
  
  theme_manuscript

# ============================================================================
# 10. PANEL D
# PRAUC ACROSS RECOMBINATION RATES
# ============================================================================

safe_prauc <- function(labels,
                       scores) {
  
  valid_indices <- !is.na(scores)
  
  labels <- labels[valid_indices]
  scores <- scores[valid_indices]
  
  if(length(unique(labels)) < 2) {
    return(NA_real_)
  }
  
  if(sum(labels == 1) < 2) {
    return(NA_real_)
  }
  
  pr_result <- PRROC::pr.curve(
    
    scores.class0 = scores[labels == 1],
    
    scores.class1 = scores[labels == 0],
    
    curve = FALSE
  )
  
  return(
    pr_result$auc.integral
  )
}

prauc_results <- phylogenetic_results %>%
  
  group_by(
    recombination_rate,
    replicate_id
  ) %>%
  
  summarise(
    
    pr_auc = safe_prauc(
      labels = is_direct,
      scores = -cophenetic_distance
    ),
    
    .groups = "drop"
  )

panel_prauc <- prauc_results %>%
  
  ggplot(
    
    aes(
      x = recombination_rate,
      y = pr_auc
    )
  ) +
  
  geom_boxplot(
    fill = "grey85",
    outlier.alpha = 0.5
  ) +
  
  geom_jitter(
    width = 0.15, 
    size = 3,
    alpha = 0.6,
    color = "blue3"
  ) +
  
  ylim(0, 1) +
  
  labs(
    x = "Recombination rate",
    y = "PRAUC"
  ) +
  
  theme_manuscript +
  
  theme(
    axis.text.x = element_text(
      size = 12,
      color = "black"
      # angle = 45,
      # hjust = 1
    )
  )

# ============================================================================
# 11. COMBINE MAIN FIGURE
# ============================================================================

figure_5 <- (
  
  panel_tmrca_distribution |
    
    panel_phylogenetic_distance
  
) / (
  
  panel_ibd_phylogenetic |
    
    panel_prauc
) +
  
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
# 12. SAVE MAIN FIGURE
# ============================================================================

# ggsave(
#   
#   filename =
#     "results/figures/main/figure5_phylogenetic_failure.pdf",
#   
#   plot = figure_5,
#   
#   width = 16,
#   height = 12,
#   dpi = 600
# )

# ============================================================================

ggsave(filename = "results/figures/main/tmrca_distribution.pdf",
  plot = panel_tmrca_distribution, width = 16, height = 12, dpi = 600)

# ============================================================================

ggsave(filename = "results/figures/main/figure5_phylogenetic_failure.pdf",
       plot = panel_phylogenetic_distance, width = 16, height = 12, dpi = 600)

# ============================================================================

ggsave(filename = "results/figures/main/figure5_phylogenetic_prauc.pdf",
       plot = panel_prauc, width = 16, height = 12, dpi = 600)

# ============================================================================
# 13. SUPPLEMENTARY:
# MISCLASSIFICATION HEATMAPS
# ============================================================================

supplementary_heatmap <- phylogenetic_results %>%
  
  mutate(
    
    predicted_direct =
      cophenetic_distance <
      median(
        cophenetic_distance,
        na.rm = TRUE
      )
  ) %>%
  
  count(
    
    recombination_rate,
    replicate_id,
    is_direct,
    predicted_direct
  ) %>%
  
  ggplot(
    
    aes(
      x = factor(is_direct),
      y = factor(predicted_direct),
      fill = n
    )
  ) +
  
  geom_tile() +
  
  facet_grid(
    replicate_id ~ recombination_rate
  ) +
  
  scale_fill_viridis_c() +
  
  labs(
    x = "True direct transmission",
    y = "Predicted direct transmission",
    fill = "Pairs"
  ) +
  
  theme_manuscript

# ============================================================================
# 14. SAVE SUPPLEMENTARY FIGURE
# ============================================================================

ggsave(
  
  filename =
    "results/figures/supplementary/supplementary_figure5_heatmaps.pdf",
  
  plot = supplementary_heatmap,
  
  width = 14,
  height = 10,
  dpi = 600
)

###############################################################################
# END OF SCRIPT
###############################################################################