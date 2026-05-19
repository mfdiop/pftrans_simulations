#!/usr/bin/env Rscript

###############################################################################
# Figure 3: Effect of Recombination on True and Inferred IBD Recovery
#
# Author: YOUR NAME
# Date: YYYY-MM-DD
#
# DESCRIPTION
# ===========
# This script evaluates how recombination influences:
#
#   1. The distribution of TRUE IBD sharing
#   2. The distribution of INFERRED IBD sharing (hmmIBD)
#   3. The recovery accuracy of inferred IBD
#   4. The correlation between true and inferred pairwise IBD
#
# Biological Motivation
# ---------------------
# In malaria transmission systems, recombination reshapes genomic ancestry.
#
# At low recombination:
#   - long ancestral haplotypes persist,
#   - indirect ancestry remains difficult to distinguish,
#   - genomic similarity may not uniquely identify recent transmission.
#
# At high recombination:
#   - ancient ancestry becomes fragmented,
#   - recent transmission retains detectable shared haplotypes,
#   - direct transmission becomes more distinguishable.
#
# This script investigates whether inferred IBD accurately recovers
# true pairwise relatedness across recombination regimes.
#
# IMPORTANT ANALYTICAL NOTE
# =========================
# Correlation analyses are performed at the PAIRWISE LEVEL.
#
# Specifically:
#   correlation(true_ibd, inferred_ibd)
#
# is computed ACROSS infection pairs within each:
#
#   - replicate
#   - recombination regime
#
# This preserves biologically meaningful pair-level structure and avoids
# artifacts caused by averaging pairwise relationships before evaluating
# recovery performance.
#
# INPUT FILES
# ===========
# Expected directory structure:
#
# results/multiple_runs/inferred/
# ├── run1_rec1e09_chr1/
# │   ├── true_ibd_summary.tsv
# │   └── inferred_ibd_hmm.rds
# ├── run2_rec1e09_chr1/
# └── ...
#
# true_ibd_summary.tsv columns:
#   - Id1
#   - Id2
#   - true_ibd_prop
#
# inferred_ibd_hmm.rds columns:
#   - Id1
#   - Id2
#   - ibd
#
# OUTPUTS
# =======
# 1. Figure:
#    results/figures/main/figure3_ibd_recombination.png
#    results/figures/main/figure3_ibd_recombination.pdf
#
# 2. Pairwise correlation statistics:
#    results/tables/ibd_pairwise_correlations.tsv
#
# 3. Full merged dataset:
#    results/tables/figure3_ibd_recombination.rds
#
###############################################################################

# ============================================================================
# 1. Load Packages
# ============================================================================

library(tidyverse)
library(boot)
library(patchwork)
library(glue)

# ============================================================================
# 2. Helper Functions
# ============================================================================

# ----------------------------------------------------------------------------
# Bootstrap confidence intervals for the mean
# ----------------------------------------------------------------------------

boot_mean_ci <- function(x,
                         R = 5000,
                         conf = 0.95) {
  
  x <- x[is.finite(x)]
  
  if(length(x) < 2) {
    
    return(
      tibble(
        mean  = NA_real_,
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }
  
  boot_obj <- boot::boot(
    data = x,
    statistic = function(data, idx) {
      mean(data[idx], na.rm = TRUE)
    },
    R = R
  )
  
  ci <- boot::boot.ci(
    boot_obj,
    conf = conf,
    type = "perc"
  )
  
  tibble(
    mean  = mean(x, na.rm = TRUE),
    lower = ci$percent[4],
    upper = ci$percent[5]
  )
}

# ----------------------------------------------------------------------------
# Standardize pair ordering
#
# Example:
#   A_B == B_A
#
# Prevents duplicated pairwise relationships
# ----------------------------------------------------------------------------

standardize_pairs <- function(df,
                              id1 = "Id1",
                              id2 = "Id2") {
  
  df %>%
    rowwise() %>%
    mutate(
      pair_key = paste(sort(c( get(id1), get(id2))), collapse = "_")
    ) %>%
    ungroup() %>%
    distinct(pair_key, .keep_all = TRUE)
}

# ----------------------------------------------------------------------------
# Safe correlation function
#
# Prevents failures when insufficient observations exist
# ----------------------------------------------------------------------------

safe_cor <- function(x, y) {
  
  valid <- complete.cases(x, y)
  
  if(sum(valid) < 3) {
    return(NA_real_)
  }
  
  cor(x[valid], y[valid])
}

# ============================================================================
# 3. Define Input Directories
# ============================================================================

root <- "simulations/multiple_runs/metrics/inferred"

scenario_dirs <- list.dirs(
  root,
  recursive = TRUE,
  full.names = TRUE
)

# ============================================================================
# 4. Parse Directory Metadata
# ============================================================================

parsed_dirs <- tibble(
  path = scenario_dirs) %>%
  mutate(
    folder = basename(path),
    run_id = str_extract(folder, "(?<=run)\\d+" ),
    r_rate = str_extract(folder, "(?<=rec)[^_]+"),
    chr = str_extract(folder, "(?<=chr)\\d+")) %>%
  filter(chr == 1)

# ============================================================================
# 5. Load and Merge Pairwise IBD Data
# ============================================================================

load_true_ibd <- function(path) {
  
  file <- file.path(
    path,
    "true_ibd_summary.tsv"
  )
  
  read_tsv(
    file,
    show_col_types = FALSE
  ) %>%
    mutate(
      Id1 = paste0("tsk_", Id1),
      Id2 = paste0("tsk_", Id2)
    ) %>%
    standardize_pairs()
}

load_hmm_ibd <- function(path) {
  
  file <- file.path(
    path,
    "inferred_ibd_hmm.rds"
  )
  
  readRDS(file) %>%
    standardize_pairs()
}

# ----------------------------------------------------------------------------
# Load datasets
# ----------------------------------------------------------------------------

df_ibd <- parsed_dirs %>%
  mutate(
    true_dat = map(path, load_true_ibd),
    hmm_dat = map(path, load_hmm_ibd)
  ) %>%
  mutate(
    
    merged = map2(
      true_dat,
      hmm_dat,
      ~ inner_join(.x, .y, by = "pair_key")
      )
  ) %>%
  select(run_id, r_rate, merged
  ) %>%
  unnest(merged)

# ============================================================================
# 6. Clean and Format Data
# ============================================================================

df_ibd <- df_ibd %>%
  mutate(
    r_rate_numeric = case_when(
      r_rate == "1e09" ~ 1e-09,
      r_rate == "1e08" ~ 1e-08,
      r_rate == "1e07" ~ 1e-07,
      r_rate == "1e06" ~ 1e-06,
      TRUE ~ NA_real_
    ),
    
    r_rate_label = factor(
      r_rate_numeric,
      levels = c(
        1e-09,
        1e-08,
        1e-07,
        1e-06
      ),
      labels = c(
        "1e-09",
        "1e-08",
        "1e-07",
        "1e-06"
      )
    )
  )

# ============================================================================
# 7. Compute Pairwise Correlations Per Replicate
#
# Critical:
# Correlations are computed across PAIRS
# within each replicate and recombination regime.
# ============================================================================

cor_ibd <- df_ibd %>%
  group_by(
    run_id,
    r_rate_numeric,
    r_rate_label
  ) %>%
  summarize(
    
    pairwise_correlation = safe_cor(
      true_ibd_prop,
      hmm
    ),
    
    mean_true_ibd = mean(
      true_ibd_prop,
      na.rm = TRUE
    ),
    
    mean_inferred_ibd = mean(
      hmm,
      na.rm = TRUE
    ),
    
    n_pairs = n(),
    
    .groups = "drop"
  )

# ============================================================================
# 8. Bootstrap Confidence Intervals Across Replicates
# ============================================================================

ci_summary <- cor_ibd %>%
  group_by(
    r_rate_numeric,
    r_rate_label
  ) %>%
  summarize(
    
    true_stats = list(
      boot_mean_ci(mean_true_ibd)
    ),
    
    inferred_stats = list(
      boot_mean_ci(mean_inferred_ibd)
    ),
    
    cor_stats = list(
      boot_mean_ci(pairwise_correlation)
    ),
    
    .groups = "drop"
  ) %>%
  unnest_wider(
    true_stats,
    names_sep = "_true"
  ) %>%
  unnest_wider(
    inferred_stats,
    names_sep = "_inferred"
  ) %>%
  unnest_wider(
    cor_stats,
    names_sep = "_cor"
  )

# ============================================================================
# 9. Panel A
# Distribution of TRUE IBD
# ============================================================================

pA <- df_ibd %>%
  ggplot(
    aes(
      r_rate_label,
      true_ibd_prop,
      fill = r_rate_label)
  ) +
  
  geom_violin(
    alpha = 0.6,
    scale = "width"
  ) +
  
  geom_boxplot(
    width = 0.1,
    outlier.alpha = 0
  ) +
  
  theme_bw(base_size = 13) +
  
  labs(
    x = "Recombination rate",
    y = "True IBD proportion",
    title = "A) Distribution of true IBD"
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 15),
    plot.title.position = "plot",
    axis.title = element_text(size = 15, color = 'black', face = 'bold'),
    axis.text = element_text(size = 14, color = 'black'),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm")
  )

print(pA)

# ============================================================================
# 10. Panel B
# Distribution of INFERRED IBD
# ============================================================================

pB <- df_ibd %>%
  ggplot(
    aes(
      r_rate_label,
      hmm,
      fill = r_rate_label
    )
  ) +
  
  geom_violin(
    alpha = 0.6,
    scale = "width"
  ) +
  
  geom_boxplot(
    width = 0.1,
    outlier.alpha = 0
  ) +
  
  theme_bw() +
  
  labs(
    x = "Recombination rate",
    y = "Inferred IBD proportion",
    title = "B) Distribution of inferred IBD"
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text( face = "bold", size = 15),
    plot.title.position = "plot",
    axis.title = element_text(size = 15, color = 'black', face = 'bold'),
    axis.text = element_text(size = 14, color = 'black'),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm")
  )

# ============================================================================
# 11. Panel C
# Mean TRUE vs INFERRED IBD
# ============================================================================

pC <- ci_summary %>%
  ggplot(
    aes(r_rate_numeric)
  ) +
  
  geom_ribbon(
    aes(
      ymin = true_stats_truelower,
      ymax = true_stats_trueupper
    ),
    fill = "steelblue",
    alpha = 0.2
  ) +
  
  geom_line(
    aes(y = true_stats_truemean),
    color = "steelblue",
    linewidth = 1.2
  ) +
  
  geom_point(
    aes(y = true_stats_truemean),
    color = "steelblue",
    size = 3
  ) +
  
  geom_ribbon(
    aes(
      ymin = inferred_stats_inferredlower,
      ymax = inferred_stats_inferredupper
    ),
    fill = "firebrick",
    alpha = 0.2
  ) +
  
  geom_line(
    aes(y = inferred_stats_inferredmean),
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  geom_point(
    aes(y = inferred_stats_inferredmean),
    color = "firebrick",
    size = 3
  ) +
  
  scale_x_log10(
    breaks = c(
      1e-09,
      1e-08,
      1e-07,
      1e-06
    ),
    labels = c(
      "1e-09",
      "1e-08",
      "1e-07",
      "1e-06"
    )
  ) +
  
  theme_bw() +
  
  labs(
    x = "Recombination rate",
    y = "Mean IBD proportion",
    title = "C) Mean true and inferred IBD"
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.title.position = "plot",
    axis.title = element_text(size = 15, color = 'black', face = 'bold'),
    axis.text = element_text(size = 14, color = 'black'),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm")
  )

# ============================================================================
# 12. Panel D
# Pairwise Recovery Correlation
# ============================================================================

pD <- ci_summary %>%
  ggplot(
    aes(
      r_rate_numeric,
      cor_stats_cormean)) +
  
  geom_ribbon(
    aes(
      ymin = cor_stats_corlower,
      ymax = cor_stats_corupper
    ),
    fill = "grey70",
    alpha = 0.3) +
  
  geom_line(linewidth = 1.2) +
  
  geom_point(size = 3) +
  
  scale_x_log10(
    breaks = c(1e-09, 1e-08, 1e-07, 1e-06),
    labels = c( "1e-09", "1e-08", "1e-07", "1e-06")) +
  
  ylim(0, 1) +
  
  theme_bw(base_size = 13) +
  
  labs(
    x = "Recombination rate",
    y = "Pairwise correlation",
    title = "D) Recovery of pairwise IBD relationships"
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.title.position = "plot",
    axis.title = element_text(size = 15, color = 'black', face = 'bold'),
    axis.text = element_text(size = 14, color = 'black'),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm")
  )

# ============================================================================
# 13. Combine Figure Panels
# ============================================================================

final_figure <- (pA | pB) / (pC | pD)

# ============================================================================
# 14. Save Figure
# ============================================================================

ggsave(
  filename = "results/figures/main/figure3_ibd_recombination.pdf",
  plot = final_figure,
  width = 16, height = 10,
  dpi = 600
)

# ============================================================================
# 15. Export Summary Statistics
# ============================================================================

write_tsv(
  cor_ibd,
  "results/tables/ibd_pairwise_correlations.tsv"
)

# ============================================================================
# 16. Save Full Workspace Object
# ============================================================================

saveRDS(
  list(
    full_data = df_ibd,
    correlation_summary = cor_ibd,
    bootstrap_summary = ci_summary
  ),
  "results/rds/figure3_ibd_recombination.rds"
)

###############################################################################
# End of Script
###############################################################################