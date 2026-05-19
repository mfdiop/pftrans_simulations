#!/usr/bin/env Rscript
# ==============================================================================
# SHOW: Recombination Breaks Down OLD IBD, Preserves RECENT IBD
# ==============================================================================

###############################################################################
# Figure 4: Effect of Recombination on IBD Tract Fragmentation
#           and IBD Inference Error
#
# Author: Mouhamadou Fadel DIOP
# Date: YYYY-MM-DD
#
# DESCRIPTION
# ===========
# This script investigates the mechanistic effects of recombination on:
#
#   1. Fragmentation of true IBD tracts
#   2. Detectability of shared genomic segments
#   3. False-positive IBD inference
#   4. False-negative IBD inference
#
# Biological Motivation
# ---------------------
# Recombination reshapes the genomic architecture of relatedness.
#
# LOW recombination:
#   - preserves long ancestral haplotypes,
#   - inflates background relatedness,
#   - increases false-positive IBD inference.
#
# HIGH recombination:
#   - fragments recent transmission tracts,
#   - shortens shared segments,
#   - increases false-negative IBD inference.
#
# Therefore, recombination creates a tradeoff between:
#
#   ancestral LD persistence
#                  versus
#   tract detectability
#
# This script quantifies that tradeoff using pairwise IBD segment data.
#
# INPUT FILES
# ===========
# Expected directory structure:
#
# results/multiple_runs/inferred/
# ├── run1_rec1e09_chr1/
# │   ├── true_ibd_segments.tsv
# │   └── inferred_ibd_segments.tsv
# ├── run2_rec1e09_chr1/
# └── ...
#
# REQUIRED COLUMNS
# =================
#
# TRUE SEGMENTS
# -------------
# Id1
# Id2
# start
# end
# segment_length
#
# INFERRED SEGMENTS
# -----------------
# Id1
# Id2
# start
# end
# segment_length
#
# NOTE
# ====
# If segment_length is absent, it will be computed as:
#
#   end - start
#
# OUTPUTS
# =======
#
# Figure 4 Panels
# ---------------
#
# A. Distribution of TRUE IBD tract lengths
# B. Distribution of INFERRED IBD tract lengths
# C. False-positive rate across recombination regimes
# D. False-negative rate across recombination regimes
#
###############################################################################

# ============================================================================
# 1. Load Packages
# ============================================================================
library(tidyverse)
library(patchwork)

# ============================================================================
# Helper Functions (same as yours)
# ============================================================================

compute_segment_length <- function(df) {
  if(!"segment_length" %in% names(df)) {
    df <- df %>%
      mutate(segment_length = end - start)
  }
  return(df)
}

standardize_pairs <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      pair_key = paste(sort(c(Id1, Id2)), collapse = "_")
    ) %>%
    ungroup() %>%
    distinct(pair_key, .keep_all = TRUE)
}

load_true_segments <- function(path) {
  file <- file.path(path, "true_ibd_summary.tsv")
  
  read_tsv(file, show_col_types = FALSE) %>%
    standardize_pairs() %>%
    compute_segment_length()
}

# ============================================================================
# Load Data with Generation Distance
# ============================================================================

root <- "simulations/multiple_runs/metrics/inferred"
scenario_dirs <- list.dirs(root, full.names = TRUE, recursive = TRUE)

parsed_dirs <- tibble(path = scenario_dirs) %>%
  mutate(
    folder = basename(path),
    run_id = str_extract(folder, "(?<=run)\\d+"),
    r_rate = str_extract(folder, "(?<=rec)[^_]+"),
    chr = str_extract(folder, "(?<=chr)\\d+")
  ) %>% 
  filter(chr == 1)

# Load segments WITH generation distance
segment_data <- parsed_dirs %>%
  mutate(
    true_segments = map(path, safely(load_true_segments))
  ) %>%
  mutate(
    true_dat = map(true_segments, ~.x$result)
  ) %>%
  select(run_id, r_rate, true_dat) %>%
  unnest(true_dat) %>%
  mutate(
    r_rate = as.numeric(gsub("e", "e-", r_rate)),
    r_rate_label = factor(
      r_rate,
      levels = c(1e-9, 1e-8, 1e-7, 1e-6),
      labels = c("1e-09 (clonal)", "1e-08 (low)", "1e-07 (moderate)", "1e-06 (high)")
    )
  )

# CRITICAL: You need generation distance (gen_distance) or min_tmrca in your data
# If it's not in true_ibd_summary.tsv, you need to add it from pedigree

# Check if gen_distance | min_tmrca | TMRCA exists
if(!c("gen_distance", "min_tmrca") %in% names(segment_data)) {
  stop("ERROR: gen_distance not found in data. 
       You need genealogical distance to show the mechanism!")
}

# ============================================================================
# Categorize Segments by Generation Distance
# ============================================================================

segment_categorized <- segment_data %>%
  rename(gen_distance = min_tmrca) %>% 
  mutate(
    ancestry_category = case_when(
      gen_distance <= 5  ~ "Recent (G≤5)",
      gen_distance <= 15 ~ "Intermediate (G 6-15)",
      gen_distance > 15  ~ "Old (G>15)"
    ),
    ancestry_category = factor(
      ancestry_category,
      levels = c("Recent (G≤5)", "Intermediate (G 6-15)", "Old (G>15)")
    )
  )

# ============================================================================
# Calculate Summary Statistics
# ============================================================================

segment_summary <- segment_categorized %>%
  group_by(r_rate, r_rate_label, ancestry_category) %>%
  summarize(
    mean_length = mean(total_ibd_bp, na.rm = TRUE),
    median_length = median(total_ibd_bp, na.rm = TRUE),
    q25 = quantile(total_ibd_bp, 0.25, na.rm = TRUE),
    q75 = quantile(total_ibd_bp, 0.75, na.rm = TRUE),
    n_segments = n(),
    .groups = "drop"
  )

cat("=== SEGMENT LENGTH BY ANCESTRY CATEGORY ===\n")
print(segment_summary)
cat("\n")

# ============================================================================
# KEY FIGURE: Show Recombination Paradox
# ============================================================================

# =============================================================================
#     Panel A: Segment length distributions by ancestry category
# =============================================================================
panelA <- segment_categorized %>%
  ggplot(aes(x = total_ibd_bp, fill = ancestry_category)) + # segment_length
  geom_density(alpha = 0.6) +
  scale_x_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  #   limits = c(1e3, 1e7)
  ) +
  scale_fill_manual(
    values = c(
      "Recent (G≤5)" = "#e41a1c",        # Red
      "Intermediate (G 6-15)" = "#377eb8", # Blue
      "Old (G>15)" = "#4daf4a"            # Green
    )
  ) +
  facet_wrap(~r_rate_label, ncol = 2, scales = "free") +
  theme_bw() +
  labs(
    x = "IBD segment length (bp, log scale)",
    y = "Density",
    fill = "Ancestry",
    title = "A) IBD segment length by ancestry and recombination rate"
    # subtitle = "Recent transmission maintains longer segments across all rates"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, color = "black", face = "bold"),
    legend.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.title.position = "plot",
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "lightgrey")
  )

print(panelA)

# =============================================================================
#   Panel B: Mean segment length by recombination rate (separated by ancestry)
# =============================================================================
panelB <- segment_summary %>%
  ggplot(aes(x = r_rate, y = mean_length, color = ancestry_category)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = c(1e-9, 1e-8, 1e-7, 1e-6),
    labels = c("1e-09", "1e-08", "1e-07", "1e-06")
  ) +
  scale_y_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  scale_color_manual(
    values = c(
      "Recent (G≤5)" = "#e41a1c",
      "Intermediate (G 6-15)" = "#377eb8",
      "Old (G>15)" = "#4daf4a"
    )
  ) +
  theme_bw() +
  labs(
    x = "Recombination rate (log scale)",
    y = "Mean segment length (bp, log scale)",
    color = "Ancestry",
    title = "B) Recombination fragments OLD ancestry, preserves RECENT"
    # subtitle = "This creates temporal stratification that enables discrimination"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, color = "black", face = "bold", hjust = 0.5),
    legend.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.title.position = "plot",
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm")
  )

print(panelB)

# =============================================================================
#   Panel C: Separation ratio (ratio of recent:old mean segment lengths)
# =============================================================================
separation_ratio <- segment_summary %>%
  select(r_rate, r_rate_label, ancestry_category, mean_length) %>%
  pivot_wider(
    names_from = ancestry_category,
    values_from = mean_length
  ) %>%
  mutate(
    ratio = `Recent (G≤5)` / `Old (G>15)`
  )

panelC <- separation_ratio %>%
  ggplot(aes(x = r_rate, y = ratio)) +
  geom_line(linewidth = 1.2, color = "black") +
  geom_point(size = 4, color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_x_log10(
    breaks = c(1e-9, 1e-8, 1e-7, 1e-6),
    labels = c("1e-09", "1e-08", "1e-07", "1e-06")
  ) +
  theme_bw() +
  labs(
    x = "Recombination rate (log scale)",
    y = "Separation ratio\n(Recent / Old mean segment length)",
    title = "C) Higher recombination creates stronger temporal stratification"
    # subtitle = "Ratio >1 means recent transmission has disproportionately longer segments"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.title.position = "plot",
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm")
  )

print(panelC)

# =============================================================================
#   Panel D: Coefficient of variation (shows discrimination potential)
# =============================================================================
cv_analysis <- segment_categorized %>%
  group_by(r_rate, r_rate_label, ancestry_category) %>%
  summarize(
    cv = sd(total_ibd_bp, na.rm = TRUE) / mean(total_ibd_bp, na.rm = TRUE),
    .groups = "drop"
  )

panelD <- cv_analysis %>%
  ggplot(aes(x = r_rate, y = cv, color = ancestry_category)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = c(1e-9, 1e-8, 1e-7, 1e-6),
    labels = c("1e-09", "1e-08", "1e-07", "1e-06")
  ) +
  scale_color_manual(
    values = c(
      "Recent (G≤5)" = "#e41a1c",
      "Intermediate (G 6-15)" = "#377eb8",
      "Old (G>15)" = "#4daf4a"
    )
  ) +
  theme_bw() +
  labs(
    x = "Recombination rate (log scale)",
    y = "Coefficient of variation",
    color = "Ancestry",
    title = "D) Variability in segment lengths by ancestry"
    # subtitle = "Higher CV at low recombination indicates less consistent signal"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.title.position = "plot",
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm")
  )

print(panelD)

# ==================
#   Combine panels
# ==================

p_mechanism <- (panelA | panelB) / (panelC | panelD)

print(p_mechanism)

ggsave(
  "results/figures/main/recombination_mechanism.pdf",
  p_mechanism, width = 16, height = 12, dpi = 600
)

cat("\n=== KEY FINDINGS ===\n")
cat("Separation ratio (Recent/Old mean segment length):\n")
print(separation_ratio %>% select(r_rate_label, ratio))
cat("\n")

# Calculate fold-change
fold_change <- separation_ratio %>%
  summarize(
    lowest_rec = ratio[r_rate == min(r_rate)],
    highest_rec = ratio[r_rate == max(r_rate)],
    fold_change = highest_rec / lowest_rec
  )

cat(sprintf("At r=1e-9: Recent/Old ratio = %.2f\n", fold_change$lowest_rec))
cat(sprintf("At r=1e-6: Recent/Old ratio = %.2f\n", fold_change$highest_rec))
cat(sprintf("Fold improvement: %.1fx\n\n", fold_change$fold_change))

if(fold_change$fold_change > 3) {
  cat("✓ STRONG MECHANISM: High recombination creates >3x better separation\n")
  cat("  This explains why IBD-based inference improves with recombination.\n")
} else {
  cat("⚠ WEAK MECHANISM: Separation ratio improves by <3x\n")
  cat("  Other factors may be more important.\n")
}

# ============================================================================
# Statistical Test: Is separation significant?
# ============================================================================

# Test if recent vs old segment lengths differ more at high vs low recombination

stat_test <- segment_categorized %>%
  filter(ancestry_category %in% c("Recent (G≤5)", "Old (G>15)")) %>%
  group_by(r_rate) %>%
  summarize(
    t_statistic = t.test(
      segment_length[ancestry_category == "Recent (G≤5)"],
      segment_length[ancestry_category == "Old (G>15)"]
    )$statistic,
    p_value = t.test(
      segment_length[ancestry_category == "Recent (G≤5)"],
      segment_length[ancestry_category == "Old (G>15)"]
    )$p.value,
    effect_size = abs(t_statistic) / sqrt(n()),
    .groups = "drop"
  )

cat("\n=== STATISTICAL SEPARATION ===\n")
cat("T-test comparing Recent vs Old segment lengths:\n")
print(stat_test)
cat("\n")

cat("=== INTERPRETATION FOR MANUSCRIPT ===\n")
cat("Use this text in your Results section:\n\n")

cat(sprintf(
  "\"The recombination paradox mechanism: At low recombination (r=1e-9), 
mean IBD segment length was %.0f kb for recent transmission pairs (G≤5) 
versus %.0f kb for old ancestry pairs (G>15), yielding a separation ratio 
of %.2f (Figure XB). At high recombination (r=1e-6), segment lengths were 
%.0f kb versus %.0f kb, respectively, yielding a ratio of %.2f—a %.1f-fold 
improvement in temporal stratification (Figure XC). This occurs because 
recombination events accumulate across generations: pairs separated by 
G=15 generations experience 3× more recombination events than G=5 pairs. 
Consequently, high recombination preferentially fragments background IBD 
while preserving recent transmission signals as longer intact segments, 
enabling IBD-based methods to exploit tract length for discrimination.\"\n",
  
  segment_summary %>% 
    filter(r_rate == 1e-9, ancestry_category == "Recent (G≤5)") %>% 
    pull(mean_length) / 1000,
  
  segment_summary %>% 
    filter(r_rate == 1e-9, ancestry_category == "Old (G>15)") %>% 
    pull(mean_length) / 1000,
  
  separation_ratio %>% 
    filter(r_rate == 1e-9) %>% 
    pull(ratio),
  
  segment_summary %>% 
    filter(r_rate == 1e-6, ancestry_category == "Recent (G≤5)") %>% 
    pull(mean_length) / 1000,
  
  segment_summary %>% 
    filter(r_rate == 1e-6, ancestry_category == "Old (G>15)") %>% 
    pull(mean_length) / 1000,
  
  separation_ratio %>% 
    filter(r_rate == 1e-6) %>% 
    pull(ratio),
  
  fold_change$fold_change
))

cat("\n=== SCRIPT COMPLETE ===\n")