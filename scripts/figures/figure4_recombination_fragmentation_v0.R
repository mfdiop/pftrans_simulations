#!/usr/bin/env Rscript

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
library(glue)

# ============================================================================
# 2. Define Input Directories
# ============================================================================

root <- "simulations/multiple_runs/metrics/inferred"

scenario_dirs <- list.dirs(
  root,
  recursive = TRUE,
  full.names = TRUE
)

# ============================================================================
# 3. Parse Metadata
# ============================================================================

parsed_dirs <- tibble(
  path = scenario_dirs) %>%
  mutate(
    
    folder = basename(path),
    run_id = str_extract(
      folder,
      "(?<=run)\\d+"),
    
    r_rate = str_extract(
      folder,
      "(?<=rec)[^_]+"
    ),
    
    chr = str_extract(
      folder,
      "(?<=chr)\\d+"
    )
    
  ) %>%
  filter(chr == 1)

# ============================================================================
# 4. Helper Functions
# ============================================================================

# ----------------------------------------------------------------------------
# Standardize pair ordering
# ----------------------------------------------------------------------------

standardize_pairs <- function(df,
                              id1 = "Id1",
                              id2 = "Id2") {
  
  df %>%
    rowwise() %>%
    mutate(
      pair_key = paste(sort(c( get(id1), get(id2) )), collapse = "_")
    ) %>%
    ungroup()
}

# ----------------------------------------------------------------------------
# Compute segment length if absent
# ----------------------------------------------------------------------------

compute_segment_length <- function(df) {
  
  if(!"segment_length" %in% names(df)) {
    
    df <- df %>%
      mutate(
        segment_length = end - start
      )
  }
  
  return(df)
}

# ============================================================================
# 5. Load TRUE IBD Segments
# ============================================================================

load_true_segments <- function(path) {
  
  file <- file.path(
    path,
    "true_ibd_summary.tsv"
  )
  
  read_tsv(
    file,
    show_col_types = FALSE
  ) %>%
    standardize_pairs() #%>%
    # compute_segment_length()
}

# ============================================================================
# 6. Load INFERRED IBD Segments
# ============================================================================

load_inferred_segments <- function(path) {
  
  file <- file.path(
    path,
    "inferred_ibd_segments.tsv"
  )
  
  read_tsv(
    file,
    show_col_types = FALSE
  ) %>%
    standardize_pairs() %>%
    compute_segment_length()
}

# ============================================================================
# 7. Load All Segment Data
# ============================================================================

segment_data <- parsed_dirs %>%
  mutate(
    
    true_segments = map(
      path,
      safely(load_true_segments)
    )
    
    # inferred_segments = map(
    #   path,
    #   safely(load_inferred_segments)
    # )
    
  )

# ============================================================================
# 8. Extract Valid Segment Tables
# ============================================================================

extract_safe <- function(x) {
  
  if(is.null(x$result)) {
    return(NULL)
  }
  
  return(x$result)
}

# ----------------------------------------------------------------------------
# True segments
# ----------------------------------------------------------------------------

true_segments <- segment_data %>%
  # mutate(
  #   data = map(true_dat, extract_safe) # true_segments
  # ) %>%
  select( run_id, r_rate, true_dat
  ) %>%
  unnest(true_dat)

# ----------------------------------------------------------------------------
# Inferred segments
# ----------------------------------------------------------------------------

inferred_segments <- segment_data %>%
  mutate(
    data = map(inferred_segments, extract_safe)
  ) %>%
  select( run_id, r_rate, data
  ) %>%
  unnest(data)

# ============================================================================
# 9. Format Recombination Labels
# ============================================================================

format_rates <- function(df) {
  
  df %>%
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
}

true_segments <- format_rates(true_segments)
inferred_segments <- format_rates(inferred_segments)

# ============================================================================
# 10. Panel A
# TRUE IBD Tract Length Distributions
# ============================================================================

pA <- true_segments %>%
  ggplot(
    aes(x = max_segment_bp, # segment_length, total_ibd_bp
      fill = r_rate_label)
  ) +
  
  geom_density(alpha = 0.7) +
  
  scale_x_log10(labels = scales::label_scientific()) + # scientific annotation
  # scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + # Mathematical annotation
  
  facet_wrap( ~ r_rate_label, scales = "free_y") +
  
  theme_bw() +
  
  labs(
    x = "True IBD tract length (log scale)",
    y = "Density") +
  
  theme(
    legend.position = "none",
    plot.title = element_text( face = "bold", size = 15),
    plot.title.position = "plot",
    axis.title = element_text(size = 15, color = 'black', face = 'bold'),
    axis.text = element_text(size = 14, color = 'black'),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm"),
    strip.text = element_text(size = 15, face = 'bold')
  )

print(pA)

# ============================================================================
# 11. Panel B
# INFERRED IBD Tract Length Distributions
# ============================================================================

pB <- inferred_segments %>%
  ggplot(
    aes(
      x = max_segment_bp,
      fill = r_rate_label
    )
  ) +
  
  geom_density(
    alpha = 0.5
  ) +
  
  scale_x_log10() +
  
  facet_wrap(
    ~ r_rate_label,
    scales = "free_y"
  ) +
  
  theme_bw() +
  
  labs(
    x = "Inferred IBD tract length (log scale)",
    y = "Density",
    title = "B. Fragmentation of inferred IBD tracts"
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

# ============================================================================
# 12. Compute Pairwise Detection Status
#
# True Positive  = inferred AND true
# False Positive = inferred only
# False Negative = true only
# ============================================================================
# ----------------------------------------------------------------------------
# True pair presence
# ----------------------------------------------------------------------------

true_pairs <- true_segments %>%
  distinct( run_id, r_rate, pair_key) %>%
  mutate(true_ibd = TRUE)

# ----------------------------------------------------------------------------
# Inferred pair presence
# ----------------------------------------------------------------------------

inferred_pairs <- inferred_segments %>%
  distinct(run_id, r_rate, pair_key) %>%
  mutate(inferred_ibd = TRUE)

# ----------------------------------------------------------------------------
# Merge pair status
# ----------------------------------------------------------------------------

pair_status <- full_join(
  true_pairs,
  inferred_pairs,
  by = c(
    "run_id",
    "r_rate",
    "pair_key"
  )
) %>%
  mutate(
    
    true_ibd = replace_na(
      true_ibd,
      FALSE
    ),
    
    inferred_ibd = replace_na(
      inferred_ibd,
      FALSE
    ),
    
    classification = case_when(
      
      true_ibd & inferred_ibd
      ~ "True Positive",
      
      !true_ibd & inferred_ibd
      ~ "False Positive",
      
      true_ibd & !inferred_ibd
      ~ "False Negative",
      
      TRUE
      ~ "True Negative"
    )
  ) %>%
  format_rates()

# ============================================================================
# 13. Compute Error Rates
# ============================================================================

error_rates <- pair_status %>%
  group_by(
    r_rate_numeric,
    r_rate_label
  ) %>%
  summarize(
    
    false_positive_rate = mean(
      classification == "False Positive"
    ),
    
    false_negative_rate = mean(
      classification == "False Negative"
    ),
    
    .groups = "drop"
  )

# ============================================================================
# 14. Panel C
# False Positive Rates
# ============================================================================

pC <- error_rates %>%
  ggplot(
    aes(
      r_rate_numeric,
      false_positive_rate
    )
  ) +
  
  geom_line(
    linewidth = 1.2,
    color = "firebrick"
  ) +
  
  geom_point(
    size = 3,
    color = "firebrick"
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
    y = "False-positive rate",
    title = "C. Inflation of false positives at low recombination"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

# ============================================================================
# 15. Panel D
# False Negative Rates
# ============================================================================

pD <- error_rates %>%
  ggplot(
    aes(
      r_rate_numeric,
      false_negative_rate
    )
  ) +
  
  geom_line(
    linewidth = 1.2,
    color = "steelblue"
  ) +
  
  geom_point(
    size = 3,
    color = "steelblue"
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
    y = "False-negative rate",
    title = "D. Fragmentation-induced false negatives"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

# ============================================================================
# 16. Combine Figure Panels
# ============================================================================

final_figure <- (
  pA | pB
) / (
  pC | pD
)

# ============================================================================
# 17. Save Figure
# ============================================================================

ggsave(
  filename = "results/figures/main/figure4_recombination_fragmentation.pdf",
  plot = final_figure,
  width = 16,
  height = 10,
  dpi = 600
)

# ============================================================================
# 18. Export Summary Tables
# ============================================================================

write_tsv(
  error_rates,
  "results/tables/recombination_error_rates.tsv"
)

# ============================================================================
# 19. Save Full Workspace
# ============================================================================

saveRDS(
  list(
    true_segments = true_segments,
    inferred_segments = inferred_segments,
    pair_status = pair_status,
    error_rates = error_rates
  ),
  "results/rds/figure4_recombination_fragmentation.rds"
)

###############################################################################
# End of Script
###############################################################################