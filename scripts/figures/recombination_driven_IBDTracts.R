#!/usr/bin/env Rscript

###############################################################################
# Figure 4: Recombination-Driven Fragmentation of IBD Tracts
#
# Author: YOUR NAME
# Date: YYYY-MM-DD
#
# DESCRIPTION
# ===========
# This script investigates how recombination reshapes genomic relatedness
# by fragmenting IBD (Identity-by-Descent) tracts across transmission depths.
#
# BIOLOGICAL HYPOTHESIS
# =====================
#
# LOW recombination:
#   - preserves long ancestral haplotypes,
#   - indirect ancestry remains genomically similar,
#   - recent and ancient transmission become difficult to distinguish.
#
# HIGH recombination:
#   - fragments ancient ancestry into shorter tracts,
#   - reduces background relatedness,
#   - preserves recent transmission tracts,
#   - improves distinguishability of recent transmission.
#
# Therefore:
#
#   Ancient ancestry should fragment faster than recent transmission.
#
# This script quantifies:
#
#   1. IBD tract length distributions
#   2. Number of tracts per pair
#   3. Fragmentation indices
#   4. False-positive IBD rates
#   5. False-negative IBD rates
#
# stratified by:
#
#   - recombination regime
#   - transmission depth (TMRCA)
#
# INPUT FILES
# ===========
#
# Expected directory structure:
#
# simulations/multiple_runs/metrics/inferred/
# ├── rep1/
# │   └── run1_rec1e09_chr1/
# │       ├── true_ibd_segments.tsv
# │       ├── inferred_ibd_segments.tsv
# │       └── true_ibd_summary.tsv
# ├── rep2/
# └── ...
#
# REQUIRED FILES
# ==============
#
# true_ibd_segments.tsv
# ---------------------
# Id1
# Id2
# start
# end
#
# inferred_ibd_segments.tsv
# -------------------------
# Id1
# Id2
# start
# end
#
# true_ibd_summary.tsv
# --------------------
# pair
# min_tmrca
#
# OUTPUTS
# =======
#
# Figure:
# results/figures/main/figure4_ibd_fragmentation.png
#
# Tables:
# results/tables/fragmentation_summary.tsv
# results/tables/error_rates.tsv
#
###############################################################################

# ============================================================================
# 1. Load Packages
# ============================================================================

library(tidyverse)
library(patchwork)
library(glue)

# ============================================================================
# 2. Helper Functions
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
      pair_key = paste(
        sort(c(
          get(id1),
          get(id2)
        )),
        collapse = "_"
      )
    ) %>%
    ungroup()
}

# ----------------------------------------------------------------------------
# Compute segment lengths
# ----------------------------------------------------------------------------

compute_segment_length <- function(df) {
  
  df %>%
    mutate(
      segment_length = end - start
    )
}

# ----------------------------------------------------------------------------
# Format recombination labels
# ----------------------------------------------------------------------------

format_recombination <- function(df) {
  
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

# ----------------------------------------------------------------------------
# Transmission class annotation
# ----------------------------------------------------------------------------

annotate_tmrca <- function(df) {
  
  df %>%
    mutate(
      
      transmission_class = case_when(
        
        min_tmrca == 1
        ~ "Direct",
        
        min_tmrca <= 2
        ~ "Near-direct",
        
        TRUE
        ~ "Indirect"
      ),
      
      transmission_class = factor(
        transmission_class,
        levels = c(
          "Direct",
          "Near-direct",
          "Indirect"
        )
      )
    )
}

# ============================================================================
# 3. Define Directories
# ============================================================================

root <- "simulations/multiple_runs/metrics/inferred"

rep_dirs <- list.dirs(
  root,
  recursive = TRUE,
  full.names = TRUE
)

# ============================================================================
# 4. Parse Metadata
# ============================================================================

parsed_dirs <- tibble(
  path = rep_dirs) %>%
  mutate(
    
    folder = basename(path),
    
    run_id = str_extract(
      folder,
      "(?<=run)\\d+"
    ),
    
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
# 5. Load TMRCA Metadata
# ============================================================================

load_tmrca <- function(path) {
  
  file <- file.path(
    path,
    "true_ibd_summary.tsv"
  )
  
  read_tsv(
    file,
    show_col_types = FALSE
  ) %>%
    separate(
      pair,
      into = c("id1", "id2"),
      sep = "_"
    ) %>%
    mutate(
      Id1 = paste0("tsk_", id1),
      Id2 = paste0("tsk_", id2)
    ) %>%
    select(
      Id1,
      Id2,
      min_tmrca
    ) %>%
    standardize_pairs() %>%
    annotate_tmrca()
}

# ============================================================================
# 6. Load TRUE Segments
# ============================================================================

load_true_segments <- function(path) {
  
  file <- file.path(
    path,
    "true_ibd_summary.tsv"
  )
  
  read_tsv(file, show_col_types = FALSE) %>%
    standardize_pairs() %>% 
    annotate_tmrca()
    # compute_segment_length()
}


# ============================================================================
# 7. Load INFERRED Segments
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
# 8. Load All Data
# ============================================================================

full_data <- parsed_dirs %>%
  mutate(
    
    # tmrca = map(path, safely(load_tmrca)),
    
    true_segments = map( path, safely(load_true_segments))
    
    # inferred_segments = map(
    #   path,
    #   safely(load_inferred_segments)
    # )
  )

# ============================================================================
# 9. Safe Extraction Helper
# ============================================================================

extract_safe <- function(x) {
  
  if(is.null(x$result)) {
    return(NULL)
  }
  
  return(x$result)
}

# ============================================================================
# 10. Extract TMRCA Metadata
# ============================================================================

tmrca_data <- full_data %>%
  mutate(
    data = map(tmrca, extract_safe)
  ) %>%
  select(
    run_id,
    r_rate,
    data
  ) %>%
  unnest(data)

# ============================================================================
# 11. Extract TRUE Segments
# ============================================================================

true_segments <- full_data %>%
  mutate(
    data = map(true_segments, extract_safe)
  ) %>%
  select(
    run_id,
    r_rate,
    data
  ) %>%
  unnest(data)

# ============================================================================
# 12. Extract INFERRED Segments
# ============================================================================

inferred_segments <- full_data %>%
  mutate(
    data = map(inferred_segments, extract_safe)
  ) %>%
  select(
    run_id,
    r_rate,
    data
  ) %>%
  unnest(data)

# ============================================================================
# 13. Add Recombination Metadata
# ============================================================================

true_segments <- annotate_tmrca(true_segments)

true_segments <- format_recombination(true_segments)

# tmrca_data <- format_recombination(tmrca_data)
# 
# inferred_segments <- format_recombination(inferred_segments)

# ============================================================================
# 14. Join Transmission Classes to TRUE Segments
# ============================================================================

true_segments <- true_segments %>%
  left_join(
    tmrca_data %>%
      select(
        pair_key,
        transmission_class,
        min_tmrca
      ),
    by = "pair_key"
  )

# ============================================================================
# 15. TRUE IBD TRACT LENGTH DISTRIBUTIONS
#
# Expected:
# LOW recombination:
#   indirect ancestry retains long tracts
#
# HIGH recombination:
#   indirect ancestry fragments
#
# DIRECT transmission:
#   retains longer tracts
# ============================================================================

pA <- true_segments %>%
  ggplot(
    aes(
      x = total_ibd_bp,
      fill = transmission_class
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
    x = "True IBD tract length (log scale)",
    y = "Density",
    fill = "Transmission class",
    title = "A. Recombination fragments ancient ancestry"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

# ============================================================================
# 16. NUMBER OF TRACTS PER PAIR
#
# High recombination should increase tract fragmentation:
# more tracts but shorter tracts.
# ============================================================================

tract_counts <- true_segments %>%
  group_by(
    pair_key,
    r_rate_label,
    transmission_class
  ) %>%
  summarize(
    n_tracts = n(),
    .groups = "drop"
  )

pB <- tract_counts %>%
  ggplot(
    aes(
      x = r_rate_label,
      y = n_tracts,
      fill = transmission_class
    )
  ) +
  
  geom_violin(
    alpha = 0.5,
    scale = "width"
  ) +
  
  geom_boxplot(
    width = 0.1,
    outlier.alpha = 0
  ) +
  
  theme_bw() +
  
  labs(
    x = "Recombination rate",
    y = "Number of tracts per pair",
    fill = "Transmission class",
    title = "B. High recombination increases tract fragmentation"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

# ============================================================================
# 17. FRAGMENTATION INDEX
#
# Fragmentation Index:
#
# number_of_tracts / total_tract_length
#
# High values indicate highly fragmented ancestry.
# ============================================================================

fragmentation_index <- true_segments %>%
  group_by(
    pair_key,
    r_rate_label,
    transmission_class
  ) %>%
  summarize(
    
    total_length = sum(
      total_ibd_bp,
      na.rm = TRUE
    ),
    
    n_tracts = n(),
    
    fragmentation_index =
      n_tracts / total_length,
    
    .groups = "drop"
  )

pC <- fragmentation_index %>%
  ggplot(
    aes(
      x = r_rate_label,
      y = fragmentation_index,
      fill = transmission_class
    )
  ) +
  
  geom_violin(
    alpha = 0.5,
    scale = "width"
  ) +
  
  geom_boxplot(
    width = 0.1,
    outlier.alpha = 0
  ) +
  
  scale_y_log10() +
  
  theme_bw() +
  
  labs(
    x = "Recombination rate",
    y = "Fragmentation index (log scale)",
    fill = "Transmission class",
    title = "C. Fragmentation of indirect ancestry under recombination"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

# ============================================================================
# 18. FALSE POSITIVE / FALSE NEGATIVE ANALYSIS
# ============================================================================

# ----------------------------------------------------------------------------
# True pair presence
# ----------------------------------------------------------------------------

true_pairs <- true_segments %>%
  distinct(
    run_id,
    r_rate,
    pair_key
  ) %>%
  mutate(
    true_ibd = TRUE
  )

# ----------------------------------------------------------------------------
# Inferred pair presence
# ----------------------------------------------------------------------------

inferred_pairs <- inferred_segments %>%
  distinct(
    run_id,
    r_rate,
    pair_key
  ) %>%
  mutate(
    inferred_ibd = TRUE
  )

# ----------------------------------------------------------------------------
# Merge detection status
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
  format_recombination()

# ----------------------------------------------------------------------------
# Compute error rates
# ----------------------------------------------------------------------------

error_rates <- pair_status %>%
  group_by(
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
# 19. FALSE POSITIVE VS FALSE NEGATIVE
#
# Expected:
#
# LOW recombination:
#   inflated false positives
#
# HIGH recombination:
#   elevated false negatives
# ============================================================================

error_long <- error_rates %>%
  pivot_longer(
    cols = c(
      false_positive_rate,
      false_negative_rate
    ),
    names_to = "error_type",
    values_to = "rate"
  )

pD <- error_long %>%
  ggplot(
    aes(
      x = r_rate_label,
      y = rate,
      color = error_type,
      group = error_type
    )
  ) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  geom_point(
    size = 3
  ) +
  
  theme_bw() +
  
  labs(
    x = "Recombination rate",
    y = "Error rate",
    color = "Error type",
    title = "D. Tradeoff between false positives and false negatives"
  ) +
  
  scale_color_manual(
    values = c(
      "false_positive_rate" = "firebrick",
      "false_negative_rate" = "steelblue"
    ),
    labels = c(
      "False Positive",
      "False Negative"
    )
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

# ============================================================================
# 20. Combine Final Figure
# ============================================================================

final_figure <- (
  pA | pB
) / (
  pC | pD
)

# ============================================================================
# 21. Save Figure
# ============================================================================

ggsave(
  filename =
    "results/figures/main/figure4_ibd_fragmentation.png",
  
  plot = final_figure,
  
  width = 18,
  height = 12,
  dpi = 600
)

# ============================================================================
# 22. Export Summary Tables
# ============================================================================

write_tsv(
  fragmentation_index,
  "results/tables/fragmentation_summary.tsv"
)

write_tsv(
  error_rates,
  "results/tables/error_rates.tsv"
)

# ============================================================================
# 23. Save Workspace
# ============================================================================

saveRDS(
  list(
    true_segments = true_segments,
    inferred_segments = inferred_segments,
    fragmentation_index = fragmentation_index,
    error_rates = error_rates
  ),
  
  "results/rds/figure4_ibd_fragmentation.rds"
)

###############################################################################
# END OF SCRIPT
###############################################################################