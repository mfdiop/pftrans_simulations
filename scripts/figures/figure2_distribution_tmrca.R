
###############################################################
# Figure 2: Overlap of TMRCA Distributions Across
# Recombination Regimes
#
# Objective:
# ----------
# This script evaluates whether direct and indirect
# transmission pairs exhibit overlapping distributions
# of true TMRCA (time to most recent common ancestor)
# across multiple recombination regimes.
#
# Biological Interpretation:
# --------------------------
# If TMRCA distributions overlap substantially,
# genomic similarity alone may not uniquely identify
# direct transmission events.
#
# This figure therefore quantifies the emergence of
# identifiability limits in genomic transmission inference.
#
# Main Components:
# ----------------
# 1. Load phylogenetic trees and inferred pairwise data
# 2. Compute pairwise cophenetic distances
# 3. Merge true transmission metadata
# 4. Generate density distributions across recombination rates
# 5. Zoom into overlapping regions
# 6. Quantify overlap using overlap coefficient (OVL)
#
# Output:
# -------
# Figure 2A:
#   Global TMRCA density distributions
#
# Figure 2B:
#   Zoomed overlap region
#
# Figure 2C:
#   Overlap coefficient across recombination regimes
#
# Author: YOUR NAME
# Date: YYYY-MM-DD
###############################################################

# ============================================================
# Load Libraries
# ============================================================

library(tidyverse)
library(ape)
library(phangorn)
library(pROC)
library(patchwork)


load_multiple_phylo <- function(rec_rate, rep_id) {
  
  true_path <- glue::glue("simulations/multiple_runs/metrics/inferred/rep{rep_id}/run{rep_id}_{rec_rate}_chr1/true_ibd_summary.tsv")
  
  # Load data
  df_true <- read_tsv(true_path, show_col_types = FALSE)

  # Merge with true link info
  df <- df_true %>%
    select(-Id1, -Id2) %>% 
    separate(., pair, into = c("id1", "id2"), sep = "_") %>% 
    mutate(id1 = paste0("tsk_", id1),
           id2 = paste0("tsk_", id2),
           pair_key = paste(id1, id2, sep = "_")) %>%
    mutate( recombination_rate = rec_rate,
            replicate = rep_id)
  
  return(df)
}

# ============================================================
# Load All Replicates Across All Recombination Rates
# ============================================================

df_multiple <- map_dfr(
  recombination_scenarios,
  ~ map_dfr(1:5, \(x) load_multiple_phylo(.x, x))
)


# ============================================================
# Clean and Annotate Data
# ============================================================
df <- df_multiple %>% 
  rename(tmrca_true = min_tmrca) %>% 
  mutate(
    # Strict direct transmission definition
    # transmission_class = case_when(
    #   tmrca_true == 1 ~ "Direct",
    #   tmrca_true <= 3 ~ "Near-direct",
    #   tmrca_true >= 4 ~ "Indirect"
    # ),
    transmission_class = case_when(
      tmrca_true <= 5 ~ "Direct",
      tmrca_true > 5 ~ "Indirect"
    ),
    
    direct_binary = if_else(
      transmission_class == "Direct",
      "Direct",
      "Indirect"
    ),
    
    recombination_rate = factor(
      recombination_rate,
      levels = c(
        "rec1e09",
        "rec1e08",
        "rec1e07",
        "rec1e06"),
      # labels = c(
      #   expression(10^-9),
      #   expression(10^-8),
      #   expression(10^-7),
      #   expression(10^-6))
      labels = c("1e-09", "1e-08", "1e-07", "1e-06")
    )
  )


# ============================================================
# Figure 2A:
# Global Density Distributions
# ============================================================


figure2A <- df %>%
  # filter(recombination_rate == expression(10^-9)) %>%
  ggplot(aes(x = tmrca_true, fill = factor(direct_binary, levels = c("Direct", "Indirect")))) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ recombination_rate, scales = "free") +
  scale_fill_manual(values = c("Direct" = "blue3", "Indirect" = "grey60"),
                    labels = c("Direct", "Indirect")) +
  labs(
    x = "True TMRCA (generations)", y = "Density",
    fill = "Transmission Class") +
  scale_x_log10() +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    axis.text = element_text(color = "black", size = 15),
    legend.text = element_text(size = 14, color = "black", face = "bold"),
    legend.title = element_text(size = 13, color = "black", face = "bold"),
    strip.text = element_text(size = 15, face = 'bold'))

# ============================================================
# Figure 2B:
# Zoomed Overlap Region
#
# Focus on low-TMRCA region where ambiguity occurs
# ============================================================
inset <- df %>%
  ggplot(aes(x = tmrca_true, fill = factor(direct_binary, levels = c("Direct", "Indirect")))) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ recombination_rate, scales = "free") +
  scale_fill_manual(values = c("Direct" = "blue3", "Indirect" = "grey60"),
                    labels = c("Direct", "Indirect")) +
  coord_cartesian( xlim = c(0, 100)) +
  # scale_x_log10() +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_text(color = "black", size = 15),
    strip.text = element_text(size = 15, face = 'bold'))


ggsave("results/figures/main/figure2_tmrca_distribution.png", 
       plot = figure2, width = 10, height = 6, dpi = 600)


# ============================================================
# Compute Density Overlap Coefficient (OVL)
#
# OVL ranges:
#   0 = no overlap
#   1 = complete overlap
# ============================================================

compute_overlap <- function(x1, x2) {
  
  d1 <- density(x1, na.rm = TRUE)
  d2 <- density(x2, na.rm = TRUE)
  
  # Common x-grid
  x_common <- seq(
    max(min(d1$x), min(d2$x)),
    min(max(d1$x), max(d2$x)),
    length.out = 1000
  )
  
  y1 <- approx(d1$x, d1$y, xout = x_common)$y
  y2 <- approx(d2$x, d2$y, xout = x_common)$y
  
  overlap <- sum(pmin(y1, y2)) *
    diff(x_common[1:2])
  
  return(overlap)
}

# ============================================================
# Calculate Overlap Statistics
# ============================================================

overlap_stats <- df_phylo_all %>%
  filter(
    transmission_class %in% c("Direct", "Indirect")
  ) %>%
  group_by(recombination_rate) %>%
  summarise(
    
    overlap_coefficient = compute_overlap(
      tmrca_true[transmission_class == "Direct"],
      tmrca_true[transmission_class == "Indirect"]
    )
    
  )

print(overlap_stats)

# ============================================================
# Figure 2C:
# Overlap Coefficient by Recombination Rate
# ============================================================

figure2C <- overlap_stats %>%
  ggplot(
    aes(
      x = recombination_rate,
      y = overlap_coefficient,
      group = 1)) +
  
  geom_line(linewidth = 1.2) +
  
  geom_point(size = 4) +
  
  ylim(0, 1) +
  
  labs(
    x = "Recombination Rate",
    y = "Overlap Coefficient (OVL)") +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.title = element_text(face = "bold"))












