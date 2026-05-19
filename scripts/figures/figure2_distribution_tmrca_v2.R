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
library(cowplot)
library(ape)
library(phangorn)
library(glue)
library(patchwork)

# ============================================================
# Define Recombination Scenarios
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
# Function to Create Main + Inset Panel
# ============================================================

create_tmrca_panel <- function(data_subset, rec_label, logp = FALSE) {
  
  # ----------------------------------------------------------
  # Main plot
  # ----------------------------------------------------------
  
  p_main <- ggplot(
    data_subset,
    aes(
      x = tmrca_true,
      fill = factor(direct_binary, levels = c("Direct", "Indirect"))) ) +
    
    geom_density(alpha = 0.5) +
    
    scale_fill_manual(
      values = c(
        "Direct" = "blue3",
        "Indirect" = "grey70"
      )) +
    
    theme_bw(base_size = 13) +
    
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 15, color = "black", face = "bold"),
      axis.text  = element_text(color = "black", size = 13),
      axis.line  = element_line(linewidth = 0.8, color = "black", lineend = "square"),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.22, "cm"),
      # legend.text = element_text(size = 14, color = "black", face = "bold"),
      # legend.title = element_text(size = 13, color = "black", face = "bold")
      legend.position = "none"
    )
  
  if(logp) {
    p_main <- p_main + 
      scale_x_log10() +
      labs( title = rec_label,
        x = "True TMRCA (log scale)", y = "Density")
  }
  else{
    p_main <- p_main +
      labs( title = rec_label,
            x = "True TMRCA (generations)", y = "Density")
  }
  
  # ----------------------------------------------------------
  # Inset zoom plot
  # ----------------------------------------------------------
  
  p_inset <- ggplot(
    data_subset,
    aes(
      x = tmrca_true,
      fill = factor(direct_binary, levels = c("Direct", "Indirect")) )
  ) +
    
    geom_density(alpha = 0.5) +
    
    coord_cartesian(xlim = c(0, 100)) +
    
    scale_fill_manual(
      values = c(
        "Direct" = "blue3",
        "Indirect" = "grey70"
      )
    ) +
    
    labs( x = NULL, y = NULL) +
    
    theme_bw(base_size = 10) +
    
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 8),
      axis.title = element_blank(),
      plot.background = element_rect(color = "black", linewidth = 0.5)
    )
  
  # ----------------------------------------------------------
  # Combine main + inset
  # ----------------------------------------------------------
  
  combined_plot <- ggdraw() +
    
    draw_plot(p_main) +
    
    draw_plot(
      p_inset,
      x = 0.52,
      y = 0.45,
      width = 0.42,
      height = 0.42
    )
  
  return(combined_plot)
}

# ===================================
# Generate all recombination panels
# ===================================

panel_list <- df %>%
  split(.$recombination_rate) %>%
  imap(
    ~ create_tmrca_panel(.x, .y, logp = FALSE)
  )

# ===================================
# Combine panels
# ===================================
final_figure <- plot_grid(
  plotlist = panel_list,
  ncol = 2
  )

print(final_figure)

# ============================================================
# Save Figure
# ============================================================

ggsave(
  filename = "results/figures/main/figure2_tmrca_overlap.pdf",
  plot = final_figure,
  width = 10,
  height = 8,
  dpi = 600
)

# ====================
# Extract legend
# ====================
legend_plot <- ggplot(df, aes( x = tmrca_true, fill = factor(direct_binary, levels = c("Direct", "Indirect")))) +
  
  geom_density(alpha = 0.5) +
  
  scale_fill_manual( values = c( "Direct" = "blue3", "Indirect" = "grey70" ) ) +
  labs(fill = "Transmission Class") +
  theme_bw() +
  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14, color = "black", face = "bold"),
    legend.title = element_text(size = 13, color = "black", face = "bold")
  )

# shared_legend <- get_legend(legend_plot)
ggsave(
  filename = "results/figures/main/figure2_legend.pdf",
  plot = legend_plot,
  width = 5,
  dpi = 600
)

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

overlap_stats <- df %>%
  # filter(transmission_class %in% c("Direct", "Indirect")) %>%
  group_by(recombination_rate) %>%
  summarise(
    
    overlap_coefficient = compute_overlap(
      tmrca_true[transmission_class == "Direct"],
      tmrca_true[transmission_class == "Indirect"]
    )
    
  )

print(overlap_stats)

# ============================================================
# Export Overlap Statistics
# ============================================================

write_tsv(
  overlap_stats,
  "results/tables/tmrca_overlap_statistics.tsv"
)

###############################################################
# End of Script
###############################################################

