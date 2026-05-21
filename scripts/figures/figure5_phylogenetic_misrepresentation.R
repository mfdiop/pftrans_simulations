

# Goal of Figure 5
# Show how phylogenetic metrics (cophenetic distance, patristic distance, tree clustering) fail to recover direct transmission events across baseline, recombination, and migration scenarios.
# 
# Panels
# A. Distribution of true TMRCAs for direct vs non-direct transmission pairs
# B. Distribution of phylogenetic distances (e.g., cophenetic) for direct vs non-direct pairs
# C. Scatter: true IBD fraction vs phylogenetic distance
# D. AUROC/PR curves for phylogenetic distance predicting “direct transmission”
# E. Supplementary: Misclassification heatmaps per method/scenario (separate)


library(tidyverse)
library(ape)
library(phangorn)
library(pROC)
library(patchwork)

scenarios <- "multiple_runs"

# Function to load one replicate
load_replicate_phylo <- function(scenario, rep_id) {
  
  tree_path <- glue::glue("simulations/{scenario}/phylo_output/replicate{rep_id}/run{rep_id}_chr1_modelfinder.treefile")
  true_path <- glue::glue("simulations/{scenario}/inferred/replicate{rep_id}/true_ibd_summary.tsv")
  ibd_path  <- glue::glue("simulations/{scenario}/inferred/replicate{rep_id}/inferred_ibd_hmm.tsv")
  
  # Load data
  tree <- read.tree(tree_path)
  df_true <- read_tsv(true_path, show_col_types = FALSE)
  df_ibd  <- read_tsv(ibd_path, show_col_types = FALSE) %>% 
    rename(., id1 = Id1, id2 = Id2)
  
  # Compute cophenetic distances
  phylo_dist <- cophenetic(tree) |> as.data.frame() |> tibble::rownames_to_column("id1")
  
  phylo_long <- phylo_dist %>%
    pivot_longer(-id1, names_to = "id2", values_to = "cophenetic_dist") %>%
    mutate(pair_key = paste(id1, id2, sep = "_"))
  
  # Merge with true link info
  df <- df_true %>%
    select(-Id1, -Id2) %>% 
    separate(., pair, into = c("id1", "id2"), sep = "_") %>% 
    mutate(id1 = paste0("tsk_", id1),
           id2 = paste0("tsk_", id2),
           pair_key = paste(id1, id2, sep = "_")) %>%
    left_join(df_ibd, by = c("id1", "id2")) %>%
    # left_join(phylo_long, by = "pair_key") %>%
    left_join(phylo_long, by = c("id1", "id2")) %>%
    mutate(
      scenario = scenario,
      replicate = rep_id)
  
  return(df)
}

load_multiple_phylo <- function(scenario, rep_id) {
  
  tree_path <- glue::glue("simulations/{scenario}/metrics/phylo_results/rep{rep_id}/run{rep_id}_rec1e06_chr1.treefile")
  true_path <- glue::glue("simulations/{scenario}/metrics/inferred/rep{rep_id}/run{rep_id}_rec1e06_chr1/true_ibd_summary.tsv")
  ibd_path  <- glue::glue("simulations/{scenario}/metrics/inferred/rep{rep_id}/run{rep_id}_rec1e06_chr1/inferred_ibd_hmm.tsv")
  
  # Load data
  tree <- read.tree(tree_path)
  df_true <- read_tsv(true_path, show_col_types = FALSE)
  df_ibd  <- read_tsv(ibd_path, show_col_types = FALSE) %>% 
    rename(., id1 = Id1, id2 = Id2)
  
  # Compute cophenetic distances
  phylo_dist <- cophenetic(tree) |> as.data.frame() |> tibble::rownames_to_column("id1")
  
  phylo_long <- phylo_dist %>%
    pivot_longer(-id1, names_to = "id2", values_to = "cophenetic_dist") %>%
    mutate(pair_key = paste(id1, id2, sep = "_"))
  
  # Merge with true link info
  df <- df_true %>%
    select(-Id1, -Id2) %>% 
    separate(., pair, into = c("id1", "id2"), sep = "_") %>% 
    mutate(id1 = paste0("tsk_", id1),
           id2 = paste0("tsk_", id2),
           pair_key = paste(id1, id2, sep = "_")) %>%
    left_join(df_ibd, by = c("id1", "id2")) %>%
    # left_join(phylo_long, by = "pair_key") %>%
    left_join(phylo_long, by = c("id1", "id2")) %>%
    mutate(
      scenario = scenario,
      replicate = rep_id)
  
  return(df)
}


# scenarios <- "single_run"
# # Load all replicates
# df_single <- map_dfr(
#   scenarios,
#   ~ map_dfr(1:5, \(x) load_replicate_phylo(.x, x))
# )

df_multiple <- map_dfr(
  scenarios,
  ~ map_dfr(1:5, \(x) load_multiple_phylo(.x, x))
)


df_phylo <- df_multiple %>% # scenarios <- "multiple_runs"
  mutate(scenario = recode(scenario, "multiple_runs" = "Recombination\nSweep"),
         direct = if_else(min_tmrca<=5, 1, 0),
         transmission_class = case_when(
           min_tmrca == 1 ~ "Direct (1 gen)",
           min_tmrca <= 5 ~ "Near-direct (≤5 gen)",
           min_tmrca > 5 ~ "Indirect (>5 gen)"),
         is_direct = min_tmrca == 1) %>%
  rename(tmrca_true = min_tmrca)
 
# fig5A <- df_phylo %>%
#   ggplot(aes(x = tmrca_true, fill = factor(direct))) +
#   geom_density(alpha = 0.6) +
#   facet_wrap(~ scenario,
#              labeller = as_labeller(c(baseline = "Baseline", 
#                                       recombination_sweep = "Recombination\nsweep"))) +
#   scale_fill_manual(values = c("0" = "grey60", "1" = "blue3"),
#                     # values = c("0" = "grey60", "1" = "red3"),
#                     labels = c("Non-direct", "Direct")) +
#   labs(
#     x = "True TMRCA (generations)",
#     y = "Density",
#     fill = "Transmission",
#     # title = "Figure 5A. TMRCA distributions for direct vs non-direct pairs"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 18, color = "black", face = "bold"),
#     axis.text = element_text(color = "black", size = 15),
#     legend.text = element_text(size = 14, color = "black", face = "bold"),
#     legend.title = element_text(size = 13, color = "black", face = "bold"),
#     strip.text = element_text(size = 15, face = 'bold'))
# 
# ggsave("simulations/main/figure5_tmrca_distribution.png", 
#        plot = fig5A, width = 10, height = 6, dpi = 600)
#
# fig5A <- df_phylo %>%
#   ggplot(aes(x = tmrca_true, fill = factor(transmission_class))) +
#   geom_density(alpha = 0.4) +
#   facet_wrap(~ scenario,
#              labeller = as_labeller(c(baseline = "Baseline", 
#                                       recombination_sweep = "Recombination\nsweep"))) +
#   scale_fill_manual(values = c("0" = "grey60", "1" = "red3"),
#                     labels = c("Non-direct", "Direct")) +
#   labs(
#     x = "True TMRCA (generations)",
#     y = "Density",
#     fill = "Transmission",
#     title = "Figure 5A. TMRCA distributions for direct vs non-direct pairs"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 14, color = "black", face = "bold"),
#     axis.text = element_text(color = "black", size = 12),
#     legend.text = element_text(size = 11, color = "black", face = "bold"),
#     legend.title = element_text(size = 13, color = "black", face = "bold"),
#     strip.text = element_text(size = 12, face = 'bold'))

# ========================================
#   Visualise distance distribution
# ========================================
figure1 <- df_phylo %>%
  mutate(direct = factor(direct)) %>% # , levels = c("Direct", "Indirect")
  ggplot(aes(x = cophenetic_dist, fill = direct)) +
  geom_density(alpha = 0.6) +
  # facet_wrap(~ scenario,
  #            labeller = as_labeller(c(baseline = "Baseline", 
  #                                     recombination_sweep = "Recombination\nsweep"))) +
  scale_fill_manual(values = c("0" = "grey40", "1" = "blue3"),
                    labels = c("Indirect", "Direct")) +
  labs(
    x = "Cophenetic distance",
    y = "Density",
    fill = "Transmission Class"
    # title = "Figure 5B. Phylogenetic distances for direct vs non-direct pairs"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    axis.text = element_text(color = "black", size = 14),
    legend.text = element_text(size = 14, color = "black", face = "bold"),
    legend.title = element_text(size = 13, color = "black", face = "bold"),
    strip.text = element_text(size = 15, face = 'bold'))

ggsave("simulations/main/figure5_phylogenetic_distances.png", 
       plot = fig5B, width = 10, height = 6, dpi = 600)


fig5C <- df_phylo %>%
  ggplot(aes(x = true_ibd_prop, y = cophenetic_dist, color = factor(direct))) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("0" = "grey40", "1" = "blue"),
                     labels = c("Indirect", "Direct")) +
  labs(
    x = "True IBD fraction",
    y = "Cophenetic distance",
    color = "Transmission"
    # title = "Figure 5C. Relationship between IBD and phylogenetic distance"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    axis.text = element_text(color = "black", size = 14),
    legend.text = element_text(size = 15, color = "black", face = "bold"),
    legend.title = element_text(size = 13, color = "black", face = "bold"),
    strip.text = element_text(size = 15, face = 'bold'))

ggsave("simulations/main/figure5_ibd_phylo.png", 
       plot = fig5C, width = 10, height = 6, dpi = 600)

# safe_auc <- function(y, score) {
#   if (length(unique(y[!is.na(y)])) < 2) return(NA_real_)
#   as.numeric(pROC::roc(y, score, quiet = TRUE)$auc)
# }
# 
# df_auc <- df_phylo %>%
#   group_by(scenario, replicate) %>%
#   summarize(
#     auc = safe_auc(direct, -cophenetic_dist),
#     .groups = "drop"
#   )

safe_prauc <- function(y, score) {
  y <- y[!is.na(score)]
  score <- score[!is.na(score)]
  
  if (length(unique(y)) < 2) return(NA_real_)
  if (sum(y == 1) < 2) return(NA_real_)
  
  pr <- PRROC::pr.curve(
    scores.class0 = score[y == 1],
    scores.class1 = score[y == 0],
    curve = FALSE
  )
  
  pr$auc.integral
}

df_prauc <- df_phylo %>%
  group_by(scenario, replicate) %>%
  summarize(
    pr_auc = safe_prauc(direct, -cophenetic_dist),
    .groups = "drop"
  )

# Summarize across replicates
df_prauc_summary <- df_prauc %>%
  group_by(scenario) %>%
  summarize(
    median_prauc = median(pr_auc, na.rm = TRUE),
    iqr_lo = quantile(pr_auc, 0.25, na.rm = TRUE),
    iqr_hi = quantile(pr_auc, 0.75, na.rm = TRUE),
    n_na = sum(is.na(pr_auc)),
    .groups = "drop"
  )


fig5D <- df_prauc %>%
  ggplot(aes(x = scenario, y = pr_auc)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_x_discrete(labels = c( baseline = "Baseline",
                               recombination_sweep = "Recombination\nsweep")) +
  ylim(0, 1) +
  labs(
    x = "",
    y = "PRAUC"
    # title = "Figure 5D. AUROC of phylogenetic distance predicting direct transmission"
    ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(color = "black", size = 13))

ggsave("simulations/main/figure5_prauc_phylo.png", 
       plot = fig5D, width = 10, height = 6, dpi = 600)


# Save combined plot
figure5 <- (fig5A / fig5B | fig5C / fig5D) +
  plot_layout(heights = c(1, 1, 1, 1))

figure5

ggsave("simulations/main/figure5_phylogenetic_misrepresentation.pdf", 
       plot = figure5, width = 12, height = 10, dpi = 600)

# =========================================================
# 7. SUPPLEMENTARY FIGURE 5E — Misclassification heatmaps
# =========================================================

supp5E <- df_phylo %>%
  mutate(predicted_direct = cophenetic_dist < median(cophenetic_dist)) %>%
  count(scenario, replicate, direct, predicted_direct) %>%
  ggplot(aes(x = factor(direct), y = factor(predicted_direct), fill = n)) +
  geom_tile() +
  facet_grid(replicate ~ scenario) +
  scale_fill_viridis_c() +
  labs(
    x = "True",
    y = "Predicted",
    title = "Supplementary Figure 5E. Misclassification matrices per replicate"
  ) +
  theme_bw()




