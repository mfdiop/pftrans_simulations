#!/usr/bin/env Rscript


library(tidyverse)
library(patchwork)

#------------------------------------------------------------
# Helper: convert symmetric matrix to long vector of upper triangle
#------------------------------------------------------------
upper_vec <- function(mat) {
  mat[upper.tri(mat, diag = FALSE)]
}

standardize_pairs <- function(df, id_cols) {
  df_std <- df %>%
    rowwise() %>%
    mutate(
      pair_key = paste(sort(c(!!sym(id_cols[1]), !!sym(id_cols[2]))), collapse = "_")
    ) %>%
    ungroup() %>%
    distinct(pair_key, .keep_all = TRUE)  # Remove duplicates
  
  return(df_std)
}


load_true_ibd <- function(path, ids) {
  file <- file.path(path, "true_ibd_summary.tsv")
  x <- read_tsv(file, show_col_types = FALSE) %>% 
    mutate(Id1 = paste0("tsk_", Id1),
           Id2= paste0("tsk_", Id2))
  # select(pair = pair, true_ibd = true_ibd_prop)
  x <- standardize_pairs(x, ids)
  return(x)
}


load_hmm_ibd <- function(path, ids) {
  file <- file.path(path, "inferred_ibd_hmm.rds")
  x <- readRDS(file)
  # x <- x$fract[, c("sample1", "sample2", "fract_sites_IBD")]
  colnames(x) <- c("Id1", "Id2", "ibd")
  x <- standardize_pairs(x, ids)
  return(x)
}

"simulations/multiple_runs/inferred/rep1/run1_rec1e09_chr1/inferred_ibd_hmm.tsv"

#------------------------------------------------------------
# Load all replicate matrices and assemble tidy dataset
#------------------------------------------------------------
# base_dir <- "simulations/single_run/inferred"
base_dir <- "results/single_run/inferred/"
# output_dir <- "simulations/single_run/figures"
output_dir <- "results/single_run/inferred/"
reps <- 1:20

# df_ibd <- map_dfr(reps, function(r) {
#   
#   truth_mat <- read_tsv(file.path(base_dir, paste0("replicate", r), "true_ibd_summary.tsv"))
#   hmm_mat   <- readRDS(file.path(base_dir, paste0("replicate", r), "inferred_ibd_hmm.rds"))
#   
#   hmm_mat <- hmm_mat$fract[, c("sample1", "sample2", "fract_sites_IBD")]
#   
#   tibble(
#     replicate = r,
#     true_ibd = truth_mat,
#     hmm_ibd  = hmm_mat
#   )
# })


scenario_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)

id_cols <- c("Id1", "Id2")

parsed_dirs <- tibble(
  path = scenario_dirs) %>%
  mutate(folder = basename(path))

parsed_dirs_sorted <- parsed_dirs %>%
  mutate(rep_num = readr::parse_number(folder)) %>%   # extract numeric part
  arrange(rep_num) %>% 
  drop_na() %>% 
  filter(rep_num == 1) %>% 
  select(-rep_num)


df_ibd <- parsed_dirs_sorted %>%
  mutate(
    true_dat = map(path, load_true_ibd, id_cols),
    hmm_dat  = map(path, load_hmm_ibd, id_cols)) %>%
  mutate(merged = map2(true_dat, hmm_dat, ~full_join(.x, .y, by = "pair_key"))) %>%
  unnest(merged) %>% 
  drop_na()

#------------------------------------------------------------
# Plot 1: Histogram of true pairwise IBD distribution
#------------------------------------------------------------
p_true <- df_ibd %>%
  ggplot(aes(x = true_ibd_prop)) +
  geom_histogram(bins = 50, color = "black", fill = "steelblue", alpha = 0.7) +
  geom_vline(aes(xintercept = median(ibd)), color = "red", size = 1) +
  geom_vline(aes(xintercept = quantile(ibd, 0.25)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(ibd, 0.75)), linetype = "dashed") +
  labs(
    title = "A) True Pairwise IBD Fraction Distribution",
    x = "True IBD fraction", y = "Count") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(size = 16, color = 'black', face = 'bold'),
        axis.title = element_text(size = 15, color = 'black', face = 'bold'),
        axis.text = element_text(size = 14, color = 'black'),
        axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
        axis.ticks = element_line(color = 'black', linewidth = .7),
        axis.ticks.length = unit(.20, "cm"))

print(p_true)

#------------------------------------------------------------
# Plot 2: Violin plot of inferred IBD distribution (hmmIBD)
#------------------------------------------------------------
p_hmm <- df_ibd %>%
  ggplot(aes(x = factor(1), y = ibd)) +  # Keep x as ibd and y as constant factor
  geom_violin(fill = "orange", alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", size = 3, color = "red") +
  labs(
    title = "B) Inferred IBD Fraction (hmmIBD)",
    y = "Inferred IBD fraction", x = "") +  # Adjust labels accordingly
  coord_flip() +  # Flip coordinates
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(size = 14, color = 'black'),  # Show x-axis text again
        plot.title = element_text(size = 16, color = 'black', face = 'bold'),
        plot.title.position = "plot",
        axis.title = element_text(size = 15, color = 'black', face = 'bold'),
        axis.text.y = element_blank(),  # Hide y-axis text because it's a constant factor
        axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
        axis.ticks = element_line(color = 'black', linewidth = .7),
        axis.ticks.length = unit(.20, "cm"))


#------------------------------------------------------------
# Combine the panels
#------------------------------------------------------------
figure_baseline_ibd <- p_true / p_hmm +
  plot_annotation(title = "Baseline Scenario: True vs Inferred IBD Distributions")

ggsave("results/single_run/figure2_baseline_ibd.png", 
       figure_baseline_ibd, , width = 16, height = 10)
