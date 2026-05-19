#!/usr/bin/env Rscript

# r_rate        # numeric recombination rate (e.g., 1e-9, 1e-8, 1e-7, 1e-6)
# replicate     # replicate ID (1…5 or 1…20)
# pair_id       # infection pair
# true_ibd      # pedigree-derived total IBD fraction for each pair
# hmm_ibd       # inferred IBD fraction for same pair

# 1. Load Packages
library(tidyverse)
library(boot)
library(patchwork)

# ==========================================
# 2. Helper Function for Bootstrapped CIs
# =========================================
boot_mean_ci <- function(x, R = 10000, conf = 0.95) {
  
  boot_obj <- boot::boot(
    data = x,
    statistic = function(data, idx) mean(data[idx], na.rm = TRUE),
    R = R
  )
  
  ci <- boot::boot.ci(boot_obj, conf = conf, type = "perc")
  
  tibble(
    mean = mean(x, na.rm = TRUE),
    lower = ci$percent[4],
    upper = ci$percent[5]
  )
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
  
  x <- standardize_pairs(x, ids)
  return(x)
}

safe_cor <- function(x, y) {
  if (sum(complete.cases(x, y)) < 2) return(NA_real_)
  cor(x, y, use = "complete.obs")
}


# ======================================================
# ==================== MAIN ============================
root <- "simulations/multiple_runs/metrics/inferred"
# root <- "results/multiple_runs/inferred/"
scenario_dirs <- list.dirs(root, full.names = TRUE, recursive = TRUE)

rec_rate <- c(1e-09, 1e-08, 1e-07, 1e-06)

results <- list()

id_cols <- c("Id1", "Id2")

parsed_dirs <- tibble(
  path = scenario_dirs) %>%
  mutate(
    folder = basename(path),
    run_id = str_extract(folder, "(?<=run)\\d+"),
    r_rate = str_extract(folder, "(?<=rec)[^_]+"),
    chr    = str_extract(folder, "(?<=chr)\\d+")
    # r_rate = parse_number(r_rate)  # Converts "1e-08" → numeric
  ) %>% 
  filter(chr == 1)


df_ibd <- parsed_dirs %>%
  mutate(
    true_dat = map(path, load_true_ibd, id_cols),
    hmm_dat  = map(path, load_hmm_ibd, id_cols)) %>%
  mutate(
    merged = map2(true_dat, hmm_dat, ~full_join(.x, .y, by = "pair_key"))
  ) %>%
  select(run_id, r_rate, chr, merged) %>%
  unnest(merged)
  
# =============================================
# 3. Prepare Summary Statistics Per Replicate
# 
# We compute mean IBD per recombination rate, 
# separately for each replicate. 
# Then we combine replicates for CIs.

summ_ibd <- df_ibd %>%
  mutate(r_rate = as.numeric(gsub("e", "e-", r_rate))) %>% 
  drop_na() %>% 
  group_by(run_id, r_rate, pair_key) %>%
  summarize(
    mean_true = mean(true_ibd_prop, na.rm = TRUE),
    mean_hmm  = mean(ibd, na.rm = TRUE),
    .groups = "drop")

# =============================================
# 4. Compute Bootstrapped CIs Across Replicates

ci_recomb <- summ_ibd %>%
  group_by(r_rate) %>%
  summarize(
    boot_true = list(boot_mean_ci(mean_true)),
    boot_hmm  = list(boot_mean_ci(mean_hmm)),
    .groups = "drop") %>%
  unnest_wider(boot_true, names_sep = "_true") %>%
  unnest_wider(boot_hmm,  names_sep = "_hmm")

# ======================================================
# 5. Panel A — Distribution of true IBD vs recombination
pA <- df_ibd %>%
  # mutate(r_rate = as.numeric(gsub("e", "e-", r_rate))) %>% 
  ggplot(aes(factor(r_rate), true_ibd_prop, fill = factor(r_rate))) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.1, outlier.alpha = 0) +
  # scale_x_discrete(labels = rec_rate) +
  theme_bw() +
  labs(
    x = "Recombination rate", y = "True IBD proportion",
    title = "A. True IBD distribution across recombination rates") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, color = 'black', face = 'bold'),
        axis.title = element_text(size = 15, color = 'black', face = 'bold'),
        axis.text = element_text(size = 14, color = 'black'),
        axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
        axis.ticks = element_line(color = 'black', linewidth = .7),
        axis.ticks.length = unit(.20, "cm")
  )

# ===========================================================
# 6. Panel B — Distribution of inferred IBD vs recombination
pB <- df_ibd %>%
  ggplot(aes(factor(r_rate), hmm, fill = factor(r_rate))) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.1, outlier.alpha = 0) +
  scale_x_discrete(labels = rate_label) +
  theme_bw() +
  labs(
    x = "Recombination rate", y = "Inferred IBD proportion (hmmIBD)",
    title = "B. Inferred IBD distribution across recombination rates") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, color = 'black', face = 'bold'),
        axis.title = element_text(size = 15, color = 'black', face = 'bold'),
        axis.text = element_text(size = 14, color = 'black'),
        axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
        axis.ticks = element_line(color = 'black', linewidth = .7),
        axis.ticks.length = unit(.20, "cm")
  )

# ========================================================
# 7. Panel C — Decay of mean IBD (true & inferred) with CI
pC <- ci_recomb %>%
  ggplot(aes(r_rate)) +
  
  # true IBD
  geom_ribbon(
    aes(ymin = lower_true, ymax = upper_true),
    fill = "steelblue", alpha = 0.2) +
  geom_line(
    aes(y = mean_true),
    color = "steelblue", linewidth = 1.2) +
  
  # inferred IBD
  geom_ribbon(
    aes(ymin = lower_hmm, ymax = upper_hmm),
    fill = "firebrick", alpha = 0.2) +
  geom_line(
    aes(y = mean_hmm),
    color = "firebrick", linewidth = 1.2) +
  
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Recombination rate (log scale)", y = "Mean IBD proportion",
    title = "C. Decay of true and inferred IBD across recombination rates") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, color = 'black', face = 'bold'),
        axis.title = element_text(size = 14, color = 'black', face = 'bold'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
        axis.ticks = element_line(color = 'black', linewidth = .7),
        axis.ticks.length = unit(.20, "cm")
  )
  

# =======================================================
# 8. Panel D — Correlation between true and inferred IBD

pD <- ci_recomb %>%
  ggplot(aes(r_rate)) +

  geom_ribbon(
    aes(ymin = lower_cor, ymax = upper_cor),
    fill = "grey70", alpha = 0.3) +
  geom_line(
    aes(y = mean_cor),
    color = "black", linewidth = 1.2) +

  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Recombination rate (log scale)", y = "Correlation (true vs inferred IBD)",
    title = "D. Correlation between true and inferred IBD") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, color = 'black', face = 'bold'),
        axis.title = element_text(size = 14, color = 'black', face = 'bold'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
        axis.ticks = element_line(color = 'black', linewidth = .7),
        axis.ticks.length = unit(.20, "cm")
  )

results <- list(ibd = df_ibd, 
                summary = summ_ibd, 
                boostrap = ci_recomb)

saveRDS("results/multiple_runs/figure2.rds")
  
# ===================
# 9. Combine Panels
# ===================
p_fig2 <- (pA | pB) / (pC | pD)

print(p_fig2)

# Save combined plot
# ggsave("simulations/main/figure2_ibd_recombination.png", p_fig2, width = 16, height = 10)
ggsave("results/multiple_runs/figure2_ibd_recombination.png", p_fig2, width = 16, height = 10)

# ======================================================
# Option 2. Replace undefined correlations with NA and 
# Compute correlation at the run × recombination rate level

# cor_ibd <- df_ibd %>%
#   group_by(run_id, r_rate) %>%
#   summarize(
#     cor_true_hmm = safe_cor(true_ibd_prop, hmm),
#     n = n(),
#     .groups = "drop"
#   )

# ====================================================
# Extended Structure for Figure 2
# Panels
# 
# A. LD decay across recombination rates
# 
# B. True IBD proportion distributions
# 
# C. Inferred IBD (hmmIBD) distributions
# 
# D. Method recovery accuracy (true vs inferred IBD, with bootstrap CI)
# 
# E. Directionality error rate (mis-oriented pairs: Is A→B predicted when true is B→A?)
























