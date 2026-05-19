
# After line 120 (after summ_ibd is created), ADD:

# Calculate correlation between true and inferred IBD per replicate
cor_per_replicate <- df_ibd %>%
  mutate(r_rate = as.numeric(gsub("e", "e-", r_rate))) %>%
  drop_na(true_ibd_prop, hmm) %>%  # Use 'ibd' or 'hmm' - check your column name
  group_by(run_id, r_rate) %>%
  summarize(
    correlation = safe_cor(true_ibd_prop, hmm),
    n_pairs = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(correlation))  # Remove replicates with undefined correlation

# Summarize correlation across replicates with bootstrap CI
ci_correlation <- cor_per_replicate %>%
  group_by(r_rate) %>%
  summarize(
    boot_cor = list(boot_mean_ci(correlation)),
    .groups = "drop"
  ) %>%
  unnest_wider(boot_cor, names_sep = "_cor") %>%
  rename(
    mean_cor = boot_cor_cormean,
    lower_cor = boot_cor_corlower,
    upper_cor = boot_cor_corupper
  )


# ============================================================
# DIAGNOSTIC 1: Raw scatterplots per recombination rate
# ============================================================

library(cowplot)

p_scatter <- df_ibd %>%
  mutate(r_rate = as.numeric(gsub("e", "e-", r_rate))) %>%
  drop_na(true_ibd_prop, hmm) %>%
  ggplot(aes(x = true_ibd_prop, y = hmm)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~r_rate, scales = "free", ncol = 2,
             labeller = labeller(r_rate = function(x) paste0("r = ", x))) +
  theme_bw() +
  labs(
    x = "True IBD proportion (from tree sequence)",
    y = "Inferred IBD proportion (hmmIBD)",
    title = "True vs. Inferred IBD - Raw Data"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave("diagnostics/true_vs_inferred_scatter.png", p_scatter, 
       width = 10, height = 10)

# ============================================================
# DIAGNOSTIC 2: Why doesn't red line track blue in Panel C?
# ============================================================

# Calculate bias (inferred - true) per recombination rate
bias_analysis <- df_ibd %>%
  mutate(r_rate = as.numeric(gsub("e", "e-", r_rate))) %>%
  drop_na(true_ibd_prop, hmm) %>%
  mutate(bias = hmm - true_ibd_prop) %>%
  group_by(r_rate) %>%
  summarize(
    mean_bias = mean(bias, na.rm = TRUE),
    median_bias = median(bias, na.rm = TRUE),
    sd_bias = sd(bias, na.rm = TRUE),
    .groups = "drop"
  )

print("=== BIAS ANALYSIS ===")
print(bias_analysis)

# Plot bias distribution
p_bias <- df_ibd %>%
  mutate(r_rate = as.numeric(gsub("e", "e-", r_rate))) %>%
  drop_na(true_ibd_prop, hmm) %>%
  mutate(bias = hmm - true_ibd_prop) %>%
  ggplot(aes(x = factor(r_rate), y = bias, fill = factor(r_rate))) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    x = "Recombination rate",
    y = "Bias (Inferred - True IBD)",
    title = "hmmIBD Bias Across Recombination Rates"
  ) +
  theme(legend.position = "none")

ggsave("diagnostics/hmm_ibd_bias.png", p_bias, width = 10, height = 6)

# ============================================================
# DIAGNOSTIC 3: Check variance at each recombination rate
# ============================================================

variance_check <- df_ibd %>%
  mutate(r_rate = as.numeric(gsub("e", "e-", r_rate))) %>%
  drop_na(true_ibd_prop, hmm) %>%
  group_by(r_rate) %>%
  summarize(
    var_true = var(true_ibd_prop, na.rm = TRUE),
    var_inferred = var(hmm, na.rm = TRUE),
    mean_true = mean(true_ibd_prop, na.rm = TRUE),
    mean_inferred = mean(hmm, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print("=== VARIANCE CHECK ===")
print(variance_check)

# ============================================================
# DIAGNOSTIC 4: Correlation per replicate distribution
# ============================================================

p_cor_dist <- cor_per_replicate %>%
  ggplot(aes(x = factor(r_rate), y = correlation, fill = factor(r_rate))) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    x = "Recombination rate",
    y = "Correlation (true vs inferred IBD)",
    title = "Correlation Distribution Across Replicates"
  ) +
  theme(legend.position = "none")

ggsave("diagnostics/correlation_distribution.png", p_cor_dist, 
       width = 10, height = 6)
# Merge correlation into ci_recomb for Panel D
ci_recomb <- ci_recomb %>%
  left_join(ci_correlation, by = "r_rate")