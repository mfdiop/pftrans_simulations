
# 4A: Performance across scenarios (F1 or PRAUC)

library(tidyverse)
library(patchwork)

custom_theme <- theme(
  axis.title = element_text(size = 14, face = "bold"),
  axis.text = element_text(size = 12),
  legend.text = element_text(size = 13, face = "bold"),
  legend.title = element_blank())

# Load data
single_results <- readRDS("simulations/single_run/evaluation_results.rds")

scenario_name <- "baseline"

baseline_summary <- imap_dfr(
  single_results,
  ~ extract_metrics(
    rep_obj = .x,
    scenario = scenario_name,
    run_id = .y))

baseline_summary <- baseline_summary %>% 
  add_column(rate = 6.6666667e-7) %>% 
  relocate(rate, , .after = run_id)

recomb <- readxl::read_xlsx("simulations/multiple_runs/recombination_evaluation/tables/evaluation_metrics.xlsx")

recomb <- recomb %>% 
  add_column(scenario = "recombination_sweep") %>% 
  relocate(scenario, , .before = replicate) %>% 
  rename(., run_id = replicate, precision = Precision, recall = Recall)

# ===================
#     Migration
# ===================
# scenario_name <- "structured_models"
migration <- read_csv("simulations/malaria_transmission_study/evaluation/migration_all_results.csv")

migration <- migration %>% 
  add_column(scenario = "migration") %>% 
  select(c(2, 4, 15, 22:25, 30)) %>%   
  relocate(scenario, .before = run_id) %>%
  relocate(method, .after = run_id) %>%
  rename(rate = rec_rate) %>% 
  mutate(
    precision = if_else(
      TP + FP == 0,
      NA_real_,
      TP / (TP + FP)),
    recall = if_else(
      TP + FN == 0,
      NA_real_,
      TP / (TP + FN)),
    f1 = if_else(
      is.na(precision) | is.na(recall) | (precision + recall) == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)))


# ================= COMBINE DATA FRAMES =====================
df_perf <- rbind.data.frame(baseline_summary, recomb, migration) %>%
  mutate(method = recode(method, "phylo" = "Phylo"),
         scenario = factor(scenario, levels = c("baseline", "recombination_sweep", "migration")))


df_perf_clean <- df_perf %>% 
  select(-run_id) %>%
  mutate(
    rate_num = as.numeric(rate),
    rate_clean = case_when(
      rate_num == 6.6666667e-07 ~ 6.67e-07,
      rate_num == 1e+06 ~ 1e-06,
      rate_num == 1e+07 ~ 1e-07,
      rate_num == 1e+08 ~ 1e-08,
      rate_num == 1e+09 ~ 1e-09
    # rate_clean = if_else(rate_num == 6.6666667e-07, 6.67e-07, rate_num)
    ))

fig4D <- df_perf %>%
  filter(scenario != "migration") %>%
  drop_na() %>% 
  ggplot(aes(x = rate, y = f1, color = method, group = method)) +
  stat_summary(fun = median, geom = "point", size = 3, show.legend = FALSE) +
  stat_summary(fun = median, geom = "line", linewidth = 1.5) +
  scale_x_log10(
    breaks = c(1e-09, 1e-08, 1e-07, 6.666667e-07, 1e-06),
    labels = scales::label_scientific()) + 
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), limits = c(0, 0.6)) +
  labs( x = "Recombination rate", y = "Precision–recall balance (F1 score)"  ) +
  custom_theme + theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

ggsave("simulations/main/figure4_f1_rcombinations.png", 
       plot = fig4D, width = 10, height = 6, dpi = 600)

# ====================================
# === CHECKING F1 DISTRIBUTION =======
df_perf_clean <- df_perf_clean %>%
  mutate( rate_log10 = log10(rate_clean))

distinct(df_perf_clean, rate_clean) %>% arrange(rate_clean)

# 1. Recall vs recombination (this directly addresses identifiability)
df_perf_clean %>%
  ggplot(aes(x = rate_clean, y = recall, color = method, group = method)) +
  stat_summary(fun = median, geom = "point", size = 2) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  scale_x_log10(
    breaks = c(1e-09, 1e-08, 1e-07, 6.67e-07, 1e-06),
    labels = scales::label_scientific()
  ) +
  labs(
    x = "Recombination rate",
    y = "Recall (TP / (TP + FN))",
    title = "Recall declines with increasing recombination"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_blank())

# 2. Number of predicted links (TP + FP)
df_perf_clean %>%
  mutate(predicted_links = TP + FP) %>%
  group_by(method, rate) %>% 
  summarise(N = sum(predicted_links)) %>%
  ggplot(aes(x = rate_clean, y = predicted_links, color = method, group = method)) +
  stat_summary(fun = median, geom = "point", size = 2) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  scale_x_log10(
    breaks = c(1e-09, 1e-08, 1e-07, 6.67e-07, 1e-06),
    labels = scales::label_scientific()) +
  labs(
    x = "Recombination rate",
    y = "Number of predicted transmission links",
    title = "Inference collapses to few predicted links at high recombination"
  ) +
  theme_bw()

# 3. PR curves with collapse highlighted
prcurve_dt <- read_tsv("simulations/multiple_runs/recombination_evaluation/tables/pr_curve.tsv")
  
prcurve_dt <- prcurve_dt %>% 
  add_column(scenario = "recombination_sweep") %>% 
  relocate(scenario, , .before = method) %>% 
  relocate(rate, , .after = method)
  
  
df_curve <- pr_baseline %>% 
  add_column(rate = 6.6666667e-7) %>% 
  select(-run_id) %>% 
  relocate(rate, , .after = method) %>% 
  rbind.data.frame(., prcurve_dt) %>%
  mutate(method = recode(method, "phylo" = "Phylo"))
  
df_curve %>%
  ggplot(aes(recall, precision, group = method, color = method)) + # interaction(run_id, rate)
  geom_path(alpha = 0.1) +
  facet_wrap(~ scenario) +
  labs(
    x = "Recall",
    y = "Precision",
    title = "Precision–Recall curves across recombination regimes") +
  coord_equal() +
  theme_bw()

# (b) Highlight collapse regime explicitly
df_curve %>%
  mutate(
    collapse = rate_num >= 6.67e-07
  ) %>%
  ggplot(aes(recall, precision, color = method)) +
  geom_path(
    aes(group = interaction(run_id, rate_num)),
    alpha = 0.15
  ) +
  geom_point(
    data = ~ filter(.x, collapse),
    size = 0.6,
    alpha = 0.6
  ) +
  facet_wrap(~ scenario) +
  labs(
    x = "Recall",
    y = "Precision",
    title = "Collapse of PR curves at high recombination"
  ) +
  coord_equal() +
  theme_bw()


# 4. TP and FN separately (this makes the story unavoidable)
# (a) True positives
df_perf_clean %>%
  ggplot(aes(x = rate_num, y = TP, color = method, group = method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_x_log10(
    breaks = c(1e-09, 1e-08, 1e-07, 6.67e-07, 1e-06),
    labels = scales::label_scientific() ) +
  labs(
    x = "Recombination rate",
    y = "True positives",
    title = "True positives collapse with recombination") +
  theme_bw()

# (b) False negatives
df_perf_clean %>%
  ggplot(aes(x = rate_num, y = FN, color = method, group = method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_x_log10(
    breaks = c(1e-09, 1e-08, 1e-07, 6.67e-07, 1e-06),
    labels = scales::label_scientific() ) +
  labs(
    x = "Recombination rate",
    y = "False negatives",
    title = "Missed transmission links increase with recombination" ) + theme_bw()


fig4E <- df_perf_all %>%
  select(method, scenario, replicate, TP, FP, FN) %>%
  pivot_longer(cols = c(TP, FP, FN), names_to = "error_type", values_to = "count") %>%
  ggplot(aes(x = scenario, y = count, fill = error_type)) +
  stat_summary(fun = median, geom = "bar", position = "stack") +
  facet_wrap(~ method) +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Count (median per replicate)",
    title = "E. Error mode decomposition")



fig4F <- df_pairs %>%
  ggplot(aes(x = true_IBD, y = inferred_IBD, color = method)) +
  geom_bin2d() +
  scale_fill_viridis_c() +
  facet_wrap(~ scenario) +
  theme_bw() +
  labs(
    x = "True pairwise IBD fraction",
    y = "Inferred IBD fraction",
    title = "F. Calibration of methods")

figure4 <- fig4A + fig4B + fig4C + fig4D + fig4E + plot_layout(ncol = 2)


# ========================================
# SUPPLEMENTARY
# ========================================
library(tidyverse)

supp4A <- df_perf_all %>%
  ggplot(aes(x = as.factor(replicate), y = F1, color = method)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_line(aes(group = method), alpha = 0.3) +
  facet_grid(scenario ~ method) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Replicate",
    y = "F1 Score",
    title = "Supplementary Figure 4A. Replicate-level performance across scenarios"
  )


supp4B <- df_pairs %>%
  ggplot(aes(x = true_IBD, y = inferred_IBD)) +
  geom_point(alpha = 0.15) +
  facet_grid(replicate ~ method) +
  theme_bw() +
  labs(
    x = "True IBD fraction",
    y = "Inferred IBD fraction",
    title = "Supplementary Figure 4B. Calibration scatter per replicate"
  )



supp4C <- df_perf_all %>%
  pivot_longer(
    cols = c(TP, FP, FN),
    names_to = "error_type",
    values_to = "count"
  ) %>%
  ggplot(aes(x = replicate, y = count, fill = error_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(scenario ~ method) +
  theme_bw() +
  labs(
    x = "Replicate",
    y = "Error Counts",
    title = "Supplementary Figure 4C. Error decomposition per replicate"
  )


supp4D <- df_perf_all %>%
  ggplot(aes(x = replicate, y = cor_true_inferred, color = method)) +
  geom_point(size = 2) +
  geom_line(aes(group = method), alpha = 0.3) +
  facet_wrap(~ scenario) +
  theme_bw() +
  labs(
    x = "Replicate",
    y = "Spearman ρ",
    title = "Supplementary Figure 4D. Correlation between true and inferred relatedness"
  )


library(patchwork)

supp_fig4 <- (
  supp4A / supp4B / supp4C / supp4D
) + plot_layout(heights = c(1, 1, 1, 1))


