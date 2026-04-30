
library(tidyverse)
library(patchwork)

# We need this following information
# scenario              # baseline, migration, factorial
# run_id                # replicate
# method                # hmmIBD, true_phasing, IBS, phylogeny
# TP, FP, FN            # transmission link classification
# precision             # TP/(TP+FP)
# recall                # TP/(TP+FN)
# F1                    # harmonic mean
# threshold             # IBD or IBS cutoff used (for PR curves)

# For the single run simulation, this was already computed, so we load this list from the 
# code: plot_aggregated_confusion_matrix.R
# all_results <- readRDS("simulations/single_run/evaluation_results.rds")
# all_confusion <- plot_all_confusion_matrices(all_results, "optimal_youden")

# ======================= FUNCTIONS ===============================
base_theme <- theme(
  axis.title = element_text(size = 14, color = "black", face = "bold"),
  # axis.text.x = element_text(color = "black", size = 14), # angle = 30, hjust = 1, 
  axis.text = element_text(color = "black", size = 14),
  axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
  axis.ticks = element_line(color = 'black', linewidth = .7),
  axis.ticks.length = unit(.22, "cm"))

sc.colors <- c("bisque4", "darkolivegreen4", "aquamarine4" )

# Transform list to a dataframe across replicates for each metrics
extract_metrics <- function(rep_obj, scenario, run_id) {
  rep_obj$evaluations$optimal_youden$metrics_summary %>%
    as_tibble() %>%
    transmute(
      scenario = scenario,
      run_id   = run_id,
      method = Metric,
      TP = TP,
      FP = FP,
      TN = TN,
      FN = FN,
      precision = Precision,
      recall    = Recall,
      f1        = F1_Score
    )
}

extract_pr_curves <- function(rep_obj, scenario, run_id) {
  
  methods <- names(rep_obj$curve_data)
  
  map_dfr(methods, function(m) {
    
    curve_df <- rep_obj$curve_data[[m]]$pr$curve %>%
      as_tibble()
    
    names(curve_df) <- c("recall", "precision", "threshold")
    
    curve_df %>%
      transmute(
        scenario  = scenario_name,
        run_id    = run_id,
        method    = m,
        recall    = recall,
        precision = precision,
        threshold = if ("threshold" %in% names(curve_df))
          threshold
        else
          NA_real_
      )
  })
}

# ========================================================

single_results <- readRDS("simulations/single_run/evaluation_results.rds")

scenario_name <- "baseline"

baseline_summary <- imap_dfr(
  single_results,
  ~ extract_metrics(
    rep_obj = .x,
    scenario = scenario_name,
    run_id = .y
  )
)

baseline_summary <- baseline_summary %>% 
  add_column(rate = 6.667e-7) %>% 
  relocate(rate, , .after = run_id)
  
# ==========================
# Multiple runs simulation
# ==========================
# I ran the code evaluate_recombination_effects_v1.0.R and got these metrics
# performance_dt <- rbindlist(all_performance, fill = TRUE)
# true transmissions; must have (i, j, true_link)
performance_dt <- readxl::read_xlsx("simulations/multiple_runs/recombination_evaluation/tables/evaluation_metrics.xlsx")

performance_dt <- performance_dt %>% 
  # select(-rate_num) %>% 
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
  select(c(2, 4, 15, 22:25, 30)) %>%   # -ends_with("label"), -replicate_id, -sampling_prop, -note, -n_pairs, -n_positive
  relocate(scenario, .before = run_id) %>%
  relocate(method, .after = run_id) %>%
  rename(rate = rec_rate) %>% 
  mutate(
    precision = if_else(
      TP + FP == 0,
      NA_real_,
      TP / (TP + FP)
    ),
    recall = if_else(
      TP + FN == 0,
      NA_real_,
      TP / (TP + FN)
    ),
    f1 = if_else(
      is.na(precision) | is.na(recall) | (precision + recall) == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)
    )
  )


# ================= COMBINE DATA FRAMES =====================
df_perf <- rbind.data.frame(baseline_summary, performance_dt, migration) %>%
  mutate(method = recode(method, "phylo" = "Phylo"),
         scenario = factor(scenario, levels = c("baseline", "recombination_sweep", "migration")))

# ============================
# 1. Panel A: Sensitivity
# ============================

pA <- df_perf %>% 
  group_by(scenario, method) %>%
  summarise(
    sensitivity = median(recall, na.rm = TRUE),
    lower = quantile(recall, 0.25, na.rm = TRUE),
    upper = quantile(recall, 0.75, na.rm = TRUE),
    .groups = "drop") %>%
  ggplot(aes(method, sensitivity, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,
                position = position_dodge(width = 0.7)) +
  scale_fill_manual(labels = c(recombination_sweep = "recombination\nsweep"),
                    values = c("#30859b", "#d7d7d7", "#f79645")) +
  labs(x = "", # Inference method
    y = "Sensitivity (TP / (TP + FN))",
    fill = "Scenario") +
  theme_classic() + base_theme +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = "black", face = "bold"))


# ============================
# 2. Panel B: False Discovery Rate
# ============================

pB <- df_perf %>%
  group_by(scenario, method) %>%
  summarise(
    FDR = median(FP / pmax(TP + FP, 1), na.rm = TRUE),
    lower = quantile(FP / pmax(TP + FP, 1), 0.25, na.rm = TRUE),
    upper = quantile(FP / pmax(TP + FP, 1), 0.75, na.rm = TRUE),
    .groups = "drop") %>%
  ggplot(aes(method, FDR, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2,
                position = position_dodge(width = 0.7))  +
  scale_fill_discrete(labels = c(recombination_sweep = "recombination\nsweep"),
                      palette = c("#30859b", "#d7d7d7", "#f79645")) +
  labs( x = "", y = "False Discovery Rate (FP / (TP + FP))", fill = "Scenario") +
  theme_classic() + base_theme +
  theme(
    # legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 15, colour = "black", face = "bold")
    )


# ===========================================================
# ======= COMBINE PRECISION-RECALL DATAFRAMES ===============
# ===========================================================

pr_baseline <- imap_dfr(
  single_results,
  ~ extract_pr_curves(
    rep_obj = .x,
    scenario = scenario_name,
    run_id  = .y
  )
)

prcurve_dt <- read_tsv("simulations/multiple_runs/recombination_evaluation/tables/pr_curve.tsv")

prcurve_dt <- prcurve_dt %>% 
  select(-rate) %>%
  add_column(scenario = "recombination_sweep") %>% 
  relocate(scenario, , .before = replicate) %>% 
  rename(run_id = replicate)

colnames(prcurve_dt)[4:6] <- c("recall", "precision", "threshold")

migration_dt <- read_tsv("simulations/malaria_transmission_study/migration_prcurves.tsv")

migration_dt <- migration_dt %>% 
  select(c(run_id, method, recall, precision, threshold)) %>%
  add_column(scenario = "migration") %>% 
  relocate(scenario, , .before = run_id)

df_curve <- pr_baseline %>% 
  # select(-run_id) %>%
  rbind.data.frame(., prcurve_dt) %>%
  rbind.data.frame(., migration_dt) %>%
  mutate(method = recode(method, "phylo" = "Phylo"))

# ======================================
# 3. Panel C: Precision–Recall Curves
# ======================================
METHOD_COLS <- c("IBD"="#E07B39", "IBS"="#4878CF", "Phylo"="#6ACC65")

df_curve_binned <- df_curve %>%
  mutate(recall_bin = round(recall, 3)) %>%   # controls smoothness
  group_by(scenario, method, recall_bin) %>%
  summarise(
    precision = mean(precision),
    .groups = "drop")


pC <- df_curve_binned %>% 
  mutate(scenario = factor(scenario, levels = c("baseline", "recombination_sweep", "migration"))) %>%
  ggplot(aes(recall_bin, precision, color = method)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = METHOD_COLS) +
  facet_wrap(
    ~ scenario,
    labeller = as_labeller(c(
      baseline = "Baseline", migration = "Migration",
      recombination_sweep = "Recombination\nsweep"))) +
  labs(x = "Recall", y = "Precision", color = "Method") +
  coord_equal() + theme_bw() + base_theme +
  theme(legend.position = "none",
    legend.text = element_text(size = 15, color = "black", face = "bold"),
    legend.title = element_blank(),
    strip.text = element_text(size = 15, color = "black", face = "bold"))

# # ===========================================
# # 4. Panel D: Method Ranking Heatmap
# # ===========================================
# 
# df_rank <- df_perf %>%
#   mutate(method = recode(method, "phylo" = "Phylo")) %>% 
#   group_by(scenario, method) %>%
#   summarise(F1_med = median(f1, na.rm = TRUE),
#             .groups = "drop") %>%
#   group_by(scenario) %>%
#   mutate(rank = dense_rank(desc(F1_med)))
# 
# pD <- df_rank %>%
#   mutate(scenario = factor(scenario, levels = c("baseline", "recombination_sweep", "migration")),
#   rank = factor(rank)) %>% 
#   ggplot(aes(method, scenario, fill = rank)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = rank), size = 8, fontface = 'bold', color = "white") +
#   scale_fill_manual( values = c("1" = "#1a9850", "2" = "orange", "3" = "darkred"), name = "Rank\n(1 = Best)") +
#   # scale_fill_viridis_c(direction = -1) +
#   scale_y_discrete(labels = c(baseline = "Baseline", migration = "Migration",
#     recombination_sweep = "Recombination\nsweep")) +
#   labs(x = "", y = "", fill = "Rank\n(1 = Best)") +
#   theme_minimal() + 
#   theme(
#     axis.text        = element_text(size = 16, color = "black", face = "bold"),
#     axis.title.x     = element_text(size = 16, color = "black", face = "bold", margin = margin(t = 8)),
#     plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
#     legend.text = element_text(size = 15, color = "black", face = "bold"),
#     legend.title = element_text(size = 15, color = "black", face = "bold", hjust = 0.5),
#     panel.grid       = element_blank()
#   )

# ============================
# Assemble Figure 3
# ============================

Figure3 <- pA + pB + #plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 16))

Figure3 <- ggguides::collect_legends(Figure3, position = "bottom")

print(Figure3)

ggsave("results/supplementary/figureSx_method_comparison.pdf", plot = Figure3, 
       width = 18, height = 10, dpi = 600) 


# ===========================================================================
#                   Precision-Recall Balance (F1 Scores)
# ===========================================================================

df_perf_adj <- df_perf %>%
  group_by(method) %>%
  mutate(
    f1 = ifelse(rate == 6.667e-07,
                median(f1[rate == 6.667e-07], na.rm = TRUE),
                f1) ) %>% ungroup()

anchor_points <- df_perf %>%
  filter(rate == 6.667e-07) %>%
  group_by(method) %>%
  summarise(
    rate = 6.667e-07,
    f1 = median(f1))

fig4D <- df_perf_adj %>%
  ggplot(aes(x = rate, y = f1, color = method,
             linetype = scenario, group = interaction(method, scenario))) +
  
  stat_summary(fun = median, geom = "point", size = 3, show.legend = FALSE) +
  
  stat_summary(fun = median, geom = "line", linewidth = 1.5) +
  
  # geom_segment(data = anchor_points, aes(x = rate, xend = rate, y = 0, yend = f1, color = method),
  #              linetype = "solid", inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_vline(xintercept = 6.667e-07, linetype = "dotted", linewidth = 1) +
  
  scale_x_log10(
    breaks = c(1e-09, 1e-08, 1e-07, 6.67e-07, 1e-06),
    labels = scales::label_scientific()) +
  
  scale_y_continuous( breaks = seq(0, 0.6, 0.1), limits = c(0, 0.6)) +
  
  scale_color_manual(values = METHOD_COLS) +
  
  guides(linetype = guide_legend( override.aes = list(color = "black"))) +
  
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +  # adjust as needed https://ggplot2.tidyverse.org/reference/scale_linetype.html
  
  labs( x = "Recombination rate", y = "Precision–recall balance (F1 score)") +
  
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(color = "black", size = 14),
    axis.line = element_line(linewidth = 1, colour = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm") )

Figure4 <- pC / fig4D + 
  # plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 16))

print(Figure4)

ggsave("results/main/figure4_f1_recombinations.pdf", 
       plot = Figure4, width = 12, height = 10, dpi = 600)

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
  summarise(N = sum(predicted_links)) %>% arrange(N)
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









# ===============================================
# Panel A–C: Precision–Recall curves per scenari
# ===============================================

pC <- df_curve %>%
  ggplot(aes(recall, precision, color = method)) +
  geom_path(alpha = 0.7, linewidth = 1) +
  facet_wrap(~ scenario, labeller = as_labeller(c(baseline = "Baseline", recombination_sweep = "Recombination\nsweep"))) +
  labs(x = "Recall", y = "Precision", color = "Method") +
  coord_equal() + theme_bw() + base_theme +
  theme(
    legend.text = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_text(size = 14, color = "black", face = "bold", hjust = 0),
    strip.text = element_text(size = 15, color = "black", face = 'bold'))

ggsave("simulations/main/figureS3_prcurve_scenarios.pdf", plot = pC, 
       width = 12, height = 10, dpi = 600)

p_f1 <- df_perf %>%
  mutate(method = recode(method, "phylo" = "Phylo")) %>%
  ggplot(aes(x = scenario, y = f1, color = method, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, position = position_dodge(width = 0.75)) +
  geom_jitter(alpha = 0.6, size = 1.6, position = position_dodge(width = 0.75), show.legend = FALSE) +
  scale_x_discrete(labels = c(baseline = "Baseline", mutation = "Mutation",
                              recombination_sweep = "Recombination\nsweep")) +
  theme_bw(base_size = 13) + base_theme +
  labs(x = "", y = "F1-score") +
  theme(
    legend.text = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_blank())

ggsave("simulations/main/figureS4_transmission_detection_accuracy.pdf", 
       plot = p_f1, width = 12, height = 10, dpi = 600) 

p_fpr <- df_perf %>%
  mutate(fp_rate = FP / (TP + FP),
         method = recode(method, "phylo" = "Phylo")) %>%
  ggplot(aes(x = scenario, y = fp_rate, fill = method)) +
  geom_boxplot(alpha = 0.6) +
  scale_x_discrete(labels = c(baseline = "Baseline", recombination_sweep = "Recombination\nsweep")) +
  theme_bw(base_size = 13) + base_theme +
  labs(x = "", y = "False Positive Link Rate") +
  theme(
    legend.text = element_text(size = 13, color = "black", face = "bold"),
    legend.title = element_blank()
    )

ggsave("simulations/main/figureS5_error_rate.pdf", 
       plot = p_fpr, width = 12, height = 10, dpi = 600) 

# Combine plots
figure3 <- 
  p_pr /                # top: PR curves
  (p_f1 | p_fpr) +      # bottom: F1 and FP-rate side by side
  plot_layout(heights = c(2.2, 1.5)) +
  plot_annotation(title = "Figure 3. Performance of genomic methods for identifying transmission links", 
                  theme = theme(plot.title = element_text(size = 18)),
                  tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14))

figure3
