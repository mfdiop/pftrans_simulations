
# ============================================================ 
# Figure 3: Comparative recoverability across inference layers
# with migration as a continuous gradient across the full factorial
# AUTHOR: [You]
# DATE:   Sys.Date()
# ============================================================

message("▶ Starting identifiability figure generation")

# ============================================================

# 0. Libraries

# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# ============================================================

# 1. Load and prepare data

# ============================================================

INPUT_FILE  <- "simulations/malaria_transmission_study/evaluation/identifiability_results.csv"
OUTPUT_DIR  <- "results/main"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

message("▶ Reading input file: ", INPUT_FILE)

df <- read_csv(INPUT_FILE, show_col_types = FALSE)

message("  Rows read: ", nrow(df))

# Derived metrics

df <- df %>%
  mutate(
    sensitivity = TP / (TP + FN),
    FDR = if_else((TP + FP) > 0,
                  FP / (TP + FP),
                  NA_real_))

# ============================================================

# 2. Global constants

# ============================================================

METHODS <- c("IBD", "IBS", "Phylo")
METHOD_COLORS <- c(
  IBD   = "#E07B39",
  IBS   = "#4878CF",
  Phylo = "#6ACC65"
)

REC_VALS   <- sort(unique(df$rec_rate))
REC_LABELS <- as.character(REC_VALS)

MIG_VALS   <- sort(unique(df$migration))
MIG_LABELS <- as.character(MIG_VALS)

G_VALS     <- sort(unique(df$G_threshold))

# ============================================================
# FIGURE 3
# Comparative recoverability across genomic inference layers
# ============================================================

message("\n▶ Building Figure 3")

g25 <- df %>% filter(G_threshold == 25)

# ------------------------------------------------------------
# Aggregation
# ------------------------------------------------------------

agg_auprc <- g25 %>%
  group_by(rec_rate, migration, method) %>%
  summarise(
    mean  = mean(auprc),
    q25   = quantile(auprc, 0.25),
    q75   = quantile(auprc, 0.75),
    ci_lo = quantile(auprc, 0.025),
    ci_hi = quantile(auprc, 0.975),
    .groups = "drop" )

agg_sens <- g25 %>%
  group_by(rec_rate, migration, method) %>%
  summarise(
    mean  = mean(sensitivity_at_90spec),
    q25   = quantile(sensitivity_at_90spec, 0.25),
    q75   = quantile(sensitivity_at_90spec, 0.75),
    ci_lo = quantile(sensitivity_at_90spec, 0.025),
    ci_hi = quantile(sensitivity_at_90spec, 0.975),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Helper: grouped bar + CI
# ------------------------------------------------------------

plot_grouped_ci <- function(dat, title, ylab, show_xlab = FALSE, ref = NULL) {  
  
  ggplot(dat, aes(factor(migration), mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.75),
             width = 0.65, alpha = 0.80) +
    geom_errorbar(
      aes(ymin = ci_lo, ymax = ci_hi),
      position = position_dodge(width = 0.75),
      width = 0.18, linewidth = 0.6 ) +
    scale_fill_manual(values = METHOD_COLORS) +
    scale_y_continuous(limits = c(0, 0.85), expand = c(0,0)) + 
    labs(title = title, x = if (show_xlab) "Migration rate" else NULL, y = ylab) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text = element_text(size = 14, colour = "black"),
      axis.title = element_text(size = 16, colour = "black", face = "bold"),
      axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
      axis.ticks = element_line(color = 'black', linewidth = .7),
      axis.ticks.length = unit(.22, "cm")) +
    {if (!is.null(ref)) geom_hline(yintercept = ref, linetype = "dashed", colour = "grey30")}
}

# ------------------------------------------------------------
# Row 0 & 1: AUPRC + Sensitivity panels
# ------------------------------------------------------------

p_auprc <- map2(
  REC_VALS, REC_LABELS,
  ~ plot_grouped_ci(
    agg_auprc %>% filter(rec_rate == .x),
    title = paste("r =", .y),
    ylab  = "Mean AUPRC (95% CI)",
    ref   = 0.5))

# ------------------------------------------------------------
# Row 2: Gradient line plot
# ------------------------------------------------------------

p_grad <- agg_auprc %>%
  ggplot(aes(log10(rec_rate), mean,
             colour = method,
             group = interaction(method, migration),
             linetype = factor(migration))) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = method),
              alpha = 0.06, colour = NA) +
  scale_colour_manual(values = METHOD_COLORS) +
  scale_fill_manual(values = METHOD_COLORS) +
  scale_x_continuous(breaks = log10(REC_VALS), labels = REC_LABELS) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey30") +
  labs(title = "AUPRC gradient across recombination and migration (G ≤ 25)",
    x = "Recombination rate", y = "Mean AUPRC (95% CI)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.25, "cm"),
    legend.position = "bottom",
    legend.text = element_text(size = 14, colour = "black"),
    legend.title = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))

# ------------------------------------------------------------
# Assemble Figure 3
# ------------------------------------------------------------

fig3 <- (wrap_plots(p_auprc, nrow = 1) / p_grad) +
  plot_annotation(tag_levels = "A") & # title = "Figure 3. Comparative recoverability of transmission links across\n", 
  theme(plot.tag = element_text(face = "bold", size = 16))

print(fig3)

ggsave(file.path(OUTPUT_DIR, "Figure6_migration_redesigned.pdf"),
       fig3, width = 16, height = 16, dpi = 300)

# Supplementary
p_sens <- map2(
  REC_VALS, REC_LABELS,
  ~ plot_grouped_ci(
    agg_sens %>% filter(rec_rate == .x),
    title = paste("r =", .y),
    ylab  = "Sensitivity @ 90% spec",
    ref   = 0.5,
    show_xlab = TRUE))

figS3 <- wrap_plots(p_sens,  nrow = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("results/supplementary/FigSxxx_migration_sensitivity.pdf",
       figS3, width = 16, height = 16, dpi = 600)


message("✔ Figure 3 saved")
