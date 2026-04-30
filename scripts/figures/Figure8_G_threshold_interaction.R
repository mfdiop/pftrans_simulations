
# ============================================================
# Figure 8: G-threshold × recombination × method interaction
# G threshold × recombination × method interaction
# how surveillance temporal window interacts with biology to define identifiability limits
# AUTHOR: [You]
# DATE:   Sys.Date()
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
G_LABELS   <- paste0("G=", G_VALS)

message("\n▶ Building Figure 8")

agg_G <- df %>%
  group_by(G_threshold, rec_rate, method) %>%
  summarise(
    mean  = mean(auprc, na.rm = TRUE),
    ci_lo = quantile(auprc, 0.025, na.rm = TRUE),
    ci_hi = quantile(auprc, 0.975, na.rm = TRUE),
    .groups = "drop")

# ------------------------------------------------------------
# Panel A: heatmaps
# ------------------------------------------------------------

plot_G_heat <- function(meth) {
  
  agg_G %>%
    filter(method == meth) %>%
    ggplot(aes(factor(rec_rate), factor(G_threshold), fill = mean)) +
    geom_tile(color = "white") +
    # geom_text(aes(label  = sprintf("%.2f", mean), colour = mean > 0.62), size = 4) +
    geom_text(aes(label  = sprintf("%.2f", mean), colour = ifelse(mean > 0.62, "white", "black")), size = 4) +
    scale_colour_identity(guide = "none", breaks = c(TRUE, FALSE), labels = c("white", "black")) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdYlGn"), limits = c(0.20, 0.75)) +
    scale_x_discrete(labels = REC_LABELS, expand = c(0, 0)) + scale_y_discrete(labels = G_LABELS, expand = c(0, 0)) +
    labs(title = meth, x = NULL, y = "G threshold (generations)") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 14, colour = "black", angle = 20),
      axis.text = element_text(size = 14, colour = "black"),
      axis.title = element_text(size = 14, colour = "black", face = "bold"),
      axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
      axis.ticks = element_line(color = 'black', linewidth = .7),
      axis.ticks.length = unit(.22, "cm"),
      plot.title = element_text(face = "bold", colour = METHOD_COLORS[meth], hjust = 0.5))
}

p_heat <- map(METHODS, plot_G_heat)

# ------------------------------------------------------------
# Panel B: G-gain curve
# ------------------------------------------------------------

p_gain <- agg_G %>%
  ggplot(aes(G_threshold, mean,
             colour = method,
             linetype = factor(rec_rate))) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = METHOD_COLORS) +
  scale_x_continuous(breaks = G_VALS) +
  labs(
    title = "G-threshold gain curve",
    x = "G threshold (generations)",
    y = "Mean AUPRC") +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 14, colour = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm"),
    plot.title = element_text(face = "bold", colour = "black", hjust = 0.5))

# ------------------------------------------------------------
# Assemble Figure 8
# ------------------------------------------------------------

fig8 <- wrap_plots(p_heat, nrow = 1) / p_gain +
  plot_annotation(tag_levels = "A") + #  title = "Figure 8. Effect of transmission generation threshold (G)\n", 
  theme(plot.tag = element_text(face = "bold", size = 16))

print(fig8)

ggsave(file.path(OUTPUT_DIR, "Figure8_G_threshold_interaction.pdf"),
       fig8, width = 14, height = 12, dpi = 300)

message("✔ Figure 8 saved")

# ============================================================

# DONE

# ============================================================

message("\n✓ All figures generated successfully")
message("  Output directory: ", normalizePath(OUTPUT_DIR))


