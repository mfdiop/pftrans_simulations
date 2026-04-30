
library(tidyverse)
library(patchwork)

df  <- read_csv("simulations/malaria_transmission_study/evaluation/identifiability_results.csv")

g25 <- df %>%
  filter(G_threshold == 25)

REC_VALS   <- sort(unique(g25$rec_rate))
REC_LOG    <- log10(REC_VALS)
REC_LABELS <- c("1.0e-9", "1.0e-8", "1.0e-7", "1.0e-6")

MIG_VALS   <- sort(unique(g25$migration))
MIG_LABELS <- as.character(MIG_VALS)

SAMP_VALS  <- sort(unique(g25$sample_size))

# Panel A data: best AUPRC across methods
best <- g25 %>%
  group_by(rec_rate, migration, sample_size) %>%
  summarise(best_auprc = max(auprc), .groups = "drop") %>%
  mutate(rec_log = log10(rec_rate))

# Panel B data: IBD vs next-best gap
gap_df <- g25 %>%
  group_by(rec_rate, migration, method) %>%
  summarise(auprc = mean(auprc), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = auprc) %>%
  mutate(
    next_best = pmax(IBS, Phylo, na.rm = TRUE),
    gap = round(IBD - next_best, 3))

# Panel A plotting function (phase diagram)
plot_phase <- function(mig_val) {
  
  sub <- best %>% filter(migration == mig_val)
  
  ggplot(sub, aes(rec_log, factor(sample_size), fill = best_auprc)) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = sprintf("%.2f", best_auprc)),
      size = 4,
      colour = ifelse(sub$best_auprc > 0.78, "white", "black")) +
    geom_contour(
      aes(z = best_auprc),
      breaks = 0.80,
      linetype = "dashed",
      linewidth = 0.8,
      colour = "black") +
    scale_fill_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "RdYlGn")),
      # colours = RColorBrewer::brewer.pal(11, "RdYlGn"),
      limits = c(0, 1),
      name = "Best AUPRC\n(max over methods)") +
    scale_x_continuous(breaks = REC_LOG, labels = REC_LABELS) +
    # scale_y_continuous(breaks = SAMP_VALS) +
    scale_y_discrete(labels = c("100", "200", "400")) +
    labs(
      title = paste("Migration =", mig_val),
      x = "Recombination rates", y = "Sample sizes") +
    theme_minimal(base_size = 14) +
    theme(
      legend.text = element_text(size = 14, colour = "black"),
      legend.title = element_text(size = 16, colour = "black", hjust = .5),
      axis.text.x = element_text(size = 14, angle = 20, hjust = 1, vjust = 2, colour = "black"),
      axis.text.y = element_text(size = 14, colour = "black"),
      axis.title = element_text(size = 16, colour = "black", face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
}

# Build Panel A (three plots)
pA <- (plot_phase(MIG_VALS[1]) |
         plot_phase(MIG_VALS[2]) |
         plot_phase(MIG_VALS[3])) +
  plot_layout(guides = "collect")

print(pA)

# Panel B: IBD advantage heatmap

pB <- ggplot(gap_df,
             aes(factor(rec_rate, levels = REC_VALS),
                 factor(migration, levels = MIG_VALS),
                 fill = gap)) +
  geom_tile(color = "white") +
  geom_text( aes(label = sprintf("%.3f", gap)),
    size = 7, colour = ifelse(gap_df$gap > 0.12, "white", "black")) +
  # scale_fill_gradient(low = "#fff7bc", high = "#7f0000", limits = c(0, 0.20), name = "AUPRC gap") +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlGn")),
                       limits = c(0, 0.20), name = "AUPRC gap") +
  scale_x_discrete(labels = REC_LABELS) +
  scale_y_discrete(labels = MIG_LABELS) +
  labs(
    title = "IBD advantage over next-best method\n(AUPRC gap = IBD − max(IBS, Phylo), G ≤ 25)",
    x = "Recombination rate",
    y = "Migration rate" ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    legend.text = element_text(size = 14, colour = "black"),
    legend.title = element_text(size = 16, colour = "black", hjust = .5),
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))

# Final layout (exact analogue of GridSpec)
# final_fig <- (pA) / (pB) +
#   plot_layout(heights = c(1.4, 1)) +
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(face = "bold", size = 16))
# 
# print(final_fig)

# Save outputs
ggsave("results/main/Figure7_detectability.pdf", pA,
       width = 14, height = 9, dpi = 600)

ggsave("results/supplementary/FigSxxx_detectability.pdf", pB,
       width = 14, height = 9, dpi = 600)




