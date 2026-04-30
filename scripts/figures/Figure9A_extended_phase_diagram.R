
# ══════════════════════════════════════════════════════════════════════════
# Figure generation script for all missing/priority figures
# ══════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(patchwork)
library(scales)
library(reshape2)

options(warn = -1)

# ── load data ──────────────────────────────────────────────────────────────
df <- read_csv("simulations/malaria_transmission_study/evaluation/identifiability_results.csv")

REC_ORDER   <- c(1e-9, 1e-8, 1e-7, 1e-6)
REC_LABELS  <- c("1e-09", "1e-08", "1e-07", "1e-06")
MIG_ORDER   <- c(0.001, 0.010, 0.050)
MIG_LABELS  <- c("0.001", "0.010", "0.050")
METHOD_COLS <- c("IBD"="#E07B39", "IBS"="#4878CF", "Phylo"="#6ACC65")
METHOD_LIST <- c("IBD", "IBS", "Phylo")

OUTPUT <- "results/main"

theme_base <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_line(linetype="dashed"),
    legend.text = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    axis.text = element_text(color = "black", size = 14),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm"),
    strip.text = element_text(size = 16, color = "black"))

# ══════════════════════════════════════════════════════════════════════════
# PRIORITY 1 — Migration × recombination AUPRC gradient
# ══════════════════════════════════════════════════════════════════════════

g25 <- df %>% filter(G_threshold == 25)

# ══════════════════════════════════════════════════════════════════════════
# PRIORITY 2 — Extended phase diagram
# ══════════════════════════════════════════════════════════════════════════

best_auprc <- g25 %>%
  group_by(rec_rate, migration, sample_size) %>%
  summarise(best_auprc=max(auprc), .groups="drop") %>%
  mutate(log_rec=log10(rec_rate))

# Transform the data to fill the gap in sample size
library(akima)  # for interpolation
library(zoo)  # for na.locf

best_auprc_interp <- best_auprc %>%
  tidyr::complete(
    migration, log_rec,
    sample_size = seq(100, 400, 50))   # fill missing sizes
  # group_by(migration) %>%
  # mutate(best_auprc = approx(sample_size, best_auprc, xout = sample_size)$y) %>%
  # ungroup()

# Assuming your dataframe is called df
df_filled <- best_auprc_interp %>%
  arrange(migration, log_rec, sample_size) %>%   # make sure data is sorted
  group_by(migration, log_rec) %>%              # fill within each migration × rec group
  mutate(rec_rate = zoo::na.locf(rec_rate),
         best_auprc = zoo::na.locf(best_auprc)) %>% # carry last observation forward
  ungroup()

# Make plot
p2 <- df_filled %>%
  ggplot(aes(x=log_rec, y=sample_size, fill=best_auprc)) +
  geom_tile() +
  geom_contour(aes(z = best_auprc), breaks = c(0.80), color = "black",
               linetype = "dashed", linewidth = .7) +
  facet_wrap(~migration, nrow=1, labeller = labeller(migration=function(x) paste("Migration =", x))) +
  scale_fill_gradientn(colours=rev(RColorBrewer::brewer.pal(11,"RdYlGn")), limits=c(0.45,1.0)) +
  scale_x_continuous(breaks=log10(REC_ORDER), labels=REC_LABELS, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs( x="Recombination rate (log10)", y="Sample size", fill="Best AUPRC") +
  theme_base

print(p2)

combined <- (p2 / (pC + pD)) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(file.path(OUTPUT,"Figure9_validation.pdf"),
       combined, width = 14, height = 12, dpi = 600)
