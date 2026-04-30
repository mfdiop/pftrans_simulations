
library(tidyverse)
library(patchwork)

mig_df <- read_csv("simulations/malaria_transmission_study/evaluation/migration_all_results.csv")

xx <- mig_df %>% 
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

base_theme <- theme(
  axis.title = element_text(size = 15, color = "black", face = "bold"),
  # axis.text.x = element_text(color = "black", size = 14), # angle = 30, hjust = 1, 
  axis.text = element_text(color = "black", size = 14),
  axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
  axis.ticks = element_line(color = 'black', linewidth = .7),
  axis.ticks.length = unit(.25, "cm"))

fig6b <- ggplot(xx, aes(x = migration, y = f1, color = method, group = method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.5) +
  stat_summary(fun = median, geom = "point", size = 5, show.legend = FALSE) +
  # scale_x_continuous(expand = c(0,0), breaks = seq(0, 0.05, 0.01), labels = as.character(seq(0, 0.05, 0.01))) +
  labs(x = "Migration Rate",
       y = "Precision–recall balance (F1 score)" ) +
  theme_classic() + base_theme +
  theme(legend.text = element_text(size = 15, colour = "black", face = "bold"),
        legend.title = element_blank())

print(fig6b)

ggsave("results/main/figure6_migration_F1Scores.png", 
       plot = fig6b, width = 10, height = 6, dpi = 600)
