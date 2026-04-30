
# sim_design_figure.R
# Produce a 2-panel "Simulation design" figure for manuscript (workflow + scenarios)
# Output: figures/simulation_design.pdf and .png
#
# Usage: source("sim_design_figure.R"); make_sim_design_figure()
#
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(grid)
  library(gridExtra)
  library(data.table)
})

make_sim_design_figure <- function(outdir = "figures",
                                   outfile_base = "simulation_design",
                                   width = 12, height = 7, dpi = 300) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------- Panel A: Workflow schematic ----------
  # We'll place boxes at fixed coordinates and draw arrows between them.
  
  boxes <- data.table(
    id = c("SLiM_sim", "tskit_ts", "tskitIBD", "inferred_metrics", "evaluation"),
    label = c("SLiM\nforward simulation\n(demography, selection)",
              "Tree sequence\n(tskit)\n Add mutation (msprime)\nGenerate VCF data", 
              "IBD truth\n(tskitIBD)\n(segment start/end, TMRCA)",
              "Inferred metrics\n .IBS\n .IBD \nPhylogenetic Patristic distances",
              "Evaluation Metrics\nAUPR, AUROC, precision@K,\nSpearman, confusion matrices"),
    x = c(0.13, 0.38, 0.63, 0.63, 0.38),
    y = c(0.85, 0.85, 0.85, 0.50, 0.50),
    w = 0.22, h = 0.18
  )
  
  # arrows: from -> to as indices in boxes$id
  arrows <- data.table(
    from = c("SLiM_sim", "tskit_ts", "tskitIBD", "inferred_metrics"),
    to   = c("tskit_ts", "tskitIBD", "inferred_metrics", "evaluation")
  )
  
  pA <- ggplot() + theme_void() + xlim(0,1) + ylim(0,1)
  
  # draw boxes
  for (i in seq_len(nrow(boxes))) {
    b <- boxes[i]
    pA <- pA +
      annotate("rect", xmin = b$x - b$w/2, xmax = b$x + b$w/2,
               ymin = b$y - b$h/2, ymax = b$y + b$h/2,
               fill = "#f7fbee", color = "#2b8cbe", 
               linewidth = 1, alpha = 0.75) +
      annotate("text", x = b$x, y = b$y, label = b$label, size = 4.5, lineheight = 0.95)
  }
  
  print(pA)
  
  # # draw arrows between box centers
  # for (i in seq_len(nrow(arrows))) {
  #   from_id <- arrows[i, from]
  #   to_id   <- arrows[i, to]
  #   f <- boxes[id == from_id]
  #   t <- boxes[id == to_id]
  # 
  #   # compute start and end just outside rectangles
  #   start_x <- f$x + ifelse(t$x > f$x, f$w/2, -f$w/2)
  #   start_y <- f$y
  #   end_x <- t$x - ifelse(t$x > f$x, t$w/2, -t$w/2)
  #   end_y <- t$y
  #   pA <- pA + geom_segment(x = start_x, y = start_y, xend = end_x, yend = end_y,
  #                           arrow = arrow(length = unit(0.02, "npc")),
  #                           size = 0.7, colour = "#4d4d4d")
  # }
  # 
  # arrow_df <- do.call(rbind, lapply(seq_len(nrow(arrows)), function(i) {
  #   from_id <- arrows[i, from]
  #   to_id   <- arrows[i, to]
  #   
  #   f <- boxes[id == from_id]
  #   t <- boxes[id == to_id]
  #   
  #   data.frame(
  #     start_x = f$x + ifelse(t$x > f$x, f$w/2, -f$w/2),
  #     start_y = f$y,
  #     end_x   = t$x - ifelse(t$x > f$x, t$w/2, -t$w/2),
  #     end_y   = t$y
  #   )
  # }))
  
  arrow_df <- do.call(rbind, lapply(seq_len(nrow(arrows)), function(i) {
    from_id <- arrows[i, from]
    to_id   <- arrows[i, to]
    
    f <- boxes[id == from_id]
    t <- boxes[id == to_id]
    
    dx <- t$x - f$x
    dy <- t$y - f$y
    
    # choose dominant direction
    if (abs(dx) > abs(dy)) {
      # horizontal arrow
      start_x <- f$x + ifelse(dx > 0, f$w/2, -f$w/2)
      start_y <- f$y
      end_x   <- t$x - ifelse(dx > 0, t$w/2, -t$w/2)
      end_y   <- t$y
      
    } else {
      # vertical arrow
      start_x <- f$x
      start_y <- f$y + ifelse(dy > 0, f$h/2, -f$h/2)
      end_x   <- t$x
      end_y   <- t$y - ifelse(dy > 0, t$h/2, -t$h/2)
    }
    
    data.frame(start_x, start_y, end_x, end_y)
  }))
  
  
  pA <- pA +
    geom_segment(
      data = arrow_df,
      aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
      arrow = arrow(length = unit(0.02, "npc")),
      size = 0.7,
      colour = "#4d4d4d"
    )
  
  
  pA <- pA + ggtitle("A) Simulation & analysis workflow") +
    theme(plot.title = element_text(hjust = 0.05, vjust = -5, size = 18, face = "bold"))
  
  print(pA)
  
  ggsave("simulations/main/figure1_simulation_framework.png", plot = pA,
         width = width, height = height, dpi = dpi)
  
  # ---------- Panel B: Scenario panel ----------
  # We'll show 3 scenario boxes with bullet lists and a small param table.
  
  # scenarios text (from your specs)
  scen_text <- list(
    scenario1 = c(
      "Single panmictic population",
      "Fixed parameters",
      "Replicates: 5",
      "Samples per replicate: 1000"
    ),
    scenario2 = c(
      "Single population, varied recombination",
      "Recombination rates: 1e-9, 1e-8, 1e-7, 1e-6",
      "Replicates: 20 (per rate)",
      "Samples per replicate: 1000"
    ),
    scenario3 = c(
      "Multiple populations (epidemiological factors)",
      "Varied: migration, µ, r, Ne, sample size",
      "Bottlenecks: 2500, 3500, 5000",
      "Selection: single-locus selection",
      "Replicates: assorted (scenarios)"
    )
  )
  
  # make a data.frame for plotting scenario boxes
  scen_df <- data.table(
    id = c("S1","S2","S3"),
    title = c("Scenario 1", "Scenario 2", "Scenario 3"),
    x = c(0.2, 0.6, 0.9),
    y = c(0.75, 0.75, 0.75),
    w = c(0.32, 0.36, 0.36),
    h = c(0.38, 0.52, 0.62)
  )
  
  # parameter table as grob
  param_dt <- data.table(
    parameter = c("Population size (Ne)", "Generations simulated", "Mutation rate (µ)", "Genome length",
                  "Sampling per replicate", "Recombination rates", "Replicates (per scenario)"),
    value     = c("50,000", "200", "1e-8", "Chromosome 1 length", "1000", "1e-9,1e-8,1e-7,1e-6", "5 (S1), 20 (S2), variable (S3)")
  )
  
  # # build scenario panels with ggplot
  # pB_base <- ggplot() + xlim(0,1) + ylim(0,1) + theme_void()
  # 
  # # draw scenario boxes with bullet text
  # pB <- pB_base
  # 
  # # draw scenario rectangles
  # for (i in 1:nrow(scen_df)) {
  #   r <- scen_df[i]
  #   # rectangle
  #   pB <- pB + annotate("rect", xmin = r$x - r$w/2, xmax = r$x + r$w/2,
  #                       ymin = r$y - r$h/2, ymax = r$y + r$h/2,
  #                       fill = "#fff7fb", color = "#7b3294", size = 0.6)
  #   # title
  #   pB <- pB + annotate("text", x = r$x, y = r$y + r$h/2 - 0.04, 
  #                       label = r$title, size = 4.2, fontface = "bold")
  #   # bullets
  #   bullets <- scen_text[[paste0("scenario", i)]]
  #   
  #   # compute lines y positions
  #   nlines <- length(bullets)
  #   ys <- seq(r$y + 0.08, r$y - r$h/2 + 0.08, length.out = nlines)
  #   
  #   for (j in seq_along(bullets)) {
  #     pB <- pB + annotate("text", x = r$x - r$w/2 + 0.03, y = ys[j],
  #                         label = paste0("\u2022  ", bullets[j]), 
  #                         hjust = 0, size = 3.2, family = "sans")
  #   }
  # }
  # 
  # print(pB)
  
  # Build rectangle data for geom_rect (xmin/xmax/ymin/ymax)
  rect_df <- copy(scen_df)
  rect_df[, `:=`(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2)]
  
  # Build title data for geom_text
  title_df <- rect_df[, .(x = x, y = ymax - 0.04, label = title)]
  
  # Build bullets data: compute y positions for each bullet inside each box
  bullet_rows <- list()
  for (i in seq_len(nrow(rect_df))) {
    box <- rect_df[i]
    bullets <- scen_text[[paste0("scenario", i)]]
    n <- length(bullets)
    # top of text area just below the title; bottom leave margin
    top <- box$y + box$h/2 - 0.08
    bottom <- box$y - box$h/2 + 0.06
    ys <- seq(top, bottom, length.out = n)
    for (j in seq_along(bullets)) {
      bullet_rows[[length(bullet_rows)+1]] <- data.table(
        x = box$x - box$w/2 + 0.03,
        y = ys[j],
        label = bullets[j],
        box_id = box$id
      )
    }
  }
  bullet_df <- rbindlist(bullet_rows)
  
  # Now build the ggplot panel
  pB <- ggplot() +
    # draw rects
    geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#fff7fb", color = "#7b3294", size = 0.6) +
    # titles
    geom_text(data = title_df, aes(x = x, y = y, label = label), fontface = "bold", size = 4.2) +
    # bullets (use \u2022 and left-justified)
    geom_text(data = bullet_df, aes(x = x, y = y, label = paste0("\u2022  ", label)),
              hjust = 0, size = 3.2, family = "sans") +
    xlim(0,1) + ylim(0,1) +
    theme_void() +
    ggtitle("B) Simulation scenarios and core parameters") +
    theme(plot.title = element_text(hjust = 0.02, vjust = -4, size = 14, face = "bold"))
  
  # show
  print(pB)
  
  # parameter table as a ggplot table (right-lower)
  # convert to table grob
  tbl <- tableGrob(param_dt, rows = NULL, theme = ttheme_minimal(core=list(fg_params = list(cex=0.9)),
                                                                 colhead=list(fg_params=list(cex=1))))
  # place table grob as annotation on pB using annotation_custom later
  
  pB <- pB + ggtitle("B) Simulation scenarios and core parameters") +
    theme(plot.title = element_text(hjust = 0.02, vjust = -4, size = 14, face = "bold"))
  
  # assemble panels: left (pA) and right (pB + table)
  # convert table grob to patchwork element using cowplot::ggdraw
  table_plot <- ggdraw() + draw_grob(tbl)
  # arrange B: put scenario panel above and table below in right column
  right_col <- plot_grid(pB, table_plot, ncol = 1, rel_heights = c(2.3, 1))
  
  # full figure: put A and right_col side by side
  full <- plot_grid(pA, right_col, ncol = 2, rel_widths = c(1.1, 1))
  
  # save outputs
  pdf_file <- file.path(outdir, paste0(outfile_base, ".pdf"))
  png_file <- file.path(outdir, paste0(outfile_base, ".png"))
  
  ggsave(pdf_file, full, width = width, height = height)
  ggsave(png_file, full, width = width, height = height, dpi = dpi)
  
  message("Saved: ", pdf_file)
  message("Saved: ", png_file)
  invisible(list(gg = full, pdf = pdf_file, png = png_file))
}

# Run it
make_sim_design_figure()


# ========================================================================
# ========================================================================


###############################################
##   FIGURE 1 — Pipeline + Scenario Diagram  ##
###############################################

## Required packages
library(tidyverse)
library(ggforce)
library(ggtext)
library(patchwork)


plot_boxes <- function(df) {
  ggplot() +
    geom_rect(
      data = df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#2b8cbe", color = "black", linewidth = 0.6) +
    geom_text(
      data = df,
      aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2, label = label),
      size = 4, vjust = 0.5, hjust = 0.5) +
    theme_void()
}

plot_scenarios <- function(df) {
  ggplot() +
    geom_rect(
      data = df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#2b8cbe", color = "black", linewidth = 0.6) + 
    geom_text(
      data = df,
      aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2, label = label),
      # fill = NA, label.color = NA,
      size = 4, vjust = 0.5, hjust = 0.5, fontface = "bold"
    ) +
    theme_void() +
    coord_fixed()
}


###############################################
## 1. Define positions for scenario boxes     ##
###############################################

scenario_boxes <- tribble(
  ~name, ~label, ~xmin, ~xmax, ~ymin, ~ymax,
  
  "baseline", 
  "Scenario1\nBaseline Biological Parameters
   μ = 1×10⁻⁸
   Ne: 10,000 → 1,000
   Selfing = 95%
   Drug pressure (s = 0.3)",
  -0.1, 0.35, 0.75, 0.98,
  
  "recombination", 
  "Scenario2\nRecombination Sweep
   r = 1×10⁻⁹ → 1×10⁻⁶",
  -0.1, 0.35, 0.42, 0.66,
  
  "multipop",
  "Scenario3\nMulti-Population Structure
   2–4 subpops
   Migration = 0.001–0.05
   + Full Factorial Exploration",
  -0.1, 0.35, 0.10, 0.38
)

###############################################
## 2. Main pipeline boxes                     ##
###############################################

pipeline_boxes <- tribble(
  ~name, ~label, ~xmin, ~xmax, ~ymin, ~ymax,
  
  "slim", 
  "SLiM Simulation\nPedigree + Tree Sequence",
  0.7, 1.35, 0.8, 0.85,
  
  "tskit", 
  "Tree-Seq Processing\ntskit / pyslim",
  1.7, 2.3, 0.8, 0.85,
  
  "mutation", 
  "Recapitate - Mutation Overlay\nmsprime",
  2.7, 3.3, 0.8, 0.85,
  
  "genome", 
  "Genome Output\n VCF",
  3.7, 4.3, 0.8, 0.85
)

###############################################
## 3. Inference layer                         ##
###############################################

inference_boxes <- tribble(
  ~name, ~label, ~xmin, ~xmax, ~ymin, ~ymax,
  
  "true_ibd", 
  "True IBD\n(Pedigree-derived)
  tskitIBD",
  3.6, 5.2, 0.05, 0.25,
  
  "hmmibd",
  "Inferred IBD\n(hmmIBD)",
  5.7, 6.7, 0.05, 0.25, # 4.0, 4.6, 0.3, 0.5,
  
  "ibs",
  "IBS Distances\n",
  7.2, 8.2, 0.05, 0.25, #0.45, 0.65,
  
  "phylo",
  "ML Phylogeny\n",
  8.8, 9.8, 0.05, 0.25 #0.7, 0.9
)

###############################################
## 4. Evaluation box                          ##
###############################################

evaluation_box <- tribble(
  ~name, ~label, ~xmin, ~xmax, ~ymin, ~ymax,
  "eval",
  "Evaluation Metrics
  AUC-PR
  F1, 
  ROC",
  5, 6.3, 0.45, 0.85
)

###############################################
## 5. Define arrows                           ##
###############################################

## Scenario → SLiM arrows
scenario_arrows <- tribble(
  ~x, ~y, ~xend, ~yend,
  0.35, 0.86, 0.7, 0.65,   # baseline
  0.35, 0.56, 0.7, 0.65,   # recombination
  0.35, 0.26, 0.7, 0.65    # multipop/full factorial
)

## Pipeline arrows
pipeline_arrows <- tribble(
  ~x, ~y, ~xend, ~yend,
  1.3, 0.65, 1.7, 0.65,
  2.3, 0.65, 2.7, 0.65,
  3.3, 0.65, 3.7, 0.65,
  
  ## genome → inference branches
  4.3, 0.65, 3.9, 0.20,   # true IBD
  4.3, 0.65, 4.3, 0.40,   # hmmIBD
  4.3, 0.65, 4.7, 0.55,   # IBS
  4.3, 0.65, 5.1, 0.80,   # phylo
  
  ## inference → evaluation
  3.9, 0.20, 5.7, 0.65,
  4.3, 0.40, 5.7, 0.65,
  4.7, 0.55, 5.7, 0.65,
  5.1, 0.80, 5.7, 0.65
)

###############################################
## 6. Plot                                     ##
###############################################

p <- plot_scenarios(scenario_boxes) +
  plot_boxes(pipeline_boxes) +
  plot_boxes(inference_boxes) +
  plot_boxes(evaluation_box) +
  
  geom_link(
    data = scenario_arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey30",
    linewidth = 0.8,
    n = 200) +
  
  geom_link(
    data = pipeline_arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey30",
    linewidth = 0.8,
    n = 200) +
  theme_void() +
  coord_fixed()

print(p)

###############################################
## 7. Save                                     ##
###############################################

ggsave("figure1_pipeline.pdf", p, width = 13, height = 6)
ggsave("figure1_pipeline.svg", p, width = 13, height = 6)
