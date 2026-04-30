
# Real-Data Validation Panel 
# =========================== 
# Gambian P. falciparum WGS data | n = 160 samples | 2014–2015 
# Epidemiological proxy : same-household, same-village pairs 
# Methods evaluated : IBD (hmmIBD) · IBS (allele sharing) · Phylo (IQ-TREE patristic) 
# Pipeline 
# -------- 
# 1. Parse all data files and cross-reference sample IDs 
# 2. Compute phylogenetic patristic (tip-to-tip) distances — pure Python, no biopython 
# 3. Build pair-level feature + label dataframe 
# 4. Compute AUPRC with 2 000-replicate bootstrap 95 % CIs 
# 5. Produce 4-panel publication figure + summary CSV


library(tidyverse)
library(readxl)
library(ape)
library(PRROC)
library(boot)
library(patchwork)
library(scales)

# ═══════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ═══════════════════════════════════════════════════════════════════
meta <- read_excel("results/tables/GamMetadata_Final_imputemissingdate.xlsx")

ibd_raw <- read_tsv("results/tables/ibd_hmm.tsv")   # p1 p2 hmm
ibs_raw <- read_tsv("results/tables/ibs.tsv")       # p1 p2 ibs
tree <- read.tree("results/tables/iqtree_boots10k.contree")

geno_ids <- unique(c(ibd_raw$p1, ibd_raw$p2))

meta_geno <- meta %>%
  filter(SampleID %in% geno_ids)

cat("Samples with genomic data:", nrow(meta_geno), "\n")
cat("COIL=1:", sum(meta_geno$COIL == 1), "\n")
cat("COIL>=2:", sum(meta_geno$COIL >= 2), "\n")
cat("Total pairs:", nrow(ibd_raw), "\n")

meta_lookup <- meta_geno %>%
  select(SampleID, VillageCode, HHCode, CompoundCode, COIL, VisitDate)

# ═══════════════════════════════════════════════════════════════════
# 2. COMPUTE PATRISTIC DISTANCES (ape)
# ═══════════════════════════════════════════════════════════════════

patristic_matrix <- cophenetic.phylo(tree)

# Convert tip labels like "SPT35443_L_1" → SampleID only
clean_id <- function(x) str_split(x, "_", simplify = TRUE)[,1]

tip_map <- tibble(
  tip = rownames(patristic_matrix),
  SampleID = clean_id(tip)
)

patristic_df <- as.data.frame(as.table(patristic_matrix)) %>%
  rename(tip1 = Var1, tip2 = Var2, dist = Freq) %>%
  left_join(tip_map, by = c("tip1" = "tip")) %>%
  rename(s1 = SampleID) %>%
  left_join(tip_map, by = c("tip2" = "tip")) %>%
  rename(s2 = SampleID) %>%
  filter(s1 < s2) %>%
  select(s1, s2, dist)

# ═══════════════════════════════════════════════════════════════════
# 3. BUILD PAIR-LEVEL DATAFRAME
# ═══════════════════════════════════════════════════════════════════

make_key <- function(a,b) map2_chr(a,b, ~ paste(sort(c(.x,.y)), collapse="_"))

ibd_raw <- ibd_raw %>%
  mutate(key = make_key(p1,p2))

ibs_raw <- ibs_raw %>%
  mutate(key = make_key(p1,p2))

maxd <- max(patristic_df$dist, na.rm = TRUE)

patristic_df <- patristic_df %>%
  mutate(key = make_key(s1,s2),
         phylo_neg = 1 - (dist / maxd))


merged <- ibd_raw %>%
  left_join(ibs_raw %>% select(key, ibs), by="key") %>%
  left_join(patristic_df %>% select(key, phylo_neg), by="key") %>%
  separate(key, into=c("s1","s2"), sep="_") %>%
  left_join(meta_lookup, by=c("s1"="SampleID")) %>%
  rename_with(~paste0(.,"_1"), c(VillageCode,HHCode,CompoundCode,COIL)) %>%
  left_join(meta_lookup, by=c("s2"="SampleID")) %>%
  rename_with(~paste0(.,"_2"), c(VillageCode,HHCode,CompoundCode,COIL)) %>%
  mutate(same_hh = as.integer(VillageCode_1==VillageCode_2 & HHCode_1==HHCode_2),
         both_mono = as.integer(COIL_1==1 & COIL_2==1) ) %>%
  drop_na(hmm, ibs, phylo_neg)

df_mono <- merged %>% filter(both_mono==1)

prev_all  <- round(mean(merged$same_hh), 3)
prev_mono <- round(mean(df_mono$same_hh), 3)

# ═══════════════════════════════════════════════════════════════════
# 4. AUPRC + BOOTSTRAP CI
# ═══════════════════════════════════════════════════════════════════
auprc_boot <- function(data, indices){
  d <- data[indices,]
  pr <- pr.curve(scores.class0 = d$score[d$label==1],
                 scores.class1 = d$score[d$label==0],
                 curve = FALSE)
  return(pr$auc.integral)
}

compute_auprc <- function(df, score_col){
  tmp <- tibble(
    label = df$same_hh,
    score = df[[score_col]]
  )
  pr <- pr.curve(scores.class0 = tmp$score[tmp$label==1],
                 scores.class1 = tmp$score[tmp$label==0],
                 curve = TRUE)
  
  b <- boot(tmp, auprc_boot, R=2000)
  ci <- boot.ci(b, type="perc")$percent[4:5]
  
  list(
    curve = pr$curve,
    auc   = pr$auc.integral,
    lo    = ci[1],
    hi    = ci[2]
  )
}

# Compute Results
methods <- c("hmm","ibs","phylo_neg")

res_all  <- map(methods, ~ compute_auprc(merged, .x))
res_mono <- map(methods, ~ compute_auprc(df_mono, .x))

names(res_all)  <- c("IBD","IBS","Phylo")
names(res_mono) <- names(res_all)

# ═══════════════════════════════════════════════════════════════════
# 5. SIMULATED UPPER BOUNDS
# ═══════════════════════════════════════════════════════════════════
# SIM <- tibble(
#   Method = c("IBD","IBS","Phylo"),
#   mean = c(0.609,0.529,0.63),
#   lo   = c(0.65,0.60,0.58),
#   hi   = c(0.75,0.70,0.68)
# )

# ----------------------
# Aggregation utilities for plotting
# ----------------------
aggregate_results <- function(all_results) {
  rows <- list(); k <- 1
  for (rep in names(all_results)) {
    df <- NULL
    rep_res <- all_results[[rep]]$curve_data
    if (is.null(rep_res)) next
    for (m in names(rep_res)) {
      auc <- round(rep_res[[m]]$auc_pr, 3)
      df <- rbind.data.frame(df, data.table(replicate = rep, method = m, aupr = auc))
    }
    
    rows[[k]] <- df
    k <- k + 1
  }
  
  return(rows)
}

single_results <- readRDS("simulations/single_run/evaluation_results.rds")

xx <- aggregate_results(single_results)

df <- rbindlist(xx)

SIM <- df %>% 
  group_by(method) %>% 
  summarise(mean = mean(aupr),
            q25   = quantile(aupr, 0.25),
            q75   = quantile(aupr, 0.75),
            lo = quantile(aupr, 0.025),
            hi = quantile(aupr, 0.975),
            .groups = "drop"
            ) %>% 
  mutate(method = recode(method, "phylo" = "Phylo"))

writexl::write_xlsx(SIM, "results/tables/simdata_auprc_summary.xlsx")

# ═══════════════════════════════════════════════════════════════════
# 6. FIGURE PANELS
# ═══════════════════════════════════════════════════════════════════
METHOD_COLORS <- c(
  IBD   = "#E07B39",
  IBS   = "#4878CF",
  Phylo = "#6ACC65"
)

plot_pr <- function(res_list, prev){
  map_df(names(res_list), function(m){
    curve <- as.data.frame(res_list[[m]]$curve)
    colnames(curve) <- c("Recall","Precision","Threshold")
    curve$Method <- m
    curve}) %>%
    ggplot(aes(Recall, Precision, color=Method)) +
    geom_line(size=1.5) +
    geom_hline(yintercept=prev, linetype="dashed") +
    scale_color_manual(values = METHOD_COLORS) +
    coord_cartesian(ylim=c(0, prev*12)) +
    theme_bw() + 
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black", face = "bold"),
      axis.title = element_text(size = 16, color = "black", face = "bold"),
      axis.text = element_text(color = "black", size = 14),
      axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
      axis.ticks = element_line(color = 'black', linewidth = .7),
      axis.ticks.length = unit(.22, "cm"))
}

pA <- plot_pr(res_all, prev_all) +
  ggtitle("A  PR — All pairs")

pB <- plot_pr(res_mono, prev_mono) +
  ggtitle("B  PR — Monoclonal only")

real_df <- tibble(
  method = names(res_all),
  auc = map_dbl(res_all, "auc"),
  lo  = map_dbl(res_all, "lo"),
  hi  = map_dbl(res_all, "hi"))

writexl::write_xlsx(real_df, "results/tables/realdata_auprc_summary.xlsx")

# C_IBD   = "#C0392B"    # deep red
# C_IBS   = "#2471A3"    # steel blue
# C_PHYLO = "#1E8449"    # forest green
# C_RAND  = "#95A5A6"    # slate grey

xx <- real_df %>%
  left_join(SIM, by="method") %>%
  mutate(method = factor(method, levels = c("Phylo", "IBS","IBD")),
         mid_y = (lo.y + hi.y) / 2)

pC <- xx %>% 
  ggplot(aes(y = method)) +
  
  geom_rect(aes(xmin = lo.y, xmax = hi.y,
                ymin = as.numeric(factor(method))+0.1,
                ymax = as.numeric(factor(method))+0.4,
                fill = method, alpha = .5)) + # fill="grey85"
  
  # 🔹 Central vertical line inside rectangle
  geom_segment(aes( x = mid_y, xend = mid_y, 
                    y = as.numeric(method) + 0.1,
                    yend = as.numeric(method) + 0.4,
                    color = method),
               linewidth = 1.2) +
  
  geom_point(aes(x = auc, color = method), size = 5) +
  geom_errorbarh(aes(xmin = lo.x, xmax = hi.x, color = method), height = 0.2) +
  geom_vline(xintercept = prev_all, linetype = "dashed") +
  scale_fill_manual(values = METHOD_COLORS) +
  scale_color_manual(values = METHOD_COLORS) +
  xlab("AUPRC") + ylab(NULL) +
  theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 16, color = "black", face = "bold"),
      axis.text = element_text(color = "black", size = 14),
      axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
      axis.ticks = element_line(color = 'black', linewidth = .7),
      axis.ticks.length = unit(.22, "cm")) +
  
  # After your existing geoms
  
  annotate("text", x = max(xx$hi.y)*0.9, y = 3.0, label = "Simulated",
           size = 3, color = "gray60", fontface = "bold") +

  annotate("text", x = max(xx$hi.x) * 1.12,
           y = 2.87, label = "Observed",
           size = 3, color = "gray60", fontface = "bold")

print(pC)

same_ibd = merged %>% count(same_hh) %>% filter(same_hh == 1) %>% pull(n)
diff_ibd = merged %>% count(same_hh) %>% filter(same_hh == 0) %>% pull(n)
med_same = merged %>% group_by(same_hh) %>% summarise(median = median(hmm), .groups = "drop") %>% filter(same_hh == 1) %>% pull(median)
med_diff = merged %>% group_by(same_hh) %>% summarise(median = median(hmm), .groups = "drop") %>% filter(same_hh == 0) %>% pull(median)

pD <- merged %>%
  ggplot(aes(hmm, fill = factor(same_hh))) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, position="identity") +
  geom_vline(xintercept = med_same, linetype = "dashed", color = "darkred", linewidth = 1) +
  geom_vline(xintercept = med_diff, linetype = "dashed", color = "gray30", linewidth = 1) +
  
  scale_fill_manual(values = c("0" = "#AAAAAA", "1" = "#E07B39"),
                    labels = c(paste0("Different-HH (n = ", diff_ibd, ")"), 
                               paste0("Same-HH (n = ", same_ibd, ")"))) +
  theme_bw() +
  labs(fill = NULL, x = "IBD proportion") +
  theme(
    legend.text = element_text(size = 10, color = "black", face = "bold"),
    legend.position = c(0.78, 0.85),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    axis.text = element_text(color = "black", size = 14),
    axis.line = element_line(linewidth = 1, colour = 'black', lineend = "square"),
    axis.ticks = element_line(color = 'black', linewidth = .7),
    axis.ticks.length = unit(.22, "cm"))

final_plot <- (pA | pB) / (pC | pD)

ggsave("realdata_validation.pdf",
       final_plot,
       width=14, height=10)


# ═══════════════════════════════════════════════════════════════════
# 7. EXPORT SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════

summary_table <- real_df %>%
  left_join(SIM, by="Method") %>%
  mutate(
    Random_baseline = prev_all,
    Lift_over_random = auc / prev_all,
    Below_sim_mean = auc < mean
  )

write_csv(summary_table, "realdata_auprc_summary.csv")















