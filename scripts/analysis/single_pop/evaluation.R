

library(tidyverse)
library(PRROC)


# -----------------------------------------
# 1. Load true IBD pairs (from tskitIBD)
# -----------------------------------------
ibd_data <- read.delim("true.ibd", sep = "\t")

ibd_data <- ibd_data %>% 
  mutate(Id1 = paste0("tsk_", Id1), Id2 = paste0("tsk_", Id2))


# -----------------------------------------------------------
# 2. Load observed relatedness metrics (computed from VCF)
# -----------------------------------------------------------
ibs_long <- read_tsv("single_pop/out_true_ibd/inferred_ibs.tsv")
# expected columns: Id1, Id2, metric (e.g., IBS, genetic_distance, etc.)

label_data %>%
  group_by(is_true) %>%
  summarise(
    mean_metric = mean(IBS),
    median_metric = median(IBS),
    .groups = 'drop')

# is_true mean_metric median_metric
# <int>       <dbl>         <dbl>
#   1       0       0.767          0.76
# 2       1       0.761          0.76

# Your metric values are essentially identical for true and false pairs:
#   
# False pairs (is_true=0): mean = 0.767, median = 0.76
# True pairs (is_true=1): mean = 0.761, median = 0.76

# The difference is only 0.006 (0.6%), which is negligible. Your metric cannot distinguish 
# between recently transmitted pairs and unrelated pairs.
# 

# ------- End of IBS precision -------

# -------- IBD evaluation --------
# 1. HMM
compute_pr(ibd_data, inferred_ibd_hmm, "hmm")

# -------------------
ibd_label %>%
  group_by(is_true) %>%
  summarise(
    mean_metric = mean(hmm),
    median_metric = median(hmm),
    .groups = 'drop'
  )

true_ibd <- ibd_data %>% select(Id1, Id2)
obs_ibd <- inferred_ibd_hmm %>% select(Id1, Id2)
all_pairs <- ibd_data %>% select(Id1, Id2) %>% as.data.frame()

auc <- compute_pr(true_ibd, obs_ibd, all_pairs)
print(auc)


# -----------  perform ROC curve analysis --------------
roc_ibs <- compute_roc(label_data, metric_col = "IBS")
roc_result <- compute_roc(ibd_label, metric_col = "hmm")

# ------------- Correlations ---------------
library(stats)

# Create upper triangular mask (excluding diagonal)
mask <- upper.tri(true_ibd_matrix, diag = FALSE)

# Extract upper triangular values
true_vals <- true_ibd_matrix[mask]
inferred_vals <- inferred_ibd_matrix[mask]

# Calculate Spearman correlation
correlation <- cor.test(true_vals, inferred_vals, method = "spearman")

# Print result
cat(sprintf("Spearman correlation: %.3f\n", correlation$estimate))

# Just get the correlation coefficient without the test:
# If you only need the correlation value (faster)
correlation_value <- cor(true_vals, inferred_vals, method = "spearman")
cat(sprintf("Spearman correlation: %.3f\n", correlation_value))

# To also get the p-value:
correlation <- cor.test(true_vals, inferred_vals, method = "spearman")
cat(sprintf("Spearman correlation: %.3f (p-value: %.3e)\n", 
            correlation$estimate, 
            correlation$p.value))


overlap_score <- topk_overlap(true_ibd_matrix, inferred_ibd_matrix, k = 100)
cat(sprintf("Top-k overlap: %.3f\n", overlap_score))




# Usage
ibs_roc_result <- compute_roc_ggplot(label_data, metric_col = "IBS", direction = "<")
hmm_roc_result <- compute_roc_ggplot(ibd_label, metric_col = "hmm", direction = "<")


# Usage example
# Assuming you have different metrics
obs_ibs <- obs_rel  # IBS metric
obs_distance <- obs_rel %>% mutate(metric = 1 - metric)  # Distance metric

results <- compare_roc_curves(
  data_list = list(obs_ibs, obs_distance),
  metric_names = c("IBS", "Genetic Distance")
)

# Usage
roc_analysis <- perform_roc_analysis(
  true_ibd = ibd_data,
  obs_rel = obs_rel,
  g = 5,
  metric_col = "metric"
)


# ROC vs PR Curve Comparison
comparison <- compare_roc_pr(obs_rel, metric_col = "metric")











