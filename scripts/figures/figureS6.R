
# Figure 6 demonstrates a fundamental modelling insight:
#   
#  Even when pairwise relatedness appears strong or clustering signals look tight, 
#  the population-level structure inferred by genomic methods does not match the true transmission clusters, 
#  especially under recombination and spatial mixing.
# 
# This figure is typically structured as a network-level comparison:
#   
#   Panel A — True transmission network
# 
# • Nodes = infections
# • Edges = true direct lineage links
# • Colored by population / village / deme
# • Layout using ggraph::layout_with_fr or tree-like layout
# 
# Panel B — Inferred IBD network
# 
# • Nodes = samples
# • Edges only if IBD > threshold (e.g., 0.05 or optimized threshold)
# • Shows missing edges (false negatives) and new edges (false positives)
# 
# Panel C — Cluster assignment comparison
# 
# • True cluster vs IBD-based cluster (or phylogeny-based cluster)
# • Visualized as an alluvial plot or confusion heatmap
# 
# Panel D — Cluster recovery metrics
# 
# • Adjusted Rand Index (ARI)
# • Normalized Mutual Information (NMI)
# • Misclassification rate (1 – ARI)
# 
# Supplementary
# 
# Per replicate → networks + cluster heatmaps
# Per scenario → distribution of ARI/NMI
# Effect of recombination rate on cluster recovery

library(tidyverse)
library(igraph)
library(ggraph)
library(patchwork)
library(aricode)  # for ARI, NMI

scenarios <- c("baseline", "recombination", "migration_fullfactorial")

load_rep_links <- function(scenario, rep_id) {
  
  true_path <- glue::glue("simulations/{scenario}/true/rep{rep_id}/true_transmissions.tsv")
  ibd_path  <- glue::glue("simulations/{scenario}/inferred/rep{rep_id}/hmm_ibd.tsv")
  meta_path <- glue::glue("simulations/{scenario}/metadata/rep{rep_id}/metadata.tsv")
  
  df_true <- read_tsv(true_path, show_col_types = FALSE)
  df_ibd  <- read_tsv(ibd_path, show_col_types = FALSE)
  df_meta <- read_tsv(meta_path, show_col_types = FALSE)
  
  df <- df_true %>%
    mutate(true_edge = 1, pair_key = paste(id1, id2, sep = "_")) %>%
    full_join(df_ibd %>% mutate(pair_key = paste(id1, id2, sep = "_")),
              by = c("pair_key", "id1", "id2")) %>%
    left_join(df_meta %>% rename(node = id), by = c("id1" = "node")) %>%
    rename(village1 = village, deme1 = deme) %>%
    left_join(df_meta %>% rename(node = id), by = c("id2" = "node")) %>%
    rename(village2 = village, deme2 = deme) %>%
    mutate(scenario = scenario, replicate = rep_id)
  
  return(df)
}

df_links_all <- map_dfr(
  scenarios,
  ~ map_dfr(1:5, \(r) load_rep_links(.x, r))
)


make_true_graph <- function(df) {
  g <- graph_from_data_frame(
    df %>% 
      filter(true_edge == 1) %>% 
      select(id1, id2),
    directed = TRUE,
    vertices = tibble(id = unique(c(df$id1, df$id2)))
  )
  return(g)
}

graph <- df_links_all %>%
  filter(scenario == "baseline", replicate == 1, total_ibd_prop >0.5) %>%   # choose one reference panel
  # mutate(true_edge = as.integer(direct == 1)) %>% 
  mutate(true_edge = as.integer(direct)) %>% 
  select(-c(1:6,11, 14, ))
  # make_true_graph() %>%
  # as_tbl_graph()


# graph <- graph %>%
#   activate(edges) %>%
#   mutate(edge_direct = true_edge == 1)
# 
# graph <- graph %>%
#   activate(nodes) %>%
#   mutate(
#     transmission_class = ifelse(
#       centrality_degree(
#         mode = "all",
#         weights = as.numeric(.E()[true_edge == 1]))>0,
#       "Direct",
#       "No_direct"))

my_nodes <- graph %>% 
  tidygraph::as_tbl_graph(., directed = F) %>% 
  activate(nodes) %>%
  dplyr::mutate(community = as.factor(tidygraph::group_label_prop(weights = ibd))) %>%
  as_tibble()

my_edges <- ibd_mle_long %>%
  select(-c(6:8, 12:14, 17, 21)) %>% 
  tidygraph::as_tbl_graph(., directed = F) %>% 
  activate(edges) %>%
  as_tibble() %>% 
  select(1:3)


fig6A <-   
  ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.3, color = "grey40") +
  # geom_node_point(size = 3, aes(color = factor(id))) + # later color by village
  geom_node_point(size = 3, aes(color = factor(transmission_class))) + # later color by village
  theme_void() +
  labs(title = "Figure 6A. True transmission network")


ibd_threshold <- 0.05


make_ibd_graph <- function(df, thr) {
  g <- graph_from_data_frame(
    df %>% filter(total_ibd_prop >= thr) %>% select(id1, id2),
    directed = FALSE
  )
  return(g)
}

fig6B <- df_links_all %>%
  filter(scenario == "baseline", replicate == 1) %>%
  make_ibd_graph(thr = ibd_threshold) %>%
  ggraph(layout = "fr") +
  geom_edge_link(alpha = 0.3, color = "dodgerblue3") +
  geom_node_point(size = 3, color = "black") +
  theme_void() +
  labs(title = "Figure 6B. Inferred IBD network (thresholded)")


compute_clusters <- function(df, thr) {
  g_true <- make_true_graph(df)
  g_ibd  <- make_ibd_graph(df, thr)
  
  tibble(
    id = V(g_true)$name,
    cluster_true = components(g_true)$membership[match(V(g_true)$name, names(components(g_true)$membership))],
    cluster_ibd = components(g_ibd)$membership[match(V(g_true)$name, names(components(g_ibd)$membership))]
  )
}

df_clusters <- df_links_all %>%
  filter(scenario == "baseline", replicate == 1) %>%
  compute_clusters(thr = ibd_threshold)


fig6C <- df_clusters %>%
  count(cluster_true, cluster_ibd) %>%
  ggplot(aes(x = factor(cluster_true), y = factor(cluster_ibd), fill = n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    x = "True cluster",
    y = "IBD-inferred cluster",
    title = "Figure 6C. Cluster assignment comparison"
  ) +
  theme_bw()



df_cluster_perf <- df_links_all %>%
  group_by(scenario, replicate) %>%
  summarize(
    ARI = aricode::ARI(
      compute_clusters(cur_data_all(), thr = ibd_threshold)$cluster_true,
      compute_clusters(cur_data_all(), thr = ibd_threshold)$cluster_ibd
    ),
    NMI = aricode::NMI(
      compute_clusters(cur_data_all(), thr = ibd_threshold)$cluster_true,
      compute_clusters(cur_data_all(), thr = ibd_threshold)$cluster_ibd
    ),
    .groups = "drop"
  )


fig6D <- df_cluster_perf %>%
  ggplot(aes(x = scenario, y = ARI)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  ylim(0, 1) +
  labs(
    x = "Scenario",
    y = "Adjusted Rand Index",
    title = "Figure 6D. Cluster recovery accuracy"
  ) +
  theme_bw()



supp6 <- df_links_all %>%
  filter(scenario == "baseline") %>%
  group_by(replicate) %>%
  group_map(~ {
    make_ibd_graph(.x, thr = ibd_threshold) %>%
      ggraph(layout = "fr") +
      geom_edge_link(alpha = 0.3) +
      geom_node_point(size = 2) +
      theme_void() +
      labs(title = paste("Replicate", unique(.x$replicate)))
  })


figure6 <- (fig6A | fig6B) /
  (fig6C | fig6D)

figure6










































