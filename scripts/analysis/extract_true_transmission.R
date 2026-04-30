# ==========================================================
# Extract true transmission (who-infected-whom) from SLiM .trees file
# ==========================================================

# --- Load dependencies ---
suppressPackageStartupMessages({
  library(optparse)
  library(reticulate)
  
  # --- Helpers to safely convert Python scalars/arrays via reticulate ---
  py_as_integer <- function(x) {
    xr <- reticulate::py_to_r(x)
    if (length(xr) == 0) return(NA_integer_)
    as.integer(xr[[1]])
  }
  py_as_numeric <- function(x) {
    xr <- reticulate::py_to_r(x)
    if (length(xr) == 0) return(NA_real_)
    as.numeric(xr[[1]])
  }
  py_len <- function(x) {
    length(reticulate::py_to_r(x))
  }
  
  # --- Bind reticulate to the correct Python BEFORE any import() calls ---
  # Priority: RETICULATE_PYTHON > CONDA_DEFAULT_ENV > "slim-msprime"
  try({
    py <- Sys.getenv("RETICULATE_PYTHON", "")
    if (nzchar(py)) {
      use_python(py, required = TRUE)
    } else {
      env <- Sys.getenv("CONDA_DEFAULT_ENV", "")
      if (nzchar(env)) {
        use_condaenv(env, required = FALSE)
      } else {
        use_condaenv("slim-msprime", required = FALSE)
      }
    }
  }, silent = TRUE)
  
  library(igraph)
  library(ggraph)
  library(tidyverse)
})


# --- Initialize Python environment ---
py_require(c("pyslim", "msprime", "tskit", "pandas"))
tskit <- import("tskit")
pyslim <- import("pyslim")

# --- Path to your SLiM output ---
trees_file <- "sim_with_transmissions_realistic.trees"

# --- Load the tree sequence ---
ts <- tskit$load(trees_file)

cat("Loaded tree sequence with", ts$num_individuals, "individuals.\n")

# --- Extract individual metadata ---
inds <- ts$individuals()
inds_list <- reticulate::iterate(inds)

meta_list <- lapply(inds_list, function(ind) {
  data.frame(
    id = ind$id,
    p1 = ind$metadata$pedigree_p1,
    p2 = ind$metadata$pedigree_p2,
    child = ind$metadata$pedigree_id,
    deme = ind$metadata$subpopulation,
    time = ind$time
  )
})

meta_df <- bind_rows(meta_list)

# --- Extract parent-child relationships from nodes ---
edges <- data.frame()

for (i in 1:py_as_integer(inds$length)) {
  child <- inds[i-1]$metadata$pedigree_id
  p1 <- inds[i-1]$metadata$pedigree_p1
  p2 <- inds[i-1]$metadata$pedigree_p2
  
  if (p1 >= 0 & p2 >= 0) edges <- rbind(edges, data.frame(p1 = p1, p2 = p2, child = child))
}


cat("Extracted", nrow(edges), "parent-child links.\n")

# --- Join with metadata ---
edges_annot <- edges %>%
  left_join(., meta_df)

# --- Create igraph object ---
g <- graph_from_data_frame(edges_annot, directed = TRUE)

# --- Compute generation lag (infection interval) ---
# edges_annot <- edges_annot %>%
#   mutate(time_diff = child_time - parent_time)

summary(edges_annot$time_diff)

# --- Visualize network ---
set.seed(123)
ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = as.factor(deme), alpha = 0.6),
                 arrow = arrow(length = unit(2, "mm")),
                 show.legend = TRUE) +
  geom_node_point(aes(color = as.factor(V(g)$deme)), size = 3) +
  theme_minimal() +
  labs(title = "True Transmission Network (SLiM Simulation)",
       color = "Deme / Subpopulation",
       subtitle = "Edges show parent→child infection links") +
  theme(legend.position = "bottom")
