
# FUNCTION DOCUMENTATION TEMPLATE:
# 
# phylo_long(treefile, method = "linear")
#   Purpose: Calculate phylogenetic distances from tree file
#   Inputs: 
#     - treefile: Path to phylogenetic tree file
#     - method: Distance calculation method ("linear" or "gaussian")
#   Output: Dataframe with pairwise phylogenetic distances
#
# merge_metrics(ibd_df, ibs_df, phylo_df)
#   Purpose: Merge different inference metrics into unified dataframe
#   Inputs: Dataframes from IBD, IBS, and phylogenetic analyses
#   Output: Merged dataframe with all inference metrics
#
# define_ground_truth(true_links, gen_cutoff)
#   Purpose: Classify relationships as true/false based on generational distance
#   Inputs:
#     - true_links: Ground truth relationship data
#     - gen_cutoff: Maximum generations to consider as "related"
#   Output: Dataframe with binary ground truth labels
#
# plot_scatter() / plot_hex()
#   Purpose: Visualization functions for relationship inference results
#   Input: Merged metrics dataframe
#   Output: ggplot2 visualization objects
#
# compute_correlations_v1()
#   Purpose: Calculate correlation statistics between inference methods
#   Inputs: 
#     - data: Merged metrics dataframe
#     - method1, method2: Names of methods to compare
#     - transform_for_pearson: Transformation to apply ("logit", "log", etc.)
#   Output: List containing Pearson and Spearman correlation results

compute_phylo_distances <- function(treefile) {
  tree <- ape::read.tree(treefile)
  D <- ape::cophenetic.phylo(tree)  # matrix of patristic distances
  return(D)
}

matrix_to_long_efficient <- function(matrix_data, value_name = "value", id_cols = c("id1", "id2")) {
  n <- nrow(matrix_data)
  
  # Get upper triangle indices
  upper_tri_indices <- which(upper.tri(matrix_data), arr.ind = TRUE)
  
  # Extract values and create dataframe
  unique_pairs <- data.frame(
    id1 = rownames(matrix_data)[upper_tri_indices[, 1]],
    id2 = colnames(matrix_data)[upper_tri_indices[, 2]],
    value = matrix_data[upper_tri_indices]
  )
  
  names(unique_pairs) <- c(id_cols, value_name)
  
  cat("Converted", n, "×", n, "matrix to", nrow(unique_pairs), "unique pairs\n")
  return(unique_pairs)
}

distance_to_long <- function(distmat, metric_name = "dist") {
  # Input validation
  if (!is.matrix(distmat)) {
    stop("Input must be a matrix")
  }
  
  n <- nrow(distmat)
  if (n != ncol(distmat)) {
    stop("Input matrix must be square")
  }
  
  # Get only upper triangle (excluding diagonal) to avoid duplicates
  upper_tri_indices <- which(upper.tri(distmat), arr.ind = TRUE)
  
  # Create long format dataframe with unique pairs only
  long_df <- data.frame(
    id1 = rownames(distmat)[upper_tri_indices[, 1]],
    id2 = colnames(distmat)[upper_tri_indices[, 2]], 
    value = distmat[upper_tri_indices],
    metric = metric_name,
    stringsAsFactors = FALSE
  )
  
  # Verify we have the correct number of pairs
  expected_pairs <- n * (n - 1) / 2
  actual_pairs <- nrow(long_df)
  
  cat("Matrix dimension:", n, "x", n, "\n")
  cat("Expected unique pairs:", expected_pairs, "\n")
  cat("Actual unique pairs:", actual_pairs, "\n")
  
  if (actual_pairs != expected_pairs) {
    warning("Unexpected number of pairs. Check matrix structure.")
  }
  
  return(long_df)
}


patristic_dist_to_similarity <- function(distvec) {
  # Input validation
  if (!is.numeric(distvec)) {
    stop("Input must be a numeric vector")
  }
  
  if (any(distvec < 0, na.rm = TRUE)) {
    warning("Negative distances found. Taking absolute values.")
    distvec <- abs(distvec)
  }
  
  sim <- 1 / (1 + distvec)
  sim[sim < 0] <- 0
  sim[is.na(sim)] <- 0  # Handle NA values
  
  return(sim)
}


phylo_dist_to_similarity <- function(distvec) {
  # Input validation
  if (!is.numeric(distvec)) {
    stop("Input must be a numeric vector")
  }
  
  if (any(distvec < 0, na.rm = TRUE)) {
    warning("Negative distances found. Taking absolute values.")
    distvec <- abs(distvec)
  }
  
  maxd <- max(distvec, na.rm = TRUE)
  
  if (maxd == 0) {
    warning("All distances are zero. Returning vector of ones.")
    return(rep(1, length(distvec)))
  }
  
  sim <- 1 - (distvec / maxd)
  sim[sim < 0] <- 0
  sim[is.na(sim)] <- 0  # Handle NA values
  
  return(sim)
}

phylo_dist_gaussian <- function(distvec, sigma = NULL) {
  # Input validation
  if (!is.numeric(distvec)) {
    stop("Input must be a numeric vector")
  }
  
  if (any(distvec < 0, na.rm = TRUE)) {
    warning("Negative distances found. Taking absolute values.")
    distvec <- abs(distvec)
  }
  
  if (is.null(sigma)) {
    sigma <- median(distvec, na.rm = TRUE)
  }
  
  if (sigma <= 0) {
    warning("Sigma is <= 0. Using default sigma = 1")
    sigma <- 1
  }
  
  similarity <- exp(-(distvec^2) / (2 * sigma^2))
  similarity[is.na(similarity)] <- 0  # Handle NA values
  
  return(similarity)
}

#' Transform patristic distances into bounded similarity (0..1)
#'
#' @param d numeric vector of patristic distances
#' @param method "exp" or "inv" or "rank" (rank returns scaled ranks 0..1)
#' @param alpha positive numeric for exp transform (if NULL it's auto-estimated)
#' @param beta positive numeric for inverse transform (if NULL set beta = 1/median(d))
#' @return numeric vector similarity (0..1)
transform_phylo_to_similarity <- function(d,
                                          method = c("exp", "inv", "rank"),
                                          alpha = NULL,
                                          beta = NULL) {
  method <- match.arg(method)
  d <- as.numeric(d)
  if (method == "exp") {
    if (is.null(alpha)) {
      med <- median(d[is.finite(d) & !is.na(d)])
      # place median distance at sim ~ 0.5 -> alpha = log(2)/med
      alpha <- ifelse(med > 0, log(2) / med, 1)
    }
    sim <- exp(-alpha * d)
  } else if (method == "inv") {
    if (is.null(beta)) {
      med <- median(d[is.finite(d) & !is.na(d)])
      beta <- ifelse(med > 0, 1 / med, 1)
    }
    sim <- 1 / (1 + beta * d)
  } else if (method == "rank") {
    # normalized rank: highest similarity for smallest distance
    r <- rank(d, ties.method = "average", na.last = "keep")
    sim <- 1 - ((r - 1) / (max(r, na.rm = TRUE) - 1 + 1e-12))  # scale to 0..1
  } else sim <- 1 / (1 + d)
  
  sim[is.na(sim)] <- 0
  return(sim)
}


phylo_long <- function(treefile, method = "linear") {
  # Read tree and get distance matrix
  dist_matrix <- compute_phylo_distances(treefile)
  
  n_individuals <- nrow(dist_matrix)
  cat("Processing", n_individuals, "individuals\n")
  
  # Convert to long format with unique pairs only
  # df <- distance_to_long(dist_matrix, metric_name = "patristic")
  df <- matrix_to_long_efficient(dist_matrix, value_name = "phylo")
  
  if(method == "patristic"){
    df <- df %>% mutate(phylo = patristic_dist_to_similarity(phylo))
  }
  else if (method == "linear") {
    df <- df %>% mutate(phylo = phylo_dist_to_similarity(phylo))
  } 
  else if (method == "gaussian") {
    df <- df %>% mutate(phylo = phylo_dist_gaussian(phylo))
  }
  
  # # Apply similarity transformation
  # if (method == "linear") {
  #   phylo_long$value <- phylo_dist_to_similarity(phylo_long$value)
  #   phylo_long$metric <- "phylo_similarity_linear"
  # } else if (method == "gaussian") {
  #   phylo_long$value <- phylo_dist_gaussian(phylo_long$value)
  #   phylo_long$metric <- "phylo_similarity_gaussian"
  # }
  
  # Verify final pair count
  expected_pairs <- n_individuals * (n_individuals - 1) / 2
  actual_pairs <- nrow(df)
  
  cat("Final unique pairs:", actual_pairs, "/ Expected:", expected_pairs, "\n")
  
  if (actual_pairs != expected_pairs) {
    stop("Mismatch in pair counts!")
  }
  
  return(df)
}
