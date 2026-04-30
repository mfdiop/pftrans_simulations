

# Compute Pearson & Spearman between two columns, with optional logit transform for proportions
compute_correlations_v1 <- function(df, var1, var2, transform_for_pearson = c("none","logit","sqrt")) {
  transform_for_pearson <- match.arg(transform_for_pearson)
  
  x <- df[[var1]]; y <- df[[var2]]
  # optional transform
  if (transform_for_pearson == "logit") {
    eps <- 1e-9
    x_t <- qlogis(pmin(pmax(x, eps), 1 - eps))
    y_t <- qlogis(pmin(pmax(y, eps), 1 - eps))
  } else if (transform_for_pearson == "sqrt") {
    x_t <- sqrt(x); y_t <- sqrt(y)
  } else {
    x_t <- x; y_t <- y
  }
  pear <- cor.test(x_t, y_t, method = "pearson")
  spear <- cor.test(x, y, method = "spearman", exact = FALSE)
  return(list(pearson = pear, spearman = spear))
}

# Mantel test: requires distance matrices
# Provide pairwise matrices (square) or build from pairwise list
mantel_test_from_pairwise <- function(pairs_df, value_col, ids = NULL, nperm = 999) {
  # if ids null infer unique ids from id1/id2
  if (is.null(ids)) ids <- unique(c(pairs_df$id1, pairs_df$id2))
  n <- length(ids)
  # build distance/similarity matrix; we convert similarity->distance: d = 1-s
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(ids, ids))
  for (r in seq_len(nrow(pairs_df))) {
    i <- as.character(pairs_df$id1[r]); j <- as.character(pairs_df$id2[r])
    if (i %in% ids & j %in% ids) {
      v <- pairs_df[[value_col]][r]
      mat[i, j] <- v; mat[j, i] <- v
    }
  }
  # convert similarity to distance if necessary (assume value_col is similarity in 0..1)
  dist_mat <- as.dist(1 - mat)  # if column is distance already, skip this step
  # simple fallback if diag zero not set
  res <- vegan::mantel(dist_mat, dist_mat, permutations = nperm, method = "pearson")
  return(res)
}

## 8. Correlation tests

### Pairwise correlation
compute_correlations_v2 <- function(df) {
  list(
    pearson_IBD_IBS = cor(df$IBD, df$IBS, method = "pearson", use = "complete.obs"),
    spearman_IBD_IBS = cor(df$IBD, df$IBS, method = "spearman", use = "complete.obs"),
    pearson_IBD_phylo = cor(df$IBD, df$phylo_sim, method = "pearson", use = "complete.obs"),
    spearman_IBD_phylo = cor(df$IBD, df$phylo_sim, method = "spearman", use = "complete.obs")
  )
}

mantel_compare <- function(df_wide) {
  # df_wide is a distance/similarity matrix, not the long df
  mantel(as.dist(df_wide$IBD_dist),
         as.dist(df_wide$IBS_dist),
         permutations = 10000)
}
