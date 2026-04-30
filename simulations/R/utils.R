
#' Source all R scripts in the R directory
#'
#' @param exclude Character vector of script names to exclude
#' @export
source_all_scripts <- function(exclude = c(normalizePath("simulations/R/utils.R", winslash = "\\"))) {
  r_files <- list.files("simulations/R", pattern = "\\.[Rr]$", full.names = TRUE)
  r_files <- r_files[!basename(r_files) %in% exclude]
  
  for (file in r_files) {
    cat("Sourcing:", file, "\n")
    source(file)
  }
  
  invisible(TRUE)
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

#' Merge multiple relationship inference metrics with proper handling of different pairs
#'
#' @param ground_truth_df Dataframe with ground truth pairs
#' @param ibd_df IBD inference results
#' @param ibs_df IBS inference results  
#' @param phylo_df Phylogenetic distance results
#' @param merge_method How to handle different pair sets: 
#'        "intersection" (only pairs present in all), 
#'        "ground_truth" (only pairs with ground truth),
#'        "union" (all pairs with NAs)
#' @param id_cols Column names for individual IDs (default: c("id1", "id2"))
#'
#' @return Merged dataframe with consistent pair set
merge_metrics <- function(ibd_df, ibs_df, phylo_df) {
  colnames(ibd_df)[1:2] <- c('id1', 'id2')
  colnames(ibs_df)[1:2] <- c('id1', 'id2')
  
  ibd_df2  <- ibd_df  %>% rename(IBD = hmm) %>% mutate(IBD = round(IBD, 3))
  ibs_df2  <- ibs_df  #%>% rename(IBS = ibs_prop)
  phylo_df2 <- phylo_df %>% select(id1, id2, phylo) %>% 
    mutate(phylo = round(phylo, 3))
  
  merged <- ibd_df2 %>%
    full_join(ibs_df2,  by = c("id1", "id2")) %>%
    full_join(phylo_df2, by = c("id1", "id2"))
  
  merged
}


add_truth <- function(metric_df, true_links) {
  metric_df %>%
    left_join(true_links, by = c("id1", "id2"))
}

define_ground_truth <- function(df_ground, gen_cutoff = 1) {
  df_ground$true_link <- df_ground$generations <= gen_cutoff
  df_ground$true_link <- as.numeric(df_ground$true_link) # ifelse(string_vector == "Yes", 1, 0)
  df_ground
}


# Create a function to standardize the order
standardize_pair <- function(df, col1, col2) {
  df <- df %>%
    mutate(
      temp_min = pmin({{col1}}, {{col2}}),
      temp_max = pmax({{col1}}, {{col2}})
    ) %>%
    select(-c({{col1}}, {{col2}})) %>%
    rename(id1 = temp_min, id2 = temp_max)
  
  return(df)
}

# Standardize pair representation across all dataframes
# standardize_pairs <- function(df, source_name) {
#   if (!all(id_cols %in% names(df))) {
#     stop("ID columns not found in ", source_name)
#   }
#   
#   df <- df %>%
#     rowwise() %>%
#     mutate(
#       pair_key = paste(sort(c(!!sym(id_cols[1]), !!sym(id_cols[2]))), collapse = "_")
#     ) %>%
#     ungroup() %>%
#     distinct(pair_key, .keep_all = TRUE)  # Remove duplicates
#   
#   return(df)
# }

# STANDARDIZE PAIRS ========================================================

standardize_pairs <- function(df, source_name, id_cols) {
  df_std <- df %>%
    rowwise() %>%
    mutate(
      pair_key = paste(sort(c(!!sym(id_cols[1]), !!sym(id_cols[2]))), collapse = "_")
    ) %>%
    ungroup() %>%
    distinct(pair_key, .keep_all = TRUE)  # Remove duplicates
  
  cat("✓ ", source_name, ": ", nrow(df_std), " unique pairs\n", sep = "")
  return(df_std)
}

# # Find rows where the pair appears in both directions
# find_reversed_pairs <- function(df, id_cols = c("id1", "id2")) {
#   
#   # Create standardized pair key for comparison
#   df <- df %>%
#     rowwise() %>%
#     mutate(
#       pair_key_forward = paste(!!sym(id_cols[1]), !!sym(id_cols[2]), sep = "_"),
#       pair_key_reversed = paste(!!sym(id_cols[2]), !!sym(id_cols[1]), sep = "_")
#     ) %>%
#     ungroup()
#   
#   # Find pairs that exist in both directions
#   reversed_pairs <- df %>%
#     filter(pair_key_forward %in% df$pair_key_reversed) %>%
#     arrange(pair_key_forward)
#   
#   return(reversed_pairs)
# }
# 
# # Usage
# reversed_rows <- find_reversed_pairs(ibs_df, id_cols = c("Id1", "Id2"))


#' Merge multiple relationship inference metrics with proper handling of different pairs
#'
#' @param ground_truth_df Dataframe with ground truth pairs
#' @param ibd_df IBD inference results
#' @param ibs_df IBS inference results  
#' @param phylo_df Phylogenetic distance results
#' @param merge_method How to handle different pair sets: 
#'        "intersection" (only pairs present in all), 
#'        "ground_truth" (only pairs with ground truth),
#'        "union" (all pairs with NAs)
#' @param id_cols Column names for individual IDs (default: c("id1", "id2"))
#'
#' @return Merged dataframe with consistent pair set
#' 
#' @examples
#' # Basic usage
#' merged <- merge_metrics_robust(
#'   ground_truth_df = df_ground,
#'   ibd_df = ibd_df,
#'   ibs_df = ibs_df,
#'   phylo_df = phylo_df,
#'   merge_method = "ground_truth"
#' )
#' 
#' # See available merge methods
#' merge_metrics_robust(show_methods = TRUE)

merge_metrics_with <- function(df_ground, ibd_df, ibs_df, phylo_df,
                               merge_method = "ground_truth",
                               id_cols = c("id1", "id2"),
                               show_methods = FALSE) {
  
  # HELP: Show available methods if requested
  if (show_methods) {
    cat("Available merge methods:\n")
    cat("• 'intersection': Only pairs present in ALL datasets (most conservative)\n")
    cat("• 'ground_truth': Only pairs with ground truth labels (recommended for evaluation)\n") 
    cat("• 'union': All pairs from any dataset (with NAs for missing data)\n")
    cat("\nUsage example:\n")
    cat('merge_metrics_robust(ground_truth_df, ibd_df, ibs_df, phylo_df, merge_method = "ground_truth")\n')
    return(invisible(NULL))
  }
  
  # INPUT VALIDATION ==========================================================
  
  # Check if required dataframes are provided
  if (is.null(df_ground)) {
    stop("ERROR: 'ground_truth_df' is required but missing.\n",
         "Please provide a dataframe with ground truth labels.")
  }
  
  if (is.null(ibd_df)) {
    stop("ERROR: 'ibd_df' is required but missing.\n", 
         "Please provide IBD inference results.")
  }
  
  if (is.null(ibs_df)) {
    stop("ERROR: 'ibs_df' is required but missing.\n",
         "Please provide IBS inference results.") 
  }
  
  if (is.null(phylo_df)) {
    stop("ERROR: 'phylo_df' is required but missing.\n",
         "Please provide phylogenetic distance results.")
  }
  
  # Validate merge_method
  valid_methods <- c("intersection", "ground_truth", "union")
  if (!merge_method %in% valid_methods) {
    stop("ERROR: Invalid merge_method '", merge_method, "'\n",
         "Available methods: ", paste(valid_methods, collapse = ", "), "\n",
         "Use merge_metrics_with(show_methods = TRUE) to see descriptions.")
  }
  
  # Validate id_cols
  if (length(id_cols) != 2) {
    stop("ERROR: 'id_cols' must contain exactly 2 column names (e.g., c('id1', 'id2'))")
  }
  
  if (!is.character(id_cols)) {
    stop("ERROR: 'id_cols' must be a character vector (e.g., c('id1', 'id2'))")
  }
  
  colnames(ibd_df)[1:2] <- c('id1', 'id2')
  colnames(ibs_df)[1:2] <- c('id1', 'id2')
  
  ibd_df <- ibd_df %>%
    rename(IBD = hmm) %>%
    mutate(IBD = round(IBD, 3))
  
  phylo_df <- phylo_df %>%
    select(id1, id2, phylo) %>%
    mutate(phylo = round(phylo, 3),
           id1 = as.character(id1),
           id2 = as.character(id2))
  
  # Function to check dataframe structure with helpful errors
  check_dataframe_structure <- function(df, df_name, id_cols) {
    if (!is.data.frame(df)) {
      stop("ERROR: '", df_name, "' must be a dataframe, but got: ", class(df)[1])
    }
    
    if (nrow(df) == 0) {
      stop("ERROR: '", df_name, "' is empty (0 rows)")
    }
    
    # Check for ID columns with helpful suggestions
    missing_ids <- setdiff(id_cols, names(df))
    if (length(missing_ids) > 0) {
      # Suggest possible column names
      possible_cols <- names(df)[grepl("id|ID|Id|sample|indiv", names(df), ignore.case = TRUE)]
      suggestion <- ""
      if (length(possible_cols) > 0) {
        suggestion <- paste("\nPossible ID columns in your data:", paste(possible_cols, collapse = ", "))
      }
      
      stop("ERROR: Missing required columns in '", df_name, "': ", paste(missing_ids, collapse = ", "),
           "\nYour dataframe columns: ", paste(names(df), collapse = ", "),
           suggestion,
           "\nYou can specify different column names using: id_cols = c('your_id1_col', 'your_id2_col')")
    }
    
    # Check for duplicate pairs
    pair_check <- df %>%
      rowwise() %>%
      mutate(pair_check = paste(sort(c(!!sym(id_cols[1]), !!sym(id_cols[2]))), collapse = "_")) %>%
      ungroup()
    
    dup_count <- sum(duplicated(pair_check$pair_check))
    if (dup_count > 0) {
      warning("Found ", dup_count, " duplicate pairs in '", df_name, "'. Keeping only first occurrence of each pair.")
    }
    
    return(TRUE)
  }
  
  # Validate all dataframes
  cat("Validating input dataframes...\n")
  check_dataframe_structure(df_ground, "df_ground", id_cols)
  check_dataframe_structure(ibd_df, "ibd_df", id_cols)
  check_dataframe_structure(ibs_df, "ibs_df", id_cols) 
  check_dataframe_structure(phylo_df, "phylo_df", id_cols)
  cat("✓ All dataframes validated successfully\n")
  
  # STANDARDIZE PAIRS ========================================================
  
  # # Apply to both dataframes
  # ground_standardized <- standardize_pair(df_ground, id1, id2)
  # ibd_standardized <- standardize_pair(ibd_df1, id1, id2)
  # ibs_standardized <- standardize_pair(ibs_df, id1, id2)
  # phylo_standardized <- standardize_pair(phylo_df1, id1, id2)
  
  # df <- ground_standardized |>
  #   dplyr::left_join(ibd_standardized, by = c("id1", "id2")) |>
  #   dplyr::left_join(ibs_standardized, by = c("id1", "id2")) |>
  #   dplyr::left_join(phylo_standardized, by = c("id1", "id2"))
  
  # Apply standardization
  cat("\nStandardizing pairs across datasets...\n")
  true_std <- standardize_pairs(df_ground, "true_link", id_cols)
  ibd_std <- standardize_pairs(ibd_df, "IBD", id_cols)
  ibs_std <- standardize_pairs(ibs_df, "IBS", id_cols) 
  phylo_std <- standardize_pairs(phylo_df, "phylo", id_cols)
  
  # Get all unique pair keys
  all_pairs <- unique(c(
    true_std$pair_key,
    ibd_std$pair_key, 
    ibs_std$pair_key,
    phylo_std$pair_key
  ))
  
  cat("==============================\n")
  cat("Pair counts:\n")
  cat("==============================\n")
  cat("   Ground truth:", length(unique(true_std$pair_key)), "\n")
  cat("   IBD:", length(unique(ibd_std$pair_key)), "\n")
  cat("   IBS:", length(unique(ibs_std$pair_key)), "\n")
  cat("   Phylo:", length(unique(phylo_std$pair_key)), "\n")
  cat("   All unique:", length(all_pairs), "\n")
  
  # Determine which pairs to keep based on merge method
  if (merge_method == "intersection") {
    keep_pairs <- Reduce(intersect, list(
      true_std$pair_key,
      ibd_std$pair_key,
      ibs_std$pair_key,
      phylo_std$pair_key
    ))
    if (length(keep_pairs) == 0) {
      stop("ERROR: No pairs found in the intersection of all datasets.\n",
           "Consider using merge_method = 'ground_truth' or 'union' instead.")
    }
    
    cat("Using intersection:", length(keep_pairs), "pairs\n")
    
  } else if (merge_method == "ground_truth") {
    keep_pairs <- true_std$pair_key
    cat("Using ground truth pairs:", length(keep_pairs), "pairs\n")
    
  } else if (merge_method == "union") {
    keep_pairs <- all_pairs
    cat("Using union:", length(keep_pairs), "pairs\n")
  }
  
  # Create base dataframe with all pairs to keep
  base_df <- data.frame(pair_key = keep_pairs)
    # separate(pair_key, into = id_cols, sep = "(?<=[0-9])_", extra = "merge")
    # separate(pair_key, into = id_cols, sep = "_tsk_", remove = FALSE)
  
  # Merge each metric dataframe
  # merged <- base_df %>%
  #   left_join(select(true_std, pair_key, generations, ibd, true_link = is_true), 
  #             by = "pair_key") %>%
  #   left_join(select(ibd_std, pair_key, IBD), by = "pair_key") %>%
  #   left_join(select(ibs_std, pair_key, IBS), by = "pair_key") %>%
  #   left_join(select(phylo_std, pair_key, phylo_dist), by = "pair_key")
  
  # Helper function for safe joining with progress reporting
  safe_join <- function(left_df, right_df, right_name, by = "pair_key") {
    cat("  Merging", right_name, "... ")
    
    # Select only the metric columns from right_df (excluding id columns)
    right_cols_to_join <- setdiff(names(right_df), c("id1", "id2", "pair_key"))
    
    result <- left_df %>%
      left_join(select(right_df, all_of(c(by, right_cols_to_join))), by = by)
    
    cat("done\n")
    return(result)
  }
  
  merged <- base_df %>%
    safe_join(true_std, "ground truth") %>%
    safe_join(ibd_std, "IBD") %>%
    safe_join(ibs_std, "IBS") %>%
    safe_join(phylo_std, "phylo")
  
  # FINAL VALIDATION AND REPORT ==============================================
  
  cat("\nMerge completed successfully!\n")
  cat("Final dataset: ", nrow(merged), " pairs\n\n", sep = "")
  
  cat("Missing data summary:\n")
  metrics <- c("true_link" = "Ground truth labels",
               "IBD" = "IBD inference", 
               "IBS" = "IBS inference",
               "phylo" = "Phylogenetic distance")
  
  for (metric in names(metrics)) {
    if (metric %in% names(merged)) {
      na_count <- sum(is.na(merged[[metric]]))
      na_pct <- round(mean(is.na(merged[[metric]])) * 100, 1)
      status <- ifelse(na_count == 0, "✓ Complete", paste("⚠ Missing:", na_count, "(", na_pct, "%)"))
      cat("• ", metrics[metric], ": ", status, "\n", sep = "")
    } else {
      cat("• ", metrics[metric], ": ✗ Not found in data\n", sep = "")
    }
  }
  
  # Warn if high missingness
  total_missing <- sum(is.na(merged$IBD) | is.na(merged$IBS) | is.na(merged$phylo))
  if (total_missing > 0.5 * nrow(merged)) {
    warning("\nHigh rate of missing data (>50%). Consider using merge_method = 'intersection' for more complete data.")
  }
  
  cat("\nTo see available merge methods: merge_metrics_with(show_methods = TRUE)\n")
  
  # Report missing data
  cat("==================================\n")
  cat("\nMissing data after merge:\n")
  cat("==================================\n")
  cat("   IBD:", sum(is.na(merged$IBD)), "(", round(mean(is.na(merged$IBD)) * 100, 1), "%)\n")
  cat("   IBS:", sum(is.na(merged$IBS)), "(", round(mean(is.na(merged$IBS)) * 100, 1), "%)\n")
  cat("   Phylo:", sum(is.na(merged$phylo)), "(", round(mean(is.na(merged$phylo)) * 100, 1), "%)\n")
  cat("   Ground truth:", sum(is.na(merged$true_link)), "(", round(mean(is.na(merged$true_link)) * 100, 1), "%)\n")
  
  return(merged)
}


# prepare symmetric matrices:
build_sym_matrix <- function(pairs_df, value_col, ids) {
  n <- length(ids)
  m <- matrix(NA, nrow = n, ncol = n, dimnames = list(ids, ids))
  for (r in seq_len(nrow(pairs_df))) {
    i <- as.character(pairs_df$id1[r]); j <- as.character(pairs_df$id2[r])
    if (i %in% ids & j %in% ids) {
      m[i,j] <- pairs_df[[value_col]][r]; m[j,i] <- m[i,j]
    }
  }
  diag(m) <- 0
  m[is.na(m)] <- 0
  return(m)
}


















