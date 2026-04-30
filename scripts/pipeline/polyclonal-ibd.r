#' Extensions for handling polyclonal infections in IBD estimation
#'
#' Three approaches:
#' 1. Within-host diversity filtering (Fws)
#' 2. Probabilistic IBD with COI
#' 3. Dominant clone IBD
#'
#' @param vcf_data VCF data object
#' @param method One of "fws_filter", "probabilistic", "dominant_clone"

# ============================================================================
# Approach 1: Filter by Within-Host Diversity (Fws)
# ============================================================================

#' Calculate within-host diversity (Fws statistic)
#'
#' Fws measures within-host diversity
#' Fws = 1 - (Het_obs / Het_exp)
#' Fws ≈ 1: monoclonal (low heterozygosity)
#' Fws < 0.95: likely polyclonal
#'
#' @param genotypes Numeric genotype matrix (variants x samples)
#' @param af Allele frequencies
#' @return Vector of Fws values per sample

calculate_fws <- function(genotypes, af) {
  n_samples <- ncol(genotypes)
  fws <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    gt <- genotypes[, i]
    valid <- !is.na(gt)
    
    # Observed heterozygosity
    het_obs <- mean(gt[valid] == 1)  # Assuming 0=hom ref, 1=het, 2=hom alt
    
    # Expected heterozygosity under random mating
    het_exp <- mean(2 * af[valid] * (1 - af[valid]))
    
    # Fws statistic
    fws[i] <- 1 - (het_obs / het_exp)
  }
  
  return(fws)
}


#' Filter samples by Fws (monoclonal only)
#'
#' @param vcf_data VCF data object
#' @param fws_threshold Minimum Fws to keep (default: 0.95)
#' @return List with filtered data and Fws values

filter_monoclonal_samples <- function(vcf_data, fws_threshold = 0.95) {
  
  # Calculate allele frequencies
  af <- rowMeans(vcf_data$genotypes_numeric, na.rm = TRUE)
  
  # Calculate Fws for each sample
  fws <- calculate_fws(vcf_data$genotypes_numeric, af)
  
  # Identify monoclonal samples
  monoclonal <- fws >= fws_threshold
  
  cat(sprintf("Fws filtering: %d/%d samples are monoclonal (Fws >= %.2f)\n",
              sum(monoclonal), length(monoclonal), fws_threshold))
  
  # Filter VCF data
  vcf_filtered <- vcf_data
  vcf_filtered$genotypes_numeric <- vcf_data$genotypes_numeric[, monoclonal]
  vcf_filtered$n_samples <- sum(monoclonal)
  vcf_filtered$metadata$sample_id <- vcf_data$metadata$sample_id[monoclonal]
  
  return(list(
    vcf_data = vcf_filtered,
    fws_values = fws,
    monoclonal_samples = vcf_data$metadata$sample_id[monoclonal],
    polyclonal_samples = vcf_data$metadata$sample_id[!monoclonal]
  ))
}


# ============================================================================
# Approach 2: Probabilistic IBD with Complexity of Infection (COI)
# ============================================================================

#' Estimate Complexity of Infection (COI) from heterozygosity
#'
#' Simple COI estimation based on heterozygosity patterns
#' More sophisticated methods exist (e.g., THE REAL McCOIL, COIL)
#'
#' @param genotypes Genotype vector for one sample
#' @param af Allele frequencies
#' @return Estimated COI (integer)

estimate_coi_simple <- function(genotypes, af) {
  valid <- !is.na(genotypes)
  
  # Heterozygous sites
  n_het <- sum(genotypes[valid] == 1)
  n_total <- sum(valid)
  
  het_rate <- n_het / n_total
  
  # Expected het rate for COI clones
  # E[het] ≈ 1 - (1/COI) for equal mixture
  # Solve: COI ≈ 1 / (1 - het_rate)
  
  if (het_rate < 0.05) {
    return(1)  # Monoclonal
  } else if (het_rate < 0.3) {
    return(2)
  } else if (het_rate < 0.5) {
    return(3)
  } else {
    return(4)  # 4+ clones
  }
}


#' Calculate probabilistic IBD for polyclonal samples
#'
#' For samples with COI > 1, calculate probability that at least one
#' clone pair shares IBD
#'
#' @param gt1 Genotype vector sample 1
#' @param gt2 Genotype vector sample 2
#' @param coi1 COI for sample 1
#' @param coi2 COI for sample 2
#' @param af Allele frequencies
#' @return Probabilistic IBD score

calculate_probabilistic_ibd <- function(gt1, gt2, coi1, coi2, af) {
  
  valid <- !is.na(gt1) & !is.na(gt2)
  gt1 <- gt1[valid]
  gt2 <- gt2[valid]
  af <- af[valid]
  
  # For each site, calculate P(at least one clone pair matches)
  prob_match <- numeric(length(gt1))
  
  for (i in seq_along(gt1)) {
    g1 <- gt1[i]
    g2 <- gt2[i]
    p <- af[i]
    
    # If both samples monoclonal (COI=1)
    if (coi1 == 1 && coi2 == 1) {
      prob_match[i] <- as.numeric(g1 == g2)
      
    # If one is polyclonal
    } else {
      # P(match) depends on genotypes observed
      
      # Both homozygous same allele
      if ((g1 == 0 && g2 == 0) || (g1 == 2 && g2 == 2)) {
        prob_match[i] <- 1.0  # Definitely share this allele
        
      # Both heterozygous  
      } else if (g1 == 1 && g2 == 1) {
        # Both have 0 and 1 alleles in mixture
        # P(at least one clone from each matches) is high
        prob_match[i] <- 0.75  # Simplified
        
      # One het, one hom
      } else if ((g1 == 1 && g2 %in% c(0, 2)) || 
                 (g2 == 1 && g1 %in% c(0, 2))) {
        # Het sample has both alleles, hom has only one
        # P(het sample has clone matching hom)
        prob_match[i] <- 0.5
        
      # Different homozygous
      } else {
        prob_match[i] <- 0.0
      }
    }
  }
  
  # Weight by allele frequency informativeness
  weights <- 2 * af * (1 - af)
  
  # Weighted average
  ibd_score <- sum(prob_match * weights) / sum(weights)
  
  return(ibd_score)
}


#' IBD estimation accounting for COI
#'
#' @param vcf_data VCF data object
#' @param min_ibd_threshold Minimum IBD score to report
#' @return IBD results with COI information

estimate_ibd_with_coi <- function(vcf_data, 
                                  min_ibd_threshold = 0.3) {
  
  n_samples <- vcf_data$n_samples
  sample_ids <- vcf_data$metadata$sample_id
  genotypes <- vcf_data$genotypes_numeric
  
  # Calculate allele frequencies
  af <- rowMeans(genotypes, na.rm = TRUE)
  
  # Estimate COI for each sample
  cat("Estimating COI for each sample...\n")
  coi <- sapply(1:n_samples, function(i) {
    estimate_coi_simple(genotypes[, i], af)
  })
  
  cat(sprintf("COI distribution: %s\n", 
              paste(table(coi), collapse = ", ")))
  
  # Initialize results
  ibd_matrix <- matrix(0, n_samples, n_samples)
  rownames(ibd_matrix) <- colnames(ibd_matrix) <- sample_ids
  diag(ibd_matrix) <- 1
  
  results <- list()
  pair_count <- 0
  
  # Calculate pairwise IBD
  cat("Calculating pairwise probabilistic IBD...\n")
  for (i in 1:(n_samples - 1)) {
    if (i %% 10 == 0) {
      cat(sprintf("  Processing sample %d/%d\n", i, n_samples))
    }
    
    for (j in (i + 1):n_samples) {
      
      gt1 <- genotypes[, i]
      gt2 <- genotypes[, j]
      
      # Calculate probabilistic IBD
      ibd_score <- calculate_probabilistic_ibd(
        gt1, gt2, coi[i], coi[j], af
      )
      
      if (ibd_score >= min_ibd_threshold) {
        ibd_matrix[i, j] <- ibd_score
        ibd_matrix[j, i] <- ibd_score
        
        pair_count <- pair_count + 1
        results[[pair_count]] <- data.frame(
          sample1 = sample_ids[i],
          sample2 = sample_ids[j],
          coi1 = coi[i],
          coi2 = coi[j],
          ibd_score = ibd_score
        )
      }
    }
  }
  
  results_df <- if (pair_count > 0) {
    do.call(rbind, results)
  } else {
    data.frame()
  }
  
  return(list(
    ibd_matrix = ibd_matrix,
    ibd_pairs = results_df,
    coi_estimates = data.frame(
      sample_id = sample_ids,
      coi = coi
    )
  ))
}


# ============================================================================
# Approach 3: Dominant Clone IBD
# ============================================================================

#' Infer dominant clone genotypes from polyclonal sample
#'
#' For polyclonal samples, attempt to infer the dominant (major) clone
#' by assuming homozygous calls represent the dominant clone
#'
#' @param genotypes Genotype matrix
#' @return Inferred haploid genotypes for dominant clones

infer_dominant_clone <- function(genotypes) {
  
  # Strategy: 
  # - Keep homozygous calls (0 or 2) as is
  # - For heterozygous calls (1), set to NA (ambiguous)
  # - This gives us partial genotype of dominant clone
  
  dominant <- genotypes
  dominant[genotypes == 1] <- NA  # Het sites are ambiguous
  dominant[genotypes == 2] <- 1   # Convert to haploid (0/1)
  
  return(dominant)
}


#' IBD estimation using dominant clone inference
#'
#' @param vcf_data VCF data object
#' @param min_overlap Minimum fraction of non-NA overlap required
#' @return IBD results

estimate_ibd_dominant_clone <- function(vcf_data,
                                        min_overlap = 0.5) {
  
  cat("Inferring dominant clones...\n")
  
  # Infer dominant clone genotypes
  dominant_genotypes <- infer_dominant_clone(vcf_data$genotypes_numeric)
  
  # Calculate IBD on dominant clones only
  vcf_dominant <- vcf_data
  vcf_dominant$genotypes_numeric <- dominant_genotypes
  
  # Use the robust IBD function
  ibd_results <- estimate_ibd_robust(
    vcf_dominant,
    min_segment_length = 1e6,  # More lenient for partial data
    min_lod = 2,               # Lower threshold
    min_snps = 15,
    max_missing = 1 - min_overlap  # Allow more missing
  )
  
  return(ibd_results)
}


# ============================================================================
# Wrapper function for all approaches
# ============================================================================

#' Estimate IBD with polyclonal handling
#'
#' @param vcf_data VCF data object  
#' @param method One of "fws_filter", "probabilistic", "dominant_clone"
#' @param ... Additional parameters passed to specific methods
#' @export

estimate_ibd_polyclonal <- function(vcf_data, 
                                    method = c("fws_filter", 
                                             "probabilistic", 
                                             "dominant_clone"),
                                    ...) {
  
  method <- match.arg(method)
  
  cat(sprintf("=== IBD Estimation with Polyclonal Handling ===\n"))
  cat(sprintf("Method: %s\n\n", method))
  
  result <- switch(method,
    fws_filter = {
      filtered <- filter_monoclonal_samples(vcf_data, ...)
      ibd <- estimate_ibd_robust(filtered$vcf_data, ...)
      list(
        ibd_results = ibd,
        fws_values = filtered$fws_values,
        monoclonal_samples = filtered$monoclonal_samples,
        polyclonal_samples = filtered$polyclonal_samples,
        method = "fws_filter"
      )
    },
    
    probabilistic = {
      ibd <- estimate_ibd_with_coi(vcf_data, ...)
      list(
        ibd_results = ibd,
        method = "probabilistic"
      )
    },
    
    dominant_clone = {
      ibd <- estimate_ibd_dominant_clone(vcf_data, ...)
      list(
        ibd_results = ibd,
        method = "dominant_clone"
      )
    }
  )
  
  return(result)
}


# ============================================================================
# Usage Examples
# ============================================================================

#' Example usage:
#' 
#' # Approach 1: Filter to monoclonal only (most reliable)
#' result1 <- estimate_ibd_polyclonal(vcf_data, method = "fws_filter", 
#'                                    fws_threshold = 0.95)
#' 
#' # Approach 2: Probabilistic IBD (keeps all samples)
#' result2 <- estimate_ibd_polyclonal(vcf_data, method = "probabilistic",
#'                                    min_ibd_threshold = 0.3)
#' 
#' # Approach 3: Dominant clone (partial information)
#' result3 <- estimate_ibd_polyclonal(vcf_data, method = "dominant_clone",
#'                                    min_overlap = 0.5)
