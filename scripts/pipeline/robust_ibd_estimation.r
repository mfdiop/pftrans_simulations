#' Robust IBD Estimation from VCF Data
#'
#' Estimates Identity by Descent (IBD) between samples using allele frequency
#' weighting and Hidden Markov Model for segment detection
#'
#' @param vcf_data List containing genotypes_numeric matrix, metadata with positions
#' @param min_segment_length Minimum segment length in base pairs (default: 2e6 for 2Mb)
#' @param min_lod Minimum LOD score for IBD segment (default: 3)
#' @param min_snps Minimum number of SNPs in a segment (default: 20)
#' @param error_rate Genotyping error rate (default: 0.01)
#' @param min_maf Minimum minor allele frequency for SNP inclusion (default: 0.01)
#' @param max_missing Maximum missing data proportion per pair (default: 0.3)
#' @return List with ibd_matrix, ibd_segments, and quality metrics

estimate_ibd_robust <- function(vcf_data,
                                min_segment_length = 2e6,  # 2 Mb
                                min_lod = 3,
                                min_snps = 20,
                                error_rate = 0.01,
                                min_maf = 0.01,
                                max_missing = 0.3) {
  
  # Extract data
  n_samples <- vcf_data$n_samples
  sample_ids <- vcf_data$metadata$sample_id
  positions <- vcf_data$metadata$position
  chromosomes <- vcf_data$metadata$chromosome
  genotypes <- vcf_data$genotypes_numeric
  
  cat("Starting robust IBD estimation...\n")
  cat(sprintf("Samples: %d, SNPs: %d\n", n_samples, nrow(genotypes)))
  
  # Step 1: Calculate allele frequencies and filter SNPs
  cat("Calculating allele frequencies...\n")
  af <- rowMeans(genotypes, na.rm = TRUE)
  maf <- pmin(af, 1 - af)
  
  # Filter SNPs by MAF
  snp_filter <- maf >= min_maf & maf <= (1 - min_maf)
  cat(sprintf("Filtered to %d SNPs (MAF >= %.3f)\n", sum(snp_filter), min_maf))
  
  genotypes <- genotypes[snp_filter, ]
  positions <- positions[snp_filter]
  chromosomes <- chromosomes[snp_filter]
  af <- af[snp_filter]
  
  # Calculate SNP weights based on informativeness
  # Weight = 2*p*(1-p) for biallelic markers (heterozygosity)
  snp_weights <- 2 * af * (1 - af)
  
  # Step 2: Initialize output structures
  ibd_matrix <- matrix(0, n_samples, n_samples)
  rownames(ibd_matrix) <- colnames(ibd_matrix) <- sample_ids
  diag(ibd_matrix) <- 1  # Perfect IBD with self
  
  ibd_segments <- list()
  segment_count <- 0
  
  # Quality metrics
  quality_metrics <- data.frame(
    sample1 = character(),
    sample2 = character(),
    n_valid_snps = integer(),
    missing_rate = numeric(),
    mean_lod = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Step 3: Calculate pairwise IBD
  cat("Calculating pairwise IBD...\n")
  n_comparisons <- n_samples * (n_samples - 1) / 2
  comparison_count <- 0
  
  for (i in 1:(n_samples - 1)) {
    if (i %% 5 == 0) {
      cat(sprintf("Processing sample %d/%d (%.1f%% complete)\n", 
                  i, n_samples, 100 * comparison_count / n_comparisons))
    }
    
    for (j in (i + 1):n_samples) {
      comparison_count <- comparison_count + 1
      
      # Get genotypes for pair
      gt1 <- genotypes[, i]
      gt2 <- genotypes[, j]
      
      # Identify valid positions
      valid <- !is.na(gt1) & !is.na(gt2)
      n_valid <- sum(valid)
      missing_rate <- 1 - (n_valid / length(valid))
      
      # Skip if too much missing data
      if (missing_rate > max_missing) {
        next
      }
      
      # Step 4: Calculate LOD scores for IBD at each SNP
      lod_scores <- calculate_lod_scores(
        gt1[valid], gt2[valid], 
        af[valid], 
        error_rate
      )
      
      # Step 5: Run HMM to identify IBD segments
      segments <- detect_ibd_segments_hmm(
        lod_scores = lod_scores,
        positions = positions[valid],
        chromosomes = chromosomes[valid],
        min_segment_length = min_segment_length,
        min_lod = min_lod,
        min_snps = min_snps
      )
      
      # Step 6: Calculate genome-wide IBD proportion
      if (nrow(segments) > 0) {
        # Total length of IBD segments
        total_ibd_length <- sum(segments$length_bp)
        
        # Genome length covered by valid SNPs
        genome_length <- max(positions[valid]) - min(positions[valid])
        
        # IBD proportion
        ibd_proportion <- total_ibd_length / genome_length
        
        # Store in matrix
        ibd_matrix[i, j] <- ibd_proportion
        ibd_matrix[j, i] <- ibd_proportion
        
        # Store segments
        for (k in 1:nrow(segments)) {
          segment_count <- segment_count + 1
          ibd_segments[[segment_count]] <- list(
            sample1 = sample_ids[i],
            sample2 = sample_ids[j],
            chromosome = segments$chromosome[k],
            start_pos = segments$start_pos[k],
            end_pos = segments$end_pos[k],
            length_bp = segments$length_bp[k],
            n_snps = segments$n_snps[k],
            mean_lod = segments$mean_lod[k],
            ibd_proportion = ibd_proportion
          )
        }
        
        # Store quality metrics
        quality_metrics <- rbind(quality_metrics, data.frame(
          sample1 = sample_ids[i],
          sample2 = sample_ids[j],
          n_valid_snps = n_valid,
          missing_rate = missing_rate,
          mean_lod = mean(lod_scores, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  cat(sprintf("\nCompleted! Found %d IBD segments in %d sample pairs\n",
              segment_count, nrow(quality_metrics)))
  
  # Convert segments list to data frame
  if (segment_count > 0) {
    ibd_segments_df <- do.call(rbind, lapply(ibd_segments, as.data.frame))
  } else {
    ibd_segments_df <- data.frame()
  }
  
  return(list(
    ibd_matrix = ibd_matrix,
    ibd_segments = ibd_segments_df,
    quality_metrics = quality_metrics,
    parameters = list(
      min_segment_length = min_segment_length,
      min_lod = min_lod,
      min_snps = min_snps,
      error_rate = error_rate,
      min_maf = min_maf,
      n_snps_used = nrow(genotypes)
    )
  ))
}


#' Calculate LOD scores for IBD at each SNP position
#'
#' @param gt1 Genotype vector for sample 1
#' @param gt2 Genotype vector for sample 2
#' @param af Allele frequencies
#' @param error_rate Genotyping error rate
#' @return Vector of LOD scores

calculate_lod_scores <- function(gt1, gt2, af, error_rate) {
  
  # Check if genotypes match (IBS)
  ibs <- (gt1 == gt2)
  
  # For matching genotypes, calculate LOD score
  # LOD = log10(P(data|IBD) / P(data|not IBD))
  
  # P(match | IBD) = 1 - error_rate (should match if truly IBD)
  # P(match | not IBD) = p^2 + (1-p)^2 for homozygous match by chance
  #                     + 2*p*(1-p) for heterozygous match
  
  lod_scores <- numeric(length(gt1))
  
  for (i in seq_along(gt1)) {
    if (ibs[i]) {
      # Probability of match if IBD
      p_match_ibd <- 1 - error_rate
      
      # Probability of match by chance (depends on genotype)
      # For haploid or homozygous calls (0 or 2):
      if (gt1[i] %in% c(0, 2)) {
        # Homozygous match probability
        p_match_no_ibd <- af[i]^2 + (1 - af[i])^2
      } else {
        # Heterozygous match probability
        p_match_no_ibd <- 2 * af[i] * (1 - af[i])
      }
      
      # Avoid division by zero
      p_match_no_ibd <- max(p_match_no_ibd, 1e-10)
      
      # LOD score
      lod_scores[i] <- log10(p_match_ibd / p_match_no_ibd)
      
    } else {
      # Mismatch - evidence against IBD
      p_mismatch_ibd <- error_rate
      p_mismatch_no_ibd <- 1 - (af[i]^2 + (1 - af[i])^2)
      p_mismatch_no_ibd <- max(p_mismatch_no_ibd, 1e-10)
      
      lod_scores[i] <- log10(p_mismatch_ibd / p_mismatch_no_ibd)
    }
  }
  
  return(lod_scores)
}


#' Detect IBD segments using Hidden Markov Model approach
#'
#' @param lod_scores Vector of LOD scores at each SNP
#' @param positions Physical positions of SNPs
#' @param chromosomes Chromosome identifiers
#' @param min_segment_length Minimum segment length in bp
#' @param min_lod Minimum mean LOD score for segment
#' @param min_snps Minimum number of SNPs
#' @return Data frame of detected IBD segments

detect_ibd_segments_hmm <- function(lod_scores, positions, chromosomes,
                                   min_segment_length, min_lod, min_snps) {
  
  # Simple HMM: classify each SNP as IBD or not based on LOD score
  # Use a threshold approach (simplified HMM)
  
  # State: IBD if LOD > 0, not IBD if LOD <= 0
  ibd_state <- lod_scores > 0
  
  # Apply smoothing: require consecutive SNPs
  # Use a sliding window to smooth noise
  window_size <- 5
  if (length(ibd_state) >= window_size) {
    smoothed_state <- stats::filter(
      as.numeric(ibd_state), 
      rep(1/window_size, window_size), 
      sides = 2
    )
    # Classify as IBD if >50% of window is IBD
    ibd_state <- smoothed_state > 0.5
    ibd_state[is.na(ibd_state)] <- FALSE
  }
  
  # Find runs of IBD state by chromosome
  segments <- data.frame()
  
  for (chrom in unique(chromosomes)) {
    chrom_idx <- which(chromosomes == chrom)
    
    if (length(chrom_idx) < min_snps) next
    
    # Get IBD state for this chromosome
    chrom_ibd <- ibd_state[chrom_idx]
    chrom_pos <- positions[chrom_idx]
    chrom_lod <- lod_scores[chrom_idx]
    
    # Find runs of TRUE (IBD segments)
    rle_result <- rle(chrom_ibd)
    
    # Extract IBD segments
    end_positions <- cumsum(rle_result$lengths)
    start_positions <- c(1, end_positions[-length(end_positions)] + 1)
    
    for (k in which(rle_result$values)) {
      start_idx <- start_positions[k]
      end_idx <- end_positions[k]
      
      # Segment properties
      seg_start <- chrom_pos[start_idx]
      seg_end <- chrom_pos[end_idx]
      seg_length <- seg_end - seg_start
      seg_n_snps <- end_idx - start_idx + 1
      seg_mean_lod <- mean(chrom_lod[start_idx:end_idx], na.rm = TRUE)
      
      # Filter by criteria
      if (seg_length >= min_segment_length && 
          seg_n_snps >= min_snps && 
          seg_mean_lod >= min_lod) {
        
        segments <- rbind(segments, data.frame(
          chromosome = chrom,
          start_pos = seg_start,
          end_pos = seg_end,
          length_bp = seg_length,
          n_snps = seg_n_snps,
          mean_lod = seg_mean_lod
        ))
      }
    }
  }
  
  return(segments)
}


#' Visualize IBD matrix as heatmap
#'
#' @param ibd_results Output from estimate_ibd_robust
#' @export

plot_ibd_matrix <- function(ibd_results) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("pheatmap package needed for visualization. Install with: install.packages('pheatmap')")
    return(invisible(NULL))
  }
  
  pheatmap::pheatmap(
    ibd_results$ibd_matrix,
    color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
    main = "Pairwise IBD Proportions",
    display_numbers = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    breaks = seq(0, 1, length.out = 101)
  )
}


#' Summary statistics for IBD analysis
#'
#' @param ibd_results Output from estimate_ibd_robust
#' @export

summarize_ibd <- function(ibd_results) {
  cat("=== IBD Analysis Summary ===\n\n")
  
  cat("Parameters:\n")
  cat(sprintf("  Min segment length: %.2f Mb\n", 
              ibd_results$parameters$min_segment_length / 1e6))
  cat(sprintf("  Min LOD score: %.2f\n", ibd_results$parameters$min_lod))
  cat(sprintf("  Min SNPs per segment: %d\n", ibd_results$parameters$min_snps))
  cat(sprintf("  Genotyping error rate: %.4f\n", ibd_results$parameters$error_rate))
  cat(sprintf("  SNPs used: %d\n\n", ibd_results$parameters$n_snps_used))
  
  cat("Results:\n")
  cat(sprintf("  Sample pairs analyzed: %d\n", nrow(ibd_results$quality_metrics)))
  cat(sprintf("  IBD segments detected: %d\n", nrow(ibd_results$ibd_segments)))
  
  if (nrow(ibd_results$ibd_segments) > 0) {
    cat(sprintf("  Mean segment length: %.2f Mb\n", 
                mean(ibd_results$ibd_segments$length_bp) / 1e6))
    cat(sprintf("  Mean LOD score: %.2f\n", 
                mean(ibd_results$ibd_segments$mean_lod)))
    cat(sprintf("  Mean IBD proportion: %.4f\n", 
                mean(ibd_results$ibd_segments$ibd_proportion)))
  }
  
  cat("\nQuality metrics:\n")
  cat(sprintf("  Mean valid SNPs per pair: %.0f\n", 
              mean(ibd_results$quality_metrics$n_valid_snps)))
  cat(sprintf("  Mean missing rate: %.4f\n", 
              mean(ibd_results$quality_metrics$missing_rate)))
}
