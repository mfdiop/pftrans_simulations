

# -------------------------------------------------
#               Helper functions#
# ------------------------------------------------
# Calculate allele frequencies from genotype vector
# x: numeric vector of genotypes (0,1)
freq <- function(x){

    n <- length(x)
    ones  <- sum(1*(x==1), na.rm = T) / n 
    zeros <- sum(1*(x==0), na.rm = T) / n 
  
    return (data.frame(ref = zeros, alt = ones, row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = default.stringsAsFactors()))
}

# ---------------------------------------------
# Function to create symmetric matrix from pairwise data
create_symmetric_matrix <- function(df, id_col1 = "Id1", id_col2 = "Id2", 
                                    value_col = "true_ibd_prop") {
  
  # Get all unique IDs
  all_ids <- unique(c(df[[id_col1]], df[[id_col2]]))
  all_ids <- sort(all_ids)
  
  # Create an empty matrix
  n <- length(all_ids)
  mat <- matrix(0, nrow = n, ncol = n)
  rownames(mat) <- all_ids
  colnames(mat) <- all_ids
  
  # Fill in the matrix (both upper and lower triangles)
  for (i in 1:nrow(df)) {
    id1 <- as.character(df[[id_col1]][i])
    id2 <- as.character(df[[id_col2]][i])
    value <- df[[value_col]][i]
    
    mat[id1, id2] <- value
    mat[id2, id1] <- value  # Make it symmetric
  }
  
  # Set diagonal to 1 (or genome_length if using total_ibd_bp)
  diag(mat) <- ifelse(value_col == "true_ibd_prop", 1, genome_length)
  
  return(mat)
}


symmetric_matrix <- function(ibd_summary){
  # Create symmetric data by duplicating rows with swapped IDs
  ibd_symmetric <- ibd_summary %>%
    # Keep original pairs
    select(Id1, Id2, true_ibd_prop) %>%
    # Add reversed pairs
    bind_rows(
      ibd_summary %>%
        select(Id1 = Id2, Id2 = Id1, true_ibd_prop)
    ) %>%
    # Add diagonal (self-pairs)
    bind_rows(
      data.frame(
        Id1 = unique(c(ibd_summary$Id1, ibd_summary$Id2)),
        Id2 = unique(c(ibd_summary$Id1, ibd_summary$Id2)),
        true_ibd_prop = 1
      )
    ) %>%
    distinct()
  
  # Convert to matrix
  ibd_matrix <- ibd_symmetric %>%
    pivot_wider(names_from = Id2, values_from = true_ibd_prop, values_fill = 0) %>%
    column_to_rownames("Id1") %>%
    as.matrix()
  
  # Ensure rows and columns are in the same order
  ibd_matrix <- ibd_matrix[order(as.numeric(rownames(ibd_matrix))), 
                           order(as.numeric(colnames(ibd_matrix)))]
}

# ---------------------------------------------
# Read and standardize IBD calls from a tskibd/hmmIBD output file
read_ibd <- function(fn) {
  # Read IBD file and standardize column names
  df <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  rename_if <- function(idx, newname) {
    if (!is.na(idx)) names(df)[idx] <<- newname
  }
  
  map <- list(
    id1    = match("id1",    tolower(names(df))),
    id2    = match("id2",    tolower(names(df))),
    chrom  = match("chrom",  tolower(names(df))),
    start  = match("start",  tolower(names(df))),
    end    = match("end",    tolower(names(df))),
    tmrca  = match("tmrca",  tolower(names(df))),
    cm     = match("cm",     tolower(names(df)))
  )
  
  #   rename_if <- function(idx, to) if (!is.na(idx)) names(df)[idx] <<- to
  rename_if(map$id1, "Id1"); rename_if(map$id2, "Id2")
  rename_if(map$chrom, "Chrom"); rename_if(map$start, "Start"); rename_if(map$end, "End")
  
  if (!is.na(map$tmrca)) names(df)[map$tmrca] <- "Tmrca"
  if (!is.na(map$cm))    names(df)[map$cm]    <- "Cm"
  
  # If still missing essentials, try positional 1:4
  req <- c("Id1","Id2","Start","End")
  if (!all(req %in% names(df)) && ncol(df) >= 4) names(df)[1:4] <- req
  
  if (!all(req %in% names(df))) stop("Cannot standardize IBD columns in ", fn)
  df[df$End > df$Start, , drop=FALSE]
}

# ---------------------------------------------
# Compute true IBD per chromosome using tskibd 
# --------------------------------------------
run_tskibd <- function(tree, chr) {
  # Run tskibd <chr> <bp_per_cm> <sample_window_bp> <mincm> <tree>
  
  bp_per_cm <- 1e6  # 1 cM per Mb
  sample_window_bp <- 1e7  # 10 Mb window
  mincm <- 0.0  # no minimum cM

  tmp <- sprintf("tskibd_%s_%s", chr, as.integer(runif(1, 1, 1e9)))
  dir.create(tmp, showWarnings=FALSE)

  owd <- getwd(); on.exit({ setwd(owd); unlink(tmp, recursive=TRUE) }, add=TRUE)
  setwd(tmp)

  cmd <- c(as.character(chr), as.character(bp_per_cm),
           as.character(sample_window_bp), as.character(mincm), normalizePath(file.path(owd, tree)))

  out <- tryCatch(system2(tskibd_bin, args=cmd, stdout=TRUE, stderr=TRUE),
                  warning=function(w) w, error=function(e) e)

  ibd_path <- file.path(getwd(), sprintf("%s.ibd", chr))

  if (!file.exists(ibd_path)) {
    stop("tskibd did not produce ", ibd_path, "\nOutput:\n", paste(out, collapse="\n"))
  }

  read_ibd(ibd_path)
}

# ------------------- main loop -------------------
# run_tskibd_all: run tskibd on all chromosomes and write per-chromosome truth files
# Args:
# trees: vector of tree file paths, one per chromosome
# chrnos: vector of chromosome numbers, one per chromosome
# outdir: output directory for truth files
# max_tmrca: optional vector of max TMRCA thresholds for filtering
# Returns: nothing; writes files to outdir
# -----------------------------------------------

run_tskibd_all <- function(trees, chrnos, outdir, max_tmrca = 10){
  all_list <- list()
  
  # Compute true IBD across chromosome using tskibd
  for (i in seq_along(chrnos)) {
    chr <- chrnos[i]; tree <- trees[i]
    message(sprintf("[chr %d] tskibd on %s ...", chr, tree))

    df <- run_tskibd(tree, chr)
    
    # Add Chrom if missing; keep standard columns
    if (!"Chrom" %in% names(df)) df$Chrom <- chr
    keep <- intersect(c("Id1","Id2","Chrom","Start","End","Tmrca","Cm"), names(df))
    df <- df[, keep, drop=FALSE]
    
    # Write per-chromosome truth
    out_fn <- file.path(outdir, sprintf("%d.true_ibd.tsv", chr))
    write_tsv(df, out_fn)
    message("  wrote: ", out_fn)
    
    # Optional: TMRCA-filtered truth
    if (length(max_tmrca) > 0 && "Tmrca" %in% names(df)) {
      for (m in max_tmrca) {
        sub <- df[df$Tmrca <= m, , drop=FALSE]
        of <- file.path(outdir, sprintf("%d.true_ibd_maxtmrca_%d.tsv", chr, m))
        write_tsv(sub, of)
        message("  wrote: ", of)
      }
    }
    all_list[[i]] <- df
  }
  
  # Merge across chromosomes
  all_df <- do.call(rbind, all_list)
  merged_fn <- file.path(outdir, "true_ibd.all.tsv")
  write_tsv(all_df, merged_fn)
}

# ----- End of script -----

# ---------------------------------------------------
# Run Identity-by-state (IBS) analysis from the VCF
# ---------------------------------------------------
# Example usage:
#      # vcf_file <- "simulated_data.vcf"
#      #      # ibs_matrix <- run_ibs(vcf_file)
#      #      # View IBS matrix
#      #      # head(ibs_matrix)
#      ---------------------------------------------------

# Function to convert VCF genotypes to 0/1/2
convert_genotypes <- function(gt_matrix) {
  result <- matrix(0, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
  
  for (i in 1:nrow(gt_matrix)) {
    for (j in 1:ncol(gt_matrix)) {
      gt <- gt_matrix[i, j]
      
      # Extract just the genotype part (before any :)
      gt_clean <- strsplit(gt, ":")[[1]][1]
      
      # Convert to numeric count of alternate alleles
      if (grepl("\\|", gt_clean)) {  # Phased
        alleles <- as.numeric(strsplit(gt_clean, "\\|")[[1]])
      } else if (grepl("/", gt_clean)) {  # Unphased
        alleles <- as.numeric(strsplit(gt_clean, "/")[[1]])
      } else {
        alleles <- as.numeric(gt_clean)
      }
      
      result[i, j] <- sum(alleles)  # Count of alternate alleles (0,1,2)
    }
  }
  return(result)
}


run_ibs <- function(vcf_file, SNPRelate = FALSE) {
  
  # Read VCF file
  # vcf <- vcfR::read.vcfR(vcf_file)
  vcf <- read.table(vcf_file, sep = "\t", comment.char = "#",
                    check.names = FALSE, header = FALSE,
                    stringsAsFactors = FALSE)

  names(vcf) <- c("#CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT",
                  paste0("tsk", seq_along(10:ncol(vcf))))
  
  # Extract genotype matrix (convert to numeric alternative allele counts)
  # gt <- vcfR::extract.gt(vcf, element = "GT")
  gt <- vcf[, 10:ncol(vcf)]
  
  # Convert GT to numeric alternative allele counts
  # genotypes <- matrix(as.numeric(gsub("\\|", "", gsub("/", "", gt))),
  #                     nrow = nrow(gt), ncol = ncol(gt))
  
  # Convert genotypes
  genotypes <- convert_genotypes(gt)
  
  # Assuming genotypes is a matrix 
  # where rows are variants and columns are samples
  n_variants <- nrow(genotypes)
  n_samples <- ncol(genotypes)

  cat("Processing", n_samples, "individuals...\n")
  cat("Expected unique pairs:", (n_samples * (n_samples - 1)) / 2, "\n")
  
  # Calculate pairwise matching (IBS)
  # OPTION 1: Calculate only unique pairs (recommended)
  unique_pairs <- data.frame()
  
  for (i in 1:(n_samples - 1)) {
    for (j in (i + 1):n_samples) {  # Only j > i to avoid duplicates
      matches <- sum(genotypes[, i] == genotypes[, j], na.rm = TRUE)
      ibs_value <- round(matches / n_variants, 4)
      
      unique_pairs <- rbind(unique_pairs, 
                           data.frame(id1 = colnames(gt)[i],
                                      id2 = colnames(gt)[j],
                                      IBS = ibs_value))
    }
  }
  
  cat("Actual unique pairs calculated:", nrow(unique_pairs), "\n")
  
  # ----- Alternative manual IBS calculation -----
  # More efficient vectorized approach
  ibs_matrix1 <- matrix(0, nrow = n_samples, ncol = n_samples)
  for (i in 1:n_samples) {
    matches <- colSums(genotypes == genotypes[, i], na.rm = TRUE)
    ibs_matrix1[i, ] <- round(matches / n_variants, 2)
  }
  
  # Add sample names
  colnames(ibs_matrix1) = rownames(ibs_matrix1) <- colnames(gt)

  # ==========================================================
  # Calculate IBS for unique pairs only
  unique_pairs <- data.frame()
  
  for (i in 1:(n_samples - 1)) {
    # Vectorized comparison for all j > i
    sample_i <- geno_mat[, i]
    
    # Compare sample i with all samples j > i
    for (j in (i + 1):n_samples) {
      matches <- sum(sample_i == geno_mat[, j], na.rm = TRUE)
      ibs_value <- matches / nrow(geno_mat)
      
      unique_pairs <- rbind(unique_pairs,
                           data.frame(id1 = colnames(gt)[i],
                                      id2 = colnames(gt)[j],
                                      IBS = round(ibs_value, 4)))
    }
  }
  
  # More efficient approach with SNPRelate:
  if(SNPRelate){
    # Load required libraries
    library(SNPRelate)
    
    # Convert VCF to GDS format (efficient for large files)
    gds_file <- gsub(".vcf.gz", ".gds", vcf_file)
    snpgdsVCF2GDS(vcf_file, gds_file, ) 
    
    # Open GDS file
    genofile <- snpgdsOpen(gds_file)
    
    # Calculate IBS matrix
    ibs <- snpgdsIBS(genofile, num.thread = 2, autosome.only = FALSE)
    ibs_matrix2 <- round(ibs$ibs, 2)
    
    # Close file
    snpgdsClose(genofile)
    
    rm(gds_file)
  }
  
  return(unique_pairs)
}

# run_ibs("1_1.vcf")
# ----- End of script -----

# ---------------------------------------------------
#      Run IBD methods to infer IBD estimates from VCF
# ---------------------------------------------------
# Run hmmIBD to infer IBD segments from the VCF file 
# # Example usage:
#       # vcf_file <- "simulated_data.vcf"
#       # out_prefix <- "hmmibd_output"
#       # hmmibd_results <- run_hmmibd(vcf_file, out_prefix, min_cm=0.5)
#       # View results
#       # head(hmmibd_results)
# --------------------------------------------------------
run_hmmibd <- function(vcf_file, out_prefix, min_cm = 0.0) {
  cmd <- c(hmmibd_bin, "--vcf", vcf_file, "--out", out_prefix, "--min_cm", as.character(min_cm))
  out <- tryCatch(system2(cmd[1], args=cmd[-1], stdout=TRUE, stderr=TRUE), warning=function(w) w, error=function(e) e)
  
  hmmibd_ibd_path <- paste0(out_prefix, ".ibd")
  if (!file.exists(hmmibd_ibd_path)) {
    stop("hmmIBD did not produce ", hmmibd_ibd_path, "\nOutput:\n", paste(out, collapse="\n"))
  }
  
  read_ibd(hmmibd_ibd_path)          
}

# ----- End of script -----

#------------------
# RUN MLE IBD
#------------------
run_mle <- function(vcf_file){
  # MIPanalyzer::inbreeding_mle function to estimate pairwise IBD from VCF genotype data
  # Args:
  #   vcf: vcfR object
  # Returns:
  #   ibd.matrix: pairwise IBD matrix
  # Example usage:
  #    vcf_file <- "simulated_data.vcf" 
  #    vcf <- vcfR::read.vcfR(vcf_file)
  #    ibd_matrix <- run_mle(vcf)

  ibd.matrix <- MIPanalyzer::inbreeding_mle(x = vcf,
                                            f = seq(0.01, 0.99, 0.01),
                                            ignore_het = FALSE,
                                            report_progress = FALSE)
  
  return(ibd.matrix)
}


#------------------
# ESTIMATE HMM IBD
#------------------
run_hmmibdr <- function(vcf_file, recomb = 6.666667e-07,
                        gt_file = "gt_hmmIBD.txt",
                        af_file = "af_hmmIBD.txt"){
  
  # load required libraries
  require(vcfR)
  require(hmmibdr)
  require(tidyverse)

  #----------
  # Part 1
  #----------
  # Read VCF file
  vcf <- read.vcfR(vcf_file)
  
  # Extract genotype matrix (convert to numeric alternative allele counts)
  gt <- extract.gt(vcf, element = "GT")
  
  # # Split genotypes and sum alleles
  # genotypes <- apply(gt, c(1, 2), function(x) {
  #   if (is.na(x)) return(NA)
  #   # Split by | or / and sum the alleles
  #   alleles <- as.numeric(unlist(strsplit(x, "[|/]")))
  #   sum(alleles)
  # })
  # 
  # # Faster vectorized alternative:
  # # Replace separators with empty string and convert
  # gt_clean <- gsub("[|/]", "", gt)
  # 
  # # Sum the digits (works because 0+0=0, 0+1=1, 1+1=2)
  # genotypes <- matrix(
  #   sapply(gt_clean, function(x) {
  #     if (is.na(x) || x == ".") return(NA)
  #     sum(as.numeric(strsplit(x, "")[[1]]))
  #   }),
  #   nrow = nrow(gt), 
  #   ncol = ncol(gt)
  # )
  
  # Most efficient approach:
  
  # Convert GT to numeric alternative allele counts
  genotypes <- matrix(
    as.numeric(gsub("0\\|0|0/0", "0", 
                    gsub("0\\|1|1\\|0|0/1|1/0", "1",
                         gsub("1\\|1|1/1", "1", gt)))),
    nrow = nrow(gt),
    ncol = ncol(gt)
  ) %>% as.data.frame()

  names(genotypes) <- colnames(gt)
    
  #------------------------------
  # Part 2: format genotypes
  #------------------------------
  gtmat <- dplyr::bind_cols("chrom" = as.numeric(vcfR::getCHROM(vcf)), 
                            "pos" = vcfR::getPOS(vcf), 
                            genotypes) 
  
  # Save genotype matrix for hmmIBD
  readr::write_tsv(x = gtmat, file = gt_file, col_names = T) 
  
  #-------------------------------------
  # Part 3: Estimate Allele Frequencies
  #-------------------------------------
  
  # Calculate allele frequencies for each variant
  afreqs <- apply(genotypes, 1, freq)
  afreqs <- do.call(rbind, afreqs)

  # Transpose and convert to data frame
  afmat <- dplyr::bind_cols("chrom" = vcfR::getCHROM(vcf),
                            "pos" = vcfR::getPOS(vcf),
                            afreqs) %>%
    dplyr::mutate(chrom = as.numeric(factor(chrom)))
  
  # Save allele frequency file for hmmIBD
  readr::write_tsv(afmat, file = af_file, col_names = F)
  
  #----------------
  # Running hmmIBD
  #----------------
  tf <- tempfile(pattern = "output")
  out <- hmmibdr::hmm_ibd(input_file = gt_file,
                          allele_freqs =  af_file,
                          max_fit_iterations = 20,
                          rec_rate = recomb, # 1e-2 note the small recombo rate relative to what would be expected in malaria
                          output_file = tf)

  # hmmIBD tidy
  ibd_hmm <- tibble::tibble(
      p1 = out$fract$sample1,
      p2 = out$fract$sample2,
      hmm = out$fract$fract_sites_IBD)
  
  return(ibd_hmm)
}

# ----- End of script -----

# Create inputs for isoRelate
isorelate_prep <- function(vcf_file, output_dir){
  # Convert VCF to isoRelate format and save to file
  # Args:
  #   vcf_file: path to VCF file
  #   output_prefix: prefix for output files
  # Returns:
  #   None; writes files to disk
  
  # Load required libraries
  library(tidyverse)
  library(vcfR)
  
  # Read VCF file
  vcf <- read.vcfR(vcf_file)

  gt <- extract.gt(vcf, element = "GT")
  sample_ids <- colnames(gt)

  # Convert GT to numeric alternative allele counts
  genotypes <- matrix(
    as.numeric(gsub("0\\|0|0/0", "1", 
                    gsub("0\\|1|1\\|0|0/1|1/0", "2",
                         gsub("1\\|1|1/1", "2", gt)))),
    nrow = nrow(gt),
    ncol = ncol(gt)
  ) %>% as.data.frame()

  #names(genotypes) <- sample_ids
   
  #......................
  # part 2
  #......................
  gtmat <- dplyr::bind_cols("chrom" = as.numeric(vcfR::getCHROM(vcf)), 
                            "pos" = vcfR::getPOS(vcf), 
                            genotypes)
   
  # Duplicate SNPs to have two alleles
  gtmat <- gtmat[rep(seq_len(nrow(gtmat)), each = 2), -c(1,2)]
   
  # Create PED file format
  ped <- data.frame(family_id = sample_ids, isoate_id = sample_ids,
                    paternal_id = 0, maternal_id = 0, 
                    moi = 1, phenotype = 0, t(gtmat))
   
  # Create MAP file format
  map <- vcf@fix[,1:2] %>%
      tibble::as_tibble() %>%
      type_convert() %>% 
      dplyr::mutate(CHROM = paste0("Pf3D7_0", CHROM, "_v3"),
                    snp_id = paste(CHROM, POS, sep = ":"),
                    pos_cm = POS/17200,
                    pos_bp = as.integer(POS)) %>%
      dplyr::select(-POS) %>% 
      as.data.frame()
      
  # Create list from PED and MAP
  ped.map <- list(ped, map)

  # Save to files
#  saveRDS(ped.map, file = file.path(output_dir, "_isorelate_pedmap.rds"))
   
  return(ped.map)

}

run_isorelate <- function(pedmap){
  # Run isoRelate to infer IBD segments from the VCF file 
  # Args:
  #   vcf_file: path to VCF file
  # Returns:
  #   inferred_ibd_iso: data frame of inferred IBD segments
  
  # Load required libraries
  library(isoRelate)

  # reformat and filter genotypes
  my_genotypes <- getGenotypes(ped.map = pedmap,
                              reference.ped.map = NULL,
                              maf = 0.01,
                              isolate.max.missing = 0.1,
                              snp.max.missing = 0.1,
                              chromosomes = NULL,
                              input.map.distance = "cM",
                              reference.map.distance = "cM")

  # estimate parameters
  my_parameters <- getIBDparameters(ped.genotypes = my_genotypes, 
                                    number.cores = 1)

  head(my_parameters)

  # infer IBD
  my_ibd <- getIBDsegments(ped.genotypes = my_genotypes,
                          parameters = my_parameters, 
                          number.cores = 1, 
                          minimum.snps = 20, 
                          minimum.length.bp = 50000,
                          error = 0.001)

  head(my_ibd)

  # get a summary of IBD segments
  getIBDsummary(ped.genotypes = my_genotypes, 
                ibd.segments = my_ibd)

  
  return(list(my_genotypes = my_genotypes,
              my_parameters = my_parameters,
              my_ibd = my_ibd))
} 

extract_related_pairs <- function (ped.genotypes, ibd.segments, prop = 1) 
{
   pedigree <- ped.genotypes[["pedigree"]]
   genotypes <- ped.genotypes[["genotypes"]]
   
   genome.length <- 0
   for (chr in unique(as.character(genotypes[, "chr"]))) {
      genotypes.0 <- genotypes[genotypes[, "chr"] == chr, ]
      genome.length <- genome.length + max(genotypes.0[, "pos_bp"]) - 
         min(genotypes.0[, "pos_bp"])
   }
   ibd.pairs <- paste(ibd.segments[, "fid1"], 
                      ibd.segments[, "iid1"],
                      ibd.segments[, "fid2"],
                      ibd.segments[, "iid2"], 
                      sep = "/")
   
   isolate.pairs <- isolatePairs(pedigree[, 1], pedigree[, 2])
   isolate.pairs <- paste(isolate.pairs[, 1], 
                          isolate.pairs[, 2], 
                          isolate.pairs[, 3], 
                          isolate.pairs[, 4], sep = "/")
   highly.related <- NULL
   for (i in 1:length(unique(ibd.pairs))) {
      ibd.segments.0 <- ibd.segments[ibd.pairs == unique(ibd.pairs)[i], ]
      # if (sum(ibd.segments.0[, "length_bp"])/genome.length >= prop) 
      #    highly.related <- c(highly.related, unique(ibd.pairs)[i])
      ibd <- sum(ibd.segments.0[, "length_bp"])/genome.length
      if(i == 1){
         highly.related <- tibble(fid1 = unique(ibd.segments.0$fid1),
                                  iid1 = unique(ibd.segments.0$iid1),
                                  fid2 = unique(ibd.segments.0$fid2),
                                  iid2 = unique(ibd.segments.0$iid2),
                                  ibd)
      }
      else{
         highly.related <- rbind(highly.related, 
                                 tibble(fid1 = unique(ibd.segments.0$fid1),
                                        iid1 = unique(ibd.segments.0$iid1),
                                        fid2 = unique(ibd.segments.0$fid2),
                                        iid2 = unique(ibd.segments.0$iid2),
                                        ibd))
      }
      
   }
   
   related.pairs <- highly.related %>% 
      mutate(related.pairs = paste(fid1, iid1, fid2, iid2, sep = "/")) %>% 
      pull(related.pairs)
   
   ibd.interval <- ibd.segments[ibd.pairs %in% related.pairs, ]
   
   if (nrow(ibd.interval) == 0) 
      stop("no pairs share >= ", prop * 100, "% of their genome IBD")
   ibd.melt <- data.frame(from = paste(ibd.interval[, "fid1"], ibd.interval[, "iid1"], sep = "/"), 
                          to = paste(ibd.interval[, "fid2"], ibd.interval[, "iid2"], sep = "/"), 
                          value = ibd.interval[, "length_bp"])
   ibd.melt <- ibd.melt[!duplicated(paste(ibd.melt[, "from"], ibd.melt[, "to"], sep = "/")), ]
   nodes <- unique(c(as.character(ibd.melt[, 1]), as.character(ibd.melt[, 2])))
   
   return(highly.related)
}

# ----- End of script -----

# ---------------------------------------