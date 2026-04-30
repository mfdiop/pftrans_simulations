#!/usr/bin/env Rscript

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
  library(tidyverse)
  library(stringr)
})

have_gg <- requireNamespace("ggplot2", quietly = TRUE)

save_plot_png <- function(plot_obj, filename, width=7, height=5, res=120) {
  if (have_gg && inherits(plot_obj, "ggplot")) {
    ggplot2::ggsave(filename, plot = plot_obj, width = width, height = height, dpi = res* (96/72))
  } else {
    # plot_obj is a function that draws the base plot
    png(filename, width = width, height = height, units = "in", res = res)
    on.exit(dev.off(), add = TRUE)
    plot_obj()
  }
}

# load useful functions
source('wrapper_codes/01_useful_functions.R', chdir = TRUE)

# ----- CLI -----
option_list <- list(
  
  make_option("--tree", type = "character", help = "Required."),
  make_option("--vcf", type = "character", help = "Required."),
  make_option("--true_ibd", type = "character", help = "Required."),
  make_option("--ibs", action = "store_true", default = FALSE,
              help = "Compute IBS estimates"),
  make_option("--hmm", action = "store_true", default = FALSE,
              help = "Compute IBS estimates"),
  make_option("--isorelate", action = "store_true", default = FALSE,
              help = "Compute IBS estimates"),
  make_option("--rec_rate", type = "double",  default = 0.01 / 15000),
  make_option("--outdir", type = "character", default = "post_simulation"),
  make_option("--test", action = "store_true", default = FALSE,
              help = "Use smaller nsam and seqlen for a quick test")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Check if the length of the arguments vector is zero
if (length(opt) == 0) {
  stop("Error: No arguments provided to the script. At least 3 arguments must be provided", call. = FALSE)
} else {
  message("Arguments provided:")
  print(opt)
}


# setup output directory

if (!dir.exists(opt$outdir)) dir.create(opt$outdir)

# ----- Python interop (tskit) -----
py_require(c("pyslim", "msprime", "tskit", "pandas"))
tskit    <- import("tskit", convert = FALSE)

# ----- Begin script -----
# (a) Extract true IBD from .trees files using tskit/tskitibd
# 1) Load trees and extract IBD segments using tskitibd
# ts <- tskit$load("single_pop/single_run/chr1_1.trees")
ts <- tskit$load(opt$tree)

# Initial sample size
n_sam <- py_as_integer(ts$num_samples) # ts$sample_size

# 2) Compute genome length (assuming 1 chromosome of 15,000 bp)
genome_length <- py_as_integer(ts$sequence_length)

# 3) Read IBD data from tskibd output file
# ibd_data <- read_tsv("single_pop/out_true_ibd/chr1_1.ibd")
true_ibd <- read_tsv(opt$true_ibd)

# Aggregate total shared length per pair of individuals per chromosome
ibd_summary <- true_ibd %>%
  mutate(segment_length = End - Start) %>%
  group_by(Id1, Id2) %>%
  summarise(total_ibd_bp = sum(segment_length), 
            total_ibd_prop = total_ibd_bp / genome_length,
            n_segments = n(),
            max_segment_bp = max(segment_length),
            max_segment_prop = max_segment_bp / genome_length,
            mean_segment = mean(segment_length),
            min_tmrca = min(Tmrca, na.rm = TRUE),
            .groups = 'drop') %>%
  ungroup() %>% 
  mutate(true_ibd_prop = total_ibd_bp / genome_length) %>%   # Normalize by genome length
  mutate(pair = paste(pmin(Id1, Id2), pmax(Id1, Id2), sep = "_"))


# Save summary to file
write_tsv(ibd_summary, file.path(opt$outdir, "true_ibd_summary.tsv"))

# 4) Create a true IBD proportion symmetric matrix
ibd_matrix <- create_symmetric_matrix(ibd_summary, 
                                      id_col1 = "Id1", 
                                      id_col2 = "Id2",
                                      value_col = "true_ibd_prop")

# # Or for other metrics:
# ibd_bp_matrix <- create_symmetric_matrix(ibd_summary, 
#                                          value_col = "total_ibd_bp")
# 
# ibd_nseg_matrix <- create_symmetric_matrix(ibd_summary, 
#                                            value_col = "n_segments")

ibd_long <- as.data.frame(as.table(ibd_matrix))
colnames(ibd_long) <- c("Id1", "Id2", "IBD")

write_tsv(ibd_long, file.path(opt$outdir, "true_ibd_summary_long.tsv"))

# # Visualize the matrix as a heatmap
# library(pheatmap)
# 
# pheatmap(ibd_matrix,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          color = colorRampPalette(c("white", "darkred", "darkblue"))(50),
#          main = "Pairwise IBD Proportions")
# 
# # Or using ggplot2
# 
# ggplot(ibd_long, aes(x = Id1, y = Id2, fill = IBD_prop)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "darkblue") +
#   theme_minimal() +
#   labs(title = "Pairwise IBD Proportions",
#        fill = "IBD Proportion") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


# ----- End of true IBD summary -----
# ---------------------------------------

# Compute inferred methods and compare to true IBD from simulated data

# ----- IBS computation -----
# (a) Compute IBS matrix from VCF
# vcf_file <- "single_pop/single_run/chr1_1.vcf.gz"
vcf_file <- opt$vcf
ibs_matrix <- run_ibs(vcf_file)


# Save IBS matrix to file
save(ibs_matrix, file = file.path(opt$outdir, "ibs_matrix.RData"))
saveRDS(ibs_matrix, file = file.path(opt$outdir, "ibs_matrix.rds"))

# Transform the wide IBS matrix to long format
# Or using ggplot2
ibs_long <- as.data.frame(as.table(ibs_matrix))
colnames(ibs_long) <- c("Id1", "Id2", "IBS")

write_tsv(ibs_long, file.path(opt$outdir, "inferred_ibs.tsv"))

# ----- IBD inference and comparison -----
# (b) Compute inferred IBD matrices using different tools
# Run hmmIBD
inferred_ibd_hmm <- run_hmmibdr(vcf_file, recomb = opt$rec_rate)
colnames(inferred_ibd_hmm) <- c("Id1", "Id2", "hmm")
write_tsv(inferred_ibd_hmm, file.path(opt$outdir, "inferred_ibd_hmm.tsv"))

# Run isoRelate
# Prepare data for isoRelate
pedmap <- isorelate_prep(vcf_file,  output_dir = opt$outdir)  

# Run isoRelate
results <- run_isorelate(pedmap)

write_rds(results, file = "single_pop/out_true_ibd/inferred_ibd_iso.rds")

# Extract related pairs
inferred_ibd_iso <- extract_related_pairs(ped.genotypes = resuts$my_genotypes, 
                                       ibd.segments = resuts$my_ibd, 
                                       prop = 0.9)

names(inferred_ibd_iso) <- c("fid1", "Id1", "fid2", "Id2", "ibd")

write_tsv(inferred_ibd_iso, file.path(opt$outdir, "inferred_ibd_iso.tsv"))
# Run other IBD tools similarly and save their outputs...








# (c) Compare “inferred” vs “true” matrices
# 
# You’ll have three comparable matrices:
#   true_ibd_matrix (from tskitIBD)
#   inferred_ibd_matrix (from IBD tool)
#   ibs_matrix (from allele matching)
#   

# compute:
  
# Correlation (Spearman or Pearson) between true and inferred IBD matrices.


# Mean absolute deviation or RMSE of inferred vs true shared proportions.

# Rank-based recovery accuracy: fraction of true top-k related pairs that appear among top-k inferred ones (important for identifying transmission pairs).




