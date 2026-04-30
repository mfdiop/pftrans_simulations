#!/usr/bin/env Rscript

# Spatial population simulation analysis pipeline
# - runs SLiM,pyslim, msprime and tskit for spatial simulation
# - uses reticulate to bind Python and R

suppressPackageStartupMessages({
  library(optparse)
  library(reticulate)
  library(tidyverse)
  
  # bind python
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
  py_to_r <- function(x) {
    reticulate::py_to_r(x)
  }
  
  # --- Bind reticulate to the correct Python BEFORE any import() calls ---
  # Priority: RETICULATE_PYTHON > CONDA_DEFAULT_ENV > "slim-msprime"
  py_available <- reticulate::py_config()
  message("Using Python: ", py_available$python)
  
  try({
    py <- Sys.getenv("RETICULATE_PYTHON", "")
    if (nzchar(py)) use_python(py, required = TRUE)
    else{
      env <- Sys.getenv("CONDA_DEFAULT_ENV", "")
      if (nzchar(env)) {
        use_condaenv(env, required = FALSE)
      } else {
        use_condaenv("slim-msprime", required = FALSE)
      }
    }
  }, silent = TRUE)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  
})

option_list <- list(
  make_option("--slim_bin", type="character", default="slim.exe"),
  make_option("--slim_script", type="character", default="wrapper_codes/spatial_pop_improved.slim"),
  make_option("--outdir", type="character", default="spatial_analysis"),
  make_option("--chrno", type="integer", default=1),
  make_option("--genome_set_id", type="integer", default=1),
  make_option("--r", type="double", default=1e-8),
  make_option("--mu", type="double", default=1e-8),
  make_option("--ancestral_Ne", type="integer", default=1000),
  make_option("--remember_gen", type="integer", default=1000),
  make_option("--seed", type="integer", default=42)
)

opt <- parse_args(OptionParser(option_list=option_list))

# validate options and prepare output dir structure
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
run_prefix <- file.path(opt$outdir, sprintf("%d_%d", opt$genome_set_id, opt$chrno))
dir.create(run_prefix, recursive = TRUE, showWarnings = FALSE)

# build SLiM command: pass constants via -d
slim_bin <- Sys.which(opt$slim_bin)
if (nchar(slim_bin) == 0) stop("SLiM binary not found on PATH (use --slim_bin)")
if (!file.exists(opt$slim_script)) stop("SLiM script not found: ", opt$slim_script)

# set -d args; ensure numeric values are formatted safely
slim_params <- list(
  # L = opt$seqlen,
  # selpos = opt$selpos,
  # num_origins = opt$num_origins,
  # N = opt$N,
  # h = opt$h,
  # s = opt$s,
  # g_sel_start = opt$g_sel_start,
  r = opt$r
  # outid = opt$chrno
  # max_restart = 100,
  # sim_relatedness = opt$sim_relatedness,
  # N0 = opt$N0,
  # g_ne_change_start = opt$g_ne_change_start
)

# dargs <- c(
#   sprintf("-d OUTID=%d", opt$chrno),
#   sprintf("-d REMEMBER_T_GEN=%d", opt$remember_gen),
#   sprintf("-d TOTAL_GEN=%d", 2 * opt$remember_gen), # keep same structure as script default
#   sprintf("-d W=35") # adjust if needed
# )

# Construct command line args for SLiM
slim_args <- paste(unlist(Map(function(k,v) sprintf("-d %s=%s", k, as.character(v)),
                              names(slim_params), slim_params)), collapse = " ")
# set seed if not provided
# seed_arg  <- if (!is.null(slim_seed)) sprintf("-seed %s", slim_seed) else ""

cmd_args <- sprintf("%s %s", slim_args, shQuote(normalizePath(opt$slim_script)))
message(">> simulate chrom with id ", opt$chrno) 
message(">> ", cmd_args)

# Write the command to a file for reproducibility
cmd_file <- file.path(run_prefix, sprintf("slim_cmd_%d.txt", opt$chrno))
writeLines(paste0(slim_bin, " ", cmd_args), cmd_file)

stdout_f <- file.path(run_prefix, "slim_stdout.txt")
stderr_f <- file.path(run_prefix, "slim_stderr.txt")

message("Running SLiM...")
res <- system2(slim_bin, args = cmd_args, stdout = stdout_f, stderr = stderr_f, wait = TRUE)

if (res != 0) stop(sprintf("SLiM failed (exit %s); see %s and %s", res, stdout_f, stderr_f))
message("SLiM completed; stdout in ", stdout_f)

# --------------------------------------------------
# ----- Python interop (tskit/msprime/pyslim) -----
#  --------------------------------------------------
# install required Python packages if needed
py_require(c("tskit", "pyslim", "msprime", "pandas"), action = "add")

# import Python modules
np       <- import("numpy", convert = TRUE)
tskit    <- import("tskit", convert = FALSE)   # keep as Py objects
pyslim   <- import("pyslim", convert = FALSE)
msprime  <- import("msprime", convert = FALSE)
builtins <- import_builtins()

# ensure helper module is importable and dir is on sys.path
reticulate::py_run_string("import sys")

# add current dir to front of sys.path
reticulate::py_run_string(sprintf("sys.path.insert(0, '%s')", normalizePath(".", winslash="/"))) 

# import helpers
helpers <- import_from_path(module = "tskit_spatial_helpers", 
                            path = "wrapper_codes")

# If you want to import a python file as a package
# Step 1 – Add an __init__.py file
# Step 2 – Import correctly in R (no .py extension):
# helpers <- import(module = "wrapper_codes.tskit_spatial_helpers")

py_run_string("print('Python environment active')")

# ------------------------------------------------------------------------------
# find trees output (script writes spatial_sim_out_<OUTID>.trees)
trees_path <- file.path(getwd(), sprintf("spatial_sim_out_%d.trees", opt$chrno))

if (!file.exists(trees_path)) {
  stop("Trees file not found: ", trees_path)
}

# Read SLiM stdout for SAMPLE_MAP and PEDIGREE lines
stdout_lines <- readLines(stdout_f, warn = FALSE)
sample_map_df <- helpers$parse_sample_map_lines(stdout_lines) # returns pandas DataFrame object

# convert to R data.frame
sample_map_r <- py_to_r(sample_map_df)

write_tsv(sample_map_r, file.path(run_prefix, "sample_map.tsv"))

# recapitate and mutate via Python helper
recap_out <- file.path(run_prefix, "spatial_sim.recap.trees")

mts <- helpers$recapitate_and_mutate(trees_path, 
                                     recap_out_path = recap_out,
                                     recombination_rate = opt$r,
                                     ancestral_Ne = opt$ancestral_Ne,
                                     mu = opt$mu,
                                     random_seed = opt$seed)

py_run_string("print('Recapitation and mutation done')")
py_run_string(sprintf("print('Recap output: %s')", recap_out))
# py_run_string("print(The tree sequence now has {mts.num_trees} trees, and {mts.num_sites} sites and {mts.num_mutations} mutations.)")

# mts is a Python tskit tree sequence object accessible via reticulate
# For convenience, assign to R variable via pyobject
ts_obj <- mts

# Now sample groups using helper: modern (time=0) and ancient (time = REMEMBER_T_GEN)
groups <- helpers$sample_groups(ts_obj, alive_time = 0, ancient_time = opt$remember_gen, W = 35, w = 5, rng_seed = opt$seed)

# convert groups into R lists for quick inspect
groups_r <- py_to_r(groups)

# write groups summary
groups_summary <- lapply(names(groups_r), function(k) length(groups_r[[k]]))
writeLines(paste(names(groups_summary), " group has ", unlist(groups_summary), sep=":\t"),
           con = file.path(run_prefix, "groups_summary.txt"))

# Optional pre-check
py_run_string("
# check that all IDs are integers
def check_integer_ids(groups):
   for k, v in groups.items():
     bad = [x for x in v if not float(x).is_integer()]
   if bad:
     print(f'⚠️ Non-integer IDs in group {k}: {bad[:5]}' + ('...' if len(bad) > 5 else ''))
   
 ")

# run the check
main <- import_main()
main$check_integer_ids(groups)

# compute group-level divergence matrix
group_div_df <- helpers$divergence_between_groups(ts_obj, groups)

# convert pandas df to R and save
group_div_r <- py_to_r(group_div_df)
write.csv(group_div_r, file.path(run_prefix, "group_divergence.csv"), row.names = TRUE)

# pick individuals for pairwise divergence & geographic distances:
# Combine all groups into a vector of unique individuals

# Get all group values from the Python dict
all_groups <- py_to_r(groups$values())

# Flatten and extract unique individuals
ind_list <- unique(unlist(all_groups))

# simpler in R: call a small py snippet to flatten
py_run_string("
def flatten_groups(d):
    res = []
    for k in d:
        res.extend(list(d[k]))
    return list(dict.fromkeys(res))
flat_inds = flatten_groups(groups)
")

flat_inds <- py$flat_inds

# write metadata and vcf for these individuals
vcf_out <- file.path(run_prefix, "sampled_individuals.vcf")
meta_out <- file.path(run_prefix, "sampled_individuals_metadata.tsv")

# Compute and write VCF and metadata
helpers$write_vcf_and_metadata(ts_obj, flat_inds, vcf_out, meta_out)

# compute pairwise divergence and geo distances
pairs, divs, geogs <- helpers$pairwise_individual_divergence_and_geo(ts_obj, flat_inds)

# convert to R
# reticulate returns tuple; access via py object
py_pairs <- py$pairs
py_divs <- py$divs
py_geogs <- py$geogs
pairs_r <- py_to_r(py_pairs)
divs_r <- py_to_r(py_divs)
geogs_r <- py_to_r(py_geogs)

# Build dataframe for scatter plot
n_pairs <- length(divs_r)
df_pairs <- data.frame(
  pair_index = seq_len(n_pairs),
  divergence = as.numeric(divs_r),
  geog_dist = as.numeric(geogs_r)
)


# plot divergence vs distance
p <- ggplot(df_pairs, aes(x = geog_dist, y = divergence)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(x = "Geographic distance", y = "Genetic divergence (differences per site)")
ggsave(file.path(run_prefix, "divergence_vs_distance.png"), p, width = 6, height = 6, dpi = 300)

# compute IBD fraction for flat_inds
ibd_df <- helpers$pairwise_ibd_fraction(ts_obj, flat_inds)

# save IBD table (pandas df -> csv)
py_run_string("import pandas as _pd; \
_ibd = ibd_df; \
_ibd.to_csv('ibd_fraction.csv')")

# move to run_prefix
file.rename("ibd_fraction.csv", file.path(run_prefix, "ibd_fraction.csv"))

cat("Done. outputs in:", run_prefix, "\n")
