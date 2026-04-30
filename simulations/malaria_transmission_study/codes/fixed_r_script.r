#!/usr/bin/env Rscript
  
suppressPackageStartupMessages({
  library(optparse)
})

# Set non-interactive graphics device for HPC environments (no X11)
options(bitmapType = "cairo")

# ----- CLI (parse early to get options) -----
option_list <- list(
  make_option("--chrno", type = "integer", default = 1),
  make_option("--seqlen", type = "integer", default = 100 * 15000),
  make_option("--selpos", type = "integer", default = as.integer(0.33 * 100 * 15000)),
  make_option("--num_origins", type = "integer", default = 1),
  make_option("--N", type = "integer", default = 10000),
  make_option("--h", type = "double",  default = 0.5),
  make_option("--s", type = "double",  default = 0.3),
  make_option("--g_sel_start", type = "integer", default = 80),
  make_option("--r", type = "double",  default = 0.01 / 15000),
  make_option("--sim_relatedness", type = "integer", default = 0),
  make_option("--g_ne_change_start", type = "integer", default = 200),
  make_option("--N0", type = "integer", default = 1000),
  make_option("--u", type = "double",  default = 1e-8),
  make_option("--nsam", type = "integer", default = 1000),
  make_option("--outdir", type = "character", default = "."),
  make_option("--test", action = "store_true", default = FALSE,
              help = "Use smaller nsam and seqlen for a quick test"),
  make_option("--genome_set_id", type = "integer", help = "Required."),
  make_option("--slim_script", type = "character", default = "my_test/single_pop.slim"),
  make_option("--slim_bin", type = "character", default = "slim") # Changed default from slim.exe
)
  
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Print start message
message(sprintf("Start simulation using %s", opt$slim_script))

# Validate required arguments
if (is.null(opt$genome_set_id)) {
  stop("--genome_set_id is required")
}

# ----- Setup Python/reticulate BEFORE loading reticulate -----
# Clear any cached Python configurations
Sys.unsetenv("RETICULATE_PYTHON_FALLBACK")

# Try to find Python in this order:
# 1. RETICULATE_PYTHON environment variable
# 2. Active conda environment
# 3. System python3
python_path <- NULL

# Check for explicit RETICULATE_PYTHON
if (nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
  python_path <- Sys.getenv("RETICULATE_PYTHON")
  message(sprintf("Using Python from RETICULATE_PYTHON: %s", python_path))
}

# Check for active conda environment
if (is.null(python_path) && nzchar(Sys.getenv("CONDA_PREFIX", ""))) {
  conda_prefix <- Sys.getenv("CONDA_PREFIX")
  python_candidates <- c(
    file.path(conda_prefix, "bin", "python"),
    file.path(conda_prefix, "bin", "python3")
  )
  for (py in python_candidates) {
    if (file.exists(py)) {
      python_path <- py
      message(sprintf("Using Python from active conda environment: %s", python_path))
      break
    }
  }
}

# Check for system python3
if (is.null(python_path)) {
  system_python <- Sys.which("python3")
  if (nzchar(system_python)) {
    python_path <- as.character(system_python)
    message(sprintf("Using system Python: %s", python_path))
  }
}

# Check for system python (fallback)
if (is.null(python_path)) {
  system_python <- Sys.which("python")
  if (nzchar(system_python)) {
    python_path <- as.character(system_python)
    message(sprintf("Using system Python: %s", python_path))
  }
}

# Validate Python path
if (is.null(python_path) || !file.exists(python_path)) {
  stop(sprintf(
    "Could not find valid Python installation.\n",
    "  Checked RETICULATE_PYTHON: %s\n",
    "  Checked CONDA_PREFIX: %s\n",
    "  Checked system python3: %s\n",
    "Please ensure Python is installed and accessible, or set RETICULATE_PYTHON explicitly."
  ))
}

# Set RETICULATE_PYTHON before loading reticulate
Sys.setenv(RETICULATE_PYTHON = python_path)

# Now load reticulate
suppressPackageStartupMessages({
  library(reticulate)
})

# Verify Python is working
tryCatch({
  py_config()
  message("Python configuration successful")
}, error = function(e) {
  message("Python configuration:")
  message(capture.output(py_config()))
  stop(sprintf("Failed to configure Python: %s", e$message))
})

# Load other required packages
suppressPackageStartupMessages({
  library(tidyverse)
#  library(dplyr)
#  library(stringr)
#  library(tidyr)
})

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

have_gg <- requireNamespace("ggplot2", quietly = TRUE)

save_plot_png <- function(plot_obj, filename, width=7, height=5, res=120) {
  # Ensure directory exists
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  if (have_gg && inherits(plot_obj, "ggplot")) {
    # Use cairo device for ggplot2
    tryCatch({
      ggplot2::ggsave(filename, plot = plot_obj, width = width, height = height, 
                      dpi = res * (96/72), device = "png", type = "cairo")
    }, error = function(e) {
      # Fallback without specifying type
      ggplot2::ggsave(filename, plot = plot_obj, width = width, height = height, 
                      dpi = res * (96/72), device = "png")
    })
  } else {
    # plot_obj is a function that draws the base plot
    # Use cairo png device explicitly to avoid X11
    tryCatch({
      png(filename, width = width * res, height = height * res, 
          units = "px", res = res, type = "cairo")
      on.exit(dev.off(), add = TRUE)
      plot_obj()
    }, error = function(e) {
      # If cairo fails, try cairo-png
      tryCatch({
        png(filename, width = width * res, height = height * res, 
            units = "px", res = res, type = "cairo-png")
        on.exit(dev.off(), add = TRUE)
        plot_obj()
      }, error = function(e2) {
        # Last resort: try without specifying type
        message("Warning: Using default PNG device, may fail on headless systems")
        png(filename, width = width * res, height = height * res, 
            units = "px", res = res)
        on.exit(dev.off(), add = TRUE)
        plot_obj()
      })
    })
  }
}

# Validate SLiM binary and script path early
slim_bin <- opt$slim_bin
slim_bin_found <- Sys.which(slim_bin)
if (nchar(slim_bin_found) == 0) {
  stop(sprintf("SLiM binary '%s' not found on PATH. Set --slim_bin / install SLiM.", slim_bin))
}
message(sprintf("Found SLiM binary: %s", slim_bin_found))

# Resolve script path (absolute or relative to cwd)
slim_script_path <- tryCatch(
  normalizePath(opt$slim_script, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)

if (is.na(slim_script_path)) {
  stop(sprintf("SLiM script not found: %s\n  getwd(): %s\n  tip: pass --slim_script with an absolute path.",
               opt$slim_script, getwd()))
}
message(sprintf("Found SLiM script: %s", slim_script_path))

if (isTRUE(opt$test)) {
  opt$nsam   <- 50L
  opt$seqlen <- as.integer(20 * (0.01 / opt$r))
  opt$selpos <- as.integer(opt$seqlen %/% 3L)
  message("Running in TEST mode with reduced parameters")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ----- Run SLiM -----
run_slim <- function(idx, slim_seed) {
  slim_params <- list(
    L = opt$seqlen,
    selpos = opt$selpos,
    num_origins = opt$num_origins,
    N = opt$N,
    h = opt$h,
    s = opt$s,
    g_sel_start = opt$g_sel_start,
    r = opt$r,
    outid = idx,
    max_restart = 100,
    sim_relatedness = opt$sim_relatedness,
    N0 = opt$N0,
    g_ne_change_start = opt$g_ne_change_start
  )
  
  # Construct command line args for SLiM
  slim_args <- paste(unlist(Map(function(k,v) sprintf("-d %s=%s", k, as.character(v)),
                                names(slim_params), slim_params)), collapse = " ")
  
  seed_arg  <- if (!is.null(slim_seed)) sprintf("-seed %s", slim_seed) else ""
  
  cmd <- sprintf("%s %s %s %s", shQuote(slim_bin_found), slim_args, seed_arg, shQuote(slim_script_path))
  message(">> simulate chrom with id ", idx)
  message(">> ", cmd)

  # Write the command to a file for reproducibility
  cmd_file <- file.path(opt$outdir, sprintf("slim_cmd_%d.txt", idx))
  writeLines(cmd, cmd_file)
  
  # Capture stdout/stderr to files
  out_file <- file.path(opt$outdir, sprintf("slim_stdout_%d.txt", idx))
  err_file <- file.path(opt$outdir, sprintf("slim_stderr_%d.txt", idx))
  
  status <- tryCatch(
    {
      system2("bash", c("-lc", shQuote(cmd)), stdout = out_file, stderr = err_file, wait = TRUE)
    },
    error = function(e) {
      attr(e, "status") <- 127
      e
    }
  )
  
  if (inherits(status, "error") || (!is.null(status) && status != 0)) {
    err_msg <- paste0("SLiM failed; see logs:\n  cmd: ", cmd_file, "\n  out: ", out_file, "\n  err: ", err_file)
    # Also echo last few lines to the console for convenience
    if (file.exists(err_file)) {
      tail_err <- tryCatch(tail(readLines(err_file, warn = FALSE), 15), error = function(e) character())
      if (length(tail_err)) message(paste(tail_err, collapse = "\n"))
    }
    stop(err_msg)
  }
  
  out_txt <- if (file.exists(out_file)) paste(readLines(out_file, warn = FALSE), collapse = "\n") else ""
  list(stdout = out_txt,
       tree_fn = sprintf("tmp_slim_out_single_pop_%s.trees", idx))
}

parse_slim_stdout <- function(slim_stdout) {
  lines <- strsplit(slim_stdout, "\n", fixed = TRUE)[[1]]
  ne_lines  <- c("GEN\tNE")
  daf_lines <- c("GEN\tDAF")
  restart_count <- 0L
  for (ln in lines) {
    if (startsWith(ln, "restart_count")) {
      parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(parts) >= 2) restart_count <- as.integer(parts[2])
    } else if (startsWith(ln, "True_Ne")) {
      ne_lines <- c(ne_lines, sub("^True_Ne\t", "", ln))
    } else if (startsWith(ln, "in")) {
      message(ln)
    } else if (startsWith(ln, "DAF")) {
      daf_lines <- c(daf_lines, sub("^DAF\t", "", ln))
    }
  }
  ne_df  <- suppressMessages(readr::read_tsv(paste0(ne_lines, collapse = "\n"), show_col_types = FALSE))
  daf_df <- suppressMessages(readr::read_tsv(paste0(daf_lines, collapse = "\n"), show_col_types = FALSE))
  ne_df  <- ne_df %>% distinct(GEN, .keep_all = TRUE)
  daf_df <- daf_df %>% distinct(GEN, .keep_all = TRUE)
  list(ne = ne_df, daf = daf_df, restart_count = restart_count)
}

# ----- Python interop (tskit/msprime/pyslim) -----
message("Importing Python modules...")
tryCatch({
  np       <- import("numpy", convert = TRUE)
  tskit    <- import("tskit", convert = FALSE)   # keep as Py objects
  pyslim   <- import("pyslim", convert = FALSE)
  msprime  <- import("msprime", convert = FALSE)
  builtins <- import_builtins()
  message("Python modules imported successfully")
}, error = function(e) {
  stop(sprintf(
    "Failed to import required Python modules.\n",
    "Error: %s\n",
    "Please ensure the following packages are installed in your Python environment:\n",
    "  - numpy\n",
    "  - tskit\n",
    "  - pyslim\n",
    "  - msprime\n",
    "  - pandas\n",
    "Install with: pip install numpy tskit pyslim msprime pandas\n",
    "Or: conda install -c conda-forge numpy tskit pyslim msprime pandas"
  ), e$message)
})

# Silence harmless time-units warning when mixing SLiM ticks with msprime generations
py_run_string("import warnings, msprime; warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)")

# helpers in Python space
py_run_string("
import gzip
import numpy as np
import pandas as pd

def write_peudo_homozygous_vcf(ts_mutated, chrno, out_vcf):
    gt_list = []
    pos_list = []
    ref_list = []
    alt_list = []
    for v in ts_mutated.variants():
        if len(v.alleles) != 2:
            continue
        gt_list.append(v.genotypes)
        i = int(v.alleles[1])
        a = 'ATGC'[i % 4]
        r = 'ATGC'[(i + 1) % 4]
        ref_list.append(r)
        alt_list.append(a)
        pos_list.append(int(v.position))

    header = f\"\"\"##fileformat=VCFv4.2
##source=tskit
##FILTER=<ID=PASS,Description=\"All filters passed\">
##contig=<ID={chrno},length={int(ts_mutated.sequence_length)}>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
\"\"\"
    df1 = pd.DataFrame({
        '#CHROM': chrno,
        'POS': pos_list,
        'ID': '.',
        'REF': ref_list,
        'ALT': alt_list,
        'QUAL': '.',
        'FILTER': 'PASS',
        'INFO': '.',
        'FORMAT': 'GT'
    })
    df2 = pd.DataFrame(gt_list)
    df2.columns = [f'tsk_{n}' for n in df2.columns]
    df2 = df2.astype(str)
    df2 = df2 + '|' + df2
    df = pd.concat([df1.reset_index(drop=True), df2.reset_index(drop=True)], axis=1)
    with gzip.open(out_vcf, 'wt') as f:
        f.write(header)
        df.to_csv(f, sep='\\t', header=True, index=False)

def variant_alt_frequencies(ts):
    n = ts.num_samples
    freqs = []
    for v in ts.variants():
        g = v.genotypes
        alt = np.sum(g==1)
        if alt>0 and alt<n:
            freqs.append(alt/n)
    return freqs

def site_frequency_spectrum(ts):
    n = ts.num_samples
    sfs = np.zeros(n+1, dtype=int)
    for v in ts.variants():
        g = v.genotypes
        alt = np.sum(g==1)
        sfs[alt]+=1
    return sfs
")

# ----- Orchestration -----
message("Starting SLiM simulation...")
# 1) SLiM
slim <- run_slim(idx = opt$chrno, slim_seed = opt$chrno + opt$genome_set_id * 14L)
parsed <- parse_slim_stdout(slim$stdout)
message("SLiM simulation completed")

# 2) Load trees, simplify, recapitate, mutate
message("Loading tree sequence...")
ts <- tskit$load(slim$tree_fn)

# choose half the sample size in individuals
alive_inds <- pyslim$individuals_alive_at(ts, as.integer(0))

# convert to R vector of ints
alive_idx <- as.integer(reticulate::py_to_r(alive_inds))
if (length(alive_idx) < opt$nsam %/% 2) {
  stop(sprintf('Not enough alive individuals to sample. Found %d, need %d', 
               length(alive_idx), opt$nsam %/% 2))
}

set.seed(123)
keep_inds <- sample(alive_idx, size = opt$nsam %/% 2, replace = FALSE)
keep_nodes <- c()

for (i in keep_inds) {
  nodes <- ts$individual(as.integer(i))$nodes
  keep_nodes <- c(keep_nodes, as.integer(reticulate::py_to_r(nodes)))
}

message(sprintf("Simplifying tree sequence (keeping %d nodes)...", length(keep_nodes)))
sts <- ts$simplify(keep_nodes, keep_input_roots = TRUE)

# recapitate
message("Recapitating...")
rts <- pyslim$recapitate(sts,
                         ancestral_Ne = py_as_numeric(opt$N),
                         recombination_rate = py_as_numeric(opt$r),
                         random_seed = as.integer(opt$chrno + opt$genome_set_id * 14L))

cur_sample_nodes <- as.integer(reticulate::py_to_r(rts$samples()))
sts2 <- rts$simplify(cur_sample_nodes)

# delete existing sites then mutate
nsites <- reticulate::py_to_r(sts2$num_sites)
if (nsites > 0) {
  message(sprintf("Removing %d existing mutation sites...", nsites))
  site_list <- as.list(0:(nsites - 1L))
  ts_nosites <- sts2$delete_sites(site_list)
} else {
  ts_nosites <- sts2
}

message("Adding neutral mutations...")
mts <- msprime$sim_mutations(ts_nosites,
                             rate = py_as_numeric(opt$u),
                             model = msprime$SLiMMutationModel(type = as.integer(0)),
                             keep = TRUE)

# ----- Outputs -----
message("Writing output files...")
prefix <- file.path(opt$outdir, sprintf('%d_%d', opt$genome_set_id, opt$chrno))
ofn_slim_restart_count <- sprintf('%s.restart_count', prefix)
ofn_true_ne <- sprintf('%s.true_ne', prefix)
ofn_daf <- sprintf('%s.daf', prefix)
ofn_tree <- sprintf('%s.trees', prefix)
ofn_vcf <- sprintf('%s.vcf.gz', prefix)
ofn_sfs <- sprintf('%s.sfs.tsv', prefix)
ofn_true_ne_png <- sprintf('%s.true_ne.png', prefix)
ofn_daf_png <- sprintf('%s.daf.png', prefix)
ofn_sfs_png <- sprintf('%s.sfs.png', prefix)

writeLines(as.character(parsed$restart_count), con = ofn_slim_restart_count)
readr::write_tsv(parsed$ne,  ofn_true_ne)
readr::write_tsv(parsed$daf, ofn_daf)

# dump trees pre-mutation simplification
sts2$dump(ofn_tree)

# call Python helper to write the VCF
main <- import_main()
main$write_peudo_homozygous_vcf(mts, py_as_integer(opt$chrno), ofn_vcf)

# ----- Compute SFS via Python -----
freqs <- reticulate::py_to_r(main$variant_alt_frequencies(mts))
sfs_vec <- as.integer(reticulate::py_to_r(main$site_frequency_spectrum(mts)))
n_samp <- py_as_integer(mts$num_samples)
# Drop 0 and n bins; keep 1..n-1
if (length(sfs_vec) >= 2) {
  k <- seq_along(sfs_vec) - 0L
  keep <- which(k > 0 & k < n_samp)
  sfs_df <- tibble::tibble(k = k[keep], count = sfs_vec[keep], maf = k[keep] / n_samp)
  readr::write_tsv(sfs_df, ofn_sfs)
} else {
  sfs_df <- tibble::tibble(k = integer(), count = integer(), maf = numeric())
}

# ----- PLOTS -----
message("Generating plots...")

# Wrapper function to safely generate plots
safe_plot <- function(plot_func, filename, plot_name) {
  tryCatch({
    plot_func()
    message(sprintf("  ✓ Generated %s", plot_name))
  }, error = function(e) {
    message(sprintf("  ✗ Failed to generate %s: %s", plot_name, e$message))
    message("  Continuing without this plot...")
  })
}

# 1) True Ne over generations
safe_plot(function() {
  if (have_gg) {
    p_ne <- ggplot2::ggplot(parsed$ne, ggplot2::aes(x = GEN, y = NE)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = "True Ne over generations",
                    x = "Generation", y = "Ne")
    save_plot_png(p_ne, ofn_true_ne_png, width=7, height=5, res=150)
  } else {
    save_plot_png(function(){
      plot(parsed$ne$GEN, parsed$ne$NE, type="l", xlab="Generation", ylab="Ne", main="True Ne over generations")
    }, ofn_true_ne_png, width=7, height=5, res=150)
  }
}, ofn_true_ne_png, "Ne plot")

# 2) DAF (from SLiM stdout) over generations
safe_plot(function() {
  if (have_gg) {
    p_daf <- ggplot2::ggplot(parsed$daf, ggplot2::aes(x = GEN, y = DAF)) +
      ggplot2::geom_line() +
      ggplot2::labs(title = "Derived Allele Frequency over generations (SLiM)",
                    x = "Generation", y = "DAF")
    save_plot_png(p_daf, ofn_daf_png, width=7, height=5, res=150)
  } else {
    save_plot_png(function(){
      plot(parsed$daf$GEN, parsed$daf$DAF, type="l", xlab="Generation", ylab="DAF", main="DAF over generations (SLiM)")
    }, ofn_daf_png, width=7, height=5, res=150)
  }
}, ofn_daf_png, "DAF plot")

# 3) Site Frequency Spectrum from mutated ts
if (nrow(sfs_df) > 0) {
  safe_plot(function() {
    if (have_gg) {
      p_sfs <- ggplot2::ggplot(sfs_df, ggplot2::aes(x = maf, y = count)) +
        ggplot2::geom_col() +
        ggplot2::labs(title = "Site Frequency Spectrum (minor allele freq)",
                      x = "Minor allele frequency", y = "Number of sites")
      save_plot_png(p_sfs, ofn_sfs_png, width=7, height=5, res=150)
    } else {
      save_plot_png(function(){
        barplot(height = sfs_df$count, names.arg = round(sfs_df$maf, 3), las = 2,
                xlab = "Minor allele frequency", ylab = "Number of sites", main = "Site Frequency Spectrum")
      }, ofn_sfs_png, width=10, height=5, res=150)
    }
  }, ofn_sfs_png, "SFS plot")
}

message("\n=== SIMULATION COMPLETE ===")
cat(sprintf('\nOutput files:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n',
            ofn_slim_restart_count, ofn_true_ne, ofn_daf, ofn_tree, ofn_vcf, ofn_sfs,
            ofn_true_ne_png, ofn_daf_png, ofn_sfs_png))
