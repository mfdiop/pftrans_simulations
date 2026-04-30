#!/usr/bin/env Rscript

# ================================================================
# true_ibd_from_trees.R
# ------------------------------------------------
# Compute TRUE IBD from tree-sequence files (.trees) using tskibd.
# - Runs tskibd per chromosome (in a temp dir to avoid clutter)
# - Standardizes output: Id1, Id2, Chrom, Start, End, (optional Tmrca, Cm)
# - Writes per-chromosome TSVs and a merged "true_ibd.all.tsv"
# - Optional: write TMRCA-filtered truth (if column present)
#
# Requires: "tskibd" binary on PATH (or pass --tskibd /path/to/tskibd)
# No R package dependencies beyond base R.
# ================================================================

# ----------------------------- CLI ------------------------------
args <- commandArgs(trailingOnly = TRUE)

usage <- "
Usage:
  Rscript true_ibd_from_trees.R \\
    --trees out/1.trees,out/2.trees,out/3.trees \\
    --chrnos 1,2,3 \\
    --r 6.66666666666667e-07 \\
    --mincm 2 \\
    --outdir out_true_ibd \\
    [--tskibd tskibd] [--max_tmrca 300,1000,3000]

Notes:
  - r is Morgans per base (e.g., 0.01/15000 ≈ 6.6667e-07).
  - mincm is the minimum IBD length in cM reported by tskibd.
  - If --max_tmrca is provided AND the tskibd output has a Tmrca column,
    extra filtered files with suffix _maxtmrca_<X>.tsv are written.
"

kv <- list()
if (length(args) == 0) { cat(usage, "\n"); quit(status = 1) }
for (i in seq(1, length(args), by = 2)) {
  k <- sub("^--", "", args[i]); v <- if (i + 1 <= length(args)) args[i + 1] else ""
  kv[[k]] <- v
}

need <- c("trees","chrnos","r","mincm","outdir")
miss <- need[!need %in% names(kv)]
if (length(miss) > 0) { cat("Missing:", paste(miss, collapse=", "), "\n\n", usage); quit(status = 1) }

split_commas <- function(x) unlist(strsplit(x, "\\s*,\\s*"))
trees  <- split_commas(kv[["trees"]])
chrnos <- as.integer(split_commas(kv[["chrnos"]]))
r      <- as.numeric(kv[["r"]])               # Morgans/bp
mincm  <- as.numeric(kv[["mincm"]])
outdir <- kv[["outdir"]]
tskibd_bin <- if ("tskibd" %in% names(kv)) kv[["tskibd"]] else "tskibd"
max_tmrca  <- if ("max_tmrca" %in% names(kv)) as.integer(split_commas(kv[["max_tmrca"]])) else integer(0)

if (length(trees) != length(chrnos)) stop("trees and chrnos lengths must match")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bp_per_cm <- as.integer(0.01 / r)            # e.g., 15000 bp per cM if r=0.01/15000
if (bp_per_cm <= 0L) stop("bp_per_cm <= 0; check r value")
sample_window <- max(1L, as.integer(0.01 * bp_per_cm))  # ≈ 1% of a cM in bp

message(sprintf("bp_per_cm=%d, sample_window(bp)=%d, mincm=%.3f cM", bp_per_cm, sample_window, mincm))

# ------------------- helpers -------------------
read_ibd <- function(fn) {
  # Try headered first, fallback to headerless
  df <- tryCatch(read.table(fn, header=TRUE, sep="\t", check.names=FALSE, comment.char=""),
                 error=function(e) read.table(fn, header=FALSE, sep="\t", check.names=FALSE, comment.char=""))
  # Standardize essential columns (best-effort)
  nms <- tolower(names(df))
  map <- list(
    id1   = which(nms %in% c("id1","sample1","ind1","i1"))[1],
    id2   = which(nms %in% c("id2","sample2","ind2","i2"))[1],
    chrom = which(nms %in% c("chrom","chr","chromosome"))[1],
    start = which(nms %in% c("start","begin","pos_start"))[1],
    end   = which(nms %in% c("end","stop","pos_end"))[1],
    tmrca = which(nms %in% c("tmrca","tmrca_gen"))[1],
    cm    = which(nms %in% c("cm","length_cm","len_cm"))[1]
  )
  rename_if <- function(idx, to) if (!is.na(idx)) names(df)[idx] <<- to
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

write_tsv <- function(df, fn) write.table(df, fn, sep="\t", row.names=FALSE, quote=FALSE)

run_tskibd_one <- function(tree, chr) {
  # Run tskibd <chr> <bp_per_cm> <sample_window_bp> <mincm> <tree>
  tmp <- sprintf("tskibd_%s_%s", chr, as.integer(runif(1,1,1e9)))
  dir.create(tmp, showWarnings=FALSE)
  owd <- getwd(); on.exit({ setwd(owd); unlink(tmp, recursive=TRUE) }, add=TRUE)
  setwd(tmp)
  cmd <- c(as.character(chr), as.character(bp_per_cm),
           as.character(sample_window), as.character(mincm), normalizePath(file.path(owd, tree)))
  out <- tryCatch(system2(tskibd_bin, args=cmd, stdout=TRUE, stderr=TRUE),
                  warning=function(w) w, error=function(e) e)
  ibd_path <- file.path(getwd(), sprintf("%s.ibd", chr))
  if (!file.exists(ibd_path)) {
    stop("tskibd did not produce ", ibd_path, "\nOutput:\n", paste(out, collapse="\n"))
  }
  read_ibd(ibd_path)
}

# ------------------- main loop -------------------
all_list <- list()

for (i in seq_along(chrnos)) {
  chr <- chrnos[i]; tree <- trees[i]
  message(sprintf("[chr %d] tskibd on %s ...", chr, tree))
  df <- run_tskibd_one(tree, chr)
  
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
message("\nDone. Merged truth: ", merged_fn, "\n")
