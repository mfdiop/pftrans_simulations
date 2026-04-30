#!/usr/bin/env Rscript
# run_spatial_full_pipeline.R
# Full spatial simulation pipeline: SLiM -> pyslim/tskit -> msprime -> analyses (VCF, IBD, Fst, diversity, divergence vs distance)

# I modified the script: "slim_msprime_wrapper.R" from David by adding the spatial SLiM simulation and python functions to estimate some stats using tskit (VCF, IBD, Fst, diversity, divergence vs distance)

suppressPackageStartupMessages({
  library(optparse)
  library(reticulate)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# --- small helpers for reticulate conversions
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

# ---- CLI ----
option_list <- list(
  make_option("--slim_bin", type="character", default="slim"),
  make_option("--slim_script", type="character", default="spatial_pop_improved.slim"),
  make_option("--outdir", type="character", default="analysis_spatial"),
  make_option("--chrno", type="integer", default=1),
  make_option("--genome_set_id", type="integer", default=1),
  make_option("--recomb", type="double", default=6.67e-7),
  make_option("--mu", type="double", default=1e-8),
  make_option("--ancestral_Ne", type="integer", default=1000),
  make_option("--remember_gen", type="integer", default=1000),
  make_option("--nsam", type="integer", default=200),
  make_option("--seed", type="integer", default=42)
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
run_prefix <- file.path(opt$outdir, sprintf("%d_%d", opt$genome_set_id, opt$chrno))
dir.create(run_prefix, recursive = TRUE, showWarnings = FALSE)

# Bind reticulate to python env (tries env vars, then fallback)
try({
  py <- Sys.getenv("RETICULATE_PYTHON", "")
  if (nzchar(py)) use_python(py, required = TRUE) else use_condaenv("slim-msprime", required = FALSE)
}, silent = TRUE)

# ensure required python packages will be present; install if desired
py_require(c("pyslim", "msprime", "tskit", "numpy", "pandas"), action = "ignore")

# import python modules
pyslim <- import("pyslim", convert = FALSE)
msprime <- import("msprime", convert = FALSE)
tskit <- import("tskit", convert = FALSE)
np <- import("numpy", convert = TRUE)
pd <- import("pandas", convert = FALSE)

# ---- write inline Python helper module to file and import it ----
py_module_path <- file.path(run_prefix, "tskit_spatial_helpers_inline.py")
py_code <- '
import tskit, pyslim, msprime, numpy as np, pandas as pd, math, gzip

def _to_int(x):
    """Ensure an input is an integer (handles numpy float types from reticulate)."""
    try:
        return int(x)
    except Exception:
        return int(float(x))

def recapitate_and_mutate(ts_path, recap_out_path, recombination_rate=1e-8, ancestral_Ne=1000, mu=1e-8, random_seed=None):
    ts = tskit.load(ts_path)
    r_ts = pyslim.recapitate(ts, recombination_rate=float(recombination_rate), ancestral_Ne=float(ancestral_Ne), random_seed=random_seed)
    mts = msprime.sim_mutations(r_ts, rate=float(mu), model=msprime.SLiMMutationModel(type=0), keep=True, random_seed=random_seed)
    mts.dump(recap_out_path)
    return mts

def parse_sample_map_lines(lines):
    rows = []
    for ln in lines:
        if not isinstance(ln, str):
            ln = str(ln)
        if ln.startswith("SAMPLE_MAP\\t"):
            parts = ln.strip().split("\\t")
            if len(parts) >= 7:
                gen = int(parts[1]); ind = int(parts[2])
                nodeA = int(parts[3]); nodeB = int(parts[4])
                x = float(parts[5]); y = float(parts[6])
                rows.append({"generation": gen, "individual": ind, "nodeA": nodeA, "nodeB": nodeB, "x": x, "y": y})
    return pd.DataFrame(rows)

def sample_groups(ts, alive_time=0, ancient_time=None, W=35, w=5, rng_seed=23):
    rng = np.random.default_rng(rng_seed)
    alive = pyslim.individuals_alive_at(ts, alive_time)
    alive = [int(a) for a in alive]
    locs = np.asarray(ts.individuals_location)[alive, :]

    groups = {}
    groups["topleft"] = [alive[i] for i in range(len(alive)) if locs[i,0] < w and locs[i,1] < w]
    groups["topright"] = [alive[i] for i in range(len(alive)) if locs[i,0] < w and locs[i,1] > W - w]
    groups["bottomleft"] = [alive[i] for i in range(len(alive)) if locs[i,0] > W - w and locs[i,1] < w]
    groups["bottomright"] = [alive[i] for i in range(len(alive)) if locs[i,0] > W - w and locs[i,1] > W - w]
    groups["center"] = [alive[i] for i in range(len(alive)) if (abs(locs[i,0] - W/2) < w/2) and (abs(locs[i,1] - W/2) < w/2)]
    if ancient_time is not None:
        old_ones = pyslim.individuals_alive_at(ts, ancient_time)
        old_ones = [int(x) for x in old_ones]
        if len(old_ones) > 0:
            k = min(5, len(old_ones))
            groups["ancient"] = list(rng.choice(old_ones, size=k, replace=False))
        else:
            groups["ancient"] = []
    else:
        groups["ancient"] = []
    return groups

def group_sampled_nodes(ts, groups):
    # returns list of lists of node ids (integers) per group
    sampled_nodes = []
    for k in list(groups.keys()):
        nodes_of_group = []
        for ind in groups[k]:
            ind_int = int(ind)
            nodes = list(ts.individual(ind_int).nodes)
            nodes_of_group.extend([int(x) for x in nodes])
        sampled_nodes.append(nodes_of_group)
    return sampled_nodes

def divergence_between_groups(ts, groups):
    # groups: dict group->list(individual ids)
    sampled_nodes = group_sampled_nodes(ts, groups)
    div = ts.divergence(sampled_nodes)
    return pd.DataFrame(div, index=list(groups.keys()), columns=list(groups.keys()))

def pairwise_individual_divergence_and_geo(ts, ind_ids):
    ind_ids = [int(i) for i in ind_ids]
    ind_nodes = [list(ts.individual(i).nodes) for i in ind_ids]
    pairs = [(i,j) for i in range(len(ind_ids)) for j in range(i, len(ind_ids))]
    div = ts.divergence(ind_nodes, indexes=pairs)
    # geographic distances
    locs = np.asarray(ts.individuals_location)
    geog = np.zeros(len(pairs))
    for k,(i,j) in enumerate(pairs):
        xi, yi = locs[ind_ids[i],0], locs[ind_ids[i],1]
        xj, yj = locs[ind_ids[j],0], locs[ind_ids[j],1]
        geog[k] = math.sqrt((xi-xj)**2 + (yi-yj)**2)
    return pairs, div, geog

def write_vcf_and_metadata(ts, ind_ids, vcf_path, meta_path):
    ind_ids = [int(i) for i in ind_ids]
    # nodes for individuals (diploid -> two nodes per individual)
    flat_nodes = []
    names = []
    meta_rows = []
    for ind in ind_ids:
        nodes = list(ts.individual(ind).nodes)
        if len(nodes) == 0:
            continue
        # ensure two nodes
        flat_nodes.extend([int(x) for x in nodes])
        names.append(f"tsk_{ind}")
        birth_time = ts.node(nodes[0]).time if len(nodes)>0 else None
        loc = ts.individual(ind).location
        meta_rows.append({"vcf_label": f"tsk_{ind}", "tskit_id": ind, "birth_time_ago": birth_time, "x": loc[0] if loc is not None else None, "y": loc[1] if loc is not None else None})
    # simplify tree to these nodes (keeps metadata)
    if len(flat_nodes) == 0:
        raise ValueError("No sample nodes found for writing VCF")
    ts_sub = ts.simplify(flat_nodes, keep_input_roots=True)
    # write VCF
    with open(vcf_path, "w") as f:
        ts_sub.write_vcf(f, ploidy=2, individual_names=names)
    pd.DataFrame(meta_rows).to_csv(meta_path, sep="\\t", index=False)
    return

def pairwise_ibd_fraction(ts, ind_ids):
    ind_ids = [int(i) for i in ind_ids]
    node_lists = [list(ts.individual(i).nodes) for i in ind_ids]
    n = len(ind_ids)
    L = ts.sequence_length
    ibd = np.zeros((n,n), dtype=float)
    for tree in ts.trees():
        interval = tree.interval
        for i in range(n):
            for j in range(n):
                shared = False
                for a in node_lists[i]:
                    for b in node_lists[j]:
                        if tree.mrca(int(a), int(b)) != tskit.NULL:
                            shared = True
                            break
                    if shared:
                        break
                if shared:
                    ibd[i,j] += interval
    ibd = ibd / L
    return pd.DataFrame(ibd, index=ind_ids, columns=ind_ids)

def nucleotide_diversity(ts, ind_ids):
    ind_ids = [int(i) for i in ind_ids]
    node_lists = [[int(x) for x in ts.individual(i).nodes] for i in ind_ids]
    # use tskit.diversity on each individual's nodes? Prefer average pairwise within group: flatten group's nodes
    # Here we compute per-group pi by calling ts.diversity on the flattened nodes
    samples = [n for pair in node_lists for n in pair]
    if len(samples) == 0:
        return 0.0
    return ts.diversity(samples)

def hudson_fst_matrix(ts, groups):
    # compute pairwise Hudson Fst between groups
    group_order = list(groups.keys())
    m = len(group_order)
    fst_mat = np.zeros((m,m))
    for i in range(m):
        for j in range(i,m):
            a = groups[group_order[i]]
            b = groups[group_order[j]]
            # get nodes for each group
            anodes = []
            bnodes = []
            for ind in a:
                anodes.extend(list(ts.individual(int(ind)).nodes))
            for ind in b:
                bnodes.extend(list(ts.individual(int(ind)).nodes))
            if len(anodes) == 0 or len(bnodes) == 0:
                fst_val = np.nan
            else:
                # ts.allele_frequency will be used per variant via iterator over variants
                # to compute Hudson numerator and denominator sums
                num = 0.0
                den = 0.0
                for v in ts.variants():
                    g = v.genotypes
                    # compute alt counts in node lists (note: v.genotypes aligns to ts.samples())
                    # We map sample index across ts.subset? Simpler: compute allele counts by indexing
                    # Create a boolean mask for which sample indices correspond to our nodes
                    # Here we retrieve sample indices via ts.samples().index?
                    # Simpler approach: use allele counts via v.alleles and ts.genotypes? For speed we fallback to get alt freq by scanning g.
                    # Build mask arrays (convert samples to python list of ints)
                    samples = list(ts.samples())
                    idx_map = {s:i for i,s in enumerate(samples)}
                    # get indices for anodes/bnodes relative to samples
                    a_idx = [idx_map[int(n)] for n in anodes if int(n) in idx_map]
                    b_idx = [idx_map[int(n)] for n in bnodes if int(n) in idx_map]
                    if len(a_idx)==0 or len(b_idx)==0:
                        continue
                    alt_a = sum([1 for ii in a_idx if g[ii]==1])
                    alt_b = sum([1 for ii in b_idx if g[ii]==1])
                    na = len(a_idx); nb = len(b_idx)
                    pa = alt_a/na; pb = alt_b/nb
                    # Hudson components
                    num += (pa - pb)**2 - (pa*(1-pa)/(na-1) + pb*(1-pb)/(nb-1))
                    den += (pa*(1-pb) + pb*(1-pa))
                if den == 0:
                    fst_val = np.nan
                else:
                    fst_val = num/den
            fst_mat[i,j] = fst_val
            fst_mat[j,i] = fst_val
    return pd.DataFrame(fst_mat, index=group_order, columns=group_order)
'
# write file
writeLines(py_code, con = py_module_path)

# import the inline module
helpers <- import_from_path("tskit_spatial_helpers_inline", path = run_prefix)

# ---- Run SLiM ----
slim_bin <- Sys.which(opt$slim_bin)
if (nchar(slim_bin) == 0) stop("SLiM binary not found; set --slim_bin or install SLiM")

# pass runtime -d args if required; assume spatial SLiM script handles passed constants OUTID, REMEMBER_T_GEN, TOTAL_GEN
slim_args <- c(sprintf("-d OUTID=%d", opt$chrno),
               sprintf("-d REMEMBER_T_GEN=%d", opt$remember_gen),
               sprintf("-d TOTAL_GEN=%d", opt$remember_gen * 2),
               sprintf("-d W=%d", 35)
)
cmd_args <- c(slim_args, normalizePath(opt$slim_script))

stdout_f <- file.path(run_prefix, "slim_stdout.txt")
stderr_f <- file.path(run_prefix, "slim_stderr.txt")

message("Running SLiM...")
res <- system2(slim_bin, args = cmd_args, stdout = stdout_f, stderr = stderr_f, wait = TRUE)
if (res != 0) stop(sprintf("SLiM failed (exit %s). See %s & %s", res, stdout_f, stderr_f))

# ---- Read SLiM outputs: parse SAMPLE_MAP, PEDIGREE, DAF, True_Ne from stdout ----
stdout_lines <- readLines(stdout_f, warn = FALSE)
# extract DAF and True_Ne and restart_counter if present
ne_lines  <- grep("^True_Ne\\t", stdout_lines, value = TRUE)
daf_lines <- grep("^DAF\\t", stdout_lines, value = TRUE)
ped_lines <- grep("^PEDIGREE\\t", stdout_lines, value = TRUE)
sample_map_lines <- grep("^SAMPLE_MAP\\t", stdout_lines, value = TRUE)

# write those logs
writeLines(ne_lines, con = file.path(run_prefix, "slim_true_ne.log"))
writeLines(daf_lines, con = file.path(run_prefix, "slim_daf.log"))
writeLines(ped_lines, con = file.path(run_prefix, "slim_pedigree.log"))
writeLines(sample_map_lines, con = file.path(run_prefix, "slim_sample_map.log"))

# ---- Load trees (produced by SLiM spatial script) ----
# expect file named spatial_sim_out_<OUTID>.trees or tmp_slim_out_single_pop_... adapt as needed
trees_candidates <- c(
  file.path(getwd(), sprintf("spatial_sim_out_%d.trees", opt$chrno)),
  file.path(getwd(), sprintf("tmp_slim_out_single_pop_%s.trees", opt$chrno)),
  file.path(getwd(), "spatial_sim.trees")
)
trees_file <- trees_candidates[file.exists(trees_candidates)][1]
if (is.na(trees_file) || length(trees_file)==0) stop("No trees file found in working directory; ensure SLiM wrote a .trees file")

message("Loading trees: ", trees_file)
# use tskit via retained python module (helpers will load as needed)
# recapitate + mutate using helper
recap_out <- file.path(run_prefix, "spatial_sim.recap.trees")
mts <- helpers$recapitate_and_mutate(trees_file, recap_out,
                                     recombination_rate = opt$recomb,
                                     ancestral_Ne = opt$ancestral_Ne,
                                     mu = opt$mu,
                                     random_seed = as.integer(opt$seed))

# ---- Subsampling & groups (use helper sample_groups) ----
ts_obj <- mts  # python tskit tree sequence object via reticulate
groups_py <- helpers$sample_groups(ts_obj, as.integer(0), as.integer(opt$remember_gen), W = 35, w = 5, rng_seed = as.integer(opt$seed))

# convert group summary to R and save
groups_r <- py_to_r(groups_py)
groups_summary <- sapply(names(groups_r), function(k) length(groups_r[[k]]))
write.csv(data.frame(group = names(groups_summary), size = as.integer(groups_summary)), file = file.path(run_prefix, "groups_summary.csv"), row.names = FALSE)

# ---- compute group-level divergence matrix ----
group_div_df_py <- helpers$divergence_between_groups(ts_obj, groups_py)
group_div_df <- py_to_r(group_div_df_py)
write.csv(group_div_df, file = file.path(run_prefix, "group_divergence.csv"), row.names = TRUE)

# ---- flatten group individuals into unique list for per-individual analyses ----
# flatten in Python to avoid reticulate type issues
py_run_string("def flatten_groups(d):\n    res = []\n    for k in d:\n        for v in d[k]:\n            if v not in res:\n                res.append(int(v))\n    return res\nflat_inds = flatten_groups(groups_py)\n")
flat_inds <- py$flat_inds
# write sampled individuals metadata using helper (VCF + meta)
vcf_out <- file.path(run_prefix, "sampled_individuals.vcf")
meta_out <- file.path(run_prefix, "sampled_individuals_metadata.tsv")
helpers$write_vcf_and_metadata(ts_obj, flat_inds, vcf_out, meta_out)

# ---- pairwise divergence vs geographic distance (individual-level) ----
py_pairs_div_geo <- helpers$pairwise_individual_divergence_and_geo(ts_obj, flat_inds)
# py returns (pairs, divs, geogs) as tuple accessible via py object names
pairs_py <- py_pairs_div_geo[[1]]
divs_py <- py_pairs_div_geo[[2]]
geogs_py <- py_pairs_div_geo[[3]]

pairs_r <- py_to_r(pairs_py)
divs_r <- py_to_r(divs_py)
geogs_r <- py_to_r(geogs_py)

# Build dataframe for plotting
# pairs is list of tuples; but we only need geog vs div
df_pairs <- data.frame(geog = as.numeric(geogs_r), divergence = as.numeric(divs_r))
p1 <- ggplot(df_pairs, aes(x = geog, y = divergence)) +
  geom_point(alpha = 0.6) + theme_minimal() + labs(x = "Geographic distance", y = "Genetic divergence (diffs/site)")
ggsave(filename = file.path(run_prefix, "divergence_vs_distance.png"), plot = p1, width = 6, height = 6, dpi = 300)

# ---- pairwise IBD ----
ibd_df_py <- helpers$pairwise_ibd_fraction(ts_obj, flat_inds)
ibd_df <- py_to_r(ibd_df_py)
write.csv(ibd_df, file = file.path(run_prefix, "pairwise_ibd_fraction.csv"), row.names = TRUE)

# ---- nucleotide diversity (pi) per group ----
diversity_list <- lapply(names(groups_r), function(g) {
  inds <- groups_r[[g]]
  py_res <- helpers$nucleotide_diversity(ts_obj, inds)
  as.numeric(py_to_r(py_res))
})
diversity_df <- data.frame(group = names(groups_r), pi = unlist(diversity_list))
write.csv(diversity_df, file = file.path(run_prefix, "group_nucleotide_diversity.csv"), row.names = FALSE)

# ---- Hudson Fst between groups ----
fst_py <- helpers$hudson_fst_matrix(ts_obj, groups_py)
fst_df <- py_to_r(fst_py)
write.csv(fst_df, file = file.path(run_prefix, "group_hudson_fst.csv"), row.names = TRUE)

# ---- SFS (from mts) via inline helper previously available in earlier scripts ----
# compute sfs vector in R by calling tskit directly (via python)
sfs_py <- py_eval("[(v.num_alleles) for v in []]", convert = TRUE) # dummy-safe default
# We can compute SFS via simple python snippet:
py_run_string("
def site_frequency_spectrum_py(ts):
    n = ts.num_samples
    sfs = np.zeros(n+1, dtype=int)
    for v in ts.variants():
        alt = np.sum(v.genotypes == 1)
        sfs[alt] += 1
    return sfs
sfs_vec = site_frequency_spectrum_py(mts)
")
sfs_vec <- py$sfs_vec
sfs_r <- py_to_r(sfs_vec)
sfs_df <- data.frame(k = seq_along(sfs_r)-1, count = as.integer(sfs_r))
write.csv(sfs_df, file = file.path(run_prefix, "sfs.csv"), row.names = FALSE)

# ---- True Ne and DAF logs (already parsed from SLiM stdout) ----
# Parse True_Ne lines into a table (GEN, NE)
if (length(ne_lines) > 0) {
  ne_parsed <- do.call(rbind, lapply(ne_lines, function(ln) {
    parts <- strsplit(ln, "\\t", fixed = TRUE)[[1]]
    if (length(parts) >= 3) {
      data.frame(GEN = as.integer(parts[2]), NE = as.integer(parts[3]))
    } else NULL
  }))
  if (!is.null(ne_parsed)) write.csv(ne_parsed, file = file.path(run_prefix, "true_ne.csv"), row.names = FALSE)
}
if (length(daf_lines) > 0) {
  daf_parsed <- do.call(rbind, lapply(daf_lines, function(ln) {
    parts <- strsplit(ln, "\\t", fixed = TRUE)[[1]]
    if (length(parts) >= 3) {
      data.frame(GEN = as.integer(parts[2]), DAF = as.numeric(parts[3]))
    } else NULL
  }))
  if (!is.null(daf_parsed)) write.csv(daf_parsed, file = file.path(run_prefix, "daf.csv"), row.names = FALSE)
}

# ---- Final output summary ----
cat("Outputs written to: ", run_prefix, "\n")
cat("Files:\n")
cat(list.files(run_prefix, full.names = TRUE), sep = "\n")





import msprime, pyslim, tskit, numpy as np, pandas as pd, matplotlib.pyplot as plt

# 1. Load SLiM tree sequence
slim_ts = pyslim.load("spatial_pf_sim_1.trees")

# 2. Recapitate (add ancestral history using coalescent)
recap_ts = pyslim.recapitate(slim_ts, recombination_rate=7.4e-7, ancestral_Ne=5000)

# 3. Add neutral mutations
mut_ts = msprime.sim_mutations(recap_ts, rate=1.5e-9, model=msprime.SLiMMutationModel(type=0), keep=True)
mut_ts.dump("spatial_pf_sim_1.recap.trees")

# 4. Sample individuals by region
W = 40
alive = pyslim.individuals_alive_at(mut_ts, 0)
locs = mut_ts.individuals_location[alive, :]
groups = {
  'topleft' : alive[(locs[:,0]<10) & (locs[:,1]<10)],
  'topright': alive[(locs[:,0]>30) & (locs[:,1]<10)],
  'bottomleft': alive[(locs[:,0]<10) & (locs[:,1]>30)],
  'bottomright': alive[(locs[:,0]>30) & (locs[:,1]>30)],
  'center': alive[(np.abs(locs[:,0]-20)<5) & (np.abs(locs[:,1]-20)<5)]
}

# 5. Compute pairwise group divergences
sampled_nodes = []
for k in groups:
  nodes = []
for ind in groups[k]:
  nodes.extend(mut_ts.individual(ind).nodes)
sampled_nodes.append(nodes)

div_matrix = mut_ts.divergence(sampled_nodes)
div_df = pd.DataFrame(div_matrix, index=groups.keys(), columns=groups.keys())
div_df.to_csv("divergence_matrix.csv")
print(div_df)

# 6. Visualization of spatial structure
fig, ax = plt.subplots(figsize=(6,6), dpi=300)
for k, inds in groups.items():
  ax.scatter(locs[inds,0], locs[inds,1], label=k, s=15)
ax.legend()
ax.set_title("Spatial structure of simulated P. falciparum population")
plt.savefig("spatial_pf_population.png", dpi=300)
# 7. Save sampled individuals to VCF
# flatten individuals to nodes  

