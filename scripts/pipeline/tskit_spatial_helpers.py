# tskit_spatial_helpers.py
import tskit
import pyslim
import msprime
import numpy as np
import pandas as pd
from math import sqrt

def recapitate_and_mutate(ts_path, recap_out_path, recombination_rate, ancestral_Ne, mu, random_seed=None): 
    ts = tskit.load(ts_path)
    r_ts = pyslim.recapitate(ts, recombination_rate=recombination_rate, ancestral_Ne=ancestral_Ne, random_seed=random_seed)
    mts = msprime.sim_mutations(r_ts, rate=mu, model=msprime.SLiMMutationModel(type=0), keep=True, random_seed=random_seed)
    mts.dump(recap_out_path)
    return mts

def parse_sample_map_lines(lines):
    # expects lines like: SAMPLE_MAP\t<generation>\t<indID>\t<nodeA>\t<nodeB>\t<x>\t<y>
    rows = []
    for ln in lines:
        if not ln.startswith("SAMPLE_MAP\t"):
            continue
        parts = ln.strip().split("\t")
        # ensure proper length
        if len(parts) < 7:
            continue
        gen = int(parts[1]); ind = int(parts[2])
        nodeA = int(parts[3]); nodeB = int(parts[4])
        x = float(parts[5]); y = float(parts[6])
        rows.append({"generation": gen, "individual": ind, "nodeA": nodeA, "nodeB": nodeB, "x": x, "y": y})
    return pd.DataFrame(rows)

def get_individual_locations(ts):
    # returns array of shape (num_individuals, 2) of locations; missing become NaN
    locs = np.asarray(ts.individuals_location, dtype=float)
    return locs

def sample_groups(ts, alive_time=0, ancient_time=None, W=35, w=5, rng_seed=23):
    """Return dict of groups -> list(individual_ids). Alive_time=0 for modern; ancient_time an int generation ago to sample ancient individuals.
    Groups: topleft, topright, bottomleft, bottomright, center, ancient
    """
    rng = np.random.default_rng(rng_seed)
    alive = pyslim.individuals_alive_at(ts, alive_time)
    locs = ts.individuals_location[alive, :]

    groups = {}
    groups['topleft'] = alive[np.logical_and(locs[:,0] < w, locs[:,1] < w)]
    groups['topright'] = alive[np.logical_and(locs[:,0] < w, locs[:,1] > W - w)]
    groups['bottomleft'] = alive[np.logical_and(locs[:,0] > W - w, locs[:,1] < w)]
    groups['bottomright'] = alive[np.logical_and(locs[:,0] > W - w, locs[:,1] > W - w)]
    groups['center'] = alive[np.logical_and(np.abs(locs[:,0] - W/2) < w/2, np.abs(locs[:,1] - W/2) < w/2)]
    if ancient_time is not None:
        old_ones = pyslim.individuals_alive_at(ts, ancient_time)
        if len(old_ones) > 0:
            choose = min(5, len(old_ones))
            groups['ancient'] = rng.choice(old_ones, size=choose, replace=False)
        else:
            groups['ancient'] = np.array([], dtype=int)
    else:
        groups['ancient'] = np.array([], dtype=int)

    return groups

def divergence_between_groups(ts, groups):
    """Compute pairwise mean sequence divergence between and within groups."""
    # produce group-level divergence matrix using ts.divergence
    group_order = list(groups.keys())
    sampled_nodes = []
    for k in group_order:
        nodes_of_group = []
        for ind in groups[k]:
	    # ensure IDs are integer
            ind = int(ind)
            nodes_of_group.extend(list(ts.individual(ind).nodes))
        sampled_nodes.append(nodes_of_group)

    # convert to list of lists of node ids
    div = ts.divergence(sampled_nodes)
    # return DataFrame
    df = pd.DataFrame(div, index=group_order, columns=group_order)
    return df

def pairwise_individual_divergence_and_geo(ts, ind_ids):
    # ind_ids: list/array of individual ids to compare
    ind_nodes = [list(ts.individual(i).nodes) for i in ind_ids]
    # compute pairwise divergence for all unique pairs (i<=j)
    n = len(ind_ids)
    pairs = [(i,j) for i in range(n) for j in range(i, n)]
    div = ts.divergence(ind_nodes, indexes=pairs)
    # geographic distances
    locs = ts.individuals_location
    geog = np.zeros(len(pairs))
    for k, (i,j) in enumerate(pairs):
        xi, yi = locs[ind_ids[i], 0], locs[ind_ids[i], 1]
        xj, yj = locs[ind_ids[j], 0], locs[ind_ids[j], 1]
        geog[k] = sqrt((xi - xj)**2 + (yi - yj)**2)
    return pairs, div, geog

def write_vcf_and_metadata(ts, ind_ids, vcf_path, meta_path):
    # ts.write_vcf requires sample nodes list (nodes per sample). We'll write individuals as "samples"
    sample_nodes = []
    sample_names = []
    metadata_rows = []
    for ind in ind_ids:
        nodes = list(ts.individual(ind).nodes)
        # represent diploid individual using the two node ids
        sample_nodes.append(nodes)
        sample_names.append(f"tsk_{ind}")
        birth_time = ts.node(nodes[0]).time if len(nodes)>0 else np.nan
        metadata_rows.append({
            "vcf_label": f"tsk_{ind}",
            "tskit_id": ind,
            "birth_time_ago": birth_time,
            "x": ts.individual(ind).location[0] if ts.individual(ind).location is not None else np.nan,
            "y": ts.individual(ind).location[1] if ts.individual(ind).location is not None else np.nan
        })
    # flatten node list to tskit.samples order: requires sample node ids
    flat_nodes = [n for pair in sample_nodes for n in pair]
    # tskit.write_vcf accepts individuals argument as sequence of node ids per sample only in python API >= rosetta; 
    # simpler: use ts.samples() selection if necessary: here we will write VCF for chosen node ids via low-level writer
    ts_subset = ts.simplify(flat_nodes, keep_input_roots=True)
    # now write VCF
    with open(vcf_path, "w") as f:
        ts_subset.write_vcf(f, ploidy=2, individual_names=sample_names)
    # write metadata tsv
    pd.DataFrame(metadata_rows).to_csv(meta_path, sep="\t", index=False)
    return

def pairwise_ibd_fraction(ts, ind_ids):
    # compute pairwise IBD fraction genome-wide per individual pair (diploid individuals)
    # We'll use nodes per individual and sum tree intervals where MRCA != NULL for pair's haplotypes.
    node_lists = [list(ts.individual(i).nodes) for i in ind_ids]
    samples = [n for lst in node_lists for n in lst]
    # mapping sample index pairs backup to individuals
    L = ts.sequence_length
    n = len(ind_ids)
    ibd = np.zeros((n,n), dtype=float)
    # for each tree interval check if any haplotype pairs coalesce recently; simpler: treat individuals as merged haplotypes
    for tree in ts.trees():
        interval = tree.interval
        for i in range(n):
            for j in range(n):
                # check if any haplotype of i and any haplotype of j share MRCA at this tree
                shared = False
                for a in node_lists[i]:
                    for b in node_lists[j]:
                        if tree.mrca(a,b) != tskit.NULL:
                            shared = True
                            break
                    if shared:
                        break
                if shared:
                    ibd[i,j] += interval
    ibd = ibd / L
    df = pd.DataFrame(ibd, index=ind_ids, columns=ind_ids)
    return df
