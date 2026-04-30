#!/usr/bin/env python3
"""
Pf_spatial_simulation.py

Comprehensive pipeline to simulate spatial Plasmodium falciparum populations
with SLiM, extract multigenerational pedigrees, recapitate+mutate with msprime,
compute geographic distances, weight pedigree edges by transmission probability,
convert pedigree to directed igraph and compute centrality/bottleneck measures,
and overlay simulated metadata (village/deme, sex, age) on the pedigree.

Usage example:
    python Pf_spatial_simulation.py --outdir out --generations 200 --n_demes 5 --n_samples 100 --test

Requirements:
    - SLiM (callable as 'slim' on PATH)
    - Python packages: pyslim, tskit, msprime, numpy, pandas, igraph, scipy (optional)
"""

import argparse
import os
import json
import math
import subprocess
import tempfile
import random
from typing import List, Tuple, Dict, Any

import numpy as np
import pandas as pd

import pyslim
import tskit
import msprime
import igraph as ig

# --------------------------
# Helper: write spatial SLiM script
# --------------------------
def write_spatial_slim_script(
    script_path: str,
    trees_output: str,
    deme_coords: List[Tuple[float, float]],
    deme_sizes: List[int],
    generations: int,
    base_migration: float = 0.01,
    decay_scale: float = 10.0,
    seed: int = 12345,
    genome_length: int = 5_000_000
) -> str:
    """
    Write a SLiM script implementing a discrete-deme spatial model.

    Returns the path to the script (script_path). Also writes a .meta.json
    containing coords and the computed migration matrix.
    """
    if len(deme_coords) != len(deme_sizes):
        raise ValueError("deme_coords and deme_sizes must have same length")

    Ndeme = len(deme_coords)
    coords = np.array(deme_coords, dtype=float)

    # Distance matrix
    dist_mat = np.zeros((Ndeme, Ndeme), dtype=float)
    for i in range(Ndeme):
        for j in range(Ndeme):
            dx = coords[i, 0] - coords[j, 0]
            dy = coords[i, 1] - coords[j, 1]
            dist_mat[i, j] = math.hypot(dx, dy)

    # Migration matrix (exponential kernel)
    mig_mat = np.zeros_like(dist_mat)
    for i in range(Ndeme):
        for j in range(Ndeme):
            if i == j:
                mig_mat[i, j] = 0.0
            else:
                mig_mat[i, j] = base_migration * math.exp(-dist_mat[i, j] / decay_scale)

    # Save meta JSON
    meta = {
        "deme_coords": deme_coords,
        "deme_sizes": deme_sizes,
        "generations": generations,
        "base_migration": base_migration,
        "decay_scale": decay_scale,
        "dist_matrix": dist_mat.tolist(),
        "migration_matrix": mig_mat.tolist()
    }
    meta_path = os.path.splitext(script_path)[0] + ".meta.json"
    with open(meta_path, "w") as mf:
        json.dump(meta, mf, indent=2)
    print(f"[INFO] Wrote SLiM metadata to {meta_path}")

    # Build SLiM script (keeps pedigrees)
    lines = []
    lines.append("// Spatial discrete-deme SLiM script auto-generated")
    lines.append(f"// Ndeme = {Ndeme}; base_migration = {base_migration}; decay_scale = {decay_scale}")
    lines.append("initialize() {")
    lines.append("    initializeSLiMOptions(keepPedigrees=T);")
    lines.append("    initializeTreeSeq();")
    lines.append(f"    initializeMutationRate(0.0);  // neutral mutations added later with msprime")
    lines.append("    initializeMutationType('m1', 0.5, 'f', 0.0);")
    lines.append("    initializeGenomicElementType('g1', m1, 1.0);")
    lines.append(f"    initializeGenomicElement(g1, 0, {genome_length - 1});")
    lines.append("}")

    # Create demes at generation 1
    lines.append("1 {")
    lines.append(f"    setSeed({seed});")
    for i, size in enumerate(deme_sizes):
        lines.append(f"    sim.addSubpop(paste('p', {i}, sep=''), {int(size)});")
    lines.append("    catn('Created subpopulations: ' + asString(sim.subpopulations));")
    lines.append("}")

    # Set migration rates at generation 1 late()
    lines.append("1 late() {")
    for i in range(Ndeme):
        targets = []
        rates = []
        for j in range(Ndeme):
            if i == j:
                continue
            targets.append(f"sim.subpopulations[{j}]")
            rates.append(f"{mig_mat[i, j]:.8g}")
        if targets:
            targets_str = ", ".join(targets)
            rates_str = ", ".join(rates)
            # Set migration rates for source subpop i
            lines.append(f"    sim.subpopulations[{i}].setMigrationRates(c({targets_str}), c({rates_str}));")
            lines.append(f"    catn('Subpop {i} migration rates set.');")
    lines.append("}")

    # Run until 'generations' then output tree sequence
    lines.append(f"{generations} late() {{")
    lines.append(f"    sim.treeSeqOutput('{trees_output}');")
    lines.append("    catn('Wrote tree sequence and finishing simulation');")
    lines.append("    sim.simulationFinished();")
    lines.append("}")

    script_text = "\n".join(lines)
    with open(script_path, "w") as fh:
        fh.write(script_text)
    print(f"[INFO] Wrote SLiM script to {script_path}")
    return script_path

# --------------------------
# Pedigree extraction (multi-generation)
# --------------------------
def extract_full_pedigree(ts: tskit.TreeSequence, sampled_inds: List[int], max_generations: int = 10) -> pd.DataFrame:
    """
    Extract multigenerational pedigree relationships up to max_generations for each sampled individual.

    Returns a DataFrame with columns:
      child_id, ancestor_id, generations_apart, time_child, time_ancestor
    """
    pedigree_links = []
    visited = set()

    sampled_inds = [int(x) for x in sampled_inds]

    for child_id in sampled_inds:
        queue = [(child_id, 0)]
        while queue:
            current_id, depth = queue.pop(0)
            if depth >= max_generations:
                continue
            parents = pyslim.get_individual_parents(ts, current_id)
            if not parents:
                continue
            for parent_id in parents:
                if parent_id is None or parent_id == tskit.NULL:
                    continue
                key = (child_id, parent_id, depth + 1)
                if key in visited:
                    continue
                visited.add(key)
                try:
                    child = ts.individual(current_id)
                    parent = ts.individual(parent_id)
                    pedigree_links.append({
                        "child_id": int(child.id),
                        "ancestor_id": int(parent.id),
                        "generations_apart": int(depth + 1),
                        "time_child": float(child.time),
                        "time_ancestor": float(parent.time)
                    })
                    queue.append((parent_id, depth + 1))
                except Exception:
                    # Some tree sequences may not expose these individuals as individuals (skip)
                    continue

    pedigree_df = pd.DataFrame(pedigree_links)
    if pedigree_df.shape[0] > 0:
        pedigree_df = pedigree_df.sort_values(["child_id", "generations_apart"]).reset_index(drop=True)
    print(f"[INFO] Extracted {len(pedigree_df)} pedigree links (up to {max_generations} generations).")
    return pedigree_df

# --------------------------
# Direct parent-child pairs
# --------------------------
def extract_direct_transmissions(ts: tskit.TreeSequence) -> pd.DataFrame:
    """
    Extract direct parent->child pairs for all individuals in ts (if metadata exists).
    """
    pairs = []
    # iterate individuals in ts
    for ind in ts.individuals():
        parents = pyslim.get_individual_parents(ts, ind.id)
        if not parents:
            continue
        for p in parents:
            if p is None or p == tskit.NULL:
                continue
            try:
                # times via associated nodes (if present)
                parent_nodes = ts.individual(p).nodes
                child_nodes = ts.individual(ind.id).nodes
                time_parent = float(ts.node(parent_nodes[0]).time) if len(parent_nodes) > 0 else np.nan
                time_child = float(ts.node(child_nodes[0]).time) if len(child_nodes) > 0 else np.nan
                pairs.append({
                    "parent": int(p),
                    "child": int(ind.id),
                    "time_parent": time_parent,
                    "time_child": time_child
                })
            except Exception:
                continue
    df = pd.DataFrame(pairs)
    print(f"[INFO] Extracted {len(df)} direct transmission pairs.")
    return df

# --------------------------
# Recapitate and add neutral mutations (msprime)
# --------------------------
def recapitate_and_mutate(ts: tskit.TreeSequence, recombination_rate: float, mutation_rate: float, Ne: int) -> tskit.TreeSequence:
    """
    Recapitate forward-time tree sequence and add neutral mutations using msprime.
    """
    print("[INFO] Recapitate (pyslim.recapitate) ...")
    ts_recap = pyslim.recapitate(ts, recombination_rate=recombination_rate, ancestral_Ne=Ne)
    print("[INFO] Simulating neutral mutations (msprime.sim_mutations) ...")
    ts_mut = msprime.sim_mutations(ts_recap, rate=mutation_rate, random_seed=None, keep=True)
    print(f"[INFO] Added {ts_mut.num_mutations} mutations.")
    return ts_mut

# --------------------------
# Simplify to sampled individuals
# --------------------------
def simplify_to_sampled(ts: tskit.TreeSequence, sampled_inds: List[int]) -> tskit.TreeSequence:
    """
    Given a tree sequence (likely recapitated+mutated), simplify to nodes corresponding to sampled individuals.
    Each sampled individual -> include all nodes associated with that individual.
    """
    sampled_nodes = []
    for ind in sampled_inds:
        try:
            ind_obj = ts.individual(ind)
            sampled_nodes.extend([n for n in ind_obj.nodes if n != tskit.NULL])
        except Exception:
            # if individual not present in ts (rare), skip
            continue
    if len(sampled_nodes) == 0:
        raise ValueError("No nodes found for sampled individuals.")
    ts_simpl = ts.simplify(samples=sampled_nodes, keep_input_roots=True)
    print(f"[INFO] Simplified to {len(sampled_nodes)} sample nodes.")
    return ts_simpl

# --------------------------
# VCF export helper
# --------------------------
def write_vcf(ts: tskit.TreeSequence, output_path: str, contig_id: str = "chr1"):
    with open(output_path, "w") as fh:
        ts.write_vcf(fh, contig_id=contig_id)
    print(f"[INFO] Wrote VCF to {output_path}")

# --------------------------
# Map sampled individuals to deme (population) and compute pairwise distances
# --------------------------
def compute_pairwise_sample_distances(ts: tskit.TreeSequence, sampled_inds: List[int], meta_json_path: str) -> pd.DataFrame:
    """
    Compute pairwise Euclidean distances between sampled individuals using deme coordinates from meta_json.

    Returns DataFrame with columns: sample_i, sample_j, deme_i, deme_j, distance_km
    """
    # Load metadata
    with open(meta_json_path, "r") as mf:
        meta = json.load(mf)

    deme_coords = meta["deme_coords"]
    # Map sampled individual -> deme by using the population id of one of its nodes
    mapping = []
    for ind in sampled_inds:
        try:
            ind_obj = ts.individual(ind)
            nodes = [n for n in ind_obj.nodes if n != tskit.NULL]
            if len(nodes) == 0:
                # fallback: unknown
                mapping.append((ind, None, None))
                continue
            node0 = ts.node(nodes[0])
            pop_id = node0.population
            # pop_id is integer index (0..Ndeme-1)
            if pop_id is None or pop_id == tskit.NULL:
                mapping.append((ind, None, None))
            else:
                coords = tuple(deme_coords[int(pop_id)])
                mapping.append((ind, int(pop_id), coords))
        except Exception:
            mapping.append((ind, None, None))

    # Create distance DataFrame
    records = []
    for i in range(len(mapping)):
        for j in range(i+1, len(mapping)):
            ind_i, pop_i, coord_i = mapping[i]
            ind_j, pop_j, coord_j = mapping[j]
            if coord_i is None or coord_j is None:
                dist = np.nan
            else:
                dx = coord_i[0] - coord_j[0]
                dy = coord_i[1] - coord_j[1]
                dist = math.hypot(dx, dy)
            records.append({
                "sample_i": int(ind_i),
                "sample_j": int(ind_j),
                "deme_i": pop_i,
                "deme_j": pop_j,
                "distance": float(dist) if dist is not None else np.nan
            })
    df = pd.DataFrame(records)
    print(f"[INFO] Computed pairwise distances for {len(sampled_inds)} samples; pairs: {len(df)}")
    return df

# --------------------------
# Weight pedigree links by transmission probability (exponential decay by generations)
# --------------------------
def weight_pedigree_by_probability(pedigree_df: pd.DataFrame, lambda_gen: float = 0.7) -> pd.DataFrame:
    """
    Given pedigree_df with 'generations_apart' column, add 'transmission_prob'
    computed as exp(-lambda_gen * generations_apart). This is a simple, tunable model.
    """
    df = pedigree_df.copy()
    if "generations_apart" not in df.columns:
        raise ValueError("pedigree_df must contain 'generations_apart' column")
    df["transmission_prob"] = np.exp(-lambda_gen * df["generations_apart"].astype(float))
    # Optionally normalize probabilities per child so they sum to 1 across ancestors of a child
    df["prob_norm_child"] = df.groupby("child_id")["transmission_prob"].transform(lambda x: x / (x.sum() if x.sum() > 0 else 1.0))
    print(f"[INFO] Weighted {len(df)} pedigree links by exponential decay (lambda={lambda_gen}).")
    return df

# --------------------------
# Convert pedigree to igraph directed network and compute centrality
# --------------------------
def pedigree_to_igraph(pedigree_df: pd.DataFrame, metadata_df: pd.DataFrame = None) -> ig.Graph:
    """
    Convert pedigree dataframe to igraph directed graph.
    Nodes: union of child_id and ancestor_id
    Edges: directed ancestor -> child (ancestral -> descendant)
    Attach edge attribute 'generations_apart', 'transmission_prob', 'prob_norm_child'.
    Attach node attributes from metadata_df if provided (metadata_df indexed by node id).
    """
    if pedigree_df.shape[0] == 0:
        raise ValueError("Empty pedigree_df")

    # Build node list
    node_ids = pd.unique(pedigree_df[["child_id", "ancestor_id"]].values.ravel())
    node_id_to_idx = {int(n): idx for idx, n in enumerate(node_ids)}

    # Build edges list (from ancestor -> child)
    edges = []
    gens = []
    tprobs = []
    prob_norms = []
    for _, row in pedigree_df.iterrows():
        a = int(row["ancestor_id"])
        c = int(row["child_id"])
        edges.append((node_id_to_idx[a], node_id_to_idx[c]))  # ancestor -> child
        gens.append(int(row["generations_apart"]))
        tprobs.append(float(row.get("transmission_prob", 1.0)))
        prob_norms.append(float(row.get("prob_norm_child", tprobs[-1] if tprobs else 1.0)))

    g = ig.Graph(directed=True)
    g.add_vertices(len(node_ids))
    g.add_edges(edges)

    # Set vertex attributes: original node id
    g.vs["node_id"] = list(node_ids.astype(int))
    # Add metadata attributes if provided
    if metadata_df is not None and not metadata_df.empty:
        # metadata_df expected indexed by node id or contain column 'node_id'
        if "node_id" in metadata_df.columns:
            md = metadata_df.set_index("node_id")
        else:
            md = metadata_df.copy().set_index(metadata_df.index)
        # Attach attributes where available
        for attr in md.columns:
            vals = []
            for nid in g.vs["node_id"]:
                if int(nid) in md.index:
                    vals.append(md.loc[int(nid), attr])
                else:
                    vals.append(None)
            g.vs[attr] = vals

    # Edge attributes
    g.es["generations_apart"] = gens
    g.es["transmission_prob"] = tprobs
    g.es["prob_norm_child"] = prob_norms

    # Compute centrality measures and attach as vertex attributes
    g.vs["in_degree"] = g.indegree()
    g.vs["out_degree"] = g.outdegree()
    try:
        # betweenness is expensive on very large graphs
        g.vs["betweenness"] = g.betweenness(directed=True)
        g.vs["closeness"] = g.closeness(mode="OUT")
    except Exception:
        g.vs["betweenness"] = [0.0] * g.vcount()
        g.vs["closeness"] = [0.0] * g.vcount()

    print(f"[INFO] Converted pedigree to igraph: vertices={g.vcount()}, edges={g.ecount()}")
    return g

# --------------------------
# Simulate and overlay metadata (deme/village, sex, age) on sampled nodes
# --------------------------
def simulate_sample_metadata(sampled_inds: List[int], ts: tskit.TreeSequence, meta_json_path: str,
                             sex_ratio: float = 0.5, age_mean: float = 20.0, age_sd: float = 10.0) -> pd.DataFrame:
    """
    Simulate simple metadata for sampled individuals: village/deme (from node.population),
    sex (Bernoulli), age (normal truncated at 0), disease status optional.

    Returns pandas DataFrame with columns: node_id, deme, x, y, sex, age
    """
    with open(meta_json_path, "r") as mf:
        meta = json.load(mf)
    deme_coords = meta["deme_coords"]

    records = []
    for ind in sampled_inds:
        try:
            ind_obj = ts.individual(ind)
            nodes = [n for n in ind_obj.nodes if n != tskit.NULL]
            if len(nodes) == 0:
                deme = None
                coord = (None, None)
            else:
                node0 = ts.node(nodes[0])
                pop_id = int(node0.population)
                deme = pop_id
                coord = tuple(deme_coords[pop_id])
        except Exception:
            deme = None
            coord = (None, None)
        sex = "M" if random.random() < sex_ratio else "F"
        age = max(0.1, random.gauss(age_mean, age_sd))
        records.append({
            "node_id": int(ind),
            "deme": deme,
            "x": coord[0],
            "y": coord[1],
            "sex": sex,
            "age": float(age)
        })
    df = pd.DataFrame(records)
    print(f"[INFO] Simulated metadata for {len(df)} sampled individuals.")
    return df

# --------------------------
# Save igraph to GraphML and CSVs
# --------------------------
def save_graph_and_tables(g: ig.Graph, node_meta_df: pd.DataFrame, outdir: str, basename: str):
    os.makedirs(outdir, exist_ok=True)
    graph_path = os.path.join(outdir, f"{basename}.graphml")
    g.write_graphml(graph_path)
    print(f"[INFO] Saved igraph GraphML to {graph_path}")

    # Node table
    node_df = pd.DataFrame({
        "node_index": list(range(g.vcount())),
        "node_id": g.vs["node_id"],
        "in_degree": g.vs["in_degree"],
        "out_degree": g.vs["out_degree"],
        "betweenness": g.vs.get("betweenness", [None]*g.vcount()),
        "closeness": g.vs.get("closeness", [None]*g.vcount())
    })
    # Merge with node_meta_df if provided
    if node_meta_df is not None:
        node_df = node_df.merge(node_meta_df, left_on="node_id", right_on="node_id", how="left")
    node_table_path = os.path.join(outdir, f"{basename}_nodes.csv")
    node_df.to_csv(node_table_path, index=False)
    print(f"[INFO] Saved node attributes to {node_table_path}")

    # Edge table
    edges = []
    for e in g.es:
        src = g.vs[e.tuple[0]]["node_id"]
        tgt = g.vs[e.tuple[1]]["node_id"]
        edges.append({
            "source_node_id": int(src),
            "target_node_id": int(tgt),
            "generations_apart": int(e["generations_apart"]) if "generations_apart" in e.attribute_names() else None,
            "transmission_prob": float(e["transmission_prob"]) if "transmission_prob" in e.attribute_names() else None,
            "prob_norm_child": float(e["prob_norm_child"]) if "prob_norm_child" in e.attribute_names() else None
        })
    edge_df = pd.DataFrame(edges)
    edge_table_path = os.path.join(outdir, f"{basename}_edges.csv")
    edge_df.to_csv(edge_table_path, index=False)
    print(f"[INFO] Saved edge table to {edge_table_path}")

# --------------------------
# Main runner that ties everything together
# --------------------------
def run_spatial_pf_pipeline(
    outdir: str,
    n_demes: int = 5,
    deme_size: int = 200,
    generations: int = 200,
    n_samples: int = 100,
    base_migration: float = 0.02,
    decay_scale: float = 5.0,
    recombination_rate: float = 6.67e-7,
    mutation_rate: float = 1e-8,
    Ne: int = 10000,
    lambda_gen: float = 0.7,
    max_pedigree_generations: int = 8,
    test: bool = False
):
    os.makedirs(outdir, exist_ok=True)

    # Build simple deme layout (grid-ish) for reproducibility
    # Place demes on a roughly circular layout to avoid colocation
    theta = np.linspace(0, 2*math.pi, n_demes, endpoint=False)
    radius = max(1.0, n_demes / 2.0)
    deme_coords = [(radius * math.cos(t), radius * math.sin(t)) for t in theta]
    deme_sizes = [deme_size] * n_demes

    script_path = os.path.join(outdir, "spatial_model.slim")
    trees_path = os.path.join(outdir, "spatial_model.trees")

    write_spatial_slim_script(
        script_path=script_path,
        trees_output=trees_path,
        deme_coords=deme_coords,
        deme_sizes=deme_sizes,
        generations=generations if not test else min(50, generations),
        base_migration=base_migration,
        decay_scale=decay_scale,
        seed=42,
        genome_length=5_000_000
    )

    # Run SLiM
    print("[INFO] Running SLiM ...")
    slim_cmd = ["slim", script_path]
    result = subprocess.run(slim_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print("[ERROR] SLiM failed. stderr (first 2000 chars):")
        print(result.stderr[:2000])
        raise RuntimeError("SLiM run failed")
    else:
        print("[INFO] SLiM finished successfully.")

    # Load original forward-time tree sequence
    ts_orig = pyslim.load(trees_path)
    print(f"[INFO] Loaded forward tree sequence: trees={ts_orig.num_trees}, individuals={ts_orig.num_individuals}")

    # Sample alive individuals at present (time 0)
    alive_inds = pyslim.individuals_alive_at(ts_orig, 0)
    if len(alive_inds) == 0:
        raise RuntimeError("No alive individuals found at time 0 (check SLiM timing).")
    n_sample = min(n_samples, len(alive_inds))
    sampled_inds = list(np.random.choice(alive_inds, size=n_sample, replace=False))
    print(f"[INFO] Selected {n_sample} sampled individuals (IDs: {sampled_inds[:10]}...)")

    # Extract pedigree (multi-generation)
    pedigree_df = extract_full_pedigree(ts_orig, sampled_inds, max_generations=max_pedigree_generations)
    pedigree_file = os.path.join(outdir, "pedigree_multigen.tsv")
    pedigree_df.to_csv(pedigree_file, sep="\t", index=False)
    print(f"[INFO] Saved pedigree (multi-gen) to {pedigree_file}")

    # Direct parent-child pairs
    direct_pairs_df = extract_direct_transmissions(ts_orig)
    direct_pairs_file = os.path.join(outdir, "direct_transmissions.tsv")
    direct_pairs_df.to_csv(direct_pairs_file, sep="\t", index=False)
    print(f"[INFO] Saved direct transmission pairs to {direct_pairs_file}")

    # Recapitate and add neutral mutations
    ts_mut = recapitate_and_mutate(ts_orig, recombination_rate=recombination_rate, mutation_rate=mutation_rate, Ne=Ne)

    # Simplify to sampled individuals (use same individual ids)
    ts_simpl = simplify_to_sampled(ts_mut, sampled_inds)

    # Write VCF
    vcf_path = os.path.join(outdir, "simulated_samples.vcf")
    write_vcf(ts_simpl, vcf_path)

    # Compute simple pairwise genetic distances and save
    genetic_dist = np.array([])  # default empty
    try:
        genetic_dist = compute_simple_genetic_distance(ts_simpl)
    except Exception:
        # fallback using our simple loop if compute_simple_genetic_distance not defined
        pass

    # Compute pairwise geographic distances between sampled individuals
    meta_json_path = os.path.splitext(script_path)[0] + ".meta.json"
    distances_df = compute_pairwise_sample_distances(ts_orig, sampled_inds, meta_json_path)
    distances_file = os.path.join(outdir, "pairwise_sample_distances.tsv")
    distances_df.to_csv(distances_file, sep="\t", index=False)
    print(f"[INFO] Saved pairwise distances to {distances_file}")

    # Weight pedigree by transmission probability
    weighted_pedigree = weight_pedigree_by_probability(pedigree_df, lambda_gen=lambda_gen)
    weighted_pedigree_file = os.path.join(outdir, "weighted_pedigree.tsv")
    weighted_pedigree.to_csv(weighted_pedigree_file, sep="\t", index=False)
    print(f"[INFO] Saved weighted pedigree to {weighted_pedigree_file}")

    # Simulate / attach sample metadata (village/deme, sex, age)
    sample_meta_df = simulate_sample_metadata(sampled_inds, ts_orig, meta_json_path)
    sample_meta_file = os.path.join(outdir, "sample_metadata.tsv")
    sample_meta_df.to_csv(sample_meta_file, sep="\t", index=False)
    print(f"[INFO] Saved sample metadata to {sample_meta_file}")

    # Convert to igraph and compute centrality / bottleneck metrics
    g = pedigree_to_igraph(weighted_pedigree, metadata_df=sample_meta_df)
    graph_base = os.path.join(outdir, "pedigree_graph")
    save_graph_and_tables(g, sample_meta_df, outdir, "pedigree_graph")
    print("[INFO] Graph + tables saved.")

    print("[INFO] Pipeline complete. Outputs in:", outdir)
    return {
        "ts_original": ts_orig,
        "ts_simpl": ts_simpl,
        "pedigree_df": pedigree_df,
        "weighted_pedigree": weighted_pedigree,
        "sample_meta_df": sample_meta_df,
        "distances_df": distances_df,
        "igraph": g,
        "vcf": vcf_path
    }

# --------------------------
# Minimal genetic distance computation (simple proportion of differing sites)
# --------------------------
def compute_simple_genetic_distance(ts: tskit.TreeSequence) -> np.ndarray:
    samples = list(ts.samples())
    n = len(samples)
    if n == 0:
        return np.zeros((0, 0))
    dist = np.zeros((n, n), dtype=float)
    sites = list(ts.variants())
    num_sites = len(sites)
    if num_sites == 0:
        return dist
    for var in sites:
        gts = var.genotypes
        for i in range(n):
            for j in range(i+1, n):
                if gts[i] != gts[j]:
                    dist[i, j] += 1
                    dist[j, i] += 1
    dist = dist / max(1, num_sites)
    print(f"[INFO] Computed simple genetic distance matrix (n={n}, sites={num_sites})")
    return dist

# --------------------------
# CLI
# --------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Spatial Pf simulation pipeline (SLiM + pyslim + msprime + igraph).")
    p.add_argument("--outdir", default="out", help="Output directory")
    p.add_argument("--n_demes", type=int, default=5, help="Number of demes (villages)")
    p.add_argument("--deme_size", type=int, default=200, help="Individuals per deme")
    p.add_argument("--generations", type=int, default=200, help="Forward generations")
    p.add_argument("--n_samples", type=int, default=100, help="Number of final-sample individuals")
    p.add_argument("--base_migration", type=float, default=0.02, help="Base migration rate m0")
    p.add_argument("--decay_scale", type=float, default=5.0, help="Distance decay scale delta")
    p.add_argument("--recombination_rate", type=float, default=6.67e-7, help="Recombination rate per bp per gen")
    p.add_argument("--mutation_rate", type=float, default=1e-8, help="Mutation rate per bp per gen")
    p.add_argument("--Ne", type=int, default=10000, help="Ancestral Ne for recapitation")
    p.add_argument("--lambda_gen", type=float, default=0.7, help="Transmission prob decay per generation")
    p.add_argument("--max_pedigree_generations", type=int, default=8, help="Max generations to trace in pedigree")
    p.add_argument("--test", action="store_true", help="Run a short test (shorter generations)")
    return p.parse_args()

def main():
    args = parse_args()
    if args.test:
        args.generations = min(50, args.generations)
        args.n_samples = min(30, args.n_samples)
    res = run_spatial_pf_pipeline(
        outdir=args.outdir,
        n_demes=args.n_demes,
        deme_size=args.deme_size,
        generations=args.generations,
        n_samples=args.n_samples,
        base_migration=args.base_migration,
        decay_scale=args.decay_scale,
        recombination_rate=args.recombination_rate,
        mutation_rate=args.mutation_rate,
        Ne=args.Ne,
        lambda_gen=args.lambda_gen,
        max_pedigree_generations=args.max_pedigree_generations,
        test=args.test
    )
    # Optionally print summary
    print("\nPipeline produced the following outputs:")
    for k in ["vcf", "pedigree_df", "weighted_pedigree", "sample_meta_df", "distances_df"]:
        if k in res:
            print(" -", k)

if __name__ == "__main__":
    main()
