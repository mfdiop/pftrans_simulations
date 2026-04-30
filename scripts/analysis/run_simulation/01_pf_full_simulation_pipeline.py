#!/usr/bin/env python3

"""
    01_pf_full_simulation_pipeline.py

    Single-population SLiM -> pyslim -> msprime pipeline that:
        - runs SLiM forward simulation (keepPedigrees=T),
        - extracts true parent->child transmission pairs (CSV),
        - recapitate + add neutral mutations (msprime),
        - simplifies to sampled individuals and writes VCF.

Usage:
    python 01_pf_full_simulation_pipeline.py --n_samples 50 --generations 200

Requirements:
    - SLiM installed and on PATH (callable as 'slim')
    - Python packages: pyslim, tskit, msprime, numpy, pandas
"""

import argparse
import subprocess
import os
import sys
import random
import math
import json
import igraph as ig
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Dict, Any, Tuple

import tskit
import msprime
import pyslim
import numpy as np
import pandas as pd

# -------------------------------------------------------------
# Parameter dataclass
# -------------------------------------------------------------

@dataclass
class PfSimParams:
    """Plasmodium falciparum simulation parameters"""
    # Genome parameters
    n_chromosomes: int = 14
    chr_lengths: List[int] = None             # bp, will use realistic sizes if None

    # Population parameters
    Ne: int = 10000                           # Ancient Effective population size
    N0: int = 1000                            # Effective population size at sampling time
    n_samples: int = 100                      # Number of samples to collect

    # Evolutionary parameters
    mutation_rate: float = 1.0e-8             # per bp per generation
    recombination_rate: float = 6.67e-7       # per bp per generation

    # Simulation control
    generations: int = 200                    # Forward simulation generations
    g_sel_start: int = 80
    num_origins: int = 1
    h: float = 0.5
    s: float = 0.0
    sim_relatedness: int = 0

    # Life cycle / sampling
    selfing_rate: float = 0.95                # High inbreeding (haploid with occasional outcrossing)

    # Random seed for reproducibility (if None, will be randomly generated)
    seed: int = None 

    outdir: str = "out"

    def __post_init__(self):
        if self.chr_lengths is None:
            self.chr_lengths = [
                640851, 947102, 1067971, 1200490, 1343557,
                1418242, 1445207, 1472805, 1541735, 1687656,
                2038340, 2271494, 2925236, 3291936
            ]

        # Generate random seed if not provided
        if self.seed is None:
            self.seed = np.random.randint(1, 2**31 - 1)

# -------------------------------------------------------------
# Main simulation class
# -------------------------------------------------------------
class PfSimulation:
    def __init__(self, params: PfSimParams):
        self.params = params
        self.tree_sequences: Dict[int, tskit.TreeSequence] = {}
        self.vcf_paths: Dict[int, str] = {}
        self.ground_truth: Dict[str, Any] = {}
        os.makedirs(self.params.outdir, exist_ok=True)

    # -------------------------
    # SLiM script writer
    # -------------------------
    
    def write_slim_script(self, chr_idx: int, output_file: str, rec_rate: float, trans_name: str, pop_params: int) -> str: # tracking_output: str, 

        """
        Build a SLiM script for one chromosome. Returns the script text.
        tracking_output: path where any SLiM pedigree / logs can be saved
        """

        L = self.params.chr_lengths[chr_idx]

        slim_script = f"""

initialize()
{{
    initializeSLiMOptions(keepPedigrees=T); 
	initializeTreeSeq();                              // F, timeUnit="generations"
	initializeMutationRate(0.0);

	defineConstant("selpos", asInteger({L / 3}));                                 // selection position in bp
	defineConstant("num_origins", {self.params.num_origins});                     // how many genomes contains the selected mutation when selection starts.
	defineConstant("h", {self.params.h});                                         // dominant coefficient
	defineConstant("s", {self.params.s});                                         // selection coefficient
	defineConstant("g_sel_start", {self.params.g_sel_start});                     // time of selected mutation being introduced (generations ago --BACKWARD)
	defineConstant("outid", 1); // idx
	defineConstant("sim_relatedness", {pop_params});             // whether simulate high relatedness

    //defineConstant("TRACK_INTERVAL", max(1, {self.params.generations // 100}) );
	defineConstant("N0", {self.params.N0}); // the effective population size at sampling time
	defineConstant("g_ne_change_start", {self.params.generations});         // Ne change time (generations ago -- BACKWARD)
	defineConstant("slim_total_generations", // time of simulation ended -- forward
		max({self.params.g_sel_start}, {self.params.generations + 1}) );

	initializeMutationType("m1", 0.5, "f", 0.0);                                // neutral
	initializeMutationType("m2", h, "f", s);                                   // balanced
	initializeGenomicElementType("g1", m1, 1);
	initializeGenomicElement(g1, 0, {L - 1});
	initializeRecombinationRate({rec_rate});

	// define global
	defineGlobal("restart_counter", 1);
	defineGlobal("max_restart", 100); // max number of restart
        
}}

// Trick: only run modifyChild in the last twenty generations so its much faster
s0:10 modifyChild()
{{
	if (!sim_relatedness)
		return T; // normal simulation or old generations, always return T

	// use this to let 50% of candidate pass without checking pedigree relatedness
	// so the Ne etimates are more stable
	if (child.index % 2 == 0)
		return T;
	prob = parent1.relatedness(parent2);
	if ((prob > 1.0 / 32) & (prob < 1.0 / 4))
		return T;

	// for unrelated samples return True with a small probability
	else
		return sample(c(T, F), size=1, replace=T, weights=c(1, 999));
}}

1 early()
{{
	sim.addSubpop("p1", {self.params.Ne});
	community.rescheduleScriptBlock(s0, slim_total_generations - 40 + 1);
	community.rescheduleScriptBlock(s1, slim_total_generations - g_ne_change_start + 1);
	community.rescheduleScriptBlock(s2, slim_total_generations - g_sel_start - 1, slim_total_generations); // minus 1 so that it allows the s2 code block the save the state
	community.rescheduleScriptBlock(s3, slim_total_generations + 1, slim_total_generations + 1);
	print(slim_total_generations);

}}

// control population size
s1 300: early()
{{
	t = slim_total_generations - sim.cycle; // generation ago
	Nt = ({self.params.Ne} / N0)^(t / g_ne_change_start) * N0; // calculate Nt
	p1.setSubpopulationSize(asInteger(Nt)); // set new population size 
}}

// condition on selection establishment (not lost)
s2 450: late()
{{
    //cat(sim.subpopulations.individuals[0].pedigreeID); 

	if (sim.cycle == slim_total_generations - g_sel_start - 1 & s != 0.0)
	{{
		sim.treeSeqOutput(paste("state_single_pop_", outid, ".trees", sep=''));
		print(c('saved state:', paste("state_single_pop_", outid, ".trees", sep='')));
		sample(p1.haplosomes, num_origins).addNewDrawnMutation(m2, selpos);
	}}
	else if (sim.cycle >= slim_total_generations - g_sel_start & s != 0)
	{{
		mut = sim.mutationsOfType(m2);
		fixed = sum(sim.substitutions.mutationType == m2);
		need_restart = 0;
		if (fixed)
		{{
			print("selected mutation fixed");
			catn(c("DAF", slim_total_generations - sim.cycle, 1.0), sep='\t');
			community.deregisterScriptBlock(self);
		}}
		else if ((mut.size() != 1) & restart_counter < max_restart)
		{{
			print("selected mutation lost; restarting...");
			sim.readFromPopulationFile(paste("state_single_pop_", outid, ".trees", sep=''));
			setSeed(rdunif(1, 0, asInteger(2^62 - 1)));
			sample(p1.haplosomes, num_origins).addNewDrawnMutation(m2, selpos);
			restart_counter = restart_counter + 1;
		}}
		else
		{{
			catn(c("DAF", slim_total_generations - sim.cycle, sim.mutationFrequencies(p1, mut)), sep='\t');
		}}
	}}
}}

s3 500 late()
{{
    // at end of simulation, write tree sequence with pedigrees
    indv = sim.subpopulations.individuals;
    //catn(indv[0].pedigreeParentIDs);

    // assign tag values to be preserved
    inds = sortBy(sim.subpopulations.individuals, "pedigreeID"); 
    tags = rdunif(size(inds), 0, 100000); 
    inds.tag = tags; 

    // record tag values and pedigree IDs in metadata
    metadataDict = Dictionary("tags", tags, "ids", inds.pedigreeID); 

    // write forward-time tree sequence with pedigrees included
	//sim.treeSeqOutput("{output_file}", simplify=T); Comment this line to explore the one below
    sim.treeSeqOutput("{output_file}", simplify=F, metadata=metadataDict);
    catn("Wrote tree sequence to {output_file}");
	sim.simulationFinished();
	catn(c("restart_counter", restart_counter), sep='\t');

}}

late()
{{
	if (sim.cycle < slim_total_generations)
		catn(c('True_Ne', slim_total_generations - sim.cycle - 1, p1.individualCount), sep='\t');

        // Added this line to remember individuals for pedigree tracking (could not track pedigree before)
        //sim.treeSeqRememberIndividuals(sim.subpopulations.individuals); // remember all individuals for pedigree tracking

    if (sim.cycle >= (slim_total_generations - 10))
        sim.treeSeqRememberIndividuals(sim.subpopulations.individuals);
}}

"""
        return slim_script

    # -------------------------
    # Function To Run SLiM
    # -------------------------

    def run_slim_simulation(self, chr_idx: int, trans_name: str, rec_rate: float, pop_name: str, pop_params: int, output_dir: str = None) -> str:

        """
        Writes and runs the SLiM script for one chromosome.
        Returns the path to the generated tree sequence (.trees).
        """

        if output_dir is None:
            output_dir = self.params.outdir

        os.makedirs(output_dir, exist_ok=True)

        script_file = f"pf_chr{chr_idx+1}.slim"
        output_file = f"pf_chr{chr_idx+1}_{trans_name}_{pop_name}.trees"
        out_file = os.path.join(output_dir, f"pf_chr{chr_idx+1}_output.log")
        err_file = os.path.join(output_dir, f"pf_chr{chr_idx+1}_error.log")
        #tracking_file = os.path.join(output_dir, f"pf_chr{chr_idx+1}_{trans_name}_{pop_name}_pedigree.txt")
 
        print(f"\n--- Writing SLiM script -> {script_file}")
        script_text = self.write_slim_script(chr_idx, output_file, rec_rate, trans_name, pop_params) 
        with open(script_file, 'w') as fh:
            fh.write(script_text)

        # -------------------------
        # run SLiM
        # -------------------------

        # Run SLiM: adjust the command if your SLiM executable is named differently or requires full path
        slim_cmd = ["slim", script_file]
        print("Running SLiM:", " ".join(slim_cmd))
        #result = subprocess.run(slim_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Run the command and capture both stdout and stderr
        with open(out_file, "w") as out, open(err_file, "w") as err:
          result = subprocess.run(f"slim.exe {script_file}", stdout=out, stderr=err, text=True)

        if result.returncode != 0:
            print(" [ERROR] SLiM failed. stderr (first 2000 chars):")
            print(result.stderr[:2000])
            raise RuntimeError(f"\n SLiM failed (exit {result.returncode}). See output above.")
        else:
            print("SLiM completed successfully.\n")

        # SLiM writes the .trees file defined in script; return its path
        if not os.path.exists(output_file):
            raise FileNotFoundError(f"Expected SLiM output not found: {output_file}\n")

        return output_file

    # -------------------------------------------
    #   Define True Transmission Pairs
    # -------------------------------------------
    def _nodes_from_individual_or_nodes(self, ts: tskit.TreeSequence, id_or_node: int) -> List[int]:
       """
       Helper: if id_or_node is an individual id that exists in ts.individuals,
       return the list of node ids for that individual; otherwise assume it's a node id
       and return it as a single-element list.
       """
       # Try to interpret as individual id first
       try:
           ind = ts.individual(id_or_node)
           # If the individual exists, return its nodes (may be empty)
           nodes = [n for n in ind.nodes if n != tskit.NULL]
           if len(nodes) > 0:
               return nodes
       except Exception:
           # not an individual id (or absent); fallthrough to treat as node id
           pass

       # fallback: treat id_or_node as node id
       return [int(id_or_node)]

    # -------------------------------------------
    #   Define Pair Weighted MRCA Time
    # -------------------------------------------
    def pair_weighted_mrca_time(self, ts: tskit.TreeSequence, node_u: int, node_v: int) -> float:
       """
       Compute length-weighted average MRCA time across the tree sequence for node pair (u, v).
       Returns a float (time units used by ts; often generations).
       If no MRCA found in any tree, returns np.inf.
       """
       total_len = 0.0
       weighted_time = 0.0
       for tree in ts.trees():
           # interval length (in same coordinate units as tree positions)
           interval = tree.interval
           interval_len = len(interval)
           if interval_len <= 0:
               continue
           # get MRCA node for this tree
           try:
               mrca = tree.get_mrca(node_u, node_v)
           except Exception:
               # if nodes not present in this tree, skip
               continue
           if mrca == tskit.NULL:
               continue

           # node time from ts (global node table)
           time = float(ts.node(mrca).time)
           weighted_time += time * interval_len
           total_len += interval_len

           if total_len == 0:
               return float("inf")

       return weighted_time / total_len

    # -------------------------------------------
    #   Recent Common Ancestor
    # -------------------------------------------

    def recent_common_ancestor(self, ts: tskit.TreeSequence, id1: int, id2: int, k: float = 5.0) -> Tuple[bool, float]:

        """
        Decide whether two samples (individual id or node id) share a 'recent' common ancestor
        within `k` generations (or time units used in ts.node.time). Returns (is_recent, mrca_time).

        Behavior:
          - If id1/id2 are individual IDs (have .nodes), we consider all node pairs (n1 in id1.nodes, n2 in id2.nodes).
          For polyclonal/multi-node individuals we take the MINIMUM weighted MRCA time across node pairs
          (i.e., conservative: treat as recent if any pair is recent).
          - If id1/id2 are node IDs, they are used directly.
          - MRCA time is computed as the length-weighted average MRCA time across trees:
          weighted_time = sum_tree [ time(mrca_in_tree) * tree_interval_length ] / total_genome_length
          - If no MRCA found (shouldn't happen for connected trees), returns (False, inf).

        Parameters
        ----------
        ts : tskit.TreeSequence
             The (possibly recapitated+mutated) tree sequence.
             id1, id2 : int
             Individual id or node id. If individual id, function will use its node(s).
        k : float
             Threshold in same time units as ts.node.time (usually generations). If weighted_mrca_time <= k, returns True.

        Returns
        -------
        (is_recent: bool, weighted_mrca_time: float)
        """
        # Get node lists for each input (handles individual IDs and node IDs)
        nodes1 = self._nodes_from_individual_or_nodes(ts, id1)
        nodes2 = self._nodes_from_individual_or_nodes(ts, id2)

        # If either has no nodes, cannot compute
        if len(nodes1) == 0 or len(nodes2) == 0:
            return False, float("inf")

        # For multiple node pairs (polyclonal), compute weighted MRCA for each node pair and take min
        best_time = float("inf")

        for n1 in nodes1:
           for n2 in nodes2:
               wt = self.pair_weighted_mrca_time(ts, int(n1), int(n2))
               if wt < best_time:
                   best_time = wt

        # If best_time is inf => no MRCA found
        if not np.isfinite(best_time):
            return False, best_time

        # Determine recentness
        is_recent = (best_time <= float(k))
        return bool(is_recent), float(best_time)
    

    # -----------------------------------------------
    # extract parent->child pairs from SLiM metadata
    # -----------------------------------------------
    def extract_transmission_pairs_from_metadata(self, ts: tskit.TreeSequence) -> pd.DataFrame:
       """
        Use individual.metadata fields 'pedigree_id', 'pedigree_p1', 'pedigree_p2' to build parent->child list.
        Returns DataFrame with columns:
           parent_individual, child_individual, parent_time, child_time, parent_pedigree_id, child_pedigree_id
        Note: parent/child times are taken from the first node associated with the individual (if present).
       """
       # --------------------------------------------
       # Build mapping: pedigree_id -> individual.id
       # --------------------------------------------
       print("[INFO] Build mapping: pedigree_id -> individual.id.")

       ped_to_ind: Dict[int, int] = {}
       for ind in ts.individuals():
            meta = ind.metadata
            # SLiM writes 'pedigree_id' for individuals; skip if not present
            #print("[INFO] SLiM writes 'pedigree_id' for individuals.")
            pid = meta.get("pedigree_id", None)
            if pid is not None:
                ped_to_ind[int(pid)] = int(ind.id)

       rows = []
       for ind in ts.individuals():
            meta = ind.metadata
            #child_ind = int(ind.id)
            child_ind = meta["pedigree_id"]
            child_pid = meta.get("pedigree_id", None)
            # parents stored as 'pedigree_p1' and 'pedigree_p2' (may be None)
            for key in ("pedigree_p1", "pedigree_p2"):
                parent_pid = meta.get(key, None)
                if parent_pid is None:
                    continue
                parent_pid = int(parent_pid)
                parent_ind = ped_to_ind.get(parent_pid, None)
                if parent_ind is None:
                    # parent pid not in map (founder or not present)
                    continue
            # Determine times using node tables (if individual has nodes)
            parent_time = None
            child_time = None
            try:
                parent_nodes = ts.individual(parent_ind).nodes
                if len(parent_nodes) > 0 and parent_nodes[0] != tskit.NULL:
                    parent_time = float(ts.node(parent_nodes[0]).time)
            except Exception:
                parent_time = None
            try:
                child_nodes = ts.individual(child_ind).nodes
                if len(child_nodes) > 0 and child_nodes[0] != tskit.NULL:
                    child_time = float(ts.node(child_nodes[0]).time)
            except Exception:
                child_time = None

            rows.append({
                "parent_individual": parent_ind,
                "child_individual": child_ind,
                "parent_time": parent_time,
                "child_time": child_time,
                "parent_pedigree_id": parent_pid,
                "child_pedigree_id": child_pid if child_pid is not None else None
            })

       df = pd.DataFrame(rows)
       if df.shape[0] == 0:
            print("[WARN] No pedigree links found in metadata. Check that keepPedigrees=T and that SLiM wrote pedigree metadata.")
       else:
            # compute generation difference if times present
            def gen_diff(row):
                if pd.isna(row.parent_time) or pd.isna(row.child_time):
                    return None
                return row.parent_time - row.child_time
            df["generation_diff"] = df.apply(gen_diff, axis=1)

       return df

    def extract_transmission_pairs_from_metadata_v2(self, ts: tskit.TreeSequence) -> pd.DataFrame:
        """
        Extract true parent→child relationships from SLiM pedigree metadata.
        Returns DataFrame with:
            parent_individual, child_individual, parent_time, child_time,
            parent_pedigree_id, child_pedigree_id, generation_diff
        """

        """
        # Build pedigree_id -> individual.id map
        print("[INFO] Building pedigree_id → individual.id map...") 
        ped_to_ind = {}
        for ind in ts.individuals():
            meta = ind.metadata
            pid = meta.get("pedigree_id")
            if pid is not None:
                ped_to_ind[int(pid)] = int(ind.id)

        # Build parent->child pairs
        rows = []
        for child in ts.individuals():
            meta = child.metadata
            child_pid = meta.get("pedigree_id")
            if child_pid is None:
                continue
            child_id = int(child.id)

            # loop over both parents
            for key in ("pedigree_p1", "pedigree_p2"):
                parent_pid = meta.get(key)
                if parent_pid is None:
                    continue
                parent_pid = int(parent_pid)
                parent_ind = ped_to_ind.get(parent_pid)
                if parent_ind is None:
                    continue  # founder or missing

                # get times from first node
                def get_time(ind_id):
                    nodes = ts.individual(ind_id).nodes
                    if len(nodes) > 0 and nodes[0] != tskit.NULL:
                        return float(ts.node(nodes[0]).time)
                    return None

                parent_time = get_time(parent_ind)
                child_time = get_time(child_id)

                # Append row to DataFrame 
                rows.append({
                    "parent_individual": parent_ind,
                    "child_individual": child_id,
                    "parent_time": parent_time,
                    "child_time": child_time,
                    "parent_pedigree_id": parent_pid,
                    "child_pedigree_id": child_pid,
                    "generation_diff": parent_time - child_time if (parent_time and child_time) else None
                })

        # Create DataFrame from rows 
        df = pd.DataFrame(rows)
        if df.empty:
            print("[WARN] No pedigree links found. Check if SLiM wrote pedigree metadata (keepPedigrees=T).")
        else:
            print(f"[INFO] Extracted {len(df)} transmission pairs.")

        # Return DataFrame
        return df
        
        # ----------------------------------------------------------
        # extract parent->child pairs from SLiM metadata (modified)
        # ----------------------------------------------------------
        # Prepare output filenames
        edges_csv = f"{output_prefix}.edges.csv"
        pairs_csv = f"{output_prefix}.pairs.csv"

        # Safety checks
        if ts.num_edges == 0:
            print("[WARN] No edges in tree sequence. Returning empty DataFrame.")
            df_empty = pd.DataFrame(columns=[
                "parent_node","child_node","parent_individual","child_individual",
                "left","right","span_bp","parent_time","child_time","gen_diff",
                "parent_pedigree_id","child_pedigree_id","parent_x","parent_y","child_x","child_y"
            ])
            df_empty.to_csv(edges_csv, index=False)
            df_empty.to_csv(pairs_csv, index=False)
            
            return df_empty

        nodes_table = ts.tables.nodes
        edges_table = ts.tables.edges

        raw_rows = []

        # Iterate edges directly (each edge is a parent->child over [left,right))
        for e in edges_table:
            parent_node = int(e.parent)
            child_node = int(e.child)
            left = float(e.left)
            right = float(e.right)
            span_bp = right - left

            # node times (generations; 0 = present)
            parent_time = float(nodes_table[parent_node].time) if parent_node != tskit.NULL else None
            child_time = float(nodes_table[child_node].time) if child_node != tskit.NULL else None
            gen_diff = None
            
            if parent_time is not None and child_time is not None:
                gen_diff = parent_time - child_time

            # map nodes -> individuals (may be tskit.NULL)
            parent_ind = ts.node(parent_node).individual if parent_node != tskit.NULL else tskit.NULL
            child_ind = ts.node(child_node).individual if child_node != tskit.NULL else tskit.NULL

        """

    # -----------------------------------------------
    # extract parent->child pairs from SLiM metadata
    # -----------------------------------------------

    def extract_transmission_network(self, ts: tskit.TreeSequence, output="true_transmissions", plot_network=True):

        """
        Extract parent->child transmission links from SLiM tree sequence.
        Saves a CSV file with columns:
              parent_node, child_node, left, right, parent_time, child_time, gen_diff
        Optionally plots the transmission network as a directed graph.

        Produces two CSVs:
            1. {output}.csv           - raw edges per genomic segment
            2. {output}_summary.csv   - one row per unique parent–child pair

        Also plots:
            - Transmission network graph
            - Histogram of generation differences

        Parameters
        ----------
        ts : tskit.TreeSequence
                The SLiM tree sequence with pedigree information.
        output : str
                Base name for output CSV file (will add .csv).
        plot_network : bool
                Whether to plot the transmission network graph.
        Returns
        -------
        pd.DataFrame
                Summary DataFrame of parent–child relationships with generation differences.
        """
        output_csv = f"{output}.csv"
        output_summary = f"{output}_summary.csv"

        # Extract nodes and edges
        # ------------------------
        nodes = ts.tables.nodes
        edges = ts.tables.edges
        
        # --------------------------------------------------------
        # Extract parent-child relationships with generation times
        # --------------------------------------------------------
        records = []
        for e in edges:
            parent_time = nodes[e.parent].time   # time of parent node
            child_time = nodes[e.child].time     # time of child node
            gen_diff = parent_time - child_time  # positive means parent older
            
            # Append record for this edge 
            records.append({
                "parent_node": e.parent,
                "child_node": e.child,
                "left": e.left,
                "right": e.right,
                "parent_time": parent_time,
                "child_time": child_time,
                "gen_diff": gen_diff
            })
        
        # Save raw edges to CSV
        df = pd.DataFrame(records)
        df.to_csv(output_csv, index=False)
        print(f"[INFO] Saved {len(df)} raw edges to {output_csv}")
        #print(f" Saved transmission links to {output_csv}")
        
        # Aggregate by parent-child pair (there may be multiple edges covering different genomic regions)

        df_summary = df.groupby(["parent_node", "child_node"]).agg({
            "left": "min",
            "right": "max",
            "parent_time": "first",
            "child_time": "first",
            "gen_diff": "first"
        }).reset_index()
        
        df_summary.to_csv(output_summary, index=False)
        print(f"\n[INFO] Saved {len(df_summary)} unique transmission links to {output_summary}")
        
        # Optional: visualize as directed graph
        if plot_network and len(df_summary) > 0:
            g = ig.Graph(directed=True)
            vertices = list(set(df_summary["parent_node"]).union(df_summary["child_node"]))
            g.add_vertices(vertices)
            g.add_edges(list(zip(df_summary["parent_node"], df_summary["child_node"])))
            
            # Layout by time (younger individuals lower)
#            times = {n: nodes[n].time for n in g.vs["name"]}
#            y_coords = [-times[n] for n in g.vs["name"]]  # invert time for visual clarity
#            layout = [(i % 20, y_coords[i]) for i in range(len(g.vs))]

            times = {int(n): nodes[int(n)].time for n in vertices}
            y_coords = [-times[int(n)] for n in vertices]
            layout = [(i % 20, y_coords[i]) for i in range(len(vertices))]

            fig, ax = plt.subplots(figsize=(12, 8))
            
            ig.plot(
                g,
                layout=layout,
                vertex_size=5,
                vertex_label=None,
                edge_arrow_size=0.5,
                bbox=(1000, 600),
                margin=50,
                #target=None,
                target=ax,  # draw on matplotlib axis
                #output=output_file,
                dpi=600,
                vertex_color="lightblue",
                edge_color="gray"

            )

            # Plot using matplotlib for better control
            plt.title("True Transmission Network (Parent → Child)")
            # Save the plot
            plt.savefig("true_transmission_network.png", dpi=600, bbox_inches="tight")

            # Plot distribution of generation differences
            plt.figure(figsize=(10, 6))
            plt.hist(df_summary["gen_diff"].dropna(), bins=30, color='skyblue', edgecolor='black')
            plt.set_title("Distribution of Generation Differences (Parent - Child)") 
            plt.set_xlabel("Generation Difference")
            plt.set_ylabel("Count")
            plt.savefig("distribution_generation_differences.png", dpi=600, bbox_inches="tight")

            plt.close('all')
        
            #print(f"Saved transmission network plot to {output_file}")
        return df_summary

        # Example usage:
        # df = extract_transmission_network("example_output.trees")

    # --------------------------------------------------
    # Plot ts edges with labels
    # --------------------------------------------------
    def plot_ts_edges_with_labels(self, ts: tskit.TreeSequence, output_file: str, label_by="spatial", ax=None):
        """
        Plot tree sequence edges with node positions labeled by:
        - spatial (x,y) if available in metadata
        - individual IDs otherwise

        Parameters
        ----------
        ts : tskit.TreeSequence
            Tree sequence object.
        label_by : str, optional
            One of {"spatial", "individual"}.
        ax : matplotlib.axes.Axes, optional
            Axis to plot on.
        """

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))

        # --- extract node table for coordinates ---
        nodes = ts.tables.nodes
        node_positions = np.zeros((ts.num_nodes, 2))  # default 2D positions

        for node_id, node in enumerate(nodes):
            # Each node may belong to an individual; get metadata from that
            if node.individual != tskit.NULL:
                ind = ts.individual(node.individual)
                meta = {}
                try:
                    meta = json.loads(ind.metadata.decode())
                except Exception:
                    meta = ind.metadata

                # Fill spatial position if available
                if label_by == "spatial" and "x" in meta and "y" in meta:
                    node_positions[node_id, 0] = meta["x"]
                    node_positions[node_id, 1] = meta["y"]
                elif label_by == "individual":
                    node_positions[node_id, 0] = node_id
                    node_positions[node_id, 1] = 0
            else:
                # fallback: use genomic position or random offset
                node_positions[node_id, :] = [node_id, np.random.rand() * 0.1]

        # --- draw edges ---
        for edge in ts.tables.edges:
            parent_pos = node_positions[edge.parent]
            child_pos = node_positions[edge.child]
            ax.plot(
                [parent_pos[0], child_pos[0]],
                [parent_pos[1], child_pos[1]],
                color="gray", alpha=0.5
            )

        # --- label nodes ---
        for node_id in range(ts.num_nodes):
            label = f"N{node_id}"
            if label_by == "individual" and ts.node(node_id).individual != tskit.NULL:
                label = f"I{ts.node(node_id).individual}"
            ax.text(
                node_positions[node_id, 0],
                node_positions[node_id, 1],
                label, fontsize=8, color="blue"
            )

        ax.set_xlabel("X position")
        ax.set_ylabel("Y position")
        ax.set_title(f"Tree Edges Labeled by {label_by}")
        plt.tight_layout()
        
        return ax

    # --------------------------------------------------
    # Add neutral mutations via recapitation + msprime
    # --------------------------------------------------
    def add_mutations_msprime(self, ts: tskit.TreeSequence, rec_rate: float, pop_params: int, random_seed: int) -> tskit.TreeSequence:
        """
        Recapitate and add neutral mutations using msprime.
        Returns a mutated TreeSequence.
        """
        print("Recapitating with msprime/pyslim...")
        print("[INFO] Recapitate (adding coalescent history before forward sim)...")
        ts_recap = pyslim.recapitate(ts, 
                                     recombination_rate=rec_rate, 
                                     ancestral_Ne=self.params.Ne, 
                                     random_seed=random_seed)

        print("Adding neutral mutations...")
        ts_mut = msprime.sim_mutations(ts_recap, 
                                       rate=self.params.mutation_rate, 
                                       model=msprime.SLiMMutationModel(type=pop_params, 
                                                                       next_id=pyslim.next_slim_mutation_id(ts_recap),
                                                                       ),
                                       random_seed=random_seed, 
                                       keep=True)
        
        print(f"[INFO] Mutations added: {ts_mut.num_mutations}")

        return ts_mut

    # ---------------------------------------------------------------------------------
    # Simplify to sampled individuals (keeps only nodes ancestral to sampled nodes)
    # ---------------------------------------------------------------------------------
    def simplify_to_samples(self, ts_mut: tskit.TreeSequence, sampled_inds: List[int]) -> tskit.TreeSequence:
        """
        Simplify a tree sequence to only include sampled individuals' nodes.
        Expects sampled_inds to be individual ids in the original ts.
        """
        # Convert individuals -> nodes
        sampled_nodes = []
        for ind in sampled_inds:
            try:
               ind_obj = ts_mut.individual(ind)
               # sometimes after recapitation some individuals persist; but if not, fallback: search by pedigree mapping / append associated nodes (haploid -> usually 1 node)
               nodes = [n for n in ind_obj.nodes if n != tskit.NULL]
               sampled_nodes.extend(nodes)
            except Exception:
               # If individual not present in ts_mut.individuals(), try find nodes by matching metadata pedigree_id (rare)
               pass

        if len(sampled_nodes) == 0:
            raise RuntimeError("No sampled nodes found in mutated tree sequence. Sampling mapping failed.")
        ts_simpl = ts_mut.simplify(samples=sampled_nodes, keep_input_roots=True)        

        print(f"Simplified to {len(sampled_nodes)} nodes (samples).")
        return ts_simpl

    # -------------------------
    # VCF export
    # -------------------------
    def extract_vcf(self, ts_simpl: tskit.TreeSequence, chr_idx: int, output_file: str):
        """Export VCF file containing variants for the simplified tree sequence."""

        print(f"Writing VCF to {output_file}...")
        with open(output_file, 'w') as f:
            ts_simpl.write_vcf(f, contig_id=f"chr{chr_idx+1}")
        print(f"Wrote VCF: {output_file}")
        self.vcf_paths[chr_idx] = output_file

    # ---------------------------------------------------
    # Simple genetic distance (pairwise SNP differences)
    # ---------------------------------------------------
    def compute_ground_truth_stats(self, ts_simpl: tskit.TreeSequence) -> np.ndarray:
        """
        Compute a simple pairwise genetic distance matrix (proportion of differing sites).
        Uses ts_simpl.samples() order for nodes.
        """
        samples = list(ts_simpl.samples())
        n = len(samples)
        if n == 0:
            return np.zeros((0, 0))
        dist = np.zeros((n, n), dtype=float)
        sites = list(ts_simpl.variants())
        num_sites = len(sites)
        if num_sites == 0:
            return dist
        for var in sites:
            gts = var.genotypes  # array of 0/1 for samples
            for i in range(n):
                for j in range(i+1, n):
                    if gts[i] != gts[j]:
                        dist[i, j] += 1
                        dist[j, i] += 1

        dist = dist / max(1, num_sites)
        return dist
    
    # -------------------------
    # Summary report
    # -------------------------
    def generate_summary_report(self, output_dir: str, rec_rate: float, scenario_name: str, master_report_path: str = None):
        report_file = os.path.join(output_dir, f"{scenario_name}_summary.txt")

        lines = [
        "="*60,
        f"PLASMODIUM FALCIPARUM SIMULATION SUMMARY — {scenario_name.upper()}",
        "="*60,
        f"Effective population size (Ne): {self.params.Ne}",
        f"Number of samples requested: {self.params.n_samples}",
        f"Generations simulated: {self.params.generations}",
        f"Mutation rate: {self.params.mutation_rate:.2e}",
        f"Recombination rate: {rec_rate:.2e}",
        f"Output directory: {output_dir}",
        "\n"
        ]
        
        with open(report_file, 'w') as fh:
#            fh.write("\n".join(lines))
            for line in lines:
                fh.write(line + "\n")
                print(line)        
    
        # Append to master report if provided
        if master_report_path:
            with open(master_report_path, 'a') as f:
                f.write("\n".join(lines))

        print("\n" + "="*30)
        print(f"✓ Wrote summary report for {scenario_name} to {report_file}")
        print("\n" + "="*30)

    # ------------------------------------------
    # Main pipeline for a single recomb rate
    # ------------------------------------------
    def run_full_simulation(self, n_chromosomes: int = None):

        """
           Execute one simulation: SLiM -> load original ts -> extract pedigree CSV -> recapitate+mutate -> simplify->VCF
        """

        if n_chromosomes is None:
            n_chromosomes = self.params.n_chromosomes

        outdir = self.params.outdir
        os.makedirs(outdir, exist_ok=True)

        # Master summary report
        master_report = os.path.join(outdir, "all_scenarios_summary.txt")

        # Transmission levels (as recombination proxies)
        transmission_levels = {
           'low': 1.0e-7,
           'moderate': 6.67e-7,
           'high': 2.0e-6
        }

        sampling_levels = [self.params.n_samples]  # you can expand the list as needed

        population_structures = {
           'panmictic': 0,
           'structured': 1
        }

        # Loop chromosomes (you can simulate a single chromosome for speed)
        for chr_idx in range(n_chromosomes):
            L = self.params.chr_lengths[chr_idx]
            random_seed = ((self.params.seed + (chr_idx + 1) * L) % (2**32 - 1)) + 1  # different seed per chromosome

            print("\n" + "="*40)
            print(f"Chromosome {chr_idx+1}/{n_chromosomes} (length {L:,} bp)")
            print("="*40)

            for trans_name, rec_rate in transmission_levels.items():
                for n_samples in sampling_levels:
                    for pop_name, pop_params in population_structures.items():
                        scenario_name = f"{trans_name}_n{n_samples}_{pop_name}"
                        scenario_dir = os.path.join(outdir, scenario_name)
                        os.makedirs(scenario_dir, exist_ok=True)

                        print("\n" + "="*30)
                        print(f"\nRunning scenario: {scenario_name}")
                        print(f"  rec_rate = {rec_rate:.2e}")
                        print(f"  n_samples = {n_samples}")
                        print(f"  pop = {pop_name}")
                        print(f"  output_dir = {scenario_dir}")
                        print(f"  seed = {random_seed}")
                        print("\n" + "="*30)

                        # ----------------------------------------------------------------------------------
                        # 1) Run SLiM - this produces original forward-time tree sequence (no recapitation)
                        # ----------------------------------------------------------------------------------
                        trees_file = self.run_slim_simulation(chr_idx, trans_name, rec_rate, pop_name, pop_params, scenario_dir)

                        # ----------------------------------------------------------------------------------
                        # 2) Load original tree sequence (forward-time, with pedigree metadata) (forward-time, keepPedigrees=T)
                        # ----------------------------------------------------------------------------------
                        ts_original = tskit.load(trees_file)
                        print(f"[INFO] Loaded original SLiM tree sequence: trees={ts_original.num_trees}, individuals={ts_original.num_individuals}\n")

                        # -----------------------------------------------------------------------------------------
                        # 3) Identify/choose alive sampled individuals at time 0 (the "present" of the simulation)
                        # -----------------------------------------------------------------------------------------
                        import random
                        random.seed(1234)  # for reproducibility

                        alive_inds = pyslim.individuals_alive_at(ts_original, 0)  # generation 0 = end of sim
                        if len(alive_inds) == 0:
                            raise RuntimeError("No alive individuals found at time 0 - check your SLiM timings.")
                        
                        print(f"[INFO] Found {len(alive_inds)} alive individuals at time 0.\n")
                        # Sample individuals to keep
                        if len(alive_inds) < n_samples // 2:
                            raise RuntimeError(' Not enough alive individuals to sample')

                        #n_sample = min(n_samples, len(alive_inds)) # cannot sample more than alive
                        n_sample = len(alive_inds) // 2  # sample up to half of alive individuals to ensure diversity

                        sampled_inds = list(np.random.choice(alive_inds, size=n_sample, replace=False))
                        print(f"[INFO] Sampling {n_sample} individuals (IDs: {sampled_inds[:10]} ... )\n")

                        # ---------------------------------------------
                        # 4) Extract transmission pairs (parent->child)
                        # ---------------------------------------------
                        print("[INFO] Extract transmission pairs ........ \n")
                        # Find the nodes corresponding to individuals still alive
                        alive_nodes = [n for ind in sampled_inds for n in ts_original.individual(ind).nodes]

                        # Simplify to those nodes (no mutations yet)
                        ts_simplified = ts_original.simplify(samples=alive_nodes, keep_input_roots=True)
                        print(f"[INFO] Simplified TS to sampled individuals: trees={ts_simplified.num_trees}, individuals={ts_simplified.num_individuals}\n")

                        # Extract nodes and edges from the simplified TS
                        print("[INFO] Extracting transmission pairs from metadata ...\n")
                        pedigree_csv = os.path.join(scenario_dir, f"pf_chr{chr_idx+1}_{scenario_name}_true_transmissions")

                        df_pairs = self.extract_transmission_pairs_from_metadata(ts_simplified) 
                        
                        if df_pairs.shape[0] > 0:
                            df_pairs.to_csv(f"{pedigree_csv}.csv", index=False)
                            print(f"[INFO] Saved transmission pairs CSV: {pedigree_csv}\n")
                        else:
                            print("[WARN] No transmission pairs extracted; saved empty DataFrame\n")
                            df_pairs.to_csv(f"{pedigree_csv}.csv", index=False)

                        # ------------------------------------------------
                        # Test second version of transmission pair extraction
                        # ------------------------------------------------
                        print("[INFO] Testing second version of transmission pair extraction ...\n")
                        ground_truth = os.path.join(scenario_dir, f"pf_chr{chr_idx+1}_{scenario_name}_true_transmissions_v2")
                        df_pairs_v2 = self.extract_transmission_pairs_from_metadata_v2(ts_simplified)

                        df_v2 = pd.DataFrame(df_pairs_v2)
                        print(df_v2.head())

                        # Save the results of v2
                        if df_v2.shape[0] > 0:
                            df_v2.to_excel(f"{ground_truth}.xlsx", index=False)
                            print(f"[INFO] Saved transmission pairs CSV (v2): {ground_truth}\n")
                        else:
                            print("[WARN] No transmission pairs extracted in v2; saved empty DataFrame")
                            df_v2.to_excel(f"{ground_truth}.xlsx", index=False)

                        # -------------------------------------------
                        # Optional: extract full transmission network
                        # -------------------------------------------
                        output = os.path.join(scenario_dir, f"pf_chr{chr_idx+1}_{scenario_name}_transmission_network")

                        df = self.extract_transmission_network(ts_simplified, output=output, plot_network=False)

                        # -------------------------------------------
                        # Optional: export nodes and edges tables
                        # -------------------------------------------
                        nodes_csv = os.path.join(scenario_dir, f"pf_chr{chr_idx+1}_{scenario_name}_nodes")
                        nodes = ts_simplified.tables.nodes

                        nodes_df = pd.DataFrame({
                            "node_id": list(range(nodes.num_rows)),
                            "time": nodes.time,
                            "individual": nodes.individual,
                            "flags": nodes.flags,
                            "population": nodes.population,
                            "metadata": [nodes.metadata[i] for i in range(nodes.num_rows)]
                        })

                        print(nodes_df.head())

                        nodes_df.to_csv(f"{nodes_csv}.csv", index=False)
                        nodes_df.to_excel(f"{nodes_csv}.xlsx", index=False)
                        print("\n" + "="*60)
                        print(f"Saved nodes table to {nodes_csv}")

                        #  -------------------------------------------
                        # Export edges table
                        edges_csv  =os.path.join(scenario_dir, f"pf_chr{chr_idx+1}_{scenario_name}_edges")
                        edges = ts_simplified.tables.edges
                        edges_df = pd.DataFrame({
                            "left": [e.left for e in edges],    
                            "right": [e.right for e in edges],
                            "parent": [e.parent for e in edges],
                            "child": [e.child for e in edges]
                        })

                        edges_df.to_csv(f"{edges_csv}.csv", index=False)

                        edges_df.to_excel(f"{edges_csv}.xlsx", index=False)

                        print("\n" + "="*60)
                        print(f"Saved edges table to {edges_csv}")
                        print("\n" + "="*60)

                        # -------------------------------------------
                        # Optional: plot tree edges with labels 
                        # plot by spatial coordinates

                        #self.plot_ts_edges_with_labels(ts_simplified, label_by="spatial")
                        #plt.show()

                        # -------------------------------------------
                        # Optional: extract multi-generation pedigree for sampled individuals
                        # -------------------------------------------

                        # Extract multi-generation pedigree for sampled individuals
                        #pedigree_df = self.extract_full_pedigree(ts_original, sampled_inds, max_generations=10)
                        #pedigree_file = os.path.join(scenario_dir, f"pedigree_chr{chr_idx+1}_{scenario_name}.tsv")
                        #pedigree_df.to_csv(pedigree_file, sep='\t', index=False)
                        #print(f"Saved pedigree (multi-gen) to {pedigree_file}")

                        # Extract and export transmission links
                        #self.export_transmissions_to_csv(ts, out_csv="true_transmission_pairs.csv", max_gen=5)

                        # Also save direct transmission pairs if desired
                        #direct_df = self.extract_direct_transmissions(ts_original)
                        #direct_file = os.path.join(scenario_dir, f"direct_transmissions_chr{chr_idx+1}_{scenario_name}.tsv")
                        #direct_df.to_csv(direct_file, sep='\t', index=False)
                        #print(f"Saved direct transmission pairs to {direct_file}")

                        # --------------------------------------------------------------------
                        # 5) Now Recapitate + add mutations + simplify to sampled individuals
                        # --------------------------------------------------------------------
                        ts_mut = self.add_mutations_msprime(ts_simplified, rec_rate, pop_params, random_seed) # ts_original

                        # ---------------------------------------------
                        # 6) Simplify to sampled individuals + export VCF
                        # Simplify to sampled individuals: use the same sampled_inds (note these are individual ids in original ts)
                        ts_simpl = self.simplify_to_samples(ts_mut, sampled_inds)

                        # -------------------------------------------------
                        # 7) Store outputs: tree sequences, VCF, genetic distance matrix
                        # Store tree sequence for this chromosome / scenario
                        self.tree_sequences[(chr_idx, scenario_name)] = ts_simpl

                        # Export VCF for downstream tools
                        vcf_path = os.path.join(scenario_dir, f"pf_chr{chr_idx+1}_{scenario_name}.vcf")
                        self.extract_vcf(ts_simpl, chr_idx, vcf_path)

                        # Compute simple ground truth genetic distance matrix
                        genetic_dist = self.compute_ground_truth_stats(ts_simpl)
                        mat_file = os.path.join(scenario_dir, f"genetic_distance_chr{chr_idx+1}_{scenario_name}")

                        # Save matrix as .npy
                        np.save(f"{mat_file}.npy", genetic_dist)
                        print(f"\n Saved genetic distance matrix to {mat_file}.npy")

                        # Also save as CSV for easy viewing 
                        np.savetxt(f"{mat_file}.csv", genetic_dist, delimiter=',', fmt='%.2f')
                        print(f"\n Saved genetic distance matrix to {mat_file}.csv ")

                        # Save mapping of sampled individual -> sample node indices for traceability
                        sample_map = {
                            "sampled_individuals": [int(i) for i in sampled_inds],  # ensure list of Python ints
                            "scenario": str(scenario_name),
                            "chr_idx": int(chr_idx),
                            "trees_file": str(trees_file),
                            "vcf_path": str(vcf_path),
                            "pedigree_csv": str(pedigree_csv)
                            }

                        sm_file = os.path.join(scenario_dir, f"sample_map_chr{chr_idx+1}_{scenario_name}.json")
                        
                        with open(sm_file, 'w') as fh:
                            json.dump(sample_map, fh, indent=2)

                        print("\n" + "="*30)
                        print(f"\n Saved sample map to {sm_file}")

                        # Generate overall summary report                    
                        #self.generate_summary_report(self.params.outdir)

                        self.generate_summary_report(scenario_dir, rec_rate, scenario_name, master_report_path=master_report)
                        
                        # Append scenario summary to master CSV
                        summary_data = {
                            "scenario": scenario_name,
                            "rec_rate": rec_rate,
                            "Ne": self.params.Ne,
                            "n_samples": self.params.n_samples,
                            "generations": self.params.generations,
                            "mutation_rate": self.params.mutation_rate,
                            "ts_orig": ts_original,
                            "ts_simpl": ts_simpl,
                            "vcf_path": vcf_path,
                            "pedigree_csv": pedigree_csv,
                            "sample_map": sample_map
                        }

                        summary_csv = os.path.join(self.params.outdir, "scenario_summary_table.csv")
                        df_summary = pd.DataFrame([summary_data])

                        if not os.path.exists(summary_csv):
                            df_summary.to_csv(summary_csv, index=False)
                        else:
                            df_summary.to_csv(summary_csv, index=False, mode='a', header=False)

                        print(f"\n Appended scenario summary to {summary_csv}")

            # End loops
            print("\n All simulations complete.")

# -------------------------------------------------------------
# CLI entrypoint
# -------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Single-pop simulation producing VCF + true transmission pairs with pedigree extraction")
    parser.add_argument("--Ne", type=int, default=10000, help="Ancestral Ne for recapitation")
    parser.add_argument("--N0", type=int, default=1000, help="Population size at end")
    parser.add_argument("--n_samples", type=int, default=100, help="Number of sampled individuals")
    parser.add_argument("--generations", type=int, default=200, help="Forward generations")
    parser.add_argument("--n_chromosomes", type=int, default=1, help="Number of chromosomes to simulate")
    parser.add_argument("--mutation_rate", type=float, default=1e-8, help="Mutation rate per bp per gen")
    parser.add_argument("--recombination_rate", type=float, default=6.67e-7)
    #parser.add_argument("--seed", type=int, help="Random seed")
    parser.add_argument("--outdir", type=str, default="out", help="Output directory")
    parser.add_argument("--test", action="store_true", help="run small test  (smaller genome/gens)")
    args = parser.parse_args()

    # adjust small test case if requested
    if args.test:
        args.n_chromosomes = 1
        args.n_samples = min(20, args.n_samples)
        args.generations = min(100, args.generations)
        print("[INFO] Running small test case with reduced parameters.")

    # Set up simulation parameters data class
    params = PfSimParams(
        N0=args.N0,
        #seed=args.seed,
        n_chromosomes=args.n_chromosomes,
        Ne=args.Ne,
        n_samples=args.n_samples,
        generations=args.generations,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        outdir=args.outdir
    )

    # Run full simulation pipeline
    sim = PfSimulation(params) # create simulation object
    sim.run_full_simulation(n_chromosomes = params.n_chromosomes) # can adjust n_chromosomes here

    print("[DONE] All runs finished.")

if __name__ == "__main__":
    main()
