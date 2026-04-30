
"""
Plasmodium falciparum simulation pipeline with SLiM, pyslim, ts-kit and msprime
Generates VCF, IBD segments, and ground truth data for analysis
Author: [Your Name] 
Date: [Date]

Single-population SLiM -> pyslim -> msprime pipeline that:
 - runs SLiM forward simulation (keepPedigrees=T),
 - Vary recombination rates,
 - Include unrelated individuals,
 - extracts true parent->child transmission pairs (CSV),
 - recapitate + add neutral mutations (msprime),
 - simplifies to sampled individuals and writes VCF.

Usage:
    python3 codes/01_run_simulation_replicates.py --outdir out --n_samples 50 --generations 200

Requirements:
    - SLiM installed and on PATH (callable as 'slim')
    - Python packages: pyslim, tskit, msprime, numpy, pandas

"""

# Import necessary libraries
import argparse
import tskit
import msprime
import pyslim
import numpy as np
import pandas as pd
import subprocess, re
from dataclasses import dataclass
from typing import List, Tuple
import allel
import os
import sys
import shutil

# ===================================================
# Define simulation parameters using dataclass
# ===================================================

@dataclass
class PfSimParams:
    """Plasmodium falciparum simulation parameters"""
    # Genome parameters
    n_chromosomes: int = 14                   # Number of chromosomes
    chr_lengths: List[int] = None             # bp, will use realistic sizes if None
    num_origins: int = 1                      # Number of genomes containing the selected mutation when selection starts.
    
    # Population parameters
    Ne: int = 10000                           # Ancient Effective population size
    N0: int = 1000                            # Effective population size at sampling time
    n_samples: int = 100                      # Number of samples to collect
    
    # Evolutionary parameters
    mutation_rate: float = 1.0e-8             # per bp per generation (Pf specific)  2.5e-9
    h: float = 0.5                            # Dominant coefficient
    s: float = 0.3                            # Selection coefficient
    g_sel_start: int = 80                     # Generations ago when selection starts (backward)
    sim_relatedness: int = 0                  # Whether to simulate high relatedness (1) or not (0)
    
    # Life cycle parameters
    generations: int = 200                    # Forward simulation generations
    g_ne_change_start: int = 200              # Forward simulation generations
    run_id: int = 1

    selfing_rate: float = 0.95                # High inbreeding (haploid with occasional outcrossing)

    outdir: str = "out"                       # Output directory
    
    # Chromosome lengths (if None, use realistic Pf sizes)
    def __post_init__(self):
        if self.chr_lengths is None:
            # Realistic Pf chromosome lengths (approximate, in bp)
            self.chr_lengths = [
                640851, 947102, 1067971, 1200490, 1343557,
                1418242, 1445207, 1472805, 1541735, 1687656,
                2038340, 2271494, 2925236, 3291936
            ]

# ===================================================
#                Main simulation class
# ==================================================
class PfSimulation:
    """
        Main simulation class for P. falciparum using SLiM and msprime
    """

    # Convert rate to clean filename format
    def format_rec_rate(self, rec_rate: float):
        """Convert to format: 1e09, 1e08, 1e07, 1e06 (no minus sign)"""
        exponent = int(abs(round(np.log10(rec_rate))))
        return f"1e{exponent:02d}"
#	  Convert "1e-02" to "1e02"
#         exponent = int(exp_str.split('e')[-1])  # Extract -02
#         positive_exponent = abs(exponent)  # Convert to 02
#         return f"1e{positive_exponent:02d}"
    
    def __init__(self, params: PfSimParams):
        self.params = params
        self.tree_sequences = {}
        self.vcf_data = {}
        self.ground_truth = {}
        
    def write_slim_script(self, chr_idx: int, rec_rate: float, prefix: str, output_file: str): # , tracking_output: str
        """
        Generate SLiM script for one chromosome based on parameters
        """
        # Set sequence legth based on chromosome size
        L = self.params.chr_lengths[chr_idx]

        slim_script = f"""

initialize()
{{
        initializeSLiMOptions(keepPedigrees=T);                                    // keep pedigrees for relatedness calculation
  	initializeTreeSeq(timeUnit="generations");                                 // enable tree sequence recording
  	initializeMutationRate(0.0);                                               // no neutral mutations during SLiM
  
        defineConstant("L", {L});                                                  // chromosome length in bp  
  	defineConstant("selpos", asInteger({L / 3}));                              // selection position in bp
  	defineConstant("num_origins", {self.params.num_origins});                  // how many genomes contains the selected mutation when selection starts.
  	defineConstant("N", {self.params.Ne});                                     // ancient effective population size
  	defineConstant("h", {self.params.h});                                         // dominant coefficient
  	defineConstant("s", {self.params.s});                                         // selection coefficient
  	defineConstant("g_sel_start", {self.params.g_sel_start});                     // time of selected mutation being introduced (generations ago --BACKWARD)
  	defineConstant("r", {rec_rate});                                              // recombinantion rate
        defineConstant("outid", '{prefix}');
  	defineConstant("sim_relatedness", {self.params.sim_relatedness});             // whether simulate high relatedness
  
  	defineConstant("N0", {self.params.N0});                                       // the effective population size at sampling time
  	defineConstant("g_ne_change_start", {self.params.g_ne_change_start});         // Ne change time (generations ago -- BACKWARD)
  	defineConstant("slim_total_generations",                                      // time of simulation ended -- forward
  		max({self.params.g_sel_start}, {self.params.g_ne_change_start + 1}) );
  
  	initializeMutationType("m1", 0.5, "f", 0.0);                                  // neutral mutations
  	initializeMutationType("m2", h, "f", s);                                      // balanced selection mutation
  	initializeGenomicElementType("g1", m1, 1);                                    // neutral genomic element
  	initializeGenomicElement(g1, 0, {L - 1});                                     // whole chromosome is neutral
  	initializeRecombinationRate(r);                                               // recombination rate 
  
  	// define global
  	defineGlobal("restart_counter", 1);                                           // restart counter
  	defineGlobal("max_restart", 100);                                             // max number of restart

}}


// -----------------------------
// Relatedness-Controlled Mating
// -----------------------------
//reproduction(NULL, "F") {{
//    subpop.addCrossed(parent1, parent2);
//}}

// Run modifyChild only in late generations (to control cost)
//late() {{
//    if (sim.generation >= sim.generation - 20)
//        sim.modifyChild("relatedness_control");
//}}

// Relatedness filtering to control diversity
s0:10 modifyChild() 
{{
    if (!exists("sim_relatedness") | !sim_relatedness)
        return T; // skip control if not enabled

    // Retrieve parental relatedness (pedigree-based)
    rel = child.parent1.relatedness(child.parent2);

    // Allow ~50% matings to occur freely for Ne stability
    if (child.index % 2 == 0)
        return T;

    // Accept if parents are unrelated enough (rel < 1/64 ≈ distant)
    if (rel < 1.0 / 64.0)
        return T;

    // Accept if moderately related (1/64 ≤ rel ≤ 1/8) with some probability
    if (rel >= 1.0 / 64.0 & rel <= 1.0 / 8.0)
        return (runif(1) < 0.5);  // 50% chance to accept moderate relatedness

    // Reject highly related matings (rel > 1/8) most of the time
    return (runif(1) < 0.02);  // 2% chance for inbred matings
}}

1 early()
{{
  	sim.addSubpop("p1", N);                                                 // initial population size
//        sim_relatedness = T;  // enable relatedness filtering
  	community.rescheduleScriptBlock(s0, slim_total_generations - 50 + 1);                  // only run modifyChild in the last 50 generations
  	community.rescheduleScriptBlock(s1, slim_total_generations - g_ne_change_start + 1);  // reschedule population size control
  	community.rescheduleScriptBlock(s2, slim_total_generations - g_sel_start - 1, slim_total_generations); // minus 1 so that it allows the s2 code block the save the state
  	community.rescheduleScriptBlock(s3, slim_total_generations + 1, slim_total_generations + 1);           // end simulation after slim_total_generations
  	print(slim_total_generations);

}}

// control population size over time 
s1 300: early()
{{
  	t = slim_total_generations - sim.cycle;                           // generation ago (backward)
  	Nt = (N / N0)^(t / g_ne_change_start) * N0;       // calculate Nt based on exponential growth 
  	p1.setSubpopulationSize(asInteger(Nt));                         // set new population size 
}}

// condition on selection establishment (not lost) 
s2 450: late()
{{
    // introduce selected mutation at selpos 
  	if (sim.cycle == slim_total_generations - g_sel_start - 1 & s != 0.0)
  	{{
  		sim.treeSeqOutput(paste("state_", outid, ".trees", sep=''));          // save state before introducing selected mutation         
  		print(c('saved state:', paste("state_", outid, ".trees", sep='')));   // print message
  		sample(p1.haplosomes, num_origins).addNewDrawnMutation(m2, selpos);              // add selected mutation to num_origins random haplotypes
  	}}
  	else if (sim.cycle >= slim_total_generations - g_sel_start & s != 0)                 // track selected mutation frequency
  	{{
  		mut = sim.mutationsOfType(m2);                                                   // get all selected mutations 
  		fixed = sum(sim.substitutions.mutationType == m2);                               // check if selected mutation is fixed
  		need_restart = 0;                                                                // flag for restart
  		if (fixed)
  		{{
  			print("selected mutation fixed");                                            // print message 
  			catn(c("DAF", slim_total_generations - sim.cycle, 1.0), sep='\t');
  			community.deregisterScriptBlock(self);                                      // stop further execution
  		}}
  		else if ((mut.size() != 1) & restart_counter < max_restart)                     // check if selected mutation is lost
  		{{
  			print("selected mutation lost; restarting...");
  			sim.readFromPopulationFile(paste("state_", outid, ".trees", sep=''));  // reload saved state
  			setSeed(rdunif(1, 0, asInteger(2^62 - 1)));
  			sample(p1.haplosomes, num_origins).addNewDrawnMutation(m2, selpos);              // reintroduce selected mutation
  			restart_counter = restart_counter + 1;                                           // increment restart counter
  		}}
		  else
		  {{
			  catn(c("DAF", slim_total_generations - sim.cycle, sim.mutationFrequencies(p1, mut)), sep='\t'); // print derived allele frequency
		  }}
	  }}
}}

s3 500 late()
{{
  	sim.simulationFinished();                                                          // end simulation
  	catn(c("restart_counter", restart_counter), sep='\t');                             // print number of restarts
  	sim.treeSeqOutput("{output_file}");                                                // save final tree sequence
}}

late()
{{
  	if (sim.cycle < slim_total_generations)                                           // print true Ne every TRACK_INTERVAL generations
  		catn(c('True_Ne', slim_total_generations - sim.cycle - 1, p1.individualCount), sep='\t');   // print true Ne
}}

"""
        return slim_script
    
    # ========================================================================
    #                       Run SLiM simulation 
    # =======================================================================
    def run_slim_simulation(self, chr_idx: int, rec_rate: float, prefix: str, output_dir: str = "pf_simulation_output") -> str:
        """
            Run SLiM for one chromosome and return output tree sequence file path
        """
        print("="*30)
        print(f"Simulating chromosome {chr_idx + 1}/{self.params.n_chromosomes}")
        print(f"Length: {self.params.chr_lengths[chr_idx]:,} bp")
        print("="*30)

        # Start Chromosome by 1
        chr = chr_idx + 1

        # Prepare file names
        script_file = os.path.join(output_dir, f"{prefix}_chr{chr}.slim")         # write SLiM script to this file
        output_file = f"{prefix}_chr{chr}.trees"                                  # output tree sequence file
        # tracking_output = f"pf_chr{chr}.txt"
        out_file = os.path.join(output_dir, f"{prefix}_chr{chr}.out")             # SLiM stdout log
        err_file = os.path.join(output_dir, f"{prefix}_chr{chr}.error")           # SLiM stderr log
        
        # Write SLiM script to file
        with open(script_file, 'w') as f:
            f.write(self.write_slim_script(chr_idx, rec_rate, prefix, output_file)) # , tracking_output
        
        # Run SLiM
        #os.system(f"slim.exe {script_file}")

        SLIM_EXEC = "/mnt/scratch/fadel/PhD/Objective1/methods_evaluation/build/slim "
	# Run the command and capture both stdout and stderr
        with open(out_file, "w") as out, open(err_file, "w") as err:
            # with shell=True:
            result = subprocess.run(f"{SLIM_EXEC} {script_file}", shell=True, stdout=out, stderr=err, text=True)
        
	    # Optional: check the exit code and print a message
        # Check the exit code
        if result.returncode == 0:
            print("SLiM finished successfully.")
        else:
            print(f"SLiM failed with code {result.returncode}. Check {err_file}.")
            sys.exit(1)  # Exit on failure
        
        # Return output tree sequence file
        return f"{output_file}"
    
    # =======================================================================
    #   Simplify (extract alived) to sampled individuals to simulate haploid
    # =======================================================================
    def simplify_to_samples(self, ts: tskit.TreeSequence) -> tskit.TreeSequence:
        """Simplify tree sequence to sampled individuals"""
        
        print(" Simplifying to sampled alived individuals ...")
        
        # Get all alive individuals
        alive_inds = pyslim.individuals_alive_at(ts, 0)

	    # choose half of the sample size in individuals
        alive_idx = len(alive_inds)
        
        if alive_idx < self.params.n_samples // 2:
            #raise RuntimeError('Not enough alive individuals to sample') # More specific error
            #raise Exception('Not enough alive individuals to sample')    # General error or generic exception
            sys.exit('Not enough alive individuals to sample')            # Exit program with message

        # Determine number of haploid samples to draw
        n_sampled = self.params.n_samples // 2                             # Haploid samples, so half the individuals
        # n_sample = min(self.params.n_samples, alive_idx)                # In case alive individuals are less than requested samples
        print(f" Alive individuals: {alive_idx}, Sampling {n_sampled} individuals")

        # Randomly sample n haploid individuals from alive individuals
        sampled_inds = np.random.choice(alive_inds, size=n_sampled, replace=False)
        
        # Get nodes (Pf is haploid, so 1 node per individual) for sampled individuals
        sampled_nodes = []
        for ind in sampled_inds:
            sampled_nodes.extend(ts.individual(ind).nodes)
        
        # Simplify tree sequence to sampled nodes
        ts_simplified = ts.simplify(sampled_nodes, keep_input_roots=True)
        
        # Summary of simplified tree sequence
        print(f"  Simplified to {n_sampled} samples ({len(sampled_nodes)} nodes)")
        
        # Return simplified tree sequence
        return ts_simplified
    
    # End of simplify_to_samples
    
    # ========================================================================
    #   Post-simulation processing with tskit, pyslim and msprime
    # ========================================================================
    def recapitate_add_mutations(self, ts: tskit.TreeSequence, recombination_rate: float) -> tskit.TreeSequence:

        """Add neutral mutations using msprime"""
        # ===============================================================================================
        # Recapitate tree sequence with given recombination rate and Ne (add ancient history if needed)
        # ===============================================================================================
        print("Recapitate tree sequence with pyslim ...")
        
        ts_recap = pyslim.recapitate(
            ts,
            recombination_rate=recombination_rate,
            ancestral_Ne=self.params.Ne,
            random_seed=42
        )
        
        # Summary of recapitated tree sequence
        print(f"  Recapitated TS: {ts_recap.num_trees:,} trees, {ts_recap.num_samples} samples")

        # =======================================
        # Simplify recapitated tree sequence
        cur_sample_nodes = list(ts_recap.samples())
        sts2 = ts_recap.simplify(cur_sample_nodes)

        # delete existing sites then mutate
        nsites = sts2.num_sites
        if nsites > 0:
            site_list = list(range(nsites))
            ts_nosites = sts2.delete_sites(site_list)
        else:
            ts_nosites = sts2

        # ========================================
        # Add mutations with msprime
        # ========================================
        print("Adding mutations with msprime ...")
        ts_mutated = msprime.sim_mutations(            # ts_recap,
            ts_nosites,
            rate=self.params.mutation_rate,
            model = msprime.SLiMMutationModel(type=0),
            random_seed=42,
            keep=True
        )
        
        # Summary of mutated tree sequence
        print(f"  Total mutations: {ts_mutated.num_mutations:,}")
        print(f"  Segregating sites: {len(ts_mutated.tables.sites):,}")
        
        return ts_mutated, ts_recap 
    
    # End of add_mutations_msprime

    # ========================================================================
    #                     Export VCF file
    # ========================================================================
    def extract_vcf(self, ts: tskit.TreeSequence, chr_idx: int, output_file: str):
        """Export VCF file"""

        print(f"Writing VCF to {output_file}...")
        print("\n" + "="*20)
        
        with open(output_file, 'w') as f:
            ts.write_vcf(
                f,
                contig_id=f"chr{chr_idx + 1}",
                individuals=range(ts.num_individuals)
            )
        
        # Summary
        print(f"\n  Written {ts.num_mutations:,} variants")
        print("\n" + "="*20)

    # End of extract_vcf

    # ========================================================================
    #              Compute ground truth statistics
    # ========================================================================
    def compute_ground_truth_stats(self, ts: tskit.TreeSequence, chr_idx: int) -> dict:
        """Compute ground truth population genetics statistics"""
        
        print("\n Computing ground truth statistics...")
        
        stats = {
            'chr': chr_idx + 1,
            'length': int(ts.sequence_length),
            'n_samples': ts.num_samples,
            'n_sites': ts.num_sites,
            'n_mutations': ts.num_mutations,
        }
        
        # Genetic diversity (pi)
        stats['pi'] = ts.diversity(mode='site')
        
        # Tajima's D
        stats['tajimas_d'] = ts.Tajimas_D(mode='site')
        
        # Segregating sites
        stats['segregating_sites'] = ts.segregating_sites(mode='site')


        # Compute SFS from final tree sequence
        # This is the "true" SFS after simplification to samples
        sfs_from_ts = np.zeros(ts.num_samples + 1, dtype=int)
        
        for var in ts.variants():
            ac = np.sum(var.genotypes)  # Allele count
            
            # Ensure allele count is within bounds (account for numerical issues)
            if ac >= 0 and ac <= ts.num_samples:
                sfs_from_ts[int(ac)] += 1
            else:
                # Count variants outside expected range (shouldn't happen but safeguard)
                print(f"  Warning: Allele count {ac} out of bounds for {ts.num_samples} samples")
        
        stats['sfs_from_samples'] = sfs_from_ts[1:-1].tolist()  # Exclude 0 and n (non-segregating)

        # Pairwise relatedness (genetic distance)
        samples = ts.samples()
        n = len(samples)
        genetic_dist = np.zeros((n, n))
        
        for var in ts.variants():
            gts = var.genotypes
            for i in range(n):
                for j in range(i+1, n):
                    if gts[i] != gts[j]:
                        genetic_dist[i, j] += 1
                        genetic_dist[j, i] += 1
        
        # Normalize by number of sites
        if ts.num_sites > 0:
            genetic_dist /= ts.num_sites
        
        stats['mean_genetic_distance'] = np.mean(genetic_dist[np.triu_indices(n, k=1)])
        stats['median_genetic_distance'] = np.median(genetic_dist[np.triu_indices(n, k=1)])
        
        print(f"  π = {stats['pi']:.6f}")
        print(f"  Tajima's D = {stats['tajimas_d']:.4f}")
        print(f"  Segregating Sites = {stats['segregating_sites']:.4f}")
        print(f"  Mean genetic distance = {stats['mean_genetic_distance']:.3f}")
        
        return stats, genetic_dist
    
    # ================================================================================================
    #                    Generate summary report
    # ================================================================================================
    
    def generate_summary_report(self, rec_rate: float, run_id: int, prefix: str, output_dir: str):
        """
        Generate summary report for the simulation run
        
        Args:
            rec_rate: Recombination rate used
            run_id: Run identifier
            prefix: File name prefix
            output_dir: Output directory path
        """
        report_file = os.path.join(output_dir, f"{prefix}_simulation_summary.txt")
        
        with open(report_file, 'w') as f:
            f.write("="*60 + "\n")
            f.write("PLASMODIUM FALCIPARUM SIMULATION SUMMARY\n")
            f.write("="*60 + "\n\n")
            
            f.write("SIMULATION PARAMETERS\n")
            f.write("-"*30 + "\n")
            f.write(f"Run ID: {run_id}\n")
            f.write(f"Effective population size (Ne): {self.params.Ne:,}\n")
            f.write(f"Current population size (N0): {self.params.N0:,}\n")
            f.write(f"Number of samples: {self.params.n_samples}\n")
            f.write(f"Generations simulated: {self.params.generations}\n")
            f.write(f"Selfing rate: {self.params.selfing_rate}\n")
            f.write(f"Mutation rate: {self.params.mutation_rate:.2e} per bp per gen\n")
            f.write(f"Recombination rate: {rec_rate:.2e} per bp per gen\n")
            f.write(f"Selection coefficient (s): {self.params.s}\n")
            f.write(f"Dominance coefficient (h): {self.params.h}\n")
            f.write(f"Selection start generation: {self.params.g_sel_start}\n")
            f.write(f"Ne change start generation: {self.params.g_ne_change_start}\n")
            f.write(f"Chromosomes simulated: {len(self.tree_sequences)}\n")
            f.write(f"Total genome length: {sum(self.params.chr_lengths[:len(self.tree_sequences)]):,} bp\n\n")
            
            f.write("GENETIC DIVERSITY STATISTICS\n")
            f.write("-"*30 + "\n")
            
            # Iterate through ground truth data
            for chr_key, data in sorted(self.ground_truth.items()):
                stats = data['stats']
                f.write(f"\nChromosome {stats['chr']}:\n")
                f.write(f"  Length: {stats['length']:,} bp\n")
                f.write(f"  Samples: {stats['n_samples']}\n")
                f.write(f"  Segregating sites: {stats['n_sites']:,}\n")
                f.write(f"  Total mutations: {stats['n_mutations']:,}\n")
                f.write(f"  Nucleotide diversity: {stats['pi']:.6f}\n")
                f.write(f"  Tajima's D: {stats['tajimas_d']:.4f}\n")
                f.write(f"  Mean genetic distance: {stats['mean_genetic_distance']:.6f}\n")
                f.write(f"  Median genetic distance: {stats['median_genetic_distance']:.6f}\n")
            
            # Calculate genome-wide averages
            all_stats = [data['stats'] for data in self.ground_truth.values()]
            f.write("\nGENOME-WIDE AVERAGES\n")
            f.write("-"*30 + "\n")
            
            total_sites = sum(s['n_sites'] for s in all_stats)
            total_mutations = sum(s['n_mutations'] for s in all_stats)
            avg_pi = np.mean([s['pi'] for s in all_stats])
            avg_tajimas_d = np.mean([s['tajimas_d'] for s in all_stats])
            avg_genetic_dist = np.mean([s['mean_genetic_distance'] for s in all_stats])
            
            f.write(f"Total segregating sites: {total_sites:,}\n")
            f.write(f"Total mutations: {total_mutations:,}\n")
            f.write(f"Average pi: {avg_pi:.6f}\n")
            f.write(f"Average Tajima's D: {avg_tajimas_d:.4f}\n")
            f.write(f"Average genetic distance: {avg_genetic_dist:.6f}\n\n")
                
            f.write("OUTPUT FILES\n")
            f.write("-"*30 + "\n")
            f.write(f"VCF files: {prefix}_chr*.vcf\n")
            f.write(f"Tree sequences: {prefix}_chr*_mutated.trees\n")
            f.write(f"Ground truth statistics: {prefix}_chr*_ground_truth.txt\n")
            f.write(f"Genetic distance matrices: {prefix}_chr*_genetic_distance.txt\n")
            f.write(f"Summary report: {report_file}\n\n")
            
            f.write("="*60 + "\n")
            f.write("END OF REPORT\n")
            f.write("="*60 + "\n")
        
        print(f"\n✓ Summary report saved: {report_file}")

    # ========================================================================
    #
    # ========================================================================
    def run_full_simulation(self, rec_rate: float, run_id: int, n_chromosomes: int = None):
        """
        Run complete simulation pipeline for all chromosomes
            1. Run SLiM simulation
            2. Post-process with pyslim and msprime
            3. Export VCF and compute ground truth statistics
            4. Save combined results
            5. Generate summary report

        """

        # Determine number of chromosomes to simulate
        if n_chromosomes is None:
            n_chromosomes = self.params.n_chromosomes

        # Assign output directory to a new variable
        output_dir = self.params.outdir

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        rec = self.format_rec_rate(rec_rate)
        prefix = f"run{run_id}_rec{rec}"

        # Print simulation parameters
        print("\n" + "="*80)
        print("     PLASMODIUM FALCIPARUM SIMULATION PIPELINE")
        print("="*80)
        print(f"Parameters:")
        print(f"  Output file name = {prefix}")
        print(f"  Ne = {self.params.Ne:,}")
        print(f"  Samples = {self.params.n_samples}")
        print(f"  Generations = {self.params.generations}")
        print(f"  Selfing rate = {self.params.selfing_rate}")
        print(f"  Mutation rate = {self.params.mutation_rate:.2e}")
        print(f"  Recombination rate = {rec_rate}")
        print(f"  Chromosomes to simulate = {n_chromosomes}")
        print("="*30)
        
        all_stats = []
        genetic_distance_matrices = {}
        
        for chr_idx in range(n_chromosomes):
            # Run SLiM
            trees_file = self.run_slim_simulation(chr_idx, rec_rate, prefix, output_dir)
            
            # Load tree sequence
            ts = tskit.load(trees_file)
            print(f"Loaded tree sequence: {ts.num_trees:,} trees, {ts.num_samples} samples")
            
            # Simplify to sampled individuals
            ts_simplified = self.simplify_to_samples(ts)
            
            # Add mutations with msprime
            ts_mutate, ts_recap = self.recapitate_add_mutations(ts_simplified, recombination_rate=rec_rate)

            # Save Recapitate tree
            recap_ts_file = os.path.join(output_dir, f"{prefix}_chr{chr_idx+1}_recap.trees")
            ts_recap.dump(recap_ts_file)
            
            # Store tree sequence
            self.tree_sequences[chr_idx] = ts_mutate

            # save mutated tree sequence
            mutated_ts_file = os.path.join(output_dir, f"{prefix}_chr{chr_idx+1}_mutated.trees")
            ts_mutate.dump(mutated_ts_file)
            
            # Export VCF
            vcf_file = os.path.join(output_dir, f"{prefix}_chr{chr_idx+1}.vcf")
            self.extract_vcf(ts_mutate, chr_idx, vcf_file)
            
            # Ground truth statistics
            stats, genetic_dist = self.compute_ground_truth_stats(ts_mutate, chr_idx)
            self.ground_truth[chr_idx+1] = {
                'stats': stats,
                'genetic_distance_matrix': genetic_dist
            }
            
            all_stats.append(stats)

            # Store genetic distance matrix separately
            #genetic_distance_matrices[chr_idx+1] = genetic_dist
            
            # Save Genetic Distance MAtrix
            dist_file = os.path.join(output_dir, f"{prefix}_chr{chr_idx+1}_genetic_distance.txt")
            np.savetxt(dist_file, genetic_dist, fmt='%.6f')
            print(f"Saved genetic distance matrix: {prefix}_chr{chr_idx+1}_genetic_distance.txt")

            # Move generated tree to output directory
            destination_path = os.path.join(output_dir, os.path.basename(trees_file))
            try:
            # Move the file
                shutil.move(trees_file, destination_path)
                print(f"File '{trees_file}' successfully moved to '{destination_path}'")
            except FileNotFoundError:
                print(f"Error: Source file '{trees_file}' not found.")
            except Exception as e:
                print(f"An error occurred: {e}")

        # ======================
        # Save combined results
        # ======================
        print("\n" + "="*30)
        print("SAVING COMBINED RESULTS")
        print("="*30)
                
        # Ground truth statistics
        print("\n Create Output directory for stats")
        stats_output = os.path.join(output_dir, f"{prefix}_chr{chr_idx+1}_ground_truth.txt")

        stats_df = pd.DataFrame(all_stats)
        stats_df.to_csv(stats_output, sep='\t', index=False)
        print(f"Saved ground truth stats: {stats_output}")

        # Save genetic distance matrices
        #for chr_idx, gen_dist in genetic_distance_matrices.items():
        #    dist_file = os.path.join(output_dir, f"{prefix}_chr{chr_idx+1}_genetic_distance.txt")
        #    np.savetxt(dist_file, gen_dist, fmt='%.6f')
        #    print(f"Saved genetic distance matrix: {prefix}_chr{chr_idx+1}_genetic_distance.txt")

        # Generate simulation summary report
        print("\n" + "="*30)
        print("GENERATING SUMMARY REPORT")
        print("="*30)
        self.generate_summary_report(rec_rate=rec_rate, run_id=run_id, prefix=prefix, output_dir=output_dir)
        
        print("\n" + "="*80)
        print(f"✓ SIMULATION COMPLETE - RUN {run_id}")
        print(f"✓ All results saved to: {output_dir}")
        print("="*80)

# ===================================================
#                Main function
# ===================================================
def main():
    """
    Parse command-line arguments and run simulation pipeline with specified parameters
    1. Parse arguments
    2. Initialize PfSimulation with parameters
    3. Run full simulation pipeline
    4. Print completion message 
    """

    parser = argparse.ArgumentParser(description="Run SLiM with dynamic parameters")

    parser.add_argument("--Ne", type=int, default=10000, help="Effective Population Size")
    parser.add_argument("--n_samples", type=int, default=100, help="Number of samples to collect")
    parser.add_argument("--generations", type=int, default=1000, help="Forward simulation generations")
    parser.add_argument("--n_chromosomes", type=int, default=14, help="Number of chromosomes to simulate")
    parser.add_argument("--selfing_rate", type=float, default=0.95, help="Selfing rate (high inbreeding)")
    parser.add_argument("--num_origins", type=int, default=1, help="Number of genomes containing the selected mutation")
    parser.add_argument("--h", type=float, default=0.5, help="Dominant coefficient")
    parser.add_argument("--s", type=float, default=0.3, help="Selection coefficient")
    parser.add_argument("--g_sel_start", type=int, default=80, help="Time of selected mutation being introduced")
    parser.add_argument("--sim_relatedness", type=int, default=0, help="Simulate high relatedness")
    parser.add_argument("--g_ne_change_start", type=int, default=200, help="Ne change time (generations ago)")
    parser.add_argument("--N0", type=int, default=1000, help="Effective population size at sampling time")
    parser.add_argument("--mutation_rate", type=float, default=1e-8, help="Mutation rate")
    parser.add_argument("--rec_rate", nargs="+", type=float, default=[1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6], 
                        help="Recombination rates to sweep")
    parser.add_argument("--run_id", type=int, default=1, help="Run ID")
    parser.add_argument("--outdir", type=str, default="out", help="Output directory")
    parser.add_argument("--test", action="store_true", default=False,
                        help="Use smaller nsam and seqlen for a quick test")
    
    # Check if no arguments provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Test mode - run and exit
    if args.test:
        print("Running in test mode with smaller parameters...")
        params = PfSimParams(
            Ne=1000,
            N0=100,
            n_samples=20,
            generations=100,
            n_chromosomes=1,
            selfing_rate=args.selfing_rate,
            num_origins=args.num_origins,
            h=args.h,
            s=args.s,
            g_sel_start=args.g_sel_start,
            sim_relatedness=args.sim_relatedness,
            g_ne_change_start=args.g_ne_change_start,
            mutation_rate=args.mutation_rate,
            outdir="test_output"
        )
        sim = PfSimulation(params)
        sim.run_full_simulation(rec_rate=6.66666667e-7, run_id=1, n_chromosomes=1)
        print("\nTest simulation complete! Check 'test_output' directory for results.")
        return  # FIXED: Exit after test mode
    
    # Main simulation loop
    for i, r in enumerate(args.rec_rate):  # FIXED: rec_rates -> rec_rate
        print("\n" + "="*60)
        print(f"[RUN {i+1}/{len(args.rec_rate)}] recombination rate = {r}")
        
        # Print parsed arguments
        print("\nParsed arguments:")
        for arg, value in vars(args).items():
            print(f"  {arg}: {value}")

        # Initialize simulation parameters
        params = PfSimParams(
            Ne=args.Ne,
            n_samples=args.n_samples,
            generations=args.generations,
            n_chromosomes=args.n_chromosomes,
            selfing_rate=args.selfing_rate,  # Added missing parameter
            num_origins=args.num_origins,
            h=args.h,
            s=args.s,
            g_sel_start=args.g_sel_start,
            sim_relatedness=args.sim_relatedness,
            g_ne_change_start=args.g_ne_change_start,
            N0=args.N0,
            mutation_rate=args.mutation_rate,
            outdir=args.outdir
        )

        sim = PfSimulation(params)
        sim.run_full_simulation(rec_rate=r, run_id=i+1, n_chromosomes=args.n_chromosomes)

        print(f"\nSimulation {i+1} complete!")

    print(f"\n[DONE] All {len(args.rec_rate)} runs finished. Check '{args.outdir}' directory for results.")


if __name__ == "__main__":
    main()
