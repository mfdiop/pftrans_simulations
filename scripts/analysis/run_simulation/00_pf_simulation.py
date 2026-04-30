"""
Plasmodium falciparum simulation pipeline with SLiM and msprime
Generates VCF, IBD segments, and ground truth data
"""

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

@dataclass
class PfSimParams:
    """Plasmodium falciparum simulation parameters"""
    # Genome parameters
    n_chromosomes: int = 14
    chr_lengths: List[int] = None             # bp, will use realistic sizes
    num_origins: int = 1
    
    # Population parameters
    Ne: int = 10000                           # Ancient Effective population size
    N0: int = 1000                            # Effective population size at sampling time
    n_samples: int = 100                      # Number of samples to collect
    
    # Evolutionary parameters
    mutation_rate: float = 1.0e-8             # per bp per generation (Pf specific)  2.5e-9
    recombination_rate: float = 6.67e-7       # per bp per generation      7.4e-8
    h: float = 0.5
    s: float = 0.3
    g_sel_start: int = 80
    sim_relatedness: int = 0
    
    # Life cycle parameters
    generations: int = 200                    # Forward simulation generations
    g_ne_change_start: int = 200              # Forward simulation generations

    selfing_rate: float = 0.95                # High inbreeding (haploid with occasional outcrossing)

    outdir: str = "out"
    
    def __post_init__(self):
        if self.chr_lengths is None:
            # Realistic Pf chromosome lengths (approximate, in bp)
            self.chr_lengths = [
                640851, 947102, 1067971, 1200490, 1343557,
                1418242, 1445207, 1472805, 1541735, 1687656,
                2038340, 2271494, 2925236, 3291936
            ]

    def __recom_init__(self):
        if self.recombination_rate is None:
            # Adjust recombination rates to simulate High, Moderate and Low transmission areas
            self.recombination_rate = [3e-9, 1e-8, 5e-8, 6.67e-7, 2e-7, 1e-6]

class PfSimulation:
    """Main simulation class for P. falciparum"""
    
    def __init__(self, params: PfSimParams):
        self.params = params
        self.tree_sequences = {}
        self.vcf_data = {}
        self.ground_truth = {}
        
    def write_slim_script(self, chr_idx: int, output_file: str, trans_name: str, rec_rate: float, pop_name: str, pop_params: int):
        """Generate SLiM script for one chromosome"""
        
        # Set sequence length based on chromosome size
        L = self.params.chr_lengths[chr_idx]
        
        slim_script = f"""

initialize()
{{
        initializeSLiMOptions(keepPedigrees=T);
	initializeTreeSeq();
	initializeMutationRate(0.0);

	defineConstant("selpos", asInteger({L / 3}));                                 // selection position in bp
	defineConstant("num_origins", {self.params.num_origins});                     // how many genomes contains the selected mutation when selection starts.
	defineConstant("h", {self.params.h});                                         // dominant coefficient
	defineConstant("s", {self.params.s});                                         // selection coefficient
	defineConstant("g_sel_start", {self.params.g_sel_start});                     // time of selected mutation being introduced (generations ago --BACKWARD)
	defineConstant("outid", 1); // idx
	defineConstant("sim_relatedness", {pop_params});             // whether simulate high relatedness

        defineConstant("TRACK_INTERVAL", max(1, {self.params.generations // 100}) );
	defineConstant("N0", {self.params.N0}); // the effective population size at sampling time
	defineConstant("g_ne_change_start", {self.params.g_ne_change_start});         // Ne change time (generations ago -- BACKWARD)
	defineConstant("slim_total_generations", // time of simulation ended -- forward
		max({self.params.g_sel_start}, {self.params.g_ne_change_start + 1}) );

	initializeMutationType("m1", 0.5, "f", 0.0);                                // neutral
	initializeMutationType("m2", h, "f", s);                                   // balanced
	initializeGenomicElementType("g1", m1, 1);
	initializeGenomicElement(g1, 0, {L - 1});
	initializeRecombinationRate({rec_rate});

	// define global
	defineGlobal("restart_counter", 1);
	defineGlobal("max_restart", 100); // max number of restart
        defineConstant("TRACK_FILE", "{tracking_output}");

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
	sim.simulationFinished();
	catn(c("restart_counter", restart_counter), sep='\t');
	sim.treeSeqOutput("{output_file}");
}}

late()
{{
	if (sim.cycle < slim_total_generations)
		catn(c('True_Ne', slim_total_generations - sim.cycle - 1, p1.individualCount), sep='\t');
}}
"""
        return slim_script
    
    def run_slim_simulation(self, chr_idx: int, trans_name: str, rec_rate: float, pop_name: str, pop_params: int, output_dir: str = "pf_simulation_output") -> str:
       
        """Run SLiM for one chromosome"""
        print(f"\n{'='*30}")
        print(f" Simulating chromosome {chr_idx + 1}/{self.params.n_chromosomes}")
        print(f" Length: {self.params.chr_lengths[chr_idx]:,} bp")
        print(f" Transmission = {trans_name} for a rec rate = {rec_rate}")
        print(f" Nbre samples = {n_samples}")
        print(f" Population structure = {pop_name}")
        print(f"{'='*30}")

        # Start Chromosome by 1
        chr_idx = chr_idx + 1

        script_file = f"pf_chr{chr_idx}.slim"
        output_file = f"pf_chr{chr_idx}_{trans_name}_{pop_name}.trees"
        out_file = os.path.join(output_dir, f"pf_chr{chr_idx}_{trans_name}_{pop_name}.out")
        err_file = os.path.join(output_dir, f"pf_chr{chr_idx}_{trans_name}_{pop_name}.err")
        
        # Write SLiM script
        with open(script_file, 'w') as f:
            f.write(self.write_slim_script(chr_idx, output_file, trans_name, rec_rate, pop_name, pop_params))

	# Run the command and capture both stdout and stderr
        with open(out_file, "w") as out, open(err_file, "w") as err:
          result = subprocess.run(f"slim.exe {script_file}", stdout=out, stderr=err, text=True)
        
	# Optional: check the exit code
        if result.returncode == 0:
           print("SLiM finished successfully.")
        else:
           print(f"SLiM failed with code {result.returncode}. Check {err_file}.")
        
        return f"{output_file}"

        # Load the tree sequence
        #ts_original = tskit.load({output_file})
        
        #return f"{output_file}", ts_original

# =================================================================================
    def extract_transmission_network(ts: tskit.TreeSequence, sampled_individuals: list) -> pd.DataFrame:
     """Extract the transmission network from the original SLiM tree sequence.
    
     Parameters
     ----------
     ts : tskit.TreeSequence
         The original SLiM tree sequence (without recapitation).
     sampled_individuals : list
         List of sampled individual IDs (as in the original tree sequence).
    
     Returns
     -------
     pd.DataFrame
         A DataFrame with columns: child, parent, time_child, time_parent.
     """
    # We will traverse the pedigree of the sampled individuals.
    # We are interested in the parent-child relationships that are in the ancestry of the sampled individuals.
    
    # Get all individuals in the tree sequence
    all_individuals = [ts.individual(i) for i in range(ts.num_individuals)]
    
    # We create a set of sampled individual IDs (the ones we are interested in)
    sampled_set = set(sampled_individuals)
    
    # We will collect the parent-child relationships that are in the ancestry of the sampled individuals.
    transmission_events = []
    
    # For each sampled individual, traverse up the pedigree until we hit the root (or until we have gone beyond the generations we care about)
    
    for ind_id in sampled_individuals:
        ind = ts.individual(ind_id)
        # In SLiM, each individual has a parent? In a haploid model, each individual has one parent.
        # The parent is stored in the metadata of the individual? Actually, in pyslim, we can use the pedigree information.
        # We can use the `pyslim` function to get the parents.
        parents = pyslim.get_individual_parents(ts, ind_id)
        for parent_id in parents:
            if parent_id != tskit.NULL:
                # Get the parent individual
                parent = ts.individual(parent_id)
                # Record the transmission event
                transmission_events.append({
                    'child': ind_id,
                    'parent': parent_id,
                    'time_child': ind.time,
                    'time_parent': parent.time
                })
    
    # We can also include the parent-child relationships that are between the sampled individuals and their ancestors.
    # But note: we are only interested in the direct parent-child relationships that are in the ancestry of the sampled individuals.
    
    return pd.DataFrame(transmission_events)

    def extract_transmission_network_v1(self, ts: tskit.TreeSequence) -> pd.DataFrame:
    """Extract true transmission pairs from SLiM pedigree"""
    transmission_pairs = []
    
    for ind in ts.individuals():
        ind_id = ind.id
        metadata = ind.metadata
        # SLiM stores parent information in metadata
        if hasattr(metadata, 'parents') and metadata.parents:
            for parent_id in metadata.parents:
                if parent_id != tskit.NULL:
                    transmission_pairs.append({
                        'donor_id': parent_id,
                        'recipient_id': ind_id,
                        'transmission_time': ind.time,
                        'generations_apart': 1  # Direct transmission
                    })
    
    return pd.DataFrame(transmission_pairs)

# =================================================================================
    
    def add_mutations_msprime(self, ts: tskit.TreeSequence) -> tskit.TreeSequence:
        """Add neutral mutations using msprime"""
        print("Adding mutations with msprime...")
        
        # 1. Recapitate (add ancient history if needed)
        ts_recap = pyslim.recapitate(
            ts,
            recombination_rate=self.params.recombination_rate,
            ancestral_Ne=self.params.Ne
        )
        
        # 2. Add mutations
        ts_mutated = msprime.sim_mutations(
            ts_recap,
            rate=self.params.mutation_rate,
            model = msprime.SLiMMutationModel(type=0),
            random_seed=None,
            keep=True
        )
        
        print(f"  Total mutations: {ts_mutated.num_mutations:,}")
        print(f"  Segregating sites: {len(ts_mutated.tables.sites):,}")
        
        return ts_mutated
    
# =================================================================================
    def simplify_to_samples(self, ts: tskit.TreeSequence) -> tskit.TreeSequence:
        """Simplify tree sequence to sampled individuals"""
        # Get all alive individuals
        alive_inds = pyslim.individuals_alive_at(ts, 0)

	# choose half of the sample size in individuals
        alive_idx = len(alive_inds)

	#if alive_idx < self.params.n_samples:
        #        sys.exit('Not enough alive individuals to sample')
        
        # Randomly sample n individuals
        n_sample = min(self.params.n_samples, len(alive_inds))
        sampled_inds = np.random.choice(alive_inds, size=n_sample, replace=False)
        
        # Get nodes (Pf is haploid, so 1 node per individual)
        sampled_nodes = []
        for ind in sampled_inds:
            sampled_nodes.extend(ts.individual(ind).nodes)
        
        # Simplify
        ts_simplified = ts.simplify(sampled_nodes, keep_input_roots=True)
        
        print(f"  Simplified to {n_sample} samples ({len(sampled_nodes)} nodes)")
        
        return ts_simplified
        
# =================================================================================
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
        
        print(f"\n  Written {ts.num_mutations:,} variants")
        print("\n" + "="*20)

# =================================================================================
    def compute_ground_truth_stats(self, ts: tskit.TreeSequence, chr_idx: int) -> dict:
        """Compute ground truth population genetics statistics"""
        print("\n Computing ground truth statistics...")
         
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
                
        return genetic_dist
    
# =================================================================================
    def run_full_simulation(self, n_chromosomes: int = None): 

        """Run complete simulation pipeline"""
        if n_chromosomes is None:
            n_chromosomes = self.params.n_chromosomes
        
        # Assign output directory to a new variable
        output_dir = self.params.outdir
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Transmission intensity (via recombination rate)
        transmission_levels = {
           'low': 1.0e-7,
           'moderate': 6.67e-7, 
           'high': 2.0e-6
        }
    
        # Sampling intensities
        sampling_levels = [50, 100, 200]  # Number of samples
    
        # Population structures
        population_structures = {
           'panmictic': {'sim_relatedness': 0},
           'structured': {'sim_relatedness': 1}
        }
        
        print("\n" + "="*80)
        
        genetic_distance_matrices = {}
        transmission_network = None
        
        for chr_idx in range(n_chromosomes):

            print("PLASMODIUM FALCIPARUM SIMULATION PIPELINE")
            print("="*80)
            print(f"Parameters:")
            print(f"  Ne = {self.params.Ne:,}")
            print(f"  Samples = {self.params.n_samples}")
            print(f"  Generations = {self.params.generations}")
            print(f"  Selfing rate = {self.params.selfing_rate}")
            print(f"  Mutation rate = {self.params.mutation_rate:.2e}")

            # Generate all combinations
            for trans_name, rec_rate in transmission_levels.items():
               for n_samples in sampling_levels:
                   for pop_name, pop_params in population_structures.items():
                       print(f"  Name = {trans_name}_transmission_{n_samples}_samples_{pop_name}")
                       print(f"  Recombination Rate: rec_rate")
                       print(f"  Number of samples: n_samples")
                       print(f"  Output directory: output_dir/{trans_name}_n{n_samples}_{pop_name}")
                       print("="*30)
                       
                       # Run SLiM
                       trees_file = self.run_slim_simulation(chr_idx, trans_name, rec_rate, pop_name, pop_params, output_dir)            
                       
                       # Load the original tree sequence
                       ts_original = tskit.load(trees_file)
                       print(f" Loaded tree sequence: {ts.num_trees:,} trees, {ts.num_samples} samples\n")

                       if chr_idx == 0:
                          # Get the sampled individuals from the original tree sequence (the ones that are alive at time 0)
                          alive_inds = pyslim.individuals_alive_at(ts_original, 0)
                          n_sample = min(self.params.n_samples, len(alive_inds))
                          sampled_inds = np.random.choice(alive_inds, size=n_sample, replace=False)
                          transmission_network = self.extract_transmission_network(ts_original, sampled_inds)
                
                          # Save the transmission network
                          transmission_network_file = os.path.join(output_dir, "pf_transmission_network.txt")
                          transmission_network.to_csv(transmission_network_file, sep='\t', index=False)
                          print(f"Saved transmission network: {transmission_network_file}")

                      # Now, we proceed with the recapitation and adding mutations on the original tree sequence
                      #ts = self.add_mutations_msprime(ts_original)
                      #ts = self.simplify_to_samples(ts)  # This simplifies to the sampled individuals (the same as in transmission_network)
            
                      # Add mutations with msprime
                      ts = self.add_mutations_msprime(ts_original)
            
                      # Simplify to sampled individuals
                      ts = self.simplify_to_samples(ts)
            
                      # Store tree sequence
                      self.tree_sequences[chr_idx] = ts
            
                      # Export VCF
                      vcf_file = os.path.join(output_dir, f"pf_chr{chr_idx+1}.vcf")
                      self.extract_vcf(ts, chr_idx, vcf_file)
                        
                      # Ground truth statistics
                      genetic_dist = self.compute_ground_truth_stats(ts, chr_idx)
                      self.ground_truth[chr_idx] = {
                         'genetic_distance_matrix': genetic_dist
                      }
            
                      # Store genetic distance matrix separately
                      genetic_distance_matrices[chr_idx] = genetic_dist
                  
# ============================================================================
  
        # Save combined results
        print("\n" + "="*30)
        print("SAVING COMBINED RESULTS")
        print("="*30)
                
        # Ground truth statistics
        print("\n Create Output directory for Matrices")
      
        # Save genetic distance matrices
        for chr_idx, gen_dist in genetic_distance_matrices.items():
            dist_file = os.path.join(output_dir, f"pf_chr{chr_idx+1}_genetic_distance_matrix.txt")
            np.savetxt(dist_file, gen_dist, fmt='%.6f')
            print(f"Saved genetic distance matrix: pf_chr{chr_idx+1}_genetic_distance_matrix.txt")

# ================================================================================================
        
        # Summary report
        self.generate_summary_report(output_dir)
        
        print("\n" + "="*60)
        print("SIMULATION COMPLETE!")
        print(f"Output directory: {output_dir}")
        print("="*60)
        
        return self
    
    def generate_summary_report(self, output_dir: str):
        """Generate summary report"""
        report_file = os.path.join(output_dir, "simulation_summary.txt")
        
        with open(report_file, 'w') as f:
            f.write("="*60 + "\n")
            f.write("PLASMODIUM FALCIPARUM SIMULATION SUMMARY\n")
            f.write("="*60 + "\n\n")
            
            f.write("SIMULATION PARAMETERS\n")
            f.write("-"*15 + "\n")
            f.write(f"Effective population size (Ne): {self.params.Ne:,}\n")
            f.write(f"Number of samples: {self.params.n_samples}\n")
            f.write(f"Generations simulated: {self.params.generations}\n")
            f.write(f"Selfing rate: {self.params.selfing_rate}\n")
            f.write(f"Mutation rate: {self.params.mutation_rate:.2e} per bp per gen\n")
            f.write(f"Recombination rate: {self.params.recombination_rate:.2e} per bp per gen\n")
            f.write(f"Chromosomes simulated: {len(self.tree_sequences)}\n")
            f.write(f"Total genome length: {sum(self.params.chr_lengths[:len(self.tree_sequences)]):,} bp\n\n")
            
            f.write("GENETIC DIVERSITY STATISTICS\n")
            f.write("-"*30 + "\n")
                                  
            f.write("\nOUTPUT FILES\n")
            f.write("-"*20 + "\n")
            f.write("VCF files: pf_chr*.vcf\n")
            f.write("Ground truth statistics: pf_ground_truth_stats.txt\n")
        
        print(f"Summary report: {report_file}")


def main():
    """Example usage"""

    parser = argparse.ArgumentParser(description="Run SLiM with dynamic parameters")


    parser.add_argument("--Ne", type=int, default=10000)                                     # Ancient effective population size
    parser.add_argument("--n_samples", type=int, default=50)
    parser.add_argument("--generations", type=int, default=1000)
    parser.add_argument("--n_chromosomes", type=int, default=14)
    parser.add_argument("--selfing_rate", type=int, default=0.95)
    parser.add_argument("--num_origins", type = int, default = 1),
    parser.add_argument("--h", type = float,  default = 0.5),                                     # dominant coefficient
    parser.add_argument("--s", type = float,  default = 0.3),                                     # selection coefficient
    parser.add_argument("--g_sel_start", type = int, default = 80),                               # time of selected mutation being introduced (generations ago --BACKWARD)
    parser.add_argument("--sim_relatedness", type = int, default = 0),                            # whether simulate high relatedness
    parser.add_argument("--g_ne_change_start", type = int, default = 200),                        # Ne change time (generations ago -- BACKWARD)
    parser.add_argument("--N0", type = int, default = 1000),                                      # Effective population size at sampling time
    parser.add_argument("--mutation_rate", type = float,  default = 1e-8),                        # Mutation rate
    parser.add_argument("--outdir", type = str, default = "out"),
    parser.add_argument("--ibd_max_time", type = int, default = 50),
    parser.add_argument("--ibd_min_length", type = float, default = 20 * (0.01 / 6.67e-7)),

    parser.add_argument("--test", action = "store_true", default = False,
              help = "Use smaller nsam and seqlen for a quick test"),

    args = parser.parse_args()

    params = PfSimParams(Ne=args.Ne, n_samples=args.n_samples, generations=args.generations, n_chromosomes=args.n_chromosomes, num_origins=args.num_origins, h=args.h, s=args.s, g_sel_start=args.g_sel_start, sim_relatedness=args.sim_relatedness, g_ne_change_start=args.g_ne_change_start, N0=args.N0, mutation_rate=args.mutation_rate, ibd_max_time=args.ibd_max_time, ibd_min_length=args.ibd_min_length, outdir=args.outdir)

    sim = PfSimulation(params)
    sim.run_full_simulation()
    
    print("\nSimulation complete! Check 'pf_simulation_output' directory for results.")


if __name__ == "__main__":
    main()


