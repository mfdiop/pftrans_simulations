#!/usr/bin/env python3

"""
Multi-population malaria transmission simulation with migration
Explores the parameter space of: recombination, bottlenecks, mutation, sampling, and outbreak size
"""

import argparse
import sys
import os
import json
import itertools
from dataclasses import dataclass, asdict
from typing import List, Dict, Tuple
import gzip
import numpy as np
import pandas as pd
from pathlib import Path


@dataclass
class TransmissionParams:
    """
    Parameters defining transmission biology and population structure
    """
    # Core biological triad
    rec_rate: float  			 # Recombination rate (proxy for within-host diversity)
    bottleneck_size: int  		 # Number of parasites transmitted
    est: float  			     # Expected substitutions per transmission
    chr_length: int = 640851

    # Population structure
    n_populations: int = 3  		# Number of demes/populations
    migration_rate: float = 0.01  	# Migration rate between populations

    # Population sizes and sampling
    Ne: int = 10000  			        # Effective population size per deme
    outbreak_size: int = 2000  		   # Total number of infections
    sampling_proportion: float = 0.2  	# Proportion of infections sampled

    # Evolutionary parameters
    generations: int = 200  		# Forward simulation time
    mutation_rate: float = 1e-8  	# Per-site per-generation

    # Selection parameters (optional - can simulate neutral or selected scenarios)
    s: float = 0.0  			     # Selection coefficient (0 = neutral)
    h: float = 0.5  			     # Dominance coefficient
    g_sel_start: int = 80  		     # Generation when selection starts (0 = no selection)

    # Selfing/inbreeding
    selfing_rate: float = 0.0  		# Selfing rate (0 = random mating, 0.95 = high inbreeding)

    # Simulation metadata
    replicate_id: int = 5
    scenario_id: str = ""
    outdir: str = "sim_migration"    # "migration_sims"


class ParameterSpace:
    """
    Defines the full parameter space for the transmission simulation study
    """

    # A. Core biological triad
    RECOMBINATION_RATES = {
        'very_low': 1e-9,      # All clonal
        'low': 1e-8,           # Almost clonal
        'medium': 1e-7,        # Moderate recombination
        'high': 1e-6          # Frequent recombination
#        'very_high': 1e-5      # Very high recombination / high COI
    }

    BOTTLENECK_SIZES = {
        'tight': 1,      	# Single parasite (severe bottleneck)
        'medium': 5,     	# Few parasites
        'loose': 20      	# Many parasites (relaxed bottleneck)
    }

    EXPECTED_SUBS_PER_TRANSMISSION = {
        'none': 0.0,         	# No mutations
        'low': 0.1,     	# Very rare mutations
        'moderate': 0.5      	# Low-moderate
#        'moderate': 1.0,     	# Moderate mutation load
#        'high': 2.0          	# High mutation accumulation
    }

    # B. Sampling and scale
    SAMPLING_PROPORTIONS = [0.05, 0.1, 0.2]  # 0.3, 0.6, 0.9
    OUTBREAK_SIZES = 2000 #[2500, 3500, 5000]

    # C. Population structure
    N_POPULATIONS = 3 #[2, 3, 4]  	# Number of interconnected populations
    MIGRATION_RATES = {
        'low': 0.001,      		# Rare migration
        'medium': 0.01,    		# Moderate gene flow
        'high': 0.05       		# High connectivity
    }

    @classmethod
    def get_est_from_mutation_rate(cls, est: float, genome_length: int = 640851, scenario: str = "typical") -> float:    #23e6
        """
        Convert expected substitutions per transmission to mutation rate
        EST = mutation_rate × genome_length × generations_per_transmission
        Assume ~10-15 days per transmission cycle in mosquito + human

       Convert EST to mutation rate with detailed malaria life cycle accounting.

       Malaria transmission cycle generations:

       HUMAN PHASE:
          - Liver stage: 7-10 days = 1 generation (hepatocyte invasion → merozoite release)
          - Blood stage: Typical 10-day infection with 48h cycles = 5 generations
          (Day 0-2, 2-4, 4-6, 6-8, 8-10)

       MOSQUITO PHASE:
          - Gametocyte → Ookinetes: 1-2 generations
          - Oocyst sporogony: 8-10 generations (massive replication)
          - Sporozoite migration: 1-2 generations

        TOTAL: ~15-25 generations depending on scenario
        """
        # Malaria generation time approximation:
        # ~1 day in mosquito, parasite replication every 48h in human

        if est == 0:
            return 0

        generation_scenarios = {
            "rapid": 18,      # Short infections, efficient transmission
            "typical": 22,    # Standard estimate (1 + 5 + 10 + 6)
            "prolonged": 28,  # Long infections, multiple mosquito feeds
            "conservative": 25  # Rounded conservative estimate
        }

        if scenario not in generation_scenarios:
            raise ValueError(f"scenario must be one of {list(generation_scenarios.keys())}")

        generations_per_transmission = generation_scenarios[scenario]

        mutation_rate = est / (genome_length * generations_per_transmission)
        # Log detailed breakdown for transparency
        if scenario == "typical":
            print(f"EST breakdown: {est} subs / (23Mbp × {generations_per_transmission} gens)")
            print("  Typical generations: Liver=1, Blood=5, Mosquito=16")

        # Validate the calculated mutation rate is biologically plausible
        if mutation_rate > 1e-6:
            print(f"WARNING: Calculated mutation rate {mutation_rate:.2e} seems high")
        elif mutation_rate < 1e-10:
            print(f"WARNING: Calculated mutation rate {mutation_rate:.2e} seems very low")

        return mutation_rate

    @classmethod
    def get_core_scenarios(cls) -> List[Dict]:
        """
        Generate all 5 × 3 × 5 = 75 core biological scenarios
        """
        scenarios = []
        scenario_id = 0

        for rec_name, rec_rate in cls.RECOMBINATION_RATES.items():
            for bottle_name, bottle_size in cls.BOTTLENECK_SIZES.items():
                for est_name, est_value in cls.EXPECTED_SUBS_PER_TRANSMISSION.items():
                    scenario_id += 1
                    scenarios.append({
                        'scenario_id': f"S{scenario_id:03d}_{rec_name}_{bottle_name}_{est_name}",
                        'rec_rate': rec_rate,
                        'bottleneck_size': bottle_size,
                        'est': est_value,
                        'rec_name': rec_name,
                        'bottle_name': bottle_name,
                        'est_name': est_name
                    })

        return scenarios

    @classmethod
    def get_full_design(cls, n_replicates: int = 10) -> List[TransmissionParams]:
        """
        Generate full factorial design: core scenarios × sampling × outbreak size × replicates
        Total: 48 cores × 3 sampling × 1 outbreak sizes × n_replicates
        """
        core_scenarios = cls.get_core_scenarios()
        all_params = []

        for scenario in core_scenarios:
            for sampling_prop in cls.SAMPLING_PROPORTIONS:
                for replicate in range(1, n_replicates + 1):
                    # Calculate mutation rate from EST
                    mutation_rate = cls.get_est_from_mutation_rate(scenario['est'])

                    params = TransmissionParams(
                        rec_rate=scenario['rec_rate'],
                        bottleneck_size=scenario['bottleneck_size'],
                        est=scenario['est'],
                        sampling_proportion=sampling_prop,
                        outbreak_size=OUTBREAK_SIZES,
                        mutation_rate=mutation_rate,
                        replicate_id=replicate,
                        scenario_id=scenario['scenario_id']
                    )
                    all_params.append(params)

        return all_params

    @classmethod
    def get_migration_scenarios(cls, n_replicates: int = 10) -> List[TransmissionParams]:
        """
        Generate migration-specific scenarios varying population number and migration rate
        Subset of core scenarios with migration parameter sweep
        """
        # Select representative core scenarios (not all 75)
        representative_scenarios = [
            ('very_low', 'tight', 'moderate'),       # Low rec, tight bottleneck, moderate mutation
            ('low', 'medium', 'moderate'),           # Balanced scenario
            ('medium', 'loose', 'low'),              # High diversity scenario
            ('medium', 'medium', 'moderate'),        # High diversity scenario
            ('high', 'medium', 'low'),               # High diversity scenario
        ]

        all_params = []
        scenario_id = 0

        for rec_name, bottle_name, est_name in representative_scenarios:
            rec_rate = cls.RECOMBINATION_RATES[rec_name]
            bottle_size = cls.BOTTLENECK_SIZES[bottle_name]
            est_value = cls.EXPECTED_SUBS_PER_TRANSMISSION[est_name]
            mutation_rate = cls.get_est_from_mutation_rate(est_value)

            for mig_name, mig_rate in cls.MIGRATION_RATES.items():
                for sampling_prop in cls.SAMPLING_PROPORTIONS:
                    for replicate in range(1, n_replicates + 1):
                        scenario_id += 1
                        params = TransmissionParams(
                            rec_rate=rec_rate,
                            bottleneck_size=bottle_size,
                            est=est_value,
                            n_populations=N_POPULATIONS,
                            migration_rate=mig_rate,
                            sampling_proportion=sampling_prop,
                            outbreak_size=OUTBREAK_SIZES,
                            mutation_rate=mutation_rate,
                            replicate_id=replicate,
                            scenario_id=f"MIG{scenario_id:04d}_{rec_name}_{bottle_name}_npop{n_pops}_m{mig_name}"
                        )
                        all_params.append(params)

        return all_params


class SimulationDesign:
    """
    Manages the simulation design and execution
    """

    def __init__(self, design_type: str = "full", n_replicates: int = 10,
                 outdir: str = "sim_migration"):
        self.design_type = design_type
        self.n_replicates = n_replicates
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)

        # Generate parameter combinations
        if design_type == "full":
            self.params_list = ParameterSpace.get_full_design(n_replicates)
        elif design_type == "migration":
            self.params_list = ParameterSpace.get_migration_scenarios(n_replicates)
        elif design_type == "core":
            # Core scenarios only, no sampling/outbreak size variation
            core = ParameterSpace.get_core_scenarios()
            self.params_list = []
            for scenario in core:
                for rep in range(1, n_replicates + 1):
                    mutation_rate = ParameterSpace.get_est_from_mutation_rate(scenario['est'])
                    params = TransmissionParams(
                        rec_rate=scenario['rec_rate'],
                        bottleneck_size=scenario['bottleneck_size'],
                        est=scenario['est'],
                        mutation_rate=mutation_rate,
                        replicate_id=rep,
                        scenario_id=scenario['scenario_id']
                    )
                    self.params_list.append(params)
        else:
            raise ValueError(f"Unknown design type: {design_type}")

        print(f"Generated {len(self.params_list)} parameter combinations")

    def save_design(self, filename: str = "simulation_design.json"):
        """Save the full simulation design to JSON"""
        design_data = {
            'design_type': self.design_type,
            'n_replicates': self.n_replicates,
            'n_scenarios': len(self.params_list),
            'parameters': [asdict(p) for p in self.params_list]
        }

        output_file = self.outdir / filename
        with open(output_file, 'w') as f:
            json.dump(design_data, f, indent=2)

        print(f"Saved simulation design to {output_file}")
        return output_file

    def generate_slurm_array(self, slurm_template: str = "run_migration_sim.slurm",
                            array_size: int = 20):
        """
        Generate SLURM array job script
        Split scenarios into manageable array jobs
        """
        n_scenarios = len(self.params_list)
        n_arrays = (n_scenarios + array_size - 1) // array_size

        slurm_script = f"""#!/bin/bash

#SBATCH --job-name=malaria_migration
#SBATCH --output=logs/migration_%A_%a.out
#SBATCH --error=logs/migration_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-{n_scenarios}%{array_size}
#SBATCH --partition=main


# Run Python simulation for this array task
python3 codes/malaria_transmission_study/malaria_migration_full_pipeline_v1.3.py \\
    --run-from-design migration_sims/simulation_design.json \\
    --task-id ${{SLURM_ARRAY_TASK_ID}} \\
    --outdir {self.outdir}

echo "Task ${{SLURM_ARRAY_TASK_ID}} complete"
"""

        slurm_file = self.outdir / slurm_template
        with open(slurm_file, 'w') as f:
            f.write(slurm_script)

        os.chmod(slurm_file, 0o755)
        print(f"Generated SLURM script: {slurm_file}")
        print(f"Submit with: sbatch {slurm_file}")

        return slurm_file

    def get_summary_stats(self) -> Dict:
        """Generate summary statistics of the design"""
        if not self.params_list:
            return {}

        stats = {
            'total_scenarios': len(self.params_list),
            'n_replicates': self.n_replicates,
            'rec_rates': sorted(set(p.rec_rate for p in self.params_list)),
            'bottleneck_sizes': sorted(set(p.bottleneck_size for p in self.params_list)),
            'est_values': sorted(set(p.est for p in self.params_list)),
            'sampling_props': sorted(set(p.sampling_proportion for p in self.params_list)),
            'outbreak_sizes': sorted(set(p.outbreak_size for p in self.params_list)),
            'n_populations': sorted(set(p.n_populations for p in self.params_list)),
            'migration_rates': sorted(set(p.migration_rate for p in self.params_list))
        }

        return stats

    def print_design_summary(self):
        """Print a summary of the simulation design"""
        stats = self.get_summary_stats()

        print("\n" + "="*70)
        print("SIMULATION DESIGN SUMMARY")
        print("="*70)
        print(f"Design type: {self.design_type}")
        print(f"Total parameter combinations: {stats['total_scenarios']}")
        print(f"Replicates per scenario: {self.n_replicates}")
        print(f"\nParameter ranges:")
        print(f"  Recombination rates: {len(stats['rec_rates'])} levels")
        print(f"  Bottleneck sizes: {stats['bottleneck_sizes']}")
        print(f"  EST values: {stats['est_values']}")
        print(f"  Sampling proportions: {stats['sampling_props']}")
        print(f"  Outbreak sizes: {stats['outbreak_sizes']}")
        print(f"  Number of populations: {stats['n_populations']}")
        print(f"  Migration rates: {stats['migration_rates']}")
        print("="*70 + "\n")


def run_single_scenario(params: TransmissionParams, dry_run: bool = False):
    """
    Execute a single simulation scenario

    This integrates with SLiM + msprime + downstream analysis
    """
    import subprocess
    import json
    from pathlib import Path

    print(f"\n{'='*60}")
    print(f"Running scenario: {params.scenario_id}")
    print(f"Replicate: {params.replicate_id}")
    print(f"{'='*60}")

    # Print parameters
    print("\nParameters:")
    for key, value in asdict(params).items():
        print(f"  {key}: {value}")

    if dry_run:
        print("\n[DRY RUN] Would execute simulation here")
        return

    # Create output directory
    scenario_dir = Path(params.outdir) / params.scenario_id / f"rep_{params.replicate_id:0d}"
    scenario_dir.mkdir(parents=True, exist_ok=True)

    # Save parameters
    param_file = scenario_dir / "parameters.json"
    with open(param_file, 'w') as f:
        json.dump(asdict(params), f, indent=2)

    print(f"\nOutput directory: {scenario_dir}")

    # =========================================================================
    # STEP 1: Run SLiM simulation with migration
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 1: Running SLiM forward simulation")
    print("="*60)

    slim_script = find_slim_script()  # Find appropriate SLiM script
    slim_seed = params.replicate_id + hash(params.scenario_id) % 1000000

    # Get the full path to the SLiM executable
    try:
        slim_path = subprocess.check_output(["realpath", "build/slim"], text=True).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Fallback: try which or use absolute path
        try:
            slim_path = subprocess.check_output(["which", "slim"], text=True).strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            # If all else fails, use the relative path
            slim_path = os.path.abspath("build/slim")

    # Get the absolute path to the SLiM executable
    # slim_path = os.path.abspath("build/slim")

    # Construct SLiM command
    slim_cmd = [
        slim_path,
        "-d", f"L={params.chr_length}",
        "-d", f"N={params.Ne}",
        "-d", f"n_pops={params.n_populations}",
        "-d", f"mig_rate={params.migration_rate}",
        "-d", f"bottleneck={params.bottleneck_size}",
        "-d", f"rec_rate={params.rec_rate}",
        "-d", f"mu={params.mutation_rate}",
        "-d", f"generations={params.generations}",
        "-d", f"s={params.s}",
        "-d", f"h={params.h}",
        "-d", f"selfing={params.selfing_rate}",
        "-d", f"outdir='{scenario_dir}'",
        "-d", f"scenario_id='{params.scenario_id}'",
        "-seed", str(slim_seed),
        slim_script
    ]

    print(f"Command: {' '.join(slim_cmd)}")

    output_path = scenario_dir / "slim_output.out"
    error_path = scenario_dir / "slim_output.error"

    print(f"Running SLiM, output will be saved to: {output_path}")

    with open(output_path, 'w') as out_file, open(error_path, 'w') as err_file:
        result = subprocess.run(
            slim_cmd,
            stdout=out_file,
            stderr=err_file,
            text=True
        )

    print(f"SLiM finished with return code: {result.returncode}")
    if result.returncode != 0:
        print(f"Check error log: {error_path}")

    print("✓ SLiM simulation complete")

    # =========================================================================
    # STEP 2: Post-process with msprime/pyslim
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 2: Post-processing tree sequences")
    print("="*60)

    # Call Python post-processing
    postprocess_trees(
        scenario_dir=scenario_dir,
        params=params
    )

    # =========================================================================
    # STEP 3: Generate outputs (VCF, metadata)
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 3: Generating output files")
    print("="*60)

    generate_outputs(
        scenario_dir=scenario_dir,
        params=params
    )

    # =========================================================================
    # STEP 4: Run downstream analysis
    # =========================================================================
    if params.est > 0:  # Only if mutations present
        print("\n" + "="*60)
        print("STEP 4: Running downstream analysis")
        print("="*60)

        run_downstream_analysis(
            scenario_dir=scenario_dir,
            params=params
        )

    print(f"\n✓ Scenario {params.scenario_id} replicate {params.replicate_id} complete")
    print(f"Results: {scenario_dir}")


def find_slim_script():
    """Find the appropriate SLiM script"""
    possible_paths = [
        "slim_scripts/migration_model.slim",
        "codes/malaria_transmission_study/migration_model_v1.2.slim",
        "../slim_scripts/migration_model.slim",
    ]

    for path in possible_paths:
        if Path(path).exists():
            return path

    raise FileNotFoundError(
        "SLiM script not found. Please create migration_model.slim or specify path"
    )


def postprocess_trees(scenario_dir: Path, params: TransmissionParams):
    """
    Post-process SLiM output with msprime
    - Load tree sequences
    - Parse SLiM stdout for DAF and Ne trajectories
    - Recapitate
    - Simplify to sampled individuals
    - Add neutral mutations
    - Save processed outputs and metadata
    """
    import tskit
    import pyslim
    import msprime
    import pandas as pd
    import re

    print("Loading tree sequences from SLiM...")

    # Find .trees file(s) from SLiM
    tree_files = list(scenario_dir.glob("slim_output_*.trees"))
    if not tree_files:
        raise FileNotFoundError(f"No SLiM output .trees files found in {scenario_dir}")

    # Parse SLiM stdout to extract tracked metrics
    slim_metrics = parse_slim_stdout(scenario_dir)

    for tree_file in tree_files:
        print(f"  Processing: {tree_file.name}")

        # Load
        ts = tskit.load(tree_file)

        # Add debugging before recapitation
        print(f"Number of trees: {ts.num_trees}")
        print(f"Sequence length: {ts.sequence_length}")

        # Now, I replace ts by the recapitated tree ts_recap
        # Sample individuals based on sampling proportion
        n_samples = int(params.outbreak_size * params.sampling_proportion)
        alive_individuals = pyslim.individuals_alive_at(ts, 0)

        if len(alive_individuals) < n_samples:
            print(f"  WARNING: Only {len(alive_individuals)} alive, requested {n_samples}")
            n_samples = len(alive_individuals)

        # Sample across all populations if multiple pops exist
        if params.n_populations > 1:
            # Sample proportionally from each population
            sampled_individuals = sample_across_populations(
                ts, alive_individuals, n_samples, params.n_populations
            )
        else:
            import numpy as np
            np.random.seed(params.replicate_id)
            sampled_individuals = np.random.choice(
                alive_individuals, n_samples, replace=False
            )

        # Get nodes for sampled individuals
        sample_nodes = []
        for ind_id in sampled_individuals:
            ind = ts.individual(int(ind_id))
            sample_nodes.extend(ind.nodes)

        # Simplify (Let's simplify the recapitated tree
        print("  Simplifying tree sequence...")
        sts = ts.simplify(sample_nodes, keep_input_roots=True)

        # Recapitate (I obtained trees with multiple roots so I am trying to recapitate first before simplification
        print("  Recapitating...")
        ts_recap = pyslim.recapitate(
            sts,
            ancestral_Ne=params.Ne,
            recombination_rate=params.rec_rate,
            random_seed=params.replicate_id
        )

        # Simplify again
        ts_recap = ts_recap.simplify(ts_recap.samples())

        # Remove SLiM mutations, add neutral ones
        if ts_recap.num_sites > 0:
            ts_recap = ts_recap.delete_sites(list(range(ts_recap.num_sites)))

        # Save recapitate tree sequence
        output_recap = scenario_dir / f"{tree_file.stem}_recap.trees"
        ts_recap.dump(str(output_recap))

        # Add mutations
        print(f"  Adding mutations (rate={params.mutation_rate})...")
        mts = msprime.sim_mutations(
            ts_recap,
            rate=params.mutation_rate,
            random_seed=params.replicate_id,
            keep=True
        )

        # Save recapitate tree sequence
        output_ts = scenario_dir / f"{tree_file.stem}_mutate.trees"
        mts.dump(str(output_ts))

        print(f"  ✓ Saved: {output_ts.name}")
        print(f"    Sites: {mts.num_sites}, Samples: {mts.num_samples}")

        # Save tracked metrics
        save_slim_metrics(slim_metrics, scenario_dir, tree_file.stem)


def parse_slim_stdout(scenario_dir: Path):
    """
    Parse SLiM stdout to extract:
    - DAF trajectories (global and per-population)
    - True Ne trajectories (per-population and total)
    - Restart count
    - Other metadata
    """
    import re
    import pandas as pd

    # Find SLiM log/output file
    log_files = list(scenario_dir.glob("*.log")) + list(scenario_dir.glob("slim_*.out"))

    metrics = {
        'daf_global': [],
        'daf_by_pop': {},
        'ne_by_pop': {},
        'ne_total': [],
        'restart_count': 0,
        'parameters': {}
    }

    if not log_files:
        print("  WARNING: No SLiM log file found for parsing metrics")
        return metrics

    # Parse the most recent log file
    log_file = log_files[0]

    with open(log_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')

            # Parse parameters
            if len(parts) == 3 and parts[0] == 'Parameter':
                metrics['parameters'][parts[1]] = parts[2]

            # Parse restart counter
            elif len(parts) == 2 and parts[0] == 'restart_counter':
                metrics['restart_count'] = int(parts[1])

            # Parse global DAF
            elif len(parts) == 3 and parts[0] == 'DAF_global':
                gen = int(parts[1])
                freq = float(parts[2])
                metrics['daf_global'].append({'generation': gen, 'frequency': freq})

            # Parse per-population DAF
            elif len(parts) == 3 and parts[0].startswith('DAF_pop'):
                pop_id = parts[0].replace('DAF_pop', '')
                gen = int(parts[1])
                freq = float(parts[2])
                if pop_id not in metrics['daf_by_pop']:
                    metrics['daf_by_pop'][pop_id] = []
                metrics['daf_by_pop'][pop_id].append({'generation': gen, 'frequency': freq})

            # Parse per-population Ne
            elif len(parts) == 3 and parts[0].startswith('True_Ne_pop'):
                pop_id = parts[0].replace('True_Ne_pop', '')
                gen = int(parts[1])
                ne = int(parts[2])
                if pop_id not in metrics['ne_by_pop']:
                    metrics['ne_by_pop'][pop_id] = []
                metrics['ne_by_pop'][pop_id].append({'generation': gen, 'Ne': ne})

            # Parse total Ne
            elif len(parts) == 3 and parts[0] == 'True_Ne_total':
                gen = int(parts[1])
                ne = int(parts[2])
                metrics['ne_total'].append({'generation': gen, 'Ne': ne})

    # Convert lists to DataFrames
    if metrics['daf_global']:
        metrics['daf_global'] = pd.DataFrame(metrics['daf_global'])

    for pop_id in metrics['daf_by_pop']:
        metrics['daf_by_pop'][pop_id] = pd.DataFrame(metrics['daf_by_pop'][pop_id])

    for pop_id in metrics['ne_by_pop']:
        metrics['ne_by_pop'][pop_id] = pd.DataFrame(metrics['ne_by_pop'][pop_id])

    if metrics['ne_total']:
        metrics['ne_total'] = pd.DataFrame(metrics['ne_total'])

    return metrics


def save_slim_metrics(metrics: dict, output_dir: Path, prefix: str):
    """Save parsed SLiM metrics to files"""
    import json
    import pandas as pd

    # Save restart count
    restart_file = output_dir / f"{prefix}_restart_count.txt"
    with open(restart_file, 'w') as f:
        f.write(str(metrics['restart_count']))

    # Save DAF trajectories
    if isinstance(metrics['daf_global'], pd.DataFrame) and not metrics['daf_global'].empty:
        daf_file = output_dir / f"{prefix}_daf_global.tsv"
        metrics['daf_global'].to_csv(daf_file, sep='\t', index=False)

    for pop_id, daf_df in metrics['daf_by_pop'].items():
        if isinstance(daf_df, pd.DataFrame) and not daf_df.empty:
            daf_file = output_dir / f"{prefix}_daf_pop{pop_id}.tsv"
            daf_df.to_csv(daf_file, sep='\t', index=False)

    # Save Ne trajectories
    for pop_id, ne_df in metrics['ne_by_pop'].items():
        if isinstance(ne_df, pd.DataFrame) and not ne_df.empty:
            ne_file = output_dir / f"{prefix}_ne_pop{pop_id}.tsv"
            ne_df.to_csv(ne_file, sep='\t', index=False)

    if isinstance(metrics['ne_total'], pd.DataFrame) and not metrics['ne_total'].empty:
        ne_file = output_dir / f"{prefix}_ne_total.tsv"
        metrics['ne_total'].to_csv(ne_file, sep='\t', index=False)

    # Save all metrics as JSON
    metrics_json = output_dir / f"{prefix}_slim_metrics.json"
    # Convert DataFrames to dicts for JSON serialization
    json_metrics = {
        'restart_count': metrics['restart_count'],
        'parameters': metrics['parameters'],
        'daf_global': metrics['daf_global'].to_dict('records') if isinstance(metrics['daf_global'], pd.DataFrame) else [],
        'daf_by_pop': {k: v.to_dict('records') if isinstance(v, pd.DataFrame) else []
                       for k, v in metrics['daf_by_pop'].items()},
        'ne_by_pop': {k: v.to_dict('records') if isinstance(v, pd.DataFrame) else []
                      for k, v in metrics['ne_by_pop'].items()},
        'ne_total': metrics['ne_total'].to_dict('records') if isinstance(metrics['ne_total'], pd.DataFrame) else []
    }

    with open(metrics_json, 'w') as f:
        json.dump(json_metrics, f, indent=2)

    print(f"  ✓ Saved SLiM metrics: {prefix}_slim_metrics.json")


def sample_across_populations(ts, alive_individuals, n_samples, n_pops):
    """
    Sample individuals proportionally across populations
    """
    import numpy as np

    # Get population assignment for each individual
    pop_assignments = {}
    for ind_id in alive_individuals:
        ind = ts.individual(int(ind_id))
        # Get population from first node
        pop_id = ts.node(ind.nodes[0]).population
        if pop_id not in pop_assignments:
            pop_assignments[pop_id] = []
        pop_assignments[pop_id].append(ind_id)

    # Sample proportionally from each population
    n_per_pop = n_samples // n_pops
    remainder = n_samples % n_pops

    sampled = []
    np.random.seed(42)

    for i, (pop_id, individuals) in enumerate(sorted(pop_assignments.items())):
        n_from_pop = n_per_pop + (1 if i < remainder else 0)
        n_from_pop = min(n_from_pop, len(individuals))

        sampled.extend(
            np.random.choice(individuals, n_from_pop, replace=False)
        )

    return np.array(sampled)

# Write a VCF file with pseudo-homozygous genotypes from a tskit tree sequence
def write_pseudo_homozygous_vcf_v1(ts_mutated, chrno, out_vcf):
    """
    Write a VCF file with pseudo-homozygous genotypes from a tskit tree sequence.

    Args:
        ts_mutated: tskit tree sequence with mutations
        chrno: Chromosome identifier
        out_vcf: Output VCF file path (will be gzipped)
    """
    gt_list = []
    pos_list = []
    ref_list = []
    alt_list = []

    for v in ts_mutated.variants():
        # Only process biallelic sites
        if len(v.alleles) != 2:
            continue

        # Store variant information
        pos_list.append(int(v.position))

        # Use the actual alleles from the variant
        ref_allele = v.alleles[0]
        alt_allele = v.alleles[1]

        # For pseudo-homozygous, we want to represent heterozygous sites as homozygous
        # Convert genotypes: 0->0, 1->1, but if we want all homozygous, we need to decide strategy
        genotypes = v.genotypes

        ref_list.append(ref_allele)
        alt_list.append(alt_allele)
        gt_list.append(genotypes)

    if not pos_list:  # No variants to write
        print(f"Warning: No biallelic variants found for chromosome {chrno}")
        return

    # Create VCF header
    header = f'''##fileformat=VCFv4.2
##source=tskit
##FILTER=<ID=PASS,Description=\"All filters passed\">
##contig=<ID={chrno},length={int(ts_mutated.sequence_length)}>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
'''

    # Create main DataFrame
    df1 = pd.DataFrame({
        '#CHROM': [chrno] * len(pos_list),
        'POS': pos_list,
        'ID': ['.'] * len(pos_list),
        'REF': ref_list,
        'ALT': alt_list,
        'QUAL': ['.'] * len(pos_list),
        'FILTER': ['PASS'] * len(pos_list),
        'INFO': ['.'] * len(pos_list),
        'FORMAT': ['GT'] * len(pos_list)
    })

    # Create genotype DataFrame with pseudo-homozygous encoding
    # Convert 0/1 genotypes to 0|0 and 1|1 respectively
    gt_array = []
    for gt in gt_list:
        # Make all genotypes homozygous by duplicating the allele
        pseudo_homozygous = [f"{g}|{g}" for g in gt]
        gt_array.append(pseudo_homozygous)

    df2 = pd.DataFrame(gt_array)
    df2.columns = [f'tsk_{n}' for n in df2.columns]

    # Combine DataFrames more efficiently
    df = pd.concat([df1, df2], axis=1)

    # Write to gzipped VCF
    with gzip.open(out_vcf, 'wt') as f:
        f.write(header)
        df.to_csv(f, sep='\t', header=True, index=False)

# Alternative Version with Better Performance:
def write_pseudo_homozygous_vcf_v2(ts_mutated, chrno, out_vcf):
    """More efficient version using list comprehensions."""

    # Collect all variant data in one pass
    variants_data = []
    for v in ts_mutated.variants():
        if len(v.alleles) != 2:
            continue

        # Create pseudo-homozygous genotypes
        genotypes = [f"{g}|{g}" for g in v.genotypes]

        variants_data.append({
            'pos': int(v.position),
            'ref': v.alleles[0],
            'alt': v.alleles[1],
            'genotypes': genotypes
        })

    if not variants_data:
        print(f"Warning: No biallelic variants found for chromosome {chrno}")
        return

    # Create header
    header = f'''##fileformat=VCFv4.2
##source=tskit
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID={chrno},length={int(ts_mutated.sequence_length)}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
'''

    # Build the VCF lines directly for better performance
    sample_ids = [f'tsk_{i}' for i in range(len(variants_data[0]['genotypes']))]
    header_line = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_ids

    with gzip.open(out_vcf, 'wt') as f:
        f.write(header)
        f.write('\t'.join(header_line) + '\n')

        for var in variants_data:
            row = [
                chrno,
                str(var['pos']),
                '.',
                var['ref'],
                var['alt'],
                '.',
                'PASS',
                '.',
                'GT'
            ] + var['genotypes']
            f.write('\t'.join(row) + '\n')


def generate_outputs(scenario_dir: Path, params: TransmissionParams):
    """Generate VCF and metadata files"""
    import tskit
    import gzip

    # Find processed tree files
    tree_files = list(scenario_dir.glob("*_processed.trees"))

    for tree_file in tree_files:
        print(f"  Processing: {tree_file.name}")

        ts = tskit.load(tree_file)

        # Generate VCF
        vcf_file = scenario_dir / f"{tree_file.stem}.vcf.gz"

        #with gzip.open(vcf_file, 'wt') as f:
        #    ts.write_vcf(f, contig_id="chr1")

        write_pseudo_homozygous_vcf_v1(ts, "Pf3D7_01_v3", vcf_file)

        vcf2 = scenario_dir / f"{tree_file.stem}_v2.vcf.gz"
        write_pseudo_homozygous_vcf_v2(ts, "Pf3D7_01_v3", vcf2)

        print(f"  ✓ VCF: {vcf_file.name}")

        # Generate metadata
        metadata = {
            'scenario_id': params.scenario_id,
            'replicate_id': params.replicate_id,
            'n_samples': ts.num_samples,
            'n_sites': ts.num_sites,
            'sequence_length': ts.sequence_length,
            'n_trees': ts.num_trees,
            'parameters': asdict(params)
        }

        import json
        meta_file = scenario_dir / f"{tree_file.stem}_metadata.json"
        with open(meta_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        print(f"  ✓ Metadata: {meta_file.name}")


def run_downstream_analysis(scenario_dir: Path, params: TransmissionParams):
    """
    Run downstream analysis pipelines:
    1. Phylogenetic tree
    2. Population structure
    3. IBD analysis
    4. Summary statistics
    """
    import subprocess

    vcf_files = list(scenario_dir.glob("*.vcf.gz"))

    for vcf_file in vcf_files:
        print(f"\n  Analyzing: {vcf_file.name}")

        # 1. Build phylogenetic tree
        print("    - Building phylogenetic tree...")
        run_phylogenetic_analysis(vcf_file, scenario_dir)

        # 2. Calculate summary statistics
        print("    - Computing summary statistics...")
        compute_summary_stats(vcf_file, scenario_dir)

#        # 3. Population structure (if multiple populations)
#        if params.n_populations > 1:
#            print("    - Analyzing population structure...")
#            analyze_population_structure(vcf_file, scenario_dir, params)

        # 3. Population genetics analysis (NEW)
        if params.n_populations > 1:
            print("    - Running population genetics analysis...")
            compute_population_genetics(
                vcf_file,
                scenario_dir,
                n_pops=params.n_populations
            )

        # 4. Population structure
        print("    - Analyzing population structure...")
        analyze_population_structure(vcf_file, scenario_dir, params)

def run_phylogenetic_analysis(vcf_file: Path, output_dir: Path):
    """Build phylogenetic tree from VCF"""
    # This would call your existing phylo pipeline
    # Or use a simplified version
    pass  # Implement based on your phylo script


def compute_summary_stats_v1(vcf_file: Path, output_dir: Path):
    """Compute genetic diversity statistics with proper error handling"""
    import allel
    import numpy as np

    try:
        # Load VCF with explicit parameter to handle multi-allelic sites
        callset = allel.read_vcf(str(vcf_file), fields=['calldata/GT', 'variants/POS'])

        if callset is None or 'calldata/GT' not in callset:
            print(f"  WARNING: No genotype data found in {vcf_file}")
            return

        # Convert to GenotypeArray - this handles the dimension issue
        gt = allel.GenotypeArray(callset['calldata/GT'])

        # Count number of variants and samples
        n_variants = gt.n_variants
        n_samples = gt.n_samples

        if n_variants == 0:
            print(f"  WARNING: No variants found in {vcf_file}")
            stats = {
                'n_variants': 0,
                'n_samples': n_samples,
                'mean_het': 0.0,
                'pi': 0.0,
                'segregating_sites': 0
            }
        else:
            # Calculate statistics with proper dimension handling
            het_matrix = gt.is_het()
            mean_het = float(np.mean(het_matrix))

            # For pi (nucleotide diversity), we need allele counts
            # Convert genotypes to allele counts
            ac = gt.count_alleles()

            # Calculate pi - this handles the dimensions correctly
            pi = float(allel.sequence_diversity(callset['variants/POS'], ac))

            # Alternative: mean pairwise difference if you prefer
            # This requires counting alleles per variant first
            # ac = gt.count_alleles()
            # pi = float(allel.mean_pairwise_difference(ac))

            stats = {
                'n_variants': n_variants,
                'n_samples': n_samples,
                'mean_het': mean_het,
                'pi': pi,
                'segregating_sites': int(np.sum(ac[:, 1] > 0))  # Sites with minor alleles
            }

        # Save statistics
        import json
        stats_file = output_dir / f"{vcf_file.stem}_stats_v1.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)

        print(f"  ✓ Computed stats: {n_variants} variants, π={stats['pi']:.4f}")

    except Exception as e:
        print(f"  ERROR computing stats for {vcf_file}: {str(e)}")
        # Save error information
        error_stats = {
            'error': str(e),
            'n_variants': 0,
            'n_samples': 0,
            'mean_het': 0.0,
            'pi': 0.0,
            'segregating_sites': 0
        }
        stats_file = output_dir / f"{vcf_file.stem}_stats_v1.json"
        with open(stats_file, 'w') as f:
            json.dump(error_stats, f, indent=2)

def compute_summary_stats(vcf_file: Path, output_dir: Path):
    """Compute comprehensive genetic diversity statistics"""
    import allel
    import numpy as np

    try:
        # Read VCF with specific fields to avoid memory issues
        callset = allel.read_vcf(
            str(vcf_file),
            fields=['calldata/GT', 'variants/POS', 'variants/CHROM'],
            alt_number=2  # Limit to biallelic sites
        )

        if callset is None:
            print(f"  WARNING: Could not read VCF {vcf_file}")
            return None

        gt_array = allel.GenotypeArray(callset['calldata/GT'])

        if gt_array.n_variants == 0:
            print(f"  WARNING: No variants in {vcf_file}")
            return {
                'n_variants': 0, 'n_samples': gt_array.n_samples,
                'mean_het': 0.0, 'pi': 0.0, 'tajima_d': 0.0,
                'segregating_sites': 0, 'private_alleles': 0
            }

        # Basic counts
        n_variants = gt_array.n_variants
        n_samples = gt_array.n_samples

        # Heterozygosity
        het_matrix = gt_array.is_het()
        mean_het = float(np.mean(het_matrix))

        # Allele counts for diversity statistics
        ac = gt_array.count_alleles()

        # Nucleotide diversity (pi)
        if 'variants/POS' in callset:
            pos = callset['variants/POS']
            pi = float(allel.sequence_diversity(pos, ac))
        else:
            pi = float(allel.mean_pairwise_difference(ac))

        # Segregating sites
        seg_sites = int(np.sum(ac[:, 1] > 0))  # Sites with at least one minor allele

        # Tajima's D (if enough variants)
        tajima_d = 0.0
        if seg_sites > 1:
            try:
                tajima_d = float(allel.tajima_d(ac, pos=callset['variants/POS'] if 'variants/POS' in callset else None))
            except:
                tajima_d = 0.0  # Fallback if calculation fails

        stats = {
            'n_variants': n_variants,
            'n_samples': n_samples,
            'mean_het': mean_het,
            'pi': pi,
            'tajima_d': tajima_d,
            'segregating_sites': seg_sites,
            'private_alleles': 0  # Could calculate if you have population info
        }

        # Save to file
        import json
        stats_file = output_dir / f"{vcf_file.stem}_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)

        print(f"  ✓ Stats: {n_variants} variants, π={pi:.4f}, D={tajima_d:.3f}")
        return stats

    except Exception as e:
        print(f"  ERROR in compute_summary_stats for {vcf_file}: {str(e)}")
        # Return minimal stats to avoid breaking the pipeline
        return {
            'n_variants': 0, 'n_samples': 0,
            'mean_het': 0.0, 'pi': 0.0, 'tajima_d': 0.0,
            'segregating_sites': 0, 'private_alleles': 0,
            'error': str(e)
        }


def compute_summary_stats_v2(vcf_file: Path, output_dir: Path):
    """Minimal robust version that won't break your pipeline"""
    import allel
    import numpy as np

    try:
        # Simple approach - just count variants and samples
        callset = allel.read_vcf(str(vcf_file), fields=['calldata/GT'])

        if callset is None:
            return {'n_variants': 0, 'n_samples': 0, 'pi': 0.0, 'mean_het': 0.0}

        gt = allel.GenotypeArray(callset['calldata/GT'])

        stats = {
            'n_variants': gt.n_variants,
            'n_samples': gt.n_samples,
            'pi': 0.0,  # Skip complex calculations
            'mean_het': 0.0
        }

        # Save basic stats
        import json
        stats_file = output_dir / f"{vcf_file.stem}_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)

        return stats

    except Exception as e:
        print(f"  Minimal stats failed for {vcf_file}: {e}")
        return {'n_variants': 0, 'n_samples': 0, 'pi': 0.0, 'mean_het': 0.0, 'error': str(e)}

def analyze_population_structure(vcf_file: Path, output_dir: Path, params: TransmissionParams):
    """Analyze population structure with PCA"""
    # This would run PCA, FST, etc.
    pass  # Implement population genetics analyses


def compute_population_genetics(vcf_file: Path, output_dir: Path,
                              population_assignments: Dict[str, List[int]] = None,
                              n_pops: int = None):
    """
    Compute comprehensive population genetic statistics including:
    - FST (global and pairwise)
    - PCA
    - Genetic diversity indices
    - Population structure metrics

    Args:
        vcf_file: Path to VCF file
        output_dir: Output directory for results
        population_assignments: Dict mapping pop names to sample indices
        n_pops: Number of populations (if assignments not provided)
    """
    import allel
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    import json
    import seaborn as sns

    results = {}

    try:
        print(f"  Computing population genetics for {vcf_file.name}")

        # Load VCF data
        callset = allel.read_vcf(
            str(vcf_file),
            fields=['calldata/GT', 'variants/POS', 'variants/CHROM', 'samples'],
            alt_number=2  # Biallelic sites only
        )

        if callset is None or 'calldata/GT' not in callset:
            print(f"  WARNING: No genotype data in {vcf_file}")
            return None

        # Create GenotypeArray
        gt = allel.GenotypeArray(callset['calldata/GT'])
        samples = callset['samples']
        n_samples = len(samples)

        print(f"    Loaded {gt.n_variants} variants for {n_samples} samples")

        # If population assignments not provided, create dummy assignments
        if population_assignments is None:
            if n_pops is None:
                n_pops = 3  # Default
            population_assignments = create_dummy_populations(samples, n_pops)

        # =========================================================================
        # 1. BASIC DIVERSITY STATISTICS PER POPULATION
        # =========================================================================
        print("    Computing diversity statistics...")
        diversity_stats = compute_diversity_per_population(gt, population_assignments)
        results['diversity'] = diversity_stats

        # =========================================================================
        # 2. FST ANALYSIS - Global and Pairwise
        # =========================================================================
        print("    Computing FST statistics...")
        fst_results = compute_fst_analysis(gt, population_assignments)
        results['fst'] = fst_results

        # =========================================================================
        # 3. PRINCIPAL COMPONENT ANALYSIS (PCA)
        # =========================================================================
        print("    Running PCA...")
        pca_results = compute_pca_analysis(gt, population_assignments, samples)
        results['pca'] = pca_results

        # =========================================================================
        # 4. POPULATION STRUCTURE (ADMIXTURE-like analysis)
        # =========================================================================
        print("    Computing population structure...")
        structure_results = compute_population_structure(gt, population_assignments, samples)
        results['structure'] = structure_results

        # =========================================================================
        # 5. GENETIC DISTANCE MATRIX
        # =========================================================================
        print("    Computing genetic distances...")
        distance_results = compute_genetic_distances(gt, population_assignments, samples)
        results['distances'] = distance_results

        # =========================================================================
        # SAVE ALL RESULTS
        # =========================================================================
        save_popgen_results(results, output_dir, vcf_file.stem)

        # Generate summary plots
        generate_popgen_plots(results, output_dir, vcf_file.stem)

        print(f"  ✓ Population genetics complete: {vcf_file.stem}")
        return results

    except Exception as e:
        print(f"  ERROR in population genetics for {vcf_file}: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


def create_dummy_populations(samples: np.ndarray, n_pops: int) -> Dict[str, List[int]]:
    """Create dummy population assignments if none provided"""
    pop_assignments = {}
    samples_per_pop = len(samples) // n_pops

    for i in range(n_pops):
        start_idx = i * samples_per_pop
        if i == n_pops - 1:  # Last population gets remaining samples
            end_idx = len(samples)
        else:
            end_idx = (i + 1) * samples_per_pop

        pop_name = f"pop_{i+1}"
        pop_assignments[pop_name] = list(range(start_idx, end_idx))

    return pop_assignments


def compute_diversity_per_population(gt, population_assignments: Dict) -> Dict:
    """Compute genetic diversity statistics for each population"""
    import allel
    import numpy as np

    diversity_stats = {}

    for pop_name, sample_indices in population_assignments.items():
        if len(sample_indices) == 0:
            continue

        # Extract genotypes for this population
        pop_gt = gt.take(sample_indices, axis=1)

        # Count alleles
        ac = pop_gt.count_alleles()

        # Filter for segregating sites
        seg_sites = ac[:, 1] > 0  # Sites with minor alleles

        if np.sum(seg_sites) == 0:
            # No variation in this population
            diversity_stats[pop_name] = {
                'n_samples': len(sample_indices),
                'n_variants': pop_gt.n_variants,
                'n_segregating': 0,
                'pi': 0.0,
                'tajima_d': 0.0,
                'waterson_theta': 0.0,
                'mean_heterozygosity': 0.0
            }
            continue

        # Segregating sites count
        n_seg = np.sum(seg_sites)

        # Nucleotide diversity (pi)
        pi = float(allel.sequence_diversity(np.arange(pop_gt.n_variants)[seg_sites],
                                          ac[seg_sites]))

        # Waterson's theta
        theta_w = float(allel.watterson_theta(np.arange(pop_gt.n_variants)[seg_sites],
                                            ac[seg_sites]))

        # Tajima's D
        try:
            tajima_d = float(allel.tajima_d(ac[seg_sites]))
        except:
            tajima_d = 0.0

        # Mean heterozygosity
        het_matrix = pop_gt[seg_sites].is_het()
        mean_het = float(np.mean(het_matrix)) if het_matrix.size > 0 else 0.0

        diversity_stats[pop_name] = {
            'n_samples': len(sample_indices),
            'n_variants': pop_gt.n_variants,
            'n_segregating': int(n_seg),
            'pi': pi,
            'tajima_d': tajima_d,
            'waterson_theta': theta_w,
            'mean_heterozygosity': mean_het
        }

    return diversity_stats


def compute_fst_analysis(gt, population_assignments: Dict) -> Dict:
    """Compute FST statistics - global and pairwise"""
    import allel
    import numpy as np
    from itertools import combinations

    fst_results = {}

    # Get population names and sample indices
    pop_names = list(population_assignments.keys())

    if len(pop_names) < 2:
        return {'global': 0.0, 'pairwise': {}, 'message': 'Need at least 2 populations for FST'}

    # Prepare genotype arrays for each population
    pop_genotypes = []
    for pop_name in pop_names:
        sample_indices = population_assignments[pop_name]
        if len(sample_indices) > 0:
            pop_gt = gt.take(sample_indices, axis=1)
            pop_genotypes.append(pop_gt)
        else:
            pop_genotypes.append(None)

    # Global FST (Weir & Cockerham 1984)
    try:
        # Convert to allele counts for FST calculation
        ac_list = []
        for pop_gt in pop_genotypes:
            if pop_gt is not None:
                ac = pop_gt.count_alleles()
                ac_list.append(ac)

        if len(ac_list) >= 2:
            global_fst, _, _ = allel.weir_cockerham_fst(ac_list)
            fst_results['global'] = float(np.nanmean(global_fst))
        else:
            fst_results['global'] = 0.0
    except Exception as e:
        print(f"      Global FST calculation failed: {e}")
        fst_results['global'] = 0.0

    # Pairwise FST
    fst_results['pairwise'] = {}
    for (i, j) in combinations(range(len(pop_names)), 2):
        pop1_name = pop_names[i]
        pop2_name = pop_names[j]

        if (pop_genotypes[i] is not None and pop_genotypes[j] is not None and
            len(population_assignments[pop1_name]) > 0 and
            len(population_assignments[pop2_name]) > 0):

            try:
                # Extract genotypes for this pair
                pop1_gt = pop_genotypes[i]
                pop2_gt = pop_genotypes[j]

                # Count alleles
                ac1 = pop1_gt.count_alleles()
                ac2 = pop2_gt.count_alleles()

                # Compute pairwise FST
                fst, _, _ = allel.weir_cockerham_fst([ac1, ac2])
                pairwise_fst = float(np.nanmean(fst))

                fst_results['pairwise'][f"{pop1_name}_{pop2_name}"] = pairwise_fst

            except Exception as e:
                print(f"      Pairwise FST {pop1_name}-{pop2_name} failed: {e}")
                fst_results['pairwise'][f"{pop1_name}_{pop2_name}"] = 0.0

    return fst_results


def compute_pca_analysis(gt, population_assignments: Dict, samples: np.ndarray) -> Dict:
    """Perform Principal Component Analysis"""
    import allel
    import numpy as np
    from sklearn.decomposition import PCA

    pca_results = {}

    try:
        # Convert to allele counts matrix for PCA
        # Use mean imputation for missing data
        gn = gt.to_n_alt()

        # Replace missing values with mean of the variant
        for i in range(gn.shape[0]):
            variant = gn[i]
            mask = variant == -1
            if np.any(mask):
                mean_val = np.mean(variant[~mask])
                gn[i, mask] = mean_val

        # Transpose for PCA (samples × variants)
        X = gn.T

        # Standardize the data
        X_std = (X - np.mean(X, axis=0)) / np.std(X, axis=0)

        # Perform PCA
        pca = PCA(n_components=min(10, X_std.shape[0], X_std.shape[1]))
        X_pca = pca.fit_transform(X_std)

        # Prepare results
        pca_results['explained_variance_ratio'] = pca.explained_variance_ratio_.tolist()
        pca_results['components'] = pca.components_.tolist()
        pca_results['n_components'] = pca.n_components_

        # Create sample information with population assignments
        sample_info = []
        for i, sample in enumerate(samples):
            # Find which population this sample belongs to
            pop_name = "unknown"
            for pop, indices in population_assignments.items():
                if i in indices:
                    pop_name = pop
                    break

            sample_info.append({
                'sample': sample,
                'population': pop_name,
                'pc1': float(X_pca[i, 0]) if X_pca.shape[1] > 0 else 0,
                'pc2': float(X_pca[i, 1]) if X_pca.shape[1] > 1 else 0,
                'pc3': float(X_pca[i, 2]) if X_pca.shape[1] > 2 else 0
            })

        pca_results['samples'] = sample_info
        pca_results['transformed'] = X_pca.tolist()

    except Exception as e:
        print(f"      PCA failed: {e}")
        pca_results['error'] = str(e)

    return pca_results


def compute_population_structure(gt, population_assignments: Dict, samples: np.ndarray) -> Dict:
    """Compute population structure metrics"""
    import allel
    import numpy as np

    structure_results = {}

    try:
        # Calculate allele frequency differentiation
        pop_names = list(population_assignments.keys())
        allele_freqs = {}

        for pop_name in pop_names:
            sample_indices = population_assignments[pop_name]
            if len(sample_indices) > 0:
                pop_gt = gt.take(sample_indices, axis=1)
                ac = pop_gt.count_alleles()
                # Calculate allele frequencies (avoid division by zero)
                total_alleles = ac.sum(axis=1)
                freq = np.divide(ac[:, 1], total_alleles,
                               out=np.zeros_like(ac[:, 1], dtype=float),
                               where=total_alleles!=0)
                allele_freqs[pop_name] = freq

        # Calculate population-specific F statistics
        structure_results['allele_frequency_differentiation'] = allele_freqs

        # Mean allele frequency differences between populations
        if len(pop_names) >= 2:
            freq_differences = {}
            for i in range(len(pop_names)):
                for j in range(i+1, len(pop_names)):
                    pop1, pop2 = pop_names[i], pop_names[j]
                    if pop1 in allele_freqs and pop2 in allele_freqs:
                        diff = np.mean(np.abs(allele_freqs[pop1] - allele_freqs[pop2]))
                        freq_differences[f"{pop1}_{pop2}"] = float(diff)

            structure_results['mean_allele_frequency_differences'] = freq_differences

    except Exception as e:
        print(f"      Population structure analysis failed: {e}")
        structure_results['error'] = str(e)

    return structure_results


def compute_genetic_distances(gt, population_assignments: Dict, samples: np.ndarray) -> Dict:
    """Compute genetic distance matrices"""
    import allel
    import numpy as np
    from scipy.spatial.distance import pdist, squareform

    distance_results = {}

    try:
        # Convert to allele sharing matrix
        gn = gt.to_n_alt()

        # Handle missing data
        for i in range(gn.shape[0]):
            variant = gn[i]
            mask = variant == -1
            if np.any(mask):
                mean_val = np.mean(variant[~mask])
                gn[i, mask] = mean_val

        # Calculate pairwise genetic distances (Euclidean)
        X = gn.T  # Samples × variants
        pairwise_distances = pdist(X, metric='euclidean')
        distance_matrix = squareform(pairwise_distances)

        distance_results['pairwise_distance_matrix'] = distance_matrix.tolist()
        distance_results['mean_distance'] = float(np.mean(pairwise_distances))
        distance_results['max_distance'] = float(np.max(pairwise_distances))

        # Population-level distances
        pop_names = list(population_assignments.keys())
        pop_distances = {}

        for pop1 in pop_names:
            for pop2 in pop_names:
                indices1 = population_assignments[pop1]
                indices2 = population_assignments[pop2]

                if indices1 and indices2:
                    # Extract submatrix for these populations
                    submatrix = distance_matrix[np.ix_(indices1, indices2)]
                    pop_distances[f"{pop1}_{pop2}"] = {
                        'mean': float(np.mean(submatrix)),
                        'std': float(np.std(submatrix)),
                        'min': float(np.min(submatrix)),
                        'max': float(np.max(submatrix))
                    }

        distance_results['population_distances'] = pop_distances

    except Exception as e:
        print(f"      Genetic distance calculation failed: {e}")
        distance_results['error'] = str(e)

    return distance_results


def save_popgen_results(results: Dict, output_dir: Path, prefix: str):
    """Save population genetics results to files"""

    # Save main results as JSON
    results_file = output_dir / f"{prefix}_popgen_results.json"

    # Convert numpy types to Python types for JSON serialization
    def convert_types(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_types(item) for item in obj]
        else:
            return obj

    results_serializable = convert_types(results)

    with open(results_file, 'w') as f:
        json.dump(results_serializable, f, indent=2, default=str)

    # Save diversity statistics as CSV
    if 'diversity' in results:
        diversity_df = pd.DataFrame(results['diversity']).T
        diversity_file = output_dir / f"{prefix}_diversity_stats.csv"
        diversity_df.to_csv(diversity_file)

    # Save FST matrix as CSV
    if 'fst' in results and 'pairwise' in results['fst']:
        fst_df = pd.DataFrame([results['fst']['pairwise']])
        fst_file = output_dir / f"{prefix}_fst_matrix.csv"
        fst_df.to_csv(fst_file)

    print(f"    Saved results: {results_file}")


def generate_popgen_plots(results: Dict, output_dir: Path, prefix: str):
    """Generate visualization plots for population genetics results"""
    import matplotlib.pyplot as plt
    import seaborn as sns

    try:
        plt.style.use('seaborn-v0_8')
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # 1. PCA Plot
        if 'pca' in results and 'samples' in results['pca']:
            samples_info = results['pca']['samples']
            if samples_info:
                pc1 = [s['pc1'] for s in samples_info]
                pc2 = [s['pc2'] for s in samples_info]
                populations = [s['population'] for s in samples_info]

                unique_pops = list(set(populations))
                colors = plt.cm.Set3(np.linspace(0, 1, len(unique_pops)))

                for i, pop in enumerate(unique_pops):
                    mask = [p == pop for p in populations]
                    axes[0,0].scatter(np.array(pc1)[mask], np.array(pc2)[mask],
                                    label=pop, alpha=0.7, color=colors[i])

                axes[0,0].set_xlabel('PC1')
                axes[0,0].set_ylabel('PC2')
                axes[0,0].set_title('PCA Plot')
                axes[0,0].legend()

                # Add variance explained
                if 'explained_variance_ratio' in results['pca']:
                    var_exp = results['pca']['explained_variance_ratio']
                    axes[0,0].text(0.02, 0.98, f"Var: PC1={var_exp[0]:.3f}, PC2={var_exp[1]:.3f}",
                                 transform=axes[0,0].transAxes, verticalalignment='top')

        # 2. FST Heatmap
        if 'fst' in results and 'pairwise' in results['fst']:
            pairwise_fst = results['fst']['pairwise']
            if pairwise_fst:
                # Create matrix for heatmap
                pop_pairs = list(pairwise_fst.keys())
                pops = sorted(set([p.split('_')[0] for p in pop_pairs] +
                                 [p.split('_')[1] for p in pop_pairs]))

                fst_matrix = np.zeros((len(pops), len(pops)))
                for i, pop1 in enumerate(pops):
                    for j, pop2 in enumerate(pops):
                        if i == j:
                            fst_matrix[i,j] = 0
                        else:
                            key1 = f"{pop1}_{pop2}"
                            key2 = f"{pop2}_{pop1}"
                            if key1 in pairwise_fst:
                                fst_matrix[i,j] = pairwise_fst[key1]
                            elif key2 in pairwise_fst:
                                fst_matrix[i,j] = pairwise_fst[key2]

                sns.heatmap(fst_matrix, annot=True, xticklabels=pops, yticklabels=pops,
                           ax=axes[0,1], cmap='YlOrRd')
                axes[0,1].set_title('Pairwise FST Matrix')

        # 3. Diversity Statistics
        if 'diversity' in results:
            diversity_data = []
            for pop, stats in results['diversity'].items():
                diversity_data.append({
                    'Population': pop,
                    'Pi': stats.get('pi', 0),
                    'Tajima_D': stats.get('tajima_d', 0),
                    'Theta_W': stats.get('waterson_theta', 0)
                })

            if diversity_data:
                div_df = pd.DataFrame(diversity_data)
                x = np.arange(len(div_df))
                width = 0.25

                axes[1,0].bar(x - width, div_df['Pi'], width, label='π')
                axes[1,0].bar(x, div_df['Tajima_D'], width, label="Tajima's D")
                axes[1,0].bar(x + width, div_df['Theta_W'], width, label="θW")

                axes[1,0].set_xlabel('Population')
                axes[1,0].set_ylabel('Value')
                axes[1,0].set_title('Diversity Statistics')
                axes[1,0].set_xticks(x)
                axes[1,0].set_xticklabels(div_df['Population'])
                axes[1,0].legend()

        # 4. Distance Distribution
        if 'distances' in results and 'population_distances' in results['distances']:
            pop_dists = results['distances']['population_distances']
            within_pop = [v['mean'] for k, v in pop_dists.items() if k.split('_')[0] == k.split('_')[1]]
            between_pop = [v['mean'] for k, v in pop_dists.items() if k.split('_')[0] != k.split('_')[1]]

            if within_pop and between_pop:
                axes[1,1].boxplot([within_pop, between_pop],
                                labels=['Within Pop', 'Between Pop'])
                axes[1,1].set_ylabel('Genetic Distance')
                axes[1,1].set_title('Within vs Between Population Distances')

        plt.tight_layout()
        plot_file = output_dir / f"{prefix}_popgen_plots.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"    Generated plots: {plot_file}")

    except Exception as e:
        print(f"    Plot generation failed: {e}")


def main():
    """Main entry point with flexible command-line interface"""

    parser = argparse.ArgumentParser(
        description="Multi-population malaria transmission simulation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate full design (75 core × 4 sampling × 3 outbreak × 10 reps = 9000 runs)
  python malaria_migration_simulation.py --generate-design full --n-replicates 10

  # Generate migration-focused design (smaller, focused on population structure)
  python malaria_migration_simulation.py --generate-design migration --n-replicates 20

  # Generate core biological scenarios only
  python malaria_migration_simulation.py --generate-design core --n-replicates 10

  # Run a single scenario from design file
  python malaria_migration_simulation.py --run-from-design simulation_design.json --task-id 0

  # Test a single parameter combination
  python malaria_migration_simulation.py --test-run \\
    --rec-rate 1e-7 --bottleneck 5 --est 0.5 --outbreak-size 200
        """
    )

    # Design generation options
    parser.add_argument("--generate-design", type=str, choices=['full', 'migration', 'core'],
                       help="Generate simulation design")
    parser.add_argument("--n-replicates", type=int, default=10,
                       help="Number of replicates per scenario")
    parser.add_argument("--outdir", type=str, default="migration_sims",
                       help="Output directory")

    # Run from existing design
    parser.add_argument("--run-from-design", type=str,
                       help="JSON file with simulation design")
    parser.add_argument("--task-id", type=int,
                       help="Task ID for array job (0-indexed)")

    # Test single run
    parser.add_argument("--test-run", action="store_true",
                       help="Run a single test scenario")
    parser.add_argument("--rec-rate", type=float, default=6.67e-7)
    parser.add_argument("--bottleneck", type=int, default=5)
    parser.add_argument("--est", type=float, default=0.5)
    parser.add_argument("--outbreak-size", type=int, default=200)
    parser.add_argument("--sampling-prop", type=float, default=0.1)
    parser.add_argument("--n-populations", type=int, default=3)
    parser.add_argument("--migration-rate", type=float, default=0.01)

    # General options
    parser.add_argument("--dry-run", action="store_true",
                       help="Print what would be done without executing")

    args = parser.parse_args()

    # Generate design
    if args.generate_design:
        design = SimulationDesign(
            design_type=args.generate_design,
            n_replicates=args.n_replicates,
            outdir=args.outdir
        )
        design.print_design_summary()
        design.save_design()
        design.generate_slurm_array()

        print("\nNext steps:")
        print("1. Review the design in simulation_design.json")
        print(f"2. Submit array job: sbatch {args.outdir}/run_migration_sim.slurm")
        print("3. Monitor progress: squeue -u $USER")
        return

    # Run from existing design (for SLURM array)
    if args.run_from_design:
        if args.task_id is None:
            print("Error: --task-id required when running from design")
            sys.exit(1)

        with open(args.run_from_design, 'r') as f:
            design_data = json.load(f)

        if args.task_id >= len(design_data['parameters']):
            print(f"Error: task_id {args.task_id} out of range")
            sys.exit(1)

        param_dict = design_data['parameters'][args.task_id]
        params = TransmissionParams(**param_dict)
        run_single_scenario(params, dry_run=args.dry_run)
        return

    # Test single run
    if args.test_run:
        mutation_rate = ParameterSpace.get_est_from_mutation_rate(args.est)
        params = TransmissionParams(
            rec_rate=args.rec_rate,
            bottleneck_size=args.bottleneck,
            est=args.est,
            mutation_rate=mutation_rate,
            outbreak_size=args.outbreak_size,
            sampling_proportion=args.sampling_prop,
            n_populations=args.n_populations,
            migration_rate=args.migration_rate,
            scenario_id="TEST",
            replicate_id=1,
            outdir=args.outdir
        )
        run_single_scenario(params, dry_run=args.dry_run)
        return

    # No action specified
    parser.print_help()


if __name__ == "__main__":
    main()
