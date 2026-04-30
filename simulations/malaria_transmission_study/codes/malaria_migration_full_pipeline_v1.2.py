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
    rec_rate: float  # Recombination rate (proxy for within-host diversity)
    bottleneck_size: int  # Number of parasites transmitted
    est: float  # Expected substitutions per transmission

    # Population structure
    n_populations: int = 3  # Number of demes/populations
    migration_rate: float = 0.01  # Migration rate between populations

    # Population sizes and sampling
    Ne: int = 10000  # Effective population size per deme
    outbreak_size: int = 200  # Total number of infections
    sampling_proportion: float = 0.3  # Proportion of infections sampled

    # Evolutionary parameters
    generations: int = 200  # Forward simulation time
    mutation_rate: float = 1e-8  # Per-site per-generation

    # Selection parameters (optional - can simulate neutral or selected scenarios)
    s: float = 0.0  # Selection coefficient (0 = neutral)
    h: float = 0.5  # Dominance coefficient
    g_sel_start: int = 0  # Generation when selection starts (0 = no selection)

    # Selfing/inbreeding
    selfing_rate: float = 0.0  # Selfing rate (0 = random mating, 0.95 = high inbreeding)

    # Simulation metadata
    replicate_id: int = 1
    scenario_id: str = ""
    outdir: str = "migration_sims"


class ParameterSpace:
    """
    Defines the full parameter space for the transmission simulation study
    """

    # A. Core biological triad
    RECOMBINATION_RATES = {
        'very_low': 1e-9,      # All clonal
        'low': 1e-8,           # Almost clonal
        'medium': 1e-7,        # Moderate recombination
        'high': 1e-6,          # Frequent recombination
        'very_high': 1e-5      # Very high recombination / high COI
    }

    BOTTLENECK_SIZES = {
        'tight': 1,      # Single parasite (severe bottleneck)
        'medium': 5,     # Few parasites
        'loose': 20      # Many parasites (relaxed bottleneck)
    }

    EXPECTED_SUBS_PER_TRANSMISSION = {
        'none': 0.0,         # No mutations
        'very_low': 0.1,     # Very rare mutations
        'low_mod': 0.5,      # Low-moderate
        'moderate': 1.0,     # Moderate mutation load
        'high': 2.0          # High mutation accumulation
    }

    # B. Sampling and scale
    SAMPLING_PROPORTIONS = [0.1, 0.3, 0.6, 0.9]
    OUTBREAK_SIZES = [50, 200, 800]

    # C. Population structure
    N_POPULATIONS = [3, 4, 5]  # Number of interconnected populations
    MIGRATION_RATES = {
        'low': 0.001,      # Rare migration
        'medium': 0.01,    # Moderate gene flow
        'high': 0.05       # High connectivity
    }

    @classmethod
    def get_est_from_mutation_rate(cls, est: float, genome_length: int = 23e6) -> float:
        """
        Convert expected substitutions per transmission to mutation rate
        EST = mutation_rate × genome_length × generations_per_transmission
        Assume ~10-15 days per transmission cycle in mosquito + human
        """
        # Malaria generation time approximation:
        # ~1 day in mosquito, parasite replication every 48h in human
        # Estimate ~1 transmission = ~10 parasite generations
        generations_per_transmission = 10

        if est == 0:
            return 0

        mutation_rate = est / (genome_length * generations_per_transmission)
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
        Total: 75 core × 4 sampling × 3 outbreak sizes × n_replicates
        """
        core_scenarios = cls.get_core_scenarios()
        all_params = []

        for scenario in core_scenarios:
            for sampling_prop in cls.SAMPLING_PROPORTIONS:
                for outbreak_size in cls.OUTBREAK_SIZES:
                    for replicate in range(1, n_replicates + 1):

                        # Calculate mutation rate from EST
                        mutation_rate = cls.get_est_from_mutation_rate(scenario['est'])

                        params = TransmissionParams(
                            rec_rate=scenario['rec_rate'],
                            bottleneck_size=scenario['bottleneck_size'],
                            est=scenario['est'],
                            sampling_proportion=sampling_prop,
                            outbreak_size=outbreak_size,
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
            ('low', 'tight', 'low_mod'),      # Low rec, tight bottleneck, moderate mutation
            ('medium', 'medium', 'moderate'),  # Balanced scenario
            ('high', 'loose', 'high'),         # High diversity scenario
        ]

        all_params = []
        scenario_id = 0

        for rec_name, bottle_name, est_name in representative_scenarios:
            rec_rate = cls.RECOMBINATION_RATES[rec_name]
            bottle_size = cls.BOTTLENECK_SIZES[bottle_name]
            est_value = cls.EXPECTED_SUBS_PER_TRANSMISSION[est_name]
            mutation_rate = cls.get_est_from_mutation_rate(est_value)

            for n_pops in cls.N_POPULATIONS:
                for mig_name, mig_rate in cls.MIGRATION_RATES.items():
                    for sampling_prop in cls.SAMPLING_PROPORTIONS:
                        for outbreak_size in cls.OUTBREAK_SIZES:
                            for replicate in range(1, n_replicates + 1):
                                scenario_id += 1

                                params = TransmissionParams(
                                    rec_rate=rec_rate,
                                    bottleneck_size=bottle_size,
                                    est=est_value,
                                    n_populations=n_pops,
                                    migration_rate=mig_rate,
                                    sampling_proportion=sampling_prop,
                                    outbreak_size=outbreak_size,
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
                 outdir: str = "migration_sims"):
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
                            array_size: int = 100):
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
#SBATCH --array=0-{n_scenarios-1}%{array_size}
#SBATCH --partition=standard

# Load modules and activate environment
#source ~/.bashrc
#conda activate slim-msprime

# Run Python simulation for this array task
python3 malaria_migration_simulation.py \\
    --run-from-design simulation_design.json \\
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
    scenario_dir = Path(params.outdir) / params.scenario_id / f"rep_{params.replicate_id:03d}"
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

    # Construct SLiM command
    slim_cmd = [
        "slim",
        "-d", f"L={params.Ne}",
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

    # Run SLiM
    result = subprocess.run(
        slim_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        print(f"ERROR: SLiM failed with code {result.returncode}")
        print(f"STDERR: {result.stderr}")
        raise RuntimeError("SLiM simulation failed")

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
        "codes/migration_model.slim",
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

        # Simplify
        print("  Simplifying tree sequence...")
        sts = ts.simplify(sample_nodes, keep_input_roots=True)

        # Recapitate
        print("  Recapitating...")
        rts = pyslim.recapitate(
            sts,
            ancestral_Ne=params.Ne,
            recombination_rate=params.rec_rate,
            random_seed=params.replicate_id
        )

        # Simplify again to just sampled nodes
        rts = rts.simplify(rts.samples())

        # Remove SLiM mutations, add neutral ones
        if rts.num_sites > 0:
            print(f"  Removing {rts.num_sites} SLiM mutations...")
            rts = rts.delete_sites(list(range(rts.num_sites)))

        # Add mutations
        print(f"  Adding mutations (rate={params.mutation_rate})...")
        mts = msprime.sim_mutations(
            rts,
            rate=params.mutation_rate,
            random_seed=params.replicate_id,
            keep=True
        )

        # Save processed tree sequence
        output_ts = scenario_dir / f"{tree_file.stem}_processed.trees"
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
    log_files = list(scenario_dir.glob("*.log")) + list(scenario_dir.glob("slim_*.txt"))

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

        with gzip.open(vcf_file, 'wt') as f:
            ts.write_vcf(f, contig_id="chr1")

        write_peudo_homozygous_vcf(ts, "Pf3D7_01_v3", vcf_file)

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

        # 3. Population structure (if multiple populations)
        if params.n_populations > 1:
            print("    - Analyzing population structure...")
            analyze_population_structure(vcf_file, scenario_dir, params)


def run_phylogenetic_analysis(vcf_file: Path, output_dir: Path):
    """Build phylogenetic tree from VCF"""
    # This would call your existing phylo pipeline
    # Or use a simplified version
    pass  # Implement based on your phylo script


def compute_summary_stats(vcf_file: Path, output_dir: Path):
    """Compute genetic diversity statistics"""
    import allel
    import numpy as np

    # Load VCF
    callset = allel.read_vcf(str(vcf_file))

    if callset is None:
        return

    gt = allel.GenotypeArray(callset['calldata/GT'])

    # Calculate statistics
    stats = {
        'n_variants': gt.n_variants,
        'n_samples': gt.n_samples,
        'mean_het': float(np.mean(gt.is_het())),
        'pi': float(allel.mean_pairwise_difference(gt)),
    }

    # Save
    import json
    stats_file = output_dir / f"{vcf_file.stem}_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)


def analyze_population_structure(vcf_file: Path, output_dir: Path, params: TransmissionParams):
    """Analyze population structure with PCA"""
    # This would run PCA, FST, etc.
    pass  # Implement population genetics analyses


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
    parser.add_argument("--rec-rate", type=float, default=1e-7)
    parser.add_argument("--bottleneck", type=int, default=5)
    parser.add_argument("--est", type=float, default=0.5)
    parser.add_argument("--outbreak-size", type=int, default=200)
    parser.add_argument("--sampling-prop", type=float, default=0.3)
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
