# Multi-Population Malaria Transmission Simulation Design

## Overview

This simulation framework explores how **within-host diversity**, **transmission bottlenecks**, **mutation accumulation**, **sampling strategy**, and **population structure** affect our ability to infer transmission chains in malaria outbreaks.

---

## Parameter Space Structure

### A. Core Biological Triad (Primary factors)

These three parameters form the foundation of transmission biology:

#### 1. **Within-Host Diversity** (Recombination Rate)
Proxy for complexity of infection (COI) and genetic diversity within hosts.

| Level | Rate | Interpretation |
|-------|------|----------------|
| Very Low | 1×10⁻⁹ | All clonal, no recombination |
| Low | 1×10⁻⁸ | Almost clonal, rare recombination |
| Medium | 1×10⁻⁷ | Moderate recombination |
| High | 1×10⁻⁶ | Frequent recombination, high COI |
| Very High | 1×10⁻⁵ | Very high diversity, multiple clones |

**Biological relevance**: In low-transmission settings, infections are typically clonal. In high-transmission areas, mixed-clone infections are common.

#### 2. **Transmission Bottleneck Size** (B)
Number of parasites successfully establishing infection in new host.

| Level | Size | Interpretation |
|-------|------|----------------|
| Tight | 1 | Single parasite (severe bottleneck) |
| Medium | 5 | Few parasites transmitted |
| Loose | 20 | Many parasites (relaxed bottleneck) |

**Biological relevance**: Mosquito-to-human transmission typically involves 1-10 sporozoites. Tighter bottlenecks preserve genetic signatures of transmission.

#### 3. **Expected Substitutions per Transmission** (EST)
Mutations accumulating per transmission event.

| Level | EST | Interpretation |
|-------|-----|----------------|
| None | 0.0 | No mutations (pure genealogy) |
| Very Low | 0.1 | Rare mutations |
| Low-Moderate | 0.5 | Some mutation load |
| Moderate | 1.0 | Moderate mutation signal |
| High | 2.0 | High mutation accumulation |

**Biological relevance**: 
- P. falciparum genome: ~23 Mb
- Mutation rate: ~10⁻⁹ to 10⁻⁸ per site per generation
- EST = μ × genome_length × generations_per_transmission
- Typical: 0.1-1.0 substitutions per transmission

**Core scenarios**: 5 × 3 × 5 = **75 biological scenarios**

---

### B. Sampling & Scale Factors

These parameters determine observational constraints and outbreak characteristics:

#### 4. **Sampling Proportion** (ρ)
Fraction of infections that are sequenced.

| Level | Proportion |
|-------|------------|
| Low | 0.1 (10%) |
| Medium-Low | 0.3 (30%) |
| Medium-High | 0.6 (60%) |
| High | 0.9 (90%) |

**Impact**: Lower sampling creates more missing links, affecting reconstruction accuracy.

#### 5. **Outbreak Size** (n)
Total number of infections in the outbreak.

| Level | Size | Context |
|-------|------|---------|
| Small | 50 | Village-level outbreak |
| Medium | 200 | District-level outbreak |
| Large | 800 | Regional outbreak |

**Impact**: Larger outbreaks increase complexity and computational burden.

---

### C. Population Structure (Migration scenarios)

#### 6. **Number of Populations**
Number of interconnected demes/villages.

- **Levels**: 3, 4, 5 populations

#### 7. **Migration Rate** (m)
Rate of parasite movement between populations.

| Level | Rate | Interpretation |
|-------|------|----------------|
| Low | 0.001 | Rare migration, mostly isolated |
| Medium | 0.01 | Moderate gene flow |
| High | 0.05 | High connectivity |

**Biological relevance**: Human movement patterns, mosquito flight range, imported cases.

---

## Simulation Designs

### Design 1: Full Factorial (Most comprehensive)

**Formula**: 75 core × 4 sampling × 3 outbreak sizes × N replicates

**Total scenarios**: 75 × 4 × 3 × 10 = **9,000 simulations** (with 10 replicates)

**Use case**: Complete exploration of parameter space

**Computational cost**: High (~36-72 hours with parallelization)

```bash
# Generate design
python malaria_migration_simulation.py \
  --generate-design full \
  --n-replicates 10 \
  --outdir full_design_sims
```

---

### Design 2: Migration-Focused (Population structure emphasis)

**Formula**: 3 representative core scenarios × 3 n_populations × 3 migration rates × 4 sampling × 3 outbreak sizes × N replicates

**Total scenarios**: 3 × 3 × 3 × 4 × 3 × 20 = **6,480 simulations** (with 20 replicates)

**Representative scenarios**:
1. Low recombination + Tight bottleneck + Low-moderate EST
2. Medium recombination + Medium bottleneck + Moderate EST  
3. High recombination + Loose bottleneck + High EST

**Use case**: Focus on how population structure affects transmission inference

```bash
# Generate design
python malaria_migration_simulation.py \
  --generate-design migration \
  --n-replicates 20 \
  --outdir migration_sims
```

---

### Design 3: Core Biological Only (Minimal)

**Formula**: 75 core scenarios × N replicates

**Total scenarios**: 75 × 10 = **750 simulations**

**Use case**: 
- Pilot study
- Focus on biological factors only
- Fixed sampling (30%) and outbreak size (200)

```bash
# Generate design
python malaria_migration_simulation.py \
  --generate-design core \
  --n-replicates 10 \
  --outdir core_sims
```

---

## Replication Strategy

### Why 10-20 replicates?

1. **Statistical power**: 
   - Compute mean ± SD for precision/recall metrics
   - Generate bootstrap-based 95% confidence intervals
   - Detect ~10% differences in performance

2. **Stochastic variation**:
   - Forward simulation is inherently stochastic
   - Genetic drift varies between runs
   - Transmission trees differ randomly

3. **Practical balance**:
   - 10 replicates: Minimum for stable estimates
   - 20 replicates: Better precision, especially for rare events
   - >20 replicates: Diminishing returns for computational cost

---

## Execution Workflow

### Step 1: Generate Design

```bash
# Choose your design
python malaria_migration_simulation.py \
  --generate-design [full|migration|core] \
  --n-replicates 10 \
  --outdir output_directory
```

**Outputs**:
- `simulation_design.json`: Full parameter specifications
- `run_migration_sim.sh`: SLURM array job script

### Step 2: Review Design

```bash
# Inspect the design
head -100 output_directory/simulation_design.json

# Check number of scenarios
jq '.n_scenarios' output_directory/simulation_design.json
```

### Step 3: Submit to HPC

```bash
# Submit SLURM array job
cd output_directory
sbatch run_migration_sim.sh

# Monitor progress
squeue -u $USER
watch -n 30 'squeue -u $USER'
```

### Step 4: Test Single Scenario

```bash
# Before full submission, test one scenario
python malaria_migration_simulation.py --test-run \
  --rec-rate 1e-7 \
  --bottleneck 5 \
  --est 0.5 \
  --outbreak-size 200 \
  --sampling-prop 0.3 \
  --dry-run
```

---

## Integration with Existing Code

### Connecting to your PfSimulation class

In the `run_single_scenario()` function, you'll integrate with your existing simulation framework:

```python
def run_single_scenario(params: TransmissionParams, dry_run: bool = False):
    """Execute simulation with migration"""
    
    if dry_run:
        print("[DRY RUN] Would execute simulation")
        return
    
    # Convert TransmissionParams to your PfSimParams format
    pf_params = PfSimParams(
        Ne=params.Ne,
        n_samples=int(params.outbreak_size * params.sampling_proportion),
        generations=params.generations,
        n_chromosomes=14,  # P. falciparum has 14 chromosomes
        selfing_rate=params.selfing_rate,
        num_origins=1,
        h=params.h,
        s=params.s,
        g_sel_start=params.g_sel_start,
        sim_relatedness=0,
        g_ne_change_start=0,
        N0=params.Ne,
        mutation_rate=params.mutation_rate,
        outdir=params.outdir
    )
    
    # Initialize and run simulation
    sim = PfSimulation(pf_params)
    
    # Run with migration between populations
    sim.run_migration_simulation(
        rec_rate=params.rec_rate,
        n_populations=params.n_populations,
        migration_rate=params.migration_rate,
        bottleneck_size=params.bottleneck_size,
        scenario_id=params.scenario_id,
        replicate_id=params.replicate_id
    )
```

---

## Output Organization

Recommended directory structure:

```
migration_sims/
├── simulation_design.json          # Full parameter specifications
├── run_migration_sim.sh           # SLURM submission script
├── logs/                          # SLURM output logs
│   ├── migration_12345_0.out
│   ├── migration_12345_1.out
│   └── ...
├── results/                       # Simulation outputs
│   ├── S001_very_low_tight_none/
│   │   ├── rep_001/
│   │   │   ├── vcf/
│   │   │   ├── trees/
│   │   │   ├── metadata.json
│   │   │   └── summary_stats.tsv
│   │   ├── rep_002/
│   │   └── ...
│   ├── S002_very_low_tight_very_low/
│   └── ...
└── analysis/                      # Downstream analysis
    ├── transmission_inference/
    ├── method_evaluation/
    └── figures/
```

---

## Key Research Questions Addressed

### 1. **When do transmission inference methods work?**
- Low recombination + tight bottleneck + high sampling → Good performance
- High recombination + loose bottleneck + low sampling → Poor performance

### 2. **How does population structure affect inference?**
- Can we distinguish local transmission from imports?
- How much migration confounds transmission chains?

### 3. **What's the minimum data quality needed?**
- Sampling threshold for reliable inference
- Outbreak size where complexity overwhelms methods

### 4. **Which biological factors matter most?**
- Relative importance of recombination vs. bottleneck vs. mutation
- Interaction effects between parameters

### 5. **Method-specific performance**
- Which inference methods excel in which scenarios?
- Are there universally robust methods?

---

## Computational Considerations

### Resource estimates (per scenario):

| Component | Time | Memory | Notes |
|-----------|------|--------|-------|
| SLiM forward sim | 5-30 min | 2-8 GB | Depends on Ne, generations |
| msprime recapitation | 1-5 min | 1-4 GB | Depends on sample size |
| Mutation overlay | 1-2 min | 1-2 GB | Fast with tskit |
| VCF generation | 2-10 min | 2-4 GB | Depends on variant density |

**Total per scenario**: ~10-50 minutes, 8-16 GB RAM

**Full design (9000 scenarios)**: 
- Serial: ~90-450 days
- Parallel (100 jobs): ~1-5 days

**Recommended**:
- SLURM array with `--array=0-8999%100` (max 100 simultaneous jobs)
- 16-32 GB memory per job
- 4-8 hour time limit per job

---

## Next Steps

1. **Test the framework**:
   ```bash
   python malaria_migration_simulation.py --test-run --dry-run
   ```

2. **Generate pilot design** (core only):
   ```bash
   python malaria_migration_simulation.py --generate-design core --n-replicates 5
   ```

3. **Run pilot** (small subset):
   ```bash
   # Manually run first 10 scenarios
   for i in {0..9}; do
     python malaria_migration_simulation.py --run-from-design simulation_design.json --task-id $i
   done
   ```

4. **Validate outputs**, then scale to full design

5. **Implement transmission inference methods** on simulated data

6. **Evaluate performance metrics**: precision, recall, F1, network accuracy
