# Quick Start Guide: Multi-Population Malaria Simulations

## Installation & Setup

### 1. Prerequisites
```bash
# Python 3.8+
conda create -n malaria-sim python=3.10
conda activate malaria-sim

# Install required packages
conda install -c conda-forge numpy pandas matplotlib seaborn
pip install tskit pyslim msprime

# Install SLiM
# On Ubuntu/Debian:
sudo apt-get install slim

# Or from source: https://messerlab.org/slim/
```

### 2. Download scripts
```bash
# Create project directory
mkdir malaria_transmission_study
cd malaria_transmission_study

# Copy the three main scripts:
# - malaria_migration_simulation.py
# - analyze_simulation_design.py  
# - Your existing PfSimulation class

chmod +x *.py
```

---

## Workflow: Start to Finish

### Step 1: Generate Simulation Design

Choose one of three designs based on your research goals:

#### Option A: Core Biological Scenarios (Recommended for initial exploration)
```bash
# 75 core scenarios × 10 replicates = 750 simulations
python malaria_migration_simulation.py \
  --generate-design core \
  --n-replicates 10 \
  --outdir core_sims

python3 codes/malaria_transmission_study/malaria_migration_full_pipeline_v1.3.py  \
  --generate-design full \
  --outdir sim_migration \
  --n-replicates 10

# Estimated time: 2-3 days with 50 parallel jobs
```

#### Option B: Migration-Focused Design
```bash
# Population structure emphasis
# ~6,480 simulations with 20 replicates
python malaria_migration_simulation.py \
  --generate-design migration \
  --n-replicates 20 \
  --outdir migration_sims

python3 codes/malaria_transmission_study/malaria_migration_full_pipeline_v1.3.py  \
  --generate-design migration \
  --outdir sim_migration \
  --n-replicates 10

# Estimated time: 5-7 days with 100 parallel jobs
```

#### Option C: Full Factorial Design (Most comprehensive)
```bash
# Complete parameter space
# 9,000 simulations with 10 replicates
python malaria_migration_simulation.py \
  --generate-design full \
  --n-replicates 10 \
  --outdir full_design_sims

# Estimated time: 1-2 weeks with 100 parallel jobs
```

### Step 2: Analyze the Design

```bash
# Generate summary statistics and visualizations
python analyze_simulation_design.py \
  core_sims/simulation_design.json \
  --output-dir core_sims/design_analysis

python codes/malaria_transmission_study/analyze_simulation_design.py sim_migration/simulation_design.json

# This creates:
#   - Parameter space heatmaps
#   - Batch processing summaries
#   - Scenario lists for reference
```

### Step 3: Test a Single Scenario

Before submitting thousands of jobs, test your setup:

```bash
# Dry run (print commands only)
python malaria_migration_simulation.py --test-run \
  --rec-rate 1e-7 \
  --bottleneck 5 \
  --est 0.5 \
  --outbreak-size 200 \
  --sampling-prop 0.3 \
  --n-populations 3 \
  --migration-rate 0.01 \
  --dry-run

# Actual test run
python malaria_migration_simulation.py --test-run \
  --rec-rate 1e-7 \
  --bottleneck 5 \
  --est 0.5 \
  --outbreak-size 200 \
  --sampling-prop 0.3 \
  --outdir test_output
```

### Step 4: Submit to HPC Cluster

```bash
cd core_sims

# Submit SLURM array job
sbatch run_migration_sim.sh

sbatch sim_migration/run_migration_sim.slurm

# Monitor progress
squeue -u $USER
watch -n 30 'squeue -u $USER | head -20'

# Check log files for errors
tail -f logs/migration_*.out
tail -f logs/migration_*.err
```

### Step 5: Run Specific Scenarios Manually

```bash
# Run scenario 0
python malaria_migration_simulation.py \
  --run-from-design core_sims/simulation_design.json \
  --task-id 0

# Run scenarios 0-9 (useful for debugging)
for i in {0..9}; do
  python malaria_migration_simulation.py \
    --run-from-design core_sims/simulation_design.json \
    --task-id $i
done
```

---

## Parameter Selection Guide

### When to use each design:

| Design Type | Use When | Scenarios | Time |
|-------------|----------|-----------|------|
| **Core** | Initial exploration, pilot study | 750 | 2-3 days |
| **Migration** | Focus on population structure | 6,480 | 5-7 days |
| **Full** | Comprehensive study, publication | 9,000 | 1-2 weeks |

### Adjusting parameters:

```bash
# Fewer replicates (faster, less precise)
--n-replicates 5

# More replicates (slower, more robust)
--n-replicates 20

# Custom output directory
--outdir /scratch/username/sims_run1
```

---

## Understanding Output

### Directory structure:
```
core_sims/
├── simulation_design.json       # Full parameter specifications
├── run_migration_sim.sh        # SLURM submission script
├── design_analysis/            # Design summaries and plots
│   ├── parameter_space_triad.png
│   ├── scenario_list.csv
│   └── batch_summary.csv
├── logs/                       # SLURM logs
│   ├── migration_JOBID_0.out
│   └── migration_JOBID_0.err
└── results/                    # Simulation outputs (created during runs)
    ├── S001_very_low_tight_none/
    │   ├── rep_001/
    │   │   ├── trees/
    │   │   ├── vcf/
    │   │   └── metadata.json
    │   └── rep_002/
    └── S002_very_low_tight_very_low/
```

### Key files:

1. **simulation_design.json**: Complete parameter specifications
2. **scenario_list.csv**: Unique scenarios (excluding replicates)
3. **batch_summary.csv**: Batch processing information
4. **Parameter space plots**: Visualizations of design coverage

---

## Troubleshooting

### Problem: "Too many scenarios, will take forever!"

**Solution**: Start with core design (750 scenarios), then expand:

```bash
# Start small
python malaria_migration_simulation.py \
  --generate-design core \
  --n-replicates 5 \
  --outdir pilot_sims
```

### Problem: SLURM jobs failing

**Check**:
1. Log files: `tail logs/migration_*.err`
2. Test single scenario: Use `--test-run`
3. Memory: Increase `#SBATCH --mem` in SLURM script
4. Time limit: Increase `#SBATCH --time`

### Problem: Need to run specific subset of scenarios

**Solution**: Create custom scenario list:

```python
# In Python
import json
with open('simulation_design.json', 'r') as f:
    design = json.load(f)

# Select scenarios with high recombination only
high_rec_scenarios = [p for p in design['parameters']
                     if p['rec_rate'] >= 1e-6]

# Save subset
subset_design = design.copy()
subset_design['parameters'] = high_rec_scenarios
with open('subset_design.json', 'w') as f:
    json.dump(subset_design, f, indent=2)
```

### Problem: Integration with existing code

**TODO**: Modify `run_single_scenario()` in `malaria_migration_simulation.py`:

```python
def run_single_scenario(params: TransmissionParams, dry_run=False):
    # TODO: Replace this with your PfSimulation class

    # Convert parameters
    pf_params = convert_to_pf_params(params)

    # Run simulation
    sim = YourPfSimulation(pf_params)
    sim.run_with_migration(...)
```

---

## Resource Requirements

### Per scenario estimates:

| Component | CPU Time | Memory | Disk |
|-----------|----------|--------|------|
| Forward sim (SLiM) | 10-30 min | 4-8 GB | 100 MB |
| Recapitation | 2-5 min | 2-4 GB | 50 MB |
| VCF generation | 5-10 min | 2-4 GB | 500 MB |
| **Total** | **20-45 min** | **8-16 GB** | **650 MB** |

### Cluster allocation:

```bash
# For core design (750 scenarios):
# 750 × 30 min / 60 min / 50 parallel = ~7.5 hours
# Request: 50 nodes × 16 GB × 12 hours

# For full design (9000 scenarios):
# 9000 × 30 min / 60 min / 100 parallel = ~45 hours
# Request: 100 nodes × 16 GB × 48 hours
```

---

## Next Steps After Simulation

### 1. Aggregate results
```bash
# Collect all VCF files
find results/ -name "*.vcf.gz" > vcf_list.txt

# Collect metadata
find results/ -name "metadata.json" -exec cat {} \; > all_metadata.json
```

### 2. Transmission inference
```bash
# Run your inference methods on simulated data
# Examples: TransPhylo, outbreaker2, SCOTTI, etc.
```

### 3. Evaluate performance
```python
# Calculate metrics for each scenario
# - Precision: fraction of inferred links that are correct
# - Recall: fraction of true links that were inferred
# - F1 score: harmonic mean of precision and recall
```

### 4. Meta-analysis
```python
# Identify parameter combinations where methods work well
# Generate decision rules for practitioners
# Create diagnostic plots
```

---

## Quick Reference: Command Cheatsheet

```bash
# Generate core design
python malaria_migration_simulation.py --generate-design core --n-replicates 10

# Analyze design
python analyze_simulation_design.py simulation_design.json

# Test single scenario
python malaria_migration_simulation.py --test-run --dry-run

# Submit jobs
sbatch run_migration_sim.sh

# Monitor
squeue -u $USER
sacct -j JOBID

# Run specific task
python malaria_migration_simulation.py --run-from-design simulation_design.json --task-id 0

# Check logs
tail -f logs/migration_*.err

N=$(($(wc -l < sim_migration/phylo_results/logs/missing_trees.csv) - 1))
sbatch --array=1-${N} sim_migration/codes/04_run_phylogenetic_migration_failed.slurm
squeue -u $USER
```

---

## Getting Help

1. **Check documentation**: `SIMULATION_DESIGN.md`
2. **Test mode**: Use `--test-run --dry-run`
3. **Start small**: Begin with 5-10 replicates
4. **Validate**: Run pilot before full-scale submission
5. **Monitor**: Watch first few jobs for errors

## Contact

For questions about:
- **Simulation design**: Review `SIMULATION_DESIGN.md`
- **Parameter choices**: See biological justification in design doc
- **Technical issues**: Check SLURM logs and test single scenarios
- **Integration**: Modify `run_single_scenario()` function
