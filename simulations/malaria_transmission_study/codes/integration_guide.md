# Complete Malaria Migration Simulation Framework

## Overview

This framework connects **three major components**:

```
Design Generator → SLiM Simulation → Post-Processing → Analysis
     (Python)          (SLiM)         (tskit/msprime)    (Various)
```

---

## Component 1: Parameter Design Generator

**File**: `malaria_migration_simulation.py`

**Purpose**: Generate all parameter combinations for systematic exploration

**What it does**:
1. Defines parameter space (recombination, bottleneck, mutation, sampling, etc.)
2. Creates factorial design (e.g., 9,000 scenarios)
3. Saves to `simulation_design.json`
4. Generates SLURM array job script

**Usage**:
```bash
# Generate design
python malaria_migration_simulation.py \
  --generate-design full \
  --n-replicates 10 \
  --outdir migration_sims
```

**Output**: `simulation_design.json` containing all parameter combinations

---

## Component 2: SLiM Forward Simulation

**File**: `migration_model.slim`

**Purpose**: Run forward-time evolutionary simulation with migration

**What it does**:
1. Creates multiple populations (demes)
2. Simulates migration between populations
3. Implements transmission bottlenecks
4. Tracks genealogies (tree sequences)
5. Outputs `.trees` file

**Key parameters** (passed from Python):
- `N`: Population size per deme
- `n_pops`: Number of populations
- `mig_rate`: Migration rate between populations
- `bottleneck`: Transmission bottleneck size
- `rec_rate`: Recombination rate
- `generations`: Simulation time

**Called by**: `run_single_scenario()` function in Python

---

## Component 3: Post-Processing Pipeline

**What it does**:
1. **Load SLiM output** (.trees files)
2. **Sample individuals** based on sampling proportion
3. **Recapitate** - add ancient ancestry with msprime
4. **Simplify** - reduce to sampled lineages
5. **Add mutations** - overlay neutral mutations
6. **Generate VCF** - output genotype data
7. **Create metadata** - save simulation parameters

**Tools used**:
- `tskit`: Tree sequence manipulation
- `pyslim`: SLiM-specific tree sequence operations
- `msprime`: Coalescent simulation for recapitation

---

## Component 4: Downstream Analysis

**What it does**:
1. **Phylogenetic trees** - IQ-TREE for ML phylogeny
2. **Population genetics** - PCA, FST, admixture
3. **IBD/IBS analysis** - Identity by descent
4. **Summary statistics** - diversity, heterozygosity
5. **Transmission inference** - Apply inference methods

---

## Complete Workflow

### Step 1: Setup

```bash
# Create directory structure
mkdir -p migration_study/{scripts,slim_scripts,results,logs}
cd migration_study

# Copy files
cp malaria_migration_simulation.py scripts/
cp migration_model.slim slim_scripts/
cp analyze_simulation_design.py scripts/
```

### Step 2: Generate Design

```bash
python scripts/malaria_migration_simulation.py \
  --generate-design core \
  --n-replicates 10 \
  --outdir results
```

**Output**:
- `results/simulation_design.json` - All parameters
- `results/run_migration_sim.sh` - SLURM script

### Step 3: Inspect Design

```bash
# Visualize parameter space
python scripts/analyze_simulation_design.py \
  results/simulation_design.json \
  --output-dir results/design_analysis

# Check number of scenarios
jq '.n_scenarios' results/simulation_design.json
```

### Step 4: Test Single Scenario

```bash
# Dry run
python scripts/malaria_migration_simulation.py --test-run --dry-run

# Actual test
python scripts/malaria_migration_simulation.py --test-run \
  --rec-rate 1e-7 \
  --bottleneck 5 \
  --est 0.5 \
  --outdir results/test
```

### Step 5: Submit Full Job Array

```bash
cd results
sbatch run_migration_sim.sh

# Monitor
watch -n 30 'squeue -u $USER | head -20'
```

### Step 6: Check Progress

```bash
# Count completed
find results -name "*_processed.trees" | wc -l

# Check for errors
grep -r "ERROR" logs/

# View specific log
tail -f logs/migration_12345_42.out
```

### Step 7: Analyze Results

```bash
# Generate summary
python scripts/summarize_results.py results/

# Run downstream analysis on completed scenarios
python scripts/run_downstream_analysis.py results/
```

---

## File Organization

```
migration_study/
├── scripts/
│   ├── malaria_migration_simulation.py    # Main framework
│   ├── analyze_simulation_design.py        # Design analysis
│   └── summarize_results.py                # Results aggregation
├── slim_scripts/
│   └── migration_model.slim                # SLiM simulation
├── results/
│   ├── simulation_design.json              # Parameter design
│   ├── run_migration_sim.sh                # SLURM script
│   ├── S001_very_low_tight_none/
│   │   ├── rep_001/
│   │   │   ├── parameters.json
│   │   │   ├── slim_output_*.trees
│   │   │   ├── *_processed.trees
│   │   │   ├── *.vcf.gz
│   │   │   ├── *_metadata.json
│   │   │   ├── *_stats.json
│   │   │   └── *.treefile
│   │   ├── rep_002/
│   │   └── ...
│   ├── S002_very_low_tight_very_low/
│   └── ...
└── logs/
    ├── migration_12345_0.out
    └── migration_12345_0.err
```

---

## Key Functions Explained

### 1. `run_single_scenario()`

**Purpose**: Execute complete simulation pipeline for one parameter set

**Steps**:
```python
1. Create output directory
2. Save parameters to JSON
3. Run SLiM simulation
4. Post-process tree sequences
5. Generate VCF files
6. Run downstream analysis
7. Save results
```

### 2. `postprocess_trees()`

**Purpose**: Convert SLiM output to analysis-ready format

**Steps**:
```python
1. Load .trees from SLiM
2. Sample individuals (based on sampling_proportion)
3. Simplify tree sequence
4. Recapitate (add ancient history)
5. Remove SLiM mutations
6. Add neutral mutations
7. Save processed tree sequence
```

### 3. `generate_outputs()`

**Purpose**: Create analysis files (VCF, metadata)

**Steps**:
```python
1. Load processed tree sequences
2. Write VCF with genotypes
3. Extract metadata (n_samples, n_sites, etc.)
4. Save to structured format
```

### 4. `run_downstream_analysis()`

**Purpose**: Apply analytical methods

**Steps**:
```python
1. Build phylogenetic trees
2. Calculate summary statistics
3. Analyze population structure (if migration)
4. Compute IBD/IBS
5. Run transmission inference
```

---

## Integration with Your Existing Code

### Option 1: Minimal Integration

Replace `run_single_scenario()` to call your existing R script:

```python
def run_single_scenario(params: TransmissionParams, dry_run=False):
    import subprocess
    
    # Convert params to R script arguments
    r_cmd = [
        "Rscript", "your_simulation_script.R",
        "--genome_set_id", str(params.replicate_id),
        "--N", str(params.Ne),
        "--rec_rate", str(params.rec_rate),
        "--n_populations", str(params.n_populations),
        "--migration_rate", str(params.migration_rate),
        # ... more parameters
    ]
    
    subprocess.run(r_cmd, check=True)
```

### Option 2: Python Wrapper

Keep your R code but add Python wrapper:

```python
# your_simulation_wrapper.py
import rpy2.robjects as ro

def run_simulation(params):
    # Load R script
    ro.r.source("your_simulation.R")
    
    # Call R function
    result = ro.r['run_simulation'](
        N=params.Ne,
        rec_rate=params.rec_rate,
        # ... more parameters
    )
    
    return result
```

### Option 3: Full Python Rewrite

Rewrite simulation logic in Python using tskit/msprime:

```python
def run_slim_simulation(params):
    import subprocess
    
    cmd = ["slim", "-d", f"N={params.Ne}", ...]
    subprocess.run(cmd, check=True)
    
def postprocess_with_msprime(trees_file, params):
    import tskit, pyslim, msprime
    
    ts = tskit.load(trees_file)
    # ... processing steps
    return mts
```

---

## Debugging Tips

### Check SLiM is Running

```bash
# Test SLiM manually
slim -d "N=1000" -d "n_pops=3" migration_model.slim

# Check output
ls -lh *.trees
```

### Verify Tree Sequences

```python
import tskit

ts = tskit.load("slim_output_test.trees")
print(f"Samples: {ts.num_samples}")
print(f"Trees: {ts.num_trees}")
print(f"Sites: {ts.num_sites}")
```

### Test Post-Processing

```python
from malaria_migration_simulation import postprocess_trees, TransmissionParams

params = TransmissionParams(...)
postprocess_trees(scenario_dir, params)
```

### Check VCF Generation

```bash
# Count variants
zcat output.vcf.gz | grep -v "^#" | wc -l

# Check samples
bcftools query -l output.vcf.gz
```

---

## Common Issues

### Issue 1: "No .trees files found"

**Cause**: SLiM didn't complete or output path wrong

**Fix**:
```python
# Check SLiM output path in migration_model.slim
outfile = outdir + "/slim_output_" + scenario_id + ".trees";
```

### Issue 2: "Not enough alive individuals"

**Cause**: Population crashed during simulation

**Fix**: Increase population size or reduce bottleneck severity

### Issue 3: "Module not found: pyslim"

**Cause**: Missing Python dependencies

**Fix**:
```bash
pip install tskit pyslim msprime scikit-allel
```

### Issue 4: VCF has no variants

**Cause**: Mutation rate too low or sites filtered out

**Fix**: Check `params.mutation_rate` and EST calculation

---

## Next Steps

1. **Test framework end-to-end** with small design (5 scenarios)
2. **Validate outputs** - check VCFs, trees, metadata
3. **Run pilot study** - 100 scenarios to check performance
4. **Scale up** - Run full design (9,000 scenarios)
5. **Analyze results** - Aggregate across replicates
6. **Method evaluation** - Test transmission inference methods

---

## Performance Estimates

| Component | Time (per scenario) | Memory | Notes |
|-----------|---------------------|--------|-------|
| SLiM | 5-30 min | 2-8 GB | Depends on N, generations |
| Recapitation | 1-5 min | 1-4 GB | Depends on samples |
| Mutations | 1-2 min | 1-2 GB | Fast with tskit |
| VCF generation | 2-5 min | 2-4 GB | Depends on variants |
| Phylogeny | 30-120 min | 8-16 GB | Depends on bootstrap |
| **Total** | **40-160 min** | **16-32 GB** | Per scenario |

**Full design (9,000 scenarios)**:
- Serial: ~1,000-2,400 hours (42-100 days)
- Parallel (100 jobs): ~10-24 hours

---

## Summary

The framework provides:

✅ **Systematic parameter exploration** - 9,000 biologically relevant scenarios  
✅ **Forward simulation** - SLiM with migration and bottlenecks  
✅ **Efficient processing** - tskit for tree sequence operations  
✅ **Standardized outputs** - VCF, trees, metadata  
✅ **Parallel execution** - SLURM array jobs  
✅ **Downstream analysis** - Phylogenetics, pop gen, inference  

**The missing piece you identified was the connection between design and execution** - now implemented in `run_single_scenario()` with full pipeline integration!
