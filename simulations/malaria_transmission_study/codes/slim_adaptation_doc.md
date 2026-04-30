# SLiM Migration Model: Adaptation from Single Population

## Overview

The new `migration_model.slim` adapts your existing `single_pop.slim` script to support **multiple populations with migration** while preserving all critical features.

---

## Key Features Preserved

### ✅ 1. **Parameter System with `set_default()`**
```slim
function (NULL) set_default(string k, lifs v) {
	if(!exists(k)) defineConstant(k, v);
	catn(c("Parameter", k, v), sep='\t');
	return NULL;
}
```
- Same flexible parameter system
- All parameters can be set from command line
- Defaults provided if not specified

### ✅ 2. **Relatedness-Based Mating (`modifyChild`)**
```slim
s0:10 modifyChild() {
	if (!sim_relatedness) return T;
	
	if (child.index % 2 == 0) return T;  // 50% pass without check
	
	prob = parent1.relatedness(parent2);
	if ((prob > 1.0 / 32) & (prob < 1.0 / 4)) return T;
	
	else return sample(c(T, F), size=1, replace=T, weights=c(1, 999));
}
```
- **Preserved exactly** for simulating inbreeding
- Only active in last 40 generations for performance
- Maintains Ne stability

### ✅ 3. **Dynamic Ne Changes**
```slim
s1 300: early() {
	t = slim_total_generations - sim.cycle;
	Nt = (N / N0)^(t / g_ne_change_start) * N0;
	
	// NEW: Applied to ALL populations
	for (pop in sim.subpopulations) {
		pop.setSubpopulationSize(asInteger(Nt));
	}
}
```
- Same exponential growth/decline formula
- Now applies to **all populations simultaneously**

### ✅ 4. **Selection with Restart Logic**
```slim
s2 450: late() {
	// Save state before selection
	if (sim.cycle == slim_total_generations - g_sel_start - 1 & s != 0.0) {
		sim.treeSeqOutput(state_file);
		
		// NEW: Introduce in first population only
		target_pop = sim.subpopulations[0];
		sample(target_pop.genomes, num_origins).addNewDrawnMutation(m2, selpos);
	}
	
	// Restart if lost
	else if ((mut.size() != 1) & restart_counter < max_restart) {
		sim.readFromPopulationFile(state_file);
		setSeed(rdunif(1, 0, asInteger(2^62 - 1)));
		sample(target_pop.genomes, num_origins).addNewDrawnMutation(m2, selpos);
		restart_counter = restart_counter + 1;
	}
}
```
- **Preserved restart logic** for selection establishment
- Mutation introduced in **population 1** initially
- Can spread to other populations via migration

### ✅ 5. **DAF Tracking**
```slim
// NEW: Track both global and per-population frequencies
global_freq = sim.mutationFrequencies(NULL, mut);
catn(c("DAF_global", time, global_freq), sep='\t');

for (pop in sim.subpopulations) {
	pop_freq = sim.mutationFrequencies(pop, mut);
	catn(c("DAF_pop" + pop.id, time, pop_freq), sep='\t');
}
```
- Tracks **global DAF** across all populations
- Tracks **per-population DAF** separately
- Shows mutation spread via migration

### ✅ 6. **True Ne Tracking**
```slim
late() {
	if (sim.cycle < slim_total_generations) {
		t_ago = slim_total_generations - sim.cycle - 1;
		
		// NEW: Track each population separately
		for (pop in sim.subpopulations) {
			catn(c('True_Ne_pop' + pop.id, t_ago, pop.individualCount), sep='\t');
		}
		
		// NEW: Track total Ne
		total_individuals = sum(sim.subpopulations.individualCount);
		catn(c('True_Ne_total', t_ago, total_individuals), sep='\t');
	}
}
```
- Tracks **Ne per population**
- Tracks **total Ne** across populations
- Useful for detecting population structure effects

---

## New Features Added

### 1. **Multiple Populations**
```slim
1 early() {
	// Create n_pops populations
	for (pop_id in 1:n_pops) {
		sim.addSubpop(pop_id, N);
	}
}
```
- Creates 3-5 populations (configurable)
- Each starts with same Ne

### 2. **Migration Matrix**
```slim
if (n_pops > 1) {
	for (i in 1:n_pops) {
		for (j in 1:n_pops) {
			if (i != j) {
				sim.subpopulations[i-1].setMigrationRates(j, mig_rate);
			}
		}
	}
}
```
- **Symmetric migration** between all populations
- Configurable migration rate
- Island model (all populations connected equally)

### 3. **Enhanced Output Tracking**
```slim
// Progress updates every 100 generations
100:999999 late() {
	if (sim.cycle % 100 == 0) {
		catn(c("Progress", "Generation", sim.cycle), sep='\t');
		
		for (pop in sim.subpopulations) {
			catn(c("  Pop" + pop.id + ":", pop.individualCount), sep='\t');
		}
		
		if (s != 0.0 & mut.size() > 0) {
			catn(c("  Selection frequency:", global_freq), sep='\t');
		}
	}
}
```

---

## Parameter Changes

### Original Parameters (All Preserved)
| Parameter | Default | Description |
|-----------|---------|-------------|
| `L` | 750000 | Chromosome length |
| `selpos` | L/3 | Selection position |
| `num_origins` | 1 | Initial mutation copies |
| `N` | 10000 | **Ancient Ne per deme** |
| `N0` | 1000 | **Final Ne per deme** |
| `h` | 0.5 | Dominance coefficient |
| `s` | 0.3 | Selection coefficient |
| `g_sel_start` | 80 | Selection start (gen ago) |
| `r` | 6.67e-7 | Recombination rate |
| `sim_relatedness` | F | Simulate inbreeding |
| `g_ne_change_start` | 200 | Ne change start |
| `max_restart` | 100 | Max selection restarts |

### New Migration Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_pops` | 3 | Number of populations |
| `mig_rate` | 0.01 | Migration rate (per gen) |
| `scenario_id` | "test" | Scenario identifier |
| `outdir` | "." | Output directory |

---

## Output Format Changes

### Original Output Lines
```
Parameter	N	10000
Parameter	s	0.3
True_Ne	199	9845
DAF	79	0.123
restart_counter	3
```

### New Output Lines
```
Parameter	N	10000
Parameter	n_pops	3
Parameter	mig_rate	0.01

# Per-population Ne tracking
True_Ne_pop1	199	3280
True_Ne_pop2	199	3295
True_Ne_pop3	199	3270
True_Ne_total	199	9845

# Per-population DAF tracking
DAF_global	79	0.123
DAF_pop1	79	0.245
DAF_pop2	79	0.089
DAF_pop3	79	0.034

restart_counter	3
```

---

## Python Integration

### Parsing SLiM Output

```python
def parse_slim_stdout(scenario_dir: Path):
    """Extract all tracked metrics"""
    metrics = {
        'daf_global': [],          # Global DAF trajectory
        'daf_by_pop': {},          # Per-population DAF
        'ne_by_pop': {},           # Per-population Ne
        'ne_total': [],            # Total Ne
        'restart_count': 0,
        'parameters': {}
    }
    
    # Parse tab-separated output
    for line in log_file:
        parts = line.split('\t')
        
        if parts[0] == 'DAF_global':
            metrics['daf_global'].append({
                'generation': int(parts[1]),
                'frequency': float(parts[2])
            })
        
        elif parts[0].startswith('DAF_pop'):
            pop_id = parts[0].replace('DAF_pop', '')
            metrics['daf_by_pop'][pop_id].append({
                'generation': int(parts[1]),
                'frequency': float(parts[2])
            })
        
        # ... similar for Ne tracking
    
    return metrics
```

### Saved Outputs Per Scenario

```
scenario_dir/
├── parameters.json                      # All parameters
├── slim_output_S001_rep1.trees         # Tree sequence from SLiM
├── slim_output_S001_rep1_processed.trees  # After recapitation
├── slim_output_S001_rep1_restart_count.txt
├── slim_output_S001_rep1_daf_global.tsv
├── slim_output_S001_rep1_daf_pop1.tsv
├── slim_output_S001_rep1_daf_pop2.tsv
├── slim_output_S001_rep1_daf_pop3.tsv
├── slim_output_S001_rep1_ne_pop1.tsv
├── slim_output_S001_rep1_ne_pop2.tsv
├── slim_output_S001_rep1_ne_pop3.tsv
├── slim_output_S001_rep1_ne_total.tsv
└── slim_output_S001_rep1_slim_metrics.json  # All metrics in JSON
```

---

## Usage Example

### Command Line Execution
```bash
slim \
  -d "N=10000" \
  -d "N0=1000" \
  -d "n_pops=3" \
  -d "mig_rate=0.01" \
  -d "s=0.3" \
  -d "h=0.5" \
  -d "g_sel_start=80" \
  -d "r=6.67e-7" \
  -d "sim_relatedness=F" \
  -d "outid=1" \
  -d "scenario_id='S001'" \
  -d "outdir='results/S001/rep_001'" \
  -seed 12345 \
  migration_model.slim
```

### From Python
```python
import subprocess

cmd = [
    "slim",
    "-d", f"N={params.Ne}",
    "-d", f"N0={params.N0}",
    "-d", f"n_pops={params.n_populations}",
    "-d", f"mig_rate={params.migration_rate}",
    "-d", f"s={params.s}",
    "-d", f"h={params.h}",
    "-d", f"g_sel_start={params.g_sel_start}",
    "-d", f"r={params.rec_rate}",
    "-d", f"sim_relatedness={'T' if params.selfing_rate > 0 else 'F'}",
    "-d", f"outid={params.replicate_id}",
    "-d", f"scenario_id='{params.scenario_id}'",
    "-d", f"outdir='{scenario_dir}'",
    "-seed", str(seed),
    "migration_model.slim"
]

subprocess.run(cmd, check=True, capture_output=True)
```

---

## Biological Interpretation

### Single Population (Original)
```
Population 1: [====== 10,000 individuals ======]
              ↓ selection at gen 80
              ↓ selected allele rises
              ↓ Ne declines to 1,000
```

### Multi-Population with Migration (New)
```
Pop 1: [===== 10k =====]  ← mutation introduced here
       ↓ ↘ ↙              ← migration
Pop 2: [===== 10k =====]  ← mutation spreads here
       ↙ ↓ ↘
Pop 3: [===== 10k =====]  ← mutation spreads here

Time →
  Selection starts (gen 80)
  Mutation rises in Pop 1 (DAF_pop1 increases)
  Spreads to Pop 2 & 3 via migration
  All Ne decline to 1,000
  Track frequencies separately per population
```

### Key Differences
1. **Mutation origin**: Starts in Pop 1, spreads via migration
2. **DAF heterogeneity**: Different frequencies in different populations
3. **Population structure**: FST between populations
4. **Sample distribution**: Can sample from multiple populations

---

## Testing the Script

### Test 1: Single Population (Backward Compatibility)
```bash
# Should behave identically to original
slim \
  -d "n_pops=1" \
  -d "N=10000" \
  -d "s=0.3" \
  migration_model.slim
```

### Test 2: Multiple Populations, No Selection
```bash
# Test migration without selection
slim \
  -d "n_pops=3" \
  -d "mig_rate=0.01" \
  -d "s=0.0" \
  migration_model.slim
```

### Test 3: Selection with Migration
```bash
# Full model
slim \
  -d "n_pops=3" \
  -d "mig_rate=0.01" \
  -d "s=0.3" \
  -d "g_sel_start=80" \
  migration_model.slim
```

---

## Validation

### Check Output Files
```bash
# Tree sequence created
ls -lh slim_output_*.trees

# Check it's valid
python3 -c "import tskit; ts = tskit.load('slim_output_test_1.trees'); print(f'Samples: {ts.num_samples}, Populations: {ts.num_populations}')"
```

### Verify DAF Tracking
```bash
# Extract DAF lines
grep "^DAF" slim_output.log

# Should see both global and per-population
```

### Check Ne Tracking
```bash
# Extract Ne lines
grep "^True_Ne" slim_output.log

# Should see per-population and total
```

---

## Summary of Adaptations

| Feature | Original | Migration Model | Status |
|---------|----------|-----------------|--------|
| `set_default()` | ✓ | ✓ | Preserved |
| `modifyChild()` relatedness | ✓ | ✓ | Preserved exactly |
| Dynamic Ne changes | ✓ | ✓ (all pops) | Enhanced |
| Selection + restart | ✓ | ✓ (pop 1) | Enhanced |
| DAF tracking | Single | Global + per-pop | Enhanced |
| Ne tracking | Single | Per-pop + total | Enhanced |
| Multiple populations | ✗ | ✓ | **NEW** |
| Migration | ✗ | ✓ | **NEW** |
| Output parsing | Simple | Structured | Enhanced |

**All critical features preserved, multi-population functionality added!**
