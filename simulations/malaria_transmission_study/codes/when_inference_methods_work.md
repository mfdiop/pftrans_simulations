
When do transmission inference methods work, and what biological/sampling factors matter most?

The framework is designed to answer: When do transmission inference methods work, and what biological/sampling factors matter most?

I'll help you design a comprehensive simulation framework for studying malaria transmission with migration between populations. This will build on your existing code structure while incorporating the complex parameter space you've outlined.

Now let me create a comprehensive documentation document explaining the simulation design:

Now, let me create a helper script to analyze and visualize the simulation design:

Perfect! Now let me create a quick reference guide for using the system:

## Summary

I've created a comprehensive simulation framework for your multi-population malaria transmission study. Here's what you have:

### 📦 **Three Main Components:**

1. **`malaria_migration_simulation.py`** - Core simulation framework
   - Defines complete parameter space (5×3×5 core biological scenarios)
   - Three design types: full (9,000 runs), migration-focused (6,480 runs), core (750 runs)
   - Generates SLURM array jobs automatically
   - Integrates with your existing code

2. **`analyze_simulation_design.py`** - Design analysis tool
   - Generates summary statistics
   - Creates visualizations of parameter space
   - Produces batch processing summaries
   - Exports scenario lists

3. **Documentation**
   - `SIMULATION_DESIGN.md` - Complete design rationale
   - `QUICK_START.md` - Step-by-step usage guide

### 🎯 **Key Features:**

**Parameter Space:**
- **Recombination rates**: 5 levels (1e-9 to 1e-5) - proxy for within-host diversity
- **Bottleneck sizes**: 3 levels (1, 5, 20 parasites)
- **Expected substitutions**: 5 levels (0 to 2.0 per transmission)
- **Sampling proportions**: 4 levels (10%, 30%, 60%, 90%)
- **Outbreak sizes**: 3 levels (50, 200, 800 infections)
- **Population structure**: 3-5 populations with variable migration

**Design Flexibility:**
- Start small with 750 core scenarios
- Scale to 9,000 comprehensive scenarios
- 10-20 replicates per scenario for robust statistics
- Easy SLURM integration with automatic array job generation

### 🚀 **Getting Started:**

```bash
# 1. Generate design
python malaria_migration_simulation.py \
  --generate-design core --n-replicates 10 --outdir sims

# 2. Analyze design  
python analyze_simulation_design.py sims/simulation_design.json

# 3. Test
python malaria_migration_simulation.py --test-run --dry-run

# 4. Submit
cd sims && sbatch run_migration_sim.sh
```

The framework is designed to answer: **When do transmission inference methods work, and what biological/sampling factors matter most?**