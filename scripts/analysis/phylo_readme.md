# Phylogenetic Analysis Pipeline

## Overview

This pipeline performs phylogenetic analysis on VCF files from evolutionary simulations:
1. Filters VCF (optional)
2. Converts VCF to FASTA
3. Runs multiple sequence alignment (optional)
4. Selects best substitution model (ModelFinder)
5. Builds phylogenetic tree (IQ-TREE with ultrafast bootstrap)

---

## Key Improvements in This Version

### 1. **Fixed Critical Bugs**
- ✅ Fixed unclosed `if` statement (missing `fi`)
- ✅ Fixed VCF to FASTA logic (append-then-pop → collect-then-append)
- ✅ Corrected bcftools MAF filtering syntax
- ✅ Proper error handling with try-catch in Python
- ✅ Color codes disabled in log files (terminal detection)

### 2. **Command-Line Arguments**
All parameters can now be overridden from command line:
```bash
sbatch script.slurm --bootstrap 2000 --min-maf 0.05 --use-iupac
```

### 3. **Better Error Handling**
- Robust genotype parsing (handles `.`, `./.`, `.|.`, haploid)
- Proper pipeline status checking (`${PIPESTATUS[0]}`)
- Detailed error messages with context

### 4. **Enhanced Reporting**
- Statistics on filtered sites
- Summary reports for each analysis
- Validation script to check all results

---

## Installation

### Required Software

```bash
# IQ-TREE (any version)
conda install -c bioconda iqtree

# Or download from: http://www.iqtree.org/

# Python with BioPython
pip install biopython

# Optional: VCF filtering tools
conda install -c bioconda bcftools vcftools

# Optional: Alignment tools
conda install -c bioconda mafft muscle
```

---

## Usage

### Basic Usage (Default Parameters)

```bash
# Submit array job
sbatch 04_run_phylogenetic_batch_improved.slurm
```

### Custom Parameters

```bash
# More stringent filtering
sbatch 04_run_phylogenetic_batch_improved.slurm \
  --min-maf 0.05 \
  --max-missing 0.9

# More bootstrap replicates
sbatch 04_run_phylogenetic_batch_improved.slurm \
  --bootstrap 2000

# Use IUPAC codes for heterozygotes
sbatch 04_run_phylogenetic_batch_improved.slurm \
  --use-iupac

# Run alignment before tree building
sbatch 04_run_phylogenetic_batch_improved.slurm \
  --run-alignment \
  --alignment-tool mafft

# Custom directories
sbatch 04_run_phylogenetic_batch_improved.slurm \
  --input-dir /path/to/vcfs \
  --output-dir /path/to/output

# Combine multiple options
sbatch 04_run_phylogenetic_batch_improved.slurm \
  --bootstrap 2000 \
  --min-maf 0.05 \
  --max-missing 0.95 \
  --threads 16 \
  --use-iupac
```

### All Available Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input-dir` | string | `results/multiple_runs` | Input VCF directory |
| `--output-dir` | string | `results/multiple_runs/phylo_results` | Output directory |
| `--bootstrap` | int | `1000` | Bootstrap replicates |
| `--threads` | int | `${SLURM_CPUS_PER_TASK}` | Number of threads |
| `--model-set` | string | `DNA` | Model set to test |
| `--min-maf` | float | `0.0` | Minimum minor allele frequency |
| `--max-missing` | float | `1.0` | Maximum missing data (0-1) |
| `--use-iupac` | flag | `false` | Use IUPAC codes for heterozygotes |
| `--min-samples` | int | `10` | Minimum non-missing samples per site |
| `--run-alignment` | flag | `false` | Run multiple sequence alignment |
| `--alignment-tool` | string | `mafft` | Alignment tool (`mafft` or `muscle`) |

---

## Checking Results

After jobs complete, validate results:

```bash
# Check all results
bash 05_check_phylo_results.sh

# Check specific directory
bash 05_check_phylo_results.sh /path/to/phylo_results
```

This generates:
- Summary of complete/incomplete jobs
- Error analysis from log files
- File size statistics
- Recommendations for next steps

---

## Output Structure

```
results/multiple_runs/phylo_results/
├── rep1/
│   ├── run1_rec1e09_chr1/
│   │   ├── run1_rec1e09_chr1.fasta              # Converted sequences
│   │   ├── run1_rec1e09_chr1.treefile           # Phylogenetic tree
│   │   ├── run1_rec1e09_chr1.iqtree             # Full IQ-TREE report
│   │   ├── run1_rec1e09_chr1_modelfinder.log    # Model selection log
│   │   ├── run1_rec1e09_chr1.log                # Execution log
│   │   └── run1_rec1e09_chr1_report.txt         # Summary report
│   ├── run1_rec1e08_chr1/
│   └── ...
├── rep2/
├── ...
├── vcf_to_fasta_batch.py                        # Conversion script
└── analysis_summary.txt                          # Overall summary
```

---

## Understanding the Output

### 1. Tree File (`.treefile`)
Newick format phylogenetic tree with:
- Branch lengths (substitutions per site)
- Bootstrap support values at nodes
- Sample names at tips

```newick
((sample1:0.001,sample2:0.002)100:0.005,(sample3:0.003,sample4:0.001)95:0.004);
```

### 2. IQ-TREE Report (`.iqtree`)
Comprehensive report including:
- Model selection results
- Log-likelihood
- Tree length
- Bootstrap proportions
- Execution time and parameters

### 3. FASTA File (`.fasta`)
Aligned sequences used for tree building:
```
>sample1
ATCGATCGATCG...
>sample2
ATCGATCGATCG...
```

---

## Troubleshooting

### Problem: "No valid sites found in VCF"

**Cause**: All sites filtered out or VCF is empty

**Solutions**:
```bash
# Check VCF has data
zcat your_file.vcf.gz | grep -v "^#" | head

# Reduce filtering stringency
sbatch script.slurm --min-maf 0.0 --max-missing 1.0

# Lower minimum samples requirement
sbatch script.slurm --min-samples 5
```

### Problem: "Not enough sequences for phylogenetic analysis"

**Cause**: Fewer than 4 samples in VCF

**Solution**: Check your VCF file has enough samples
```bash
# Count samples in VCF
bcftools query -l your_file.vcf | wc -l
```

### Problem: IQ-TREE fails with "No informative sites"

**Cause**: All samples have identical sequences

**Solutions**:
- Check that VCF contains variants
- Verify mutation rate in simulation was non-zero
- Check that sites weren't over-filtered

### Problem: Jobs timing out

**Cause**: Bootstrap or dataset too large

**Solutions**:
```bash
# Increase time limit in SBATCH header
#SBATCH --time=120:00:00

# Reduce bootstrap replicates
sbatch script.slurm --bootstrap 500

# Use more threads
#SBATCH --cpus-per-task=16
sbatch script.slurm --threads 16
```

### Problem: Memory errors

**Cause**: Large datasets

**Solutions**:
```bash
# Increase memory in SBATCH header
#SBATCH --mem=64G

# Filter more stringently
sbatch script.slurm --min-maf 0.05 --max-missing 0.9
```

---

## Best Practices

### 1. **Filtering VCF**
```bash
# For high-quality phylogenies, filter variants:
sbatch script.slurm \
  --min-maf 0.05 \
  --max-missing 0.9
```

### 2. **Bootstrap Replicates**
- **Standard**: 1000 (good for most purposes)
- **Publication**: 2000+ (more robust support values)
- **Quick test**: 100 (fast but less reliable)

### 3. **IUPAC Codes**
- **Use `--use-iupac`** if heterozygotes are real (diploid organisms)
- **Don't use** if pseudo-haploidizing makes biological sense (malaria)

### 4. **Alignment**
- **Skip** (`--run-alignment` off) for SNP data from VCF (already aligned)
- **Use** if you have unaligned sequences from other sources

### 5. **Model Selection**
- `DNA` (default): Tests all DNA models
- `GTR`: Only GTR family models
- `JC,HKY,GTR`: Test specific models

---

## Example Workflows

### Workflow 1: Quick Test Run
```bash
# Test on subset of data
sbatch --array=1-5 script.slurm \
  --bootstrap 100 \
  --threads 4
```

### Workflow 2: High-Quality Analysis
```bash
# Full analysis with stringent filtering
sbatch script.slurm \
  --bootstrap 2000 \
  --min-maf 0.05 \
  --max-missing 0.9 \
  --threads 16
```

### Workflow 3: Diploid Data
```bash
# Preserve heterozygous information
sbatch script.slurm \
  --use-iupac \
  --bootstrap 1000
```

### Workflow 4: Sequential Mode (No SLURM)
```bash
# Run locally without SLURM array
bash script.slurm \
  --bootstrap 500 \
  --threads 8
```

---

## Downstream Analysis

### Visualize Trees in R

```r
library(ape)
library(phytools)

# Read tree
tree <- read.tree("results/.../run1_rec1e09_chr1.treefile")

# Plot
plot(tree, cex=0.8)
nodelabels(tree$node.label, cex=0.6, bg="lightblue")

# Extract distances
dist_matrix <- cophenetic.phylo(tree)
```

### Visualize with FigTree
1. Download: http://tree.bio.ed.ac.uk/software/figtree/
2. Open `.treefile`
3. Display bootstrap values: Node Labels → Display
4. Color branches, adjust layout, export figures

### Compare Trees Across Replicates

```r
# Compare topologies
library(TreeDist)

trees <- lapply(tree_files, read.tree)
robinson_foulds <- RobinsonFoulds(trees[[1]], trees[[2]])
```

---

## Performance Optimization

### Time Estimates (per job)

| Dataset Size | Samples | Sites | Bootstrap | Time |
|--------------|---------|-------|-----------|------|
| Small | 50 | 1K | 1000 | 5-10 min |
| Medium | 100 | 10K | 1000 | 30-60 min |
| Large | 200 | 50K | 1000 | 2-4 hours |
| Very Large | 500 | 100K | 2000 | 8-12 hours |

### Resource Recommendations

```bash
# Small datasets
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00

# Medium datasets (recommended)
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# Large datasets
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
```

---

## Citation

If using this pipeline, please cite:

**IQ-TREE**:
- Nguyen et al. (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol, 32:268-274.

**ModelFinder**:
- Kalyaanamoorthy et al. (2017) ModelFinder: Fast model selection for accurate phylogenetic estimates. Nat Methods, 14:587-589.

**UFBoot**:
- Hoang et al. (2018) UFBoot2: Improving the ultrafast bootstrap approximation. Mol Biol Evol, 35:518-522.

---

## Support

For issues:
1. Check log files in `phylo_logs/`
2. Run validation script: `bash 05_check_phylo_results.sh`
3. Review individual job logs in output directories
4. Check IQ-TREE documentation: http://www.iqtree.org/doc/

---

## Changelog

### Version 2.0 (2025-11-20)
- ✅ Fixed all bash syntax errors (unclosed if/fi)
- ✅ Fixed VCF to FASTA conversion logic
- ✅ Added command-line argument parsing
- ✅ Improved error handling throughout
- ✅ Added terminal detection for colored output
- ✅ Enhanced Python script with better stats
- ✅ Created validation/checking script
- ✅ Comprehensive documentation

### Version 1.0 (2025-11-05)
- Initial version
