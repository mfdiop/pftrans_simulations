#!/bin/bash
#SBATCH --job-name=relationship_pipeline
#SBATCH --output=logs/pipeline_%A_%a.out
#SBATCH --error=logs/pipeline_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-60%10  # 60 combinations, 10 concurrent jobs

# Set parameters
BASE_DIR="results/multiple_runs"
OUTDIR="results_multi_chrom"
REPLICATES="1 2 3 4 5"
REC_RATES="1e-9 1e-8 1e-7 1e-6" 
CHROMOSOMES="1 2 3"

# Calculate total combinations
TOTAL_COMBINATIONS=$((5 * 3 * 3))  # 5 reps * 4 rates * 3 chroms

# Run R script
Rscript -e "
source('multi_chromosome_pipeline.R')
run_multi_chromosome_pipeline(
  base_dir = '$BASE_DIR',
  replicates = c($REPLICATES),
  recombination_rates = c($REC_RATES),
  chromosomes = c($CHROMOSOMES),
  outdir = '$OUTDIR',
  slurm_job_id = $SLURM_ARRAY_TASK_ID
)
"

echo "Job $SLURM_ARRAY_TASK_ID completed"