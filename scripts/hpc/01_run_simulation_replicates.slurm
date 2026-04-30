#!/bin/bash

#SBATCH --job-name=pf_sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=1-60%4                    # 4 rec_rates × 5 replicates = 20 jobs

###############################################################################
# Script: 01_run_simulation_replicates.slurm
# Purpose: Run P. falciparum simulations across multiple recombination rates
#          with multiple replicates per rate
# Author: [Your Name]
# Date: 2025-11-05
#
# Description:
#   Runs SLiM-msprime simulation pipeline across different recombination rates
#   with multiple replicates. Each array job processes one replicate of one
#   recombination rate.
#
# Usage:
#   # Edit configuration section below, then:
#   sbatch codes/01_run_simulation_replicates.slurm
#
#   # Or run single replicate/rec_rate combination:
#   bash codes/01_run_simulation_replicates.slurm <rec_rate_index> <replicate_number>
#
#   Parallel execution (recommended)
#   sbatch codes/01_run_simulation_replicates.slurm
#
#   Single replicate/rec_rate
#   bash codes/01_run_simulation_replicates.slurm 1 1  # rec_rate_index=1, replicate=1
#
#   Sequential (all combinations)
#   bash codes/01_run_simulation_replicates.slurm
#
#   With custom parameters
#   N_CHROMOSOMES=14 N_SAMPLES=200 sbatch run_simulation_replicates.slurm
#
# Output structure:
#   out/
#   ├── rep1/
#   │   ├── run1_rec1.0e-09_chr*.trees
#   │   ├── run1_rec1.0e-09_chr*.vcf
#   │   └── ...
#   ├── rep2/
#   │   └── ...
#   └── rep5/
#       └── ...
###############################################################################

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

###############################################################################
# CONFIGURATION
###############################################################################

# Python script location
#PYTHON_SCRIPT="codes/01_run_simulation_replicates.py"
PYTHON_SCRIPT="codes/01_run_simulation_replicates_for_unrelated_indvs.py"
# Check if Python script exists
if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
    echo "ERROR: Python script not found: ${PYTHON_SCRIPT}"
    exit 1
fi

# Base output directory
#BASE_OUTDIR="results/scenario"

# Create directories
#mkdir -p "${BASE_OUTDIR}" logs

# Recombination rates to simulate
declare -a REC_RATES=(
    "1.0e-09"
    "1.0e-08"
    "1.0e-07"
    "1.0e-06"
)

# Number of replicates per recombination rate
N_REPLICATES=5

# Simulation parameters (can be overridden by command line)
NE=${NE:-10000}                          # Effective population size
N0=${N0:-1000}                           # Current population size
N_SAMPLES=${N_SAMPLES:-1000}              # Number of samples
GENERATIONS=${GENERATIONS:-1000}         # Forward simulation generations
N_CHROMOSOMES=${N_CHROMOSOMES:-3}        # Number of chromosomes (use 3 for testing, 14 for full)
SELFING_RATE=${SELFING_RATE:-0.95}       # Selfing rate
NUM_ORIGINS=${NUM_ORIGINS:-1}            # Number of mutation origins
H=${H:-0.5}                              # Dominance coefficient
S=${S:-0.3}                              # Selection coefficient
G_SEL_START=${G_SEL_START:-80}           # Selection start generation
SIM_RELATEDNESS=${SIM_RELATEDNESS:-0}    # Simulate relatedness
G_NE_CHANGE=${G_NE_CHANGE:-200}          # Ne change start generation
MUTATION_RATE=${MUTATION_RATE:-1e-8}     # Mutation rate
BASE_OUTDIR=${BASE_OUTDIR:-"out"}

mkdir -p "$BASE_OUTDIR" logs

# Environment setup (uncomment and modify as needed)
module load python/3.12.2
# source activate slim-env
export PATH="/mnt/scratch/fadel/PhD/Objective1/methods_evaluation/build/slim:$PATH"

###############################################################################
# FUNCTIONS
###############################################################################

# Function to log messages with timestamp
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Function to run simulation for one replicate and recombination rate
run_simulation() {
    local rec_rate=$1
    local replicate=$2
    local run_id=$3
    
    # Create output directory for this replicate
    local output_dir="${BASE_OUTDIR}/rep${replicate}"
    mkdir -p "${output_dir}"
    
    log_message "=========================================="
    log_message "Starting simulation"
    log_message "=========================================="
    log_message "Replicate: ${replicate}"
    log_message "Recombination rate: ${rec_rate}"
    log_message "Run ID: ${run_id}"
    log_message "Output directory: ${output_dir}"
    log_message "=========================================="
    log_message "Parameters:"
    log_message "  Ne: ${NE}"
    log_message "  N0: ${N0}"
    log_message "  N_samples: ${N_SAMPLES}"
    log_message "  Generations: ${GENERATIONS}"
    log_message "  N_chromosomes: ${N_CHROMOSOMES}"
    log_message "  Mutation_rate: ${MUTATION_RATE}"
    log_message "  Selfing_rate: ${SELFING_RATE}"
    log_message "  Selection (s): ${S}"
    log_message "  Dominance (h): ${H}"
    log_message "=========================================="
    
    # Build Python command
    python_cmd="python ${PYTHON_SCRIPT}"
    python_cmd+=" --Ne ${NE}"
    python_cmd+=" --N0 ${N0}"
    python_cmd+=" --n_samples ${N_SAMPLES}"
    python_cmd+=" --generations ${GENERATIONS}"
    python_cmd+=" --n_chromosomes ${N_CHROMOSOMES}"
    python_cmd+=" --selfing_rate ${SELFING_RATE}"
    python_cmd+=" --num_origins ${NUM_ORIGINS}"
    python_cmd+=" --h ${H}"
    python_cmd+=" --s ${S}"
    python_cmd+=" --g_sel_start ${G_SEL_START}"
    python_cmd+=" --sim_relatedness ${SIM_RELATEDNESS}"
    python_cmd+=" --g_ne_change_start ${G_NE_CHANGE}"
    python_cmd+=" --mutation_rate ${MUTATION_RATE}"
    python_cmd+=" --rec_rate ${rec_rate}"
    python_cmd+=" --run_id ${run_id}"
    python_cmd+=" --outdir ${output_dir}"
    
    log_message "Command: ${python_cmd}"
    log_message ""
    
    # Run Python script
    if eval "${python_cmd}"; then
        log_message "=========================================="
        log_message "SUCCESS: Simulation completed"
        log_message "Replicate: ${replicate}, Rec rate: ${rec_rate}"
        log_message "Results saved to: ${output_dir}"
        log_message "=========================================="
        return 0
    else
        log_message "=========================================="
        log_message "ERROR: Simulation failed"
        log_message "Replicate: ${replicate}, Rec rate: ${rec_rate}"
        log_message "=========================================="
        return 1
    fi
}

###############################################################################
# MAIN EXECUTION
###############################################################################

log_message "=========================================="
log_message "    P. falciparum Simulation Pipeline     "
log_message "=========================================="

# Calculate total number of jobs
TOTAL_JOBS=$((${#REC_RATES[@]} * 3 * N_REPLICATES))

# Parse command line arguments (for standalone mode)
if [[ $# -ge 2 ]]; then
    # Standalone mode with explicit arguments
    REC_RATE_IDX=$1
    REPLICATE=$2
    
    # Validate indices
    if [[ ${REC_RATE_IDX} -lt 1 || ${REC_RATE_IDX} -gt ${#REC_RATES[@]} ]]; then
        echo "ERROR: Invalid recombination rate index: ${REC_RATE_IDX}"
        echo "Valid range: 1-${#REC_RATES[@]}"
        exit 1
    fi
    
    if [[ ${REPLICATE} -lt 1 || ${REPLICATE} -gt ${N_REPLICATES} ]]; then
        echo "ERROR: Invalid replicate number: ${REPLICATE}"
        echo "Valid range: 1-${N_REPLICATES}"
        exit 1
    fi
    
    # Get recombination rate
    REC_RATE="${REC_RATES[$((REC_RATE_IDX - 1))]}"
    
    # Calculate run_id (unique across all combinations)
    RUN_ID=$(( (REC_RATE_IDX - 1) * N_REPLICATES + REPLICATE ))
    
    log_message "Running in standalone mode"
    log_message "Rec rate index: ${REC_RATE_IDX}, Replicate: ${REPLICATE}"
    
    run_simulation "${REC_RATE}" "${REPLICATE}" "${RUN_ID}"
    
elif [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # SLURM array job mode
    log_message "Running as SLURM array job: Task ${SLURM_ARRAY_TASK_ID}"
    
    if [[ ${SLURM_ARRAY_TASK_ID} -lt 1 || ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_JOBS} ]]; then
        log_message "ERROR: Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
        log_message "Valid range: 1-${TOTAL_JOBS}"
        exit 1
    fi
    
    # Calculate which combination to process
    # Array layout: rec_rate1_rep1, rec_rate1_rep2, ..., rec_rate2_rep1, ...
    TASK_ID=${SLURM_ARRAY_TASK_ID}
    
    # Calculate rec_rate index and replicate number
    REC_RATE_IDX=$(( (TASK_ID - 1) / N_REPLICATES + 1 ))
    REPLICATE=$(( (TASK_ID - 1) % N_REPLICATES + 1 ))
    
    # Get recombination rate
    REC_RATE="${REC_RATES[$((REC_RATE_IDX - 1))]}"
    
    # Use task ID as run_id for uniqueness
    RUN_ID=${TASK_ID}
    
    log_message "Processing combination ${TASK_ID}/${TOTAL_JOBS}"
    log_message "Rec rate: ${REC_RATE} (index ${REC_RATE_IDX})"
    log_message "Replicate: ${REPLICATE}"
    
    # Run simulation for this combination
    run_simulation "${REC_RATE}" "${REPLICATE}" "${RUN_ID}"
    
else
    # Sequential mode - process all combinations
    log_message "Running in sequential mode (processing all combinations)"
    log_message "Total jobs: ${TOTAL_JOBS}"
    
    CURRENT_JOB=0
    FAILED_JOBS=0
    RUN_ID=1
    
    # Loop through all recombination rates
    for rec_rate_idx in "${!REC_RATES[@]}"; do
        REC_RATE="${REC_RATES[$rec_rate_idx]}"
        
        log_message ""
        log_message "=========================================="
        log_message "Recombination rate: ${REC_RATE}"
        log_message "=========================================="
        
        # Loop through all replicates
        for rep in $(seq 1 ${N_REPLICATES}); do
            CURRENT_JOB=$((CURRENT_JOB + 1))
            
            log_message ""
            log_message "Job ${CURRENT_JOB}/${TOTAL_JOBS}: Rep ${rep}, Rec ${REC_RATE}"
            
            # Run simulation
            if ! run_simulation "${REC_RATE}" "${rep}" "${RUN_ID}"; then
                FAILED_JOBS=$((FAILED_JOBS + 1))
                log_message "WARNING: Job ${CURRENT_JOB} failed, continuing..."
            fi
            
            RUN_ID=$((RUN_ID + 1))
            log_message ""
        done
    done
    
    log_message "=========================================="
    log_message "Sequential processing complete"
    log_message "Successful jobs: $((TOTAL_JOBS - FAILED_JOBS))/${TOTAL_JOBS}"
    log_message "Failed jobs: ${FAILED_JOBS}/${TOTAL_JOBS}"
    log_message "=========================================="
    
    if [[ ${FAILED_JOBS} -gt 0 ]]; then
        exit 1
    fi
fi

log_message "=========================================="
log_message "Simulation pipeline completed"
log_message "=========================================="

exit 0
