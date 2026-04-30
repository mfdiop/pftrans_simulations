#!/bin/bash
#SBATCH --job-name=single_pop_sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-20
#SBATCH --partition=main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mdiop@mrc.gm

################################################################################
# SLURM Script for Single Population Evolutionary Simulation
################################################################################
#
# PURPOSE:
#   This script runs forward-time evolutionary simulations using SLiM and 
#   processes the output using tskit/msprime/pyslim via an R wrapper script.
#   It executes multiple replicates (job array) to enable statistical analysis
#   with confidence intervals across replicate simulations.
#
# WORKFLOW:
#   1. Environment Setup - Load required modules and conda environments
#   2. Parameter Parsing - Read simulation parameters from command line
#   3. Directory Creation - Set up output directories for this run
#   4. Validation - Check for required files (R script, SLiM script)
#   5. Simulation Execution - Run R script which:
#      a. Executes SLiM forward simulation
#      b. Loads tree sequences from SLiM output
#      c. Recapitates (adds ancient ancestry) using msprime
#      d. Mutates tree sequences
#      e. Generates VCF, tree files, and summary statistics
#   6. Results Collection - Aggregate outputs across replicates
#   7. Cleanup - Remove temporary files if requested
#
# INPUT PARAMETERS:
#   --chrno         : Chromosome number/ID (default: 1)
#   --seqlen        : Sequence length in bp (default: 1,500,000)
#   --selpos        : Position of selected site (default: 0.33*seqlen)
#   --num_origins   : Number of independent mutation origins (default: 1)
#   --N             : Population size (default: 10,000)
#   --h             : Dominance coefficient (default: 0.5)
#   --s             : Selection coefficient (default: 0.3)
#   --g_sel_start   : Generation when selection starts (default: 80)
#   --r             : Recombination rate per bp (default: 0.01/15000)
#   --sim_relatedness : Flag for simulating relatedness (default: 0)
#   --g_ne_change_start : Generation for Ne change (default: 200)
#   --N0            : Initial/ancestral population size (default: 1000)
#   --u             : Mutation rate per bp per generation (default: 1e-8)
#   --nsam          : Number of sampled haploid genomes (default: 1000)
#   --test          : Run in test mode (smaller parameters)
#   --slim_script   : Path to SLiM script file
#   --slim_bin      : SLiM executable name (default: slim or slim.exe)
#   --r_script      : Path to R wrapper script (default: auto-detect)
#   --outdir        : Base output directory (default: ./results)
#   --cleanup       : Remove intermediate files after success (default: 0)
#
# OUTPUT FILES (per replicate, per chromosome):
#   {genome_set_id}_{chrno}.restart_count - Number of SLiM restarts
#   {genome_set_id}_{chrno}.true_ne       - Effective population size over time
#   {genome_set_id}_{chrno}.daf           - Derived allele frequency trajectory
#   {genome_set_id}_{chrno}.trees         - Tree sequence file (tskit format)
#   {genome_set_id}_{chrno}.vcf.gz        - Genotype data in VCF format
#   {genome_set_id}_{chrno}.sfs.tsv       - Site frequency spectrum
#   {genome_set_id}_{chrno}.*.png         - Visualization plots
#
# SLURM ARRAY USAGE:
#   The script uses SLURM job arrays to run multiple replicates in parallel.
#   Each array task (SLURM_ARRAY_TASK_ID) represents one replicate simulation.
#   The genome_set_id is set to the array task ID, ensuring unique outputs.
#
# REQUIREMENTS:
#   - R (with packages: optparse, reticulate, readr, dplyr, stringr, tidyr)
#   - Python environment with: numpy, pandas, tskit, pyslim, msprime
#   - SLiM (forward evolutionary simulator)
#   - Sufficient disk space for tree sequences and VCF files
#
# EXAMPLE USAGE:
#   # Submit job with default parameters
#   sbatch codes/single_run/01_slim_msprime_wrapper.slurm
#
#   # Submit with custom parameters
#   sbatch codes/single_run/01_slim_msprime_wrapper.slurm --N 50000 --s 0.5 --nsam 2000
#
#   # Run test mode with 5 replicates
#   sbatch --array=1-5 codes/single_run/01_slim_msprime_wrapper.slurm --test
#
#   # Specify custom scripts and output directory
#   sbatch codes/single_run/01_slim_msprime_wrapper.slurm \
#     --slim_script /path/to/custom.slim \
#     --r_script /path/to/custom.R \
#     --outdir /scratch/project/sims
#
# TROUBLESHOOTING:
#   - Check log files in logs/sim_{JOBID}_{ARRAYID}.{out,err}
#   - Verify conda environment has all required Python packages
#   - Ensure SLiM binary is in PATH or specify with --slim_bin
#   - For memory errors, increase --mem in SBATCH header
#   - For timeout, increase --time in SBATCH header
#
# AUTHOR: Mouhamadou Fadel DIOP
# DATE: $(date '+%Y-%m-%d')
# VERSION: 1.0
################################################################################

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# Function to print colored messages
log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

log_step() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$*${NC}"
    echo -e "${BLUE}========================================${NC}"
}

log_success() {
    echo -e "${MAGENTA}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

################################################################################
# Error Handler
################################################################################
error_exit() {
    log_error "$1"
    log_error "Job failed. Check log files for details."
    exit 1
}

trap 'error_exit "Script failed at line $LINENO"' ERR

################################################################################
# Print Job Information
################################################################################
log_step "SLURM Job Information"
log_info "Job ID: ${SLURM_JOB_ID:-N/A}"
log_info "Array Task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
log_info "Node: ${SLURM_NODELIST:-N/A}"
log_info "CPUs: ${SLURM_CPUS_PER_TASK:-N/A}"
log_info "Memory: ${SLURM_MEM_PER_NODE:-N/A}MB"
log_info "Working Directory: $(pwd)"

################################################################################
# Default Parameters
################################################################################
CHRNO=1
SEQLEN=1500000             # 100 * 15000
SELPOS=495000              # as.integer(0.33 * 100 * 15000)
NUM_ORIGINS=1
N=10000
H=0.5
S=0.3
G_SEL_START=80
R=6.66666666667            # 0.01 / 15000
SIM_RELATEDNESS=0
G_NE_CHANGE_START=200
N0=1000
U=1e-8
NSAM=1000
OUTDIR="./results"
TEST=0
SLIM_SCRIPT="codes/single_run/single_pop.slim"
SLIM_BIN="build/slim"
R_SCRIPT="codes/single_run/01_slim_msprime_wrapper.R"
CLEANUP=0

################################################################################
# PARSE COMMAND LINE ARGUMENTS
################################################################################
log_step "Parsing Command Line Arguments"

usage() {
    cat << EOF
Usage: $0 --chrno 1 [OPTIONS]

Required Arguments:
  --chrno INT              Input VCF file (can be gzipped)
  --seqlen INT             Output file prefix (default: phylo)
  --selpos FLOAT           Use IUPAC for heterozygous genotypes
  --num_origins INT        Minimum number of samples to create a fasta file (default: 10)
  --h INT                  Number of threads (default: 4)
  --s INT                  Number of bootstrap replicates (default: 1000)
  --g_sel_start INT        Run multiple sequence alignment
  --r FLOAT                Alignment tool: mafft or muscle (default: mafft)
  --sim_relatedness INT    Minimum minor allele frequency (default: 0.0)
  --g_ne_change_start INT  Maximum missing data proportion (default: 1.0)
  --N0 INT      
  --nsam DIR  
  --u FLOAT              
  --slim_script STR        
  --slim_bin STR   
  --r_script     
  --outdir DIR             Output directory (default: "out")
  --help                  Show this help message

Examples:
  # Basic usage
  bash $0 --N 50000 --s 0.5 --nsam 2000

  # With filtering and alignment
  bash $0 --N 50000 --s 0.5 --nsam 2000

  # Custom output and threading
  bash $0 --test

EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --chrno)
            CHRNO="$2"
            shift 2
            ;;
        --seqlen)
            SEQLEN="$2"
            shift 2
            ;;
        --selpos)
            SELPOS="$2"
            shift 2
            ;;
        --num_origins)
            NUM_ORIGINS="$2"
            shift 2
            ;;
        --N)
            N="$2"
            shift 2
            ;;
        --h)
            H="$2"
            shift 2
            ;;
        --s)
            S="$2"
            shift 2
            ;;
        --g_sel_start)
            G_SEL_START="$2"
            shift 2
            ;;
        --r)
            R="$2"
            shift 2
            ;;
        --sim_relatedness)
            SIM_RELATEDNESS="$2"
            shift 2
            ;;
        --g_ne_change_start)
            G_NE_CHANGE_START="$2"
            shift 2
            ;;
        --N0)
            N0="$2"
            shift 2
            ;;
        --u)
            U="$2"
            shift 2
            ;;
        --nsam)
            NSAM="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --test)
            TEST=1
            shift
            ;;
        --slim_script)
            SLIM_SCRIPT="$2"
            shift 2
            ;;
        --slim_bin)
            SLIM_BIN="$2"
            shift 2
            ;;
        --r_script)
            R_SCRIPT="$2"
            shift 2
            ;;
        --cleanup)
            CLEANUP=1
            shift
            ;;
        *)
            log_error "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

################################################################################
# Environment Setup
################################################################################
log_step "Setting Up Environment"

# Load required modules (adjust for your HPC environment)
if command -v module &> /dev/null; then
    log_info "Loading environment modules..."
    # Uncomment and modify these lines based on your cluster
    # module purge
    # module load r/4.3.3
    # module load python/3.12
    # module load slim/4.0
fi

# Activate conda environment
if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
    log_info "Conda environment already active: $CONDA_DEFAULT_ENV"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    log_info "Attempting to activate conda environment: slim-msprime"
    conda activate slim-msprime || log_warn "Could not activate slim-msprime environment"
else
    log_warn "Conda not found. Ensure Python packages are available."
fi

# Set environment variables for reticulate
if [ -n "${CONDA_PREFIX:-}" ]; then
    export RETICULATE_PYTHON="${CONDA_PREFIX}/bin/python"
    log_info "Set RETICULATE_PYTHON to: $RETICULATE_PYTHON"
fi

################################################################################
# Set genome_set_id from SLURM array task ID
################################################################################
GENOME_SET_ID=${SLURM_ARRAY_TASK_ID:-1}
log_info "Genome Set ID (replicate): $GENOME_SET_ID"

################################################################################
# Create Output Directories
################################################################################
log_step "Creating Output Directories"

# Create replicate-specific output directory
REPLICATE_DIR="${OUTDIR}/replicate_${GENOME_SET_ID}"
mkdir -p "$REPLICATE_DIR" || error_exit "Failed to create output directory: $REPLICATE_DIR"
mkdir -p logs || log_warn "Could not create logs directory"

log_info "Output directory: $REPLICATE_DIR"

################################################################################
# Locate R Script
################################################################################
log_step "Locating R Script"

if [ -z "$R_SCRIPT" ]; then
    # Try to find the R script in common locations
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    POSSIBLE_PATHS=(
        "${SCRIPT_DIR}/run_simulation.R"
        "${SCRIPT_DIR}/simulation.R"
        "${SCRIPT_DIR}/../scripts/run_simulation.R"
        "./run_simulation.R"
    )
    
    for path in "${POSSIBLE_PATHS[@]}"; do
        if [ -f "$path" ]; then
            R_SCRIPT="$path"
            log_info "Found R script: $R_SCRIPT"
            break
        fi
    done
    
    if [ -z "$R_SCRIPT" ]; then
        error_exit "R script not found. Please specify with --r_script"
    fi
else
    if [ ! -f "$R_SCRIPT" ]; then
        error_exit "R script not found: $R_SCRIPT"
    fi
    log_info "Using R script: $R_SCRIPT"
fi

################################################################################
# Validate Required Files
################################################################################
log_step "Validating Required Files and Executables"

# Check SLiM binary
if ! command -v "$SLIM_BIN" &> /dev/null; then
    error_exit "SLiM binary '$SLIM_BIN' not found in PATH"
fi
log_info "SLiM binary found: $(which $SLIM_BIN)"

# Check SLiM script
if [ ! -f "$SLIM_SCRIPT" ]; then
    error_exit "SLiM script not found: $SLIM_SCRIPT"
fi
log_info "SLiM script found: $SLIM_SCRIPT"

# Check R executable
if ! command -v Rscript &> /dev/null; then
    error_exit "Rscript not found in PATH"
fi
log_info "R executable found: $(which Rscript)"

################################################################################
# Print Simulation Parameters
################################################################################
log_step "Simulation Parameters"
log_info "Replicate (genome_set_id): $GENOME_SET_ID"
log_info "Chromosome: $CHRNO"
log_info "Sequence length: $SEQLEN bp"
log_info "Selection position: $SELPOS"
log_info "Number of origins: $NUM_ORIGINS"
log_info "Population size (N): $N"
log_info "Dominance coefficient (h): $H"
log_info "Selection coefficient (s): $S"
log_info "Selection start generation: $G_SEL_START"
log_info "Recombination rate (r): $R"
log_info "Simulate relatedness: $SIM_RELATEDNESS"
log_info "Ne change start generation: $G_NE_CHANGE_START"
log_info "Ancestral population size (N0): $N0"
log_info "Mutation rate (u): $U"
log_info "Sample size (nsam): $NSAM"
log_info "Test mode: $TEST"
log_info "Output directory: $REPLICATE_DIR"

################################################################################
# Construct R Command
################################################################################
log_step "Constructing R Command"

R_CMD="Rscript $R_SCRIPT \
  --chrno $CHRNO \
  --seqlen $SEQLEN \
  --selpos $SELPOS \
  --num_origins $NUM_ORIGINS \
  --N $N \
  --h $H \
  --s $S \
  --g_sel_start $G_SEL_START \
  --r $R \
  --sim_relatedness $SIM_RELATEDNESS \
  --g_ne_change_start $G_NE_CHANGE_START \
  --N0 $N0 \
  --u $U \
  --nsam $NSAM \
  --outdir $REPLICATE_DIR \
  --genome_set_id $GENOME_SET_ID \
  --slim_script $SLIM_SCRIPT \
  --slim_bin $SLIM_BIN"

if [ "$TEST" -eq 1 ]; then
    R_CMD="$R_CMD --test --seqlen 300000 --selpos 100000"
fi

log_info "R command:"
echo "$R_CMD" | sed 's/--/\n  --/g'

################################################################################
# Run Simulation
################################################################################
log_step "Running Simulation (Replicate $GENOME_SET_ID)"

START_TIME=$(date +%s)
log_info "Start time: $(date '+%Y-%m-%d %H:%M:%S')"

# Run R script and capture exit status
if eval "$R_CMD"; then
    EXIT_STATUS=0
    log_success "Simulation completed successfully!"
else
    EXIT_STATUS=$?
    error_exit "Simulation failed with exit status $EXIT_STATUS"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s"

################################################################################
# Verify Output Files
################################################################################
log_step "Verifying Output Files"

EXPECTED_FILES=(
    "${GENOME_SET_ID}_${CHRNO}.restart_count"
    "${GENOME_SET_ID}_${CHRNO}.true_ne"
    "${GENOME_SET_ID}_${CHRNO}.daf"
    "${GENOME_SET_ID}_${CHRNO}.trees"
    "${GENOME_SET_ID}_${CHRNO}.vcf.gz"
    "${GENOME_SET_ID}_${CHRNO}.sfs.tsv"
)

MISSING_FILES=0
for file in "${EXPECTED_FILES[@]}"; do
    if [ -f "$REPLICATE_DIR/$file" ]; then
        SIZE=$(du -h "$REPLICATE_DIR/$file" | cut -f1)
        log_info "✓ $file ($SIZE)"
    else
        log_warn "✗ Missing: $file"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
done

if [ $MISSING_FILES -gt 0 ]; then
    log_warn "$MISSING_FILES expected output files are missing"
else
    log_success "All expected output files generated"
fi

################################################################################
# Cleanup Temporary Files
################################################################################
if [ "$CLEANUP" -eq 1 ]; then
    log_step "Cleaning Up Temporary Files"
    
    # Remove SLiM intermediate files
    rm -f "$REPLICATE_DIR"/tmp_slim_out_*.trees 2>/dev/null || true
    rm -f "$REPLICATE_DIR"/slim_stdout_*.txt 2>/dev/null || true
    rm -f "$REPLICATE_DIR"/slim_stderr_*.txt 2>/dev/null || true
    rm -f "$REPLICATE_DIR"/slim_cmd_*.txt 2>/dev/null || true
    
    log_info "Temporary files removed"
fi

################################################################################
# Generate Summary Statistics
################################################################################
log_step "Generating Summary Statistics"

if [ -f "$REPLICATE_DIR/${GENOME_SET_ID}_${CHRNO}.sfs.tsv" ]; then
    TOTAL_VARIANTS=$(awk 'NR>1 {sum+=$2} END {print sum}' "$REPLICATE_DIR/${GENOME_SET_ID}_${CHRNO}.sfs.tsv")
    log_info "Total variants: $TOTAL_VARIANTS"
fi

if [ -f "$REPLICATE_DIR/${GENOME_SET_ID}_${CHRNO}.restart_count" ]; then
    RESTART_COUNT=$(cat "$REPLICATE_DIR/${GENOME_SET_ID}_${CHRNO}.restart_count")
    log_info "SLiM restart count: $RESTART_COUNT"
fi

################################################################################
# Final Summary
################################################################################
log_step "Simulation Complete"
log_success "Replicate $GENOME_SET_ID finished successfully"
log_info "Output directory: $REPLICATE_DIR"
log_info "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"

# Print resource usage if available
if command -v sacct &> /dev/null && [ -n "${SLURM_JOB_ID:-}" ]; then
    log_info "Resource usage:"
    sacct -j "${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}" --format=JobID,MaxRSS,MaxVMSize,CPUTime,Elapsed
fi

log_step "Job Completed Successfully"

exit 0