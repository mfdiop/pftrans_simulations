#!/bin/bash

#SBATCH --job-name=ibd_infer_migration
#SBATCH --output=logs/ibd_migration_%A_%a.out
#SBATCH --error=logs/ibd_migration_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-432%10

###############################################################################
# Script: codes/malaria_transmission_study/03_inferred_metrics_migration.slurm
# Purpose: Run IBD inference methods (IBS, hmmIBD, isoRelate) for migration scenarios
# Author: Mouhamadou Fadel DIOP
# Date: $(date +%Y-%m-%d)
#
# Description:
#   Computes inferred IBD using multiple methods and compares with true IBD
#   from tskibd. Processes all migration scenarios and chromosomes.
#
# Input structure:
#   migration_sims/<scenario_folder>/rep_001/slim_output_<scenario>_1_processed.{trees,vcf.gz}
#   results/malaria_transmission/<scenario_folder>/chr{1,2,3}/  (true IBD from tskibd)
#
# Output structure:
#   results/malaria_transmission/<scenario_folder>/inferred/
#
# Usage:
#   sbatch codes/malaria_transmission_study/03_inferred_metrics_migration.slurm
###############################################################################

set -euo pipefail

###############################################################################
# CONFIGURATION
###############################################################################

# Directories
BASE_SIM_DIR="migration_sims"                            # Simulation output directory
BASE_TRUE_IBD_DIR="results/malaria_transmission"         # True IBD from tskibd
BASE_OUTPUT_DIR="results/malaria_transmission"           # Output directory (same as true IBD)

LOG_DIR="logs"

# Create directories
mkdir -p "${LOG_DIR}"

# R script location
R_SCRIPT="codes/malaria_transmission_study/03_inferred_metrics.R"

# Check if R script exists
if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    echo "Please update R_SCRIPT path or create the R script"
    exit 1
fi

# Methods to run (can be overridden via environment variables)
RUN_IBS=${RUN_IBS:-true}
RUN_HMM=${RUN_HMM:-true}
RUN_ISORELATE=${RUN_ISORELATE:-true}

# Chromosomes to process
#declare -a CHROMOSOMES=(1 2 3)
CHROMOSOMES=1

###############################################################################
# FUNCTIONS
###############################################################################

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Parse folder name to extract parameters
parse_folder_name() {
    local folder_name=$1
    
    # Pattern: MIG{scenario_id:04d}_{rec_name}_{bottle_name}_npop{n_pops}_m{mig_name}
    if [[ "${folder_name}" =~ ^MIG([0-9]{4})_([^_]+)_([^_]+)_npop([0-9]+)_m([^_]+)$ ]]; then
        echo "${BASH_REMATCH[1]}"  # scenario_id
        echo "${BASH_REMATCH[2]}"  # rec_name
        echo "${BASH_REMATCH[3]}"  # bottle_name
        echo "${BASH_REMATCH[4]}"  # n_pops
        echo "${BASH_REMATCH[5]}"  # mig_name
    else
        echo "ERROR: Folder name '${folder_name}' doesn't match expected pattern" >&2
        return 1
    fi
}

# Get all scenario folders
get_scenario_folders() {
    find "${BASE_SIM_DIR}" -maxdepth 1 -type d -name "MIG[0-9][0-9][0-9][0-9]_*" | sort | xargs -I {} basename {}
}

# Find tree and VCF files for a scenario
find_simulation_files() {
    local scenario_folder=$1
    local chr_num=$2
    
    # Base directory for this scenario
    local base_dir="${BASE_SIM_DIR}/${scenario_folder}/rep_001"
    
    # Look for files with pattern: slim_output_{scenario}_{chr}_processed.{trees,vcf.gz}
    # The pattern is: slim_output_MIG0002_low_tight_npop2_mlow_1_processed.{trees,vcf.gz}
    
    # Construct the base filename without extension
    local base_file="${base_dir}/slim_output_${scenario_folder}_${chr_num}_processed"
    
    local tree_file="${base_file}.trees"
    local vcf_file="${base_file}.vcf.gz"
    
    # Alternative patterns to try if not found
    if [[ ! -f "${tree_file}" ]] || [[ ! -f "${vcf_file}" ]]; then
        log_message "WARNING: Files not found with standard pattern, searching alternatives..."
        
        # Try to find any matching files
        local found_tree=$(find "${base_dir}" -name "*${chr_num}*.trees" -type f | head -1)
        local found_vcf=$(find "${base_dir}" -name "*${chr_num}*.vcf.gz" -type f | head -1)
        
        [[ -n "${found_tree}" ]] && tree_file="${found_tree}"
        [[ -n "${found_vcf}" ]] && vcf_file="${found_vcf}"
    fi
    
    echo "${tree_file}"
    echo "${vcf_file}"
}

# Run IBD analysis for a single scenario and chromosome
run_ibd_analysis() {
    local scenario_folder=$1
    local chr_num=$2
    
    # Parse folder name for logging
    local params
    params=$(parse_folder_name "${scenario_folder}")
    if [[ $? -ne 0 ]]; then
        log_message "ERROR: Failed to parse folder name: ${scenario_folder}"
        return 1
    fi
    
    # Read parsed parameters
    read -r scenario_id rec_name bottle_name n_pops mig_name <<< "${params}"
    
    # Find simulation files
    local files
    files=$(find_simulation_files "${scenario_folder}" "${chr_num}")
    
    local tree_file=$(echo "${files}" | sed -n '1p')
    local vcf_file=$(echo "${files}" | sed -n '2p')
    
    # True IBD directory (from previous tskibd script)
    local true_ibd_dir="${BASE_TRUE_IBD_DIR}/${scenario_folder}/chr${chr_num}"
    local true_ibd_file="${true_ibd_dir}/chr${chr_num}.ibd"  # Default tskibd output name
    
    # Alternative true IBD file patterns
    if [[ ! -f "${true_ibd_file}" ]]; then
        local alt_ibd=$(find "${true_ibd_dir}" -name "*.ibd" -type f | head -1)
        [[ -n "${alt_ibd}" ]] && true_ibd_file="${alt_ibd}"
    fi
    
    # Output directory
    local output_dir="${BASE_OUTPUT_DIR}/${scenario_folder}/inferred/chr${chr_num}"
    
    # Check if input files exist
    local missing_files=()
    [[ ! -f "${tree_file}" ]] && missing_files+=("${tree_file}")
    [[ ! -f "${vcf_file}" ]] && missing_files+=("${vcf_file}")
    [[ ! -f "${true_ibd_file}" ]] && missing_files+=("${true_ibd_file}")
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        log_message "ERROR: Missing input files for ${scenario_folder} chr${chr_num}:"
        for file in "${missing_files[@]}"; do
            log_message "  ✗ ${file}"
        done
        
        # List available files for debugging
        log_message "Available files in ${BASE_SIM_DIR}/${scenario_folder}/rep_001/:"
        ls -la "${BASE_SIM_DIR}/${scenario_folder}/rep_001/" 2>/dev/null || log_message "  (directory not found)"
        
        log_message "Available files in ${true_ibd_dir}:"
        ls -la "${true_ibd_dir}" 2>/dev/null || log_message "  (directory not found)"
        
        return 1
    fi
    
    # Create output directory
    mkdir -p "${output_dir}"
    
    # Log parameters
    log_message "=========================================================================="
    log_message "Processing: ${scenario_folder} chr${chr_num}"
    log_message "=========================================================================="
    log_message "Scenario details:"
    log_message "  ID: ${scenario_id}"
    log_message "  Rec rate: ${rec_name}"
    log_message "  Bottleneck: ${bottle_name}"
    log_message "  Populations: ${n_pops}"
    log_message "  Migration: ${mig_name}"
    log_message "Input files:"
    log_message "  ✓ Tree: ${tree_file}"
    log_message "  ✓ VCF: ${vcf_file}"
    log_message "  ✓ True IBD: ${true_ibd_file}"
    log_message "Parameters:"
    log_message "  Chromosome: ${chr_num}"
    log_message "  Output: ${output_dir}"
    log_message "Methods:"
    log_message "  IBS: ${RUN_IBS}"
    log_message "  HMM-IBD: ${RUN_HMM}"
    log_message "  isoRelate: ${RUN_ISORELATE}"
    
    # Build R command
    local r_cmd="Rscript ${R_SCRIPT}"
    r_cmd+=" --tree '${tree_file}'"
    r_cmd+=" --vcf '${vcf_file}'"
    r_cmd+=" --true_ibd '${true_ibd_file}'"
    r_cmd+=" --outdir '${output_dir}'"
    r_cmd+=" --scenario '${scenario_folder}'"
    r_cmd+=" --chromosome ${chr_num}"
    
    # Add method flags
    [[ "${RUN_IBS}" == "true" ]] && r_cmd+=" --ibs"
    [[ "${RUN_HMM}" == "true" ]] && r_cmd+=" --hmm"
    [[ "${RUN_ISORELATE}" == "true" ]] && r_cmd+=" --isorelate"
    
    # Run R script
    log_message "Running R analysis..."
    log_message "Command: ${r_cmd}"
    log_message ""
    
    # Create a wrapper script for better error handling
    local wrapper_script="${output_dir}/run_analysis.sh"
    cat > "${wrapper_script}" << EOF
#!/bin/bash
set -euo pipefail
cd "$(pwd)"
${r_cmd}
EOF
    
    chmod +x "${wrapper_script}"
    
    if "${wrapper_script}" >> "${output_dir}/analysis.log" 2>&1; then
        log_message "✓ SUCCESS: Analysis completed"
        
        # Check output files
        local output_files=()
        [[ -f "${output_dir}/true_ibd_summary.tsv" ]] && output_files+=("true_ibd_summary.tsv")
        [[ -f "${output_dir}/inferred_ibs.tsv" ]] && output_files+=("inferred_ibs.tsv")
        [[ -f "${output_dir}/inferred_ibd_hmm.tsv" ]] && output_files+=("inferred_ibd_hmm.tsv")
        [[ -f "${output_dir}/inferred_ibd_iso.tsv" ]] && output_files+=("inferred_ibd_iso.tsv")
        
        log_message "  Generated ${#output_files[@]} output file(s):"
        for file in "${output_files[@]}"; do
            log_message "    ✓ ${file}"
        done
        log_message "  Log: ${output_dir}/analysis.log"
    else
        log_message "✗ ERROR: Analysis failed"
        log_message "  Check log: ${output_dir}/analysis.log"
        
        # Show last few lines of log for debugging
        if [[ -f "${output_dir}/analysis.log" ]]; then
            log_message "  Last 20 lines of log:"
            tail -20 "${output_dir}/analysis.log" | while IFS= read -r line; do
                log_message "    ${line}"
            done
        fi
        
        return 1
    fi
    
    log_message "=========================================="
    return 0
}

###############################################################################
# MAIN EXECUTION
###############################################################################

log_message "=========================================="
log_message "IBD Inference Pipeline - Migration Scenarios"
log_message "=========================================="

log_message "Configuration:"
log_message "  R script: ${R_SCRIPT}"
log_message "  Base simulation dir: ${BASE_SIM_DIR}"
log_message "  True IBD dir: ${BASE_TRUE_IBD_DIR}"
log_message "  Output dir: ${BASE_OUTPUT_DIR}"
log_message "=========================================="

# Get all scenario folders
scenario_folders=($(get_scenario_folders))
N_SCENARIOS=${#scenario_folders[@]}
N_CHRS=${#CHROMOSOMES[@]}
TOTAL_JOBS=$((N_SCENARIOS * N_CHRS))

log_message "Found ${N_SCENARIOS} scenario folders"
log_message "Processing ${N_CHRS} chromosomes per scenario"
log_message "Total jobs: ${TOTAL_JOBS}"

# Verify expected count
if [[ ${N_SCENARIOS} -ne 432 ]]; then
    log_message "WARNING: Expected 432 scenarios, found ${N_SCENARIOS}"
    log_message "  This might be OK if some simulations are still running or failed"
fi

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # SLURM array mode
    TASK_ID=${SLURM_ARRAY_TASK_ID}
    
    if [[ ${TASK_ID} -lt 1 || ${TASK_ID} -gt ${TOTAL_JOBS} ]]; then
        log_message "ERROR: Invalid task ID ${TASK_ID} (valid: 1-${TOTAL_JOBS})"
        exit 1
    fi
    
    # Calculate scenario and chromosome indices
    SCENARIO_IDX=$(( (TASK_ID - 1) / N_CHRS ))
    CHR_IDX=$(( (TASK_ID - 1) % N_CHRS ))
    
    # Get actual values
    SCENARIO_FOLDER="${scenario_folders[$SCENARIO_IDX]}"
    CHR_NUM="${CHROMOSOMES[$CHR_IDX]}"
    
    log_message "SLURM array job: Task ${TASK_ID}/${TOTAL_JOBS}"
    log_message "  Scenario: ${SCENARIO_FOLDER} (index: ${SCENARIO_IDX})"
    log_message "  Chromosome: ${CHR_NUM}"
    log_message ""
    
    # Run analysis
    run_ibd_analysis "${SCENARIO_FOLDER}" "${CHR_NUM}"
    
else
    # Sequential mode (for testing)
    log_message "Running in sequential mode (testing first scenario)"
    
    if [[ ${N_SCENARIOS} -eq 0 ]]; then
        log_message "ERROR: No scenario folders found in ${BASE_SIM_DIR}"
        exit 1
    fi
    
    # Test with first scenario and first chromosome
    TEST_SCENARIO="${scenario_folders[0]}"
    TEST_CHR="${CHROMOSOMES[0]}"
    
    log_message "Testing with: ${TEST_SCENARIO} chr${TEST_CHR}"
    
    # Run with all methods
    RUN_IBS=true RUN_HMM=true RUN_ISORELATE=true \
        run_ibd_analysis "${TEST_SCENARIO}" "${TEST_CHR}"
    
    exit_code=$?
    if [[ ${exit_code} -eq 0 ]]; then
        log_message "✓ Test successful"
        log_message "To run all jobs, submit with: sbatch --array=1-${TOTAL_JOBS} $0"
        log_message ""
        log_message "To run specific methods on all jobs, use:"
        log_message "  RUN_IBS=true RUN_HMM=false RUN_ISORELATE=false sbatch --array=1-${TOTAL_JOBS} $0"
    else
        log_message "✗ Test failed, check configuration and file paths"
        log_message "Common issues:"
        log_message "  1. R script not found at ${R_SCRIPT}"
        log_message "  2. Simulation files not in expected location"
        log_message "  3. True IBD files not generated yet (run tskibd first)"
        log_message "  4. File naming patterns don't match"
    fi
fi

log_message "=========================================="
log_message "✓ IBD inference pipeline completed"
log_message "=========================================="

exit 0