#!/bin/bash

#SBATCH --job-name=tskibd_migration
#SBATCH --output=logs/tskibd_%A_%a.out
#SBATCH --error=logs/tskibd_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-432%10  
                 
###############################################################################
# Script: 02_run_tskibd_migration.slurm
# Purpose: Run tskibd on tree sequences from migration simulations
# Author: Mouhamadou Fadel DIOP
# Date: $(date +%Y-%m-%d)
#
# Description:
#   Processes tree sequences from migration simulation scenarios (432 folders)
#   Runs tskibd to extract true IBD segments with parameters scaled by rec rate
#
# Usage:
#   sbatch 02_run_tskibd_migration.slurm
###############################################################################

set -euo pipefail

###############################################################################
# CONFIGURATION
###############################################################################

# Directories
BASE_INPUT_DIR="migration_sims"                    # Base input directory
BASE_OUTPUT_DIR="results/malaria_transmission"     # Base output directory
LOG_DIR="logs"

# Create directories
mkdir -p "${BASE_OUTPUT_DIR}" "${LOG_DIR}"

# MISSING INFORMATION 1: Please provide the mapping from rec_name to actual recombination rate values
# Example format (uncomment and modify):
# declare -A REC_RATE_MAP=(
#     ["low"]="1.0e-09"
#     ["medium"]="1.0e-08"
#     ["high"]="1.0e-07"
#     ["very high"]="1.0e-06"
# )
declare -A REC_RATE_MAP

# Chromosomes to process
#declare -a CHROMOSOMES=(1 2 3)

# tskibd parameters
MIN_CM=2                                # Minimum segment length in cM

# Reference values for P. falciparum
REF_REC_RATE=6.666667e-07       # 1 cM / 15,000 bp
REF_BP_PER_CM=15000
REF_SAMPLING_WINDOW=150

# tskibd executable
TSKIBD_EXEC="tskibd/build/tskibd"

# Check if tskibd exists
if [[ ! -f "${TSKIBD_EXEC}" ]]; then
    echo "ERROR: tskibd executable not found: ${TSKIBD_EXEC}"
    echo "Please compile tskibd or update TSKIBD_EXEC path"
    exit 1
fi

###############################################################################
# FUNCTIONS
###############################################################################

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Extract recombination rate from rec_name
get_rec_rate() {
    local rec_name=$1
    if [[ -v REC_RATE_MAP["${rec_name}"] ]]; then
        echo "${REC_RATE_MAP[${rec_name}]}"
    else
        echo "ERROR: rec_name '${rec_name}' not found in REC_RATE_MAP" >&2
        return 1
    fi
}

# Calculate bp_per_cm from recombination rate
# Formula: bp_per_cm = 0.01 / rec_rate
calculate_bp_per_cm() {
    local rec_rate=$1
    python3 -c "print(int(0.01 / ${rec_rate}))"
}

# Calculate sampling window scaled by rec rate
calculate_sampling_window() {
    local bp_per_cm=$1
    python3 -c "print(int((${bp_per_cm} / ${REF_BP_PER_CM}) * ${REF_SAMPLING_WINDOW}))"
}

# Parse folder name to extract parameters
# Allowed parameter keys
REC_KEYS=("low" "medium" "high" "very_high")
BOTTLENECK_KEYS=("tight" "medium" "loose")
MIG_KEYS=("low" "medium" "high")

in_array() {
    local needle="$1"; shift
    local item
    for item in "$@"; do
        [[ "$item" == "$needle" ]] && return 0
    done
    return 1
}

parse_folder_name() {
    local input="$1"

    # Remove trailing slash (if any)
    input="${input%/}"

    # Extract basename (remove path)
    local folder_name
    folder_name="$(basename "$input")"

    local regex='^MIG([0-9]{4})_([a-z_]+)_([a-z]+)_npop([0-9]+)_m([a-z]+)$'

    if [[ ! "$folder_name" =~ $regex ]]; then
        echo "ERROR: '$folder_name' does not match expected naming scheme" >&2
        return 1
    fi

    local scenario_id="${BASH_REMATCH[1]}"
    local rec_name="${BASH_REMATCH[2]}"
    local bottle_name="${BASH_REMATCH[3]}"
    local n_pops="${BASH_REMATCH[4]}"
    local mig_name="${BASH_REMATCH[5]}"

    # Validate semantics
    if ! in_array "$rec_name" "${REC_KEYS[@]}"; then
        echo "ERROR: Invalid recombination key '$rec_name'" >&2
        return 1
    fi

    if ! in_array "$bottle_name" "${BOTTLENECK_KEYS[@]}"; then
        echo "ERROR: Invalid bottleneck key '$bottle_name'" >&2
        return 1
    fi

    if ! in_array "$mig_name" "${MIG_KEYS[@]}"; then
        echo "ERROR: Invalid migration key '$mig_name'" >&2
        return 1
    fi

    printf '%s %s %s %s %s\n' "$scenario_id" "$rec_name" "$bottle_name" "$n_pops" "$mig_name"

}


# Get all scenario folders
get_scenario_folders() {
    find "${BASE_INPUT_DIR}" -maxdepth 1 -type d -name "MIG[0-9][0-9][0-9][0-9]_*" | sort | xargs -I {} basename {}
}

# Run tskibd for a single scenario and chromosome
run_tskibd() {
    local folder_name=$1
    local chr_num=$2
    
    # Parse folder name
    local params
    params=$(parse_folder_name "${folder_name}")
    if [[ $? -ne 0 ]]; then
        log_message "ERROR: Failed to parse folder name: ${folder_name}"
        return 1
    fi
    
    # Read parsed parameters
    read -r scenario_id rec_name bottle_name n_pops mig_name <<< "${params}"
#    read -r scenario_id rec_name bottle_name n_pops mig_name \
#    <<< "$(parse_folder_name "$folder_name")"

#    mapfile -t params < <(parse_folder_name "$folder_name")

#    scenario_id="${params[0]}"
#    rec_name="${params[1]}"
#    bottle_name="${params[2]}"
#    n_pops="${params[3]}"
#    mig_name="${params[4]}"


    
    # Get recombination rate
    local rec_rate
    rec_rate=$(get_rec_rate "${rec_name}")
    if [[ $? -ne 0 ]]; then
        log_message "ERROR: ${rec_rate}"
        return 1
    fi
    
    # Construct file paths
    local input_file="${BASE_INPUT_DIR}/${folder_name}/rep001/slim_output_${folder_name}_1_processed.trees"
    local output_dir="${BASE_OUTPUT_DIR}/${folder_name}/chr1"
    
    # MISSING INFORMATION 2: Verify the actual filename pattern
    # The input file pattern above assumes: run1_rec{rec_name}_chr{chr_num}_recap.trees
    # Please confirm if this matches your actual filenames
    
    # Check if input file exists
    if [[ ! -f "${input_file}" ]]; then
        log_message "WARNING: Input file not found: ${input_file}"
        log_message "  Trying alternative patterns..."
        
        # Try alternative patterns
        local alt_patterns=(
            "${BASE_INPUT_DIR}/${folder_name}/rep001/slim_output_${folder_name}_1_processed.trees"
            "${BASE_INPUT_DIR}/${folder_name}/rep001/chr${chr_num}_recap.trees"
            "${BASE_INPUT_DIR}/${folder_name}/rep001/chr${chr_num}.trees"
            "${BASE_INPUT_DIR}/${folder_name}/chr${chr_num}_recap.trees"
            "${BASE_INPUT_DIR}/${folder_name}/chr${chr_num}.trees"
        )
        
        for alt_pattern in "${alt_patterns[@]}"; do
            if [[ -f "${alt_pattern}" ]]; then
                input_file="${alt_pattern}"
                log_message "  Found alternative: ${input_file}"
                break
            fi
        done
        
        if [[ ! -f "${input_file}" ]]; then
            log_message "ERROR: No tree sequence file found for ${folder_name} chr${chr_num}"
            return 1
        fi
    fi
    
    # Create output directory
    mkdir -p "${output_dir}"
    
    # Use reference parameters (or calculate based on rec_rate if needed)
    # Uncomment the following lines if you want to scale parameters by rec_rate:
    # local bp_per_cm=$(calculate_bp_per_cm "${rec_rate}")
    # local sampling_window=$(calculate_sampling_window "${bp_per_cm}")
    
    # For consistency, using reference parameters
    local bp_per_cm=${REF_BP_PER_CM}
    local sampling_window=${REF_SAMPLING_WINDOW}
    
    # Get sequence length and validate parameters
    local seq_length=$(python3 -c "import tskit; ts = tskit.load('${input_file}'); print(int(ts.sequence_length))" 2>/dev/null || echo "0")
    
    if [[ ${seq_length} -eq 0 ]]; then
        log_message "ERROR: Could not read sequence length from ${input_file}"
        return 1
    fi
    
    # Log parameters
    log_message "=========================================================================="
    log_message "Processing: ${folder_name} chr${chr_num}"
    log_message "  Scenario ID: ${scenario_id}"
    log_message "  Rec name: ${rec_name} (rate: ${rec_rate})"
    log_message "  Bottleneck: ${bottle_name}"
    log_message "  Populations: ${n_pops}"
    log_message "  Migration: ${mig_name}"
    log_message "  Input: ${input_file}"
    log_message "  Seq length: ${seq_length} bp"
    log_message "  bp_per_cm: ${bp_per_cm}"
    log_message "  sampling_window: ${sampling_window}"
    log_message "  min_cm: ${MIN_CM}"
    log_message "  Output: ${output_dir}"
    
    # Run tskibd
    log_message "Running tskibd..."
    "${TSKIBD_EXEC}" "${chr_num}" "${bp_per_cm}" "${sampling_window}" "${MIN_CM}" \
                     "${input_file}" "${output_dir}" 2>&1 | tee "${output_dir}/tskibd.log"
    
    local exit_code=$?
    
    if [[ ${exit_code} -eq 0 ]]; then
        log_message "✓ SUCCESS: tskibd completed"
        
        # Check if output files were created
        local ibd_files=$(find "${output_dir}" -name "*.ibd" | wc -l)
        if [[ ${ibd_files} -eq 0 ]]; then
            log_message "WARNING: No .ibd files generated (might be no IBD segments)"
        else
            log_message "  Generated ${ibd_files} .ibd file(s)"
        fi
    else
        log_message "✗ ERROR: tskibd failed (exit code: ${exit_code})"
        return ${exit_code}
    fi
    
    log_message "=========================================="
    return 0
}

###############################################################################
# MAIN EXECUTION
###############################################################################

log_message "=========================================="
log_message "tskibd Migration Simulations Processing"
log_message "=========================================="

# Get all scenario folders
scenario_folders=($(get_scenario_folders))
N_SCENARIOS=${#scenario_folders[@]}
N_CHRS=${#CHROMOSOMES[@]}
TOTAL_JOBS=$((N_SCENARIOS * N_CHRS))

log_message "Configuration:"
log_message "  Scenarios: ${N_SCENARIOS}"
log_message "  Chromosomes: ${N_CHRS}"
log_message "  Total jobs: ${TOTAL_JOBS}"
log_message "  Array task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
log_message "=========================================="

# MISSING INFORMATION 3: Check if we have the right number of folders
if [[ ${N_SCENARIOS} -ne 432 ]]; then
    log_message "WARNING: Expected 432 scenarios, found ${N_SCENARIOS}"
    log_message "  This might be OK if some simulations are still running"
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
    FOLDER_NAME="${scenario_folders[$SCENARIO_IDX]}"
    CHR_NUM="${CHROMOSOMES[$CHR_IDX]}"
    
    log_message "SLURM array job: Task ${TASK_ID}/${TOTAL_JOBS}"
    log_message "  Scenario: ${FOLDER_NAME} (index: ${SCENARIO_IDX})"
    log_message "  Chromosome: ${CHR_NUM}"
    
    # Run tskibd
    run_tskibd "${FOLDER_NAME}" "${CHR_NUM}"
    
else
    # Sequential mode (for testing)
    log_message "Running in sequential mode (testing first scenario)"
    
    if [[ ${N_SCENARIOS} -eq 0 ]]; then
        log_message "ERROR: No scenario folders found in ${BASE_INPUT_DIR}"
        exit 1
    fi
    
    # Test with first scenario
    TEST_FOLDER="${scenario_folders[0]}"
    TEST_CHR="${CHROMOSOMES[0]}"
    
    log_message "Testing with: ${TEST_FOLDER} chr${TEST_CHR}"
    run_tskibd "${TEST_FOLDER}" "${TEST_CHR}"
    
    exit_code=$?
    if [[ ${exit_code} -eq 0 ]]; then
        log_message "✓ Test successful"
        log_message "To run all jobs, submit with: sbatch --array=1-${TOTAL_JOBS} $0"
    else
        log_message "✗ Test failed, check configuration"
    fi
fi

log_message "=========================================="
log_message "✓ Processing completed"
log_message "=========================================="

exit 0