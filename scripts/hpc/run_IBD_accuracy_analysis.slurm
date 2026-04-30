#!/bin/bash
#SBATCH --job-name=ibd_accuracy
#SBATCH --output=logs/accuracy_%j.out
#SBATCH --error=logs/accuracy_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

################################################################################
# Script: run_accuracy_analysis.sh
# Purpose: Run IBD accuracy analysis and generate plots
# Usage: 
#   sbatch run_accuracy_analysis.sh
#   bash run_accuracy_analysis.sh
################################################################################

set -euo pipefail

# Configuration
R_SCRIPT="analyze_ibd_accuracy.R"
INPUT_DIR="ibd_results"
OUTPUT_DIR="ibd_accuracy"
MIN_LENGTH=0
OVERLAP_THRESHOLD=0.5
METHODS="hmm,iso,ibs"

# Create log directory
mkdir -p logs

# Check if R script exists
if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    exit 1
fi

# Check if input directory exists
if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "ERROR: Input directory not found: ${INPUT_DIR}"
    exit 1
fi

echo "=========================================="
echo "IBD Accuracy Analysis"
echo "=========================================="
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Methods: ${METHODS}"
echo "Min segment length: ${MIN_LENGTH} bp"
echo "Overlap threshold: ${OVERLAP_THRESHOLD}"
echo "=========================================="

# Run R script
Rscript "${R_SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --min_length "${MIN_LENGTH}" \
    --overlap_threshold "${OVERLAP_THRESHOLD}" \
    --methods "${METHODS}"

exit_code=$?

if [[ ${exit_code} -eq 0 ]]; then
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "Results saved to: ${OUTPUT_DIR}"
    echo "=========================================="
else
    echo "=========================================="
    echo "ERROR: Analysis failed with exit code ${exit_code}"
    echo "=========================================="
    exit ${exit_code}
fi