#!/bin/bash

###############################################################################
# Script: 05_check_phylo_results.sh
# Purpose: Validate and summarize phylogenetic analysis results
# Usage: bash 05_check_phylo_results.sh [output_dir]
###############################################################################

set -euo pipefail

# Colors
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    RED='\033[0;31m'
    BLUE='\033[0;34m'
    NC='\033[0m'
else
    GREEN=''
    YELLOW=''
    RED=''
    BLUE=''
    NC=''
fi

log_info() { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$*${NC}"
    echo -e "${BLUE}========================================${NC}"
}

# Default directory
PHYLO_DIR="${1:-results/multiple_runs/phylo_results}"

if [[ ! -d "${PHYLO_DIR}" ]]; then
    log_error "Directory not found: ${PHYLO_DIR}"
    exit 1
fi

log_step "Checking Phylogenetic Analysis Results"
log_info "Directory: ${PHYLO_DIR}"

# Initialize counters
TOTAL_EXPECTED=0
TOTAL_FOUND=0
MISSING_TREES=0
MISSING_FASTA=0
FAILED_JOBS=0

# Output file
SUMMARY_FILE="${PHYLO_DIR}/analysis_summary.txt"

# Find all expected outputs
log_info "Scanning for output files..."

# Count replicates, rec_rates, chromosomes
N_REPLICATES=$(find "${PHYLO_DIR}" -maxdepth 1 -type d -name "rep*" | wc -l)
log_info "Found ${N_REPLICATES} replicate directories"

# Create detailed report
{
    echo "=================================================================================="
    echo "PHYLOGENETIC ANALYSIS SUMMARY"
    echo "=================================================================================="
    echo "Date: $(date)"
    echo "Directory: ${PHYLO_DIR}"
    echo ""
    echo "OVERVIEW"
    echo "--------"
    
    # Check each replicate
    for rep_dir in "${PHYLO_DIR}"/rep*/; do
        if [[ ! -d "${rep_dir}" ]]; then
            continue
        fi
        
        rep_name=$(basename "${rep_dir}")
        echo ""
        echo "=== ${rep_name} ==="
        
        # Count subdirectories (each is a job)
        n_jobs=$(find "${rep_dir}" -maxdepth 1 -mindepth 1 -type d | wc -l)
        echo "  Expected jobs: ${n_jobs}"
        
        # Check for required files in each job
        complete=0
        incomplete=0
        
        for job_dir in "${rep_dir}"/*/; do
            if [[ ! -d "${job_dir}" ]]; then
                continue
            fi
            
            job_name=$(basename "${job_dir}")
            TOTAL_EXPECTED=$((TOTAL_EXPECTED + 1))
            
            # Check for required files
            has_tree=false
            has_fasta=false
            has_iqtree=false
            
            if ls "${job_dir}"/*.treefile &>/dev/null; then
                has_tree=true
            fi
            
            if ls "${job_dir}"/*.fasta &>/dev/null; then
                has_fasta=true
            fi
            
            if ls "${job_dir}"/*.iqtree &>/dev/null; then
                has_iqtree=true
            fi
            
            if ${has_tree} && ${has_fasta} && ${has_iqtree}; then
                complete=$((complete + 1))
                TOTAL_FOUND=$((TOTAL_FOUND + 1))
            else
                incomplete=$((incomplete + 1))
                echo "    INCOMPLETE: ${job_name}"
                ${has_tree} || echo "      - Missing: treefile"
                ${has_fasta} || echo "      - Missing: FASTA"
                ${has_iqtree} || echo "      - Missing: IQ-TREE report"
                
                ${has_tree} || MISSING_TREES=$((MISSING_TREES + 1))
                ${has_fasta} || MISSING_FASTA=$((MISSING_FASTA + 1))
            fi
        done
        
        echo "  Complete: ${complete}"
        echo "  Incomplete: ${incomplete}"
    done
    
    echo ""
    echo "=================================================================================="
    echo "SUMMARY STATISTICS"
    echo "=================================================================================="
    echo "Total expected jobs: ${TOTAL_EXPECTED}"
    echo "Complete jobs: ${TOTAL_FOUND}"
    echo "Incomplete jobs: $((TOTAL_EXPECTED - TOTAL_FOUND))"
    echo "Success rate: $(echo "scale=2; 100 * ${TOTAL_FOUND} / ${TOTAL_EXPECTED}" | bc)%"
    echo ""
    echo "Missing files:"
    echo "  - Tree files: ${MISSING_TREES}"
    echo "  - FASTA files: ${MISSING_FASTA}"
    echo ""
    
    # Check log files for errors
    echo "ERROR ANALYSIS"
    echo "--------------"
    if find "${PHYLO_DIR}" -name "*.log" -exec grep -l "ERROR\|Error\|error" {} \; | head -n 10 | grep -q .; then
        echo "Jobs with errors found in logs:"
        find "${PHYLO_DIR}" -name "*.log" -exec grep -l "ERROR\|Error\|error" {} \; | head -n 10 | while read -r logfile; do
            echo "  - ${logfile}"
        done
        echo ""
        echo "Check individual log files for details"
    else
        echo "No errors found in log files (checked first 10)"
    fi
    
    echo ""
    echo "=================================================================================="
    echo "FILE SIZE STATISTICS"
    echo "=================================================================================="
    
    # Tree file sizes
    if find "${PHYLO_DIR}" -name "*.treefile" | head -n 1 | grep -q .; then
        echo "Tree files:"
        find "${PHYLO_DIR}" -name "*.treefile" -exec du -h {} \; | \
            awk '{sum+=$1; count++} END {print "  Total: " count " files"}'
        
        # Average size
        echo -n "  Average size: "
        find "${PHYLO_DIR}" -name "*.treefile" -exec du -b {} \; | \
            awk '{sum+=$1; count++} END {printf "%.2f KB\n", sum/count/1024}'
    fi
    
    # FASTA file sizes
    if find "${PHYLO_DIR}" -name "*.fasta" | head -n 1 | grep -q .; then
        echo "FASTA files:"
        find "${PHYLO_DIR}" -name "*.fasta" -exec du -h {} \; | \
            awk '{sum+=$1; count++} END {print "  Total: " count " files"}'
        
        echo -n "  Average size: "
        find "${PHYLO_DIR}" -name "*.fasta" -exec du -b {} \; | \
            awk '{sum+=$1; count++} END {printf "%.2f KB\n", sum/count/1024}'
    fi
    
    echo ""
    echo "Total disk usage:"
    du -sh "${PHYLO_DIR}"
    
    echo ""
    echo "=================================================================================="
    echo "RECOMMENDATIONS"
    echo "=================================================================================="
    
    if [[ $((TOTAL_EXPECTED - TOTAL_FOUND)) -gt 0 ]]; then
        echo "✗ Some jobs incomplete or failed"
        echo "  Action: Check SLURM logs in phylo_logs/ for error messages"
        echo "  Resubmit failed jobs with: sbatch --array=X,Y,Z script.slurm"
    else
        echo "✓ All jobs completed successfully"
    fi
    
    echo ""
    echo "Next steps:"
    echo "  1. View trees with FigTree, iTOL, or R (ape, phytools)"
    echo "  2. Compare tree topologies across replicates"
    echo "  3. Analyze bootstrap support values"
    echo "  4. Extract phylogenetic distances for transmission inference"
    
    echo ""
    echo "=================================================================================="
    
} | tee "${SUMMARY_FILE}"

log_info ""
log_info "Summary saved to: ${SUMMARY_FILE}"

# Exit with appropriate code
if [[ ${TOTAL_FOUND} -eq ${TOTAL_EXPECTED} ]]; then
    log_step "✓ All phylogenetic analyses completed successfully"
    exit 0
else
    log_warn "Some analyses incomplete or failed"
    exit 1
fi
