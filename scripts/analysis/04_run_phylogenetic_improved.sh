#!/bin/bash

#SBATCH --job-name=phylo_batch
#SBATCH --output=phylo_logs/phylo_%A_%a.out
#SBATCH --error=phylo_logs/phylo_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-60%5

###############################################################################
# Script: 04_run_phylogenetic_batch.slurm
# Purpose: Run IQ-TREE phylogenetic analysis for each VCF file
# Author: Improved Version
# Date: 2025-11-20
#
# Description:
#   Converts VCF to FASTA and runs IQ-TREE with ModelFinder
#   Each rep/run/rec_rate/chr gets its own phylogenetic tree
#
# Usage:
#   sbatch codes/multiple_runs/04_run_phylogenetic_batch.slurm
#   sbatch codes/multiple_runs/04_run_phylogenetic_batch.slurm --bootstrap 2000 --min-maf 0.05
#
# Command-line options (override defaults):
#   --input-dir DIR        : Input directory (default: results/multiple_runs)
#   --output-dir DIR       : Output directory (default: results/multiple_runs/phylo_results)
#   --bootstrap N          : Bootstrap replicates (default: 1000)
#   --threads N            : Number of threads (default: ${SLURM_CPUS_PER_TASK})
#   --model-set SET        : Model set to test (default: DNA)
#   --min-maf FLOAT        : Minimum MAF for filtering (default: 0.0)
#   --max-missing FLOAT    : Maximum missing data (default: 1.0)
#   --use-iupac           : Use IUPAC codes for heterozygotes
#   --min-samples N        : Minimum non-missing samples per site (default: 10)
#   --run-alignment        : Run MSA before tree building
#   --alignment-tool TOOL  : Alignment tool: mafft or muscle (default: mafft)
###############################################################################

set -euo pipefail

# Detect if output is a terminal (for color support)
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

log_info() { echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"; }
log_step() { 
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$*${NC}" 
    echo -e "${BLUE}========================================${NC}"
}

###############################################################################
# DEFAULT CONFIGURATION
###############################################################################

# Directories
BASE_INPUT_DIR="results/multiple_runs"
BASE_OUTPUT_DIR="results/multiple_runs/phylo_results"
LOG_DIR="phylo_logs"

# Simulation parameters
N_REPLICATES=5
declare -a REC_RATES=("1.0e-09" "1.0e-08" "1.0e-07" "1.0e-06")
declare -a CHROMOSOMES=(1 2 3)

# Phylogenetic analysis parameters
THREADS=${SLURM_CPUS_PER_TASK:-8}
BOOTSTRAP=1000
MODEL_SET="DNA"
MIN_MAF=0.0
MAX_MISSING=1.0
USE_IUPAC="false"
MIN_SAMPLES=10

# Alignment
RUN_ALIGNMENT="false"
ALIGNMENT_TOOL="mafft"

###############################################################################
# PARSE COMMAND LINE ARGUMENTS
###############################################################################

parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --input-dir)
                BASE_INPUT_DIR="$2"
                shift 2
                ;;
            --output-dir)
                BASE_OUTPUT_DIR="$2"
                shift 2
                ;;
            --bootstrap)
                BOOTSTRAP="$2"
                shift 2
                ;;
            --threads)
                THREADS="$2"
                shift 2
                ;;
            --model-set)
                MODEL_SET="$2"
                shift 2
                ;;
            --min-maf)
                MIN_MAF="$2"
                shift 2
                ;;
            --max-missing)
                MAX_MISSING="$2"
                shift 2
                ;;
            --use-iupac)
                USE_IUPAC="true"
                shift
                ;;
            --min-samples)
                MIN_SAMPLES="$2"
                shift 2
                ;;
            --run-alignment)
                RUN_ALIGNMENT="true"
                shift
                ;;
            --alignment-tool)
                ALIGNMENT_TOOL="$2"
                shift 2
                ;;
            --help)
                cat << 'HELPEOF'
Usage: sbatch script.slurm [OPTIONS]

Options:
  --input-dir DIR        Input directory
  --output-dir DIR       Output directory
  --bootstrap N          Bootstrap replicates (default: 1000)
  --threads N            Number of threads
  --model-set SET        Model set (DNA, GTR, etc.)
  --min-maf FLOAT        Minimum minor allele frequency
  --max-missing FLOAT    Maximum missing data proportion
  --use-iupac           Use IUPAC codes for heterozygotes
  --min-samples N        Minimum samples per site
  --run-alignment        Run multiple sequence alignment
  --alignment-tool TOOL  Alignment tool (mafft/muscle)
  --help                 Show this help

Example:
  sbatch script.slurm --bootstrap 2000 --min-maf 0.05
HELPEOF
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                exit 1
                ;;
        esac
    done
}

# Parse arguments
parse_args "$@"

# Create directories
mkdir -p "${BASE_OUTPUT_DIR}" "${LOG_DIR}"

log_step "Configuration"
log_info "Input directory: ${BASE_INPUT_DIR}"
log_info "Output directory: ${BASE_OUTPUT_DIR}"
log_info "Bootstrap: ${BOOTSTRAP}"
log_info "Threads: ${THREADS}"
log_info "Model set: ${MODEL_SET}"
log_info "Min MAF: ${MIN_MAF}"
log_info "Max missing: ${MAX_MISSING}"
log_info "Use IUPAC: ${USE_IUPAC}"
log_info "Min samples: ${MIN_SAMPLES}"
log_info "Run alignment: ${RUN_ALIGNMENT}"

################################################################################
# CHECK DEPENDENCIES
################################################################################

log_step "Checking dependencies"

check_command() {
    if command -v "$1" &> /dev/null; then
        log_info "✓ Found: $1 ($(command -v "$1"))"
        return 0
    else
        log_warn "✗ Not found: $1"
        return 1
    fi
}

MISSING_DEPS=()

# Get full path for iqtree3 if installed locally
if [ -f "iqtree-3.0.1-Linux/bin/iqtree3" ]; then
    iqtree3=$(realpath "iqtree-3.0.1-Linux/bin/iqtree3")
else
    iqtree3=""
fi

# Required tools
if ! check_command "iqtree2" && ! check_command "iqtree"; then
    if [ -n "${iqtree3}" ] && [ -x "${iqtree3}" ]; then
        log_info "✓ Found: iqtree3 (${iqtree3})"
    else
        MISSING_DEPS+=("iqtree/iqtree2")
    fi
fi

check_command "python3" || check_command "python" || MISSING_DEPS+=("python")

# Optional but recommended
check_command "bcftools" || check_command "vcftools" || log_warn "Neither bcftools nor vcftools found (VCF filtering disabled)"

# Check alignment tool if needed
if [[ "${RUN_ALIGNMENT}" == "true" ]]; then
    if [[ "${ALIGNMENT_TOOL}" == "mafft" ]]; then
        check_command "mafft" || MISSING_DEPS+=("mafft")
    elif [[ "${ALIGNMENT_TOOL}" == "muscle" ]]; then
        check_command "muscle" || MISSING_DEPS+=("muscle")
    fi
fi

if [[ ${#MISSING_DEPS[@]} -gt 0 ]]; then
    log_error "Missing required dependencies: ${MISSING_DEPS[*]}"
    exit 1
fi

# Check Python modules
log_info "Checking Python modules..."
if ! python3 -c "import Bio" 2>/dev/null; then
    log_error "BioPython not found. Install with: pip install biopython"
    exit 1
fi
log_info "✓ BioPython found"

# Python script path
PYTHON_SCRIPT="${BASE_OUTPUT_DIR}/vcf_to_fasta_batch.py"

###############################################################################
# FUNCTIONS
###############################################################################

format_rec_rate() {
    local rec_rate=$1
    python3 -c "import math; exp=int(abs(round(math.log10(${rec_rate})))); print(f'1e{exp:02d}')"
}

# Create improved VCF to FASTA conversion script
create_vcf2fasta_script() {
    cat > "${PYTHON_SCRIPT}" << 'EOFPYTHON'
#!/usr/bin/env python3
"""
VCF to FASTA converter for phylogenetic analysis
Handles haploid and diploid genotypes with robust error handling
"""
import sys
import gzip
from collections import defaultdict

def open_file(filename):
    """Open regular or gzipped file"""
    return gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r')

def vcf_to_fasta(vcf_file, output_file, use_iupac=False, min_samples=0):
    """
    Convert VCF to FASTA format
    
    Args:
        vcf_file: Input VCF file (can be .gz)
        output_file: Output FASTA file
        use_iupac: Use IUPAC codes for heterozygotes (default: False, uses first allele)
        min_samples: Minimum non-missing samples required per site
    """
    print(f"Converting {vcf_file} to {output_file}")
    
    # Parse header for sample names
    samples = []
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                break
    
    if not samples:
        raise ValueError("No samples found in VCF header")
    
    print(f"  Found {len(samples)} samples")
    
    # IUPAC ambiguity codes for heterozygotes
    iupac = {
        ('A','G'):'R', ('G','A'):'R', ('C','T'):'Y', ('T','C'):'Y',
        ('G','C'):'S', ('C','G'):'S', ('A','T'):'W', ('T','A'):'W',
        ('G','T'):'K', ('T','G'):'K', ('A','C'):'M', ('C','A'):'M'
    }
    
    # Initialize storage
    sequences = {s: [] for s in samples}
    n_sites = 0
    n_skipped_multiallelic = 0
    n_skipped_indel = 0
    n_skipped_low_coverage = 0
    n_missing_total = 0
    
    # Process variants
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            
            ref, alt = fields[3], fields[4]
            
            # Skip multi-allelic sites
            if ',' in alt:
                n_skipped_multiallelic += 1
                continue
            
            # Skip indels (only keep SNPs)
            if len(ref) != 1 or len(alt) != 1:
                n_skipped_indel += 1
                continue
            
            # Skip invalid bases
            if not (set(ref).issubset('ACGT') and set(alt).issubset('ACGT')):
                continue
            
            genotypes = fields[9:]
            non_missing = 0
            site_bases = []  # Collect all bases for this site first
            
            for sample, gt_field in zip(samples, genotypes):
                gt = gt_field.split(':')[0]
                
                # Parse genotype robustly
                if '/' in gt:
                    alleles = gt.split('/')
                elif '|' in gt:
                    alleles = gt.split('|')
                else:
                    # Handle edge cases: '.', './.', '.|.', or single values
                    if gt in ['.', './.', '.|.']:
                        alleles = ['.', '.']
                    else:
                        # Assume haploid or unusual format - treat as homozygous
                        alleles = [gt, gt]
                
                # Convert genotype to nucleotide
                try:
                    if '.' in alleles:
                        # Missing data
                        site_bases.append('N')
                        n_missing_total += 1
                    elif len(set(alleles)) == 1:
                        # Homozygous
                        idx = int(alleles[0])
                        site_bases.append(ref if idx == 0 else alt)
                        non_missing += 1
                    else:
                        # Heterozygous
                        non_missing += 1
                        if use_iupac:
                            # Use IUPAC ambiguity codes
                            a1 = ref if alleles[0] == '0' else alt
                            a2 = ref if alleles[1] == '0' else alt
                            site_bases.append(iupac.get((a1, a2), 'N'))
                        else:
                            # Pseudo-haploidize: use first allele
                            idx = int(alleles[0])
                            site_bases.append(ref if idx == 0 else alt)
                except (ValueError, IndexError) as e:
                    # Malformed genotype - treat as missing
                    site_bases.append('N')
                    n_missing_total += 1
            
            # Only add site if enough non-missing samples
            if non_missing >= int(min_samples):
                for sample, base in zip(samples, site_bases):
                    sequences[sample].append(base)
                n_sites += 1
            else:
                n_skipped_low_coverage += 1
            
            # Progress indicator
            if n_sites % 10000 == 0 and n_sites > 0:
                print(f"  Processed {n_sites:,} sites", flush=True)
    
    if n_sites == 0:
        raise ValueError("No valid sites found in VCF after filtering")
    
    # Write FASTA
    with open(output_file, 'w') as out:
        for sample in samples:
            seq = ''.join(sequences[sample])
            out.write(f">{sample}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + '\n')
    
    # Report statistics
    print(f"  ✓ Wrote {len(samples)} sequences with {n_sites:,} sites")
    print(f"  Statistics:")
    print(f"    - Sites retained: {n_sites:,}")
    print(f"    - Skipped (multiallelic): {n_skipped_multiallelic:,}")
    print(f"    - Skipped (indels): {n_skipped_indel:,}")
    print(f"    - Skipped (low coverage): {n_skipped_low_coverage:,}")
    print(f"    - Missing genotypes: {n_missing_total:,}")
    
    return n_sites

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python vcf_to_fasta.py <input.vcf> <output.fasta> <use_iupac> <min_samples>")
        print("  use_iupac: 'true' or 'false'")
        print("  min_samples: minimum number of non-missing samples per site")
        sys.exit(1)
    
    try:
        n_sites = vcf_to_fasta(
            sys.argv[1], 
            sys.argv[2], 
            sys.argv[3].lower() == 'true', 
            int(sys.argv[4])
        )
        print(f"\n✓ Conversion successful: {n_sites} sites written")
        sys.exit(0)
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        sys.exit(1)
EOFPYTHON
    chmod +x "${PYTHON_SCRIPT}"
}

# Run phylogenetic analysis for single VCF
run_phylo_analysis() {
    local replicate=$1
    local run_id=$2
    local rec_rate=$3
    local chr_num=$4
    
    local rec_formatted
    rec_formatted=$(format_rec_rate "${rec_rate}")
    
    # File paths
    local VCF_FILE="${BASE_INPUT_DIR}/rep${replicate}/run${run_id}_rec${rec_formatted}_chr${chr_num}.vcf"
    local OUTPUT_DIR="${BASE_OUTPUT_DIR}/rep${replicate}/run${run_id}_rec${rec_formatted}_chr${chr_num}"
    local PREFIX="run${run_id}_rec${rec_formatted}_chr${chr_num}"
    
    # Check input exists
    if [[ ! -f "${VCF_FILE}" ]]; then
        log_error "VCF not found: ${VCF_FILE}"
        return 1
    fi
    
    mkdir -p "${OUTPUT_DIR}"
    
    # Output files
    local FASTA_FILE="${OUTPUT_DIR}/${PREFIX}.fasta"
    local ALIGNED_FASTA="${OUTPUT_DIR}/${PREFIX}_aligned.fasta"
    local FILTERED_VCF="${OUTPUT_DIR}/${PREFIX}_filtered.vcf"
    local LOG_FILE="${OUTPUT_DIR}/${PREFIX}.log"
    
    # Redirect all output to log
    exec > >(tee -a "${LOG_FILE}") 2>&1
    
    log_step "PHYLOGENETIC ANALYSIS PIPELINE"
    log_info "VCF: ${VCF_FILE}"
    log_info "Output: ${OUTPUT_DIR}"
    log_info "Prefix: ${PREFIX}"
    log_info "Parameters:"
    log_info "  - IUPAC codes: ${USE_IUPAC}"
    log_info "  - Min samples: ${MIN_SAMPLES}"
    log_info "  - Threads: ${THREADS}"
    log_info "  - Bootstrap: ${BOOTSTRAP}"
    log_info "  - Run alignment: ${RUN_ALIGNMENT}"
    log_info "  - Min MAF: ${MIN_MAF}"
    log_info "  - Max missing: ${MAX_MISSING}"
    
    ################################################################################
    # STEP 1: FILTER VCF (OPTIONAL)
    ################################################################################
    
    local VCF_TO_USE="${VCF_FILE}"
    
    if (( $(echo "${MIN_MAF} > 0" | bc -l) )) || (( $(echo "${MAX_MISSING} < 1.0" | bc -l) )); then
        log_step "STEP 1: Filtering VCF"
        log_info "Filters: MAF>=${MIN_MAF}, missing<=${MAX_MISSING}"
        
        if command -v bcftools &> /dev/null; then
            # Calculate max allele frequency
            local MAX_AF
            MAX_AF=$(echo "1.0 - ${MIN_MAF}" | bc -l)
            
            bcftools view \
                --min-af "${MIN_MAF}" \
                --max-af "${MAX_AF}" \
                --threads "${THREADS}" \
                "${VCF_FILE}" | \
            bcftools view \
                --max-missing "${MAX_MISSING}" \
                -O v \
                -o "${FILTERED_VCF}"
                
            VCF_TO_USE="${FILTERED_VCF}"
            log_info "✓ Filtered VCF: ${FILTERED_VCF}"
            
        elif command -v vcftools &> /dev/null; then
            vcftools \
                --vcf "${VCF_FILE}" \
                --maf "${MIN_MAF}" \
                --max-missing "${MAX_MISSING}" \
                --recode \
                --recode-INFO-all \
                --out "${OUTPUT_DIR}/${PREFIX}_filtered"
            
            mv "${OUTPUT_DIR}/${PREFIX}_filtered.recode.vcf" "${FILTERED_VCF}"
            VCF_TO_USE="${FILTERED_VCF}"
            log_info "✓ Filtered VCF: ${FILTERED_VCF}"
        else
            log_warn "No filtering tool available, using original VCF"
        fi
    else
        log_info "STEP 1: Skipping VCF filtering (no filters specified)"
    fi
    
    ################################################################################
    # STEP 2: VCF TO FASTA CONVERSION
    ################################################################################
    
    log_step "STEP 2: Converting VCF to FASTA"
    
    if ! python3 "${PYTHON_SCRIPT}" "${VCF_TO_USE}" "${FASTA_FILE}" "${USE_IUPAC}" "${MIN_SAMPLES}"; then
        log_error "VCF to FASTA conversion failed"
        return 1
    fi
    
    if [[ ! -f "${FASTA_FILE}" ]]; then
        log_error "FASTA file not created"
        return 1
    fi
    
    # Check for sufficient sequences
    local n_seqs
    n_seqs=$(grep -c "^>" "${FASTA_FILE}" || echo "0")
    if [[ ${n_seqs} -lt 4 ]]; then
        log_warn "Insufficient sequences (${n_seqs}) for phylogenetic analysis, skipping"
        return 0
    fi
    
    log_info "✓ FASTA created: ${n_seqs} sequences"
    
    ################################################################################
    # STEP 3: MULTIPLE SEQUENCE ALIGNMENT (OPTIONAL)
    ################################################################################
    
    local FASTA_FOR_TREE="${FASTA_FILE}"
    
    if [[ "${RUN_ALIGNMENT}" == "true" ]]; then
        log_step "STEP 3: Running multiple sequence alignment"
        log_info "Tool: ${ALIGNMENT_TOOL}"
        
        if [[ "${ALIGNMENT_TOOL}" == "mafft" ]]; then
            mafft \
                --auto \
                --thread "${THREADS}" \
                "${FASTA_FILE}" > "${ALIGNED_FASTA}"
        elif [[ "${ALIGNMENT_TOOL}" == "muscle" ]]; then
            muscle \
                -in "${FASTA_FILE}" \
                -out "${ALIGNED_FASTA}" \
                -threads "${THREADS}"
        else
            log_error "Unknown alignment tool: ${ALIGNMENT_TOOL}"
            return 1
        fi
        
        if [[ ! -f "${ALIGNED_FASTA}" ]]; then
            log_error "Alignment failed"
            return 1
        fi
        
        log_info "✓ Aligned FASTA: ${ALIGNED_FASTA}"
        FASTA_FOR_TREE="${ALIGNED_FASTA}"
    else
        log_info "STEP 3: Skipping alignment (SNPs are pre-aligned)"
    fi
    
    ################################################################################
    # STEP 4: MODEL SELECTION
    ################################################################################
    
    log_step "STEP 4: Model selection with ModelFinder"
    log_info "Testing model set: ${MODEL_SET}"
    
    # Determine IQ-TREE command
    local IQTREE_CMD
    if command -v iqtree2 &> /dev/null; then
        IQTREE_CMD="iqtree2"
    elif command -v iqtree &> /dev/null; then
        IQTREE_CMD="iqtree"
    elif [ -n "${iqtree3}" ] && [ -x "${iqtree3}" ]; then
        IQTREE_CMD="${iqtree3}"
    else
        log_error "IQ-TREE not found"
        return 1
    fi
    
    log_info "Using: ${IQTREE_CMD}"
    
    # Run ModelFinder
    "${IQTREE_CMD}" \
        -s "${FASTA_FOR_TREE}" \
        -m MFP \
        -mset "${MODEL_SET}" \
        -nt "${THREADS}" \
        -pre "${OUTPUT_DIR}/${PREFIX}_modelfinder" \
        -quiet 2>&1 | tee -a "${LOG_FILE}"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        log_error "ModelFinder failed"
        return 1
    fi
    
    # Extract best model
    local MODEL_LOG="${OUTPUT_DIR}/${PREFIX}_modelfinder.log"
    local BEST_MODEL
    BEST_MODEL=$(grep "Best-fit model:" "${MODEL_LOG}" 2>/dev/null | head -n 1 | awk '{print $3}')
    
    if [[ -z "${BEST_MODEL}" ]]; then
        log_warn "Could not extract best model, defaulting to GTR+G"
        BEST_MODEL="GTR+G"
    fi
    
    log_info "✓ Best model: ${BEST_MODEL}"
    
    ################################################################################
    # STEP 5: BUILD PHYLOGENETIC TREE
    ################################################################################
    
    log_step "STEP 5: Building phylogenetic tree"
    log_info "Model: ${BEST_MODEL}"
    log_info "Bootstrap: ${BOOTSTRAP}"
    
    "${IQTREE_CMD}" \
        -s "${FASTA_FOR_TREE}" \
        -m "${BEST_MODEL}" \
        -bb "${BOOTSTRAP}" \
        -nt "${THREADS}" \
        -pre "${OUTPUT_DIR}/${PREFIX}" \
        -redo \
        -quiet 2>&1 | tee -a "${LOG_FILE}"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        log_error "Tree building failed"
        return 1
    fi
    
    log_info "✓ Tree: ${OUTPUT_DIR}/${PREFIX}.treefile"
    
    ################################################################################
    # STEP 6: SUMMARY REPORT
    ################################################################################
    
    log_step "STEP 6: Generating summary report"
    
    cat > "${OUTPUT_DIR}/${PREFIX}_report.txt" << EOFREPORT
================================================================================
PHYLOGENETIC ANALYSIS REPORT
================================================================================
Generated: $(date)

INPUT
-----
VCF file: ${VCF_FILE}
Replicate: ${replicate}
Run ID: ${run_id}
Recombination rate: ${rec_rate}
Chromosome: ${chr_num}

FILTERING
---------
Min MAF: ${MIN_MAF}
Max missing: ${MAX_MISSING}

ANALYSIS PARAMETERS
-------------------
Sequences: ${n_seqs}
Bootstrap replicates: ${BOOTSTRAP}
Model tested: ${MODEL_SET}
Best model: ${BEST_MODEL}
Threads: ${THREADS}
IUPAC codes: ${USE_IUPAC}
Minimum samples: ${MIN_SAMPLES}

OUTPUT FILES
------------
- FASTA: ${FASTA_FILE}
- Tree: ${OUTPUT_DIR}/${PREFIX}.treefile
- IQ-TREE report: ${OUTPUT_DIR}/${PREFIX}.iqtree
- ModelFinder log: ${MODEL_LOG}
- Full log: ${LOG_FILE}

================================================================================
View tree: ${OUTPUT_DIR}/${PREFIX}.treefile
Full report: ${OUTPUT_DIR}/${PREFIX}.iqtree
================================================================================
EOFREPORT
    
    log_info "✓ Analysis complete"
    
    return 0
}

###############################################################################
# MAIN EXECUTION
###############################################################################

log_step "Phylogenetic Analysis Batch Processing"

# Create VCF to FASTA conversion script
create_vcf2fasta_script
log_info "Created VCF to FASTA converter"

# Calculate total jobs
N_REC_RATES=${#REC_RATES[@]}
N_CHRS=${#CHROMOSOMES[@]}
TOTAL_JOBS=$((N_REPLICATES * N_REC_RATES * N_CHRS))

log_info "Job configuration:"
log_info "  Replicates: ${N_REPLICATES}"
log_info "  Rec rates: ${N_REC_RATES}"
log_info "  Chromosomes: ${N_CHRS}"
log_info "  Total jobs: ${TOTAL_JOBS}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # SLURM array mode
    TASK_ID=${SLURM_ARRAY_TASK_ID}
    
    if [[ ${TASK_ID} -lt 1 || ${TASK_ID} -gt ${TOTAL_JOBS} ]]; then
        log_error "Invalid task ID ${TASK_ID} (valid range: 1-${TOTAL_JOBS})"
        exit 1
    fi
    
    # Calculate which replicate/rec_rate/chr this task corresponds to
    JOBS_PER_REP=$((N_REC_RATES * N_CHRS))
    
    REPLICATE=$(( (TASK_ID - 1) / JOBS_PER_REP + 1 ))
    RUN_ID=${REPLICATE}
    
    POS_IN_REP=$(( (TASK_ID - 1) % JOBS_PER_REP ))
    REC_RATE_IDX=$(( POS_IN_REP / N_CHRS ))
    CHR_IDX=$(( POS_IN_REP % N_CHRS ))
    
    REC_RATE="${REC_RATES[$REC_RATE_IDX]}"
    CHR_NUM="${CHROMOSOMES[$CHR_IDX]}"
    
    log_info "SLURM array job: Task ${TASK_ID}/${TOTAL_JOBS}"
    log_info "  Replicate: ${REPLICATE}"
    log_info "  Run: ${RUN_ID}"
    log_info "  Rec rate: ${REC_RATE}"
    log_info "  Chromosome: ${CHR_NUM}"
    
    # Run analysis
    if run_phylo_analysis "${REPLICATE}" "${RUN_ID}" "${REC_RATE}" "${CHR_NUM}"; then
        log_step "✓ Task ${TASK_ID} completed successfully"
        exit 0
    else
        log_error "Task ${TASK_ID} failed"
        exit 1
    fi
    
else
    # Sequential mode (no SLURM array)
    log_info "Running in sequential mode"
    
    CURRENT_JOB=0
    FAILED_JOBS=0
    SUCCESS_JOBS=0
    
    for replicate in $(seq 1 ${N_REPLICATES}); do
        run_id=${replicate}
        
        log_step "Replicate ${replicate} (run${run_id})"
        
        for rec_rate in "${REC_RATES[@]}"; do
            for chr_num in "${CHROMOSOMES[@]}"; do
                CURRENT_JOB=$((CURRENT_JOB + 1))
                
                log_info ""
                log_info "Job ${CURRENT_JOB}/${TOTAL_JOBS}"
                
                if run_phylo_analysis "${replicate}" "${run_id}" "${rec_rate}" "${chr_num}"; then
                    SUCCESS_JOBS=$((SUCCESS_JOBS + 1))
                else
                    FAILED_JOBS=$((FAILED_JOBS + 1))
                    log_warn "Job ${CURRENT_JOB} failed, continuing..."
                fi
            done
        done
    done
    
    # Final summary
    log_step "Sequential processing complete"
    log_info "Results:"
    log_info "  ✓ Successful: ${SUCCESS_JOBS}/${TOTAL_JOBS}"
    log_info "  ✗ Failed: ${FAILED_JOBS}/${TOTAL_JOBS}"
    
    if [[ ${FAILED_JOBS} -gt 0 ]]; then
        log_error "Some jobs failed"
        exit 1
    else
        log_step "✓ All jobs completed successfully"
        exit 0
    fi
fi