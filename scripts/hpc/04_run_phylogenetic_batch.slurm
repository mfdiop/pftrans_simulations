#!/bin/bash

#SBATCH --job-name=phylo_batch
#SBATCH --output=phylo_logs/phylo_%A_%a.out
#SBATCH --error=phylo_logs/phylo_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-60%5                    # 5 reps × 4 rec_rates × 3 chrs = 60 jobs, max 5 concurrent

###############################################################################
# Script: 04_run_phylogenetic_batch.slurm
# Purpose: Run IQ-TREE phylogenetic analysis for each VCF file
# Author: [Your Name]
# Date: 2025-11-05
#
# Description:
#   Converts each VCF to FASTA and runs IQ-TREE with ModelFinder
#   Each rep/run/rec_rate/chr gets its own phylogenetic tree
#
# Input structure:
#   results/vary_rec_rate/rep1/run1_rec1e09_chr1.vcf
#   results/vary_rec_rate/rep1/run1_rec1e08_chr1.vcf
#   ...
#   results/vary_rec_rate/rep5/run5_rec1e06_chr3.vcf
#
# Output structure:
#   phylo_results/rep1/
#     ├── run1_rec1e09_chr1/
#     │   ├── phylo.fasta
#     │   ├── phylo.treefile
#     │   ├── phylo.iqtree
#     │   └── phylo_report.txt
#     ├── run1_rec1e08_chr1/
#     └── ...
#
# Usage:
#   sbatch codes/multiple_runs/04_run_phylogenetic_batch.slurm
###############################################################################

set -euo pipefail

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'  # No Color

log_info() { echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*";}
log_step() {
	echo ""
	echo -e "${BLUE}========================================${NC}"
	echo -e "${BLUE}$*${NC}"
	echo -e "${BLUE}========================================${NC}"
}

###############################################################################
# CONFIGURATION
###############################################################################

# Directories
BASE_INPUT_DIR="results/multiple_runs"                    # Input VCFs (rep1-5)
BASE_OUTPUT_DIR="results/multiple_runs/phylo_results"         # Output directory

LOG_DIR="logs"

mkdir -p "${BASE_OUTPUT_DIR}" "${LOG_DIR}"

# Simulation parameters
N_REPLICATES=5
declare -a REC_RATES=("1.0e-09" "1.0e-08" "1.0e-07" "1.0e-06")
declare -a CHROMOSOMES=(1 2 3)

# Phylogenetic analysis parameters
THREADS=${SLURM_CPUS_PER_TASK:-8}
BOOTSTRAP=1000
MODEL_SET="DNA"                         # DNA models for SNP data
MIN_MAF=0.0
MAX_MISSING=1.0
USE_IUPAC=false
MIN_SAMPLES=10

# Alignment
RUN_ALIGNMENT=false
ALIGNMENT_TOOL="mafft"  # Options: mafft, muscle, none

# IQ-TREE executable
IQTREE_CMD=$(command -v iqtree2 || command -v iqtree || echo "iqtree-3.0.1-Linux/bin/iqtree3")

if [[ ! -x "${IQTREE_CMD}" ]] && [[ ! -f "${IQTREE_CMD}" ]]; then
    log_error "IQ-TREE not found. Checked: iqtree2, iqtree, ${IQTREE_CMD}"
    exit 1
fi

log_info "Using IQ-TREE: ${IQTREE_CMD}"

################################################################################
# CHECK DEPENDENCIES
################################################################################

log_step "STEP 0: Checking dependencies"

check_command() {
    if command -v "$1" &> /dev/null; then
        log_info "✓ Found: $1 ($(command -v $1))"
        return 0
    else
        log_warn "✗ Not found: $1"
        return 1
    fi
}

MISSING_DEPS=()

# Get full/absolute path of .....
iqtree3=$(realpath iqtree-3.0.1-Linux/bin/iqtree3)   # readlink -f program

# Required tools
check_command "iqtree2" || check_command "iqtree" || check_command "${iqtree3}" || MISSING_DEPS+=("iqtree/iqtree2")
check_command "python3" || check_command "python" || MISSING_DEPS+=("python")

# Optional tools
check_command "bcftools" || check_command "vcftools" || log_warn "Neither bcftools nor vcftools found (needed for filtering)"

# Check alignment tool
if [[ "${RUN_ALIGNMENT}" == "true" ]]; then
    if [[ "${ALIGNMENT_TOOL}" == "mafft" ]]; then
        check_command "mafft" || MISSING_DEPS+=("mafft")
    elif [[ "${ALIGNMENT_TOOL}" == "muscle" ]]; then
        check_command "muscle" || MISSING_DEPS+=("muscle")
    fi
fi

if [[ ${#MISSING_DEPS[@]} -gt 0 ]]; then
    log_error "Missing required dependencies: ${MISSING_DEPS[*]}"
    log_error "Please install missing tools and try again"
    exit 1
fi


# Check Python modules
log_info "Checking Python modules..."
python3 << 'EOF'
import sys
missing = []
try:
    import Bio
    print("  ✓ BioPython found")
except ImportError:
    missing.append("biopython")
    print("  ✗ BioPython not found")

if missing:
    print(f"\nMissing Python modules: {', '.join(missing)}")
    print("Install with: pip install " + " ".join(missing))
    sys.exit(1)
EOF

if [[ $? -ne 0 ]]; then
    log_error "Python dependency check failed"
    exit 1
fi

# Python script for VCF to FASTA conversion
PYTHON_SCRIPT="${BASE_OUTPUT_DIR}/vcf_to_fasta_batch.py"

###############################################################################
# FUNCTIONS
###############################################################################

format_rec_rate() {
    local rec_rate=$1
    python3 -c "import math; exp=int(abs(round(math.log10(${rec_rate})))); print(f'1e{exp:02d}')"
}

# Create VCF to FASTA conversion script
create_vcf2fasta_script() {
    cat > "${PYTHON_SCRIPT}" << 'EOFPYTHON'
#!/usr/bin/env python3
import sys, gzip
from collections import defaultdict

def open_file(f):
    return gzip.open(f, 'rt') if f.endswith('.gz') else open(f, 'r')
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

    # Parse header
    samples = []
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                break

    if not samples:
        raise ValueError("No samples found in VCF")

    print(f"  Found {len(samples)} samples")

    # IUPAC ambiguity codes for heterozygotes
    iupac = {
        ('A','G'):'R', ('G','A'):'R', ('C','T'):'Y', ('T','C'):'Y',
        ('G','C'):'S', ('C','G'):'S', ('A','T'):'W', ('T','A'):'W',
        ('G','T'):'K', ('T','G'):'K', ('A','C'):'M', ('C','A'):'M'
    }

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

            # Skip indels
            if len(ref) != 1 or len(alt) != 1:
                n_skipped_indel += 1
                continue

            # Skip invalid bases
            if not (set(ref).issubset('ACGT') and set(alt).issubset('ACGT')):
                continue

            genotypes = fields[9:]
            non_missing = 0
            site_bases = []

            for sample, gt_field in zip(samples, genotypes):
                gt = gt_field.split(':')[0]

                # Parse genotype more robustly
                if '/' in gt:
                    alleles = gt.split('/')
                elif '|' in gt:
                    alleles = gt.split('|')
                else:
                    # Handle edge cases
                    if gt in ['.', './.', '.|.']:
                        alleles = ['.', '.']
                    else:
                        # Assume haploid or unusual format
                        alleles = [gt, gt]

                # Convert to nucleotide
                try:
                    if '.' in alleles:
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

            if n_sites % 10000 == 0 and n_sites > 0:
                print(f"  Processed {n_sites:,} sites")

    if n_sites == 0:
        raise ValueError("No valid sites found in VCF")

    # Write FASTA
    with open(output_file, 'w') as out:
        for sample in samples:
            seq = ''.join(sequences[sample])
            out.write(f">{sample}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + '\n')

    # Report statistics
    print(f"  ✓ Wrote {len(samples)} sequences, {n_sites:,} sites")
    print(f"  Statistics:")
    print(f"    - Sites kept: {n_sites:,}")
    print(f"    - Skipped (multiallelic): {n_skipped_multiallelic:,}")
    print(f"    - Skipped (indels): {n_skipped_indel:,}")
    print(f"    - Skipped (low coverage): {n_skipped_low_coverage:,}")
    print(f"    - Missing genotypes: {n_missing_total:,}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python vcf_to_fasta.py <input.vcf> <output.fasta> <use_iupac> <min_samples>")
        sys.exit(1)

    vcf_to_fasta(sys.argv[1], sys.argv[2],
                 sys.argv[3].lower() == 'true',
                 sys.argv[4])
EOFPYTHON
    chmod +x "${PYTHON_SCRIPT}"
}

# Run phylogenetic analysis for single VCF
run_phylo_analysis() {
    local replicate=$1
    local run_id=$2
    local rec_rate=$3
    local chr_num=$4

    local rec_formatted=$(format_rec_rate "${rec_rate}")

    # File paths
    local VCF_FILE="${BASE_INPUT_DIR}/rep${replicate}/run${run_id}_rec${rec_formatted}_chr${chr_num}.vcf"
    local OUTPUT_DIR="${BASE_OUTPUT_DIR}/rep${replicate}/run${run_id}_rec${rec_formatted}_chr${chr_num}"

    # Check input
    if [[ ! -f "${VCF_FILE}" ]]; then
        log_error "VCF not found: ${VCF_FILE}"
        return 1
    fi

    mkdir -p "${OUTPUT_DIR}"

    local PREFIX="run${run_id}_rec${rec_formatted}_chr${chr_num}"
    #local fasta_file="${output_dir}/${prefix}.fasta"
    #local log_file="${output_dir}/${prefix}.log"

    # Set file paths
    local FASTA_FILE="${OUTPUT_DIR}/${PREFIX}.fasta"
    local ALIGNED_FASTA="${OUTPUT_DIR}/${PREFIX}_aligned.fasta"
    local FILTERED_VCF="${OUTPUT_DIR}/${PREFIX}_filtered.vcf"
    local MODEL_FILE="${OUTPUT_DIR}/${PREFIX}.iqtree"
    local TREE_FILE="${OUTPUT_DIR}/${PREFIX}.treefile"
    local LOG_FILE="${OUTPUT_DIR}/${PREFIX}.log"

    # Redirect all output to log file
#    exec > >(tee -a "${LOG_FILE}") 2>&1

    log_step "PHYLOGENETIC ANALYSIS PIPELINE"
    log_info "   VCF file: ${VCF_FILE}"
    log_info "   Output directory: ${OUTPUT_DIR}"
    log_info "   Prefix: ${PREFIX}"
    log_info "   IUPAC: ${USE_IUPAC}"
    log_info "   Minimum samples: ${MIN_SAMPLES}"
    log_info "   Threads: ${THREADS}"
    log_info "   Bootstrap replicates: ${BOOTSTRAP}"
    log_info "   Run alignment: ${RUN_ALIGNMENT}"
    log_info "   Min MAF: ${MIN_MAF}"
    log_info "   Max missing: ${MAX_MISSING}"


#    log_step "Processing rep${replicate}/run${run_id}_rec${rec_formatted}_chr${chr_num}"
#    log_info "VCF: ${vcf_file}"
#    log_info "Output: ${output_dir}"

    # Redirect to log
#    exec > >(tee -a "${log_file}") 2>&1
    exec > >(tee -a "${LOG_FILE}") 2>&1

    ################################################################################
    # STEP 1: FILTER VCF (OPTIONAL)
    ################################################################################

    if (( $(echo "${MIN_MAF} > 0" | bc -l) )) || (( $(echo "${MAX_MISSING} < 1.0" | bc -l) )); then
        log_step "STEP 1: Filtering VCF"

        log_info "Applying filters:"
        log_info "  - Minimum MAF: ${MIN_MAF}"
        log_info "  - Maximum missing data: ${MAX_MISSING}"

        if command -v bcftools &> /dev/null; then
            bcftools view \
                --min-af "${MIN_MAF}" \
                --max-af "$(echo "1.0 - ${MIN_MAF}" | bc -l)" \
                --threads "${THREADS}" \
                "${VCF_FILE}" | \
            bcftools view \
                --max-missing "${MAX_MISSING}" \
                -O v \
                -o "${FILTERED_VCF}"
       elif command -v vcftools &> /dev/null; then
            vcftools \
                --gzvcf "${VCF_FILE}" \
                --maf "${MIN_MAF}" \
                --max-missing "${MAX_MISSING}" \
                --recode \
                --recode-INFO-all \
                --out "${OUTPUT_DIR}/${PREFIX}_filtered"
           mv "${OUTPUT_DIR}/${PREFIX}_filtered.recode.vcf" "${FILTERED_VCF}"
       else
           log_warn "No filtering tool available, using original VCF"
           FILTERED_VCF="${VCF_FILE}"
       fi

       VCF_TO_USE="${FILTERED_VCF}"
       log_info "Filtered VCF saved: ${FILTERED_VCF}"
    else
       log_info "STEP 1: Skipping VCF filtering (no filters specified)"
       VCF_TO_USE="${VCF_FILE}"
    fi

    #####################################################################
    # Step 2: VCF to FASTA
    #####################################################################
    log_step "Step 1: Converting VCF to FASTA..."
    log_info "Creating Python script for VCF to FASTA conversion..."

    if ! python3 "${PYTHON_SCRIPT}" "${VCF_TO_USE}" "${FASTA_FILE}" "${USE_IUPAC}" "${MIN_SAMPLES}"; then
        log_error "VCF to FASTA conversion failed"
        return 1
    fi

    if [[ ! -f "${FASTA_FILE}" ]]; then
        log_error "FASTA file not created"
        return 1
    fi

    # Check FASTA has enough sequences
    n_seqs=$(grep -c "^>" "${FASTA_FILE}" || echo "0")
    if [[ ${n_seqs} -lt 4 ]]; then
        log_warn "Not enough sequences (${n_seqs}) for phylogenetic analysis, skipping"
        return 0
    fi

    log_info "✓ FASTA created: ${n_seqs} sequences"

    ################################################################################
    # STEP 3: MULTIPLE SEQUENCE ALIGNMENT (OPTIONAL)
    ################################################################################

    if [[ "${RUN_ALIGNMENT}" == "true" ]]; then
        log_step "STEP 3: Running multiple sequence alignment"

        log_info "Alignment tool: ${ALIGNMENT_TOOL}"

        if [[ "${ALIGNMENT_TOOL}" == "mafft" ]]; then
            log_info "Running MAFFT..."
            mafft \
                --auto \
                --thread "${THREADS}" \
                "${FASTA_FILE}" > "${ALIGNED_FASTA}"
        elif [[ "${ALIGNMENT_TOOL}" == "muscle" ]]; then
            log_info "Running MUSCLE..."
            muscle \
                -in "${FASTA_FILE}" \
                -out "${ALIGNED_FASTA}" \
                -threads "${THREADS}"
        else
            log_error "Unknown alignment tool: ${ALIGNMENT_TOOL}"
            exit 1
        fi

        if [[ ! -f "${ALIGNED_FASTA}" ]]; then
            log_error "Alignment failed!"
            exit 1
        fi

        log_info "✓ Aligned FASTA: ${ALIGNED_FASTA}"
        FASTA_FOR_TREE="${ALIGNED_FASTA}"
    else
        log_info "STEP 3: Skipping alignment (SNPs from VCF are already aligned)"
        FASTA_FOR_TREE="${FASTA_FILE}"
    fi


    ################################################################################
    # STEP 4: MODEL SELECTION WITH MODELFINDER
    ################################################################################

    log_step "STEP 4: Running ModelFinder to identify best substitution model"

    log_info "Testing model set: ${MODEL_SET}"
    log_info "This may take some time..."

    # Determine IQ-TREE command (iqtree2 or iqtree)
    if command -v iqtree2 &> /dev/null; then
        IQTREE_CMD="iqtree2"
    elif command -v iqtree &> /dev/null; then    # FIXED: Added condition
        IQTREE_CMD="iqtree"
    else
        IQTREE_CMD="${iqtree3}"
    fi

    # Run ModelFinder
    log_info "Running ModelFinder..."

    "${IQTREE_CMD}" \
        -s "${FASTA_FOR_TREE}" \
        -m MFP \
        -mset "${MODEL_SET}" \
        -nt "${THREADS}" \
        -pre "${OUTPUT_DIR}/${PREFIX}_modelfinder" \
        -quiet 2>&1 | tee -a "${LOG_FILE}"

    if [[ $? -ne 0 ]]; then
        log_error "ModelFinder failed"
        return 1
    fi

    # Extract best model
    local MODEL_LOG="${OUTPUT_DIR}/${PREFIX}_modelfinder.log"
    local BEST_MODEL=$(grep "Best-fit model:" "${MODEL_LOG}" | head -n 1 | awk '{print $3}')

    if [[ -z "${BEST_MODEL}" ]]; then
        log_warn "Could not extract best model, using GTR+G"
        BEST_MODEL="GTR+G"
    fi

    log_info "✓ ModelFinder completed"
    log_info ""
    log_info "=========================================="
    log_info "BEST MODEL IDENTIFIED: ${BEST_MODEL}"
    log_info "=========================================="
    log_info ""
    log_info "Model selection results saved to: ${MODEL_LOG}"

    # Display top models
    log_info "Top 5 models according to BIC:"
    grep -A 10 "List of models sorted by BIC scores:" "${MODEL_LOG}" | head -n 15

    ################################################################################
    # STEP 5: BUILD PHYLOGENETIC TREE WITH IQ-TREE
    ################################################################################

    log_step "STEP 5: Building phylogenetic tree with IQ-TREE"

    log_info "Using model: ${BEST_MODEL}"
    log_info "Bootstrap replicates: ${BOOTSTRAP}"
    log_info "Threads: ${THREADS}"

    # Run IQ-TREE with best model

    log_info "Step 3: Building phylogenetic tree..."
    "${IQTREE_CMD}" \
        -s "${FASTA_FOR_TREE}" \
        -m "${BEST_MODEL}" \
        -bb "${BOOTSTRAP}" \
        -nt "${THREADS}" \
        -pre "${OUTPUT_DIR}/${PREFIX}" \
        -redo \
        -quiet 2>&1 | tee -a "${LOG_FILE}"

    if [[ $? -ne 0 ]]; then
        log_error "IQ-TREE failed"
        return 1
    fi

    log_info "✓ Tree constructed: ${OUTPUT_DIR}/${PREFIX}.treefile"

    ################################################################################
    # STEP 6: GENERATE SUMMARY REPORT
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
Run: ${run_id}
Rec rate: ${rec_rate}
Chromosome: ${chr_num}
Output directory: ${OUTPUT_DIR}
Prefix: ${PREFIX}

FILTERING PARAMETERS
--------------------
Minimum MAF: ${MIN_MAF}
Maximum missing data: ${MAX_MISSING}

PARAMETERS
----------
Bootstrap replicates: ${BOOTSTRAP}
Model set tested: ${MODEL_SET}
Best model: ${BEST_MODEL}
Threads: ${THREADS}
Number of sequences: ${n_seqs}

OUTPUT FILES
------------
- FASTA: ${OUTPUT_DIR}/${PREFIX}.fasta
- Final Tree: ${OUTPUT_DIR}/${PREFIX}.treefile
- IQ-TREE report: ${OUTPUT_DIR}/${PREFIX}.iqtree
- ModelFinder log: ${OUTPUT_DIR}/${PREFIX}_modelfinder.log
- Full log: ${OUTPUT_DIR}/${PREFIX}.log

TREE FILES
----------
- Newick format: ${OUTPUT_DIR}/${PREFIX}.treefile
- With branch lengths: ${OUTPUT_DIR}/${PREFIX}.treefile
- Bootstrap support values included

MODEL SELECTION DETAILS
-----------------------
$(grep -A 5 "Best-fit model:" "${MODEL_LOG}" || echo "See ${MODEL_LOG} for details")

================================================================================
To view the tree:
  - Use FigTree, iTOL, or other tree viewer
  - Tree file: ${OUTPUT_DIR}/${PREFIX}.treefile

For more details, see:
  - Full IQ-TREE report: ${OUTPUT_DIR}/${PREFIX}.iqtree
  - ModelFinder report: ${OUTPUT_DIR}/${PREFIX}_modelfinder.iqtree

================================================================================
EOFREPORT

    log_info "✓ Analysis complete"
    log_info "Tree file: ${OUTPUT_DIR}/${PREFIX}.treefile"

    return 0
}

###############################################################################
# MAIN EXECUTION
###############################################################################

log_step "Phylogenetic Analysis Batch Processing"

# Create conversion script
create_vcf2fasta_script
log_info "Created VCF to FASTA conversion script"

# Calculate total jobs
N_REC_RATES=${#REC_RATES[@]}
N_CHRS=${#CHROMOSOMES[@]}
TOTAL_JOBS=$((N_REPLICATES * N_REC_RATES * N_CHRS))

log_info "Configuration:"
log_info "  Replicates: ${N_REPLICATES}"
log_info "  Rec rates: ${N_REC_RATES}"
log_info "  Chromosomes: ${N_CHRS}"
log_info "  Total jobs: ${TOTAL_JOBS}"
log_info "  Bootstrap: ${BOOTSTRAP}"
log_info "  Threads: ${THREADS}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # SLURM array mode
    TASK_ID=${SLURM_ARRAY_TASK_ID}

    if [[ ${TASK_ID} -lt 1 || ${TASK_ID} -gt ${TOTAL_JOBS} ]]; then
        log_error "Invalid task ID ${TASK_ID} (valid: 1-${TOTAL_JOBS})"
        exit 1
    fi

    # Task layout: rep1 has all rec_rates × chrs, rep2 has all, etc.
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

    run_phylo_analysis "${REPLICATE}" "${RUN_ID}" "${REC_RATE}" "${CHR_NUM}"

else
    # Sequential mode
    log_info "Running in sequential mode"

    CURRENT_JOB=0
    FAILED_JOBS=0

    for replicate in $(seq 1 ${N_REPLICATES}); do
        run_id=${replicate}

        log_info ""
        log_step "Replicate ${replicate} (run${run_id})"

        for rec_rate in "${REC_RATES[@]}"; do
            for chr_num in "${CHROMOSOMES[@]}"; do
                CURRENT_JOB=$((CURRENT_JOB + 1))

                log_info ""
                log_info "Job ${CURRENT_JOB}/${TOTAL_JOBS}"

                if ! run_phylo_analysis "${replicate}" "${run_id}" "${rec_rate}" "${chr_num}"; then
                    FAILED_JOBS=$((FAILED_JOBS + 1))
                    log_warn "Job ${CURRENT_JOB} failed, continuing..."
                fi
            done
        done
    done

    log_info ""
    log_step "Sequential processing complete"
    log_info "Successful: $((TOTAL_JOBS - FAILED_JOBS))/${TOTAL_JOBS}"
    log_info "Failed: ${FAILED_JOBS}/${TOTAL_JOBS}"

    [[ ${FAILED_JOBS} -gt 0 ]] && exit 1
fi

log_step "✓ Phylogenetic analysis batch completed"


log_info ""
log_info "=========================================="
log_info "SUMMARY"
log_info "=========================================="
log_info "✓ VCF to FASTA conversion: SUCCESS"
if [[ "${RUN_ALIGNMENT}" == "true" ]]; then
    log_info "✓ Multiple sequence alignment: SUCCESS"
fi
log_info "✓ Model selection (ModelFinder): SUCCESS"
log_info "  → Best model: ${BEST_MODEL}"
log_info "✓ Phylogenetic tree construction: SUCCESS"
log_info "  → Bootstrap replicates: ${BOOTSTRAP}"
log_info ""
log_info "Output files:"
log_info "  • Tree file: ${OUTPUT_DIR}/${PREFIX}.treefile"
log_info "  • Full report: ${OUTPUT_DIR}/${PREFIX}.iqtree"
log_info "  • Summary: ${REPORT_FILE}"
log_info "=========================================="


exit 0
