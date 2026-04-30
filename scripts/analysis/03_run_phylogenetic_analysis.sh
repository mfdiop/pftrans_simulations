#!/bin/bash

#SBATCH --job-name=phylo_iqtree
#SBATCH --output=logs/phylo_%j.out
#SBATCH --error=logs/phylo_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

################################################################################
# Script: run_phylogenetic_analysis.sh
# Purpose: Complete phylogenetic analysis pipeline from VCF to IQ-TREE
# Author: [Your Name]
# Date: 2025-11-05
#
# Pipeline steps:
#   1. Convert VCF to FASTA format
#   2. (Optional) Multiple sequence alignment
#   3. Run ModelFinder to identify best substitution model
#   4. Run IQ-TREE with best model
#   5. Generate summary report
#
# Requirements:
#   - vcftools or bcftools
#   - Python with BioPython
#   - IQ-TREE (http://www.iqtree.org/)
#   - (Optional) MAFFT or MUSCLE for alignment
#
# Usage:
#   sbatch run_phylogenetic_analysis.sh --vcf input.vcf.gz
#   bash run_phylogenetic_analysis.sh --vcf input.vcf.gz --threads 8
################################################################################

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

################################################################################
# DEFAULT CONFIGURATION
################################################################################

# Input/Output
VCF_FILE=""
OUTPUT_DIR="results/phylo_output"
PREFIX="phylo"

# Computational parameters
THREADS=${SLURM_CPUS_PER_TASK:-4}
MEMORY="32G"

# Convert VCF to Fasta
USE_IUPAC=false
MIN_SAMPLES=10

# Analysis options
RUN_ALIGNMENT=false
ALIGNMENT_TOOL="mafft"  # Options: mafft, muscle, none
MIN_MAF=0.0             # Minimum minor allele frequency filter
MAX_MISSING=1.0         # Maximum proportion of missing data per site
BOOTSTRAP=1000          # Number of bootstrap replicates
MODEL_SET="ALL"         # Model set for ModelFinder (ALL, DNA, or specific models)

# IQ-TREE options
IQTREE_EXTRA=""         # Additional IQ-TREE parameters

################################################################################
# PARSE COMMAND LINE ARGUMENTS
################################################################################

usage() {
    cat << EOF
Usage: $0 --vcf <input.vcf> [OPTIONS]

Required Arguments:
  --vcf FILE              Input VCF file (can be gzipped)

Optional Arguments:
  --output-dir DIR        Output directory (default: phylo_output)
  --prefix STRING         Output file prefix (default: phylo)
  --use-iupac             Use IUPAC for heterozygous genotypes
  --min-samples INT       Minimum number of samples to create a fasta file (default: 10)
  --threads INT           Number of threads (default: 4)
  --bootstrap INT         Number of bootstrap replicates (default: 1000)
  --align                 Run multiple sequence alignment
  --aligner TOOL          Alignment tool: mafft or muscle (default: mafft)
  --min-maf FLOAT         Minimum minor allele frequency (default: 0.0)
  --max-missing FLOAT     Maximum missing data proportion (default: 1.0)
  --model-set STRING      Model set for testing: ALL, DNA (default: ALL)
  --help                  Show this help message

Examples:
  # Basic usage
  bash $0 --vcf input.vcf.gz

  # With filtering and alignment
  bash $0 --vcf input.vcf.gz --min-maf 0.05 --max-missing 0.2 --align

  # Custom output and threading
  bash $0 --vcf input.vcf.gz --output-dir results --threads 16 --bootstrap 10000

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf)
            VCF_FILE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --prefix)
            PREFIX="$2"
            shift 2
            ;;
        --use-iupac)
            USE_IUPAC=true
            shift 2
            ;;
        --min-samples)
            MIN_SAMPLES="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --bootstrap)
            BOOTSTRAP="$2"
            shift 2
            ;;
        --align)
            RUN_ALIGNMENT=true
            shift
            ;;
        --aligner)
            ALIGNMENT_TOOL="$2"
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
        --model-set)
            MODEL_SET="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "${VCF_FILE}" ]]; then
    log_error "Error: --vcf argument is required"
    usage
    exit 1
fi

if [[ ! -f "${VCF_FILE}" ]]; then
    log_error "Error: VCF file not found: ${VCF_FILE}"
    exit 1
fi

################################################################################
# SETUP
################################################################################

# Create output directory
mkdir -p "${OUTPUT_DIR}" logs

# Set file paths
FASTA_FILE="${OUTPUT_DIR}/${PREFIX}.fasta"
ALIGNED_FASTA="${OUTPUT_DIR}/${PREFIX}_aligned.fasta"
FILTERED_VCF="${OUTPUT_DIR}/${PREFIX}_filtered.vcf"
MODEL_FILE="${OUTPUT_DIR}/${PREFIX}.iqtree"
TREE_FILE="${OUTPUT_DIR}/${PREFIX}.treefile"
LOG_FILE="${OUTPUT_DIR}/${PREFIX}.log"

# Redirect all output to log file
exec > >(tee -a "${LOG_FILE}") 2>&1

log_step "PHYLOGENETIC ANALYSIS PIPELINE"
log_info "VCF file: ${VCF_FILE}"
log_info "Output directory: ${OUTPUT_DIR}"
log_info "Prefix: ${PREFIX}"
log_info "IUPAC: ${USE_IUPAC}"
log_info "Minimum samples: ${MIN_SAMPLES}"
log_info "Threads: ${THREADS}"
log_info "Bootstrap replicates: ${BOOTSTRAP}"
log_info "Run alignment: ${RUN_ALIGNMENT}"
log_info "Min MAF: ${MIN_MAF}"
log_info "Max missing: ${MAX_MISSING}"

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

################################################################################
# STEP 2: CONVERT VCF TO FASTA
################################################################################

log_step "STEP 2: Converting VCF to FASTA"

log_info "Creating Python script for VCF to FASTA conversion..."

cat > "${OUTPUT_DIR}/vcf_to_fasta.py" << 'EOFPYTHON'
#!/usr/bin/env python3
"""
Convert VCF file to FASTA format for phylogenetic analysis
Handles both phased and unphased genotypes

Robust VCF to FASTA converter for phylogenetic analysis
Handles common VCF format issues and provides detailed diagnostics

Usage:
    python vcf_to_fasta_robust.py input.vcf output.fasta
    python vcf_to_fasta_robust.py input.vcf.gz output.fasta --min-samples 10
"""

import sys
import gzip
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def open_file(filename):
    """Open file, handling gzip compression"""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')

def parse_vcf_header(vcf_file):
    """Extract sample names from VCF header"""
    print(f"Parsing VCF header from: {vcf_file}")
    
    with open_file(vcf_file) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('##'):
                continue  # Skip meta-information lines
            
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                
                if len(fields) < 10:
                    raise ValueError(
                        f"Invalid VCF header at line {line_num}: "
                        f"expected at least 10 columns, got {len(fields)}\n"
                        f"Line: {line.strip()}"
                    )
                
                samples = fields[9:]
                
                if len(samples) == 0:
                    raise ValueError("No samples found in VCF header")
                
                print(f"  ✓ Found {len(samples)} samples")
                return samples
    
    raise ValueError("No header line (#CHROM) found in VCF")

def validate_vcf(vcf_file):
    """Validate VCF format and return statistics"""
    print("\nValidating VCF format...")
    
    stats = {
        'total_lines': 0,
        'header_lines': 0,
        'data_lines': 0,
        'empty_lines': 0,
        'malformed_lines': 0
    }
    
    with open_file(vcf_file) as f:
        for line in f:
            stats['total_lines'] += 1
            
            if not line.strip():
                stats['empty_lines'] += 1
                continue
            
            if line.startswith('#'):
                stats['header_lines'] += 1
                continue
            
            fields = line.strip().split('\t')
            if len(fields) >= 8:  # Minimum VCF columns
                stats['data_lines'] += 1
            else:
                stats['malformed_lines'] += 1
    
    print(f"  Total lines: {stats['total_lines']:,}")
    print(f"  Header lines: {stats['header_lines']:,}")
    print(f"  Data lines: {stats['data_lines']:,}")
    print(f"  Empty lines: {stats['empty_lines']:,}")
    print(f"  Malformed lines: {stats['malformed_lines']:,}")
    
    if stats['data_lines'] == 0:
        raise ValueError("No data lines found in VCF!")
    
    if stats['malformed_lines'] > 0:
        print(f"  ⚠ WARNING: {stats['malformed_lines']} malformed lines will be skipped")
    
    print("  ✓ VCF validation passed")
    return stats

def vcf_to_fasta(vcf_file, output_file, use_iupac=True, min_samples=0):
    """
    Convert VCF to FASTA format with robust error handling
    
    Args:
        vcf_file: Input VCF file path
        output_file: Output FASTA file path
        use_iupac: If True, use IUPAC codes for heterozygous sites
        min_samples: Minimum number of non-missing samples required per site
    """

    print("="*70)
    print(f"Reading VCF: {vcf_file}")
    print("VCF TO FASTA CONVERSION")
    print("="*70)
    
    # Validate VCF
    validate_vcf(vcf_file)

    # Get sample names
    try:
        samples = parse_vcf_header(vcf_file)
    except ValueError as e:
        print(f"\n✗ ERROR: {e}")
        sys.exit(1)
    
    n_samples = len(samples)
    print(f"\nProcessing {n_samples} samples...")
    
    # Initialize sequences for each sample
    sequences = {sample: [] for sample in samples}
    n_sites = 0
    n_missing = 0
    
    # IUPAC codes for heterozygous sites
    iupac_codes = {
        ('A', 'G'): 'R', ('G', 'A'): 'R',
        ('C', 'T'): 'Y', ('T', 'C'): 'Y',
        ('G', 'C'): 'S', ('C', 'G'): 'S',
        ('A', 'T'): 'W', ('T', 'A'): 'W',
        ('G', 'T'): 'K', ('T', 'G'): 'K',
        ('A', 'C'): 'M', ('C', 'A'): 'M',
    }

    # Statistics
    stats = defaultdict(int)
    valid_bases = set('ACGT')
    
    # Read VCF and extract genotypes
    with open_file(vcf_file) as f:
        for line in f:
            # Skip headers
            if line.startswith('#'):
                continue
            
            # Skip empty lines
            if not line.strip():
                stats['empty_lines'] += 1
                continue
            
            fields = line.strip().split('\t')
            
            # Validate line format
            if len(fields) < 10:
                stats['malformed_lines'] += 1
                if stats['malformed_lines'] <= 5:  # Show first 5 errors
                    print(f"  ⚠ Line {line_num}: only {len(fields)} fields (expected ≥10)")
                continue

            # Extract variant information
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            
            # Skip multi-allelic sites
            if ',' in alt:
                stats['multiallelic'] += 1
                continue
            
            # Validate REF and ALT
            if not (set(ref).issubset(valid_bases) and set(alt).issubset(valid_bases)):
                stats['invalid_bases'] += 1
                continue
            
            # Skip indels (focus on SNPs)
            if len(ref) != 1 or len(alt) != 1:
                stats['indels'] += 1
                continue

            # Get genotypes for all samples
            genotypes = fields[9:]
            
            if len(genotypes) != n_samples:
                stats['sample_mismatch'] += 1
                continue
            
            # Count non-missing genotypes
            non_missing = 0
            
            # Process each sample's genotype
            site_valid = True
            for sample, gt_field in zip(samples, genotypes):
                gt = gt_field.split(':')[0]  # Extract GT field
                
                # Parse genotype
                if '/' in gt:
                    alleles = gt.split('/')
                elif '|' in gt:
                    alleles = gt.split('|')
                else:
                    # Handle cases like "0" or "1" without separator
                    if gt in ['.', './.', '.|.']:
                        alleles = ['.', '.']
                    else:
                        alleles = [gt, gt]
                
                # Convert to nucleotides
                if '.' in alleles or alleles == ['.']:
                    # Missing data
                    sequences[sample].append('N')
                    n_missing += 1
                    stats['missing_gt'] += 1
                elif len(set(alleles)) == 1:
                    # Homozygous
                    try:
                        allele_idx = int(alleles[0])
                        if allele_idx == 0:
                            sequences[sample].append(ref)
                        elif allele_idx == 1:
                            sequences[sample].append(alt)
                        else:
                            sequences[sample].append('N')
                            stats['invalid_gt'] += 1
                        non_missing += 1
                    except ValueError:
                        sequences[sample].append('N')
                        stats['invalid_gt'] += 1
                else:
                    # Heterozygous
                    non_missing += 1
                    if use_iupac:
                        try:
                            a1 = ref if alleles[0] == '0' else alt
                            a2 = ref if alleles[1] == '0' else alt
                            iupac = iupac_codes.get((a1, a2), 'N')
                            sequences[sample].append(iupac)
                        except (ValueError, KeyError):
                            sequences[sample].append('N')
                            stats['invalid_gt'] += 1
                    else:
                        # Use first allele only
                        try:
                            allele_idx = int(alleles[0])
                            if allele_idx == 0:
                                sequences[sample].append(ref)
                            else:
                                sequences[sample].append(alt)
                        except ValueError:
                            sequences[sample].append('N')
                            stats['invalid_gt'] += 1
            
            # Check if site meets minimum sample requirement
            if non_missing < int(min_samples):
                stats['insufficient_samples'] += 1
                # Remove this site from all sequences
                for sample in samples:
                    sequences[sample].pop()
                continue
            
            stats['sites_processed'] += 1
            
            if stats['sites_processed'] % 10000 == 0:
                print(f"  Processed {stats['sites_processed']:,} sites...")
    
    #print(f"\nProcessed {n_sites:,} sites")
    #print(f"Missing data: {n_missing:,} genotypes ({n_missing/(n_sites*n_samples)*100:.2f}%)")

    # Print statistics
    print("\n" + "="*70)
    print("CONVERSION STATISTICS")
    print("="*70)
    print(f"Sites processed: {stats['sites_processed']:,}")
    print(f"\nSkipped:")
    print(f"  - Malformed lines: {stats['malformed_lines']:,}")
    print(f"  - Multi-allelic: {stats['multiallelic']:,}")
    print(f"  - Invalid bases: {stats['invalid_bases']:,}")
    print(f"  - Indels: {stats['indels']:,}")
    print(f"  - Sample mismatch: {stats['sample_mismatch']:,}")
    print(f"  - Insufficient samples: {stats['insufficient_samples']:,}")
    print(f"\nGenotypes:")
    print(f"  - Missing (N): {stats['missing_gt']:,}")
    print(f"  - Invalid: {stats['invalid_gt']:,}")
    
    if stats['sites_processed'] == 0:
        print("\n✗ ERROR: No valid sites found in VCF!")
        sys.exit(1)
    
    # Write FASTA file
    print(f"Writing FASTA: {output_file}")
    records = []
    for sample in samples:
        seq = ''.join(sequences[sample])
        record = SeqRecord(Seq(seq), id=sample, description="")
        records.append(record)
    
    with open(output_file, 'w') as out:
        SeqIO.write(records, out, "fasta")
    
    print(f"✓ FASTA file created: {output_file}")
    print(f"  Samples: {n_samples}")
    print(f"  Sequence length: {len(seq):,} bp")


#    print(f"\nWriting FASTA file: {output_file}")
#    
#    try:
#        with open(output_file, 'w') as out:
#            for sample in samples:
#                seq = ''.join(sequences[sample])
#                out.write(f">{sample}\n")
#                
#                # Write sequence in 80-character lines
#                for i in range(0, len(seq), 80):
#                    out.write(seq[i:i+80] + '\n')
#        
#        print(f"✓ FASTA file created successfully!")
#        print(f"  - Samples: {n_samples}")
#        print(f"  - Sequence length: {len(seq):,} bp")
#        print("="*70)
#        
#    except IOError as e:
#        print(f"\n✗ ERROR writing output file: {e}")
#        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python vcf_to_fasta.py <input.vcf> <output.fasta> <use_iupac> <min_samples>")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    output_file = sys.argv[2]
    use_iupac = sys.argv[3]
    min_samples = sys.argv[4]
    
    vcf_to_fasta(vcf_file, output_file, use_iupac, min_samples)
EOFPYTHON

chmod +x "${OUTPUT_DIR}/vcf_to_fasta.py"

# Run conversion
log_info "Converting VCF to FASTA..."
python3 "${OUTPUT_DIR}/vcf_to_fasta.py" "${VCF_TO_USE}" "${FASTA_FILE}" "${USE_IUPAC}" "${MIN_SAMPLES}"

if [[ ! -f "${FASTA_FILE}" ]]; then
    log_error "FASTA conversion failed!"
    exit 1
fi

log_info "✓ FASTA file created: ${FASTA_FILE}"

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
${IQTREE_CMD} \
    -s "${FASTA_FOR_TREE}" \
    -m MFP \
    -mset "${MODEL_SET}" \
    -nt "${THREADS}" \
    -pre "${OUTPUT_DIR}/${PREFIX}_modelfinder"

if [[ $? -ne 0 ]]; then
    log_error "ModelFinder failed!"
    exit 1
fi

# Extract best model
MODEL_LOG="${OUTPUT_DIR}/${PREFIX}_modelfinder.log"
BEST_MODEL=$(grep "Best-fit model:" "${MODEL_LOG}" | head -n 1 | awk '{print $3}')

if [[ -z "${BEST_MODEL}" ]]; then
    log_error "Could not extract best model from ModelFinder output"
    exit 1
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
${IQTREE_CMD} \
    -s "${FASTA_FOR_TREE}" \
    -m "${BEST_MODEL}" \
    -bb "${BOOTSTRAP}" \
    -nt "${THREADS}" \
    -pre "${OUTPUT_DIR}/${PREFIX}" \
    ${IQTREE_EXTRA}

if [[ $? -ne 0 ]]; then
    log_error "IQ-TREE phylogeny reconstruction failed!"
    exit 1
fi

log_info "✓ Phylogenetic tree constructed"

################################################################################
# STEP 6: GENERATE SUMMARY REPORT
################################################################################

log_step "STEP 6: Generating summary report"

REPORT_FILE="${OUTPUT_DIR}/${PREFIX}_report.txt"

cat > "${REPORT_FILE}" << EOFREPORT
================================================================================
PHYLOGENETIC ANALYSIS REPORT
================================================================================
Generated: $(date)
================================================================================

INPUT DATA
----------
VCF file: ${VCF_FILE}
Output directory: ${OUTPUT_DIR}
Prefix: ${PREFIX}

FILTERING PARAMETERS
--------------------
Minimum MAF: ${MIN_MAF}
Maximum missing data: ${MAX_MISSING}

ANALYSIS PARAMETERS
-------------------
Alignment performed: ${RUN_ALIGNMENT}
$(if [[ "${RUN_ALIGNMENT}" == "true" ]]; then echo "Alignment tool: ${ALIGNMENT_TOOL}"; fi)
Model set tested: ${MODEL_SET}
Bootstrap replicates: ${BOOTSTRAP}
Threads used: ${THREADS}

RESULTS
-------
Best substitution model: ${BEST_MODEL}

OUTPUT FILES
------------
1. FASTA file: ${FASTA_FILE}
$(if [[ "${RUN_ALIGNMENT}" == "true" ]]; then echo "2. Aligned FASTA: ${ALIGNED_FASTA}"; fi)
3. ModelFinder log: ${OUTPUT_DIR}/${PREFIX}_modelfinder.log
4. Final tree: ${OUTPUT_DIR}/${PREFIX}.treefile
5. IQ-TREE log: ${OUTPUT_DIR}/${PREFIX}.iqtree
6. Full log: ${OUTPUT_DIR}/${PREFIX}.log

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

log_info "Report saved: ${REPORT_FILE}"

################################################################################
# FINAL SUMMARY
################################################################################

log_step "ANALYSIS COMPLETE"

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
