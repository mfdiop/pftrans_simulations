#!/usr/bin/env python3
import argparse
import sys
from collections import Counter
import math

# Import pysam with a helpful error if it's missing
try:
    import pysam
except Exception as e:
    sys.stderr.write("pysam is required. Install with: pip install pysam\n")
    raise

# Import numpy with a helpful error if it's missing
try:
    import numpy as np
except Exception as e:
    sys.stderr.write("numpy is required. Install with: pip install numpy\n")
    raise

"""
compute_ibs_matrix.py

Compute pairwise IBS (Identity-By-State) from a VCF file.

Outputs a tab-delimited matrix of pairwise IBS proportions (0.0 - 1.0),
where IBS at a site is the number of shared alleles (0,1,2) divided by 2,
and the reported value is the average across all non-missing sites for the pair.

Usage:
    python compute_ibs_matrix.py --vcf input.vcf.gz --out ibs_matrix.tsv

Requires:
    pysam
"""



def parse_args():
    p = argparse.ArgumentParser(description="Compute pairwise IBS from a VCF")
    p.add_argument("--vcf", "-v", required=True, help="Input VCF (can be bgzipped .vcf.gz)")
    p.add_argument("--out", "-o", required=True, help="Output TSV file (matrix)")
    p.add_argument("--min-af", type=float, default=0.0, help="Skip sites with minor allele frequency < MIN_AF (default: 0.0)")
    p.add_argument("--biallelic-only", action="store_true", default=True, help="Only use biallelic sites (skip multiallelic)")
    
    return p.parse_args()


def gt_to_alleles(gt):
    """
    Convert pysam GT tuple to a list of alleles (ints).
    Returns None if genotype is missing (any allele is None).
    """
    if gt is None:
        return None
    # gt may be like (0, 1) for diploid, (0,) for haploid
    if any(a is None for a in gt):
        return None
    return list(gt)


def shared_alleles_count(a1, a2):
    """
    Count number of shared alleles between two genotypes a1 and a2.
    a1 and a2 are lists of allele integers, assumed non-missing.
    Works for diploid/haploid by multiset intersection.
    Returns integer 0..min(len(a1),len(a2))* number of matches (for diploid up to 2).
    """
    c1 = Counter(a1)
    c2 = Counter(a2)
    shared = 0
    for allele in c1:
        shared += min(c1[allele], c2.get(allele, 0))
    return shared


def compute_ibs_matrix(vcf_path, biallelic_only=False, min_af=0.0):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    n = len(samples)
    if n < 2:
        raise ValueError("VCF must contain at least two samples")

    # Accumulators: sum of shared alleles (0..2 per site) and number of compared sites
    shared_sum = np.zeros((n, n), dtype=np.float64)
    compared = np.zeros((n, n), dtype=np.int64)

    for rec in vcf:
        # Optionally skip multiallelic sites
        if biallelic_only and len(rec.alleles) != 2:
            continue

        # compute allele counts across samples to get MAF if requested
        if min_af > 0.0:
            allele_counts = Counter()
            total_alleles = 0
            for s in samples:
                gt = rec.samples[s].get("GT")
                alleles = gt_to_alleles(gt)
                if alleles is None:
                    continue
                for a in alleles:
                    allele_counts[a] += 1
                    total_alleles += 1
            if total_alleles == 0:
                continue
            # compute minor allele frequency ignoring missing
            freqs = [allele_counts[a] / total_alleles for a in allele_counts]
            maf = min(freqs) if freqs else 0.0
            if maf < min_af:
                continue

        # Collect genotypes for all samples; skip missing
        genos = []
        present_idx = []
        for i, s in enumerate(samples):
            gt = rec.samples[s].get("GT")
            alleles = gt_to_alleles(gt)
            if alleles is None:
                genos.append(None)
            else:
                genos.append(alleles)
                present_idx.append(i)

        # If fewer than 2 samples have genotype, skip
        if len(present_idx) < 2:
            continue

        # Update pairwise shared allele counts for present samples
        # For each pair i <= j, compute shared and update symmetric entries
        for a_i_idx_idx, i in enumerate(present_idx):
            g_i = genos[i]
            for j in present_idx[a_i_idx_idx:]:
                g_j = genos[j]
                shared = shared_alleles_count(g_i, g_j)  # 0..2 typically
                shared_sum[i, j] += shared
                compared[i, j] += 1
                if i != j:
                    shared_sum[j, i] += shared
                    compared[j, i] += 1

    # Build IBS proportion matrix: average shared alleles / 2 (normalize to 0..1)
    ibs = np.full((n, n), np.nan, dtype=np.float64)
    nonzero = compared > 0
    ibs[nonzero] = (shared_sum[nonzero] / compared[nonzero]) / 2.0
    # Diagonal: for completeness set to 1.0 if compared>0 (self vs self)
    for i in range(n):
        if compared[i, i] > 0:
            ibs[i, i] = 1.0

    return samples, ibs, compared


def write_matrix(out_path, samples, matrix):
    n = len(samples)
    with open(out_path, "w") as fo:
        fo.write("#samples\t" + "\t".join(samples) + "\n")
        for i, s in enumerate(samples):
            row = []
            for j in range(n):
                v = matrix[i, j]
                if v is None or (isinstance(v, float) and math.isnan(v)):
                    row.append("NA")
                else:
                    row.append(f"{v:.6f}")
            fo.write(s + "\t" + "\t".join(row) + "\n")


def main():
    args = parse_args()
    samples, ibs, compared = compute_ibs_matrix(
        args.vcf, biallelic_only=args.biallelic_only, min_af=args.min_af
    )
    write_matrix(args.out, samples, ibs)
    sys.stdout.write(f"Wrote IBS matrix to {args.out}\n")


if __name__ == "__main__":
    main()