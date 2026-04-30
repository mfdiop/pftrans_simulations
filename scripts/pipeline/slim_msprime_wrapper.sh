#!/usr/bin/env bash

# Chromosome lengths (based on Pf3D7 reference)
declare -A CHR_LEN=(
  [1]=640851 [2]=947102 [3]=1060597 [4]=1203928 [5]=1343531 [6]=1418587
  [7]=1684441 [8]=1472805 [9]=1545076 [10]=1687656 [11]=2038348
  [12]=2271718 [13]=2925725 [14]=3291872
)

# Example selection positions (modify as needed)
declare -A SEL_POS=(
  [1]=114473 [2]=798866 [3]=139452 [4]=1100039 [5]=702042 [6]=851783
  [7]=900000 [8]=1317660 [9]=900000 [10]=1000000 [11]=1293674
  [12]=1300000 [13]=1500000 [14]=1996254
)

# Loop over chromosomes
for i in {1..14}; do
  seqlen=${CHR_LEN[$i]}
  selpos=${SEL_POS[$i]}
  echo "Running chromosome $i: length=$seqlen, selection site=$selpos"
  
  Rscript wrapper_codes/slim_msprime_wrapper.R \
    --chrno $i \
    --seqlen $seqlen \
    --selpos $selpos \
    --outdir "results/chr${i}_sim.trees"
done
