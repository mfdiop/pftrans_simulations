#!/bin/bash
#SBATCH --job-name=figure2
#SBATCH --output=logs/figure2.out
#SBATCH --error=logs/figure2.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mdiop@mrc.gm

Rscript codes/figure2_ibd_distribution_baseline.R

Rscript codes/figure2_impact_recombination_ibd_accuracy.R

echo "Job Completed Successfully"

exit 0