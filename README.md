# Benchmarking Methods for Malaria Transmission Inference — Simulation Study

PhD Thesis — Objective 1 | London School of Hygiene and Tropical Medicine

## Overview

This project benchmarks genomic methods (IBD, IBS, phylogenetics) for inferring malaria parasite transmission using forward-time simulations (SLiM + msprime).

## Project Structure

```
pftrans_simulations/
├── scripts/
│   ├── slim/         # SLiM forward-time simulation scripts
│   ├── pipeline/     # Wrapper scripts for running the full pipeline
│   ├── analysis/     # Post-simulation analysis scripts (R, Python)
│   ├── figures/      # Figure generation scripts
│   └── hpc/          # HPC/SLURM job submission scripts
│
├── simulations/      # Simulation scenarios (data NOT tracked by git — see .gitignore)
│   ├── single_run/            # Baseline single-replicate scenario
│   ├── multiple_runs/         # Multi-replicate scenarios (recombination, genomic)
│   └── malaria_transmission_study/  # Migration scenario
│
├── results/
│   ├── figures/
│   │   ├── main/              # Main manuscript figures
│   │   ├── supplementary/     # Current supplementary figures
│   │   └── supplementary_v1/  # Previous version of supplementary figures
│   └── tables/                # Summary tables and real data inputs
│
├── manuscript/       # Manuscript Rmd, LaTeX, bibliography, and submission files
│   └── notes/        # Working notes
│
├── docs/             # Code reviews, analysis guides, thesis package notes
│
├── simulation_design.json   # Simulation parameter design
├── environment.yml          # Conda environment specification
└── pftrans_simulations.Rproj  # RStudio project file
```

## Simulation Scenarios

| Folder | Description |
|--------|-------------|
| `simulations/single_run/` | Baseline scenario: single replicate, single population |
| `simulations/multiple_runs/` | Multi-replicate: recombination rate sweep, genomic evaluation |
| `simulations/malaria_transmission_study/` | Migration scenario: two-population model |

## Quick Start

See `scripts/analysis/quick_start.md` for setup and running instructions.

## Dependencies

- SLiM (forward-time simulator)
- msprime / tskit (Python)
- R packages: see `environment.yml`
- Conda environment: `conda env create -f environment.yml`

## Large Data

Simulation output data (~10GB+) is excluded from this repository (see `.gitignore`).
Data is stored on the LSHTM HPC cluster and locally via OneDrive.

## Reference

Diop MF. *Benchmarking genomic metrics for malaria parasite transmission inference*. PhD Thesis, LSHTM.
