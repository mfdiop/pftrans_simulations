#!/usr/bin/env python3

"""
Compute patristic distances from IQ-TREE phylogenetic trees
Run BEFORE evaluation script to generate phylo_results.csv for each scenario
"""

from Bio import Phylo
import pandas as pd
from pathlib import Path
import itertools


SCENARIO_PATTERN = "MIG*_*_*_*_*"


def compute_patristic_distances(treefile):
    """
    Compute all pairwise patristic distances from a phylogenetic tree
    """

    tree = Phylo.read(treefile, "newick")

    terminals = tree.get_terminals()
    terminal_names = [t.name for t in terminals]

    results = []
    for id1, id2 in itertools.combinations(terminal_names, 2):
        distance = tree.distance(id1, id2)
        results.append({
            "Id1": id1,
            "Id2": id2,
            "patristic_distance": distance
        })

    return pd.DataFrame(results)


def process_all_scenarios(phylo_dir, outdir="phylo_results"):

    phylo_dir = Path(phylo_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    scenario_dirs = sorted(phylo_dir.glob(SCENARIO_PATTERN))
    print(f"\nFound {len(scenario_dirs)} scenarios")

    missing_trees = []
    failed_scenarios = []
    processed = []

    for idx, scenario_dir in enumerate(scenario_dirs, 1):

        scenario_id = scenario_dir.name

        print("\n" + "=" * 60)
        print(f"Scenario: {scenario_id}")
        print("=" * 60)

        if idx % 50 == 0:
            print(f"Progress: {idx}/{len(scenario_dirs)}")

        # Expect exactly one treefile per scenario
        treefile = scenario_dir / f"{scenario_id}.treefile"

        if not treefile.exists():
            print(f"  ✗ Missing treefile: {treefile}")
            missing_trees.append(scenario_id)
            continue

        print(f"  ✓ Found treefile")

        try:
            distances_df = compute_patristic_distances(treefile)

            scenario_outdir = outdir / scenario_id
            scenario_outdir.mkdir(parents=True, exist_ok=True)

            output_file = scenario_outdir / "phylo_results.csv"

            distances_df.to_csv(output_file, index=False)

            processed.append(scenario_id)
            print(f"  ✓ Saved: {output_file}")

        except Exception as e:
            print(f"  ✗ Error processing {scenario_id}: {e}")
            failed_scenarios.append({
                "scenario": scenario_id,
                "error": str(e)
            })

    # ---- Save logs AFTER loop ----
    logs_dir = Path(phylo_dir.parent) / "logs"
    logs_dir.mkdir(exist_ok=True)

    pd.DataFrame({"scenario": processed}).to_csv(
        logs_dir / "processed_files.csv", index=False
    )

    pd.DataFrame({"scenario": missing_trees}).to_csv(
        logs_dir / "missing_trees.csv", index=False
    )

    pd.DataFrame(failed_scenarios).to_csv(
        logs_dir / "failed_scenarios.csv", index=False
    )

    print("\n" + "=" * 60)
    print("Processing complete")
    print(f"  ✓ Success: {len(processed)}")
    print(f"  ⚠ Missing trees: {len(missing_trees)}")
    print(f"  ✗ Failed: {len(failed_scenarios)}")
    print(f"  Logs saved to: {logs_dir}")
    print("=" * 60)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phylo-dir",
        default="sim_migration/phylo_results",
        help="Directory containing scenario folders"
    )
    parser.add_argument(
        "--outdir",
        default="phylo_results",
        help="Output directory"
    )

    args = parser.parse_args()

    process_all_scenarios(args.phylo_dir, args.outdir)
