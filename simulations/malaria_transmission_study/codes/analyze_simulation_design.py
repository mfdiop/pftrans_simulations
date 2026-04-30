#!/usr/bin/env python3
"""
Analyze and visualize simulation design
Generate summary statistics and plots of parameter space coverage
"""

import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)


def load_design(design_file: str) -> pd.DataFrame:
    """Load simulation design from JSON and convert to DataFrame"""
    with open(design_file, 'r') as f:
        design = json.load(f)
    
    df = pd.DataFrame(design['parameters'])
    
    # Add categorical labels for better plotting
    df['rec_category'] = pd.cut(
        df['rec_rate'], 
        bins=[0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1],
        labels=['very_low', 'low', 'medium', 'high', 'very_high', 'ultra_high']
    )
    
    df['bottle_category'] = df['bottleneck_size'].map({
        1: 'tight', 
        5: 'medium', 
        20: 'loose'
    })
    
    df['est_category'] = df['est'].map({
        0.0: 'none',
        0.1: 'very_low',
        0.5: 'low_mod',
        1.0: 'moderate',
        2.0: 'high'
    })
    
    return df


def print_summary(df: pd.DataFrame):
    """Print comprehensive summary of simulation design"""
    
    print("\n" + "="*70)
    print("SIMULATION DESIGN ANALYSIS")
    print("="*70)
    
    print(f"\nTotal scenarios: {len(df)}")
    print(f"Unique parameter combinations: {len(df.drop('replicate_id', axis=1).drop_duplicates())}")
    print(f"Replicates per combination: {df['replicate_id'].max()}")
    
    print("\n" + "-"*70)
    print("PARAMETER DISTRIBUTIONS")
    print("-"*70)
    
    # Recombination rates
    print("\n1. Recombination Rate:")
    rec_counts = df['rec_category'].value_counts().sort_index()
    for cat, count in rec_counts.items():
        pct = count / len(df) * 100
        print(f"   {cat:12s}: {count:5d} ({pct:5.1f}%)")
    
    # Bottleneck sizes
    print("\n2. Bottleneck Size:")
    bottle_counts = df['bottle_category'].value_counts().sort_index()
    for cat, count in bottle_counts.items():
        pct = count / len(df) * 100
        print(f"   {cat:12s}: {count:5d} ({pct:5.1f}%)")
    
    # EST
    print("\n3. Expected Substitutions per Transmission:")
    est_counts = df['est_category'].value_counts().sort_index()
    for cat, count in est_counts.items():
        pct = count / len(df) * 100
        est_val = df[df['est_category']==cat]['est'].iloc[0]
        print(f"   {cat:12s} (EST={est_val}): {count:5d} ({pct:5.1f}%)")
    
    # Sampling proportions
    print("\n4. Sampling Proportion:")
    samp_counts = df['sampling_proportion'].value_counts().sort_index()
    for prop, count in samp_counts.items():
        pct = count / len(df) * 100
        print(f"   {prop:5.1f}: {count:5d} ({pct:5.1f}%)")
    
    # Outbreak sizes
    print("\n5. Outbreak Size:")
    outbreak_counts = df['outbreak_size'].value_counts().sort_index()
    for size, count in outbreak_counts.items():
        pct = count / len(df) * 100
        print(f"   {size:5d}: {count:5d} ({pct:5.1f}%)")
    
    # Population structure
    if 'n_populations' in df.columns and df['n_populations'].nunique() > 1:
        print("\n6. Number of Populations:")
        npop_counts = df['n_populations'].value_counts().sort_index()
        for npop, count in npop_counts.items():
            pct = count / len(df) * 100
            print(f"   {npop:5d}: {count:5d} ({pct:5.1f}%)")
        
        print("\n7. Migration Rate:")
        mig_counts = df['migration_rate'].value_counts().sort_index()
        for rate, count in mig_counts.items():
            pct = count / len(df) * 100
            print(f"   {rate:8.4f}: {count:5d} ({pct:5.1f}%)")
    
    print("\n" + "="*70)


def plot_parameter_space(df: pd.DataFrame, output_dir: Path):
    """Generate visualization of parameter space coverage"""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Core triad heatmap
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Recombination vs Bottleneck
    pivot1 = df.groupby(['rec_category', 'bottle_category']).size().unstack(fill_value=0)
    sns.heatmap(pivot1, annot=True, fmt='d', cmap='YlOrRd', ax=axes[0], cbar_kws={'label': 'Count'})
    axes[0].set_title('Recombination × Bottleneck')
    axes[0].set_xlabel('Bottleneck Size')
    axes[0].set_ylabel('Recombination Rate')
    
    # Recombination vs EST
    pivot2 = df.groupby(['rec_category', 'est_category']).size().unstack(fill_value=0)
    sns.heatmap(pivot2, annot=True, fmt='d', cmap='YlGnBu', ax=axes[1], cbar_kws={'label': 'Count'})
    axes[1].set_title('Recombination × EST')
    axes[1].set_xlabel('EST Level')
    axes[1].set_ylabel('Recombination Rate')
    
    # Bottleneck vs EST
    pivot3 = df.groupby(['bottle_category', 'est_category']).size().unstack(fill_value=0)
    sns.heatmap(pivot3, annot=True, fmt='d', cmap='BuPu', ax=axes[2], cbar_kws={'label': 'Count'})
    axes[2].set_title('Bottleneck × EST')
    axes[2].set_xlabel('EST Level')
    axes[2].set_ylabel('Bottleneck Size')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'parameter_space_triad.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'parameter_space_triad.png'}")
    
    # 2. Sampling and outbreak size
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Distribution by sampling proportion
    samp_data = df.groupby(['sampling_proportion', 'outbreak_size']).size().unstack(fill_value=0)
    samp_data.plot(kind='bar', ax=axes[0], color=['#1f77b4', '#ff7f0e', '#2ca02c'])
    axes[0].set_title('Scenarios by Sampling Proportion')
    axes[0].set_xlabel('Sampling Proportion')
    axes[0].set_ylabel('Number of Scenarios')
    axes[0].legend(title='Outbreak Size', loc='upper right')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=0)
    
    # Distribution by outbreak size
    outbreak_data = df.groupby(['outbreak_size', 'sampling_proportion']).size().unstack(fill_value=0)
    outbreak_data.plot(kind='bar', ax=axes[1], color=['#d62728', '#9467bd', '#8c564b', '#e377c2'])
    axes[1].set_title('Scenarios by Outbreak Size')
    axes[1].set_xlabel('Outbreak Size')
    axes[1].set_ylabel('Number of Scenarios')
    axes[1].legend(title='Sampling Prop', loc='upper right')
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'sampling_outbreak_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'sampling_outbreak_distribution.png'}")
    
    # 3. Migration scenarios (if present)
    if 'migration_rate' in df.columns and df['migration_rate'].nunique() > 1:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Population number distribution
        npop_data = df.groupby('n_populations').size()
        npop_data.plot(kind='bar', ax=axes[0], color='steelblue')
        axes[0].set_title('Scenarios by Number of Populations')
        axes[0].set_xlabel('Number of Populations')
        axes[0].set_ylabel('Number of Scenarios')
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=0)
        
        # Migration rate distribution
        mig_data = df.groupby('migration_rate').size()
        mig_data.plot(kind='bar', ax=axes[1], color='coral')
        axes[1].set_title('Scenarios by Migration Rate')
        axes[1].set_xlabel('Migration Rate')
        axes[1].set_ylabel('Number of Scenarios')
        axes[1].set_xticklabels([f'{x:.4f}' for x in mig_data.index], rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'migration_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {output_dir / 'migration_distribution.png'}")
    
    # 4. Complete parameter combination overview
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create a comprehensive combination identifier
    df['combination'] = (df['rec_category'].astype(str) + '_' + 
                        df['bottle_category'].astype(str) + '_' + 
                        df['est_category'].astype(str))
    
    combo_counts = df['combination'].value_counts().head(30)
    combo_counts.plot(kind='barh', ax=ax, color='teal')
    ax.set_title('Top 30 Parameter Combinations (by frequency)')
    ax.set_xlabel('Number of Scenarios')
    ax.set_ylabel('Parameter Combination')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'parameter_combinations.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'parameter_combinations.png'}")


def generate_batch_summary(df: pd.DataFrame, batch_size: int = 100) -> pd.DataFrame:
    """
    Generate summary for batch processing (useful for SLURM array jobs)
    """
    n_scenarios = len(df)
    n_batches = (n_scenarios + batch_size - 1) // batch_size
    
    batch_summary = []
    for batch_id in range(n_batches):
        start_idx = batch_id * batch_size
        end_idx = min((batch_id + 1) * batch_size, n_scenarios)
        
        batch_df = df.iloc[start_idx:end_idx]
        
        batch_summary.append({
            'batch_id': batch_id,
            'start_task': start_idx,
            'end_task': end_idx - 1,
            'n_tasks': end_idx - start_idx,
            'unique_scenarios': len(batch_df['scenario_id'].unique()),
            'replicates': batch_df['replicate_id'].nunique()
        })
    
    summary_df = pd.DataFrame(batch_summary)
    return summary_df


def export_scenario_list(df: pd.DataFrame, output_file: Path):
    """Export simplified scenario list for reference"""
    
    # Get unique scenarios (excluding replicates)
    unique_scenarios = df.drop('replicate_id', axis=1).drop_duplicates()
    
    # Select key columns
    columns = ['scenario_id', 'rec_rate', 'bottleneck_size', 'est', 
               'sampling_proportion', 'outbreak_size', 'n_populations', 'migration_rate']
    
    export_df = unique_scenarios[[col for col in columns if col in unique_scenarios.columns]]
    export_df = export_df.sort_values('scenario_id')
    
    export_df.to_csv(output_file, index=False)
    print(f"\nExported scenario list to: {output_file}")
    print(f"Unique scenarios: {len(export_df)}")


def main():
    parser = argparse.ArgumentParser(description="Analyze simulation design")
    parser.add_argument("design_file", type=str, help="Path to simulation_design.json")
    parser.add_argument("--output-dir", type=str, default="design_analysis",
                       help="Output directory for plots and summaries")
    parser.add_argument("--batch-size", type=int, default=100,
                       help="Batch size for SLURM array jobs")
    parser.add_argument("--no-plots", action="store_true",
                       help="Skip generating plots")
    
    args = parser.parse_args()
    
    # Load design
    print(f"Loading design from: {args.design_file}")
    df = load_design(args.design_file)
    
    # Print summary
    print_summary(df)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    if not args.no_plots:
        print("\nGenerating visualizations...")
        try:
            plot_parameter_space(df, output_dir)
        except Exception as e:
            print(f"Warning: Could not generate plots: {e}")
    
    # Generate batch summary
    print("\nGenerating batch processing summary...")
    batch_summary = generate_batch_summary(df, args.batch_size)
    batch_file = output_dir / "batch_summary.csv"
    batch_summary.to_csv(batch_file, index=False)
    print(f"Saved batch summary to: {batch_file}")
    print(f"Total batches: {len(batch_summary)}")
    
    # Export scenario list
    scenario_file = output_dir / "scenario_list.csv"
    export_scenario_list(df, scenario_file)
    
    # Save full dataframe
    full_file = output_dir / "full_design.csv"
    df.to_csv(full_file, index=False)
    print(f"\nSaved full design to: {full_file}")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nResults saved to: {output_dir}/")
    print("\nGenerated files:")
    print(f"  - {batch_file.name}")
    print(f"  - {scenario_file.name}")
    print(f"  - {full_file.name}")
    if not args.no_plots:
        print(f"  - parameter_space_triad.png")
        print(f"  - sampling_outbreak_distribution.png")
        print(f"  - parameter_combinations.png")


if __name__ == "__main__":
    main()
