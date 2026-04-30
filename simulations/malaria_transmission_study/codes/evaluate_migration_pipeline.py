#!/usr/bin/env python3

"""
Evaluation script for migration scenarios with tree-sequence ground truth
Adapted for: true_ibd from tree sequences + generation distance thresholds
"""

import numpy as np
import pandas as pd
# from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    confusion_matrix,
    adjusted_rand_score,
    normalized_mutual_info_score
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef
)
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob

# ============================================================================
# CONFIGURATION
# ============================================================================

# Generation thresholds to test
G_THRESHOLDS = [5, 10, 15, 25]  # Test multiple definitions

# Identifiability criteria
AUROC_THRESHOLD = 0.80
SENSITIVITY_AT_90SPEC_THRESHOLD = 0.60

# Folder pattern
SCENARIO_PATTERN = "MIG*_*_*_*_*"  # Matches your naming

# ============================================================================
#                        1. CHECK FILE EXTENSION
# ============================================================================
def check_extension_pathlib(filename):

    """
    Checks if the file has the target extension using pathlib.
    The check is case-insensitive.
    , target_extension
    """
    # Get the suffix (extension with the dot) and convert to lowercase
    file_extension = Path(filename).suffix.lower()
    # Ensure the target extension is also in lowercase and includes the dot
    # if not target_extension.startswith('.'):
    #    target_extension = '.' + target_extension

    return file_extension

# ============================================================================
#           1. LOAD GROUND TRUTH (From Tree Sequence)
# ============================================================================

def load_ground_truth(scenario_dir):
    """
    Load ground truth from tree sequence analysis

    Expected file: ground_truth.csv with columns:
    - id1, id2: sample pair
    - gen_distance: generations to MRCA
    - true_ibd_prop: proportion of genome sharing IBD (from tree sequence)
    - population_i, population_j: spatial info (for RQ3)
    """
    truth_file = scenario_dir / "true_ibd_summary.tsv"  # Or your actual filename

    if not truth_file.exists():
        return None

    # Get file EXTENSION
    file_extension = check_extension_pathlib(truth_file)

    if file_extension == ".csv"
        df = pd.read_csv(truth_file)
    elif file_extension == ".tsv"
        # Read the TSV file by specifying the separator as a tab character ('\t')
        df = pd.read_csv(truth_file, sep='\t')
    else:
        print(f"'{truth_file}' does not have an allowed extension.")

    # ALTERNATIVE
#    allowed_extensions = ('.csv', '.tsv', '.xlsx')
#    if filename.lower().endswith(allowed_extensions):
#        print(f"'{filename}' has an allowed extension.")
#    else:
#        print(f"'{filename}' does not have an allowed extension.")


    # Ensure proper columns
    required = ['Id1', 'Id2', 'min_tmrca', 'total_ibd_prop']
    if not all(col in df.columns for col in required):
        raise ValueError(f"Ground truth missing required columns: {required}")

    # Create binary labels for each threshold
    for G in G_THRESHOLDS:
        df[f'related_G{G}'] = (df['min_tmrca'] <= G).astype(int)

    return df


# ============================================================================
# 2. LOAD INFERENCE RESULTS
# ============================================================================

def load_inference_results(scenario_dir):
    """
    Load IBD, IBS, and phylogenetic distance predictions

    Expected files:
    - ibd_results.csv: id1, id2, ibd_proportion (or ibd_score)
    - ibs_results.csv: id1, id2, ibs_similarity
    - phylo_results.csv: id1, id2, patristic_distance
    """
    results = {}

    # IBD
    ibd_file = scenario_dir / "inferred_ibd_hmm.tsv"  # Adjust to your filename
    if ibd_file.exists():
        results['IBD'] = pd.read_csv(ibd_file, sep='\t')

    # IBS
    ibs_file = scenario_dir / "inferred_ibs.tsv"
    if ibs_file.exists():
        results['IBS'] = pd.read_csv(ibs_file, sep='\t')

    # Phylogenetic
    phylo_file = scenario_dir / "phylo_results.csv"
    if phylo_file.exists():
        phylo = pd.read_csv(phylo_file)
        # Convert distance to similarity (higher = more related)
        phylo['similarity'] = 1 / (1 + phylo['patristic_distance'])
        results['Phylo'] = phylo

    return results


# ============================================================================
# 3. CALCULATE SENSITIVITY AT SPECIFICITY
# ============================================================================

def calculate_sensitivity_at_specificity(y_true, y_score, target_spec=0.90):
    """
    Calculate sensitivity at a given specificity threshold
    """
    fpr, tpr, thresholds = roc_curve(y_true, y_score)

    # Specificity = 1 - FPR
    specificity = 1 - fpr

    # Find closest point to target specificity
    idx = np.argmin(np.abs(specificity - target_spec))

    sensitivity = tpr[idx]
    threshold = thresholds[idx] if idx < len(thresholds) else np.nan

    return sensitivity, threshold


# ============================================================================
# 4. EVALUATE SINGLE SCENARIO
# ============================================================================

def evaluate_scenario(scenario_dir, scenario_id):
    """
    Complete evaluation for one scenario
    """
    # Load data
    truth = load_ground_truth(scenario_dir)
    if truth is None:
        print(f"  ⚠️  No ground truth for {scenario_id}")
        return None

    inference = load_inference_results(scenario_dir)
    if not inference:
        print(f"  ⚠️  No inference results for {scenario_id}")
        return None

    # Extract scenario parameters from folder name
    # Pattern: MIG001_verylow_tight_moderate_mlow
    # Meaning: runid, rec_rate, bottleneck, est, migration
    params = parse_scenario_id(scenario_id)

    results = []

    # Evaluate each method × threshold combination
    for G in G_THRESHOLDS:
        y_true = truth[f'related_G{G}'].values

        for method_name, method_df in inference.items():
            # Merge with ground truth
            merged = truth.merge(
                method_df[['Id1', 'Id2', method_df.columns[-1]]],  # Last column is the score
                on=['Id1', 'Id2'],
                how='inner'
            )

            if len(merged) == 0:
                continue

            # Get predictions
            score_col = merged.columns[-1]  # Last column after merge
            y_score = merged[score_col].values
            y_true_merged = merged[f'related_G{G}'].values

            # Calculate metrics
            try:
                auroc = roc_auc_score(y_true_merged, y_score)
                sens_at_90spec, threshold_90 = calculate_sensitivity_at_specificity(
                    y_true_merged, y_score, target_spec=0.90
                )

                # Identifiable?
                identifiable = (auroc >= AUROC_THRESHOLD and
                               sens_at_90spec >= SENSITIVITY_AT_90SPEC_THRESHOLD)

                # Store result
                results.append({
                    'scenario_id': scenario_id,
                    **params,
                    'G_threshold': G,
                    'method': method_name,
                    'n_pairs': len(merged),
                    'n_positive': sum(y_true_merged),
                    'auroc': auroc,
                    'sensitivity_at_90spec': sens_at_90spec,
                    'threshold_at_90spec': threshold_90,
                    'identifiable': identifiable
                })

            except Exception as e:
                print(f"  ⚠️  Error evaluating {method_name} at G={G}: {e}")
                continue

    return pd.DataFrame(results)


# ============================================================================
# 5. PARSE SCENARIO PARAMETERS
# ============================================================================

def parse_scenario_id(scenario_id):
    """
    Extract parameters from folder name
    Pattern: MIG001_verylow_tight_moderate_mlow

    Returns dict with:
    - run_id: 1
    - rec_rate_label: 'verylow'
    - rec_rate: 1e-9 (actual value)
    - bottleneck_label: 'tight'
    - bottleneck: 1 (actual value)
    - est_label: 'moderate'
    - est: 0.3 (actual value)
    - migration_label: 'mlow'
    - migration: 0.001 (actual value)
    """
    parts = scenario_id.replace('MIG', '').split('_')

    # Mapping labels to values (adjust based on your actual encoding)
    rec_map = {'verylow': 1e-9, 'low': 1e-8, 'medium': 1e-7, 'high': 1e-6}
    bottleneck_map = {'tight': 1, 'medium': 5, 'loose': 20}
    est_map = {'low': 0.145, 'moderate': 0.3}
    migration_map = {'mlow': 0.001, 'mmedium': 0.01, 'mhigh': 0.05}

    return {
        'run_id': int(parts[0]),
        'rec_rate_label': parts[1],
        'rec_rate': rec_map.get(parts[1], np.nan),
        'bottleneck_label': parts[2],
        'bottleneck': bottleneck_map.get(parts[2], np.nan),
        'est_label': parts[3],
        'est': est_map.get(parts[3], np.nan),
        'migration_label': parts[4],
        'migration': migration_map.get(parts[4], np.nan)
    }


def classify_identifiability(auc, sensitivity_at_90spec):
    """
    Classify scenario as identifiable or not based on thresholds

    Returns:
    --------
    category : str
        'High confidence', 'Acceptable', 'Marginal', or 'Not identifiable'
    """
    if auc >= 0.90 and sensitivity_at_90spec >= 0.70:
        return "IDENTIFIABLE - High confidence"
    elif auc >= 0.80 and sensitivity_at_90spec >= 0.60:
        return "IDENTIFIABLE - Acceptable"
    elif auc >= 0.70:
        return "MARGINAL - Context dependent"
    else:
        return "NOT IDENTIFIABLE"


# ============================================================================
# 6. PARAMETER EFFECT ANALYSIS (RQ2)
# ============================================================================

def analyze_parameter_effects(results_df):
    """
    Fit regression model to quantify parameter effects on AUC

    Returns:
    --------
    effects : pd.DataFrame
        Standardized coefficients for each method
    """
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler

    # Prepare predictors (log-transform rates)
    X = results_df[['rec_rate', 'bottleneck', 'est', 'migration', 'sampling']].copy()
    X['rec_rate'] = np.log10(X['rec_rate'])
    X['migration'] = np.log10(X['migration'])

    # Standardize for effect size comparison
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    effects = []
    for method in ['ibd', 'ibs', 'phylo']:
        # Subset to method
        method_data = results_df[results_df['method'] == method]
        y = method_data['auc'].values

        # Fit model
        model = LinearRegression()
        model.fit(X_scaled, y)

        # Store coefficients
        for param, coef in zip(X.columns, model.coef_):
            effects.append({
                'method': method,
                'parameter': param,
                'coefficient': coef,
                'effect_size': abs(coef)
            })

    effects_df = pd.DataFrame(effects)
    return effects_df


# ============================================================================
# 7. VISUALIZATION
# ============================================================================

def plot_roc_curves(scenario_dir, methods=['ibd', 'ibs', 'phylo']):
    """
    Plot ROC curves for all methods on same scenario
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    for method in methods:
        # Load data
        truth = load_ground_truth(scenario_dir)
        inference = load_inference_results(scenario_dir, method)
        data = merge_data(truth, inference, method)

        y_true = data['is_direct_transmission'].values
        y_score = data['distance_metric'].values

        # Invert phylogenetic distance
        if method == 'phylo':
            y_score = -1 * y_score

        # Calculate ROC
        fpr, tpr, _ = roc_curve(y_true, y_score)
        auc = roc_auc_score(y_true, y_score)

        # Plot
        ax.plot(fpr, tpr, label=f'{method.upper()} (AUC={auc:.2f})', linewidth=2)

    # Format
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='Random')
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title(f'ROC Curves - {scenario_dir}', fontsize=14)
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    return fig


def plot_performance_heatmap(results_df, metric='auc'):
    """
    Heatmap of method performance across recombination rates
    """
    # Pivot for heatmap
    pivot = results_df.pivot_table(
        index='method',
        columns='rec_rate',
        values=metric,
        aggfunc='mean'
    )

    # Plot
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdYlGn',
                vmin=0.5, vmax=1.0, ax=ax, cbar_kws={'label': metric.upper()})
    ax.set_xlabel('Recombination Rate', fontsize=12)
    ax.set_ylabel('Method', fontsize=12)
    ax.set_title(f'Method Performance Heatmap ({metric.upper()})', fontsize=14)

    plt.tight_layout()
    return fig


def plot_identifiability_phase_diagram(results_df):
    """
    Phase diagram showing identifiable parameter space
    """
    # Get best AUC across methods for each parameter combination
    best_auc = results_df.groupby(['rec_rate', 'sampling'])['auc'].max().reset_index()

    # Pivot for contour plot
    pivot = best_auc.pivot(index='sampling', columns='rec_rate', values='auc')

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    X = np.log10(pivot.columns.values)  # Log recombination rate
    Y = pivot.index.values  # Sampling proportion
    Z = pivot.values

    # Contour plot
    contour = ax.contourf(X, Y, Z, levels=np.linspace(0.5, 1.0, 11),
                          cmap='RdYlGn', extend='both')

    # Add identifiability threshold line
    ax.contour(X, Y, Z, levels=[0.80], colors='black', linewidths=2,
              linestyles='dashed')

    # Format
    ax.set_xlabel('Recombination Rate (log10)', fontsize=12)
    ax.set_ylabel('Sampling Proportion', fontsize=12)
    ax.set_title('Identifiability Phase Diagram\n(Dashed line = AUC 0.80 threshold)',
                fontsize=14)

    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Best AUC', fontsize=12)

    plt.tight_layout()
    return fig


# ============================================================================
# 8. SUMMARY STATISTICS
# ============================================================================

def summarize_identifiability(results_df):
    """
    Generate summary statistics for identifiability
    """
    summary = {
        'overall': {},
        'by_method': {},
        'by_recombination': {},
        'by_sampling': {}
    }

    # Overall
    summary['overall'] = {
        'total_scenarios': len(results_df['scenario'].unique()),
        'pct_identifiable': (results_df['auc'] >= 0.80).mean() * 100,
        'mean_auc': results_df['auc'].mean(),
        'median_auc': results_df['auc'].median()
    }

    # By method
    for method in results_df['method'].unique():
        method_data = results_df[results_df['method'] == method]
        summary['by_method'][method] = {
            'pct_identifiable': (method_data['auc'] >= 0.80).mean() * 100,
            'mean_auc': method_data['auc'].mean()
        }

    # By recombination rate
    for rec_rate in results_df['rec_rate'].unique():
        rec_data = results_df[results_df['rec_rate'] == rec_rate]
        summary['by_recombination'][rec_rate] = {
            'pct_identifiable': (rec_data['auc'] >= 0.80).mean() * 100,
            'mean_auc': rec_data['auc'].mean(),
            'best_method': rec_data.loc[rec_data['auc'].idxmax(), 'method']
        }

    # By sampling proportion
    for samp in results_df['sampling'].unique():
        samp_data = results_df[results_df['sampling'] == samp]
        summary['by_sampling'][samp] = {
            'pct_identifiable': (samp_data['auc'] >= 0.80).mean() * 100,
            'mean_auc': samp_data['auc'].mean()
        }

    return summary

# ============================================================================
# 6. SANITY CHECK ON PILOT SCENARIOS
# ============================================================================

def run_sanity_check(base_dir):
    """
    Test on 3 representative scenarios before full run
    """
    print("\n" + "="*80)
    print("SANITY CHECK: Testing 3 pilot scenarios")
    print("="*80)

    # Select 3 scenarios (you adjust based on your actual folders)
    pilot_scenarios = [
        'MIG001_verylow_tight_moderate_mlow',  # Low recomb, low migration
        'MIG841_high_medium_low_mhigh',          # High recomb, high migration
        'MIG601_medium_medium_moderate_mmedium' # Moderate everything
    ]

    for scenario_id in pilot_scenarios:
        scenario_dir = Path(base_dir) / scenario_id

        if not scenario_dir.exists():
            print(f"\n⚠️  Scenario not found: {scenario_id}")
            continue

        print(f"\n{'='*80}")
        print(f"Evaluating: {scenario_id}")
        print(f"{'='*80}")

        # results_df = evaluate_scenario(scenario_dir, scenario_id)
        results_df = evaluate_scenario(base_dir, scenario_id)

        if results_df is not None:
            print("\nResults summary:")
            print(results_df.groupby(['method', 'G_threshold'])['auroc'].agg(['mean', 'count']))

            # Check hypothesis for low recombination scenario
            if 'verylow' in scenario_id:
                print("\n✓ Expected: IBS > IBD for low recombination")
                ibs_auc = results_df[results_df['method'] == 'IBS']['auroc'].mean()
                ibd_auc = results_df[results_df['method'] == 'IBD']['auroc'].mean()
                print(f"  IBS AUC: {ibs_auc:.3f}")
                print(f"  IBD AUC: {ibd_auc:.3f}")
                if ibs_auc > ibd_auc:
                    print("  ✓ PASS: IBS outperforms IBD")
                else:
                    print("  ✗ FAIL: IBD outperforms IBS (unexpected!)")

            # Check for high recombination
            if 'high' in scenario_id and 'rec' in scenario_id:
                print("\n✓ Expected: IBD > IBS for high recombination")
                ibs_auc = results_df[results_df['method'] == 'IBS']['auroc'].mean()
                ibd_auc = results_df[results_df['method'] == 'IBD']['auroc'].mean()
                print(f"  IBD AUC: {ibd_auc:.3f}")
                print(f"  IBS AUC: {ibs_auc:.3f}")
                if ibd_auc > ibs_auc:
                    print("  ✓ PASS: IBD outperforms IBS")
                else:
                    print("  ✗ FAIL: IBS outperforms IBD (unexpected!)")


# ============================================================================
# 7. FULL ANALYSIS (All 900 scenarios)
# ============================================================================

def run_full_analysis(base_dir, output_file='migration_results.csv'):
    """
    Analyze all scenarios and compile results
    """
    print("\n" + "="*80)
    print("FULL ANALYSIS: Processing all scenarios")
    print("="*80)

    scenario_dirs = sorted(Path(base_dir).glob(SCENARIO_PATTERN))
    print(f"\nFound {len(scenario_dirs)} scenarios")

    all_results = []

    for i, scenario_dir in enumerate(scenario_dirs, 1):
        scenario_id = scenario_dir.name

        if i % 50 == 0:
            print(f"  Progress: {i}/{len(scenario_dirs)} scenarios processed")

        results_df = evaluate_scenario(scenario_dir, scenario_id)

        if results_df is not None:
            all_results.append(results_df)

    # Combine all results
    final_df = pd.concat(all_results, ignore_index=True)

    # Save
    final_df.to_csv(output_file, index=False)
    print(f"\n✓ Results saved to: {output_file}")

    return final_df


# ============================================================================
# 8. ANSWER RESEARCH QUESTIONS
# ============================================================================

def answer_research_questions(results_df):
    """
    Generate summaries addressing each RQ
    """
    print("\n" + "="*80)
    print("ANSWERING RESEARCH QUESTIONS")
    print("="*80)

    # RQ1: Identifiability conditions
    print("\n" + "-"*80)
    print("RQ1: Under what conditions is inference identifiable?")
    print("-"*80)

    # Group by recombination rate
    rq1_summary = results_df.groupby(['rec_rate', 'method', 'G_threshold']).agg({
        'identifiable': 'mean',
        'auroc': 'mean'
    }).round(3)

    print("\n% Identifiable by recombination rate:")
    print(rq1_summary)

    # RQ2: Effect of biological processes (requires regression - see separate function)
    print("\n" + "-"*80)
    print("RQ2: Which biological processes most erode identifiability?")
    print("-"*80)
    print("  (Run variance decomposition separately)")

    # RQ3: Migration effects (if population info available)
    print("\n" + "-"*80)
    print("RQ3: Migration rate effects on importation detection")
    print("-"*80)

    migration_summary = results_df.groupby(['migration', 'method']).agg({
        'auroc': 'mean'
    }).round(3)
    print(migration_summary)

    # RQ4: Failure taxonomy
    print("\n" + "-"*80)
    print("RQ4: Methodological vs. fundamental failure")
    print("-"*80)

    # Find scenarios where at least one method fails
    failed = results_df[results_df['auroc'] < 0.70]

    # Calculate AUC range per scenario
    auc_range = results_df.groupby(['scenario_id', 'G_threshold']).agg({
        'auroc': lambda x: x.max() - x.min()
    }).reset_index()

    auc_range.columns = ['scenario_id', 'G_threshold', 'auc_range']

    # Classify failures
    auc_range['failure_type'] = 'SUCCESS'
    auc_range.loc[auc_range['auc_range'] < 0.10, 'failure_type'] = 'FUNDAMENTAL'
    auc_range.loc[auc_range['auc_range'] >= 0.15, 'failure_type'] = 'METHODOLOGICAL'


    if auc_max < 0.70:
        if auc_range < 0.10:
            failure_type = "FUNDAMENTAL (all methods fail equally)"
        else:
            failure_type = "MIXED (some methods work marginally)"
    elif auc_min < 0.70:
        if auc_range >= 0.15:
            failure_type = "METHODOLOGICAL (method choice critical)"
        else:
            failure_type = "MIXED (methods differ moderately)"
    else:
        failure_type = "SUCCESS (all methods work)"

    print("\nFailure type distribution:")
    print(auc_range['failure_type'].value_counts())


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Evaluate migration scenarios")
    parser.add_argument("--base-dir", required=True, help="Directory containing MIG* folders")
    parser.add_argument("--sanity-check", action="store_true", help="Run sanity check only")
    parser.add_argument("--full", action="store_true", help="Run full analysis on all 900")
    parser.add_argument("--output", default="migration_results.csv", help="Output file")

    args = parser.parse_args()

    if args.sanity_check:
        run_sanity_check(args.base_dir)

    if args.full:
        results = run_full_analysis(args.base_dir, args.output)
        answer_research_questions(results)

        # results_df.to_csv("identifiability_results.csv", index=False)

        # 4. Parameter effects
        effects = analyze_parameter_effects(results)
        print(effects.sort_values('effect_size', ascending=False))

        # 5. Generate plots
        # fig = plot_roc_curves("sim_migration/MIG001")
        # fig.savefig("roc_curves_MIG001.png", dpi=300)

        # 6. Summary
        summary = summarize_identifiability(results)
        print(json.dumps(summary, indent=2))

    if not args.sanity_check and not args.full:
        print("Please specify --sanity-check or --full")
        parser.print_help()


# ==================================
#               USAGE
# ==================================
# 1. Sanity Check First (2-3 scenarios) This will test 3 representative scenarios and verify your hypothesis.
# python evaluate_migration_pipeline.py \
#   --base-dir /path/to/sim_migration \
#   --sanity-check
#
# # 2. Full Analysis (900 scenarios)
# python evaluate_migration_pipeline.py \
#   --base-dir /path/to/sim_migration \
#   --full \
#   --output migration_all_results.csv
