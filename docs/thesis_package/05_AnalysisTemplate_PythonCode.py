#!/usr/bin/env python3
"""
Template for calculating identifiability metrics
For PhD Objective 1: Transmission inference evaluation
"""

import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_auc_score, 
    roc_curve, 
    confusion_matrix,
    adjusted_rand_score,
    normalized_mutual_info_score
)
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================================
# 1. LOAD DATA
# ============================================================================

def load_ground_truth(scenario_dir):
    """
    Load ground truth transmission relationships from simulation
    
    Expected format:
    sample_i, sample_j, is_direct_transmission, transmission_distance, 
    same_population, generation_gap
    """
    truth = pd.read_csv(f"{scenario_dir}/ground_truth_pairs.csv")
    return truth


def load_inference_results(scenario_dir, method='ibd'):
    """
    Load inference method results
    
    Expected format:
    sample_i, sample_j, distance_metric
    
    For IBD: higher values = more related (0-1)
    For IBS: higher values = more similar (0-1)
    For phylo: lower values = more related (patristic distance)
    """
    results = pd.read_csv(f"{scenario_dir}/{method}_results.csv")
    return results


def merge_data(ground_truth, inference_results, method='ibd'):
    """
    Merge ground truth with inference results
    Ensure same pairs in both datasets
    """
    merged = pd.merge(
        ground_truth,
        inference_results,
        on=['sample_i', 'sample_j'],
        how='inner'
    )
    
    # Verify no missing pairs
    n_expected = len(ground_truth)
    n_merged = len(merged)
    if n_merged < n_expected:
        print(f"WARNING: {n_expected - n_merged} pairs missing from {method} results")
    
    return merged


# ============================================================================
# 2. PRIMARY METRIC: ROC-AUC
# ============================================================================

def calculate_auc(y_true, y_score, higher_is_related=True):
    """
    Calculate ROC-AUC for transmission inference
    
    Parameters:
    -----------
    y_true : array-like
        Binary labels (1 = direct transmission, 0 = unrelated)
    y_score : array-like
        Distance/similarity scores from inference method
    higher_is_related : bool
        True if high scores = more related (IBD, IBS)
        False if low scores = more related (phylogenetic distance)
    
    Returns:
    --------
    auc : float
        Area under ROC curve
    """
    # Invert scores if lower values indicate relatedness
    if not higher_is_related:
        y_score = -1 * y_score
    
    auc = roc_auc_score(y_true, y_score)
    return auc


def calculate_sensitivity_at_specificity(y_true, y_score, target_spec=0.90, 
                                         higher_is_related=True):
    """
    Calculate sensitivity at a given specificity threshold
    
    Returns:
    --------
    sensitivity : float
        True positive rate at target specificity
    threshold : float
        Distance threshold achieving target specificity
    """
    if not higher_is_related:
        y_score = -1 * y_score
    
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    
    # Find closest point to target specificity (spec = 1 - fpr)
    target_fpr = 1 - target_spec
    idx = np.argmin(np.abs(fpr - target_fpr))
    
    sensitivity = tpr[idx]
    threshold = thresholds[idx]
    
    return sensitivity, threshold


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
# 3. EVALUATE SINGLE SCENARIO
# ============================================================================

def evaluate_scenario(scenario_dir, method='ibd', higher_is_related=True):
    """
    Complete evaluation pipeline for one scenario
    
    Returns:
    --------
    results : dict
        Dictionary with all metrics
    """
    # Load data
    truth = load_ground_truth(scenario_dir)
    inference = load_inference_results(scenario_dir, method)
    data = merge_data(truth, inference, method)
    
    # Extract labels and scores
    y_true = data['is_direct_transmission'].values
    y_score = data['distance_metric'].values
    
    # Calculate metrics
    auc = calculate_auc(y_true, y_score, higher_is_related)
    sens, thresh = calculate_sensitivity_at_specificity(
        y_true, y_score, target_spec=0.90, higher_is_related=higher_is_related
    )
    
    # Classify identifiability
    identifiability = classify_identifiability(auc, sens)
    
    # Store results
    results = {
        'scenario': scenario_dir,
        'method': method,
        'n_pairs': len(data),
        'n_transmission_pairs': sum(y_true),
        'auc': auc,
        'sensitivity_at_90spec': sens,
        'threshold_at_90spec': thresh,
        'identifiable': identifiability
    }
    
    return results


# ============================================================================
# 4. METHOD COMPARISON (RQ4)
# ============================================================================

def compare_methods(scenario_dir):
    """
    Compare IBD, IBS, and phylogenetic distance on same scenario
    
    Returns:
    --------
    comparison : pd.DataFrame
        Results for all three methods
    failure_type : str
        'Fundamental', 'Methodological', or 'Mixed'
    """
    # Evaluate each method
    results_ibd = evaluate_scenario(scenario_dir, 'ibd', higher_is_related=True)
    results_ibs = evaluate_scenario(scenario_dir, 'ibs', higher_is_related=True)
    results_phylo = evaluate_scenario(scenario_dir, 'phylo', higher_is_related=False)
    
    comparison = pd.DataFrame([results_ibd, results_ibs, results_phylo])
    
    # Determine failure type
    aucs = comparison['auc'].values
    auc_max = aucs.max()
    auc_min = aucs.min()
    auc_range = auc_max - auc_min
    
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
    
    comparison['failure_type'] = failure_type
    
    return comparison, failure_type


# ============================================================================
# 5. BATCH ANALYSIS ACROSS SCENARIOS
# ============================================================================

def analyze_all_scenarios(scenario_list, methods=['ibd', 'ibs', 'phylo']):
    """
    Analyze all scenarios and compile results
    
    Parameters:
    -----------
    scenario_list : list
        List of scenario directories
    methods : list
        Methods to evaluate
    
    Returns:
    --------
    results_df : pd.DataFrame
        Compiled results across all scenarios
    """
    all_results = []
    
    for scenario_dir in scenario_list:
        # Extract scenario parameters from directory name or metadata
        params = extract_scenario_params(scenario_dir)
        
        for method in methods:
            # Evaluate scenario
            res = evaluate_scenario(scenario_dir, method, 
                                   higher_is_related=(method != 'phylo'))
            
            # Add parameter info
            res.update(params)
            all_results.append(res)
    
    results_df = pd.DataFrame(all_results)
    return results_df


def extract_scenario_params(scenario_dir):
    """
    Extract parameter values from scenario directory or metadata file
    
    Returns:
    --------
    params : dict
        Dictionary with parameter values
    """
    # Example: read from scenario metadata JSON
    import json
    meta_file = f"{scenario_dir}/scenario_metadata.json"
    
    with open(meta_file, 'r') as f:
        params = json.load(f)
    
    return {
        'rec_rate': params['rec_rate'],
        'bottleneck': params['bottleneck_size'],
        'est': params['est'],
        'migration': params['migration_rate'],
        'sampling': params['sampling_proportion']
    }


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
# 9. MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    # Example usage
    
    # 1. Evaluate single scenario
    print("Evaluating single scenario...")
    results = evaluate_scenario(
        scenario_dir="sim_migration/MIG001",
        method='ibd',
        higher_is_related=True
    )
    print(results)
    
    # 2. Compare methods
    print("\nComparing methods...")
    comparison, failure_type = compare_methods("sim_migration/MIG001")
    print(comparison)
    print(f"Failure type: {failure_type}")
    
    # 3. Batch analysis (example with placeholder list)
    # scenario_list = [f"sim_migration/MIG{i:03d}" for i in range(1, 46)]
    # results_df = analyze_all_scenarios(scenario_list)
    # results_df.to_csv("identifiability_results.csv", index=False)
    
    # 4. Parameter effects
    # effects = analyze_parameter_effects(results_df)
    # print(effects.sort_values('effect_size', ascending=False))
    
    # 5. Generate plots
    # fig = plot_roc_curves("sim_migration/MIG001")
    # fig.savefig("roc_curves_MIG001.png", dpi=300)
    
    # 6. Summary
    # summary = summarize_identifiability(results_df)
    # print(json.dumps(summary, indent=2))
