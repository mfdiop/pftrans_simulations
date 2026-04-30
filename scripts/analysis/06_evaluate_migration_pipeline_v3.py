#!/usr/bin/env python3

"""
Evaluation script for migration scenarios with tree-sequence ground truth
Adapted for: true_ibd from tree sequences + generation distance thresholds
"""

import numpy as np
import pandas as pd
import os
import json
# from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import (
    roc_auc_score, roc_curve,
    precision_recall_curve, auc,
    average_precision_score,
    confusion_matrix, f1_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef
)
from sklearn.utils import resample
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
import pyreadr

# ============================================================================
# CONFIGURATION
# ============================================================================

# Generation thresholds to test
G_THRESHOLDS = [1, 3, 5, 10, 15, 25]  # Test multiple definitions
N_FULL = 2000
REPS_PER_SAMPLE = 20
SAMPLING_PROPORTIONS = [0.05, 0.1, 0.2]
SAMPLE_SIZES = [100, 200, 400]

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

def extract_replicate_id(scenario_id):
    # MIG001_verylow_tight_moderate_mlow → 1
    return int(scenario_id[3:6])

def infer_sampling_proportion(scenario_id):
    rep_id = extract_replicate_id(scenario_id)
    sample_block = (rep_id - 1) // 20 + 1
    # L'astuce : (1-1)%3=0, (2-1)%3=1, (3-1)%3=2, (4-1)%3=0 ...
    list_index = (sample_block - 1) % len(SAMPLING_PROPORTIONS)
    return SAMPLING_PROPORTIONS[list_index]


def infer_sample_size(scenario_id):
    # p = infer_sampling_proportion(scenario_id)
    # return int(p * N_FULL)
    rep_id = extract_replicate_id(scenario_id)
    sample_block = (rep_id - 1) // 20 + 1
    # L'astuce : (1-1)%3=0, (2-1)%3=1, (3-1)%3=2, (4-1)%3=0 ...
    list_index = (sample_block - 1) % len(SAMPLE_SIZES)
    return SAMPLE_SIZES[list_index]

def precision_at_k(y_true, y_score, k):
    order = np.argsort(y_score)[::-1]
    top_k = order[:k]
    return np.mean(y_true[top_k])


def recall_at_fixed_precision(y_true, y_score, target_precision=0.8):
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    valid = precision >= target_precision
    if not np.any(valid):
        return 0.0
    return np.max(recall[valid])

# ============================================================================
#
# ============================================================================
def canonical_pair(a, b):
    """
    Canonical pair key (ensures id1 <= id2).
    Returns None if either value is missing.
    """
    if pd.isna(a) or pd.isna(b):
        return None

    a_str = str(a)
    b_str = str(b)

    if a_str <= b_str:
        return f"{a_str}--{b_str}"
    else:
        return f"{b_str}--{a_str}"

# ============================================================================
# Vectorized Version for DataFrames (Recommended for Large Data)
# ============================================================================
def canonical_pair_series(df, col1="Id1", col2="Id2"):
    mask = df[col1].notna() & df[col2].notna()

    a = df.loc[mask, col1].astype(str)
    b = df.loc[mask, col2].astype(str)

    key = np.where(a <= b, a + "--" + b, b + "--" + a)

    result = pd.Series(None, index=df.index)
    result.loc[mask] = key

    return result

# ==========================================
# Create New column to have Id1 <= Id2
# ==========================================
# # Ensure true_ibd is a pandas DataFrame
# true_ibd = pd.DataFrame(true_ibd)
#
# # Add pair_key column only if it does not already exist
# if "pair_key" not in true_ibd.columns:
#     mask = true_ibd["Id1"].notna() & true_ibd["Id2"].notna()
#
#     id1 = true_ibd.loc[mask, "Id1"].astype(str)
#     id2 = true_ibd.loc[mask, "Id2"].astype(str)
#
#     pair_key = np.where(
#         id1 <= id2,
#         id1 + "--" + id2,
#         id2 + "--" + id1
#     )
#
#     true_ibd["pair_key"] = None
#     true_ibd.loc[mask, "pair_key"] = pair_key

def add_pair_key(df, col1="Id1", col2="Id2"):
    if "pair_key" in df.columns:
        return df

    mask = df[col1].notna() & df[col2].notna()
    a = df.loc[mask, col1].astype(str)
    b = df.loc[mask, col2].astype(str)

    df["pair_key"] = None
    df.loc[mask, "pair_key"] = np.where(
        a <= b,
        a + "--" + b,
        b + "--" + a
    )
    return df

# ============================================================================
#           1. LOAD GROUND TRUTH (From Tree Sequence)
# ============================================================================

def load_ground_truth(scenario_dir, scenario_id):
    """
    Load ground truth from tree sequence analysis

    Expected file: ground_truth.csv with columns:
    - id1, id2: sample pair
    - gen_distance: generations to MRCA
    - true_ibd_prop: proportion of genome sharing IBD (from tree sequence)
    - population_i, population_j: spatial info (for RQ3)
    """
    print(f" Formatting true IBD in {scenario_dir}  ...\n")
    scenario_dir = Path(scenario_dir) / scenario_id
    truth_file = f"{scenario_dir}/true_ibd_summary.tsv"  # Or your actual filename

    if not os.path.exists(truth_file):
        print(f"{truth_file} does not exist ...")
        return None

    print(f"{truth_file} exists ...\n")

    # Get file EXTENSION
    file_extension = check_extension_pathlib(truth_file)

    if file_extension == ".csv":
        df = pd.read_csv(truth_file)
    elif file_extension == ".tsv":
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
    required = ['Id1', 'Id2', 'total_ibd_prop', 'min_tmrca'] # 'gen_distance', 'true_ibd_prop']
    if not all(col in df.columns for col in required):
        raise ValueError(f"Ground truth missing required columns: {required}")

    # Modify the sample IDs to match the inferred metrics
    # Option 1: Modification des colonnes id1 et id2
    df['Id1'] = 'tsk_' + df['Id1'].astype(str)
    df['Id2'] = 'tsk_' + df['Id2'].astype(str)

    # Option 2 : La méthode .assign() (Style "Pipe")
    # df = df.assign(
    #    Id1 = lambda x: 'tsk_' + x['Id1'].astype(str),
    #    Id2 = lambda x: 'tsk_' + x['Id2'].astype(str)
    # )

    df = add_pair_key(df)

    # Create binary labels for each threshold
    for G in G_THRESHOLDS:
        df[f'related_G{G}'] = (df['min_tmrca'] <= G).astype(int)

    return df


# ============================================================================
# 2. LOAD INFERENCE RESULTS
# ============================================================================

def load_inference_results(scenario_dir, scenario_id):
    """
    Load IBD, IBS, and phylogenetic distance predictions

    Expected files:
    - ibd_results.csv: id1, id2, ibd_proportion (or ibd_score)
    - ibs_results.csv: id1, id2, ibs_similarity
    - phylo_results.csv: id1, id2, patristic_distance
    """
    results = {}
    scenario_dir = Path(scenario_dir) / scenario_id

    # --------------------------------------------------
    # 1. LOAD INFERRED IBD File
    # --------------------------------------------------
    ibd_file = f"{scenario_dir}/inferred_ibd_hmm.tsv"  # Adjust to your filename
    if os.path.exists(ibd_file):
        inferred_ibd = pd.read_csv(ibd_file, sep='\t')
        inferred_ibd = add_pair_key(inferred_ibd)
        results['IBD'] = inferred_ibd

    # --------------------------------------------------
    # 1. IBS File existence check
    # --------------------------------------------------
    ibs_file = f"{scenario_dir}/ibs_matrix.rds"
    if not os.path.exists(ibs_file):
        return None

    # --------------------------------------------------
    # 2. Safe readRDS equivalent
    # --------------------------------------------------
    try:
        result = pyreadr.read_r(ibs_file)
    except Exception as e:
        print(f"readRDS failed: {ibs_file}")
        return None

    if len(result.keys()) == 0:
        return None

    # Extract object (RDS contains single object)
    obj = list(result.values())[0]

    # --------------------------------------------------
    # 3. Convert to matrix
    # --------------------------------------------------
    if isinstance(obj, pd.DataFrame):
        mat = obj.copy()
    else:
        try:
            mat = pd.DataFrame(obj)
        except:
            return None

    # --------------------------------------------------
    # 4. Handle rownames
    # --------------------------------------------------
    ids = mat.index.astype(str)

    if ids.isnull().any() or len(ids) == 0:
        # Maybe first column contains IDs (like in R fallback)
        if isinstance(obj, pd.DataFrame):
            ids = obj.iloc[:, 0].astype(str).values
            mat = obj.iloc[:, 1:].copy()
            mat.index = ids
            mat.columns = ids
        else:
            return None

    ids = mat.index.astype(str).values
    mat_values = mat.values

    # --------------------------------------------------
    # 5. Convert upper triangle to long format
    # --------------------------------------------------
    n = len(ids)

    if mat_values.shape[0] != mat_values.shape[1]:
        print("Matrix is not square.")
        return None

    # Get upper triangle indices (excluding diagonal)
    i_idx, j_idx = np.triu_indices(n, k=1)

    long_df = pd.DataFrame({
        "Id1": ids[i_idx],
        "Id2": ids[j_idx],
        "IBS": mat_values[i_idx, j_idx].astype(float)
    })

    long_df = add_pair_key(long_df)

    results['IBS'] = long_df

    # --------------------------------------------------
    # 1. LOAD INFERRED PHYLOGENETIC DISTANCE File
    # --------------------------------------------------
    phylo_file = f"{scenario_dir}/phylo_results.csv"
    if os.path.exists(phylo_file):
        # Read phylogenetic distance file
        phylo = pd.read_csv(phylo_file)

        # Convert distance to similarity (higher = more related)
        # phylo['similarity'] = 1 / (1 + phylo['patristic_distance'])
        max_dist = phylo["patristic_distance"].max(skipna=True)
        phylo["similarity"] = 1 - (phylo["patristic_distance"] / max_dist)
        phylo = add_pair_key(phylo)
        results['Phylo'] = phylo
    else:
        raise ValueError(f"{phylo_file} does not exist ...")

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

# ============================================================================
# 4. EVALUATE SINGLE SCENARIO
# ============================================================================

def evaluate_scenario(scenario_dir, scenario_id, output):
    """
    Complete evaluation for one scenario
    """
    # Load data
    truth = load_ground_truth(scenario_dir, scenario_id)
    if truth is None:
        print(f"  ⚠️  No ground truth for {scenario_id}")
        return None

    inference = load_inference_results(scenario_dir, scenario_id)
    if not inference:
        print(f"  ⚠️  No inference results for {scenario_id}")
        return None

    # Extract scenario parameters from folder name
    # Pattern: MIG001_verylow_tight_moderate_mlow
    # Meaning: runid, rec_rate, bottleneck, est, migration
    print("\n Extract scenario parameters from {scenario_id}\n")
    params = parse_scenario_id(scenario_id)
    prefix = scenario_id.split('_')[0]

    # Extract sample size
    sampling_prop = infer_sampling_proportion(scenario_id)
    sample_size = infer_sample_size(scenario_id)

    results = []
    pr_curves = []  # NEW: Store PR curves

    # Evaluate each method × threshold combination
    print("\n Evaluate each method × threshold combination\n")
    for G in G_THRESHOLDS:

        for method_name, method_df in inference.items():
            # Merge with ground truth
            print(method_df.columns)
            print(f"\n Method Name: {method_name}")
            print("\n" + "="*50 + "\n")

            # Standardize the score column
            if method_name == "IBD":
                score_col = "hmm"
            elif method_name == "IBS":
                score_col = "IBS"
            elif method_name == "Phylo":
                score_col = "similarity"

            method_df = method_df.rename(columns={score_col: "score"})
            print(" Truth duplicate pair_key:", truth['pair_key'].duplicated().sum())
            print(" Method duplicate pair_key:", method_df['pair_key'].duplicated().sum())

            print("\n  Merge with ground truth ...\n")
            # Ensure uniqueness
            assert truth['pair_key'].is_unique
            assert method_df['pair_key'].is_unique

            merged = truth.merge(
                method_df[['score', 'pair_key']],  # Last column is the score 'Id1', 'Id2',
                on="pair_key",   #['Id1', 'Id2'],
                how='inner'
            )

            if merged.empty:
                print(f"{method_name}: no merged pairs")
                continue

            # Get predictions
            # Class diagnostic
            y_true = merged[f'related_G{G}'].values
            y_score = merged['score'].values

            n_pos = np.sum(y_true)
            n_neg = len(y_true) - n_pos
            n_total = len(y_true)
            prevalence = n_pos / n_total if n_total > 0 else np.nan

            if n_pos == 0 or n_neg == 0:
                print(f"{method_name}: single class at G={G}")
                results.append({
                    'scenario_id': scenario_id,
                    **params,
                    'replicate_id': extract_replicate_id(scenario_id),
                    'sampling_prop': sampling_prop,
                    'sample_size': sample_size,
                    'G_threshold': G,
                    'method': method_name,
                    'n_pairs': len(merged),
                    'n_positive': sum(y_true),
                    'prevalence': prevalence,
                    'auroc': np.nan, 'auc': np.nan,
                    'auprc': np.nan,
                    'TP': np.nan, 'FP': np.nan,
                    'TN': np.nan, 'FN': np.nan,
                    'sensitivity_at_90spec': np.nan,
                    'threshold_at_90spec': np.nan,
                    'identifiable': False,

                    # NEW: Metrics at optimal F1 threshold
                    'precision_opt': np.nan, 'recall_opt': np.nan,
                    'F1_score': np.nan, 'threshold_opt': np.nan,

                    # NEW: Metrics at 90% specificity threshold
                    'precision_90spec': np.nan,
                    'recall_90spec': np.nan,
                    'F1_90spec': np.nan,
                    'threshold_90spec': np.nan,

                    'note': 'single_class'
                })
                continue

            if len(np.unique(y_score)) < 2:
                # print(f"{method_name}: constant scores at G={G}")
                print(f"  ⚠️  Only one class present for {method_name} at G={G}. Skipping.")
                continue

            # Calculate metrics
            try:
                auroc = roc_auc_score(y_true, y_score)
                sens_at_90spec, threshold_90 = calculate_sensitivity_at_specificity(
                    y_true, y_score, target_spec=0.90)

                # ============================================
                #   Compute AUCPR (Precision-Recall AUC)
                # ============================================
                auprc = average_precision_score(y_true, y_score)
                precision, recall, _ = precision_recall_curve(y_true, y_score)
                auc_pr = auc(recall, precision)

                # ============================================================
                # NEW: Calculate Precision, Recall, F1 at optimal threshold
                # ============================================================
                # Option A: Use threshold that maximizes F1 score
                precision_curve, recall_curve, thresholds_curve = precision_recall_curve(y_true, y_score)


                # Identifiable?
                # identifiable = (auroc >= AUROC_THRESHOLD and sens_at_90spec >= SENSITIVITY_AT_90SPEC_THRESHOLD)
                identifiable = (
                    not np.isnan(auroc) and
                    auroc >= AUROC_THRESHOLD and
                    auprc >= prevalence * 2 and   # must beat random baseline
                    sens_at_90spec >= SENSITIVITY_AT_90SPEC_THRESHOLD)

                # ==============================================
                #  Compute Confusion Matrix Components
                # ==============================================
                # 1. Fixed threshold (e.g., 0.5)
                threshold = 0.5

                merged['pred'] = (merged['score'] >= threshold).astype(int)
                tn, fp, fn, tp = confusion_matrix(
                    y_true, merged['pred']
                ).ravel()

                # Calculate F1 for each threshold
                # (skip last point where precision/recall can be undefined)
                f1_scores = 2 * (precision_curve[:-1] * recall_curve[:-1]) / (precision_curve[:-1] + recall_curve[:-1] + 1e-10)

                # Find threshold that maximizes F1
                best_f1_idx = np.argmax(f1_scores)
                optimal_threshold = thresholds_curve[best_f1_idx]

                # Make binary predictions at optimal threshold
                y_pred = (y_score >= optimal_threshold).astype(int)

                # Calculate metrics at this threshold
                precision = precision_score(y_true, y_pred, zero_division=0)
                recall = recall_score(y_true, y_pred, zero_division=0)
                F1 = f1_score(y_true, y_pred, zero_division=0)

                # ============================================================
                # ALTERNATIVE: Calculate at threshold_90 (90% specificity)
                # ============================================================
                y_pred_90spec = (y_score >= threshold_90).astype(int)
                precision_90spec = precision_score(y_true, y_pred_90spec, zero_division=0)
                recall_90spec = recall_score(y_true, y_pred_90spec, zero_division=0)
                F1_90spec = f1_score(y_true, y_pred_90spec, zero_division=0)

                # # 2. F1-max threshold
                # best_f1 = 0
                # best_threshold = 0
                #
                # for t in np.linspace(0,1,100):
                #     pred = (merged['score'] >= t).astype(int)
                #     f1 = f1_score(merged['label'], pred)
                #     if f1 > best_f1:
                #         best_f1 = f1
                #         best_threshold = t

                results.append({
                     'scenario_id': scenario_id,
                     **params,
                     'replicate_id': extract_replicate_id(scenario_id),
                     'sampling_prop': sampling_prop,
                     'sample_size': sample_size,
                     'G_threshold': G,
                     'method': method_name,
                     'n_pairs': len(merged),
                     'n_positive': sum(y_true),
                     'prevalence': prevalence,
                     'auroc': auroc, 'auc': auc_pr,
                     'auprc': auprc,
                     'TP': tp, 'FP': fp,
                     'TN': tn, 'FN': fn,
                     'sensitivity_at_90spec': sens_at_90spec,
                     'threshold_at_90spec': threshold_90,

                     # NEW: Metrics at optimal F1 threshold
                     'precision_opt': precision, 'recall_opt': recall,
                     'F1_score': F1, 'threshold_opt': optimal_threshold,

                     # NEW: Metrics at 90% specificity threshold
                     'precision_90spec': precision_90spec,
                     'recall_90spec': recall_90spec,
                     'F1_90spec': F1_90spec,
                     'threshold_90spec': threshold_90,
                     'identifiable': identifiable,
                     'note': 'multiple_classes'
                })


                # ============================================================
                # NEW: Save PR curve data
                # ============================================================
                pr_curve_data = pd.DataFrame({
                    'scenario_id': scenario_id, **params,
                    'replicate_id': extract_replicate_id(scenario_id),
                    'sample_size': sample_size,
                    'method': method_name, 'G_threshold': G,
                    'precision': precision_curve[:-1],
                    'recall': recall_curve[:-1],
                    'threshold': thresholds_curve
                })
                pr_curves.append(pr_curve_data)

            except Exception as e:
                print(f"  ⚠️  Error evaluating {method_name} at G={G}: {e}")
                continue
    # Save PR curves to separate file
    if pr_curves:
        pr_curves_df = pd.concat(pr_curves, ignore_index=True)
        pr_output = f"{output}/{prefix}_prcurves.csv"
        pr_curves_df.to_csv(pr_output, index=False)

    return pd.DataFrame(results)

# ============================================================================
# 7. VISUALIZATION
# ============================================================================

# ==============================
#           ROC curve
# ==============================
def plot_roc_curves(scenario_dir, scenario_id, methods=None):
    """
    Plot ROC curves for all methods on same scenario
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # scenario_dir = Path(scenario_dir)
    # scenario_id = scenario_dir.name

    # Load data
    truth = load_ground_truth(scenario_dir, scenario_id)
    if truth is None:
        print(f"  ⚠️  No ground truth for {scenario_id}")
        return None

    inference = load_inference_results(scenario_dir, scenario_id)
    if not inference:
        print(f"  ⚠️  No inference results for {scenario_id}")
        return None

    G = 25
    for method_name, method_df in inference.items():
        # Merge with ground truth
        print(method_df.columns)
        print(f"\n Method Name: {method_name}")
        print("\n" + "="*50 + "\n")

        # Standardize the score column
        if method_name == "IBD":
            score_col = "hmm"
        elif method_name == "IBS":
            score_col = "IBS"
        elif method_name == "Phylo":
            score_col = "similarity"

        method_df = method_df.rename(columns={score_col: "score"})
        print(" Truth duplicate pair_key:", truth['pair_key'].duplicated().sum())
        print(" Method duplicate pair_key:", method_df['pair_key'].duplicated().sum())

        print("\n  Merge with ground truth ...\n")
        # Ensure uniqueness
        assert truth['pair_key'].is_unique
        assert method_df['pair_key'].is_unique

        merged = truth.merge(
            method_df[['score', 'pair_key']],
            on="pair_key", how='inner')

        if merged.empty:
            print(f"{method_name}: no merged pairs")
            continue

        # Get predictions
        # Class diagnostic
        y_true = merged[f'related_G{G}'].values
        y_score = merged['score'].values

        n_pos = np.sum(y_true)
        n_neg = len(y_true) - n_pos

        if n_pos == 0 or n_neg == 0:
            print(f"{method_name}: single class at G={G}")
            continue

        if len(np.unique(y_score)) < 2:
            print(f"{method_name}: constant scores at G={G}")
            continue

        # Calculate metrics
        try:
            # Calculate ROC
            fpr, tpr, _ = roc_curve(y_true, y_score)
            auc = roc_auc_score(y_true, y_score)

            # Plot
            ax.plot(fpr, tpr, label=f'{method_name.upper()} (AUC={auc:.2f})', linewidth=2)

        except Exception as e:
            print(f"  ⚠️  Error evaluating {method_name} at G={G}: {e}")
            continue

    # Format
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='Random')
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title(f'ROC Curves - {scenario_dir} at G={G}', fontsize=14)
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    return fig

# ==================================
#   AUCPR curve with bootstrap CIs
# ==================================
def pr_curve_with_ci(y_true, y_score, n_boot=1000, seed=42):
    """
    Compute mean PR curve with bootstrap confidence intervals
    """
    rng = np.random.RandomState(seed)

    precisions = []
    recalls = []

    # Reference recall grid
    recall_grid = np.linspace(0, 1, 200)

    for _ in range(n_boot):
        idx = resample(
            np.arange(len(y_true)),
            replace=True,
            random_state=rng)

        p, r, _ = precision_recall_curve(y_true[idx], y_score[idx])

        # Interpolate precision onto common recall grid
        p_interp = np.interp(recall_grid, r[::-1], p[::-1])
        precisions.append(p_interp)

    precisions = np.array(precisions)

    return {
        "recall": recall_grid,
        "precision_mean": precisions.mean(axis=0),
        "precision_low": np.percentile(precisions, 2.5, axis=0),
        "precision_high": np.percentile(precisions, 97.5, axis=0)
    }


# ==========================================================
# Scenario-level PR curves (single scenario, all methods)
# ==========================================================
def plot_pr_curves_single_scenario(scenario_dir, scenario_id):

    #  Define G threshold
    G = 25
    # Load data
    truth = load_ground_truth(scenario_dir, scenario_id)
    if truth is None:
        print(f"  ⚠️  No ground truth for {scenario_id}")
        return None

    inference = load_inference_results(scenario_dir, scenario_id)
    if not inference:
        print(f"  ⚠️  No inference results for {scenario_id}")
        return None

    fig, ax = plt.subplots(figsize=(7, 6))

    for method_name, method_df in inference.items():
        # Merge with ground truth
        print(method_df.columns)
        print(f"\n Method Name: {method_name}")
        print("\n" + "="*50 + "\n")

        # Standardize the score column
        if method_name == "IBD":
            score_col = "hmm"
        elif method_name == "IBS":
            score_col = "IBS"
        elif method_name == "Phylo":
            score_col = "similarity"

        method_df = method_df.rename(columns={score_col: "score"})
        print(" Truth duplicate pair_key:", truth['pair_key'].duplicated().sum())
        print(" Method duplicate pair_key:", method_df['pair_key'].duplicated().sum())

        print("\n  Merge with ground truth ...\n")
        # Ensure uniqueness
        assert truth['pair_key'].is_unique
        assert method_df['pair_key'].is_unique

        merged = truth.merge(
            method_df[['score', 'pair_key']],
            on="pair_key", how='inner')

        if merged.empty:
            print(f"{method_name}: no merged pairs")
            continue

        # Get predictions
        y_true = merged[f'related_G{G}'].values
        y_score = merged['score'].values

        n_pos = np.sum(y_true)
        n_neg = len(y_true) - n_pos

        if n_pos == 0 or n_neg == 0:
            print(f"{method_name}: single class at G={G}")
            continue

        if len(np.unique(y_score)) < 2:
            print(f"{method_name}: constant scores at G={G}")
            continue

        # Calculate metrics
        try:
            pr = pr_curve_with_ci(y_true, y_score)
            auprc = average_precision_score(y_true, y_score)

            ax.plot( pr["recall"], pr["precision_mean"],
                label=f'{method_name.upper()} (AUPRC={auprc:.2f})',
                linewidth=2)

            ax.fill_between(pr["recall"], pr["precision_low"], pr["precision_high"], alpha=0.2)

        except Exception as e:
            print(f"  ⚠️  Error evaluating {method_name} at G={G}: {e}")
            continue

    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    return fig

# Scenario-level aggregation (this is the key innovation)
# Now we answer: Across all epidemiological scenarios, which method is consistently better?

# ============================================================
# def aggregate_pr_across_scenarios(scenario_dirs, method):
# ============================================================
def aggregate_pr_across_scenarios(y_true, y_score):
    all_precisions = []

    recall_grid = np.linspace(0, 1, 200)

    for scen in scenario_dirs:
        truth = load_ground_truth(scen)
        inf = load_inference_results(scen, method)

        data = truth.merge(
            inf, on=['sample_i', 'sample_j'], how='inner'
        )

        y_true = data['is_direct_transmission'].values

        if method == 'ibd':
            y_score = data['ibd_relatedness'].values
        elif method == 'ibs':
            y_score = -data['ibs_distance'].values
        elif method == 'phylo':
            y_score = -data['phylo_distance'].values

        p, r, _ = precision_recall_curve(y_true, y_score)
        p_interp = np.interp(recall_grid, r[::-1], p[::-1])
        all_precisions.append(p_interp)

    all_precisions = np.array(all_precisions)

    return {
        "recall": recall_grid,
        "precision_mean": all_precisions.mean(axis=0),
        "precision_low": np.percentile(all_precisions, 2.5, axis=0),
        "precision_high": np.percentile(all_precisions, 97.5, axis=0)
    }

# ============================================================
# def aggregate_pr_across_scenarios(scenario_dirs, method):
# ============================================================
def aggregate_evaluation_results(df, metrics=("auprc",),
    group_cols=("scenario_id", "rec_rate_label", "bottleneck_label",
                "est_label", "migration_label", "method",
                "sampling_prop", "sample_size", "G_threshold"),
    ci=0.95, min_replicates=5
):
    """
    Aggregate evaluation results across replicates.

    Parameters
    ----------
    df : pd.DataFrame
        Output of evaluate_scenario() concatenated across scenarios
    metrics : tuple
        Metrics to aggregate (e.g. 'auprc', 'auroc')
    group_cols : tuple
        Columns defining aggregation groups
    ci : float
        Confidence interval level (default 95%)
    min_replicates : int
        Minimum replicates required to compute summary

    Returns
    -------
    pd.DataFrame
        Aggregated metrics with mean, std, CI bounds, and n_reps
    """

    z = {
        0.90: 1.645,
        0.95: 1.96,
        0.99: 2.576
    }[ci]

    agg_rows = []

    for group_keys, gdf in df.groupby(list(group_cols)):

        n_reps = gdf["replicate_id"].nunique()

        if n_reps < min_replicates:
            continue

        row = dict(zip(group_cols, group_keys))
        row["n_replicates"] = n_reps

        for m in metrics:
            vals = gdf[m].dropna()

            if len(vals) == 0:
                row[f"{m}_mean"] = np.nan
                row[f"{m}_median"] = np.nan
                row[f"{m}_ci_low"] = np.nan
                row[f"{m}_ci_high"] = np.nan
                continue

            mean = vals.mean()
            median = vals.median()
            sd = vals.std(ddof=1)
            se = sd / np.sqrt(len(vals))

            row[f"{m}_mean"] = mean
            row[f"{m}_median"] = median
            row[f"{m}_ci_low"] = mean - z * se
            row[f"{m}_ci_high"] = mean + z * se

        agg_rows.append(row)

    return pd.DataFrame(agg_rows)

# ===================================================
# Main manuscript panel (clean, minimal, powerful)
# ===================================================
def detect_ibd_collapse(agg_df, delta=0.05, random_baseline=0.02):
    """
    Detect scenarios where IBD collapses relative to IBS.

    Returns a DataFrame with collapse flags.
    """

    ibd = agg_df[agg_df["method"] == "IBD"]
    ibs = agg_df[agg_df["method"] == "IBS"]

    merged = ibd.merge(
        ibs,
        on=["sample_size"],
        suffixes=("_ibd", "_ibs")
    )

    merged["ibd_collapse"] = (
        (merged["auprc_mean_ibd"] + delta < merged["auprc_mean_ibs"]) &
        (merged["auprc_ci_high_ibd"] < random_baseline)
    )

    return merged


def plot_pr_main_panel(pr_main, collapse_df):
    """
    Main manuscript PR panel with confidence bands and IBD failure annotations.
    """

    fig, ax = plt.subplots(figsize=(7, 5))

    methods = ["IBS", "IBD", "Phylo"]
    colors = {
        "IBS": "#2c7fb8",
        "IBD": "#d7191c",
        "Phylo": "#31a354"
    }

    for method in methods:
        df = pr_main[pr_main["method"] == method].sort_values("sample_size")

        ax.plot(
            df["sample_size"],
            df["auprc_mean"],
            label=method,
            linewidth=2,
            color=colors[method]
        )

        ax.fill_between(
            df["sample_size"],
            df["auprc_ci_low"],
            df["auprc_ci_high"],
            alpha=0.25,
            color=colors[method]
        )

    # --- Failure annotations (IBD collapse) ---
    collapse_points = collapse_df[collapse_df["ibd_collapse"]]

    for _, row in collapse_points.iterrows():
        ax.annotate(
            "IBD collapse",
            xy=(row["sample_size"], row["auprc_mean_ibd"]),
            xytext=(row["sample_size"], row["auprc_mean_ibd"] + 0.1),
            arrowprops=dict(arrowstyle="->", lw=1),
            fontsize=9,
            color="black",
            ha="center"
        )

    ax.set_xlabel("Sample size")
    ax.set_ylabel("AUPRC")
    ax.set_title("Precision–Recall Performance vs Sampling Density")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)

    return fig

# ============================================================================================
# Supplementary multi-G panel: This is not for interpretation — it’s for robustness evidence.
# ============================================================================================
def plot_pr_by_G(pr_by_G):
    """
    Supplementary figure: PR across G thresholds.
    """

    fig, axes = plt.subplots(1, len(pr_by_G["G_threshold"].unique()),
                             figsize=(16, 4), sharey=True)

    for ax, G in zip(axes, sorted(pr_by_G["G_threshold"].unique())):
        dfG = pr_by_G[pr_by_G["G_threshold"] == G]

        for method in ["IBS", "IBD", "Phylo"]:
            dfm = dfG[dfG["method"] == method].sort_values("sample_size")

            ax.plot(
                dfm["sample_size"],
                dfm["auprc_mean"],
                label=method,
                linewidth=2
            )

            ax.fill_between(
                dfm["sample_size"],
                dfm["auprc_ci_low"],
                dfm["auprc_ci_high"],
                alpha=0.2
            )

        ax.set_title(f"G = {G}")
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("AUPRC")
    axes[0].legend(frameon=False)

    return fig

# =========================================================
# =========================================================
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

def analyze_parameter_effects(results_df, outcome_col="auprc", min_samples=10):
    """
    Fit regression model to quantify parameter effects on AUC

    Quantify standardized parameter effects on performance (AUPRC by default).

    This function:
    - Fits separate regression models per method
    - Drops failed runs (NaN outcomes)
    - Standardizes predictors for effect-size comparison
    - Returns tidy coefficient table for plotting

    Parameters
    ----------
    results_df : pd.DataFrame
        Output of evaluation pipeline (long format)
    outcome_col : str
        Performance metric to model (default: 'auprc')
    methods : tuple
        Methods to analyze
    min_samples : int
        Minimum rows required to fit a model

    Returns
    -------
    pd.DataFrame
        Columns: method, parameter, coefficient, abs_effect
    """
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler

    # Prepare predictors (log-transform rates)
    # X = results_df[['rec_rate', 'bottleneck', 'est', 'migration', 'sample_size']].copy() # , 'sampling'
    # X['rec_rate'] = np.log10(X['rec_rate'])
    # X['migration'] = np.log10(X['migration'])

    # Standardize for effect size comparison
    # scaler = StandardScaler()
    # X_scaled = scaler.fit_transform(X)

    effects = []
    for method in ['IBD', 'IBS', 'Phylo']:
        # Subset to method
        # method_data = results_df[results_df['method'] == method]
        method_data = results_df[results_df['method'] == method].copy()

        # Drop failed evaluations
        method_data = method_data.dropna(subset=[outcome_col])

        if method_data.shape[0] < min_samples:
            print(f"⚠️  Skipping {method}: insufficient valid samples ({method_data.shape[0]})")
            continue

        # Prepare predictors
        X = method_data[['rec_rate', 'bottleneck', 'est', 'migration', 'sample_size']].copy()

        # Guard against log(0)
        X['rec_rate'] = np.log10(X['rec_rate'].clip(lower=1e-12))
        X['migration'] = np.log10(X['migration'].clip(lower=1e-12))

        # Outcome
        if 'auprc' in method_data:
            y = method_data['auprc'].values
        elif 'auc' in method_data:
            y = method_data['auc'].values
        else:
            raise ValueError("Expected column 'auprc' not found")

        # Standardize predictors
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Fit regression
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


# =================================
#       PERFORMANCE HEATMAP
# =================================
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

# =======================================
#       DIAGRAM PHASE IDENTIFIABILITY
# =======================================
def plot_identifiability_phase_diagram(results_df):
    """
    Phase diagram showing identifiable parameter space
    """
    # Get best AUC across methods for each parameter combination
    best_auc = results_df.groupby(['rec_rate', 'sample_size'])['auc'].max().reset_index() # 'G_threshold',

    # Pivot for contour plot
    pivot = best_auc.pivot(index='sample_size', columns='rec_rate', values='auc') # index='sampling'

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
        'total_scenarios': len(results_df['scenario_id'].unique()),
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
    for samp in results_df['sample_size'].unique():
        samp_data = results_df[results_df['sample_size'] == samp]
        summary['by_sampling'][samp] = {
            'pct_identifiable': (samp_data['auc'] >= 0.80).mean() * 100,
            'mean_auc': samp_data['auc'].mean()
        }

    return summary

def summarize_identifiability_v1(results_df):
    """
    Generate summary statistics for identifiability
    JSON-safe output
    """

    summary = {
        "overall": {},
        "by_method": {},
        "by_recombination": {},
        "by_sampling": {}
    }

    # ---- Overall ----
    auc_vals = results_df["auc"].dropna()

    summary["overall"] = {
        "total_scenarios": int(results_df["scenario_id"].nunique()),
        "pct_identifiable": float((auc_vals >= 0.80).mean() * 100) if len(auc_vals) else None,
        "mean_auc": float(auc_vals.mean()) if len(auc_vals) else None,
        "median_auc": float(auc_vals.median()) if len(auc_vals) else None
    }

    # ---- By method ----
    for method in results_df["method"].unique():
        method_data = results_df[results_df["method"] == method]
        auc_m = method_data["auc"].dropna()

        summary["by_method"][str(method)] = {
            "pct_identifiable": float((auc_m >= 0.80).mean() * 100) if len(auc_m) else None,
            "mean_auc": float(auc_m.mean()) if len(auc_m) else None,
            "n_runs": int(len(auc_m))
        }

    # ---- By recombination rate ----
    for rec_rate in results_df["rec_rate"].unique():
        rec_data = results_df[results_df["rec_rate"] == rec_rate]
        auc_r = rec_data["auc"].dropna()

        if len(auc_r) > 0:
            best_row = rec_data.loc[auc_r.idxmax()]
            best_method = best_row["method"]
        else:
            best_method = None

        summary["by_recombination"][str(rec_rate)] = {
            "pct_identifiable": float((auc_r >= 0.80).mean() * 100) if len(auc_r) else None,
            "mean_auc": float(auc_r.mean()) if len(auc_r) else None,
            "best_method": best_method
        }

    # ---- By sampling size ----
    for samp in results_df["sample_size"].unique():
        samp_data = results_df[results_df["sample_size"] == samp]
        auc_s = samp_data["auc"].dropna()

        summary["by_sampling"][str(samp)] = {
            "pct_identifiable": float((auc_s >= 0.80).mean() * 100) if len(auc_s) else None,
            "mean_auc": float(auc_s.mean()) if len(auc_s) else None,
            "n_runs": int(len(auc_s))
        }

    return summary

def summary_to_dataframes(summary):
    """
    Convert summarize_identifiability() output into tidy DataFrames.
    Returns a dict of DataFrames.
    """

    dfs = {}

    # --- Overall ---
    dfs["overall"] = (
        pd.DataFrame([summary["overall"]])
        .assign(scope="overall")
    )

    # --- By method ---
    dfs["by_method"] = (
        pd.DataFrame.from_dict(summary["by_method"], orient="index")
        .reset_index()
        .rename(columns={"index": "method"})
    )

    # --- By recombination ---
    dfs["by_recombination"] = (
        pd.DataFrame.from_dict(summary["by_recombination"], orient="index")
        .reset_index()
        .rename(columns={"index": "rec_rate"})
    )

    # --- By sampling ---
    dfs["by_sampling"] = (
        pd.DataFrame.from_dict(summary["by_sampling"], orient="index")
        .reset_index()
        .rename(columns={"index": "sample_size"})
    )

    return dfs

# ============================================================================
# 6. SANITY CHECK ON PILOT SCENARIOS
# ============================================================================

def run_sanity_check(base_dir):
    """
    Test on 3 representative scenarios before full run
    """
    print("\n" + "="*80)
    print("    SANITY CHECK: Testing 3 pilot scenarios")
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
        else:
            print(f"\n {scenario_dir} found ..........")

        print(f"\n{'='*80}")
        print(f"   Evaluating: {scenario_dir} ID: {scenario_id}")
        print(f"{'='*80}")

#        results_df = evaluate_scenario(scenario_dir, scenario_id)
        results_df = evaluate_scenario(Path(base_dir), scenario_id)

        if results_df is None:
            print("⚠️ evaluate_scenario returned None")
            continue

        if results_df.empty:
            print("⚠️ Results dataframe is empty")
            continue

        if results_df is not None:

            print("\nRaw results preview:")
            print(results_df.head())
            print("\nAvailable methods:", results_df['method'].unique())
            print("\nResults summary:")
            print("\n" + "="* 40)
            print(results_df.groupby(['method', 'G_threshold'])['auroc'].agg(['mean', 'count']))
            print("\nResults summary AUCPR:")
            print("\n" + "="* 40)
            print(results_df.groupby(['method', 'G_threshold'])['auc'].agg(['mean', 'count']))

            # Inspect class balance if stored
            if 'prevalence' in results_df.columns:
                print("Mean prevalence:", results_df['prevalence'].mean())

            # Check hypothesis for low recombination scenario
            if 'verylow' in scenario_id:
                print("\n✓ Expected: IBS > IBD for low recombination")
                ibs_auc = results_df[results_df['method'] == 'IBS']['auroc'].mean()
                ibd_auc = results_df[results_df['method'] == 'IBD']['auroc'].mean()
                dist_auc = results_df[results_df['method'] == 'Phylo']['auroc'].mean()
                print(f"  IBS AUC: {ibs_auc:.3f}")
                print(f"  IBD AUC: {ibd_auc:.3f}")
                print(f"  Phylo AUC: {dist_auc:.3f}")

                ibs_aucpr = results_df[results_df['method'] == 'IBS']['auc'].mean()
                ibd_aucpr = results_df[results_df['method'] == 'IBD']['auc'].mean()
                dist_aucpr = results_df[results_df['method'] == 'Phylo']['auc'].mean()
                print(f"  IBS AUPRC: {ibs_aucpr:.3f}")
                print(f"  IBD AUPRC: {ibd_aucpr:.3f}")
                print(f"  Phylo AUPRC: {dist_aucpr:.3f}")

                if ibs_auc > ibd_auc:
                    print("  ✓ PASS: IBS outperforms IBD")
                else:
                    print("  ✗ FAIL: IBD outperforms IBS (unexpected!)")

            # Check for high recombination
            if 'high' in scenario_id:
                print("\n✓ Expected: IBD > IBS for high recombination")
                ibs_auc = results_df[results_df['method'] == 'IBS']['auroc'].mean()
                ibd_auc = results_df[results_df['method'] == 'IBD']['auroc'].mean()
                dist_auc = results_df[results_df['method'] == 'Phylo']['auroc'].mean()
                print(f"  IBS AUC: {ibs_auc:.3f}")
                print(f"  IBD AUC: {ibd_auc:.3f}")
                print(f"  Phylo AUC: {dist_auc:.3f}")

                ibs_aucpr = results_df[results_df['method'] == 'IBS']['auc'].mean()
                ibd_aucpr = results_df[results_df['method'] == 'IBD']['auc'].mean()
                dist_aucpr = results_df[results_df['method'] == 'Phylo']['auc'].mean()
                print(f"  IBS AUPRC: {ibs_aucpr:.3f}")
                print(f"  IBD AUPRC: {ibd_aucpr:.3f}")
                print(f"  Phylo AUPRC: {dist_aucpr:.3f}")

                if ibd_auc > ibs_auc:
                    print("  ✓ PASS: IBD outperforms IBS")
                else:
                    print("  ✗ FAIL: IBS outperforms IBD (unexpected!)")

            # Check for medium recombination
            if 'medium' in scenario_id: #  and 'rec' in scenario_id
                print("\n✓ Expected: IBD > IBS for high recombination")
                ibs_auc = results_df[results_df['method'] == 'IBS']['auroc'].mean()
                ibd_auc = results_df[results_df['method'] == 'IBD']['auroc'].mean()
                dist_auc = results_df[results_df['method'] == 'Phylo']['auroc'].mean()
                print(f"  IBS AUC: {ibs_auc:.3f}")
                print(f"  IBD AUC: {ibd_auc:.3f}")
                print(f"  Phylo AUC: {dist_auc:.3f}")

                ibs_aucpr = results_df[results_df['method'] == 'IBS']['auc'].mean()
                ibd_aucpr = results_df[results_df['method'] == 'IBD']['auc'].mean()
                dist_aucpr = results_df[results_df['method'] == 'Phylo']['auc'].mean()
                print(f"  IBS AUPRC: {ibs_aucpr:.3f}")
                print(f"  IBD AUPRC: {ibd_aucpr:.3f}")
                print(f"  Phylo AUPRC: {dist_aucpr:.3f}")

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

        prefix = scenario_id.split('_')[0]
        output = Path(output_file).parent
        outdir = output / "plots"
        outdir.mkdir(parents=True, exist_ok=True)

        if i % 50 == 0:
            print(f"  Progress: {i}/{len(scenario_dirs)} scenarios processed")

        results_df = evaluate_scenario(Path(base_dir), scenario_id, output)

        if results_df is not None:
            all_results.append(results_df)

        # 5. Generate plots
        # scenario = Path(base_dir) / scenario_id
        fig_roc = plot_roc_curves(Path(base_dir), scenario_id)
        fig_pr = plot_pr_curves_single_scenario(Path(base_dir), scenario_id)

        fig_roc.savefig(f"{outdir}/{prefix}_auroc.png", dpi=300)
        fig_pr.savefig(f"{outdir}/{prefix}_aupr.png", dpi=300)

        if results_df is not None:

            print("\nRaw results preview:")
            print(results_df.head())
            print("\nAvailable methods:", results_df['method'].unique())
            print("\nResults summary:")
            print("\n" + "="* 40)
            print(results_df.groupby(['method', 'G_threshold'])['auroc'].agg(['mean', 'count']))
            print("\nResults summary AUCPR:")
            print("\n" + "="* 40)
            print(results_df.groupby(['method', 'G_threshold'])['auc'].agg(['mean', 'count']))

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
        'identifiable': 'mean', 'auroc': 'mean',
        'auc': 'mean', 'auprc': 'mean'
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
        'auroc': 'mean',
        'auc': 'mean'
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
        'auroc': lambda x: x.max() - x.min(),
        'auc': lambda x: x.max() - x.min()
    }).reset_index()

    auc_range.columns = ['scenario_id', 'G_threshold', 'auroc_range', 'auc_range']

    # Classify failures
    auc_range['failure_type'] = 'SUCCESS'
    auc_range.loc[auc_range['auc_range'] < 0.10, 'failure_type'] = 'FUNDAMENTAL'
    auc_range.loc[auc_range['auc_range'] >= 0.15, 'failure_type'] = 'METHODOLOGICAL'


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

    # 1. Sanity Check
    if args.sanity_check:
        run_sanity_check(args.base_dir)

    # 2. Run full evaluation across folders
    if args.full:
#        results = run_full_analysis(args.base_dir, args.output)

        output = Path(args.output)
        identifiability = output.parent / "identifiability_results.csv"
#        results.to_csv(identifiability, index=False)

        results = pd.read_csv(identifiability)

        # Aggregate for PR main panel
#        pr_main = aggregate_evaluation_results(results, metrics=("auprc",),
#            group_cols=("rec_rate_label", "method", "sample_size"))

        # Aggregate for supplement (multi-G comparison)
#        pr_by_G = aggregate_evaluation_results(
#            results,
#            metrics=("auprc",),
#            group_cols=("rec_rate_label", "method", "sample_size", "G_threshold")
#        )

#        pr_main_df = output.parent / "summary_results.csv"
#        pr_main.to_csv(pr_main_df, index=False)
#        pr_by_G_df = output.parent / "summary_results_byG.csv"
#        pr_by_G.to_csv(pr_by_G_df, index=False)

#        collapse_df = detect_ibd_collapse(pr_main)
#        main_plot = plot_pr_main_panel(pr_main, collapse_df)
#        supp = plot_pr_by_G(pr_by_G)

        # 3. Answer Research Questions
#        answer_research_questions(results)

#        heatmap = plot_performance_heatmap(results)
#        diagram = plot_identifiability_phase_diagram(results)

        out_file = output.parent
#        main_plot.savefig(f"{out_file}/main_figure.png", dpi=300)
#        supp.savefig(f"{out_file}/supp_figure.png", dpi=300)
#        heatmap.savefig(f"{out_file}/performance_heatmap.png", dpi=300)
#        diagram.savefig(f"{out_file}/identifiability.png", dpi=300)

        # 4. Parameter effects
        #effects = analyze_parameter_effects(results)
        #print(effects.sort_values('effect_size', ascending=False))
        #effects.to_csv(f"{out_file}/effect_size.csv", index=False)

        # 6. Summary
        summary = summarize_identifiability_v1(results)
        print(json.dumps(summary, indent=2))
        dfs = summary_to_dataframes(summary)

        for name, df in dfs.items():
            df.to_csv(out_file / f"{name}.csv", index=False)

        overall_latex = dfs["overall"].rename(columns={"total_scenarios": "Total scenarios", "pct_identifiable": "\\% identifiable", "mean_auc": "Mean AUROC","median_auc": "Median AUROC"})

        print(overall_latex.to_latex(index=False, float_format="%.2f",
              caption="Overall identifiability performance across all scenarios.",
              label="tab:overall_identifiability",
              escape=False
             )
        )

    if not args.sanity_check and not args.full:
        print("\n" + "="*80)
        print("  Please specify --sanity-check or --full \n")
        parser.print_help()
        print("\n" + "="*80)


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
