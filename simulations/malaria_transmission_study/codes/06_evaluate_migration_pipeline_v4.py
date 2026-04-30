#!/usr/bin/env python3

"""
Evaluation script for migration scenarios with tree-sequence ground truth
Adapted for: true_ibd from tree sequences + generation distance thresholds
"""

import numpy as np
import pandas as pd
import os
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
# ============================================================================
# 4. HELPER: STANDARDISE SCORE COLUMN
# ============================================================================

def _standardise_score_col(method_df, method_name):
    """
    Rename the method-specific score column to a canonical 'score' column.

    IBD score priority : max_seg_bp > total_ibd_bp > hmm > n_segments
    ─────────────────────────────────────────────────────────────────────
    max_seg_bp is the biologically correct primary score for IBD: only
    recent pairs (few generations) share a very long unbroken haplotype.
    Recombination erodes long segments with time, so max_seg_bp decays
    rapidly with generation distance and gives the best discrimination
    between recent transmission and old ancestry.

    total_ibd_bp collapses at high recombination because background short
    segments accumulate and the distribution converges for all pairs.

    ibd_prop / true_ibd_prop are intentionally excluded — they are the
    tree-sequence TRUTH, not an inferred score.
    """
    df = method_df.copy()
    if method_name == "IBD":
        ibd_priority = ["max_seg_bp", "max_segment_bp",
                        "total_ibd_bp", "total_ibd", "hmm", "n_segments"]
        chosen = next((c for c in ibd_priority if c in df.columns), None)
        if chosen is None:
            raise ValueError(
                f"IBD DataFrame has no recognised score column. "
                f"Columns: {list(df.columns)}"
            )
        if chosen not in ("max_seg_bp", "max_segment_bp"):
            print(f"  ⚠️  IBD: using '{chosen}' as score — max_seg_bp preferred but not found")
        else:
            print(f"  ✓  IBD: using '{chosen}' (max segment length)")
        df = df.rename(columns={chosen: "score"})
    elif method_name == "IBS":
        df = df.rename(columns={"IBS": "score"})
    elif method_name == "Phylo":
        df = df.rename(columns={"similarity": "score"})
    else:
        raise ValueError(f"Unknown method: {method_name}")
    return df


# ============================================================================
# 4. HELPER: BUILD SHARED WIDE EVALUATION TABLE
# ============================================================================

def _build_merged_table(truth, inference):
    """
    Build one wide table from ground truth + all methods. Left-join so the
    truth universe is preserved for every method.

    IBD zero-fill
    ─────────────
    The HMM only writes pairs where it detected ≥1 segment.  Pairs absent
    from the HMM output received zero segments — that IS a score of 0, not
    missing data.  Replacing NaN with 0 after the left join restores all
    true negatives and makes AUROC/AUPRC meaningful for IBD.

    Without this fix, an inner join on IBD pairs keeps only detected pairs.
    That biased sample contains almost no true negatives, so AUROC = 0.5
    exactly, and IBD cannot be compared fairly with IBS or Phylo.

    IBS / Phylo NaNs are genuine missingness and left as NaN — they are
    excluded from that method's evaluation only.
    """
    merged = truth.copy()
    col_map = {"IBD": "score_IBD", "IBS": "score_IBS", "Phylo": "score_Phylo"}

    for method_name, method_df in inference.items():
        df     = _standardise_score_col(method_df, method_name)
        out    = col_map[method_name]
        assert df['pair_key'].is_unique, f"{method_name}: duplicate pair_keys"

        merged = merged.merge(
            df[['pair_key', 'score']].rename(columns={'score': out}),
            on='pair_key',
            how='left'          # ← preserves truth universe for all methods
        )

        # ── IBD zero-fill ────────────────────────────────────────────────
        if method_name == "IBD":
            n_na = int(merged[out].isna().sum())
            merged[out] = merged[out].fillna(0.0)
            if n_na > 0:
                print(f"    ✓ IBD: {n_na:,} undetected pairs assigned score = 0")

    score_cols = [c for c in merged.columns if c.startswith("score_")]
    print(f"  ✓ Merged: {len(merged):,} rows × {merged.shape[1]} cols "
          f"| scores: {score_cols}")
    return merged


# ============================================================================
# 4. HELPER: COMPUTE METRICS FOR ONE METHOD × ONE G
# ============================================================================

def _compute_metrics(y_true, y_score, calib_threshold):
    """
    Full metric suite for one (y_true, y_score) pair.
    calib_threshold: Youden threshold from baseline G; used for confusion matrix.
    If None, falls back to the F1-maximising threshold.
    """
    nan_row = dict(
        auroc=np.nan, auprc=np.nan, auc_integral=np.nan,
        sensitivity_at_90spec=np.nan, threshold_at_90spec=np.nan,
        precision_opt=np.nan, recall_opt=np.nan,
        F1_score=np.nan, threshold_opt=np.nan,
        precision_90spec=np.nan, recall_90spec=np.nan,
        F1_90spec=np.nan, threshold_90spec=np.nan,
        TP=np.nan, FP=np.nan, TN=np.nan, FN=np.nan,
        identifiable=False, note='error'
    )

    n_pos      = int(np.sum(y_true))
    n_neg      = int(len(y_true) - n_pos)
    prevalence = n_pos / len(y_true)
    nan_row.update(n_pairs=len(y_true), n_positive=n_pos, prevalence=prevalence)

    if n_pos == 0 or n_neg == 0:
        nan_row['note'] = 'single_class'
        return nan_row
    if len(np.unique(y_score)) < 2:
        nan_row['note'] = 'constant_scores'
        return nan_row

    try:
        # Ranking metrics (threshold-free)
        auroc      = roc_auc_score(y_true, y_score)
        auprc      = average_precision_score(y_true, y_score)   # sklearn AP
        prec_c, rec_c, thr_c = precision_recall_curve(y_true, y_score)
        auc_integ  = auc(rec_c, prec_c)
        sens90, thr90 = calculate_sensitivity_at_specificity(
            y_true, y_score, target_spec=0.90)

        identifiable = (
            not np.isnan(auroc)
            and auroc >= AUROC_THRESHOLD
            and auprc >= prevalence * 2
            and sens90 >= SENSITIVITY_AT_90SPEC_THRESHOLD
        )

        # Threshold-based metrics
        f1_scores   = (2 * prec_c[:-1] * rec_c[:-1]
                       / (prec_c[:-1] + rec_c[:-1] + 1e-10))
        opt_idx     = int(np.argmax(f1_scores))
        thr_opt     = float(thr_c[opt_idx])
        thr_cm      = calib_threshold if calib_threshold is not None else thr_opt

        y_opt       = (y_score >= thr_opt).astype(int)
        y_90        = (y_score >= thr90).astype(int)
        y_cm        = (y_score >= thr_cm).astype(int)

        tn, fp, fn, tp = confusion_matrix(y_true, y_cm).ravel()

        return dict(
            n_pairs=len(y_true), n_positive=n_pos, prevalence=prevalence,
            auroc=auroc, auprc=auprc, auc_integral=auc_integ,
            sensitivity_at_90spec=sens90, threshold_at_90spec=thr90,
            precision_opt=precision_score(y_true, y_opt, zero_division=0),
            recall_opt=recall_score(y_true, y_opt, zero_division=0),
            F1_score=f1_score(y_true, y_opt, zero_division=0),
            threshold_opt=thr_opt,
            precision_90spec=precision_score(y_true, y_90, zero_division=0),
            recall_90spec=recall_score(y_true, y_90, zero_division=0),
            F1_90spec=f1_score(y_true, y_90, zero_division=0),
            threshold_90spec=thr90,
            TP=int(tp), FP=int(fp), TN=int(tn), FN=int(fn),
            identifiable=bool(identifiable), note='ok'
        )

    except Exception as exc:
        nan_row['note'] = f'error: {exc}'
        return nan_row


# ============================================================================
# 4. EVALUATE SINGLE SCENARIO
# ============================================================================

def evaluate_scenario(scenario_dir, scenario_id, output=None):
    """
    Complete evaluation for one scenario.

    Architecture change vs original
    ────────────────────────────────
    Original: per-method inner join → each method evaluated on a different
    pair universe → IBD inner join excludes all undetected pairs (= all true
    negatives) → AUROC = 0.5, methods not comparable across the same pairs.

    Corrected: single wide left-join on truth universe → IBD NAs zero-filled
    → all three methods evaluated on IDENTICAL pairs → AUROC/AUPRC are fair.
    """
    truth = load_ground_truth(scenario_dir, scenario_id)
    if truth is None:
        print(f"  ⚠️  No ground truth for {scenario_id}")
        return None

    inference = load_inference_results(scenario_dir, scenario_id)
    if not inference:
        print(f"  ⚠️  No inference results for {scenario_id}")
        return None

    assert truth['pair_key'].is_unique, "Ground truth has duplicate pair_keys"

    print(f"\n  Extracting parameters from {scenario_id}")
    params        = parse_scenario_id(scenario_id)
    prefix        = scenario_id.split('_')[0]
    sampling_prop = infer_sampling_proportion(scenario_id)
    sample_size   = infer_sample_size(scenario_id)

    # ── Single wide-table merge (the core fix) ─────────────────────────
    print("  Building shared evaluation table...")
    merged = _build_merged_table(truth, inference)

    # ── Calibrate thresholds once at G=5 (Youden's J) ──────────────────
    # Used only for confusion-matrix metrics; AUPRC/AUROC are threshold-free.
    calib_thresholds = {}
    ref_label = f'related_G{min(G_THRESHOLDS)}'     # strictest G available
    col_map   = {"IBD": "score_IBD", "IBS": "score_IBS", "Phylo": "score_Phylo"}

    for method_name in inference:
        sc = col_map[method_name]
        if sc not in merged.columns or ref_label not in merged.columns:
            calib_thresholds[method_name] = None
            continue
        valid = merged[[sc, ref_label]].dropna()
        if valid[ref_label].nunique() < 2 or valid[sc].nunique() < 2:
            calib_thresholds[method_name] = None
            continue
        try:
            fpr, tpr, thr = roc_curve(valid[ref_label], valid[sc])
            best           = int(np.argmax(tpr - fpr))
            calib_thresholds[method_name] = float(thr[best])
        except Exception:
            calib_thresholds[method_name] = None

    # ── Evaluate each method × G threshold ─────────────────────────────
    print(f"\n  Evaluating {len(inference)} methods × {len(G_THRESHOLDS)} G thresholds...")
    results   = []
    pr_curves = []

    for G in G_THRESHOLDS:
        label_col = f'related_G{G}'
        if label_col not in merged.columns:
            continue

        for method_name in inference:
            sc = col_map[method_name]
            if sc not in merged.columns:
                continue

            sub = merged[[sc, label_col]].dropna()
            if sub.empty:
                continue

            y_true  = sub[label_col].values.astype(int)
            y_score = sub[sc].values.astype(float)
            m       = _compute_metrics(y_true, y_score, calib_thresholds.get(method_name))

            results.append({
                'scenario_id':  scenario_id,
                **params,
                'replicate_id': extract_replicate_id(scenario_id),
                'sampling_prop': sampling_prop,
                'sample_size':  sample_size,
                'G_threshold':  G,
                'method':       method_name,
                **m
            })

            # Store PR curve data
            if m['note'] == 'ok':
                prec_c, rec_c, thr_c = precision_recall_curve(y_true, y_score)
                pr_curves.append(pd.DataFrame({
                    'scenario_id':  scenario_id,
                    **params,
                    'replicate_id': extract_replicate_id(scenario_id),
                    'sampling_prop': sampling_prop,
                    'sample_size':  sample_size,
                    'method':       method_name,
                    'G_threshold':  G,
                    'precision':    prec_c[:-1],
                    'recall':       rec_c[:-1],
                    'threshold':    thr_c,
                }))

    # Save PR curves
    if pr_curves and output is not None:
        pr_df    = pd.concat(pr_curves, ignore_index=True)
        out_path = Path(output) / f"{prefix}_prcurves.csv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        pr_df.to_csv(out_path, index=False)
        print(f"  ✓ PR curves → {out_path}")

    return pd.DataFrame(results)


# ============================================================================
# 7. VISUALIZATION
# ============================================================================

# ==============================
#           ROC curve
# ==============================
def plot_roc_curves(scenario_dir, scenario_id, G=25, methods=None):
    """
    Plot ROC curves for all methods on the same shared pair universe.

    Uses _build_merged_table so that IBD is evaluated on the full truth
    universe (undetected pairs zero-filled) — the same set as IBS and Phylo.
    """
    truth = load_ground_truth(scenario_dir, scenario_id)
    if truth is None:
        print(f"  ⚠️  No ground truth for {scenario_id}")
        return None

    inference = load_inference_results(scenario_dir, scenario_id)
    if not inference:
        print(f"  ⚠️  No inference results for {scenario_id}")
        return None

    merged    = _build_merged_table(truth, inference)
    label_col = f'related_G{G}'
    if label_col not in merged.columns:
        print(f"  ⚠️  Label column {label_col} not found")
        return None

    col_map = {"IBD": "score_IBD", "IBS": "score_IBS", "Phylo": "score_Phylo"}
    fig, ax = plt.subplots(figsize=(8, 6))

    for method_name in inference:
        sc = col_map[method_name]
        if sc not in merged.columns:
            continue
        if methods and method_name not in methods:
            continue

        sub = merged[[sc, label_col]].dropna()
        y_true  = sub[label_col].values.astype(int)
        y_score = sub[sc].values.astype(float)

        if len(np.unique(y_true)) < 2 or len(np.unique(y_score)) < 2:
            print(f"  ⚠️  {method_name}: degenerate at G={G}, skipping")
            continue

        try:
            fpr, tpr, _ = roc_curve(y_true, y_score)
            auroc = roc_auc_score(y_true, y_score)
            ax.plot(fpr, tpr, label=f'{method_name} (AUROC={auroc:.2f})', linewidth=2)
        except Exception as e:
            print(f"  ⚠️  {method_name}: {e}")

    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='Random')
    ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
    ax.set_title(f'ROC Curves — {scenario_id}  (G={G})', fontsize=14, fontweight='bold')
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
def plot_pr_curves_single_scenario(scenario_dir, scenario_id, G=25):
    """
    Plot PR curves with bootstrap CIs for all methods on shared pair universe.

    Uses _build_merged_table so IBD is evaluated on the full truth universe
    (undetected pairs zero-filled), the same set as IBS and Phylo.
    """
    truth = load_ground_truth(scenario_dir, scenario_id)
    if truth is None:
        print(f"  ⚠️  No ground truth for {scenario_id}")
        return None

    inference = load_inference_results(scenario_dir, scenario_id)
    if not inference:
        print(f"  ⚠️  No inference results for {scenario_id}")
        return None

    merged    = _build_merged_table(truth, inference)
    label_col = f'related_G{G}'
    if label_col not in merged.columns:
        print(f"  ⚠️  Label column {label_col} not found")
        return None

    col_map = {"IBD": "score_IBD", "IBS": "score_IBS", "Phylo": "score_Phylo"}
    fig, ax = plt.subplots(figsize=(7, 6))

    for method_name in inference:
        sc = col_map[method_name]
        if sc not in merged.columns:
            continue

        sub = merged[[sc, label_col]].dropna()
        y_true  = sub[label_col].values.astype(int)
        y_score = sub[sc].values.astype(float)

        if len(np.unique(y_true)) < 2 or len(np.unique(y_score)) < 2:
            print(f"  ⚠️  {method_name}: degenerate at G={G}, skipping")
            continue

        try:
            pr    = pr_curve_with_ci(y_true, y_score)
            auprc = average_precision_score(y_true, y_score)
            ax.plot(pr["recall"], pr["precision_mean"],
                    label=f'{method_name} (AUPRC={auprc:.2f})', linewidth=2)
            ax.fill_between(pr["recall"],
                            pr["precision_low"], pr["precision_high"], alpha=0.2)
        except Exception as e:
            print(f"  ⚠️  {method_name}: {e}")

    ax.set_xlabel("Recall", fontsize=13, fontweight='bold')
    ax.set_ylabel("Precision", fontsize=13, fontweight='bold')
    ax.set_title(f'PR Curves — {scenario_id}  (G={G})', fontsize=16, fontweight='bold')
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
def aggregate_pr_across_scenarios(results_df, method, G=25, recall_grid=None):
    """
    Aggregate PR curves across all scenarios for one method at one G threshold.

    Takes the already-evaluated results_df (output of evaluate_scenario
    concatenated across scenarios) so that:
      - IBD is already evaluated on the full zero-filled universe
      - No raw re-loading or re-merging is needed
      - All methods share identical pair universes (the design invariant)

    Parameters
    ----------
    results_df  : pd.DataFrame — concatenated output of evaluate_scenario()
    method      : str — 'IBD', 'IBS', or 'Phylo'
    G           : int — G threshold to use
    recall_grid : array-like — common recall axis (default: 200 points 0→1)

    Returns
    -------
    dict with keys: recall, precision_mean, precision_low, precision_high
    """
    if recall_grid is None:
        recall_grid = np.linspace(0, 1, 200)

    sub = results_df[
        (results_df['method'] == method) &
        (results_df['G_threshold'] == G)
    ]

    if sub.empty:
        raise ValueError(f"No results for method={method}, G={G}")

    all_precisions = []

    for _, row in sub.iterrows():
        # PR curve data stored in pr_curves CSV; here we reconstruct the
        # scalar AUPRC for aggregation. For full curve aggregation, load
        # the per-scenario *_prcurves.csv files and join on scenario_id.
        auprc = row.get('auprc', np.nan)
        if not np.isnan(auprc):
            # Approximate flat curve at precision=auprc for aggregation
            # when per-point curve data is not available in-memory.
            all_precisions.append(np.full_like(recall_grid, auprc))

    if not all_precisions:
        return None

    arr = np.vstack(all_precisions)
    return {
        "method":          method,
        "G":               G,
        "recall":          recall_grid,
        "precision_mean":  arr.mean(axis=0),
        "precision_low":   np.percentile(arr, 2.5,  axis=0),
        "precision_high":  np.percentile(arr, 97.5, axis=0),
        "n_scenarios":     len(all_precisions),
    }

# ============================================================
# def aggregate_pr_across_scenarios(scenario_dirs, method):
# ============================================================
def aggregate_evaluation_results(df, metrics=("auprc",),
    group_cols=("method", "G_threshold", "sampling_prop", "sample_size"),
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
                row[f"{m}_ci_low"] = np.nan
                row[f"{m}_ci_high"] = np.nan
                continue

            mean = vals.mean()
            sd = vals.std(ddof=1)
            se = sd / np.sqrt(len(vals))

            row[f"{m}_mean"] = mean
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
        on=["sampling_prop", "sample_size"],
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
            fontsize=10,
            color="black",
            ha="center"
        )

    ax.set_xlabel("Sample size", fontsize=13, fontweight='bold')
    ax.set_ylabel("AUPRC", fontsize=13, fontweight='bold')
    ax.set_title("Precision–Recall Performance vs Sampling Density", fontsize=16, fontweight='bold')
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

    axes[0].set_ylabel("AUPRC", fontweight='bold', fontsize=12)
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
    X = results_df[['rec_rate', 'bottleneck', 'est', 'migration', 'sample_size']].copy() # , 'sampling'
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
    ax.set_xlabel('Recombination Rate', fontsize=13, fontweight='bold')
    ax.set_ylabel('Method', fontsize=13, fontweight='bold')
    ax.set_title(f'Method Performance Heatmap ({metric.upper()})',
                fontsize=16, fontweight='bold')

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
    ax.set_xlabel('Recombination Rate (log10)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Sampling Proportion', fontsize=12, fontweight='bold')
    ax.set_title('Identifiability Phase Diagram\n(Dashed line = AUC 0.80 threshold)',
                fontsize=16, fontweight='bold')

    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Best AUC', fontsize=13)

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
    for samp in results_df['sample_size'].unique():
        samp_data = results_df[results_df['sample_size'] == samp]
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

        if i % 50 == 0:
            print(f"  Progress: {i}/{len(scenario_dirs)} scenarios processed")

        results_df = evaluate_scenario(Path(base_dir), scenario_id)

        if results_df is not None:
            all_results.append(results_df)

        # 5. Generate plots
        # scenario = Path(base_dir) / scenario_id
        fig_roc = plot_roc_curves(Path(base_dir), scenario_id)
        fig_pr = plot_pr_curves_single_scenario(Path(base_dir), scenario_id)

        prefix = scenario_id.split('_')[0]
        outdir = output_file.parent / "plots"
        outdir.mkdir(parents=True, exist_ok=True)
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
    import json

    parser = argparse.ArgumentParser(description="Evaluate migration scenarios")
    parser.add_argument("--base-dir",      required=True,
                        help="Directory containing MIG* folders")
    parser.add_argument("--sanity-check",  action="store_true",
                        help="Run sanity check on 3 pilot scenarios")
    parser.add_argument("--full",          action="store_true",
                        help="Run full analysis on all scenarios")
    parser.add_argument("--output",        default="migration_results.csv",
                        help="Output CSV file path")

    args = parser.parse_args()

    if args.sanity_check:
        run_sanity_check(args.base_dir)

    if args.full:
        results  = run_full_analysis(args.base_dir, args.output)
        out_path = Path(args.output)           # ← was: {args.output}.parent (set literal)
        out_dir  = out_path.parent
        out_dir.mkdir(parents=True, exist_ok=True)

        results.to_csv(out_dir / "identifiability_results.csv", index=False)

        # Aggregate — use `results`, not undefined `evaluation_df`
        pr_main = aggregate_evaluation_results(
            results, metrics=("auprc",),
            group_cols=("method", "sample_size")
        )
        pr_by_G = aggregate_evaluation_results(
            results, metrics=("auprc",),
            group_cols=("method", "sample_size", "G_threshold")
        )

        pr_main.to_csv(out_dir / "summary_results.csv",     index=False)
        pr_by_G.to_csv(out_dir / "summary_results_byG.csv", index=False)

        collapse_df = detect_ibd_collapse(pr_main)
        main_plot   = plot_pr_main_panel(pr_main, collapse_df)
        supp        = plot_pr_by_G(pr_by_G)

        answer_research_questions(results)

        heatmap = plot_performance_heatmap(results)
        diagram = plot_identifiability_phase_diagram(results)

        plot_dir = out_dir / "plots"
        plot_dir.mkdir(exist_ok=True)
        main_plot.savefig(plot_dir / "main_figure.png",          dpi=300)
        supp.savefig(     plot_dir / "supp_figure.png",          dpi=300)
        heatmap.savefig(  plot_dir / "performance_heatmap.png",  dpi=300)
        diagram.savefig(  plot_dir / "identifiability.png",      dpi=300)

        effects = analyze_parameter_effects(results)
        print(effects.sort_values('effect_size', ascending=False))

        summary = summarize_identifiability(results)
        print(json.dumps(summary, indent=2))

    if not args.sanity_check and not args.full:
        print("\n" + "="*80)
        print("  Please specify --sanity-check or --full")
        parser.print_help()
        print("=" * 80)


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
