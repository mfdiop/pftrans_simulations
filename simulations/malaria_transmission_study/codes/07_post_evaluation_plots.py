#!/usr/bin/env python3

"""
Evaluation script for migration scenarios with tree-sequence ground truth
Adapted for: true_ibd from tree sequences + generation distance thresholds
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve

def precision_at_k(y_true, y_score, k):
    order = np.argsort(y_score)[::-1]
    top_k = order[:k]
    return np.mean(y_true[top_k])

Choose multiple K (e.g. 50, 100, 500)

Plot Precision@K vs K

Or show a table for key K values

def recall_at_fixed_precision(y_true, y_score, target_precision=0.8):
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    valid = precision >= target_precision
    if not np.any(valid):
        return 0.0
    return np.max(recall[valid])


# ====================================================
# Panel A — AUPRC difference vs recombination
# Goal: Show IBD collapsing as recombination increases
# ====================================================
def plot_auprc_diff_vs_recomb(agg_df):
    """
    Panel A: AUPRC difference (IBD - IBS) vs recombination
    """

    # Pivot to wide
    wide = agg_df.pivot_table(
        index=["scenario_id", "recomb_rate"],
        columns="method",
        values="auprc",
        aggfunc="mean"
    ).reset_index()

    # Keep only scenarios with both methods
    wide = wide.dropna(subset=["IBD", "IBS"])

    wide["delta_auprc"] = wide["IBD"] - wide["IBS"]

    fig, ax = plt.subplots(figsize=(6, 4))

    sns.regplot(
        data=wide,
        x="recomb_rate",
        y="delta_auprc",
        scatter_kws=dict(alpha=0.4),
        line_kws=dict(color="black"),
        ax=ax
    )

    ax.axhline(0, linestyle="--", color="gray", alpha=0.6)
    ax.set_xlabel("Recombination rate")
    ax.set_ylabel("Δ AUPRC (IBD − IBS)")
    ax.set_title("IBD collapses as recombination increases")

    return fig

# ====================================================
# Panel B — Precision@K difference vs prevalence
# Goal: Show IBS robustness when positives are rare
# ====================================================
def plot_precision_at_k_vs_prevalence(agg_df, K=100):
    """
    Panel B: Precision@K difference vs prevalence
    """

    df = agg_df[agg_df["K"] == K].copy()

    wide = df.pivot_table(
        index=["scenario_id", "prevalence"],
        columns="method",
        values="precision_at_k",
        aggfunc="mean"
    ).reset_index()

    wide = wide.dropna(subset=["IBS", "IBD"])
    wide["delta_p_at_k"] = wide["IBS"] - wide["IBD"]

    fig, ax = plt.subplots(figsize=(6, 4))

    sns.scatterplot(
        data=wide,
        x="prevalence",
        y="delta_p_at_k",
        alpha=0.6,
        ax=ax
    )

    sns.lineplot(
        data=wide,
        x="prevalence",
        y="delta_p_at_k",
        estimator="mean",
        ci=95,
        color="black",
        ax=ax
    )

    ax.axhline(0, linestyle="--", color="gray", alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("Prevalence (log scale)")
    ax.set_ylabel(f"Δ Precision@{K} (IBS − IBD)")
    ax.set_title("IBS dominates when true links are rare")

    return fig

# ====================================================
# Panel C — Heatmap (method × regime)
# Goal: Instant visual summary: who wins where
# ====================================================
def assign_regime(row):
    if row["recomb_rate"] < 1e-8 and row["prevalence"] > 0.05:
        return "Low recomb\nHigh prev"
    if row["recomb_rate"] >= 1e-8 and row["prevalence"] <= 0.05:
        return "High recomb\nLow prev"
    return "Mixed"

#
def plot_method_regime_heatmap(agg_df):
    """
    Panel C: Heatmap summarizing method performance by regime
    """

    df = agg_df.copy()
    df["regime"] = df.apply(assign_regime, axis=1)

    summary = (
        df.groupby(["method", "regime"])["auprc"]
          .mean()
          .reset_index()
          .pivot(index="method", columns="regime", values="auprc")
    )

    fig, ax = plt.subplots(figsize=(6, 4))

    sns.heatmap(
        summary,
        annot=True,
        fmt=".2f",
        cmap="RdYlGn",
        center=summary.mean().mean(),
        ax=ax
    )

    ax.set_title("Who wins where: method × regime")
    ax.set_xlabel("Epidemiological regime")
    ax.set_ylabel("Method")

    return fig

# ==========================================================
# Panel D — standardized effect sizes (method × parameter)
# Output of analyze_parameter_effects()
# ==========================================================
def plot_effect_size_heatmap(effects_df):
    pivot = effects_df.pivot(
        index='parameter',
        columns='method',
        values='effect_size')

    sns.heatmap(
        pivot,
        cmap='viridis',
        annot=True,
        fmt=".2f")

    plt.title("Standardized parameter effects on AUPRC")

# ======================================================
def plot_sample_size_scaling(agg_df, metric="auprc"):
    sns.lineplot(data=agg_df,
        x='sample_size', y=f'{metric}_mean',
        hue='method', marker='o'
    )
    plt.ylabel("AUPRC")
    plt.xlabel("Sample size")

# ====================================================
#
# ====================================================
def aggregate_evaluation_results(df, metrics=("auprc",),
    group_cols=("method", "rec_rate_label", "sample_size"),
    ci=0.95, min_replicates=3):

    z_lookup = {0.90: 1.645, 0.95: 1.96, 0.99: 2.576}
    z = z_lookup[ci]

    required_cols = set(group_cols) | set(metrics) | {"replicate_id"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

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
                row.update({
                    f"{m}_mean": np.nan,
                    f"{m}_median": np.nan,
                    f"{m}_ci_low": np.nan,
                    f"{m}_ci_high": np.nan
                })
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

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--phylo-dir", default="sim_migration/phylo_results",
                        help="Directory containing scenario folders")
    parser.add_argument("--outdir", default="phylo_results", help="Output directory")

    args = parser.parse_args()

    results_df = pd.read_csv(file, sep=None, engine="python")
