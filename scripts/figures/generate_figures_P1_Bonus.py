"""
Figure generation script for all missing/priority figures.
Figures produced:
  Priority 1  → Fig_P1_migration_auprc_gradient.pdf
  Priority 2  → Fig_P2_phase_diagram_extended.pdf
  Priority 3  → Fig_P3_bottleneck_auprc.pdf
  Priority 4  → Fig_P4_ibd_method_gap_heatmap.pdf  (IBD vs next-best gap)
  Missing 1   → Fig_M1_sensitivity_at_90spec_migration.pdf  (proxy for IBD decay)
  Missing 2   → Fig_M2_identifiable_fraction_migration.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib import cm
import warnings
warnings.filterwarnings("ignore")

# ── load data ──────────────────────────────────────────────────────────────
df = pd.read_csv("/mnt/user-data/uploads/identifiability_results.csv")

# Ordered labels for axes
REC_ORDER   = [1e-9, 1e-8, 1e-7, 1e-6]
REC_LABELS  = ["1e-09", "1e-08", "1e-07", "1e-06"]
MIG_ORDER   = [0.001, 0.010, 0.050]
MIG_LABELS  = ["0.001", "0.010", "0.050"]
METHOD_COLS = {"IBD": "#E07B39", "IBS": "#4878CF", "Phylo": "#6ACC65"}
METHOD_LIST = ["IBD", "IBS", "Phylo"]

# ── shared style ───────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})

OUTPUT = "/mnt/user-data/outputs"

# ══════════════════════════════════════════════════════════════════════════
# PRIORITY 1 — Migration × recombination AUPRC gradient  (Fig_P1)
# Line plot: x = recombination rate (log), y = mean AUPRC ± 95 CI
# Facets = migration rate (3 panels), lines = method
# G_threshold = 25 (broadest, most interpretable)
# ══════════════════════════════════════════════════════════════════════════
g25 = df[df["G_threshold"] == 25].copy()

agg = (g25.groupby(["migration", "rec_rate", "method"])["auprc"]
         .agg(mean="mean", std="std", n="count")
         .reset_index())
agg["ci"] = 1.96 * agg["std"] / np.sqrt(agg["n"])

fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
fig.suptitle("AUPRC by Recombination Rate and Migration Rate (G ≤ 25 generations)",
             fontsize=11, y=1.02)

for ax, mig in zip(axes, MIG_ORDER):
    sub = agg[agg["migration"] == mig]
    for method in METHOD_LIST:
        m = sub[sub["method"] == method].sort_values("rec_rate")
        ax.plot(np.log10(m["rec_rate"]), m["mean"],
                color=METHOD_COLS[method], label=method,
                marker="o", linewidth=2, markersize=5)
        ax.fill_between(np.log10(m["rec_rate"]),
                        m["mean"] - m["ci"], m["mean"] + m["ci"],
                        color=METHOD_COLS[method], alpha=0.15)
    ax.axhline(0.6, color="gray", linestyle="--", linewidth=1,
               label="AUPRC = 0.6" if mig == MIG_ORDER[0] else "")
    ax.set_title(f"Migration = {mig}", fontsize=10)
    ax.set_xlabel("Recombination rate (log₁₀)", fontsize=9)
    ax.set_xticks([-9, -8, -7, -6])
    ax.set_xticklabels(["1e-9", "1e-8", "1e-7", "1e-6"], fontsize=8)

axes[0].set_ylabel("Mean AUPRC (± 95% CI)", fontsize=9)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles[:4], labels[:4], loc="lower center",
           ncol=4, frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, -0.08))
plt.tight_layout()
plt.savefig(f"{OUTPUT}/Fig_P1_migration_auprc_gradient.pdf",
            bbox_inches="tight", dpi=300)
plt.close()
print("Saved Fig_P1")

# ══════════════════════════════════════════════════════════════════════════
# PRIORITY 2 — Extended phase diagram: recombination × sample size
#              faceted by migration rate (3 panels)
#              colour = best AUPRC (max over methods) at G=25
# ══════════════════════════════════════════════════════════════════════════
best_auprc = (g25.groupby(["rec_rate", "migration", "sample_size"])["auprc"]
                 .max()
                 .reset_index()
                 .rename(columns={"auprc": "best_auprc"}))

rec_vals  = np.log10(sorted(df["rec_rate"].unique()))
samp_vals = sorted(df["sample_size"].unique())

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
fig.suptitle("Identifiability Phase Diagram by Migration Rate\n"
             "(Best AUPRC across methods, G ≤ 25, dashed = AUC 0.80 threshold)",
             fontsize=11, y=1.03)

cmap   = plt.cm.RdYlGn
levels = np.linspace(0.45, 1.0, 200)
norm   = plt.Normalize(vmin=0.45, vmax=1.0)

for ax, mig in zip(axes, MIG_ORDER):
    sub = best_auprc[best_auprc["migration"] == mig]
    Z = np.zeros((len(samp_vals), len(rec_vals)))
    for i, sv in enumerate(samp_vals):
        for j, rv in enumerate(10**rec_vals):
            row = sub[(sub["sample_size"] == sv) &
                      (np.isclose(sub["rec_rate"], rv))]
            Z[i, j] = row["best_auprc"].values[0] if len(row) else np.nan

    cf = ax.contourf(rec_vals, samp_vals, Z,
                     levels=levels, cmap=cmap, norm=norm)
    # Draw 0.80 contour line
    try:
        ax.contour(rec_vals, samp_vals, Z,
                   levels=[0.80], colors="black",
                   linestyles="--", linewidths=1.5)
    except Exception:
        pass

    ax.set_title(f"Migration = {mig}", fontsize=10)
    ax.set_xlabel("Recombination rate (log₁₀)", fontsize=9)
    ax.set_xticks([-9, -8, -7, -6])
    ax.set_xticklabels(["1e-9", "1e-8", "1e-7", "1e-6"], fontsize=8)

axes[0].set_ylabel("Sample size", fontsize=9)
cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                    ax=axes, orientation="vertical",
                    fraction=0.015, pad=0.02)
cbar.set_label("Best AUPRC", fontsize=9)
plt.savefig(f"{OUTPUT}/Fig_P2_phase_diagram_extended.pdf",
            bbox_inches="tight", dpi=300)
plt.close()
print("Saved Fig_P2")

# ══════════════════════════════════════════════════════════════════════════
# PRIORITY 3 — Bottleneck effect on AUPRC
# Design note: bottleneck=5 & bottleneck=20 both appear at rec_rate=medium
#              bottleneck=1 only at verylow rec_rate
#              → compare 5 vs 20 at medium rec_rate; show verylow separately
# Facets = rec_rate (medium only for direct comparison, verylow for context)
# x = bottleneck, y = AUPRC, lines = method, colour = method
# G = 25 for consistency
# ══════════════════════════════════════════════════════════════════════════
agg_bn = (g25.groupby(["rec_rate_label", "rec_rate", "bottleneck", "method"])["auprc"]
             .agg(mean="mean", std="std", n="count")
             .reset_index())
agg_bn["ci"] = 1.96 * agg_bn["std"] / np.sqrt(agg_bn["n"])

# Use rec_rates that have >1 bottleneck value (medium has 5 & 20)
# Also show verylow (bottleneck=1) as context bar
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
fig.suptitle("Effect of Bottleneck Size on AUPRC by Method (G ≤ 25 generations)",
             fontsize=11, y=1.02)

for ax, (label, title) in zip(
        axes,
        [("medium", "Recombination rate = 1×10⁻⁷ (medium)"),
         ("verylow", "Recombination rate = 1×10⁻⁹ (very low)")]):
    sub = agg_bn[agg_bn["rec_rate_label"] == label]
    bns = sorted(sub["bottleneck"].unique())
    x = np.arange(len(bns))
    width = 0.25
    for k, method in enumerate(METHOD_LIST):
        m = sub[sub["method"] == method].sort_values("bottleneck")
        ax.bar(x + k * width, m["mean"], width,
               color=METHOD_COLS[method], alpha=0.85,
               label=method, yerr=m["ci"], capsize=3,
               error_kw={"elinewidth": 1})
    ax.set_xticks(x + width)
    ax.set_xticklabels([f"Ne bottleneck\n= {b}×1000" for b in bns], fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Effective population size at bottleneck", fontsize=9)
    ax.axhline(0.6, color="gray", linestyle="--", linewidth=1)

axes[0].set_ylabel("Mean AUPRC (± 95% CI)", fontsize=9)
handles, labels_ = axes[0].get_legend_handles_labels()
fig.legend(handles[:3], labels_[:3], loc="lower center",
           ncol=3, frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, -0.08))
plt.tight_layout()
plt.savefig(f"{OUTPUT}/Fig_P3_bottleneck_auprc.pdf",
            bbox_inches="tight", dpi=300)
plt.close()
print("Saved Fig_P3")

# ══════════════════════════════════════════════════════════════════════════
# PRIORITY 4 — IBD vs next-best method gap heatmap
# For each (rec_rate, migration) cell: IBD_auprc - max(IBS_auprc, Phylo_auprc)
# Colour = gap magnitude; red = IBD much better; white = convergence
# G = 25, averaged across sample sizes and replicates
# ══════════════════════════════════════════════════════════════════════════
pivot = (g25.groupby(["rec_rate", "migration", "method"])["auprc"]
            .mean()
            .unstack("method")
            .reset_index())
pivot["next_best"] = pivot[["IBS", "Phylo"]].max(axis=1)
pivot["gap"]       = pivot["IBD"] - pivot["next_best"]

gap_mat = pivot.pivot_table(index="migration", columns="rec_rate",
                            values="gap")

fig, ax = plt.subplots(figsize=(7, 4))
im = ax.imshow(gap_mat.values, cmap="RdYlGn_r",
               vmin=0, vmax=0.15, aspect="auto")
ax.set_xticks(range(len(REC_ORDER)))
ax.set_xticklabels(REC_LABELS, fontsize=9)
ax.set_yticks(range(len(MIG_ORDER)))
ax.set_yticklabels(MIG_LABELS, fontsize=9)
ax.set_xlabel("Recombination rate", fontsize=10)
ax.set_ylabel("Migration rate", fontsize=10)
ax.set_title("IBD advantage over next-best method\n(AUPRC gap, G ≤ 25)",
             fontsize=11)
for i in range(len(MIG_ORDER)):
    for j in range(len(REC_ORDER)):
        ax.text(j, i, f"{gap_mat.values[i, j]:.3f}",
                ha="center", va="center", fontsize=9,
                color="black" if gap_mat.values[i, j] < 0.10 else "white")
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("AUPRC gap (IBD − next best)", fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUTPUT}/Fig_P4_ibd_method_gap_heatmap.pdf",
            bbox_inches="tight", dpi=300)
plt.close()
print("Saved Fig_P4")

# ══════════════════════════════════════════════════════════════════════════
# MISSING 1 — Sensitivity at 90% specificity as proxy for IBD signal
#             under migration (closest available metric to IBD decay plot)
# x = recombination rate (log), y = mean sensitivity@90spec
# Facets = migration rate, lines = method, G = 25
# ══════════════════════════════════════════════════════════════════════════
agg_sens = (g25.groupby(["migration", "rec_rate", "method"])["sensitivity_at_90spec"]
               .agg(mean="mean", std="std", n="count")
               .reset_index())
agg_sens["ci"] = 1.96 * agg_sens["std"] / np.sqrt(agg_sens["n"])

fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
fig.suptitle("Sensitivity at 90% Specificity by Recombination and Migration Rate (G ≤ 25)",
             fontsize=11, y=1.02)

for ax, mig in zip(axes, MIG_ORDER):
    sub = agg_sens[agg_sens["migration"] == mig]
    for method in METHOD_LIST:
        m = sub[sub["method"] == method].sort_values("rec_rate")
        ax.plot(np.log10(m["rec_rate"]), m["mean"],
                color=METHOD_COLS[method], label=method,
                marker="o", linewidth=2, markersize=5)
        ax.fill_between(np.log10(m["rec_rate"]),
                        m["mean"] - m["ci"], m["mean"] + m["ci"],
                        color=METHOD_COLS[method], alpha=0.15)
    ax.set_title(f"Migration = {mig}", fontsize=10)
    ax.set_xlabel("Recombination rate (log₁₀)", fontsize=9)
    ax.set_xticks([-9, -8, -7, -6])
    ax.set_xticklabels(["1e-9", "1e-8", "1e-7", "1e-6"], fontsize=8)

axes[0].set_ylabel("Sensitivity at 90% specificity (± 95% CI)", fontsize=9)
handles, labels_ = axes[0].get_legend_handles_labels()
fig.legend(handles[:3], labels_[:3], loc="lower center",
           ncol=3, frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, -0.08))
plt.tight_layout()
plt.savefig(f"{OUTPUT}/Fig_M1_sensitivity_90spec_migration.pdf",
            bbox_inches="tight", dpi=300)
plt.close()
print("Saved Fig_M1")

# ══════════════════════════════════════════════════════════════════════════
# MISSING 2 — Identifiable fraction by migration x rec_rate x method
# Heatmap: rows = method, columns = rec_rate, facets = migration
# colour = fraction of replicates where identifiable == True (G=25)
# ══════════════════════════════════════════════════════════════════════════
ident_frac = (g25.groupby(["migration", "rec_rate", "method"])["identifiable"]
                 .mean()
                 .reset_index()
                 .rename(columns={"identifiable": "frac_identifiable"}))

fig, axes = plt.subplots(1, 3, figsize=(13, 3.5), sharey=True)
fig.suptitle("Fraction of Replicates Classified as Identifiable (G ≤ 25)",
             fontsize=11, y=1.04)

for ax, mig in zip(axes, MIG_ORDER):
    sub = ident_frac[ident_frac["migration"] == mig]
    mat = sub.pivot_table(index="method", columns="rec_rate",
                          values="frac_identifiable")
    mat = mat.reindex(METHOD_LIST)
    im = ax.imshow(mat.values, cmap="RdYlGn",
                   vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(REC_ORDER)))
    ax.set_xticklabels(REC_LABELS, fontsize=8, rotation=30)
    ax.set_yticks(range(len(METHOD_LIST)))
    ax.set_yticklabels(METHOD_LIST, fontsize=9)
    ax.set_title(f"Migration = {mig}", fontsize=10)
    ax.set_xlabel("Recombination rate", fontsize=9)
    for i in range(len(METHOD_LIST)):
        for j in range(len(REC_ORDER)):
            val = mat.values[i, j]
            ax.text(j, i, f"{val:.2f}",
                    ha="center", va="center", fontsize=8,
                    color="black" if val < 0.75 else "white")

cbar = fig.colorbar(im, ax=axes, orientation="vertical",
                    fraction=0.015, pad=0.02)
cbar.set_label("Fraction identifiable", fontsize=9)
plt.savefig(f"{OUTPUT}/Fig_M2_identifiable_fraction_migration.pdf",
            bbox_inches="tight", dpi=300)
plt.close()
print("Saved Fig_M2")

# ══════════════════════════════════════════════════════════════════════════
# BONUS — G threshold × rec_rate interaction for all methods
# Shows how temporal window definition interacts with recombination
# Heatmap per method: rows = G, columns = rec_rate, colour = mean AUPRC
# averaged across migration and sample_size
# ══════════════════════════════════════════════════════════════════════════
agg_G = (df.groupby(["G_threshold", "rec_rate", "method"])["auprc"]
           .mean()
           .reset_index())

G_ORDER = [1, 3, 5, 10, 15, 25]

fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
fig.suptitle("AUPRC by G Threshold and Recombination Rate\n"
             "(mean across migration and sample size)",
             fontsize=11, y=1.04)

for ax, method in zip(axes, METHOD_LIST):
    sub = agg_G[agg_G["method"] == method]
    mat = sub.pivot_table(index="G_threshold", columns="rec_rate",
                          values="auprc")
    mat = mat.reindex(G_ORDER)
    im = ax.imshow(mat.values, cmap="RdYlGn",
                   vmin=0.25, vmax=0.75, aspect="auto")
    ax.set_xticks(range(len(REC_ORDER)))
    ax.set_xticklabels(REC_LABELS, fontsize=8, rotation=30)
    ax.set_yticks(range(len(G_ORDER)))
    ax.set_yticklabels([f"G={g}" for g in G_ORDER], fontsize=9)
    ax.set_title(method, fontsize=11, fontweight="bold")
    ax.set_xlabel("Recombination rate", fontsize=9)
    for i in range(len(G_ORDER)):
        for j in range(len(REC_ORDER)):
            val = mat.values[i, j]
            ax.text(j, i, f"{val:.2f}",
                    ha="center", va="center", fontsize=7.5,
                    color="black" if val < 0.60 else "white")

cbar = fig.colorbar(im, ax=axes, orientation="vertical",
                    fraction=0.015, pad=0.02)
cbar.set_label("Mean AUPRC", fontsize=9)
axes[0].set_ylabel("G threshold (generations)", fontsize=9)
plt.savefig(f"{OUTPUT}/Fig_BONUS_G_recrate_interaction.pdf",
            bbox_inches="tight", dpi=300)
plt.close()
print("Saved Fig_BONUS_G_recrate_interaction")

print("\nAll figures saved to /mnt/user-data/outputs/")
