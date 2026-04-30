"""
Figure 7 — Surveillance constraints on transmission link recoverability
  Panel A: AUPRC vs sample size (bottom axis) and sampling proportion (top axis)
           lines = method, shading = 95% CI
           averaged across rec_rate, migration, bottleneck
  Panel B: AUPRC vs G threshold
           lines = method x rec_rate
           shows how temporal breadth of transmission definition
           interacts with biology

Figure 8 — Joint identifiability limits across biological and surveillance space
  Three contourf panels (one per migration rate)
  x = recombination rate (log10), y = sample size / sampling proportion
  colour = best AUPRC across methods
  dashed contour = AUPRC 0.80 threshold
  This is the synthesis figure: where is identifiability achievable?
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib import cm
import warnings
warnings.filterwarnings("ignore")

# ── data ───────────────────────────────────────────────────────────────────
df  = pd.read_csv("/mnt/user-data/uploads/identifiability_results.csv")

METHOD_COLORS = {"IBD": "#E07B39", "IBS": "#4878CF", "Phylo": "#6ACC65"}
METHOD_LIST   = ["IBD", "IBS", "Phylo"]
REC_VALS      = sorted(df["rec_rate"].unique())
REC_LABELS    = ["1×10⁻⁹", "1×10⁻⁸", "1×10⁻⁷", "1×10⁻⁶"]
MIG_VALS      = sorted(df["migration"].unique())
MIG_LABELS    = ["0.001", "0.010", "0.050"]
SAMP_VALS     = sorted(df["sample_size"].unique())   # [100, 200, 400]
PROP_VALS     = sorted(df["sampling_prop"].unique()) # [0.05, 0.10, 0.20]
G_VALS        = sorted(df["G_threshold"].unique())   # [1,3,5,10,15,25]

REC_LS = {1e-9: "-", 1e-8: "--", 1e-7: "-.", 1e-6: ":"}

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         10,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.grid":         True,
    "grid.alpha":        0.25,
    "grid.linestyle":    "--",
})

OUTPUT = "/mnt/user-data/outputs"

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 7
# ══════════════════════════════════════════════════════════════════════════

fig7, (ax7a, ax7b) = plt.subplots(1, 2, figsize=(14, 5.5))
fig7.subplots_adjust(wspace=0.35)

# ── Panel A: AUPRC vs sample size (dual x-axis: sample size + prop) ───────
# Average across ALL biological parameters to isolate surveillance effect
agg_samp = (df.groupby(["sample_size", "sampling_prop", "method"])["auprc"]
              .agg(mean="mean",
                   ci_lo=lambda x: np.percentile(x, 2.5),
                   ci_hi=lambda x: np.percentile(x, 97.5))
              .reset_index()
              .sort_values("sample_size"))

for method in METHOD_LIST:
    sub = agg_samp[agg_samp["method"] == method]
    ax7a.plot(sub["sample_size"], sub["mean"],
              color=METHOD_COLORS[method],
              linewidth=2.5, marker="o", markersize=7,
              label=method, zorder=3)
    ax7a.fill_between(sub["sample_size"],
                      sub["ci_lo"], sub["ci_hi"],
                      color=METHOD_COLORS[method], alpha=0.15)

ax7a.axhline(0.5, color="gray", linestyle="--",
             linewidth=1, alpha=0.6, label="AUPRC = 0.5")
ax7a.set_xlabel("Sample size  (absolute)", fontsize=10)
ax7a.set_ylabel("Mean AUPRC (95% CI)\naveraged across biological parameters", fontsize=9.5)
ax7a.set_xticks(SAMP_VALS)
ax7a.set_xticklabels(SAMP_VALS, fontsize=9)
ax7a.set_ylim(0.35, 0.75)
ax7a.set_title("A   AUPRC as a function of sampling intensity", fontsize=10,
               loc="left", fontweight="bold")

# Add secondary x-axis for sampling proportion
ax7a_top = ax7a.twiny()
ax7a_top.set_xlim(ax7a.get_xlim())
ax7a_top.set_xticks(SAMP_VALS)
ax7a_top.set_xticklabels([f"{p:.0%}" for p in PROP_VALS], fontsize=9)
ax7a_top.set_xlabel("Sampling proportion  (% of population)", fontsize=9.5)
ax7a_top.spines["top"].set_visible(True)
ax7a_top.spines["right"].set_visible(False)

ax7a.legend(fontsize=9, frameon=False, loc="lower right")

# ── Panel B: AUPRC vs G threshold, by method x rec_rate ──────────────────
# Average across migration, sample_size, bottleneck
agg_G = (df.groupby(["G_threshold", "rec_rate", "method"])["auprc"]
           .agg(mean="mean",
                ci_lo=lambda x: np.percentile(x, 2.5),
                ci_hi=lambda x: np.percentile(x, 97.5))
           .reset_index()
           .sort_values("G_threshold"))

for method in METHOD_LIST:
    for rv in REC_VALS:
        sub = agg_G[(agg_G["method"] == method) &
                    (np.isclose(agg_G["rec_rate"], rv))]
        ax7b.plot(sub["G_threshold"], sub["mean"],
                  color=METHOD_COLORS[method],
                  linestyle=REC_LS[rv],
                  linewidth=1.8, marker="o", markersize=4,
                  alpha=0.85)
        # shading only for IBD to avoid clutter
        if method == "IBD":
            ax7b.fill_between(sub["G_threshold"],
                              sub["ci_lo"], sub["ci_hi"],
                              color=METHOD_COLORS[method], alpha=0.07)

ax7b.axhline(0.5, color="gray", linestyle="--",
             linewidth=1, alpha=0.6)
ax7b.set_xlabel("G threshold (maximum generations since common ancestor)", fontsize=9.5)
ax7b.set_ylabel("Mean AUPRC (95% CI)", fontsize=9.5)
ax7b.set_xticks(G_VALS)
ax7b.set_xticklabels(G_VALS, fontsize=9)
ax7b.set_ylim(0.30, 0.80)
ax7b.set_title("B   AUPRC as a function of transmission window (G threshold)",
               fontsize=10, loc="left", fontweight="bold")

# legend for Panel B
method_h = [Line2D([0], [0], color=METHOD_COLORS[m],
                    linewidth=2, marker="o", markersize=4, label=m)
             for m in METHOD_LIST]
rec_h = [Line2D([0], [0], color="gray",
                 linestyle=REC_LS[rv], linewidth=2,
                 label=f"r = {REC_LABELS[i]}")
          for i, rv in enumerate(REC_VALS)]
ref_h = Line2D([0], [0], color="gray", linestyle="--",
                linewidth=1, alpha=0.6, label="AUPRC = 0.5")
ax7b.legend(handles=method_h + rec_h + [ref_h],
            fontsize=8, frameon=False,
            loc="upper left", ncol=1,
            bbox_to_anchor=(0.01, 0.99))

fig7.suptitle("Figure 7.  Surveillance constraints on transmission link recoverability:\n"
              "effect of sampling intensity and temporal transmission window definition",
              fontsize=11, y=1.02)

plt.savefig(f"{OUTPUT}/Figure7_surveillance_constraints.pdf",
            bbox_inches="tight", dpi=300)
plt.savefig(f"{OUTPUT}/Figure7_surveillance_constraints.png",
            bbox_inches="tight", dpi=180)
plt.close()
print("Figure 7 saved.")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 8 — Joint identifiability limits
# Three contourf panels (migration), x=rec_rate, y=sample_size
# Colour = best AUPRC across methods at G=25
# Dashed = AUPRC 0.80 contour
# Right panel: best AUPRC at G=1 (narrow window) vs G=25 (broad window)
#              to show how G interacts with the joint space
# ══════════════════════════════════════════════════════════════════════════

g25 = df[df["G_threshold"] == 25].copy()
g1  = df[df["G_threshold"] == 1].copy()

best25 = (g25.groupby(["rec_rate", "migration", "sample_size"])["auprc"]
             .max().reset_index().rename(columns={"auprc": "best_auprc"}))
best1  = (g1.groupby(["rec_rate", "migration", "sample_size"])["auprc"]
            .max().reset_index().rename(columns={"auprc": "best_auprc"}))

REC_LOG  = np.log10(REC_VALS)
cmap8    = plt.cm.RdYlGn
norm8    = Normalize(vmin=0.45, vmax=1.0)
levels8  = np.linspace(0.45, 1.0, 300)

# Layout: 2 rows x 3 cols
# Row 0: G=25, one panel per migration (3 panels)
# Row 1: G=1,  one panel per migration (3 panels)
# Right: shared colorbar
fig8 = plt.figure(figsize=(16, 10))
gs8  = gridspec.GridSpec(2, 3, figure=fig8, hspace=0.42, wspace=0.28)

axes_top = [fig8.add_subplot(gs8[0, k]) for k in range(3)]
axes_bot = [fig8.add_subplot(gs8[1, k]) for k in range(3)]

def phase_panel(ax, data, mig, title, show_ylabel, show_xlabel):
    sub = data[np.isclose(data["migration"], mig)]
    Z   = np.full((len(SAMP_VALS), len(REC_VALS)), np.nan)
    for i, sv in enumerate(SAMP_VALS):
        for j, rv in enumerate(REC_VALS):
            row = sub[(sub["sample_size"] == sv) &
                      (np.isclose(sub["rec_rate"], rv))]
            if len(row):
                Z[i, j] = row["best_auprc"].values[0]

    cf = ax.contourf(REC_LOG, SAMP_VALS, Z,
                     levels=levels8, cmap=cmap8, norm=norm8)

    # annotate cells
    for i, sv in enumerate(SAMP_VALS):
        for j, rl in enumerate(REC_LOG):
            val = Z[i, j]
            if not np.isnan(val):
                ax.text(rl, sv, f"{val:.2f}",
                        ha="center", va="center",
                        fontsize=8.5,
                        color="white" if val > 0.78 else "black")

    # 0.80 contour
    try:
        ax.contour(REC_LOG, SAMP_VALS, Z,
                   levels=[0.80],
                   colors="black", linestyles="--", linewidths=2)
    except Exception:
        pass

    ax.set_title(title, fontsize=9.5)
    ax.set_xticks(REC_LOG)
    if show_xlabel:
        ax.set_xticklabels(REC_LABELS, fontsize=8, rotation=20)
        ax.set_xlabel("Recombination rate", fontsize=9)
    else:
        ax.set_xticklabels([])
    ax.set_yticks(SAMP_VALS)
    if show_ylabel:
        ax.set_yticklabels(SAMP_VALS, fontsize=8.5)
        ax.set_ylabel("Sample size", fontsize=9)
    else:
        ax.set_yticklabels([])

    # add sampling proportion on right y-axis for leftmost panels
    if show_ylabel:
        ax2 = ax.twinx() if not show_ylabel else None

    return cf

# Row 0: G=25
for k, (ax, mig) in enumerate(zip(axes_top, MIG_VALS)):
    cf = phase_panel(ax, best25, mig,
                     title=f"G ≤ 25 generations  |  Migration = {mig}",
                     show_ylabel=(k == 0),
                     show_xlabel=False)

# Row 1: G=1
for k, (ax, mig) in enumerate(zip(axes_bot, MIG_VALS)):
    cf = phase_panel(ax, best1, mig,
                     title=f"G = 1 generation  |  Migration = {mig}",
                     show_ylabel=(k == 0),
                     show_xlabel=True)

# Row labels
fig8.text(0.005, 0.75, "Broad\ntransmission\nwindow\n(G ≤ 25)",
          fontsize=9, va="center", ha="left",
          rotation=0, color="#555555",
          bbox=dict(boxstyle="round,pad=0.3", fc="#f0f0f0", ec="none"))
fig8.text(0.005, 0.28, "Narrow\ntransmission\nwindow\n(G = 1)",
          fontsize=9, va="center", ha="left",
          rotation=0, color="#555555",
          bbox=dict(boxstyle="round,pad=0.3", fc="#f0f0f0", ec="none"))

# Panel labels
for k, ax in enumerate(axes_top):
    ax.text(-0.10, 1.07, chr(65 + k),
            transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")
for k, ax in enumerate(axes_bot):
    ax.text(-0.10, 1.07, chr(68 + k),
            transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")

# Dashed contour legend
contour_h = Line2D([0], [0], color="black", linestyle="--",
                    linewidth=2, label="AUPRC = 0.80 threshold")
fig8.legend(handles=[contour_h], loc="lower center",
            fontsize=9, frameon=False,
            bbox_to_anchor=(0.5, -0.02), ncol=1)

# Shared colorbar
cbar8 = fig8.colorbar(
    cm.ScalarMappable(norm=norm8, cmap=cmap8),
    ax=list(axes_top) + list(axes_bot),
    orientation="vertical",
    fraction=0.012, pad=0.03, shrink=0.85
)
cbar8.set_label("Best AUPRC  (max over IBD, IBS, Phylo)", fontsize=9.5)
cbar8.ax.tick_params(labelsize=8.5)

# Sampling proportion secondary tick annotation on colorbar side
fig8.text(0.97, 0.75, "Sampling\nproportion:", fontsize=8, ha="center",
          color="#555555")
for sv, pv in zip(SAMP_VALS, PROP_VALS):
    # locate approximate y position
    pass  # handled by dual-axis in phase_panel if needed

fig8.suptitle("Figure 8.  Joint identifiability limits across recombination rate,\n"
              "migration rate, and surveillance sampling intensity",
              fontsize=11, y=1.01)

plt.savefig(f"{OUTPUT}/Figure8_joint_identifiability.pdf",
            bbox_inches="tight", dpi=300)
plt.savefig(f"{OUTPUT}/Figure8_joint_identifiability.png",
            bbox_inches="tight", dpi=180)
plt.close()
print("Figure 8 saved.")
print("\nAll done.")
