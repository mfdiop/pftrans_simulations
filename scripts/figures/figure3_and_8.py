"""
Figure 3 (redesigned): Comparative recoverability across genomic inference layers
               with migration as a continuous gradient across the full factorial

Figure 8 (new):        G threshold × recombination × method interaction
                       — how surveillance temporal window interacts with
                         biology to define identifiability limits
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, to_rgba
from matplotlib import cm
import warnings
warnings.filterwarnings("ignore")

# ── data ───────────────────────────────────────────────────────────────────
df = pd.read_csv("/mnt/user-data/uploads/identifiability_results.csv")
df["sensitivity"] = df["TP"] / (df["TP"] + df["FN"])
df["FDR"]         = np.where(
    (df["TP"] + df["FP"]) > 0,
    df["FP"] / (df["TP"] + df["FP"]),
    np.nan
)

# ── shared style ───────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "font.size":          9.5,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.alpha":         0.25,
    "grid.linestyle":     "--",
    "axes.labelsize":     9.5,
    "xtick.labelsize":    8.5,
    "ytick.labelsize":    8.5,
})

METHOD_COLORS = {"IBD": "#E07B39", "IBS": "#4878CF", "Phylo": "#6ACC65"}
METHOD_LIST   = ["IBD", "IBS", "Phylo"]

REC_VALS      = sorted(df["rec_rate"].unique())
REC_LABELS    = ["1×10⁻⁹", "1×10⁻⁸", "1×10⁻⁷", "1×10⁻⁶"]
MIG_VALS      = sorted(df["migration"].unique())
MIG_LABELS    = ["0.001", "0.010", "0.050"]

# Migration colormap for gradient (light → dark blue)
MIG_CMAP   = cm.Blues
MIG_NORM   = Normalize(vmin=-0.5, vmax=len(MIG_VALS) - 0.3)
MIG_COLORS = {m: MIG_CMAP(MIG_NORM(i + 0.7)) for i, m in enumerate(MIG_VALS)}

OUTPUT = "/mnt/user-data/outputs"

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3 (redesigned)
# Layout:
#   Row 0 (3 panels): AUPRC distributions — violin + jitter
#                     x = migration rate (3 positions), colour = method
#                     facets = recombination rate (4 panels in row)  ← 4 panels
#   Row 1 (3 panels): Sensitivity@90spec same structure
#   Row 2 (1 panel):  AUPRC gradient lines — x=rec_rate, colour=method,
#                     linetype=migration, averaged across sample_size/bottleneck
#   Row 3 (1 panel):  Proportion of replicates where IBD > IBS > Phylo
#                     shown as stacked bar or line, x=rec_rate, facet=migration
# ══════════════════════════════════════════════════════════════════════════

g25 = df[df["G_threshold"] == 25].copy()

# ── pre-aggregate ──────────────────────────────────────────────────────────
agg_auprc = (g25.groupby(["rec_rate", "migration", "method"])
               .agg(mean=("auprc", "mean"),
                    median=("auprc", "median"),
                    q25=("auprc", lambda x: np.percentile(x, 25)),
                    q75=("auprc", lambda x: np.percentile(x, 75)),
                    ci_lo=("auprc", lambda x: np.percentile(x, 2.5)),
                    ci_hi=("auprc", lambda x: np.percentile(x, 97.5)))
               .reset_index())

agg_sens = (g25.groupby(["rec_rate", "migration", "method"])
               .agg(mean=("sensitivity_at_90spec", "mean"),
                    q25=("sensitivity_at_90spec", lambda x: np.percentile(x, 25)),
                    q75=("sensitivity_at_90spec", lambda x: np.percentile(x, 75)),
                    ci_lo=("sensitivity_at_90spec", lambda x: np.percentile(x, 2.5)),
                    ci_hi=("sensitivity_at_90spec", lambda x: np.percentile(x, 97.5)))
               .reset_index())

# ── Figure 3 layout ────────────────────────────────────────────────────────
fig3 = plt.figure(figsize=(16, 13))
gs3  = gridspec.GridSpec(3, 4,
                          figure=fig3,
                          hspace=0.52, wspace=0.28,
                          height_ratios=[1, 1, 1.1])

# Row 0: AUPRC by method x migration, faceted by rec_rate
axes_auprc = [fig3.add_subplot(gs3[0, k]) for k in range(4)]
# Row 1: Sensitivity@90spec
axes_sens  = [fig3.add_subplot(gs3[1, k]) for k in range(4)]
# Row 2: gradient line plot spanning all 4 columns
ax_grad    = fig3.add_subplot(gs3[2, :])

panel_labels = list("ABCDEFGHIJ")
label_idx = 0

def add_panel_label(ax, label):
    ax.text(-0.14, 1.08, label, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")

# ── helper: grouped bar with CI ──────────────────────────────────────────
def grouped_bar_ci(ax, data_agg, metric_col, title, ylabel, show_xlab,
                   show_ref=False):
    n_mig     = len(MIG_VALS)
    n_method  = len(METHOD_LIST)
    group_w   = 0.8
    bar_w     = group_w / n_method
    x_pos     = np.arange(n_mig)

    for k, method in enumerate(METHOD_LIST):
        sub = (data_agg[data_agg["method"] == method]
               .groupby("migration", as_index=False)[["migration","mean","ci_lo","ci_hi"]]
               .mean()
               .sort_values("migration"))
        sub = sub.set_index("migration").reindex(MIG_VALS).reset_index()
        offset = (k - (n_method - 1) / 2) * bar_w
        bars = ax.bar(x_pos + offset,
                      sub["mean"],
                      bar_w * 0.88,
                      color=METHOD_COLORS[method],
                      alpha=0.85,
                      label=method)
        ax.errorbar(x_pos + offset,
                    sub["mean"],
                    yerr=[sub["mean"] - sub["ci_lo"],
                          sub["ci_hi"] - sub["mean"]],
                    fmt="none",
                    color="black",
                    capsize=2.5,
                    linewidth=0.8,
                    alpha=0.7)

    if show_ref:
        ax.axhline(0.5, color="gray", linestyle="--",
                   linewidth=1, alpha=0.6, zorder=0)

    ax.set_xticks(x_pos)
    if show_xlab:
        ax.set_xticklabels(MIG_LABELS, fontsize=8)
        ax.set_xlabel("Migration rate", fontsize=8.5)
    else:
        ax.set_xticklabels([])
    ax.set_title(title, fontsize=9, pad=4)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=8.5)
    ax.set_ylim(0, 0.85)

# ── Row 0: AUPRC ──────────────────────────────────────────────────────────
for k, (ax, rv) in enumerate(zip(axes_auprc, REC_VALS)):
    sub = agg_auprc[np.isclose(agg_auprc["rec_rate"], rv)]
    grouped_bar_ci(ax, sub, "mean",
                   title=f"r = {REC_LABELS[k]}",
                   ylabel="Mean AUPRC (95% CI)" if k == 0 else "",
                   show_xlab=False,
                   show_ref=True)
    add_panel_label(ax, panel_labels[label_idx]); label_idx += 1

# ── Row 1: Sensitivity@90spec ─────────────────────────────────────────────
for k, (ax, rv) in enumerate(zip(axes_sens, REC_VALS)):
    sub = agg_sens[np.isclose(agg_sens["rec_rate"], rv)]
    grouped_bar_ci(ax, sub, "mean",
                   title=f"r = {REC_LABELS[k]}",
                   ylabel="Sensitivity @ 90% spec (95% CI)" if k == 0 else "",
                   show_xlab=True,
                   show_ref=False)
    add_panel_label(ax, panel_labels[label_idx]); label_idx += 1

# ── Row 2: gradient line plot ─────────────────────────────────────────────
# x = recombination rate, y = mean AUPRC
# solid/dashed/dotted line = migration rate, colour = method
LINESTYLS = {0.001: "-", 0.010: "--", 0.050: ":"}
LINEWIDTHS = {0.001: 2.2, 0.010: 1.8, 0.050: 1.5}

for method in METHOD_LIST:
    for mig in MIG_VALS:
        sub = (agg_auprc[(agg_auprc["method"] == method) &
                          (np.isclose(agg_auprc["migration"], mig))]
               .sort_values("rec_rate"))
        ax_grad.plot(np.log10(sub["rec_rate"]),
                     sub["mean"],
                     color=METHOD_COLORS[method],
                     linestyle=LINESTYLS[mig],
                     linewidth=LINEWIDTHS[mig],
                     marker="o", markersize=5,
                     alpha=0.9)
        ax_grad.fill_between(np.log10(sub["rec_rate"]),
                              sub["ci_lo"], sub["ci_hi"],
                              color=METHOD_COLORS[method],
                              alpha=0.06)

ax_grad.axhline(0.5, color="gray", linestyle="--",
                linewidth=1, alpha=0.5, label="AUPRC = 0.5")
ax_grad.set_xticks(np.log10(REC_VALS))
ax_grad.set_xticklabels(REC_LABELS, fontsize=9)
ax_grad.set_xlabel("Recombination rate", fontsize=10)
ax_grad.set_ylabel("Mean AUPRC (95% CI)", fontsize=10)
ax_grad.set_title("AUPRC gradient across recombination rate, migration rate, and inference layer  (G ≤ 25)",
                  fontsize=10)
add_panel_label(ax_grad, panel_labels[label_idx]); label_idx += 1

# ── legends ────────────────────────────────────────────────────────────────
method_handles = [mpatches.Patch(color=METHOD_COLORS[m],
                                  label=m, alpha=0.85)
                  for m in METHOD_LIST]
mig_handles = [Line2D([0], [0],
                       color="gray",
                       linestyle=LINESTYLS[m],
                       linewidth=2,
                       label=f"Migration = {m}")
               for m in MIG_VALS]
ref_handle = Line2D([0], [0], color="gray", linestyle="--",
                     linewidth=1, alpha=0.6, label="AUPRC = 0.5")

fig3.legend(handles=method_handles + mig_handles + [ref_handle],
            loc="lower center", ncol=7,
            frameon=False, fontsize=9,
            bbox_to_anchor=(0.5, -0.03))

fig3.suptitle("Figure 3.  Comparative recoverability of transmission links across\n"
              "genomic inference layers, recombination rates, and migration rates",
              fontsize=11, y=1.01)

plt.savefig(f"{OUTPUT}/Figure3_redesigned.pdf",
            bbox_inches="tight", dpi=300)
plt.savefig(f"{OUTPUT}/Figure3_redesigned.png",
            bbox_inches="tight", dpi=180)
plt.close()
print("Figure 3 saved.")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 8 — G threshold × recombination × method interaction
#
# This is the appropriate Figure 8 because:
# 1. It is the only remaining major result not yet in main text
# 2. It directly addresses a surveillance constraint dimension:
#    how does the temporal breadth of the transmission definition
#    interact with biology to set identifiability limits?
# 3. It closes the narrative: biology (rec, migration) constrains the ceiling;
#    surveillance design (G, sample size) determines where within that ceiling
#    you operate
#
# Layout:
#   Panel A (left 3): heatmap per method — rows=G, cols=rec_rate
#                     colour = mean AUPRC (averaged across migration, sample_size)
#   Panel B (right):  line plot — x=G, y=mean AUPRC,
#                     colour=method, linetype=rec_rate
#                     Shows the G-threshold gain curve per method and rec_rate
# ══════════════════════════════════════════════════════════════════════════

G_VALS      = sorted(df["G_threshold"].unique())   # [1,3,5,10,15,25]
G_LABELS    = [f"G={g}" for g in G_VALS]

# Mean AUPRC by G x rec_rate x method (averaged across migration, sample_size)
agg_G = (df.groupby(["G_threshold", "rec_rate", "method"])["auprc"]
           .agg(mean="mean",
                ci_lo=lambda x: np.percentile(x, 2.5),
                ci_hi=lambda x: np.percentile(x, 97.5))
           .reset_index())

# Colormap for rec_rate in line plot
REC_CMAP   = cm.plasma
REC_NORM   = Normalize(vmin=-0.5, vmax=len(REC_VALS) - 0.3)
REC_COLORS = {rv: REC_CMAP(REC_NORM(i + 0.5)) for i, rv in enumerate(REC_VALS)}
REC_LS     = {1e-9: "-", 1e-8: "--", 1e-7: "-.", 1e-6: ":"}

fig8 = plt.figure(figsize=(16, 7))
gs8  = gridspec.GridSpec(1, 4,
                          figure=fig8,
                          wspace=0.32)

axes_heat = [fig8.add_subplot(gs8[0, k]) for k in range(3)]
ax_line   = fig8.add_subplot(gs8[0, 3])

# ── Panel A: heatmaps ──────────────────────────────────────────────────────
cmap_h = plt.cm.RdYlGn
norm_h = Normalize(vmin=0.35, vmax=0.75)

for ax, method in zip(axes_heat, METHOD_LIST):
    sub  = agg_G[agg_G["method"] == method]
    mat  = sub.pivot_table(index="G_threshold",
                            columns="rec_rate",
                            values="mean")
    mat  = mat.reindex(G_VALS)

    im = ax.imshow(mat.values,
                   cmap=cmap_h, norm=norm_h,
                   aspect="auto")

    # annotate cells
    for i in range(len(G_VALS)):
        for j in range(len(REC_VALS)):
            val = mat.values[i, j]
            ax.text(j, i, f"{val:.2f}",
                    ha="center", va="center",
                    fontsize=8.5,
                    color="white" if val > 0.62 else "black")

    ax.set_xticks(range(len(REC_VALS)))
    ax.set_xticklabels(REC_LABELS, fontsize=8, rotation=25)
    ax.set_yticks(range(len(G_VALS)))
    ax.set_yticklabels(G_LABELS, fontsize=9)
    ax.set_title(method, fontsize=11, fontweight="bold",
                 color=METHOD_COLORS[method])
    ax.set_xlabel("Recombination rate", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes_heat[0].set_ylabel("G threshold (generations)", fontsize=9.5)
axes_heat[0].text(-0.18, 1.06, "A", transform=axes_heat[0].transAxes,
                   fontsize=13, fontweight="bold")

cbar8 = fig8.colorbar(cm.ScalarMappable(norm=norm_h, cmap=cmap_h),
                       ax=axes_heat,
                       orientation="vertical",
                       fraction=0.018, pad=0.03, shrink=0.85)
cbar8.set_label("Mean AUPRC", fontsize=9)
cbar8.ax.tick_params(labelsize=8)

# ── Panel B: G-gain curve ─────────────────────────────────────────────────
ax_line.text(-0.18, 1.06, "B", transform=ax_line.transAxes,
              fontsize=13, fontweight="bold")

for method in METHOD_LIST:
    for rv in REC_VALS:
        sub = (agg_G[(agg_G["method"] == method) &
                      (np.isclose(agg_G["rec_rate"], rv))]
               .sort_values("G_threshold"))
        ax_line.plot(sub["G_threshold"],
                     sub["mean"],
                     color=METHOD_COLORS[method],
                     linestyle=REC_LS[rv],
                     linewidth=1.8,
                     marker="o", markersize=4,
                     alpha=0.85)

ax_line.axhline(0.5, color="gray", linestyle="--",
                linewidth=1, alpha=0.5)
ax_line.set_xlabel("G threshold (generations)", fontsize=9.5)
ax_line.set_ylabel("Mean AUPRC", fontsize=9.5)
ax_line.set_title("G-threshold gain curve\nby method and recombination rate",
                   fontsize=9.5)
ax_line.set_xticks(G_VALS)
ax_line.set_xticklabels(G_VALS, fontsize=8.5)
ax_line.set_ylim(0.30, 0.80)

# ── legend for Panel B ────────────────────────────────────────────────────
method_h = [Line2D([0], [0], color=METHOD_COLORS[m],
                    linewidth=2, marker="o", markersize=4,
                    label=m)
             for m in METHOD_LIST]
rec_h = [Line2D([0], [0], color="gray",
                 linestyle=REC_LS[rv],
                 linewidth=2,
                 label=f"r = {REC_LABELS[i]}")
          for i, rv in enumerate(REC_VALS)]
ref_h = Line2D([0], [0], color="gray", linestyle="--",
                linewidth=1, alpha=0.6, label="AUPRC = 0.5")

ax_line.legend(handles=method_h + rec_h + [ref_h],
               fontsize=7.5, frameon=False,
               loc="upper left",
               ncol=1)

fig8.suptitle("Figure 8.  Effect of transmission generation threshold (G) on AUPRC\n"
              "by inference layer and recombination rate",
              fontsize=11, y=1.02)

plt.savefig(f"{OUTPUT}/Figure8_G_threshold_interaction.pdf",
            bbox_inches="tight", dpi=300)
plt.savefig(f"{OUTPUT}/Figure8_G_threshold_interaction.png",
            bbox_inches="tight", dpi=180)
plt.close()
print("Figure 8 saved.")
print("\nDone. Both figures written to /mnt/user-data/outputs/")
