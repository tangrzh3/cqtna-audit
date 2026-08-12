"""Fig E — spatial compartment attribution of TPI1 (Thrane ST, 4 patients x 2 sections)"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "figures")
ST_DIR = r"D:/Downloads/ST-Melanoma-Datasets_1"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": .8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150,
})
C_TUM, C_LYM, C_STR, C_GLY = "#C4453C", "#3B7DD8", "#8C8C8C", "#E8A33D"

spots = pd.read_csv(f"{MR}/23a_ST_spot_level.tsv", sep="\t")
mt = pd.read_csv(f"{MR}/23d_ST_meta_correlations.tsv", sep="\t")
mk = pd.read_csv(f"{MR}/23h_ST_TPI1_vs_markers.tsv", sep="\t")

comp_cols = ["tumor", "lymphoid", "stroma"]
z = spots.groupby("section")[comp_cols].transform(lambda s: (s - s.mean()) / s.std())
spots["region"] = z.idxmax(axis=1)
spots["margin"] = z.max(axis=1) - z.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
spots.loc[spots.margin < 0.25, "region"] = "mixed"

fig = plt.figure(figsize=(10.6, 7.0))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=.42, wspace=.40,
                      left=.095, right=.985, top=.93, bottom=.10)

# ---------------------------------------------------------------- panel a-c
SEC = "mel2_rep1"
d = spots[spots.section == SEC]
for i, (col, cmap, ttl) in enumerate([
        ("tumor", "Reds", "Melanocyte / tumour score"),
        ("lymphoid", "Blues", "Lymphoid score"),
        ("TPI1", "viridis", "TPI1  (ρ=+.14 tumour, −.08 lymphoid)")]):
    ax = fig.add_subplot(gs[0, i])
    v = d[col]
    sc = ax.scatter(d.x, d.y, c=v, cmap=cmap, s=26, edgecolors="none",
                    vmin=np.nanpercentile(v, 2), vmax=np.nanpercentile(v, 98))
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(ttl, fontsize=9, pad=4)
    cb = fig.colorbar(sc, ax=ax, fraction=.045, pad=.02)
    cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
    if i == 0:
        ax.text(-.02, 1.14, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")
        ax.text(0, -.06, f"section {SEC}  ({len(d)} spots)", transform=ax.transAxes,
                fontsize=7.5, color="#555")

# ---------------------------------------------------------------- panel d
ax = fig.add_subplot(gs[1, 0])
genes = ["TPI1", "HLA-C", "SMC2", "SPSB2", "ZFYVE19", "KIAA0040"]
comps = [("tumor", C_TUM, "tumour"), ("lymphoid", C_LYM, "lymphoid")]
ypos = np.arange(len(genes))[::-1]
for j, (comp, col, lab) in enumerate(comps):
    off = (j - .5) * .30
    sub = mt[mt.compartment == comp].set_index("gene").reindex(genes)
    ax.errorbar(sub.rho_meta, ypos + off,
                xerr=[sub.rho_meta - sub.ci_lo, sub.ci_hi - sub.rho_meta],
                fmt="o", ms=5, lw=1.4, color=col, capsize=2.5, label=lab)
ax.axvline(0, color="k", lw=.7, ls="--", alpha=.6)
ax.set_yticks(ypos); ax.set_yticklabels(genes, fontsize=8.5)
ax.set_xlabel("Spearman ρ with compartment score\n(meta of 8 sections, 2,317 spots)", fontsize=8)
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.text(-.24, 1.10, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")
ax.set_title("TPI1 tracks tumour; HLA-C tracks lymphoid", fontsize=8.5, pad=6)

# ---------------------------------------------------------------- panel e
ax = fig.add_subplot(gs[1, 1])
order = {"glycolysis": 0, "melanocyte": 1, "T cell": 2}
cols = {"glycolysis": C_GLY, "melanocyte": C_TUM, "T cell": C_LYM}
mk = mk.sort_values(["marker_class", "rho_meta"],
                    key=lambda s: s.map(order) if s.name == "marker_class" else s,
                    ascending=[True, True])
y = np.arange(len(mk))
ax.barh(y, mk.rho_meta, color=[cols[c] for c in mk.marker_class], height=.7)
ax.errorbar(mk.rho_meta, y, xerr=[mk.rho_meta - mk.ci_lo, mk.ci_hi - mk.rho_meta],
            fmt="none", ecolor="#333", lw=.8, capsize=1.6)
ax.axvline(0, color="k", lw=.7)
ax.set_yticks(y); ax.set_yticklabels(mk.marker, fontsize=7.2)
ax.set_ylim(-.8, len(mk) - .2)
ax.set_xlabel("Spearman ρ with TPI1", fontsize=8)
ax.legend(handles=[Line2D([], [], color=cols[k], lw=6, label=k)
                   for k in ["glycolysis", "melanocyte", "T cell"]],
          frameon=False, fontsize=7.5, loc="upper left")
ax.text(-.26, 1.10, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")
ax.set_title("TPI1 vs individual marker genes", fontsize=8.5, pad=6)

# ---------------------------------------------------------------- panel f
ax = fig.add_subplot(gs[1, 2])
regions = ["tumor", "lymphoid", "stroma"]
rcol = {"tumor": C_TUM, "lymphoid": C_LYM, "stroma": C_STR}
w, gap = .34, .46
for gi, gene in enumerate(["TPI1", "HLA-C"]):
    v = spots.groupby("section")[gene].transform(lambda s: s - s.mean())
    for ri, reg in enumerate(regions):
        vals = v[spots.region == reg].dropna()
        pos = gi * (len(regions) * w + gap) + ri * w
        bp = ax.boxplot(vals, positions=[pos], widths=w * .82, patch_artist=True,
                        showfliers=False, medianprops=dict(color="k", lw=1.1),
                        whiskerprops=dict(lw=.8), capprops=dict(lw=.8),
                        boxprops=dict(lw=.7))
        bp["boxes"][0].set_facecolor(rcol[reg]); bp["boxes"][0].set_alpha(.8)
ax.axhline(0, color="k", lw=.6, ls=":", alpha=.6)
ax.set_xticks([w * 1, len(regions) * w + gap + w * 1])
ax.set_xticklabels(["TPI1", "HLA-C\n(positive control)"], fontsize=8.5)
ax.set_ylabel("expression, section-centred", fontsize=8)
ax.legend(handles=[Line2D([], [], color=rcol[r], lw=6,
                          label=f"{r} (n={(spots.region == r).sum()})") for r in regions],
          frameon=False, fontsize=7.5, loc="lower right")
ax.text(-.24, 1.10, "d", transform=ax.transAxes, fontsize=13, fontweight="bold")
ax.set_title("Expression by assigned region", fontsize=8.5, pad=16)
ylo, yhi = ax.get_ylim()
ax.set_ylim(ylo, yhi + .55)
for x0, x1, lab in [(0, w, "tumour > lymphoid\nP=2.5e-7"),
                    (len(regions) * w + gap, len(regions) * w + gap + w,
                     "tumour < lymphoid\nP=1.9e-6")]:
    yb = yhi + .05
    ax.plot([x0, x0, x1, x1], [yb, yb + .07, yb + .07, yb], color="k", lw=.8)
    ax.text((x0 + x1) / 2, yb + .10, lab, ha="center", va="bottom", fontsize=6.8)

fig.savefig(os.path.join(OUT, "FigE_spatial_compartment.pdf"))
fig.savefig(os.path.join(OUT, "FigE_spatial_compartment.png"), dpi=300)
print("FigE_spatial_compartment ok")
