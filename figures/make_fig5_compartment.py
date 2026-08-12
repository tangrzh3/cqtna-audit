"""Fig 5 — Compartment attribution（按 Step 68/69/70 重做，取代原空转四联图）

面板：
  a  各细胞类型的 TPI1 表达与检出率（GSE115978，全部细胞）
  b  Mal vs CD4 的患者内配对：未匹配 / 深度匹配 / 独立队列 三组并列
  c  预设阳性对照（MLANA、PTPRC）
  d  空间语境（Thrane 跨切片 meta ρ）+ TCGA bulk 生存演示

设计要点：主证据是细胞层面（b），空转退居语境（d）。
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "figures")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": .8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150,
})
C_MAL, C_CD4 = "#C4453C", "#3B7DD8"
C_TUM, C_LYM = "#C4453C", "#3B7DD8"

desc = pd.read_csv(f"{MR}/68a_celltype_descriptive.tsv", sep="\t")
p115 = pd.read_csv(f"{MR}/68b_paired_sample_level.tsv", sep="\t")
pdm  = pd.read_csv(f"{MR}/69a_GSE115978_depth_matched_samples.tsv", sep="\t")
p720 = pd.read_csv(f"{MR}/69b_GSE72056_paired.tsv", sep="\t")
st   = pd.read_csv(f"{MR}/23d_ST_meta_correlations.tsv", sep="\t")

fig = plt.figure(figsize=(13.5, 8.2))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1.15, 1.25, 1],
                      hspace=.46, wspace=.36, left=.105, right=.985,
                      top=.90, bottom=.09)

# ---------------------------------------------------------------- panel a
ax = fig.add_subplot(gs[0, 0])
d = desc[desc.type != "?"].sort_values("TPI1_mean")
cols = ["#C4453C" if t == "Mal" else ("#3B7DD8" if t == "T.CD4" else "#B8BCC2")
        for t in d.type]
ax.barh(np.arange(len(d)), d.TPI1_mean, color=cols, zorder=3)
for i, (v, det, n) in enumerate(zip(d.TPI1_mean, d.TPI1_detect, d.n)):
    ax.text(v + .08, i, f"{det:.0%}", va="center", fontsize=6.8, color="#555")
ax.set_yticks(np.arange(len(d)))
ax.set_yticklabels([f"{t} ({n:,})" for t, n in zip(d.type, d.n)], fontsize=8)
ax.set_xlabel("TPI1, mean log-normalised expression")
ax.set_xlim(0, 7.4)
ax.text(.98, .04, "% = fraction of cells\nwith TPI1 detected", transform=ax.transAxes,
        ha="right", fontsize=7, color="#666")
ax.set_title("a  TPI1 by cell type (GSE115978)", fontsize=9.5, loc="left", pad=6)

# ---------------------------------------------------------------- panel b
ax = fig.add_subplot(gs[0, 1:])
blocks = [("GSE115978\nall cells", p115.TPI1_Mal, p115["TPI1_T.CD4"]),
          ("GSE115978\ndepth-matched", pdm.TPI1_mal, pdm.TPI1_cd4),
          ("GSE72056\nindependent cohort", p720.TPI1_Mal, p720["TPI1_T.CD4"])]
xs = [0, 1, 2]
for k, (lab, mal, cd4) in enumerate(blocks):
    xm, xc = k * 1.0 - .13, k * 1.0 + .13
    for a, b in zip(mal, cd4):
        ax.plot([xm, xc], [a, b], color="#BBB", lw=.7, zorder=1)
    ax.scatter([xm] * len(mal), mal, s=26, color=C_MAL, zorder=3,
               edgecolors="white", linewidths=.5)
    ax.scatter([xc] * len(cd4), cd4, s=26, color=C_CD4, zorder=3,
               edgecolors="white", linewidths=.5)
    n_up = int((np.array(mal) > np.array(cd4)).sum())
    ax.text(k, ax.get_ylim()[1], "", ha="center")
    lbl = f"{n_up}/{len(mal)} patients" if k == 0 else f"{n_up}/{len(mal)}"
    ax.text(k, -0.7, lbl, ha="center", fontsize=8, color="#333", weight="bold")
ax.set_xticks(xs)
ax.set_xticklabels([b[0] for b in blocks], fontsize=8)
ax.set_ylabel("TPI1, mean log-normalised expression")
ax.set_ylim(-1.4, 8.2)
ax.text(-.42, -0.55, "patients with\nMal > CD4:", ha="left", fontsize=7.5, color="#333")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, color=C_MAL, label="Malignant"),
                   Line2D([], [], marker="o", ls="", ms=6, color=C_CD4, label="CD4⁺ T")],
          frameon=False, fontsize=8, loc="upper right", ncol=2)
ax.set_title("b  Paired within patients — the effect survives depth matching and replicates",
             fontsize=9.5, loc="left", pad=6)

# ---------------------------------------------------------------- panel c
ax = fig.add_subplot(gs[1, 0])
pcs = [("MLANA", p115.MLANA_Mal, p115["MLANA_T.CD4"]),
       ("PTPRC", p115.PTPRC_Mal, p115["PTPRC_T.CD4"])]
for k, (lab, mal, cd4) in enumerate(pcs):
    xm, xc = k * 1.0 - .13, k * 1.0 + .13
    for a, b in zip(mal, cd4):
        ax.plot([xm, xc], [a, b], color="#BBB", lw=.7, zorder=1)
    ax.scatter([xm] * len(mal), mal, s=24, color=C_MAL, zorder=3,
               edgecolors="white", linewidths=.5)
    ax.scatter([xc] * len(cd4), cd4, s=24, color=C_CD4, zorder=3,
               edgecolors="white", linewidths=.5)
ax.set_xticks([0, 1]); ax.set_xticklabels(["MLANA\n(expect Mal high)",
                                           "PTPRC\n(expect CD4 high)"], fontsize=8)
ax.set_ylabel("Mean log-normalised expression")
ax.set_title("c  Pre-specified positive controls", fontsize=9.5, loc="left", pad=6)

# ---------------------------------------------------------------- panel d
ax = fig.add_subplot(gs[1, 1])
sub = st[st.gene.isin(["TPI1", "HLA-C"])].copy()
order = ["glycolysis", "tumor", "lymphoid", "myeloid"]
sub = sub[sub.compartment.isin(order)]
w = .36
for k, (gene, c) in enumerate([("TPI1", C_TUM), ("HLA-C", "#2E7D5B")]):
    g = sub[sub.gene == gene].set_index("compartment").reindex(order)
    ax.barh(np.arange(len(order)) + (k - .5) * w, g.rho_meta, height=w,
            color=c, zorder=3, label=gene)
ax.axvline(0, lw=.8, color="#666")
ax.set_yticks(np.arange(len(order)))
ax.set_yticklabels(["Glycolysis\nmodule", "Tumour", "Lymphoid", "Myeloid"], fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("Spearman ρ across 8 sections (meta)")
ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(.02, .30))
ax.set_title("d  Spatial context (Thrane)", fontsize=9.5, loc="left", pad=6)

# ---------------------------------------------------------------- panel e
ax = fig.add_subplot(gs[1, 2])
labels = ["TPI1", "TPI1\n+ tumour", "TPI1\n+ tumour\n+ immune", "Immune\nscore"]
hr = [1.232, 1.214, 1.134, 0.702]
pv = [1.61e-3, 3.57e-3, 6.24e-2, 5.62e-7]
cols = ["#C4453C" if p < .05 else "#B8BCC2" for p in pv]
cols[-1] = "#2E7D5B"
ax.bar(np.arange(4), hr, color=cols, zorder=3)
ax.axhline(1, ls="--", lw=.9, color="#666", zorder=2)
for i, (h, p) in enumerate(zip(hr, pv)):
    ax.text(i, h + .012, f"P={p:.0e}".replace("e-0", "e−"), ha="center", fontsize=6.8)
ax.set_xticks(np.arange(4)); ax.set_xticklabels(labels, fontsize=7.2)
ax.set_ylim(0.6, 1.33)
ax.set_ylabel("Hazard ratio, TCGA-SKCM")
ax.set_title("e  Bulk survival is composition", fontsize=9.5, loc="left", pad=6)

fig.suptitle("A gene nominated with CD4-specific instruments is measured, in tissue, as tumour",
             fontsize=11.5, y=.965)
fig.savefig(os.path.join(OUT, "Fig5_compartment.pdf"))
fig.savefig(os.path.join(OUT, "Fig5_compartment.png"), dpi=300)
print("Fig5_compartment ok")
