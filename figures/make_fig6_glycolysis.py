"""Fig 6 -- TPI1 / CD4 glycolysis: genetic window, pathway shift, ICB response

Panel a  cis-eQTL strength across CD4 activation (TPI1 pulse vs SPSB2 constitutive)
Panel b  glycolytic enzymes, responder vs non-responder (corrected, Step 26)
Panel c  TPI1 and the glycolysis module by response, both timepoints
All numbers from 27a (eQTL) and 32b/32c/32e (corrected single-cell layer).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "figures")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": .8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150,
})
C_R, C_NR = "#3B7DD8", "#C4453C"          # responder / non-responder
C_TPI1, C_SPSB2 = "#C4453C", "#7A8794"

tc = pd.read_csv(f"{MR}/27a_TPI1_eqtl_timecourse.tsv", sep="\t")
gg = pd.read_csv(f"{MR}/32e_glyco_gene_by_gene_corrected.tsv", sep="\t")
sm = pd.read_csv(f"{MR}/32b_sample_level_corrected.tsv", sep="\t")
fam = pd.read_csv(f"{MR}/32c_response_family_corrected.tsv", sep="\t")

fig = plt.figure(figsize=(11.2, 7.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1, 1.1],
                      hspace=.45, wspace=.30, left=.085, right=.98,
                      top=.92, bottom=.09)

# ----------------------------------------------------------------- panel a
ax = fig.add_subplot(gs[0, 0])
stages = ["0h resting", "16h pre-division", "40h post-1st-division", "5d effector"]
xs = np.arange(4)
for gene, col, mk in [("TPI1", C_TPI1, "o"), ("SPSB2", C_SPSB2, "s")]:
    for lin, ls in [("Naive", "-"), ("Memory", "--")]:
        d = tc[(tc.gene == gene) & (tc.profile.str.contains(lin))]
        d = d.set_index("stage").reindex(stages)
        ax.plot(xs, -np.log10(d.top_pval.values), ls, marker=mk, ms=5,
                color=col, lw=1.6, alpha=.95 if lin == "Naive" else .55,
                label=f"{gene} {lin}")
ax.axhline(-np.log10(5e-8), color="k", lw=.8, ls=":", alpha=.7)
ax.text(3.02, -np.log10(5e-8) + .4, "p = 5e-8", fontsize=7, ha="right", color="#444")
ax.set_xticks(xs)
ax.set_xticklabels(["0h\nresting", "16h\npre-division", "40h\npost-division", "5d\neffector"],
                   fontsize=7.5)
ax.set_ylabel("strongest cis-eQTL,  $-\\log_{10}P$", fontsize=8.5)
ax.set_title("TPI1's eQTL is an activation-window pulse", fontsize=9, pad=6)
ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper right")
ax.text(-.20, 1.10, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")

# ----------------------------------------------------------------- panel b
ax = fig.add_subplot(gs[0, 1])
d = gg[gg.timepoint == "Post"].copy().sort_values("diff")
d = d[d.detection >= 5]
y = np.arange(len(d))
# bake per-bar transparency into RGBA: solid = FDR<0.05, faded = not significant
_rgba = matplotlib.colors.to_rgba
cols = [_rgba("#C4453C" if v < 0 else "#3B7DD8", 1.0 if f < .05 else .35)
        for v, f in zip(d["diff"], d.FDR)]
ax.barh(y, d["diff"], color=cols, height=.72)
ax.axvline(0, color="k", lw=.8)
ax.set_yticks(y)
ax.set_yticklabels([f"{g}" + (" *" if f < .05 else "")
                    for g, f in zip(d.gene, d.FDR)], fontsize=7.2)
for i, (g, f) in enumerate(zip(d.gene, d.FDR)):
    if g == "TPI1":
        ax.get_yticklabels()[i].set_fontweight("bold")
ax.set_xlabel("median difference, responder − non-responder\n(log$_2$(TPM+1), post-treatment)",
              fontsize=8)
ax.set_title("Glycolytic enzymes are coordinately higher in non-responders",
             fontsize=9, pad=6)
nb = gg[gg.timepoint == "Post"]
# the binomial P that used to be printed here (6.1e-5) was withdrawn during
# revision: the enzymes are correlated, so a binomial against 1/2 is
# anticonservative by 2-3 orders of magnitude. The count is descriptive; the
# correlation-preserving permutation value is the one reported (manuscript 2.7, 3.3).
ax.text(.03, .95, f"20/22 genes higher in non-responders\n"
                  f"(descriptive count; permutation P = 0.088)   "
                  f"(* FDR < 0.05, n = {(nb.FDR < .05).sum()})",
        transform=ax.transAxes, fontsize=7.3, va="top",
        bbox=dict(fc="white", ec="#ccc", lw=.6, pad=3))
ax.text(-.24, 1.10, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")

# ----------------------------------------------------------------- panel c
for j, (var, lab) in enumerate([("TPI1", "TPI1"),
                                ("GlycoScore", "Glycolysis module\n(TPI1 excluded)")]):
    ax = fig.add_subplot(gs[1, j])
    xpos = {"Pre": 0, "Post": 1}
    for tp in ["Pre", "Post"]:
        for resp, col, off in [("Responder", C_R, -.17), ("Non-responder", C_NR, .17)]:
            v = sm[(sm.timepoint == tp) & (sm.response == resp)][var].values
            x = xpos[tp] + off
            ax.scatter(np.random.default_rng(0).normal(x, .035, len(v)), v,
                       s=17, color=col, alpha=.75, edgecolors="none", zorder=3)
            ax.plot([x - .09, x + .09], [np.median(v)] * 2, color=col, lw=2.2, zorder=4)
    for tp in ["Pre", "Post"]:
        f = fam[(fam.timepoint == tp) & (fam.variable == var)].iloc[0]
        ax.text(xpos[tp], ax.get_ylim()[1], f"FDR = {f.FDR:.3g}", ha="center",
                va="bottom", fontsize=7.6,
                fontweight="bold" if f.FDR < .05 else "normal")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Pre\n(9 R / 10 NR)", f"Post\n(8 R / 20 NR)"], fontsize=8)
    ax.set_xlim(-.5, 1.5)
    ax.set_ylabel(f"{lab}\nmean across CD4$^+$ T cells", fontsize=8)
    ax.margins(y=.16)
    if j == 0:
        ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C_R, label="Responder"),
                           Line2D([], [], marker="o", ls="", color=C_NR, label="Non-responder")],
                  frameon=False, fontsize=7.5, loc="lower left")
        ax.text(-.20, 1.12, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")

fig.savefig(os.path.join(OUT, "Fig6_glycolysis_TPI1.pdf"))
fig.savefig(os.path.join(OUT, "Fig6_glycolysis_TPI1.png"), dpi=300)
print("Fig6_glycolysis_TPI1 ok")
