"""Fig — FinnGen release natural power experiment (Step 59).

面板在看到观测结果之前定死（预注册做法）：
  a  命中数 vs 病例数，叠加登记的预测区间
  b  已知位点 vs 新位点的恢复率轨迹，叠加各自的登记预测带
  c  观测 vs 预测的逐点比对（对角线 = 完全吻合）

配色沿用 FigA/Fig4：已知位点红、新位点蓝。
预测带来自 58a（登记于任何 R8–R11 数据被查看之前）。
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
C_KNOWN, C_NOVEL, C_HITS = "#C4453C", "#3B7DD8", "#333333"

obs = pd.read_csv(f"{MR}/59a_release_trajectory.tsv", sep="\t")
obs = obs[obs["mode"] == "common"].sort_values("cases").reset_index(drop=True)
pred = pd.read_csv(f"{MR}/58a_finngen_reference_predictions.tsv", sep="\t")
m = obs.merge(pred, on=["release", "cases"], suffixes=("", "_p"))
x = m.cases / 1000


def band(ax, lo, hi, c, label):
    ax.fill_between(x, lo, hi, color=c, alpha=.13, lw=0, zorder=1, label=label)


fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))

# ------------------------------------------------------------------ panel a
ax = axes[0]
band(ax, m.pred_hits_lo, m.pred_hits_hi, C_HITS, "registered prediction (5–95%)")
ax.plot(x, m.pred_hits, ls="--", lw=1.2, color=C_HITS, alpha=.6, zorder=2,
        label="predicted mean")
ax.plot(x, m.n_hits, "-o", ms=5, lw=1.8, color=C_HITS, zorder=4, label="observed")
for xi, yi, r in zip(x, m.n_hits, m.release):
    ax.annotate(r, (xi, yi), textcoords="offset points", xytext=(0, 7),
                ha="center", fontsize=7, color="#555")
ax.set_xlabel("Melanoma cases (thousands)")
ax.set_ylabel("Discoveries at FDR < 0.05")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
ax.set_title("a  Discoveries rise with case number", fontsize=9.5, loc="left", pad=6)

# ------------------------------------------------------------------ panel b
ax = axes[1]
band(ax, m.pred_recovery_known_lo, m.pred_recovery_known_hi, C_KNOWN, None)
band(ax, m.pred_recovery_novel_lo, m.pred_recovery_novel_hi, C_NOVEL, None)
ax.plot(x, m.pred_recovery_known, ls="--", lw=1.1, color=C_KNOWN, alpha=.6)
ax.plot(x, m.pred_recovery_novel, ls="--", lw=1.1, color=C_NOVEL, alpha=.6)
ax.plot(x, m.recovery_known, "-o", ms=5, lw=1.9, color=C_KNOWN, zorder=4,
        label="Known pigmentation / naevus loci")
ax.plot(x, m.recovery_novel, "-o", ms=5, lw=1.9, color=C_NOVEL, zorder=4,
        label="Novel loci")
ax.set_ylim(-.03, 1.05)
ax.set_xlabel("Melanoma cases (thousands)")
ax.set_ylabel("Recovery of the R12 reference list")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
ax.set_title("b  Known loci recover early; novel loci do not\n"
             "(dashed = registered prediction)", fontsize=9.5, loc="left", pad=6)

# ------------------------------------------------------------------ panel c
ax = axes[2]
pairs = [("recovery_known", "pred_recovery_known", C_KNOWN, "known-locus recovery"),
         ("recovery_novel", "pred_recovery_novel", C_NOVEL, "novel-locus recovery"),
         ("recovery_main", "pred_recovery_main", C_HITS, "overall recovery")]
for o, p, c, lab in pairs:
    ax.scatter(m[p], m[o], s=42, color=c, alpha=.85, edgecolors="white",
               linewidths=.6, zorder=3, label=lab)
ax.plot([0, 1], [0, 1], ls=":", lw=1, color="#888", zorder=1)
ax.text(.97, .05, "above the line =\nmore stable than predicted\n(expected: releases are nested)",
        transform=ax.transAxes, ha="right", fontsize=7, color="#666")
ax.set_xlim(-.03, 1.05); ax.set_ylim(-.03, 1.05)
ax.set_xlabel("Predicted (registered before observation)")
ax.set_ylabel("Observed")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
ax.set_title("c  Observed vs registered prediction", fontsize=9.5, loc="left", pad=6)

fig.suptitle("FinnGen sequential releases: candidate-list stability as outcome power grows",
             fontsize=11, y=1.0)
fig.tight_layout(rect=[0, 0, 1, .95])
fig.savefig(os.path.join(OUT, "FigG_release_trajectory.pdf"))
fig.savefig(os.path.join(OUT, "FigG_release_trajectory.png"), dpi=300)
print("FigG_release_trajectory ok")
