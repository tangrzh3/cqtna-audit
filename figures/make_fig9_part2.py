"""Fig 9 — Part II: instrument availability, the glycolytic state, and patient response.

面板（在 Part II 承重命题改为"工具变量可得性"之后重新设计）：
  a  同一变异在各时点的 beta 与 95% CI + p 值 —— 显示 16h 的 p 最小是因为 SE 小，
     而 naive 的 beta 在 0h 反而最大。这是 Part II 最重要的一张图：它把
     "16h 窗口"从"效应增强"更正为"可检测性"。
  b  锁定 signature 与背景的 lead-SE 比值随时点变化 —— 精度混杂的直接证据
  c  患者侧：锁定 16 / 去 TPI1 15 基因在两个时点的方向一致性

数据：step62 的 beta/se 核查、63b_precision_by_profile.tsv、63c_patient_locked.tsv
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
C_N, C_M = "#C4453C", "#3B7DD8"          # naive / memory
C_LOCK, C_BG = "#2E7D5B", "#7A8794"

# ------------------------------------------------------------------ 数据
# 同一变异 12:6867132 在各 profile 的 beta/se/p（Step 62 核查所得）
TP = ["0h", "16h", "40h", "5d"]
BETA = {"Naive":  [-0.427, -0.177, -0.085, +0.017],
        "Memory": [-0.151, -0.193, -0.049, +0.033]}
SE = {"Naive":  [0.107, 0.021, 0.040, 0.079],
      "Memory": [0.086, 0.029, 0.030, 0.077]}
PV = {"Naive":  [1.51e-4, 1.74e-12, 3.50e-2, 8.27e-1],
      "Memory": [8.48e-2, 4.22e-9, 1.11e-1, 6.75e-1]}

prec = pd.read_csv(f"{MR}/63b_precision_by_profile.tsv", sep="\t")
pat = pd.read_csv(f"{MR}/63c_patient_locked.tsv", sep="\t")

fig = plt.figure(figsize=(13, 4.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=.30,
                      left=.06, right=.985, top=.86, bottom=.14)

# ------------------------------------------------------------------ panel a
ax = fig.add_subplot(gs[0, 0])
x = np.arange(4)
for lin, c, off in [("Naive", C_N, -.09), ("Memory", C_M, +.09)]:
    b = np.array(BETA[lin]); s = np.array(SE[lin])
    ax.errorbar(x + off, b, yerr=1.96 * s, fmt="o", ms=5.5, lw=1.4,
                capsize=3.5, color=c, label=lin, zorder=4)
ax.axhline(0, ls=":", lw=.8, color="#888", zorder=1)
ax.axvspan(0.5, 1.5, color="#F0F0F0", zorder=0)
ax.text(1, .16, "the only timepoint\nwith an instrument", ha="center",
        fontsize=7.5, color="#666")
# p 值标注
for lin, c, off, va in [("Naive", C_N, -.09, "top"), ("Memory", C_M, +.09, "bottom")]:
    for i, p in enumerate(PV[lin]):
        txt = f"{p:.0e}".replace("e-0", "e−").replace("e-", "e−")
        ax.annotate(txt, (x[i] + off, BETA[lin][i] - 1.96 * SE[lin][i]),
                    textcoords="offset points", xytext=(0, -9),
                    ha="center", fontsize=6.2, color=c)
ax.set_xticks(x); ax.set_xticklabels(TP)
ax.set_xlabel("Time after anti-CD3/CD28 stimulation")
ax.set_ylabel("cis-eQTL effect on TPI1  (β ± 95% CI)")
ax.set_ylim(-0.68, 0.26)
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_title("a  The 16 h signal is precision, not amplitude",
             fontsize=9.5, loc="left", pad=6)

# ------------------------------------------------------------------ panel b
ax = fig.add_subplot(gs[0, 1])
for lin, c in [("Naive", C_N), ("Memory", C_M)]:
    d = prec[prec.lineage == lin].set_index("timepoint").loc[TP]
    ax.plot(np.arange(4), d.ratio.values, "-o", ms=5, lw=1.6, color=c, label=lin)
ax.axhline(1, ls=":", lw=.8, color="#888")
ax.text(3.02, 1.005, "no advantage\nover background", ha="right", va="bottom",
        fontsize=7, color="#777")
ax.set_xticks(np.arange(4)); ax.set_xticklabels(TP)
ax.set_ylim(0.5, 1.05)
ax.set_xlabel("Time after stimulation")
ax.set_ylabel("Lead-variant SE, locked signature / background")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_title("b  Measurability, not regulation, tracks activation",
             fontsize=9.5, loc="left", pad=6)

# ------------------------------------------------------------------ panel c
ax = fig.add_subplot(gs[0, 2])
sets = ["all22", "locked16", "locked15_noTPI1"]
labels = ["all 22", "locked 16", "locked 15\n(TPI1 removed)"]
w = 0.36
for j, (tpt, c) in enumerate([("Pre", "#7A5AA8"), ("Post", "#E8A33D")]):
    d = pat[pat.timepoint == tpt].set_index("set").loc[sets]
    frac = d.higher_in_NR / d.n
    ax.bar(np.arange(3) + (j - .5) * w, frac, width=w, color=c,
           label=f"{tpt}-treatment", zorder=3)
    for i, (f, k, n, p) in enumerate(zip(frac, d.higher_in_NR, d.n, d.binom_p)):
        ax.text(i + (j - .5) * w, f + .015, f"{k}/{n}\nP={p:.0e}".replace("e-0", "e−"),
                ha="center", fontsize=6.4, color="#333")
ax.axhline(.5, ls="--", lw=.9, color="#666", zorder=2)
ax.text(2.45, .515, "chance", ha="right", fontsize=7, color="#555")
ax.set_xticks(np.arange(3)); ax.set_xticklabels(labels, fontsize=7.8)
ax.set_ylim(0, 1.12)
ax.set_ylabel("Fraction of enzymes higher in non-responders")
ax.legend(frameon=False, fontsize=7.5, loc="lower center", ncol=2)
ax.set_title("c  Patient CD4 cells: not driven by TPI1",
             fontsize=9.5, loc="left", pad=6)

fig.suptitle("Instrument availability, the glycolytic CD4 state, and checkpoint-blockade response",
             fontsize=11, y=.985)
fig.savefig(os.path.join(OUT, "Fig9_part2.pdf"))
fig.savefig(os.path.join(OUT, "Fig9_part2.png"), dpi=300)
print("Fig9_part2 ok")
