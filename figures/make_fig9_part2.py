"""Fig 9 — Part II: instrument availability, the glycolytic state, and patient response.

面板（在 Part II 承重命题改为"工具变量可得性"之后重新设计）：
  a  同一变异在各时点的 beta 与 95% CI + p 值 —— 显示 16h 的 p 最小是因为 SE 小，
     而 naive 的 beta 在 0h 反而最大。这是 Part II 最重要的一张图：它把
     "16h 窗口"从"效应增强"更正为"可检测性"。
  b  锁定 signature 与背景的 lead-SE 比值随时点变化 —— 精度混杂的直接证据
  c  患者侧：锁定 16 / 去 TPI1 15 基因在两个时点的方向一致性
  d  ★ 新增（回应 R3 / B3）：multiome 残差轴的家族组成 —— 该轴确实是一个可定义的
     细胞状态，且与活化强度正交（43c 中全部生物学模块均未超过随机地板）
  e  ★ 新增：该轴的染色质特征（活化不变 peak 内的 motif 富集）
     ⚠ 这与 Step 65 被推翻的"工具变量经 AP-1 motif disruption 起作用"**不是同一件事**：
     那条是关于单个变异是否破坏 motif（经验 p=1.0，已撤回）；这里是该轴的
     染色质总体特征，两者不得混为一谈。

数据：step62 的 beta/se 核查、63b、63c、44b_family_composition.tsv、
      43c_residual_axis_verdict.tsv、48b_motif_axis.tsv
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

fig = plt.figure(figsize=(19.5, 4.6))
gs = fig.add_gridspec(1, 5, width_ratios=[1.25, 1, 1, 1.05, 1.15], wspace=.34,
                      left=.042, right=.99, top=.86, bottom=.16)

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
    # NB: the binomial P values that used to be printed here were WITHDRAWN during
    # revision -- the enzymes are correlated (effective n ~ 11.7 of 16), so a
    # binomial against 1/2 is anticonservative by 2-3 orders of magnitude. Counts
    # are descriptive; inference is the correlation-preserving permutation
    # reported in the text (manuscript 2.5, 3.3). Do not restore the P values.
    for i, (f, k, n) in enumerate(zip(frac, d.higher_in_NR, d.n)):
        ax.text(i + (j - .5) * w, f + .015, "{}/{}".format(k, n),
                ha="center", fontsize=6.4, color="#333")
ax.axhline(.5, ls="--", lw=.9, color="#666", zorder=2)
ax.text(2.45, .515, "chance", ha="right", fontsize=7, color="#555")
ax.set_xticks(np.arange(3)); ax.set_xticklabels(labels, fontsize=7.8)
ax.set_ylim(0, 1.12)
ax.set_ylabel("Fraction of enzymes higher in non-responders")
ax.legend(frameon=False, fontsize=7.5, loc="lower center", ncol=2)
ax.set_title("c  Patient CD4 cells: not driven by TPI1",
             fontsize=9.5, loc="left", pad=6)

# ------------------------------------------------------------------ panel d
# the residual axis is a definable state, and it is orthogonal to activation
fam = pd.read_csv(f"{MR}/44b_family_composition.tsv", sep="	")
ver = pd.read_csv(f"{MR}/43c_residual_axis_verdict.tsv", sep="	")
up = (fam[fam.direction == "up"].sort_values("enrichment", ascending=False)
        .head(5).iloc[::-1])
ax = fig.add_subplot(gs[0, 3])
lbl = [f.replace("_", " ") for f in up.family]
cols = [C_LOCK if "glyco" in f else C_BG for f in up.family]
ax.barh(np.arange(len(up)), up.enrichment.values, color=cols, height=.62)
for i, (e, n) in enumerate(zip(up.enrichment.values, up.n_in_top.values)):
    ax.text(e + .5, i, f"{e:.1f}x  (n={n})", va="center", fontsize=7.4)
ax.axvline(1, color="#999", lw=.8, ls=":")
ax.set_yticks(np.arange(len(up))); ax.set_yticklabels(lbl, fontsize=7.8)
ax.set_xlim(0, max(up.enrichment) * 1.42)
ax.set_xlabel("Enrichment in the top of the residual axis")
# ⚠ 2026-08-15 更正：原注释写的是"all 11 biological modules within the matched
# null (9/11)"，自相矛盾且不实。43c 里 OXPHOS 与 Proliferation 两个模块**超过**
# 零模型，而 Proliferation 是 step43 预先写死的**阳性对照**（它不超过则该检验
# 功效不足、零模型无信息）。要说明"该轴不是活化强度的换个说法"，该引的是
# Activation 模块自己贴在零模型上这一条。正文与 v2 图注已同步改。
nnull = int((~ver.exceeds_null).sum())
act = float(ver.loc[ver.module == "Activation", "smd"].iloc[0])
over = ", ".join(ver.loc[ver.exceeds_null, "module"])
ax.text(.98, .04,
        "Activation sits at the matched null (SMD = {:+.3f})\n"
        "{}/{} modules within it; {} exceed\n"
        "(Proliferation is the pre-specified positive control)"
        .format(act, nnull, len(ver), over),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=6.6,
        color="#444", bbox=dict(fc="white", ec="#DDD", lw=.6, pad=2.5))
ax.set_title("d  A definable state, orthogonal to activation",
             fontsize=9.5, loc="left", pad=6)

# ------------------------------------------------------------------ panel e
# chromatin signature of the axis. NB this is NOT the withdrawn Step 65 claim
# (that the instrument acts by disrupting an AP-1 motif; empirical p = 1.0).
mot = pd.read_csv(f"{MR}/48b_motif_axis.tsv", sep="	")
top = mot[mot.contrast == "glyco_high"].sort_values("odds", ascending=False).head(6).iloc[::-1]
ax = fig.add_subplot(gs[0, 4])
ax.barh(np.arange(len(top)), top.odds.values, color="#8E6BB3", height=.62)
for i, (o, q) in enumerate(zip(top.odds.values, top.FDR.values)):
    ax.text(o + .03, i, f"FDR {q:.0e}", va="center", fontsize=7.0)
ax.axvline(1, color="#999", lw=.8, ls=":")
ax.set_yticks(np.arange(len(top)))
ax.set_yticklabels([m[:14] for m in top.motif], fontsize=7.6)
ax.set_xlim(0, max(top.odds) * 1.30)
ax.set_xlabel("Motif odds ratio, axis-high vs background peaks")
ax.text(.98, .04,
        "axis-level chromatin;\nnot the withdrawn\nmotif-disruption claim",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.0,
        color="#444", bbox=dict(fc="white", ec="#DDD", lw=.6, pad=2.5))
ax.set_title("e  Chromatin signature of the same axis",
             fontsize=9.5, loc="left", pad=6)

fig.suptitle("Instrument availability, the glycolytic CD4 state, and checkpoint-blockade response",
             fontsize=11, y=.985)
fig.savefig(os.path.join(OUT, "Fig9_part2.pdf"))
fig.savefig(os.path.join(OUT, "Fig9_part2.png"), dpi=300)
print("Fig9_part2 ok")
