"""GB Fig 6 — the activation window in which TPI1 is instrumentable:
effect size against precision across the eight profiles.

这张图承担的是本文对自己早先读法的**更正**：
16 h 不是"效应被放大"的时点，而是"效应恰好可测"的时点。
效应量在静息时最大（Naive 0h，β = −0.427），只是那里噪声大得多；
唯一越过全基因组显著线的是两个 16 h profile。

  a  同一变异 12:6867132 在八个 profile 上的 β ± 95% CI。
  b  |β| 对 SE 的平面：全基因组显著线在这张平面上是一条过原点的直线
     （|β| = 5.45 × SE）。0h 的点在**纵轴上最高**却落在线下——
     "最大的效应"与"唯一可用的工具变量"不是同一个时点，一眼可见。

⚠ 这张图与已发表的 genotype × pseudotime 交互检验[8]**不矛盾也不重复**：
   那个检验问的是"效应是否随伪时间变化"（三个 lead 变异全部为阴性），
   本图问的是"效应在哪个时点测得准到能当工具变量"。两者是不同的问题，
   图注与正文都必须保持这个区分。

数据源：120a_tpi1_instrument_window.tsv（step120 从八个 profile 的 parquet 重算，
        逐位复现了此前写死在 make_fig9_part2.py 里的数）。
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
TP = ["0h", "16h", "40h", "5d"]
Z_GWS = 5.4513104                        # |z| at P = 5e-8, two-sided

d = pd.read_csv(f"{MR}/120a_tpi1_instrument_window.tsv", sep="\t")
d["timepoint"] = pd.Categorical(d.timepoint, categories=TP, ordered=True)
d = d.sort_values(["lineage", "timepoint"])


def sub(lin):
    return d[d.lineage == lin].set_index("timepoint").loc[TP]


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    x = np.arange(4)
    ax.axvspan(.5, 1.5, color="#F0F0F0", zorder=0)
    ax.axhline(0, ls=":", lw=.8, color="#888", zorder=1)

    for lin, c, off in (("Naive", C_N, -.10), ("Memory", C_M, +.10)):
        s = sub(lin)
        ax.errorbar(x + off, s.beta, yerr=1.96 * s.se, fmt="o", ms=5.5, lw=1.4,
                    capsize=3.5, color=c, label=lin, zorder=4)
        for i, (b, se, p, inst) in enumerate(zip(s.beta, s.se, s.pval,
                                                 s.has_instrument)):
            ax.annotate(f"{p:.0e}".replace("e-0", "e−").replace("e-", "e−"),
                        (x[i] + off, b - 1.96 * se), textcoords="offset points",
                        xytext=(0, -9), ha="center", fontsize=6.3,
                        color=c, fontweight="bold" if inst else "normal")

    ax.text(1, .175, "the only timepoint with\na usable instrument",
            ha="center", fontsize=7.4, color="#666", linespacing=1.35)
    s = sub("Naive")
    ax.annotate("largest effect — and the noisiest",
                (0 - .10, s.beta.iloc[0]), textcoords="offset points",
                xytext=(34, -6), fontsize=7, color=C_N,
                arrowprops=dict(arrowstyle="-", lw=.6, color="#CCC"))

    ax.set_xticks(x); ax.set_xticklabels(TP)
    ax.set_xlim(-.55, 3.55); ax.set_ylim(-.72, .30)
    ax.set_xlabel("Time after anti-CD3/CD28 stimulation")
    ax.set_ylabel("cis-eQTL effect on TPI1  (β ± 95% CI)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("a  The same variant, all eight profiles",
                 fontsize=9.5, loc="left", pad=6)


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    xmax = .125
    ymax = .62
    xs = np.linspace(0, xmax, 50)
    ax.fill_between(xs, Z_GWS * xs, ymax, color="#C4453C", alpha=.055, lw=0, zorder=0)
    ax.plot(xs, Z_GWS * xs, ls="--", lw=1.1, color="#8B2E28", zorder=2)
    ax.set_xlim(0, xmax); ax.set_ylim(0, ymax)

    # 斜线标签的角度必须按显示坐标算，按数据坐标算会随长宽比跑偏
    x0, x1 = .088, .098
    (px0, py0), (px1, py1) = ax.transData.transform(
        [(x0, Z_GWS * x0), (x1, Z_GWS * x1)])
    ax.text((x0 + x1) / 2, Z_GWS * (x0 + x1) / 2 + .012, "$P$ = 5×10⁻⁸",
            ha="center", va="bottom", fontsize=7.4, color="#8B2E28",
            rotation=np.degrees(np.arctan2(py1 - py0, px1 - px0)),
            rotation_mode="anchor")
    ax.text(.0055, .128, "instrumentable", fontsize=7.6, color="#8B2E28",
            va="center")

    for lin, c, mk in (("Naive", C_N, "o"), ("Memory", C_M, "s")):
        s = sub(lin)
        ax.scatter(s.se, s.beta.abs(), s=58, marker=mk, color=c,
                   edgecolors="white", linewidths=.8, zorder=4, label=lin)
        for tp, se, ab, inst in zip(TP, s.se, s.beta.abs(), s.has_instrument):
            ax.annotate(tp, (se, ab), textcoords="offset points",
                        xytext=(-7, 7) if lin == "Naive" else (8, -10),
                        ha="right" if lin == "Naive" else "left",
                        fontsize=7.2, color=c,
                        fontweight="bold" if inst else "normal")

    ax.set_xlabel("Standard error of the eQTL effect  (imprecision →)")
    ax.set_ylabel("|β|  of the eQTL effect on TPI1")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.text(.035, .965,
            "The point highest on the vertical axis — naive 0 h — is the\n"
            "largest effect in the series, and it is not instrumentable.\n"
            "Both 16 h points cross the line on precision alone.",
            transform=ax.transAxes, ha="left", va="top", fontsize=7,
            color="#444", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#DDD", lw=.6))
    ax.set_title("b  Effect size against precision", fontsize=9.5, loc="left", pad=6)


def main():
    fig = plt.figure(figsize=(12.4, 6.0))
    gs = fig.add_gridspec(1, 2, wspace=.24, left=.068, right=.985,
                          top=.875, bottom=.245)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    fig.suptitle("Time-specific significance is not a dynamic genetic effect",
                 fontsize=11.5, y=.962)
    fig.text(.068, .015,
             "What the window is not.  In the published genotype × pseudotime "
             "interaction test [8], none of the three TPI1 lead variants shows an "
             "interaction — memory cells, linear $P$ = 0.974 and\n"
             "quadratic 0.992; naive cells, 0.707 and 0.891 for the variant 1,063 bp "
             "from our instrument, 0.378 and 0.605 for the other — while TPI1 "
             "expression is among the most strongly\n"
             "pseudotime-dependent genes in that dataset (Moran's I = 0.664).  "
             "What moves across the activation window is measurability, not regulation.",
             ha="left", va="bottom", fontsize=7.2, color="#666", linespacing=1.6)
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "Fig6_tpi1_window" + ext), **kw)
    plt.close(fig)

    big = d.loc[d.beta.abs().idxmax()]
    tight = d.loc[d.se.idxmin()]
    inst = d[d.has_instrument]
    print(f"Fig6_tpi1_window ok  |  instrumentable in "
          f"{'; '.join(inst.lineage + ' ' + inst.timepoint.astype(str))}  |  "
          f"largest |beta| {big.lineage} {big.timepoint} ({big.beta:+.3f}, "
          f"se {big.se:.3f})  |  smallest se {tight.lineage} {tight.timepoint}")


if __name__ == "__main__":
    main()
