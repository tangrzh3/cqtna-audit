"""GB Fig 3 — power and list stability by locus class, and the same stratification
against full-power |z|.

GB 的图注比全文源的 FigF 多要了一件事：**"the same stratification against
full-power |z|"**。FigF 只有 a/b/c 三个面板（功效曲线、五癌种、按位点类别的恢复率），
条件化在效应量上的那一步没有图。本脚本复用 FigF 的三个面板函数，补上 d、e：

  d  两类候选的全功效 |z| 分布：新位点候选**按构造**贴着检出阈值
     （全部落在 3.73–4.55），已知位点候选可达 15.99。
  e  |z| 匹配后残余的类别差：原始差 47.6（按位点）/ 63.0（按基因）个百分点，
     匹配后只剩 +5.2 [−1.6, +12.0] / +6.4 [+0.9, +12.0]。
     按位点那一格跨过 0，**故不主张残余可归因于类别本身**。

⚠ `FigF_power_stability` 保持不变（它是全文源八图方案里的 Fig 4）。本脚本另存文件。
⚠ 正文的"仅用 |z| 的模型重现 68–80% 的差距"**没有落盘的表**，
   故图上不画这个数，只画 101a/101b 里有的量。

数据源：53a · 53b · 54a · 54b · 55a（经 FigF 的面板函数）· 101a_zdist.tsv · 101b_matched.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_fig_power_stability as F          # noqa: E402  面板 a/b/c 直接复用

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "figures")
C_KNOWN, C_NOVEL = F.C_KNOWN, F.C_NOVEL

zd = pd.read_csv(f"{MR}/101a_zdist.tsv", sep="\t")
mt = pd.read_csv(f"{MR}/101b_matched.tsv", sep="\t")

UNIT_LAB = {"locus": "by independent locus", "gene": "by gene"}


# ------------------------------------------------------------------ panel d
def panel_d(ax):
    rows = [("locus", "known"), ("locus", "novel"),
            ("gene", "known"), ("gene", "novel")]
    ypos = [3.4, 2.6, 1.2, .4]

    for (unit, cls), y in zip(rows, ypos):
        r = zd[(zd.unit == unit) & (zd.cls == cls)].iloc[0]
        c = C_KNOWN if cls == "known" else C_NOVEL
        ax.plot([r.lo, r.hi], [y, y], lw=1.1, color=c, zorder=3)
        for e in (r.lo, r.hi):
            ax.plot([e, e], [y - .09, y + .09], lw=1.1, color=c, zorder=3)
        ax.add_patch(plt.Rectangle((r.q1, y - .17), r.q3 - r.q1, .34,
                                   facecolor=c, alpha=.28, edgecolor=c,
                                   lw=1.1, zorder=4))
        ax.plot([r["median"], r["median"]], [y - .19, y + .19], lw=2.1,
                color=c, zorder=5)
        ax.text(r.hi + .45, y, f"n = {int(r.n)}", va="center", fontsize=7,
                color=c)

    ax.axvspan(3.7339154127190537, 4.546214923018536, color=C_NOVEL, alpha=.07,
               zorder=0)
    ax.annotate("every novel-locus candidate lies in this band (3.73–4.55)",
                (4.55, 1.9), textcoords="offset points", xytext=(30, 0),
                ha="left", va="center", fontsize=7, color=C_NOVEL,
                arrowprops=dict(arrowstyle="-", lw=.6, color=C_NOVEL))

    ax.set_yticks(ypos)
    ax.set_yticklabels(["known", "novel", "known", "novel"], fontsize=8)
    for y, lab in ((3.0, "by independent\nlocus"), (.8, "by gene")):
        ax.text(-.155, y, lab, transform=ax.get_yaxis_transform(),
                ha="center", va="center", fontsize=7.4, color="#666",
                linespacing=1.35)
    ax.set_xlim(2.6, 18.6); ax.set_ylim(-.15, 4.35)
    ax.set_xlabel("Full-power |z| of the candidate")
    ax.tick_params(axis="y", length=0)
    ax.set_title("d  Novel-locus candidates sit against the detection threshold "
                 "by construction", fontsize=9.5, loc="left", pad=6)


# ------------------------------------------------------------------ panel e
def panel_e(ax):
    for k, unit in enumerate(("locus", "gene")):
        r = mt[mt.unit == unit].iloc[0]
        y = 1 - k
        raw, res = 100 * r.gap_observed, 100 * r.delta_matched
        lo, hi = 100 * r.ci_lo, 100 * r.ci_hi

        ax.barh(y + .17, raw, height=.30, color=C_KNOWN, alpha=.80, zorder=3)
        ax.text(raw + 1.4, y + .17, f"{raw:.1f} pp", va="center", fontsize=7.6,
                color=C_KNOWN)
        ax.barh(y - .17, res, height=.30, color="#8A8A8A", alpha=.70, zorder=3)
        ax.errorbar(res, y - .17, xerr=[[res - lo], [hi - res]], fmt="none",
                    ecolor="#555", elinewidth=1.1, capsize=3, zorder=5)
        ax.text(hi + 1.4, y - .17, f"{res:+.1f} pp  [{lo:+.1f}, {hi:+.1f}]",
                va="center", fontsize=7.6, color="#555")
        ax.text(-1.6, y, UNIT_LAB[unit], ha="right", va="center", fontsize=8.2)

    ax.axvline(0, lw=.9, color="#333", zorder=2)
    ax.set_yticks([])
    ax.set_xlim(-14, 88); ax.set_ylim(-1.42, 2.05)
    ax.set_xlabel("Known-minus-novel recovery difference (percentage points)")
    ax.spines["left"].set_visible(False)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=C_KNOWN, alpha=.80),
                       plt.Rectangle((0, 0), 1, 1, fc="#8A8A8A", alpha=.70)],
              labels=["raw difference", "residual after matching on |z| (95% CI)"],
              frameon=False, fontsize=7.2, loc="upper right")
    ax.text(37, -.92,
            "By independent locus the residual interval crosses zero, and only one known\n"
            "locus falls inside the novel |z| range — so the remainder is not separable\n"
            "from the failure of matching, and we do not attribute it to the category.\n"
            "What is retained is threshold proximity (Supplementary S28).",
            ha="center", va="center", fontsize=6.9, color="#666", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#DDD", lw=.6))
    ax.set_title("e  Conditioning on effect size absorbs most of the class gap",
                 fontsize=9.5, loc="left", pad=6)


def main():
    fig = plt.figure(figsize=(13.0, 12.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.12, 1.02], hspace=.44,
                          wspace=.34, left=.095, right=.965, top=.925, bottom=.055)
    F.panel_a(fig.add_subplot(gs[0, 0]))
    F.panel_b(fig.add_subplot(gs[0, 1]))
    F.panel_c(fig.add_subplot(gs[1, :]))
    panel_d(fig.add_subplot(gs[2, 0]))
    panel_e(fig.add_subplot(gs[2, 1]))
    fig.suptitle("Outcome GWAS power sets the reproducibility of the candidate list — "
                 "and the class gap is mostly threshold proximity",
                 fontsize=11.5, y=.962)
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "Fig3_power_stability" + ext), **kw)
    plt.close(fig)

    for _, r in mt.iterrows():
        print(f"Fig3  {r.unit:<6} raw gap {100*r.gap_observed:.1f} pp  ->  "
              f"matched residual {100*r.delta_matched:+.1f} pp "
              f"[{100*r.ci_lo:+.1f}, {100*r.ci_hi:+.1f}]  "
              f"criterion_met={r.criterion_met}")


if __name__ == "__main__":
    main()
