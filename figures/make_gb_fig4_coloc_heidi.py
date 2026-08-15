"""GB Fig 4 — colocalisation versus SMR/HEIDI, with in-sample fine-mapping.

GB 的图注比 FigB 多要了 **"with in-sample fine-mapping"** 一件事。
FigB 只有 coloc PP.H4 对 HEIDI 的散点，没有精细定位那一面。

  a  coloc 判为"不同因果变异"的记录里，HEIDI 有多少没能拒绝同质性。
     253/291 = 86.9%。VPS9D1-AS1 P_HEIDI = 0.649、CDK10 0.086/0.081。
  b  MC1R 区的精细定位：FinnGen 官方**样本内 LD** 给出 3 个高纯度 credible set
     （log₁₀BF 45.6 / 36.6 / 19.5）——该区确实带多个独立因果信号，
     正是"邻近基因的 eQTL 被其中之一 tag 上却并不共享它"的那种构型。
     同一区域我们用 1000G EUR **代理 LD** 跑 susie_rss 得到 10（FinnGen）+ 9（meta）
     = 19 个，**过度拆分**。故代理 LD 的计数**只能单向解读**：
     恰好 1 个 credible set 是"单信号"的保守证据；多个**不是**多信号的证据。

⚠ 样本内计数只有 MC1R 这一个区有（来自 FinnGen 的已发布精细定位），
   其余三区没有，图上不得凭空补。

数据源：16_coloc_meta_results.tsv · 15_SMR_meta_results.tsv · 60a_susie_outcome_signals.tsv
        MC1R 的样本内 3 个 credible set 与其 log₁₀BF 来自 FinnGen 已发布结果，
        登记在 step60a_extract_regions.py 的文件头（正文 Methods 与 [23,25] 同源）。
"""
import math
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_figures as M                      # noqa: E402  复用 CAT / rd / f

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "figures")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": .8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150,
})

# FinnGen 已发布的样本内精细定位（MC1R 区，3 个高纯度 credible set）
INSAMPLE_MC1R_LOG10BF = [45.6, 36.6, 19.5]
C_PROXY, C_INSAMPLE = "#8A8A8A", "#2E7D5B"

su = pd.read_csv(f"{MR}/60a_susie_outcome_signals.tsv", sep="\t")


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    col = {r["exposure"]: r for r in M.rd("16_coloc_meta_results.tsv")}
    smr = {r["gene"] + "|" + r["profile"]: r for r in M.rd("15_SMR_meta_results.tsv")}

    dat = []
    for e, c in col.items():
        s = smr.get(e)
        if not s:
            continue
        h4, ph = M.f(c["PP.H4"]), M.f(s["p_HEIDI"])
        if h4 is None or ph is None or ph <= 0:
            continue
        dat.append((h4, -math.log10(ph), M.CAT.get(c["category"], ("#999", ""))[0],
                    c["SYMBOL"]))

    hi = max(d[1] for d in dat)
    ax.axvspan(0, .5, color="#FCEEED", zorder=0)
    ax.axhline(-math.log10(.05), ls="--", lw=.8, c="#666", zorder=2)
    ax.axvline(.7, ls="--", lw=.8, c="#666", zorder=2)
    ax.scatter([d[0] for d in dat], [d[1] for d in dat], s=20,
               c=[d[2] for d in dat], alpha=.8, edgecolors="none", zorder=3)

    n_lo = [d for d in dat if d[0] < .2]
    n_pass = [d for d in n_lo if d[1] < -math.log10(.05)]
    ax.text(.035, hi * .97,
            "coloc: distinct causal variants (PP.H4 < 0.2)\n"
            "HEIDI: homogeneity not rejected\n"
            f"{len(n_pass)}/{len(n_lo)} = {len(n_pass)/len(n_lo):.1%} pass HEIDI",
            fontsize=8.4, va="top", color="#8B2E28", linespacing=1.45)

    tag = {"VPS9D1-AS1", "CDK10", "SPATA33", "CHMP1A", "CTU2"}
    done = set()
    for h4, y, c, s in dat:
        if h4 < .15 and y < -math.log10(.05) and s in tag and s not in done:
            done.add(s)
            ax.annotate(s, (h4, y), textcoords="offset points", xytext=(5, 2),
                        fontsize=6.6, style="italic", color="#8B2E28")

    ax.text(.985, .74,
            "Evidence tiers that accept MR plus SMR without\n"
            "colocalisation would have reported MC1R linkage\n"
            "spillover as CD4⁺-mediated immune targets.",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.9,
            color="#666", linespacing=1.45,
            bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#DDD", lw=.6))

    ax.set_xlabel("coloc PP.H4  (posterior of a shared causal variant)")
    ax.set_ylabel(r"$-\log_{10}(P_{\rm HEIDI})$   (lower = HEIDI passes)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=5, color=c, label=l)
                       for c, l in M.CAT.values()],
              frameon=False, fontsize=7.2, loc="upper right")
    ax.set_title("a  HEIDI fails to reject the LD-confounded signals coloc identifies",
                 fontsize=9.5, loc="left", pad=6)


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    regions = ["MC1R", "PARP1", "ZFYVE19", "TPI1"]
    y = range(len(regions))[::-1]
    h = .26

    for k, reg in enumerate(regions):
        yy = list(y)[k]
        fg = su[(su.region == reg) & (su.source == "finngen")].iloc[0]
        mt = su[(su.region == reg) & (su.source == "meta")].iloc[0]

        for off, r, lab in ((h + .03, fg, "FinnGen"), (-.03, mt, "meta")):
            ax.barh(yy + off, r.n_cs_highpurity, height=h, color=C_PROXY,
                    alpha=.7, zorder=3)
            ax.text(r.n_cs_highpurity + .22, yy + off, f"{int(r.n_cs_highpurity)}",
                    va="center", fontsize=7.4, color="#555")
            ax.text(-.35, yy + off, lab, ha="right", va="center", fontsize=6.6,
                    color="#999")

        if reg == "MC1R":
            n_in = len(INSAMPLE_MC1R_LOG10BF)
            ax.barh(yy - h - .09, n_in, height=h, color=C_INSAMPLE, alpha=.85,
                    zorder=3)
            ax.text(n_in + .22, yy - h - .09, f"{n_in}", va="center",
                    fontsize=7.4, fontweight="bold", color=C_INSAMPLE)
            ax.text(-.35, yy - h - .09, "in-sample", ha="right", va="center",
                    fontsize=6.6, color=C_INSAMPLE)
            ax.annotate("FinnGen's published in-sample fine-mapping:\n"
                        "three high-purity credible sets, log₁₀BF "
                        + " / ".join(f"{v}" for v in INSAMPLE_MC1R_LOG10BF) +
                        "\n— the region genuinely carries several independent signals",
                        (n_in + .6, yy - h - .09), textcoords="offset points",
                        xytext=(52, -2), fontsize=7, color=C_INSAMPLE,
                        va="center", linespacing=1.4,
                        arrowprops=dict(arrowstyle="-", lw=.6, color=C_INSAMPLE))
        else:
            ax.text(.35, yy - h - .09, "no in-sample release for this region",
                    va="center", fontsize=6.5, color="#BBB", fontstyle="italic")

    ax.set_yticks(list(y))
    ax.set_yticklabels(regions, fontsize=8.6)
    ax.set_xlim(-2.6, 20); ax.set_ylim(-.72, 3.62)
    ax.set_xlabel("High-purity credible sets in the 1 Mb window")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=C_PROXY, alpha=.7),
                       plt.Rectangle((0, 0), 1, 1, fc=C_INSAMPLE, alpha=.85)],
              labels=["susie_rss with a 1000G EUR proxy LD matrix",
                      "FinnGen in-sample LD (published)"],
              frameon=False, fontsize=7.2, loc="lower right")
    ax.text(.985, .30,
            "Proxy LD recovers 19 credible sets at MC1R against 3 in-sample, so it\n"
            "over-splits. These counts are therefore read in one direction only:\n"
            "exactly one credible set is conservative evidence of a single signal;\n"
            "several is not evidence of multiple signals. The eQTL side was not\n"
            "fine-mapped at all — 85–100 donors without in-sample LD.",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.9,
            color="#666", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#DDD", lw=.6))
    ax.set_title("b  In-sample fine-mapping, and what proxy LD does to it",
                 fontsize=9.5, loc="left", pad=6)


def main():
    fig = plt.figure(figsize=(13.6, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.18], wspace=.20,
                          left=.065, right=.985, top=.855, bottom=.105)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    fig.suptitle("Colocalisation and SMR/HEIDI disagree where the region carries "
                 "more than one causal signal", fontsize=11.5, y=.955)
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "Fig4_coloc_heidi" + ext), **kw)
    plt.close(fig)

    mc = su[su.region == "MC1R"]
    print(f"Fig4_coloc_heidi ok  |  MC1R proxy-LD credible sets "
          f"{int(mc[mc.source=='finngen'].n_cs_highpurity.iloc[0])} (FinnGen) + "
          f"{int(mc[mc.source=='meta'].n_cs_highpurity.iloc[0])} (meta) = "
          f"{int(mc.n_cs_highpurity.sum())}  vs  "
          f"{len(INSAMPLE_MC1R_LOG10BF)} in-sample")


if __name__ == "__main__":
    main()
