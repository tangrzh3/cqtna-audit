"""GB Fig 1 — locus attribution under both outcomes, by bounded locus.

图注要求三件事，缺一不可：
  · **两个**结局（FinnGen R12 轮 与 meta 轮），不是只画 meta；
  · **按有界位点**计数，不是按基因记录——正文报的 4.96× 是位点级的数
  · 位点划分 = 非递归固定锚定窗口 1000 kb（S36 冻结），不再是单连锁
    （记录级为 5.26×，Step 30 的自设纪律是正文一律用位点级，见 FIGURES_plan.md §2）；
  · 端到端对照：MC1R 区给出全研究最强关联 P = 4×10⁻³⁷，PARP1 方向复现。

面板：
  a  镜像 Manhattan，**每个有界位点一个点**（取该位点内最强记录）。
     上 = FinnGen R12 轮，下 = meta 轮。颜色 = 位点类别。
  b  位点级归属：两个结局各自的"显著位点中已知的比例"对上各自的背景比例。

⚠ `FigA_locus_attribution`（`make_figures.py::fig_A`）是**按记录**的单结局 Manhattan，
   是全文源八图方案里的 Fig 1，与本图并存，**不要互相覆盖**。

数据源：06_locus_annotation.tsv（FinnGen 轮）· 13_meta_locus_annotation.tsv（meta 轮）
        位点划分用与 step85/step119 同一条 1 Mb 单连锁规则，脚本内现算，不外读。
"""
import os
from math import lgamma, exp

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
LOCUS_KB = 1000
# 按位点取最强记录后，两个结局各自**只有一个**位点（MC1R 区）超过 -log10 P = 8，
# 其余 550 余个全在 0–6。沿用 fig_A 的 CAP=30 会把整张图压到中线上，
# 故此处截到 8，超出者画成三角并单独标注。
CAP = 8.0

CAT = {
    "潜在新位点":      ("#3B7DD8", "Novel locus"),
    "痣数目通路":      ("#E8A33D", "Naevus-count locus"),
    "色素/发色通路":   ("#C4453C", "Pigmentation locus"),
    "黑色素瘤已知位点": ("#7A5AA8", "Known melanoma locus"),
}
NOVEL = "潜在新位点"
CHRLEN = {str(i): l for i, l in enumerate(
    [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
     145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
     101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
     50818468], 1)}


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    hi = min(r1, c1)
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, hi + 1)), 1.0)


def assign_loci(df):
    """1 Mb single-linkage clustering of instrument positions within a chromosome."""
    # Non-recursive fixed-anchor partition, the frozen main analysis
    # (manuscript/PREREG_locus_partition.md, S36). The first unassigned variant
    # on a chromosome becomes an anchor and claims everything within LOCUS_KB of
    # it; the next unassigned variant becomes the next anchor. The span is
    # bounded by the window whatever the density, unlike the single-linkage rule
    # this figure used to draw, which chains a dense resource into blocks of tens
    # of megabases. Identical to cqtna:::cq_assign_loci, ">" boundary included.
    out = {}
    for ch, sub in df.groupby("chr"):
        sub = sub.sort_values("pos")
        lid, anchor = 0, None
        for _, r in sub.iterrows():
            if anchor is None or r.pos - anchor > LOCUS_KB * 1000:
                lid += 1
                anchor = r.pos
            out[(ch, r.pos)] = f"{ch}_{lid}"
    return out


def load(fn):
    d = pd.read_csv(f"{MR}/{fn}", sep="\t")
    d[["chr", "pos"]] = d.SNP.str.split(":", expand=True)
    d["pos"] = d.pos.astype(int)
    d = d[d.chr.isin(CHRLEN)].copy()
    d["FDR"] = pd.to_numeric(d.FDR, errors="coerce")
    d["pval"] = pd.to_numeric(d.pval, errors="coerce")
    snp = d.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp)
    d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
    d["known"] = d.category != NOVEL
    return d


FG = load("06_locus_annotation.tsv")
MT = load("13_meta_locus_annotation.tsv")

OFF, CUM = {}, 0
for c in map(str, range(1, 23)):
    OFF[c] = CUM
    CUM += CHRLEN[c]


def per_locus(d):
    """每个有界位点取最强记录，Manhattan 按位点画而不是按记录。"""
    d = d.copy()
    d["y"] = -np.log10(d.pval.clip(lower=1e-320))
    idx = d.groupby("locus").y.idxmax()
    top = d.loc[idx].copy()
    top["x"] = [OFF[c] + p for c, p in zip(top.chr, top.pos)]
    return top


def counts(d):
    bg = d.groupby("locus").agg(known=("known", "any"))
    sg = d[d.FDR < 0.05].groupby("locus").agg(known=("known", "any"))
    BG_T, BG_K, S_T, S_K = len(bg), int(bg.known.sum()), len(sg), int(sg.known.sum())
    fold = (S_K / S_T) / (BG_K / BG_T)
    p = fisher_greater(S_K, S_T - S_K, BG_K - S_K, (BG_T - BG_K) - (S_T - S_K))
    return dict(bg_loci=BG_T, bg_known=BG_K, sig_loci=S_T, sig_known=S_K,
                bg_pct=100 * BG_K / BG_T, sig_pct=100 * S_K / S_T,
                fold=fold, p=p)


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    fg, mt = per_locus(FG), per_locus(MT)
    fg_thr = -np.log10(FG.loc[FG.FDR < 0.05, "pval"].max())
    mt_thr = -np.log10(MT.loc[MT.FDR < 0.05, "pval"].max())

    for c in map(str, range(1, 23)):
        if int(c) % 2 == 0:
            ax.axvspan(OFF[c], OFF[c] + CHRLEN[c], color="#F4F4F4", zorder=0)
    ax.axhline(0, lw=.9, color="#333", zorder=5)

    for d, sign in ((fg, 1), (mt, -1)):
        for cat, (col, _) in CAT.items():
            s = d[d.category == cat]
            if not len(s):
                continue
            yy = np.minimum(s.y, CAP)
            over = s.y > CAP
            ax.scatter(s.x[~over], sign * yy[~over], s=26, c=col, alpha=.85,
                       edgecolors="white", linewidths=.4, zorder=3)
            if over.any():
                ax.scatter(s.x[over], sign * yy[over], s=52, c=col,
                           marker="^" if sign > 0 else "v",
                           edgecolors="black", linewidths=.45, zorder=4)

    for thr, sign in ((fg_thr, 1), (mt_thr, -1)):
        ax.axhline(sign * thr, ls="--", lw=.8, color="#666", zorder=2)
    ax.text(CUM * .998, fg_thr + .35, "FDR = 0.05", ha="right", fontsize=7, color="#555")
    ax.text(CUM * .998, -mt_thr - .95, "FDR = 0.05", ha="right", fontsize=7, color="#555")

    # 端到端对照：全研究最强关联在 MC1R 区；PARP1 复现已发表方向
    mc1r = fg.loc[fg.y.idxmax()]
    ax.annotate("MC1R region\n$P$ = 4×10⁻³⁷ — strongest in the study\n"
                "(and $P$ below float precision under the meta)",
                (mc1r.x, CAP), textcoords="offset points", xytext=(-16, 14),
                ha="right", fontsize=7, color="#8B2E28", linespacing=1.35,
                arrowprops=dict(arrowstyle="-", lw=.6, color="#BBB"))
    for d, sign in ((fg, 1), (mt, -1)):
        pp = d[d.SYMBOL == "PARP1"]
        if len(pp):
            r = pp.iloc[0]
            ax.annotate("PARP1", (r.x, sign * min(r.y, CAP)),
                        textcoords="offset points", xytext=(9, 3 * sign),
                        fontsize=6.8, style="italic", color="#555")

    # Counted from the data rather than typed in. The hardcoded version said
    # "7 independent loci", which was both the wrong word and the single-linkage
    # count -- the fixed-anchor partition this figure now draws gives 8.
    fg_n = int((FG.FDR < 0.05).sum())
    mt_sig = MT[MT.FDR < 0.05]
    mt_n, mt_loci = int(len(mt_sig)), mt_sig.locus.nunique()
    mt_novel = int((~mt_sig.known).sum())
    ax.text(CUM * .006, CAP + 3.6,
            f"FinnGen R12 round — {fg_n} significant records, "
            "every one on a known pigmentation or naevus locus",
            fontsize=7.8, color="#333")
    ax.text(CUM * .006, -CAP - 3.6,
            f"Meta outcome — {mt_n} significant records across {mt_loci} "
            f"bounded loci, {mt_novel} records on novel loci",
            fontsize=7.8, color="#333", va="top")

    ax.set_xticks([OFF[c] + CHRLEN[c] / 2 for c in map(str, range(1, 23))])
    ax.set_xticklabels(range(1, 23), fontsize=7)
    ax.set_xlim(0, CUM); ax.set_ylim(-CAP - 4.6, CAP + 4.6)
    ax.set_yticks([-8, -4, 0, 4, 8])
    ax.set_yticklabels([8, 4, 0, 4, 8], fontsize=7.5)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P_{\rm MR})$   of the strongest record per locus")
    ax.spines["left"].set_bounds(-CAP, CAP)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=5.5, color=c,
                              markeredgecolor="white", label=l)
                       for c, l in CAT.values()] +
                      [Line2D([], [], marker="^", ls="", ms=6, color="#777",
                              markeredgecolor="black", label="beyond the axis limit")],
              frameon=False, fontsize=7.2, loc="upper left", ncol=5,
              bbox_to_anchor=(0, -.105), columnspacing=1.4, handletextpad=.4)
    ax.set_title("a  One point per bounded locus, both outcomes",
                 fontsize=9.5, loc="left", pad=16)


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    fg, mt = counts(FG), counts(MT)

    for k, (c, lab) in enumerate(((fg, "FinnGen R12\nround"),
                                  (mt, "Meta\noutcome"))):
        x = k * 1.5
        ax.bar(x, c["sig_pct"], width=.62, color="#C4453C", alpha=.85, zorder=3)
        ax.bar(x + .70, c["bg_pct"], width=.62, color="#BBBBBB", alpha=.75, zorder=3)
        ax.text(x, c["sig_pct"] + 3.5, f"{c['sig_pct']:.1f}%", ha="center",
                fontsize=9, fontweight="bold", color="#C4453C")
        ax.text(x, c["sig_pct"] + 10.5, f"{c['sig_known']}/{c['sig_loci']}\nloci",
                ha="center", fontsize=7, color="#8B2E28", linespacing=1.3)
        ax.text(x + .70, c["bg_pct"] + 3.5, f"{c['bg_pct']:.1f}%", ha="center",
                fontsize=8, color="#777")
        ax.text(x + .70, c["bg_pct"] + 9.0,
                f"{c['bg_known']}/{c['bg_loci']}\nloci", ha="center",
                fontsize=6.6, color="#999", linespacing=1.3)
        ax.text(x + .35, -13,
                f"{c['fold']:.2f}×\nP = {c['p']:.3f}" if c["p"] >= 1e-3 else
                f"{c['fold']:.2f}×\nP = {c['p']:.0e}",
                ha="center", va="top", fontsize=8.2, color="#333", linespacing=1.35)
        ax.text(x + .35, -30, lab, ha="center", va="top", fontsize=8.4, color="#333")

    ax.set_xlim(-.65, 2.9); ax.set_ylim(-42, 128)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Loci carrying a known melanoma,\nnaevus or pigmentation lead SNP  (%)")
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_bounds(0, 100)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="#C4453C", alpha=.85),
                       plt.Rectangle((0, 0), 1, 1, fc="#BBBBBB", alpha=.75)],
              labels=["FDR < 0.05 loci", "all testable loci (background)"],
              frameon=False, fontsize=7.2, loc="upper center",
              bbox_to_anchor=(.5, 1.02))
    ax.text(.5, -.30,
            "Counted by independent locus, not by gene record.\n"
            "The record-level figure is larger (5.26×) and is\n"
            "reported only in the supplement.",
            transform=ax.transAxes, ha="center", va="top", fontsize=6.9,
            color="#666", linespacing=1.45)
    ax.set_title("b  Attribution, by locus", fontsize=9.5, loc="left", pad=16)


def main():
    fig = plt.figure(figsize=(14.4, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.85, 1], wspace=.20,
                          left=.055, right=.985, top=.865, bottom=.115)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    fig.suptitle("Significant signal sits on loci already known for the outcome",
                 fontsize=11.5, y=.955)
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "Fig1_locus_attribution" + ext), **kw)
    plt.close(fig)

    fg, mt = counts(FG), counts(MT)
    for lab, c in (("FinnGen R12", fg), ("meta", mt)):
        print(f"Fig1  {lab:<12} {c['sig_known']}/{c['sig_loci']} sig loci known "
              f"({c['sig_pct']:.1f}%), background {c['bg_known']}/{c['bg_loci']} "
              f"({c['bg_pct']:.1f}%)  ->  {c['fold']:.2f}x, P = {c['p']:.4g}")


if __name__ == "__main__":
    main()
