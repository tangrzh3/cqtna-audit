"""GB Fig 2 — generality: the attribution finding under one change at a time.

正文口径（MANUSCRIPT_GB.md §"Significant signal sits on loci already known"
之后的一整节）。图注要求覆盖六件事：
  五个嵌套 release · 迁移检验 · 第二疾病 · 非癌结局 · 第二暴露资源 · 双轴交叉网格 + 错配位点对照

三个面板分摊：
  a  同一资源内改变功效：R8–R12 五个嵌套 release，加 R13 迁移检验
     （R13 **不是第六个功效层**——endpoint 换了，控制组排除规则也换了，
       故用竖线隔开、单独着色，绝不连成一条曲线）
  b  双轴交叉网格：三疾病 × 二资源 = 六格。第二疾病(HCC)与非癌结局(RA)
     是这张网格的两列，第二暴露资源是它的第二行——三件事在同一张图里读完。
     HCC 低功效层是**功效敏感性**，不是第七、八格，故置于分隔线之下。
  c  错配位点对照：把每格改用**错的那个疾病**的已知位点名单打分。
     没有这一栏，b 只能说明"显著位点落在某些位点上"，不能说明落在**它自己疾病**的位点上。

数据源（全部来自 step119 重建的权威表，**不要读 94d/94e**——
94d 里混着预注册判为作废的 lung/colorectal 两行，94e 汇总的是含那两行的旧网格）：
  119a_grid_main.tsv · 119c_mismatch_controls.tsv
  59a_release_trajectory.tsv · 58a_finngen_reference_predictions.tsv
  99a_r13_trajectory.tsv · 99c_r13_attribution.tsv · 108a_ra_attribution.tsv
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
C_TRANSFER = "#8A6BBE"
C_MATCH, C_MISMATCH = "#C4453C", "#9A9A9A"

grid = pd.read_csv(f"{MR}/119a_grid_main.tsv", sep="\t")
mism = pd.read_csv(f"{MR}/119c_mismatch_controls.tsv", sep="\t")
rel = pd.read_csv(f"{MR}/59a_release_trajectory.tsv", sep="\t")
rel = rel[rel["mode"] == "own"].sort_values("cases").reset_index(drop=True)
pred = pd.read_csv(f"{MR}/58a_finngen_reference_predictions.tsv", sep="\t")
r13 = pd.read_csv(f"{MR}/99a_r13_trajectory.tsv", sep="\t")
r13 = r13[(r13["mode"] == "own") & (r13.release == "R13")].iloc[0]
r13a = pd.read_csv(f"{MR}/99c_r13_attribution.tsv", sep="\t")
ra = pd.read_csv(f"{MR}/108a_ra_attribution.tsv", sep="\t")

DIS_LAB = {"melanoma": "Melanoma", "HCC_high": "HCC\n(3,748 cases)",
           "RA": "Rheumatoid\narthritis", "HCC_low": "HCC\n(947 cases)"}
EXP_LAB = {"Soskic_CD4": "CD4⁺ T activation\ntime course",
           "eQTLGen_blood": "eQTLGen\nwhole blood"}


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    m = rel.merge(pred, on=["release", "cases"], suffixes=("", "_p"))
    x = m.cases / 1000

    ax.fill_between(x, m.pred_hits_lo, m.pred_hits_hi, color=C_HITS, alpha=.12,
                    lw=0, zorder=1, label="registered prediction (5–95%)")
    ax.plot(x, m.n_hits, "-o", ms=6, lw=1.9, color=C_HITS, zorder=4,
            label="discoveries at FDR < 0.05")
    for xi, yi, r in zip(x, m.n_hits, m.release):
        ax.annotate(r, (xi, yi), textcoords="offset points", xytext=(0, 8),
                    ha="center", fontsize=7, color="#555")

    # R13 不是第六个功效层：endpoint 与控制组排除规则都换了
    xr = r13.cases / 1000
    ax.axvline((5.753 + xr) / 2, ls="--", lw=1, color="#BBB", zorder=1)
    ax.scatter([xr], [r13.n_hits], marker="D", s=62, color=C_TRANSFER,
               edgecolors="white", linewidths=.8, zorder=5)
    ax.annotate("R13", (xr, r13.n_hits), textcoords="offset points",
                xytext=(0, 9), ha="center", fontsize=7, color=C_TRANSFER)
    ax.text((5.753 + xr) / 2 - .07, 7.4,
            "endpoint redefined →\ntransfer test, not a\nsixth power level",
            fontsize=6.9, color=C_TRANSFER, va="top", ha="right",
            linespacing=1.35)

    ax.text(.035, .965,
            "No novel-locus gene reaches FDR < 0.05\n"
            "in any release. The same five genes are\n"
            "the significant list at three consecutive\n"
            "releases without change.",
            transform=ax.transAxes, va="top", fontsize=7.3, color="#444",
            linespacing=1.45,
            bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#DDD", lw=.6))

    ax.set_ylim(0, 14.5)
    ax.set_xlabel("Melanoma cases (thousands)")
    ax.set_ylabel("Discoveries at FDR < 0.05")
    ax.legend(frameon=False, fontsize=7.2, loc="lower left")
    ax.set_title("a  One resource, five nested power levels — and a transfer test",
                 fontsize=9.5, loc="left", pad=6)


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    diseases = ["melanoma", "HCC_high", "RA"]
    resources = ["Soskic_CD4", "eQTLGen_blood"]

    def cell(e, d):
        return grid[(grid.exposure == e) & (grid.disease == d)].iloc[0]

    for j, d in enumerate(diseases):
        for i, e in enumerate(resources):
            r = cell(e, d)
            sig = r.fisher_p < 0.05
            ax.add_patch(plt.Rectangle((j - .44, i - .38), .88, .76,
                                       facecolor=C_KNOWN,
                                       alpha=.06 + .16 * min(r.fold / 9, 1),
                                       edgecolor=C_KNOWN if sig else "#CCC",
                                       lw=1.6 if sig else .9,
                                       ls="-" if sig else "--", zorder=2))
            ax.text(j, i + .13, f"{r.fold:.2f}×", ha="center", va="center",
                    fontsize=13.5, fontweight="bold",
                    color=C_KNOWN if sig else "#8A8A8A", zorder=4)
            pf = (f"P = {r.fisher_p:.3f}" if r.fisher_p >= 1e-3
                  else f"P = {r.fisher_p:.0e}".replace("e-", "×10⁻"))
            ax.text(j, i - .10, pf, ha="center", va="center", fontsize=7.6,
                    color="#444", zorder=4)
            ax.text(j, i - .26, f"{r.sig_known}/{r.sig_loci} loci",
                    ha="center", va="center", fontsize=6.8, color="#888", zorder=4)

    # 功效敏感性行，明确隔开
    ax.axhline(-.62, ls="-", lw=.9, color="#888")
    for j, d in enumerate(["HCC_low"]):
        for i, e in enumerate(resources):
            r = grid[(grid.exposure == e) & (grid.disease == d)].iloc[0]
            ax.text(1 + (i - .5) * .96, -.90, f"{r.fold:.2f}×  (P = {r.fisher_p:.3f})",
                    ha="center", va="center", fontsize=8, color="#777")
            ax.text(1 + (i - .5) * .96, -1.11, EXP_LAB[e].replace("\n", " "),
                    ha="center", va="center", fontsize=6.6, color="#AAA")
    ax.text(-.44, -.90, "Power sensitivity\nHCC at 947 cases", ha="left",
            va="center", fontsize=7.2, color="#777", linespacing=1.35)
    ax.text(2.62, -1.40,
            "not a seventh and eighth cell — the two HCC "
            "levels come from one resource (S20 §9)",
            ha="right", va="center", fontsize=6.6, color="#AAA")

    ax.set_xticks(range(3))
    ax.set_xticklabels([DIS_LAB[d] for d in diseases], fontsize=8.4)
    ax.set_yticks(range(2))
    ax.set_yticklabels([EXP_LAB[e] for e in resources], fontsize=8.4)
    ax.set_xlim(-.62, 2.62); ax.set_ylim(-1.52, 1.72)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(1, 1.44, "All six enrich; four reach P < 0.05.  The two that do not "
                     "are the same cell —\nHCC at its higher power — on numerators "
                     "as small as one of two loci.",
            ha="center", va="bottom", fontsize=7.3, color="#444", linespacing=1.4)
    ax.set_title("b  Both axes crossed: three diseases × two exposure resources",
                 fontsize=9.5, loc="left", pad=8)


# ------------------------------------------------------------------ panel c
def panel_c(ax):
    nc1 = mism[mism.control == "NC1"].iloc[0]
    nc2 = mism[mism.control == "NC2"].iloc[0]
    ra_m = ra[ra.cell == "N1 Soskic x RA"].iloc[0]
    ra_nc = ra[ra.cell == "NC Soskic x RA scored with melanoma list"].iloc[0]
    r13_m = r13a[r13a.cell == "C1 Soskic x melanoma R13"].iloc[0]
    r13_nc = r13a[r13a.cell == "NC C1 R13 scored with HCC list"].iloc[0]

    items = [
        ("eQTLGen ×\nmelanoma", nc1.matched_fold, nc1.matched_p,
         nc1.fold, nc1.fisher_p, "HCC list"),
        ("CD4⁺ ×\nmelanoma (R13)", r13_m.fold, r13_m.fisher_p,
         r13_nc.fold, r13_nc.fisher_p, "HCC list"),
        ("CD4⁺ ×\nHCC-low", nc2.matched_fold, nc2.matched_p,
         nc2.fold, nc2.fisher_p, "melanoma list"),
        ("CD4⁺ ×\nRA", ra_m.fold, ra_m.fisher_p,
         ra_nc.fold, ra_nc.fisher_p, "melanoma list"),
    ]
    y = np.arange(len(items))[::-1]
    h = .34
    for k, (lab, mf, mp, xf, xp, which) in enumerate(items):
        yy = y[k]
        ax.barh(yy + h / 2 + .02, mf, height=h, color=C_MATCH, alpha=.85, zorder=3)
        ax.barh(yy - h / 2 - .02, xf, height=h, color=C_MISMATCH, alpha=.65, zorder=3)
        ax.text(mf + .3, yy + h / 2 + .02,
                f"{mf:.2f}×  " + (f"P = {mp:.3f}" if mp >= 1e-3
                                  else f"P = {mp:.0e}".replace("e-", "×10⁻")),
                va="center", fontsize=7.3, color=C_MATCH)
        ax.text(max(xf, 0) + .3, yy - h / 2 - .02,
                f"{xf:.2f}×  P = {xp:.2f}   ({which})",
                va="center", fontsize=7.3, color="#777")

    ax.axvline(1, ls=":", lw=1, color="#888", zorder=2)
    ax.text(1.15, -.66, "1× = no enrichment", fontsize=6.8, color="#888")
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items], fontsize=7.8)
    ax.set_xlim(0, 26); ax.set_ylim(-.85, 3.62)
    ax.set_xlabel("Fold enrichment on known outcome loci")
    ax.tick_params(axis="y", length=0)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=C_MATCH, alpha=.85),
                       plt.Rectangle((0, 0), 1, 1, fc=C_MISMATCH, alpha=.65)],
              labels=["scored against the outcome's own known loci",
                      "scored against the wrong disease's list"],
              frameon=False, fontsize=7.2, loc="lower right")
    ax.text(.995, .615,
            "A mismatched list is a narrower control than it looks:\n"
            "it rules out enrichment on loci indiscriminately dense\n"
            "across diseases, not a density that is itself\n"
            "disease-specific. RA is the least clean of the four —\n"
            "RA and melanoma share immune loci across the MHC.",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.7,
            color="#666", linespacing=1.45,
            bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#DDD", lw=.6))
    ax.set_title("c  The same cells scored against the wrong disease's known loci",
                 fontsize=9.5, loc="left", pad=6)


def main():
    fig = plt.figure(figsize=(15.2, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.02], hspace=.34, wspace=.24,
                          left=.075, right=.975, top=.895, bottom=.065)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, :]))
    fig.suptitle("Changing one thing at a time: outcome power, disease, "
                 "and exposure resource",
                 fontsize=11.5, y=.955)
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "Fig2_generality" + ext), **kw)
    plt.close(fig)

    mc = grid[grid.status == "main"]
    print(f"Fig2_generality ok  |  main grid {len(mc)} cells, "
          f"fold>1 {int((mc.fold>1).sum())}, P<0.05 {int((mc.fisher_p<0.05).sum())}  |  "
          f"mismatch NC1 {mism[mism.control=='NC1'].iloc[0].fold}x, "
          f"NC2 {mism[mism.control=='NC2'].iloc[0].fold}x")


if __name__ == "__main__":
    main()
