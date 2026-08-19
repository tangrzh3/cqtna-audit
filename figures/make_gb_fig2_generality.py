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
  123d_fixed_anchor_full_grid.tsv（全网格 + 错配对照，冻结分区）
  126a_offgrid_attribution.tsv（迁移检验等非网格格子，同一分区）
  59a_release_trajectory.tsv · 58a_finngen_reference_predictions.tsv · 99a_r13_trajectory.tsv

⚠ 不要再读 119a/119c/108a/99c —— 它们是单连锁分区上的数，
与本图其余部分不同单位。
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

# One table for the whole grid and its mismatched controls: 123d is every
# cell recomputed on the frozen partition (S36). The figure used to read
# 119a/119c for the grid and 108a/99c for the RA and transfer cells, which
# meant three panels drawn on two different statistical units.
grid = pd.read_csv(f"{MR}/123d_fixed_anchor_full_grid.tsv", sep="\t")
grid = grid[grid.analysis == "main"].copy()
grid = grid[~grid.cell.str.contains("MHC")]        # post-hoc rows are not the grid
grid["disease"] = grid.cell.str.split(" x ").str[0]
grid["exposure"] = grid.cell.str.split(" x ").str[1]
grid["void"] = grid.control.eq("FAILED")
offgrid = pd.read_csv(f"{MR}/126a_offgrid_attribution.tsv", sep="\t")
rel = pd.read_csv(f"{MR}/59a_release_trajectory.tsv", sep="\t")
rel = rel[rel["mode"] == "own"].sort_values("cases").reset_index(drop=True)
pred = pd.read_csv(f"{MR}/58a_finngen_reference_predictions.tsv", sep="\t")
r13 = pd.read_csv(f"{MR}/99a_r13_trajectory.tsv", sep="\t")
r13 = r13[(r13["mode"] == "own") & (r13.release == "R13")].iloc[0]

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
            void = bool(r.void)
            sig = (r.fisher_p < 0.05) and not void
            # A void cell is not a weak cell. It is drawn grey and hatched, with
            # its own fold withheld, because its mismatched-list control also
            # enriched and the pre-registration voids it on that alone.
            ax.add_patch(plt.Rectangle((j - .44, i - .38), .88, .76,
                                       facecolor="#EDEDED" if void else C_KNOWN,
                                       alpha=1 if void else .06 + .16 * min(r.fold / 11, 1),
                                       edgecolor="#999" if void else (C_KNOWN if sig else "#CCC"),
                                       hatch="////" if void else None,
                                       lw=1.2 if void else (1.6 if sig else .9),
                                       ls="-" if (sig or void) else "--", zorder=2))
            if void:
                # The hatching says "excluded"; the labels have to stay readable
                # over it, so they sit on their own opaque strip.
                ax.text(j, i + .15, "VOID", ha="center", va="center",
                        fontsize=12.5, fontweight="bold", color="#5A5A5A",
                        zorder=5,
                        bbox=dict(boxstyle="square,pad=0.16", fc="#EDEDED",
                                  ec="none"))
                ax.text(j, i - .16,
                        "mismatched control also\n"
                        f"enriches: {r.mismatch_fold:.2f}×, P = "
                        + f"{r.mismatch_p:.0e}".replace("e-", "×10⁻"),
                        ha="center", va="center", fontsize=6.8, color="#5A5A5A",
                        linespacing=1.4, zorder=5,
                        bbox=dict(boxstyle="square,pad=0.22", fc="#EDEDED",
                                  ec="none"))
                continue
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
    ax.text(1, 1.44, "One cell is void on its own mismatched-list control; the "
                     "other five all enrich, four\nreaching P < 0.05. The one that "
                     "does not is on a numerator of one of two loci.",
            ha="center", va="bottom", fontsize=7.3, color="#444", linespacing=1.4)
    ax.set_title("b  Both axes crossed: three diseases × two exposure resources",
                 fontsize=9.5, loc="left", pad=8)


# ------------------------------------------------------------------ panel c
def panel_c(ax):
    # Every row comes from the one grid table, so the matched and mismatched
    # bars of a row are two scorings of the same loci under one partition.
    def g(cell):
        return grid[grid.cell == cell].iloc[0]

    eq_mel = g("melanoma x eQTLGen_blood")
    hcc_lo = g("HCC_low x Soskic_CD4")
    ra_cd4 = g("RA x Soskic_CD4")
    ra_eq = g("RA x eQTLGen_blood")
    r13_m = offgrid[offgrid.cell == "C1 Soskic x melanoma R13"].iloc[0]

    items = [
        ("eQTLGen ×\nmelanoma", eq_mel.fold, eq_mel.fisher_p,
         eq_mel.mismatch_fold, eq_mel.mismatch_p, "HCC list", False),
        ("CD4⁺ ×\nmelanoma (R13)", r13_m.fold, r13_m.fisher_p,
         r13_m.mismatch_fold, r13_m.mismatch_p, "HCC list", False),
        ("CD4⁺ ×\nHCC-low", hcc_lo.fold, hcc_lo.fisher_p,
         hcc_lo.mismatch_fold, hcc_lo.mismatch_p, "melanoma list", False),
        ("CD4⁺ ×\nRA", ra_cd4.fold, ra_cd4.fisher_p,
         ra_cd4.mismatch_fold, ra_cd4.mismatch_p, "melanoma list", False),
        ("eQTLGen ×\nRA  (VOID)", ra_eq.fold, ra_eq.fisher_p,
         ra_eq.mismatch_fold, ra_eq.mismatch_p, "melanoma list", True),
    ]
    y = np.arange(len(items))[::-1]
    h = .34
    for k, (lab, mf, mp, xf, xp, which, void) in enumerate(items):
        yy = y[k]
        ax.barh(yy + h / 2 + .02, mf, height=h,
                color="#BFBFBF" if void else C_MATCH, alpha=.85,
                hatch="////" if void else None, zorder=3)
        ax.barh(yy - h / 2 - .02, xf, height=h, color=C_MISMATCH,
                alpha=.9 if void else .65,
                edgecolor="#C4453C" if void else "none",
                lw=1.2 if void else 0, zorder=3)
        ax.text(mf + .3, yy + h / 2 + .02,
                f"{mf:.2f}×  " + (f"P = {mp:.3f}" if mp >= 1e-3
                                  else f"P = {mp:.0e}".replace("e-", "×10⁻"))
                + ("   — not reportable" if void else ""),
                va="center", fontsize=7.3, color="#8A8A8A" if void else C_MATCH)
        ax.text(max(xf, 0) + .3, yy - h / 2 - .02,
                f"{xf:.2f}×  " + (f"P = {xp:.4f}" if xp < .01
                                   else f"P = {xp:.2f}") + f"   ({which})",
                va="center", fontsize=7.3,
                color="#C4453C" if void else "#777")

    ax.axvline(1, ls=":", lw=1, color="#888", zorder=2)
    ax.text(1.15, -.66, "1× = no enrichment", fontsize=6.8, color="#888")
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items], fontsize=7.8)
    ax.set_xlim(0, 26); ax.set_ylim(-.85, 4.62)
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
            "disease-specific. In RA it is not orthogonal at all —\n"
            "RA and melanoma share immune loci, and on whole blood\n"
            "the wrong list enriches outright, voiding that cell.",
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

    mc = grid[grid.role == "main"]
    ok = mc[~mc.void]
    print(f"Fig2_generality ok  |  {len(mc)} registered cells, "
          f"{int(mc.void.sum())} void on the mismatched control  |  "
          f"of the {len(ok)} usable: fold>1 {int((ok.fold>1).sum())}, "
          f"P<0.05 {int((ok.fisher_p<0.05).sum())}  |  "
          f"void: {', '.join(mc.cell[mc.void]) or 'none'}")


if __name__ == "__main__":
    main()
