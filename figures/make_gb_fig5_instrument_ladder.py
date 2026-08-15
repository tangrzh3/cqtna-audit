"""GB Fig 5 — instrument availability across the glycolytic pathway, at three levels.

正文口径（MANUSCRIPT_GB.md §"What can be interrogated, stated at the level it applies to"）：
28 个糖酵解基因中 **3 个在本暴露资源里可工具化**（某活化状态下有全基因组显著 cis-eQTL），
**2 个在与本结局 harmonise 后仍可分析**，**1 个给出名义关联**。

三个数不是同一类陈述，图必须把它们分开画，不能压成"一个可用工具变量"：
  · 第一级只关乎通路 + 暴露数据，换任何结局都不变；
  · 第二级加了"与某个特定结局 GWAS harmonise"这一条件（SLC2A1 栽在这里，
    换个结局未必栽 → 不等于 SLC2A1 没有工具变量）；
  · 第三级才加上结局自己的显著性。

面板：
  a  28 基因 × 8 活化状态的最强 cis-eQTL P；虚线 = 5e-8。资源里根本没测到的基因单列。
  b  三级阶梯，逐级标出掉队的基因与掉队的原因。

数据源：75a_glycolysis_instrument_matrix.tsv（第一级）、
        34a_glycolysis_MR_records.tsv / 34b_glycolysis_MR_by_gene.tsv（第二、三级）。
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "figures")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": .8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150,
})
C_INST, C_MISS, C_DROP = "#C4453C", "#BBBBBB", "#3B7DD8"
GWS = 5e-8

PROFILES = ["Naive_0h", "Naive_16h", "Naive_40h", "Naive_5d",
            "Memory_0h", "Memory_16h", "Memory_40h", "Memory_5d"]
PROF_LAB = ["0 h", "16 h", "40 h", "5 d"] * 2

m = pd.read_csv(f"{MR}/75a_glycolysis_instrument_matrix.tsv", sep="\t")
rec = pd.read_csv(f"{MR}/34a_glycolysis_MR_records.tsv", sep="\t")
byg = pd.read_csv(f"{MR}/34b_glycolysis_MR_by_gene.tsv", sep="\t")

N_PATHWAY = len(m)
INSTRUMENTED = m.loc[m.has_instrument == 1, "SYMBOL"].tolist()
ANALYSABLE = byg.SYMBOL.tolist()
NOMINAL = byg.loc[byg.best_p < 0.05, "SYMBOL"].tolist()
# 在资源里完全没有 cis-eQTL 记录的基因（整行空白）——不是"没通过阈值"，是没测到
UNTESTED = m.loc[m[PROFILES].isna().all(axis=1), "SYMBOL"].tolist()


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    d = m.copy()
    d["minp"] = d[PROFILES].min(axis=1)
    # 有数据的按最强 P 排序在上，没测到的沉底
    tested = d[~d.SYMBOL.isin(UNTESTED)].sort_values("minp")
    untested = d[d.SYMBOL.isin(UNTESTED)].sort_values("SYMBOL")
    d = pd.concat([tested, untested], ignore_index=True)

    mat = -np.log10(d[PROFILES].astype(float).values)
    ax.imshow(np.ma.masked_invalid(mat), aspect="auto", cmap="RdPu",
              vmin=0, vmax=12, interpolation="nearest")

    # 没测到的格子画成灰底斜纹，避免与"测了但不显著"混为一谈
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isnan(mat[i, j]):
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, facecolor="#F0F0F0",
                                           edgecolor="white", lw=.5, hatch="///"))
    # 达到全基因组显著的格子描红框
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if not np.isnan(mat[i, j]) and mat[i, j] >= -np.log10(GWS):
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, facecolor="none",
                                           edgecolor=C_INST, lw=1.8, zorder=5))

    ax.axvline(3.5, color="#333", lw=1.2)
    ax.set_xticks(range(8)); ax.set_xticklabels(PROF_LAB, fontsize=7.5)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(
        [f"$\\bf{{{s}}}$" if s in INSTRUMENTED else s for s in d.SYMBOL],
        fontsize=7, fontstyle="italic")
    ax.text(1.5, -1.35, "Naive", ha="center", fontsize=8.5, color="#333")
    ax.text(5.5, -1.35, "Memory", ha="center", fontsize=8.5, color="#333")
    ax.set_ylim(len(d) - .5, -.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("a  Strongest cis-eQTL per gene and activation profile",
                 fontsize=9.5, loc="left", pad=22)

    cb = ax.figure.colorbar(ax.images[0], ax=ax, fraction=.035, pad=.02)
    cb.set_label(r"$-\log_{10}(P_{\rm eQTL})$", fontsize=8, labelpad=2)
    cb.ax.axhline(-np.log10(GWS), color=C_INST, lw=1.5)
    cb.ax.text(-.35, -np.log10(GWS), "5×10⁻⁸ ", va="center", ha="right", fontsize=7,
               color=C_INST, transform=cb.ax.get_yaxis_transform())
    cb.ax.tick_params(labelsize=7)

    ax.plot([], [], "s", ms=8, mfc="none", mec=C_INST, mew=1.8,
            label="genome-wide significant cis-eQTL")
    ax.plot([], [], "s", ms=8, mfc="#F0F0F0", mec="#CCC",
            label=f"no cis-eQTL record in the resource ({len(UNTESTED)} genes)")
    ax.legend(frameon=False, fontsize=7.2, loc="upper left",
              bbox_to_anchor=(0, -.055), ncol=2, handletextpad=.6,
              columnspacing=1.6)


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    tpi1 = rec[rec.SYMBOL == "TPI1"].iloc[0]
    eno1 = rec[rec.SYMBOL == "ENO1"].iloc[0]

    steps = [
        (N_PATHWAY, "Glycolytic genes in the pathway definition",
         "a property of the pathway", "#555555"),
        (len(INSTRUMENTED), "Instrumentable in this exposure resource",
         "a property of the pathway and the exposure data;\n"
         "unchanged whatever outcome were analysed", C_INST),
        (len(ANALYSABLE), "Still analysable against this outcome",
         "adds harmonisation with one particular outcome GWAS", "#B5793A"),
        (len(NOMINAL), "Yields a nominal association",
         "adds the outcome's own significance", C_DROP),
    ]
    drops = [
        (f"{N_PATHWAY - len(INSTRUMENTED)} genes carry no genome-wide significant\n"
         f"cis-eQTL in any of the eight profiles — including PGAM1, the\n"
         f"functionally strongest enzyme at baseline\n"
         f"({len(UNTESTED)} of them have no cis-eQTL record at all)"),
        ("SLC2A1 fails harmonisation with this outcome.\n"
         "That is an outcome-side failure: it is not a statement\n"
         "that SLC2A1 lacks an instrument"),
        (f"ENO1 is analysable but not associated\n"
         f"(P = {eno1.pval:.2f}, FDR = {eno1.FDR:.2f})"),
    ]
    genes_at = [None, INSTRUMENTED, ANALYSABLE, NOMINAL]

    ymax = 4.0
    for k, (n, title, kind, col) in enumerate(steps):
        y = ymax - k * 1.15
        w = 2.05 * np.sqrt(n / N_PATHWAY) + .55
        ax.add_patch(plt.Rectangle((-w / 2, y - .27), w, .54, facecolor=col,
                                   alpha=.16, edgecolor=col, lw=1.3, zorder=3))
        ax.text(0, y, f"{n}", ha="center", va="center", fontsize=15,
                fontweight="bold", color=col, zorder=4)
        ax.text(w / 2 + .18, y + .11, title, ha="left", va="center",
                fontsize=8.6, color="#222")
        ax.text(w / 2 + .18, y - .13, kind, ha="left", va="top",
                fontsize=7.2, color="#777", linespacing=1.35)
        if genes_at[k]:
            ax.text(-w / 2 - .18, y, ", ".join(genes_at[k]), ha="right", va="center",
                    fontsize=7.4, fontstyle="italic", color=col)
        if k < len(drops):
            ax.add_patch(FancyArrowPatch((0, y - .30), (0, y - .86),
                                         arrowstyle="-|>", mutation_scale=11,
                                         lw=1.1, color="#999", zorder=2))
            ax.text(-.30, y - .58, drops[k], ha="right", va="center",
                    fontsize=6.9, color="#666", linespacing=1.4)

    ax.text(.95, ymax - 3 * 1.15 - .62,
            "Compressing these into “one usable instrument” would attribute\n"
            "two outcome-side limitations to the biology of the pathway.",
            ha="center", va="top", fontsize=7.6, color="#333",
            bbox=dict(boxstyle="round,pad=0.4", fc="#F7F7F7", ec="#DDD", lw=.6))

    ax.set_xlim(-4.5, 6.4); ax.set_ylim(-.5, ymax + .62)
    ax.axis("off")
    ax.set_title("b  The same pathway, counted at three levels",
                 fontsize=9.5, loc="left", pad=8)


def main():
    fig = plt.figure(figsize=(13.2, 6.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.32], wspace=.06,
                          left=.065, right=.995, top=.87, bottom=.115)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    fig.suptitle("What a pathway-level question can actually be asked of: "
                 "instrument availability at three levels",
                 fontsize=11, y=.965)
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "Fig5_instrument_ladder" + ext), **kw)
    plt.close(fig)
    print(f"Fig5_instrument_ladder ok  |  pathway {N_PATHWAY} -> "
          f"instrumentable {len(INSTRUMENTED)} {INSTRUMENTED} -> "
          f"analysable {len(ANALYSABLE)} {ANALYSABLE} -> "
          f"nominal {len(NOMINAL)} {NOMINAL}  |  untested {len(UNTESTED)}")


if __name__ == "__main__":
    main()
