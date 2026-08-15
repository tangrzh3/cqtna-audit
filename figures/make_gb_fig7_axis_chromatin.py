"""GB Fig 7 — the CD4⁺ metabolic axis and its chromatin signature.

⚠⚠ 与正文现有措辞不符的一处，画图时按**表**画，不按正文画：

   GB §"What the instrumented gene marks" 与 v2 的 Fig9 图注都写着
   "**none of the eleven biological modules tested exceeds a composition-matched
   null**"。`43c_residual_axis_verdict.tsv` 说的不是这个：11 个生物学模块里
   **有 2 个超过**——**Proliferation**（|SMD| = 0.118）与 **OXPHOS**（0.212），
   零模型阈值为 50 个随机模块 |SMD| 的 95 分位数 = 0.0720。

   而且 Proliferation 是 `step43_glyco_vs_activation.py` 里**预先写死的阳性对照**：
   脚本明写"若阳性对照不超过零模型，则该检验功效不足、其零模型无信息"。
   也就是说 **Proliferation 超过零模型是检验成立的必要条件**，不是反例。

   正文真正要说的那件事仍然成立，而且比原句更强：
   **Activation 模块自己就贴在零模型上（SMD = −0.003）**，
   即这条轴不是"活化强度"的换个说法；超过零模型的两个是阳性对照与 OXPHOS，
   后者与"合成代谢状态"的读法一致而非相左。

   → 已记入 FIGURES_GB_mapping.md §五，正文须改。本图 panel b 照表画全部 11 个模块。

面板：
  a  残差轴顶端的家族组成（44b）：糖酵解 21.3 倍、核糖体蛋白 10.3 倍
  b  组成匹配零模型：11 个生物学模块的 |SMD| 对 95 分位阈值（43b/43c）
  c  轴高 vs 背景 peak 的 motif 富集（48b）
  d  全部 peak 对**活化不变** peak 的 motif odds（50b）：AP-1 在扣掉活化后仍在

⚠ panel c/d 是**轴层面的染色质刻画**，与修订中已撤回的
   "工具变量经破坏 AP-1 motif 起作用"（Step 65，经验 P = 1.0）**不是同一件事**。
"""
import os

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
RED, BLUE, PURPLE = "#C4453C", "#3B7DD8", "#8E6BB3"
C_NULL, C_OVER, C_CTRL = "#BBBBBB", "#C4453C", "#2E7D5B"

fam = pd.read_csv(f"{MR}/44b_family_composition.tsv", sep="\t")
split = pd.read_csv(f"{MR}/43b_residual_split_modules.tsv", sep="\t")
ver = pd.read_csv(f"{MR}/43c_residual_axis_verdict.tsv", sep="\t")
mot_axis = pd.read_csv(f"{MR}/48b_motif_axis.tsv", sep="\t")
mot_cmp = pd.read_csv(f"{MR}/50b_motif_comparison_all_vs_invariant.tsv", sep="\t")

# 零模型阈值：50 个随机模块 |SMD| 的 95 分位（与 step43 同式重算，已核对逐行一致）
_agg = split.groupby(["module", "kind"]).smd.mean().reset_index()
NULL_THR = float(np.percentile(_agg[_agg.kind == "random"].smd.abs(), 95))
N_RANDOM = int((_agg.kind == "random").sum())
POSITIVE_CONTROL = "Proliferation"       # step43 里预先写死的阳性对照


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    f = fam[fam.n_in_top > 0].sort_values("enrichment").tail(9)
    cols = [RED if d.lower().startswith("up") else BLUE for d in f.direction]
    y = np.arange(len(f))
    ax.barh(y, f.enrichment, color=cols, height=.66, zorder=3)
    for i, (e, n) in enumerate(zip(f.enrichment, f.n_in_top)):
        ax.text(e + .45, i, f"{e:.1f}×  (n = {int(n)})", va="center", fontsize=7.3,
                color="#444")
    ax.axvline(1, ls="--", lw=1, color="#666", zorder=2)
    ax.set_yticks(y)
    # interferon 等家族在两端各出现一次，标签必须带方向，否则看起来是重复行
    dup = f.family.duplicated(keep=False)
    ax.set_yticklabels(
        [fm.replace("_", " ") + (f" ({d})" if k else "")
         for fm, d, k in zip(f.family, f.direction, dup)], fontsize=8)
    ax.set_xlim(0, f.enrichment.max() * 1.42)
    ax.set_xlabel("Enrichment among the axis-defining genes (fold over background)")
    ax.tick_params(axis="y", length=0)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=8, color=RED,
                              label="glycolysis-high end"),
                       Line2D([], [], marker="s", ls="", ms=8, color=BLUE,
                              label="glycolysis-low end")],
              frameon=False, fontsize=7.6, loc="lower right")
    ax.set_title("a  What the axis is made of — no pre-specified modules",
                 fontsize=9.5, loc="left", pad=6)


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    v = ver.sort_values("smd", key=np.abs).reset_index(drop=True)
    y = np.arange(len(v))
    for i, r in v.iterrows():
        if r.module == POSITIVE_CONTROL:
            c = C_CTRL
        elif r.exceeds_null:
            c = C_OVER
        else:
            c = C_NULL
        ax.barh(i, abs(r.smd), color=c, height=.62,
                alpha=.9 if r.exceeds_null else .65, zorder=3)

    ax.axvline(NULL_THR, ls="--", lw=1.2, color="#333", zorder=4)
    ax.text(.155, 5.4,
            f"composition-matched null:\n95th percentile of {N_RANDOM} random\n"
            f"modules, |SMD| = {NULL_THR:.3f}",
            fontsize=6.9, color="#333", va="center", linespacing=1.45,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#DDD", lw=.6))

    ax.set_yticks(y)
    ax.set_yticklabels([m.replace("_", " ") for m in v.module], fontsize=8)
    for i, r in v.iterrows():
        if r.module == POSITIVE_CONTROL:
            ax.text(abs(r.smd) + .006, i, "pre-specified positive control —\n"
                                          "must exceed, or the null is uninformative",
                    va="center", fontsize=6.7, color=C_CTRL, linespacing=1.35)
        elif r.exceeds_null:
            ax.text(abs(r.smd) + .006, i, "exceeds", va="center", fontsize=7,
                    color=C_OVER)
    act = v[v.module == "Activation"].iloc[0]
    ax.annotate("Activation itself sits at the null (SMD = %.3f) —\n"
                "the axis is not a restatement of activation strength" % act.smd,
                (abs(act.smd), v[v.module == "Activation"].index[0]),
                textcoords="offset points", xytext=(108, 4), fontsize=7,
                color="#333", va="center", linespacing=1.4,
                arrowprops=dict(arrowstyle="-", lw=.6, color="#BBB"))

    ax.set_xlim(0, .30)
    ax.set_xlabel("|standardised mean difference| between the axis ends")
    ax.tick_params(axis="y", length=0)
    ax.set_title("b  Eleven biological modules against a composition-matched null",
                 fontsize=9.5, loc="left", pad=6)


# ------------------------------------------------------------------ panel c
def panel_c(ax):
    top = (mot_axis[mot_axis.contrast == "glyco_high"]
           .sort_values("odds", ascending=False).head(8).iloc[::-1])
    y = np.arange(len(top))
    ax.barh(y, top.odds.values, color=PURPLE, height=.66, zorder=3)
    for i, (o, q) in enumerate(zip(top.odds.values, top.FDR.values)):
        ax.text(o + .015, i, f"FDR {q:.0e}".replace("e-0", "e−"), va="center",
                fontsize=6.9, color="#555")
    ax.axvline(1, ls=":", lw=1, color="#999", zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([m[:16] for m in top.motif], fontsize=7.6)
    ax.set_xlim(0, top.odds.max() * 1.28)
    ax.set_xlabel("Motif odds ratio, axis-high versus background peaks")
    ax.tick_params(axis="y", length=0)
    ax.set_title("c  A concordant chromatin signature", fontsize=9.5,
                 loc="left", pad=6)


# ------------------------------------------------------------------ panel d
def panel_d(ax):
    m = mot_cmp.dropna(subset=["odds_all", "odds_inv"]).copy()
    ap1 = m.motif.str.contains("JUN|FOS|BATF|JDP2", case=False, na=False)
    irf = m.motif.str.contains("IRF|STAT", case=False, na=False)

    ax.scatter(m.odds_all[~(ap1 | irf)], m.odds_inv[~(ap1 | irf)], s=13,
               color="#D5D8DC", zorder=2, label="other motifs")
    ax.scatter(m.odds_all[irf], m.odds_inv[irf], s=34, color=BLUE, zorder=4,
               edgecolors="white", linewidths=.4, label="IRF / STAT family")
    ax.scatter(m.odds_all[ap1], m.odds_inv[ap1], s=38, color=RED, zorder=5,
               edgecolors="white", linewidths=.4, label="AP-1 family")

    lim = [m[["odds_all", "odds_inv"]].min().min() * .95,
           m[["odds_all", "odds_inv"]].max().max() * 1.05]
    ax.plot(lim, lim, ls=":", lw=1, color="#888", zorder=1)
    ax.axhline(1, lw=.7, color="#BBB"); ax.axvline(1, lw=.7, color="#BBB")
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel("Odds ratio, all peaks")
    ax.set_ylabel("Odds ratio, activation-invariant peaks")
    ax.text(.97, .05, "below the diagonal = weakened\nonce activation is held constant",
            transform=ax.transAxes, ha="right", fontsize=7, color="#666",
            linespacing=1.4)
    ax.legend(frameon=False, fontsize=7.6, loc="upper left")
    ax.set_title("d  AP-1 enrichment persists in activation-invariant peaks",
                 fontsize=9.5, loc="left", pad=6)


def main():
    fig = plt.figure(figsize=(13.6, 9.6))
    gs = fig.add_gridspec(2, 2, hspace=.30, wspace=.34,
                          left=.105, right=.985, top=.905, bottom=.115)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, 0]))
    panel_d(fig.add_subplot(gs[1, 1]))
    fig.suptitle("What the instrumented gene marks: a definable anabolic state, "
                 "not a restatement of activation strength", fontsize=11.5, y=.955)
    fig.text(.105, .022,
             "Panels c and d characterise the chromatin of the axis.  They are NOT "
             "the claim, withdrawn during revision, that the instrument acts by "
             "disrupting an AP-1 motif — that test\nreturned an empirical "
             "$P$ = 1.0 and is reported nowhere in this paper.  Rebuilt without TPI1, "
             "which had entered both the defining score and the enrichment family, "
             "the glycolysis\nenrichment is 21.0-fold on 14 genes against 21.3 on 15.",
             ha="left", va="bottom", fontsize=7.1, color="#666", linespacing=1.6)
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "Fig7_axis_chromatin" + ext), **kw)
    plt.close(fig)

    ex = ver[ver.exceeds_null]
    gly = fam[(fam.direction == "up") & (fam.family == "glycolysis")].iloc[0]
    rib = fam[(fam.direction == "up") & (fam.family == "ribosomal_protein")].iloc[0]
    print(f"Fig7_axis_chromatin ok  |  glycolysis {gly.enrichment:.1f}x "
          f"({int(gly.n_in_top)} genes), ribosomal {rib.enrichment:.1f}x  |  "
          f"null threshold {NULL_THR:.4f} ({N_RANDOM} random modules)")
    print(f"  modules exceeding the null: {len(ex)} of {len(ver)} — "
          f"{', '.join(ex.module)}  "
          f"(Proliferation is the pre-specified positive control)")
    print(f"  *** the manuscript says 'none of the eleven ... exceeds' — that is "
          f"wrong; see FIGURES_GB_mapping.md §五")


if __name__ == "__main__":
    main()
