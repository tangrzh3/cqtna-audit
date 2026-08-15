"""GB Fig 9 — patients across three cohorts.

正文承重的说法是**"两个 arm 都在发现队列显著，但没有任何一个 arm 被第二次确认，
且两次失败是镜像的"**：
  · 治疗前 arm 在**同病种**队列量级吻合（+0.55 对 +0.74）但功效不足，在跨病种队列消失；
  · 治疗后 arm 在**同病种**队列翻向，却在跨病种队列以完整量级重现（+0.874，P = 0.043）。

因此这张图必须同时画出三件事，少一件就会被读成"复制成功"或"复制失败"：
  a  三个队列 × 两个 arm 的效应量与单侧置换 P；
  b  ★ 设计上的 P 下限。跨病种队列的**治疗前** arm 是 4 名应答者对 **2** 名非应答者，
     无论效应量多大都**不可能**达到 P < 0.05——这是在读任何表达值之前就算出来并
     预注册的（S21）。把它和"检验了但不显著"画成同一种东西是错的。

口径（三个队列各自的主定义，逐行记明来源，不得混用）：
  发现队列   74a_patient_level_inference.tsv   locked16 | all samples
  同病种复制 79a_pozniak_authoritative.tsv     CD4_Tcells（作者注释，主）· locked16 · all
  跨病种     87b_arm_results.tsv               tissue T · cd4_def nonTreg · sig16
             （87b 的 nonTreg 行才是正文引的 +0.874；all 行为 +0.869）
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
C_PRE, C_POST = "#7A5AA8", "#E8A33D"
C_FLOOR = "#B03A3A"

disc = pd.read_csv(f"{MR}/74a_patient_level_inference.tsv", sep="\t")
poz = pd.read_csv(f"{MR}/79a_pozniak_authoritative.tsv", sep="\t")
hcc = pd.read_csv(f"{MR}/87b_arm_results.tsv", sep="\t")


def d_disc(tp):
    r = disc[(disc.timepoint == tp) & (disc.genes == 16) &
             (disc.subset == "locked16 | all samples")].iloc[0]
    return dict(delta=r.score_diff, lo=r.ci_lo, hi=r.ci_hi, p=r.perm_p_score,
                n_R=int(r.n_R), n_NR=int(r.n_NR), floor=np.nan)


def d_poz(tp):
    r = poz[(poz.definition.str.startswith("CD4_Tcells (")) &
            (poz.timepoint == tp) & (poz.genes == "locked16") &
            (poz.subset == "all")].iloc[0]
    return dict(delta=r.score_diff, lo=r.ci_lo, hi=r.ci_hi, p=r.perm_p,
                n_R=int(r.n_R), n_NR=int(r.n_NR), floor=np.nan)


def d_hcc(hyp):
    r = hcc[(hcc.hypothesis == hyp) & (hcc.tissue == "T") &
            (hcc.cd4_def == "nonTreg") & (hcc.signature == "sig16")].iloc[0]
    return dict(delta=r.delta_NR_minus_R, lo=np.nan, hi=np.nan,
                p=r.p_perm_onesided, n_R=int(r.n_R), n_NR=int(r.n_NR),
                floor=r.min_possible_p)


ROWS = [
    ("Discovery\nmelanoma, anti-PD-1/CTLA-4", "post", d_disc("Post")),
    ("Discovery\nmelanoma, anti-PD-1/CTLA-4", "pre", d_disc("Pre")),
    ("Same-disease replication\nmelanoma, anti-PD-1", "post", d_poz("OT")),
    ("Same-disease replication\nmelanoma, anti-PD-1", "pre", d_poz("BT")),
    ("Cross-disease\nHCC, anti-PD-1 + TKI", "post", d_hcc("H1")),
    ("Cross-disease\nHCC, anti-PD-1 + TKI", "pre", d_hcc("H2")),
]


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    y = np.arange(len(ROWS))[::-1]
    ax.axvline(0, lw=.9, color="#333", zorder=2)
    for k in (1.5, 3.5):
        ax.axhline(k, lw=.7, color="#E4E4E4", zorder=0)

    for k, (coh, arm, r) in enumerate(ROWS):
        yy = y[k]
        c = C_POST if arm == "post" else C_PRE
        sig = r["p"] < 0.05
        if np.isfinite(r["lo"]):
            ax.plot([r["lo"], r["hi"]], [yy, yy], lw=1.3, color=c, zorder=3)
            for e in (r["lo"], r["hi"]):
                ax.plot([e, e], [yy - .10, yy + .10], lw=1.3, color=c, zorder=3)
        ax.scatter([r["delta"]], [yy], s=78 if sig else 58,
                   marker="o" if np.isfinite(r["lo"]) else "D",
                   color=c if sig else "white", edgecolors=c,
                   linewidths=1.6, zorder=5)
        ax.text(-1.30, yy, f"{arm}-treatment", ha="right", va="center",
                fontsize=8, color=c)
        ax.text(2.55, yy, f"Δ = {r['delta']:+.2f}", ha="right", va="center",
                fontsize=8, color="#333")
        ptxt = f"P = {r['p']:.4f}".rstrip("0") if r["p"] < .01 else f"P = {r['p']:.3f}"
        ax.text(4.05, yy, ptxt, ha="right", va="center", fontsize=8,
                color=c if sig else "#888",
                fontweight="bold" if sig else "normal")
        ax.text(5.35, yy, f"{r['n_R']} R / {r['n_NR']} NR", ha="right",
                va="center", fontsize=7.2, color="#999")
        if np.isfinite(r["floor"]) and r["floor"] > 0.05:
            ax.text(5.50, yy, "◀ cannot reach P < 0.05", ha="left", va="center",
                    fontsize=7.2, color=C_FLOOR)

    for k, lab in ((0.5, ROWS[4][0]), (2.5, ROWS[2][0]), (4.5, ROWS[0][0])):
        ax.text(-3.45, k, lab, ha="left", va="center", fontsize=8.2,
                color="#333", linespacing=1.4)

    ax.set_yticks([])
    ax.set_xlim(-3.5, 8.2); ax.set_ylim(-.7, 5.7)
    ax.set_xticks([-1, 0, 1])
    ax.set_xlabel("Composite signature score, non-responders minus responders")
    ax.xaxis.set_label_coords(.21, -.075)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(-1.2, 1.4)
    ax.set_title("a  The same pre-specified test in three cohorts",
                 fontsize=9.5, loc="left", pad=8)


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    items = [(coh.split("\n")[0], arm, r) for coh, arm, r in ROWS
             if np.isfinite(r["floor"])]
    y = np.arange(len(items))[::-1]

    for k, (coh, arm, r) in enumerate(items):
        yy = y[k]
        blocked = r["floor"] > 0.05
        ax.barh(yy, r["floor"], height=.46,
                color=C_FLOOR if blocked else "#BBBBBB",
                alpha=.85 if blocked else .6, zorder=3)
        ax.text(r["floor"] + .003, yy, f"{r['floor']:.3f}", va="center",
                fontsize=7.6, color=C_FLOOR if blocked else "#777")
        ax.text(-.004, yy, f"{arm}-treatment   {r['n_R']} R / {r['n_NR']} NR",
                ha="right", va="center", fontsize=7.8,
                color=C_FLOOR if blocked else "#555")

    ax.axvline(.05, ls="--", lw=1, color="#333", zorder=4)
    ax.text(.052, len(items) - .42, "P = 0.05", fontsize=7.4, color="#333")
    ax.set_yticks([])
    ax.set_xlim(-.030, .098); ax.set_ylim(-2.35, len(items) - .28)
    ax.set_xlabel("Smallest one-sided permutation P the design can produce")
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(0, .09)
    ax.text(.985, .02,
            "The cross-disease pre-treatment arm has four responders against two\n"
            "non-responders, so no effect size whatever could have reached P < 0.05.\n"
            "This was computed from the sample structure and registered before any\n"
            "expression value was read (Supplementary S21) — it is the same criterion\n"
            "by which another dataset was excluded from the replication search.",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.9,
            color="#666", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#DDD", lw=.6))
    ax.set_title("b  What the cross-disease design could have shown at all",
                 fontsize=9.5, loc="left", pad=6)


def main():
    fig = plt.figure(figsize=(14.6, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.95, 1], wspace=.10,
                          left=.105, right=.985, top=.805, bottom=.105)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    fig.suptitle("Present in each cohort, confirmed in none twice",
                 fontsize=11.5, y=.962)
    fig.text(.5, .885,
             "Both arms significant on discovery · no arm confirmed a second time · "
             "the two follow-ups fail in mirror image",
             ha="center", va="bottom", fontsize=8.4, color="#444")
    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "Fig9_patients" + ext), **kw)
    plt.close(fig)

    for coh, arm, r in ROWS:
        print(f"Fig9  {coh.splitlines()[0]:<26} {arm:<5} "
              f"delta {r['delta']:+.4f}  P {r['p']:.4f}  "
              f"{r['n_R']}R/{r['n_NR']}NR"
              + (f"  floor {r['floor']:.4f}" if np.isfinite(r['floor']) else ""))


if __name__ == "__main__":
    main()
