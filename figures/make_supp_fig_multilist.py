"""Supplementary figure for S39 — every cell against every reference list.

The table this accompanies is read row by row, and row by row it says the
outcome's own list wins. What it cannot show is by how much, and the margin is
the whole point: the six cancer cells span 1.82- to 7.63-fold over their best
comparator, and the two rheumatoid arthritis cells sit at 1.28 and 1.32.

One row per cell. The filled diamond is the outcome's own list; open circles are
the four unrelated lists, filled where that list also enriches at P < 0.05. The
bracket spans own fold to the largest comparator fold, so the margin is what the
eye measures.

An earlier version of this figure drew that bracket to the largest *significant*
comparator, which made it infinite for cells where no comparator reached
P < 0.05, and the caption called that "an unbounded margin". That was wrong: a
comparator that misses significance still has a point estimate -- prostate is
2.73-fold on melanoma x CD4. The bracket now runs to the largest comparator fold
whether or not it is significant, and every margin is finite.

Data: 130a_multilist_main.tsv, 130c_multilist_verdict.tsv (step130).
Output: figures/FigS_multilist_control.pdf / .png
"""
import math
import os

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

C_OWN, C_CTL, C_HIT, C_RA = "#C4453C", "#B8B8B8", "#3B7DD8", "#C4453C"

main = pd.read_csv(f"{MR}/130a_multilist_main.tsv", sep="\t")
vd = pd.read_csv(f"{MR}/130c_multilist_verdict.tsv", sep="\t")

# cancer cells first, RA last, so the margin collapse reads top to bottom
ORDER = ["melanoma x Soskic_CD4", "melanoma x eQTLGen_blood",
         "HCC_high x Soskic_CD4", "HCC_high x eQTLGen_blood",
         "HCC_low x Soskic_CD4", "HCC_low x eQTLGen_blood",
         "RA x Soskic_CD4", "RA x eQTLGen_blood"]
LAB = {"melanoma x Soskic_CD4": "Melanoma × CD4⁺",
       "melanoma x eQTLGen_blood": "Melanoma × blood",
       "HCC_high x Soskic_CD4": "HCC × CD4⁺",
       "HCC_high x eQTLGen_blood": "HCC × blood",
       "HCC_low x Soskic_CD4": "HCC-low × CD4⁺",
       "HCC_low x eQTLGen_blood": "HCC-low × blood",
       "RA x Soskic_CD4": "RA × CD4⁺",
       "RA x eQTLGen_blood": "RA × blood  (void)"}
FLOOR = 0.055     # where a 0.00 fold is drawn on a log axis


def main_fig():
    # Title, caption, plot and legend each get their own band. Packing them by
    # eye put the legend on top of the bottom row and the caption through the
    # title, which is how the first draft of this figure read.
    fig, ax = plt.subplots(figsize=(11.2, 6.6))
    fig.subplots_adjust(left=.145, right=.985, top=.775, bottom=.185)

    for i, cell in enumerate(ORDER):
        y = len(ORDER) - 1 - i
        s = main[main.cell == cell]
        v = vd[vd.cell == cell].iloc[0]
        own = s[s.relation == "own outcome"].iloc[0]
        ctl = s[s.relation == "unrelated control"]

        ax.axhline(y, color="#EEE", lw=8, zorder=0)

        # The margin the argument rests on, measured to the largest comparator
        # whether or not it reached significance. Restricting it to significant
        # comparators is what produced the "unbounded margin" this figure used
        # to claim.
        fmax = float(v.F_max_all)
        if fmax > 0:
            ax.plot([fmax, v.F_own], [y, y], color="#666", lw=1.6, zorder=2)
            ax.text((fmax * v.F_own) ** .5, y + .26,
                    f"×{v.F_own / fmax:.2f}", ha="center", va="bottom",
                    fontsize=7.6, color="#222", fontweight="bold",
                    bbox=dict(boxstyle="square,pad=0.12", fc="white",
                              ec="none"))

        placed = []          # x positions of labels already drawn on this row
        for _, r in ctl.sort_values("fold").iterrows():
            enr = pd.notna(r.fisher_p) and r.fisher_p < 0.05 and r.fold > 1
            x = r.fold if r.fold > 0 else FLOOR
            ax.plot(x, y, "o", ms=8 if enr else 6,
                    mfc=C_HIT if enr else "white",
                    mec=C_HIT if enr else C_CTL, mew=1.6, zorder=4)
            if not enr:
                continue
            # drop to a second line when the previous label is within ~15% on
            # the log axis, rather than letting two names overprint
            crowded = any(abs(math.log10(x) - math.log10(px)) < .15
                          for px in placed)
            ax.text(x, y - (.52 if crowded else .28), r.list_name,
                    ha="center", va="top", fontsize=6.5, color=C_HIT)
            placed.append(x)

        ax.plot(own.fold, y, "D", ms=9, mfc=C_OWN, mec="white", mew=1.1,
                zorder=6)

    ax.axvline(1, ls=":", lw=1, color="#999", zorder=1)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels([LAB[c] for c in reversed(ORDER)], fontsize=8.6)
    ax.set_xscale("log")
    ax.set_xlim(.045, 30)
    ax.set_xticks([FLOOR, 1, 2, 5, 10, 20])
    ax.set_xticklabels(["0", "1", "2", "5", "10", "20"], fontsize=8.4)
    ax.set_ylim(-.95, len(ORDER) - .25)
    ax.set_xlabel("Fold enrichment on a disease's known loci "
                  "(log scale; 1× = no enrichment)", fontsize=8.8)
    ax.tick_params(axis="y", length=0)

    handles = [
        Line2D([], [], lw=0, marker="D", ms=8, mfc=C_OWN, mec="white",
               label="the outcome's OWN known loci"),
        Line2D([], [], lw=0, marker="o", ms=6, mfc="white", mec=C_CTL, mew=1.6,
               label="an unrelated disease's list (n.s.)"),
        Line2D([], [], lw=0, marker="o", ms=8, mfc=C_HIT, mec=C_HIT,
               label="an unrelated list that also enriches (P < 0.05)"),
        Line2D([], [], color="#666", lw=1.6,
               label="margin: own fold ÷ largest comparator fold"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7.7, ncol=2,
              loc="upper center", bbox_to_anchor=(.5, -.115),
              handletextpad=.7, columnspacing=2.6)

    fig.suptitle("Every cell beats its comparators, but rheumatoid arthritis beats "
                 "them by the\nnarrowest margin in the grid",
                 fontsize=11.4, x=.145, ha="left", y=.985, va="top",
                 linespacing=1.35)
    fig.text(.145, .885,
             "Each row is one cell, scored against every reference list in the "
             "study with \u226530 placeable lead SNPs. A post-hoc diagnostic (S39), "
             "not a registered test, and\nnot used to upgrade any cell. Cancer "
             "cells span 1.82\u2013 to 7.63-fold over their largest comparator; the "
             "two RA cells sit at 1.28 and 1.32. Of the four comparators\n"
             "available for RA, two leave the cell standing \u2014 prostate (1.04\u00d7, "
             "P = 0.59) and melanoma (1.71\u00d7, P = 0.10), the one this study "
             "designated \u2014 while HCC (3.05\u00d7,\nP = 0.0028) and colorectal "
             "(1.77\u00d7, P = 0.0045) would each have voided it. Whether the cell "
             "survives depends on which comparator was named.",
             fontsize=7.35, color="#444", va="top", linespacing=1.65)

    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "FigS_multilist_control" + ext), **kw)
    plt.close(fig)

    for c in ORDER:
        v = vd[vd.cell == c].iloc[0]
        print(f"  {c:<26} k={v.k_enriching}  own {v.F_own:5.2f}  "
              f"largest comparator {v.F_max_all:5.2f} ({v.best_rival_any})  "
              f"margin {v.F_own / v.F_max_all:.2f}")


if __name__ == "__main__":
    main_fig()
