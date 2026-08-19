"""Supplementary figure for S36 — the mismatched-list control across windows.

The main grid reports each negative control at one setting. Sweeping the setting
shows something the grid cannot: how much of "the control is clean" is a property
of the data and how much is a property of the window.

Every cell is drawn separately, because they do not fall into two tidy groups.
Three cells have literally no mismatch enrichment at any width. Melanoma on whole
blood does, and its fold trajectory is almost identical to RA on CD4+ T cells --
they part company on whether the enrichment reaches significance, not on its
size. Rheumatoid arthritis on whole blood fails everywhere and is already void.

Two panels, one per convention the paper calls "1 Mb":
  a  KNOWN_KB  how near a lead SNP a locus must be to count as known
  b  LOCUS_KB  the width of a locus (under fixed anchoring this IS the max span)

Read the y axis as "how much the WRONG disease's list enriches". 1x is no
enrichment; a filled marker is a setting at which that cell fails its own
pre-registered negative control.

Data: 128a_known_kb_sweep.tsv, 128b_locus_kb_sweep.tsv (step128).
Output: figures/FigS_window_control.pdf / .png
"""
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

# cell -> (label, colour, linewidth, linestyle, z, x-dodge)
# The dodge exists because two series nearly coincide at the tightest window and
# the marker underneath is the one that carries information.
STYLE = {
    "RA x Soskic_CD4":          ("RA × CD4⁺",            "#C4453C", 2.6, "-",  6, -0.055),
    "RA x eQTLGen_blood":       ("RA × blood (void)",    "#8A6BBE", 1.8, "--", 5, 0.0),
    "melanoma x eQTLGen_blood": ("Melanoma × blood",     "#3B7DD8", 1.8, "-",  4, 0.055),
    # These three sit at exactly 0.00 at every setting, so they plot as one
    # line. Drawing three overlapping lines and listing three legend entries
    # would promise a distinction the figure cannot show.
    "melanoma x Soskic_CD4":    ("Melanoma × CD4⁺, HCC × CD4⁺,\nHCC × blood — 0.00 throughout",
                                 "#5FA37A", 1.6, "-",  3, 0.0),
}
HIDDEN_AT_FLOOR = ("HCC_high x Soskic_CD4", "HCC_high x eQTLGen_blood")
FLOOR = 0.028          # where a 0.00 fold is drawn on a log axis

kn = pd.read_csv(f"{MR}/128a_known_kb_sweep.tsv", sep="\t")
lk = pd.read_csv(f"{MR}/128b_locus_kb_sweep.tsv", sep="\t")

# Collapsing three cells into one line is only honest while they agree.
# Check it rather than assume it; if a rerun ever separates them, stop.
for _t in (kn, lk):
    _f = _t[_t.cell.isin(("melanoma x Soskic_CD4",) + HIDDEN_AT_FLOOR)].mismatch_fold
    assert (_f == 0).all(), (
        "the three floor cells are no longer all 0.00; give them separate lines")


def panel(ax, tab, key, xlabel, title):
    xs = sorted(tab[key].unique())
    xi = {v: i for i, v in enumerate(xs)}
    ax.axhline(1, ls=":", lw=1, color="#999", zorder=1)
    ax.axvline(xi[1000], color="#333", lw=.9, alpha=.3, zorder=1)

    for cell, (lab, col, lw, ls, z, dx) in STYLE.items():
        s = tab[tab.cell == cell].sort_values(key)
        if s.empty:
            continue
        x = [xi[v] + dx for v in s[key]]
        y = [v if v > 0 else FLOOR for v in s.mismatch_fold]
        ax.plot(x, y, color=col, lw=lw, ls=ls, zorder=z, alpha=.95,
                solid_capstyle="round")
        for xx, yy, failed in zip(x, y, s.control == "FAILED"):
            ax.plot(xx, yy, "o", ms=8 if failed else 5,
                    mfc=col if failed else "white", mec=col,
                    mew=1.9 if failed else 1.4, zorder=z + 10)

    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels([f"{v:,}" for v in xs], fontsize=8.4)
    ax.set_xlim(-.4, len(xs) - .5)
    ax.set_yscale("log")
    ax.set_ylim(.02, 7)
    ax.set_yticks([FLOOR, 0.5, 1, 2, 4])
    ax.set_yticklabels(["0", "0.5", "1", "2", "4"], fontsize=8.4)
    ax.set_xlabel(xlabel, fontsize=8.8)
    ax.set_title(title, fontsize=9.4, loc="left", pad=18)
    ax.text(xi[1000] - .06, 5.6, "registered\nmain analysis", fontsize=6.9,
            color="#555", ha="right", va="top", linespacing=1.4)


def main():
    fig = plt.figure(figsize=(11.6, 4.9))
    gs = fig.add_gridspec(1, 2, wspace=.14, left=.085, right=.615,
                          top=.775, bottom=.135)
    ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    panel(ax1, kn, "known_kb", "Known-locus window (kb)",
          "a  How near a lead SNP a locus must be")
    panel(ax2, lk, "locus_kb", "Locus width (kb)",
          "b  How wide a locus may be")
    ax1.set_ylabel("Enrichment of the WRONG disease's known loci\n"
                   "(mismatched-list negative control)",
                   fontsize=8.5, linespacing=1.5)
    ax2.tick_params(labelleft=False)

    handles = [Line2D([], [], color=c, lw=lw, ls=ls, marker="o", ms=5,
                      mfc="white", mec=c, label=lab)
               for lab, c, lw, ls, _, _dx in STYLE.values()]
    handles.append(Line2D([], [], color="#444", lw=0, marker="o", ms=8,
                          mfc="#444", mec="#444",
                          label="filled = fails its control here"))
    fig.legend(handles=handles, frameon=False, fontsize=7.7,
               loc="upper left", bbox_to_anchor=(.632, .80),
               handlelength=2.4, labelspacing=.78)

    fig.text(.632, .445,
             "Three cells show no mismatch enrichment at any width. The\n"
             "other three rise as the window narrows, and two of them —\n"
             "melanoma on whole blood and RA on CD4⁺ T cells — trace\n"
             "almost the same curve. They differ in whether that\n"
             "enrichment reaches significance, not in its size.\n\n"
             "RA × CD4⁺ is clean at exactly one setting: the registered\n"
             "one, at the coarse end. A wide window is where the control\n"
             "discriminates least — both lists hit a great deal and the\n"
             "ratio is diluted towards 1 — so the cell clears its negative\n"
             "control where that control is weakest. Hence the non-cancer\n"
             "result is reported as descriptive.",
             fontsize=7.05, color="#444", va="top", linespacing=1.55)

    fig.suptitle("A negative control that depends on the window: rheumatoid "
                 "arthritis versus cancer",
                 fontsize=11.2, y=.955, x=.085, ha="left")

    for ext, kw in ((".pdf", {}), (".png", {"dpi": 300})):
        fig.savefig(os.path.join(OUT, "FigS_window_control" + ext), **kw)
    plt.close(fig)

    both = pd.concat([kn, lk])
    for cell in STYLE:
        s = both[both.cell == cell]
        print(f"  {cell:<26} fails {int((s.control == 'FAILED').sum())}/8"
              f"   fold {s.mismatch_fold.min():.2f}-{s.mismatch_fold.max():.2f}")


if __name__ == "__main__":
    main()
