"""Fig — power-stability curves (Step 53-55).

Three panels:
  (a) melanoma recovery curve + EMPIRICAL calibration at FinnGen power.
      NB not 'external' and not 'independent': FinnGen contributes to the meta
      outcome, so this is a separately observed round, not an independent one.
  (b) five cancers, recovery vs fraction of observed cases
  (c) *** recovery stratified by locus class: known pigmentation/naevus vs novel

预设的表述纪律（跑之前写死，勿事后改）:
  - 曲线只画到观测功效为止，不外推（两次外推尝试均失败，见 FINDINGS Step 53）
  - (c) 的百分比基数很小（已知 6 基因 / 新位点 4 基因），画 CI 带，不报小数点后一位
"""
import csv, os
import numpy as np
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

C_KNOWN, C_NOVEL = "#C4453C", "#3B7DD8"      # 与 FigA 的位点类别配色一致
C_HITS, C_JAC = "#333333", "#2E7D5B"
TRAIT_COL = {"Melanoma": "#C4453C", "Lung": "#3B7DD8", "Colorectal": "#E8A33D",
             "Breast": "#7A5AA8", "Prostate": "#2E7D5B"}

FINNGEN_CASES = 5753          # 独立测得的校准点
REAL_HITS, REAL_JAC = 10, 0.4545


def rd(fp):
    with open(os.path.join(MR, fp), encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def col(rows, k, cast=float):
    return np.array([cast(r[k]) for r in rows])


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300)
    plt.close(fig)
    print(name, "ok", flush=True)


# ------------------------------------------------------------------ panel a
def panel_a(ax):
    r = rd("53a_power_stability_curve.tsv")
    cases, hits = col(r, "cases"), col(r, "mean_hits")
    jac, jlo, jhi = col(r, "jaccard_genes"), col(r, "j_lo"), col(r, "j_hi")
    cal = rd("53b_calibration.tsv")[0]

    ax.plot(cases / 1000, hits, "-o", ms=3.5, lw=1.4, color=C_HITS, zorder=4,
            label="FDR<0.05 discoveries")
    ax.set_xlabel("Melanoma cases (thousands)")
    ax.set_ylabel("Mean discoveries at FDR < 0.05", color=C_HITS)
    ax.set_ylim(0, 24)

    ax2 = ax.twinx()
    ax2.spines["top"].set_visible(False)
    ax2.fill_between(cases / 1000, jlo, jhi, color=C_JAC, alpha=.14, lw=0, zorder=1)
    ax2.plot(cases / 1000, jac, "-", lw=1.4, color=C_JAC, zorder=3,
             label="Jaccard vs full-power gene list")
    ax2.set_ylabel("Gene-list Jaccard vs full-power list", color=C_JAC)
    ax2.set_ylim(0, 1.05)
    ax2.tick_params(axis="y", colors=C_JAC)
    ax.tick_params(axis="y", colors=C_HITS)

    # 外部校准：把 meta 降到 FinnGen 功效，必须重现独立测得的真实值
    # 两个真实值在两套坐标下几乎重叠，故左右各偏移 0.2k 病例以便区分
    ax.axvline(FINNGEN_CASES / 1000, ls=":", lw=1, color="#888", zorder=2)
    ax.scatter([FINNGEN_CASES / 1000 - .22], [REAL_HITS], marker="*", s=170,
               facecolor="white", edgecolor=C_HITS, lw=1.2, zorder=6)
    ax2.scatter([FINNGEN_CASES / 1000 + .22], [REAL_JAC], marker="*", s=170,
                facecolor="white", edgecolor=C_JAC, lw=1.2, zorder=6)
    ax.annotate("★ = separately observed FinnGen round\n     (a component of the meta outcome,\n      not an independent study)",
                (FINNGEN_CASES / 1000, 8.4), textcoords="offset points",
                xytext=(10, -34), fontsize=7, color="#444",
                arrowprops=dict(arrowstyle="-", lw=.5, color="#AAA"))
    ax.text(.03, .97,
            "Calibration at 5,753 cases\n"
            f"hits  {float(cal['sim_hits']):.1f} sim  vs  {int(float(cal['real_hits']))} real\n"
            f"Jaccard  {float(cal['sim_jaccard']):.2f} "
            f"[{float(cal['sim_lo']):.2f}, {float(cal['sim_hi']):.2f}]  vs  "
            f"{float(cal['real_jaccard']):.2f} real",
            transform=ax.transAxes, va="top", fontsize=7.2, color="#444",
            bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#DDD", lw=.6))
    ax.set_title("a  Down-sampling reproduces a separately observed round",
                 fontsize=9.5, loc="left", pad=6)


# ------------------------------------------------------------------ panel b
def panel_b(ax):
    rows = rd("54a_multitrait_power_curves.tsv")
    thr = {r["trait"]: r for r in rd("54b_power_thresholds.tsv") if r["target"] == "0.5"}
    for t, c in TRAIT_COL.items():
        d = [r for r in rows if r["trait"] == t]
        if not d:
            continue
        x, y = col(d, "frac"), col(d, "sensitivity")
        n = int(float(thr[t]["obs_cases"])) if t in thr else 0
        ax.plot(x, y, "-o", ms=3, lw=1.3, color=c,
                label=f"{t} ({n:,} cases)")
    ax.axhline(.5, ls="--", lw=.8, color="#666")
    ax.text(1.0, .515, "50% of the full-power list recovered",
            ha="right", fontsize=7, color="#555")
    ax.set_xlim(0.05, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Fraction of the study's own observed cases")
    ax.set_ylabel("Recovery of the full-power gene list")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.set_title("b  The same curve in five cancers", fontsize=9.5, loc="left", pad=6)


# ------------------------------------------------------------------ panel c
def panel_c(ax):
    r = rd("55a_recovery_by_locus_class.tsv")
    cases = col(r, "cases") / 1000
    kn, klo, khi = col(r, "known_recovery"), col(r, "known_lo"), col(r, "known_hi")
    nv, nlo, nhi = col(r, "novel_recovery"), col(r, "novel_lo"), col(r, "novel_hi")

    # 同类研究结局常见规模
    ax.axvspan(3, 8, color="#EDEDED", zorder=0)
    ax.text(7.9, .04, "case numbers typical of\nstudies using this framework",
            ha="right", va="bottom", fontsize=7.5, color="#777", zorder=1)

    for y, lo, hi, c, lab in [(kn, klo, khi, C_KNOWN, "Known pigmentation / naevus loci (6 genes)"),
                              (nv, nlo, nhi, C_NOVEL, "Novel loci (4 genes)")]:
        ax.fill_between(cases, lo, hi, color=c, alpha=.10, lw=0, zorder=2)
        ax.plot(cases, lo, lw=.6, ls="-", color=c, alpha=.45, zorder=2)
        ax.plot(cases, hi, lw=.6, ls="-", color=c, alpha=.45, zorder=2)
        ax.plot(cases, y, "-o", ms=4, lw=1.9, color=c, label=lab, zorder=4)

    ax.axhline(.5, ls="--", lw=.8, color="#666", zorder=3)
    # 达到 50% 恢复所需病例数：已知 2,506 / 新位点 8,771 → 3.5 倍
    for x, c in [(2.506, C_KNOWN), (8.771, C_NOVEL)]:
        ax.plot([x, x], [0, .5], ls=":", lw=1.1, color=c, zorder=3)
    ax.annotate("", xy=(2.506, .5), xytext=(8.771, .5),
                arrowprops=dict(arrowstyle="<->", lw=1, color="#333"))
    ax.text(5.64, .525, "3.5× more cases to recover half the list",
            ha="center", fontsize=8, color="#222", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"))
    ax.text(2.506, -.06, "2,506", ha="center", fontsize=7.5, color=C_KNOWN)
    ax.text(8.771, -.06, "8,771", ha="center", fontsize=7.5, color=C_NOVEL)

    ax.annotate("46% vs 0.4%\nat 10% of observed power", (1.253, .23),
                textcoords="offset points", xytext=(16, -18), fontsize=7.5,
                color="#444", zorder=6,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"),
                arrowprops=dict(arrowstyle="-", lw=.5, color="#AAA"))

    ax.set_xlim(0.8, 11.8)
    ax.set_ylim(0, 1.04)
    ax.set_xlabel("Melanoma cases (thousands)")
    ax.set_ylabel("Recovery of the full-power gene list")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_title("c  The reproducible part of the candidate list is the part that is not a discovery",
                 fontsize=9.5, loc="left", pad=6)


def main():
    fig = plt.figure(figsize=(11, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=.42, wspace=.32,
                          left=.07, right=.93, top=.93, bottom=.09)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, :]))
    fig.suptitle("Outcome GWAS power sets the reproducibility of the candidate list",
                 fontsize=11, y=.975)
    save(fig, "FigF_power_stability")


if __name__ == "__main__":
    main()
