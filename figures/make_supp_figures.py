"""Supplementary figures S2–S8.

每张图对应正文里一个"我们自己查出来的问题"或一个必须交代的对照，
故补充材料不是堆砌，而是把正文中被压缩的证据展开。

S2  Steiger：TwoSampleMR 的 R² 越界与有界形式
S3  四种估计量对比（含 weighted mode 为何不适用）
S4  痣阴性对照与共病 MR
S5  谱系纯度对照：效应量随对照层级塌缩
S6  两处数据处理错误：修正前后
S7  糖酵解轴：家族构成 + motif（全部 peak vs 活化不变 peak）
S8  外周血 ICC：随机模块地板与上界对照
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
RED, BLUE, GREEN, GREY = "#C4453C", "#3B7DD8", "#2E7D5B", "#8B9096"


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300)
    plt.close(fig)
    print(name, "ok", flush=True)


# ================================================================== S2
def s2():
    d = pd.read_csv(f"{MR}/09_steiger_filtering.tsv", sep="\t")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    a = ax[0]
    a.hist(d["rsq.exposure"].dropna(), bins=60, color=GREY, zorder=3)
    a.axvline(1, ls="--", lw=1.1, color=RED, zorder=4)
    n_over = int((d["rsq.exposure"] > 1).sum()); n_08 = int((d["rsq.exposure"] > .8).sum())
    a.text(.97, .93, f"values > 1: {n_over}\nvalues > 0.8: {n_08}\nmaximum {d['rsq.exposure'].max():.3f}",
           transform=a.transAxes, ha="right", va="top", fontsize=8, color=RED)
    a.set_xlabel(r"$R^2_{\rm exposure}$, TwoSampleMR formula for SD units")
    a.set_ylabel("Records")
    a.set_title("a  The published formula is unbounded", fontsize=9.5, loc="left", pad=6)

    b = ax[1]
    v = d["rsq.exposure.bounded"].dropna()
    b.hist(v, bins=60, color=GREEN, zorder=3)
    b.text(.97, .93, f"range {v.min():.3f}–{v.max():.3f}\nmedian {v.median():.3f}\nno value exceeds 1",
           transform=b.transAxes, ha="right", va="top", fontsize=8, color=GREEN)
    b.set_xlabel(r"$R^2 = F/(F+N-2)$   (bounded form, reported)")
    b.set_title("b  Bounded form; direction unchanged", fontsize=9.5, loc="left", pad=6)
    fig.suptitle("Fig S2 | Steiger filtering is reported as procedure, not evidence",
                 fontsize=10.5, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, .94]); save(fig, "FigS2_steiger")


# ================================================================== S3
def s3():
    d = pd.read_csv(f"{MR}/18_sensitivity_meta.tsv", sep="\t")
    d = d[d.method != "method"].copy()
    d["pval"] = pd.to_numeric(d.pval, errors="coerce")
    d["FDR"] = pd.to_numeric(d.FDR, errors="coerce")
    d["nsnp"] = pd.to_numeric(d.nsnp, errors="coerce")
    meth = ["Inverse variance weighted",
            "Inverse variance weighted (multiplicative random effects)",
            "Weighted median"]
    lab = ["IVW", "IVW-MRE", "Weighted\nmedian", "Weighted\nmode"]
    n_p = [int((d[d.method == m].pval < .05).sum()) for m in meth] + [1]
    n_f = [int((d[d.method == m].FDR < .05).sum()) for m in meth] + [0]
    n_t = [int((d.method == m).sum()) for m in meth] + [226]

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    x = np.arange(4); w = .38
    ax[0].bar(x - w/2, n_p, w, color=GREY, label="nominal P < 0.05", zorder=3)
    ax[0].bar(x + w/2, n_f, w, color=RED, label="FDR < 0.05", zorder=3)
    for i, (p, f, t) in enumerate(zip(n_p, n_f, n_t)):
        ax[0].text(i - w/2, p + 3, str(p), ha="center", fontsize=7.5)
        ax[0].text(i + w/2, f + 3, str(f), ha="center", fontsize=7.5)
        ax[0].text(i, -14, f"of {t}", ha="center", fontsize=7, color="#777")
    ax[0].set_xticks(x); ax[0].set_xticklabels(lab, fontsize=8)
    ax[0].set_ylabel("Exposures"); ax[0].set_ylim(-22, max(n_p) * 1.2)
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].annotate("1 of 226 significant\n→ no discriminating power,\nreported as inapplicable",
                   (3, 8), xytext=(2.35, 95), fontsize=7.5, color=RED,
                   arrowprops=dict(arrowstyle="->", lw=.9, color=RED))
    ax[0].set_title("a  Weighted mode is not usable at this instrument count",
                    fontsize=9.5, loc="left", pad=6)

    nn = d[d.method == "Inverse variance weighted"].nsnp.dropna()
    ax[1].hist(nn, bins=np.arange(1, nn.max() + 2) - .5, color=BLUE, zorder=3)
    ax[1].axvline(nn.median(), ls="--", lw=1.1, color=RED)
    ax[1].text(nn.median() + .4, ax[1].get_ylim()[1] * .9,
               f"median {nn.median():.0f} instruments", fontsize=8, color=RED)
    ax[1].set_xlabel("Independent instruments per exposure (r² < 0.1)")
    ax[1].set_ylabel("Exposures"); ax[1].set_xlim(0, 20)
    ax[1].set_title("b  Why: the instrument count", fontsize=9.5, loc="left", pad=6)
    fig.suptitle("Fig S3 | Multi-instrument sensitivity analysis, and the estimator we excluded",
                 fontsize=10.5, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, .93]); save(fig, "FigS3_estimators")


# ================================================================== S4
def s4():
    d = pd.read_csv(f"{MR}/20_comorbidity_MR.tsv", sep="\t")
    nv = d[d.pheno == "CD2_BENIGN_MELANOCYTIC"].copy()
    nv["pval"] = pd.to_numeric(nv.pval, errors="coerce")
    nv = nv.dropna(subset=["OR", "pval"]).sort_values(["group", "OR"])
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    y = np.arange(len(nv))
    cols = [RED if g == "阳性对照" else BLUE for g in nv.group]
    ax.errorbar(nv.OR, y, xerr=[nv.OR - nv.OR_L, nv.OR_U - nv.OR], fmt="o",
                ms=5, lw=1.2, capsize=3, ecolor="#BBB", zorder=3,
                mfc="none", mec="none")
    ax.scatter(nv.OR, y, s=42, c=cols, zorder=4, edgecolors="white", linewidths=.6)
    ax.axvline(1, ls="--", lw=1, color="#666")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g}   P={p:.2g}" for g, p in zip(nv.gene, nv.pval)], fontsize=8)
    ax.set_xlabel("Odds ratio for melanocytic naevi (95% CI)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=7, color=RED,
                              label="Known pigmentation locus (positive control)"),
                       Line2D([], [], marker="o", ls="", ms=7, color=BLUE,
                              label="Candidate")],
              frameon=False, fontsize=8, loc="lower right")
    ax.set_title("Fig S4 | Naevi as a pigmentation-pathway negative control\n"
                 "Test sensitivity is established by the positive controls, which do affect naevi",
                 fontsize=9.5, loc="left", pad=8)
    fig.tight_layout(); save(fig, "FigS4_naevus_control")


# ================================================================== S5
def s5():
    a = pd.read_csv(f"{MR}/29a_purity_control_module_tests.tsv", sep="\t")
    g = pd.read_csv(f"{MR}/29g_CD8ness_matched_modules.tsv", sep="\t")
    mods = ["CytotoxicScore", "HelperRegScore", "EffectorClean", "HelperClean"]
    levels = [("No control\n(depth only)", "all matched pairs", a),
              ("CD8-negative\npairs", "CD8-free pairs", a),
              ("CD8ness + depth\nmatched", None, g)]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    x = np.arange(len(mods)); w = .26
    for k, (lab, subset, src) in enumerate(levels):
        vals, ps = [], []
        for m in mods:
            if subset is None:
                r = src[src.module == m]
            else:
                r = src[(src.module == m) & (src.set == subset)]
            vals.append(float(r.mean_delta.iloc[0]) if len(r) else np.nan)
            ps.append(float(r.p.iloc[0]) if len(r) else np.nan)
        cols = [[GREY, "#D9A7A4", "#F0D9D8"][k] if p >= .05 else [RED, "#D9736C", "#E8A8A3"][k]
                for p in ps]
        ax.bar(x + (k - 1) * w, vals, w, color=cols, zorder=3,
               label=lab, edgecolor="white", linewidth=.6)
        for i, (v, p) in enumerate(zip(vals, ps)):
            ax.text(i + (k - 1) * w, v + (.012 if v >= 0 else -.028),
                    "*" if p < .05 else "ns", ha="center", fontsize=7.5,
                    color="#333" if p < .05 else "#999")
    ax.axhline(0, lw=.9, color="#444")
    ax.set_xticks(x); ax.set_xticklabels(
        ["Cytotoxic\n(contains CD8A/B)", "Helper/Reg\n(contains CD4)",
         "Effector\n(lineage-free)", "Helper\n(lineage-free)"], fontsize=8)
    ax.set_ylabel("Mean Δ, glycolysis-high vs -low")
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center",
              bbox_to_anchor=(.5, 1.02))
    ax.set_ylim(-.45, .40)
    ax.set_title("Fig S5 | A score-based split also splits lineage purity\n"
                 "Effects collapse 75–100% once cells are matched on lineage as well as depth",
                 fontsize=9.5, loc="left", pad=26)
    fig.tight_layout(); save(fig, "FigS5_purity_control")


# ================================================================== S6
def s6():
    rows = [("TPI1, pre-treatment", 0.0128, 0.0076), ("TPI1, post-treatment", 0.0019, 0.0057),
            ("Proliferation, pre", 0.0101, 0.0101), ("Exhaustion, pre", 0.0128, 0.0172),
            ("SMC2, post", 0.029, 0.0047), ("KIAA0040, post", 0.017, 0.070),
            ("Glycolysis module, post", 0.0094, 0.00169)]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    a = ax[0]
    y = np.arange(len(rows))
    for i, (lab, before, after) in enumerate(rows):
        a.plot([before, after], [i, i], color="#CCC", lw=1.2, zorder=1)
        a.scatter(before, i, s=44, color=GREY, zorder=3, edgecolors="white", linewidths=.6)
        a.scatter(after, i, s=44, color=RED if after < .05 else "#D8C7C6",
                  zorder=4, edgecolors="white", linewidths=.6)
    a.axvline(.05, ls="--", lw=1, color="#666")
    a.set_xscale("log")
    a.set_xticks([1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
    a.set_xticklabels(["0.001", "0.003", "0.01", "0.03", "0.1"], fontsize=8)
    a.set_xlim(8e-4, 1.2e-1)
    a.set_yticks(y)
    a.set_yticklabels([r[0] for r in rows], fontsize=8)
    a.set_xlabel("P value (log scale)")
    a.legend(handles=[Line2D([], [], marker="o", ls="", ms=7, color=GREY, label="as first analysed"),
                      Line2D([], [], marker="o", ls="", ms=7, color=RED, label="after correction")],
             frameon=False, fontsize=8, loc="lower right")
    a.set_title("a  Every patient-level number changed", fontsize=9.5, loc="left", pad=6)

    b = ax[1]
    b.axis("off")
    b.text(0, 1.06, "b  Error 1 — re-normalisation", fontsize=9.5, weight="bold",
           transform=b.transAxes)
    b.text(0, .97,
           "log₂(TPM+1) values passed through NormalizeData().\n"
           "Diagnostic: back-solved per-cell sums.\n"
           "   correct object   9.94×10⁵  ≈ 10⁶  ✓\n"
           "   affected object  no clean solution  ✗\n"
           "Spearman between versions: 0.874 (not 1.0),\n"
           "so it is not a monotone rescaling.\n"
           "Cell-level rank tests unaffected; patient-level means are not.",
           fontsize=8, va="top", transform=b.transAxes, family="DejaVu Sans")
    b.text(0, .47, "   Error 2 — merged response labels", fontsize=9.5, weight="bold",
           transform=b.transAxes)
    b.text(0, .39,
           "Patients with two post-treatment biopsies can carry\n"
           "different outcomes:\n"
           "   Post_P1   responder        Post_P1_2  non-responder\n"
           "   Post_P5   non-responder    Post_P5_2  responder\n"
           "These had been collapsed into one label per patient.\n"
           "Corrected samples: 9R/10NR pre, 8R/20NR post.",
           fontsize=8, va="top", transform=b.transAxes)
    b.text(0, .02, "The criterion for error 1 is whether the SOURCE file is already\n"
                    "normalised — not whether a derived object matches it. Applying the\n"
                    "same test to the second cohort gives the opposite answer.",
           fontsize=7.6, va="top", transform=b.transAxes, color=RED)
    fig.suptitle("Fig S6 | Two processing errors we found in our own completed analyses",
                 fontsize=10.5, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, .93]); save(fig, "FigS6_data_errors")


# ================================================================== S7
def s7():
    fam = pd.read_csv(f"{MR}/44b_family_composition.tsv", sep="\t")
    mot = pd.read_csv(f"{MR}/50b_motif_comparison_all_vs_invariant.tsv", sep="\t")
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    a = ax[0]
    f = fam.sort_values("enrichment", ascending=True).tail(10)
    cols = [RED if d.lower().startswith("up") else BLUE for d in f.direction]
    a.barh(np.arange(len(f)), f.enrichment, color=cols, zorder=3)
    a.set_yticks(np.arange(len(f)))
    a.set_yticklabels([f"{fm}  ({p:.1f}%)" for fm, p in zip(f.family, f.pct_in_top)], fontsize=8)
    a.axvline(1, ls="--", lw=1, color="#666")
    a.set_xlabel("Enrichment in the top 300 axis genes (fold over background)")
    a.legend(handles=[Line2D([], [], marker="s", ls="", ms=8, color=RED, label="glycolysis-high end"),
                      Line2D([], [], marker="s", ls="", ms=8, color=BLUE, label="glycolysis-low end")],
             frameon=False, fontsize=8, loc="lower right")
    a.set_title("a  What the axis is made of (no pre-specified modules)",
                fontsize=9.5, loc="left", pad=6)

    b = ax[1]
    m = mot.dropna(subset=["odds_all", "odds_inv"]).copy()
    ap1 = m.motif.str.contains("JUN|FOS|BATF|JDP2", case=False, na=False)
    irf = m.motif.str.contains("IRF|STAT", case=False, na=False)
    b.scatter(m.odds_all[~(ap1 | irf)], m.odds_inv[~(ap1 | irf)], s=14,
              color="#D5D8DC", zorder=2, label="other motifs")
    b.scatter(m.odds_all[ap1], m.odds_inv[ap1], s=34, color=RED, zorder=4,
              edgecolors="white", linewidths=.4, label="AP-1 family")
    b.scatter(m.odds_all[irf], m.odds_inv[irf], s=34, color=BLUE, zorder=4,
              edgecolors="white", linewidths=.4, label="IRF / STAT family")
    lim = [m[["odds_all", "odds_inv"]].min().min() * .95,
           m[["odds_all", "odds_inv"]].max().max() * 1.05]
    b.plot(lim, lim, ls=":", lw=1, color="#888", zorder=1)
    b.axhline(1, lw=.7, color="#BBB"); b.axvline(1, lw=.7, color="#BBB")
    b.set_xlabel("Odds ratio, all peaks")
    b.set_ylabel("Odds ratio, activation-invariant peaks")
    b.text(.97, .05, "below the diagonal =\nweakened once activation\nis held constant",
           transform=b.transAxes, ha="right", fontsize=7.2, color="#666")
    b.legend(frameon=False, fontsize=8, loc="upper left")
    b.set_title("b  Motif enrichment before and after restricting to\n"
                "activation-invariant peaks", fontsize=9.5, loc="left", pad=6)
    fig.suptitle("Fig S7 | The metabolic axis and its chromatin signature",
                 fontsize=10.5, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, .92]); save(fig, "FigS7_axis_motif")


# ================================================================== S8
def s8():
    d = pd.read_csv(f"{MR}/52a_icc_all.tsv", sep="\t")
    rnd = d[d.kind == "random_floor"].icc.values
    fig, ax = plt.subplots(figsize=(8.4, 4.3))
    ax.hist(rnd, bins=22, color=GREY, alpha=.85, zorder=3, label="random gene modules (10 donors)")
    p95 = np.percentile(rnd, 95)
    ax.axvline(p95, ls="--", lw=1.1, color="#666", zorder=4)
    ax.text(p95 + .004, ax.get_ylim()[1] * .82, "95th pct", fontsize=7.5, color="#666")
    for v, lab, c in [(0.1042, "Glycolysis, all 10 donors", RED),
                      (0.0285, "Glycolysis, 8 patients only", BLUE)]:
        ax.axvline(v, lw=2, color=c, zorder=5)
        ax.text(v + .004, ax.get_ylim()[1] * (.98 if c == RED else .62), lab,
                fontsize=8, color=c)
    ax.axvline(0.0403, ls=":", lw=1.4, color=BLUE, zorder=4)
    ax.text(0.0403 + .004, ax.get_ylim()[1] * .45,
            "95th pct of the floor\nrecomputed within patients", fontsize=7.5, color=BLUE)
    ax.set_xlabel("Intraclass correlation across donors")
    ax.set_ylabel("Random modules")
    ax.text(.98, .55,
            "Upper-bound control: XIST ICC = 0.796\n"
            "(and it correctly separates the two female\n"
            "healthy donors from eight male patients)",
            transform=ax.transAxes, ha="right", fontsize=7.6, color=GREEN)
    ax.legend(frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1, .78))
    ax.set_title("Fig S8 | The precondition test that closed a planned analysis\n"
                 "Within a clinically homogeneous patient group, the programme is\n"
                 "indistinguishable from random gene modules",
                 fontsize=9.5, loc="left", pad=8)
    fig.tight_layout(); save(fig, "FigS8_icc")


for f in (s2, s3, s4, s5, s6, s7, s8):
    try:
        f()
    except Exception as e:
        print(f"[FAIL] {f.__name__}: {type(e).__name__}: {e}", flush=True)
