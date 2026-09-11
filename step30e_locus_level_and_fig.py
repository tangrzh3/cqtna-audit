#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 30e  (1) redo the locus-attribution enrichment counting INDEPENDENT LOCI
              rather than gene records, and (2) draw the cross-cancer figure.

WHY LOCI, NOT RECORDS
Step 30d counted FDR<0.05 gene records. Genes at the same locus are not
independent: lung's apparent 2.83x enrichment came largely from FADS1 +
TMEM258 (one chr11 cluster) plus HLA-DPA1. Melanoma's came from several
genuinely separate loci (MC1R, PARP1, ...). Counting loci removes that
inflation and is the number that belongs in the manuscript.

Loci are defined by 1 Mb single-linkage clustering of instrument positions
within a chromosome. A locus counts as "known" if any of its instruments is
annotated to a Landi 2020 melanoma/nevus/pigmentation lead SNP.

Output: 36e (table), figures/Fig7_crosscancer.{png,pdf}
"""
import sys
import os
import csv
import collections
from math import lgamma, exp
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
OUT = os.path.join(MR, "figures")
NOVEL = "潜在新位点"
LOCUS_KB = 1000

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.linewidth": .8, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    hi = min(r1, c1)
    denom = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - denom)
                   for x in range(a, hi + 1)), 1.0)


# ------------------------------------------------------------------ loci
ann = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
ann[["chrom", "bp"]] = ann.SNP.str.split(":", expand=True)
ann["bp"] = ann.bp.astype(int)

snp = (ann.groupby("SNP")
          .agg(chrom=("chrom", "first"), bp=("bp", "first"),
               known=("category", lambda s: (s != NOVEL).any()))
          .reset_index())

locus_of = {}
for ch, sub in snp.groupby("chrom"):
    sub = sub.sort_values("bp")
    lid, prev = 0, None
    for _, r in sub.iterrows():
        if prev is not None and r.bp - prev > LOCUS_KB * 1000:
            lid += 1
        locus_of[r.SNP] = f"{ch}_{lid}"
        prev = r.bp
snp["locus"] = snp.SNP.map(locus_of)
loc = snp.groupby("locus").agg(known=("known", "any")).reset_index()
BG_TOT, BG_KNOWN = len(loc), int(loc.known.sum())
print(f"instrument SNPs: {len(snp):,} -> independent loci ({LOCUS_KB} kb): {BG_TOT:,}")
print(f"background known-locus share: {BG_KNOWN}/{BG_TOT} = {100*BG_KNOWN/BG_TOT:.1f}%\n")

ann["locus"] = ann.SNP.map(locus_of)
mel = ann[ann.FDR < 0.05]
rows = [dict(cancer="Melanoma", n_loci=mel.locus.nunique(),
             n_known=mel[mel.category != NOVEL].locus.nunique())]

cc = pd.read_csv(f"{MR}/36a_crosscancer_MR_all.tsv", sep="\t")
cc["locus"] = cc.SNP.map(locus_of)
for c in ["Lung", "Colorectal", "Pancreas", "Breast", "Prostate"]:
    s = cc[(cc.cancer == c) & (cc.FDR < 0.05)]
    rows.append(dict(cancer=c, n_loci=s.locus.nunique(),
                     n_known=s[s.category != NOVEL].locus.nunique()))

res = pd.DataFrame(rows)
res["pct_known"] = (100 * res.n_known / res.n_loci).round(1)
res["bg_pct"] = round(100 * BG_KNOWN / BG_TOT, 1)
res["enrichment"] = (res.pct_known / res.bg_pct).round(2)
res["fisher_p"] = [fisher_greater(int(r.n_known), int(r.n_loci - r.n_known),
                                  BG_KNOWN - int(r.n_known),
                                  (BG_TOT - BG_KNOWN) - int(r.n_loci - r.n_known))
                   if r.n_loci else np.nan for r in res.itertuples()]
res.to_csv(f"{MR}/36e_locus_level_attribution.tsv", sep="\t", index=False)
print("=== enrichment of known melanoma loci, counted by INDEPENDENT LOCUS ===")
print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

# ------------------------------------------------------------------ figure
CANC = ["Melanoma", "Lung", "Colorectal", "Pancreas", "Breast", "Prostate"]
GENES = ["TPI1", "SMC2", "ZFYVE19", "SPSB2", "HLA-C", "KIAA0040"]

mel_best = (ann.sort_values("pval").groupby("SYMBOL")
               .first().reset_index()[["SYMBOL", "OR", "pval"]])
mat_or = pd.DataFrame(index=GENES, columns=CANC, dtype=float)
mat_p = pd.DataFrame(index=GENES, columns=CANC, dtype=float)
for g in GENES:
    r = mel_best[mel_best.SYMBOL == g]
    if len(r):
        mat_or.loc[g, "Melanoma"] = r.OR.iloc[0]
        mat_p.loc[g, "Melanoma"] = r.pval.iloc[0]
for c in CANC[1:]:
    sub = cc[cc.cancer == c].sort_values("pval").groupby("symbol").first()
    for g in GENES:
        if g in sub.index:
            mat_or.loc[g, c] = sub.loc[g, "OR"]
            mat_p.loc[g, c] = sub.loc[g, "pval"]

fig = plt.figure(figsize=(11.4, 4.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.35, .78], wspace=.42,
                      left=.075, right=.98, top=.86, bottom=.20)

# panel a -- locus attribution
ax = fig.add_subplot(gs[0, 0])
d = res.set_index("cancer").reindex(CANC).dropna(subset=["n_loci"])
d = d[d.n_loci > 0]
y = np.arange(len(d))[::-1]
cols = ["#C4453C" if i == "Melanoma" else "#8C99A6" for i in d.index]
ax.barh(y, d.enrichment, color=cols, height=.65)
ax.axvline(1, color="k", lw=.9, ls="--", alpha=.7)
ax.set_yticks(y); ax.set_yticklabels(d.index, fontsize=8.5)
for yy, (e, p, n) in zip(y, zip(d.enrichment, d.fisher_p, d.n_loci)):
    ax.text(e + .12, yy, f"P={p:.2g}  (n={int(n)})", va="center", fontsize=7)
ax.set_xlabel("enrichment of known melanoma loci\n(vs %.1f%% background, by locus)"
              % d.bg_pct.iloc[0], fontsize=8)
ax.set_xlim(0, max(d.enrichment) * 1.55)
# NOT "melanoma-specific": counted by independent locus, prostate reaches the
# same nominal significance as melanoma and neither survives correction for the
# six outcomes tested. The honest statement is that melanoma has the largest
# point estimate on a small number of loci.
ax.set_title("Known-locus enrichment, counted by\nindependent locus", fontsize=9, pad=6)
ax.text(-.30, 1.14, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")

# panel b -- candidate OR heatmap
ax = fig.add_subplot(gs[0, 1])
M = np.log2(mat_or.values.astype(float))
v = np.nanmax(np.abs(M))
im = ax.imshow(M, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
ax.set_xticks(range(len(CANC)))
ax.set_xticklabels(CANC, rotation=35, ha="right", fontsize=8)
ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES, fontsize=8.5)
for i in range(len(GENES)):
    for j in range(len(CANC)):
        p = mat_p.iloc[i, j]
        if np.isnan(p):
            continue
        star = "***" if p < 1e-3 else "**" if p < .01 else "*" if p < .05 else ""
        if star:
            ax.text(j, i, star, ha="center", va="center", fontsize=9,
                    color="white" if abs(M[i, j]) > v * .55 else "black")
ax.set_title("Candidate effects across cancers\n(* P<.05  ** P<.01  *** P<.001)",
             fontsize=9, pad=6)
cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.02)
cb.set_label("log$_2$ OR", fontsize=8); cb.ax.tick_params(labelsize=7)
cb.outline.set_visible(False)
ax.text(-.22, 1.14, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")

# panel c -- hit count vs outcome case count
# Case counts from the FinnGen R12 manifest. Melanoma is the FinnGen+Rashkin
# meta (12,530), not a FinnGen endpoint, so it is drawn open to mark that it
# comes from a different outcome dataset.
CASES = {"Melanoma": 12530, "Lung": 9639, "Colorectal": 11790,
         "Pancreas": 3139, "Breast": 24270, "Prostate": 20368}
ax = fig.add_subplot(gs[0, 2])
d2 = res.set_index("cancer").reindex(CANC)
xs = np.array([CASES[c] for c in d2.index], float)
ys = d2.n_loci.values.astype(float)
for x, yv, c in zip(xs, ys, d2.index):
    is_mel = c == "Melanoma"
    ax.scatter(x / 1000, yv, s=52, zorder=3,
               facecolors="none" if is_mel else "#8C99A6",
               edgecolors="#C4453C" if is_mel else "none", linewidths=1.8)
    ax.annotate(c, (x / 1000, yv), textcoords="offset points",
                xytext=(5, 4), fontsize=7)
rho = pd.Series(ys).corr(pd.Series(xs), method="spearman")
ax.set_xlabel("outcome cases (thousands)", fontsize=8)
ax.set_ylabel("independent loci at FDR < 0.05", fontsize=8)
ax.set_title("Hit count tracks outcome size\n" +
             f"Spearman $\\rho$ = {rho:.2f} (n = 6)", fontsize=9, pad=6)
ax.margins(.16)
ax.text(-.38, 1.14, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")
print(f"\npanel c: Spearman rho(cases, loci) = {rho:.3f} over {len(xs)} outcomes")
print("  melanoma is the FinnGen+Rashkin meta, not a FinnGen endpoint")

fig.savefig(os.path.join(OUT, "Fig7_crosscancer.pdf"))
fig.savefig(os.path.join(OUT, "Fig7_crosscancer.png"), dpi=300)
print("\nwritten: 36e, figures/Fig7_crosscancer")
