#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 49  Is the chromatin axis just activation again?

Step 48 found glycolysis-high peaks enriched for AP-1 (BATF::JUN odds 2.22) --
almost the same motif ranking as the activation positive control (BATF::JUN
odds 2.64). Cells were matched on an RNA activation score and glycolysis is
nearly orthogonal to it (R2 = 0.024, Step 43), but an 11-gene RNA module is a
crude proxy for chromatin activation state, so the resemblance has to be tested
rather than explained away.

Three checks, on the same consensus reference:
  1. across peaks, how strongly does the axis log fold change correlate with the
     activation log fold change?
  2. how much do the top axis-up peaks overlap the top activation-opening peaks?
  3. restricted to peaks that do NOT change with activation, does the axis still
     separate them -- i.e. is there an axis signal outside the activation
     programme at all?

Output: 49a
"""
import numpy as np
import pandas as pd
from scipy import stats

MR = r"D:/R_ex/MR"
N_TOP = 2000

ax = pd.read_csv(f"{MR}/47a_axis_differential_peaks.tsv.gz", sep="\t")
ctl = pd.read_csv(f"{MR}/47b_activation_differential_peaks.tsv.gz", sep="\t")
d = ax.merge(ctl[["peak", "rest_cpm", "act_cpm", "lfc"]], on="peak",
             suffixes=("", "_act"))
d = d.rename(columns={"lfc": "act_lfc"})
d = d[d.consistent & (d.total >= 100)]
print(f"peaks tested: {len(d):,}")

# ---- 1. correlation ------------------------------------------------------
r_p = stats.pearsonr(d.mean_lfc, d.act_lfc)
r_s = stats.spearmanr(d.mean_lfc, d.act_lfc)
print("\n" + "=" * 72)
print("1  axis log2FC vs activation log2FC, across peaks")
print("=" * 72)
print(f"  Pearson  r = {r_p.statistic:+.3f}  (P = {r_p.pvalue:.3g})")
print(f"  Spearman r = {r_s.statistic:+.3f}")
print(f"  variance of the axis explained by activation: {r_p.statistic**2:.3f}")

# ---- 2. overlap of top sets ---------------------------------------------
top_ax = set(d.nlargest(N_TOP, "mean_lfc").peak)
top_ct = set(d.nlargest(N_TOP, "act_lfc").peak)
ov = len(top_ax & top_ct)
exp = N_TOP * N_TOP / len(d)
print("\n" + "=" * 72)
print(f"2  overlap of the top {N_TOP} peaks in each contrast")
print("=" * 72)
print(f"  observed {ov}, expected by chance {exp:.0f}, fold {ov/exp:.2f}")
tbl = [[ov, N_TOP - ov], [N_TOP - ov, len(d) - 2 * N_TOP + ov]]
print(f"  Fisher P = {stats.fisher_exact(tbl, alternative='greater')[1]:.3g}")

# ---- 3. axis signal among activation-invariant peaks ---------------------
print("\n" + "=" * 72)
print("3  axis signal restricted to peaks that do not change with activation")
print("=" * 72)
q = d.act_lfc.abs().quantile(0.5)
inv = d[d.act_lfc.abs() <= q]
print(f"  activation-invariant peaks (|act log2FC| <= {q:.3f}): {len(inv):,}")
print(f"  their axis log2FC: mean {inv.mean_lfc.mean():+.4f}, "
      f"sd {inv.mean_lfc.std():.4f}, "
      f"range {inv.mean_lfc.min():+.3f} to {inv.mean_lfc.max():+.3f}")
print(f"  axis log2FC spread in all peaks:      sd {d.mean_lfc.std():.4f}")
print(f"  -> spread retained: {100*inv.mean_lfc.std()/d.mean_lfc.std():.0f}%")

out = inv.nlargest(N_TOP, "mean_lfc")[["peak"]]
out.to_csv(f"{MR}/49a_axis_up_activation_invariant.tsv", sep="\t", index=False)
inv.nsmallest(N_TOP, "mean_lfc")[["peak"]].to_csv(
    f"{MR}/49b_axis_down_activation_invariant.tsv", sep="\t", index=False)
inv[["peak", "mean_lfc", "act_lfc", "total"]].to_csv(
    f"{MR}/49c_activation_invariant_peaks.tsv.gz", sep="\t", index=False,
    compression="gzip")

print("\n" + "=" * 72)
print("verdict")
print("=" * 72)
r2 = r_p.statistic ** 2
if r2 > 0.25 or ov / exp > 3:
    print(f"  the chromatin axis overlaps the activation programme substantially")
    print(f"  (R2 = {r2:.3f}, top-set overlap {ov/exp:.1f}x). The AP-1 signal")
    print(f"  cannot be attributed to glycolysis independently of activation")
    print(f"  from these data alone; re-run motif enrichment using")
    print(f"  activation-invariant peaks only (49a/49b written for that).")
else:
    print(f"  the chromatin axis is largely distinct from the activation")
    print(f"  programme (R2 = {r2:.3f}, top-set overlap {ov/exp:.1f}x)")
print("\nwritten: 49a-49c")
