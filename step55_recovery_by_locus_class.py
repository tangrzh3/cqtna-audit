#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 55  Recovery curves split by locus class.

Step 54 showed melanoma's candidate list is the most robust to loss of outcome
power of the six diseases -- 60% recovered at half power, against 19-24% for
lung, breast and prostate. That looks reassuring until it is put beside finding
1: melanoma's list is dominated by known pigmentation and naevus loci, whose
effects are enormous (MC1R region, P to 4e-37). Large effects survive
down-sampling. So the apparent robustness may belong entirely to the part of the
list that finding 1 says is not the immune signal being sought.

This tests that directly: recovery is recomputed separately for genes at known
pigmentation/naevus/melanoma loci and for genes at novel loci.

PREDICTION, stated before running
  known loci recover quickly, novel loci slowly. If so, the stable part of the
  list is the part that is not the finding, and the part that is the finding is
  the least reproducible.

Output: 55a
"""
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
rng = np.random.default_rng(1)
N_REP = 400
FDR_THR = 0.05
NOVEL = "潜在新位点"
META_CASE, META_CTRL = 12530, 789099


def n_eff(nc, nk):
    return 4.0 / (1.0 / nc + 1.0 / nk)


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def wald_p(bo, seo, be):
    return 2 * stats.norm.sf(np.abs((bo / be) / (seo / np.abs(be))))


d = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
d = d[["exposure", "SYMBOL", "category", "beta_exposure", "beta_outcome",
       "se_outcome"]].dropna()
d = d[(d.beta_exposure != 0) & (d.se_outcome > 0)].reset_index(drop=True)
d["is_novel"] = d.category == NOVEL

p_full = wald_p(d.beta_outcome.values, d.se_outcome.values, d.beta_exposure.values)
hit_full = bh(p_full) < FDR_THR
full_known = set(d.SYMBOL[hit_full & ~d.is_novel])
full_novel = set(d.SYMBOL[hit_full & d.is_novel])
print(f"exposures: {len(d):,}")
print(f"full-power list: {hit_full.sum()} records | "
      f"known-locus genes {len(full_known)}, novel-locus genes {len(full_novel)}")
print(f"  known: {sorted(full_known)}")
print(f"  novel: {sorted(full_novel)}")

NE_OBS = n_eff(META_CASE, META_CTRL)
rows = []
for frac in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
    nc = int(META_CASE * frac)
    ne = n_eff(nc, META_CTRL)
    scale = np.sqrt(NE_OBS / ne)
    se_s = d.se_outcome.values * scale
    extra = np.sqrt(np.maximum(se_s ** 2 - d.se_outcome.values ** 2, 0))
    rk, rn = [], []
    for _ in range(N_REP):
        bo = d.beta_outcome.values + rng.normal(0, extra)
        h = bh(wald_p(bo, se_s, d.beta_exposure.values)) < FDR_THR
        gk = set(d.SYMBOL[h & ~d.is_novel])
        gn = set(d.SYMBOL[h & d.is_novel])
        rk.append(len(gk & full_known) / len(full_known) if full_known else np.nan)
        rn.append(len(gn & full_novel) / len(full_novel) if full_novel else np.nan)
    rows.append(dict(frac=frac, cases=nc,
                     known_recovery=np.nanmean(rk),
                     known_lo=np.nanpercentile(rk, 2.5),
                     known_hi=np.nanpercentile(rk, 97.5),
                     novel_recovery=np.nanmean(rn),
                     novel_lo=np.nanpercentile(rn, 2.5),
                     novel_hi=np.nanpercentile(rn, 97.5)))
    print(f"  {nc:>6,} cases ({frac:.0%})   known {100*np.nanmean(rk):5.1f}%   "
          f"novel {100*np.nanmean(rn):5.1f}%")

res = pd.DataFrame(rows)
res.to_csv(f"{MR}/55a_recovery_by_locus_class.tsv", sep="\t", index=False)

print("\n" + "=" * 74)
print("verdict")
print("=" * 74)
half = res[res.frac == 0.5].iloc[0]
print(f"  at half power: known-locus genes recover {100*half.known_recovery:.0f}%, "
      f"novel-locus genes {100*half.novel_recovery:.0f}%")
gap = half.known_recovery - half.novel_recovery
print(f"  gap {100*gap:+.0f} percentage points")
for target in (0.5, 0.8):
    k = res[res.known_recovery >= target]
    n = res[res.novel_recovery >= target]
    ck = int(k.cases.min()) if len(k) else None
    cn = int(n.cases.min()) if len(n) else None
    print(f"  cases to recover {target:.0%}:  known "
          f"{f'{ck:,}' if ck else '>11,277'}   novel "
          f"{f'{cn:,}' if cn else '>11,277'}")
print("\nwritten: 55a")
