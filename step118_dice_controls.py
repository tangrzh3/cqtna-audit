#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 118 -- the three controls S35 §3 registered, run on the one cell that qualified.

Only DICE (Schmiedel_2018, CD4 anti-CD3/CD28 4 h) cleared both gates in step 116,
at 3.36-fold with a one-sided Fisher P of 0.112. The controls were registered to
guard against an artefactual positive; the cell is not significant, so they cannot
rescue or overturn anything. They are run because they were registered, and
because the attenuation pattern is reportable in its own right.

  control 1  mismatched list: score the same loci against HCC's known loci
  control 2  strength-frequency matched permutation: eQTL-p decile x outcome
             EAF quintile, 10,000 resamples
  control 3  density matched permutation: records-per-locus quintile x
             genes-per-locus tertile, 10,000 resamples, plus the span proxy

Every stratum relaxation is counted and printed, as S30 requires.

Inputs:  116b_mr_Schmiedel_2018.tsv, landi2020_known_loci_grch38.csv,
         84a_hcc_known_loci_grch38.csv
Outputs: 118a_dice_controls.tsv, 118_console.log
"""
import os
import io
import sys
from math import exp, lgamma

import numpy as np
import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
SRC = f"{MR}/116b_mr_Schmiedel_2018.tsv"
KNOWN_MEL = f"{MR}/landi2020_known_loci_grch38.csv"
KNOWN_HCC = f"{MR}/84a_hcc_known_loci_grch38.csv"
KNOWN_KB = 1000
N_PERM = 10000
SEED = 118


class Tee:
    def __init__(self, p):
        self.f = io.open(p, "w", encoding="utf-8")

    def write(self, s):
        enc = sys.__stdout__.encoding or "utf-8"
        sys.__stdout__.write(s.encode(enc, "replace").decode(enc, "replace"))
        self.f.write(s)

    def flush(self):
        sys.__stdout__.flush()
        self.f.flush()


sys.stdout = Tee(f"{MR}/118_console.log")


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    hi = min(r1, c1)
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, hi + 1)), 1.0)


def known_map(path, chr_col="chr", pos_col="pos"):
    k = pd.read_csv(path)
    k.columns = [c.lower() for c in k.columns]
    return {str(c): np.sort(s[pos_col].values) for c, s in k.groupby(chr_col)}


def is_known(kn, ch, pos):
    arr = kn.get(str(ch))
    if arr is None or not len(arr):
        return False
    i = np.searchsorted(arr, pos)
    return any(0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000
               for j in (i - 1, i))


m = pd.read_csv(SRC, sep="\t")
m["chr"] = m["chr"].astype(str)
print("=" * 78)
print("Step 118 -- S35 registered controls, DICE cell")
print("=" * 78)
print(f"records {len(m)}   significant records {(m.fdr < 0.05).sum()}")

# locus-level table with the covariates the permutations need
loc = (m.groupby("locus")
         .agg(chr=("chr", "first"), pos=("pos", "min"),
              n_records=("gene_id", "size"), n_genes=("gene_id", "nunique"),
              span=("pos", lambda s: int(s.max() - s.min())),
              min_pexp=("pval_exp", "min"), eaf=("eaf", "median"),
              known=("known", "any"), sig=("fdr", lambda s: bool((s < 0.05).any())))
         .reset_index())
BG_T = len(loc)
BG_K = int(loc.known.sum())
sig = loc[loc.sig]
S_T, S_K = len(sig), int(sig.known.sum())
fold0 = (S_K / S_T) / (BG_K / BG_T)
print(f"background loci {BG_T} ({BG_K} known)   significant {S_T} ({S_K} known)"
      f"   fold {fold0:.2f}")

rows = []

# ---- control 1: mismatched list -------------------------------------------
print("\n" + "=" * 78)
print("control 1 -- mismatched list (HCC known loci on a melanoma outcome)")
print("=" * 78)
kn_hcc = known_map(KNOWN_HCC)
loc["known_mis"] = [is_known(kn_hcc, c, p) for c, p in zip(loc.chr, loc.pos)]
BGm = int(loc.known_mis.sum())
Sm = int(loc.loc[loc.sig, "known_mis"].sum())
if BGm == 0:
    print("  no background locus matches the mismatched list -> fold undefined")
    foldm = pm = float("nan")
else:
    foldm = (Sm / S_T) / (BGm / BG_T)
    pm = fisher_greater(Sm, S_T - Sm, BGm - Sm, BG_T - S_T - (BGm - Sm))
    print(f"  background {BGm}/{BG_T} known under HCC list; significant {Sm}/{S_T}")
    print(f"  fold {foldm:.2f}   one-sided Fisher P {pm:.4g}")
rows.append(dict(control="mismatched_list", fold=foldm, p=pm,
                 note=f"HCC list; bg {BGm}/{BG_T}, sig {Sm}/{S_T}"))

# ---- permutation helper ----------------------------------------------------
rng = np.random.default_rng(SEED)


def stratified_perm(strata, label):
    """resample S_T background loci matching the significant loci's strata"""
    bg = loc[~loc.sig].reset_index(drop=True)
    tgt = loc[loc.sig].reset_index(drop=True)
    relax = 0
    draws = np.empty(N_PERM, dtype=int)
    pools = {}
    for s in bg[strata].unique():
        pools[s] = bg.index[bg[strata] == s].to_numpy()
    allidx = bg.index.to_numpy()
    for it in range(N_PERM):
        got = 0
        for s in tgt[strata]:
            pool = pools.get(s, allidx)
            if len(pool) == 0:
                pool = allidx
                if it == 0:
                    relax += 1
            j = rng.choice(pool)
            got += bool(bg.known.iloc[j])
        draws[it] = got
    exp_k = draws.mean()
    fold = (S_K / exp_k) if exp_k > 0 else float("inf")
    p = (1 + int((draws >= S_K).sum())) / (N_PERM + 1)
    print(f"  {label}: expected known {exp_k:.2f}, observed {S_K}"
          f"  -> fold {fold:.2f}, empirical P {p:.4g}"
          f"   (stratum relaxations {relax})")
    return fold, p, relax


# ---- control 2: strength x frequency ---------------------------------------
print("\n" + "=" * 78)
print("control 2 -- eQTL-p decile x outcome EAF quintile")
print("=" * 78)
loc["pdec"] = pd.qcut(loc.min_pexp.rank(method="first"), 10, labels=False)
loc["eafq"] = pd.qcut(loc.eaf.rank(method="first"), 5, labels=False)
loc["strat2"] = loc.pdec.astype(str) + "_" + loc.eafq.astype(str)
f2, p2, r2 = stratified_perm("strat2", "strength x frequency")
rows.append(dict(control="strength_frequency", fold=f2, p=p2,
                 note=f"decile x quintile; relaxations {r2}"))

# ---- control 3: density ----------------------------------------------------
print("\n" + "=" * 78)
print("control 3 -- density (records quintile x genes tertile), and span proxy")
print("=" * 78)
loc["recq"] = pd.qcut(loc.n_records.rank(method="first"), 5, labels=False)
loc["gent"] = pd.qcut(loc.n_genes.rank(method="first"), 3, labels=False,
                      duplicates="drop")
loc["strat3"] = loc.recq.astype(str) + "_" + loc.gent.astype(str)
f3, p3, r3 = stratified_perm("strat3", "density")
rows.append(dict(control="density", fold=f3, p=p3,
                 note=f"records quintile x genes tertile; relaxations {r3}"))

loc["spanq"] = pd.qcut(loc.span.rank(method="first"), 5, labels=False)
loc["strat3b"] = loc.spanq.astype(str) + "_" + loc.gent.astype(str)
f3b, p3b, r3b = stratified_perm("strat3b", "density (span proxy)")
rows.append(dict(control="density_span_proxy", fold=f3b, p=p3b,
                 note=f"span quintile x genes tertile; relaxations {r3b}"))

out = pd.DataFrame([dict(control="primary_fisher", fold=fold0,
                         p=fisher_greater(S_K, S_T - S_K, BG_K - S_K,
                                          BG_T - S_T - (BG_K - S_K)),
                         note=f"bg {BG_K}/{BG_T}, sig {S_K}/{S_T}")] + rows)
out.to_csv(f"{MR}/118a_dice_controls.tsv", sep="\t", index=False)

print("\n" + "=" * 78)
print(out.to_string(index=False))
print("=" * 78)
print("\n⚠ Read with the numerator in mind: 2 known of 6 significant loci. These")
print("controls cannot make a non-significant cell significant, and are reported")
print("because S35 §3 registered them, not as support for a positive result.")
