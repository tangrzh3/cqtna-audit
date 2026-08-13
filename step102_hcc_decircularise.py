#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 102 -- de-circularise the HCC known-locus reference.

Executes manuscript/PREREG_hcc_decircularisation.md (S29).

The concern, raised in review: the HCC known-locus list was built partly from the
outcome GWAS's own Table 1, so asking whether MR hits computed against that GWAS
land on "known" loci is partly circular. Provenance was resolved before this
script was written (S29 section 0): of 73 loci, exactly 3 rest on GCST90809296
alone -- rs4089 and rs7628416, which no other GWAS Catalog study reports, and
rs9277534, which is not in the Catalog at all. The other 7 paper-sourced entries
carry heavy independent support (APOE 1,399 studies, TM6SF2 1,073, PNPLA3 225),
so dropping all 10 would over-correct; that variant is kept as a bound, not as
the primary.

Three references are scored and all three reported (S29 section 3):
  F  full list, as published                                    73 loci
  P  primary: drop loci resting on the outcome GWAS alone       70 loci
  C  conservative bound: drop every paper-sourced entry         63 loci

Machinery is identical to step94c/step85; only the reference list changes.

Output: 102a_hcc_decirc.tsv
"""
from math import lgamma, exp

import numpy as np
import pandas as pd

MR = r"D:/R_ex/MR"
KNOWN_KB = 1000

# resting on GCST90809296 alone -- see S29 section 0.1 and hcc/known_locus_independence.json
OUTCOME_ONLY = {"rs4089", "rs7628416", "rs9277534"}


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, min(r1, c1) + 1)), 1.0)


def flagger(known):
    KN = {c: np.sort(s.pos.values) for c, s in known.astype({"chr": str}).groupby("chr")}

    def is_known(ch, pos):
        arr = KN.get(str(ch))
        if arr is None or not len(arr):
            return False
        i = np.searchsorted(arr, pos)
        return any(0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000
                   for j in (i - 1, i))
    return is_known


def score(d, is_known, label, ref):
    d = d.copy()
    d["known"] = [is_known(c, p) for c, p in zip(d.chr, d.pos)]
    bg = d.groupby("locus").agg(known=("known", "any"))
    BT, BK = len(bg), int(bg.known.sum())
    sig = d[d.fdr < 0.05]
    sg = sig.groupby("locus").agg(known=("known", "any"))
    ST, SK = len(sg), int(sg.known.sum())
    fold = ((SK / ST) / (BK / BT)) if ST and BK else np.nan
    p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)) if ST else np.nan
    print(f"  {label:9s} {ref:26s} background {BK:3d}/{BT} = {100*BK/BT:4.1f}%   "
          f"significant {SK}/{ST}   fold {fold:6.2f}   P = {p:.4g}")
    return dict(cell=label, reference=ref, n_known_list=None,
                bg_known=BK, bg_loci=BT, sig_known=SK, sig_loci=ST,
                pct_known=round(100 * SK / ST, 1) if ST else np.nan,
                fold=round(fold, 3) if ST and BK else np.nan, fisher_p=p)


def main():
    full = pd.read_csv(f"{MR}/84a_hcc_known_loci_grch38.csv")
    refs = {
        "F full (published)": full,
        "P drop outcome-only": full[~full.rsid.isin(OUTCOME_ONLY)],
        "C drop all paper-sourced": full[full.source != "paper_table1"],
    }
    print("reference lists:")
    for name, k in refs.items():
        print(f"  {name:26s} {len(k):3d} loci")
    dropped = sorted(set(full.rsid) & OUTCOME_ONLY)
    print(f"  outcome-only loci dropped in P: {dropped}")

    rows = []
    for cell, fn in (("HCC-high", "85a_HCC_high_annotated.tsv"),
                     ("HCC-low", "85a_HCC_low_annotated.tsv")):
        d = pd.read_csv(f"{MR}/{fn}", sep="\t")
        d["chr"] = d["chr"].astype(str)
        print(f"\n{cell}  ({len(d):,} MR records)")
        for name, k in refs.items():
            r = score(d, flagger(k), cell, name)
            r["n_known_list"] = len(k)
            rows.append(r)

    out = pd.DataFrame(rows)
    out.to_csv(f"{MR}/102a_hcc_decirc.tsv", sep="\t", index=False)

    print("\n" + "=" * 92)
    print("VERDICT (S29 section 5)")
    print("=" * 92)
    for cell in ("HCC-high", "HCC-low"):
        s = out[out.cell == cell].set_index("reference")
        f = s.loc["F full (published)"]
        p = s.loc["P drop outcome-only"]
        c = s.loc["C drop all paper-sourced"]
        print(f"  {cell}:  F {f.fold:.2f} (P={f.fisher_p:.4g})   "
              f"P {p.fold:.2f} (P={p.fisher_p:.4g})   "
              f"C {c.fold:.2f} (P={c.fisher_p:.4g})")
        if p.fold > 1 and np.sign(p.fold - 1) == np.sign(f.fold - 1):
            print(f"           -> primary reference keeps fold > 1, same direction as "
                  f"published: proposition D holds for this cell")
        else:
            print(f"           -> ** fold <= 1 or direction reversed under the primary "
                  f"reference: interpretation B, HCC generalisation withdrawn **")
    print("\n  Reminder (S29 section 0.3): the three de-circularised loci sit on "
          "chr3, chr4 and\n  chr6, none within 1 Mb of a significant HCC locus "
          "(chr19 ~19.2 Mb, chr22 ~43.9 Mb),\n  so this test could only leave the "
          "numerator unchanged. Its informative content is\n  the size of the "
          "circularity (3 of 73), not the robustness it appears to demonstrate.")
    print("\nwrote 102a")


if __name__ == "__main__":
    main()
