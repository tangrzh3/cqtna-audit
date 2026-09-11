"""Step 122 -- re-score every published cell under the corrected locus-level
convention, and report what moved.

External review of the R package found that module A mixed three conventions:
the background counted a locus as known if ANY of its records was near a known
lead SNP, the numerator counted only the SIGNIFICANT records, and the gene
labels came from per-record flags. The same mismatch is in this project's own
scoring machinery (step85 / step94c / step119 / step121), because the R package
was ported from it.

Corrected convention, applied everywhere:
  a locus is known if ANY record at that locus lies within KNOWN_KB of a known
  lead SNP; the significant loci inherit that status; the distance reported for
  a locus is the minimum over the same record set the flag used.

This script recomputes the four cells that have per-record tables and prints an
old-vs-new comparison, so nothing is quietly carried over.

Output: 122a_convention_recheck.tsv
"""
import os
import sys
from math import lgamma, exp

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
KB = [100, 250, 500, 1000]


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, min(r1, c1) + 1)), 1.0)


def assign_loci(df, locus_kb):
    out = {}
    for ch, sub in df.groupby("chr"):
        sub = sub.sort_values("pos")
        lid, prev = 0, None
        for _, r in sub.iterrows():
            if prev is not None and r.pos - prev > locus_kb * 1000:
                lid += 1
            out[(ch, r.pos)] = f"{ch}_{lid}"
            prev = r.pos
    return out


def known_map(path):
    k = pd.read_csv(f"{MR}/{path}")
    k["chr"] = k["chr"].astype(str)
    return {c: np.sort(s.pos.values) for c, s in k.groupby("chr")}


def nearest_bp(by_chr, ch, pos):
    arr = by_chr.get(str(ch))
    if arr is None or not len(arr):
        return np.inf
    i = np.searchsorted(arr, pos)
    return min((abs(int(arr[j]) - int(pos))
                for j in (i - 1, i) if 0 <= j < len(arr)), default=np.inf)


def score(d, kmap, locus_kb, known_kb, convention):
    d = d.copy()
    snp = d.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp, locus_kb)
    d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
    d["dist"] = [nearest_bp(kmap, c, p) for c, p in zip(d.chr, d.pos)]
    d["near"] = d.dist <= known_kb * 1000
    sig = d[d.fdr < 0.05]

    bg = d.groupby("locus").near.any()
    if convention == "old":
        # 分母走全部记录，分子只走显著记录 —— 就是被查出的那个不一致
        sg = sig.groupby("locus").near.any()
    else:
        sg = bg[sig.locus.unique()]
    BT, BK, ST, SK = len(bg), int(bg.sum()), len(sg), int(sg.sum())
    fold = ((SK / ST) / (BK / BT)) if ST and BK else np.nan
    p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)) if ST else np.nan
    if convention == "old":
        dist = sig.groupby("locus").dist.min()
    else:
        dist = d.groupby("locus").dist.min()[sig.locus.unique()]
    dist = sorted(v for v in dist.values if np.isfinite(v))
    return dict(bg_loci=BT, bg_known=BK, sig_loci=ST, sig_known=SK,
                fold=round(fold, 2) if ST and BK else None, fisher_p=p,
                distances_kb=[round(v / 1000) for v in dist])


def main():
    mel = known_map("landi2020_known_loci_grch38.csv")
    hcc = known_map("84a_hcc_known_loci_grch38.csv")
    cells = {}

    d = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
    d[["chr", "pos"]] = d.SNP.str.split(":", expand=True)
    d = d.rename(columns={"FDR": "fdr"})
    cells["melanoma x Soskic_CD4"] = (d[["chr", "pos", "fdr"]], mel)
    d = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    cells["melanoma x eQTLGen_blood"] = (d[["chr", "pos", "fdr"]], mel)
    for tag in ("high", "low"):
        d = pd.read_csv(f"{MR}/85a_HCC_{tag}_annotated.tsv", sep="\t")
        cells[f"HCC_{tag} x Soskic_CD4"] = (d[["chr", "pos", "fdr"]], hcc)
    for _, (d, _) in cells.items():
        d["chr"] = d["chr"].astype(str)
        d["pos"] = d["pos"].astype(int)
        d["fdr"] = pd.to_numeric(d["fdr"], errors="coerce")

    rows = []
    print("=" * 86)
    print("Step 122 -- old (mixed) vs corrected (locus-level) convention")
    print("=" * 86)
    for name, (d, kmap) in cells.items():
        print(f"\n--- {name} " + "-" * (66 - len(name)))
        for kb in KB:
            o = score(d, kmap, 1000, kb, "old")
            n = score(d, kmap, 1000, kb, "new")
            same = (o["sig_known"] == n["sig_known"] and o["fold"] == n["fold"])
            print(f"  known_kb {kb:>4}  old {o['sig_known']}/{o['sig_loci']} "
                  f"fold {o['fold']}  P {o['fisher_p']:.4g}   ->   "
                  f"new {n['sig_known']}/{n['sig_loci']} fold {n['fold']} "
                  f"P {n['fisher_p']:.4g}   {'unchanged' if same else '*** CHANGED ***'}")
            rows.append(dict(cell=name, known_kb=kb,
                             old_sig_known=o["sig_known"], old_fold=o["fold"],
                             old_fisher_p=o["fisher_p"],
                             new_sig_known=n["sig_known"], new_fold=n["fold"],
                             new_fisher_p=n["fisher_p"],
                             changed=not same))
        o = score(d, kmap, 1000, 1000, "old")
        n = score(d, kmap, 1000, 1000, "new")
        print(f"  distances (kb)  old {o['distances_kb']}")
        print(f"                  new {n['distances_kb']}"
              f"   {'unchanged' if o['distances_kb'] == n['distances_kb'] else '*** CHANGED ***'}")
        rows.append(dict(cell=name, known_kb="distances",
                         old_sig_known=";".join(map(str, o["distances_kb"])),
                         old_fold=None, old_fisher_p=None,
                         new_sig_known=";".join(map(str, n["distances_kb"])),
                         new_fold=None, new_fisher_p=None,
                         changed=o["distances_kb"] != n["distances_kb"]))

    out = pd.DataFrame(rows)
    out.to_csv(f"{MR}/122a_convention_recheck.tsv", sep="\t", index=False)
    ch = out[out.changed]
    print("\n" + "=" * 86)
    print(f"rows compared: {len(out)};  changed: {len(ch)}")
    if len(ch):
        for _, r in ch.iterrows():
            print(f"  CHANGED  {r.cell}  at {r.known_kb}")
    print("wrote 122a_convention_recheck.tsv")


if __name__ == "__main__":
    main()
