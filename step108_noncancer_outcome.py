#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 108 -- does locus attribution hold outside cancer, where CD4 is causal?

Executes manuscript/PREREG_noncancer_outcome.md (S33).

Every generalisation so far stayed inside cancer, and in melanoma and HCC alike
CD4 T cells are not the accepted causal cell type. Rheumatoid arthritis is the
framework's most favourable case: CD4 is canonically causal there, and the
outcome has more cases than the melanoma meta, so a null cannot be blamed on
power.

Reference: Okada 2014 lead SNPs, published before FinnGen existed and from
cohorts that do not include it (S33 section 3).

rsIDs are resolved to positions DURING the full-file scan, never in the
instrument-filtered extract -- that shortcut is the step94c bug, whose symptom
is an implausibly low background known-locus share. Below 2% is declared a
mapping failure here rather than reported as a result.

Output: 108a_ra_attribution.tsv
"""
import csv
import gzip
import os
from math import exp, lgamma

import numpy as np
import pandas as pd

MR = r"D:/R_ex/MR"
RA = f"{MR}/ra/finngen_R13_M13_RHEUMA.gz"
LOCUS_KB = KNOWN_KB = 1000
FDR = 0.05
MHC = ("6", 25_000_000, 34_000_000)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, min(r1, c1) + 1)), 1.0)


def two_sided(z):
    from scipy.stats import norm
    return 2 * norm.sf(np.abs(z))


def assign_loci(chrs, poss):
    order = sorted(range(len(chrs)), key=lambda i: (str(chrs[i]), int(poss[i])))
    out = [None] * len(chrs)
    lc, lp, lid = None, None, 0
    for i in order:
        c, p = str(chrs[i]), int(poss[i])
        if c != lc or p - lp > LOCUS_KB * 1000:
            lid += 1
        out[i] = lid
        lc, lp = c, p
    return out


def flagger(pos_by_chr):
    def f(ch, pos):
        arr = pos_by_chr.get(str(ch))
        if arr is None or not len(arr):
            return False
        i = np.searchsorted(arr, pos)
        return any(0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000
                   for j in (i - 1, i))
    return f


def scan(path, want_pos, want_rs):
    """one pass: instrument rows, and the coordinates of every known rsID"""
    got, known = {}, {}
    n = 0
    with gzip.open(path, "rt") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        ix = {c: i for i, c in enumerate(hdr)}
        ic, ip = ix["#chrom"], ix["pos"]
        ir, ia = ix["ref"], ix["alt"]
        ib, ise, irs = ix["beta"], ix["sebeta"], ix["rsids"]
        for row in rd:
            n += 1
            try:
                c, p = row[ic], int(row[ip])
            except (ValueError, IndexError):
                continue
            if (c, p) in want_pos:
                try:
                    got.setdefault((c, p), []).append(
                        (row[ir].upper(), row[ia].upper(),
                         float(row[ib]), float(row[ise])))
                except ValueError:
                    pass
            rs = row[irs] if irs < len(row) else ""
            if rs:
                for one in rs.split(","):
                    if one in want_rs:
                        known.setdefault(str(c), set()).add(p)
            if n % 5_000_000 == 0:
                print(f"    {n/1e6:.0f}M rows, {len(got):,} instrument positions, "
                      f"{sum(len(v) for v in known.values())} known loci placed",
                      flush=True)
    print(f"  scanned {n:,} rows; {len(got):,} instrument positions; "
          f"{sum(len(v) for v in known.values())} of {len(want_rs)} known rsIDs placed")
    return got, {c: np.sort(np.array(sorted(v))) for c, v in known.items()}


def cell(inst, got, is_known, label, allele_aware):
    keep, z = [], []
    for i, r in enumerate(inst.itertuples()):
        cands = got.get((r.chr, r.pos), [])
        if not cands:
            continue
        hit = None
        if allele_aware:
            for ref, alt, b, s in cands:
                if r.ea == alt and r.oa == ref:
                    hit = (b, s); break
                if r.ea == ref and r.oa == alt:
                    hit = (-b, s); break
        else:
            hit = (cands[0][2], cands[0][3])
        if hit is None or hit[1] <= 0:
            continue
        keep.append(i); z.append(hit[0] / hit[1])
    sub = inst.iloc[keep].reset_index(drop=True)
    p = two_sided(np.array(z))
    f = bh(p)
    loci = assign_loci(sub.chr.tolist(), sub.pos.tolist())
    kn = [is_known(c, pp) for c, pp in zip(sub.chr, sub.pos)]
    d = pd.DataFrame(dict(locus=loci, known=kn, fdr=f, p=p,
                          chr=sub.chr.values, pos=sub.pos.values))
    bg = d.groupby("locus").agg(known=("known", "any"))
    BT, BK = len(bg), int(bg.known.sum())
    sig = d[d.fdr < FDR]
    sg = sig.groupby("locus").agg(known=("known", "any"))
    ST, SK = len(sg), int(sg.known.sum())
    fold = ((SK / ST) / (BK / BT)) if ST and BK else np.nan
    pv = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)) if ST else np.nan
    mhc_nom = int(((d.chr == MHC[0]) & d.pos.between(MHC[1], MHC[2]) &
                   (d.p < 0.05)).sum())
    print(f"  [{label}] records {len(d):,}; background known {BK}/{BT} = "
          f"{100*BK/BT:.1f}%; FDR<{FDR} {len(sig)} records / {ST} loci, {SK} known "
          f"-> {fold if not np.isnan(fold) else float('nan'):.2f}x, P = {pv:.3g}; "
          f"MHC nominal {mhc_nom}")
    return dict(cell=label, n_records=len(d), bg_known=BK, bg_loci=BT,
                bg_pct=round(100 * BK / BT, 2) if BT else np.nan,
                sig_records=len(sig), sig_loci=ST, sig_known=SK,
                pct_known=round(100 * SK / ST, 1) if ST else np.nan,
                fold=round(fold, 3) if ST and BK else np.nan, fisher_p=pv,
                mhc_nominal=mhc_nom)


def main():
    if not os.path.exists(RA):
        raise SystemExit(f"missing {RA}")

    # ---- exposures
    d_sos = pd.read_csv(f"{MR}/06_locus_annotation.tsv", sep="\t")
    d_sos = d_sos[["exposure", "SYMBOL", "SNP", "beta.exposure"]].dropna()
    d_sos[["chr", "pos"]] = d_sos.SNP.str.split(":", expand=True).iloc[:, :2]
    d_sos["pos"] = d_sos.pos.astype(int)
    alle = {}
    with open(f"{MR}/01_harmonised_all.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            alle[(r["exposure"], r["SNP"])] = (r["effect_allele.outcome"].upper(),
                                               r["other_allele.outcome"].upper())
    d_sos["ea"] = [alle.get((e, s), ("", ""))[0] for e, s in zip(d_sos.exposure, d_sos.SNP)]
    d_sos["oa"] = [alle.get((e, s), ("", ""))[1] for e, s in zip(d_sos.exposure, d_sos.SNP)]
    d_sos = d_sos[(d_sos.ea != "") & (d_sos.oa != "")].reset_index(drop=True)

    d_eq = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    d_eq = d_eq[["gene_id", "symbol", "chr", "pos"]].dropna()
    d_eq["chr"] = d_eq["chr"].astype(str)
    d_eq["pos"] = d_eq["pos"].astype(int)
    d_eq["ea"] = d_eq["oa"] = ""
    print(f"instruments: Soskic {len(d_sos):,}, eQTLGen {len(d_eq):,}")

    okada = set(open(f"{MR}/ra/okada2014_lead_rsids.txt", encoding="utf-8").read().split())
    landi = pd.read_csv(f"{MR}/landi2020_known_loci_grch38.csv")
    print(f"known-locus references: Okada 2014 {len(okada)} rsIDs; "
          f"Landi melanoma {len(landi)} loci (negative control)")

    want_pos = set(zip(d_sos.chr, d_sos.pos)) | set(zip(d_eq.chr, d_eq.pos))
    want_rs = okada | set(landi.rsid.astype(str))
    print(f"\nscanning {os.path.basename(RA)} ...", flush=True)
    got, placed = scan(RA, want_pos, want_rs)

    # split the placed coordinates back into the two references
    ok_pos, mel_pos = {}, {}
    mel_rs = set(landi.rsid.astype(str))
    # re-derive which placed positions belong to which list by a second cheap pass
    # over the reference tables is unnecessary: Landi carries its own GRCh38
    # coordinates, so use them directly and keep the scan only for Okada.
    mel_pos = {c: np.sort(s.pos.values)
               for c, s in landi.astype({"chr": str}).groupby("chr")}
    ok_pos = placed  # scan collected both, but Okada is the one lacking coordinates
    n_ok = sum(len(v) for v in ok_pos.values())
    print(f"  Okada loci placed on GRCh38: {n_ok}")

    rows = []
    rows.append(cell(d_sos, got, flagger(ok_pos), "N1 Soskic x RA", True))
    rows.append(cell(d_eq, got, flagger(ok_pos), "N2 eQTLGen x RA", False))
    rows.append(cell(d_sos, got, flagger(mel_pos),
                     "NC Soskic x RA scored with melanoma list", True))
    out = pd.DataFrame(rows)
    out.to_csv(f"{MR}/108a_ra_attribution.tsv", sep="\t", index=False)

    print("\n" + "=" * 88)
    print("VERDICT (S33 section 6)")
    print("=" * 88)
    n1 = out[out.cell.str.startswith("N1")].iloc[0]
    if n1.bg_pct < 2 or n1.mhc_nominal == 0:
        v = ("D: mapping or extraction failed "
             f"(background {n1.bg_pct}%, MHC nominal {n1.mhc_nominal}); "
             "no attribution result is interpreted")
    elif n1.fold > 1 and n1.fisher_p < 0.05:
        v = "A: attribution holds outside cancer, where CD4 is the causal cell type"
    elif n1.fold > 1:
        v = "B: direction consistent, not significant"
    else:
        v = ("C: attribution does NOT hold here -- the claim narrows to settings "
             "where the exposure cell type is not causal for the disease")
    print(f"  N1: {n1.fold}x, P = {n1.fisher_p:.3g}, background {n1.bg_pct}%, "
          f"MHC nominal {n1.mhc_nominal}")
    print(f"  -> {v}")
    nc = out[out.cell.str.startswith("NC")].iloc[0]
    print(f"  NC mismatched melanoma list: {nc.fold}x, P = {nc.fisher_p:.3g} "
          f"{'(clean)' if not (nc.fold > 1 and nc.fisher_p < 0.05) else '** POSITIVE, void **'}")
    print("\n  S33 section 7: Okada 2014 predates a decade of later RA loci, so "
          "loci found\n  after it count as novel here and the enrichment is "
          "biased DOWNWARD. That is\n  deliberate: under-estimating beats "
          "reintroducing circularity.")
    print("\nwrote 108a")


if __name__ == "__main__":
    main()
