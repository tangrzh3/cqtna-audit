"""Step 94d -- the missing cell: eQTLGen exposure against both HCC outcomes.

This completes a 2 x 2 grid (melanoma, HCC) x (Soskic CD4, eQTLGen whole blood)
in which every cell is scored against its OWN disease's known-locus list, both of
which carry their own GRCh38 coordinates (Landi 157 for melanoma, the 73-locus
list built in step 84 for HCC). The five FinnGen cancers are deferred: their
GWAS Catalog rsID lists could not be placed on GRCh38 reliably (see step94b2).

Output: 94f_eqtlgen_hcc.tsv
"""
import csv
import gzip
import io
from math import lgamma, exp

import numpy as np
import pandas as pd
from scipy.stats import norm

MR = r"D:/R_ex/MR"
LOCUS_KB = KNOWN_KB = 1000


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, min(r1, c1) + 1)), 1.0)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    a = np.empty(n)
    a[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(a, 0, 1)


def assign_loci(df):
    out = {}
    for ch, sub in df.groupby("chr"):
        sub = sub.sort_values("pos")
        lid, prev = 0, None
        for _, r in sub.iterrows():
            if prev is not None and r.pos - prev > LOCUS_KB * 1000:
                lid += 1
            out[(ch, r.pos)] = f"{ch}_{lid}"
            prev = r.pos
    return out


inst = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
inst["chr"] = inst["chr"].astype(str)
want = set(zip(inst.chr, inst.pos))
print(f"eQTLGen instrument positions: {len(want):,}")

known = pd.read_csv(f"{MR}/84a_hcc_known_loci_grch38.csv")
known["chr"] = known["chr"].astype(str)
KN = {c: np.sort(s.pos.values) for c, s in known.groupby("chr")}
print(f"HCC known loci with GRCh38 coordinates: {len(known)}")


def scan(path, cols, label):
    got = {}
    op = gzip.open(path, "rt")
    rd = csv.reader(op, delimiter="\t")
    hdr = next(rd)
    ix = {c: i for i, c in enumerate(hdr)}
    n = 0
    for row in rd:
        n += 1
        try:
            key = (row[ix[cols["chr"]]], int(row[ix[cols["pos"]]]))
        except (ValueError, IndexError):
            continue
        if key in want and key not in got:
            try:
                got[key] = (float(row[ix[cols["beta"]]]), float(row[ix[cols["se"]]]))
            except ValueError:
                pass
        if n % 5_000_000 == 0:
            print(f"    {label}: {n:,} rows, {len(got):,} matched", end="\r")
    op.close()
    print(f"\n  {label}: {n:,} rows scanned, {len(got):,} instruments matched")
    return got


# 2026-08-18: only summary rows used to survive this script, so the cell
# could not be recomputed under a different locus partition. The
# standardised per-record table is now persisted; 94f is left untouched.
RECORDS = []


def cell(got, label):
    rows = []
    for r in inst.itertuples():
        m = got.get((r.chr, r.pos))
        if m is None:
            continue
        bo, so = m
        if so <= 0:
            continue
        z = bo / so
        rows.append(dict(chr=r.chr, pos=r.pos, gene_id=r.gene_id,
                         p_mr=2 * norm.sf(abs(z))))
    d = pd.DataFrame(rows)
    d["fdr"] = bh(d.p_mr.values)
    RECORDS.append(pd.DataFrame(dict(
        cell=label, record_id=[label + "|" + str(i) for i in range(len(d))],
        gene=d.gene_id.astype(str).values, chr=d.chr.values,
        pos=d.pos.values, p=d.p_mr.values)))

    def is_known(ch, pos):
        arr = KN.get(str(ch))
        if arr is None or not len(arr):
            return False
        i = np.searchsorted(arr, pos)
        return any(0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000
                   for j in (i - 1, i))

    snp = d.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp)
    d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
    d["known"] = [is_known(c, p) for c, p in zip(d.chr, d.pos)]
    bg = d.groupby("locus").agg(known=("known", "any"))
    BT, BK = len(bg), int(bg.known.sum())
    sig = d[d.fdr < .05]
    sg = sig.groupby("locus").agg(known=("known", "any"))
    ST, SK = len(sg), int(sg.known.sum())
    fold = ((SK / ST) / (BK / BT)) if ST and BK else np.nan
    p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)) if ST else np.nan
    print(f"  [{label}] records {len(d):,}; background known {BK}/{BT} = {100*BK/BT:.1f}%; "
          f"FDR<0.05 {len(sig)} records / {ST} loci, {SK} known "
          f"({100*SK/ST if ST else float('nan'):.1f}%) -> {fold:.2f}x, P = {p:.3g}")
    return dict(exposure="eQTLGen_blood", disease=label, bg_loci=BT, bg_known=BK,
                sig_records=len(sig), sig_loci=ST, sig_known=SK,
                pct_known=round(100 * SK / ST, 1) if ST else np.nan,
                fold=round(fold, 2) if ST else np.nan, fisher_p=p,
                nominal_hits=int((d.p_mr < .05).sum()))


out = []
g = scan(f"{MR}/hcc/finngen_R12_C3_HEPATOCELLU_CARC_EXALLC.gz",
         dict(chr="#chrom", pos="pos", beta="beta", se="sebeta"), "HCC_low")
out.append(cell(g, "HCC_low"))
g = scan(f"{MR}/hcc/GCST90809296.h.tsv.gz",
         dict(chr="chromosome", pos="base_pair_location", beta="beta",
              se="standard_error"), "HCC_high")
out.append(cell(g, "HCC_high"))
pd.DataFrame(out).to_csv(f"{MR}/94f_eqtlgen_hcc.tsv", sep="\t", index=False)
rec = pd.concat(RECORDS, ignore_index=True)
rec.to_csv(f"{MR}/123c_eqtlgen_hcc_records.tsv.gz", sep="\t",
           index=False, compression="gzip")
print("wrote 123c_eqtlgen_hcc_records.tsv.gz: " + format(len(rec), ",") +
      " records across " + str(rec.cell.nunique()) + " cells")
print("\nwrote 94f")
