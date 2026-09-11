"""Step 92 -- swap the EXPOSURE resource, hold the outcome fixed.

Executes manuscript/PREREG_exposure_resource.md. The outcome file is byte-for-byte
the one used throughout the paper (meta_melanoma_final.tsv.gz); the only thing
changed is where the instruments come from: eQTLGen whole blood (n = 31,684)
instead of Soskic CD4+ T cells (n = 85-100).

Design points fixed in the pre-registration and not chosen here:
  * one instrument per gene = the cis-eQTL with the smallest P, then P < 5e-8
  * known-locus reference = the same Landi 2020 list, 1 Mb window
  * counting unit = independent locus (1 Mb single linkage), not record
  * enrichment = one-sided Fisher against all instrument loci, PLUS the
    decile x quintile matched-background version of step85e
  * the primary quantity is the SHARE of significant loci that are known, not
    the count -- eQTLGen has ~300x the eQTL sample size, so more hits are
    expected and would neither support nor refute the proposition

Build note: eQTLGen 2019 is GRCh37, the outcome is GRCh38. We match on rsID and
take all positions from the outcome file, so no liftover is involved.

Effect-size note: under one instrument the Wald z is beta_out/se_out up to sign,
so the MR P value does not depend on the eQTL effect magnitude at all -- which is
this paper's own central identity. We still record an exposure beta (from Z and
the outcome allele frequency) for direction and for the OR, and flag that this is
an approximation, because nothing downstream of the P value depends on it.

Output: 92a (instruments), 92b (MR records), 92c (locus attribution), 92d (summary)
"""
import sys
import csv
import gzip
import io
import os
from itertools import islice
from math import lgamma, exp, sqrt

import numpy as np
import pandas as pd
from scipy.stats import norm

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
EQTLGEN = f"{MR}/eqtlgen/cis-eQTLsFDR0.05.txt.gz"
OUTCOME = f"{MR}/meta_melanoma_final.tsv.gz"
KNOWN = f"{MR}/landi2020_known_loci_grch38.csv"
P_EXP = 5e-8
LOCUS_KB = 1000
KNOWN_KB = 1000


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    hi = min(r1, c1)
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, hi + 1)), 1.0)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


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


def main():
    # ---------------------------------------------- 1. eQTLGen: best SNP per gene
    best = {}
    n_rows = 0
    with gzip.open(EQTLGEN, "rt") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ix = {c: i for i, c in enumerate(hdr)}
        for row in rd:
            n_rows += 1
            p = float(row[ix["Pvalue"]])
            if p >= P_EXP:
                continue
            g = row[ix["Gene"]]
            cur = best.get(g)
            if cur is None or p < cur[0]:
                best[g] = (p, row[ix["SNP"]], row[ix["AssessedAllele"]],
                           row[ix["OtherAllele"]], float(row[ix["Zscore"]]),
                           row[ix["GeneSymbol"]], int(row[ix["NrSamples"]]))
            if n_rows % 2_000_000 == 0:
                print(f"  eQTLGen rows {n_rows:,}, genes so far {len(best):,}", end="\r")
    print(f"\neQTLGen rows scanned: {n_rows:,}")
    print(f"genes with a cis-eQTL at P < {P_EXP:g}: {len(best):,}")
    inst = pd.DataFrame(
        [dict(gene_id=g, pval_exp=v[0], rsid=v[1], ea=v[2], oa=v[3], z_exp=v[4],
              symbol=v[5], n_exp=v[6]) for g, v in best.items()])
    inst.to_csv(f"{MR}/92a_eqtlgen_instruments.tsv", sep="\t", index=False)

    # ---------------------------------------------- 2. outcome lookup by rsID
    want = set(inst.rsid)
    got = {}
    with gzip.open(OUTCOME, "rt") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        ix = {c: i for i, c in enumerate(hdr)}
        for k, row in enumerate(rd):
            rs = row[ix["rsids"]]
            if not rs:
                continue
            for r in rs.split(","):
                if r in want and r not in got:
                    try:
                        got[r] = (row[ix["#chrom"]], int(row[ix["pos"]]),
                                  row[ix["ref"]], row[ix["alt"]],
                                  float(row[ix["beta"]]), float(row[ix["sebeta"]]),
                                  float(row[ix["af_alt"]]))
                    except ValueError:
                        pass
            if k % 5_000_000 == 0:
                print(f"  outcome rows {k:,}, matched {len(got):,}", end="\r")
    print(f"\noutcome SNPs matched: {len(got):,} of {len(want):,} instruments")

    # ---------------------------------------------- 3. harmonise + Wald ratio
    rows, drop = [], {"nomatch": 0, "harm": 0, "palin": 0}
    for r in inst.itertuples():
        m = got.get(r.rsid)
        if m is None:
            drop["nomatch"] += 1
            continue
        ch, pos, ref, alt, bo, so, af = m
        ea, oa = r.ea.upper(), r.oa.upper()
        if {ea, oa} == {alt.upper(), ref.upper()}:
            sign = 1.0 if ea == alt.upper() else -1.0
        else:
            drop["harm"] += 1
            continue
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        if comp.get(ea) == oa:
            drop["palin"] += 1
            continue
        z = r.z_exp * sign
        f_ = af if 0 < af < 1 else 0.5
        be = z / sqrt(2 * f_ * (1 - f_) * (r.n_exp + z * z))   # approximation, see docstring
        if be == 0:
            continue
        b_mr = bo / be
        se_mr = abs(so / be)
        zz = b_mr / se_mr
        rows.append(dict(gene_id=r.gene_id, symbol=r.symbol, rsid=r.rsid,
                         chr=str(ch), pos=pos, pval_exp=r.pval_exp, n_exp=r.n_exp,
                         beta_exp=be, beta_out=bo, se_out=so, eaf=af,
                         b_mr=b_mr, se_mr=se_mr, z=zz,
                         p_mr=2 * norm.sf(abs(zz))))
    d = pd.DataFrame(rows)
    d["fdr"] = bh(d.p_mr.values)
    d.to_csv(f"{MR}/92b_mr_eqtlgen_melanoma.tsv", sep="\t", index=False)
    print(f"MR records: {len(d):,}  dropped: {drop}")
    print(f"FDR<0.05: {(d.fdr < .05).sum():,}   nominal p<0.05: {(d.p_mr < .05).sum():,}")

    # ---------------------------------------------- 4. locus attribution
    known = pd.read_csv(KNOWN)
    known["chr"] = known["chr"].astype(str)
    kn = {c: np.sort(s.pos.values) for c, s in known.groupby("chr")}
    print(f"known melanoma/naevus/pigmentation loci: {len(known)}")

    def is_known(ch, pos):
        arr = kn.get(str(ch))
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
    BG_T, BG_K = len(bg), int(bg.known.sum())
    sig = d[d.fdr < 0.05]
    sg = sig.groupby("locus").agg(known=("known", "any"))
    S_T, S_K = len(sg), int(sg.known.sum())
    fold = (S_K / S_T) / (BG_K / BG_T) if S_T and BG_K else float("nan")
    pf = fisher_greater(S_K, S_T - S_K, BG_K - S_K, (BG_T - BG_K) - (S_T - S_K))
    print(f"\n=== E1 locus attribution (eQTLGen exposure) ===")
    print(f"  background: {BG_K}/{BG_T} = {100*BG_K/BG_T:.1f}% of instrument loci are known")
    print(f"  FDR<0.05:   {len(sig):,} records across {S_T} independent loci, "
          f"{S_K} known = {100*S_K/S_T:.1f}%  -> {fold:.2f}x, one-sided P = {pf:.3g}")
    print(f"  (Soskic exposure, same outcome: 3/7 = 42.9%, 4.09x, P = 0.028)")

    d.to_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t", index=False)
    pd.DataFrame([dict(exposure="eQTLGen_wholeblood", n_instruments=len(d),
                       bg_loci=BG_T, bg_known=BG_K, sig_records=len(sig),
                       sig_loci=S_T, sig_known=S_K, pct_known=100 * S_K / S_T,
                       fold=fold, fisher_p=pf)]).to_csv(
        f"{MR}/92d_summary.tsv", sep="\t", index=False)

    # ---------------------------------------------- 5. E2: overlap with Soskic
    mel = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
    sos_novel = set(mel[(mel.FDR < 0.05) & (mel.category == "潜在新位点")].SYMBOL.dropna())
    eq_novel = set(sig[~sig.known].symbol.dropna())
    print(f"\n=== E2 novel-locus gene overlap ===")
    print(f"  Soskic novel-locus genes:   {len(sos_novel)} {sorted(sos_novel)}")
    print(f"  eQTLGen novel-locus genes:  {len(eq_novel)}")
    print(f"  INTERSECTION: {sorted(sos_novel & eq_novel)} (n={len(sos_novel & eq_novel)})")

    # ---------------------------------------------- 6. PC1: MC1R region present?
    mc1r = d[(d.chr == "16") & (d.pos.between(89_900_000, 90_100_000))]
    print(f"\n=== PC1 MC1R region (chr16:89.9-90.1 Mb) ===")
    print(f"  instruments in region: {len(mc1r)};  min MR P = "
          f"{mc1r.p_mr.min() if len(mc1r) else float('nan'):.3g}")
    print("\nwrote 92a-92d")


if __name__ == "__main__":
    main()
