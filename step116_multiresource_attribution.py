#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 116 -- S35: locus attribution across three further exposure resources.

Executes manuscript/PREREG_multiresource_attribution.md with the corrections
registered in its §9bis. The outcome is byte-for-byte the file used throughout
the paper; the only thing that changes between cells is where the instruments
come from.

Cells (S35 §0, as corrected):
  Randolph_2021  QTD000588  CD4 T, influenza 6 h        n=89   STIMULATED
  Schmiedel_2018 QTD000484  CD4 T, anti-CD3/CD28 4 h    n=89   STIMULATED
  Nathan_2022    QTD000666  CD4+ activated (cell state) n=147  CONTRAST, unstimulated

Design points fixed in the pre-registration and not chosen here:
  * one instrument per gene = the lead cis-eQTL, then nominal P < 5e-8
  * known-locus reference = the same Landi 2020 list, 1 Mb window
  * counting unit = independent locus (1 Mb single linkage)
  * enrichment = one-sided Fisher against that resource's own instrument loci
  * gates: >= 200 instrument records, >= 5 significant loci, else the cell is
    reported as underpowered and does NOT enter the omnibus

Note on matching: the eQTL Catalogue is GRCh38, as is the outcome, so variants
match on chr:pos with allele consistency -- the same code path as the Soskic
exposure, not the rsID workaround step92 needed for GRCh37 eQTLGen.

Note on the exposure beta: the permuted files carry beta but no standard error.
Nothing downstream needs one: under a single instrument the Wald z is
beta_out/se_out, which is this paper's own central identity, and F is recovered
from the nominal p-value as a chi-square deviate. Both are recorded per cell.

Outputs: 116a_instruments_<res>.tsv, 116b_mr_<res>.tsv, 116c_attribution.tsv,
         116_console.log
"""
import os
import csv
import gzip
import io
import sys
import urllib.request
from math import exp, lgamma

import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
OUTCOME = f"{MR}/meta_melanoma_final.tsv.gz"
KNOWN = f"{MR}/landi2020_known_loci_grch38.csv"
FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats"
UA = {"User-Agent": "Mozilla/5.0 (s35)"}

P_EXP = 5e-8
LOCUS_KB = 1000
KNOWN_KB = 1000
MIN_RECORDS = 200      # S35 §3.2 gate 2
MIN_SIG_LOCI = 5       # S35 §3.2 gate 3

CELLS = [
    ("Randolph_2021", "QTS000036", "QTD000588", "stimulated"),
    ("Schmiedel_2018", "QTS000026", "QTD000484", "stimulated"),
    ("Nathan_2022", "QTS000040", "QTD000666", "contrast_unstimulated"),
]


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


sys.stdout = Tee(f"{MR}/116_console.log")


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
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
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


print("=" * 78)
print("Step 116 -- S35 multi-resource locus attribution")
print("=" * 78)

# ---------------------------------------------- 1. instruments per resource
inst_by_res = {}
for res, qts, qtd, role in CELLS:
    url = f"{FTP}/{qts}/{qtd}/{qtd}.permuted.tsv.gz"
    raw = urllib.request.urlopen(
        urllib.request.Request(url, headers=UA), timeout=120).read()
    txt = gzip.decompress(raw).decode("utf-8", "replace")
    d = pd.read_csv(io.StringIO(txt), sep="\t")
    n_genes = len(d)
    d = d[pd.to_numeric(d.pvalue, errors="coerce") < P_EXP].copy()

    # variant is chr<CH>_<POS>_<REF>_<ALT> on GRCh38
    parts = d.variant.astype(str).str.split("_")
    d["chr"] = parts.str[0].str.replace("chr", "", regex=False)
    d["pos"] = pd.to_numeric(parts.str[1], errors="coerce")
    d["ref"] = parts.str[2].str.upper()
    d["alt"] = parts.str[3].str.upper()
    d = d.dropna(subset=["pos"])
    d["pos"] = d["pos"].astype(int)
    d["F"] = chi2.isf(d.pvalue.astype(float), 1)   # z^2 from the nominal p
    d = d.rename(columns={"molecular_trait_id": "gene_id",
                          "pvalue": "pval_exp", "beta": "beta_exp"})
    d = d[["gene_id", "chr", "pos", "ref", "alt", "pval_exp", "beta_exp", "F"]]
    d["resource"] = res
    inst_by_res[res] = d
    d.to_csv(f"{MR}/116a_instruments_{res}.tsv", sep="\t", index=False)
    print(f"{res:<16} genes tested {n_genes:>6}   "
          f"instruments at P<5e-8 {len(d):>5}   min F {d.F.min():.1f}"
          f"   [{role}]")

want = set()
for d in inst_by_res.values():
    want |= set(zip(d.chr.astype(str), d.pos.astype(int)))
print(f"\nunique instrument positions to look up: {len(want):,}")

# ---------------------------------------------- 2. one pass over the outcome
got = {}
with gzip.open(OUTCOME, "rt") as f:
    rd = csv.reader(f, delimiter="\t")
    hdr = next(rd)
    ix = {c: i for i, c in enumerate(hdr)}
    for k, row in enumerate(rd):
        try:
            key = (row[ix["#chrom"]], int(row[ix["pos"]]))
        except (ValueError, IndexError):
            continue
        if key in want and key not in got:
            try:
                got[key] = (row[ix["ref"]].upper(), row[ix["alt"]].upper(),
                            float(row[ix["beta"]]), float(row[ix["sebeta"]]),
                            float(row[ix["af_alt"]]))
            except ValueError:
                pass
        if k % 5_000_000 == 0:
            print(f"  outcome rows {k:,}, matched {len(got):,}", end="\r")
print(f"\noutcome positions matched: {len(got):,} of {len(want):,}")

# ---------------------------------------------- 3. known loci
kn_raw = pd.read_csv(KNOWN)
kn = {str(c): np.sort(s.pos.values) for c, s in kn_raw.groupby("chr")}
print(f"known melanoma loci: {len(kn_raw)}")


def is_known(ch, pos):
    arr = kn.get(str(ch))
    if arr is None or not len(arr):
        return False
    i = np.searchsorted(arr, pos)
    return any(0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000
               for j in (i - 1, i))


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
summary = []

for res, qts, qtd, role in CELLS:
    d = inst_by_res[res]
    rows, drop = [], {"nomatch": 0, "harm": 0, "palin": 0}
    for r in d.itertuples():
        m = got.get((str(r.chr), int(r.pos)))
        if m is None:
            drop["nomatch"] += 1
            continue
        ref, alt, bo, so, af = m
        ea, oa = r.alt, r.ref            # catalogue effect allele is alt
        if {ea, oa} != {ref, alt}:
            drop["harm"] += 1
            continue
        if COMP.get(ea) == oa:
            drop["palin"] += 1
            continue
        sign = 1.0 if ea == alt else -1.0
        be = float(r.beta_exp) * sign
        if be == 0:
            drop["harm"] += 1
            continue
        z = bo / so                       # the paper's identity
        p = 2 * norm.sf(abs(z))
        rows.append(dict(gene_id=r.gene_id, chr=str(r.chr), pos=int(r.pos),
                         pval_exp=r.pval_exp, beta_exp=be, F=r.F,
                         beta_out=bo, se_out=so, eaf=af,
                         b_mr=bo / be, se_mr=so / abs(be), z=z, p_mr=p))
    m = pd.DataFrame(rows)
    print("\n" + "=" * 78)
    print(f"{res}  [{role}]")
    print("=" * 78)
    print(f"  harmonised records: {len(m)}   dropped "
          f"nomatch={drop['nomatch']} allele={drop['harm']} palindromic={drop['palin']}")

    if len(m) == 0:
        summary.append(dict(resource=res, role=role, n_records=0,
                            gate_records=False, verdict="no records"))
        continue

    m["fdr"] = bh(m.p_mr)
    snp = m.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp)
    m["locus"] = [loci[(c, p)] for c, p in zip(m.chr, m.pos)]
    m["known"] = [is_known(c, p) for c, p in zip(m.chr, m.pos)]
    m.to_csv(f"{MR}/116b_mr_{res}.tsv", sep="\t", index=False)

    bg = m.groupby("locus").agg(known=("known", "any"))
    BG_T, BG_K = len(bg), int(bg.known.sum())
    sig = m[m.fdr < 0.05]
    sg = sig.groupby("locus").agg(known=("known", "any"))
    S_T, S_K = len(sg), int(sg.known.sum())

    gate_rec = len(m) >= MIN_RECORDS
    gate_loci = S_T >= MIN_SIG_LOCI
    print(f"  background loci {BG_T} ({BG_K} known, {100*BG_K/BG_T:.1f}%)")
    print(f"  significant records {len(sig)}   significant loci {S_T} "
          f"({S_K} known)")
    print(f"  gate records >= {MIN_RECORDS}: {gate_rec}    "
          f"gate significant loci >= {MIN_SIG_LOCI}: {gate_loci}")

    fold = pv = None
    if gate_rec and gate_loci and S_T:
        share = S_K / S_T
        bgshare = BG_K / BG_T
        fold = share / bgshare if bgshare else float("nan")
        pv = fisher_greater(S_K, S_T - S_K, BG_K - S_K, BG_T - S_T - (BG_K - S_K))
        print(f"  known share {100*share:.1f}% vs background {100*bgshare:.1f}%"
              f"  -> fold {fold:.2f}   one-sided Fisher P {pv:.4g}")
        verdict = "in omnibus"
    else:
        verdict = "UNDERPOWERED -- direction only, not in omnibus"
        print(f"  -> {verdict}")

    summary.append(dict(resource=res, role=role, n_records=len(m),
                        bg_loci=BG_T, bg_known=BG_K, sig_loci=S_T,
                        sig_known=S_K, fold=fold, fisher_p=pv,
                        gate_records=gate_rec, gate_sig_loci=gate_loci,
                        verdict=verdict))

s = pd.DataFrame(summary)
s.to_csv(f"{MR}/116c_attribution.tsv", sep="\t", index=False)
print("\n" + "=" * 78)
print(s.to_string(index=False))
print("=" * 78)
print("\nwrote 116a/116b/116c. Omnibus (CMH) is step 117 and runs only over the")
print("cells that cleared both gates.")
