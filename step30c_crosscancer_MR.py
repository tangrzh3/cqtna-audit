#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 30c  Cross-cancer MR: the same CD4 dynamic-eQTL instruments against five
          other FinnGen R12 cancers.

TWO PURPOSES

 (1) MAIN LINE -- a second, independent test of finding 4.
     Finding 4 showed the candidate list turns over completely when the OUTCOME
     GWAS POWER changes (FinnGen 5,753 cases -> meta 12,530 cases; zero overlap).
     Here the outcome DISEASE changes instead, with the exposure side held
     identical. If lists again barely overlap, instability is demonstrated along
     a second axis. Note this is a weaker claim than (1) on its own -- different
     cancers genuinely have different biology, so non-overlap is expected to
     some degree. The informative quantities are therefore:
        - how much of each list is pigmentation/nevus loci (locus attribution
          repeated per cancer -- melanoma should be the outlier)
        - whether the SAME shared-biology genes recur, or whether each cancer
          produces its own idiosyncratic list

 (2) APPLICATION -- is TPI1's effect melanoma-specific or pan-cancer?
     Pan-cancer would support an immune-surveillance mechanism; melanoma-only
     would support specificity. Both are usable.

Allele handling reuses the harmonised melanoma alignment: effect_allele.outcome
in 01_harmonised_all.tsv is already oriented to the exposure effect allele, and
all FinnGen R12 releases share genome build and ref/alt convention, so the same
orientation applies to every endpoint.

Wald ratio only (one instrument per exposure), as in the main analysis.
Output: 36a-36d
"""
import csv
import math
import os
import collections
from statistics import NormalDist

MR = r"D:/R_ex/MR"
EXT = os.path.join(MR, "cancer_extracts")
ND = NormalDist()

CANCERS = {
    "C3_BRONCHUS_LUNG_EXALLC": "Lung",
    "C3_COLORECTAL_EXALLC": "Colorectal",
    "C3_PANCREAS_EXALLC": "Pancreas",
    "C3_BREAST_EXALLC": "Breast",
    "C3_PROSTATE_EXALLC": "Prostate",
}
CAND = {"ZFYVE19", "SMC2", "SPSB2", "TPI1", "HLA-C", "KIAA0040"}


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    adj, prev = [0.0] * n, 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        q = ps[i] * n / (rank + 1)
        prev = min(prev, q)
        adj[i] = min(prev, 1.0)
    return adj


# ---------------------------------------------------------------- exposure side
alle = {}
with open(os.path.join(MR, "01_harmonised_all.tsv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        alle[(r["exposure"], r["SNP"])] = (r["effect_allele.outcome"].upper(),
                                           r["other_allele.outcome"].upper())

expo = []
with open(os.path.join(MR, "13_meta_locus_annotation.tsv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        a = alle.get((r["exposure"], r["SNP"]))
        if a is None:
            continue
        try:
            be = float(r["beta_exposure"])
        except ValueError:
            continue
        if be == 0:
            continue
        expo.append(dict(exposure=r["exposure"], SNP=r["SNP"], symbol=r["SYMBOL"],
                         cell_type=r["cell_type"], timepoint=r["timepoint"],
                         beta_exposure=be, ea_out=a[0], oa_out=a[1],
                         category=r["category"],
                         mel_p=float(r["pval"]) if r["pval"] else None,
                         mel_FDR=float(r["FDR"]) if r["FDR"] else None,
                         mel_OR=float(r["OR"]) if r["OR"] else None))
print(f"exposures with usable instruments: {len(expo):,} "
      f"({len({e['symbol'] for e in expo}):,} genes)")

# ---------------------------------------------------------------- per cancer
rows = []
for endpoint, label in CANCERS.items():
    path = os.path.join(EXT, f"{endpoint}.tsv")
    if not os.path.exists(path):
        print(f"  [missing] {label}")
        continue
    out = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                out[r["SNP"]] = (r["ref"].upper(), r["alt"].upper(),
                                 float(r["beta"]), float(r["sebeta"]))
            except ValueError:
                continue

    recs = []
    for e in expo:
        g = out.get(e["SNP"])
        if g is None:
            continue
        ref, alt, bo, seo = g
        if e["ea_out"] == alt:
            b_out = bo
        elif e["ea_out"] == ref:
            b_out = -bo
        else:
            continue
        if seo <= 0:
            continue
        b = b_out / e["beta_exposure"]
        se = seo / abs(e["beta_exposure"])
        z = b / se
        p = 2 * ND.cdf(-abs(z)) if abs(z) < 37 else 0.0
        recs.append(dict(cancer=label, endpoint=endpoint, **{
            k: e[k] for k in ("exposure", "symbol", "cell_type", "timepoint",
                              "SNP", "category", "mel_p", "mel_FDR", "mel_OR")},
            beta=b, se=se, OR=math.exp(b), pval=p))
    fdr = bh([r["pval"] for r in recs])
    for r, q in zip(recs, fdr):
        r["FDR"] = q
    rows.extend(recs)
    n_sig = sum(1 for r in recs if r["FDR"] < 0.05)
    print(f"  {label:<12} matched {len(recs):,}/{len(expo):,} exposures | "
          f"FDR<0.05: {n_sig}")

if not rows:
    raise SystemExit("no cancer extracts found -- run step30b first")

with open(os.path.join(MR, "36a_crosscancer_MR_all.tsv"), "w", newline="",
          encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# ---------------------------------------------------------------- list turnover
mel = {r["exposure"] for r in expo if r["mel_FDR"] is not None and r["mel_FDR"] < 0.05}
mel_g = {r["symbol"] for r in expo if r["mel_FDR"] is not None and r["mel_FDR"] < 0.05}
print(f"\nmelanoma FDR<0.05: {len(mel)} records / {len(mel_g)} genes")

by = collections.defaultdict(list)
for r in rows:
    by[r["cancer"]].append(r)

print("\n" + "=" * 78)
print("LIST TURNOVER  (FDR<0.05 gene sets, melanoma vs each cancer)")
print("=" * 78)
summ = []
for c, rs in by.items():
    sig = {r["symbol"] for r in rs if r["FDR"] < 0.05}
    inter = mel_g & sig
    jac = len(inter) / len(mel_g | sig) if (mel_g | sig) else 0
    # locus attribution of this cancer's own list
    cats = collections.Counter(r["category"] for r in rs if r["FDR"] < 0.05)
    tot = sum(cats.values())
    pig = sum(v for k, v in cats.items() if k != "潜在新位点")
    summ.append(dict(cancer=c, n_genes=len(sig), overlap_with_melanoma=len(inter),
                     jaccard=round(jac, 3),
                     pct_known_pigment_or_nevus=round(100 * pig / tot, 1) if tot else 0,
                     shared_genes=",".join(sorted(inter)) if inter else "-"))
    print(f"{c:<12} genes={len(sig):<5} overlap with melanoma={len(inter):<4} "
          f"Jaccard={jac:.3f}  known-locus share={100*pig/tot if tot else 0:.0f}%")
    if inter:
        print(f"             shared: {', '.join(sorted(inter))}")

with open(os.path.join(MR, "36b_list_turnover.tsv"), "w", newline="",
          encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(summ[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(summ)

# ---------------------------------------------------------------- candidates
print("\n" + "=" * 78)
print("CANDIDATE GENES ACROSS CANCERS  (best record per gene per cancer)")
print("=" * 78)
best = {}
for r in rows:
    if r["symbol"] not in CAND:
        continue
    k = (r["symbol"], r["cancer"])
    if k not in best or r["pval"] < best[k]["pval"]:
        best[k] = r
hdr = f"{'gene':<10}" + "".join(f"{c:>22}" for c in CANCERS.values())
print(hdr); print("-" * len(hdr))
for g in sorted(CAND):
    line = f"{g:<10}"
    for c in CANCERS.values():
        r = best.get((g, c))
        if r is None:
            cellstr = "--"
        else:
            cellstr = "{:.3f} (p={:.2g})".format(r["OR"], r["pval"])
        line += "{:>22}".format(cellstr)
    print(line)

with open(os.path.join(MR, "36c_candidates_across_cancers.tsv"), "w", newline="",
          encoding="utf-8") as fo:
    vals = list(best.values())
    w = csv.DictWriter(fo, fieldnames=list(vals[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(vals)

print("\nmelanoma reference (meta):")
for g in sorted(CAND):
    e = [x for x in expo if x["symbol"] == g]
    if e:
        b = min(e, key=lambda x: x["mel_p"] if x["mel_p"] is not None else 1)
        print(f"  {g:<10} OR={b['mel_OR']:.3f}  p={b['mel_p']:.2g}  FDR={b['mel_FDR']:.3g}")

print("\nwritten: 36a-36c")
