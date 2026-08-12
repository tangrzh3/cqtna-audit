#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 38  Independent replication of the activation-window eQTL in DICE.

DICE (Schmiedel 2018) via the eQTL Catalogue: 15 sorted immune populations,
n = 82-91 donors, comparable in design and size to Soskic (~90). It has
CD4_T-cell_naive and CD4_T-cell_anti-CD3-CD28 (4 h) side by side.

Why this works where the chromatin test failed (Step 35-36): eQTLs are defined
per gene, so TPI1 and SPSB2 -- 20 kb apart and inseparable at ATAC peak
resolution -- are cleanly separable. The SPSB2 internal control is usable.

PRE-SPECIFIED, fixed before the data were seen
  D1  TPI1 cis-eQTL stronger in anti-CD3-CD28 than in naive CD4
      SPSB2 comparable in both (constitutive)          <- internal control
  D2  the variant is an eQTL in CD4 but not broadly across other cell types

DECLARED ASYMMETRY
  DICE stimulates 4 h; the Soskic window is 16 h and our time course has nothing
  between 0 h and 16 h, so TPI1's eQTL may not have risen by 4 h.
  A POSITIVE replicates. A NULL is ambiguous and must NOT be explained away as
  "4 h is too early" after the fact. Recorded before running.

FAILURE CRITERIA
  - SPSB2 also strongly activation-dependent -> internal control fails
  - TPI1 equally strong in naive -> activation dependence not replicated

METHOD NOTE
  Associations are queried by rsid. Querying by gene returns records ordered by
  POSITION and capped at 1000 per page, so a "top p" taken from one page is
  meaningless -- an earlier version of this script made that mistake.

Output: 38a, 38b
"""
import json
import time
import urllib.request
import urllib.parse
import pandas as pd

MR = r"D:/R_ex/MR"
API = "https://www.ebi.ac.uk/eqtl/api/v2"
UA = {"User-Agent": "Mozilla/5.0"}

INSTR = [("TPI1", "ENSG00000111669", "rs12302749", 6867132),
         ("SPSB2", "ENSG00000111671", "rs5446", 6847298),
         ("SPSB2", "ENSG00000111671", "rs2071069", 6869846),
         ("SPSB2", "ENSG00000111671", "rs12831467", 6886522)]


def get(path, **params):
    url = f"{API}/{path}?" + urllib.parse.urlencode(params)
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as ex:
            if attempt == 5:
                print(f"    [fail] {path} {params}: {type(ex).__name__} "
                      f"{getattr(ex,'code','')}", flush=True)
                return None
            time.sleep(3 * (attempt + 1))
    return None


ds = get("datasets/", size=1000)
dice = [d for d in ds if d["study_label"] == "Schmiedel_2018"
        and d["quant_method"] == "ge"]
print(f"DICE gene-expression datasets: {len(dice)}\n")

rows, other = [], []
for gene, ensg, rsid, pos in INSTR:
    for d in dice:
        a = get(f"datasets/{d['dataset_id']}/associations", rsid=rsid, size=200)
        time.sleep(0.4)
        if not a:
            continue
        df = pd.DataFrame(a)
        hit = df[df.molecular_trait_id == ensg]
        rec = dict(gene=gene, rsid=rsid, position=pos,
                   dataset=d["dataset_id"], cell=d["sample_group"],
                   condition=d["condition_label"], n=int(d["sample_size"]),
                   n_genes_tested=len(df))
        if len(hit):
            h = hit.iloc[0]
            rec.update(pvalue=float(h.pvalue), beta=float(h.beta),
                       se=float(h.se), maf=float(h.maf),
                       median_tpm=float(h.median_tpm) if h.median_tpm else None)
        rows.append(rec)
        # what else does this variant regulate here?
        for _, r2 in df.iterrows():
            if float(r2.pvalue) < 1e-4:
                other.append(dict(rsid=rsid, cell=d["sample_group"],
                                  condition=d["condition_label"],
                                  gene_id=r2.molecular_trait_id,
                                  pvalue=float(r2.pvalue), beta=float(r2.beta)))
        print(f"  {gene:<6}{rsid:<12}{d['sample_group']:<28}"
              f"{d['condition_label']:<18}"
              f"p={rec.get('pvalue', float('nan')):.3g}", flush=True)

res = pd.DataFrame(rows)
res.to_csv(f"{MR}/38a_DICE_instrument_eQTL.tsv", sep="\t", index=False)
pd.DataFrame(other).to_csv(f"{MR}/38b_DICE_variant_other_genes.tsv",
                           sep="\t", index=False)

print("\n" + "=" * 90)
print("D1  activation dependence in CD4  (SPSB2 = internal control)")
print("=" * 90)
cd4 = res[res.cell.str.startswith("CD4_T-cell")].copy()
piv = cd4.pivot_table(index=["gene", "rsid"], columns="condition",
                      values="pvalue", aggfunc="first")
pb = cd4.pivot_table(index=["gene", "rsid"], columns="condition",
                     values="beta", aggfunc="first")
for idx in piv.index:
    nv = piv.loc[idx].get("naive")
    st = piv.loc[idx].get("anti-CD3-CD28_4h")
    bnv = pb.loc[idx].get("naive")
    bst = pb.loc[idx].get("anti-CD3-CD28_4h")
    print(f"  {idx[0]:<7}{idx[1]:<12} naive p={nv:<12.3g} beta={bnv:+.3f}   "
          f"stim4h p={st:<12.3g} beta={bst:+.3f}")

print("\n" + "=" * 90)
print("D2  cell-type specificity")
print("=" * 90)
for gene in ("TPI1", "SPSB2"):
    g = res[res.gene == gene].sort_values("pvalue")
    print(f"\n{gene}:")
    print(g[["rsid", "cell", "condition", "n", "pvalue", "beta", "maf"]]
          .head(18).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    sig = (g.pvalue < 0.05).sum()
    print(f"  nominal p<0.05 in {sig}/{len(g)} cell-type x instrument combinations")

print("\nwritten: 38a, 38b")
