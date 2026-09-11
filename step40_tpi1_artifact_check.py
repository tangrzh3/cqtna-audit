#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 40  Is our TPI1 cis-eQTL a real signal, or an artefact of including a gene
         that the canonical pipeline excludes?

THE CONCERN (Step 39b)
TPI1 is not quantified in ANY stimulated T-cell eQTL dataset in the eQTL
Catalogue (8/8 missing, across Schmiedel_2018, Nathan_2022 and Cytoimmgen),
while SPSB2 is present in 7/8. Cytoimmgen is the SAME experiment as our
exposure data (identical 16 h sample sizes: Naive 99, Memory 95), so the
canonical processing of our own source excludes the gene our whole application
example rests on. TPI1 has processed pseudogenes (TPI1P1/P2/P3), and low
mappability is the obvious candidate reason.

CHECKS
  A  signal shape -- a real promoter cis-eQTL should peak at the TSS and decay;
     a mapping artefact tends to be diffuse or displaced
  B  locus specificity -- does the same variant look like an eQTL for many
     neighbouring genes (hub or artefact) or mainly TPI1?
  C  pseudogenes -- are TPI1P* present in our data, and do they carry a
     correlated signal at the same variants?
  D  *** pipeline validation *** -- for SPSB2, which IS in both, do our numbers
     agree with the catalogue's Cytoimmgen values? If they do, our pipeline is
     sound and TPI1's absence there is their filtering choice, not our error.
"""
import sys
import glob
import json
import os
import time
import urllib.request
import urllib.parse
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

PARQ = r"D:/Downloads/CD4_eqtl_step1_clean"
MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
API = "https://www.ebi.ac.uk/eqtl/api/v2"
UA = {"User-Agent": "Mozilla/5.0"}

TPI1, SPSB2 = "ENSG00000111669", "ENSG00000111671"
INSTR_TPI1 = 6867132
NAIVE16 = os.path.join(PARQ, "CD4_Naive_stim_16h_step1_clean.parquet")


def api(path, **p):
    url = f"{API}/{path}?" + urllib.parse.urlencode(p)
    for a in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                        timeout=120) as r:
                return json.load(r)
        except Exception:
            if a == 4:
                return None
            time.sleep(2 * (a + 1))


# ======================================================== A. signal shape
print("=" * 82)
print("A  TPI1 cis-eQTL signal shape, Naive 16h")
print("=" * 82)
t = pq.read_table(NAIVE16, columns=["gene_id", "variant_id", "chr", "pos", "beta",
                                    "se", "pval", "eaf", "ma_count", "tss_distance"],
                  filters=[("gene_id", "==", TPI1)]).to_pandas()
t["nlp"] = -np.log10(t.pval.clip(lower=1e-320))
print(f"  variants in cis window: {len(t):,}")
print(f"  p<5e-8: {(t.pval < 5e-8).sum()}   p<1e-5: {(t.pval < 1e-5).sum()}")
top = t.loc[t.pval.idxmin()]
print(f"  top: chr{top.chr}:{top.pos}  p={top.pval:.3g}  beta={top.beta:+.3f}  "
      f"eaf={top.eaf:.3f}  ma_count={top.ma_count}  tss_dist={top.tss_distance}")

sig = t[t.pval < 5e-8]
if len(sig):
    print(f"  genome-wide significant variants span "
          f"{sig.pos.min():,}-{sig.pos.max():,} ({(sig.pos.max()-sig.pos.min())/1000:.0f} kb)")
    print(f"  their tss_distance range: {sig.tss_distance.min():,} to "
          f"{sig.tss_distance.max():,}")
print("\n  -log10 p by distance from TSS (binned):")
t["bin"] = pd.cut(t.tss_distance, [-1e6, -2e5, -5e4, -1e4, -2e3, 2e3, 1e4, 5e4, 2e5, 1e6])
print(t.groupby("bin", observed=True).agg(n=("nlp", "size"), max_nlp=("nlp", "max"),
                                          median_nlp=("nlp", "median")
                                          ).to_string(float_format=lambda v: f"{v:.2f}"))

# ======================================================== B. locus specificity
print("\n" + "=" * 82)
print("B  what else does the TPI1 instrument associate with, in our own data")
print("=" * 82)
allg = pq.read_table(NAIVE16, columns=["gene_id", "pos", "pval", "beta"],
                     filters=[("pos", "==", INSTR_TPI1)]).to_pandas()
allg = allg.sort_values("pval")
print(f"  genes tested against chr12:{INSTR_TPI1} in Naive 16h: {len(allg)}")
print(allg.head(12).to_string(index=False, float_format=lambda v: f"{v:.3g}"))

# ======================================================== C. pseudogenes
print("\n" + "=" * 82)
print("C  TPI1 pseudogenes in our data")
print("=" * 82)
PSEUDO = {"TPI1P1": "ENSG00000229314", "TPI1P2": "ENSG00000232112",
          "TPI1P3": "ENSG00000226232"}
ids = pq.read_table(NAIVE16, columns=["gene_id"]).to_pandas().gene_id.unique()
print(f"  genes in this profile: {len(ids):,}")
for nm, e in PSEUDO.items():
    print(f"    {nm} ({e}): {'PRESENT' if e in set(ids) else 'absent'}")
near = pq.read_table(NAIVE16, columns=["gene_id", "pos", "pval"],
                     filters=[("pos", "==", INSTR_TPI1)]).to_pandas()
print(f"  (the instrument is tested against {len(near)} genes; a multi-mapping "
      f"artefact would typically act in trans, not as a clean cis peak)")

# ======================================================== D. pipeline check
print("\n" + "=" * 82)
print("D  pipeline validation on SPSB2 -- ours vs the catalogue's Cytoimmgen")
print("=" * 82)
ds = api("datasets/", size=1000)
cy = [d for d in ds if d["study_label"] == "Cytoimmgen"
      and d["sample_group"] == "combined_CD4_Naive_STIM_16H"
      and d["quant_method"] == "ge"]
if not cy:
    print("  Cytoimmgen Naive 16h not found")
else:
    dsid = cy[0]["dataset_id"]
    print(f"  catalogue dataset {dsid} (n={cy[0]['sample_size']}) vs our "
          f"CD4_Naive_stim_16h")
    ours = pq.read_table(NAIVE16,
                         columns=["gene_id", "pos", "beta", "se", "pval", "eaf"],
                         filters=[("gene_id", "==", SPSB2)]).to_pandas()
    rows = []
    for rsid, pos in (("rs5446", 6847298), ("rs2071069", 6869846),
                      ("rs12831467", 6886522), ("rs12302749", 6867132)):
        a = api(f"datasets/{dsid}/associations", rsid=rsid, size=200)
        time.sleep(0.4)
        cat = None
        if a:
            df = pd.DataFrame(a)
            hit = df[df.molecular_trait_id == SPSB2]
            if len(hit):
                cat = hit.iloc[0]
        o = ours[ours.pos == pos]
        rows.append(dict(rsid=rsid, pos=pos,
                         our_p=float(o.pval.iloc[0]) if len(o) else None,
                         our_beta=float(o.beta.iloc[0]) if len(o) else None,
                         our_eaf=float(o.eaf.iloc[0]) if len(o) else None,
                         cat_p=float(cat.pvalue) if cat is not None else None,
                         cat_beta=float(cat.beta) if cat is not None else None,
                         cat_maf=float(cat.maf) if cat is not None else None))
    cmp = pd.DataFrame(rows)
    print(cmp.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    cmp.to_csv(f"{MR}/40a_pipeline_validation_SPSB2.tsv", sep="\t", index=False)
    ok = cmp.dropna(subset=["our_beta", "cat_beta"])
    if len(ok):
        same_dir = (np.sign(ok.our_beta) == np.sign(ok.cat_beta)).sum()
        print(f"\n  direction agrees in {same_dir}/{len(ok)} variants")
        print("  -> if betas and p values track closely, our pipeline reproduces "
              "the canonical one\n     and TPI1's absence there is their gene "
              "filter, not our error")
print("\nwritten: 40a")
