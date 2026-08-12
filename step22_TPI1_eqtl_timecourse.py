#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 22  When does the TPI1 cis-eQTL appear across CD4 activation?

Motivation: TPI1's instrument sits in Naive 16h. In Soskic's design 16h is the
PRE-DIVISION activation window -- which is also when T cells run the glycolytic
switch. If the TPI1 eQTL is specifically detectable there, the genetic
regulatory window coincides with the metabolic switch window.

Same treatment Step 13 applied to ZFYVE19 / SMC2 (where eQTLs appeared only
after the first division, at 40h).

*** WORDING DISCIPLINE (Step 13 already stepped on this) ***
A profile with no genome-wide significant eQTL has NO INSTRUMENT -- that is not
evidence that the causal effect is absent there. The only claim available is
"the dynamic data supply an instrument that resting-state data do not".

Output: 27a_TPI1_eqtl_timecourse.tsv
"""
import glob
import os
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyarrow.compute as pc

PARQ = "D:/Downloads/CD4_eqtl_step1_clean"
MR = "D:/R_ex/MR"

GENES = {
    "TPI1":    "ENSG00000111669",
    "SPSB2":   "ENSG00000111671",   # co-locus partner, for contrast
    "ZFYVE19": "ENSG00000166140",   # Step 13 reference: 40h gene
    "SMC2":    "ENSG00000136824",   # Step 13 reference: 40h gene
}
ORDER = ["CD4_Naive_uns_0h", "CD4_Naive_stim_16h", "CD4_Naive_stim_40h",
         "CD4_Naive_stim_5d", "CD4_Memory_uns_0h", "CD4_Memory_stim_16h",
         "CD4_Memory_stim_40h", "CD4_Memory_stim_5d"]
STAGE = {"uns_0h": "0h resting", "stim_16h": "16h pre-division",
         "stim_40h": "40h post-1st-division", "stim_5d": "5d effector"}

rows = []
for path in sorted(glob.glob(os.path.join(PARQ, "*.parquet"))):
    prof = os.path.basename(path).replace("_step1_clean.parquet", "")
    tbl = pq.read_table(
        path,
        columns=["gene_id", "variant_id", "chr", "pos", "beta", "se", "pval",
                 "eaf", "ma_count", "tss_distance"],
        filters=[("gene_id", "in", list(GENES.values()))])
    df = tbl.to_pandas()
    if df.empty:
        continue
    for sym, gid in GENES.items():
        d = df[df.gene_id == gid]
        if d.empty:
            continue
        top = d.loc[d.pval.idxmin()]
        stage = next((v for k, v in STAGE.items() if prof.endswith(k)), prof)
        rows.append(dict(
            gene=sym, profile=prof, stage=stage,
            n_cis_variants=len(d),
            top_variant=f"{top.chr}:{top.pos}",
            top_pval=top.pval, top_beta=top.beta, top_se=top.se,
            top_eaf=top.eaf, tss_distance=top.tss_distance,
            n_p5e8=int((d.pval < 5e-8).sum()),
            n_p1e5=int((d.pval < 1e-5).sum()),
            has_instrument=bool((d.pval < 5e-8).any())))
    print(f"{prof}: done", flush=True)

res = pd.DataFrame(rows)
res["profile"] = pd.Categorical(res.profile, categories=ORDER, ordered=True)
res = res.sort_values(["gene", "profile"])
res.to_csv(f"{MR}/27a_TPI1_eqtl_timecourse.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
for sym in GENES:
    d = res[res.gene == sym]
    if d.empty:
        continue
    print(f"\n=== {sym} ===")
    print(d[["profile", "stage", "top_variant", "top_pval", "top_beta",
             "n_p5e8", "has_instrument"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3g}"))

print("\n=== strongest cis-eQTL p by gene x stage (Naive) ===")
piv = res[res.profile.astype(str).str.contains("Naive")].pivot_table(
    index="gene", columns="stage", values="top_pval", aggfunc="min")
cols = [c for c in ["0h resting", "16h pre-division",
                    "40h post-1st-division", "5d effector"] if c in piv.columns]
print(piv[cols].to_string(float_format=lambda v: f"{v:.2e}"))

print("\nwritten: 27a_TPI1_eqtl_timecourse.tsv")
