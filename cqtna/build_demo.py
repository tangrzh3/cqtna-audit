#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build the CQTNA demo inputs from this audit's own tables.

The demo is the paper auditing itself. That is the point: a tool for detecting
selective emphasis should be run first on the study proposing it, and its output
should be checkable against what the paper reports.

Writes cqtna/demo/.
"""
import json
import os

import pandas as pd

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "cqtna", "demo")
NOVEL = "潜在新位点"
os.makedirs(OUT, exist_ok=True)


def mr_table(path, out, p_col, gene_col):
    d = pd.read_csv(path, sep="\t")
    d = d[d[p_col].notna()].copy()
    d[["chr", "pos"]] = d.SNP.str.split(":", expand=True).iloc[:, :2]
    t = pd.DataFrame(dict(
        record_id=[f"r{i}" for i in range(1, len(d) + 1)],
        gene=d[gene_col].fillna(d.gene_id).values,
        chr=d.chr.values, pos=d.pos.astype(int).values,
        p=d[p_col].values,
        exposure_profile=d.exposure.values if "exposure" in d.columns else ""))
    t.to_csv(os.path.join(OUT, out), sep="\t", index=False)
    print(f"  {out}: {len(t):,} records")


# ---- the two outcome rounds: meta, and the FinnGen-only round
mr_table(f"{MR}/13_meta_locus_annotation.tsv", "mr_meta.tsv", "pval", "SYMBOL")
mr_table(f"{MR}/06_locus_annotation.tsv", "mr_finngen.tsv", "pval", "SYMBOL")

# ---- known-locus lists: the outcome's own, and a mismatched one
mel = pd.read_csv(f"{MR}/landi2020_known_loci_grch38.csv")
mel[["chr", "pos"]].assign(rsid=mel.rsid).to_csv(
    os.path.join(OUT, "known_melanoma.tsv"), sep="\t", index=False)
hcc = pd.read_csv(f"{MR}/84a_hcc_known_loci_grch38.csv")
hcc[["chr", "pos"]].assign(rsid=hcc.rsid).to_csv(
    os.path.join(OUT, "known_hcc_mismatched.tsv"), sep="\t", index=False)
print(f"  known_melanoma.tsv: {len(mel)} loci")
print(f"  known_hcc_mismatched.tsv: {len(hcc)} loci (negative control)")

# ---- instrument attrition over the glycolytic pathway
GLYCO28 = ["SLC2A1", "SLC2A3", "SLC2A4", "HK1", "HK2", "HK3", "GPI", "PFKL",
           "PFKM", "PFKP", "PFKFB1", "PFKFB2", "PFKFB3", "PFKFB4", "ALDOA",
           "ALDOB", "ALDOC", "TPI1", "GAPDH", "GAPDHS", "PGK1", "PGK2", "PGAM1",
           "PGAM2", "ENO1", "ENO2", "PKM", "LDHA"]
inst = pd.DataFrame(dict(gene=GLYCO28))
inst["instrumentable"] = inst.gene.isin(["TPI1", "ENO1", "SLC2A1"]).astype(int)
inst["analysable"] = inst.gene.isin(["TPI1", "ENO1"]).astype(int)
inst["associated"] = inst.gene.isin(["TPI1"]).astype(int)
inst.to_csv(os.path.join(OUT, "instruments.tsv"), sep="\t", index=False)
print(f"  instruments.tsv: {len(inst)} pathway genes, "
      f"{inst.instrumentable.sum()}/{inst.analysable.sum()}/{inst.associated.sum()}")

# ---- compartment expression for the nominated gene
expr = pd.DataFrame([
    dict(gene="TPI1", cell_type="CD4_T", mean_expression=1.00),
    dict(gene="TPI1", cell_type="Malignant", mean_expression=6.70),
    dict(gene="TPI1", cell_type="CD8_T", mean_expression=1.15),
    dict(gene="TPI1", cell_type="Macrophage", mean_expression=2.10),
])
expr.to_csv(os.path.join(OUT, "expression.tsv"), sep="\t", index=False)
print("  expression.tsv: compartment ratios for the nominated gene")

cfg = dict(mr_results="mr_meta.tsv",
           mr_results_alt="mr_finngen.tsv",
           known_loci="known_melanoma.tsv",
           known_loci_mismatched="known_hcc_mismatched.tsv",
           instrument_table="instruments.tsv",
           expression_table="expression.tsv",
           target_cell_type="CD4_T",
           fdr_threshold=0.05, locus_window_kb=1000, known_window_kb=1000)
json.dump(cfg, open(os.path.join(OUT, "demo_config.json"), "w", encoding="utf-8"),
          indent=1)
print("  demo_config.json")
