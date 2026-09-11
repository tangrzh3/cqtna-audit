#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 42  Finding 7, demonstrated positively: mRNA-based lineage calls are the
         weak link, and a score-based split inside a population also splits it
         by cell-type purity.

DOGMA-seq (GSE166188) measures 210 surface proteins alongside RNA in the same
nuclei, 16 h anti-CD3/CD28, two lysis buffers as replicates. Protein gives a
lineage call that does not depend on transcript dropout, so it can serve as the
reference against which mRNA-based calls are judged.

STEPS
  1. how well does CD4/CD8 mRNA detection track surface protein in the same cell?
  2. take protein-defined CD4 cells; split them by a gene score exactly as in
     Step 20; does the split create a gradient in PROTEIN CD8 -- i.e. is the
     purity artefact real even when purity is measured independently?
  3. repeat with mRNA-defined CD4 cells; is the gradient larger?

Output: 42a-42c
"""
import sys
import gzip
import os
import numpy as np
import pandas as pd

D = r"D:/Downloads/GSE166188"
MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")

SAMPLES = {
    "LLL_stim": dict(
        gex_prefix="GSM5065528_LLL_STIM_GExp_ATAC_filtered",
        adt_prefix="GSM5065529_LLL_stim_ADT_allCounts"),
    "DIG_stim": dict(
        gex_prefix="GSM5065534_DIG_STIM_GExp_ATAC_filtered",
        adt_prefix="GSM5065535_DIG_stim_ADT_allCounts"),
}
GLYCO = ["SLC2A1", "SLC2A3", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
         "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3"]


def read_lines(p):
    with gzip.open(p, "rt") as fh:
        return [l.rstrip("\n") for l in fh]


def read_mtx(path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                nr, nc, nnz = (int(x) for x in line.split())
                break
    m = pd.read_csv(path, sep=r"\s+", comment="%", skiprows=1, header=None,
                    names=["r", "c", "v"],
                    dtype={"r": np.int32, "c": np.int32, "v": np.int32}, engine="c")
    if len(m) == nnz + 1:
        m = m.iloc[1:]
    return m, nr, nc


def norm_bc(b):
    return b.split("-")[0].strip()


out_rows, split_rows = [], []
for name, cfg in SAMPLES.items():
    print("=" * 76)
    print(name)
    print("=" * 76)

    # ---- protein ---------------------------------------------------------
    prot = read_lines(os.path.join(D, cfg["adt_prefix"] + ".proteins.txt.gz"))
    abc = [norm_bc(b) for b in read_lines(os.path.join(D, cfg["adt_prefix"] + ".barcodes.txt.gz"))]
    am, anr, anc = read_mtx(os.path.join(D, cfg["adt_prefix"] + ".mtx.gz"))
    print(f"  ADT: {len(prot)} proteins x {len(abc):,} barcodes (mtx {anr}x{anc})")
    pidx = {p: i + 1 for i, p in enumerate(prot)}
    # exact tag names -- a prefix match would also catch CD40/CD44/CD45/CD47/CD48
    # and CD80/CD81/CD83/CD86, which are unrelated markers
    CD4_TAGS = {"CD4", "CD4-1", "CD4-2"}
    CD8_TAGS = {"CD8", "CD8a", "CD8A", "CD8-1", "CD8-2", "CD8b"}
    cd4p = [p for p in prot if p in CD4_TAGS]
    cd8p = [p for p in prot if p in CD8_TAGS]
    print(f"    CD4 tags: {cd4p}   CD8 tags: {cd8p}")
    if not cd4p or not cd8p:
        print("    *** CD4 or CD8 antibody not found; cannot proceed ***")
        continue

    def adt_vec(tag):
        r = pidx[tag]
        s = am[am.r == r].set_index("c").v
        return s.reindex(range(1, len(abc) + 1), fill_value=0).values

    tot_adt = am.groupby("c").v.sum().reindex(range(1, len(abc) + 1), fill_value=0).values
    p_cd4 = np.sum([adt_vec(t) for t in cd4p], axis=0)
    p_cd8 = np.sum([adt_vec(t) for t in cd8p], axis=0)
    clr = lambda x: np.log1p(1e4 * x / np.maximum(tot_adt, 1))
    prot_df = pd.DataFrame(dict(bc=abc, adt_total=tot_adt,
                                CD4_prot=clr(p_cd4), CD8_prot=clr(p_cd8)))

    # ---- RNA -------------------------------------------------------------
    feats = [l.split("\t") for l in read_lines(os.path.join(D, cfg["gex_prefix"] + "_features.tsv.gz"))]
    gbc = [norm_bc(b) for b in read_lines(os.path.join(D, cfg["gex_prefix"] + "_barcodes.tsv.gz"))]
    gm, gnr, gnc = read_mtx(os.path.join(D, cfg["gex_prefix"] + "_matrix.mtx.gz"))
    # features columns: 0 = Ensembl id, 1 = SYMBOL, 2 = kind. Use the symbol.
    kinds = pd.Series([f[2] if len(f) > 2 else "?" for f in feats])
    names = pd.Series([f[1] if len(f) > 1 else f[0] for f in feats])
    gene_rows = {}
    for i, (n, k) in enumerate(zip(names, kinds), start=1):
        if k == "Gene Expression":
            gene_rows.setdefault(n, i)
    print(f"  GEX: {gnr:,} features x {len(gbc):,} barcodes "
          f"({(kinds=='Gene Expression').sum():,} genes)")

    gene_row_set = set(gene_rows.values())
    rna_tot = gm[gm.r.isin(gene_row_set)].groupby("c").v.sum()
    rna_tot = rna_tot.reindex(range(1, len(gbc) + 1), fill_value=0).values

    def rna_vec(g):
        r = gene_rows.get(g)
        if r is None:
            return np.zeros(len(gbc))
        s = gm[gm.r == r].set_index("c").v
        return s.reindex(range(1, len(gbc) + 1), fill_value=0).values

    rna_df = pd.DataFrame(dict(bc=gbc, rna_total=rna_tot,
                               CD4_rna=rna_vec("CD4"), CD8A_rna=rna_vec("CD8A"),
                               CD8B_rna=rna_vec("CD8B")))
    glyco_present = [g for g in GLYCO if g in gene_rows]
    gl = np.vstack([rna_vec(g) for g in glyco_present])
    cp = 1e4 * gl / np.maximum(rna_tot, 1)
    z = (cp - cp.mean(axis=1, keepdims=True)) / np.maximum(cp.std(axis=1, keepdims=True), 1e-9)
    rna_df["Glyco"] = z.mean(axis=0)
    rna_df["nFeature"] = gm[gm.r.isin(gene_row_set)].groupby("c").size().reindex(
        range(1, len(gbc) + 1), fill_value=0).values

    # ---- join ------------------------------------------------------------
    d = rna_df.merge(prot_df, on="bc", how="inner")
    d = d[(d.rna_total >= 200) & (d.adt_total >= 100)]
    print(f"  cells with both modalities after QC: {len(d):,}")
    if len(d) < 200:
        print("    too few joint cells; skipping")
        continue

    # ---- 1. concordance --------------------------------------------------
    d["prot_CD4"] = d.CD4_prot > d.CD4_prot.median()
    d["prot_CD8"] = d.CD8_prot > d.CD8_prot.median()
    t_cells = d[d.prot_CD4 != d.prot_CD8]
    pc4 = t_cells[t_cells.prot_CD4]
    pc8 = t_cells[t_cells.prot_CD8]
    print(f"\n  protein-defined: CD4 {len(pc4):,}  CD8 {len(pc8):,}")
    print(f"    within protein CD4 cells:  CD4 mRNA detected "
          f"{100*np.mean(pc4.CD4_rna>0):.1f}%   CD8A mRNA detected "
          f"{100*np.mean(pc4.CD8A_rna>0):.1f}%")
    print(f"    within protein CD8 cells:  CD4 mRNA detected "
          f"{100*np.mean(pc8.CD4_rna>0):.1f}%   CD8A mRNA detected "
          f"{100*np.mean(pc8.CD8A_rna>0):.1f}%")
    out_rows.append(dict(sample=name, n_joint=len(d),
                         n_prot_CD4=len(pc4), n_prot_CD8=len(pc8),
                         CD4mRNA_in_protCD4=100*np.mean(pc4.CD4_rna > 0),
                         CD8AmRNA_in_protCD4=100*np.mean(pc4.CD8A_rna > 0),
                         CD4mRNA_in_protCD8=100*np.mean(pc8.CD4_rna > 0),
                         CD8AmRNA_in_protCD8=100*np.mean(pc8.CD8A_rna > 0)))

    # ---- 2/3. does a score split create a purity gradient? ---------------
    for label, pop in (("protein-defined CD4", pc4),
                       ("mRNA-defined CD4", d[(d.CD4_rna > 0) &
                                              (d.CD8A_rna == 0) & (d.CD8B_rna == 0)])):
        if len(pop) < 100:
            print(f"\n  {label}: only {len(pop)} cells, skipped")
            continue
        q40, q60 = pop.Glyco.quantile([.4, .6])
        hi = pop[pop.Glyco > q60]
        lo = pop[pop.Glyco <= q40]
        dc8 = hi.CD8_prot.mean() - lo.CD8_prot.mean()
        pooled = np.sqrt((hi.CD8_prot.var() + lo.CD8_prot.var()) / 2)
        smd = dc8 / pooled if pooled > 0 else np.nan
        print(f"\n  {label} (n={len(pop):,}) split by glycolysis score:")
        print(f"    hi {len(hi):,} / lo {len(lo):,}")
        print(f"    CD8 PROTEIN  hi {hi.CD8_prot.mean():.3f} vs lo "
              f"{lo.CD8_prot.mean():.3f}   SMD = {smd:+.3f}")
        print(f"    CD4 PROTEIN  hi {hi.CD4_prot.mean():.3f} vs lo "
              f"{lo.CD4_prot.mean():.3f}")
        print(f"    nFeature     hi {hi.nFeature.mean():.0f} vs lo "
              f"{lo.nFeature.mean():.0f}")
        split_rows.append(dict(sample=name, population=label, n=len(pop),
                               n_hi=len(hi), n_lo=len(lo),
                               CD8prot_hi=hi.CD8_prot.mean(),
                               CD8prot_lo=lo.CD8_prot.mean(),
                               CD8prot_SMD=smd,
                               CD4prot_hi=hi.CD4_prot.mean(),
                               CD4prot_lo=lo.CD4_prot.mean(),
                               nFeature_hi=hi.nFeature.mean(),
                               nFeature_lo=lo.nFeature.mean()))

pd.DataFrame(out_rows).to_csv(f"{MR}/42a_protein_vs_mRNA_concordance.tsv",
                              sep="\t", index=False)
pd.DataFrame(split_rows).to_csv(f"{MR}/42b_score_split_purity_gradient.tsv",
                                sep="\t", index=False)
print("\nwritten: 42a, 42b")
