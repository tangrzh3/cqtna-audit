#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 35  Does chromatin accessibility at glycolytic genes open on TCR activation?

GSE282266: sc-multiome (RNA + ATAC), purified naive/memory human CD4+ T cells,
TCR activation time course rest / 2.5h / 5h / 15h, 4 sets each.

PRE-SPECIFIED PREDICTIONS (fixed before looking at the data, from Step 22/29)
  1. TPI1  TSS-proximal accessibility RISES rest -> 15h
     SPSB2 TSS-proximal accessibility stays FLAT (constitutive eQTL at the same
     locus, 20 kb away)
  2. Glycolytic genes as a class gain more accessibility on activation than
     expression-matched background genes
  3. Sanity: TPI1 RNA and the glycolysis module rise on activation. If they do
     not, the premise fails and nothing downstream is interpretable.

DESIGN NOTE -- why TSS windows and not the instrument SNPs.
No instrument SNP falls inside a called peak, and the peak nearest TPI1's
instrument (92 bp away) is also the peak nearest one of SPSB2's instruments
(1,744 bp). At this peak resolution the two genes' instrument neighbourhoods
cannot be separated, so accessibility is anchored on TSS windows instead, where
the two genes are ~20 kb apart and separable. Windows are fixed at +/-2 kb
(primary) and +/-10 kb (sensitivity) and are not adjusted afterwards.

Peak sets differ between samples, so windows are defined in genomic coordinates
and summed over whatever peaks a given sample calls inside them; the number of
peaks per window per sample is reported so that "no peak called" is never
silently read as "no signal".

Usage:  python step35_atac_activation.py [--test]   (--test = one sample only)
"""
import gzip
import os
import sys
import time
import numpy as np
import pandas as pd

D = r"D:/Downloads/GSE282266"
MR = r"D:/R_ex/MR"
TPS = ["rest", "act_2.5", "act_5", "act_15"]
SETS = [1, 2, 3, 4]
WIN_PRIMARY, WIN_SENS = 2000, 10000

GLYCO = ["SLC2A1", "SLC2A3", "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP",
         "ALDOA", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2",
         "PKM", "LDHA", "PGM1", "PFKFB3", "PFKFB4"]
FOCUS = ["TPI1", "SPSB2"]
ACT_MARKERS = ["IL2RA", "CD69", "MYC", "IFNG", "IL2"]      # activation controls
REST_MARKERS = ["SELL", "CCR7", "TCF7", "KLF2"]            # should fall


def read_features(path):
    """CellRanger-ARC features.tsv.gz -> (genes DataFrame, peaks DataFrame)."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            p = line.rstrip("\n").split("\t")
            rows.append((i, p[0], p[2] if len(p) > 2 else "?",
                         p[3] if len(p) > 3 else "", p[4] if len(p) > 4 else "",
                         p[5] if len(p) > 5 else ""))
    df = pd.DataFrame(rows, columns=["row", "name", "kind", "chrom", "start", "end"])
    genes = df[df.kind == "Gene Expression"].copy()
    peaks = df[df.kind == "Peaks"].copy()
    for t in (genes, peaks):
        t["start"] = pd.to_numeric(t.start, errors="coerce")
        t["end"] = pd.to_numeric(t.end, errors="coerce")
    return genes, peaks


def load_sample(tp, s, want_genes, verbose=False):
    base = os.path.join(D, f"GSE282266_{tp}_set_{s}")
    genes, peaks = read_features(base + "_features.tsv.gz")
    t0 = time.time()
    m = pd.read_csv(base + "_matrix.mtx.gz", sep=r"\s+", comment="%",
                    skiprows=1, header=None, names=["r", "c", "v"],
                    dtype={"r": np.int32, "c": np.int32, "v": np.int32},
                    engine="c")
    if verbose:
        print(f"    matrix {len(m):,} nonzero entries in {time.time()-t0:.0f}s",
              flush=True)

    gene_rows = set(genes.row)
    peak_rows = set(peaks.row)
    is_gene = m.r.isin(gene_rows)
    rna_tot = m[is_gene].groupby("c").v.sum()
    atac_tot = m[~is_gene].groupby("c").v.sum()
    ncell = max(int(m.c.max()), 1)

    # RNA for requested genes
    gmap = genes.set_index("name").row
    out_rna = {}
    for g in want_genes:
        if g not in gmap.index:
            continue
        rr = gmap.loc[g]
        rr = int(rr) if np.isscalar(rr) or isinstance(rr, (int, np.integer)) else int(rr.iloc[0])
        sub = m[m.r == rr].set_index("c").v
        out_rna[g] = sub.reindex(range(1, ncell + 1), fill_value=0).values

    # ATAC in TSS windows
    out_atac, nwin = {}, {}
    pk = peaks.dropna(subset=["start", "end"])
    for g in want_genes:
        gi = genes[genes.name == g]
        if not len(gi):
            continue
        gi = gi.iloc[0]
        if pd.isna(gi.start):
            continue
        tss, ch = int(gi.start), gi.chrom
        for tag, w in (("w2k", WIN_PRIMARY), ("w10k", WIN_SENS)):
            sel = pk[(pk.chrom == ch) & (pk.end >= tss - w) & (pk.start <= tss + w)]
            nwin[(g, tag)] = len(sel)
            if not len(sel):
                out_atac[(g, tag)] = np.zeros(ncell)
                continue
            sub = m[m.r.isin(set(sel.row))].groupby("c").v.sum()
            out_atac[(g, tag)] = sub.reindex(range(1, ncell + 1), fill_value=0).values
    return dict(ncell=ncell,
                rna_tot=rna_tot.reindex(range(1, ncell + 1), fill_value=0).values,
                atac_tot=atac_tot.reindex(range(1, ncell + 1), fill_value=0).values,
                rna=out_rna, atac=out_atac, npeaks=nwin,
                n_genes=len(genes), n_peaks=len(peaks))


WANT = sorted(set(GLYCO + FOCUS + ACT_MARKERS + REST_MARKERS))

if "--test" in sys.argv:
    print("=== single-sample feasibility test: rest_set_1 ===", flush=True)
    t0 = time.time()
    r = load_sample("rest", 1, WANT, verbose=True)
    print(f"  cells={r['ncell']:,}  genes={r['n_genes']:,}  peaks={r['n_peaks']:,}")
    print(f"  median RNA/cell={np.median(r['rna_tot']):.0f}  "
          f"median ATAC/cell={np.median(r['atac_tot']):.0f}")
    print("\n  peaks found in each TSS window:")
    for g in FOCUS + ["GAPDH", "LDHA"]:
        print(f"    {g:<8} +/-2kb: {r['npeaks'].get((g,'w2k'),'-')}   "
              f"+/-10kb: {r['npeaks'].get((g,'w10k'),'-')}")
    print("\n  mean CP10K in this sample:")
    for g in FOCUS + ["GAPDH"]:
        if g in r["rna"]:
            cpm = 1e4 * r["rna"][g] / np.maximum(r["rna_tot"], 1)
            print(f"    {g:<8} RNA {cpm.mean():8.2f}   detected in "
                  f"{100*np.mean(r['rna'][g]>0):5.1f}% of cells")
    print(f"\n  elapsed {time.time()-t0:.0f}s for one sample "
          f"-> ~{16*(time.time()-t0)/60:.0f} min for all 16")
    sys.exit(0)

print("run with --test first to check feasibility")
