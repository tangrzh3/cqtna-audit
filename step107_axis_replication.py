#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 107 -- does the metabolic axis, and TPI1's place on it, replicate?

Executes manuscript/PREREG_axis_replication.md (S32).

Step 106 rebuilt the axis without TPI1 and found TPI1 at rank 19 of 7,653 on it.
Those four sets are one study. This repeats the construction in GSE166188
DOGMA-seq, where CD4 is called by surface protein rather than by transcript, and
asks the same held-out question.

The axis is again defined by a glycolysis score that does not contain TPI1, so
TPI1's rank here is a prediction being tested, not a quantity being fitted.

Registered criterion: TPI1 in the top 5% with the same direction in both arms.

Output: 107a_replication.tsv
"""
import gzip
import os

import numpy as np
import pandas as pd

D = r"D:/Downloads/GSE166188"
MR = r"D:/R_ex/MR"
MIN_DETECT = 0.05
TOP_N = 300
CRIT_PCT = 5.0

SAMPLES = {
    "LLL_stim": dict(gex="GSM5065528_LLL_STIM_GExp_ATAC_filtered",
                     adt="GSM5065529_LLL_stim_ADT_allCounts"),
    "DIG_stim": dict(gex="GSM5065534_DIG_STIM_GExp_ATAC_filtered",
                     adt="GSM5065535_DIG_stim_ADT_allCounts"),
}
# TPI1 absent by construction (S32 section 3)
GLYCO = ["SLC2A1", "SLC2A3", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
         "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3"]
ACTIV = ["IL2RA", "CD69", "TNFRSF9", "TNFRSF4", "ICOS", "BATF", "NFKB1",
         "REL", "IRF4", "MYC", "IL2"]
CD4_TAGS = {"CD4", "CD4-1", "CD4-2"}
CD8_TAGS = {"CD8", "CD8a", "CD8A", "CD8-1", "CD8-2", "CD8b"}


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


def score(g, rows_for, names, ncell, tot):
    z = []
    for nm in names:
        r = rows_for.get(nm)
        if r is None:
            continue
        v = g[g.r == r].set_index("c").v.reindex(range(1, ncell + 1),
                                                 fill_value=0).values
        cp = 1e4 * v / np.maximum(tot, 1)
        s = cp.std()
        if s > 0:
            z.append((cp - cp.mean()) / s)
    return np.mean(z, axis=0) if z else np.zeros(ncell)


def run(name, cfg):
    print("=" * 78)
    print(name)
    # ---------------- protein: gate CD4 by antibody, not transcript
    prot = read_lines(os.path.join(D, cfg["adt"] + ".proteins.txt.gz"))
    abc = [norm_bc(b) for b in read_lines(os.path.join(D, cfg["adt"] + ".barcodes.txt.gz"))]
    am, anr, anc = read_mtx(os.path.join(D, cfg["adt"] + ".mtx.gz"))
    pidx = {p: i + 1 for i, p in enumerate(prot)}
    cd4p = [p for p in prot if p in CD4_TAGS]
    cd8p = [p for p in prot if p in CD8_TAGS]
    if not cd4p or not cd8p:
        print("  CD4/CD8 antibody missing -> interpretation D")
        return None
    # ⚠ The ADT matrix is transposed relative to the RNA one: rows are BARCODES
    # (711,378) and columns are PROTEINS (210). Indexing it the other way round
    # -- as the RNA convention would suggest, and as step42 does -- silently
    # returns one barcode's profile in place of one protein's, and gated 0 cells
    # on the first two runs. Asserted rather than assumed.
    assert anr == len(abc) and anc == len(prot), (
        f"unexpected ADT orientation: mtx {anr}x{anc}, "
        f"{len(abc)} barcodes, {len(prot)} proteins")
    tot_adt = am.groupby("r").v.sum().reindex(range(1, len(abc) + 1), fill_value=0).values

    def adt(tag):
        return am[am.c == pidx[tag]].set_index("r").v.reindex(
            range(1, len(abc) + 1), fill_value=0).values
    clr = lambda x: np.log1p(1e4 * x / np.maximum(tot_adt, 1))
    p4 = clr(np.sum([adt(t) for t in cd4p], axis=0))
    p8 = clr(np.sum([adt(t) for t in cd8p], axis=0))
    pdf = pd.DataFrame(dict(bc=abc, CD4=p4, CD8=p8))
    print(f"  ADT barcodes: {len(abc):,} (unfiltered droplets)")

    # ---------------- RNA
    feat = []
    with gzip.open(os.path.join(D, cfg["gex"] + "_features.tsv.gz"), "rt") as fh:
        for i, line in enumerate(fh, start=1):
            p = line.rstrip("\n").split("\t")
            feat.append((i, p[1] if len(p) > 1 else p[0], p[2] if len(p) > 2 else "?"))
    fdf = pd.DataFrame(feat, columns=["row", "name", "kind"])
    genes = fdf[fdf.kind.str.contains("Gene Expression", na=False)]
    if genes.empty:
        genes = fdf[~fdf.kind.str.contains("Peaks", na=False)]
    gbc = [norm_bc(b) for b in read_lines(os.path.join(D, cfg["gex"] + "_barcodes.tsv.gz"))]
    gm, _, _ = read_mtx(os.path.join(D, cfg["gex"] + "_matrix.mtx.gz"))
    g = gm[gm.r.isin(set(genes.row))]
    del gm
    ncell = len(gbc)
    # The ADT file lists every droplet while the RNA file is cell-filtered, so
    # the CD4/CD8 medians must be taken over REAL CELLS. Computing them over all
    # 711k droplets is dominated by empties and gates almost everything out --
    # 6 barcodes survived on the first run, which is what interpretation D caught.
    gset = pd.DataFrame(dict(bc=gbc, col=np.arange(1, ncell + 1)))
    d = gset.merge(pdf, on="bc", how="inner")
    print(f"  RNA barcodes {ncell:,}; with ADT: {len(d):,}")
    if len(d) < 300:
        print("  too few barcodes shared between RNA and ADT -> interpretation D")
        return None
    gate = (d.CD4 > d.CD4.median()) & (d.CD8 <= d.CD8.median())
    cells = d.col[gate].values
    print(f"  protein-gated CD4 cells: {len(cells):,} of {len(d):,} "
          f"({100*len(cells)/len(d):.1f}%)")
    if len(cells) < 300:
        print("  too few CD4 cells -> interpretation D")
        return None
    g = g[g.c.isin(set(cells))]

    tot = g.groupby("c").v.sum().reindex(range(1, ncell + 1), fill_value=0).values
    nfeat = g.groupby("c").size().reindex(range(1, ncell + 1), fill_value=0).values
    rows_for = dict(zip(genes.name, genes.row))
    r2n = dict(zip(genes.row, genes.name))

    gly = score(g, rows_for, GLYCO, ncell, tot)
    act = score(g, rows_for, ACTIV, ncell, tot)
    mask = np.zeros(ncell, bool)
    mask[cells - 1] = True
    keep = mask & (nfeat >= 200)
    X = np.column_stack([np.ones(keep.sum()), act[keep], np.log1p(nfeat[keep])])
    beta, *_ = np.linalg.lstsq(X, gly[keep], rcond=None)
    resid = np.full(ncell, np.nan)
    resid[keep] = gly[keep] - X @ beta

    cd = pd.DataFrame(dict(cell=np.arange(1, ncell + 1), act=act, nfeat=nfeat,
                           resid=resid)).dropna()
    cd["abin"] = pd.qcut(cd.act, 10, labels=False, duplicates="drop")
    cd["dbin"] = pd.qcut(np.log1p(cd.nfeat), 5, labels=False, duplicates="drop")
    hi, lo = [], []
    for _, z in cd.groupby(["abin", "dbin"]):
        if len(z) < 20:
            continue
        q40, q60 = z.resid.quantile([.4, .6])
        h, l = z[z.resid > q60], z[z.resid <= q40]
        n = min(len(h), len(l))
        if n:
            hi.append(h.nsmallest(n, "resid"))
            lo.append(l.nlargest(n, "resid"))
    if not hi:
        print("  no matched hi/lo strata -> interpretation D")
        return None
    hi = pd.concat(hi).cell.values
    lo = pd.concat(lo).cell.values
    print(f"  matched {len(hi):,} hi / {len(lo):,} lo")

    det = g.groupby("r").size() / len(cells)
    gh = g[g.c.isin(set(hi))].groupby("r").v.sum()
    gl = g[g.c.isin(set(lo))].groupby("r").v.sum()
    df = pd.DataFrame({"hi": gh, "lo": gl}).fillna(0)
    df["gene"] = [r2n[i] for i in df.index]
    df["detect"] = [det.get(i, 0) for i in df.index]
    df = df[df.detect >= MIN_DETECT]
    th, tl = gh.sum(), gl.sum()
    df["lfc"] = np.log2((1e6 * df.hi / th + 1) / (1e6 * df.lo / tl + 1))
    df = df.groupby("gene", as_index=False).lfc.mean().sort_values(
        "lfc", ascending=False).reset_index(drop=True)
    n = len(df)
    print(f"  genes passing detection filter: {n:,}")

    res = dict(arm=name, n_cd4=len(cells), n_hi=len(hi), n_genes=n)
    for gname in ("TPI1", "PGAM1", "GAPDH", "ENO1", "PKM", "LDHA"):
        hit = df.index[df.gene == gname]
        if len(hit):
            r = int(hit[0]) + 1
            res[f"{gname}_rank"] = r
            res[f"{gname}_pct"] = round(100 * r / n, 3)
            res[f"{gname}_lfc"] = round(float(df.lfc.iloc[hit[0]]), 4)
        else:
            res[f"{gname}_rank"] = np.nan
            res[f"{gname}_pct"] = np.nan
            res[f"{gname}_lfc"] = np.nan
    top = set(df.gene.head(TOP_N))
    n_top = sum(x in top for x in GLYCO)
    n_bg = sum(x in set(df.gene) for x in GLYCO)
    res["glyco_in_top"] = n_top
    res["glyco_detected"] = n_bg
    res["glyco_enrichment"] = round((n_top / TOP_N) / (n_bg / n), 2) if n_bg else np.nan
    print(f"  TPI1 rank {res['TPI1_rank']} of {n:,} (top {res['TPI1_pct']}%), "
          f"lfc {res['TPI1_lfc']:+}")
    print(f"  glycolysis (TPI1-free) in top {TOP_N}: {n_top}/{n_bg} -> "
          f"{res['glyco_enrichment']}x")
    del g
    return res


def main():
    rows = [r for r in (run(k, v) for k, v in SAMPLES.items()) if r]
    if not rows:
        print("\nno arm produced a result -> interpretation D")
        return
    out = pd.DataFrame(rows)
    out.to_csv(f"{MR}/107a_replication.tsv", sep="\t", index=False)

    print("\n" + "=" * 78)
    print("VERDICT (S32 section 5)")
    print("=" * 78)
    ok = [(r["arm"], r["TPI1_pct"], r["TPI1_lfc"]) for r in rows]
    met = [a for a, p, l in ok if not np.isnan(p) and p <= CRIT_PCT and l > 0]
    for a, p, l in ok:
        print(f"  {a}: TPI1 top {p}%  lfc {l:+}  "
              f"{'meets' if a in met else 'does not meet'} the 5% criterion")
    if len(met) == len(rows):
        v = "A: replicates in both arms"
    elif met:
        v = "B: replicates in one arm only"
    else:
        v = "C: does not replicate -- state membership downgrades to discovery-only"
    print(f"  -> {v}")
    print("\n  ranks of the enzymes MR cannot instrument (discovery: all above TPI1)")
    for r in rows:
        s = "  ".join(f"{g} {r[f'{g}_rank']}" for g in
                      ("PGAM1", "GAPDH", "ENO1", "PKM", "LDHA")
                      if not (isinstance(r[f"{g}_rank"], float)
                              and np.isnan(r[f"{g}_rank"])))
        print(f"    {r['arm']}: TPI1 {r['TPI1_rank']}   {s}")
    print("\n  S32 section 2: LLL and DIG are lysis buffers, i.e. technical "
          "replicates of the\n  same material. This replicates across study, "
          "platform and lineage calling,\n  NOT across donors.")
    print("\nwrote 107a")


if __name__ == "__main__":
    main()
