#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 106 -- define the metabolic axis without TPI1, then ask where TPI1 sits.

Raised in review: the axis is built from a glycolysis score that contains TPI1,
and its glycolysis enrichment counts TPI1 among the enriched genes; the paper
then says TPI1 marks that axis. Part of that is self-fulfilling.

There are two circular couplings in step44 and this removes both:
  1. `gly = module(..., GLYCO, ...)` -- GLYCO contains TPI1, so the residual
     glycolysis score that DEFINES the axis is partly TPI1's own expression.
  2. the glycolysis family regex used for the enrichment contains TPI1, so TPI1
     is one of the 15 genes generating the 21.3-fold figure.

Everything else is held identical to step44 -- same detection filter, same
activation and depth matching, same 4/4 direction rule, same families, same
random-free procedure -- so the only difference is TPI1's participation.

TPI1 then becomes a held-out test rather than part of the construction: its rank
in the axis derived without it is the quantity of interest.

Outputs: 106a_axis_noTPI1_ranked.tsv, 106b_family_noTPI1.tsv, 106c_comparison.tsv
"""
import gzip
import os
import re
import time

import numpy as np
import pandas as pd

D = r"D:/Downloads/GSE282266"
MR = r"D:/R_ex/MR"
SETS = [1, 2, 3, 4]
TP = "act_15"
MIN_DETECT = 0.05
TOP_N = 300

GLYCO_FULL = ["SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA",
              "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3"]
GLYCO = [g for g in GLYCO_FULL if g != "TPI1"]          # <- coupling 1 removed
ACTIV = ["IL2RA", "CD69", "TNFRSF9", "TNFRSF4", "ICOS", "BATF", "NFKB1",
         "REL", "IRF4", "MYC", "IL2"]

# coupling 2 removed: TPI1 dropped from the glycolysis family pattern
FAM = {
    "ribosomal_protein": r"^(RPL|RPS)\d",
    "mito_encoded": r"^MT-",
    "OXPHOS_nuclear": r"^(NDUF|SDH[ABCD]|UQCR|COX\d|ATP5)",
    "glycolysis_noTPI1": r"^(HK[123]|GPI|PFK|ALDO|GAPDH|PGK1|PGAM1|ENO[123]|PKM|LDHA|SLC2A)",
    "histone": r"^(H1-|H2A|H2B|H3-|H4-|HIST)",
    "HLA": r"^HLA-",
}


def load(sample):
    base = os.path.join(D, f"GSE282266_{sample}")
    rows = []
    with gzip.open(base + "_features.tsv.gz", "rt") as fh:
        for i, line in enumerate(fh, start=1):
            p = line.rstrip("\n").split("\t")
            rows.append((i, p[0], p[2] if len(p) > 2 else "?"))
    feat = pd.DataFrame(rows, columns=["row", "name", "kind"])
    genes = feat[feat.kind == "Gene Expression"]
    with gzip.open(base + "_matrix.mtx.gz", "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                nrow, ncol, nnz = (int(x) for x in line.split())
                break
    m = pd.read_csv(base + "_matrix.mtx.gz", sep=r"\s+", comment="%", skiprows=1,
                    header=None, names=["r", "c", "v"],
                    dtype={"r": np.int32, "c": np.int32, "v": np.int32}, engine="c")
    if len(m) == nnz + 1:
        m = m.iloc[1:]
    g = m[m.r.isin(set(genes.row))]
    del m
    return genes, g, int(g.c.max())


def module(g, genes, gene_list, ncell, tot):
    n2r = dict(zip(genes.name, genes.row))
    rows = [n2r[x] for x in gene_list if x in n2r]
    sub = g[g.r.isin(set(rows))]
    z = []
    for r in rows:
        v = sub[sub.r == r].set_index("c").v.reindex(range(1, ncell + 1),
                                                     fill_value=0).values
        cp = 1e4 * v / np.maximum(tot, 1)
        s = cp.std()
        z.append((cp - cp.mean()) / s if s > 0 else np.zeros(ncell))
    return np.mean(z, axis=0) if z else np.zeros(ncell)


def main():
    print(f"glycolysis score genes: {len(GLYCO)} (TPI1 removed from "
          f"{len(GLYCO_FULL)})")
    per_sample, detect_frames = [], []
    for s in SETS:
        t0 = time.time()
        genes, g, ncell = load(f"{TP}_set_{s}")
        tot = g.groupby("c").v.sum().reindex(range(1, ncell + 1), fill_value=0).values
        nfeat = g.groupby("c").size().reindex(range(1, ncell + 1), fill_value=0).values
        det = g.groupby("r").size() / ncell
        r2n = dict(zip(genes.row, genes.name))
        detect_frames.append(pd.Series({r2n[r]: v for r, v in det.items()},
                                       name=f"set{s}"))

        gly = module(g, genes, GLYCO, ncell, tot)
        act = module(g, genes, ACTIV, ncell, tot)
        keep = nfeat >= 200
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
            hi.append(h.nsmallest(n, "resid"))
            lo.append(l.nlargest(n, "resid"))
        hi = pd.concat(hi).cell.values
        lo = pd.concat(lo).cell.values

        gh = g[g.c.isin(set(hi))].groupby("r").v.sum()
        gl = g[g.c.isin(set(lo))].groupby("r").v.sum()
        th, tl = gh.sum(), gl.sum()
        df = pd.DataFrame({"hi": gh, "lo": gl}).fillna(0)
        df["gene"] = [r2n[i] for i in df.index]
        df["lfc"] = np.log2((1e6 * df.hi / th + 1) / (1e6 * df.lo / tl + 1))
        per_sample.append(df.set_index("gene")["lfc"].rename(f"set{s}"))
        print(f"  set_{s}: {ncell:,} nuclei, matched {len(hi):,}/{len(lo):,}, "
              f"{time.time()-t0:.0f}s", flush=True)
        del g

    lfc = pd.concat(per_sample, axis=1)
    det = pd.concat(detect_frames, axis=1).fillna(0)
    ok = det.min(axis=1) >= MIN_DETECT
    lfc = lfc[lfc.index.isin(det[ok].index)].dropna()
    lfc["mean_lfc"] = lfc[[f"set{s}" for s in SETS]].mean(axis=1)
    lfc["n_pos"] = (lfc[[f"set{s}" for s in SETS]] > 0).sum(axis=1)
    lfc["consistent"] = lfc.n_pos.isin([0, 4])
    lfc = lfc.sort_values("mean_lfc", ascending=False)
    lfc.to_csv(f"{MR}/106a_axis_noTPI1_ranked.tsv", sep="\t")
    print(f"\ngenes passing detection filter in all sets: {len(lfc):,}")

    # ------------------------------------------------ family enrichment
    cons = lfc[lfc.consistent]
    up = cons.head(TOP_N).index
    rows = []
    for fam, pat in FAM.items():
        rx = re.compile(pat)
        n_top = sum(bool(rx.match(x)) for x in up)
        n_bg = sum(bool(rx.match(x)) for x in lfc.index)
        pct_top = 100 * n_top / len(up)
        pct_bg = 100 * n_bg / len(lfc)
        rows.append(dict(family=fam, n_in_top=n_top, pct_in_top=round(pct_top, 2),
                         pct_background=round(pct_bg, 3),
                         enrichment=round(pct_top / pct_bg, 2) if pct_bg else np.nan))
    fam_df = pd.DataFrame(rows)
    fam_df.to_csv(f"{MR}/106b_family_noTPI1.tsv", sep="\t", index=False)
    print("\nfamily composition of the top", TOP_N, "(axis defined without TPI1)")
    print(fam_df.to_string(index=False))

    # ------------------------------------------------ TPI1 as a held-out test
    n = len(lfc)
    out = []
    if "TPI1" in lfc.index:
        r = lfc.index.get_loc("TPI1") + 1
        row = lfc.loc["TPI1"]
        out.append(dict(gene="TPI1", rank=r, of=n, percentile=round(100 * r / n, 2),
                        mean_lfc=round(float(row.mean_lfc), 4),
                        n_sets_positive=int(row.n_pos),
                        consistent=bool(row.consistent)))
        print(f"\nHELD-OUT TEST: TPI1 ranks {r} of {n:,} "
              f"(top {100*r/n:.2f}%) on an axis built without it; "
              f"mean lfc {row.mean_lfc:+.4f}, positive in {int(row.n_pos)}/4 sets")
    else:
        print("\nTPI1 did not pass the detection filter in all sets")

    # the enzymes MR cannot see, for comparison
    for gname in ("PGAM1", "GAPDH", "ENO1", "PKM", "LDHA"):
        if gname in lfc.index:
            r = lfc.index.get_loc(gname) + 1
            row = lfc.loc[gname]
            out.append(dict(gene=gname, rank=r, of=n,
                            percentile=round(100 * r / n, 2),
                            mean_lfc=round(float(row.mean_lfc), 4),
                            n_sets_positive=int(row.n_pos),
                            consistent=bool(row.consistent)))
    cmp_df = pd.DataFrame(out)
    cmp_df.to_csv(f"{MR}/106c_comparison.tsv", sep="\t", index=False)
    print("\nranks of the glycolytic enzymes on the TPI1-free axis")
    print(cmp_df.to_string(index=False))

    pub = pd.read_csv(f"{MR}/44b_family_composition.tsv", sep="\t")
    g_pub = pub[(pub.direction == "up") & (pub.family == "glycolysis")].iloc[0]
    g_new = fam_df[fam_df.family == "glycolysis_noTPI1"].iloc[0]
    print(f"\nglycolysis enrichment: published {g_pub.enrichment:.1f}-fold "
          f"({int(g_pub.n_in_top)} genes, TPI1 included in both score and set) -> "
          f"{g_new.enrichment:.1f}-fold ({int(g_new.n_in_top)} genes, TPI1 in neither)")
    print("\nwrote 106a / 106b / 106c")


if __name__ == "__main__":
    main()
