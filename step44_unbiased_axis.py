#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 44  Unbiased characterisation of the residual glycolysis axis.

Step 43 tested pre-defined immune modules and found nothing. Step 43b showed why
that was uninformative: the modules were built from genes this assay cannot
measure (IL5 0.6%, IL17A 0.6%, RORC 1.5%, FOXP3 4.1%) describing cell states a
15 h anti-CD3/CD28 culture of purified CD4 cells cannot produce. Testing for
identities that cannot exist returns a null for reasons unrelated to the axis.

Here nothing is pre-defined. After matching cells on activation score and depth,
cells are split on residual glycolysis and EVERY detectable gene is ranked. What
the axis is, is then read off the ranking rather than assumed.

PRE-SPECIFIED
  - detectability filter fixed in advance: a gene must be detected in >=5% of
    nuclei in every set. This is the filter whose absence invalidated Step 43.
  - the unit is the sample: 4 sets give 4 paired hi/lo pseudobulk observations.
    A gene is called only if all 4 sets agree in direction.
  - no gene set is consulted until the ranking exists.

Output: 44a (ranked genes), 44b (family composition), 44c (top lists)
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

GLYCO = ["SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA",
         "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3"]
ACTIV = ["IL2RA", "CD69", "TNFRSF9", "TNFRSF4", "ICOS", "BATF", "NFKB1",
         "REL", "IRF4", "MYC", "IL2"]


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
        v = sub[sub.r == r].set_index("c").v.reindex(range(1, ncell + 1), fill_value=0).values
        cp = 1e4 * v / np.maximum(tot, 1)
        s = cp.std()
        z.append((cp - cp.mean()) / s if s > 0 else np.zeros(ncell))
    return np.mean(z, axis=0) if z else np.zeros(ncell)


per_sample, detect_frames = [], []
for s in SETS:
    t0 = time.time()
    genes, g, ncell = load(f"{TP}_set_{s}")
    tot = g.groupby("c").v.sum().reindex(range(1, ncell + 1), fill_value=0).values
    nfeat = g.groupby("c").size().reindex(range(1, ncell + 1), fill_value=0).values
    det = g.groupby("r").size() / ncell
    r2n = dict(zip(genes.row, genes.name))
    detect_frames.append(pd.Series({r2n[r]: v for r, v in det.items()}, name=f"set{s}"))

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
        hi.append(h.nsmallest(n, "resid")); lo.append(l.nlargest(n, "resid"))
    hi = pd.concat(hi).cell.values
    lo = pd.concat(lo).cell.values

    hset, lset = set(hi), set(lo)
    gh = g[g.c.isin(hset)].groupby("r").v.sum()
    gl = g[g.c.isin(lset)].groupby("r").v.sum()
    th, tl = gh.sum(), gl.sum()
    df = pd.DataFrame({"hi": gh, "lo": gl}).fillna(0)
    df["gene"] = [r2n[i] for i in df.index]
    df[f"lfc"] = np.log2((1e6 * df.hi / th + 1) / (1e6 * df.lo / tl + 1))
    per_sample.append(df.set_index("gene")["lfc"].rename(f"set{s}"))
    print(f"  {TP}_set_{s}: {ncell:,} nuclei, matched {len(hi):,}/{len(lo):,}, "
          f"{time.time()-t0:.0f}s", flush=True)
    del g

lfc = pd.concat(per_sample, axis=1)
det = pd.concat(detect_frames, axis=1).fillna(0)
ok = det.min(axis=1) >= MIN_DETECT
lfc = lfc[lfc.index.isin(det[ok].index)].dropna()
print(f"\ngenes passing the >={MIN_DETECT:.0%} detection filter in all sets: "
      f"{len(lfc):,}")

lfc["mean_lfc"] = lfc[[f"set{s}" for s in SETS]].mean(axis=1)
lfc["n_pos"] = (lfc[[f"set{s}" for s in SETS]] > 0).sum(axis=1)
lfc["consistent"] = lfc.n_pos.isin([0, 4])
lfc["detect"] = det.loc[lfc.index].min(axis=1)
lfc = lfc.sort_values("mean_lfc", ascending=False)
lfc.to_csv(f"{MR}/44a_residual_axis_ranked_genes.tsv", sep="\t")

cons = lfc[lfc.consistent]
print(f"consistent in direction across all 4 sets: {len(cons):,}")
print("\n" + "=" * 74)
print("TOP 40 UP with residual glycolysis (4/4 consistent)")
print("=" * 74)
print(", ".join(cons.head(40).index))
print("\n" + "=" * 74)
print("TOP 40 DOWN (4/4 consistent)")
print("=" * 74)
print(", ".join(cons.tail(40).index[::-1]))

# ---- family composition, by name pattern (no external database needed) ----
FAM = {
    "ribosomal_protein": r"^(RPL|RPS)\d",
    "mito_encoded": r"^MT-",
    "OXPHOS_nuclear": r"^(NDUF|SDH[ABCD]|UQCR|COX\d|ATP5)",
    "glycolysis": r"^(HK[123]|GPI|PFK|ALDO|TPI1|GAPDH|PGK1|PGAM1|ENO[123]|PKM|LDHA|SLC2A)",
    "histone": r"^(H1-|H2A|H2B|H3-|H4-|HIST)",
    "HLA": r"^HLA-",
    "interferon": r"^(IFI|ISG|OAS|MX[12]|IFIT)",
    "cell_cycle": r"^(MKI67|TOP2A|CCN[ABE]|CDK[124]|CDC20|UBE2C|RRM2|TYMS|PCNA)",
    "heat_shock": r"^(HSP|DNAJ|HSPA)",
    "lncRNA_AC/AL": r"^(AC\d|AL\d|LINC)",
}
n_top = 300
fam_rows = []
for tag, direction, sub in (("up", 1, cons.head(n_top)), ("down", -1, cons.tail(n_top))):
    for fam, pat in FAM.items():
        k = sub.index.str.match(pat).sum()
        bg = lfc.index.str.match(pat).sum()
        fam_rows.append(dict(direction=tag, family=fam, n_in_top=int(k),
                             pct_in_top=100 * k / len(sub),
                             pct_background=100 * bg / len(lfc)))
fam = pd.DataFrame(fam_rows)
fam["enrichment"] = fam.pct_in_top / fam.pct_background.replace(0, np.nan)
fam.to_csv(f"{MR}/44b_family_composition.tsv", sep="\t", index=False)
print("\n" + "=" * 74)
print(f"family composition of the top {n_top} genes each way")
print("=" * 74)
print(fam.pivot(index="family", columns="direction",
                values=["pct_in_top", "enrichment"]).round(2).to_string())

pd.DataFrame({"up": pd.Series(cons.head(200).index),
              "down": pd.Series(cons.tail(200).index[::-1])}
             ).to_csv(f"{MR}/44c_top_lists.tsv", sep="\t", index=False)
print("\nwritten: 44a-44c")
