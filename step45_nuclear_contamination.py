#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 45  Is the residual glycolysis axis a nuclear/cytoplasmic composition
         artefact rather than a cell state?

Step 44 found the axis runs: ribosomal proteins (12% of the top 300, 10x
enriched) and glycolytic enzymes UP, nuclear-retained lncRNAs (7x enriched) and
T-cell identity genes DOWN. In single-NUCLEUS data that is precisely the
signature of variable cytoplasmic carry-over: abundant cytoplasmic mRNAs rise
together while nuclear-retained transcripts fall.

TEST
  Nuclear-retention index = fraction of a nucleus's counts in MALAT1 and NEAT1,
  the canonical nuclear-retained lncRNAs. If the residual glycolysis axis is a
  contamination axis, it will be strongly negatively correlated with this index,
  and adjusting for it should collapse the axis.

DECISION RULE, fixed before running
  |rho| > 0.4 between residual glycolysis and the nuclear index, or a fall of
  more than half in the ribosomal-protein enrichment after adjustment, means the
  axis is substantially technical and cannot be interpreted as a cell state.
"""
import sys
import gzip
import os
import time
import numpy as np
import pandas as pd
from scipy import stats

D = r"D:/Downloads/GSE282266"
MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
SETS = [1, 2, 3, 4]
TP = "act_15"
NUCLEAR = ["MALAT1", "NEAT1"]
GLYCO = ["SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA",
         "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3"]
ACTIV = ["IL2RA", "CD69", "TNFRSF9", "TNFRSF4", "ICOS", "BATF", "NFKB1",
         "REL", "IRF4", "MYC", "IL2"]

rows_out = []
for s in SETS:
    t0 = time.time()
    base = os.path.join(D, f"GSE282266_{TP}_set_{s}")
    fr = []
    with gzip.open(base + "_features.tsv.gz", "rt") as fh:
        for i, line in enumerate(fh, start=1):
            p = line.rstrip("\n").split("\t")
            fr.append((i, p[0], p[2] if len(p) > 2 else "?"))
    feat = pd.DataFrame(fr, columns=["row", "name", "kind"])
    genes = feat[feat.kind == "Gene Expression"]
    n2r = dict(zip(genes.name, genes.row))

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
    ncell = int(g.c.max())
    tot = g.groupby("c").v.sum().reindex(range(1, ncell + 1), fill_value=0).values
    nfeat = g.groupby("c").size().reindex(range(1, ncell + 1), fill_value=0).values

    def vec(name):
        r = n2r.get(name)
        if r is None:
            return np.zeros(ncell)
        return g[g.r == r].set_index("c").v.reindex(
            range(1, ncell + 1), fill_value=0).values

    def mod(lst):
        z = []
        for x in lst:
            v = 1e4 * vec(x) / np.maximum(tot, 1)
            sd = v.std()
            z.append((v - v.mean()) / sd if sd > 0 else np.zeros(ncell))
        return np.mean(z, axis=0)

    nuc_counts = np.sum([vec(x) for x in NUCLEAR], axis=0)
    nuc_frac = nuc_counts / np.maximum(tot, 1)
    rp_rows = [n2r[n] for n in genes.name if n.startswith(("RPL", "RPS"))
               and n in n2r]
    rp = g[g.r.isin(set(rp_rows))].groupby("c").v.sum().reindex(
        range(1, ncell + 1), fill_value=0).values
    rp_frac = rp / np.maximum(tot, 1)

    gly, act = mod(GLYCO), mod(ACTIV)
    keep = nfeat >= 200
    X = np.column_stack([np.ones(keep.sum()), act[keep], np.log1p(nfeat[keep])])
    beta, *_ = np.linalg.lstsq(X, gly[keep], rcond=None)
    resid = gly[keep] - X @ beta

    nf, rf = nuc_frac[keep], rp_frac[keep]
    r_nuc = stats.spearmanr(resid, nf).statistic
    r_rp = stats.spearmanr(resid, rf).statistic
    # after additionally adjusting for the nuclear index
    X2 = np.column_stack([X, nf])
    b2, *_ = np.linalg.lstsq(X2, gly[keep], rcond=None)
    resid2 = gly[keep] - X2 @ b2
    var_drop = 1 - resid2.var() / resid.var()

    rows_out.append(dict(sample=f"set{s}", n=int(keep.sum()),
                         median_nuclear_frac=float(np.median(nf)),
                         median_RP_frac=float(np.median(rf)),
                         rho_resid_nuclear=r_nuc, rho_resid_RP=r_rp,
                         var_removed_by_nuclear_index=var_drop))
    print(f"  set{s}: n={keep.sum():,}  MALAT1+NEAT1 frac median={np.median(nf):.3f}  "
          f"RP frac median={np.median(rf):.3f}  "
          f"rho(resid,nuclear)={r_nuc:+.3f}  rho(resid,RP)={r_rp:+.3f}  "
          f"var removed={var_drop:.3f}", flush=True)
    del g

res = pd.DataFrame(rows_out)
res.to_csv(f"{MR}/45a_nuclear_contamination_check.tsv", sep="\t", index=False)
print("\n" + "=" * 76)
print("verdict")
print("=" * 76)
mn = res.rho_resid_nuclear.mean()
mr = res.rho_resid_RP.mean()
mv = res.var_removed_by_nuclear_index.mean()
print(f"  mean rho(residual glycolysis, nuclear-retention index) = {mn:+.3f}")
print(f"  mean rho(residual glycolysis, ribosomal-protein fraction) = {mr:+.3f}")
print(f"  mean variance removed by adding the nuclear index = {mv:.3f}")
if abs(mn) > 0.4 or mv > 0.5:
    print("\n  -> the axis is substantially a nuclear/cytoplasmic composition"
          " artefact;\n     it cannot be interpreted as a cell state")
else:
    print("\n  -> not explained by nuclear/cytoplasmic composition")
print("\nwritten: 45a")
