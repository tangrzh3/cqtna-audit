#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 43  THRESHOLD TEST: in activated CD4 T cells, is the glycolytic programme
         separable from activation state at all?

WHY THIS COMES FIRST
Step 20's phenotype collapsed because splitting cells by a score also split them
by cell-type purity. GSE282266 removes that confound (FACS-purified CD4, 10x,
~160k nuclei), but a new one takes its place: in an activation time course,
"glycolysis-high" cells may simply be "more activated" cells. If so, any
phenotype we describe is activation, not glycolysis, and the claim
"CD4 glycolysis marks ICB non-response" reduces to "CD4 activation marks it".
That is the same trap in different clothing, so it is tested before anything
else and before any motif or regulon work is attempted.

PRE-SPECIFIED DECISION RULE
  - variance of the glycolysis score explained by activation score + depth:
      R2 > 0.70  -> largely inseparable; the line stops here
  - after matching cells on activation score and depth, split on residual
    glycolysis: if no biological module exceeds the 95th percentile of 50
    expression-matched RANDOM modules, there is no independent axis
  - POSITIVE CONTROL for the test itself: proliferation is a genuinely distinct
    subset, so it must separate. If nothing separates, including proliferation,
    the test is underpowered and its null is uninformative -- exactly the
    standard applied in Steps 19 and 24.

Data: GSE282266 act_15 (15 h, 4 sets). 10x Multiome RNA is snRNA-seq, so
per-cell module scores are sparse; that is the main power concern and is why
the positive control matters.

Output: 43a-43c
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
rng = np.random.default_rng(1)

MODULES = {
    "Glycolysis": ["SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKL", "PFKP",
                   "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM",
                   "LDHA", "PFKFB3"],
    "Activation": ["IL2RA", "CD69", "TNFRSF9", "TNFRSF4", "ICOS", "BATF",
                   "NFKB1", "REL", "IRF4", "MYC", "IL2"],
    "Proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA", "TYMS",
                      "STMN1", "TUBA1B", "RRM2", "UBE2C"],
    "Exhaustion": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"],
    "Treg": ["FOXP3", "IKZF2", "IL2RA", "CTLA4", "TNFRSF18"],
    "Th1": ["TBX21", "IFNG", "CXCR3", "IL12RB2"],
    "Th2": ["GATA3", "IL4", "IL5", "IL13", "CCR4"],
    "Th17": ["RORC", "IL17A", "IL23R", "CCR6"],
    "Tfh": ["BCL6", "CXCR5", "PDCD1", "IL21"],
    "Cytotoxic": ["GZMB", "GZMK", "PRF1", "NKG7", "GNLY"],
    "OXPHOS": ["NDUFA4", "NDUFB2", "SDHB", "UQCRB", "COX5A", "COX7C", "ATP5F1A"],
    "IFN_response": ["ISG15", "IFI6", "MX1", "OAS1", "IFIT3", "STAT1", "IRF7"],
}
N_RANDOM = 50


def load_features(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            p = line.rstrip("\n").split("\t")
            rows.append((i, p[0], p[2] if len(p) > 2 else "?"))
    df = pd.DataFrame(rows, columns=["row", "name", "kind"])
    return df[df.kind == "Gene Expression"].copy()


def cell_scores(sample):
    base = os.path.join(D, f"GSE282266_{sample}")
    genes = load_features(base + "_features.tsv.gz")
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
    gene_rows = set(genes.row)
    g = m[m.r.isin(gene_rows)]
    ncell = int(m.c.max())
    tot = g.groupby("c").v.sum().reindex(range(1, ncell + 1), fill_value=0).values
    nfeat = g.groupby("c").size().reindex(range(1, ncell + 1), fill_value=0).values
    del m

    name2row = dict(zip(genes.name, genes.row))
    # per-gene detection, for expression-matched random modules
    det = g.groupby("r").size()
    row2name = dict(zip(genes.row, genes.name))
    detect = pd.Series({row2name[r]: n / ncell for r, n in det.items()})

    wanted = sorted({x for v in MODULES.values() for x in v})
    pool = detect[(detect > 0.02) & (detect < 0.9)].index.tolist()
    rand_sets = {}
    for i in range(N_RANDOM):
        k = len(MODULES["Glycolysis"])
        rand_sets[f"rand{i}"] = list(rng.choice(pool, k, replace=False))
    wanted += sorted({x for v in rand_sets.values() for x in v})
    wanted = sorted(set(wanted) & set(name2row))

    rows_needed = {name2row[w] for w in wanted}
    sub = g[g.r.isin(rows_needed)]
    mat = {}
    for gene in wanted:
        r = name2row[gene]
        s = sub[sub.r == r].set_index("c").v
        mat[gene] = s.reindex(range(1, ncell + 1), fill_value=0).values
    del g, sub

    cp = pd.DataFrame(mat)
    cp = 1e4 * cp.div(np.maximum(tot, 1), axis=0)
    z = (cp - cp.mean()) / cp.std().replace(0, np.nan)
    z = z.fillna(0)

    out = pd.DataFrame(dict(sample=sample, cell=np.arange(1, ncell + 1),
                            total=tot, nFeature=nfeat))
    for nm, gs in {**MODULES, **rand_sets}.items():
        gg = [x for x in gs if x in z.columns]
        out[nm] = z[gg].mean(axis=1) if gg else np.nan
    return out, list(rand_sets)


frames, randnames = [], None
for s in SETS:
    t0 = time.time()
    f, rn = cell_scores(f"{TP}_set_{s}")
    randnames = rn
    frames.append(f)
    print(f"  {TP}_set_{s}: {len(f):,} nuclei, {time.time()-t0:.0f}s", flush=True)
d = pd.concat(frames, ignore_index=True)
d = d[(d.nFeature >= 200)]
print(f"\nnuclei after QC: {len(d):,}")
d.to_csv(f"{MR}/43a_GSE282266_cell_scores.tsv.gz", sep="\t", index=False,
         compression="gzip")

# ================================================================= 1. overlap
print("\n" + "=" * 76)
print("1  how much of glycolysis is activation?")
print("=" * 76)
for s, sub in d.groupby("sample"):
    r_raw = stats.spearmanr(sub.Glycolysis, sub.Activation).statistic
    X = np.column_stack([np.ones(len(sub)), sub.Activation, np.log1p(sub.nFeature)])
    y = sub.Glycolysis.values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    r2 = 1 - resid.var() / y.var()
    print(f"  {s:<16}rho(glyco,activation)={r_raw:+.3f}   "
          f"R2(glyco ~ activation + depth)={r2:.3f}")

X = np.column_stack([np.ones(len(d)), d.Activation, np.log1p(d.nFeature)])
beta, *_ = np.linalg.lstsq(X, d.Glycolysis.values, rcond=None)
d["Glyco_resid"] = d.Glycolysis.values - X @ beta
R2_all = 1 - d.Glyco_resid.var() / d.Glycolysis.var()
print(f"\n  pooled R2 = {R2_all:.3f}   "
      f"-> {'INSEPARABLE (>0.70), line stops' if R2_all > .70 else 'residual axis remains'}")

# ============================================== 2. residual split vs random null
print("\n" + "=" * 76)
print("2  split on residual glycolysis, matched on activation and depth")
print("=" * 76)
res = []
for s, sub in d.groupby("sample"):
    sub = sub.copy()
    sub["abin"] = pd.qcut(sub.Activation, 10, labels=False, duplicates="drop")
    sub["dbin"] = pd.qcut(np.log1p(sub.nFeature), 5, labels=False, duplicates="drop")
    hi, lo = [], []
    for (a, b), z in sub.groupby(["abin", "dbin"]):
        if len(z) < 20:
            continue
        q40, q60 = z.Glyco_resid.quantile([.4, .6])
        h, l = z[z.Glyco_resid > q60], z[z.Glyco_resid <= q40]
        n = min(len(h), len(l))
        hi.append(h.nsmallest(n, "Glyco_resid")); lo.append(l.nlargest(n, "Glyco_resid"))
    if not hi:
        continue
    hi, lo = pd.concat(hi), pd.concat(lo)
    for mod in list(MODULES) + randnames:
        if mod == "Glycolysis":
            continue
        a, b = hi[mod].values, lo[mod].values
        pooled = np.sqrt((a.var() + b.var()) / 2)
        smd = (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan
        res.append(dict(sample=s, module=mod, smd=smd,
                        kind="random" if mod.startswith("rand") else "biological",
                        n_hi=len(hi), n_lo=len(lo)))
res = pd.DataFrame(res)
res.to_csv(f"{MR}/43b_residual_split_modules.tsv", sep="\t", index=False)

agg = res.groupby(["module", "kind"]).smd.mean().reset_index()
null = agg[agg.kind == "random"].smd.abs()
thr = np.percentile(null, 95)
bio = agg[agg.kind == "biological"].copy()
bio["exceeds_null"] = bio.smd.abs() > thr
print(f"  random-module |SMD| 95th percentile = {thr:.4f}  (n={len(null)} sets)")
print(f"  matched cells: {res.n_hi.iloc[0]:,} hi / {res.n_lo.iloc[0]:,} lo per sample\n")
print(bio.sort_values("smd", key=np.abs, ascending=False)
      .to_string(index=False, float_format=lambda v: f"{v:+.4f}"))

pc = bio[bio.module == "Proliferation"]
print(f"\n  POSITIVE CONTROL (Proliferation): |SMD|="
      f"{abs(pc.smd.iloc[0]):.4f}  exceeds null: {bool(pc.exceeds_null.iloc[0])}")
if not bool(pc.exceeds_null.iloc[0]):
    print("  *** positive control fails -> the test is underpowered and its "
          "null is uninformative ***")
else:
    n_ex = int(bio.exceeds_null.sum())
    print(f"  modules exceeding the random null: {n_ex} of {len(bio)}")
bio.to_csv(f"{MR}/43c_residual_axis_verdict.tsv", sep="\t", index=False)
print("\nwritten: 43a-43c")
