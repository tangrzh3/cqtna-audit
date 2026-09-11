#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 19b  Spatial follow-up: the fairest remaining test of the CD4 hypothesis.

Main analysis (23a-23f) showed TPI1 tracks the tumor/glycolysis compartment and
is NEGATIVELY correlated with the lymphoid compartment.  Before concluding, give
the CD4 hypothesis its best shot:

  (1) restrict to lymphocyte-dominant spots -- within T-cell-rich tissue, does
      more lymphocyte content mean more TPI1?
  (2) partial correlation controlling tumor score + sequencing depth
  (3) HLA-C run identically as a positive control (an immune gene must pass)
  (4) gene-by-gene correlations against individual T-cell / melanocyte markers
"""
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats

OUT = (sys.argv[1] if len(sys.argv) > 1
       else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
spots = pd.read_csv(f"{OUT}/23a_ST_spot_level.tsv", sep="\t")

# region labels reproduce the main script
comp_cols = ["tumor", "lymphoid", "stroma"]
z = spots.groupby("section")[comp_cols].transform(lambda s: (s - s.mean()) / s.std())
spots["region"] = z.idxmax(axis=1)
spots["region_margin"] = z.max(axis=1) - z.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
spots.loc[spots.region_margin < 0.25, "region"] = "mixed"


def partial_spearman(d, a, b, covars):
    """Spearman correlation of a and b after regressing both on covars (ranks)."""
    d = d[[a, b] + covars].dropna()
    if len(d) < 25:
        return np.nan, np.nan, len(d)
    R = d.rank()
    C = np.column_stack([np.ones(len(R))] + [R[c].values for c in covars])
    ra = R[a].values - C @ np.linalg.lstsq(C, R[a].values, rcond=None)[0]
    rb = R[b].values - C @ np.linalg.lstsq(C, R[b].values, rcond=None)[0]
    r, p = stats.pearsonr(ra, rb)
    return r, p, len(d)


def fisher_z_meta(rhos, ns):
    rhos, ns = np.asarray(rhos, float), np.asarray(ns, float)
    ok = np.isfinite(rhos) & (ns > 3)
    if ok.sum() == 0:
        return np.nan, np.nan, np.nan, np.nan
    zz = np.arctanh(np.clip(rhos[ok], -0.999999, 0.999999))
    w = ns[ok] - 3
    zbar = np.sum(w * zz) / np.sum(w)
    se = 1.0 / np.sqrt(np.sum(w))
    return (np.tanh(zbar), 2 * stats.norm.sf(abs(zbar / se)),
            np.tanh(zbar - 1.96 * se), np.tanh(zbar + 1.96 * se))


print("=" * 74)
print("(1)+(2)  WITHIN LYMPHOCYTE-DOMINANT SPOTS ONLY")
print("         partial Spearman  gene ~ lymphoid | tumor, depth")
print("=" * 74)

rows = []
for gene in ["TPI1", "HLA-C", "SMC2", "ZFYVE19", "SPSB2", "KIAA0040"]:
    for subset, label in [(spots.region == "lymphoid", "lymphoid-dominant spots"),
                          (slice(None), "all spots")]:
        sub = spots[subset] if not isinstance(subset, slice) else spots
        rs, ns = [], []
        for sec, d in sub.groupby("section"):
            r, p, n = partial_spearman(d, gene, "lymphoid", ["tumor", "log_total"])
            if np.isfinite(r):
                rs.append(r); ns.append(n)
        r, p, lo, hi = fisher_z_meta(rs, ns)
        rows.append(dict(gene=gene, subset=label, n_sections=len(rs), n_spots=int(np.sum(ns)),
                         rho_partial=r, ci_lo=lo, ci_hi=hi, p=p,
                         n_sections_pos=int(np.sum(np.array(rs) > 0))))
res1 = pd.DataFrame(rows)
res1.to_csv(f"{OUT}/23g_ST_within_lymphoid.tsv", sep="\t", index=False)
print(res1.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

print()
print("=" * 74)
print("(4)  CORRELATION WITH INDIVIDUAL MARKER GENES (meta over 8 sections)")
print("=" * 74)

markers = {"T cell": ["CD2", "CD3D", "CD3E", "CD247", "IL7R", "CD52", "PTPRC", "CD4"],
           "melanocyte": ["MLANA", "PMEL", "TYR", "DCT", "TYRP1", "SOX10"],
           "glycolysis": ["GAPDH", "PKM", "LDHA", "ENO1", "ALDOA"]}

spot_gene = pd.read_csv(f"{OUT}/23a_ST_spot_level.tsv", sep="\t", nrows=1).columns
rows = []
import glob, os, re
ST_DIR = "D:/Downloads/ST-Melanoma-Datasets_1"
for path in sorted(glob.glob(os.path.join(ST_DIR, "ST_mel*_counts.tsv"))):
    name = os.path.basename(path).replace("_counts.tsv", "").replace("ST_", "")
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = [str(i).split(" ")[0] for i in df.index]
    df = df[~df.index.duplicated(keep="first")]
    total = df.sum(axis=0)
    keep = (total >= 500) & ((df > 0).sum(axis=0) >= 200)
    df, total = df.loc[:, keep], total[keep]
    ln = np.log1p(df.div(total, axis=1) * 1e4)
    for cls, gs in markers.items():
        for m in gs:
            if m not in ln.index or "TPI1" not in ln.index:
                continue
            r, p = stats.spearmanr(ln.loc["TPI1"], ln.loc[m])
            rows.append(dict(section=name, marker=m, marker_class=cls, n=ln.shape[1], rho=r, p=p))
mg = pd.DataFrame(rows)
summ = []
for m, d in mg.groupby("marker"):
    r, p, lo, hi = fisher_z_meta(d.rho.values, d.n.values)
    summ.append(dict(marker=m, marker_class=d.marker_class.iloc[0], rho_meta=r,
                     ci_lo=lo, ci_hi=hi, p=p, n_sections_pos=int((d.rho > 0).sum())))
summ = pd.DataFrame(summ).sort_values(["marker_class", "rho_meta"], ascending=[True, False])
summ.to_csv(f"{OUT}/23h_ST_TPI1_vs_markers.tsv", sep="\t", index=False)
print(summ.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

print()
print("=" * 74)
print("SUMMARY")
print("=" * 74)
t = res1[(res1.gene == "TPI1") & (res1.subset == "lymphoid-dominant spots")].iloc[0]
h = res1[(res1["gene"] == "HLA-C") & (res1.subset == "lymphoid-dominant spots")].iloc[0]
print(f"  TPI1  within lymphoid spots: rho={t.rho_partial:+.3f} p={t.p:.3g} "
      f"({t.n_sections_pos}/{t.n_sections} sections positive)")
print(f"  HLA-C within lymphoid spots: rho={h.rho_partial:+.3f} p={h.p:.3g} "
      f"({h.n_sections_pos}/{h.n_sections} sections positive)  [positive control]")
print("\nwritten: 23g, 23h")
