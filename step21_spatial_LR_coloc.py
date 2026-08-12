#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 21  Spatial validation of the CellChat-nominated interactions.

WHAT THIS DOES *NOT* DO -- and why
----------------------------------
It does not look for "TPI1-high CD4 spots". It cannot: legacy ST spots are
100 um wide and hold 10-40 mixed cells, and Step 19 already showed that spot-
level TPI1 is dominated by tumour glycolysis (rho=+0.17 with the glycolysis
module, 8/8 sections; rho=-0.08 with the lymphoid compartment). Any analysis
keyed on "TPI1 in a spot" is therefore a tumour readout, not a CD4 readout.

WHAT IT DOES INSTEAD
--------------------
CellChat nominates ligand-receptor pairs on the glycolysis-high CD4 axis (and,
as a sensitivity run, the TPI1-high axis). The spatially checkable claim is not
about TPI1 or the module -- it is that the LIGAND and the RECEPTOR of that axis
occupy neighbouring tissue. So we test, for each pair, whether ligand expression
predicts receptor expression in ADJACENT spots (bivariate Moran's I on the ST
lattice).

TWO NULLS, because one is not enough
------------------------------------
(1) spatial permutation -- shuffle the receptor map. Tests "co-localised at
    all". Almost everything passes this, because tissue architecture makes
    any two abundant genes co-vary. Necessary but not sufficient.
(2) *** expression-matched random gene-pair null *** -- draw random gene pairs
    matched on mean expression decile, compute the same statistic. Tests
    "co-localised MORE than an arbitrary pair of equally abundant genes".
    This is the null that carries the claim. Report this one.

RESOLUTION CAVEAT -- must go in the manuscript
---------------------------------------------
Legacy ST spot pitch is ~200 um. This tests tissue-scale co-localisation, NOT
cell-cell contact. Interactions annotated "Cell-Cell Contact" in CellChatDB
are NOT validatable at this resolution; they are tagged and reported
separately, and must not be claimed as spatially confirmed.

Usage : python step21_spatial_LR_coloc.py [LR_file] [tag]
        defaults to 24f_Glyco_top_LR_for_spatial.tsv (from step20)
        run again with 24f_TPI1_top_LR_for_spatial.tsv for the sensitivity axis
Output: 26a/26b/26d tagged per run; 26c (compartment co-localisation) is
        run-independent and always overwritten with the same values
"""
import os
import re
import sys
import glob
import numpy as np
import pandas as pd
from scipy import stats

ST_DIR = "D:/Downloads/ST-Melanoma-Datasets_1"
MR = "D:/R_ex/MR"

# usage: python step21_spatial_LR_coloc.py [LR_file] [output_tag]
#   default is the glycolysis-module split from step20
LR_FILE = sys.argv[1] if len(sys.argv) > 1 else \
    os.path.join(MR, "24f_Glyco_top_LR_for_spatial.tsv")
if not os.path.isabs(LR_FILE):
    LR_FILE = os.path.join(MR, LR_FILE)
TAG = sys.argv[2] if len(sys.argv) > 2 else \
    (re.search(r"24f_(\w+?)_top_LR", os.path.basename(LR_FILE)).group(1)
     if re.search(r"24f_(\w+?)_top_LR", os.path.basename(LR_FILE)) else "run")
print(f"LR file : {LR_FILE}\noutput tag: {TAG}\n")

N_PERM = 999          # spatial permutations
N_MATCHED = 500       # matched random gene pairs
NEIGHBOR_DIST = 1.5   # grid units; <=1.5 gives the 8-neighbourhood
MIN_DETECT = 0.05     # gene must be detected in >=5% of spots
rng = np.random.default_rng(1)


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    adj = np.empty(n)
    adj[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def load_section(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = [str(i).split(" ")[0] for i in df.index]
    df = df[~df.index.duplicated(keep="first")]
    total = df.sum(axis=0)
    keep = (total >= 500) & ((df > 0).sum(axis=0) >= 200)
    df, total = df.loc[:, keep], total[keep]
    ln = np.log1p(df.div(total, axis=1) * 1e4)
    xy = np.array([[float(m.group(1)), float(m.group(2))]
                   for m in (re.match(r"^([\d.]+)x([\d.]+)$", str(s)) for s in df.columns)])
    return ln, xy


def neighbour_weights(xy, dist=NEIGHBOR_DIST):
    """Row-standardised binary adjacency on the ST lattice."""
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
    W = ((d > 0) & (d <= dist)).astype(float)
    rs = W.sum(1, keepdims=True)
    rs[rs == 0] = 1
    return W / rs


def zmat(ln, genes):
    """genes x spots, z-scored per gene (rows with zero variance -> zeros)."""
    X = ln.loc[genes].to_numpy(dtype=float)
    sd = X.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    return (X - X.mean(axis=1, keepdims=True)) / sd


def members(x):
    """CellChatDB encodes complexes as A_B; split into member genes."""
    return [g for g in re.split(r"[_\-]", str(x)) if g and g.lower() != "nan"]


# ==========================================================================
if not os.path.exists(LR_FILE):
    raise SystemExit(f"missing {LR_FILE} -- run step20_cellchat_GSE120575.R first")

lr = pd.read_csv(LR_FILE, sep="\t")
pairs = []
for _, r in lr.iterrows():
    for L in members(r.get("ligand")):
        for R in members(r.get("receptor")):
            pairs.append(dict(interaction_name=r["interaction_name"],
                              pathway=r.get("pathway_name"),
                              annotation=r.get("annotation"),
                              direction=r.get("direction"),
                              partner=r.get("partner"),
                              ligand=L, receptor=R))
pairs = pd.DataFrame(pairs).drop_duplicates(subset=["ligand", "receptor"])
print(f"{len(lr)} CellChat rows -> {len(pairs)} unique ligand/receptor gene pairs\n")

sections = sorted(glob.glob(os.path.join(ST_DIR, "ST_mel*_counts.tsv")))
rows = []
# many immune ligands/receptors (chemokine receptors, IFNG, ...) sit below the
# detection floor of 100 um ST. Track them explicitly -- a pair that was never
# testable is not a pair that failed validation, and the two must not be merged.
untestable = {}

for path in sections:
    name = os.path.basename(path).replace("_counts.tsv", "").replace("ST_", "")
    ln, xy = load_section(path)
    W = neighbour_weights(xy)
    n = ln.shape[1]

    detect = (ln > 0).mean(axis=1)
    testable = list(detect[detect >= MIN_DETECT].index)
    gidx = {g: i for i, g in enumerate(testable)}

    # Precompute once per section:
    #   Z   = z-scored expression        (g x n)
    #   WZ  = (W @ z) for every gene     (g x n)   -> bivariate Moran = Z[a] . WZ[b] / n
    #   WtZ = (W.T @ z) for every gene   (g x n)   -> permutation null in O(n)
    Z = zmat(ln, testable)
    WZ = Z @ W.T          # row b of WZ is (W @ z_b)
    WtZ = Z @ W           # row a of WtZ is (W.T @ z_a)

    mean_expr = ln.loc[testable].mean(axis=1).to_numpy()
    mbin = pd.qcut(pd.Series(mean_expr), 10, labels=False, duplicates="drop").to_numpy()
    bin_members = {b: np.where(mbin == b)[0] for b in np.unique(mbin)}

    print(f"{name}: {n} spots, {len(testable)} testable genes, "
          f"mean neighbours {(W > 0).sum(1).mean():.1f}")

    for _, p in pairs.iterrows():
        L, R = p.ligand, p.receptor
        if L not in gidx or R not in gidx:
            miss = [g for g in (L, R) if g not in gidx]
            untestable.setdefault((L, R), set()).update(miss)
            continue
        li, ri = gidx[L], gidx[R]
        I = float(Z[li] @ WZ[ri] / n)

        # null (1): spatial permutation of the receptor map -- O(n) per draw
        u = WtZ[li]
        zr = Z[ri]
        perm = np.array([u @ rng.permutation(zr) for _ in range(N_PERM)]) / n
        p_spatial = (1 + (perm >= I).sum()) / (1 + N_PERM)

        # null (2): expression-matched random gene pairs -- O(n) per draw
        la = bin_members[mbin[li]]
        rb = bin_members[mbin[ri]]
        aa = rng.choice(la, N_MATCHED)
        bb = rng.choice(rb, N_MATCHED)
        mn = np.einsum("ij,ij->i", Z[aa], WZ[bb]) / n
        p_matched = (1 + (mn >= I).sum()) / (1 + N_MATCHED)
        z_matched = (I - mn.mean()) / mn.std() if mn.std() > 0 else np.nan

        rows.append(dict(section=name, n_spots=n, ligand=L, receptor=R,
                         interaction_name=p.interaction_name, pathway=p.pathway,
                         annotation=p.annotation, direction=p.direction,
                         partner=p.partner, moran_I=I,
                         p_spatial=p_spatial, p_matched=p_matched,
                         z_vs_matched=z_matched,
                         null_matched_mean=mn.mean(), null_matched_sd=mn.std()))

res = pd.DataFrame(rows)
res.to_csv(f"{MR}/26a_{TAG}_LR_coloc_per_section.tsv", sep="\t", index=False)

tested = set(map(tuple, res[["ligand", "receptor"]].values)) if len(res) else set()
never = {k: v for k, v in untestable.items() if k not in tested}
if never:
    nt = pd.DataFrame([dict(ligand=L, receptor=R, below_detection=",".join(sorted(g)))
                       for (L, R), g in never.items()])
    nt.to_csv(f"{MR}/26d_{TAG}_untestable_pairs.tsv", sep="\t", index=False)
    print(f"\n=== NOT TESTABLE: {len(never)} of {len(pairs)} pairs "
          f"(gene below {MIN_DETECT:.0%} detection in every section) ===")
    print("    these are absent from the ST assay, NOT refuted by it")
    print(nt.to_string(index=False))


def fisher_combine(p):
    p = np.clip(np.asarray(p, float), 1e-12, 1)
    return stats.chi2.sf(-2 * np.log(p).sum(), 2 * len(p))


summ = []
for (L, R), d in res.groupby(["ligand", "receptor"]):
    summ.append(dict(
        ligand=L, receptor=R,
        interaction_name=d.interaction_name.iloc[0], pathway=d.pathway.iloc[0],
        annotation=d.annotation.iloc[0], direction=d.direction.iloc[0],
        partner=d.partner.iloc[0], n_sections=len(d),
        moran_I_mean=d.moran_I.mean(),
        z_vs_matched_mean=d.z_vs_matched.mean(),
        n_sections_above_matched=int((d.p_matched < 0.05).sum()),
        p_spatial_combined=fisher_combine(d.p_spatial),
        p_matched_combined=fisher_combine(d.p_matched)))
summ = pd.DataFrame(summ)

if len(summ):
    summ["FDR_matched"] = bh_fdr(summ.p_matched_combined)
    summ = summ.sort_values("p_matched_combined")
    summ.to_csv(f"{MR}/26b_{TAG}_LR_coloc_summary.tsv", sep="\t", index=False)

    contact = summ.annotation.astype(str).str.contains("Contact", case=False, na=False)
    show = ["ligand", "receptor", "pathway", "partner", "moran_I_mean",
            "z_vs_matched_mean", "n_sections_above_matched", "FDR_matched"]

    print("\n=== validated ABOVE the expression-matched null (FDR<0.05, secreted/ECM) ===")
    ok = summ[(summ.FDR_matched < 0.05) & (~contact)]
    print(ok.head(25)[show].to_string(index=False, float_format=lambda v: f"{v:.3g}")
          if len(ok) else "    none")

    print("\n=== NOT validatable at 200 um pitch (Cell-Cell Contact) ===")
    cc = summ[contact]
    print(cc.head(15)[["ligand", "receptor", "pathway", "z_vs_matched_mean"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3g}")
          if len(cc) else "    none")

    weak = summ[(summ.p_spatial_combined < 0.05) & (summ.p_matched_combined >= 0.05)]
    print(f"\n=== passes spatial permutation but NOT the matched null: "
          f"{len(weak)} of {len(summ)} pairs ===")
    print("    explained by tissue architecture alone -- do not claim these")

# ---- compartment co-localisation, as context ----------------------------
MARKERS = {
    "lymphoid": ["CD2", "CD3D", "CD3E", "CD3G", "CD247", "TRAC", "LCK", "IL7R",
                 "CD52", "CCL5", "CD8A", "CD4", "SKAP1", "CD27"],
    "myeloid": ["LYZ", "CD68", "CD14", "AIF1", "ITGAX", "TYROBP", "FCER1G"],
    "bcell": ["MS4A1", "CD79A", "CD79B", "IGHM", "BANK1"],
    "tumor": ["MLANA", "PMEL", "TYR", "DCT", "TYRP1", "SOX10", "MITF", "S100B", "PRAME"],
}
ct_rows = []
for path in sections:
    name = os.path.basename(path).replace("_counts.tsv", "").replace("ST_", "")
    ln, xy = load_section(path)
    W = neighbour_weights(xy)
    n = ln.shape[1]
    S = {}
    for k, gs in MARKERS.items():
        av = [g for g in gs if g in ln.index]
        s = zmat(ln, av).mean(axis=0)
        S[k] = (s - s.mean()) / (s.std() if s.std() > 0 else 1)
    keys = list(S)
    for i, a in enumerate(keys):
        for b in keys[i:]:
            I = float(S[a] @ (W @ S[b]) / n)
            u = W.T @ S[a]
            perm = np.array([u @ rng.permutation(S[b]) for _ in range(N_PERM)]) / n
            ct_rows.append(dict(section=name, comp_a=a, comp_b=b, moran_I=I,
                                p=(1 + (perm >= I).sum()) / (1 + N_PERM)))
ct = pd.DataFrame(ct_rows)
ct.to_csv(f"{MR}/26c_ST_compartment_colocalization.tsv", sep="\t", index=False)
print("\n=== compartment co-localisation (mean bivariate Moran's I) ===")
print(ct.groupby(["comp_a", "comp_b"])
        .agg(I=("moran_I", "mean"), n_sig=("p", lambda s: int((s < 0.05).sum())))
        .sort_values("I", ascending=False)
        .to_string(float_format=lambda v: f"{v:.3f}"))

print("\nwritten: 26a-26c")

