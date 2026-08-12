#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 19  Spatial transcriptomics (Thrane et al. 2018, legacy ST 100um spots)

Question: does TPI1 in melanoma tissue track the TUMOR compartment (Warburg
glycolysis) or the LYMPHOCYTE-INFILTRATED compartment?  This is the direct
falsification test for Limitation #2 of the manuscript.

Data: D:/Downloads/ST-Melanoma-Datasets_1/ST_mel{1-4}_rep{1,2}_counts.tsv
      rows = "SYMBOL ENSG...", columns = spots named "<x>x<y>"
"""
import os
import re
import glob
import numpy as np
import pandas as pd
from scipy import stats


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    adj[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def ols(y, X):
    """OLS with intercept. X: DataFrame. Returns (params, pvalues, r2) as Series."""
    Xm = np.column_stack([np.ones(len(X)), X.values])
    yv = np.asarray(y, float)
    beta, *_ = np.linalg.lstsq(Xm, yv, rcond=None)
    resid = yv - Xm @ beta
    dof = len(yv) - Xm.shape[1]
    s2 = resid @ resid / dof
    se = np.sqrt(np.diag(s2 * np.linalg.pinv(Xm.T @ Xm)))
    tval = beta / se
    pv = 2 * stats.t.sf(np.abs(tval), dof)
    names = ["const"] + list(X.columns)
    sst = ((yv - yv.mean()) ** 2).sum()
    r2 = 1 - (resid @ resid) / sst
    return pd.Series(beta, index=names), pd.Series(pv, index=names), r2, resid

ST_DIR = "D:/Downloads/ST-Melanoma-Datasets_1"
OUT = "D:/R_ex/MR"

MIN_SPOT_COUNTS = 500        # legacy ST spots; Thrane used ~500 as usable depth
MIN_SPOT_GENES = 200

CANDIDATES = ["TPI1", "SMC2", "ZFYVE19", "SPSB2", "HLA-C", "KIAA0040"]

# marker panels ------------------------------------------------------------
MARKERS = {
    "tumor": ["MLANA", "PMEL", "TYR", "DCT", "TYRP1", "SOX10", "MITF",
              "S100B", "MAGEA6", "PRAME"],
    "lymphoid": ["CD2", "CD3D", "CD3E", "CD3G", "CD247", "TRAC", "LCK",
                 "IL7R", "CD52", "CCL5", "CD8A", "CD4", "SKAP1", "CD27"],
    "bcell": ["MS4A1", "CD79A", "CD79B", "IGHM", "BANK1"],
    "myeloid": ["LYZ", "CD68", "CD14", "AIF1", "ITGAX", "TYROBP", "FCER1G"],
    "stroma": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1"],
    # glycolysis WITHOUT TPI1 -- used to ask whether TPI1 is separable from
    # the module it belongs to
    "glycolysis": ["GAPDH", "PKM", "LDHA", "ENO1", "ALDOA", "PGK1", "PGAM1",
                   "HK1", "HK2", "GPI", "SLC2A1"],
    "prolif": ["MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA"],
}


def load_section(path):
    """Return (lognorm DataFrame genes x spots, spot metadata DataFrame)."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    # index looks like "ANXA2 ENSG00000182718" -> take symbol
    df.index = [str(i).split(" ")[0] for i in df.index]
    df = df[~df.index.duplicated(keep="first")]

    total = df.sum(axis=0)
    ngene = (df > 0).sum(axis=0)
    keep = (total >= MIN_SPOT_COUNTS) & (ngene >= MIN_SPOT_GENES)
    df = df.loc[:, keep]
    total = total[keep]
    ngene = ngene[keep]

    # CPM + log1p
    ln = np.log1p(df.div(total, axis=1) * 1e4)

    coords = []
    for s in df.columns:
        m = re.match(r"^([\d.]+)x([\d.]+)$", str(s))
        coords.append((float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan))
    meta = pd.DataFrame(coords, index=df.columns, columns=["x", "y"])
    meta["total_counts"] = total
    meta["n_genes"] = ngene
    meta["log_total"] = np.log10(total)
    return ln, meta


def score(ln, genes):
    """Mean of per-section z-scored log-expression over available markers."""
    avail = [g for g in genes if g in ln.index]
    if not avail:
        return None, []
    sub = ln.loc[avail]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0, skipna=True), avail


def fisher_z_meta(rhos, ns):
    """Random-effects-free (fixed) Fisher z meta of Spearman rho."""
    rhos = np.asarray(rhos, float)
    ns = np.asarray(ns, float)
    ok = np.isfinite(rhos) & (ns > 3)
    z = np.arctanh(np.clip(rhos[ok], -0.999999, 0.999999))
    w = ns[ok] - 3
    zbar = np.sum(w * z) / np.sum(w)
    se = 1.0 / np.sqrt(np.sum(w))
    p = 2 * stats.norm.sf(abs(zbar / se))
    return np.tanh(zbar), p, np.tanh(zbar - 1.96 * se), np.tanh(zbar + 1.96 * se)


# ==========================================================================
sections = sorted(glob.glob(os.path.join(ST_DIR, "ST_mel*_counts.tsv")))
print(f"{len(sections)} sections\n")

spot_tables, corr_rows, reg_rows, region_rows = [], [], [], []

for path in sections:
    name = os.path.basename(path).replace("_counts.tsv", "").replace("ST_", "")
    patient = name.split("_")[0]
    ln, meta = load_section(path)

    for key, genes in MARKERS.items():
        sc, avail = score(ln, genes)
        meta[key] = sc
        if key == "tumor":
            print(f"{name}: {ln.shape[1]} spots kept, tumor markers {len(avail)}/{len(genes)}")

    for g in CANDIDATES:
        meta[g] = ln.loc[g] if g in ln.index else np.nan

    meta["section"] = name
    meta["patient"] = patient
    spot_tables.append(meta.reset_index().rename(columns={"index": "spot"}))

    # ---- per-section Spearman: candidate vs compartment scores -----------
    for g in CANDIDATES:
        if meta[g].isna().all():
            continue
        for comp in ["tumor", "lymphoid", "bcell", "myeloid", "stroma",
                     "glycolysis", "prolif"]:
            rho, p = stats.spearmanr(meta[g], meta[comp], nan_policy="omit")
            corr_rows.append(dict(section=name, patient=patient, gene=g,
                                  compartment=comp, n_spots=meta[g].notna().sum(),
                                  rho=rho, p=p))

    # ---- per-section regression: gene ~ tumor + lymphoid (+ depth) -------
    for g in CANDIDATES:
        if meta[g].isna().all():
            continue
        X = meta[["tumor", "lymphoid", "log_total"]].copy()
        X = (X - X.mean()) / X.std()
        y = (meta[g] - meta[g].mean()) / meta[g].std()
        d = pd.concat([y.rename("y"), X], axis=1).dropna()
        if len(d) < 30:
            continue
        bet, pv, r2, _ = ols(d["y"], d[["tumor", "lymphoid", "log_total"]])
        reg_rows.append(dict(section=name, patient=patient, gene=g, n=len(d),
                             beta_tumor=bet["tumor"], p_tumor=pv["tumor"],
                             beta_lymphoid=bet["lymphoid"], p_lymphoid=pv["lymphoid"],
                             beta_depth=bet["log_total"], r2=r2))

spots = pd.concat(spot_tables, ignore_index=True)
spots.to_csv(f"{OUT}/23a_ST_spot_level.tsv", sep="\t", index=False)
corr = pd.DataFrame(corr_rows)
corr.to_csv(f"{OUT}/23b_ST_correlations_per_section.tsv", sep="\t", index=False)
reg = pd.DataFrame(reg_rows)
reg.to_csv(f"{OUT}/23c_ST_regression_per_section.tsv", sep="\t", index=False)
print(f"\ntotal spots after QC: {len(spots)}")

# ---- meta-analysis across the 8 sections --------------------------------
meta_rows = []
for g in CANDIDATES:
    for comp in ["tumor", "lymphoid", "bcell", "myeloid", "stroma",
                 "glycolysis", "prolif"]:
        sub = corr[(corr.gene == g) & (corr.compartment == comp)]
        if sub.empty:
            continue
        r, p, lo, hi = fisher_z_meta(sub.rho.values, sub.n_spots.values)
        meta_rows.append(dict(gene=g, compartment=comp, n_sections=len(sub),
                              rho_meta=r, ci_lo=lo, ci_hi=hi, p=p,
                              n_sections_pos=(sub.rho > 0).sum()))
mt = pd.DataFrame(meta_rows)
mt["FDR"] = bh_fdr(mt.p)
mt.to_csv(f"{OUT}/23d_ST_meta_correlations.tsv", sep="\t", index=False)

print("\n=== meta Spearman across 8 sections ===")
print(mt.sort_values(["gene", "compartment"]).to_string(index=False,
      float_format=lambda v: f"{v:.3g}"))

print("\n=== joint regression (standardised betas, per section) ===")
if not reg.empty:
    print(reg.groupby("gene")[["beta_tumor", "beta_lymphoid", "r2"]]
             .agg(["mean", "min", "max"]).to_string(float_format=lambda v: f"{v:.3f}"))

# ---- region assignment and TPI1 by region -------------------------------
comp_cols = ["tumor", "lymphoid", "stroma"]
z = spots.groupby("section")[comp_cols].transform(lambda s: (s - s.mean()) / s.std())
spots["region"] = z.idxmax(axis=1)
spots["region_margin"] = z.max(axis=1) - z.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
spots.loc[spots.region_margin < 0.25, "region"] = "mixed"

print("\n=== spots per region ===")
print(spots.region.value_counts().to_string())

reg_stat = []
for g in CANDIDATES:
    if g not in spots.columns or spots[g].isna().all():
        continue
    # within-section centring so region contrast is not driven by section
    v = spots.groupby("section")[g].transform(lambda s: s - s.mean())
    for a, b in [("tumor", "lymphoid"), ("tumor", "stroma"), ("lymphoid", "stroma")]:
        va, vb = v[spots.region == a].dropna(), v[spots.region == b].dropna()
        if len(va) < 10 or len(vb) < 10:
            continue
        u, p = stats.mannwhitneyu(va, vb)
        reg_stat.append(dict(gene=g, group_a=a, group_b=b, n_a=len(va), n_b=len(vb),
                             median_a=va.median(), median_b=vb.median(),
                             diff=va.median() - vb.median(), p=p))
rs = pd.DataFrame(reg_stat)
if not rs.empty:
    rs["FDR"] = bh_fdr(rs.p)
    rs.to_csv(f"{OUT}/23e_ST_region_contrasts.tsv", sep="\t", index=False)
    print("\n=== candidate expression by region (section-centred) ===")
    print(rs.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

# ---- decisive test: is TPI1 separable from the glycolysis module? --------
print("\n=== TPI1 residual after removing glycolysis module ===")
res_rows = []
for name, d in spots.groupby("section"):
    d = d.dropna(subset=["TPI1", "glycolysis", "lymphoid", "tumor", "log_total"])
    if len(d) < 30:
        continue
    _, _, _, resid = ols(d["TPI1"], d[["glycolysis", "log_total"]])
    r_l, p_l = stats.spearmanr(resid, d["lymphoid"])
    r_t, p_t = stats.spearmanr(resid, d["tumor"])
    res_rows.append(dict(section=name, n=len(d), rho_resid_lymphoid=r_l, p_lymphoid=p_l,
                         rho_resid_tumor=r_t, p_tumor=p_t))
rr = pd.DataFrame(res_rows)
rr.to_csv(f"{OUT}/23f_ST_TPI1_residual.tsv", sep="\t", index=False)
print(rr.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
if not rr.empty:
    for col in ["rho_resid_lymphoid", "rho_resid_tumor"]:
        r, p, lo, hi = fisher_z_meta(rr[col].values, rr.n.values)
        print(f"  meta {col}: rho={r:.3f} [{lo:.3f},{hi:.3f}] p={p:.3g}")

print("\nwritten: 23a-23f")

