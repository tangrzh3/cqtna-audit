#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 54  Power curves across six diseases.

Step 53 produced one calibrated curve, for melanoma. One curve is a case study;
the claim that candidate-list recovery is a general function of outcome power
needs several diseases with different genetic architectures. The five FinnGen
cancers already extracted in Step 30 span 3,139 to 24,270 cases, which covers
the informative range of the curve.

Method as in Step 53: down-sample the outcome in summary space,
    se_sim   = se_obs * sqrt(N_eff_obs / N_eff_sim)
    beta_sim = beta_obs + N(0, se_sim^2 - se_obs^2)
recompute the Wald-ratio MR, re-derive the FDR<0.05 list, and measure agreement
with that trait's own observed full-power list.

SCOPE, carried over from Step 53 and not relaxed here
  - curves stop at each trait's observed power; upward extrapolation needs a
    model of the true-effect distribution that this data cannot validate (a flat
    prior over-produced hits 25-fold, a global normal empirical-Bayes prior
    under-produced them 100-fold)
  - the reference list is not truth, only a better-powered list that is itself
    unstable, so recovery fractions are agreement with a moving target

Output: 54a (all curves), 54b (summary at matched recovery levels)
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

MR = r"D:/R_ex/MR"
EXT = os.path.join(MR, "cancer_extracts")
rng = np.random.default_rng(1)
N_REP = 200
FDR_THR = 0.05

# FinnGen R12 manifest counts, verified in Step 30
TRAITS = {
    "Melanoma":   dict(file=None, ncase=12530, nctrl=789099),
    "Lung":       dict(file="C3_BRONCHUS_LUNG_EXALLC.tsv", ncase=9639, nctrl=378749),
    "Colorectal": dict(file="C3_COLORECTAL_EXALLC.tsv", ncase=11790, nctrl=378749),
    "Pancreas":   dict(file="C3_PANCREAS_EXALLC.tsv", ncase=3139, nctrl=378749),
    "Breast":     dict(file="C3_BREAST_EXALLC.tsv", ncase=24270, nctrl=222078),
    "Prostate":   dict(file="C3_PROSTATE_EXALLC.tsv", ncase=20368, nctrl=156671),
}


def n_eff(nc, nk):
    return 4.0 / (1.0 / nc + 1.0 / nk)


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def wald(bo, seo, be):
    b = bo / be; se = seo / np.abs(be)
    return 2 * stats.norm.sf(np.abs(b / se))


# ---------------------------------------------------------------- exposures
ann = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
base = ann[["exposure", "SYMBOL", "SNP", "beta_exposure", "beta_outcome",
            "se_outcome"]].dropna()
base = base[(base.beta_exposure != 0) & (base.se_outcome > 0)].reset_index(drop=True)

harm = pd.read_csv(f"{MR}/01_harmonised_all.tsv", sep="\t")
# column names contain dots, which itertuples renames -- index the columns directly
alle = dict(zip(zip(harm["exposure"], harm["SNP"]),
                zip(harm["effect_allele.outcome"].astype(str).str.upper(),
                    harm["other_allele.outcome"].astype(str).str.upper())))


def load_trait(name):
    cfg = TRAITS[name]
    if cfg["file"] is None:                      # melanoma: already harmonised
        d = base.copy()
        return d, n_eff(cfg["ncase"], cfg["nctrl"])
    ext = pd.read_csv(os.path.join(EXT, cfg["file"]), sep="\t")
    out = {r.SNP: (str(r.ref).upper(), str(r.alt).upper(), r.beta, r.sebeta)
           for r in ext.itertuples()}
    rows = []
    for r in base.itertuples():
        g = out.get(r.SNP)
        if g is None:
            continue
        ref, alt, bo, seo = g
        a = alle.get((r.exposure, r.SNP))
        if a is None or seo <= 0:
            continue
        ea = a[0]
        if ea == alt:
            b = bo
        elif ea == ref:
            b = -bo
        else:
            continue
        rows.append((r.exposure, r.SYMBOL, r.beta_exposure, b, seo))
    d = pd.DataFrame(rows, columns=["exposure", "SYMBOL", "beta_exposure",
                                    "beta_outcome", "se_outcome"])
    return d, n_eff(cfg["ncase"], cfg["nctrl"])


all_rows, summary = [], []
for name, cfg in TRAITS.items():
    d, ne_obs = load_trait(name)
    if len(d) < 500:
        print(f"  [skip] {name}: only {len(d)} usable exposures")
        continue
    p_full = wald(d.beta_outcome.values, d.se_outcome.values,
                  d.beta_exposure.values)
    full = set(d.SYMBOL[bh(p_full) < FDR_THR])
    print(f"\n{name}: {len(d):,} exposures, {cfg['ncase']:,} cases, "
          f"N_eff {ne_obs:,.0f}, full-power list {len(full)} genes")
    if len(full) < 3:
        print("  full-power list too small for a recovery curve; skipped")
        continue

    grid = sorted({int(cfg["ncase"] * f) for f in
                   (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)})
    for nc in grid:
        ne = n_eff(nc, cfg["nctrl"])
        if nc >= cfg["ncase"]:
            all_rows.append(dict(trait=name, cases=nc, frac=1.0,
                                 mean_hits=len(full), sensitivity=1.0,
                                 jaccard=1.0, j_lo=1.0, j_hi=1.0))
            continue
        scale = np.sqrt(ne_obs / ne)
        se_s = d.se_outcome.values * scale
        extra = np.sqrt(np.maximum(se_s ** 2 - d.se_outcome.values ** 2, 0))
        hits, sens, jac = [], [], []
        for _ in range(N_REP):
            bo = d.beta_outcome.values + rng.normal(0, extra)
            g = set(d.SYMBOL[bh(wald(bo, se_s, d.beta_exposure.values)) < FDR_THR])
            hits.append(len(g))
            sens.append(len(g & full) / len(full))
            jac.append(len(g & full) / len(g | full) if (g | full) else np.nan)
        jac = np.array(jac, float)
        all_rows.append(dict(trait=name, cases=nc, frac=nc / cfg["ncase"],
                             mean_hits=np.mean(hits),
                             sensitivity=np.mean(sens),
                             jaccard=np.nanmean(jac),
                             j_lo=np.nanpercentile(jac, 2.5),
                             j_hi=np.nanpercentile(jac, 97.5)))
        print(f"    {nc:>7,} cases ({nc/cfg['ncase']:.0%})  hits "
              f"{np.mean(hits):6.1f}  recovers {100*np.mean(sens):5.1f}%")

    c = pd.DataFrame([r for r in all_rows if r["trait"] == name])
    for target in (0.5, 0.8):
        hit = c[c.sensitivity >= target]
        summary.append(dict(trait=name, obs_cases=cfg["ncase"],
                            full_list_genes=len(full), target=target,
                            cases_needed=int(hit.cases.min()) if len(hit) else np.nan,
                            frac_of_observed=(hit.cases.min() / cfg["ncase"]
                                              if len(hit) else np.nan)))

curve = pd.DataFrame(all_rows)
curve.to_csv(f"{MR}/54a_multitrait_power_curves.tsv", sep="\t", index=False)
summ = pd.DataFrame(summary)
summ.to_csv(f"{MR}/54b_power_thresholds.tsv", sep="\t", index=False)

print("\n" + "=" * 78)
print("recovery as a fraction of each trait's own observed power")
print("=" * 78)
piv = curve.pivot_table(index="frac", columns="trait", values="sensitivity")
print((100 * piv).round(0).to_string())

print("\n" + "=" * 78)
print("outcome size needed to recover 50% / 80% of the full-power list")
print("=" * 78)
print(summ.pivot_table(index="trait", columns="target",
                       values=["cases_needed", "frac_of_observed"]).to_string())
print("\nwritten: 54a, 54b")
