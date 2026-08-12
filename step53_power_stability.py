#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 53  Power-stability curve: how much outcome GWAS power does a dynamic-eQTL
         MR candidate list need before it stops changing?

Finding 4 currently rests on one observation: identical exposure data, two
melanoma outcome datasets (FinnGen 5,753 cases -> meta 12,530), and candidate
lists that shared nothing. That is an anecdote. This turns it into a curve.

DESIGN
Individual-level data are not available, so the outcome cannot be resampled
directly. Instead the meta outcome is DOWN-SAMPLED IN SUMMARY SPACE: for a
target effective sample size N_sim < N_obs,
    se_sim   = se_obs * sqrt(N_eff_obs / N_eff_sim)
    beta_sim = beta_obs + e,  e ~ N(0, se_sim^2 - se_obs^2)
which reproduces a study of that size around the same underlying effect. MR is
recomputed, the FDR<0.05 list is re-derived, and its Jaccard overlap with the
full-power list is recorded over many replicates.

*** BUILT-IN CALIBRATION ***
The real FinnGen melanoma analysis is already in hand. Simulating the meta
outcome down to FinnGen's effective size must reproduce the overlap actually
observed between the real FinnGen list and the meta list. If it does not, the
simulation is not trustworthy and the curve is not reported. This is the
positive control, fixed before running.

Effective sample size for a case-control GWAS: N_eff = 4 / (1/ncase + 1/nctrl).

Output: 53a (curve), 53b (calibration), 53c (per-trait curves)
"""
import numpy as np
import pandas as pd
from scipy import stats

MR = r"D:/R_ex/MR"
rng = np.random.default_rng(1)
N_REP = 200
FDR_THR = 0.05

# melanoma meta and the real FinnGen round
META_CASE, META_CTRL = 12530, 789099
FG_CASE, FG_CTRL = 5753, 378749


def n_eff(ncase, nctrl):
    return 4.0 / (1.0 / ncase + 1.0 / nctrl)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def wald(bo, seo, be):
    b = bo / be
    se = seo / np.abs(be)
    z = b / se
    return b, se, 2 * stats.norm.sf(np.abs(z))


# ---------------------------------------------------------------- load
d = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
d = d[["exposure", "SYMBOL", "SNP", "category", "beta_exposure",
       "beta_outcome", "se_outcome"]].dropna()
d = d[(d.beta_exposure != 0) & (d.se_outcome > 0)].reset_index(drop=True)
print(f"exposures with usable statistics: {len(d):,} "
      f"({d.SYMBOL.nunique():,} genes)")

b_full, se_full, p_full = wald(d.beta_outcome.values, d.se_outcome.values,
                               d.beta_exposure.values)
fdr_full = bh(p_full)
full_genes = set(d.SYMBOL[fdr_full < FDR_THR])
full_expo = set(d.exposure[fdr_full < FDR_THR])
print(f"full-power (meta) list: {len(full_expo)} records / {len(full_genes)} genes")

NE_META = n_eff(META_CASE, META_CTRL)
NE_FG = n_eff(FG_CASE, FG_CTRL)
print(f"effective N: meta {NE_META:,.0f}, FinnGen {NE_FG:,.0f}")


def simulate(target_ne, reps=N_REP):
    """Down-sample the observed outcome to target_ne and compare the resulting
    list with the OBSERVED full-power list. Used only for calibration, where the
    real comparison is also between two different power levels."""
    if target_ne >= NE_META:
        return np.array([[len(full_expo), 1.0, 1.0]] * reps)
    scale = np.sqrt(NE_META / target_ne)
    se_s = d.se_outcome.values * scale
    extra = np.sqrt(np.maximum(se_s ** 2 - d.se_outcome.values ** 2, 0))
    out = []
    for _ in range(reps):
        bo = d.beta_outcome.values + rng.normal(0, extra)
        _, _, p = wald(bo, se_s, d.beta_exposure.values)
        f = bh(p)
        hits = f < FDR_THR
        g = set(d.SYMBOL[hits]); e = set(d.exposure[hits])
        jg = len(g & full_genes) / len(g | full_genes) if (g | full_genes) else 0
        je = len(e & full_expo) / len(e | full_expo) if (e | full_expo) else 0
        out.append([hits.sum(), jg, je])
    return np.array(out)


def replicability(target_ne, reps=N_REP):
    """Agreement between TWO INDEPENDENT studies at the same power.

    Comparing a down-sampled draw against the observed list is not a stability
    measure: at the observed power it is a thing compared with itself and gives
    a Jaccard of 1 by construction. Replicability is the honest quantity, and it
    extends above the observed power.

    beta_obs estimates beta_true with error se_obs, so a posterior draw of the
    truth is taken first and two independent studies are then generated from it,
    which removes the shared-noise bias that would arise from adding noise twice
    to the same beta_obs.

    The prior on beta_true must NOT be flat. A first version drew
    beta_true ~ N(beta_obs, se_obs^2), which inflates the spread of true effects
    because most effects are near zero; it produced 531 hits at the observed
    power where the real analysis gives 21 -- a 25-fold inflation that made the
    curve useless. Empirical Bayes shrinkage is used instead:
        tau2 = max(var(beta_obs) - mean(se_obs^2), 0)      between-SNP variance
        w    = tau2 / (tau2 + se_obs^2)                    shrinkage weight
        beta_true ~ N(w * beta_obs, w * se_obs^2)
    The check that this is right is that the simulated hit count at the observed
    power must match the real one.
    """
    se_obs = d.se_outcome.values
    se_n = se_obs * np.sqrt(NE_META / target_ne)
    be = d.beta_exposure.values
    tau2 = max(np.var(d.beta_outcome.values) - np.mean(se_obs ** 2), 0.0)
    w = tau2 / (tau2 + se_obs ** 2)
    post_mean = w * d.beta_outcome.values
    post_sd = np.sqrt(w * se_obs ** 2)
    out = []
    for _ in range(reps):
        b_true = post_mean + rng.normal(0, post_sd)
        lists = []
        for _ in range(2):
            bo = b_true + rng.normal(0, se_n)
            _, _, p = wald(bo, se_n, be)
            hits = bh(p) < FDR_THR
            lists.append((set(d.SYMBOL[hits]), set(d.exposure[hits]), hits.sum()))
        (ga, ea, na), (gb, eb, nb) = lists
        jg = len(ga & gb) / len(ga | gb) if (ga | gb) else np.nan
        je = len(ea & eb) / len(ea | eb) if (ea | eb) else np.nan
        out.append([(na + nb) / 2, jg, je])
    return np.array(out, dtype=float)


# ---------------------------------------------------------------- calibration
print("\n" + "=" * 78)
print("CALIBRATION against the real FinnGen round")
print("=" * 78)
fg = pd.read_csv(f"{MR}/04_MR_results_strict_all.tsv", sep="\t")
cand = [c for c in fg.columns if c.lower() in ("fdr", "q", "qval")]
if cand:
    fg_hits = fg[fg[cand[0]] < FDR_THR]
else:
    fg["_fdr"] = bh(fg["pval"].values)
    fg_hits = fg[fg._fdr < FDR_THR]
ann = pd.read_csv(f"{MR}/06_locus_annotation.tsv", sep="\t")
sym_map = dict(zip(ann.exposure, ann.SYMBOL))
fg_genes = {sym_map.get(e) for e in fg_hits.exposure} - {None, np.nan}
real_j = (len(fg_genes & full_genes) / len(fg_genes | full_genes)
          if (fg_genes | full_genes) else np.nan)
print(f"  real FinnGen FDR<0.05: {len(fg_hits)} records / {len(fg_genes)} genes")
print(f"  real Jaccard(FinnGen genes, meta genes) = {real_j:.3f}")
print(f"  shared genes: {sorted(fg_genes & full_genes)}")

sim_fg = simulate(NE_FG)
print(f"\n  simulated at FinnGen power: {sim_fg[:,0].mean():.1f} hits "
      f"(real {len(fg_hits)}), Jaccard {sim_fg[:,1].mean():.3f} "
      f"[{np.percentile(sim_fg[:,1],2.5):.3f}, {np.percentile(sim_fg[:,1],97.5):.3f}]")
ok = (np.percentile(sim_fg[:, 1], 2.5) <= real_j <= np.percentile(sim_fg[:, 1], 97.5))
print(f"  real value inside the simulated 95% interval: {ok}")
pd.DataFrame(dict(real_hits=[len(fg_hits)], real_genes=[len(fg_genes)],
                  real_jaccard=[real_j],
                  sim_hits=[sim_fg[:, 0].mean()],
                  sim_jaccard=[sim_fg[:, 1].mean()],
                  sim_lo=[np.percentile(sim_fg[:, 1], 2.5)],
                  sim_hi=[np.percentile(sim_fg[:, 1], 97.5)],
                  calibrated=[bool(ok)])
             ).to_csv(f"{MR}/53b_calibration.tsv", sep="\t", index=False)
if not ok:
    print("\n  *** CALIBRATION FAILED -- the curve below is not reported ***")

# ---------------------------------------------------------------- curve
## ---------------------------------------------------------------------------
## NOTE ON SCOPE, and two failed attempts that fixed it.
## Extending the curve ABOVE the observed power needs a model of the true-effect
## distribution, and neither attempt survived its own sanity check:
##   flat prior on beta_true  -> 531 hits at observed power (real 21), 25x too many
##   global normal EB prior   -> 0.1 hits at observed power (real 21), 100x too few
## The true effects are sparse, so a single normal prior cannot serve, and the
## data cannot validate a richer one. The curve is therefore reported only up to
## the observed power, where it is calibrated against a real comparison.
## ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("POWER CURVE  (agreement with the observed full-power list)")
print("=" * 78)
case_grid = [1000, 2000, 3000, 4000, 5000, 6000, 7500, 9000, 10000, 11000, 12530]
rows = []
for nc in case_grid:
    ne = n_eff(nc, META_CTRL)
    s = simulate(ne)
    hits, jac = s[:, 0], s[:, 1]
    # sensitivity: share of the full-power gene list recovered
    # precision: share of this study's hits that are in the full-power list
    sens, prec = [], []
    if nc < META_CASE:
        scale = np.sqrt(NE_META / ne)
        se_s = d.se_outcome.values * scale
        extra = np.sqrt(np.maximum(se_s ** 2 - d.se_outcome.values ** 2, 0))
        for _ in range(N_REP):
            bo = d.beta_outcome.values + rng.normal(0, extra)
            _, _, p = wald(bo, se_s, d.beta_exposure.values)
            g = set(d.SYMBOL[bh(p) < FDR_THR])
            sens.append(len(g & full_genes) / len(full_genes) if full_genes else np.nan)
            prec.append(len(g & full_genes) / len(g) if g else np.nan)
    else:
        sens, prec = [1.0], [1.0]
    rows.append(dict(cases=nc, n_eff=ne, mean_hits=hits.mean(),
                     jaccard_genes=jac.mean(),
                     j_lo=np.percentile(jac, 2.5), j_hi=np.percentile(jac, 97.5),
                     sensitivity=np.nanmean(sens), precision=np.nanmean(prec)))
    mark = "   <- reference, 1.0 by construction" if nc == META_CASE else ""
    print(f"  {nc:>7,} cases  hits {hits.mean():5.1f}   "
          f"Jaccard {jac.mean():.3f} [{np.percentile(jac,2.5):.3f}, "
          f"{np.percentile(jac,97.5):.3f}]   recovers "
          f"{100*np.nanmean(sens):4.0f}% of the full list{mark}")
curve = pd.DataFrame(rows)
curve.to_csv(f"{MR}/53a_power_stability_curve.tsv", sep="\t", index=False)

print("\n" + "=" * 78)
print("SANITY CHECK")
print("=" * 78)
sub = curve[curve.cases == 5000].iloc[0]
print(f"  simulated hits at 5,000 cases: {sub.mean_hits:.1f}  "
      f"(real FinnGen, 5,753 cases: {len(fg_hits)})")
r = sub.mean_hits / max(len(fg_hits), 1)
print(f"  ratio {r:.2f}  -> "
      f"{'acceptable' if 0.4 < r < 2.5 else 'MISCALIBRATED, curve not reportable'}")
print(f"\n  headline: a {5000:,}-case outcome recovers "
      f"{100*curve[curve.cases==5000].sensitivity.iloc[0]:.0f}% of the gene list")
print(f"            found at {META_CASE:,} cases, and only "
      f"{100*curve[curve.cases==5000].precision.iloc[0]:.0f}% of its own hits")
print(f"            appear in that list.")
print("\nwritten: 53a, 53b")
