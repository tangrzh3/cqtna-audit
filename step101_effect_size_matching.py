#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 101 -- is the known-vs-novel recovery gap a property of the class, or of |z|?

Executes manuscript/PREREG_effect_size_matching.md (S28), registered before any
matched analysis was run.

The generative model in step55 contains no class label: recovery is a function of
full-power |z| and the global BH threshold alone (S28 section 0). So this script
is not asking whether the gap survives conditioning on |z| -- it is quantifying
how much |z| accounts for, testing the one thing that is genuinely open (whether
the label adds anything among loci of comparable |z|, which is not settled by the
code because BH is global and loci are not independent), and producing the numbers
the pre-registered replacement wording needs.

Everything upstream is held identical to step55, including its column subset, row
filters, seed and replicate count, because PC requires reproducing 85.8% / 22.8%
before anything else is interpreted.

Outputs: 101a_zdist.tsv, 101b_matched.tsv, 101c_stratified.tsv, 101d_zonly_model.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

# Resolvable outside the author's machine, like the other packet scripts.
MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
rng = np.random.default_rng(1)          # identical to step55
N_REP = 400                             # identical to step55
FDR_THR = 0.05
NOVEL = "潜在新位点"
META_CASE, META_CTRL = 12530, 789099
FRAC_MAIN = 0.5                         # the power point the claim is made at
CALIPER = 0.20                          # |log|z| difference|, pre-registered
N_BOOT = 2000
N_PERM = 1000
LOCUS_KB = 1000


def n_eff(nc, nk):
    return 4.0 / (1.0 / nc + 1.0 / nk)


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def wald_p(bo, seo, be):
    return 2 * stats.norm.sf(np.abs((bo / be) / (seo / np.abs(be))))


def assign_loci(chrs, poss):
    order = sorted(range(len(chrs)), key=lambda i: (str(chrs[i]), int(poss[i])))
    out = [None] * len(chrs)
    lc, lp, lid = None, None, 0
    for i in order:
        c, p = str(chrs[i]), int(poss[i])
        if c != lc or p - lp > LOCUS_KB * 1000:
            lid += 1
        out[i] = lid
        lc, lp = c, p
    return out


# ---------------------------------------------------------------- data (as step55)
raw = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
d = raw[["exposure", "SYMBOL", "category", "beta_exposure", "beta_outcome",
         "se_outcome", "SNP"]].dropna()
d = d[(d.beta_exposure != 0) & (d.se_outcome > 0)].reset_index(drop=True)
d["is_novel"] = d.category == NOVEL
d["absz"] = np.abs(d.beta_outcome.values / d.se_outcome.values)
d[["chr", "pos"]] = d.SNP.str.split(":", expand=True).iloc[:, :2]
d["pos"] = d.pos.astype(int)
d["locus"] = assign_loci(d.chr.tolist(), d.pos.tolist())

p_full = wald_p(d.beta_outcome.values, d.se_outcome.values, d.beta_exposure.values)
hit_full = bh(p_full) < FDR_THR
print(f"records {len(d):,} | full-power FDR<0.05 records {hit_full.sum()}")


def unit_tables(unit):
    """gene-level or locus-level view: representative |z| = max over records,
    recovered if any record passes (identical to step55's set semantics)"""
    key = "SYMBOL" if unit == "gene" else "locus"
    full_known = set(d[key][hit_full & ~d.is_novel])
    full_novel = set(d[key][hit_full & d.is_novel])
    # a locus carrying any known-category record counts as known
    full_novel = full_novel - full_known
    z = d.groupby(key).absz.max()
    return key, full_known, full_novel, z


FRACS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


def full_pass(units):
    """One sweep with a single generator consumed exactly as step55 consumes it:
    created once, drawn once per replicate, replicates nested inside the power
    loop. Re-seeding per power point instead gives a different draw sequence and
    the 80% point then misses the published value by 3.5 pp -- which is Monte
    Carlo noise, but PC is a bitwise reproduction check and must not absorb it.
    Per-unit counting is free of randomness, so both unit levels are tallied
    inside this same sweep rather than in separate passes."""
    r = np.random.default_rng(1)
    curve = []
    cnt = {u: {x: 0 for x in (fk | fn)} for u, (_, fk, fn, _) in units.items()}
    for frac in FRACS:
        nc = int(META_CASE * frac)
        scale = np.sqrt(n_eff(META_CASE, META_CTRL) / n_eff(nc, META_CTRL))
        se_s = d.se_outcome.values * scale
        extra = np.sqrt(np.maximum(se_s ** 2 - d.se_outcome.values ** 2, 0))
        rk, rn = [], []
        for _ in range(N_REP):
            bo = d.beta_outcome.values + r.normal(0, extra)
            h = bh(wald_p(bo, se_s, d.beta_exposure.values)) < FDR_THR
            gk = set(d.SYMBOL[h & ~d.is_novel])
            gn = set(d.SYMBOL[h & d.is_novel])
            fk_g, fn_g = units["gene"][1], units["gene"][2]
            rk.append(len(gk & fk_g) / len(fk_g) if fk_g else np.nan)
            rn.append(len(gn & fn_g) / len(fn_g) if fn_g else np.nan)
            if frac == FRAC_MAIN:
                for uname, (key, fk, fn, _) in units.items():
                    hit = set(d[key][h])
                    for x in hit & cnt[uname].keys():
                        cnt[uname][x] += 1
        curve.append(dict(frac=frac, known=np.nanmean(rk), novel=np.nanmean(rn)))
    freq = {u: {x: c / N_REP for x, c in cnt[u].items()} for u in cnt}
    return curve, freq


# ---------------------------------------------------------------- PC
print("\n" + "=" * 88)
print("PC  reproduce 55a (gene level, step55's own seed and replicate count)")
print("=" * 88)
key_g, fk_g, fn_g, z_g = unit_tables("gene")
print(f"  full-power list: known-locus genes {len(fk_g)}, novel-locus genes {len(fn_g)}")
print(f"    known: {sorted(fk_g)}")
print(f"    novel: {sorted(fn_g)}")

UNITS = {u: unit_tables(u) for u in ("gene", "locus")}
curve, FREQ = full_pass(UNITS)

pub = pd.read_csv(f"{MR}/55a_recovery_by_locus_class.tsv", sep="\t")
rows = []
for c in curve:
    ref = pub[np.isclose(pub.frac, c["frac"])].iloc[0]
    ok = (abs(c["known"] - ref.known_recovery) < 1e-9 and
          abs(c["novel"] - ref.novel_recovery) < 1e-9)
    rows.append(dict(match=ok))
    print(f"  {c['frac']:.0%}  known {100*c['known']:5.1f}% "
          f"(pub {100*ref.known_recovery:5.1f}%)   "
          f"novel {100*c['novel']:5.1f}% (pub {100*ref.novel_recovery:5.1f}%)"
          f"   {'exact' if ok else 'MISMATCH'}")
pc_pass = all(r["match"] for r in rows)
print(f"  PC -> {'PASS (bitwise)' if pc_pass else 'FAIL'}")
if not pc_pass:
    raise SystemExit("\n  PC failed: nothing below is interpreted (S28 section 6, D).")

# ---------------------------------------------------------------- per unit
out_z, out_m, out_s, out_2 = [], [], [], []
for unit in ("gene", "locus"):
    key, fk, fn, zmap = UNITS[unit]
    freq = FREQ[unit]
    kr = np.mean([freq[u] for u in fk])
    nr = np.mean([freq[u] for u in fn])
    gap_obs = kr - nr
    print("\n" + "=" * 88)
    print(f"{unit.upper()} LEVEL   at {FRAC_MAIN:.0%} power: known {100*kr:.1f}%, "
          f"novel {100*nr:.1f}%, gap {100*gap_obs:+.1f} pp")
    print("=" * 88)

    kz = np.array([zmap[u] for u in sorted(fk)])
    nz = np.array([zmap[u] for u in sorted(fn)])

    # ---- M0 |z| distributions
    def desc(a):
        return dict(n=len(a), median=np.median(a), q1=np.percentile(a, 25),
                    q3=np.percentile(a, 75), lo=a.min(), hi=a.max())
    dk, dn = desc(kz), desc(nz)
    print(f"  M0  known |z|: n={dk['n']} median {dk['median']:.2f} "
          f"IQR {dk['q1']:.2f}-{dk['q3']:.2f} range {dk['lo']:.2f}-{dk['hi']:.2f}")
    print(f"      novel |z|: n={dn['n']} median {dn['median']:.2f} "
          f"IQR {dn['q1']:.2f}-{dn['q3']:.2f} range {dn['lo']:.2f}-{dn['hi']:.2f}")
    for cls, dd in (("known", dk), ("novel", dn)):
        out_z.append(dict(unit=unit, cls=cls, **dd))

    # ---- common support: the matched design is only meaningful if the two |z|
    # distributions overlap. S28 did not ask for this and should have.
    nlo, nhi = nz.min(), nz.max()
    in_supp = sorted(u for u in fk if nlo <= zmap[u] <= nhi)
    print(f"  ** common support: {len(in_supp)} of {len(fk)} known units lie inside "
          f"the novel |z| range [{nlo:.2f}, {nhi:.2f}] -> {in_supp}")

    # ---- M1 caliper matching on log|z|
    pairs = []
    for u in sorted(fn):
        lz = np.log(zmap[u])
        cand = [(abs(np.log(zmap[v]) - lz), v) for v in sorted(fk)
                if abs(np.log(zmap[v]) - lz) <= CALIPER]
        if cand:
            pairs.append((u, min(cand)[1], min(cand)[0]))
    print(f"  M1  matched pairs within caliper {CALIPER}: {len(pairs)} "
          f"of {len(fn)} novel units")
    if len(pairs) >= 3:
        dn_r = np.array([freq[a] for a, _, _ in pairs])
        dk_r = np.array([freq[b] for _, b, _ in pairs])
        delta = dk_r.mean() - dn_r.mean()
        bs = [np.mean(dk_r[i] - dn_r[i]) for i in
              (rng.integers(0, len(pairs), len(pairs)) for _ in range(N_BOOT))]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        met = abs(delta) < 0.15 and lo <= 0 <= hi
        print(f"      delta_matched {100*delta:+.1f} pp  [{100*lo:+.1f}, {100*hi:+.1f}]"
              f"   criterion |delta|<15pp and CI contains 0 -> "
              f"{'MET' if met else 'NOT MET'}")
        for a, b, dd in pairs:
            print(f"        {a} (|z|={zmap[a]:.2f}, rec {100*freq[a]:.0f}%)  <-  "
                  f"{b} (|z|={zmap[b]:.2f}, rec {100*freq[b]:.0f}%)")
    else:
        delta = lo = hi = np.nan
        met = None
        print(f"      fewer than 3 matched pairs -> interpretation C "
              f"(rewording still applies, S28 section 6 table)")
    out_m.append(dict(unit=unit, gap_observed=gap_obs, n_novel=len(fn),
                      n_matched=len(pairs), delta_matched=delta,
                      ci_lo=lo, ci_hi=hi, criterion_met=met))

    # ---- M2 how much of the gap does |z| alone account for
    units = sorted(fk | fn)
    x = np.log(np.array([zmap[u] for u in units]))
    y = np.array([freq[u] for u in units])
    novel_flag = np.array([u in fn for u in units], float)
    eps = 1.0 / (2 * N_REP)
    ylog = np.log(np.clip(y, eps, 1 - eps) / (1 - np.clip(y, eps, 1 - eps)))
    b_z = np.polyfit(x, ylog, 1)
    pred = 1 / (1 + np.exp(-np.polyval(b_z, x)))
    gap_pred = pred[novel_flag == 0].mean() - pred[novel_flag == 1].mean()
    expl = gap_pred / gap_obs if gap_obs else np.nan
    X = np.column_stack([np.ones_like(x), x, novel_flag])
    coef, *_ = np.linalg.lstsq(X, ylog, rcond=None)
    print(f"  M2  |z|-only model predicts a gap of {100*gap_pred:+.1f} pp against "
          f"{100*gap_obs:+.1f} pp observed -> |z| accounts for {100*expl:.0f}%")
    print(f"      logit coefficients: log|z| {coef[1]:+.2f}, is_novel {coef[2]:+.2f} "
          f"(point estimates only, n={len(units)}; no significance claimed)")
    # M2 is the number the manuscript quotes ("a model using |z| with no
    # class label reproduces 68-80% of the gap"). It was printed and never
    # written, so no figure could draw it and no audit could reconcile it.
    out_2.append(dict(unit=unit, n_units=len(units), gap_observed=gap_obs,
                      gap_predicted_zonly=gap_pred, frac_explained=expl,
                      coef_log_z=coef[1], coef_is_novel=coef[2]))

    # ---- M3 tertile-stratified
    t1, t2 = np.percentile(np.array([zmap[u] for u in units]), [33.3, 66.7])
    print(f"  M3  |z| tertiles split at {t1:.2f} and {t2:.2f}")
    for lab, sel in (("low", lambda v: v <= t1),
                     ("mid", lambda v: t1 < v <= t2),
                     ("high", lambda v: v > t2)):
        ks = [freq[u] for u in fk if sel(zmap[u])]
        ns = [freq[u] for u in fn if sel(zmap[u])]
        g = (np.mean(ks) - np.mean(ns)) if ks and ns else np.nan
        print(f"      {lab:5s} known n={len(ks):2d} "
              f"{100*np.mean(ks):5.1f}%   novel n={len(ns):2d} "
              f"{100*np.mean(ns):5.1f}%" if ks and ns else
              f"      {lab:5s} known n={len(ks):2d}   novel n={len(ns):2d}"
              f"   (one class empty, no within-stratum gap)")
        out_s.append(dict(unit=unit, stratum=lab, n_known=len(ks), n_novel=len(ns),
                          known_rec=np.mean(ks) if ks else np.nan,
                          novel_rec=np.mean(ns) if ns else np.nan,
                          gap=g))

    # ---- NC label permutation
    allf = np.array([freq[u] for u in units])
    nn = int(novel_flag.sum())
    perm = []
    for _ in range(N_PERM):
        idx = rng.permutation(len(units))
        perm.append(allf[idx[nn:]].mean() - allf[idx[:nn]].mean())
    perm = np.array(perm)
    pval = (np.sum(perm >= gap_obs) + 1) / (N_PERM + 1)
    print(f"  NC  label permutation: observed gap {100*gap_obs:+.1f} pp, "
          f"permuted mean {100*perm.mean():+.1f} pp, one-sided p = {pval:.3g}")
    print(f"      (a low p shows the classes differ; it does NOT show the label "
          f"rather than |z| is why -- S28 section 5)")

pd.DataFrame(out_z).to_csv(f"{MR}/101a_zdist.tsv", sep="\t", index=False)
pd.DataFrame(out_m).to_csv(f"{MR}/101b_matched.tsv", sep="\t", index=False)
pd.DataFrame(out_s).to_csv(f"{MR}/101c_stratified.tsv", sep="\t", index=False)
pd.DataFrame(out_2).to_csv(f"{MR}/101d_zonly_model.tsv", sep="\t", index=False)
print("\nwrote 101a / 101b / 101c / 101d")
