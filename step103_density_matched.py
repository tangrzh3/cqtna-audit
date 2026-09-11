#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 103 -- is the known-locus enrichment explained by locus density?

Executes manuscript/PREREG_density_matched.md (S30).

The mismatched-disease control rules out enrichment on loci indiscriminately
dense across diseases; it does not rule out a density that is itself
disease-specific. The existing permutation (step85e) matches on exposure eQTL-p
decile and outcome allele-frequency quintile, so instrument strength and allele
frequency are controlled and density is not. This runs the missing control.

Design fixed in S30 section 3 and not altered here:
  unit      independent locus, 1 Mb single linkage
  strata    quintile(records per locus) x tertile(unique genes per locus)
            secondary analysis: tertile(span) x tertile(genes)
  null      10,000 replicates; per replicate draw one background locus per
            significant locus from the same stratum, without replacement
  relaxation  sparse stratum -> records quintile only -> whole pool; counted
  statistic  number of drawn loci that are known; one-sided empirical p
  fold       observed known / expected known under the null

Output: 103a_density_matched.tsv
"""
import sys
import os
import numpy as np
import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
LOCUS_KB = KNOWN_KB = 1000
N_PERM = 10000
NOVEL = "潜在新位点"
rng = np.random.default_rng(103)


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


def qcut_codes(v, q):
    """rank-based quantile bins that tolerate heavy ties"""
    r = pd.Series(v).rank(method="average", pct=True)
    return np.minimum((r * q).astype(int), q - 1).values


def permute(loci_tab, sig_ids, strata_cols, label):
    """draw a matched background locus for each significant locus, N_PERM times"""
    tab = loci_tab.copy()
    tab["stratum"] = list(zip(*[tab[c] for c in strata_cols]))
    tab["fallback"] = tab[strata_cols[0]]
    obs = int(tab.loc[tab.index.isin(sig_ids), "known"].sum())
    pool_by_stratum = {k: g.index.values for k, g in tab.groupby("stratum")}
    pool_by_fallback = {k: g.index.values for k, g in tab.groupby("fallback")}
    all_ids = tab.index.values
    known = tab["known"].to_dict()

    relax_full = relax_pool = 0
    null = np.empty(N_PERM, dtype=int)
    for it in range(N_PERM):
        taken, cnt = set(), 0
        for sid in sig_ids:
            for level, pool in ((0, pool_by_stratum.get(tab.at[sid, "stratum"])),
                                (1, pool_by_fallback.get(tab.at[sid, "fallback"])),
                                (2, all_ids)):
                cand = [x for x in pool if x not in taken] if pool is not None else []
                if len(cand) >= 1:
                    if it == 0 and level == 1:
                        relax_full += 1
                    if it == 0 and level == 2:
                        relax_pool += 1
                    pick = cand[rng.integers(len(cand))]
                    taken.add(pick)
                    cnt += int(known[pick])
                    break
        null[it] = cnt
    exp = null.mean()
    fold = obs / exp if exp > 0 else np.nan
    p = (1 + int((null >= obs).sum())) / (N_PERM + 1)
    print(f"  {label:34s} obs {obs}/{len(sig_ids)}  exp {exp:5.2f}  "
          f"fold {fold:6.2f}  p = {p:.4f}   relax: {relax_full} to quintile, "
          f"{relax_pool} to full pool")
    return dict(cell=label, n_sig=len(sig_ids), obs_known=obs,
                exp_known=round(float(exp), 3),
                fold=round(float(fold), 3) if exp > 0 else np.nan,
                emp_p=round(p, 5), relax_to_fallback=relax_full,
                relax_to_pool=relax_pool)


def build(df, known_col, fdr_col):
    g = df.groupby("locus").agg(
        n_rec=(fdr_col, "size"),
        n_gene=("gene_key", "nunique"),
        span=("pos", lambda s: int(s.max() - s.min())),
        known=(known_col, "any"),
        sig=(fdr_col, lambda s: bool((s < 0.05).any())))
    g["q_rec"] = qcut_codes(g.n_rec.values, 5)
    g["t_gene"] = qcut_codes(g.n_gene.values, 3)
    g["t_span"] = qcut_codes(g.span.values, 3)
    return g


def main():
    rows = []

    # ------------------------------------------------ melanoma meta (the 4.09x)
    d = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
    d = d[d.pval.notna()].copy()
    d[["chr", "pos"]] = d.SNP.str.split(":", expand=True).iloc[:, :2]
    d["pos"] = d.pos.astype(int)
    d["locus"] = assign_loci(d.chr.tolist(), d.pos.tolist())
    d["gene_key"] = d.gene_id.fillna(d.SYMBOL)
    d["known_flag"] = d.category != NOVEL
    g = build(d, "known_flag", "FDR")
    sig = g.index[g.sig].tolist()
    print(f"melanoma-meta: {len(g)} background loci, {len(sig)} significant, "
          f"{int(g.known.sum())} known")
    rows.append(permute(g, sig, ["q_rec", "t_gene"], "melanoma-meta primary"))
    rows.append(permute(g, sig, ["t_span", "t_gene"], "melanoma-meta secondary"))

    # ------------------------------------------------ HCC, both power levels
    for cell, fn in (("HCC-high", "85a_HCC_high_annotated.tsv"),
                     ("HCC-low", "85a_HCC_low_annotated.tsv")):
        h = pd.read_csv(f"{MR}/{fn}", sep="\t")
        h["gene_key"] = h.gene_id
        g = build(h, "known", "fdr")
        sig = g.index[g.sig].tolist()
        print(f"{cell}: {len(g)} background loci, {len(sig)} significant, "
              f"{int(g.known.sum())} known")
        rows.append(permute(g, sig, ["q_rec", "t_gene"], f"{cell} primary"))

    out = pd.DataFrame(rows)
    out.to_csv(f"{MR}/103a_density_matched.tsv", sep="\t", index=False)

    print("\n" + "=" * 92)
    print("VERDICT (S30 section 5)")
    print("=" * 92)
    prim = out[out.cell.str.contains("primary")]
    for r in prim.itertuples():
        if r.relax_to_pool > 0.5 * r.n_sig:
            v = "D: stratification failed, no information"
        elif r.fold > 1 and r.emp_p < 0.05:
            v = "A: enrichment survives density matching"
        elif r.fold > 1:
            v = "B: direction holds, not significant after density matching"
        else:
            v = "C: enrichment does not survive density matching"
        print(f"  {r.cell:26s} fold {r.fold:6.2f}  p {r.emp_p:.4f}  ->  {v}")
    print("\n  S30 section 4: attenuation does not by itself mean artefact, since "
          "gene-dense\n  regions may host more genuine signal; and that caveat may "
          "not be used to rescue\n  the claim if the enrichment collapses.")
    print("\nwrote 103a")


if __name__ == "__main__":
    main()
