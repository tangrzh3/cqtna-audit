"""Step 85e -- G1 with the background actually named in the pre-registration.

WHY THIS EXISTS
PREREG_hcc_generalisation.md §4 defines G1 as "the share of FDR<0.05 signal that
falls on known HCC loci is higher than the background expectation matched on
eQTL-p decile + EAF quintile".  step85_hcc_tests.py did NOT do that: it used the
unmatched background (every instrument locus in the strict set), one-sided
Fisher.  That was not a silent shortcut -- it is exactly what the melanoma
comparator in the manuscript does (step30e_locus_level_and_fig.py -> 36e, the
4.09x / P=0.028 that Part I reports), and using a different background for HCC
would make the two tumours non-comparable.  But the prereg says what it says, so
the literal version is run here as well and BOTH are reported.  Registered as a
deviation in PREREG §10.

Design fixed before looking at the output:
  * unit = independent locus (1 Mb single linkage), same as step85 / 36e
  * each locus is represented by its lead record = smallest exposure eQTL p
  * strata = decile of exposure eQTL p x quintile of outcome EAF
    (melanoma's annotation table carries no EAF column, so melanoma is run on
     the decile-only variant; HCC is run on both so the two can be compared)
  * null = 10,000 replicates; per replicate draw, for each significant locus,
    one background locus from the same stratum without replacement.  Strata with
    too few members relax to decile-only, then to the whole pool; the number of
    relaxations is printed.
  * statistic = number of drawn loci that are known; one-sided empirical
    p = (1 + #{null >= observed}) / (N + 1)

Output: 85e_matched_background_fixed_anchor.tsv  (the frozen main partition)
        85e_matched_background.tsv  with LOCUS_METHOD=single_linkage, which
        reproduces the originally published rows and carries no inference
"""
import os

import numpy as np
import pandas as pd

MR = r"D:/R_ex/MR"
LOCUS_KB = 1000
# The partition this run uses. "fixed_centre" is the frozen main analysis;
# "single_linkage" reproduces the originally published rows and carries no
# inference. Set by the LOCUS_METHOD environment variable so both can be
# produced without editing the file.
LOCUS_METHOD = os.environ.get("LOCUS_METHOD", "fixed_centre")
N_REP = 10000
NOVEL_MEL = "潜在新位点"
rng = np.random.default_rng(85)


def assign_loci(df, chrom="chr", pos="pos", method=None):
    """Partition variants into independent loci.

    "fixed_centre" is the frozen main analysis (PREREG_locus_partition.md): the
    first unassigned variant on a chromosome becomes a centre and claims every
    variant within LOCUS_KB of it, then the next unassigned variant becomes the
    next centre. Non-recursive, so a locus spans at most LOCUS_KB whatever the
    density, and centres are chosen by position rather than by significance --
    which matters because this same partition supplies the denominator.

    "single_linkage" is the rule this script originally used and the one the
    published numbers were computed on. It is kept so those numbers stay
    reproducible; on a dense resource it chains, and it carries no inference.

    Identical to cqtna:::cq_assign_loci, including the ">" boundary.
    """
    method = method or LOCUS_METHOD
    out = {}
    for ch, sub in df.groupby(chrom):
        sub = sub.sort_values(pos)
        lid = 0
        anchor = None
        for _, r in sub.iterrows():
            if method == "single_linkage":
                if anchor is not None and r[pos] - anchor > LOCUS_KB * 1000:
                    lid += 1
                anchor = r[pos]                     # chains: anchor is the previous variant
            else:
                if anchor is None or r[pos] - anchor > LOCUS_KB * 1000:
                    lid += 1
                    anchor = r[pos]                 # bounded: anchor is the centre
            out[(ch, r[pos])] = f"{ch}_{lid}"
    return out


def qbin(v, k):
    """quantile bin index 0..k-1, ties tolerated, NaN -> -1"""
    v = pd.Series(v).astype(float)
    try:
        b = pd.qcut(v.rank(method="first"), k, labels=False)
    except ValueError:
        b = pd.Series(np.zeros(len(v)), index=v.index)
    b = b.fillna(-1).astype(int)
    return b.values


def run(name, lead, use_eaf):
    """lead: one row per locus with columns locus, known, pdec, eafq"""
    sig = lead[lead.is_sig]
    S = len(sig)
    obs = int(sig.known.sum())
    pool = lead[~lead.is_sig]
    if S == 0:
        return dict(dataset=name, matching="decile+quintile" if use_eaf else "decile",
                    n_sig_loci=0, obs_known=0, exp_known=np.nan, fold=np.nan,
                    emp_p=np.nan, relaxed=0)

    def stratum_pool(r, level):
        if level == 0 and use_eaf:
            m = (pool.pdec == r.pdec) & (pool.eafq == r.eafq)
        elif level <= 1:
            m = pool.pdec == r.pdec
        else:
            m = pd.Series(True, index=pool.index)
        return pool.index[m].to_numpy()

    cand, relaxed = [], 0
    for r in sig.itertuples():
        for level in (0, 1, 2):
            idx = stratum_pool(r, level)
            if len(idx) >= 5:
                if level > (0 if use_eaf else 1):
                    relaxed += 1
                break
        cand.append(idx)

    knownv = lead.known.to_numpy()
    counts = np.empty(N_REP, int)
    for i in range(N_REP):
        picked = set()
        c = 0
        for idx in cand:
            free = [j for j in idx if j not in picked]
            j = rng.choice(free if free else idx)
            picked.add(j)
            c += bool(knownv[lead.index.get_loc(j)])
        counts[i] = c
    exp_ = counts.mean()
    p = (1 + int((counts >= obs).sum())) / (N_REP + 1)
    return dict(dataset=name, matching="decile+quintile" if use_eaf else "decile",
                n_sig_loci=S, obs_known=obs, exp_known=round(float(exp_), 3),
                fold=round(obs / exp_, 3) if exp_ > 0 else np.nan,
                emp_p=round(p, 5), relaxed=relaxed)


rows = []

# ---------------------------------------------------------------- HCC
for tag, fn in (("HCC_high", "85a_HCC_high_annotated.tsv"),
                ("HCC_low", "85a_HCC_low_annotated.tsv")):
    d = pd.read_csv(f"{MR}/{fn}", sep="\t")
    d["chr"] = d["chr"].astype(str)
    # 85a carries a `locus` column, but it was written with the single-linkage
    # rule, so it cannot be reused when this script is asked for another
    # partition. Reassign from coordinates every time.
    snp = d.groupby(["chr", "pos"], as_index=False).size()
    loc_of = assign_loci(snp)
    d["locus"] = [loc_of[(c, p)] for c, p in zip(d.chr, d.pos)]
    d = d.sort_values("pval")
    lead = (d.groupby("locus")
              .agg(known=("known", "any"), is_sig=("fdr", lambda s: (s < .05).any()),
                   pval=("pval", "min"), eaf=("eaf_out", "first"))
              .reset_index())
    lead["pdec"] = qbin(lead.pval, 10)
    lead["eafq"] = qbin(lead.eaf.where(lead.eaf <= .5, 1 - lead.eaf), 5)
    print(f"[{tag}] loci {len(lead)}, known {int(lead.known.sum())}, "
          f"significant {int(lead.is_sig.sum())}")
    for use_eaf in (True, False):
        rows.append(run(tag, lead, use_eaf))

# ---------------------------------------------------------------- melanoma
mel = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
mel[["chr", "pos"]] = mel.SNP.str.split(":", expand=True)
mel["pos"] = mel.pos.astype(int)
snp = mel.groupby(["chr", "pos"], as_index=False).size()
loc_of = assign_loci(snp)
mel["locus"] = [loc_of[(c, p)] for c, p in zip(mel.chr, mel.pos)]
mel = mel.sort_values("pval_exposure")
lead = (mel.groupby("locus")
           .agg(known=("category", lambda s: (s != NOVEL_MEL).any()),
                is_sig=("FDR", lambda s: (s < .05).any()),
                pval=("pval_exposure", "min"))
           .reset_index())
lead["pdec"] = qbin(lead.pval, 10)
lead["eafq"] = 0
print(f"[melanoma_meta] loci {len(lead)}, known {int(lead.known.sum())}, "
      f"significant {int(lead.is_sig.sum())}")
rows.append(run("melanoma_meta", lead, False))

# ---------------------------------------------------------------- eQTLGen
# The manuscript quotes a matched-background version of the whole-blood cell
# (92e). No script in the tree produced 92e, so it is computed here instead, by
# the same machinery as the other datasets rather than by a second
# implementation -- 92c carries everything the matching needs.
eq = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
eq["chr"] = eq["chr"].astype(str)
snp = eq.groupby(["chr", "pos"], as_index=False).size()
loc_of = assign_loci(snp)
eq["locus"] = [loc_of[(c, p)] for c, p in zip(eq.chr, eq.pos)]
eq = eq.sort_values("pval_exp")
lead = (eq.groupby("locus")
          .agg(known=("known", "any"), is_sig=("fdr", lambda s: (s < .05).any()),
               pval=("pval_exp", "min"), eaf=("eaf", "first"))
          .reset_index())
lead["pdec"] = qbin(lead.pval, 10)
lead["eafq"] = qbin(lead.eaf.where(lead.eaf <= .5, 1 - lead.eaf), 5)
print(f"[eQTLGen_blood] loci {len(lead)}, known {int(lead.known.sum())}, "
      f"significant {int(lead.is_sig.sum())}")
for use_eaf in (True, False):
    rows.append(run("eQTLGen_blood", lead, use_eaf))

res = pd.DataFrame(rows)
res.insert(0, "partition", LOCUS_METHOD)
suffix = "" if LOCUS_METHOD == "single_linkage" else "_fixed_anchor"
res.to_csv(f"{MR}/85e_matched_background{suffix}.tsv", sep="\t", index=False)
print()
print(res.to_string(index=False))
print(f"\nwrote 85e_matched_background{suffix}.tsv  [{LOCUS_METHOD}]")
