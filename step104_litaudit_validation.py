#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 104 -- turn the literature count into a validated reporting audit.

Raised in review: the 152-paper survey is full-text keyword matching, which
supports "these checks were rarely reported in the accessible full texts we
surveyed" but not "current practice rarely performs these checks".

This does three things and refuses to do a fourth.

  1. PRISMA-style flow, reconstructed from the search record so every number
     between the query and the scored set is accounted for.
  2. A fixed-seed 30% subsample drawn for manual adjudication, with the evidence
     each paper would be judged on extracted for it: for papers the automated
     coder scored positive, every matching passage; for papers it scored
     negative, every sentence that pairs a prior-knowledge word with a locus word,
     which is where a false negative would hide.
  3. Precision and recall of the automated coder against that adjudication, once
     the adjudication file is filled in.

  4. Kappa is NOT computed here. It needs a second independent coder, and when
     one is available the blind file, the agreement statistics and the joint
     adjudication live in step105; this script then reads the adjudicated codes
     rather than one person's. Computing a kappa by having the same person code
     twice, or by treating the regex as a "coder", would misrepresent the design.

Focus is C1 and C3, the two criteria the paper's own claims rest on (7.9% and
0.7%); C2/C4/C5 are contextual and are left as reported counts.

Outputs: 104a_prisma.tsv, 104b_adjudication.tsv (to be filled), 104c_validation.tsv
"""
import sys
import csv
import glob
import os
import random
import re

import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
SCRATCH = os.path.join(MR, "litaudit")
SEED = 104
SUBSAMPLE_FRAC = 0.30

# The patterns MUST be the ones that produced the published scores, not a
# paraphrase of them. Importing step91b guarantees that: an earlier version of
# this script reconstructed them by hand, missed that WIDE C1 contains
# r"novel loc" -- which fires on "novel locus" in almost any GWAS paper -- and
# so extracted evidence that did not correspond to the score being audited.
import importlib.util as _ilu

_spec = _ilu.spec_from_file_location("s91b", os.path.join(MR, "step91b_rescore.py"))
_s91b = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_s91b)
WIDE_C1 = _s91b.WIDE["C1_known_locus"]
WIDE_C3 = _s91b.WIDE["C3_power_stability"]
# where a false negative on C1 would hide: prior-knowledge word near a locus word
FN_PROBE = re.compile(
    r"[^.]*\b(previously|known|established|reported|prior|literature)\b[^.]*"
    r"\b(loci|locus|GWAS|susceptibility|risk variant|risk allele)\b[^.]*\.",
    re.I)


def sentences_matching(text, patterns, width=260):
    hits = []
    for p in patterns:
        for m in re.finditer(p, text, re.I):
            a, b = max(0, m.start() - width // 2), min(len(text), m.end() + width // 2)
            hits.append(" ".join(text[a:b].split()))
    seen, out = set(), []
    for h in hits:
        k = h[:80]
        if k not in seen:
            seen.add(k)
            out.append(h)
    return out


def main():
    # ---------------------------------------------------------- 1. PRISMA flow
    a = pd.read_csv(f"{MR}/91a_search.tsv", sep="\t")
    scored = pd.read_csv(f"{MR}/91d_rescored.tsv", sep="\t")
    have_txt = {os.path.basename(p)[:-4] for p in glob.glob(f"{SCRATCH}/*.txt")}

    elig = a[a.eligible.astype(bool)]
    flow = [
        ("records returned by the four PubMed queries", len(a)),
        ("  excluded: not eQTL-instrumented MR target nomination",
         int(len(a) - len(elig))),
        ("eligible on title/abstract screen", len(elig)),
        ("  excluded: no open-access full text (no PMC identifier)",
         int(elig.pmc.isna().sum())),
        ("eligible with a PMC identifier", int(elig.pmc.notna().sum())),
        ("  excluded: PMC retrieval failed or returned no body text",
         int(elig.pmc.notna().sum()) - len(have_txt)),
        ("full text retrieved", len(have_txt)),
        ("  excluded: retrieved but not scorable", len(have_txt) - len(scored)),
        ("full text scored (the denominator reported)", len(scored)),
    ]
    pd.DataFrame(flow, columns=["stage", "n"]).to_csv(
        f"{MR}/104a_prisma.tsv", sep="\t", index=False)
    print("PRISMA flow")
    for s, n in flow:
        print(f"  {s:<52} {n:>5}")

    # ---------------------------------------------- 2. fixed-seed subsample
    rng = random.Random(SEED)
    pmids = sorted(scored.pmid.astype(str))
    k = round(len(pmids) * SUBSAMPLE_FRAC)
    sub = sorted(rng.sample(pmids, k))
    print(f"\nsubsample for manual adjudication: {k} of {len(pmids)} "
          f"({100*k/len(pmids):.0f}%), seed {SEED}")

    rows = []
    for pmid in sub:
        r = scored[scored.pmid.astype(str) == pmid].iloc[0]
        path = os.path.join(SCRATCH, f"{r.pmc}.txt")
        if not os.path.exists(path):
            continue
        text = open(path, encoding="utf-8", errors="ignore").read()
        for crit, pats, auto in (
                ("C1_known_locus", WIDE_C1, int(r.wide_C1_known_locus)),
                ("C3_power_stability", WIDE_C3, int(r.wide_C3_power_stability))):
            ev = sentences_matching(text, pats)
            if not ev and crit == "C1_known_locus":
                ev = [" ".join(m.group().split())
                      for m in list(FN_PROBE.finditer(text))[:4]]
            rows.append(dict(pmid=pmid, pmc=r.pmc, journal=r.journal, year=r.year,
                             criterion=crit, automated=auto,
                             n_passages=len(ev),
                             evidence=" ||| ".join(ev[:4]) if ev else "(none found)",
                             manual="", note=""))
    adj = pd.DataFrame(rows)
    # never clobber adjudication already entered: re-running must regenerate the
    # evidence without discarding the human codes it was extracted for
    prior_path = f"{MR}/104b_adjudication.tsv"
    if os.path.exists(prior_path):
        prior = pd.read_csv(prior_path, sep="\t", dtype={"pmid": str})
        if "manual" in prior.columns:
            keep = prior[["pmid", "criterion", "manual", "note"]]
            adj["pmid"] = adj.pmid.astype(str)
            adj = adj.drop(columns=["manual", "note"]).merge(
                keep, on=["pmid", "criterion"], how="left")
            n_kept = adj.manual.notna().sum()
            print(f"  carried forward {n_kept} existing adjudications")
    adj.to_csv(prior_path, sep="\t", index=False)
    print(f"wrote 104b with {len(adj)} rows "
          f"({int(adj.automated.sum())} scored positive by the regex)")

    # ------------------------------------------- 3. validation, if filled in
    # Once a second coder has run and the disagreements have been adjudicated
    # jointly (step105), the adjudicated codes supersede coder 1's. Precision is
    # then a two-coder quantity rather than one person's judgement.
    adj_path = f"{MR}/105c_adjudicated.tsv"
    label = "coder 1 only"
    if os.path.exists(adj_path):
        a2 = pd.read_csv(adj_path, sep="\t", dtype={"pmid": str})
        adj = adj.drop(columns=["manual"]).merge(
            a2[["pmid", "criterion", "final"]].rename(columns={"final": "manual"}),
            on=["pmid", "criterion"], how="left")
        label = "two coders, disagreements adjudicated"
    print(f"\nvalidation basis: {label}")
    filled = adj[adj.manual.astype(str).str.strip().isin(["0", "1", "0.0", "1.0"])]
    filled = filled.assign(manual=filled.manual.astype(float).astype(int))
    if len(filled) == 0:
        print("\n104b is not yet adjudicated; precision/recall will be computed "
              "on the next run once the manual column is filled.")
        return
    out = []
    for crit, g in filled.groupby("criterion"):
        auto = g.automated.astype(int).values
        man = g.manual.astype(int).values
        tp = int(((auto == 1) & (man == 1)).sum())
        fp = int(((auto == 1) & (man == 0)).sum())
        fn = int(((auto == 0) & (man == 1)).sum())
        tn = int(((auto == 0) & (man == 0)).sum())
        prec = tp / (tp + fp) if tp + fp else float("nan")
        rec = tp / (tp + fn) if tp + fn else float("nan")
        out.append(dict(criterion=crit, n=len(g), tp=tp, fp=fp, fn=fn, tn=tn,
                        precision=round(prec, 3), recall=round(rec, 3),
                        agreement=round((tp + tn) / len(g), 3)))
        print(f"  {crit}: n={len(g)} TP{tp} FP{fp} FN{fn} TN{tn}  "
              f"precision {prec:.2f}  recall {rec:.2f}")

        # corrected field-level estimate: apply the measured precision to the
        # full-corpus automated count, with a Wilson interval on the precision
        col = {"C1_known_locus": "wide_C1_known_locus",
               "C3_power_stability": "wide_C3_power_stability"}[crit]
        n_auto = int(scored[col].sum())
        n_tot = len(scored)
        z, npos = 1.96, tp + fp
        if npos:
            ph = tp / npos
            den = 1 + z * z / npos
            cen = (ph + z * z / (2 * npos)) / den
            hw = z * ((ph * (1 - ph) / npos + z * z / (4 * npos * npos)) ** .5) / den
            plo, phi = max(0, cen - hw), min(1, cen + hw)
        else:
            plo = phi = float("nan")
        strict_col = col.replace("wide_", "strict_")
        n_strict = int(scored[strict_col].sum())
        print(f"     automated (wide) {n_auto}/{n_tot} = {100*n_auto/n_tot:.1f}%   "
              f"strict {n_strict}/{n_tot} = {100*n_strict/n_tot:.1f}%")
        print(f"     precision-corrected: {n_auto*prec:.1f}/{n_tot} = "
              f"{100*n_auto*prec/n_tot:.1f}%  "
              f"[{100*n_auto*plo/n_tot:.1f}%, {100*n_auto*phi/n_tot:.1f}%]")
        if n_strict / n_tot < n_auto * plo / n_tot:
            print(f"     ** the strict count ({100*n_strict/n_tot:.1f}%) lies BELOW "
                  f"the corrected interval: reporting it understates the rate **")
        out[-1].update(n_auto_wide=n_auto, n_strict=n_strict, n_total=n_tot,
                       corrected_n=round(n_auto * prec, 1),
                       corrected_pct=round(100 * n_auto * prec / n_tot, 1),
                       corrected_lo_pct=round(100 * n_auto * plo / n_tot, 1),
                       corrected_hi_pct=round(100 * n_auto * phi / n_tot, 1))
    pd.DataFrame(out).to_csv(f"{MR}/104c_validation.tsv", sep="\t", index=False)
    print("\nwrote 104c")
    print("  NOTE: recall is estimated from a near-miss probe over the negatives, "
          "not from\n  exhaustive reading, so it is an upper bound on the "
          "automated coder's sensitivity.\n  Cohen's kappa for the two coders is "
          "in step105 and is reported BEFORE the adjudication\n  whose codes are "
          "used here.")


if __name__ == "__main__":
    main()
