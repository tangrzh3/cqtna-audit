#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 105 -- prepare a blind second-coder file, and score kappa once returned.

A valid Cohen's kappa needs the second coder to be independent of the first.
That means the file handed over must not carry the first coder's decisions, and
must not carry the automated score either, since either one would anchor the
judgement. Both columns are therefore stripped and the rows are shuffled under a
fixed seed so that the automated positives are not clustered.

Run once to produce 105a_blind_coding.tsv. Fill the `coder2` column, then run
again to get kappa in 105b_kappa.tsv.

Outputs: 105a_blind_coding.tsv, 105b_kappa.tsv
"""
import os
import random

import pandas as pd

MR = r"D:/R_ex/MR"
SEED = 105
BLIND = f"{MR}/105a_blind_coding.tsv"

RULES = """\
CODING RULES -- apply to the passages in the `evidence` column only.
Enter 1 or 0 in the `coder2` column. Leave blank if you truly cannot decide.

C1_known_locus = 1 only if the paper compares ITS OWN significant signal against
  previously reported loci for ITS OWN outcome trait. Examples that are 1: "of
  our N significant genes, k lie within known risk loci for <outcome>"; "we
  checked our candidates against the GWAS Catalog entries for <outcome>".
  Examples that are 0: any match falling inside a reference-list entry; a cited
  paper's title containing "novel loci"; a statement about another study's locus
  count; "annotated to known cell types"; naming a locus without comparing the
  paper's own results to a prior list.

C3_power_stability = 1 only if the paper reports how ITS CANDIDATE LIST depends
  on the outcome GWAS used. Examples that are 1: repeating the analysis against a
  second outcome dataset and comparing the resulting lists; down-sampling the
  outcome; explicitly discussing that more or fewer findings arise with a
  higher-powered outcome GWAS. Examples that are 0: generic "statistical power"
  of the instruments, F-statistics, the eQTL sample size, power calculations for
  detecting an effect, or a limitations sentence saying the cohort was small.

The passages are the same ones the first coder saw. "(none found)" means no
passage matched, which is normally 0.
"""


def kappa(a, b):
    a, b = list(a), list(b)
    n = len(a)
    po = sum(x == y for x, y in zip(a, b)) / n
    pe = sum((a.count(v) / n) * (b.count(v) / n) for v in (0, 1))
    return (po - pe) / (1 - pe) if pe < 1 else float("nan"), po, pe


def main():
    adj = pd.read_csv(f"{MR}/104b_adjudication.tsv", sep="\t", dtype={"pmid": str})

    if not os.path.exists(BLIND):
        out = adj[["pmid", "pmc", "journal", "year", "criterion", "evidence"]].copy()
        out = out.sample(frac=1, random_state=SEED).reset_index(drop=True)
        out.insert(0, "row", range(1, len(out) + 1))
        out["coder2"] = ""
        out["coder2_note"] = ""
        with open(BLIND, "w", encoding="utf-8", newline="") as fh:
            out.to_csv(fh, sep="\t", index=False)
        with open(f"{MR}/105_CODING_RULES.txt", "w", encoding="utf-8") as fh:
            fh.write(RULES)
        print(f"wrote {os.path.basename(BLIND)}: {len(out)} rows, shuffled, "
              f"with no automated score and no first-coder decision")
        print("wrote 105_CODING_RULES.txt")
        print("\nFill the coder2 column, then re-run this script.")
        return

    b = pd.read_csv(BLIND, sep="\t", dtype={"pmid": str})
    b = b[b.coder2.astype(str).str.strip().isin(["0", "1"])]
    if len(b) == 0:
        print(f"{os.path.basename(BLIND)} exists but coder2 is empty; nothing to score.")
        return

    m = b.merge(adj[["pmid", "criterion", "manual", "automated"]],
                on=["pmid", "criterion"], how="left")
    m = m[m.manual.notna()]
    rows = []
    print(f"scoring {len(m)} doubly-coded rows\n")
    for crit, g in m.groupby("criterion"):
        c1 = g.manual.astype(int).tolist()
        c2 = g.coder2.astype(int).tolist()
        k, po, pe = kappa(c1, c2)
        disagree = int(sum(x != y for x, y in zip(c1, c2)))
        rows.append(dict(criterion=crit, n=len(g), agreement=round(po, 3),
                         expected=round(pe, 3), kappa=round(k, 3),
                         n_disagreements=disagree))
        print(f"  {crit}: n={len(g)}  observed agreement {po:.3f}  "
              f"expected {pe:.3f}  kappa {k:.3f}  ({disagree} disagreements)")
        for r in g.itertuples():
            if int(r.manual) != int(r.coder2):
                print(f"     disagree pmid {r.pmid}: coder1={int(r.manual)} "
                      f"coder2={int(r.coder2)}")
    k_all, po_all, pe_all = kappa(m.manual.astype(int), m.coder2.astype(int))
    rows.append(dict(criterion="ALL", n=len(m), agreement=round(po_all, 3),
                     expected=round(pe_all, 3), kappa=round(k_all, 3),
                     n_disagreements=int((m.manual.astype(int) !=
                                          m.coder2.astype(int)).sum())))
    print(f"\n  overall kappa {k_all:.3f} on {len(m)} rows")
    pd.DataFrame(rows).to_csv(f"{MR}/105b_kappa.tsv", sep="\t", index=False)
    print("wrote 105b")
    print("\n  Disagreements are not to be resolved by the first coder alone. "
          "Adjudicate them\n  jointly, record which way each went, and report "
          "kappa BEFORE that reconciliation.")


if __name__ == "__main__":
    main()
