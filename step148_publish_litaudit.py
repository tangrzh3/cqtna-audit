#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 148 -- assemble the literature audit for publication.

A third-party review made the point that the audit is the paper's evidence that
this misreading is a field-level problem rather than a local one, and that it
cannot carry that weight while only the summary percentages are visible. So the
corpus, the per-paper scores, the double-coded subsample and the coding rules
all go out as one supplementary set, with every PMID.

What goes out
  S1  the 152 scorable full texts, per paper, with both scoring passes
  S2  the 46-paper two-coder subsample, per judgement, both coders and the
      adjudicated value
  S3  the coding rules verbatim, as written before any full text was read

Nothing here recomputes a score. If a number in the manuscript disagrees with a
number here, the manuscript is wrong.

Outputs: 148a_litaudit_corpus.tsv
         148b_litaudit_doublecoded.tsv
         148c_litaudit_coding_rules.txt
         148d_console.log
"""
import io
import os
import sys

import pandas as pd

MR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CQTNA_DIR") or \
     os.path.dirname(os.path.abspath(__file__))
CRIT = ["C1_known_locus", "C2_coloc", "C3_power_stability",
        "C4_instrument_count", "C5_smr_heidi"]

# ---- S1: the corpus, with both passes side by side -------------------------
auto = pd.read_csv(os.path.join(MR, "91b_fulltext_scores.tsv"), sep="\t")
strict = pd.read_csv(os.path.join(MR, "91d_rescored.tsv"), sep="\t")

corpus = auto[["pmid", "pmc", "year", "journal", "title", "chars"]].copy()
for c in CRIT:
    corpus["automated_" + c] = auto[c]
    corpus["strict_" + c] = strict["strict_" + c]
    corpus["wide_" + c] = strict["wide_" + c]
corpus = corpus.sort_values(["year", "journal", "pmid"])
corpus.to_csv(os.path.join(MR, "148a_litaudit_corpus.tsv"), sep="\t", index=False)
print("S1  corpus: %d papers, %d journals, %d-%d"
      % (len(corpus), corpus.journal.nunique(),
         corpus.year.min(), corpus.year.max()), flush=True)
for c in CRIT:
    print("      %-22s automated %3d  strict %3d  wide %3d"
          % (c, int(auto[c].sum()), int(strict["strict_" + c].sum()),
             int(strict["wide_" + c].sum())), flush=True)

# ---- S2: the double-coded subsample ----------------------------------------
blind = pd.read_csv(os.path.join(MR, "105a_blind_coding_coded.tsv"), sep="\t")
adj = pd.read_csv(os.path.join(MR, "105c_adjudicated.tsv"), sep="\t")
dbl = adj.merge(
    blind[["pmid", "criterion", "pmc", "journal", "year", "evidence",
           "coder2_note"]],
    on=["pmid", "criterion"], how="left")
dbl = dbl[["pmid", "pmc", "journal", "year", "criterion", "evidence",
           "automated", "manual", "coder2", "coder2_note", "final"]]
dbl = dbl.sort_values(["criterion", "pmid"])
dbl.to_csv(os.path.join(MR, "148b_litaudit_doublecoded.tsv"), sep="\t",
           index=False)
print("\nS2  double-coded: %d judgements over %d papers"
      % (len(dbl), dbl.pmid.nunique()), flush=True)
for c, g in dbl.groupby("criterion"):
    dis = int((g.manual != g.coder2).sum())
    print("      %-22s coder1 %d, coder2 %d, final %d, disagreements %d"
          % (c, int(g.manual.sum()), int(g.coder2.sum()), int(g["final"].sum()),
             dis), flush=True)

# ---- S3: the rules, verbatim ------------------------------------------------
src = os.path.join(MR, "105_CODING_RULES.txt")
rules = io.open(src, encoding="utf-8").read()
with io.open(os.path.join(MR, "148c_litaudit_coding_rules.txt"), "w",
             encoding="utf-8", newline="\n") as f:
    f.write("# Coding rules for the literature audit\n"
            "#\n"
            "# Reproduced verbatim from 105_CODING_RULES.txt, which was written\n"
            "# before any full text was read. Not edited for this release.\n"
            "#\n")
    f.write(rules)
print("\nS3  coding rules copied verbatim (%d characters)" % len(rules),
      flush=True)
print("\nwrote 148a / 148b / 148c", flush=True)
