#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 151 -- how much of the audited literature the identity actually reaches.

WHY
The sharpest objection to this paper is that the identity holds only for a
single-instrument Wald ratio, and that competent studies use IVW or weighted
median, so the target is small. That is an empirical claim about the literature
and it has never been measured here. This measures it on the same 152-paper
corpus the audit already uses.

WHAT IT COUNTS
For each cached full text, whether it names a Wald ratio, IVW, weighted median,
MR-Egger or mode-based estimator, and whether it says anywhere that its cis
instrument is a single variant. A paper can name several: the pipelines that
matter here are the ones whose *cis-eQTL* nomination rests on one variant per
gene, and papers frequently run IVW on trans or on a relaxed set while the cis
nomination is a Wald ratio.

WHAT IT CANNOT DO -- REVISED 2026-08-26, AND THIS IS THE POINT OF THE SCRIPT NOW
The first version of this header claimed the single-instrument count was a LOWER
bound on the identity's reach and the multi-instrument count an UPPER bound on
what escapes it. Both claims were wrong and are withdrawn.

The `single` pattern matches ANYWHERE in the full text. It is not required to
occur near a cis-eQTL mention, to describe the primary candidate list, or to
describe a Wald estimate. Twenty-three of the seventy-one papers it hits never
mention a cis-eQTL at all. Requiring the phrase to fall within 500 characters of
a cis mention leaves 13; requiring the same sentence leaves 10. A count that
falls sevenfold under a minimal proximity requirement is a phrase-hit tally of
unknown direction -- it can be inflated by irrelevant "single-cell" adjacent
prose and deflated by papers whose per-gene instrument count appears only in a
supplementary table.

So this script is a CORPUS DISCOVERY TOOL and nothing else. The interpretable
quantity -- which estimator produced each paper's primary cis nomination list --
is coded by hand as criterion C6 (step156).

Outputs: 151a_estimator_usage.tsv, 151b_console.log
"""
import io
import os
import re
import sys

MR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CQTNA_DIR") or \
     os.path.dirname(os.path.abspath(__file__))
DIRS = [os.path.join(MR, "litaudit"), os.path.join(MR, "105a_research_cache",
                                                    "fulltext")]

PAT = {
    "wald": r"\bwald(?:'s)?\s+(?:ratio|estimat|method)",
    "ivw": r"inverse[- ]variance[- ]weight|\bIVW\b",
    "wmedian": r"weighted\s+median",
    "egger": r"MR[- ]Egger",
    "mode": r"weighted\s+mode|mode[- ]based\s+estimat",
    # a cis nomination resting on one variant
    "single": (r"single[- ](?:SNP|variant|instrument)|only\s+one\s+(?:SNP|variant|instrument)"
               r"|one\s+(?:SNP|variant|instrument)\s+(?:per|for\s+each)"
               r"|a\s+single\s+cis[- ]eQTL|top\s+cis[- ]eQTL"),
    "cis": r"\bcis[- ]eQTL|\bcis[- ]QTL",
}

seen, rows = set(), []
for d in DIRS:
    if not os.path.isdir(d):
        continue
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".txt"):
            continue
        pmc = fn[:-4]
        if pmc in seen:
            continue
        seen.add(pmc)
        try:
            t = io.open(os.path.join(d, fn), encoding="utf-8",
                        errors="replace").read()
        except Exception:
            continue
        low = t.lower()
        hit = {k: bool(re.search(p, low, re.I)) for k, p in PAT.items()}
        # How close is the single-variant phrase to a cis mention? A phrase
        # matched anywhere in a full text says nothing about the cis analysis,
        # and the distance between "anywhere" and "same sentence" turned out to
        # be a factor of seven.
        near = sent = False
        for m in re.finditer(PAT["single"], low, re.I):
            w = low[max(0, m.start() - 500):m.end() + 500]
            if re.search(PAT["cis"], w, re.I):
                near = True
            a0 = low.rfind(".", 0, m.start())
            a0 = 0 if a0 < 0 else a0 + 1
            b0 = low.find(".", m.end())
            b0 = len(low) if b0 < 0 else b0
            if re.search(PAT["cis"], low[a0:b0], re.I):
                sent = True
        rows.append(dict(pmc=pmc, chars=len(t),
                         single_near_cis=int(near), single_same_sentence=int(sent),
                         **{k: int(v) for k, v in hit.items()}))

n = len(rows)
print("cached full texts scanned: %d" % n, flush=True)
if not n:
    sys.exit(1)


def c(k):
    return sum(r[k] for r in rows)


cis = [r for r in rows if r["cis"]]
print("  mention cis-eQTL/cis-QTL: %d" % len(cis), flush=True)
print(flush=True)
print("estimator named anywhere in the paper:", flush=True)
for k, lbl in [("wald", "Wald ratio"), ("ivw", "IVW"),
               ("wmedian", "weighted median"), ("egger", "MR-Egger"),
               ("mode", "mode-based")]:
    print("  %-18s %3d / %d  (%.0f%%)" % (lbl, c(k), n, 100 * c(k) / n),
          flush=True)

multi = [r for r in rows if r["ivw"] or r["wmedian"] or r["egger"] or r["mode"]]
onlyw = [r for r in rows if r["wald"] and not (r["ivw"] or r["wmedian"]
                                               or r["egger"] or r["mode"])]
single = [r for r in rows if r["single"]]
neither = [r for r in rows if not r["wald"] and not (r["ivw"] or r["wmedian"]
                                                     or r["egger"] or r["mode"])]
print(flush=True)
print("  names Wald and no multi-instrument estimator : %3d  (%.0f%%)"
      % (len(onlyw), 100 * len(onlyw) / n), flush=True)
print("  names any multi-instrument estimator         : %3d  (%.0f%%)"
      % (len(multi), 100 * len(multi) / n), flush=True)
print("  carries a single-variant phrase anywhere    : %3d  (%.0f%%)"
      % (len(single), 100 * len(single) / n), flush=True)
print("  names no estimator this scan recognises      : %3d  (%.0f%%)"
      % (len(neither), 100 * len(neither) / n), flush=True)
print(flush=True)
print("  of the papers naming a multi-instrument estimator,", flush=True)
print("  %d also carry the single-variant phrase somewhere."
      % sum(1 for r in multi if r["single"]), flush=True)

print(flush=True)
print("HOW FAST THE SINGLE-VARIANT COUNT COLLAPSES", flush=True)
print("  phrase anywhere in the full text             : %3d" % len(single),
      flush=True)
print("  ... of which never mention a cis-eQTL        : %3d"
      % sum(1 for r in single if not r["cis"]), flush=True)
print("  phrase within 500 characters of a cis mention: %3d"
      % sum(r["single_near_cis"] for r in rows), flush=True)
print("  phrase in the same sentence as a cis mention : %3d"
      % sum(r["single_same_sentence"] for r in rows), flush=True)
print("  -> sevenfold. This is not a bound in either direction.", flush=True)

with io.open(os.path.join(MR, "151a_estimator_usage.tsv"), "w",
             encoding="utf-8", newline="\n") as f:
    cols = ["pmc", "chars", "cis", "wald", "ivw", "wmedian", "egger", "mode",
            "single", "single_near_cis", "single_same_sentence"]
    f.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda x: x["pmc"]):
        f.write("\t".join(str(r[c_]) for c_ in cols) + "\n")
print("\nwrote 151a_estimator_usage.tsv", flush=True)
print("", flush=True)
print("NOTE: this is a corpus discovery tool, not an estimate. Phrase matching",
      flush=True)
print("      cannot say which estimator produced a primary cis nomination list,",
      flush=True)
print("      and these counts are description only. The interpretable quantity",
      flush=True)
print("      is coded by hand as criterion C6 -- see step156.", flush=True)