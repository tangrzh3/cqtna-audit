"""Build the Genome Biology version from the full manuscript.

GB constraints: structured abstract ~350 words; main text 6,000-8,000 by
convention; no hard display-item cap. Reviewer 3 asked for a 25-35% cut of the
main text with the process material moved to Supplementary, and for the seven
diagnostics to be the single backbone.

This script is DERIVED, like assemble.py: it never edits the source. It
  * swaps in a GB-format structured abstract (written here, in one place);
  * deletes the blocks Reviewer 3 named as belonging in Supplementary, each by an
    explicit anchor so the deletion is auditable;
  * leaves every retained number byte-identical to the full version.

Anything deleted here already exists in a Supplementary document; the mapping is
printed at the end of the run so it can be checked.

Usage: python build_gb.py
"""
import io
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "MANUSCRIPT_v2_dual_thread.md")
OUT = os.path.join(HERE, "MANUSCRIPT_GB.md")

GB_ABSTRACT = """## Abstract

**Background.** Context-specific expression quantitative trait loci (eQTLs)
combined with Mendelian randomization (MR) are widely used to nominate immune
targets in cancer. Because an instrument selected in one cell state at one
timepoint is usually a single variant, the Wald ratio has test statistic
z = β_out/se_out, so the outcome GWAS supplies every MR significance claim. We
audited what that constraint implies, using CD4⁺ T cell cis-eQTLs from eight
activation profiles against a 12,530-case melanoma meta-analysis.

**Results.** Significant signal concentrates on loci already known for the
outcome: 4.09-fold by independent locus. The pattern recurs across five nested
power levels of one GWAS resource, where no novel-locus gene reaches significance
at any case number between 2,705 and 5,753; in a second disease, where three of
four significant loci are known hepatocellular carcinoma loci including *PNPLA3*;
and when the exposure resource is replaced entirely by a whole-blood eQTL dataset
300-fold larger, which raises significant loci from 7 to 30 and leaves the
enrichment at 4.44-fold (P = 3.7×10⁻¹¹). All six disease-by-resource combinations
enrich, against a mismatched-locus negative control that does not. Raising outcome
power increased MR discoveries while lowering colocalisation support and replaced
the candidate list entirely: two lists from identical exposure data share no
genes. Down-sampling shows the loss falls unevenly by locus class — known-locus
genes are recovered about four times as often as novel-locus ones at half of
observed power. Of 28 glycolytic genes, 3 are instrumentable in this exposure
resource, 2 remain analysable against this outcome and 1 yields a nominal
association; the nominated gene is 6.7-fold higher in malignant cells than in
CD4⁺ T cells, so tissue-level validation of it measures tumour. Published
genotype × pseudotime tests show no dynamic genetic effect for that gene despite
strongly dynamic expression.

**Conclusions.** In this analysis, and within the power range we could observe,
the reproducible part of a candidate list produced by this framework was the part
that did not constitute a discovery. Nomination is a property of the outcome
GWAS; what survives it is a statement about where a pathway's regulation can be
measured well enough to yield an instrument at all.

"""

# (anchor_start, anchor_end_exclusive, what it is, where it now lives)
CUTS = [
    ("**No candidate is strong on all axes (Fig 8).**",
     "# PART II",
     "per-candidate comparison across axes (Fig 8 legend carries the summary)",
     "Supplementary S26"),
    ("**We do not claim that these instruments are confined to activated states.**",
     "**Is this gene quantified reliably at all?**",
     "retraction narrative for the activated-state claim",
     "Box 1; Supplementary S13"),
    ("**Is this gene quantified reliably at all?**",
     "## 3.2 The state",
     "TPI1 quantification-reliability argument (pseudogenes, mismapping controls)",
     "Supplementary S25"),
    ("We tested \u2463b directly. Because the exposure data are identical",
     "**The list itself turned over completely.**",
     "multiple-signal SuSiE adjudication in full",
     "Supplementary S15"),
    ("Applying the identical exposures to five other cancers",
     "This result was then tested twice more",
     "cross-cancer specificity comparison",
     "Supplementary S10"),
    ("The discovery result is otherwise heavily audited:",
     "## 3.4 What Part II licenses",
     "discovery-cohort audit trail and the TOX boundary",
     "Supplementary S19"),
    ("Three of the outcomes in that record are *not* scientific negatives",
     "**How much was chosen, and on what criteria.**",
     "stopping-rule instances and not-testable categories",
     "Supplementary S12 \u00a73"),
    ("**How much was chosen, and on what criteria.**",
     "The consequence for Part II is structural.",
     "selection denominators at four levels",
     "Supplementary S12 \u00a75"),
    ("Two earlier candidate failures belong here",
     "## 2.5 Finding",
     "two earlier candidate failures (PADI4 multi-instrument, GDI2 heterogeneity)",
     "Box 1; Supplementary S12"),
    ("The simulation is calibrated empirically",
     "Melanoma recovers 53%",
     "down-sampling calibration and the two failed extrapolations",
     "Supplementary S9; Methods"),
]


def main():
    s = io.open(SRC, encoding="utf-8").read()

    i = s.index("## Abstract")
    j = s.index("**Keywords**")
    s = s[:i] + GB_ABSTRACT + s[j:]

    removed = []
    for a, b, what, where in CUTS:
        if a not in s or b not in s:
            print(f"  [anchor missing] {what}")
            continue
        i, j = s.index(a), s.index(b)
        if j <= i:
            print(f"  [anchors out of order] {what}")
            continue
        removed.append((len(s[i:j].split()), what, where))
        s = s[:i] + s[j:]

    io.open(OUT, "w", encoding="utf-8").write(s)

    L = s.split("\n")
    def seg(a, b):
        return sum(len(x.split()) for x in L[a:b] if x.strip())
    def find(p):
        return next(k for k, l in enumerate(L) if l.startswith(p))
    body = seg(find("## 1. Introduction"), find("## 5. Methods"))
    ab = seg(find("## Abstract") + 1, find("## 1. Introduction"))
    print(f"wrote {os.path.basename(OUT)}")
    print(f"  abstract {ab} words (GB convention ~350)")
    print(f"  main text (Introduction through Discussion) {body:,} words")
    print("  moved to Supplementary:")
    for w, what, where in removed:
        print(f"    -{w:4d} w  {what}  ->  {where}")


if __name__ == "__main__":
    main()
