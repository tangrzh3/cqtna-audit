#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 30d  Make the cross-cancer locus attribution interpretable.

Step 30c reported that 31% of lung and 21% of colorectal FDR<0.05 hits sit in
"known melanoma locus" categories, versus 100% for melanoma itself. Those
percentages mean nothing without the BACKGROUND rate: what share of all 3,556
instruments falls in those categories anyway? Melanoma's 100% is only a finding
if the background is much lower, and lung's 31% is only unremarkable if it
matches background.

This is the proper test of finding 1: pigmentation/nevus loci should be enriched
among melanoma hits and NOT among other cancers' hits -- which would show the
enrichment is real melanoma biology, not an artefact of the instrument panel.

Output: 36d
"""
import csv
import os
import collections
from math import lgamma, exp, log

MR = r"D:/R_ex/MR"
NOVEL = "潜在新位点"


def fisher_greater(a, b, c, d):
    """one-sided Fisher exact, P(X >= a) for table [[a,b],[c,d]]"""
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    p = 0.0
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    denom = logC(n, c1)
    for x in range(a, hi + 1):
        p += exp(logC(r1, x) + logC(n - r1, c1 - x) - denom)
    return min(p, 1.0)


# ---------------------------------------------------------------- background
bg = collections.Counter()
mel_sig = collections.Counter()
with open(os.path.join(MR, "13_meta_locus_annotation.tsv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        bg[r["category"]] += 1
        try:
            if float(r["FDR"]) < 0.05:
                mel_sig[r["category"]] += 1
        except (ValueError, TypeError):
            pass

bg_tot = sum(bg.values())
bg_known = bg_tot - bg[NOVEL]
print("background: all instruments in the strict set")
for k, v in bg.most_common():
    print(f"  {k:<18} {v:6,}  {100*v/bg_tot:5.1f}%")
print(f"  {'known-locus total':<18} {bg_known:6,}  {100*bg_known/bg_tot:5.1f}%\n")

# ---------------------------------------------------------------- per cancer
rows = []
mel_tot = sum(mel_sig.values())
mel_known = mel_tot - mel_sig[NOVEL]
p = fisher_greater(mel_known, mel_tot - mel_known, bg_known - mel_known,
                   (bg_tot - bg_known) - (mel_tot - mel_known))
rows.append(dict(cancer="Melanoma (meta)", n_sig=mel_tot, n_known=mel_known,
                 pct_known=round(100 * mel_known / mel_tot, 1) if mel_tot else 0,
                 bg_pct_known=round(100 * bg_known / bg_tot, 1),
                 enrichment=round((mel_known / mel_tot) / (bg_known / bg_tot), 2)
                 if mel_tot else 0, fisher_p=p))

cc = os.path.join(MR, "36a_crosscancer_MR_all.tsv")
by = collections.defaultdict(lambda: collections.Counter())
genes = collections.defaultdict(list)
with open(cc, encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if float(r["FDR"]) < 0.05:
            by[r["cancer"]][r["category"]] += 1
            genes[r["cancer"]].append((r["symbol"], r["category"],
                                       float(r["OR"]), float(r["pval"])))

for c, cnt in by.items():
    tot = sum(cnt.values())
    known = tot - cnt[NOVEL]
    p = fisher_greater(known, tot - known, bg_known - known,
                       (bg_tot - bg_known) - (tot - known))
    rows.append(dict(cancer=c, n_sig=tot, n_known=known,
                     pct_known=round(100 * known / tot, 1) if tot else 0,
                     bg_pct_known=round(100 * bg_known / bg_tot, 1),
                     enrichment=round((known / tot) / (bg_known / bg_tot), 2) if tot else 0,
                     fisher_p=p))

print("=" * 88)
print("Enrichment of KNOWN MELANOMA loci among each cancer's FDR<0.05 hits")
print("=" * 88)
print(f"{'outcome':<20}{'n sig':>7}{'n known':>9}{'% known':>9}"
      f"{'background %':>14}{'enrichment':>12}{'Fisher p':>12}")
print("-" * 88)
for r in rows:
    print(f"{r['cancer']:<20}{r['n_sig']:>7}{r['n_known']:>9}{r['pct_known']:>9}"
          f"{r['bg_pct_known']:>14}{r['enrichment']:>12}{r['fisher_p']:>12.3g}")

with open(os.path.join(MR, "36d_locus_attribution_across_cancers.tsv"), "w",
          newline="", encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

print("\n" + "=" * 88)
print("FDR<0.05 genes per cancer")
print("=" * 88)
for c, gs in genes.items():
    seen, uniq = set(), []
    for s, cat, orr, pv in sorted(gs, key=lambda x: x[3]):
        if s in seen:
            continue
        seen.add(s); uniq.append((s, cat, orr, pv))
    print(f"\n{c}  ({len(uniq)} genes)")
    for s, cat, orr, pv in uniq:
        tag = "" if cat == NOVEL else f"  [{cat}]"
        print(f"   {s:<14} OR={orr:.3f}  p={pv:.2g}{tag}")

print("\nwritten: 36d")
