"""Step 127 -- does the manuscript still agree with the tables?

Two questions, both of which a person re-reading 1,200 lines will get wrong:

  1. every attribution number the manuscript quotes must exist in a table that
     was computed on the frozen partition (S36);
  2. no single-linkage number may survive in the text, except where the text
     explicitly flags it as a superseded and different test.

Question 2 is the one that matters. The failure mode this whole revision came
from is a manuscript carrying two statistical units at once, and that failure is
silent -- every individual number looks fine, and only the mixture is wrong.

Run it after any edit to MANUSCRIPT_GB.md or to any of the tables below.

  python step127_audit_manuscript_numbers.py
Exits non-zero if anything is stale or missing.
"""
import io
import re

import pandas as pd

MR = "D:/R_ex/MR"
TAB = chr(9)

txt = io.open(f"{MR}/manuscript/MANUSCRIPT_GB.md", encoding="utf-8").read()

grid = pd.read_csv(f"{MR}/123d_fixed_anchor_full_grid.tsv", sep=TAB)
main = grid[grid.analysis == "main"]
off = pd.read_csv(f"{MR}/126a_offgrid_attribution.tsv", sep=TAB)
perm = pd.read_csv(f"{MR}/126c_permutation_primary.tsv", sep=TAB)
mb = pd.read_csv(f"{MR}/85e_matched_background_fixed_anchor.tsv", sep=TAB)
ml = pd.read_csv(f"{MR}/130c_multilist_verdict.tsv", sep=TAB)

print("=" * 74)
print("1. numbers the manuscript must contain")
print("=" * 74)
want = []
for _, r in main.iterrows():
    if "MHC" in r.cell or r.control == "FAILED":
        continue
    want.append((r.cell, f"{r.fold:.2f}"))
want += [
    ("R13 transfer", f"{off[off.cell == 'C1 Soskic x melanoma R13'].fold.iloc[0]:.2f}"),
    ("Schmiedel gate", f"{off[off.cell == 'Schmiedel_2018 x melanoma'].fold.iloc[0]:.2f}"),
    ("perm melanoma CD4",
     f"{perm[perm.cell == 'melanoma x Soskic_CD4'].fold_vs_null.iloc[0]:.2f}"),
    ("perm RA CD4",
     f"{perm[perm.cell == 'RA x Soskic_CD4'].fold_vs_null.iloc[0]:.2f}"),
    ("matched bg melanoma",
     f"{mb[(mb.dataset == 'melanoma_meta')].fold.iloc[0]:.2f}"),
]

# The multi-list control (S39). The margin ratios are what the manuscript leans
# on, so they are checked as rendered rather than recomputed by eye.
_ra = ml[ml.cell == "RA x Soskic_CD4"].iloc[0]
_raq = ml[ml.cell == "RA x eQTLGen_blood"].iloc[0]
want += [
    ("S39 RA x CD4 own fold", f"{_ra.F_own:.2f}"),
    ("S39 RA x CD4 best competing", f"{_ra.F_max:.2f}"),
    ("S39 RA x CD4 margin", f"{_ra.F_own / _ra.F_max:.2f}"),
    ("S39 RA x blood margin", f"{_raq.F_own / _raq.F_max:.2f}"),
]
bad = 0
for label, v in want:
    # the manuscript rounds some folds to one decimal; accept either rendering
    hit = v in txt or f"{float(v):.1f}" in txt
    if not hit:
        bad += 1
    print(f"  {'OK ' if hit else 'MISSING'}  {v:>8}   {label}")

print()
print("=" * 74)
print("2. legacy-partition numbers that must NOT appear")
print("=" * 74)
# each entry: value, and the substrings whose presence makes an occurrence legitimate
# Each entry: the legacy value, and the substrings whose presence near an
# occurrence make it legitimate.
#
# 3.44 has two unrelated lives. It is the legacy quantile-stratified density
# fold, which the text may cite only while flagging it as a different test; and
# it is also, by coincidence, one of melanoma-on-blood's mismatch folds in the
# S36 window sweep. Matching on value alone cannot separate them, so the
# mismatch-curve context is allowed too. The check still catches a bare
# "3.44-fold" density claim, which is the thing worth catching.
legacy = {
    "4.09": ["quantile-stratified"],
    "4.44": [],
    "3.46": [],
    "3.57": [],
    "8.85": [],
    "5.11": [],
    "8.67": [],
    "4.36": [],
    "3.44": ["quantile-stratified", "mismatch curve"],
    "11.98": [],
    "6.58": [],
    "9.6-fold": [],
    "4.14": [],
}
for v, allowed in legacy.items():
    for m in re.finditer(re.escape(v), txt):
        ctx = txt[max(0, m.start() - 160):m.start() + 160]
        if any(a in ctx for a in allowed):
            print(f"  ok (flagged as superseded)  {v}")
            continue
        bad += 1
        line = txt[:m.start()].count("\n") + 1
        print(f"  STALE line {line}: {v}  ...{ctx[140:200].strip()}...")

print()
print("=" * 74)
print("audit:", "CLEAN" if bad == 0 else f"{bad} problem(s)")
raise SystemExit(0 if bad == 0 else 1)
