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
import os
import re

import pandas as pd

import sys
# Resolve from this file, or from argv[1] / CQTNA_DIR, so the packet runs
# wherever it is extracted rather than only on the author's machine.
MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR")
      or os.path.dirname(os.path.abspath(__file__)))
TAB = chr(9)

def _find(*parts):
    """Locate a file whether we were pointed at the project root, the review
    packet root, or the packet's reproduce/ subdirectory."""
    for base in (MR, os.path.join(MR, ".."), os.path.join(MR, "reproduce")):
        c = os.path.join(base, *parts)
        if os.path.exists(c):
            return c
    return os.path.join(MR, *parts)


_ms = _find("manuscript", "MANUSCRIPT_GB.md")
if not os.path.exists(_ms):
    _ms = _find("MANUSCRIPT_GB.md")
txt = io.open(_ms, encoding="utf-8").read()

grid = pd.read_csv(_find("123d_fixed_anchor_full_grid.tsv"), sep=TAB)
main = grid[grid.analysis == "main"]
off = pd.read_csv(_find("126a_offgrid_attribution.tsv"), sep=TAB)
perm = pd.read_csv(_find("126c_permutation_primary.tsv"), sep=TAB)
mb = pd.read_csv(_find("85e_matched_background_fixed_anchor.tsv"), sep=TAB)
ml = pd.read_csv(_find("130c_multilist_verdict.tsv"), sep=TAB)

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
# The margin is taken over ALL comparators, not only significant ones: the
# significance-restricted version is 0 when nothing reaches P<0.05, which is
# what produced the withdrawn "unbounded margin" claim.
_ra = ml[ml.cell == "RA x Soskic_CD4"].iloc[0]
_raq = ml[ml.cell == "RA x eQTLGen_blood"].iloc[0]
_mel = ml[ml.cell == "melanoma x Soskic_CD4"].iloc[0]
# The narrowest and widest cancer margins are COMPUTED, not named. They used to
# be hardcoded as melanoma x CD4 and HCC_low x CD4, which was right when it was
# written and stopped being right the moment the comparator lists were repaired:
# the widest moved to HCC_low x eQTLGen. An audit that names the cell it expects
# cannot notice that the ranking changed.
_cancer = ml[~ml.cell.str.startswith("RA ")]
_narrow = _cancer.margin_all.min()
_wide = _cancer.margin_all.max()
want += [
    ("S39 RA x CD4 own fold", f"{_ra.F_own:.2f}"),
    ("S39 RA x CD4 largest comparator", f"{_ra.F_max_all:.2f}"),
    ("S39 RA x CD4 margin", f"{_ra.margin_all:.2f}"),
    ("S39 RA x blood margin", f"{_raq.margin_all:.2f}"),
    ("S39 narrowest cancer margin", f"{_narrow:.2f}"),
    ("S39 widest cancer margin", f"{_wide:.2f}"),
]
bad = 0
for label, v in want:
    # The manuscript rounds some folds to one decimal, so either rendering
    # counts -- but the match must not run into a longer number. Plain
    # substring matching passed "4.84" because the text contained 4.87, and the
    # one-decimal fallback made that collision easy to hit. Require that the
    # matched digits are not followed by another digit.
    hit = any(re.search(re.escape(c) + r"(?!\d)", txt)
              for c in (v, f"{float(v):.1f}"))
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
    # Withdrawn claims. These are not numbers, but the same check applies:
    # they must not reappear in the manuscript.
    "unbounded margin": [],
    "most favourable of the four": [],
    "stronger test": [],
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
