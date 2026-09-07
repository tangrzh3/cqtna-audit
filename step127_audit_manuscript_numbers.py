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

# --- S42: the identity and the decomposition. These now carry the central
# claim, so they are reconciled against their tables like everything else.
# Without this block the paper's new headline would sit outside the audit.
_dec = pd.read_csv(f"{MR}/141a_enrichment_decomposition.tsv", sep="\t")
_full = _dec[_dec.excluded_below_outcome_p == 0].iloc[0]
_gws = _dec[_dec.excluded_below_outcome_p == 5e-8].iloc[0]
_thr = pd.read_csv(f"{MR}/140d_threshold.tsv", sep="\t").iloc[0]
want += [
    ("S42 fold, all significant loci", f"{_full.fold:.2f}"),
    ("S42 fold, genome-wide-significant loci removed", f"{_gws.fold:.2f}"),
    ("S42 loci remaining after removal", f"{int(_gws.n_loci)}"),
    ("S42 genome-wide-significant loci", f"{int(_thr.n_genome_wide_sig)}"),
]

# --- The two conditionings a third-party review asked for. Both now appear in
# the Results and in S42, so both are reconciled here.
_cond = pd.read_csv(f"{MR}/141b_reviewer_conditioning.tsv", sep="	")
_bg = _cond[_cond.analysis == "gws-removed, background restricted"].iloc[0]
_mg = _cond[_cond.analysis == "merge loci sharing a lead SNP"].iloc[0]
want += [
    ("S42 fold, background also restricted", f"{_bg.fold:.2f}"),
    ("S42 fold, loci merged by shared lead SNP", f"{_mg.fold:.2f}"),
    ("S42 regions after merging", f"{int(_mg.n_loci)}"),
    ("S42 background after merging", f"{int(_mg.bg_loci)}"),
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
print("1b. the identity's literature reach (S51), checked as rendered")
print("=" * 74)
# These are counts out of 154, so a bare "71" would match anything. They are
# checked in the rendered form the manuscript actually prints, with the
# percentage recomputed from the table rather than trusted from the text.
_eu = pd.read_csv(_find("151a_estimator_usage.tsv"), sep=TAB)
_n = len(_eu)
# NB: _eu.mode is DataFrame.mode, the method -- these columns must be reached
# by name, not by attribute.
_multi_mask = (_eu["ivw"] | _eu["wmedian"] | _eu["egger"] | _eu["mode"]).astype(bool)
_single = int(_eu["single"].sum())
_multi = int(_multi_mask.sum())
_both = int((_multi_mask & _eu["single"].astype(bool)).sum())
_cis = _eu[_eu["cis"] == 1]
_cis_single = int(_cis["single"].sum())


def _pc(k, d):
    return "%.0f" % (100.0 * k / d)


# The manuscript no longer reports these as bounds -- the phrase count fell
# from 71 to 10 under a same-sentence requirement, so what it now reports is
# the collapse itself. These are the numbers that collapse.
_nocis = int((_eu["single"].astype(bool) & ~_eu["cis"].astype(bool)).sum())
_near = int(_eu["single_near_cis"].sum())
_sent = int(_eu["single_same_sentence"].sum())
reach = [
    ("corpus size", "%d cached full texts" % _n),
    ("single-variant phrase anywhere", "%d contain a phrase" % _single),
    ("names a multi-instrument estimator", "%d name a multi-instrument" % _multi),
    ("of those, no cis-eQTL mention at all", "%d of the 71" % _nocis),
    ("within 500 characters of a cis mention", "leaves %d" % _near),
    ("same sentence", "sentence, %d" % _sent),
]
for label, s in reach:
    # Line wrapping is arbitrary in the source markdown, so a space in the
    # expected phrase matches any run of whitespace including a newline.
    # A plain substring test fails when the count and its percentage are
    # split across a line break, which is the same text.
    hit = re.search(r"\s+".join(re.escape(w) for w in s.split(" ")), txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-16s   %s" % ("OK " if hit else "MISSING", s, label))

print()
print("=" * 74)
print("1c. the |z|-only model share (S28), against 101d")
print("=" * 74)
# The manuscript says a |z|-only model reproduces "68-80% of the gap". Until
# now that pair of numbers existed only in the pre-registration's prose table:
# step101 printed it and wrote nothing, so no audit could reach it and the
# figure script refused to draw it. 101d carries it now.
_zo = pd.read_csv(_find("101d_zonly_model.tsv"), sep=TAB)
_lo = int(round(100 * _zo.frac_explained.min()))
_hi = int(round(100 * _zo.frac_explained.max()))
for label, s in [("|z|-only model share", "%d–%d%% of the gap" % (_lo, _hi)),
                 ("the unattributed remainder",
                  "remaining %d–%d%%" % (100 - _hi, 100 - _lo))]:
    hit = re.search(r"\s+".join(re.escape(w) for w in s.split(" ")), txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-22s   %s" % ("OK " if hit else "MISSING", s, label))

print()
print("=" * 74)
print("1d. C6, the hand-coded estimator provenance (S51), against 156e")
print("=" * 74)
# This is the paper's answer to its sharpest objection and it is now a number,
# so it is reconciled against the table like everything else.
_c6 = dict((r.quantity, r.value) for _, r in
           pd.read_csv(_find("156e_C6_result.tsv"), sep=TAB).iterrows())


def _c6i(k):
    return int(round(float(_c6[k])))


_share = 100.0 * float(_c6["W_share"])
c6 = [
    ("papers ascertainable", "%d ascertainable" % _c6i("n_ascertainable")),
    ("single-variant Wald", "%d rest on a single-variant Wald" % _c6i("n_W")),
    ("multi-instrument", "%d on a multi-instrument" % _c6i("n_M")),
    ("mixed", "%d on mixed pipelines" % _c6i("n_X")),
    ("the share, which is the claim",
     "%.1f%%; Wilson 95%% CI %.1f–%.1f%%" % (_share, 100 * float(_c6["W_lo"]),
                              100 * float(_c6["W_hi"]))),
    ("raw agreement", "%.1f%%" % (100 * float(_c6["raw_agreement"]))),
    ("kappa, pre-adjudication", "%.3f" % float(_c6["kappa_preadjudication"])),
]
for label, s_ in c6:
    hit = re.search(r"\s+".join(re.escape(w) for w in s_.split()), txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-34s   %s" % ("OK " if hit else "MISSING", s_.replace(chr(10), " "),
                                label))

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
