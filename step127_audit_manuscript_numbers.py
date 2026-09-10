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
import json
import math
import os
import re

import numpy as np
from scipy import stats

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
def _load_module(fname):
    """Import a deposited step script so its own functions can be reused.

    Reusing the published implementation beats reimplementing it. The scripts
    have no __main__ guard, so importing runs their analysis; stdout is
    swallowed and any SystemExit absorbed. Their outputs are deterministic and
    identical to what is committed -- verified for step85e, whose table is
    byte-identical after an import -- so the side effect is a rerun, not a
    change.
    """
    import contextlib
    import importlib.util
    try:
        spec = importlib.util.spec_from_file_location(
            "_dep_" + fname.replace(".", "_"), _find(fname))
        mod = importlib.util.module_from_spec(spec)
        # An audit must not modify result tables. Importing runs the script,
        # which rewrites its outputs -- byte-identical in content, but enough
        # to dirty the working tree on every run, and a permanently noisy
        # `git status` is how a real change goes unnoticed. Snapshot the
        # tracked files before, restore exactly what the import touched after;
        # the same approach step159's run() uses for failed steps.
        import subprocess
        def _tracked():
            try:
                # status, not diff. The imported script writes CRLF on
                # Windows while the repository normalises to LF, so `git diff`
                # correctly reports no CONTENT change and would leave the file
                # looking modified anyway. Cleaning the tree is the point.
                r = subprocess.run(["git", "status", "--porcelain"], cwd=MR,
                                   capture_output=True, text=True, timeout=60)
                if r.returncode != 0:
                    return None
                return set(ln[3:].strip() for ln in r.stdout.splitlines()
                           if ln[:2].strip() == "M")
            except Exception:
                return None
        _before = _tracked()
        _argv = sys.argv
        sys.argv = [fname, MR]
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                spec.loader.exec_module(mod)
        except SystemExit:
            pass
        finally:
            sys.argv = _argv
            _after = _tracked()
            if _before is not None and _after is not None:
                _touched = sorted(_after - _before)
                if _touched:
                    subprocess.run(["git", "checkout", "--"] + _touched,
                                   cwd=MR, capture_output=True, timeout=60)
        return mod
    except Exception:
        return None


def _num_in_text(v, decimals=(2, 3, 4), sci_below=0.01):
    """Is v printed in the manuscript, at any precision the paper uses?

    The text renders the same quantity as 0.10, 0.026 or 1.5 x 10-5 depending
    on the sentence, so one format string reports MISSING on correct numbers.
    Candidates that round to all zeros are dropped: "0.00" matches somewhere in
    almost any manuscript -- kappa 0.00 here -- so keeping it would pass every
    P below 0.005 on nothing. Returns the rendering that matched, so the output
    shows what was actually found rather than what was tried first.
    """
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return None
    if 0 < abs(v) < sci_below:
        # The paper mixes conventions below a percent: 2.9 x 10-4 in one
        # sentence, 0.0028 in another. Offer the scientific mantissa AND the
        # decimal renderings rather than guessing which the author reached for.
        # Decimals FIRST, mantissa last. Offering "5.2" ahead of "0.0052" let a
        # P value match "+5.2 percentage points" in an unrelated sentence -- the
        # fifth time in this file that a check passed on a looser rendering than
        # the text actually prints. Most precise wins, always.
        cands = [c for c in ("%.*f" % (d, v) for d in sorted(decimals, reverse=True))
                 if float(c) != 0]
        cands.append("%.1f" % (v / 10 ** math.floor(math.log10(abs(v)))))
    else:
        # Most precise FIRST. Trying 2 decimals first let 0.0495 pass on a
        # coincidental "0.05" elsewhere in the text while the sentence actually
        # prints 0.0495 -- a pass for the wrong reason, and the reported match
        # then hides which rendering was found.
        cands = [c for c in ("%.*f" % (d, v) for d in sorted(decimals, reverse=True))
                 if float(c) != 0]
    for c in cands:
        if re.search(re.escape(c) + r"(?!\d)", txt):
            return c
    return None



print("=" * 74)
print("1. numbers the manuscript must contain")
print("=" * 74)
want = []
# The grid table in Results covers five cells; the HCC paragraph quotes the two
# HCC CD4 folds at one decimal. HCC_low x eQTLGen_blood is computed and quoted
# nowhere, and requiring it here did active harm: its only appearance in the
# text was the WRONG-TABLE 13.53 in the density-matched sentence (see 1f), so
# this loop reported OK on a value that was there by mistake and MISSING the
# moment it was corrected. An audit satisfied by an error is worse than one
# that is silent. Cells the manuscript does not quote are now listed rather
# than demanded -- deliberate omission is the author call; drift is not.
GRID_UNQUOTED = []
GRID_RENDERINGS = []
for _, r in main.iterrows():
    if "MHC" in r.cell or r.control == "FAILED":
        continue
    _s2 = f"{r.fold:.2f}"
    _s1 = f"{r.fold:.1f}"
    # Both renderings get recorded: the grid table prints 15.49 while the HCC
    # paragraph writes 15.5-fold, and step160 counts those as separate values.
    GRID_RENDERINGS.append((r.cell, _s2, _s1))
    _seen = (re.search(re.escape(_s2) + r"(?!\d)", txt) is not None
             or re.search(re.escape(_s1) + r"(?!\d)", txt) is not None)
    if _seen:
        # Record the rendering the manuscript ACTUALLY uses. The grid table
        # prints 15.49 but the HCC paragraph writes 15.5-fold, and reporting
        # only the two-decimal form left the one-decimal one looking
        # unaudited when it is the same number from the same cell.
        want.append((r.cell, _s2 if re.search(re.escape(_s2)
                     + r"(?!\d)", txt) else _s1))
    else:
        GRID_UNQUOTED.append((r.cell, _s2))
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
if GRID_UNQUOTED:
    print()
    print("  main-grid cells the manuscript does not quote (not a failure,",
          "but they should not vanish unnoticed):")
    for _c, _v in GRID_UNQUOTED:
        print("      %-34s %s" % (_c, _v))

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
print("2u. TPI1's compartment difference (S30), against 68c")
print("=" * 74)
# "In melanoma single-cell data TPI1 is 2.75 log2 units higher in malignant
# cells than in the CD4+ T cells the instrument came from, in 16 of 16
# patients." This is the compartment problem stated at full strength against
# the paper's own instrument, so all three parts are checked: the difference,
# the unanimity, and the P.
_ct = pd.read_csv(_find("68c_tests.tsv"), sep=TAB)
_t1 = _ct[_ct.variable.astype(str).str.contains("TPI1", na=False)
          & ~_ct.variable.astype(str).str.contains("Glyco", na=False)]
if _t1.empty:
    print("  MISSING  no TPI1 row in 68c")
    bad += 1
else:
    r = _t1.iloc[0]
    for _lab, s2 in (("malignant minus CD4, log2", "%.2f" % r["diff"]),
                     ("patients with malignant higher", "%d" % int(r.n_mal_higher)),
                     ("patients tested", "%d" % int(r.n))):
        hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s2, _lab))

print()
print("=" * 74)
print("2v. the Cochran Q calibration (Methods), against the meta output")
print("=" * 74)
# "the exp(-Q/2) form used initially gave a significant-Q rate of 1.59% against
# 5.33% for the correct form". 5.33% against a null expectation of 5% is the
# quality-control evidence that back-deriving Rashkin's standard errors from
# OR and P worked, so the number is load-bearing for the whole meta.
#
# ⚠ It cannot be reproduced exactly from the deposited file, and the reason is
# worth recording. meta_finngen_rashkin.py counts `qp < 0.05` on the UNROUNDED
# value and only then writes `f"{qp:.4g}"`. Ninety-two records round to exactly
# "0.05" on the way out, so recomputing from the file gives 5.32% while the
# script's own count gave 5.33%. The paper is right; the stored column simply
# cannot answer the question that produced it.
_mgz = _find("meta_melanoma_finngen_rashkin.tsv.gz")
if not os.path.exists(_mgz):
    print("  --   meta output not present; rate not recomputed here")
else:
    import gzip as _gz
    _n = _k = _b = 0
    with _gz.open(_mgz, "rt", encoding="utf-8", errors="replace") as _f:
        _hdr = _f.readline().rstrip("\n").split(TAB)
        _ns, _qp = _hdr.index("n_studies"), _hdr.index("Q_pval")
        for _ln in _f:
            _c = _ln.rstrip("\n").split(TAB)
            if len(_c) <= max(_ns, _qp) or _c[_ns] != "2":
                continue
            _n += 1
            try:
                _v = float(_c[_qp])
            except Exception:
                continue
            if _v < 0.05:
                _k += 1
            elif _c[_qp] == "0.05":
                _b += 1
    _rate = 100.0 * (_k + _b) / _n if _n else 0.0
    s3 = "%.2f%%" % _rate
    hit = re.search(re.escape(s3) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   Q significance rate (%d of %d, +%d rounded to the boundary)"
          % ("OK " if hit else "MISSING", s3, _k, _n, _b))

print()
print("=" * 74)
print("2t. the MC1R region counted twice (S37), recomputed from 13")
print("=" * 74)
# "Two of the eight bounded loci are not two published regions. The windows at
# 16:88.86-89.73 Mb and 16:89.87 Mb fall within 1 Mb of an identical set of
# seven Landi lead variants, so the fixed-anchor partition counts the MC1R
# region twice." The paper reporting a flaw in its OWN partition, and the
# coordinates are what make the claim checkable at all.
#
# No table stores locus boundaries for this cell, so they are recomputed with
# step85e's assign_loci -- whose docstring states it is identical to
# cqtna:::cq_assign_loci including the ">" boundary. Reusing the deposited
# implementation rather than reimplementing the partition is the point: three
# of my apparent discrepancies in this audit were my own reconstructions being
# wrong, and a partition is exactly the kind of thing that is easy to get
# subtly wrong.
_s85 = _load_module("step85e_matched_background.py")
if _s85 is None or not hasattr(_s85, "assign_loci"):
    print("  MISSING  cannot load assign_loci from step85e")
    bad += 1
else:
    _ml = pd.read_csv(_find("13_meta_locus_annotation.tsv"), sep=TAB)
    _ml[["chr", "pos"]] = _ml.SNP.str.split(":", expand=True)
    _ml["pos"] = _ml.pos.astype(int)
    _sn = _ml.groupby(["chr", "pos"], as_index=False).size()
    _lo = _s85.assign_loci(_sn)
    _ml["locus"] = [_lo[(c, p)] for c, p in zip(_ml.chr, _ml.pos)]
    _sg = _ml.groupby("locus").agg(is_sig=("FDR", lambda s: (s < .05).any()))
    _sig = set(_sg[_sg.is_sig].index)
    print("  OK   %-8s   bounded loci reaching FDR < 0.05" % len(_sig))
    for _L in sorted(_sig):
        _w = _ml[_ml.locus == _L]
        if _w.chr.iloc[0] != "16":
            continue
        for _v in (_w.pos.min() / 1e6, _w.pos.max() / 1e6):
            s2 = "%.2f" % _v
            hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
            if not hit:
                bad += 1
            print("  %s  %-8s   chr16 window boundary, Mb"
                  % ("OK " if hit else "MISSING", s2))

print()
print("=" * 74)
print("2s. axis versus activation at the chromatin level (S49), recomputed")
print("=" * 74)
# "Because axis and activation log2 fold-changes correlate at r = 0.647 at the
# chromatin level despite near-orthogonality at the RNA level, the analysis was
# repeated restricted to activation-invariant peaks". The correlation is the
# REASON the restricted analysis exists; if it were small the extra work would
# be unmotivated, and if it were larger the axis would be harder to separate
# from activation at all. step49 prints it and lands it nowhere.
_ax = pd.read_csv(_find("47a_axis_differential_peaks.tsv.gz"), sep=TAB,
                  compression="gzip")
_ct = pd.read_csv(_find("47b_activation_differential_peaks.tsv.gz"), sep=TAB,
                  compression="gzip")
_pk = _ax.merge(_ct[["peak", "rest_cpm", "act_cpm", "lfc"]], on="peak",
                suffixes=("", "_act")).rename(columns={"lfc": "act_lfc"})
_pk = _pk[_pk.consistent & (_pk.total >= 100)]
_rr = stats.pearsonr(_pk.mean_lfc, _pk.act_lfc)
s2 = "%.3f" % _rr.statistic
hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
if not hit:
    bad += 1
print("  %s  %-8s   axis vs activation log2FC, %d peaks"
      % ("OK " if hit else "MISSING", s2, len(_pk)))

print()
print("=" * 74)
print("2r. the residual axis's glycolysis enrichment (S31), against 44b")
print("=" * 74)
# "The axis the gene marks is enriched 21.3-fold for glycolysis and is not
# specific to TPI1". The fold is the evidence that the axis is glycolytic at
# all -- the claim the leave-one-out then defends -- and it was unchecked.
_fam = pd.read_csv(_find("44b_family_composition.tsv"), sep=TAB)
_gly = _fam[(_fam.direction == "up") & (_fam.family == "glycolysis")]
if _gly.empty:
    print("  MISSING  no up/glycolysis row in 44b")
    bad += 1
else:
    s2 = "%.1f" % float(_gly.iloc[0].enrichment)
    hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   glycolysis enrichment, residual axis up-genes"
          % ("OK " if hit else "MISSING", s2))

print()
print("=" * 74)
print("2q. the identity, confirmed numerically (Abstract), against 12")
print("=" * 74)
# "Benjamini-Hochberg at 0.05 corresponds to an outcome P of 2.1e-4." This is
# the identity itself made concrete in the Abstract: the FDR threshold on the
# MR side maps to a threshold on the OUTCOME side alone, because |z| =
# |b_out|/se_out and the exposure contributes nothing to significance. Of every
# number in this paper it is the one most directly carrying the central claim,
# and it was unaudited. Recomputed from the outcome columns rather than read
# off a table, since no table stores it.
_ms = pd.read_csv(_find("12_MR_meta_strict.tsv"), sep=TAB)
_fc = next((c for c in _ms.columns if c.upper() in ("FDR", "FDR_BH")), None)
_bc = next((c for c in _ms.columns
            if "beta_out" in c.lower() or c in ("beta.outcome", "beta_outcome")), None)
_sc = next((c for c in _ms.columns
            if "se_out" in c.lower() or c in ("se.outcome", "se_outcome")), None)
if None in (_fc, _bc, _sc):
    print("  MISSING  12_MR_meta_strict lacks the outcome columns this needs")
    bad += 1
else:
    _sig = _ms[_ms[_fc] < 0.05]
    _z = (_sig[_bc] / _sig[_sc]).abs()
    _pout = 2 * stats.norm.sf(_z)
    s2 = _num_in_text(float(_pout.max()))
    ok = s2 is not None
    if not ok:
        bad += 1
    print("  %s  %-8s   outcome P at the BH boundary (%d records at FDR < 0.05)"
          % ("OK " if ok else "MISSING", s2 or "?", len(_sig)))

print()
print("=" * 74)
print("2p. the single novel HCC nomination (S44), against 85a")
print("=" * 74)
# "The single novel nomination, SUPV3L1, appears only at high power and has
# FDR = 0.98 at low power." Both halves are checked, because the claim is a
# contrast: 0.98 alone would say nothing without the high-power value being
# significant, and if the two were swapped the sentence would invert.
_SUPV = "ENSG00000156502"
for _f, _lab, _sig in (("85a_HCC_high_annotated.tsv", "high power", True),
                       ("85a_HCC_low_annotated.tsv", "low power", False)):
    _t = pd.read_csv(_find(_f), sep=TAB)
    _r = _t[_t.gene_id == _SUPV]
    if _r.empty:
        print("  MISSING  SUPV3L1 absent from %s" % _f)
        bad += 1
        continue
    _fdr = float(_r.iloc[0].fdr)
    if _sig:
        ok = _fdr < 0.05
        print("  %s  %-8s   SUPV3L1 FDR, %s (significant: %s)"
              % ("OK " if ok else "MISSING", "%.4f" % _fdr, _lab, ok))
        if not ok:
            bad += 1
    else:
        s2 = "%.2f" % _fdr
        ok = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
        if not ok:
            bad += 1
        print("  %s  %-8s   SUPV3L1 FDR, %s" % ("OK " if ok else "MISSING", s2, _lab))

print()
print("=" * 74)
print("2o. the chained locus and the discarded spatial test, against 123d and 23g")
print("=" * 74)
# "on the whole-blood resource it chains a chromosome arm into a single 30.8 Mb
# locus" -- the concrete demonstration that single-linkage is unusable on a
# dense resource, which is why the fixed-anchor partition exists at all. From
# the single_linkage rows of 123d, not the fixed_centre ones the rest of the
# audit reads.
_sl = grid[grid.partition == "single_linkage"]
if not _sl.empty:
    s2 = "%.1f" % (_sl.max_sig_span_kb.max() / 1000.0)
    hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   widest single-linkage locus, Mb" % ("OK " if hit else "MISSING", s2))

# "a test whose positive control fails is discarded rather than interpreted --
# the within-lymphoid-subset spatial test (HLA-C P = 0.23) ... discarded under
# this rule". A number quoted to justify DISCARDING a result: if it were
# actually significant the discard would look like suppression, so it is worth
# checking that the test really is null.
_sp = pd.read_csv(_find("23g_ST_within_lymphoid.tsv"), sep=TAB)
_hl = _sp[(_sp.gene == "HLA-C") & _sp.subset.astype(str).str.contains("lymphoid")]
if not _hl.empty:
    s3 = "%.2f" % float(_hl.iloc[0].p)
    hit = re.search(re.escape(s3) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   HLA-C within-lymphoid P" % ("OK " if hit else "MISSING", s3))

print()
print("=" * 74)
print("2m. the corrected literature rates and their intervals (S23), against 104c")
print("=" * 74)
# "not one paper ... compared its significant signal against previously
# reported loci for its own outcome trait (0% [0-3.4]), and an estimated 7.1%
# [2.5-16.1] report how the candidate list depends on the outcome GWAS used".
# These are the paper's headline claims about the literature, they are
# CORRECTED figures with the automated matcher's precision folded in, and not
# one of the six numbers was audited. The precisions themselves -- 0.00 and
# 0.20 -- are what justify calling the raw rates wrong, so they are checked too.
_vl = pd.read_csv(_find("104c_validation.tsv"), sep=TAB)
for _, r in _vl.iterrows():
    for _lab, _v, _fmt in (("corrected rate", r.corrected_pct, "%.1f%%"),
                           ("interval low", r.corrected_lo_pct, "%.1f"),
                           ("interval high", r.corrected_hi_pct, "%.1f"),
                           ("matcher precision", r.precision, "%.2f")):
        if pd.isna(_v):
            continue
        s2 = _fmt % _v
        # A rate of exactly zero is written "0%", not "0.0%". Accept both
        # renderings of the same value rather than reporting a correct number
        # as missing because of a trailing digit.
        _alts = [s2]
        if _fmt.endswith("%%") and float(_v) == 0:
            _alts.append("0%")
        _m = next((a for a in _alts
                   if re.search(re.escape(a) + r"(?!\d)", txt)), None)
        if _m is None:
            bad += 1
        print("  %s  %-8s   %s, %s" % ("OK " if _m else "MISSING", _m or s2,
                                       _lab, r.criterion))

print()
print("=" * 74)
print("2n. lung's leave-one-out (S41), against 139c")
print("=" * 74)
# "lung's result does not survive dropping its single strongest locus
# (P = 0.0138 to 0.0569 without the chr11 FADS1/TMEM258 cluster)". The paper
# reporting that one of its own five outcomes is fragile, with the named
# cluster; if the dropped locus or the resulting P were wrong the concession
# would be misdescribed, and nothing checked either.
_lo1 = pd.read_csv(_find("139c_leave_one_out.tsv"), sep=TAB)
_lg = _lo1[_lo1.row.astype(str).str.contains("Lung", case=False, na=False)]
if not _lg.empty:
    r = _lg.iloc[0]
    for _lab, s2 in (("P, full", _num_in_text(r.p_full) or "?"),
                     ("P, worst locus dropped", _num_in_text(r.p_worst_drop) or "?")):
        hit = s2 != "?" and re.search(re.escape(s2) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s2, _lab))
    for _g in str(r.worst_genes).split(","):
        _g = _g.strip()
        if _g and _g not in txt:
            bad += 1
            print("  MISSING  %-8s   named as the dropped cluster" % _g)
        elif _g:
            print("  OK   %-8s   named as the dropped cluster" % _g)

print()
print("=" * 74)
print("2l. the R13 transfer comparator (S22), against 126a")
print("=" * 74)
# "at 6,226 cases the significant list is the same six genes at the same two
# loci ... locus attribution is 10.0-fold (P = 9.1e-5), all four significant
# loci known". Section 1 checks the fold as 10.02; the text prints 10.0, and
# the P and the counts were unchecked. This is the pre-registered transfer
# test, so the numbers showing the list SURVIVES it are the ones a reader
# would want pinned.
_og = pd.read_csv(_find("126a_offgrid_attribution.tsv"), sep=TAB)
_r13 = _og[_og.note.astype(str).str.contains("R13 transfer", na=False)]
if _r13.empty:
    print("  MISSING  the R13 transfer comparator row is gone from 126a")
    bad += 1
else:
    r = _r13.iloc[0]
    for label, s2 in (("fold, as printed", "%.1f" % r.fold),
                      ("P", _num_in_text(r.fisher_p) or "?"),
                      ("significant loci", "%d" % int(r.sig_loci)),
                      ("of them known", "%d" % int(r.sig_known))):
        ok = s2 != "?" and re.search(re.escape(s2) + r"(?!\d)", txt) is not None
        if not ok:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if ok else "MISSING", s2, label))

print()
print("=" * 74)
print("2k. the transport margins (S41), recomputed from 138a")
print("=" * 74)
# "only melanoma exceeds every mismatched list by a margin (2.30-fold against
# the best rival, against 1.29 for lung and below 1 for the rest)". The margin
# is own fold divided by the strongest comparator's, stored nowhere, and it is
# the number the word "only" rests on -- if a rival's fold rose, the claim
# would fail before the margin looked wrong.
_tm = pd.read_csv(_find("138a_transport_main.tsv"), sep=TAB)
for _row, _own in (("Melanoma", "melanoma"), ("Lung", "lung")):
    _sub = _tm[_tm.row == _row]
    if _sub.empty:
        continue
    _o = float(_sub[_sub.list_name == _own].fold.iloc[0])
    _rivals = _sub[_sub.list_name != _own].fold.astype(float)
    if _rivals.empty or _rivals.max() <= 0:
        continue
    _mg = "%.2f" % (_o / _rivals.max())
    hit = re.search(re.escape(_mg) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s margin over best rival (%.2f / %.2f)"
          % ("OK " if hit else "MISSING", _mg, _row, _o, _rivals.max()))

print()
print("=" * 74)
print("2j. the glycolysis kill criterion (S49), recomputed from 43a")
print("=" * 74)
# "if R2(glycolysis ~ activation + depth) > 0.70 the axis is not separable from
# activation and the analysis stops. Observed R2 = 0.024." A pre-registered
# kill criterion and the value that cleared it -- if the observed R2 were
# wrong, an analysis that should have stopped would have continued.
#
# step43 PRINTS this and never lands it, the same pattern as step101's
# 68-80% before 101d existed. Recomputed here from 43a, which is gzipped;
# nothing in the audit chain read compressed tables until 2026-09-10.
_cs = pd.read_csv(_find("43a_GSE282266_cell_scores.tsv.gz"), sep=TAB,
                  compression="gzip")
_X = np.column_stack([np.ones(len(_cs)), _cs.Activation, np.log1p(_cs.nFeature)])
_beta, *_rest = np.linalg.lstsq(_X, _cs.Glycolysis.values, rcond=None)
_res = _cs.Glycolysis.values - _X @ _beta
_r2 = 1 - _res.var() / _cs.Glycolysis.var()
for label, s2 in (("observed R2", "%.3f" % _r2),
                  ("kill threshold", "0.70")):
    hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s2, label))
if _r2 > 0.70:
    print("  ⚠ R2 exceeds the pre-registered threshold; the analysis should")
    print("    have stopped. This is a kill criterion, not a diagnostic.")
    bad += 1

print()
print("=" * 74)
print("2h. HEIDI on the signals coloc calls distinct (S46), against 16 and 15")
print("=" * 74)
# "HEIDI failed to reject homogeneity for 253 (86.9%)". I could not reproduce
# this at first and nearly filed it as an open question for the author. Two
# mistakes of mine: the wrong tables (14 and 08 rather than 16 and 15), and no
# PP.H4 < 0.2 filter -- the denominator is not every HEIDI test, it is the
# records coloc calls DISTINCT causal variants, which is the whole point of the
# sentence. The definition is in figures/make_gb_fig4_coloc_heidi.py, which
# draws the same panel; it is encoded here so the number stops depending on a
# figure script nobody re-reads.
_hc = pd.read_csv(_find("16_coloc_meta_results.tsv"), sep=TAB)
_hs = pd.read_csv(_find("15_SMR_meta_results.tsv"), sep=TAB)
_hs["_k"] = _hs["gene"].astype(str) + "|" + _hs["profile"].astype(str)
_hm = _hc.merge(_hs[["_k", "p_HEIDI"]], left_on="exposure", right_on="_k")
_hm = _hm[_hm["PP.H4"].notna() & _hm["p_HEIDI"].notna() & (_hm["p_HEIDI"] > 0)]
_lo = _hm[_hm["PP.H4"] < 0.2]
_ps = _lo[_lo["p_HEIDI"] > 0.05]
for label, s2 in (("records coloc calls distinct", "%d" % len(_lo)),
                  ("of them, HEIDI does not reject", "%d" % len(_ps)),
                  ("share", "%.1f%%" % (100.0 * len(_ps) / len(_lo)))):
    hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s2, label))

print()
print("=" * 74)
print("2i. winner's curse, tested and rejected (S44), against 07 and 16")
print("=" * 74)
# "newly entering candidates had LOWER PP.H3+H4 (0.125 versus 0.202)" -- an
# explanation the paper tested and rejected. Same story as 2h: I reconstructed
# it against 14 and got 0.126 / 0.207, decided my definition was wrong, and was
# right about that but wrong about why. It is 16, and the groups are the 244
# exposures new to the meta round against the 127 carried over
# (FINDINGS_step5_pigmentation.md, the winner's-curse entry).
_w2 = pd.read_csv(_find("16_coloc_meta_results.tsv"), sep=TAB)
_prev = set(pd.read_csv(_find("07_coloc_results.tsv"),
                        sep=TAB)["exposure"].astype(str))
_isnew = ~_w2["exposure"].astype(str).isin(_prev)
for label, _sub in (("newly entering", _w2[_isnew]),
                    ("carried over", _w2[~_isnew])):
    _md = (_sub["PP.H3"] + _sub["PP.H4"]).median()
    for s2 in ("%d" % len(_sub), "%.3f" % _md):
        hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
        if s2.startswith("0.") and not hit:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if hit else "-- ", s2, label))

print()
print("=" * 74)
print("2f. the instrument-strength floor (Methods), against 01, 04 and 12")
print("=" * 74)
# "The minimum F statistic was 36.1 in either strict set and 22.2 across the
# harmonised set before instrument selection, so the F > 10 filter was never
# binding at any stage." The claim is that a filter never bound, which is only
# as good as the two minima it rests on -- and a minimum is the one summary a
# single new row can change without touching anything else.
for _f, _lab in (("04_MR_results_strict_all.tsv", "strict set"),
                 ("12_MR_meta_strict.tsv", "strict meta set"),
                 ("01_harmonised_all.tsv", "harmonised set")):
    _d = pd.read_csv(_find(_f), sep=TAB)
    _fc = next((c for c in ("F_stat", "F") if c in _d.columns), None)
    if _fc is None:
        continue
    s2 = "%.1f" % _d[_fc].min()
    hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   minimum F, %s" % ("OK " if hit else "MISSING", s2, _lab))

print()
print("=" * 74)
print("2g. C6's mixed-pipeline ceiling (S51), against 156e")
print("=" * 74)
# "at most 37.8% if every mixed pipeline is counted as partly single-variant".
# Section 1d checks the 22.2% headline and its interval; this is the upper
# bound the paper offers AGAINST itself, derived as (W + X) / ascertainable and
# stored nowhere, so no cell comparison could reach it.
_c6d = dict((r.quantity, float(r.value)) for _, r in
            pd.read_csv(_find("156e_C6_result.tsv"), sep=TAB).iterrows())
_ceil = 100.0 * (_c6d["n_W"] + _c6d["n_X"]) / _c6d["n_ascertainable"]
s3 = "%.1f%%" % _ceil
hit = re.search(re.escape(s3) + r"(?!\d)", txt) is not None
if not hit:
    bad += 1
print("  %s  %-8s   (W %d + X %d) / %d ascertainable"
      % ("OK " if hit else "MISSING", s3, int(_c6d["n_W"]), int(_c6d["n_X"]),
         int(_c6d["n_ascertainable"])))

print()
print("=" * 74)
print("2c. the two conditionings' P values (S42), against 141b")
print("=" * 74)
# Section 1 checks 141b's FOLDS -- 2.07 and 5.28 -- and never their P values.
# The Abstract leans on "1.98-fold, P = 0.41" as its statement that the
# locus-level signal is undetectable, so the P is the load-bearing half.
_cd = pd.read_csv(_find("141b_reviewer_conditioning.tsv"), sep=TAB)
for _, r in _cd.iterrows():
    m = _num_in_text(r.fisher_p)
    if m is None:
        bad += 1
    print("  %s  %-8s   P for %s" % ("OK " if m else "MISSING", m or "?", r.analysis))

print()
print("=" * 74)
print("2d. the matched-background permutation (S38), against 85e")
print("=" * 74)
# 6.55 and 5.66 for eQTLGen, 6.66 with its empirical P of 0.0016 for melanoma.
# 6.55 is the value migrated on 2026-09-09 after the lead-record tie-break was
# made deterministic (PRESPEC_lead_record_tiebreak.md), so it is exactly the
# number that must not drift back unnoticed.
_mb = pd.read_csv(_find("85e_matched_background_fixed_anchor.tsv"), sep=TAB)
for _, r in _mb.iterrows():
    if pd.isna(r.fold) or r.fold <= 0:
        continue
    _lab = "%s %s" % (r.dataset, r.matching)
    s2 = "%.2f" % r.fold
    _fq = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
    if not _fq:
        # Row not quoted at all, so its P is not quoted either. Checking the P
        # anyway found HCC_high's 0.12549 "present" because the text contains
        # 0.125 in a completely unrelated sentence about winner's curse. A P is
        # only checked when its own fold is quoted, which ties it to its row.
        print("  --   %-8s   %s (not quoted)" % (s2, _lab))
        continue
    print("  OK   %-8s   %s" % (s2, _lab))
    if r.emp_p > 1e-4:
        m = _num_in_text(r.emp_p)
        if m is None:
            bad += 1
        print("  %s  %-8s   empirical P, %s"
              % ("OK " if m else "MISSING", m or "?", _lab))

print()
print("=" * 74)
print("2e. the CD4 background share (S42), against 123d")
print("=" * 74)
# "enrichment over the 10.1% background" -- the denominator the headline
# 4.96-fold is a fold OVER. The fold was audited from the first version of this
# script; the number it is relative to was not.
_bm = main[main.cell == "melanoma x Soskic_CD4"].iloc[0]
s3 = "%.1f%%" % (100.0 * _bm.bg_known / _bm.bg_loci)
hit = re.search(re.escape(s3) + r"(?!\d)", txt) is not None
if not hit:
    bad += 1
print("  %s  %-8s   CD4 background (%d/%d)"
      % ("OK " if hit else "MISSING", s3, int(_bm.bg_known), int(_bm.bg_loci)))

print()
print("=" * 74)
print("2a. the automated matcher's own rates (S23), against 148a")
print("=" * 74)
# "Automated matching returned 7.9% and 35.5%" -- the figures the paper reports
# in order to say they were WRONG, against corrected values of 0% and 7.1%.
# Numbers quoted to be disowned still have to be the numbers that were
# produced, and these had never been checked. Also "about 59% performed",
# which is the automated C2 rate.
_cp = pd.read_csv(_find("148a_litaudit_corpus.tsv"), sep=TAB)
_n = len(_cp)
# The strict rate is quoted too -- "0.7% for the second criterion was too
# strict and falls below the corrected interval" -- and it is the number that
# makes the case that BOTH matchers were wrong, in opposite directions.
for _col, _lab in (("wide_C1_known_locus", "wide matcher, locus attribution"),
                   ("wide_C3_power_stability", "wide matcher, power stability"),
                   ("automated_C2_coloc", "automated, colocalisation"),
                   ("strict_C3_power_stability", "strict matcher, power stability")):
    if _col not in _cp.columns:
        continue
    _k = int((pd.to_numeric(_cp[_col], errors="coerce").fillna(0) > 0).sum())
    _pcv = 100.0 * _k / _n
    s = next((c for c in ("%.1f%%" % _pcv, "%.0f%%" % _pcv)
              if re.search(re.escape(c) + r"(?!\d)", txt)), None)
    if s is None:
        bad += 1
    print("  %s  %-8s   %s (%d/%d)"
          % ("OK " if s else "MISSING", s or "%.1f%%" % _pcv, _lab, _k, _n))

print()
print("=" * 74)
print("2b. the HCC power ratio (S44), against 135a")
print("=" * 74)
# "The higher-powered HCC study carries 30.3% of melanoma's effective sample
# size, so we down-sampled melanoma to match". The ratio decides the whole
# down-sampling comparison that follows, and it is a quotient of two numbers
# in one table, so nothing would have recomputed it.
_gd = pd.read_csv(_find("135a_hcc_gate_distance.tsv"), sep=TAB)
_mel = _gd[_gd.cell.str.startswith("melanoma")].n_eff_disc.iloc[0]
_hcc = _gd[_gd.cell.str.startswith("HCC_high")].n_eff_disc.iloc[0]
_r = 100.0 * _hcc / _mel
s2 = "%.1f%%" % _r
hit = re.search(re.escape(s2) + r"(?!\d)", txt) is not None
if not hit:
    bad += 1
print("  %s  %-8s   HCC effective N as share of melanoma (%d/%d)"
      % ("OK " if hit else "MISSING", s2, _hcc, _mel))

print()
print("=" * 74)
print("1y. the list that did not survive the meta round, against 04 and 12")
print("=" * 74)
# "IMPA1's MR P moved from 9.3e-4 to 0.11" -- the single example given for why
# the two lists share no genes, which is one of the paper's sharpest findings
# about its own pipeline. One gene, two tables, and neither number was checked.
_IMPA1 = "ENSG00000133731"
for _f, _lab in (("04_MR_results_strict_all.tsv", "single round"),
                 ("12_MR_meta_strict.tsv", "meta round")):
    _d = pd.read_csv(_find(_f), sep=TAB)
    _col = next((c for c in ("SYMBOL", "gene_id", "exposure") if c in _d.columns), None)
    _m = _d[_d[_col].astype(str).str.contains(_IMPA1 if _col != "SYMBOL" else "IMPA1",
                                              na=False)]
    if _m.empty:
        continue
    _v = float(_m.iloc[0]["pval"])
    s = _num_in_text(_v)
    ok = s is not None
    if not ok:
        bad += 1
    print("  %s  %-8s   IMPA1 P, %s" % ("OK " if ok else "MISSING", s or "?", _lab))

print()
print("=" * 74)
print("1x. the MC1R-region colocalisation and HEIDI (S46), against 07 and 08")
print("=" * 74)
# "PP.H3 dominant and PP.H4 at or near zero (CHMP1A 0.99; VPS9D1-AS1 0.96;
# SPATA33 0.92-1.00)" and "VPS9D1-AS1 (P_HEIDI = 0.649) and CDK10 (0.086,
# 0.081)". These carry the argument that the strongest MR signals in the paper
# point at DIFFERENT causal variants, which is a conclusion against the
# study's own headline and therefore worth pinning.
#
# The HEIDI values come from 08, the single-round table, NOT 15, the meta one --
# 15 gives VPS9D1-AS1 0.634 where the text says 0.649. Checking the wrong table
# first made the text look wrong; recording which table is right stops the next
# reader repeating that.
_cl = pd.read_csv(_find("07_coloc_results.tsv"), sep=TAB)
for _sym in ("CHMP1A", "VPS9D1-AS1", "SPATA33"):
    for _, r in _cl[_cl.SYMBOL == _sym].iterrows():
        s = "%.2f" % r["PP.H3"]
        hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   PP.H3 %s" % ("OK " if hit else "MISSING", s, _sym))
_hd = pd.read_csv(_find("08_SMR_HEIDI_results.tsv"), sep=TAB)
for _sym, _n in (("VPS9D1-AS1", 1), ("CDK10", 2)):
    _v = sorted(_hd[_hd.SYMBOL == _sym].p_HEIDI.dropna(), reverse=True)[:_n]
    for _x in _v:
        s = "%.3f" % _x
        hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   P_HEIDI %s" % ("OK " if hit else "MISSING", s, _sym))

print()
print("=" * 74)
print("1v. the voided RA whole-blood cell (S33), against 123d")
print("=" * 74)
# "3.94-fold (P = 3.3e-17), but its pre-registered mismatched-list control also
# enriches -- melanoma's known loci give 2.05-fold on the RA nominations,
# P = 2.9e-4 -- ... excluding the MHC does not clean it (1.83-fold, P = 0.0052)".
# Six numbers explaining why a cell is VOID. A paper is least likely to be
# reread where it rules against itself, which is why these get pinned.
# The MHC-excluded rows are a separate CELL, suffixed, not a filter on
# the same one -- "RA x eQTLGen_blood (MHC excluded)". Matching on the
# bare name silently found nothing and printed four items instead of six.
_g = grid[grid.cell.str.startswith("RA x eQTLGen_blood")
          & (grid.analysis == "main")]
_ag = _g[_g.region_filter == "all_genome"]
_mh = _g[_g.region_filter == "MHC_excluded"]
_items = []
if not _ag.empty:
    r = _ag.iloc[0]
    _items += [("own fold", "%.2f" % r.fold), ("own P", _num_in_text(r.fisher_p)),
               ("mismatch fold", "%.2f" % r.mismatch_fold),
               ("mismatch P", _num_in_text(r.mismatch_p))]
if not _mh.empty:
    r = _mh.iloc[0]
    _items += [("MHC-excluded mismatch fold", "%.2f" % r.mismatch_fold),
               ("MHC-excluded mismatch P", _num_in_text(r.mismatch_p))]
for label, s in _items:
    ok = s is not None and re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not ok:
        bad += 1
    print("  %s  %-8s   %s" % ("OK " if ok else "MISSING", s or "?", label))

print()
print("=" * 74)
print("1w. the attribution ceilings (S43), against 142a")
print("=" * 74)
# "the share runs from 0.101 for melanoma to 0.403 for prostate -- ceilings of
# 9.92 and 2.48". The ceiling is the reciprocal of the background share, so it
# is fully determined by a number already in the table; quoting it is a place
# an arithmetic slip would never be recomputed.
_pmb = pd.read_csv(_find("142a_power_matched.tsv"), sep=TAB)
_base = _pmb[_pmb.matched_k.isna() | (_pmb.matched_k.astype(str) == "NA")]
for _out in ("Melanoma", "Prostate"):
    _r = _base[_base.outcome == _out]
    if _r.empty:
        continue
    _pb = float(_r.iloc[0].p_bg)
    for label, s in ((_out + " share", "%.3f" % _pb),
                     (_out + " ceiling", "%.2f" % (1.0 / _pb))):
        ok = re.search(re.escape(s) + r"(?!\d)", txt) is not None
        if not ok:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if ok else "MISSING", s, label))

print()
print("=" * 74)
print("1t. the Steiger R-squared bounds (Methods), against 09")
print("=" * 74)
# "TwoSampleMR's R2 formula for SD units is unbounded and produced values above
# 1 (maximum 1.027), so the bounded form R2 = F/(F + N - 2) (range
# 0.185-0.933)". The whole point of the sentence is that one formula breaks its
# own bound, so the number that demonstrates it had better be the real maximum.
_st = pd.read_csv(_find("09_steiger_filtering.tsv"), sep=TAB)
for label, v in (("unbounded maximum", _st["rsq.exposure"].max()),
                 ("bounded minimum", _st["rsq.exposure.bounded"].min()),
                 ("bounded maximum", _st["rsq.exposure.bounded"].max())):
    s = "%.3f" % v
    hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, label))

print()
print("=" * 74)
print("1u. the melanoma cell's multi-list margin (S39), against 130a and 130c")
print("=" * 74)
# "breast at 2.15-fold (P = 0.0495), so its multi-list verdict moves from no
# enriching comparator to one; the cell's margin over its best rival widens from
# 1.82- to 2.31-fold". This is the paper reporting a verdict CHANGE against
# itself, so the numbers behind the change are the ones worth pinning.
_mlm = pd.read_csv(_find("130a_multilist_main.tsv"), sep=TAB)
_mel = _mlm[(_mlm.cell == "melanoma x Soskic_CD4") & (_mlm.role == "main")
            & (_mlm.list_name == "breast")]
if not _mel.empty:
    r = _mel.iloc[0]
    for label, s in (("breast fold", "%.2f" % r.fold),
                     ("breast P", _num_in_text(r.fisher_p) or "%.4f" % r.fisher_p)):
        hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, label))
_mlv = pd.read_csv(_find("130c_multilist_verdict.tsv"), sep=TAB)
_mv = _mlv[_mlv.cell == "melanoma x Soskic_CD4"]
if not _mv.empty:
    s = "%.2f" % _mv.iloc[0].margin_all
    hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   melanoma margin over best rival"
          % ("OK " if hit else "MISSING", s))

print()
print("=" * 74)
print("1r. every main-grid P value (S42), against 123d")
print("=" * 74)
# Section 1 checks the grid FOLDS. Their P values -- 0.0041 for HCC-low,
# 0.123 for HCC-high, 0.0015 on whole blood -- were never checked, so a fold
# could stay right while the significance attached to it drifted. Cells the
# manuscript does not quote are listed, not demanded, as in section 1.
_unq_p = []
for _, r in main.iterrows():
    if "MHC" in r.cell or r.control == "FAILED" or pd.isna(r.fisher_p):
        continue
    m = _num_in_text(r.fisher_p)
    if m:
        print("  OK   %-10s   P for %s" % (m, r.cell))
    else:
        _unq_p.append((r.cell, "%.4g" % r.fisher_p))
if _unq_p:
    print("  not quoted (not a failure):")
    for _c, _v in _unq_p:
        print("      %-32s %s" % (_c, _v))

print()
print("=" * 74)
print("1s. the RA sweeps quoted as ranges (S37), against 128a and 128b")
print("=" * 74)
# "flat across locus widths (3.70- to 3.96-fold)" and "rises as the known-locus
# window narrows (3.91- to 7.30-fold)". A range is quoted by its endpoints, so
# the endpoints are what get checked -- and a range is the easy place to widen
# a claim by a digit without anyone recomputing it.
for _f, _col, _lab in ((("128b_locus_kb_sweep.tsv"), "locus_kb", "locus width"),
                       (("128a_known_kb_sweep.tsv"), "known_kb", "known window")):
    _t = pd.read_csv(_find(_f), sep=TAB)
    _t = _t[_t.cell == "RA x Soskic_CD4"]
    if "role" in _t.columns:
        _t = _t[_t.role == "main"]
    for _v, _end in ((_t.fold.min(), "min"), (_t.fold.max(), "max")):
        s = "%.2f" % _v
        hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   RA %s %s" % ("OK " if hit else "MISSING", s, _lab, _end))

print()
print("=" * 74)
print("1p. the mismatched-list window sweep (S33/S37), against 128a")
print("=" * 74)
# "1.71-, 2.38-, 3.03- and 4.17-fold at 1000, 500, 250 and 100 kb" for RA, and
# "4.11-, 3.44-, 2.55-, 1.59-fold as the window widens" for melanoma on whole
# blood. Eight numbers arguing AGAINST the paper's own control, which is the
# last place a quiet transposition should be allowed to sit.
_sw = pd.read_csv(_find("128a_known_kb_sweep.tsv"), sep=TAB)
for _cell in ("RA x Soskic_CD4", "melanoma x eQTLGen_blood"):
    _c = _sw[(_sw.cell == _cell) & (_sw.role == "main")]
    for _, r in _c.sort_values("known_kb").iterrows():
        if pd.isna(r.mismatch_fold) or r.mismatch_fold <= 0:
            continue
        s = "%.2f" % r.mismatch_fold
        hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   %s at %d kb"
              % ("OK " if hit else "MISSING", s, _cell, int(r.known_kb)))

print()
print("=" * 74)
print("1q. the RA multi-list control (S39), against 130a")
print("=" * 74)
# "melanoma ... at 1.71-fold (P = 0.10) ... HCC 3.05 (0.0028), lung 3.32
# (1.5e-5), colorectal 1.71 (0.0018), breast 1.69 (0.0042), prostate 1.39
# (0.026) would each have voided it". Twelve numbers deciding whether the RA
# cell survives its negative control, and the decision turns on which
# comparator was nominated -- so every comparator's pair is checked, not the
# nominated one alone.
_ml = pd.read_csv(_find("130a_multilist_main.tsv"), sep=TAB)
_ra = _ml[(_ml.cell == "RA x Soskic_CD4") & (_ml.role == "main")]
for _, r in _ra.iterrows():
    if r.fold <= 0:
        continue
    # The text renders P at whichever precision suits it -- 0.10 here, 0.026
    # there, 1.5e-5 elsewhere -- so a single format string reports MISSING on
    # numbers that are perfectly correct. Any of the precisions the paper
    # actually uses counts; the digits still have to be right.
    if r.fisher_p >= 5e-5:
        # A candidate that rounds to all zeros carries no information and
        # would match "0.00" anywhere in the text -- kappa 0.00, for one -- so
        # every P below 0.005 would pass on nothing. Dropped.
        _cands = [c for c in ("%.*f" % (d, r.fisher_p) for d in (2, 3, 4))
                  if float(c) > 0]
    else:
        _cands = ["%.1f" % (r.fisher_p / 10 ** math.floor(math.log10(r.fisher_p)))]
    _match = next((c for c in _cands
                   if re.search(re.escape(c) + r"(?!\d)", txt)), None)
    if _match is None:
        bad += 1
    print("  %s  %-8s   RA vs %s (P)"
          % ("OK " if _match else "MISSING", _match or _cands[-1], r.list_name))
    _fs = "%.2f" % r.fold
    _fok = re.search(re.escape(_fs) + r"(?!\d)", txt) is not None
    if not _fok:
        bad += 1
    print("  %s  %-8s   RA vs %s (fold)"
          % ("OK " if _fok else "MISSING", _fs, r.list_name))

print()
print("=" * 74)
print("1o. the |z|-matched residual (S28), against 101b")
print("=" * 74)
# "+5.2 percentage points [-1.6, +12.0] by bounded locus and +6.4 [+0.9, +12.0]
# by gene, against raw gaps of 47.6 and 63.0". Eight numbers, none audited, and
# they are the ones that decide whether the class gap survives conditioning on
# effect size -- the residual IS the claim in that paragraph.
_mt = pd.read_csv(_find("101b_matched.tsv"), sep=TAB)
for _, r in _mt.iterrows():
    for s in ("%.1f" % (100 * r.gap_observed), "%.1f" % (100 * r.delta_matched),
              "%.1f" % (100 * r.ci_lo), "%.1f" % (100 * r.ci_hi)):
        hit = re.search(re.escape(s.lstrip("-")) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, r.unit))

print()
print("=" * 74)
print("1n. the Software section, against 150a and requirements.txt")
print("=" * 74)
# The Methods name twenty package versions. None was audited, and a version
# claim is exactly the kind of number that rots quietly: nothing recomputes it,
# so it stays as first typed while the environment moves. Checked against the
# deposited lock (R) and the pinned requirements (Python), normalising the
# lock's dashes to the dots the text uses (1.7-0 -> 1.7.0).
_SOFT_R = {
    "Seurat": "5.5.1", "Matrix": "1.7.0", "data.table": "1.16.0",
    "coloc": "5.2.3", "susieR": "0.14.2",
    "arrow": "25.0.0", "hdf5r": "1.3.12", "TFBSTools": "1.42.0",
    "JASPAR2020": "0.99.10", "motifmatchr": "1.26.0", "chromVAR": "1.26.0",
    "BSgenome.Hsapiens.UCSC.hg38": "1.4.5", "survival": "3.8.9",
    "org.Hs.eg.db": "3.19.1",
}
_SOFT_PY = {"numpy": "2.0.0", "pandas": "2.2.2", "scipy": "1.18.0",
            "pyarrow": "25.0.0", "matplotlib": "3.11.1"}
_LOCK_GAPS = []
_lockp = json.load(io.open(_find("150a_environment.lock"),
                           encoding="utf-8"))["Packages"]
for _pk, _pv in sorted(_SOFT_R.items()):
    _got = _lockp.get(_pk, {}).get("Version")
    if _got is None:
        # A DOCUMENTED gap, not an error: chromVAR and org.Hs.eg.db are named
        # in the Methods, absent from 150a, and present in the container at
        # exactly the stated versions -- installed by install_r_packages.R's
        # unpinned fallback, so right by luck rather than by record. The lock
        # cannot simply be regenerated to fix this: step150 would rebuild it
        # from the CURRENT machine, which has drifted (numpy 2.5.2 against the
        # recorded 2.0.0), and that would corrupt the very record S54 compares
        # against. Recorded in DEPOSIT_GAPS.md and reported here every run so
        # it stays visible, but not failed -- an audit permanently red on a
        # gap nobody can close today is one people stop reading.
        _LOCK_GAPS.append((_pk, _pv))
        print("  GAP      %-28s text %-9s not in 150a (see DEPOSIT_GAPS.md)"
              % (_pk, _pv))
    else:
        _ok = _got.replace("-", ".") == _pv
        if not _ok:
            bad += 1
        print("  %s  %-28s text %-9s lock %s"
              % ("OK " if _ok else "MISSING", _pk, _pv, _got))
_req = io.open(_find("container", "requirements.txt"), encoding="utf-8").read()
for _pk, _pv in sorted(_SOFT_PY.items()):
    _m = re.search(r"^%s==(\S+)$" % re.escape(_pk), _req, re.M)
    _ok = _m is not None and _m.group(1) == _pv
    if not _ok:
        bad += 1
    print("  %s  %-28s text %-9s req %s"
          % ("OK " if _ok else "MISSING", _pk, _pv,
             _m.group(1) if _m else "ABSENT"))

print()
print("=" * 74)
print("1l. the five-outcome transport result (S41), against 138b")
print("=" * 74)
# "melanoma 4.96-fold (P = 0.0097 corrected), lung 3.85-fold (0.0138),
# colorectal 3.15-fold (0.0029), breast 2.25-fold (0.0011) and prostate
# 2.04-fold (2.3e-5)". Ten numbers carrying the claim that the attribution
# recurs across outcomes, of which one was audited. Driven off the table, so an
# outcome cannot be added or dropped without this following.
_tr = pd.read_csv(_find("138b_transport_verdict.tsv"), sep=TAB)
for _, r in _tr.iterrows():
    if pd.isna(r.fold):
        continue
    _ps = "%.4f" % r.fisher_p_holm
    if _ps.startswith("0.0000"):
        # Too small for four decimals; the text writes it as "2.3 x 10-5", so
        # check the mantissa it actually prints rather than skipping the row.
        _ps = "%.1f" % (r.fisher_p_holm / 10 ** math.floor(
            math.log10(r.fisher_p_holm)))
    for s in ("%.2f" % r.fold, _ps):
        hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, r.row))

print()
print("=" * 74)
print("1m. fold versus chance-corrected agreement at ten loci (S43), against 142a")
print("=" * 74)
# "at ten loci, fold gives melanoma 5.95, lung 3.85, breast 3.10 and colorectal
# 2.52, while A gives breast 0.859, colorectal 0.707, melanoma 0.555 and lung
# 0.330". This is the sentence that says the two statistics rank the outcomes
# almost in reverse, so it is exactly the place a transposed pair would hide.
_pmx = pd.read_csv(_find("142a_power_matched.tsv"), sep=TAB)
_k10 = _pmx[_pmx.matched_k == 10]
for _, r in _k10.iterrows():
    # NaN fails every comparison, so "<= 0" let prostate through and "%.2f"
    # rendered it as the string "nan", which then matched text and reported OK.
    # A check that passes on a missing value is worse than no check.
    if pd.isna(r.fold) or pd.isna(r.A_chance_corrected) or r.fold <= 0:
        continue                          # prostate has no matched draw at k=10
    for s in ("%.2f" % r.fold, "%.3f" % r.A_chance_corrected):
        hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
        if not hit:
            bad += 1
        print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, r.outcome))

print()
print("=" * 74)
print("1j. the down-sampling calibration (S28), against 53b")
print("=" * 74)
# "10.2 simulated versus 10 real discoveries; simulated Jaccard 0.52
# [0.36, 0.73] containing the observed 0.455". The calibration is what licenses
# every recovery curve in Fig. 3, and none of its six numbers was audited.
_cal = pd.read_csv(_find("53b_calibration.tsv"), sep=TAB).iloc[0]
cals = [
    ("real discoveries", "%d" % int(_cal.real_hits)),
    ("simulated discoveries", "%.1f" % _cal.sim_hits),
    ("simulated Jaccard", "%.2f" % _cal.sim_jaccard),
    ("Jaccard lower", "%.2f" % _cal.sim_lo),
    ("Jaccard upper", "%.2f" % _cal.sim_hi),
    ("observed Jaccard", "%.3f" % _cal.real_jaccard),
]
for label, s in cals:
    hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, label))

print()
print("=" * 74)
print("1k. the meta-round colocalisation comparison, recomputed from 07 vs 14")
print("=" * 74)
# "across the 127 exposures run in both rounds, median PP.H3+H4 fell from 0.249
# to 0.207 and median PP.H4 from 0.124 to 0.095, with only 44.1% improving".
# Stored nowhere -- it is a join between the single-round and meta-round coloc
# tables -- so no audit comparing text against cells could reach it. Recomputing
# checks the claim against the coloc output itself.
#
# "Improving" means PP.H3+H4, not PP.H4: on H4 the figure is 40.2%. Established
# by testing both against the text rather than assuming which was meant, and
# recorded here so the next reader does not have to redo it.
_c1 = pd.read_csv(_find("07_coloc_results.tsv"), sep=TAB)
_c2 = pd.read_csv(_find("14_coloc_meta_results.tsv"), sep=TAB)
_keys = [k for k in ("gene_id", "exposure", "cell_type", "timepoint")
         if k in _c1.columns and k in _c2.columns]
_mg = _c1.merge(_c2, on=_keys, suffixes=("_1", "_2"))
_h1 = _mg["PP.H3_1"] + _mg["PP.H4_1"]
_h2 = _mg["PP.H3_2"] + _mg["PP.H4_2"]
cocs = [
    ("exposures in both rounds", "%d" % len(_mg)),
    ("median H3+H4 before", "%.3f" % _h1.median()),
    ("median H3+H4 after", "%.3f" % _h2.median()),
    ("median H4 before", "%.3f" % _mg["PP.H4_1"].median()),
    ("median H4 after", "%.3f" % _mg["PP.H4_2"].median()),
    ("share improving on H3+H4", "%.1f%%" % (100.0 * (_h2 > _h1).mean())),
]
for label, s in cocs:
    hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, label))

print()
print("=" * 74)
print("1g. the eQTLGen attribution sentence (S42), against 123d")
print("=" * 74)
# The sentence carrying the second exposure resource: "23 of 34 significant loci
# (67.6%) ... 6.31-fold over this resource's own 10.7% background". Only 6.31 was
# audited; the counts, the percentage and the background it is a fold OVER were
# not, so the fold could have stayed correct while its own denominator drifted.
_eq = main[main.cell == "melanoma x eQTLGen_blood"].iloc[0]
eqs = [
    ("significant loci", "%d" % int(_eq.sig_loci)),
    ("of them known", "%d" % int(_eq.sig_known)),
    ("share known", "%.1f%%" % float(_eq.pct_known)),
    ("resource background", "%.1f%%" % (100.0 * _eq.bg_known / _eq.bg_loci)),
]
for label, s in eqs:
    hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, label))

print()
print("=" * 74)
print("1h. the |z| distribution behind the class gap (S28), against 101a")
print("=" * 74)
# "their full-power |z| all lie between 3.73 and 4.55, while known-locus
# candidates reach 15.99 (medians 4.23 and 11.07)". These five numbers carry the
# argument that the known-vs-novel gap is an effect-size statement rather than a
# class one, and none of them was audited.
_zd = pd.read_csv(_find("101a_zdist.tsv"), sep=TAB)
_zg = _zd[_zd.unit == "gene"]
_zk = _zg[_zg.cls == "known"].iloc[0]
_zn = _zg[_zg.cls == "novel"].iloc[0]
zs = [
    ("novel |z| lower", "%.2f" % _zn.lo),
    ("novel |z| upper", "%.2f" % _zn.hi),
    ("known |z| upper", "%.2f" % _zk.hi),
    # .median is DataFrame.median, the method -- reach it by name.
    ("novel median", "%.2f" % _zn["median"]),
    ("known median", "%.2f" % _zk["median"]),
]
for label, s in zs:
    hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s" % ("OK " if hit else "MISSING", s, label))

print()
print("=" * 74)
print("1i. two-coder agreement (S23), recomputed from 148b")
print("=" * 74)
# "Agreement was 97.8% on each criterion with one disagreement each". This is not
# stored anywhere -- it is derived from the two coders' columns -- so it was
# unreachable by any audit that only compares text against table cells. Deriving
# it here means the claim is checked against the coding itself rather than
# against a number someone once wrote down.
_dc = pd.read_csv(_find("148b_litaudit_doublecoded.tsv"), sep=TAB)
for crit, g in _dc.groupby("criterion"):
    agree = (g["manual"] == g["coder2"])
    s = "%.1f%%" % (100.0 * agree.mean())
    hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-8s   %s, %d of %d, %d disagreement(s)"
          % ("OK " if hit else "MISSING", s, crit, int(agree.sum()), len(g),
             int((~agree).sum())))
    # Cohen's kappa alongside the raw rate, recomputed the same way. The paper
    # reports 0.00 for locus attribution as a PREVALENCE ARTEFACT -- with one
    # positive in 46 the expected agreement equals the observed -- so the zero
    # is a claim about the statistic, not a coding failure, and it should be
    # reproduced rather than taken on trust.
    _a = g["manual"].to_numpy()
    _b = g["coder2"].to_numpy()
    _po = float((_a == _b).mean())
    _cats = sorted(set(_a) | set(_b))
    _pe = sum(float((_a == c).mean()) * float((_b == c).mean()) for c in _cats)
    _k = (_po - _pe) / (1 - _pe) if _pe < 1 else 0.0
    _ks = "%.2f" % _k
    _khit = re.search(re.escape(_ks) + r"(?!\d)", txt) is not None
    if not _khit:
        bad += 1
    print("  %s  %-8s   Cohen kappa, %s"
          % ("OK " if _khit else "MISSING", _ks, crit))

print()
print("=" * 74)
print("1f. the density-matched permutation sentence (S38), against 126c")
print("=" * 74)
# The Results sentence beginning "Matching on density" quotes six folds. Two
# were audited; of the four that were not, one came from the WRONG TABLE --
# the text read 13.53 for HCC_low x eQTLGen_blood, which is 123d's unmatched
# grid value, where the density-matched permutation gives 12.09. Five matched
# and the sixth did not, through every green audit and the whole container
# acceptance, because nothing read it.
#
# The test runs from the SENTENCE outwards, not from the table: every fold the
# sentence quotes must correspond to a row of 126c. Checking the other way --
# demanding the text quote every row -- fails on HCC_high x Soskic_CD4, which
# is computed and deliberately not quoted, and an audit that fires on a
# correct omission is one people learn to silence.
_flat = re.sub(r"\s+", " ", txt)
_i = _flat.find("Matching on density")
if _i < 0:
    print("  MISSING  the 'Matching on density' sentence is gone -- if it was")
    print("           rewritten, re-anchor this check rather than deleting it.")
    bad += 1
else:
    _sent = _flat[_i:_i + 420]
    _perm = pd.read_csv(_find("126c_permutation_primary.tsv"), sep=TAB)
    _folds = set("%.2f" % r.fold_vs_null for _, r in _perm.iterrows())
    _quoted = re.findall(r"(?<![\w.-])\d+\.\d+(?![\w])", _sent)
    _quoted = [q for q in _quoted if float(q) >= 1.0]   # P values are not folds
    for q in _quoted:
        ok = q in _folds
        if not ok:
            bad += 1
        print("  %s  %-8s   quoted fold %s in 126c"
              % ("OK " if ok else "MISSING", q,
                 "found" if ok else "NOT FOUND"))
    # The sentence also quotes the empirical P for the melanoma cell. Checked
    # against the same table so a P and the fold it belongs to cannot drift
    # apart -- which is exactly how 13.53 survived: right table for five
    # numbers, wrong table for the sixth.
    _mp = _perm[_perm.cell == "melanoma x Soskic_CD4"]
    if not _mp.empty:
        _ep = "%.4f" % float(_mp.iloc[0].empirical_p)
        _eok = re.search(re.escape(_ep) + r"(?!\d)", txt) is not None
        if not _eok:
            bad += 1
        print("  %s  %-8s   empirical P, melanoma x Soskic_CD4"
              % ("OK " if _eok else "MISSING", _ep))
    _unq = sorted(f for f in _folds if f not in _quoted)
    print("  (computed but not quoted, not a failure: %s)" % ", ".join(_unq))
print()
print("=" * 74)
print("1e. the power-stratified recovery gap (S28), against 55a")
print("=" * 74)
# S54 section 2 names step101's positive control -- "reproduces 85.8% / 22.8%"
# -- as a tier-two criterion the container must satisfy. It was never audited:
# step160 measured coverage at 22 of 191 distinct decimals and these were among
# the 169 nothing looked at. They are also the numbers a reader meets first in
# the power paragraph, and they carry the claim that the recovery gap is a
# property of locus class.
_rec = pd.read_csv(_find("55a_recovery_by_locus_class.tsv"), sep=TAB)


def _recrow(frac):
    r = _rec[(_rec.frac - frac).abs() < 1e-9]
    if r.empty:
        raise SystemExit("55a has no frac=%g row" % frac)
    return r.iloc[0]


_lo, _hi = _recrow(0.1), _recrow(0.5)
recov = [
    ("known recovery at 10% power", "%.1f%%" % (100 * _lo.known_recovery)),
    ("novel recovery at 10% power", "%.1f%%" % (100 * _lo.novel_recovery)),
    ("known recovery at 50% power", "%.1f%%" % (100 * _hi.known_recovery)),
    ("novel recovery at 50% power", "%.1f%%" % (100 * _hi.novel_recovery)),
]
for label, s in recov:
    hit = re.search(re.escape(s) + r"(?!\d)", txt) is not None
    if not hit:
        bad += 1
    print("  %s  %-10s   %s" % ("OK " if hit else "MISSING", s, label))

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
