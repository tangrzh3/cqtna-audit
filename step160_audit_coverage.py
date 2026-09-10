#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 160 -- how much of the manuscript's arithmetic is actually audited.

WHY THIS EXISTS
step127 reconciles 41 numbers against 8 named tables. The main text contains
roughly 191 distinct decimal values. Nothing measured that gap, so nobody knew
it was there -- and the gap is not hypothetical: the matched-background fold
in the eQTLGen paragraph sat inside it, was produced by a rule with no unique
answer, gave three different values in three environments, and every audit ran
green throughout (S54 sections 9.8, 9.10.6). It surfaced because a portability
fix happened to let the producing script run, not because anything checked.

WHAT IT MEASURES, AND WHAT IT DOES NOT
It reports which numbers in the text are reconciled by step127 and which are
not. It does NOT claim the unreconciled ones are wrong -- only that no audit
looks at them. That distinction is the whole point: "not checked" is a fact
about this repository, "wrong" would be a claim about the science.

It deliberately does not try to guess provenance by searching the tables for a
matching value. 360 result files contain some 640,000 distinct rounded forms,
so almost any two- or three-decimal number matches something by coincidence.
A coincidence is not a reconciliation, and a coverage report built on one
would be worse than none: it would read as reassurance. The one search that
does carry signal is the negative -- a number that appears in NO table cannot
have come from one, so it is derived, quoted from a source, or unbacked, and
those are listed separately.

  python step160_audit_coverage.py [dir]

Outputs: 160a_audit_coverage.tsv   every text number, covered or not
         160b_console.log
"""
import glob
import io
import os
import re
import subprocess
import sys

MR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CQTNA_DIR") or \
     os.path.dirname(os.path.abspath(__file__))
LOG = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)



# Values that CANNOT be reconciled against a result table, with the reason.
# Curated by hand and deliberately not inferred: calling something
# un-auditable is a claim, and a script that guessed would quietly excuse
# whatever it failed to find. Anything not listed here is simply not yet
# audited, which is a different statement.
UNAUDITABLE = {
    # Thresholds and parameters -- inputs to the analysis, not outputs of it.
    "0.05": ("threshold", "the BH level itself"),
    "0.8": ("threshold", "PP.H4 > 0.8, a stated cutoff"),
    "0.5": ("threshold", "PP.H3+PP.H4 > 0.5, an additional requirement"),
    "0.70": ("threshold", "pre-registered kill criterion for the glycolysis axis"),
    "0.1": ("parameter", "--clump-r2 0.1"),
    "10.0": ("threshold", "the F > 10 filter, quoted as 10.0-fold elsewhere"),
    # Quoted from outside this study; no local table is their source.
    "45.6": ("external", "FinnGen's own published fine-mapping, log10BF"),
    "36.6": ("external", "FinnGen's own published fine-mapping, log10BF"),
    "19.5": ("external", "FinnGen's own published fine-mapping, log10BF"),
    "0.015631": ("external", "Rashkin case fraction, GWAS metadata"),
    "0.014962": ("external", "outcome prevalence, GWAS metadata"),
    "5.45": ("definitional", "|b| = 5.45 x SE is the 5e-8 boundary, not a result"),
    # Software versions, checked by step127 section 1n against the lock; these
    # tokens are fragments of a version string rather than quantities.
    "4.4": ("version", "fragment of R 4.4.1"),
    "3.19": ("version", "fragment of org.Hs.eg.db 3.19.1"),
    "3.12": ("version", "fragment of Python 3.12.4"),
    "25.0": ("version", "fragment of arrow/pyarrow 25.0.0"),
    # Arithmetic on numbers that ARE audited, carrying no independent content.
    "86.9": ("derived", "253/291, both quoted in the same sentence"),
    "3.5": ("derived", "8771/2506, both quoted in the same sentence"),
    "1.6": ("derived", "the lower end of [-1.6, +12.0], checked in step127 1o"),
}

def manuscript_numbers(path):
    """Decimal values in the body, with the line they sit on.

    Citation brackets are stripped first: "[12, 15]" is not arithmetic. Integers
    are left out on purpose -- years, reference numbers, section numbers and
    sample sizes are indistinguishable from result counts by pattern alone, and
    a report padded with those would bury the numbers that matter.
    """
    t = io.open(path, encoding="utf-8").read()
    body = t.split("## Declarations")[0]
    out = []
    for i, line in enumerate(body.split("\n"), 1):
        clean = re.sub(r'\[\d+(?:[,\-\u2013]\s*\d+)*\]', ' ', line)
        for m in re.finditer(r'(?<![\w.\-])\d+\.\d+(?![\w])', clean):
            out.append((m.group(0), i, line.strip()))
    return out


def step127_checked():
    """The renderings step127 verifies, taken from its own output.

    Parsed rather than re-derived so this cannot drift away from what step127
    actually does; if step127 gains a check, this sees it on the next run.
    """
    try:
        p = subprocess.run(["python", "step127_audit_manuscript_numbers.py", MR],
                           cwd=MR, capture_output=True, text=True, timeout=900)
    except Exception as e:
        say("  cannot run step127 (%s); coverage cannot be measured." % e)
        return None
    vals = set()
    for ln in (p.stdout or "").splitlines():
        m = re.match(r'\s+(?:OK |MISSING)\s+(.*)$', ln)
        if m:
            for tok in re.findall(r'(?<![\w.\-])\d+\.\d+(?![\w])', m.group(1)):
                vals.add(tok)
    return vals


def table_forms():
    """Every numeric cell in every result table, as rounded strings.

    Used ONLY to answer "does this value exist anywhere at all", which is
    informative when the answer is no.
    """
    forms = set()
    numeric = re.compile(r'^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$')
    for fn in sorted(glob.glob(os.path.join(MR, "*.tsv"))) + \
              sorted(glob.glob(os.path.join(MR, "*.csv"))):
        sep = "\t" if fn.endswith("tsv") else ","
        try:
            with io.open(fn, encoding="utf-8", errors="replace") as f:
                for line in f:
                    for cell in line.rstrip("\n").split(sep):
                        c = cell.strip().strip('"')
                        if c and numeric.match(c):
                            try:
                                v = float(c)
                            except Exception:
                                continue
                            for d in range(0, 5):
                                s = "%.*f" % (d, v)
                                forms.add(s)
                                if s.startswith("0."):
                                    forms.add(s[1:])
        except Exception:
            continue
    return forms


def main():
    ms = os.path.join(MR, "manuscript", "MANUSCRIPT_GB.md")
    if not os.path.exists(ms):
        say("manuscript not found: %s" % ms)
        return 1

    nums = manuscript_numbers(ms)
    uniq = sorted(set(n for n, _, _ in nums), key=lambda x: -float(x))
    say("=" * 74)
    say("audit coverage of the manuscript's arithmetic")
    say("=" * 74)
    say("  decimal values in the body : %d occurrences, %d distinct"
        % (len(nums), len(uniq)))

    checked = step127_checked()
    if checked is None:
        return 1
    say("  reconciled by step127      : %d distinct" % len(checked & set(uniq)))

    forms = table_forms()
    say("  numeric forms in 360 result files: %d (coincidence is cheap here,"
        % len(forms))
    say("    which is why presence is not treated as evidence)")
    say()

    rows = []
    unbacked, unchecked = [], []
    first_line = {}
    for n, ln, ctx in nums:
        first_line.setdefault(n, (ln, ctx))
    for n in uniq:
        ln, ctx = first_line[n]
        if n in checked:
            cov = "step127"
        elif n in UNAUDITABLE:
            cov = "unauditable:" + UNAUDITABLE[n][0]
        else:
            cov = "NOT AUDITED"
        back = "yes" if n in forms else "NO TABLE"
        rows.append((n, cov, back, ln, ctx[:120]))
        if n not in checked and n not in UNAUDITABLE:
            unchecked.append(n)
        if n not in forms:
            unbacked.append((n, ln, ctx))

    _rec = len(checked & set(uniq))
    _una = len([n for n in uniq if n in UNAUDITABLE])
    _pend = len(unchecked)
    _auditable = len(uniq) - _una
    say("  reconciled          : %d" % _rec)
    say("  un-auditable        : %d  (threshold, external, derived, version)" % _una)
    say("  auditable, not done : %d" % _pend)
    say("  COVERAGE of what CAN be reconciled: %d of %d (%.0f%%)"
        % (_rec, _auditable, 100.0 * _rec / _auditable if _auditable else 0))
    say()
    say("  Un-auditable is a curated claim, not an inference -- see UNAUDITABLE")
    say("  in this file. A value is listed there only with a stated reason, so")
    say("  the category cannot quietly absorb whatever the audit failed to find.")
    say()
    say("  %d value(s) appear in NO result table -- derived, quoted, or unbacked:"
        % len(unbacked))
    for n, ln, ctx in unbacked:
        say("    %-12s line %-5d %s" % (n, ln, ctx[:88]))

    with io.open(os.path.join(MR, "160a_audit_coverage.tsv"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write("value\tcoverage\tin_any_table\tfirst_line\tcontext\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    say()
    say("  wrote 160a_audit_coverage.tsv (%d rows)" % len(rows))
    say()
    say("  This is a coverage report, not a verdict. A row marked NOT AUDITED")
    say("  is not a row known to be wrong -- it is a row nothing looks at.")
    return 0


if __name__ == "__main__":
    rc = main()
    io.open(os.path.join(MR, "160b_console.log"), "w",
            encoding="utf-8", newline="\n").write("\n".join(LOG) + "\n")
    sys.exit(rc)
