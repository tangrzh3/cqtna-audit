#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 154 -- the cross-references, which step127 cannot see.

WHY
step127 reconciles numbers against tables. It found nothing wrong on the day a
manual sweep found six things wrong, because none of them were numbers:

  * seven supplementary items existed in the numbering and were cited nowhere;
  * five of the nine figures had no in-text callout at all;
  * none of the nine figures had a panel legend, and the figures carry
    twenty-five panels between them;
  * ten references were listed and never cited;
  * one caveat appeared verbatim in both Results and Methods;
  * the manuscript had no author block and no declarations.

Every one of them has the same cause: text moved into the supplement during
compression and its pointers stayed behind. That will happen again -- the
compression plan still has candidates in it -- so the sweep is a script now.

  python step154_crossref_audit.py [repo]

Exits non-zero if any check fails.
"""
import io
import os
import re
import sys

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(MR, "manuscript", "MANUSCRIPT_GB.md")
FIGDIR = os.path.join(MR, "figures")

BAD = []


def fail(section, msg):
    BAD.append((section, msg))
    print("  FAIL  " + msg)


def ok(msg):
    print("  ok    " + msg)


def head(n):
    print()
    print("=" * 78)
    print(n)
    print("=" * 78)


def read(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


def sections(t):
    """body = everything a reader reads before the reference list; the figure
    legends and the supplementary index are addressed separately."""
    return {
        "all": t,
        # Every place a citation can legitimately live. Methods sits AFTER
        # the figure legends in this file, so a scan that stops at
        # "## Figures" misses most citations and calls them orphans.
        "cite": t[:t.index("## References")],
        "body": t[:t.index("## Figures")],
        "figs": t[t.index("## Figures"):t.index("## Methods")],
        "supp": t[t.index("## Supplementary information"):t.index("## References")],
        "refs": t[t.index("## References"):],
    }


# ----------------------------------------------------------------- 1. supplements
def check_supplements(S):
    head("1. supplementary items: declared, cited, both")
    m = re.search(r"S(\d+)\s*[–-]\s*S(\d+)", S["supp"])
    if not m:
        fail("supp", "the supplementary index declares no S-number range")
        return
    lo, hi = int(m.group(1)), int(m.group(2))
    declared = set(range(lo, hi + 1))
    cited = set(int(x) for x in re.findall(r"\bS(\d+)\b",
                                           S["cite"]))
    missing = sorted(declared - cited)
    beyond = sorted(cited - declared)
    if missing:
        fail("supp", "declared S%d-S%d but never cited: %s"
             % (lo, hi, ", ".join("S%d" % n for n in missing)))
    else:
        ok("every item in the declared range S%d-S%d is cited" % (lo, hi))
    if beyond:
        fail("supp", "cited but outside the declared range: %s"
             % ", ".join("S%d" % n for n in beyond))
    else:
        ok("no citation points outside the declared range")


# --------------------------------------------------------------------- 2. figures
def script_panels():
    """Panel letters as the figure scripts actually draw them. A gb script may
    reuse another script's panel functions (Fig. 3 borrows a/b/c from the
    power-stability figure, Fig. 8 is the compartment figure renumbered), so
    local imports are followed one level."""
    out = {}
    if not os.path.isdir(FIGDIR):
        return out
    for fn in sorted(os.listdir(FIGDIR)):
        m = re.match(r"make_gb_fig(\d+)_.*\.py$", fn)
        if not m:
            continue
        n = int(m.group(1))
        src = read(os.path.join(FIGDIR, fn))
        texts = [src]
        for imp in re.findall(r"^import\s+(make_[a-z0-9_]+)", src, re.M):
            q = os.path.join(FIGDIR, imp + ".py")
            if os.path.exists(q):
                texts.append(read(q))
        letters = set()
        for s in texts:
            letters.update(re.findall(r"set_title\(\s*[\"']([a-h])\s\s",
                                      s))
        out[n] = letters
    return out


def check_figures(S):
    head("2. figures: callout, legend, and every panel the script draws")
    panels = script_panels()
    if not panels:
        fail("fig", "no figure scripts found under figures/")
        return
    legends = dict((int(a), b) for a, b in
                   re.findall(r"\*\*Fig\.\s*(\d+)\s*\|(.*?)(?=\n\n)",
                              S["figs"], re.S))
    for n in sorted(panels):
        calls = len(re.findall(r"Fig\.\s*%d\b" % n, S["body"]))
        if calls == 0:
            fail("fig", "Fig. %d is never referred to in the text" % n)
        if n not in legends:
            fail("fig", "Fig. %d has no legend in the Figures section" % n)
            continue
        want = panels[n]
        have = set(re.findall(r"\*\*([a-h])\*\*,", legends[n]))
        miss = sorted(want - have)
        extra = sorted(have - want)
        if miss:
            fail("fig", "Fig. %d legend does not describe panel(s) %s, which the "
                        "script draws" % (n, ", ".join(miss)))
        if extra:
            fail("fig", "Fig. %d legend describes panel(s) %s that the script "
                        "does not draw" % (n, ", ".join(extra)))
        if calls and not miss and not extra:
            ok("Fig. %d: cited %d time(s), legend covers %d panel(s)"
               % (n, calls, len(want)))


# ------------------------------------------------------------------ 3. references
def cited_refs(txt):
    out = set()
    for m in re.findall(r"\[([\d,\s–-]+)\]", txt):
        for part in m.split(","):
            part = part.strip()
            if re.match(r"^\d+$", part):
                out.add(int(part))
                continue
            r = re.match(r"^(\d+)\s*[–-]\s*(\d+)$", part)
            if r:
                out.update(range(int(r.group(1)), int(r.group(2)) + 1))
    return out


def check_references(S):
    head("3. references: every entry cited, every citation resolvable")
    entries = set(int(x) for x in re.findall(r"^(\d+)\.\s", S["refs"], re.M))
    cited = cited_refs(S["cite"])
    if not entries:
        fail("ref", "no reference entries found")
        return
    dangling = sorted(cited - entries)
    orphan = sorted(entries - cited)
    if dangling:
        fail("ref", "cited with no entry in the list: %s"
             % ", ".join(str(n) for n in dangling))
    else:
        ok("every citation resolves to an entry")
    if orphan:
        fail("ref", "listed and never cited: %s"
             % ", ".join(str(n) for n in orphan))
    else:
        ok("every one of the %d entries is cited" % len(entries))


# ------------------------------------------------------------------ 4. duplicates
def check_duplicates(S):
    head("4. text that appears twice word for word")
    flat = re.sub(r"\s+", " ", re.sub(r"[*`#>]", "", S["cite"]))
    seen, dup = {}, []
    for s in re.split(r"(?<=[.!?]) ", flat):
        s = s.strip()
        if len(s) < 60:
            continue
        k = re.sub(r"[^a-z0-9 ]", "", s.lower())
        if k in seen:
            dup.append(s)
        seen[k] = 1
    for s in dup:
        fail("dup", "repeated verbatim: %s..." % s[:90])
    if not dup:
        ok("no sentence over 60 characters appears twice")


# --------------------------------------------------------------- 5. placeholders
def check_placeholders(S):
    head("5. placeholders still bracketed")
    left = re.findall(r"⟨([^⟩]{1,220})⟩", S["all"])
    for w in left:
        print("  open  %s" % w.split("\n")[0][:70])
    if left:
        print("  (%d outstanding; step152_set_identity.py fills them)" % len(left))
    else:
        ok("none")


def main():
    if not os.path.exists(MS):
        print("manuscript not found: %s" % MS)
        return 2
    S = sections(read(MS))
    check_supplements(S)
    check_figures(S)
    check_references(S)
    check_duplicates(S)
    check_placeholders(S)
    print()
    print("=" * 78)
    if BAD:
        print("cross-reference audit: %d problem(s)" % len(BAD))
        return 1
    print("cross-reference audit: CLEAN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
