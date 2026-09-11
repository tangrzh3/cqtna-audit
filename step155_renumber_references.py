#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 155 -- put the reference list in citation order, and drop what nothing cites.

WHY
Two separate defects, one fix. The list was not in citation order at all -- the
first-appearance sequence ran 1, 42, 44, 35, 45, ... -- which the target journal
requires, and it carried entries no sentence in the manuscript cites, left behind
when compression moved their analyses into the supplement. Renumbering by hand
across forty-seven entries and thirty-odd citation groups is precisely the task
where one silently wrong number survives every proofread.

WHAT IT DOES
Reads the manuscript, expands every citation group, orders the entries by first
appearance, drops the ones never cited, rewrites both the in-text markers and
the list, and prints the full old -> new map plus what it removed.

It refuses to write if the result would not round-trip: every cited number must
have an entry, the map must be a bijection, and re-running must be a no-op.

  python step155_renumber_references.py [repo] [--apply]

Without --apply it reports and changes nothing. Run step154_crossref_audit.py
afterwards.
"""
import io
import os
import re
import sys

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
APPLY = "--apply" in sys.argv
MR = (ARGS[0] if ARGS
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(MR, "manuscript", "MANUSCRIPT_GB.md")
DASH = "–"
GROUP = re.compile(r"\[([\d,\s–-]+)\]")


def asc(x):
    """This console is not always UTF-8; author names are."""
    return x.encode("ascii", "replace").decode("ascii")


def expand(g):
    out = []
    for part in g.split(","):
        part = part.strip()
        if re.match(r"^\d+$", part):
            out.append(int(part))
            continue
        m = re.match(r"^(\d+)\s*[%s-]\s*(\d+)$" % DASH, part)
        if m:
            out.extend(range(int(m.group(1)), int(m.group(2)) + 1))
        else:
            raise ValueError("unparsable citation group: [%s]" % g)
    return out


def render(nums):
    """Collapse runs of three or more into a range, as the source style does."""
    nums = sorted(set(nums))
    parts, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j - i >= 2:
            parts.append("%d%s%d" % (nums[i], DASH, nums[j]))
        else:
            parts.extend(str(n) for n in nums[i:j + 1])
        i = j + 1
    return "[" + ",".join(parts) + "]"


def supplement_brackets():
    """Renumbering the main list silently rewrites what every bracketed number
    in the SUPPLEMENTS points at, because nothing links the two. That already
    happened once: S47 cited the pseudotime test as [8] under the old numbering,
    and [8] became a different paper. Supplements should cite by author name;
    this reports anything that still does not. Numeric intervals like
    "8.3 [7,11]" are indistinguishable from a citation and are listed too --
    saying so is better than guessing."""
    import glob
    out = []
    for pat in ("manuscript/SUPP_*.md", "manuscript/PREREG_*.md", "SUPP_*.md"):
        for f in sorted(glob.glob(os.path.join(MR, pat))):
            txt = io.open(f, encoding="utf-8", errors="replace").read()
            fenced = False
            for i, line in enumerate(txt.split(chr(10)), 1):
                # Code is not citation. A fenced block or an inline span can
                # hold sys.argv[1] or x[0,1], and reporting those buries the
                # real cases: on 2026-09-11, eight of eighteen hits were
                # argv[1]/argv[2]/argv[3] inside a paragraph about argument
                # parsing. A warning nobody can scan is a warning nobody reads.
                if line.lstrip().startswith("```"):
                    fenced = not fenced
                    continue
                if fenced or line.startswith("    "):
                    continue
                bare = re.sub(r"`[^`]*`", "", line)
                for m in re.finditer(r"\[\d{1,2}(?:,\d{1,2})*\]", bare):
                    out.append((os.path.relpath(f, MR), i, m.group(0)))
    return out


def main():
    t = io.open(MS, encoding="utf-8").read()
    cut = t.index("## References")
    body, refs = t[:cut], t[cut:]

    entries = {}
    for m in re.finditer(r"^(\d+)\.\s(.*)$", refs, re.M):
        entries[int(m.group(1))] = m.group(2)
    if not entries:
        print("no reference entries found")
        return 2

    order = []
    for m in GROUP.finditer(body):
        for n in expand(m.group(1)):
            if n not in order:
                order.append(n)

    missing = [n for n in order if n not in entries]
    if missing:
        print("refusing: cited with no entry -> %s"
              % ", ".join(str(n) for n in missing))
        return 1

    new = dict((old, i + 1) for i, old in enumerate(order))
    dropped = sorted(set(entries) - set(order))

    print("=" * 74)
    print("old -> new  (%d entries kept, %d dropped)" % (len(order), len(dropped)))
    print("=" * 74)
    for old in sorted(new):
        print("  %3d -> %3d   %s" % (old, new[old], asc(entries[old][:72])))
    if dropped:
        print()
        print("dropped, cited nowhere in the manuscript:")
        for old in dropped:
            print("  %3d          %s" % (old, asc(entries[old][:72])))
        print()
        print("  These belong to analyses the compression moved out of the main")
        print("  text. If a supplement still uses one, it carries its own list.")

    new_body = GROUP.sub(lambda m: render(new[n] for n in expand(m.group(1))), body)
    lines = ["## References", ""]
    for old in sorted(new, key=lambda k: new[k]):
        lines.append("%d. %s" % (new[old], entries[old]))
    new_refs = "\n".join(lines) + "\n"

    # round-trip: the rewritten text must renumber to itself
    check_order = []
    for m in GROUP.finditer(new_body):
        for n in expand(m.group(1)):
            if n not in check_order:
                check_order.append(n)
    if check_order != list(range(1, len(order) + 1)):
        print("\nrefusing: the rewritten text does not number 1..N in order")
        return 1

    sb = supplement_brackets()
    if sb:
        print("")
        print("[!] bracketed numerals in supplements. Renumbering the main list")
        print("    does not touch them, so any that are citations now point")
        print("    elsewhere. Cite by author name in supplements instead.")
        for f, i2, w in sb[:20]:
            print("      %-46s line %-5d %s" % (f, i2, w))
        if len(sb) > 20:
            print("      ... and %d more" % (len(sb) - 20))

    if not APPLY:
        # ⚠ The exit code has to mean something, because step159 reads
        # it as one of S54's four acceptance audits. Until 2026-09-11 this
        # returned 0 unconditionally in report mode, so "step155 reference
        # order ok" meant "the script ran" -- an out-of-order list would
        # have printed a whole non-identity map above and still exited
        # clean. Same shape as the testthat skip in S54 section 9.11.7: a
        # check whose green says less than it reads.
        #   Report mode now exits non-zero when work is outstanding, the
        # convention a --check flag follows elsewhere. Nothing is written
        # either way; --apply is still what writes.
        moved = sorted(o for o in new if new[o] != o)
        if moved or dropped:
            print("\nreference list is NOT in citation order: "
                  "%d entry(ies) would move, %d would be dropped."
                  % (len(moved), len(dropped)))
            print("re-run with --apply to write the corrected list")
            return 1
        print("\nreference list is already in citation order "
              "(%d entries, none moved, none dropped); nothing to write"
              % len(new))
        return 0
    io.open(MS, "w", encoding="utf-8", newline="\n").write(new_body + new_refs)
    print("\nwrote %s" % MS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
