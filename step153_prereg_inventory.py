#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 153 -- how many pre-registration documents are there, exactly?

WHY
The manuscript states the number five times and disagrees with itself: two
places say twelve (one of them enumerating twelve S-numbers) and three say
fourteen. There are nineteen PREREG_*.md files on disk. A count asserted five
times and enumerated once is a count that will drift again, so this reads it off
the files instead of trusting any of the five.

WHAT COUNTS AS ONE
Not the filename prefix. Fourteen of the nineteen files belong to one declared
series: each says in its own header that it is of the "same specification" as
the earlier members, naming them, and the later ones carry an ordinal (the
seventh, the eighth, ... the fourteenth). The founding document cannot make that
claim -- there was nothing yet to be the same specification as -- so membership
is also conferred by being named in someone else's claim. Counting only the
self-declarations silently drops S9, the one document every other member points
at.

The remaining five say in their own headers that they are something else: a
frozen decision record, a post-hoc exploratory analysis, a partly-frozen
estimand whose tolerance is an analysis choice, an adoption record that states
it is "not yet a pre-registration", and a transport test whose header forbids
calling it pre-registered. Those disclaimers are the point of the documents, so
the count follows them rather than the glob.

  python step153_prereg_inventory.py [repo]

Exits non-zero if the manuscript's stated count or enumerated list disagrees
with the files.
"""
import io
import os
import re
import sys

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))
MSDIR = os.path.join(MR, "manuscript")

WORDS = {"Ten": 10, "Eleven": 11, "Twelve": 12, "Thirteen": 13, "Fourteen": 14,
         "Fifteen": 15, "Sixteen": 16, "Seventeen": 17, "Eighteen": 18,
         "Nineteen": 19, "Twenty": 20}

SAME_SPEC = "同规格"
OUTSIDE = [
    ("决策记录", "frozen decision record"),
    ("事后（post-hoc）探索性", "post-hoc exploratory"),
    ("事后探索性", "post-hoc exploratory"),
    ("尚未构成", "adopted, not yet a pre-registration"),
    ("部分冻结", "partly frozen; the tolerance is an analysis choice"),
    ("不得称 pre-registered", "explicitly not to be called pre-registered"),
]


def read(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


def cjk_int(s):
    """Headers number themselves in Chinese numerals; the console is not always
    UTF-8, so render them as integers."""
    d = dict(zip("一二三四五六七八九", range(1, 10)))
    if s == "十":
        return 10
    if s.startswith("十"):
        return 10 + d.get(s[1:], 0)
    if s.endswith("十"):
        return d.get(s[:-1], 0) * 10
    if "十" in s:
        a, b = s.split("十")
        return d.get(a, 1) * 10 + d.get(b, 0)
    return d.get(s, 0)


def assembled_map():
    """S-number -> file for S9-S33, from the supplementary table in the older
    assembled build. Six of the early documents never wrote an S-number into
    their own header, and that table is where the number actually lives."""
    p = os.path.join(MSDIR, "MANUSCRIPT_assembled.md")
    out = {}
    if not os.path.exists(p):
        return out
    for m in re.finditer(r"^\|\s*S(\d+)\s*\|.*?`(PREREG_[a-z0-9_]+\.md)`",
                         read(p), re.M):
        out[m.group(2)] = int(m.group(1))
    return out


def main():
    amap = assembled_map()
    files = sorted(f for f in os.listdir(MSDIR)
                   if f.startswith("PREREG_") and f.endswith(".md"))
    rows = []
    for f in files:
        head = read(os.path.join(MSDIR, f))[:1400]
        m = re.search(r"\*\*编号\*\*：S(\d+)", head)
        num = int(m.group(1)) if m else amap.get(f)
        series = SAME_SPEC in head
        mo = re.search(r"第([一二三四五六七八九十]+)份", head)
        ordinal = cjk_int(mo.group(1)) if mo else None
        why = ""
        if not series:
            for pat, lab in OUTSIDE:
                if pat in head:
                    why = lab
                    break
            why = why or "no series claim in header"
        rows.append([f, num, series, ordinal, why, head])

    # membership by being named in another member's same-specification claim
    named = set()
    for r in rows:
        if r[2]:
            i = r[5].index(SAME_SPEC)
            named.update(int(x) for x in re.findall(r"S(\d+)", r[5][max(0, i - 240):i]))
            # early members are named by filename, not by number
            for other in rows:
                if other[0] in r[5][max(0, i - 240):i + 200] and other[1]:
                    named.add(other[1])
    for r in rows:
        if not r[2] and r[1] in named:
            r[2], r[4] = True, "the founding document, named by every other member"

    print("=" * 92)
    print("pre-registration family: %d files in manuscript/" % len(rows))
    print("=" * 92)
    print("  %-42s %-5s %-7s %s" % ("file", "S", "series", "status in its own header"))
    for f, num, series, ordinal, why, _ in sorted(
            rows, key=lambda r: (not r[2], r[1] or 999)):
        note = why if why else ("no. %d in the series" % ordinal if ordinal else "member")
        print("  %-42s %-5s %-7s %s"
              % (f, ("S%d" % num) if num else "?", "yes" if series else "no", note))

    ser = sorted(r[1] for r in rows if r[2] and r[1])
    print()
    print("registered series: %d documents -> %s"
          % (len(ser), ", ".join("S%d" % n for n in ser)))
    print("outside the series: %d" % sum(1 for r in rows if not r[2]))

    ms = os.path.join(MSDIR, "MANUSCRIPT_GB.md")
    txt = read(ms)
    bad = 0
    print()
    print("=" * 92)
    print("what MANUSCRIPT_GB.md says")
    print("=" * 92)
    pat = r"(%s)\s+pre-registrations?(?:\s+documents)?" % "|".join(
        list(WORDS) + [w.lower() for w in WORDS])
    said = re.findall(pat, txt)
    for w in said:
        ok = WORDS[w.capitalize()] == len(ser)
        bad += 0 if ok else 1
        print("  %s  \"%s pre-registration(s)\"  (files say %d)"
              % ("OK " if ok else "WRONG", w, len(ser)))
    if not said:
        bad += 1
        print("  WRONG  no count sentence found at all")

    m = re.search(r"pre-registration documents \(Supplementary ([^)]*)\)", txt, re.S)
    if m:
        listed = sorted(int(x) for x in re.findall(r"S(\d+)", m.group(1)))
        missing = [n for n in ser if n not in listed]
        extra = [n for n in listed if n not in ser]
        ok = not missing and not extra
        bad += 0 if ok else 1
        print("  %s  enumerated list has %d entries"
              % ("OK " if ok else "WRONG", len(listed)))
        if missing:
            print("         missing: %s" % ", ".join("S%d" % n for n in missing))
        if extra:
            print("         not in the series: %s" % ", ".join("S%d" % n for n in extra))
    else:
        bad += 1
        print("  WRONG  no enumerated list found")

    print()
    print("inventory:", "CLEAN" if bad == 0
          else "%d disagreement(s) between the files and the manuscript" % bad)
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
