#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 156 -- C6: which estimator produced each paper's primary cis list.

WHY
step151 tried to answer this by phrase matching and could not. Its count fell
from 71 papers to 10 when the single-variant phrase was required to occur in the
same sentence as a cis-eQTL mention, and 23 of the 71 never mention a cis-eQTL at
all. That is a measure of wording, not of design, and the paper has withdrawn it
as a bound in either direction.

C6 replaces it with a judgement two people make independently on the SAME
fixed-seed 46-paper subsample already used for C1 and C3, so the estimate is
comparable with the audit's other criteria rather than drawn from a fresh sample.
Rules are in 156_C6_CODING_RULES.txt, written before any paper was read.

USAGE
  python step156_estimator_provenance_coding.py [repo] prepare
      writes 156a_C6_coderA.tsv and 156b_C6_coderB.tsv -- blind, differently
      shuffled, carrying candidate evidence passages and no automated label.

  python step156_estimator_provenance_coding.py [repo] score
      reads both files back once their `coder` columns are filled and writes
      156c_C6_agreement.tsv, 156d_C6_disagreements.tsv and 156e_console.log.

The passages are an aid to finding the answer, not the answer. A coder who
cannot resolve a paper from them is expected to open the full text in
litaudit/<pmc>.txt, which is included in both files as `fulltext`.

Outputs: 156a, 156b (prepare); 156c, 156d, 156e (score)
"""
import io
import os
import random
import re
import sys

ARGS = [a for a in sys.argv[1:] if not a.startswith("-")]
MR = (ARGS[0] if ARGS and os.path.isdir(ARGS[0])
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))
MODE = ([a for a in ARGS if a in ("prepare", "score")] or ["prepare"])[0]

SUBSAMPLE = os.path.join(MR, "105a_blind_coding.tsv")
RULES = os.path.join(MR, "156_C6_CODING_RULES.txt")
CATS = ("W", "M", "X", "U")

# Sentences worth putting in front of a coder: they mention an instrument count,
# an estimator, or a cis analysis. Deliberately generous -- the cost of an extra
# passage is a coder's second of reading, the cost of a missing one is a wrong U.
WANTED = re.compile(
    r"(instrument|IVW|inverse[- ]variance|weighted median|MR[- ]Egger|"
    r"weighted mode|wald|cis[- ]eQTL|cis[- ]QTL|SNP per|per gene|lead SNP|"
    r"top SNP|index SNP|clump)", re.I)
CIS = re.compile(r"cis[- ]eQTL|cis[- ]QTL", re.I)


def read(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


def fulltext_path(pmc):
    for d in ("litaudit", os.path.join("105a_research_cache", "fulltext")):
        q = os.path.join(MR, d, pmc + ".txt")
        if os.path.exists(q):
            return q
    return None


def passages(pmc, limit=14):
    """Candidate evidence sentences, cis-mentioning ones first."""
    q = fulltext_path(pmc)
    if not q:
        return "(full text not cached)"
    t = re.sub(r"\s+", " ", read(q))
    sents = [s.strip() for s in re.split(r"(?<=[.!?])\s+", t)
             if 40 < len(s.strip()) < 400 and WANTED.search(s)]
    with_cis = [s for s in sents if CIS.search(s)]
    without = [s for s in sents if not CIS.search(s)]
    picked = with_cis[:limit] + without[:max(0, limit - len(with_cis[:limit]))]
    if not picked:
        return "(none found)"
    return " || ".join(picked)


def subsample():
    rows, seen = [], set()
    for line in read(SUBSAMPLE).split("\n")[1:]:
        if not line.strip():
            continue
        f = line.split("\t")
        pmid, pmc, journal, year = f[1], f[2], f[3], f[4]
        if pmc in seen:
            continue
        seen.add(pmc)
        rows.append(dict(pmid=pmid, pmc=pmc, journal=journal, year=year))
    return rows


def prepare():
    rows = subsample()
    print("papers in the fixed-seed subsample: %d" % len(rows))
    if not os.path.exists(RULES):
        print("refusing: %s is missing. The rules are written before the coding,"
              % RULES)
        print("not after it.")
        return 1
    for r in rows:
        r["evidence"] = passages(r["pmc"])
        r["fulltext"] = (os.path.relpath(fulltext_path(r["pmc"]), MR)
                         if fulltext_path(r["pmc"]) else "")
    missing = sum(1 for r in rows if r["evidence"] == "(full text not cached)")
    if missing:
        print("[!] %d paper(s) have no cached full text; they are still listed, "
              "and are a U unless the coder can obtain the text" % missing)

    cols = ["row", "pmid", "pmc", "journal", "year", "evidence", "fulltext",
            "coder", "evidence_quote", "evidence_where", "note"]
    for tag, seed, out in (("A", 1561, "156a_C6_coderA.tsv"),
                           ("B", 1562, "156b_C6_coderB.tsv")):
        shuffled = list(rows)
        random.Random(seed).shuffle(shuffled)
        with io.open(os.path.join(MR, out), "w", encoding="utf-8",
                     newline="\n") as f:
            f.write("# C6 -- coder %s. Rules: 156_C6_CODING_RULES.txt\n" % tag)
            f.write("# Enter W, M, X or U in `coder`. Quote your evidence.\n")
            f.write("# This file carries no automated label and no other "
                    "coder's judgement, by design.\n")
            f.write("\t".join(cols) + "\n")
            for i, r in enumerate(shuffled, 1):
                f.write("\t".join([str(i), r["pmid"], r["pmc"], r["journal"],
                                   r["year"], r["evidence"], r["fulltext"],
                                   "", "", "", ""]) + "\n")
        print("wrote %s (%d rows, seed %d)" % (out, len(shuffled), seed))
    print("\nBoth files are blind and differently shuffled. Fill `coder` in each,")
    print("independently, then run: python %s %s score"
          % (os.path.basename(__file__), MR))
    return 0


def load_coded(name):
    p = os.path.join(MR, name)
    if not os.path.exists(p):
        return None
    out = {}
    hdr = None
    for line in read(p).split("\n"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if hdr is None:
            hdr = f
            continue
        d = dict(zip(hdr, f))
        v = (d.get("coder") or "").strip().upper()
        out[d["pmc"]] = dict(code=v, quote=(d.get("evidence_quote") or "").strip(),
                             where=(d.get("evidence_where") or "").strip())
    return out


def kappa(a, b, cats):
    n = len(a)
    if not n:
        return float("nan")
    po = sum(1 for x, y in zip(a, b) if x == y) / float(n)
    pe = 0.0
    for c in cats:
        pe += (a.count(c) / float(n)) * (b.count(c) / float(n))
    return (po - pe) / (1 - pe) if pe < 1 else float("nan")


def wilson(k, n, z=1.96):
    if not n:
        return (float("nan"), float("nan"))
    p = k / float(n)
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * ((p * (1 - p) / n + z * z / (4.0 * n * n)) ** 0.5)
    return ((c - h) / d, (c + h) / d)


def score():
    A, B = load_coded("156a_C6_coderA.tsv"), load_coded("156b_C6_coderB.tsv")
    if A is None or B is None:
        print("refusing: run `prepare` first, then fill both coder files.")
        return 1
    both = [p for p in A if p in B and A[p]["code"] in CATS and B[p]["code"] in CATS]
    if not both:
        print("no rows are coded in both files yet -- nothing to score.")
        print("coder A filled: %d, coder B filled: %d"
              % (sum(1 for p in A if A[p]["code"] in CATS),
                 sum(1 for p in B if B[p]["code"] in CATS)))
        return 1

    ca = [A[p]["code"] for p in both]
    cb = [B[p]["code"] for p in both]
    agree = sum(1 for x, y in zip(ca, cb) if x == y)
    k = kappa(ca, cb, CATS)

    lines = []
    lines.append("C6 -- estimator provenance of the primary cis nomination list")
    lines.append("papers coded by both: %d of %d" % (len(both), len(A)))
    lines.append("raw agreement       : %d/%d = %.1f%%"
                 % (agree, len(both), 100.0 * agree / len(both)))
    lines.append("Cohen's kappa       : %.3f  (pre-adjudication, as registered)"
                 % k)
    lines.append("")
    lines.append("distribution, coder A / coder B")
    for c in CATS:
        lines.append("  %s : %3d / %3d" % (c, ca.count(c), cb.count(c)))
    # A conservative reading until adjudication: agreed rows only.
    agreed = [x for x, y in zip(ca, cb) if x == y]
    asc = [c for c in agreed if c != "U"]
    if asc:
        w = agreed.count("W")
        lo, hi = wilson(w, len(asc))
        lines.append("")
        lines.append("agreed rows only, before adjudication:")
        lines.append("  ascertainable (W+M+X): %d;  W = %d" % (len(asc), w))
        lines.append("  W share of ascertainable: %.0f%% [%.0f, %.0f]"
                     % (100.0 * w / len(asc), 100 * lo, 100 * hi))
        lines.append("  U (not ascertainable) reported separately: %d"
                     % agreed.count("U"))
    lines.append("")
    lines.append("[!] This is not the reportable estimate. Adjudicate the")
    lines.append("    disagreements in 156d jointly, then recompute. Kappa above")
    lines.append("    stays as it is -- it is defined pre-adjudication.")

    out = "\n".join(lines)
    print(out)
    io.open(os.path.join(MR, "156e_console.log"), "w", encoding="utf-8",
            newline="\n").write(out + "\n")

    with io.open(os.path.join(MR, "156c_C6_agreement.tsv"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write("pmc\tcoderA\tcoderB\tagree\tA_where\tB_where\n")
        for p in sorted(both):
            f.write("%s\t%s\t%s\t%d\t%s\t%s\n"
                    % (p, A[p]["code"], B[p]["code"],
                       int(A[p]["code"] == B[p]["code"]),
                       A[p]["where"], B[p]["where"]))
    with io.open(os.path.join(MR, "156d_C6_disagreements.tsv"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write("pmc\tcoderA\tA_quote\tcoderB\tB_quote\tadjudicated\treason\n")
        for p in sorted(both):
            if A[p]["code"] != B[p]["code"]:
                f.write("%s\t%s\t%s\t%s\t%s\t\t\n"
                        % (p, A[p]["code"], A[p]["quote"],
                           B[p]["code"], B[p]["quote"]))
    print("\nwrote 156c / 156d / 156e")
    return 0


if __name__ == "__main__":
    raise SystemExit(prepare() if MODE == "prepare" else score())
