#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 156 -- C6: which estimator produced each paper's primary cis nomination list.

WHY
step151 tried to answer this by phrase matching and failed: its count fell from
71 papers to 10 when the single-variant phrase was required to appear in the same
sentence as a cis-eQTL mention, and 23 of the 71 never mention a cis-eQTL at all.
That number is a phrase-hit tally of unknown direction, not a bound, and the
manuscript now says so. C6 replaces it with a judgement two people make
independently on the SAME 46-paper subsample already used for C1 and C3, so the
result is comparable with the criteria beside it rather than a fresh sample with
its own selection.

  python step156_estimator_provenance_coding.py [repo] prepare   # write the two files
  python step156_estimator_provenance_coding.py [repo] score     # once both are filled

PREPARE writes 156a_C6_coderA.tsv and 156b_C6_coderB.tsv. Neither carries the
step151 scan result, neither carries the other coder's file, and the two are
shuffled under different seeds so that any ordering in the corpus is not shared.
Each row offers candidate evidence passages; the rules require the coder to open
the full text when those do not settle it.

SCORE refuses to compute anything until both files are complete. It reports raw
agreement and Cohen's kappa BEFORE adjudication, writes the disagreements to
156d for joint adjudication, and -- once those are adjudicated -- the
distribution over W/M/X/U with a Wilson interval for the W share among the
ascertainable papers. U is reported separately and never redistributed.

Rules: 156_C6_CODING_RULES.txt, fixed before any paper was read for this
criterion.

Outputs: 156a_C6_coderA.tsv, 156b_C6_coderB.tsv, 156c_C6_agreement.tsv,
         156d_C6_disagreements.tsv, 156e_C6_result.tsv
"""
import io
import math
import os
import random
import re
import sys

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
MR = (ARGS[0] if ARGS and os.path.isdir(ARGS[0])
      else os.environ.get("CQTNA_DIR") or os.path.dirname(os.path.abspath(__file__)))
MODE = ([a for a in ARGS if a in ("prepare", "score")] or ["prepare"])[0]

SUBSAMPLE = os.path.join(MR, "105a_blind_coding.tsv")
FULLTEXT_DIRS = [os.path.join(MR, "litaudit"),
                 os.path.join(MR, "105a_research_cache", "fulltext")]
CODES = ["W", "M", "X", "U"]
SEED_A, SEED_B = 1561, 1562
MAX_PASSAGES = 14

# Sentences worth putting in front of a coder. Ordered: the ones that can
# actually settle the question first.
CIS = re.compile(r"cis[- ]eQTL|cis[- ]QTL", re.I)
INSTRUMENT = re.compile(
    r"instrument|\bSNP[s]?\b|variant[s]?\b|IV\b|genetic proxy|proxies", re.I)
ESTIMATOR = re.compile(
    r"wald|inverse[- ]variance|\bIVW\b|weighted median|MR[- ]Egger|weighted mode"
    r"|mode[- ]based", re.I)
COUNTY = re.compile(r"\bone\b|\bsingle\b|\bonly\b|\btop\b|\blead\b|\bnumber of\b"
                    r"|\d+\s+(?:SNP|variant|instrument)", re.I)


def read(p):
    return io.open(p, encoding="utf-8", errors="replace").read()


def subsample():
    """The 46 papers already double-coded for C1 and C3, in file order."""
    rows = [l.split("\t") for l in read(SUBSAMPLE).rstrip("\n").split("\n")]
    head = rows[0]
    ix = dict((c, i) for i, c in enumerate(head))
    seen, out = set(), []
    for r in rows[1:]:
        pmc = r[ix["pmc"]]
        if pmc in seen:
            continue
        seen.add(pmc)
        out.append(dict(pmc=pmc, pmid=r[ix["pmid"]], journal=r[ix["journal"]],
                        year=r[ix["year"]]))
    return out


def fulltext_path(pmc):
    for d in FULLTEXT_DIRS:
        p = os.path.join(d, pmc + ".txt")
        if os.path.exists(p):
            return p
    return ""


def passages(path):
    """Candidate evidence, best first. This is a reading aid and nothing else:
    the rules require the full text to be opened when these do not settle it,
    and a judgement resting on a passage that is not quoted back is treated as
    U at adjudication."""
    if not path:
        return "(no cached full text -- open the paper)"
    txt = re.sub(r"\s+", " ", read(path))
    sents = re.split(r"(?<=[.!?]) ", txt)
    scored = []
    for s in sents:
        if len(s) < 40 or len(s) > 600:
            continue
        sc = 0
        if CIS.search(s):
            sc += 4
        if ESTIMATOR.search(s):
            sc += 3
        if INSTRUMENT.search(s):
            sc += 2
        if COUNTY.search(s):
            sc += 1
        # Threshold deliberately low. Requiring "cis-eQTL" in the same sentence
        # as the instrument description is exactly the mistake that broke
        # step151: most papers describe their instruments in sentences that
        # never repeat the word. A cis mention still ranks a sentence first.
        if sc >= 3:
            scored.append((sc, s.strip()))
    scored.sort(key=lambda x: -x[0])
    keep = [s for _, s in scored[:MAX_PASSAGES]]
    if not keep:
        return "(nothing matched -- open the paper)"
    return "  ||  ".join(keep)


def write_coder_file(rows, seed, dest, who):
    r = list(rows)
    random.Random(seed).shuffle(r)
    cols = ["row", "pmc", "pmid", "journal", "year", "fulltext", "evidence",
            "coder", "evidence_quote", "evidence_where"]
    with io.open(dest, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(cols) + "\n")
        for i, x in enumerate(r, 1):
            f.write("\t".join([
                str(i), x["pmc"], x["pmid"], x["journal"], x["year"],
                x["fulltext"], x["evidence"], "", "", ""]) + "\n")
    print("wrote %s  (coder %s, %d papers, shuffle seed %d)"
          % (os.path.basename(dest), who, len(r), seed))


def prepare():
    rows = subsample()
    for x in rows:
        p = fulltext_path(x["pmc"])
        x["fulltext"] = os.path.relpath(p, MR) if p else ""
        x["evidence"] = passages(p)
    missing = sum(1 for x in rows if not x["fulltext"])
    write_coder_file(rows, SEED_A, os.path.join(MR, "156a_C6_coderA.tsv"), "A")
    write_coder_file(rows, SEED_B, os.path.join(MR, "156b_C6_coderB.tsv"), "B")
    print()
    print("%d of %d papers have no cached full text and must be opened at the"
          % (missing, len(rows)))
    print("publisher. They are not dropped: an unreadable paper is coded U only")
    print("after someone has tried to read it.")
    print()
    print("Read 156_C6_CODING_RULES.txt first. Fill `coder` with W, M, X or U,")
    print("paste the sentence you relied on into `evidence_quote`, and where it")
    print("came from into `evidence_where`. Do not open the other coder's file,")
    print("and do not open 151a_estimator_usage.tsv.")


def load_coded(path, who):
    if not os.path.exists(path):
        print("missing: %s" % os.path.basename(path))
        return None
    rows = [l.split("\t") for l in read(path).rstrip("\n").split("\n")]
    ix = dict((c, i) for i, c in enumerate(rows[0]))
    out, blank, bad = {}, [], []
    for r in rows[1:]:
        code = (r[ix["coder"]] or "").strip().upper()
        pmc = r[ix["pmc"]]
        if not code:
            blank.append(pmc)
            continue
        if code not in CODES:
            bad.append((pmc, code))
            continue
        out[pmc] = dict(code=code,
                        quote=(r[ix["evidence_quote"]] or "").strip(),
                        where=(r[ix["evidence_where"]] or "").strip())
    if blank:
        print("coder %s: %d row(s) not yet coded (%s%s)"
              % (who, len(blank), ", ".join(blank[:5]),
                 " ..." if len(blank) > 5 else ""))
    for pmc, code in bad:
        print("coder %s: %s carries %r, which is not one of W/M/X/U"
              % (who, pmc, code))
    return None if (blank or bad) else out


def kappa(a, b, keys):
    n = len(keys)
    obs = sum(1 for k in keys if a[k]["code"] == b[k]["code"]) / float(n)
    exp = 0.0
    for c in CODES:
        pa = sum(1 for k in keys if a[k]["code"] == c) / float(n)
        pb = sum(1 for k in keys if b[k]["code"] == c) / float(n)
        exp += pa * pb
    k = (obs - exp) / (1 - exp) if exp < 1 else float("nan")
    return obs, exp, k


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / float(n)
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4.0 * n * n))
    return ((c - h) / d, (c + h) / d)


def score():
    A = load_coded(os.path.join(MR, "156a_C6_coderA.tsv"), "A")
    B = load_coded(os.path.join(MR, "156b_C6_coderB.tsv"), "B")
    if A is None or B is None:
        print("\nnothing computed. Both files must be complete first --"
              " a partial kappa is worse than none.")
        return 1
    keys = sorted(set(A) & set(B))
    obs, exp, k = kappa(A, B, keys)
    print("=" * 74)
    print("C6, before adjudication")
    print("=" * 74)
    print("  papers coded by both      : %d" % len(keys))
    print("  raw agreement             : %.1f%%" % (100 * obs))
    print("  expected agreement        : %.3f" % exp)
    print("  Cohen's kappa             : %.3f" % k)
    print()
    print("  Reported beside the raw agreement on purpose. With four categories")
    print("  and one likely dominant, kappa can collapse for prevalence reasons")
    print("  while agreement is high -- the artefact already documented for C1.")

    dis = [x for x in keys if A[x]["code"] != B[x]["code"]]
    dpath = os.path.join(MR, "156d_C6_disagreements.tsv")
    prior = {}
    if os.path.exists(dpath):
        rows = [l.split("\t") for l in read(dpath).rstrip("\n").split("\n")]
        ix = dict((c, i) for i, c in enumerate(rows[0]))
        for r in rows[1:]:
            if len(r) > ix["adjudicated"]:
                v = (r[ix["adjudicated"]] or "").strip().upper()
                if v in CODES:
                    prior[r[ix["pmc"]]] = v
    with io.open(dpath, "w", encoding="utf-8", newline="\n") as f:
        f.write("pmc\tcoderA\tA_quote\tA_where\tcoderB\tB_quote\tB_where"
                "\tadjudicated\treason\n")
        for x in dis:
            f.write("\t".join([x, A[x]["code"], A[x]["quote"], A[x]["where"],
                               B[x]["code"], B[x]["quote"], B[x]["where"],
                               prior.get(x, ""), ""]) + "\n")
    print()
    print("  %d disagreement(s) written to 156d_C6_disagreements.tsv"
          % len(dis))

    # ---- sensitivity over the unresolved cells. NOT in the registered
    # reporting list; added deliberately BEFORE adjudication rather than after.
    # Putting both extremes on the record first is what stops an adjudication
    # from being steered by knowing which way the headline needs to move. The
    # coders are not shown this until they have finished.
    agreed = dict((x, A[x]["code"]) for x in keys if A[x]["code"] == B[x]["code"])
    base = dict((c, sum(1 for v in agreed.values() if v == c)) for c in CODES)
    print()
    print("=" * 74)
    print("sensitivity over the unresolved %d (unregistered; see the note in code)"
        % len(dis))
    print("=" * 74)
    print("  agreed %d: %s" % (len(agreed),
                             ", ".join("%s %d" % (c, base[c]) for c in CODES)))
    print("  %-26s %5s %5s %5s %5s  %8s  %s"
        % ("scenario", "W", "M", "X", "U", "W share", "95% CI"))
    scen = [("all unresolved -> %s" % c, {c: len(dis)}) for c in CODES]
    scen.append(("each to coder A",
                 dict((c, sum(1 for x in dis if A[x]["code"] == c)) for c in CODES)))
    scen.append(("each to coder B",
                 dict((c, sum(1 for x in dis if B[x]["code"] == c)) for c in CODES)))
    brows = []
    for name, add in scen:
        cc = dict((c, base[c] + add.get(c, 0)) for c in CODES)
        asc2 = cc["W"] + cc["M"] + cc["X"]
        lo2, hi2 = wilson(cc["W"], asc2)
        print("  %-26s %5d %5d %5d %5d  %7.1f%%  [%.1f, %.1f]"
            % (name, cc["W"], cc["M"], cc["X"], cc["U"],
               100.0 * cc["W"] / asc2 if asc2 else float("nan"),
               100 * lo2, 100 * hi2))
        brows.append(dict(scenario=name, **cc))
    with io.open(os.path.join(MR, "156f_C6_bounds.tsv"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write("scenario\tW\tM\tX\tU\n")
        for r in brows:
            f.write("%s\t%d\t%d\t%d\t%d\n"
                    % (r["scenario"], r["W"], r["M"], r["X"], r["U"]))
    print("  wrote 156f_C6_bounds.tsv")

    unresolved = [x for x in dis if x not in prior]
    if unresolved:
        print("  %d still unadjudicated. Resolve them JOINTLY -- never by one"
              % len(unresolved))
        print("  coder alone -- fill `adjudicated` and `reason`, then re-run.")
        print("  kappa above is final either way; it is defined pre-adjudication.")
        return 0

    final = dict((x, A[x]["code"]) for x in keys)
    final.update(prior)
    counts = dict((c, sum(1 for v in final.values() if v == c)) for c in CODES)
    asc = counts["W"] + counts["M"] + counts["X"]
    lo, hi = wilson(counts["W"], asc)
    print()
    print("=" * 74)
    print("C6, adjudicated")
    print("=" * 74)
    for c, lbl in [("W", "single-variant Wald"), ("M", "multi-instrument"),
                   ("X", "mixed"), ("U", "not ascertainable")]:
        print("  %s  %-22s %3d" % (c, lbl, counts[c]))
    print()
    print("  ascertainable (W+M+X)     : %d" % asc)
    if asc:
        print("  W share of ascertainable  : %.1f%%  [%.1f, %.1f] (Wilson)"
              % (100.0 * counts["W"] / asc, 100 * lo, 100 * hi))
    print("  U is reported as itself and NOT redistributed: a paper that does")
    print("  not say is a finding about the literature's reporting.")

    with io.open(os.path.join(MR, "156e_C6_result.tsv"), "w",
                 encoding="utf-8", newline="\n") as f:
        f.write("quantity\tvalue\n")
        f.write("n_coded\t%d\n" % len(keys))
        f.write("raw_agreement\t%.4f\n" % obs)
        f.write("kappa_preadjudication\t%.4f\n" % k)
        for c in CODES:
            f.write("n_%s\t%d\n" % (c, counts[c]))
        f.write("n_ascertainable\t%d\n" % asc)
        if asc:
            f.write("W_share\t%.4f\nW_lo\t%.4f\nW_hi\t%.4f\n"
                    % (counts["W"] / float(asc), lo, hi))
    print("\nwrote 156e_C6_result.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(prepare() if MODE == "prepare" else score())
