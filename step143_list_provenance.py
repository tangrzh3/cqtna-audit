#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 143 -- provenance for the lung and colorectal known-locus lists.

WHY THIS EXISTS
step94b produced the breast, pancreas and prostate rsID lists and its docstring
says the first pass "timed out on the three largest traits" -- meaning lung and
colorectal were produced by an earlier run of the same code that was not kept.
So two of the six reference lists in this study have no producing script in the
repository. For a paper whose claim is auditability that is not acceptable, and
the fix is not to assert the provenance but to test it.

This re-queries the GWAS Catalog with step94b's criteria and traits, writes to a
SEPARATE path, and diffs against the file already in the repository. It does not
overwrite anything. A list that no longer reproduces is a finding, not an error
to be papered over -- the Catalog grows, so some drift is expected and the size
and direction of the drift is what needs recording.

Criteria, unchanged from step94b / PREREG_generality_grid section 3:
  lead rsIDs at p < 5e-8, one entry per rsID keeping the smallest p.

Outputs: known_loci_<name>_refetch_<date>.csv   the re-query
         143a_provenance_diff.tsv               per list: sizes and set diffs
         143b_console.log
"""
import csv
import datetime
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request

MR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CQTNA_DIR") or \
     os.path.dirname(os.path.abspath(__file__))
API = ("https://www.ebi.ac.uk/gwas/rest/api/associations/"
       "search/findByEfoTrait?efoTrait={}&size=3000")
# The two lists with no producing script. step94b's TRAITS values are LISTS of
# Catalog labels, so the original run may have unioned more than one; which ones
# is exactly what is not recorded. Rather than guess, try each candidate label
# and every union of them, and report which combination reproduces the file.
CANDIDATES = {"lung": ["lung carcinoma", "lung adenocarcinoma",
                       "squamous cell lung carcinoma",
                       "small cell lung carcinoma",
                       "non-small cell lung carcinoma"],
              "colorectal": ["colorectal carcinoma", "colorectal cancer"]}
STAMP = datetime.date.today().isoformat()


def fetch(term, tries=4):
    url = API.format(urllib.parse.quote(term))
    last = None
    for k in range(tries):
        try:
            return json.loads(urllib.request.urlopen(url, timeout=600).read())
        except Exception as ex:                       # noqa: BLE001
            last = ex
            time.sleep(3 * (k + 1))
    raise last


rows = []
for name, terms in CANDIDATES.items():
    src = os.path.join(MR, "known_loci_%s.csv" % name)
    old = set()
    if os.path.exists(src):
        with io.open(src, encoding="utf-8") as f:
            old = {r["rsid"] for r in csv.DictReader(f)}

    per_term, errs = {}, {}
    for t in terms:
        try:
            d = fetch(t)
        except Exception as ex:                       # noqa: BLE001
            errs[t] = "%s: %s" % (type(ex).__name__, ex)
            print("[FAIL] %s / %s: %s" % (name, t, errs[t]), flush=True)
            continue
        best = {}
        for a in d.get("_embedded", {}).get("associations", []):
            p = a.get("pvalue")
            if p is None or p >= 5e-8:
                continue
            for l in a.get("loci", []):
                for rs in l.get("strongestRiskAlleles", []):
                    r = rs.get("riskAlleleName", "").split("-")[0]
                    if r.startswith("rs") and (r not in best or p < best[r]):
                        best[r] = p
        per_term[t] = best
        print("  %s / %s: %d lead rsIDs" % (name, t, len(best)), flush=True)

    # every non-empty combination of the terms that fetched
    got = [t for t in terms if t in per_term]
    combos = []
    for mask in range(1, 1 << len(got)):
        sel = [got[i] for i in range(len(got)) if mask >> i & 1]
        u = {}
        for t in sel:
            for r, pv in per_term[t].items():
                if r not in u or pv < u[r]:
                    u[r] = pv
        combos.append((sel, u))

    best_share = max((len(set(c[1]) & old) for c in combos), default=0)
    for sel, u in combos:
        new = set(u)
        interesting = (len(sel) == 1 or new == old or
                       len(new & old) == best_share)
        rows.append(dict(
            list=name, traits=" | ".join(sel), refetch_date=STAMP,
            n_in_repo=len(old), n_refetched=len(new),
            n_shared=len(old & new), n_only_repo=len(old - new),
            n_only_refetch=len(new - old),
            # Four outcomes, not two. The one that matters is old subset of new:
            # every rsID in the repository re-fetches, and the query returns
            # some the Catalog has added since. That is a reproduction plus
            # growth, and my first version labelled it "no".
            reproduces=("exact" if old == new
                        else "reproduced + catalog growth"
                        if old and not (old - new) and (new - old)
                        else "incomplete: repo has entries not re-fetched"
                        if old and not (new - old)
                        else "no"),
            fetch_error="; ".join("%s -> %s" % kv for kv in errs.items())))
        if not interesting:
            continue
        print("  %-12s %-44s repo %4d  refetch %4d  shared %4d  only-repo %3d  only-new %3d  %s"
              % (name, " | ".join(sel)[:44], len(old), len(new), len(old & new),
                 len(old - new), len(new - old), rows[-1]["reproduces"]),
              flush=True)

    # write the best combination (most shared with the repo) to its own file
    if combos:
        sel, u = max(combos, key=lambda c: len(set(c[1]) & old))
        dest = os.path.join(MR, "known_loci_%s_refetch_%s.csv" % (name, STAMP))
        with io.open(dest, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["rsid", "pval"])
            for r, pv in sorted(u.items(), key=lambda x: x[1]):
                w.writerow([r, pv])

with io.open(os.path.join(MR, "143a_provenance_diff.tsv"), "w",
             encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("\nwrote 143a_provenance_diff.tsv", flush=True)
print("NOTE: the re-query is written to its own file. Nothing was overwritten;",
      flush=True)
print("      the analysis lists in the repository are unchanged.", flush=True)
