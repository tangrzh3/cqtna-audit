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
# The two traits with no producing script. Names follow step94b's convention of
# using the Catalog's own EFO label.
TRAITS = {"lung": ["lung carcinoma"],
          "colorectal": ["colorectal carcinoma"]}
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
for name, terms in TRAITS.items():
    best = {}
    failed = None
    for t in terms:
        try:
            d = fetch(t)
        except Exception as ex:                       # noqa: BLE001
            failed = "%s: %s" % (type(ex).__name__, ex)
            print("[FAIL] %s/%s: %s" % (name, t, failed), flush=True)
            continue
        for a in d.get("_embedded", {}).get("associations", []):
            p = a.get("pvalue")
            if p is None or p >= 5e-8:
                continue
            for l in a.get("loci", []):
                for rs in l.get("strongestRiskAlleles", []):
                    r = rs.get("riskAlleleName", "").split("-")[0]
                    if r.startswith("rs") and (r not in best or p < best[r]):
                        best[r] = p
        print("  %s/%s: running total %d" % (name, t, len(best)), flush=True)

    dest = os.path.join(MR, "known_loci_%s_refetch_%s.csv" % (name, STAMP))
    with io.open(dest, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rsid", "pval"])
        for r, p in sorted(best.items(), key=lambda x: x[1]):
            w.writerow([r, p])

    src = os.path.join(MR, "known_loci_%s.csv" % name)
    old = set()
    if os.path.exists(src):
        with io.open(src, encoding="utf-8") as f:
            old = {r["rsid"] for r in csv.DictReader(f)}
    new = set(best)
    rows.append(dict(
        list=name, trait="; ".join(terms), refetch_date=STAMP,
        n_in_repo=len(old), n_refetched=len(new),
        n_shared=len(old & new), n_only_repo=len(old - new),
        n_only_refetch=len(new - old),
        reproduces=("yes" if old and not (old - new) else "no"),
        fetch_error=failed or ""))
    print("[done] %s: repo %d, refetch %d, shared %d, only-repo %d, only-new %d"
          % (name, len(old), len(new), len(old & new),
             len(old - new), len(new - old)), flush=True)

with io.open(os.path.join(MR, "143a_provenance_diff.tsv"), "w",
             encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("\nwrote 143a_provenance_diff.tsv", flush=True)
print("NOTE: the re-query is written to its own file. Nothing was overwritten;",
      flush=True)
print("      the analysis lists in the repository are unchanged.", flush=True)
