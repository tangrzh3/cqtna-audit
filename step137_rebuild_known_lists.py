#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 137 -- rebuild the GRCh38 known-locus lists, and record every failure.

WHY THIS EXISTS
step94b2 resolved rsIDs against Ensembl and, on any batch exception, did
`continue` -- dropping the whole batch with no record. Build rates:

    lung        316 ->   0   (0.0%)
    breast      766 ->   6   (0.8%)
    prostate   1175 ->  35   (3.0%)   <- cleared S39's floor of 30 anyway
    colorectal  833 -> 642  (77.1%)

Nothing reported this, because 36e -- the only consumer of a cross-cancer
known list -- does not use these files at all; it scores every cancer against
the Landi melanoma annotation carried in 36a's `category` column.

Rules are fixed in S41 (manuscript/PREREG_transport_grid.md) §5.1, before any
rebuilt list existed: P < 5e-8, Ensembl GRCh38, chromosomes 1-22 and X, dedupe
by rsID, batch failures retried then queried singly, every unresolved rsID
written out with its reason, and a row is inconclusive below 30 loci or below
a 50% build rate.

Outputs: known_loci_<disease>_grch38.csv   (rebuilt; the old file is kept as
                                            known_loci_<disease>_grch38.pre137.csv)
         137a_list_rebuild_audit.tsv       one row per disease
         137b_unresolved_rsids.tsv         one row per rsID that did not resolve
         137c_console.log
"""
import json
import os
import sys
import time
import urllib.error
import urllib.request

import pandas as pd

MR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CQTNA_DIR") or \
     os.path.dirname(os.path.abspath(__file__))
URL = "https://rest.ensembl.org/variation/homo_sapiens?pops=0"
DISEASES = ["lung", "breast", "prostate", "colorectal", "pancreas"]
P_MAX = 5e-8
KEEP_CHR = [str(i) for i in range(1, 23)] + ["X"]
MIN_LEAD = 30            # S34/S39 frozen floor, reused by S41 §5.1 rule 5
MIN_BUILD_RATE = 0.50    # S41 §5.1 rule 6


def post(ids, tries=4):
    """POST with retries. Raises only after every retry is spent."""
    last = None
    for k in range(tries):
        req = urllib.request.Request(
            URL, data=json.dumps({"ids": ids}).encode(),
            headers={"Content-Type": "application/json",
                     "Accept": "application/json"})
        try:
            return json.loads(urllib.request.urlopen(req, timeout=180).read())
        except Exception as ex:                       # noqa: BLE001
            last = ex
            time.sleep(2 * (k + 1))
    raise last


def place(rec):
    """First GRCh38 mapping on a kept chromosome, or None."""
    for m in rec.get("mappings", []):
        if m.get("assembly_name") == "GRCh38" and \
                str(m.get("seq_region_name")) in KEEP_CHR:
            return str(m["seq_region_name"]), int(m["start"])
    return None


audit, unresolved = [], []
for name in DISEASES:
    src = os.path.join(MR, "known_loci_%s.csv" % name)
    dest = os.path.join(MR, "known_loci_%s_grch38.csv" % name)
    if not os.path.exists(src):
        audit.append(dict(disease=name, n_source=0, n_pass_p=0, n_placed=0,
                          build_rate=float("nan"), verdict="no source file"))
        print("[%s] no source file" % name, flush=True)
        continue

    d = pd.read_csv(src)
    n_src = len(d)
    d = d[pd.to_numeric(d["pval"], errors="coerce") < P_MAX]
    d = d.drop_duplicates("rsid")
    rs = list(d["rsid"])
    if not rs:
        audit.append(dict(disease=name, n_source=n_src, n_pass_p=0, n_placed=0,
                          build_rate=float("nan"), verdict="inconclusive (list)"))
        print("[%s] source empty after P filter" % name, flush=True)
        continue

    rows, missing = [], []
    for i in range(0, len(rs), 190):
        chunk = rs[i:i + 190]
        try:
            res = post(chunk)
        except Exception as ex:                       # noqa: BLE001
            # batch still failing -> query singly so one bad ID cannot sink 190
            print("\n  [%s] batch at %d failed (%s); querying singly"
                  % (name, i, type(ex).__name__), flush=True)
            res = {}
            for r in chunk:
                try:
                    res.update(post([r], tries=2))
                except Exception as ex2:              # noqa: BLE001
                    missing.append((r, "request failed: %s" % type(ex2).__name__))
        for r in chunk:
            rec = res.get(r)
            if rec is None:
                missing.append((r, "not returned by Ensembl"))
                continue
            pl = place(rec)
            if pl is None:
                missing.append((r, "no GRCh38 mapping on chr 1-22,X"))
                continue
            rows.append(dict(rsid=r, chr=pl[0], pos=pl[1]))
        print("  [%s] %d/%d  placed %d  unresolved %d"
              % (name, min(i + 190, len(rs)), len(rs), len(rows), len(missing)),
              end="\r", flush=True)
        time.sleep(0.2)
    print()

    out = pd.DataFrame(rows).drop_duplicates("rsid")
    rate = len(out) / len(rs)
    verdict = "usable"
    if len(out) < MIN_LEAD:
        verdict = "inconclusive (list): below %d loci" % MIN_LEAD
    elif rate < MIN_BUILD_RATE:
        verdict = "inconclusive (list): build rate %.1f%% < %.0f%%" % (
            100 * rate, 100 * MIN_BUILD_RATE)

    if os.path.exists(dest):
        keep = dest.replace("_grch38.csv", "_grch38.pre137.csv")
        if not os.path.exists(keep):
            os.replace(dest, keep)
    out.to_csv(dest, index=False)

    audit.append(dict(disease=name, n_source=n_src, n_pass_p=len(rs),
                      n_placed=len(out), build_rate=round(rate, 4),
                      verdict=verdict))
    unresolved.extend(dict(disease=name, rsid=r, reason=why) for r, why in missing)
    print("[%s] %d of %d placed (%.1f%%) -- %s"
          % (name, len(out), len(rs), 100 * rate, verdict), flush=True)

pd.DataFrame(audit).to_csv(os.path.join(MR, "137a_list_rebuild_audit.tsv"),
                           sep="\t", index=False)
pd.DataFrame(unresolved or [dict(disease="", rsid="", reason="")]).to_csv(
    os.path.join(MR, "137b_unresolved_rsids.tsv"), sep="\t", index=False)
print("\nwrote 137a (%d rows) and 137b (%d rows)"
      % (len(audit), len(unresolved)), flush=True)
