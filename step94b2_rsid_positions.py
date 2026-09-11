"""Step 94b2 -- GRCh38 positions for the known-locus rsIDs, from Ensembl.

WHY THIS EXISTS (a bug, recorded rather than quietly fixed)
step94c mapped known-locus rsIDs to positions by looking them up in the FinnGen
*extract* -- but that extract had already been filtered down to eQTLGen instrument
positions. So a known lead SNP was only found a position if it happened to BE an
instrument. The symptom was obvious once compared across diseases: the known-locus
background came out as 2/559 loci for lung and 5/559 for colorectal, against
58/554 for melanoma and 31/549 for HCC, whose lists carried their own coordinates.
With a background that sparse the cell has essentially no power and the resulting
"0/8 known, fold = 0" for lung was uninterpretable, not negative.

Fixed here by resolving rsIDs against Ensembl GRCh38 directly, independent of any
instrument list. The lung and colorectal numbers produced before this fix are void
and are not reported.
"""
import sys
import io
import json
import os
import time
import urllib.request

import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
URL = "https://rest.ensembl.org/variation/homo_sapiens?pops=0"
DISEASES = ["lung", "colorectal", "pancreas", "breast", "prostate"]


def post(ids):
    req = urllib.request.Request(
        URL, data=json.dumps({"ids": ids}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    return json.loads(urllib.request.urlopen(req, timeout=180).read())


for name in DISEASES:
    src = f"{MR}/known_loci_{name}.csv"
    dest = f"{MR}/known_loci_{name}_grch38.csv"
    if not os.path.exists(src) or os.path.getsize(src) < 200:
        print(f"[skip] {name}: no rsID list yet", flush=True)
        continue
    if os.path.exists(dest):
        print(f"[skip] {name}: already resolved", flush=True)
        continue
    rs = list(pd.read_csv(src).rsid)
    rows = []
    for i in range(0, len(rs), 190):
        chunk = rs[i:i + 190]
        try:
            d = post(chunk)
        except Exception as ex:
            print(f"  {name} batch {i}: {type(ex).__name__}", flush=True)
            time.sleep(2)
            continue
        for r, v in d.items():
            for m in v.get("mappings", []):
                if m.get("assembly_name") == "GRCh38" and str(m.get("seq_region_name")) in \
                        [str(x) for x in range(1, 23)] + ["X", "Y"]:
                    rows.append(dict(rsid=r, chr=str(m["seq_region_name"]),
                                     pos=int(m["start"])))
                    break
        print(f"  {name}: {min(i+190, len(rs))}/{len(rs)} resolved {len(rows)}",
              end="\r", flush=True)
        time.sleep(0.2)
    print()
    pd.DataFrame(rows).drop_duplicates("rsid").to_csv(dest, index=False)
    print(f"[done] {name}: {len(rows):,} of {len(rs):,} rsIDs placed on GRCh38",
          flush=True)
print("all done", flush=True)
