"""Step 94b -- known-locus reference lists for the remaining diseases in the grid.

Same GWAS Catalog endpoint as the HCC list (step 84). Run separately because the
first pass timed out on the three largest traits; nothing about the criteria
changed (PREREG_generality_grid §3: lead SNPs at p < 5e-8, one entry per rsID).
"""
import csv
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
API = ("https://www.ebi.ac.uk/gwas/rest/api/associations/"
       "search/findByEfoTrait?efoTrait={}&size=3000")
TRAITS = {"breast": ["breast carcinoma"],
          "pancreas": ["pancreatic carcinoma"],
          "prostate": ["prostate carcinoma"]}

for name, terms in TRAITS.items():
    dest = os.path.join(MR, f"known_loci_{name}.csv")
    if os.path.exists(dest) and os.path.getsize(dest) > 200:
        print(f"[skip] {name}", flush=True)
        continue
    best = {}
    for t in terms:
        url = API.format(urllib.parse.quote(t))
        try:
            raw = urllib.request.urlopen(url, timeout=600).read()
        except Exception as ex:
            print(f"[FAIL] {name}/{t}: {type(ex).__name__}: {ex}", flush=True)
            continue
        d = json.loads(raw)
        for a in d.get("_embedded", {}).get("associations", []):
            p = a.get("pvalue")
            if p is None or p >= 5e-8:
                continue
            for l in a.get("loci", []):
                for rs in l.get("strongestRiskAlleles", []):
                    r = rs.get("riskAlleleName", "").split("-")[0]
                    if r.startswith("rs") and (r not in best or p < best[r]):
                        best[r] = p
        print(f"  {name}/{t}: running total {len(best):,}", flush=True)
        time.sleep(0.3)
    with io.open(dest, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rsid", "pval"])
        for r, p in sorted(best.items(), key=lambda x: x[1]):
            w.writerow([r, p])
    print(f"[done] {name}: {len(best):,} lead rsIDs", flush=True)
print("all done", flush=True)
