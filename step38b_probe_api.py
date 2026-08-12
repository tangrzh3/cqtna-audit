"""Probe the eQTL Catalogue v2 API to find the correct query form and confirm
whether results are paginated by position (making 'top p' meaningless)."""
import json
import time
import urllib.request
import urllib.parse

API = "https://www.ebi.ac.uk/eqtl/api/v2"
UA = {"User-Agent": "Mozilla/5.0"}
DS = "QTD000479"                     # CD4_T-cell_naive
GENE = "ENSG00000111669"             # TPI1
RSID = "rs12302749"
POS = 6867132


def try_call(label, path, **params):
    url = f"{API}/{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90) as r:
                d = json.load(r)
            n = len(d) if isinstance(d, list) else 1
            print(f"  OK   {label:<42} n={n}")
            if isinstance(d, list) and d:
                k = d[0]
                print(f"       keys: {sorted(k.keys())}")
                print(f"       first: { {kk: k[kk] for kk in list(k)[:9]} }")
            return d
        except Exception as ex:
            code = getattr(ex, "code", None)
            if attempt == 4:
                print(f"  FAIL {label:<42} {type(ex).__name__} {code}")
                return None
            time.sleep(3 * (attempt + 1))


print("=== query forms ===")
try_call("by variant rsid", f"datasets/{DS}/associations", variant_id=RSID, size=100)
time.sleep(1)
try_call("by rsid param", f"datasets/{DS}/associations", rsid=RSID, size=100)
time.sleep(1)
try_call("by molecular_trait_id", f"datasets/{DS}/associations",
         molecular_trait_id=GENE, size=100)
time.sleep(1)
try_call("by pos window", f"datasets/{DS}/associations",
         pos=f"12:{POS-5000}-{POS+5000}", size=100)
time.sleep(1)
d = try_call("gene_id + nlog10p sort", f"datasets/{DS}/associations",
             gene_id=GENE, size=100)

print("\n=== is the gene_id result ordered by position? ===")
if d:
    pos = [r.get("position") for r in d[:10]]
    print(f"  first 10 positions: {pos}")
    ps = [r.get("nlog10p") for r in d[:10]]
    print(f"  first 10 nlog10p  : {ps}")
    allpos = [r.get("position") for r in d if r.get("position")]
    if allpos:
        print(f"  position range in page: {min(allpos):,}-{max(allpos):,}")
        print(f"  monotonic? {allpos == sorted(allpos)}")
