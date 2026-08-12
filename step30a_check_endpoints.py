"""Check which FinnGen R12 cancer endpoints are reachable before committing to
multi-GB downloads. Prints size and HTTP status for each candidate."""
import urllib.request

BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{}.gz"
CAND = [
    "C3_BRONCHUS_LUNG_EXALLC", "C3_BRONCHUS_LUNG",
    "C3_COLORECTAL_EXALLC", "C3_COLORECTAL",
    "C3_PANCREAS_EXALLC", "C3_PANCREAS",
    "C3_BREAST_EXALLC", "C3_BREAST",
    "C3_PROSTATE_EXALLC", "C3_PROSTATE",
    "C3_MELANOMA_SKIN_EXALLC",          # positive control -- already used locally
]
for e in CAND:
    url = BASE.format(e)
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            mb = int(r.headers.get("Content-Length", 0)) / 1e6
            print(f"  OK   {e:<28} {mb:8.1f} MB")
    except Exception as ex:
        print(f"  FAIL {e:<28} {type(ex).__name__}: {ex}")
