"""Check what DICE (Schmiedel 2018) data is publicly downloadable.

Wanted: cis-eQTL summary statistics for CD4+ T naive (unstimulated) and
CD4+ T stimulated (anti-CD3/CD28 4h), ideally also stimulated CD8 for a
cell-type-specificity contrast.
"""
import urllib.request

CANDIDATES = [
    "https://dice-database.org/downloads",
    "https://dice-database.org/",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118165/matrix/GSE118165_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118165/suppl/",
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/",
    "https://www.ebi.ac.uk/eqtl/api/v2/datasets/?study=QTS000010",
]

for u in CANDIDATES:
    try:
        req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read(4000)
            ct = r.headers.get("Content-Type", "")
            print(f"OK   {u}\n     status={r.status} type={ct} len={r.headers.get('Content-Length')}")
            txt = body.decode("utf-8", "replace")
            snippet = " ".join(txt.split())[:400]
            print(f"     {snippet}\n")
    except Exception as ex:
        print(f"FAIL {u}\n     {type(ex).__name__}: {ex}\n")
