"""Step 94a -- stream the five FinnGen R12 cancer endpoints, keeping only the rows
at eQTLGen instrument positions.

Executes row 2 of the grid in PREREG_generality_grid.md. Mirrors step30b (which did
the same for Soskic instruments) so the two columns of the grid are extracted by
identical machinery; only the instrument list differs.

Each endpoint is ~815 MB gzipped and is never stored: it is streamed, filtered on
the fly, and only the ~12k instrument rows are written. Resumable -- an endpoint
whose extract already exists is skipped.
"""
import csv
import gzip
import io
import os
import sys
import time
import urllib.request

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "cancer_extracts_eqtlgen")
os.makedirs(OUT, exist_ok=True)
BASE = ("https://storage.googleapis.com/finngen-public-data-r12/"
        "summary_stats/release/finngen_R12_{}.gz")
CANCERS = ["C3_BRONCHUS_LUNG_EXALLC", "C3_COLORECTAL_EXALLC",
           "C3_PANCREAS_EXALLC", "C3_BREAST_EXALLC", "C3_PROSTATE_EXALLC"]

need = set()
with io.open(os.path.join(MR, "92c_locus_annotated.tsv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        need.add((r["chr"], r["pos"]))
print(f"eQTLGen instrument positions to extract: {len(need):,}", flush=True)

for endpoint in CANCERS:
    dest = os.path.join(OUT, f"{endpoint}.tsv")
    if os.path.exists(dest):
        print(f"[skip] {endpoint}", flush=True)
        continue
    url = BASE.format(endpoint)
    t0 = time.time()
    print(f"[start] {endpoint}", flush=True)
    try:
        got, n = [], 0
        with urllib.request.urlopen(url, timeout=300) as resp:
            with gzip.open(resp, "rt") as gz:
                hdr = gz.readline().rstrip("\n").split("\t")
                ix = {c: i for i, c in enumerate(hdr)}
                for line in gz:
                    n += 1
                    f = line.rstrip("\n").split("\t")
                    key = (f[ix["#chrom"]], f[ix["pos"]])
                    if key in need:
                        got.append(f)
                    if n % 5_000_000 == 0:
                        print(f"    {endpoint}: {n:,} lines, {len(got):,} kept, "
                              f"{time.time()-t0:.0f}s", flush=True)
        with io.open(dest, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(hdr)
            w.writerows(got)
        print(f"[done] {endpoint}: {len(got):,} rows kept of {n:,} "
              f"({time.time()-t0:.0f}s)", flush=True)
    except Exception as ex:
        print(f"[FAIL] {endpoint}: {type(ex).__name__}: {ex}", flush=True)
print("all done", flush=True)
