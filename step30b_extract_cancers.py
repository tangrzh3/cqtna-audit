#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 30b  Stream FinnGen R12 cancer GWAS and keep only the instrument SNPs.

~815 MB per endpoint x 5 = ~4 GB. Streaming and filtering on the fly avoids
storing any of it: only the ~2,100 instrument rows per cancer are written.
Each extract is small, so re-running the MR later is instant.

Resumable: an endpoint whose extract already exists is skipped.
"""
import csv
import gzip
import io
import os
import sys
import time
import urllib.request

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
OUT = os.path.join(MR, "cancer_extracts")
os.makedirs(OUT, exist_ok=True)
BASE = ("https://storage.googleapis.com/finngen-public-data-r12/"
        "summary_stats/release/finngen_R12_{}.gz")

CANCERS = ["C3_BRONCHUS_LUNG_EXALLC", "C3_COLORECTAL_EXALLC",
           "C3_PANCREAS_EXALLC", "C3_BREAST_EXALLC", "C3_PROSTATE_EXALLC"]

# ---------------------------------------------------------------- instruments
need = set()
with open(os.path.join(MR, "13_meta_locus_annotation.tsv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        need.add(r["SNP"])                      # "chr:pos"
print(f"instrument SNPs to extract: {len(need):,}", flush=True)

# ---------------------------------------------------------------- stream
for endpoint in CANCERS:
    dest = os.path.join(OUT, f"{endpoint}.tsv")
    if os.path.exists(dest):
        print(f"[skip] {endpoint} already extracted", flush=True)
        continue
    url = BASE.format(endpoint)
    t0 = time.time()
    print(f"[start] {endpoint}", flush=True)
    got, n = [], 0
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            with gzip.GzipFile(fileobj=resp) as gz:
                txt = io.TextIOWrapper(gz, encoding="utf-8")
                header = txt.readline().rstrip("\n").split("\t")
                ic, ip = header.index("#chrom"), header.index("pos")
                ir, ia = header.index("ref"), header.index("alt")
                ib = header.index("beta")
                ise = header.index("sebeta")
                ipv = header.index("pval")
                for line in txt:
                    n += 1
                    f = line.rstrip("\n").split("\t")
                    k = f[ic] + ":" + f[ip]
                    if k in need:
                        got.append((k, f[ir], f[ia], f[ib], f[ise], f[ipv]))
                    if n % 5_000_000 == 0:
                        print(f"    {endpoint}: {n:,} lines, {len(got):,} kept, "
                              f"{time.time()-t0:.0f}s", flush=True)
    except Exception as ex:
        print(f"[FAIL] {endpoint}: {type(ex).__name__}: {ex}", flush=True)
        continue

    tmp = dest + ".part"
    with open(tmp, "w", newline="", encoding="utf-8") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["SNP", "ref", "alt", "beta", "sebeta", "pval"])
        w.writerows(got)
    os.replace(tmp, dest)
    print(f"[done ] {endpoint}: {n:,} lines scanned, {len(got):,} SNPs kept, "
          f"{time.time()-t0:.0f}s -> {dest}", flush=True)

print("\nextraction complete")
for endpoint in CANCERS:
    p = os.path.join(OUT, f"{endpoint}.tsv")
    print(f"  {endpoint:<28} {'OK' if os.path.exists(p) else 'MISSING'}")
