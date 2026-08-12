#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 41  Fetch only what is needed from GSE166188 (DOGMA-seq).

Purpose: demonstrate finding 7 -- that splitting a cell population by a gene
score also splits it by cell-type purity, and that mRNA-based lineage calls are
the weak link. DOGMA-seq measures 210 surface proteins alongside RNA in the same
nuclei, so CD4/CD8 can be assigned by PROTEIN rather than by dropout-prone
transcripts.

Only the STIM arms are taken (16 h anti-CD3/CD28, matching the Soskic window),
both lysis buffers as technical replicates. The ATAC fragments files are 1.1-1.4
GB each and are not needed; skipping them cuts 6.3 GB to ~0.5 GB.
"""
import os
import time
import urllib.request

DEST = r"D:/Downloads/GSE166188"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5065nnn/{gsm}/suppl/{fn}"

FILES = [
    # LLL stim -- RNA/ATAC combined matrix
    ("GSM5065528", "GSM5065528_LLL_STIM_GExp_ATAC_filtered_barcodes.tsv.gz"),
    ("GSM5065528", "GSM5065528_LLL_STIM_GExp_ATAC_filtered_features.tsv.gz"),
    ("GSM5065528", "GSM5065528_LLL_STIM_GExp_ATAC_filtered_matrix.mtx.gz"),
    # LLL stim -- ADT protein
    ("GSM5065529", "GSM5065529_LLL_stim_ADT_allCounts.barcodes.txt.gz"),
    ("GSM5065529", "GSM5065529_LLL_stim_ADT_allCounts.proteins.txt.gz"),
    ("GSM5065529", "GSM5065529_LLL_stim_ADT_allCounts.mtx.gz"),
    # DIG stim -- RNA/ATAC combined matrix
    ("GSM5065534", "GSM5065534_DIG_STIM_GExp_ATAC_filtered_barcodes.tsv.gz"),
    ("GSM5065534", "GSM5065534_DIG_STIM_GExp_ATAC_filtered_features.tsv.gz"),
    ("GSM5065534", "GSM5065534_DIG_STIM_GExp_ATAC_filtered_matrix.mtx.gz"),
    # DIG stim -- ADT protein
    ("GSM5065535", "GSM5065535_DIG_stim_ADT_allCounts.barcodes.txt.gz"),
    ("GSM5065535", "GSM5065535_DIG_stim_ADT_allCounts.proteins.txt.gz"),
    ("GSM5065535", "GSM5065535_DIG_stim_ADT_allCounts.mtx.gz"),
]

os.makedirs(DEST, exist_ok=True)
for gsm, fn in FILES:
    dest = os.path.join(DEST, fn)
    url = BASE.format(gsm=gsm, fn=fn)
    try:
        head = urllib.request.urlopen(
            urllib.request.Request(url, method="HEAD"), timeout=60)
        size = int(head.headers.get("Content-Length", 0))
    except Exception as ex:
        print(f"  [head fail] {fn}: {type(ex).__name__}: {ex}", flush=True)
        continue
    if os.path.exists(dest) and os.path.getsize(dest) == size:
        print(f"  [skip] {fn}", flush=True)
        continue
    t0 = time.time()
    tmp = dest + ".part"
    try:
        with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fo:
            while True:
                c = r.read(1 << 20)
                if not c:
                    break
                fo.write(c)
        os.replace(tmp, dest)
        print(f"  [done] {fn}  {size/1e6:.1f} MB in {time.time()-t0:.0f}s", flush=True)
    except Exception as ex:
        print(f"  [FAIL] {fn}: {type(ex).__name__}: {ex}", flush=True)

print("\ncontents:")
tot = 0
for f in sorted(os.listdir(DEST)):
    s = os.path.getsize(os.path.join(DEST, f))
    tot += s
    print(f"  {f:<62}{s/1e6:9.1f} MB")
print(f"  total {tot/1e9:.2f} GB")
