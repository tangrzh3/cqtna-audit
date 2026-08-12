#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 33  Download GSE282266 -- sc-multiome (RNA + ATAC) TCR activation time
         course of naive and memory human CD4+ T cells.

Why this dataset: 0h / 2.5h / 5h / 15h on purified CD4 T cells, both lineages,
with chromatin accessibility. 15h corresponds to the Soskic 16h window in which
the TPI1 cis-eQTL appears and then disappears (Step 22). It therefore allows a
direct, internally controlled test of whether TPI1's regulatory element is
activation-induced, with SPSB2 -- constitutive at the same locus, 3-20 kb away
-- as the built-in negative control.

Run with no arguments to list files and sizes; pass --download to fetch.
Resumable: existing complete files are skipped, partial files are resumed.
"""
import os
import re
import sys
import time
import urllib.request

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282266/suppl/"
DEST = r"D:/Downloads/GSE282266"


def listing():
    idx = urllib.request.urlopen(BASE, timeout=120).read().decode("utf-8", "replace")
    files = [m for m in re.findall(r'href="([^"]+)"', idx)
             if not m.startswith(("/", "http")) and m != "../"]
    out = []
    for f in files:
        try:
            r = urllib.request.urlopen(urllib.request.Request(BASE + f, method="HEAD"),
                                       timeout=60)
            out.append((f, int(r.headers.get("Content-Length", 0))))
        except Exception:
            out.append((f, -1))
    return out


def fetch(name, size):
    os.makedirs(DEST, exist_ok=True)
    dest = os.path.join(DEST, name)
    if os.path.exists(dest) and size > 0 and os.path.getsize(dest) == size:
        print(f"  [skip] {name}", flush=True)
        return
    tmp = dest + ".part"
    have = os.path.getsize(tmp) if os.path.exists(tmp) else 0
    req = urllib.request.Request(BASE + name)
    if have:
        req.add_header("Range", f"bytes={have}-")
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "ab" if have else "wb") as fo:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fo.write(chunk)
    os.replace(tmp, dest)
    mb = os.path.getsize(dest) / 1e6
    print(f"  [done] {name}  {mb:.1f} MB in {time.time()-t0:.0f}s", flush=True)


files = listing()
total = sum(s for _, s in files if s > 0)
for f, s in files:
    print(f"  {f:<58}{s/1e6:9.1f} MB" if s > 0 else f"  {f:<58}      ?")
print(f"\n  {len(files)} files, total {total/1e9:.2f} GB")

# GSE282266_integrated_data.rds.gz is 22 GB compressed -- a Seurat object over
# 16 multiome samples, which would need far more RAM than this machine has once
# loaded. Everything needed (RNA counts, ATAC peak counts, the pooled peak set)
# is in the per-sample CellRanger-ARC matrices, so it is excluded by default.
if "--with-integrated" not in sys.argv:
    before = len(files)
    files = [(f, s) for f, s in files if "integrated_data" not in f]
    skipped = sum(s for f, s in listing() if "integrated_data" in f and s > 0)
    if before != len(files):
        print(f"  excluding integrated_data.rds.gz ({skipped/1e9:.1f} GB); "
              f"pass --with-integrated to include it")
        print(f"  -> {len(files)} files, "
              f"{sum(s for _, s in files if s > 0)/1e9:.2f} GB")

if "--download" in sys.argv:
    print(f"\ndownloading to {DEST}\n")
    for f, s in files:
        try:
            fetch(f, s)
        except Exception as ex:
            print(f"  [FAIL] {f}: {type(ex).__name__}: {ex}", flush=True)
    print("\ndownload pass complete; re-run to retry any failures")
