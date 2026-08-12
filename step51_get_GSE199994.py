#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 51  Download GSE199994 (10x multiome PBMC, 8 melanoma patients at baseline
         + 2 healthy donors).

Used for two analyses, in this order:
  C  variance decomposition -- is the CD4 glycolysis programme a person-level
     trait in circulating blood at all? This is a prerequisite for any biomarker
     claim and needs donors, not outcome, so all 10 are downloaded.
  A  the pre-specified IRC test (see manuscript/GSE199994_prespecified_design.md,
     written before any of this data was fetched).

~20 GB. Resumable: complete files are skipped, partial ones resume.
"""
import os
import time
import urllib.request

DEST = r"D:/Downloads/GSE199994"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE199nnn/GSE199994/suppl/GSE199994_{s}.tar.gz"
SAMPLES = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "HD1", "HD2"]

def head_size(url, tries=6):
    """NCBI drops connections after sustained transfer; back off and retry."""
    for a in range(tries):
        try:
            h = urllib.request.urlopen(
                urllib.request.Request(url, method="HEAD"), timeout=90)
            return int(h.headers.get("Content-Length", 0))
        except Exception as ex:
            if a == tries - 1:
                print(f"    head failed after {tries} tries: "
                      f"{type(ex).__name__}", flush=True)
                return None
            time.sleep(10 * (a + 1))
    return None


os.makedirs(DEST, exist_ok=True)
total = 0
for s in SAMPLES:
    url = BASE.format(s=s)
    dest = os.path.join(DEST, f"GSE199994_{s}.tar.gz")
    size = head_size(url)
    if size is None:
        print(f"  [head fail] {s}", flush=True)
        continue
    if os.path.exists(dest) and os.path.getsize(dest) == size:
        print(f"  [skip] {s}  {size/1e9:.2f} GB", flush=True)
        total += size
        continue
    tmp = dest + ".part"
    t0 = time.time()
    # resume across attempts: a dropped connection leaves a partial file that
    # the next attempt continues with a Range request
    for attempt in range(8):
        have = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        if have >= size > 0:
            break
        req = urllib.request.Request(url)
        if have:
            req.add_header("Range", f"bytes={have}-")
        try:
            with urllib.request.urlopen(req, timeout=600) as r, \
                 open(tmp, "ab" if have else "wb") as fo:
                while True:
                    c = r.read(1 << 22)
                    if not c:
                        break
                    fo.write(c)
        except Exception as ex:
            got = os.path.getsize(tmp) if os.path.exists(tmp) else 0
            print(f"    {s} attempt {attempt+1}: {type(ex).__name__} at "
                  f"{got/1e9:.2f}/{size/1e9:.2f} GB, retrying", flush=True)
            time.sleep(15 * (attempt + 1))
    have = os.path.getsize(tmp) if os.path.exists(tmp) else 0
    if have >= size > 0:
        os.replace(tmp, dest)
        total += have
        print(f"  [done] {s}  {have/1e9:.2f} GB in {time.time()-t0:.0f}s",
              flush=True)
    else:
        print(f"  [FAIL] {s}: stopped at {have/1e9:.2f}/{size/1e9:.2f} GB",
              flush=True)
    time.sleep(5)          # brief pause between samples

print(f"\ntotal on disk: {total/1e9:.2f} GB")
print("re-run to retry anything that failed")
