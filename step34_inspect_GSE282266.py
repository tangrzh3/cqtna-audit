#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 34  Inspect GSE282266 before any analysis.

Needs to establish, from the files themselves rather than assumption:
  - do the features files contain BOTH genes and ATAC peaks (CellRanger-ARC)?
  - what does "set_N" mean -- donor, or a pool of donors?
  - are the peak coordinates on GRCh38 (our instruments are GRCh38)?
  - is there a peak overlapping the TPI1 instrument (chr12:6,867,132) and the
    SPSB2 instruments (6,847,298 / 6,869,846 / 6,886,522)?
"""
import gzip
import os
import collections

D = r"D:/Downloads/GSE282266"
TP = [("TPI1", 12, 6867132), ("SPSB2", 12, 6847298),
      ("SPSB2", 12, 6869846), ("SPSB2", 12, 6886522)]

print("=== files ===")
for f in sorted(os.listdir(D)):
    print(f"  {f:<52}{os.path.getsize(os.path.join(D, f))/1e6:9.1f} MB")

# ---------------------------------------------------------------- features
print("\n=== features (act_15_set_1) ===")
fp = os.path.join(D, "GSE282266_act_15_set_1_features.tsv.gz")
kinds = collections.Counter()
examples = {}
n = 0
with gzip.open(fp, "rt") as fh:
    for line in fh:
        n += 1
        parts = line.rstrip("\n").split("\t")
        k = parts[2] if len(parts) > 2 else "?"
        kinds[k] += 1
        if k not in examples:
            examples[k] = parts
print(f"  total features: {n:,}")
for k, v in kinds.most_common():
    print(f"    {k:<24}{v:>10,}   e.g. {examples[k]}")

# ---------------------------------------------------------------- barcodes
print("\n=== barcodes per sample ===")
for tp in ["rest", "act_2.5", "act_5", "act_15"]:
    row = []
    for s in range(1, 5):
        p = os.path.join(D, f"GSE282266_{tp}_set_{s}_barcodes.tsv.gz")
        if not os.path.exists(p):
            row.append("--"); continue
        with gzip.open(p, "rt") as fh:
            bcs = [l.strip() for l in fh]
        row.append(f"{len(bcs):,}")
    print(f"  {tp:<10}" + "".join(f"{x:>12}" for x in row))
with gzip.open(os.path.join(D, "GSE282266_rest_set_1_barcodes.tsv.gz"), "rt") as fh:
    head = [next(fh).strip() for _ in range(3)]
print(f"  barcode format e.g. {head}")

# ---------------------------------------------------------------- peaks
print("\n=== pooled peaks ===")
pk = os.path.join(D, "GSE282266_pooled_peaks.bed.gz")
peaks12 = []
npk = 0
with gzip.open(pk, "rt") as fh:
    for line in fh:
        if line.startswith(("#", "track")):
            continue
        f = line.rstrip("\n").split("\t")
        npk += 1
        if npk <= 3:
            print(f"  e.g. {f[:4]}")
        if f[0] in ("chr12", "12"):
            peaks12.append((int(f[1]), int(f[2])))
print(f"  total peaks: {npk:,} | on chr12: {len(peaks12):,}")

print("\n=== peaks overlapping the instruments (+/- 5 kb) ===")
for sym, ch, pos in TP:
    hit = [(s, e) for s, e in peaks12 if s - 5000 <= pos <= e + 5000]
    direct = [(s, e) for s, e in peaks12 if s <= pos <= e]
    print(f"  {sym:<7} chr{ch}:{pos:<10} peaks within 5 kb: {len(hit):<3} "
          f"| directly containing the SNP: {len(direct)}")
    for s, e in hit[:4]:
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        print(f"        peak {s:,}-{e:,}  width {e-s:,}  distance {d:,}")

# ---------------------------------------------------------------- RAW.tar
raw = os.path.join(D, "GSE282266_RAW.tar")
if os.path.exists(raw):
    import tarfile
    print("\n=== RAW.tar contents (first 20) ===")
    with tarfile.open(raw) as t:
        for i, m in enumerate(t.getmembers()):
            if i >= 20:
                print("  ...")
                break
            print(f"  {m.name:<58}{m.size/1e6:8.2f} MB")
