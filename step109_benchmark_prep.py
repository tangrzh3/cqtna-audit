#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 109 -- S34 preparation: fix the nomination file before any truth set exists.

S34 §4.1 requires process blinding: the pipeline's per-locus nomination (file A)
and the external truth set (file B) must be written separately and joined only by
the scoring script, so that the locus definition, the significance threshold and
the inclusion rules cannot be adjusted after the truth is visible.

This script writes file A and nothing else about the external benchmark. It runs
before any candidate truth set has been opened, which is the strongest available
ordering: A is fixed first, and its SHA-256 is recorded so that any later edit is
detectable.

It also recomputes the EXISTING self-assembled ten-locus check under the strict
and lenient hit definitions that S34 §4 fixes. That is not part of the blinded
comparison -- the ten-locus list is our own and already published -- but S34 §4.2
requires the two checks to be reported side by side, which is only meaningful if
both are scored the same way.

Inputs:  96a_attribution_selfcheck.tsv
Outputs: 109a_fileA_nominations.tsv   (locus + named genes, NO truth columns)
         109b_selfcheck_rescored.tsv  (strict vs lenient on the existing ten)
         109_console.log
"""
import hashlib
import io
import os
import sys

import pandas as pd

MR = r"D:/R_ex/MR"
SRC = f"{MR}/96a_attribution_selfcheck.tsv"
FILE_A = f"{MR}/109a_fileA_nominations.tsv"
RESCORED = f"{MR}/109b_selfcheck_rescored.tsv"


class Tee:
    def __init__(self, path):
        self.f = io.open(path, "w", encoding="utf-8")

    def write(self, s):
        sys.__stdout__.write(s)
        self.f.write(s)

    def flush(self):
        sys.__stdout__.flush()
        self.f.flush()


sys.stdout = Tee(f"{MR}/109_console.log")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def split_genes(cell):
    """Gene lists are ';'-joined; a truncated cell ends in '...'."""
    if pd.isna(cell) or not str(cell).strip():
        return [], False
    raw = str(cell).strip()
    truncated = raw.endswith("...")
    if truncated:
        raw = raw[:-3]
    return [g for g in raw.split(";") if g], truncated


df = pd.read_csv(SRC, sep="\t", dtype=str).fillna("")

print("=" * 78)
print("Step 109 -- S34 file A, written before any truth set is opened")
print("=" * 78)
print(f"source: 96a_attribution_selfcheck.tsv  sha256 {sha256(SRC)[:16]}...")
print(f"loci with a significant nomination under the eQTLGen exposure: {len(df)}")

# ---- file A: nominations only, no truth -----------------------------------
a = df[["locus", "n_named", "named"]].copy()
a["n_named"] = a["n_named"].astype(int)
a["named_truncated"] = [split_genes(x)[1] for x in a["named"]]
a = a.sort_values("locus").reset_index(drop=True)
a.to_csv(FILE_A, sep="\t", index=False)

n_trunc = int(a["named_truncated"].sum())
print(f"\nfile A written: {FILE_A}")
print(f"  loci: {len(a)}   sha256: {sha256(FILE_A)[:16]}...")
print(f"  columns: {list(a.columns)}  <- contains no accepted-gene information")
if n_trunc:
    print(f"  !! {n_trunc} loci have a TRUNCATED gene list in the source table")
    print("     (cell ends in '...'), so file A cannot be scored strictly for")
    print("     them until the full nomination is regenerated. Affected loci:")
    for L in a.loc[a["named_truncated"], "locus"]:
        print(f"       {L}")

print("\nnomination set size distribution (all loci):")
print(f"  median {a['n_named'].median():.1f}   "
      f"IQR {a['n_named'].quantile(.25):.1f}-{a['n_named'].quantile(.75):.1f}   "
      f"max {a['n_named'].max()}")
print(f"  single-gene nominations: {(a['n_named'] == 1).sum()} of {len(a)}")

# ---- rescore the existing ten under both definitions ----------------------
print("\n" + "=" * 78)
print("existing self-assembled ten loci, rescored under S34 §4 definitions")
print("=" * 78)

rows = []
for _, r in df.iterrows():
    acc, _ = split_genes(r["accepted"])
    if not acc:
        continue
    named, trunc = split_genes(r["named"])
    lenient = any(g in named for g in acc)
    strict = (int(r["n_named"]) == 1) and len(named) == 1 and named[0] in acc
    rows.append(dict(locus=r["locus"], accepted=r["accepted"],
                     n_named=int(r["n_named"]), named=r["named"],
                     named_truncated=trunc,
                     lenient_hit=lenient, strict_hit=strict,
                     published_hit=(r["hit"] == "True")))

s = pd.DataFrame(rows).sort_values("locus").reset_index(drop=True)
s.to_csv(RESCORED, sep="\t", index=False)

n = len(s)
nl, ns = int(s["lenient_hit"].sum()), int(s["strict_hit"].sum())
print(f"evaluable loci: {n}")
print(f"  lenient (accepted gene is among the named): {nl}/{n}")
print(f"  strict  (nomination is exactly that gene):  {ns}/{n}")
agree = int((s["lenient_hit"] == s["published_hit"]).sum())
print(f"\nreproduces the published 'hit' column under the LENIENT rule: "
      f"{agree}/{n}")
print("-> the published 6/10 is the lenient rate; S34 fixes strict as primary,")
print("   so the two checks must be reported under the same definition.")

print("\nper-locus:")
print(f"  {'locus':<8}{'accepted':<16}{'n':>3}  {'lenient':<8}{'strict':<7}named")
for _, r in s.iterrows():
    print(f"  {r['locus']:<8}{r['accepted'][:15]:<16}{r['n_named']:>3}  "
          f"{str(r['lenient_hit']):<8}{str(r['strict_hit']):<7}{r['named'][:44]}")

print(f"\nwrote {RESCORED}")
print("\nNOT done here, and deliberately: no external truth set has been opened,")
print("scored against S34 §3.1, or joined to file A. That is step 110, and it")
print("requires the §9.1 admissibility table to be filled in first.")
