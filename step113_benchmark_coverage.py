#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 113 -- S34 §9.1a amendment 2, step 2: coverage without reading gene names.

Decides whether the external benchmark can run at all, by intersecting the Open
Targets gold standards with this pipeline's significant loci. Per S34 the gene
identity column must not be read at this stage, so `gold_standard_info.gene_id`
is dropped immediately after download and an assertion enforces its absence for
the rest of the script. Only coordinates, trait names and evidence classes are
used.

Registered decisions this script applies (all fixed before it was run):
  * evidence classes kept: expert curated, functional experimental, drug
    (functional observational excluded -- it rests on molecular-QTL
    colocalisation, the evidence type this paper audits)
  * primary domain: melanoma x eQTLGen significant loci, the same domain as the
    self-assembled ten, so the two are comparable
  * trait-matched primary if it yields >= 30 evaluable loci, else cross-trait
    primary with the trait-mismatch caveat as a headline limitation
  * a borderline count such as 29 counts as < 30

Outputs: 113a_goldstandard_coords.tsv  (coordinates + class + trait, NO gene_id)
         113b_coverage.tsv
         113_console.log
"""
import io
import sys
import urllib.request

import numpy as np
import pandas as pd

MR = r"D:/R_ex/MR"
UA = {"User-Agent": "Mozilla/5.0 (attribution-benchmark-survey)"}
URL = ("https://raw.githubusercontent.com/opentargets/genetics-gold-standards/"
       "master/gold_standards/processed/gwas_gold_standards.191108.tsv")

GENE_COL = "gold_standard_info.gene_id"
KEEP_CLASSES = {"expert curated", "functional experimental", "drug"}
EXCLUDE_CLASS = "functional observational"
LOCUS_KB = 1000          # same 1 Mb single-linkage as the rest of the paper
MATCH_KB = 1000          # a gold-standard sentinel counts for a locus within 1 Mb


class Tee:
    def __init__(self, p):
        self.f = io.open(p, "w", encoding="utf-8")

    def write(self, s):
        enc = sys.__stdout__.encoding or "utf-8"
        sys.__stdout__.write(s.encode(enc, "replace").decode(enc, "replace"))
        self.f.write(s)

    def flush(self):
        sys.__stdout__.flush()
        self.f.flush()


sys.stdout = Tee(f"{MR}/113_console.log")

print("=" * 78)
print("Step 113 -- can the external benchmark run? (gene names NOT read)")
print("=" * 78)

raw = urllib.request.urlopen(
    urllib.request.Request(URL, headers=UA), timeout=90
).read().decode("utf-8", "replace")
gs = pd.read_csv(io.StringIO(raw), sep="\t", dtype=str)
print(f"downloaded gold standards: {len(gs)} rows, {gs.shape[1]} columns")

# ---- drop gene identity immediately and enforce it ------------------------
assert GENE_COL in gs.columns, "expected gene column missing -- check the file"
gs = gs.drop(columns=[GENE_COL])
assert GENE_COL not in gs.columns
assert not any("gene" in c.lower() for c in gs.columns), \
    "a gene-bearing column survived the drop"
print(f"dropped {GENE_COL}; remaining columns carry no gene identity")

CH = "sentinel_variant.locus_GRCh38.chromosome"
PO = "sentinel_variant.locus_GRCh38.position"
CL = "gold_standard_info.evidence.class"
TR = "trait_info.reported_trait_name"
TS = "trait_info.standard_trait_name"

gs = gs[[CH, PO, CL, TR, TS, "gold_standard_info.highest_confidence"]].copy()
gs.columns = ["chr", "pos", "cls", "trait_reported", "trait_std", "confidence"]
gs = gs.dropna(subset=["chr", "pos"])
gs["pos"] = pd.to_numeric(gs["pos"], errors="coerce")
gs = gs.dropna(subset=["pos"])
gs["pos"] = gs["pos"].astype(int)
gs["chr"] = gs["chr"].astype(str).str.replace("chr", "", regex=False)
gs.to_csv(f"{MR}/113a_goldstandard_coords.tsv", sep="\t", index=False)
print(f"usable GRCh38 sentinels: {len(gs)}")

print("\nevidence classes present (an entry may list several, '|'-joined):")
flat = gs.cls.fillna("").str.split("|").explode().str.strip()
for c, n in flat.value_counts().items():
    mark = "  <- EXCLUDED" if c == EXCLUDE_CLASS else ""
    print(f"  {c:<28}{n:>6}{mark}")

# an entry qualifies if ANY of its evidence classes is in the kept set
def qualifies(cell):
    if not isinstance(cell, str):
        return False
    return any(c.strip() in KEEP_CLASSES for c in cell.split("|"))


def only_excluded(cell):
    if not isinstance(cell, str):
        return False
    cs = {c.strip() for c in cell.split("|")}
    return cs == {EXCLUDE_CLASS}


gs["keep"] = gs.cls.map(qualifies)
print(f"\nentries with at least one qualifying class: {int(gs.keep.sum())}"
      f" of {len(gs)}")
print(f"entries resting ONLY on {EXCLUDE_CLASS}: "
      f"{int(gs.cls.map(only_excluded).sum())}")

# ---- our significant loci -------------------------------------------------
d = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
sig = d[d.fdr < 0.05].copy()
sig["chr"] = sig["chr"].astype(str)
loci = (sig.groupby("locus")
           .agg(chr=("chr", "first"), lo=("pos", "min"), hi=("pos", "max"),
                n_genes=("symbol", "nunique"), known=("known", "any"))
           .reset_index())
print(f"\npipeline significant loci (melanoma x eQTLGen): {len(loci)}")
print(f"  of which known-locus: {int(loci.known.sum())}")


def intersect(gsub, lsub):
    """count pipeline loci with >=1 gold-standard sentinel within MATCH_KB"""
    hits = []
    for _, L in lsub.iterrows():
        g = gsub[gsub.chr == L.chr]
        if not len(g):
            continue
        near = g[(g.pos >= L.lo - MATCH_KB * 1000) &
                 (g.pos <= L.hi + MATCH_KB * 1000)]
        if len(near):
            hits.append(dict(locus=L.locus, chr=L.chr, n_sentinels=len(near),
                             n_genes_named=L.n_genes, known=L.known))
    return pd.DataFrame(hits)


MEL = ("melanoma", "naevus", "nevus", "pigment", "skin cancer", "tanning",
       "hair colo", "freckl", "sunburn")


def is_mel(row):
    blob = f"{row.trait_reported} {row.trait_std}".lower()
    return any(k in blob for k in MEL)


gs_keep = gs[gs.keep]
gs_mel = gs_keep[gs_keep.apply(is_mel, axis=1)]

print("\n" + "=" * 78)
print("coverage")
print("=" * 78)

rows = []
for label, gsub in [("trait-matched (melanoma-like)", gs_mel),
                    ("cross-trait (all qualifying)", gs_keep),
                    ("cross-trait, no class filter", gs)]:
    inter = intersect(gsub, loci)
    n = len(inter)
    print(f"  {label:<32} sentinels={len(gsub):>5}   evaluable loci={n:>3}")
    rows.append(dict(version=label, n_sentinels=len(gsub), n_loci=n))
    if label.startswith("cross-trait (all"):
        inter.to_csv(f"{MR}/113b_coverage.tsv", sep="\t", index=False)

pd.DataFrame(rows).to_csv(f"{MR}/113c_coverage_summary.tsv", sep="\t",
                          index=False)

n_mel = rows[0]["n_loci"]
n_cross = rows[1]["n_loci"]

print("\n" + "=" * 78)
print("registered decision (S34 §3.3 and §9.1a amendment 3)")
print("=" * 78)
print(f"  trait-matched evaluable loci: {n_mel}  (floor is 30)")
if n_mel >= 30:
    print("  -> trait-matched is the PRIMARY analysis")
elif n_cross >= 30:
    print("  -> trait-matched below floor; CROSS-TRAIT becomes primary,")
    print("     and the trait-mismatch caveat is a headline limitation:")
    print("     an accepted gene for another trait need not be the causal gene")
    print("     for this one, so a miss may be real trait specificity.")
    print("     The cross-trait rate is a LOWER BOUND, not an error rate.")
else:
    print("  -> both below floor: S34 §5 branch E. The check is NOT run,")
    print("     and the ten self-assembled loci keep their existing caveat.")
    print("     Report the intersection counts and the reason, nothing else.")
print("=" * 78)
print("\nNo gene name has been read. Scoring is step 114 and only proceeds if")
print("the branch above allows it.")
