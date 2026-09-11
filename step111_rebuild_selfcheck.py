#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 111 -- rebuild the self-administered attribution check from its inputs.

96a_attribution_selfcheck.tsv is cited by Supplementary S27 and underlies the
published "six of ten" result, but no script in the repository or its history
produces it, and three of its gene lists are stored truncated with a trailing
ellipsis, so the published count could not be verified from the deposited table.

This script reconstructs it from 92c_locus_annotated.tsv, the eQTLGen-exposure MR
output with the locus assignment already applied. The reconstruction rule was
recovered by testing candidate rules against the stored table:

    a locus's NAMED genes are the FDR < 0.05 records at that locus whose own
    instrument lies within 1 Mb of a known melanoma/naevus/pigmentation lead
    (the `known` flag), grouped by the 1 Mb single-linkage locus id

That rule reproduces n_named for 20 of 20 loci; taking all significant genes
regardless of the known flag reproduces only 19 of 20, differing at 22_1, where
TOM1 is significant but sits outside any known lead. The 20 loci are exactly the
known-locus subset of the 30 significant loci reported for this exposure.

The accepted-gene column cannot be regenerated: the ten-locus list was assembled
by hand and was hardcoded in the lost script. It is carried over verbatim from
the stored table, which is legitimate because that column is short and none of it
is truncated. This is a recovery, not a regeneration, and is labelled as such.

Inputs:  92c_locus_annotated.tsv, 96a_attribution_selfcheck.tsv
Outputs: 111a_selfcheck_rebuilt.tsv, 111_console.log
"""
import os
import io
import sys

import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")


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


sys.stdout = Tee(f"{MR}/111_console.log")

d = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
old = pd.read_csv(f"{MR}/96a_attribution_selfcheck.tsv", sep="\t",
                  dtype=str).fillna("")

print("=" * 78)
print("Step 111 -- rebuilding 96a from 92c under the recovered rule")
print("=" * 78)

sig = d[d.fdr < 0.05]
all_loci = sig.groupby("locus").agg(known=("known", "any"))
print(f"significant loci (FDR < 0.05): {len(all_loci)}")
print(f"  of which known-locus:        {int(all_loci.known.sum())}"
      f"   <- the denominator of this check")

named = (sig[sig.known].groupby("locus")["symbol"]
         .apply(lambda s: sorted(set(s.dropna()))))

rows = []
acc_map = dict(zip(old.locus, old.accepted))
for locus, genes in named.items():
    acc_raw = acc_map.get(locus, "")
    accepted = [g for g in acc_raw.split(";") if g]
    lenient = any(g in genes for g in accepted) if accepted else None
    strict = (len(genes) == 1 and genes[0] in accepted) if accepted else None
    rows.append(dict(locus=locus, n_named=len(genes), named=";".join(genes),
                     accepted=acc_raw,
                     lenient_hit="" if lenient is None else lenient,
                     strict_hit="" if strict is None else strict))

new = pd.DataFrame(rows).sort_values("locus").reset_index(drop=True)
new.to_csv(f"{MR}/111a_selfcheck_rebuilt.tsv", sep="\t", index=False)

# ---- verify against the stored table where it is verifiable ---------------
print("\n" + "=" * 78)
print("verification against the stored 96a")
print("=" * 78)
o = old.set_index("locus")
n = new.set_index("locus")
print(f"locus sets identical: {set(o.index) == set(n.index)}")

n_ok = sum(int(o.loc[L, "n_named"]) == int(n.loc[L, "n_named"]) for L in o.index)
print(f"n_named agrees: {n_ok}/{len(o)}")

trunc = [L for L in o.index if o.loc[L, "named"].endswith("...")]
gene_ok = gene_bad = 0
for L in o.index:
    if L in trunc:
        continue
    if sorted(o.loc[L, "named"].split(";")) == sorted(n.loc[L, "named"].split(";")):
        gene_ok += 1
    else:
        gene_bad += 1
        print(f"  gene-list mismatch at {L}")
        print(f"    stored: {o.loc[L,'named']}")
        print(f"    rebuilt: {n.loc[L,'named']}")
print(f"gene lists agree on the {len(o)-len(trunc)} untruncated loci: "
      f"{gene_ok}/{len(o)-len(trunc)}")
print(f"\nthe {len(trunc)} loci stored truncated are now complete:")
for L in trunc:
    print(f"  {L:<8} stored {o.loc[L,'named'][:52]}")
    print(f"  {'':<8} full   {n.loc[L,'named']}")

# ---- the published count, now verifiable ----------------------------------
print("\n" + "=" * 78)
print("the ten self-assembled loci, both definitions, from complete gene lists")
print("=" * 78)
ev = new[new.accepted != ""].copy()
nl = int((ev.lenient_hit == True).sum())
ns = int((ev.strict_hit == True).sum())
print(f"evaluable loci: {len(ev)}")
print(f"  lenient (accepted gene among those named): {nl}/{len(ev)}")
print(f"  strict  (named alone and correctly):       {ns}/{len(ev)}")

pub = int((old.hit == "True").sum())
print(f"\npublished 'hit' column: {pub}/10")
print(f"lenient rebuilt:        {nl}/10  -> "
      f"{'CONFIRMED' if nl == pub else 'DOES NOT MATCH'}")

print("\nper-locus:")
print(f"  {'locus':<8}{'accepted':<14}{'n':>3}  {'len':<6}{'str':<6}named")
for _, r in ev.iterrows():
    print(f"  {r.locus:<8}{r.accepted[:13]:<14}{r.n_named:>3}  "
          f"{str(r.lenient_hit):<6}{str(r.strict_hit):<6}{r.named[:46]}")

print(f"\nwrote 111a_selfcheck_rebuilt.tsv")
print("\nNote: `accepted` is carried over from the stored table, not regenerated.")
print("The hand-assembled ten-locus list has no machine-readable provenance and")
print("that remains true after this rebuild -- which is exactly why S34 exists.")
