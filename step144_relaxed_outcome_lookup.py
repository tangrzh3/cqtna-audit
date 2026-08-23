#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 144 -- outcome-side statistics for the relaxed instrument set.

The benchmark in step145 needs, for every relaxed instrument, the meta outcome
P value at that position, so each method's nominated loci can be split by
whether the outcome GWAS already reaches 5e-8 there. 213 of the 1,158 relaxed
positions are in 13_meta; the other 945 are not, so they come from the meta
summary statistics directly.

Output: 144a_relaxed_outcome_lookup.tsv   chr, pos, outcome beta/se/p, source
        144b_console.log
"""
import gzip
import io
import os
import sys

import pandas as pd

MR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CQTNA_DIR") or \
     os.path.dirname(os.path.abspath(__file__))
SUMSTATS = os.path.join(MR, "meta_melanoma_final.tsv.gz")

rel = pd.read_csv(os.path.join(MR, "17_relaxed_instruments_meta.tsv"), sep="\t")
rel["key"] = rel.chr.astype(str) + ":" + rel.pos.astype(int).astype(str)
want = set(rel.key)
print("relaxed instrument positions: %d unique" % len(want), flush=True)

# what 13_meta already carries
meta = pd.read_csv(os.path.join(MR, "13_meta_locus_annotation.tsv"), sep="\t")
have = {}
for _, r in meta.drop_duplicates("SNP").iterrows():
    if r.SNP in want:
        have[r.SNP] = (float(r.beta_outcome), float(r.se_outcome),
                       float(r.pval), "13_meta")
print("  found in 13_meta: %d" % len(have), flush=True)

missing = want - set(have)
print("  to read from meta sumstats: %d" % len(missing), flush=True)

# one streaming pass; the file is ~765 MB gzipped
if missing:
    if not os.path.exists(SUMSTATS):
        print("  MISSING FILE: %s" % SUMSTATS, flush=True)
    else:
        n = 0
        with gzip.open(SUMSTATS, "rt") as f:
            header = f.readline().lstrip("#").rstrip("\n").split("\t")
            ic, ip = header.index("chrom"), header.index("pos")
            ipv, ib, ise = (header.index("pval"), header.index("beta"),
                            header.index("sebeta"))
            for line in f:
                n += 1
                if n % 5000000 == 0:
                    print("    %d lines, %d of %d found"
                          % (n, len(have) - (len(want) - len(missing)) +
                             (len(want) - len(missing)), len(want)), flush=True)
                p = line.split("\t")
                k = p[ic] + ":" + p[ip]
                if k in missing:
                    try:
                        have[k] = (float(p[ib]), float(p[ise]), float(p[ipv]),
                                   "meta_sumstats")
                    except ValueError:
                        have[k] = (float("nan"), float("nan"), float("nan"),
                                   "meta_sumstats_unparsed")
                    missing.discard(k)
                    if not missing:
                        break
        print("  read %d lines; still unresolved: %d" % (n, len(missing)),
              flush=True)

rows = []
for k in sorted(want):
    b, se, p, src = have.get(k, (float("nan"), float("nan"), float("nan"),
                                 "not found"))
    c, pos = k.split(":")
    rows.append(dict(key=k, chr=c, pos=int(pos), beta_outcome=b,
                     se_outcome=se, pval_outcome=p, source=src))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(MR, "144a_relaxed_outcome_lookup.tsv"), sep="\t",
           index=False)
print("\nwrote 144a_relaxed_outcome_lookup.tsv (%d rows)" % len(out), flush=True)
print(out.source.value_counts().to_string(), flush=True)
