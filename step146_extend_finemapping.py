#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 146 -- extend the outcome-side fine-mapping to the rest of module H's loci.

WHY
step60a-c fine-mapped four regions chosen for a different question (MC1R as the
positive control, PARP1 and ZFYVE19 as the coloc false-positive and
false-negative cases, TPI1 as the worked example). Module H splits the eight
significant bounded loci into three the outcome GWAS already reaches 5e-8 on and
five it does not, and asks what the second group contains. Four of those five
have never been fine-mapped: SMC2, KANSL1, MDM4 and KIAA0040.

This extracts their regional summary statistics and builds their LD, reusing
step60a/60b's conventions unchanged: 500 kb half-window, the same 1000G EUR
panel of 525 unrelated samples, the same plink2 call.

READING RULE -- already frozen, and tightened relative to registration in
PREREG_power_trajectory.md section 9.1 after FinnGen's official in-sample
fine-mapping became available:
    >= 2 credible sets   NOT admissible (proxy LD can split spuriously)
    exactly 1, stable in L   admissible, and conservative
    0                    no distinguishable outcome-side signal at this power
Nothing here relaxes it.

Outputs: regions/{name}_{finngen,meta}.tsv and regions/{name}_ld.*
         146a_new_regions.tsv, 146b_console.log
"""
import gzip
import os
import subprocess
import sys

import pandas as pd

MR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CQTNA_DIR") or \
     os.path.dirname(os.path.abspath(__file__))
REG = os.path.join(MR, "regions")
PLINK = os.path.join(MR, "bin", "plink2.exe")
PFILE = os.path.join(MR, "ref", "all_hg38")
os.makedirs(REG, exist_ok=True)

HALF = 500_000                    # unchanged from step60a
# Centres are the lead significant record at each locus, taken from 140b.
NEW = {
    "SMC2":     ("9", 104_125_300),
    "KANSL1":   ("17", 45_592_652),
    "MDM4":     ("1", 203_769_637),
    "KIAA0040": ("1", 175_143_931),
}
SOURCES = {
    "finngen": os.path.join(MR, "finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz"),
    "meta": os.path.join(MR, "meta_melanoma_final.tsv.gz"),
}

# ---- 1. regional summary statistics, one pass per source -------------------
for src, path in SOURCES.items():
    if not os.path.exists(path):
        print("[missing] %s" % path, flush=True)
        continue
    want = {n: (c, p - HALF, p + HALF) for n, (c, p) in NEW.items()}
    buf = {n: [] for n in want}
    with gzip.open(path, "rt") as f:
        head = f.readline().lstrip("#").rstrip("\n").split("\t")
        ic, ip = head.index("chrom"), head.index("pos")
        n = 0
        for line in f:
            n += 1
            p = line.rstrip("\n").split("\t")
            ch, po = p[ic], int(p[ip])
            for name, (c, lo, hi) in want.items():
                if ch == c and lo <= po <= hi:
                    buf[name].append(p)
                    break
    for name in want:
        d = pd.DataFrame(buf[name], columns=head)
        d.to_csv(os.path.join(REG, "%s_%s.tsv" % (name, src)),
                 sep="\t", index=False)
        print("  %-9s %-8s %6d variants" % (name, src, len(d)), flush=True)

# ---- 2. LD, same panel and same call as step60b ----------------------------
keep = os.path.join(REG, "eur.keep")
if not os.path.exists(keep):
    with open(os.path.join(MR, "ref", "1kg_eur.fam")) as fh, open(keep, "w") as fo:
        for line in fh:
            f = line.split()
            fo.write("%s\t%s\n" % (f[0], f[1]))

def run(cmd, name):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print("  [plink2 failed] %s: %s" % (name, " ".join(cmd[-4:])), flush=True)
        print(p.stdout[-1000:], flush=True)
        print(p.stderr[-1000:], flush=True)
        return False
    return True


rows = []
for name, (c, pos) in NEW.items():
    out = os.path.join(REG, "%s_ld" % name)
    if os.path.exists(out + ".unphased.vcor1"):
        print("  [have LD] %s" % name, flush=True)
    else:
        # Two stages, exactly as step60b: the panel's .pvar is zst-compressed, so
        # `vzs` is required, and the variant IDs must be set to chr:pos:ref:alt
        # because that is the key step60c aligns the summary statistics on.
        tmp = os.path.join(REG, "%s_tmp" % name)
        ok = run([PLINK, "--pfile", PFILE, "vzs", "--keep", keep,
                  "--chr", c, "--from-bp", str(pos - HALF),
                  "--to-bp", str(pos + HALF),
                  "--snps-only", "--max-alleles", "2", "--maf", "0.01",
                  "--set-all-var-ids", "@:#:$r:$a",
                  "--new-id-max-allele-len", "60", "missing",
                  "--rm-dup", "exclude-all",
                  "--make-pgen", "--out", tmp], name)
        if ok:
            ok = run([PLINK, "--pfile", tmp, "--r-unphased", "square",
                      "ref-based", "--out", out], name)
        print("  [%s] %s" % ("LD built" if ok else "LD FAILED", name), flush=True)
    nvar = 0
    vf = out + ".unphased.vcor1.vars"
    if os.path.exists(vf):
        nvar = sum(1 for _ in open(vf))
    rows.append(dict(region=name, chrom=c, centre=pos,
                     half_window=HALF, n_ld_variants=nvar,
                     genome_wide_sig="FALSE"))

pd.DataFrame(rows).to_csv(os.path.join(MR, "146a_new_regions.tsv"),
                          sep="\t", index=False)
print("\nwrote 146a_new_regions.tsv", flush=True)
