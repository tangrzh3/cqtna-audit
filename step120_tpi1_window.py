"""Step 120 -- recover the TPI1 instrument's effect across all eight activation
profiles, so the activation-window panel stops depending on hardcoded numbers.

Why this exists:
`figures/make_fig9_part2.py` panel a carries BETA/SE/PV for variant 12:6867132
across the eight profiles as literals, annotated "(Step 62 核查所得)". There is no
step62 script in the repository and no table holding those values. That panel is
what GB promotes to a main figure of its own (Fig 6), and it is what corrects the
paper's earlier reading of the 16 h window from "the effect is amplified" to
"the effect is merely measurable there" -- so it should not rest on literals.

Same class of problem as the 96a rebuild (Step 111) and the mismatched-locus
controls (Step 119): the number is almost certainly right, but nothing regenerates it.

What is recomputed: the cis-eQTL statistics for TPI1 (ENSG00000111669) at the
instrument variant 12:6867132 in each of the eight profiles, read straight from the
cleaned per-profile parquet files. The recomputed values are compared against the
literals currently in the figure script and any disagreement is printed.

Output: 120a_tpi1_instrument_window.tsv
"""
import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyarrow.compute as pc

MR = r"D:/R_ex/MR"
PARQ = r"D:/Downloads/CD4_eqtl_step1_clean"

GENE = "ENSG00000111669"          # TPI1
CHROM, POS = "12", 6867132        # the instrument selected at Naive 16 h

PROFILES = [("Naive", "0h", "CD4_Naive_uns_0h"),
            ("Naive", "16h", "CD4_Naive_stim_16h"),
            ("Naive", "40h", "CD4_Naive_stim_40h"),
            ("Naive", "5d", "CD4_Naive_stim_5d"),
            ("Memory", "0h", "CD4_Memory_uns_0h"),
            ("Memory", "16h", "CD4_Memory_stim_16h"),
            ("Memory", "40h", "CD4_Memory_stim_40h"),
            ("Memory", "5d", "CD4_Memory_stim_5d")]

# 现写死在 figures/make_fig9_part2.py 面板 a 里的数，用来比对
LITERAL = {
    ("Naive", "0h"):   (-0.427, 0.107, 1.51e-4),
    ("Naive", "16h"):  (-0.177, 0.021, 1.74e-12),
    ("Naive", "40h"):  (-0.085, 0.040, 3.50e-2),
    ("Naive", "5d"):   (+0.017, 0.079, 8.27e-1),
    ("Memory", "0h"):  (-0.151, 0.086, 8.48e-2),
    ("Memory", "16h"): (-0.193, 0.029, 4.22e-9),
    ("Memory", "40h"): (-0.049, 0.030, 1.11e-1),
    ("Memory", "5d"):  (+0.033, 0.077, 6.75e-1),
}
GWS = 5e-8


def read_profile(fname):
    f = pq.ParquetFile(os.path.join(PARQ, fname + "_step1_clean.parquet"))
    keep = []
    for batch in f.iter_batches(columns=["gene_id", "chr", "pos", "beta", "se",
                                         "pval", "eaf", "effect_allele",
                                         "other_allele", "ma_count"]):
        t = batch.filter(pc.and_(pc.equal(batch.column("gene_id"), GENE),
                                 pc.equal(batch.column("pos"), POS)))
        if t.num_rows:
            keep.append(t.to_pandas())
    return pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()


def main():
    rows = []
    print(f"TPI1 ({GENE}) at {CHROM}:{POS} across the eight activation profiles")
    print("-" * 92)
    print(f"{'profile':<14}{'beta':>10}{'se':>10}{'pval':>13}"
          f"{'  vs literal beta/se/p':<28}{'GWS':>6}")
    print("-" * 92)

    for lineage, tp, fname in PROFILES:
        d = read_profile(fname)
        if not len(d):
            print(f"{lineage+' '+tp:<14}  variant not present in this profile")
            rows.append(dict(lineage=lineage, timepoint=tp, beta=np.nan,
                             se=np.nan, pval=np.nan, eaf=np.nan,
                             has_instrument=False, agrees_with_figure_literal=False))
            continue
        r = d.iloc[0]
        lb, ls_, lp = LITERAL[(lineage, tp)]
        ok = (abs(r.beta - lb) <= 0.0015 and abs(r.se - ls_) <= 0.0015
              and abs(np.log10(r.pval) - np.log10(lp)) <= 0.05)
        print(f"{lineage+' '+tp:<14}{r.beta:>10.4f}{r.se:>10.4f}{r.pval:>13.3e}"
              f"   {lb:+.3f} / {ls_:.3f} / {lp:.2e}  "
              f"{'ok' if ok else '*** DIFFERS ***':<12}"
              f"{'yes' if r.pval < GWS else '':>6}")
        rows.append(dict(lineage=lineage, timepoint=tp, beta=float(r.beta),
                         se=float(r.se), pval=float(r.pval), eaf=float(r.eaf),
                         effect_allele=r.effect_allele, other_allele=r.other_allele,
                         has_instrument=bool(r.pval < GWS),
                         agrees_with_figure_literal=bool(ok)))

    d = pd.DataFrame(rows)
    d["abs_beta"] = d.beta.abs()
    d.to_csv(f"{MR}/120a_tpi1_instrument_window.tsv", sep="\t", index=False)

    inst = d[d.has_instrument]
    biggest = d.loc[d.abs_beta.idxmax()]
    tightest = d.loc[d.se.idxmin()]
    print("-" * 92)
    print(f"  genome-wide significant in {len(inst)}/8 profiles: "
          f"{', '.join(inst.lineage + ' ' + inst.timepoint)}")
    print(f"  largest |beta|  : {biggest.lineage} {biggest.timepoint} "
          f"({biggest.beta:+.3f}, se {biggest.se:.3f})")
    print(f"  smallest se     : {tightest.lineage} {tightest.timepoint} "
          f"({tightest.beta:+.3f}, se {tightest.se:.3f})")
    print(f"  -> the window is a statement about precision, not amplitude")
    print(f"  all eight agree with the figure's literals: {bool(d.agrees_with_figure_literal.all())}")
    print("\nwrote 120a_tpi1_instrument_window.tsv")


if __name__ == "__main__":
    main()
