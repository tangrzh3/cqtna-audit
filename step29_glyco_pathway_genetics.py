#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 29  Pathway-level genetic evidence for CD4 glycolysis, without instruments.

Step 28 established the structural problem: of 28 glycolytic genes, only TPI1
and ENO1 carry a genome-wide significant dynamic-eQTL instrument, and ENO1 is
null (p=0.29). Single-gene MR therefore cannot go pathway-level. Two questions
that do not require instruments:

  ANALYSIS 1  Is melanoma GWAS signal enriched at the CD4 cis-eQTLs of
              glycolytic genes as a class, relative to expression- and
              eQTL-strength-matched background genes?
              -> if yes, the pathway claim survives without instruments.

  ANALYSIS 2  Is the activation-window (16h) pattern seen for TPI1 a property
              of glycolytic genes as a class, or is TPI1 an isolated case?
              -> tests "genetic control of the glycolytic switch concentrates
                 in the activation window".

Both use the full cis windows in CD4_eqtl_step1_clean (8 profiles x ~8M rows),
not the instrument subset, so genes without instruments are included.
Outputs 35b-35f.
"""
import sys
import os
import glob
import gzip
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import stats

PARQ = "D:/Downloads/CD4_eqtl_step1_clean"
MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
GWAS = os.path.join(MR, "meta_melanoma_final.tsv.gz")
rng = np.random.default_rng(1)

ORDER = ["CD4_Naive_uns_0h", "CD4_Naive_stim_16h", "CD4_Naive_stim_40h",
         "CD4_Naive_stim_5d", "CD4_Memory_uns_0h", "CD4_Memory_stim_16h",
         "CD4_Memory_stim_40h", "CD4_Memory_stim_5d"]
STAGE = {"uns_0h": "0h", "stim_16h": "16h", "stim_40h": "40h", "stim_5d": "5d"}

gmap = pd.read_csv(f"{MR}/35a_glyco_ensembl_map.tsv", sep="\t")
GLYCO_ENSG = set(gmap.ENSEMBL)
ENSG2SYM = dict(zip(gmap.ENSEMBL, gmap.SYMBOL))
print(f"glycolytic genes: {gmap.SYMBOL.nunique()} symbols -> {len(GLYCO_ENSG)} ENSG\n")

# ==========================================================================
# pass 1: top cis-eQTL per gene per profile, for ALL genes
# ==========================================================================
cache = f"{MR}/35b_top_eqtl_per_gene_profile.tsv"
if os.path.exists(cache):
    top = pd.read_csv(cache, sep="\t")
    print(f"loaded cached {cache}: {len(top):,} rows")
else:
    rows = []
    for path in sorted(glob.glob(os.path.join(PARQ, "*.parquet"))):
        prof = os.path.basename(path).replace("_step1_clean.parquet", "")
        df = pq.read_table(path, columns=["gene_id", "variant_id", "chr", "pos",
                                          "eaf", "pval"]).to_pandas()
        idx = df.groupby("gene_id", sort=False)["pval"].idxmin()
        t = df.loc[idx].copy()
        t["profile"] = prof
        rows.append(t)
        print(f"  {prof}: {len(df):,} rows -> {len(t):,} genes", flush=True)
        del df
    top = pd.concat(rows, ignore_index=True)
    top.to_csv(cache, sep="\t", index=False)
    print(f"\nwritten {cache}: {len(top):,} rows")

top["is_glyco"] = top.gene_id.isin(GLYCO_ENSG)
top["stage"] = [next((v for k, v in STAGE.items() if p.endswith(k)), p) for p in top.profile]
top["lineage"] = np.where(top.profile.str.contains("Naive"), "Naive", "Memory")
print(f"\ngenes covered: {top.gene_id.nunique():,} | glycolytic present: "
      f"{top.loc[top.is_glyco, 'gene_id'].nunique()} of {len(GLYCO_ENSG)}")

# ==========================================================================
# ANALYSIS 2 -- is the 16h peak a class property?
# ==========================================================================
print("\n" + "=" * 70)
print("ANALYSIS 2  activation-window bias of the strongest eQTL")
print("=" * 70)

best = top.loc[top.groupby(["gene_id", "lineage"])["pval"].idxmin()]
tab = (best.groupby(["is_glyco", "stage"]).size().unstack(fill_value=0))
frac = tab.div(tab.sum(axis=1), axis=0) * 100
print("\npercentage of genes whose strongest eQTL falls in each stage:")
print(frac.round(1).to_string())

res2 = []
for st in ["0h", "16h", "40h", "5d"]:
    a = int(tab.loc[True, st]) if st in tab.columns else 0
    b = int(tab.loc[True].sum() - a)
    c = int(tab.loc[False, st]) if st in tab.columns else 0
    d = int(tab.loc[False].sum() - c)
    odds, p = stats.fisher_exact([[a, b], [c, d]])
    res2.append(dict(stage=st, glyco_n=a, glyco_total=a + b,
                     glyco_pct=100 * a / (a + b) if a + b else np.nan,
                     bg_n=c, bg_total=c + d, bg_pct=100 * c / (c + d),
                     odds_ratio=odds, p=p))
res2 = pd.DataFrame(res2)
res2.to_csv(f"{MR}/35c_activation_window_class_test.tsv", sep="\t", index=False)
print("\nFisher tests, glycolytic vs all other genes:")
print(res2.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

# ==========================================================================
# ANALYSIS 1 -- melanoma GWAS signal at glycolytic cis-eQTLs
# ==========================================================================
print("\n" + "=" * 70)
print("ANALYSIS 1  melanoma GWAS signal at CD4 cis-eQTLs of glycolytic genes")
print("=" * 70)

EQTL_P = 1e-4          # keep plausible regulatory variants only
sel = top[top.pval < EQTL_P].copy()
sel["snp"] = sel.chr.astype(str) + ":" + sel.pos.astype(str)
print(f"gene x profile pairs with top eQTL p<{EQTL_P:g}: {len(sel):,} "
      f"({sel.loc[sel.is_glyco].shape[0]} glycolytic)")

need = set(sel.snp)
print(f"unique SNPs to look up in the melanoma meta GWAS: {len(need):,}")

gcache = f"{MR}/35d_gwas_lookup.tsv"
if os.path.exists(gcache):
    gw = pd.read_csv(gcache, sep="\t")
    print(f"loaded cached {gcache}: {len(gw):,} rows")
else:
    keep = []
    with gzip.open(GWAS, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        # meta_melanoma_final.tsv.gz header:
        # #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ...
        cc, cp = ci["#chrom"], ci["pos"]
        cb, cs = ci["beta"], ci["sebeta"]
        cm = ci["mlogp"]
        for n, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            s = f[cc] + ":" + f[cp]
            if s in need:
                try:
                    b, se = float(f[cb]), float(f[cs])
                    z = abs(b / se) if se > 0 else np.nan
                except Exception:
                    z = np.nan
                keep.append((s, z, float(f[cm]) if cm is not None and f[cm] else np.nan))
            if n % 5_000_000 == 0 and n:
                print(f"    scanned {n:,} lines, matched {len(keep):,}", flush=True)
    gw = pd.DataFrame(keep, columns=["snp", "absz", "mlogp"]).drop_duplicates("snp")
    gw.to_csv(gcache, sep="\t", index=False)
    print(f"written {gcache}: {len(gw):,} rows")

sel = sel.merge(gw, on="snp", how="left")
ok = sel.dropna(subset=["absz"])
print(f"\nmatched in GWAS: {len(ok):,} of {len(sel):,} "
      f"({ok.is_glyco.sum()} glycolytic)")

# match background on eQTL-p decile and EAF decile
ok = ok.copy()
ok["p_bin"] = pd.qcut(np.log10(ok.pval), 10, labels=False, duplicates="drop")
ok["eaf_bin"] = pd.qcut(ok.eaf, 5, labels=False, duplicates="drop")
gl = ok[ok.is_glyco]
bg = ok[~ok.is_glyco]

obs = gl.absz.mean()
null = []
for _ in range(5000):
    draw = []
    for (pb, eb), sub in gl.groupby(["p_bin", "eaf_bin"]):
        pool = bg[(bg.p_bin == pb) & (bg.eaf_bin == eb)]
        if len(pool) == 0:
            pool = bg[bg.p_bin == pb]
        if len(pool) == 0:
            continue
        # random_state=rng, not None. Until 2026-09-11 this read
        # random_state=None, so the 5000 matched draws used numpy's global
        # unseeded state and the three summary columns of
        # 35e_pathway_enrichment.tsv came out different on every run --
        # measured across two runs of the SAME container image on the SAME
        # machine, which leaves nothing but the draw itself to explain it.
        #   The seeded generator on line 38 existed the whole time and was
        # never reached from here, which is worse than having no seed at all:
        # it tells anyone reading the file that this permutation reproduces.
        # S54 section 1 rules randomness out as an explanation for a moved
        # number precisely on that assumption.
        #   The spread was ordinary Monte Carlo error for 5000 draws -- the
        # two empirical P values were 0.5363 and 0.5303 against a standard
        # error of 0.0071, 0.85 SE apart, and neither is near significance --
        # so no verdict ever moved. What was broken was reproducibility, not
        # the result. Authorised by the author 2026-09-11; registered in S54
        # section 9.12.4, and the value it now fixes on is one arbitrary draw
        # replacing another, not a correction of a wrong number.
        draw.append(pool.absz.sample(len(sub), replace=True,
                                     random_state=rng).values)
    if draw:
        null.append(np.concatenate(draw).mean())
null = np.array(null)
p_emp = (1 + (null >= obs).sum()) / (1 + len(null))

print(f"\nmean |z| in the melanoma meta GWAS:")
print(f"  glycolytic cis-eQTLs : {obs:.4f}  (n = {len(gl)})")
print(f"  matched background   : {null.mean():.4f} +/- {null.std():.4f}")
print(f"  empirical P (one-sided, 5000 matched draws) = {p_emp:.4g}")
u, p_mw = stats.mannwhitneyu(gl.absz, bg.absz, alternative="greater")
print(f"  unmatched Mann-Whitney P = {p_mw:.4g}  (context only; not the test)")

pd.DataFrame([dict(n_glyco=len(gl), n_background=len(bg), obs_mean_absz=obs,
                   null_mean=null.mean(), null_sd=null.std(),
                   p_matched_permutation=p_emp, p_mannwhitney_unmatched=p_mw)]
             ).to_csv(f"{MR}/35e_pathway_enrichment.tsv", sep="\t", index=False)

per_gene = (gl.assign(symbol=gl.gene_id.map(ENSG2SYM))
              .groupby("symbol")
              .agg(n_profiles=("absz", "size"), max_absz=("absz", "max"),
                   mean_absz=("absz", "mean"), best_eqtl_p=("pval", "min"))
              .sort_values("max_absz", ascending=False))
per_gene.to_csv(f"{MR}/35f_glyco_gene_gwas_signal.tsv", sep="\t")
print("\nper-gene melanoma signal at its CD4 cis-eQTL:")
print(per_gene.to_string(float_format=lambda v: f"{v:.3g}"))
print("\nwritten: 35b-35f")
