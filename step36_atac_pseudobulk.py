#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 36  GSE282266: does chromatin accessibility at glycolytic genes open on
         TCR activation, and does TPI1 behave differently from SPSB2?

Design fixed in advance (Step 35 header); nothing below is tuned after seeing
results. Sample-level pseudobulk is the unit: 4 sets x 4 timepoints = 16
samples, so n = 4 per timepoint. This mirrors the Step 26 lesson that cells are
not replicates.

NOTE ON MODALITY: 10x Multiome profiles NUCLEI. The RNA is snRNA-seq -- MALAT1
dominates, cytoplasmic mRNAs such as GAPDH are depleted. Per-cell detection
rates are therefore not comparable to scRNA-seq and are never reported here;
RNA is summed to pseudobulk. ATAC is unaffected and is the primary readout.

PRE-SPECIFIED
  T1  TPI1 TSS accessibility rises rest -> 15h; SPSB2 stays flat
  T2  glycolytic genes gain more accessibility than baseline-matched background
  T3  positive control: activation markers rise, quiescence markers fall
      If T3 fails, the activation itself did not work and T1/T2 are void.

Windows fixed at +/-2 kb (primary) and +/-10 kb (sensitivity).
Output: 37a-37d
"""
import sys
import gzip
import os
import time
import numpy as np
import pandas as pd
from scipy import stats

D = r"D:/Downloads/GSE282266"
MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
TPS = ["rest", "act_2.5", "act_5", "act_15"]
SETS = [1, 2, 3, 4]
WINDOWS = {"w2k": 2000, "w10k": 10000}

GLYCO = ["SLC2A1", "SLC2A3", "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP",
         "ALDOA", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2",
         "PKM", "LDHA", "PGM1", "PFKFB3", "PFKFB4"]
ACT_UP = ["IL2RA", "CD69", "MYC", "IFNG", "IL2", "TNFRSF9", "BATF"]
ACT_DOWN = ["SELL", "CCR7", "TCF7", "KLF2", "LEF1"]
FOCUS = ["TPI1", "SPSB2"]


def load_features(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            p = line.rstrip("\n").split("\t")
            rows.append((i, p[0], p[2] if len(p) > 2 else "?",
                         p[3] if len(p) > 3 else "",
                         p[4] if len(p) > 4 else "", p[5] if len(p) > 5 else ""))
    df = pd.DataFrame(rows, columns=["row", "name", "kind", "chrom", "start", "end"])
    df["start"] = pd.to_numeric(df.start, errors="coerce")
    df["end"] = pd.to_numeric(df.end, errors="coerce")
    return df[df.kind == "Gene Expression"].copy(), df[df.kind == "Peaks"].copy()


def gene_window_sums(genes, peaks, rowtot, win):
    """Sum peak pseudobulk counts inside each gene's TSS +/- win."""
    out = np.zeros(len(genes))
    npk = np.zeros(len(genes), dtype=int)
    gi = genes.reset_index(drop=True)
    for ch, gsub in gi.groupby("chrom"):
        psub = peaks[peaks.chrom == ch].dropna(subset=["start", "end"])
        if not len(psub):
            continue
        ps = psub.start.values.astype(np.int64)
        pe = psub.end.values.astype(np.int64)
        pr = psub.row.values
        order = np.argsort(ps)
        ps, pe, pr = ps[order], pe[order], pr[order]
        vals = rowtot.reindex(pr).fillna(0).values
        cum = np.concatenate([[0.0], np.cumsum(vals)])
        for idx, tss in zip(gsub.index, gsub.start.values):
            if np.isnan(tss):
                continue
            lo, hi = tss - win, tss + win
            # peaks whose start <= hi and end >= lo; approximate with start-sorted
            i0 = np.searchsorted(ps, lo - 10000, side="left")
            i1 = np.searchsorted(ps, hi, side="right")
            if i1 <= i0:
                continue
            sel = (pe[i0:i1] >= lo)
            if not sel.any():
                continue
            out[idx] = vals[i0:i1][sel].sum()
            npk[idx] = int(sel.sum())
    gi["atac"] = out
    gi["npeaks"] = npk
    return gi


recs = []
t_start = time.time()
for tp in TPS:
    for s in SETS:
        base = os.path.join(D, f"GSE282266_{tp}_set_{s}")
        if not os.path.exists(base + "_matrix.mtx.gz"):
            print(f"  [missing] {tp}_set_{s}", flush=True)
            continue
        t0 = time.time()
        genes, peaks = load_features(base + "_features.tsv.gz")
        with gzip.open(base + "_matrix.mtx.gz", "rt") as fh:
            for line in fh:
                if not line.startswith("%"):
                    nrow, ncol, nnz = (int(x) for x in line.split())
                    break
        m = pd.read_csv(base + "_matrix.mtx.gz", sep=r"\s+", comment="%",
                        skiprows=1, header=None, names=["r", "c", "v"],
                        dtype={"r": np.int32, "c": np.int32, "v": np.int32},
                        engine="c")
        if len(m) == nnz + 1:
            m = m.iloc[1:]
        rowtot = m.groupby("r").v.sum()
        ncell = int(m.c.max())
        del m

        gene_tot = rowtot.reindex(genes.row).fillna(0).sum()
        peak_tot = rowtot.reindex(peaks.row).fillna(0).sum()
        genes = genes.reset_index(drop=True)
        genes["rna"] = rowtot.reindex(genes.row).fillna(0).values

        rec = dict(timepoint=tp, set=s, ncell=ncell, n_genes=len(genes),
                   n_peaks=len(peaks), gene_tot=gene_tot, peak_tot=peak_tot)
        per_gene = genes[["name", "chrom", "start", "rna"]].copy()
        for tag, w in WINDOWS.items():
            gw = gene_window_sums(genes, peaks, rowtot, w)
            per_gene[f"atac_{tag}"] = gw.atac.values
            per_gene[f"npeaks_{tag}"] = gw.npeaks.values
        per_gene["timepoint"] = tp
        per_gene["set"] = s
        per_gene["gene_tot"] = gene_tot
        per_gene["peak_tot"] = peak_tot
        recs.append(per_gene)
        print(f"  {tp}_set_{s}: {ncell:,} nuclei, {len(peaks):,} peaks, "
              f"{time.time()-t0:.0f}s", flush=True)

allg = pd.concat(recs, ignore_index=True)
allg.to_csv(f"{MR}/37a_GSE282266_pseudobulk.tsv.gz", sep="\t", index=False,
            compression="gzip")
print(f"\ntotal {time.time()-t_start:.0f}s | rows {len(allg):,}")

# ---------------------------------------------------------------- normalise
allg["rna_cpm"] = 1e6 * allg.rna / allg.gene_tot
for tag in WINDOWS:
    allg[f"acc_{tag}"] = 1e6 * allg[f"atac_{tag}"] / allg.peak_tot

# ---------------------------------------------------------------- T3 control
print("\n" + "=" * 74)
print("T3  activation positive control (pseudobulk CPM, n=4 per timepoint)")
print("=" * 74)
ctrl = []
for g in ACT_UP + ACT_DOWN:
    sub = allg[allg.name == g]
    if not len(sub):
        continue
    means = sub.groupby("timepoint").rna_cpm.mean().reindex(TPS)
    r = stats.spearmanr(np.arange(len(TPS)).repeat(
        [sub[sub.timepoint == t].shape[0] for t in TPS]),
        np.concatenate([sub[sub.timepoint == t].rna_cpm.values for t in TPS]))
    ctrl.append(dict(gene=g, expected="up" if g in ACT_UP else "down",
                     **{t: round(means[t], 2) for t in TPS},
                     rho=r.statistic, p=r.pvalue))
ctrl = pd.DataFrame(ctrl)
print(ctrl.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
ok_up = ((ctrl.expected == "up") & (ctrl.rho > 0) & (ctrl.p < .05)).sum()
ok_dn = ((ctrl.expected == "down") & (ctrl.rho < 0) & (ctrl.p < .05)).sum()
print(f"\n  markers behaving as expected: up {ok_up}/{(ctrl.expected=='up').sum()}, "
      f"down {ok_dn}/{(ctrl.expected=='down').sum()}")
ctrl.to_csv(f"{MR}/37b_activation_control.tsv", sep="\t", index=False)
if ok_up + ok_dn < 4:
    print("  *** POSITIVE CONTROL FAILED -- T1/T2 below are not interpretable ***")

# ---------------------------------------------------------------- T1
print("\n" + "=" * 74)
print("T1  TPI1 vs SPSB2 TSS accessibility across activation")
print("=" * 74)
t1 = []
for g in FOCUS:
    sub = allg[allg.name == g]
    for tag in WINDOWS:
        means = sub.groupby("timepoint")[f"acc_{tag}"].mean().reindex(TPS)
        x = np.concatenate([[i] * sub[sub.timepoint == t].shape[0]
                            for i, t in enumerate(TPS)])
        y = np.concatenate([sub[sub.timepoint == t][f"acc_{tag}"].values for t in TPS])
        r = stats.spearmanr(x, y)
        lfc = np.log2((means["act_15"] + 1) / (means["rest"] + 1))
        t1.append(dict(gene=g, window=tag,
                       **{t: round(means[t], 1) for t in TPS},
                       log2FC_15h_vs_rest=lfc, rho=r.statistic, p=r.pvalue,
                       npeaks=int(sub[f"npeaks_{tag}"].median())))
t1 = pd.DataFrame(t1)
print(t1.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
t1.to_csv(f"{MR}/37c_TPI1_vs_SPSB2_accessibility.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- T2
print("\n" + "=" * 74)
print("T2  glycolytic genes vs baseline-matched background")
print("=" * 74)
res2 = []
for tag in WINDOWS:
    piv = allg.pivot_table(index="name", columns="timepoint",
                           values=f"acc_{tag}", aggfunc="mean")
    piv = piv.reindex(columns=TPS).dropna()
    piv = piv[piv["rest"] > 0]
    piv["lfc"] = np.log2((piv["act_15"] + 1) / (piv["rest"] + 1))
    piv["bin"] = pd.qcut(piv["rest"], 20, labels=False, duplicates="drop")
    gl = piv[piv.index.isin(GLYCO)]
    bg = piv[~piv.index.isin(GLYCO)]
    obs = gl.lfc.mean()
    rng = np.random.default_rng(1)
    null = []
    for _ in range(5000):
        draw = [bg[bg.bin == b].lfc.sample(1, random_state=None).values[0]
                for b in gl.bin if (bg.bin == b).any()]
        null.append(np.mean(draw))
    null = np.array(null)
    p = (1 + (null >= obs).sum()) / (1 + len(null))
    print(f"  {tag}:  glycolytic mean log2FC = {obs:+.4f} (n={len(gl)})  |  "
          f"matched background {null.mean():+.4f} +/- {null.std():.4f}  |  P = {p:.4g}")
    res2.append(dict(window=tag, n_glyco=len(gl), obs_lfc=obs,
                     null_mean=null.mean(), null_sd=null.std(), p=p))
    top = gl.sort_values("lfc", ascending=False)
    print("     per gene (log2FC 15h vs rest): " +
          ", ".join(f"{i}={v:+.2f}" for i, v in top.lfc.items()))
pd.DataFrame(res2).to_csv(f"{MR}/37d_glyco_class_accessibility.tsv", sep="\t",
                          index=False)
print("\nwritten: 37a-37d")
