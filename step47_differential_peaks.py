#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 47 (v3)  Differential chromatin accessibility along the residual glycolysis
              axis, on a consensus peak reference.

TWO EARLIER FAILURES, both self-inflicted, both recorded so they are not repeated
  v1  matched peaks between samples by NAME. Each sample has its own peak calls
      (123k-134k), so exact coordinate strings never coincide: intersection was
      empty. HANDOFF note 15c warned against exactly this, and was written before
      the script that violated it.
  v2  mapped peaks by interval but called reset_index(drop=True) on the
      per-chromosome reference, so rc.index returned WITHIN-chromosome positions.
      Every chromosome's peaks were mapped onto chr1's global indices; the give-
      away was that all top peaks, in every contrast, were on chr1.

v3 fixes both: the global index is preserved explicitly, and the reference is a
consensus built by merging peak intervals across all eight samples rather than
the supplied pooled file, which retains only 20,676 peaks and captured 12% of
each sample's calls.

Per-sample peak counts are cached keyed by the SAMPLE's own peak names, so the
reference can be changed later without re-parsing any matrix.

Output: 47a (axis), 47b (activation control)
"""
import gzip
import os
import time
import numpy as np
import pandas as pd

D = r"D:/Downloads/GSE282266"
MR = r"D:/R_ex/MR"
CACHE = os.path.join(MR, "cache_47")
os.makedirs(CACHE, exist_ok=True)
SETS = [1, 2, 3, 4]
GLYCO = ["SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA",
         "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3"]
ACTIV = ["IL2RA", "CD69", "TNFRSF9", "TNFRSF4", "ICOS", "BATF", "NFKB1",
         "REL", "IRF4", "MYC", "IL2"]
SAMPLES = [f"{tp}_set_{s}" for tp in ("rest", "act_15") for s in SETS]


def read_sample(sample):
    base = os.path.join(D, f"GSE282266_{sample}")
    fr = []
    with gzip.open(base + "_features.tsv.gz", "rt") as fh:
        for i, line in enumerate(fh, start=1):
            p = line.rstrip("\n").split("\t")
            fr.append((i, p[0], p[2] if len(p) > 2 else "?"))
    feat = pd.DataFrame(fr, columns=["row", "name", "kind"])
    with gzip.open(base + "_matrix.mtx.gz", "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                nrow, ncol, nnz = (int(x) for x in line.split())
                break
    m = pd.read_csv(base + "_matrix.mtx.gz", sep=r"\s+", comment="%", skiprows=1,
                    header=None, names=["r", "c", "v"],
                    dtype={"r": np.int32, "c": np.int32, "v": np.int32}, engine="c")
    if len(m) == nnz + 1:
        m = m.iloc[1:]
    return feat, m, int(m.c.max())


def axis_split(feat, m, ncell):
    genes = feat[feat.kind == "Gene Expression"]
    n2r = dict(zip(genes.name, genes.row))
    g = m[m.r.isin(set(genes.row))]
    tot = g.groupby("c").v.sum().reindex(range(1, ncell + 1), fill_value=0).values
    nfeat = g.groupby("c").size().reindex(range(1, ncell + 1), fill_value=0).values

    def mod(lst):
        z = []
        for x in lst:
            r = n2r.get(x)
            if r is None:
                continue
            v = g[g.r == r].set_index("c").v.reindex(
                range(1, ncell + 1), fill_value=0).values
            cp = 1e4 * v / np.maximum(tot, 1)
            sd = cp.std()
            z.append((cp - cp.mean()) / sd if sd > 0 else np.zeros(ncell))
        return np.mean(z, axis=0)

    gly, act = mod(GLYCO), mod(ACTIV)
    keep = nfeat >= 200
    X = np.column_stack([np.ones(keep.sum()), act[keep], np.log1p(nfeat[keep])])
    beta, *_ = np.linalg.lstsq(X, gly[keep], rcond=None)
    cd = pd.DataFrame(dict(cell=np.arange(1, ncell + 1)[keep],
                           act=act[keep], nfeat=nfeat[keep],
                           resid=gly[keep] - X @ beta))
    cd["abin"] = pd.qcut(cd.act, 10, labels=False, duplicates="drop")
    cd["dbin"] = pd.qcut(np.log1p(cd.nfeat), 5, labels=False, duplicates="drop")
    hi, lo = [], []
    for _, z in cd.groupby(["abin", "dbin"]):
        if len(z) < 20:
            continue
        q40, q60 = z.resid.quantile([.4, .6])
        h, l = z[z.resid > q60], z[z.resid <= q40]
        n = min(len(h), len(l))
        hi.append(h.nsmallest(n, "resid")); lo.append(l.nlargest(n, "resid"))
    return set(pd.concat(hi).cell), set(pd.concat(lo).cell)


# ---------------------------------------------- pass 1: cache per-sample counts
for sample in SAMPLES:
    cf = os.path.join(CACHE, f"v3_{sample}.parquet")
    if os.path.exists(cf):
        print(f"  [cache] {sample}")
        continue
    t0 = time.time()
    feat, m, ncell = read_sample(sample)
    pk = feat[feat.kind == "Peaks"]
    peak_rows = set(pk.row)
    sub = m[m.r.isin(peak_rows)]
    out = pd.DataFrame({"peak": pk.name.values})
    out["all"] = sub.groupby("r").v.sum().reindex(pk.row).fillna(0).values
    if sample.startswith("act_15"):
        hi, lo = axis_split(feat, m, ncell)
        out["hi"] = sub[sub.c.isin(hi)].groupby("r").v.sum().reindex(pk.row).fillna(0).values
        out["lo"] = sub[sub.c.isin(lo)].groupby("r").v.sum().reindex(pk.row).fillna(0).values
        print(f"  {sample}: {len(hi):,} hi / {len(lo):,} lo cells, "
              f"{len(pk):,} peaks, {time.time()-t0:.0f}s", flush=True)
    else:
        print(f"  {sample}: {len(pk):,} peaks, {time.time()-t0:.0f}s", flush=True)
    out.to_parquet(cf)
    del m

# ---------------------------------------------- consensus reference
def parse(names):
    p = pd.Series(names).str.extract(r"^(.+):(\d+)-(\d+)$")
    p.columns = ["chrom", "start", "end"]
    p["start"] = p.start.astype(np.int64); p["end"] = p.end.astype(np.int64)
    return p

frames = {s: pd.read_parquet(os.path.join(CACHE, f"v3_{s}.parquet")) for s in SAMPLES}
allp = pd.concat([parse(f.peak) for f in frames.values()], ignore_index=True)
allp = allp.sort_values(["chrom", "start"]).reset_index(drop=True)
merged = []
for c, g in allp.groupby("chrom", sort=False):
    s_, e_ = None, None
    for a, b in zip(g.start.values, g.end.values):
        if s_ is None:
            s_, e_ = a, b
        elif a <= e_:
            e_ = max(e_, b)
        else:
            merged.append((c, s_, e_)); s_, e_ = a, b
    if s_ is not None:
        merged.append((c, s_, e_))
ref = pd.DataFrame(merged, columns=["chrom", "start", "end"])
ref["peak_id"] = ref.chrom + ":" + ref.start.astype(str) + "-" + ref.end.astype(str)
print(f"\nconsensus reference peaks (merged across 8 samples): {len(ref):,}")

# global index preserved deliberately -- this is what v2 destroyed
ref_idx = {c: (g.start.values, g.end.values, g.index.values)
           for c, g in ref.groupby("chrom")}


def to_ref(names):
    p = parse(names)
    p["mid"] = (p.start + p.end) // 2
    out = np.full(len(p), -1, dtype=np.int64)
    for c, g in p.groupby("chrom"):
        if c not in ref_idx:
            continue
        st, en, gi = ref_idx[c]
        pos = np.searchsorted(st, g.mid.values, side="right") - 1
        ok = pos >= 0
        good = np.zeros(len(g), dtype=bool)
        good[ok] = g.mid.values[ok] <= en[pos[ok]]
        idx = np.full(len(g), -1, dtype=np.int64)
        idx[good] = gi[pos[good]]          # GLOBAL index
        out[g.index.values] = idx
    return out


def project(df, col):
    j = to_ref(df.peak.values)
    acc = np.zeros(len(ref))
    ok = j >= 0
    np.add.at(acc, j[ok], df[col].values[ok])
    return acc, ok.mean()


# ---------------------------------------------- axis contrast
hi_m, lo_m = [], []
for s in SETS:
    f = frames[f"act_15_set_{s}"]
    a, frac = project(f, "hi"); b, _ = project(f, "lo")
    hi_m.append(a); lo_m.append(b)
    print(f"  act_15_set_{s}: {100*frac:.1f}% of peaks mapped to consensus")

ax = pd.DataFrame({"peak": ref.peak_id})
for i, s in enumerate(SETS):
    th, tl = hi_m[i].sum(), lo_m[i].sum()
    ax[f"lfc{s}"] = np.log2((1e6 * hi_m[i] / th + 1) / (1e6 * lo_m[i] / tl + 1))
ax["total"] = np.sum(hi_m, axis=0) + np.sum(lo_m, axis=0)
ax["mean_lfc"] = ax[[f"lfc{s}" for s in SETS]].mean(axis=1)
ax["n_pos"] = (ax[[f"lfc{s}" for s in SETS]] > 0).sum(axis=1)
ax["consistent"] = ax.n_pos.isin([0, 4])
ax = ax.sort_values("mean_lfc", ascending=False)
ax.to_csv(f"{MR}/47a_axis_differential_peaks.tsv.gz", sep="\t", index=False,
          compression="gzip")
cons = ax[ax.consistent & (ax.total >= 100)]
print(f"\nreference peaks: {len(ax):,} | consistent 4/4 and >=100 counts: {len(cons):,}")
chrs = pd.Series([p.split(":")[0] for p in cons.head(200).peak]).value_counts()
print(f"  chromosomes among top 200 up: {dict(chrs.head(8))}")
print(f"  top 5 up  : {list(cons.head(5).peak)}")
print(f"  top 5 down: {list(cons.tail(5).peak)}")

# ---------------------------------------------- activation control
tot = {}
for tp in ("rest", "act_15"):
    acc = np.zeros(len(ref))
    for s in SETS:
        a, _ = project(frames[f"{tp}_set_{s}"], "all")
        acc += a
    tot[tp] = acc
ctl = pd.DataFrame({"peak": ref.peak_id,
                    "rest_cpm": 1e6 * tot["rest"] / tot["rest"].sum(),
                    "act_cpm": 1e6 * tot["act_15"] / tot["act_15"].sum()})
ctl["lfc"] = np.log2((ctl.act_cpm + 1) / (ctl.rest_cpm + 1))
ctl = ctl.sort_values("lfc", ascending=False)
ctl.to_csv(f"{MR}/47b_activation_differential_peaks.tsv.gz", sep="\t",
           index=False, compression="gzip")
print(f"\ncontrol: peaks with counts {(ctl.rest_cpm + ctl.act_cpm > 0).sum():,}")
chrs = pd.Series([p.split(":")[0] for p in ctl.head(200).peak]).value_counts()
print(f"  chromosomes among top 200 opening: {dict(chrs.head(8))}")
print("\nwritten: 47a, 47b")
