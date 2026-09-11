#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 52  (analysis C) Is the CD4 glycolysis programme a person-level trait in
         circulating blood?

WHY THIS COMES BEFORE THE IRC TEST
Something that varies only from cell to cell within a person, with no structure
between people, cannot be a blood biomarker no matter how well it correlates
with outcome. This is a prerequisite, and unlike the IRC test it needs donors
rather than an outcome variable, so all 10 donors are usable and n is not 4v4.

MEASURE
Intraclass correlation of the per-cell glycolysis module score across donors:
ICC = between-donor variance / total variance, from a one-way random-effects
decomposition over CD4 T cells.

*** TWO-SIDED CONTROLS, fixed in advance ***
An ICC on its own is uninterpretable without knowing what this assay gives for
things that certainly are, and certainly are not, person-level traits.
  upper control  XIST and RPS4Y1 -- sex-determined, so ICC must approach 1
  lower control  50 random gene modules matched on detection rate -- these give
                 the ICC floor produced by noise plus donor batch alone
The glycolysis ICC is read against those two, never against zero.

A high ICC alone would still not prove biology: donor batch (processing day,
depth, ambient RNA) inflates every ICC, which is exactly what the random-module
floor measures. Only the margin above that floor is interpretable.

Depth is regressed out per cell before scoring, as in Steps 43-44.

Output: 52a-52c
"""
import sys
import glob
import gzip
import os
import re
import tarfile
import time
import numpy as np
import pandas as pd

D = r"D:/Downloads/GSE199994"
WORK = r"D:/Downloads/GSE199994/extracted"
MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
rng = np.random.default_rng(1)

GLYCO = ["SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA",
         "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3"]
CD4_PANEL = ["CD4", "IL7R", "CD40LG", "MAL", "LTB", "TRAT1", "ANXA1", "CCR7", "AQP3"]
CD8_PANEL = ["CD8A", "CD8B", "GZMK", "NKG7", "CCL5", "GZMB", "PRF1", "KLRD1", "GNLY"]
TCELL = ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "CD2"]
UPPER_CTRL = ["XIST", "RPS4Y1"]
N_RANDOM = 50
MIN_DETECT = 0.05


def extract(sample):
    """Unpack one patient tar into a per-sample folder, once."""
    out = os.path.join(WORK, sample)
    if os.path.isdir(out) and os.listdir(out):
        return out
    os.makedirs(out, exist_ok=True)
    tp = os.path.join(D, f"GSE199994_{sample}.tar.gz")
    with tarfile.open(tp) as t:
        members = [m for m in t.getmembers()
                   if re.search(r"(barcodes|features|matrix)\.(tsv|mtx)\.gz$",
                                m.name)]
        if not members:
            members = t.getmembers()
        t.extractall(out, members=members)
    return out


def find(folder, pat):
    hits = glob.glob(os.path.join(folder, "**", pat), recursive=True)
    return hits[0] if hits else None


def load_sample(sample):
    folder = extract(sample)
    fp = find(folder, "*features.tsv.gz") or find(folder, "*features.tsv")
    mp = find(folder, "*matrix.mtx.gz") or find(folder, "*matrix.mtx")
    if fp is None or mp is None:
        raise FileNotFoundError(f"{sample}: features/matrix not found in {folder}")
    op = gzip.open if fp.endswith(".gz") else open
    rows = []
    with op(fp, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            p = line.rstrip("\n").split("\t")
            rows.append((i, p[1] if len(p) > 1 else p[0],
                         p[2] if len(p) > 2 else "Gene Expression"))
    feat = pd.DataFrame(rows, columns=["row", "name", "kind"])
    opener = gzip.open if mp.endswith(".gz") else open
    with opener(mp, "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                nrow, ncol, nnz = (int(x) for x in line.split())
                break
    m = pd.read_csv(mp, sep=r"\s+", comment="%", skiprows=1, header=None,
                    names=["r", "c", "v"],
                    dtype={"r": np.int32, "c": np.int32, "v": np.int32}, engine="c")
    if len(m) == nnz + 1:
        m = m.iloc[1:]
    genes = feat[feat.kind == "Gene Expression"]
    g = m[m.r.isin(set(genes.row))]
    return genes, g, int(m.c.max())


def raw_values(genes, g, ncell, wanted):
    """Per-cell log1p(CP10K) for the requested genes.

    *** Standardisation is deliberately NOT done here. ***
    A first version z-scored each gene WITHIN each sample before forming module
    scores, which centres every donor at zero and therefore removes exactly the
    between-donor variance the ICC is meant to measure. The upper controls
    caught it: XIST gave an ICC of 0.0006 and RPS4Y1 0.0031, where sex-determined
    genes must approach 1. Genes are standardised once, across the pooled cells
    of all donors, after concatenation.
    """
    n2r = dict(zip(genes.name, genes.row))
    tot = g.groupby("c").v.sum().reindex(range(1, ncell + 1), fill_value=0).values
    nfeat = g.groupby("c").size().reindex(range(1, ncell + 1), fill_value=0).values
    need = sorted(set(wanted) & set(n2r))
    sub = g[g.r.isin({n2r[x] for x in need})]
    out = {}
    for x in need:
        v = sub[sub.r == n2r[x]].set_index("c").v.reindex(
            range(1, ncell + 1), fill_value=0).values
        out[x] = np.log1p(1e4 * v / np.maximum(tot, 1))
    return pd.DataFrame(out), tot, nfeat


def icc(values, donor):
    """One-way random effects ICC: between-donor variance / total."""
    df = pd.DataFrame(dict(v=values, d=donor)).dropna()
    grp = df.groupby("d").v
    k = grp.size()
    if len(k) < 3:
        return np.nan
    grand = df.v.mean()
    msb = (k * (grp.mean() - grand) ** 2).sum() / (len(k) - 1)
    msw = ((grp.transform("count") - 1) * 0).sum()
    msw = df.groupby("d").v.apply(lambda s: ((s - s.mean()) ** 2).sum()).sum()
    msw = msw / (len(df) - len(k))
    k0 = (len(df) - (k ** 2).sum() / len(df)) / (len(k) - 1)
    var_b = max((msb - msw) / k0, 0)
    return var_b / (var_b + msw) if (var_b + msw) > 0 else np.nan


SAMPLES = [f"P{i}" for i in range(1, 9)] + ["HD1", "HD2"]
frames = []
random_sets = None
for s in SAMPLES:
    tp = os.path.join(D, f"GSE199994_{s}.tar.gz")
    if not os.path.exists(tp):
        print(f"  [missing] {s}")
        continue
    t0 = time.time()
    genes, g, ncell = load_sample(s)
    det = (g.groupby("r").size() / ncell)
    r2n = dict(zip(genes.row, genes.name))
    detect = pd.Series({r2n[r]: v for r, v in det.items()})
    if random_sets is None:
        pool = detect[(detect >= MIN_DETECT) & (detect < 0.9)].index.tolist()
        random_sets = {f"rand{i}": list(rng.choice(pool, len(GLYCO), replace=False))
                       for i in range(N_RANDOM)}
    sets = {"Glyco": GLYCO, "CD4s": CD4_PANEL, "CD8s": CD8_PANEL, "Ts": TCELL,
            **{f"ctrl_{x}": [x] for x in UPPER_CTRL}, **random_sets}
    wanted = {x for v in sets.values() for x in v}
    vals, tot, nfeat = raw_values(genes, g, ncell, wanted)
    vals["donor"] = s
    vals["nfeat"] = nfeat
    vals = vals[vals.nfeat >= 200]
    frames.append(vals)
    print(f"  {s}: {ncell:,} nuclei kept {len(vals):,}, "
          f"{time.time()-t0:.0f}s", flush=True)
    del g

if not frames:
    raise SystemExit("no samples available -- run step51 first")
raw = pd.concat(frames, ignore_index=True).fillna(0.0)
gene_cols = [c for c in raw.columns if c not in ("donor", "nfeat")]

# standardise ONCE across the pooled cells of all donors
z = raw[gene_cols]
z = (z - z.mean()) / z.std().replace(0, np.nan)
z = z.fillna(0.0)

sets = {"Glyco": GLYCO, "CD4s": CD4_PANEL, "CD8s": CD8_PANEL, "Ts": TCELL,
        **{f"ctrl_{x}": [x] for x in UPPER_CTRL}, **random_sets}
d = pd.DataFrame(dict(donor=raw.donor.values, nfeat=raw.nfeat.values))
for nm, gs in sets.items():
    pres = [x for x in gs if x in z.columns]
    d[nm] = z[pres].mean(axis=1).values if pres else np.nan

# CD4 T cells: T-cell score positive and CD4 score above CD8 score
d = d[(d.Ts > 0) & (d.CD4s > d.CD8s)].reset_index(drop=True)
print(f"\nCD4 T cells across {d.donor.nunique()} donors: {len(d):,}")
print(d.donor.value_counts().to_string())

# depth-residualise every score before the ICC
X = np.column_stack([np.ones(len(d)), np.log1p(d.nfeat)])
for col in [c for c in d.columns if c not in ("donor", "nfeat")]:
    y = d[col].values
    ok = np.isfinite(y)
    if ok.sum() < 100:
        continue
    b, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    d.loc[ok, col] = y[ok] - X[ok] @ b

rows = []
for col in ["Glyco"] + [f"ctrl_{x}" for x in UPPER_CTRL] + list(random_sets):
    if col not in d.columns:
        continue
    rows.append(dict(measure=col, icc=icc(d[col].values, d.donor.values),
                     kind=("glycolysis" if col == "Glyco"
                           else "upper_control" if col.startswith("ctrl_")
                           else "random_floor")))
res = pd.DataFrame(rows)
res.to_csv(f"{MR}/52a_icc_all.tsv", sep="\t", index=False)

floor = res[res.kind == "random_floor"].icc
gly = res[res.measure == "Glyco"].icc.iloc[0]
print("\n" + "=" * 70)
print("ICC of the CD4 glycolysis module across donors")
print("=" * 70)
print(f"  glycolysis          {gly:.4f}")
for _, r in res[res.kind == "upper_control"].iterrows():
    print(f"  {r.measure:<20}{r.icc:.4f}   (upper control, expected near 1)")
print(f"  random modules      {floor.mean():.4f} +/- {floor.std():.4f}  "
      f"(95th pct {np.percentile(floor, 95):.4f})")
print(f"\n  glycolysis exceeds the random floor: "
      f"{gly > np.percentile(floor, 95)}")
print(f"  margin above floor mean: {gly - floor.mean():+.4f}")
pd.DataFrame(dict(glyco_icc=[gly], floor_mean=[floor.mean()],
                  floor_p95=[np.percentile(floor, 95)],
                  exceeds=[bool(gly > np.percentile(floor, 95))])
             ).to_csv(f"{MR}/52b_icc_verdict.tsv", sep="\t", index=False)
# the random modules are saved per cell as well: the noise floor has to be
# recomputable within any subgroup, since a floor estimated over a mixed cohort
# is not the right comparator for an ICC computed inside a homogeneous one
keep_cols = (["donor", "Glyco", "nfeat"] + [f"ctrl_{x}" for x in UPPER_CTRL]
             + list(random_sets))
d[keep_cols].to_csv(f"{MR}/52c_cell_scores.tsv.gz", sep="\t",
                    index=False, compression="gzip")

# per-donor means of the controls -- needed to interpret them. A sex-determined
# gene can only show a high ICC if the cohort actually contains both sexes, so
# an apparently weak upper control must be checked against cohort composition
# rather than taken at face value.
per_donor = d.groupby("donor")[["Glyco"] + [f"ctrl_{x}" for x in UPPER_CTRL]].mean()
per_donor["n_cells"] = d.groupby("donor").size()
per_donor.to_csv(f"{MR}/52d_per_donor_means.tsv", sep="\t")
print("\n" + "=" * 70)
print("per-donor means (controls are sex-determined; check cohort composition)")
print("=" * 70)
print(per_donor.round(3).to_string())
xi = per_donor[f"ctrl_XIST"]
print(f"\n  XIST spread across donors: {xi.min():+.3f} to {xi.max():+.3f}")
print(f"  donors above the XIST midpoint: "
      f"{(xi > (xi.min()+xi.max())/2).sum()} of {len(xi)}")
print("\nwritten: 52a-52d")
