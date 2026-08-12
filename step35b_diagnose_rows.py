"""Diagnose the feature-row mapping in GSE282266 before trusting anything.

Symptom: GAPDH appeared detected in only 18.7% of cells, which is impossible.
Checks, in order:
  1. does the features file line count match the matrix header's row count?
  2. do genes come before or after peaks in the features file?
  3. what are the highest-total rows, and are they genes or peaks?
  4. where does GAPDH's row actually sit and what is its total?
"""
import gzip
import numpy as np
import pandas as pd

BASE = r"D:/Downloads/GSE282266/GSE282266_rest_set_1"

# ---- 1. header vs features ------------------------------------------------
with gzip.open(BASE + "_matrix.mtx.gz", "rt") as fh:
    hdr = []
    for line in fh:
        hdr.append(line.rstrip("\n"))
        if not line.startswith("%"):
            break
print("matrix header lines:")
for h in hdr:
    print("   ", h[:100])
dims = hdr[-1].split()
nrow_hdr, ncol_hdr, nnz_hdr = (int(x) for x in dims)
print(f"  -> declared: {nrow_hdr:,} rows x {ncol_hdr:,} cols, {nnz_hdr:,} nonzero")

rows = []
with gzip.open(BASE + "_features.tsv.gz", "rt") as fh:
    for i, line in enumerate(fh, start=1):
        p = line.rstrip("\n").split("\t")
        rows.append((i, p[0], p[2] if len(p) > 2 else "?"))
feat = pd.DataFrame(rows, columns=["row", "name", "kind"])
print(f"  features file lines: {len(feat):,}")
print(f"  MATCH: {len(feat) == nrow_hdr}")

# ---- 2. ordering ----------------------------------------------------------
print("\nfirst 3 and last 3 features:")
for _, r in pd.concat([feat.head(3), feat.tail(3)]).iterrows():
    print(f"   row {r.row:>7}  {r.kind:<16} {r['name'][:45]}")
g = feat[feat.kind == "Gene Expression"]
p = feat[feat.kind == "Peaks"]
print(f"  Gene Expression rows {g.row.min():,}-{g.row.max():,}  (n={len(g):,})")
print(f"  Peaks rows           {p.row.min():,}-{p.row.max():,}  (n={len(p):,})")

# ---- 3/4. row totals ------------------------------------------------------
print("\nreading matrix ...")
m = pd.read_csv(BASE + "_matrix.mtx.gz", sep=r"\s+", comment="%",
                skiprows=1, header=None, names=["r", "c", "v"],
                dtype={"r": np.int32, "c": np.int32, "v": np.int32}, engine="c")
print(f"  rows parsed: {len(m):,}  (declared nnz {nnz_hdr:,})")
if len(m) == nnz_hdr + 1:
    print("  -> one extra row: the dimension line was parsed as data; dropping it")
    m = m.iloc[1:]
print(f"  r range {m.r.min()}-{m.r.max()} | c range {m.c.min()}-{m.c.max()}")

tot = m.groupby("r").v.sum()
det = m.groupby("r").size()
kind = feat.set_index("row").kind
name = feat.set_index("row")["name"]

top = tot.sort_values(ascending=False).head(15)
print("\ntop 15 rows by total count:")
for r, v in top.items():
    print(f"   row {r:>7}  {kind.get(r,'?'):<16} {str(name.get(r,'?'))[:38]:<40}"
          f" total={v:>12,}  cells={det.get(r,0):>7,} ({100*det.get(r,0)/ncol_hdr:5.1f}%)")

print("\nspecific genes:")
for gene in ["GAPDH", "ACTB", "B2M", "TMSB4X", "MALAT1", "TPI1", "SPSB2", "CD4"]:
    hit = feat[(feat.name == gene) & (feat.kind == "Gene Expression")]
    if not len(hit):
        print(f"   {gene:<8} NOT FOUND in features")
        continue
    for _, h in hit.iterrows():
        r = h.row
        print(f"   {gene:<8} row {r:>7}  total={tot.get(r,0):>10,}  "
              f"cells={det.get(r,0):>7,} ({100*det.get(r,0)/ncol_hdr:5.1f}%)")

gs = tot[tot.index.isin(set(g.row))].sum()
ps = tot[tot.index.isin(set(p.row))].sum()
print(f"\ntotal counts  genes={gs:,}  peaks={ps:,}  ratio peaks/genes={ps/max(gs,1):.2f}")
