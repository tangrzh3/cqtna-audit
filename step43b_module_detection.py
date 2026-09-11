"""Were the module genes in Step 43 detectable at all? A module built from genes
that snRNA-seq cannot measure, in a system that cannot produce the cell state,
returns a null for reasons that have nothing to do with the axis being tested."""
import gzip
import sys
import os
import numpy as np
import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")

D = r"D:/Downloads/GSE282266"
MODULES = {
    "Glycolysis": ["SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKL", "PFKP",
                   "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM",
                   "LDHA", "PFKFB3"],
    "Activation": ["IL2RA", "CD69", "TNFRSF9", "TNFRSF4", "ICOS", "BATF",
                   "NFKB1", "REL", "IRF4", "MYC", "IL2"],
    "Proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA", "TYMS",
                      "STMN1", "TUBA1B", "RRM2", "UBE2C"],
    "OXPHOS": ["NDUFA4", "NDUFB2", "SDHB", "UQCRB", "COX5A", "COX7C", "ATP5F1A"],
    "Exhaustion": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"],
    "Treg": ["FOXP3", "IKZF2", "IL2RA", "CTLA4", "TNFRSF18"],
    "Th1": ["TBX21", "IFNG", "CXCR3", "IL12RB2"],
    "Th2": ["GATA3", "IL4", "IL5", "IL13", "CCR4"],
    "Th17": ["RORC", "IL17A", "IL23R", "CCR6"],
    "Tfh": ["BCL6", "CXCR5", "PDCD1", "IL21"],
    "Cytotoxic": ["GZMB", "GZMK", "PRF1", "NKG7", "GNLY"],
    "IFN_response": ["ISG15", "IFI6", "MX1", "OAS1", "IFIT3", "STAT1", "IRF7"],
}

base = os.path.join(D, "GSE282266_act_15_set_1")
rows = []
with gzip.open(base + "_features.tsv.gz", "rt") as fh:
    for i, line in enumerate(fh, start=1):
        p = line.rstrip("\n").split("\t")
        rows.append((i, p[0], p[2] if len(p) > 2 else "?"))
feat = pd.DataFrame(rows, columns=["row", "name", "kind"])
genes = feat[feat.kind == "Gene Expression"]
name2row = dict(zip(genes.name, genes.row))

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
ncell = int(m.c.max())
det = m.groupby("r").size()

print(f"act_15_set_1: {ncell:,} nuclei\n")
out = []
for mod, gs in MODULES.items():
    print(f"{mod}")
    pcts = []
    for g in gs:
        r = name2row.get(g)
        if r is None:
            print(f"    {g:<10} NOT IN ANNOTATION")
            out.append(dict(module=mod, gene=g, pct=np.nan)); continue
        pct = 100 * det.get(r, 0) / ncell
        pcts.append(pct)
        flag = "  <-- effectively undetected" if pct < 5 else ""
        print(f"    {g:<10}{pct:6.1f}%{flag}")
        out.append(dict(module=mod, gene=g, pct=pct))
    if pcts:
        print(f"    -> median detection {np.median(pcts):.1f}%, "
              f"{sum(p < 5 for p in pcts)}/{len(pcts)} genes under 5%\n")
pd.DataFrame(out).to_csv(os.path.join(MR, "43d_module_gene_detection.tsv"),
                         sep="\t", index=False)
print("written: 43d")
