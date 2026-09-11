"""Step 87 -- Part II patient stratification in a second tumour type (GSE235863, HCC).

Executes manuscript/PREREG_hcc_part2_generalisation.md exactly. Nothing here is
chosen after seeing an expression value: the 16-gene signature, the CD4
definition and its Treg sensitivity, the sample filter, the permutation scheme,
the direction of every hypothesis, the positive controls and the reading table
are all fixed in that document, as is the fact that the pre-treatment tumour arm
cannot reach significance (4 vs 2 -> minimum one-sided p = 0.0667).

Output: 87a (sample-level scores), 87b (arm-level results), 87c (positive controls)
"""
import sys
import os
import collections
from math import comb
from itertools import combinations

import h5py
import numpy as np
import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
H5 = f"{MR}/hcc/GSE235863/GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad"
CHUNK = 20000
MIN_CD4 = 20

SIG16 = ["ALDOA", "ENO1", "GAPDH", "GPI", "HK1", "HK2", "LDHA", "PFKFB3",
         "PFKL", "PFKP", "PGAM1", "PGK1", "PKM", "SLC2A1", "SLC2A3", "TPI1"]
CTRL = ["PTPRC", "CD8A", "CD8B", "CD40LG", "IL7R", "CD4"]
RESP = {"P5": "R", "P11": "R", "P18": "R", "P27": "R",
        "P1": "NR", "P15": "NR", "P26": "NR", "P51": "NR", "P52": "NR"}
TREG = "CD4_C08_FOXP3"


def read_cat(f, k):
    g = f["obs"][k]
    if isinstance(g, h5py.Group):
        c = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
        return np.array([c[i] for i in g["codes"][:]])
    return g[:]


def main():
    f = h5py.File(H5, "r")
    genes = np.array([x.decode() for x in f["var"]["_index"][:]])
    gidx = {g: i for i, g in enumerate(genes)}
    want = SIG16 + CTRL
    missing = [g for g in want if g not in gidx]
    print(f"target genes missing: {missing if missing else 'none'} "
          f"(signature 16: {sum(g in gidx for g in SIG16)}/16)")
    cols = np.array([gidx[g] for g in want])
    col_rank = {c: i for i, c in enumerate(cols)}

    sample = read_cat(f, "sample")
    major = read_cat(f, "major_cluster")
    sub = read_cat(f, "sub_cluster")
    n_cells = len(sample)
    total_obs = f["obs"]["total_counts"][:].astype(np.float64)

    X = f["X"]
    indptr = X["indptr"][:]
    data_d, ind_d = X["data"], X["indices"]

    # ---- pass over cells, accumulating log-normalised expression per group ----
    # groups we need: (sample, cd4_nonTreg), (sample, cd4_all), and cluster-level for PCs
    acc = collections.defaultdict(lambda: np.zeros(len(cols)))
    cnt = collections.Counter()
    ptprc_pos, ptprc_tot = 0, 0
    ptprc_local = col_rank[gidx["PTPRC"]]
    raw_int = True
    lib_check = []

    for start in range(0, n_cells, CHUNK):
        stop = min(start + CHUNK, n_cells)
        lo, hi = indptr[start], indptr[stop]
        d = data_d[lo:hi].astype(np.float64)
        ix = ind_d[lo:hi]
        ptr = indptr[start:stop + 1] - lo
        if raw_int and start == 0:
            raw_int = bool(np.all(d[:5000] == np.round(d[:5000])))
        rows = np.repeat(np.arange(stop - start), np.diff(ptr))
        # library size from the data itself, cross-checked against obs.total_counts
        libs = np.bincount(rows, weights=d, minlength=stop - start)
        if start == 0:
            lib_check = (libs[:10].copy(), total_obs[start:start + 10].copy())
        keep = np.isin(ix, cols)
        if keep.any():
            r = rows[keep]
            c = np.array([col_rank[v] for v in ix[keep]])
            v = np.log1p(1e4 * d[keep] / np.maximum(libs[r], 1))
            gsamp = sample[start:stop]
            gmaj = major[start:stop]
            gsub = sub[start:stop]
            for rr, cc, vv in zip(r, c, v):
                s, m, sb = gsamp[rr], gmaj[rr], gsub[rr]
                if m == "CD4T":
                    acc[(s, "all")][cc] += vv
                    if sb != TREG:
                        acc[(s, "nonTreg")][cc] += vv
                acc[("CLUSTER", m)][cc] += vv
                if sb == "CD4_C01_CCR7":
                    acc[("CLUSTER", "CD4_naive")][cc] += vv
                if cc == ptprc_local:
                    ptprc_pos += 1
        for m in major[start:stop]:
            cnt[("CLUSTER", m)] += 1
        for s, m, sb in zip(sample[start:stop], major[start:stop], sub[start:stop]):
            if m == "CD4T":
                cnt[(s, "all")] += 1
                if sb != TREG:
                    cnt[(s, "nonTreg")] += 1
            if sb == "CD4_C01_CCR7":
                cnt[("CLUSTER", "CD4_naive")] += 1
        ptprc_tot = stop
        print(f"  cells {stop:,}/{n_cells:,}", end="\r")

    print()
    print(f"raw integer counts: {raw_int}; "
          f"library sizes from X vs obs.total_counts (first 10): "
          f"{np.allclose(lib_check[0], lib_check[1])}")

    # ---------------------------------------------------------- positive controls
    def cluster_mean(name):
        return acc[("CLUSTER", name)] / max(cnt[("CLUSTER", name)], 1)

    sig_local = [col_rank[gidx[g]] for g in SIG16]
    pc_rows = []
    my, nv = cluster_mean("Myeloid"), cluster_mean("CD4_naive")
    pc1 = float(np.mean(my[sig_local]) - np.mean(nv[sig_local]))
    pc_rows.append(dict(control="PC1 glycolysis Myeloid > naive CD4",
                        value=round(pc1, 4), pass_=bool(pc1 > 0)))
    frac = ptprc_pos / ptprc_tot
    pc_rows.append(dict(control="PC2 PTPRC detected fraction",
                        value=round(frac, 4), pass_=bool(frac >= 0.90)))
    cd4c, cd8c = cluster_mean("CD4T"), cluster_mean("CD8T")
    d_cd8 = float(np.mean([cd8c[col_rank[gidx[g]]] - cd4c[col_rank[gidx[g]]]
                           for g in ("CD8A", "CD8B")]))
    d_cd4 = float(np.mean([cd4c[col_rank[gidx[g]]] - cd8c[col_rank[gidx[g]]]
                           for g in ("CD40LG", "IL7R")]))
    pc_rows.append(dict(control="PC3a CD8A/CD8B higher in CD8T",
                        value=round(d_cd8, 4), pass_=bool(d_cd8 > 0)))
    pc_rows.append(dict(control="PC3b CD40LG/IL7R higher in CD4T",
                        value=round(d_cd4, 4), pass_=bool(d_cd4 > 0)))
    pc = pd.DataFrame(pc_rows)
    pc.to_csv(f"{MR}/87c_positive_controls.tsv", sep="\t", index=False)
    print("\n=== positive controls ===")
    print(pc.to_string(index=False))
    if not (pc_rows[0]["pass_"] and pc_rows[2]["pass_"] and pc_rows[3]["pass_"]):
        print("\nPC1 or PC3 FAILED -> reading table cell E, H1-H4 not interpreted")
        return

    # ---------------------------------------------------------- sample-level scores
    rows = []
    for (s, defn), v in acc.items():
        if s == "CLUSTER":
            continue
        n = cnt[(s, defn)]
        if n < MIN_CD4:
            continue
        pat, tp, tis = s.split("-")
        rec = dict(sample=s, patient=pat, timepoint=tp, tissue=tis,
                   response=RESP[pat], cd4_def=defn, n_cd4=n)
        for g in SIG16:
            rec[g] = v[col_rank[gidx[g]]] / n
        rows.append(rec)
    d = pd.DataFrame(rows)

    def score(sub_df, gene_set):
        m = sub_df[gene_set].to_numpy(float)
        z = (m - m.mean(0)) / np.where(m.std(0) > 0, m.std(0), 1)
        return z.mean(1)

    out = []
    for defn in ("nonTreg", "all"):
        for gset, gname in ((SIG16, "sig16"), ([g for g in SIG16 if g != "TPI1"], "sig15_noTPI1")):
            for tis, tp, hyp in (("T", "post", "H1"), ("T", "pre", "H2"),
                                 ("P", "post", "H3"), ("P", "pre", "H4")):
                sub_df = d[(d.cd4_def == defn) & (d.tissue == tis) & (d.timepoint == tp)].copy()
                if len(sub_df) < 3:
                    continue
                sub_df = sub_df.sort_values("sample").reset_index(drop=True)
                sub_df["score"] = score(sub_df, gset)
                nr = sub_df.score[sub_df.response == "NR"].to_numpy()
                r = sub_df.score[sub_df.response == "R"].to_numpy()
                delta = nr.mean() - r.mean()
                # exact permutation of the response label over patients, one-sided (NR > R)
                idx = np.arange(len(sub_df))
                lab = (sub_df.response == "NR").to_numpy()
                nnr = int(lab.sum())
                stats = []
                for combo in combinations(idx, nnr):
                    mask = np.zeros(len(sub_df), bool)
                    mask[list(combo)] = True
                    stats.append(sub_df.score[mask].mean() - sub_df.score[~mask].mean())
                stats = np.array(stats)
                p_perm = float((stats >= delta - 1e-12).sum() / len(stats))
                from scipy.stats import mannwhitneyu
                p_w = float(mannwhitneyu(nr, r, alternative="two-sided").pvalue) \
                    if len(nr) and len(r) else np.nan
                out.append(dict(hypothesis=hyp, tissue=tis, timepoint=tp, cd4_def=defn,
                                signature=gname, n_R=len(r), n_NR=len(nr),
                                delta_NR_minus_R=round(float(delta), 4),
                                p_perm_onesided=round(p_perm, 5),
                                min_possible_p=round(1 / comb(len(sub_df), nnr), 5),
                                p_wilcoxon=round(p_w, 5),
                                direction="NR>R" if delta > 0 else "R>NR"))
    res = pd.DataFrame(out).sort_values(["hypothesis", "cd4_def", "signature"])
    d.to_csv(f"{MR}/87a_sample_scores.tsv", sep="\t", index=False)
    res.to_csv(f"{MR}/87b_arm_results.tsv", sep="\t", index=False)

    print("\n=== primary (locked 16-gene signature, CD4 excluding Tregs) ===")
    print(res[(res.cd4_def == "nonTreg") & (res.signature == "sig16")].to_string(index=False))
    print("\n=== all configurations ===")
    print(res.to_string(index=False))
    print("\nwrote 87a-87c")


if __name__ == "__main__":
    main()
