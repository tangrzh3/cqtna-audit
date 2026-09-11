"""Step 119 -- rebuild the disease x resource grid into one authoritative table,
and recompute the two mismatched-locus controls that had no generating script.

Why this exists (三条，都在 FIGURES_GB_mapping.md §三 记过):

  1. `94d_grid.tsv` 里混着 PREREG_generality_grid.md §8.2 明确判为作废的两行
     (eQTLGen x lung, eQTLGen x colorectal -- 坐标映射有 bug，背景已知位点只有
     2/559 与 5/559，而黑色素瘤为 58/554)。表里没有任何作废标记，
     任何读表的脚本都会把它们当成有效格子。
  2. `94e_grid_summary.tsv` 汇总的是那个旧网格(含作废行、不含 RA)，
     与现行主网格(melanoma/HCC-high/RA x 2 资源)不是一回事。
  3. 两个错配对照 (eQTLGen x melanoma 用 HCC 名单 1.30x/P=0.411;
     Soskic x HCC-low 用 melanoma 名单 0.00x/P=1.0) 只存在于预注册 §8 的记录里，
     仓库中没有 TSV、也没有能重跑出它们的脚本 -- 与上一窗口修掉的 96a 同类问题。

本脚本不改任何原表。原表保留，另出三张新表：
  119a_grid_main.tsv        每格一行，带 status 列 (main / sensitivity / void_bad_coords)
  119b_grid_summary.tsv     只对 status=main 的六格重算 U1
  119c_mismatch_controls.tsv 重算的错配对照，并与预注册记录逐位比对

打分规则与 step85/step94c 完全一致，不得改动：
  · 位点 = 染色体内工具变量位置的 1 Mb 单连锁聚类
  · 某位点为"已知" = 该位点上任一工具变量落在某个已知 lead SNP 的 1 Mb 内
  · 倍数 = (显著位点中已知的比例) / (背景位点中已知的比例)
  · P = 单侧 Fisher
"""
import sys
import os
from math import lgamma, exp

import numpy as np
import pandas as pd

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
LOCUS_KB = 1000
KNOWN_KB = 1000


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    hi = min(r1, c1)
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, hi + 1)), 1.0)


def assign_loci(df):
    """1 Mb single-linkage clustering of instrument positions within a chromosome."""
    out = {}
    for ch, sub in df.groupby("chr"):
        sub = sub.sort_values("pos")
        lid, prev = 0, None
        for _, r in sub.iterrows():
            if prev is not None and r.pos - prev > LOCUS_KB * 1000:
                lid += 1
            out[(ch, r.pos)] = f"{ch}_{lid}"
            prev = r.pos
    return out


def known_lookup(path):
    k = pd.read_csv(f"{MR}/{path}")
    k["chr"] = k["chr"].astype(str)
    by_chr = {c: np.sort(s.pos.values) for c, s in k.groupby("chr")}

    def is_known(ch, pos):
        arr = by_chr.get(str(ch))
        if arr is None:
            return False
        i = np.searchsorted(arr, pos)
        for j in (i - 1, i):
            if 0 <= j < len(arr) and abs(int(arr[j]) - int(pos)) <= KNOWN_KB * 1000:
                return True
        return False

    return is_known, len(k)


def score(d, is_known, label):
    """d needs chr, pos, fdr. Returns the enrichment of the FDR<0.05 loci."""
    d = d.copy()
    d["chr"] = d["chr"].astype(str)
    d["pos"] = d["pos"].astype(int)
    snp = d.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp)
    d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
    d["known"] = [is_known(c, p) for c, p in zip(d.chr, d.pos)]

    bg = d.groupby("locus").agg(known=("known", "any"))
    sg = d[d.fdr < 0.05].groupby("locus").agg(known=("known", "any"))
    BG_T, BG_K = len(bg), int(bg.known.sum())
    S_T, S_K = len(sg), int(sg.known.sum())
    fold = ((S_K / S_T) / (BG_K / BG_T)) if S_T and BG_K else np.nan
    p = fisher_greater(S_K, S_T - S_K, BG_K - S_K, (BG_T - BG_K) - (S_T - S_K))
    print(f"  {label:<46} {S_K}/{S_T} sig loci known, background {BG_K}/{BG_T} "
          f"= {100*BG_K/BG_T:.1f}%  ->  {fold:.2f}x, P = {p:.4g}")
    return dict(bg_loci=BG_T, bg_known=BG_K, sig_loci=S_T, sig_known=S_K,
                pct_known=round(100 * S_K / S_T, 1) if S_T else np.nan,
                fold=round(fold, 2) if S_T else np.nan, fisher_p=p)


# ---------------------------------------------------------------- main table
def build_main():
    """每格从它自己的权威来源读，不从 94d 读 -- 94d 是要被取代的那张表。"""
    rows = []

    def add(exposure, disease, status, src, note, **kw):
        rows.append(dict(exposure=exposure, disease=disease, status=status,
                         source=src, note=note, **kw))

    # --- melanoma x Soskic: 13_meta_locus_annotation, known 编码在 category 列
    mel = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
    mel[["chr", "pos"]] = mel.SNP.str.split(":", expand=True)
    mel["pos"] = mel.pos.astype(int)
    mel = mel.rename(columns={"FDR": "fdr"})
    snp = mel.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp)
    mel["locus"] = [loci[(c, p)] for c, p in zip(mel.chr, mel.pos)]
    mel["known"] = mel.category != "潜在新位点"
    bg = mel.groupby("locus").agg(known=("known", "any"))
    sg = mel[mel.fdr < 0.05].groupby("locus").agg(known=("known", "any"))
    BG_T, BG_K, S_T, S_K = len(bg), int(bg.known.sum()), len(sg), int(sg.known.sum())
    print(f"  {'Soskic x melanoma':<46} {S_K}/{S_T} sig loci known, "
          f"background {BG_K}/{BG_T} = {100*BG_K/BG_T:.1f}%  ->  "
          f"{(S_K/S_T)/(BG_K/BG_T):.2f}x")
    add("Soskic_CD4", "melanoma", "main", "13_meta_locus_annotation.tsv", "",
        bg_loci=BG_T, bg_known=BG_K, sig_loci=S_T, sig_known=S_K,
        pct_known=round(100 * S_K / S_T, 1), fold=round((S_K / S_T) / (BG_K / BG_T), 2),
        fisher_p=fisher_greater(S_K, S_T - S_K, BG_K - S_K,
                                (BG_T - BG_K) - (S_T - S_K)))

    # --- eQTLGen x melanoma: 92c 已带 known 列
    eq = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    bg = eq.groupby("locus").agg(known=("known", "any"))
    sg = eq[eq.fdr < 0.05].groupby("locus").agg(known=("known", "any"))
    BG_T, BG_K, S_T, S_K = len(bg), int(bg.known.sum()), len(sg), int(sg.known.sum())
    print(f"  {'eQTLGen x melanoma':<46} {S_K}/{S_T} sig loci known, "
          f"background {BG_K}/{BG_T} = {100*BG_K/BG_T:.1f}%  ->  "
          f"{(S_K/S_T)/(BG_K/BG_T):.2f}x")
    add("eQTLGen_blood", "melanoma", "main", "92c_locus_annotated.tsv", "",
        bg_loci=BG_T, bg_known=BG_K, sig_loci=S_T, sig_known=S_K,
        pct_known=round(100 * S_K / S_T, 1), fold=round((S_K / S_T) / (BG_K / BG_T), 2),
        fisher_p=fisher_greater(S_K, S_T - S_K, BG_K - S_K,
                                (BG_T - BG_K) - (S_T - S_K)))

    # --- HCC x Soskic (85c) 与 HCC x eQTLGen (94f)，both power levels
    g1 = pd.read_csv(f"{MR}/85c_G1_enrichment.tsv", sep="\t")
    for _, r in g1.iterrows():
        tag = r.outcome
        add("Soskic_CD4", tag,
            "main" if tag == "HCC_high" else "sensitivity",
            "85c_G1_enrichment.tsv",
            "" if tag == "HCC_high" else
            "power sensitivity: the two HCC levels come from one resource "
            "and are not independent (S20 §9)",
            bg_loci=int(r.bg_loci), bg_known=int(r.bg_known),
            sig_loci=int(r.sig_loci), sig_known=int(r.sig_known),
            pct_known=round(100 * r.sig_known / r.sig_loci, 1),
            fold=round(r.fold, 2), fisher_p=r.fisher_p)

    f94 = pd.read_csv(f"{MR}/94f_eqtlgen_hcc.tsv", sep="\t")
    for _, r in f94.iterrows():
        tag = r.disease
        add("eQTLGen_blood", tag,
            "main" if tag == "HCC_high" else "sensitivity",
            "94f_eqtlgen_hcc.tsv",
            "" if tag == "HCC_high" else
            "power sensitivity: the two HCC levels come from one resource "
            "and are not independent (S20 §9)",
            bg_loci=int(r.bg_loci), bg_known=int(r.bg_known),
            sig_loci=int(r.sig_loci), sig_known=int(r.sig_known),
            pct_known=round(100 * r.sig_known / r.sig_loci, 1),
            fold=round(r.fold, 2), fisher_p=r.fisher_p)

    # --- RA x both resources (108a)
    ra = pd.read_csv(f"{MR}/108a_ra_attribution.tsv", sep="\t")
    for cell, expo in (("N1 Soskic x RA", "Soskic_CD4"),
                       ("N2 eQTLGen x RA", "eQTLGen_blood")):
        r = ra[ra.cell == cell].iloc[0]
        add(expo, "RA", "main", "108a_ra_attribution.tsv", "",
            bg_loci=int(r.bg_loci), bg_known=int(r.bg_known),
            sig_loci=int(r.sig_loci), sig_known=int(r.sig_known),
            pct_known=round(r.pct_known, 1), fold=round(r.fold, 2),
            fisher_p=r.fisher_p)

    # --- 作废行：原样带过来，但标明作废与理由
    void = pd.read_csv(f"{MR}/94d_grid.tsv", sep="\t")
    for _, r in void[void.disease.isin(["lung", "colorectal"])].iterrows():
        add("eQTLGen_blood", r.disease, "void_bad_coords", "94d_grid.tsv",
            "VOID -- PREREG_generality_grid.md §8.2 偏离1: known-locus rsIDs were "
            "looked up in an outcome extract already filtered to instrument "
            "positions, leaving a background of only 2/559 and 5/559 known loci "
            "against 58/554 for melanoma. MUST NOT BE CITED.",
            bg_loci=int(r.bg_loci), bg_known=int(r.bg_known),
            sig_loci=int(r.sig_loci), sig_known=int(r.sig_known),
            pct_known=r.pct_known, fold=r.fold, fisher_p=r.fisher_p)

    return pd.DataFrame(rows)


# ------------------------------------------------------- mismatched controls
def build_mismatch():
    """预注册 §8 记了这两个数，但没有 TSV、也没有能重跑的脚本。此处重算。"""
    is_known_hcc, n_hcc = known_lookup("84a_hcc_known_loci_grch38.csv")
    is_known_mel, n_mel = known_lookup("landi2020_known_loci_grch38.csv")
    print(f"  known-locus lists: melanoma {n_mel}, HCC {n_hcc}")

    rows = []

    # NC1: eQTLGen x melanoma 的格子，改用 HCC 名单打分
    eq = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    s = score(eq[["chr", "pos", "fdr"]], is_known_hcc,
              "NC1 eQTLGen x melanoma scored with HCC list")
    rows.append(dict(control="NC1", cell="eQTLGen_blood x melanoma",
                     scored_with="HCC known loci", matched_fold=4.44,
                     matched_p=3.724056802770761e-11,
                     prereg_fold=1.30, prereg_p=0.411, **s))

    # NC2: Soskic x HCC-low 的格子，改用 melanoma 名单打分
    lo = pd.read_csv(f"{MR}/85a_HCC_low_annotated.tsv", sep="\t")
    s = score(lo[["chr", "pos", "fdr"]], is_known_mel,
              "NC2 Soskic x HCC-low scored with melanoma list")
    rows.append(dict(control="NC2", cell="Soskic_CD4 x HCC_low",
                     scored_with="melanoma known loci", matched_fold=17.62,
                     matched_p=0.0031240945794450406,
                     prereg_fold=0.00, prereg_p=1.0, **s))

    d = pd.DataFrame(rows)
    d["agrees_with_prereg"] = [
        (abs(f - pf) <= 0.02 and abs(p - pp) <= 0.005)
        for f, pf, p, pp in zip(d.fold.fillna(0), d.prereg_fold,
                                d.fisher_p, d.prereg_p)]
    return d


def main():
    print("=" * 78)
    print("Step 119 -- rebuilding the grid from each cell's own authoritative source")
    print("=" * 78)
    g = build_main()
    g.to_csv(f"{MR}/119a_grid_main.tsv", sep="\t", index=False)

    main_cells = g[g.status == "main"]
    n_gt1 = int((main_cells.fold > 1).sum())
    n_sig = int((main_cells.fisher_p < 0.05).sum())
    summary = pd.DataFrame([dict(
        n_main_cells=len(main_cells), n_fold_gt1=n_gt1,
        pct_fold_gt1=round(100 * n_gt1 / len(main_cells), 1),
        n_p_lt_05=n_sig,
        n_sensitivity=int((g.status == "sensitivity").sum()),
        n_void=int((g.status == "void_bad_coords").sum()),
        diseases="melanoma;HCC_high;RA", resources="Soskic_CD4;eQTLGen_blood")])
    summary.to_csv(f"{MR}/119b_grid_summary.tsv", sep="\t", index=False)

    print()
    print(f"  main grid: {len(main_cells)} cells, fold>1 in {n_gt1}, "
          f"P<0.05 in {n_sig}")
    print(f"  (94e_grid_summary.tsv said 6 testable / 5 fold>1 / 83.3% -- that was "
          f"the OLD grid,\n   which included the two voided cells and excluded RA)")

    print()
    print("-" * 78)
    print("mismatched-locus controls (no TSV existed; recomputed here)")
    print("-" * 78)
    mm = build_mismatch()
    mm.to_csv(f"{MR}/119c_mismatch_controls.tsv", sep="\t", index=False)
    print()
    for _, r in mm.iterrows():
        mark = "AGREES" if r.agrees_with_prereg else "*** DIFFERS ***"
        print(f"  {r.control} {r.cell:<28} recomputed {r.fold}x P={r.fisher_p:.4g}"
              f"  |  prereg {r.prereg_fold}x P={r.prereg_p}  ->  {mark}")

    print()
    print("wrote 119a_grid_main.tsv, 119b_grid_summary.tsv, 119c_mismatch_controls.tsv")


if __name__ == "__main__":
    main()
