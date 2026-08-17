"""Step 121 -- the 1 Mb window sensitivity that HANDOFF v5 and v6 both carried
as outstanding, plus the continuous distance distribution.

两个"1 Mb"其实是**两个不同的约定**，此前一直被一起叫作"1 Mb 窗口"，
分开跑才说得清哪个在承重：

  KNOWN_KB  某位点算不算"已知" = 该位点上任一工具变量到最近的已知 lead SNP 的距离
            是否 ≤ 阈值。**这一个直接决定归属倍数。**
  LOCUS_KB  独立位点的划法 = 染色体内工具变量位置的单连锁聚类间距。
            这一个决定分母有多少个位点，以及显著记录被收成几个位点。

输出：
  121a  KNOWN_KB ∈ {100, 250, 500, 1000} kb，LOCUS_KB 固定 1000 —— 倍数与 Fisher P
  121b  LOCUS_KB ∈ {100, 250, 500, 1000} kb，KNOWN_KB 固定 1000 —— 位点数与倍数
  121c  ★ **连续距离分布**：每个位点到最近已知 lead SNP 的实际距离（bp），
        显著位点与背景位点分开列。HANDOFF 要的"报连续距离分布而非二分类"就是这张。

⚠ **只跑了四格**：melanoma × {Soskic, eQTLGen} 与 HCC-{high,low} × Soskic。
   RA 两格跑不了：`step108` 的 Okada rsID→GRCh38 坐标是**边扫 sumstats 边解的、
   没有落盘**，且 RA 的 per-record 注释表也没存（`108a` 只有汇总行）。
   要补 RA 必须重跑 step108 那一遍全量 sumstats 扫描。已在结尾打印说明。

打分规则与 step85/step119 完全一致；1000 kb 那一行必须重现已发表的数字，
脚本会自动核对并在对不上时报错退出。
"""
import os
import sys
from math import lgamma, exp

import numpy as np
import pandas as pd

# Windows 控制台是 GBK，⚠ 之类字符会让 print 崩掉（HANDOFF 坑 #34）
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

MR = r"D:/R_ex/MR"
KB = [100, 250, 500, 1000]
PUBLISHED = {          # 1000 kb 下必须重现的倍数（保留两位小数）
    "melanoma x Soskic_CD4": 4.09,
    "melanoma x eQTLGen_blood": 4.44,
    "HCC_high x Soskic_CD4": 8.85,
    "HCC_low x Soskic_CD4": 17.62,
}


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    hi = min(r1, c1)
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, hi + 1)), 1.0)


def assign_loci(df, locus_kb):
    out = {}
    for ch, sub in df.groupby("chr"):
        sub = sub.sort_values("pos")
        lid, prev = 0, None
        for _, r in sub.iterrows():
            if prev is not None and r.pos - prev > locus_kb * 1000:
                lid += 1
            out[(ch, r.pos)] = f"{ch}_{lid}"
            prev = r.pos
    return out


def known_map(path):
    k = pd.read_csv(f"{MR}/{path}")
    k["chr"] = k["chr"].astype(str)
    return {c: np.sort(s.pos.values) for c, s in k.groupby("chr")}, len(k)


def nearest_bp(by_chr, ch, pos):
    """到最近已知 lead SNP 的距离（bp）；该染色体上没有已知位点则 inf。"""
    arr = by_chr.get(str(ch))
    if arr is None or not len(arr):
        return np.inf
    i = np.searchsorted(arr, pos)
    best = np.inf
    for j in (i - 1, i):
        if 0 <= j < len(arr):
            best = min(best, abs(int(arr[j]) - int(pos)))
    return best


# ----------------------------------------------------------------- 四个格子
def load_cells():
    cells = {}

    mel_map, n_mel = known_map("landi2020_known_loci_grch38.csv")
    hcc_map, n_hcc = known_map("84a_hcc_known_loci_grch38.csv")
    print(f"  known-locus lists: melanoma {n_mel}, HCC {n_hcc}")

    d = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
    d[["chr", "pos"]] = d.SNP.str.split(":", expand=True)
    d = d.rename(columns={"FDR": "fdr"})
    cells["melanoma x Soskic_CD4"] = (d[["chr", "pos", "fdr"]], mel_map)

    d = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    cells["melanoma x eQTLGen_blood"] = (d[["chr", "pos", "fdr"]], mel_map)

    for tag in ("high", "low"):
        d = pd.read_csv(f"{MR}/85a_HCC_{tag}_annotated.tsv", sep="\t")
        cells[f"HCC_{tag} x Soskic_CD4"] = (d[["chr", "pos", "fdr"]], hcc_map)

    for k, (d, _) in cells.items():
        d["chr"] = d["chr"].astype(str)
        d["pos"] = d["pos"].astype(int)
        d["fdr"] = pd.to_numeric(d["fdr"], errors="coerce")
    return cells


def score(d, kmap, locus_kb, known_kb):
    d = d.copy()
    snp = d.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp, locus_kb)
    d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
    d["dist"] = [nearest_bp(kmap, c, p) for c, p in zip(d.chr, d.pos)]
    d["known"] = d.dist <= known_kb * 1000

    bg = d.groupby("locus").agg(known=("known", "any"))
    sg = d[d.fdr < 0.05].groupby("locus").agg(known=("known", "any"))
    BT, BK, ST, SK = len(bg), int(bg.known.sum()), len(sg), int(sg.known.sum())
    fold = ((SK / ST) / (BK / BT)) if ST and BK else np.nan
    p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK))
    return dict(bg_loci=BT, bg_known=BK, sig_loci=ST, sig_known=SK,
                pct_known=round(100 * SK / ST, 1) if ST else np.nan,
                fold=round(fold, 2) if ST else np.nan, fisher_p=p), d


def main():
    print("=" * 84)
    print("Step 121 -- window sensitivity for the two 1 Mb conventions")
    print("=" * 84)
    cells = load_cells()

    # ---- 121a: KNOWN_KB 扫描（LOCUS_KB 固定 1000）
    rows, checked = [], {}
    print("\n--- 121a  KNOWN_KB sweep (locus definition held at 1 Mb) " + "-" * 22)
    print(f"{'cell':<28}{'known_kb':>9}{'sig known/loci':>16}{'bg known/loci':>15}"
          f"{'fold':>8}{'P':>11}")
    for name, (d, kmap) in cells.items():
        for kb in KB:
            s, _ = score(d, kmap, 1000, kb)
            rows.append(dict(cell=name, locus_kb=1000, known_kb=kb, **s))
            sig_txt = "{}/{}".format(s["sig_known"], s["sig_loci"])
            bg_txt = "{}/{}".format(s["bg_known"], s["bg_loci"])
            print(f"{name:<28}{kb:>9}{sig_txt:>16}{bg_txt:>15}"
                  f"{s['fold']:>8}{s['fisher_p']:>11.4g}")
            if kb == 1000:
                checked[name] = s["fold"]
        print()
    pd.DataFrame(rows).to_csv(f"{MR}/121a_window_sensitivity.tsv", sep="\t", index=False)

    bad = [(k, v, PUBLISHED[k]) for k, v in checked.items()
           if k in PUBLISHED and abs(v - PUBLISHED[k]) > 0.02]
    if bad:
        raise SystemExit(f"*** 1000 kb row does not reproduce published folds: {bad}")
    print("  1000 kb reproduces every published fold exactly.\n")

    # ---- 121b: LOCUS_KB 扫描（KNOWN_KB 固定 1000）
    rows = []
    print("--- 121b  LOCUS_KB sweep (known-locus window held at 1 Mb) " + "-" * 19)
    print(f"{'cell':<28}{'locus_kb':>9}{'bg loci':>9}{'sig loci':>10}{'fold':>8}{'P':>11}")
    for name, (d, kmap) in cells.items():
        for kb in KB:
            s, _ = score(d, kmap, kb, 1000)
            rows.append(dict(cell=name, locus_kb=kb, known_kb=1000, **s))
            print(f"{name:<28}{kb:>9}{s['bg_loci']:>9}{s['sig_loci']:>10}"
                  f"{s['fold']:>8}{s['fisher_p']:>11.4g}")
        print()
    pd.DataFrame(rows).to_csv(f"{MR}/121b_locus_definition_sensitivity.tsv",
                              sep="\t", index=False)

    # ---- 121c: 连续距离分布
    dist_rows = []
    print("--- 121c  continuous distance to the nearest known lead SNP " + "-" * 18)
    for name, (d, kmap) in cells.items():
        _, ann = score(d, kmap, 1000, 1000)

        # ⚠ 距离必须按**与二分类同一套记录**来取，否则两者会打架：
        # 二分类里"显著位点是否已知" = 该位点的 **FDR<0.05 记录** 里有没有落进窗口的，
        # 所以显著位点的距离也必须只在它自己的显著记录上取 min。
        # （初版对位点内全部记录取 min，于是 HCC-high 报出 235 kb，
        #   而那 235 kb 属于一条不显著的记录——它的显著记录其实在 408 kb 外，
        #   与 121a 在 250 kb 下判 0/2 直接矛盾。）
        per_locus = ann.groupby("locus").agg(dist_all=("dist", "min"),
                                             fdr=("fdr", "min"),
                                             chr=("chr", "first"))
        sig_rec = ann[ann.fdr < 0.05]
        dist_sig = sig_rec.groupby("locus").dist.min()
        pos_sig = sig_rec.groupby("locus").pos.first()
        per_locus["significant"] = per_locus.index.isin(dist_sig.index)
        per_locus["dist_sig"] = dist_sig
        per_locus["pos_sig"] = pos_sig

        for locus, r in per_locus.iterrows():
            dd = r.dist_sig if r.significant else r.dist_all
            dist_rows.append(dict(
                cell=name, locus=locus, chr=r["chr"],
                pos=int(r.pos_sig) if r.significant and np.isfinite(r.pos_sig) else None,
                significant=bool(r.significant),
                dist_bp=None if not np.isfinite(dd) else int(dd),
                dist_bp_any_record=(None if not np.isfinite(r.dist_all)
                                    else int(r.dist_all)),
                min_fdr=r.fdr))
        sig = per_locus.loc[per_locus.significant, "dist_sig"].replace(np.inf, np.nan).dropna()
        bg = per_locus.dist_all.replace(np.inf, np.nan).dropna()
        if len(sig):
            print(f"  {name}")
            print(f"    significant loci (n={len(sig)}): median "
                  f"{sig.median()/1000:,.0f} kb, min {sig.min()/1000:,.0f} kb, "
                  f"max {sig.max()/1000:,.0f} kb")
            print(f"      distances (kb): "
                  f"{', '.join(f'{v/1000:,.0f}' for v in sorted(sig))}")
            print(f"    all testable loci (n={len(bg)}): median "
                  f"{bg.median()/1000:,.0f} kb")
    pd.DataFrame(dist_rows).to_csv(f"{MR}/121c_distance_distribution.tsv",
                                   sep="\t", index=False)

    print("\n" + "=" * 84)
    print("wrote 121a_window_sensitivity.tsv, 121b_locus_definition_sensitivity.tsv,")
    print("      121c_distance_distribution.tsv")
    print("\n⚠ RA (both resources) is NOT included. step108 resolves the Okada rsIDs to")
    print("  GRCh38 while streaming the FinnGen RA sumstats and never saves the map, and")
    print("  108a holds only summary rows -- there is no per-record RA table to re-score.")
    print("  Adding RA means re-running that full sumstats pass.")


if __name__ == "__main__":
    main()
