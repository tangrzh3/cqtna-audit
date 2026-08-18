"""复核脚本：位点分区与归属口径。一条命令复现评审包 §5b 的那张表。

    python review_convention.py                 # 四窗口 × 三口径 × 匹配/错配 全表
    python review_convention.py --loci eQTLGen  # 逐位点摊开某一格的差异
    python review_convention.py --spans         # 各格的位点跨度分布

⚠ 早先版本把 LOCUS_KB 写死为 1000，只能出单窗口的结果，
   而评审包却引用它来复现四窗口扫描——脚本与文档对不上。此处修好：
   `--table` 扫的是 **locus window**（分区），不是 known window（归属窗口）。
   两者都叫"1 Mb"，但改的是不同的东西，`step122` 扫的是后者。
"""
import argparse
import sys
from math import lgamma, exp

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

MR = r"D:/R_ex/MR"
KNOWN_KB = 1000                 # 归属窗口，全程固定
LOCUS_WINDOWS = (100, 250, 500, 1000)


def fisher_greater(a, b, c, d):
    def logC(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    den = logC(n, c1)
    return min(sum(exp(logC(r1, x) + logC(n - r1, c1 - x) - den)
                   for x in range(a, min(r1, c1) + 1)), 1.0)


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
    return {c: np.sort(s.pos.values) for c, s in k.groupby("chr")}


def nearest(kmap, ch, pos):
    arr = kmap.get(str(ch))
    if arr is None or not len(arr):
        return np.inf
    i = np.searchsorted(arr, pos)
    return min((abs(int(arr[j]) - int(pos))
                for j in (i - 1, i) if 0 <= j < len(arr)), default=np.inf)


def cells():
    mel = known_map("landi2020_known_loci_grch38.csv")
    hcc = known_map("84a_hcc_known_loci_grch38.csv")
    out = {}
    d = pd.read_csv(f"{MR}/13_meta_locus_annotation.tsv", sep="\t")
    d[["chr", "pos"]] = d.SNP.str.split(":", expand=True)
    d = d.rename(columns={"FDR": "fdr"})
    out["melanoma x Soskic_CD4"] = (d[["chr", "pos", "fdr", "SYMBOL"]], mel, hcc)
    d = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    d = d.rename(columns={"symbol": "SYMBOL"})
    out["melanoma x eQTLGen_blood"] = (d[["chr", "pos", "fdr", "SYMBOL"]], mel, hcc)
    for tag in ("high", "low"):
        d = pd.read_csv(f"{MR}/85a_HCC_{tag}_annotated.tsv", sep="\t")
        d["SYMBOL"] = d.gene_id
        out[f"HCC_{tag} x Soskic_CD4"] = (d[["chr", "pos", "fdr", "SYMBOL"]], hcc, mel)
    for _, (d, _, _) in out.items():
        d["chr"] = d["chr"].astype(str)
        d["pos"] = d["pos"].astype(int)
        d["fdr"] = pd.to_numeric(d["fdr"], errors="coerce")
    return out


def score(d, kmap, locus_kb):
    """三种口径下的 (倍数, P)，外加位点跨度。"""
    d = d.copy()
    snp = d.groupby(["chr", "pos"], as_index=False).size()
    loci = assign_loci(snp, locus_kb)
    d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
    d["near"] = [nearest(kmap, c, p) <= KNOWN_KB * 1000
                 for c, p in zip(d.chr, d.pos)]
    d["sig"] = d.fdr < 0.05
    sig_loci = d[d.sig].locus.unique()

    by_all = d.groupby("locus").near.any()
    by_sig = d[d.sig].groupby("locus").near.any()
    lead = d.loc[d.groupby("locus").fdr.idxmin()]

    def stat(sk, st, bk, bt):
        fold = (sk / st) / (bk / bt) if st and bk else float("nan")
        return fold, fisher_greater(sk, st - sk, bk - sk, (bt - bk) - (st - sk))

    res = {
        "C1": stat(int(by_all[sig_loci].sum()), len(sig_loci),
                   int(by_all.sum()), len(by_all)),
        "C2": stat(int(by_sig.sum()), len(by_sig),
                   int(by_all.sum()), len(by_all)),
        "C3": stat(int(lead[lead.locus.isin(sig_loci)].near.sum()), len(sig_loci),
                   int(lead.near.sum()), len(lead)),
    }
    span = d.groupby("locus").pos.agg(lambda v: (v.max() - v.min()) / 1000)
    return res, span, sig_loci


def cmd_table():
    """评审包 §5b 的表：locus window × 三口径 × 匹配/错配。"""
    for name, (d, matched, mismatched) in cells().items():
        print("=" * 108)
        print(name)
        print(f"{'locus_kb':>9}{'max span kb':>13}{'sig>5x win':>11}   "
              f"{'C1 匹配/错配':>24}{'C2 匹配/错配':>24}{'C3 匹配/错配':>24}")
        print("-" * 108)
        for kb in LOCUS_WINDOWS:
            m, span, sig = score(d, matched, kb)
            x, _, _ = score(d, mismatched, kb)
            wide = int((span[sig] > 5 * kb).sum())
            row = ""
            for c in ("C1", "C2", "C3"):
                row += f"  {m[c][0]:5.2f}/{x[c][0]:<5.2f} P={x[c][1]:<7.3g}"
            print(f"{kb:>9}{span[sig].max():>13,.0f}{wide:>11}{row}")
        print()
    print("错配那一列的 P 是判据：< 0.05 即该格作废（PREREG_generality_grid.md §5/§6-E）。")
    print("⚠ 'sig>5x win' 是跨度超过聚类窗口 5 倍的显著位点数——")
    print("   它在任何窗口下都不为零，所以没有哪个窗口可以被称作'不串联'。")


def cmd_spans():
    for name, (d, kmap, _) in cells().items():
        print("=" * 78)
        print(name)
        for kb in LOCUS_WINDOWS:
            _, span, sig = score(d, kmap, kb)
            print(f"  locus_kb {kb:>4}: {len(span):>5} 位点  "
                  f"跨度中位 {span.median():>7,.0f} kb  最大 {span.max():>9,.0f} kb  "
                  f"显著位点最大 {span[sig].max():>8,.0f} kb  "
                  f"超窗口5倍的显著位点 {int((span[sig] > 5 * kb).sum())}")
        print()


def cmd_loci(want, locus_kb):
    for name, (d, kmap, _) in cells().items():
        if want and want.lower() not in name.lower():
            continue
        d = d.copy()
        snp = d.groupby(["chr", "pos"], as_index=False).size()
        loci = assign_loci(snp, locus_kb)
        d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
        d["dist"] = [nearest(kmap, c, p) for c, p in zip(d.chr, d.pos)]
        d["near"] = d.dist <= KNOWN_KB * 1000
        d["sig"] = d.fdr < 0.05
        sig_loci = d[d.sig].locus.unique()
        by_all = d.groupby("locus").near.any()
        by_sig = d[d.sig].groupby("locus").near.any()
        moved = [L for L in sig_loci if bool(by_all[L]) != bool(by_sig[L])]

        print("=" * 78)
        print(f"{name}   (locus_kb = {locus_kb})")
        print(f"  显著位点 {len(sig_loci)}；C1 判已知 {int(by_all[sig_loci].sum())}，"
              f"C2 判已知 {int(by_sig.sum())}；状态不同的位点 {len(moved)}")
        for L in moved:
            sub = d[d.locus == L].sort_values("pos")
            span = (sub.pos.max() - sub.pos.min()) / 1000
            print(f"\n  --- {L}  chr{sub.chr.iloc[0]}:{sub.pos.min():,}-{sub.pos.max():,}"
                  f"  跨度 {span:,.0f} kb，{len(sub)} 条记录")
            print(f"      {'pos':>12} {'gene':<18} {'到最近已知':>12} {'≤1Mb':>6} {'显著':>6}")
            for _, r in sub.iterrows():
                dd = "inf" if not np.isfinite(r.dist) else f"{r.dist/1000:,.0f} kb"
                print(f"      {r.pos:>12,} {str(r.SYMBOL)[:18]:<18} {dd:>12} "
                      f"{'yes' if r.near else '-':>6} {'YES' if r.sig else '-':>6}")
            print("      → C1 判 KNOWN，C2 判 NOVEL。这个位点该不该算'已知位点'？")
        print()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", metavar="CELL",
                    help="逐位点摊开某一格（模糊匹配格名，如 eQTLGen）")
    ap.add_argument("--spans", action="store_true", help="只看位点跨度分布")
    ap.add_argument("--locus-kb", type=int, default=1000,
                    help="--loci 用的位点窗口（默认 1000）")
    a = ap.parse_args()
    if a.spans:
        cmd_spans()
    elif a.loci:
        cmd_loci(a.loci, a.locus_kb)
    else:
        cmd_table()


if __name__ == "__main__":
    main()
