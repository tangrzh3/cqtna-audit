"""复核脚本：口径更正到底改了哪些位点，逐个摊开。

    python review_convention.py                # 全部四格
    python review_convention.py eQTLGen        # 只看某一格

对每一个在两种口径下状态不同的位点，打印：
  · 位点坐标与跨度
  · 该位点上的全部工具变量，各自到最近已知 lead SNP 的距离，以及是否显著
  · 因此在"全部记录"口径下 known、在"仅显著记录"口径下 novel

这样复核的不是我的结论，而是**每一个具体位点该不该算 known** 这件事本身。
"""
import sys
import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

MR = r"D:/R_ex/MR"
LOCUS_KB = KNOWN_KB = 1000


def assign_loci(df):
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
    out["melanoma x Soskic_CD4"] = (d[["chr", "pos", "fdr", "SYMBOL"]], mel)
    d = pd.read_csv(f"{MR}/92c_locus_annotated.tsv", sep="\t")
    d = d.rename(columns={"symbol": "SYMBOL"})
    out["melanoma x eQTLGen_blood"] = (d[["chr", "pos", "fdr", "SYMBOL"]], mel)
    for tag in ("high", "low"):
        d = pd.read_csv(f"{MR}/85a_HCC_{tag}_annotated.tsv", sep="\t")
        d["SYMBOL"] = d.gene_id
        out[f"HCC_{tag} x Soskic_CD4"] = (d[["chr", "pos", "fdr", "SYMBOL"]], hcc)
    for _, (d, _) in out.items():
        d["chr"] = d["chr"].astype(str)
        d["pos"] = d["pos"].astype(int)
        d["fdr"] = pd.to_numeric(d["fdr"], errors="coerce")
    return out


def main():
    want = sys.argv[1] if len(sys.argv) > 1 else None
    for name, (d, kmap) in cells().items():
        if want and want.lower() not in name.lower():
            continue
        d = d.copy()
        snp = d.groupby(["chr", "pos"], as_index=False).size()
        loci = assign_loci(snp)
        d["locus"] = [loci[(c, p)] for c, p in zip(d.chr, d.pos)]
        d["dist"] = [nearest(kmap, c, p) for c, p in zip(d.chr, d.pos)]
        d["near"] = d.dist <= KNOWN_KB * 1000
        d["sig"] = d.fdr < 0.05

        sig_loci = d[d.sig].locus.unique()
        by_all = d.groupby("locus").near.any()
        by_sig = d[d.sig].groupby("locus").near.any()
        moved = [L for L in sig_loci if bool(by_all[L]) != bool(by_sig[L])]

        print("=" * 78)
        print(f"{name}")
        print(f"  significant loci: {len(sig_loci)}")
        print(f"  known under ALL-records  (adopted) : {int(by_all[sig_loci].sum())}")
        print(f"  known under SIG-records  (previous): {int(by_sig.sum())}")
        print(f"  loci that change status: {len(moved)}")
        for L in moved:
            sub = d[d.locus == L].sort_values("pos")
            span = (sub.pos.max() - sub.pos.min()) / 1000
            print(f"\n  --- locus {L}  chr{sub.chr.iloc[0]}:"
                  f"{sub.pos.min():,}-{sub.pos.max():,}  span {span:,.0f} kb, "
                  f"{len(sub)} records")
            print(f"      {'pos':>12} {'gene':<18} {'dist to known':>14} "
                  f"{'<=1Mb':>6} {'FDR<0.05':>9}")
            for _, r in sub.iterrows():
                dd = "inf" if not np.isfinite(r.dist) else f"{r.dist/1000:,.0f} kb"
                print(f"      {r.pos:>12,} {str(r.SYMBOL)[:18]:<18} {dd:>14} "
                      f"{'yes' if r.near else '-':>6} {'YES' if r.sig else '-':>9}")
            print("      -> ALL-records says KNOWN; SIG-records says NOVEL."
                  "\n         判断题：这个位点算不算'已知位点'？")
        print()


if __name__ == "__main__":
    main()
