#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 60a  为分析 B 抽取目标区域的结局侧 sumstats（FinnGen R12 与 meta 各一次单遍扫描）

设计见 manuscript/PREREG_power_trajectory.md §3。
主检验：结局侧在该区域有几个独立信号（暴露侧两轮完全相同，故变化必来自结局侧）。

*** 我们自己流程的阳性对照（比原登记的更强）***
MC1R 区已有 FinnGen 官方用**样本内 LD** 做出的答案：3 个高纯度 credible set
（cs1/cs2/cs3，log10BF 36.6/45.6/19.5）。我们用 1000G EUR 代理 LD 跑 susie_rss，
必须在同一区域复现 ≥2 个 credible set。做不到 = 代理 LD 不够用 =
按预注册把该分析记为无信息，而不是当作"未推翻结论"。

输出：regions/{region}_{finngen,meta}.tsv
"""
import gzip
import os

MR = r"D:/R_ex/MR"
OUT = os.path.join(MR, "regions")
os.makedirs(OUT, exist_ok=True)

HALF = 500_000
REGIONS = {
    # name: (chrom, 中心位置)  中心 = eQTL 峰（PARP1/ZFYVE19）或已知因果基因（MC1R）
    "PARP1":   ("1", 226_349_687),     # ④b 假阳性案例：PP.H4 0.92 -> 0.05
    "ZFYVE19": ("15", 40_832_560),     # ④b 假阴性案例：-> 0.99
    "MC1R":    ("16", 89_919_709),     # 阳性对照：官方样本内 LD 给出 3 个 credible set
    "TPI1":    ("12", 6_867_132),      # 应用示例位点
}

SOURCES = {
    "finngen": (os.path.join(MR, "finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz"),
                dict(chrom="#chrom", pos="pos", ref="ref", alt="alt",
                     beta="beta", se="sebeta", pval="pval", af="af_alt")),
    "meta":    (os.path.join(MR, "meta_melanoma_final.tsv.gz"),
                dict(chrom="#chrom", pos="pos", ref="ref", alt="alt",
                     beta="beta", se="sebeta", pval="pval", af="af_alt")),
}


def in_any(c, p):
    for name, (ch, mid) in REGIONS.items():
        if c == ch and abs(p - mid) <= HALF:
            return name
    return None


for src, (path, cols) in SOURCES.items():
    if not os.path.exists(path):
        print(f"[缺] {path}")
        continue
    buf = {k: [] for k in REGIONS}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {}
        for k, want in cols.items():
            ix[k] = header.index(want) if want in header else None
        if ix["chrom"] is None or ix["pos"] is None:
            raise SystemExit(f"{src}: 表头不认识 -> {header[:10]}")
        n = 0
        for line in fh:
            n += 1
            f = line.rstrip("\n").split("\t")
            c = f[ix["chrom"]].replace("chr", "")
            try:
                p = int(f[ix["pos"]])
            except ValueError:
                continue
            r = in_any(c, p)
            if r is None:
                continue
            row = {k: (f[i] if i is not None and i < len(f) else "")
                   for k, i in ix.items()}
            buf[r].append(row)
    for r, rows in buf.items():
        if not rows:
            print(f"  {src} {r}: 0 行")
            continue
        fp = os.path.join(OUT, f"{r}_{src}.tsv")
        with open(fp, "w", encoding="utf-8") as fo:
            fo.write("chrom\tpos\tref\talt\tbeta\tse\tpval\taf\n")
            for x in rows:
                fo.write("\t".join(x[k] for k in
                                   ("chrom", "pos", "ref", "alt", "beta",
                                    "se", "pval", "af")) + "\n")
        print(f"  {src} {r}: {len(rows):,} 行 -> {fp}", flush=True)
    print(f"[done] {src}: 扫描 {n:,} 行", flush=True)
print("完成")
