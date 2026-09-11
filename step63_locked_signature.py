#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 63  锁定三套数据共用的糖酵解 signature，并做 leave-TPI1-out 验证

起因：外部审稿指出 eQTL 侧用 24-28 个基因、multiome 轴用 16 个、患者分析用 22 个，
"三个相关但不完全相同的 glycolytic programmes"。

核查结果：三套集合是**嵌套**的 —— multiome(16) ⊂ 患者(22) ⊂ eQTL(28)，
交集恰为 multiome 的 16 个基因。故锁定集 = 该 16 个基因。

本脚本在锁定集上重跑两项关键结果，并各做一次 leave-TPI1-out：
  (A) eQTL 侧：最强 eQTL 所在时点是否回避静息态（类水平 timing）
  (B) 患者侧：ICB 应答方向一致性（二项检验）
同时输出**测量精度对照**（lead 变异 SE 中位数），因为 Step 62 已证明该结果受
差异化可测性混杂，必须与结果同行。

输出：63a_locked_signature.tsv、63b_timing_locked.tsv、63c_patient_locked.tsv
"""
import sys
import csv
import os
import numpy as np
import pyarrow.parquet as pq
from math import comb

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
P = r"D:/Downloads/CD4_eqtl_step1_clean"

PROFILES = [("CD4_Naive_uns_0h", "Naive", "0h"), ("CD4_Naive_stim_16h", "Naive", "16h"),
            ("CD4_Naive_stim_40h", "Naive", "40h"), ("CD4_Naive_stim_5d", "Naive", "5d"),
            ("CD4_Memory_uns_0h", "Memory", "0h"), ("CD4_Memory_stim_16h", "Memory", "16h"),
            ("CD4_Memory_stim_40h", "Memory", "40h"), ("CD4_Memory_stim_5d", "Memory", "5d")]

MULTIOME = ["SLC2A1", "SLC2A3", "HK1", "HK2", "GPI", "PFKL", "PFKP", "PFKFB3",
            "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA"]


def binom_p(k, n, p=0.5):
    """单侧二项检验（>=k）"""
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


# ---------------------------------------------------------------- 集合关系
eqtl_map = {r["SYMBOL"]: r["ENSEMBL"] for r in
            csv.DictReader(open(f"{MR}/35a_glyco_ensembl_map.tsv", encoding="utf-8"),
                           delimiter="\t")}
patient = sorted({r["gene"] for r in
                  csv.DictReader(open(f"{MR}/32e_glyco_gene_by_gene_corrected.tsv",
                                      encoding="utf-8"), delimiter="\t")})
locked = sorted(set(MULTIOME) & set(patient) & set(eqtl_map))
print(f"eQTL 侧 {len(eqtl_map)} | multiome {len(MULTIOME)} | 患者 {len(patient)}")
print(f"嵌套关系: multiome ⊂ 患者 ⊂ eQTL 侧 —— "
      f"{set(MULTIOME) <= set(patient) <= set(eqtl_map)}")
print(f"★ 锁定集 = {len(locked)} 个基因: {locked}\n")
with open(f"{MR}/63a_locked_signature.tsv", "w", newline="", encoding="utf-8") as fo:
    w = csv.writer(fo, delimiter="\t")
    w.writerow(["SYMBOL", "ENSEMBL", "in_eqtl", "in_multiome", "in_patient", "locked"])
    for s in sorted(set(eqtl_map)):
        w.writerow([s, eqtl_map[s], 1, int(s in MULTIOME), int(s in patient),
                    int(s in locked)])

locked_ens = {eqtl_map[s] for s in locked}
locked_ens_noTPI1 = {eqtl_map[s] for s in locked if s != "TPI1"}

# ---------------------------------------------------------------- (A) timing
best = {}          # gene -> lineage -> (p, timepoint, se)
for prof, lin, tp in PROFILES:
    t = pq.read_table(os.path.join(P, prof + "_step1_clean.parquet"),
                      columns=["gene_id", "pval", "se"]).to_pydict()
    cur = {}
    for g, p, s in zip(t["gene_id"], t["pval"], t["se"]):
        if g not in cur or p < cur[g][0]:
            cur[g] = (p, s)
    for g, (p, s) in cur.items():
        best.setdefault(g, {}).setdefault(lin, []).append((p, tp, s))

rows = []
for label, gset in [("locked16", locked_ens), ("locked15_noTPI1", locked_ens_noTPI1)]:
    n0h_g = n_g = n0h_b = n_b = 0
    for g, per_lin in best.items():
        for lin, lst in per_lin.items():
            if len(lst) < 4:          # 需四个时点齐全才可比
                continue
            tp = min(lst)[1]
            if g in gset:
                n_g += 1; n0h_g += (tp == "0h")
            elif g not in locked_ens:  # 背景不含任一糖酵解基因
                n_b += 1; n0h_b += (tp == "0h")
    fg, fb = n0h_g / n_g, n0h_b / n_b
    orat = (fg / (1 - fg)) / (fb / (1 - fb))
    # Fisher 精确检验（单侧，糖酵解在 0h 更少）
    from scipy.stats import fisher_exact
    _, pv = fisher_exact([[n0h_g, n_g - n0h_g], [n0h_b, n_b - n0h_b]],
                         alternative="less")
    print(f"[timing] {label:16} 糖酵解 {n0h_g}/{n_g} ({fg:.1%}) vs "
          f"背景 {n0h_b}/{n_b} ({fb:.1%})  OR={orat:.3f}  p={pv:.4f}")
    rows.append(dict(set=label, n_gene_lineage=n_g, n_0h=n0h_g, frac=fg,
                     bg_n=n_b, bg_0h=n0h_b, bg_frac=fb, OR=orat, p=pv))

# 测量精度对照：锁定集 vs 背景的 lead-SE 中位数比值（必须与上表同行报告）
prec = []
for prof, lin, tp in PROFILES:
    t = pq.read_table(os.path.join(P, prof + "_step1_clean.parquet"),
                      columns=["gene_id", "pval", "se"]).to_pydict()
    cur = {}
    for g, p, s in zip(t["gene_id"], t["pval"], t["se"]):
        if g not in cur or p < cur[g][0]:
            cur[g] = (p, s)
    gl = [v[1] for g, v in cur.items() if g in locked_ens]
    bg = [v[1] for g, v in cur.items() if g not in locked_ens]
    prec.append(dict(lineage=lin, timepoint=tp, se_locked=np.median(gl),
                     se_background=np.median(bg),
                     ratio=np.median(gl) / np.median(bg), n_locked=len(gl)))
    print(f"[precision] {lin:7}{tp:4} SE 中位 锁定={np.median(gl):.4f} "
          f"背景={np.median(bg):.4f} 比值={np.median(gl)/np.median(bg):.2f}")

with open(f"{MR}/63b_timing_locked.tsv", "w", newline="", encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
with open(f"{MR}/63b_precision_by_profile.tsv", "w", newline="", encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(prec[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(prec)

# ---------------------------------------------------------------- (B) 患者侧
pat = list(csv.DictReader(open(f"{MR}/32e_glyco_gene_by_gene_corrected.tsv",
                               encoding="utf-8"), delimiter="\t"))
out = []
print()
for tpt in ("Pre", "Post"):
    sub = [r for r in pat if r["timepoint"] == tpt]
    for label, gset in [("all22", set(patient)), ("locked16", set(locked)),
                        ("locked15_noTPI1", set(locked) - {"TPI1"})]:
        s = [r for r in sub if r["gene"] in gset]
        # diff = median_R - median_NR；非应答者更高 = diff < 0
        hi_nr = sum(1 for r in s if float(r["diff"]) < 0)
        n = len(s)
        p = binom_p(hi_nr, n)
        print(f"[patient] {tpt:5}{label:18} 非应答者更高 {hi_nr}/{n}  二项 p={p:.2e}")
        out.append(dict(timepoint=tpt, set=label, higher_in_NR=hi_nr, n=n, binom_p=p))

with open(f"{MR}/63c_patient_locked.tsv", "w", newline="", encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(out[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(out)
print("\n写出 63a / 63b / 63b_precision / 63c")
