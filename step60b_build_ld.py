#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 60b  为分析 B 的目标区域构建 LD 矩阵（1000G EUR 525 人，GRCh38 全面板）

⚠ 这是**代理 LD**：结局是 FinnGen（芬兰人群，有奠基者效应），参考面板是 1000G EUR。
susie_rss 对 LD 失配敏感，故：
  - 必须报告 estimate_s（失配诊断）
  - MC1R 区作为流程阳性对照（官方样本内 LD 已知有 3 个 credible set）
  - 对照不过则整段分析记为无信息（预注册 §3.3–3.4）

输出：regions/{region}_ld.vcor1 与 .vars
"""
import os
import subprocess
import sys

MR = r"D:/R_ex/MR"
REG = os.path.join(MR, "regions")
PLINK = os.path.join(MR, "bin", "plink2.exe")
PFILE = os.path.join(MR, "ref", "all_hg38")

HALF = 500_000
REGIONS = {
    "PARP1":   ("1", 226_349_687),
    "ZFYVE19": ("15", 40_832_560),
    "MC1R":    ("16", 89_919_709),
    "TPI1":    ("12", 6_867_132),
}

# EUR 去亲缘 525 人：沿用 SMR 用的同一批样本
keep = os.path.join(REG, "eur.keep")
with open(os.path.join(MR, "ref", "1kg_eur.fam")) as fh, \
     open(keep, "w") as fo:
    n = 0
    for line in fh:
        f = line.split()
        fo.write(f"{f[0]}\t{f[1]}\n")
        n += 1
print(f"EUR 样本: {n}")


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print("  [plink2 失败]", " ".join(cmd[-6:]))
        print(p.stdout[-1500:])
        print(p.stderr[-1500:])
        return False
    return True


for name, (ch, mid) in REGIONS.items():
    lo, hi = mid - HALF, mid + HALF
    tmp = os.path.join(REG, f"{name}_tmp")
    out = os.path.join(REG, f"{name}_ld")
    if os.path.exists(out + ".unphased.vcor1"):
        print(f"[skip] {name}")
        continue
    print(f"[{name}] chr{ch}:{lo}-{hi}", flush=True)
    ok = run([PLINK, "--pfile", PFILE, "vzs", "--keep", keep,
              "--chr", ch, "--from-bp", str(lo), "--to-bp", str(hi),
              "--snps-only", "--max-alleles", "2", "--maf", "0.01",
              "--set-all-var-ids", "@:#:$r:$a",
              "--new-id-max-allele-len", "60", "missing",
              "--rm-dup", "exclude-all",
              "--make-pgen", "--out", tmp])
    if not ok:
        continue
    ok = run([PLINK, "--pfile", tmp, "--r-unphased", "square", "ref-based",
              "--out", out])
    if not ok:
        continue
    nvar = sum(1 for _ in open(out + ".unphased.vcor1.vars"))
    print(f"  -> {nvar:,} 个变异的 LD 矩阵", flush=True)

print("完成")
