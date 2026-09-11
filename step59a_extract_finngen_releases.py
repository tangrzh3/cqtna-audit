#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 59a  流式下载 FinnGen R8–R11 的黑色素瘤 sumstats，只保留工具变量所在的行。

设计与判据：manuscript/PREREG_power_trajectory.md
预测已登记：58a_finngen_reference_predictions.tsv（在本脚本运行前完成）

每个 endpoint 约 0.8 GB，边下边过滤，只写出约 2,100 行。已存在的 extract 跳过（可续跑）。

⚠ 路径命名跨代不一致（第 0 步已核实）：
    R8/R9/R10/R11 : summary_stats/finngen_R{N}_C3_MELANOMA_SKIN_EXALLC.gz
    R12           : summary_stats/release/finngen_R12_...（本项目已在本地）
"""
import sys
import csv
import gzip
import io
import os
import time
import urllib.request

MR = (sys.argv[1] if len(sys.argv) > 1
      else os.environ.get("CQTNA_DIR") or r"D:/R_ex/MR")
OUT = os.path.join(MR, "release_extracts")
os.makedirs(OUT, exist_ok=True)

ENDPOINT = "C3_MELANOMA_SKIN_EXALLC"
RELEASES = ["R8", "R9", "R10", "R11"]


def url_for(rel):
    n = rel[1:]
    return (f"https://storage.googleapis.com/finngen-public-data-r{n}/"
            f"summary_stats/finngen_{rel}_{ENDPOINT}.gz")


def idx(header, *names):
    for nm in names:
        if nm in header:
            return header.index(nm)
    raise KeyError(f"none of {names} in header: {header[:12]}")


# ---------------------------------------------------------------- 工具变量
need = set()
with open(os.path.join(MR, "06_locus_annotation.tsv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        need.add(r["SNP"])                       # "chr:pos"
print(f"工具变量 SNP 数: {len(need):,}", flush=True)

# ---------------------------------------------------------------- 逐 release
for rel in RELEASES:
    dest = os.path.join(OUT, f"{rel}_{ENDPOINT}.tsv")
    if os.path.exists(dest):
        print(f"[skip] {rel} 已提取", flush=True)
        continue
    u = url_for(rel)
    t0, n, got = time.time(), 0, []
    print(f"[start] {rel}  {u}", flush=True)
    try:
        req = urllib.request.Request(u, headers={"User-Agent": "curl/8"})
        with urllib.request.urlopen(req, timeout=180) as resp:
            with gzip.GzipFile(fileobj=resp) as gz:
                txt = io.TextIOWrapper(gz, encoding="utf-8")
                header = txt.readline().rstrip("\n").split("\t")
                ic = idx(header, "#chrom", "chrom", "#CHROM")
                ip = idx(header, "pos", "POS")
                ir = idx(header, "ref", "REF")
                ia = idx(header, "alt", "ALT")
                ib = idx(header, "beta")
                ise = idx(header, "sebeta")
                ipv = idx(header, "pval")
                iaf = None
                for cand in ("af_alt", "maf", "af_alt_cases"):
                    if cand in header:
                        iaf = header.index(cand)
                        break
                for line in txt:
                    n += 1
                    f = line.rstrip("\n").split("\t")
                    k = f[ic] + ":" + f[ip]
                    if k in need:
                        got.append(dict(SNP=k, ref=f[ir], alt=f[ia],
                                        beta=f[ib], sebeta=f[ise], pval=f[ipv],
                                        af=(f[iaf] if iaf is not None else "")))
                    if n % 5_000_000 == 0:
                        print(f"   {rel} {n/1e6:.0f}M 行, 命中 {len(got):,}, "
                              f"{time.time()-t0:.0f}s", flush=True)
    except Exception as exc:
        print(f"[FAIL] {rel}: {type(exc).__name__}: {exc}", flush=True)
        continue

    with open(dest, "w", newline="", encoding="utf-8") as fo:
        w = csv.DictWriter(fo, fieldnames=["SNP", "ref", "alt", "beta",
                                           "sebeta", "pval", "af"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(got)
    print(f"[done] {rel}: {n:,} 行扫描, {len(got):,} 命中, "
          f"{time.time()-t0:.0f}s -> {dest}", flush=True)

print("全部完成", flush=True)
