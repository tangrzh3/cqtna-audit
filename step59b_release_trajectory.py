#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 59b  FinnGen release 功效轨迹：观测值 vs 已登记的预测

设计、判据、判读表：manuscript/PREREG_power_trajectory.md
预测（登记于任何 R8–R11 数据被查看之前）：58a_finngen_reference_predictions.tsv

固定不变（与主分析逐位一致）：暴露数据、工具变量、等位定向、Wald ratio、
BH-FDR 的检验家族、位点标注、1 Mb 位点合并规则。唯一系统变化的是结局侧病例数。

*** 三个过程阳性对照（不通过则不解释结果）***
 P1 重现校验：R12 那一轮必须重现 06_locus_annotation 的 FDR<0.05 = 10 条 / 6 基因
 P2 单调性：命中数与病例数的 Spearman rho > 0
 P3 已知大效应位点应在最早的 release 即出现（MC1R 簇 / PARP1）

输出：59a_release_trajectory.tsv、59b_release_lists.tsv、59c_prediction_check.tsv
"""
import csv
import math
import os
from statistics import NormalDist

import numpy as np
import pandas as pd

MR = r"D:/R_ex/MR"
EXT = os.path.join(MR, "release_extracts")
ND = NormalDist()
FDR_MAIN, FDR_STRAT = 0.05, 0.20
NOVEL = "潜在新位点"
ENDPOINT = "C3_MELANOMA_SKIN_EXALLC"

RELEASES = [("R8", 2705, 259583), ("R9", 2993, 287137),
            ("R10", 3194, 314193), ("R11", 3932, 345118),
            ("R12", 5753, 378749)]


def n_eff(nc, nk):
    return 4.0 / (1.0 / nc + 1.0 / nk)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


# ---------------------------------------------------------------- 暴露侧
d = pd.read_csv(f"{MR}/06_locus_annotation.tsv", sep="\t")
d = d[["exposure", "SYMBOL", "category", "SNP", "beta.exposure",
       "beta.outcome", "se.outcome"]].dropna()
d = d[(d["beta.exposure"] != 0) & (d["se.outcome"] > 0)].reset_index(drop=True)
d["novel"] = d["category"].eq(NOVEL)

# 位点 ID：与 step58 完全相同的规则
d[["_chr", "_pos"]] = d.SNP.str.split(":", expand=True).iloc[:, :2]
d["_pos"] = d["_pos"].astype(int)
d = d.sort_values(["_chr", "_pos"]).reset_index(drop=True)
loc, cur, lc, lp = [], 0, None, None
for c, p in zip(d._chr, d._pos):
    if c != lc or p - lp > 1_000_000:
        cur += 1
    loc.append(cur)
    lc, lp = c, p
d["locus"] = loc

# 等位定向：01_harmonised_all 中 effect_allele.outcome 已对齐到暴露效应等位
alle = {}
with open(f"{MR}/01_harmonised_all.tsv", encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        alle[(r["exposure"], r["SNP"])] = (r["effect_allele.outcome"].upper(),
                                           r["other_allele.outcome"].upper())
d["ea_out"] = [alle.get((e, s), ("", ""))[0] for e, s in zip(d.exposure, d.SNP)]
d["oa_out"] = [alle.get((e, s), ("", ""))[1] for e, s in zip(d.exposure, d.SNP)]
d = d[(d.ea_out != "") & (d.oa_out != "")].reset_index(drop=True)
print(f"暴露侧可用: {len(d):,} 条 / {d.SYMBOL.nunique():,} 基因 / "
      f"{d.locus.nunique():,} 位点")

# ---------------------------------------------------------------- 参照名单（R12）
p_r12 = np.array([2 * ND.cdf(-abs(bo / se)) if se > 0 else 1.0
                  for bo, se in zip(d["beta.outcome"] / d["beta.exposure"],
                                    d["se.outcome"] / d["beta.exposure"].abs())])
fdr_r12 = bh(p_r12)
ref_main = set(d.SYMBOL[fdr_r12 < FDR_MAIN])
rs = d[fdr_r12 < FDR_STRAT]
ref_known, ref_novel = set(rs.SYMBOL[~rs.novel]), set(rs.SYMBOL[rs.novel])
ref_known_loc, ref_novel_loc = set(rs.locus[~rs.novel]), set(rs.locus[rs.novel])
print(f"参照 FDR<{FDR_MAIN}: {(fdr_r12 < FDR_MAIN).sum()} 条 / {len(ref_main)} 基因")
print(f"分层参照: 已知 {len(ref_known)} 基因/{len(ref_known_loc)} 位点, "
      f"新位点 {len(ref_novel)} 基因/{len(ref_novel_loc)} 位点")

# ---------------------------------------------------------------- 读各 release
ext = {}
for rel, _, _ in RELEASES:
    if rel == "R12":
        ext[rel] = dict(zip(d.SNP, zip(d["beta.outcome"], d["se.outcome"])))
        continue
    fp = os.path.join(EXT, f"{rel}_{ENDPOINT}.tsv")
    if not os.path.exists(fp):
        print(f"  [缺] {rel}")
        continue
    # ⚠ 多等位位点在同一 chr:pos 上有多行。用 dict 直接覆盖会随机留下一行，
    # 可能是与工具变量不同的等位（本项目已有 indel 表示法歧义的前例）。
    # 故按位置存全部行，选等位时再匹配。
    m = {}
    n_rows = 0
    for r in csv.DictReader(open(fp, encoding="utf-8"), delimiter="\t"):
        try:
            rec = (r["ref"].upper(), r["alt"].upper(),
                   float(r["beta"]), float(r["sebeta"]))
        except ValueError:
            continue
        m.setdefault(r["SNP"], []).append(rec)
        n_rows += 1
    ext[rel] = m
    multi = sum(1 for v in m.values() if len(v) > 1)
    print(f"  {rel}: {n_rows:,} 行 / {len(m):,} 个位置"
          f"（其中 {multi} 个位置为多等位）")

avail = [r for r, _, _ in RELEASES if r in ext]
# 主分析：限定到所有可用 release 共有的变异
common = set(d.SNP)
for rel in avail:
    common &= set(ext[rel].keys())
print(f"\n所有 release 共有的工具变量: {len(common):,} / {len(d):,} "
      f"（丢弃 {len(d)-len(common):,}）")


def run(rel, snp_subset):
    """对一个 release 跑完整 MR，返回该轮的名单与统计量。"""
    sub = d[d.SNP.isin(snp_subset)].reset_index(drop=True)
    bo, se = [], []
    keep = []
    for i, row in sub.iterrows():
        if rel == "R12":
            b_out, s_out = row["beta.outcome"], row["se.outcome"]
        else:
            cands = ext[rel].get(row.SNP)
            if not cands:
                continue
            # 多等位位点：只接受等位与工具变量一致的那一行；无匹配则丢弃该记录
            hit_rec = None
            for ref, alt, b, s in cands:
                if row.ea_out == alt and row.oa_out == ref:
                    hit_rec = (b, s); break
                if row.ea_out == ref and row.oa_out == alt:
                    hit_rec = (-b, s); break
            if hit_rec is None:
                continue
            b_out, s_out = hit_rec
        if s_out <= 0:
            continue
        keep.append(i); bo.append(b_out); se.append(s_out)
    sub = sub.loc[keep].reset_index(drop=True)
    be = sub["beta.exposure"].values
    b = np.array(bo) / be
    s = np.array(se) / np.abs(be)
    z = b / s
    p = np.array([2 * ND.cdf(-abs(v)) if abs(v) < 37 else 0.0 for v in z])
    f = bh(p)
    hit = f < FDR_MAIN
    strat = f < FDR_STRAT
    return dict(
        n_tested=len(sub), n_nominal=int((p < 0.05).sum()), n_hits=int(hit.sum()),
        genes=set(sub.SYMBOL[hit]), loci=set(sub.locus[hit]),
        novel_genes=set(sub.SYMBOL[hit & sub.novel.values]),
        strat_genes=set(sub.SYMBOL[strat]), strat_loci=set(sub.locus[strat]),
        table=sub.assign(pval=p, FDR=f, OR=np.exp(b)))


def jac(a, b):
    return len(a & b) / len(a | b) if (a | b) else np.nan


rows, lists = [], []
for mode, subset in [("common", common), ("own", set(d.SNP))]:
    for rel, nc, nk in RELEASES:
        if rel not in ext:
            continue
        r = run(rel, subset)
        rec_known = len(r["strat_genes"] & ref_known) / len(ref_known)
        rec_novel = len(r["strat_genes"] & ref_novel) / len(ref_novel)
        rec_known_l = len(r["strat_loci"] & ref_known_loc) / len(ref_known_loc)
        rec_novel_l = len(r["strat_loci"] & ref_novel_loc) / len(ref_novel_loc)
        rows.append(dict(mode=mode, release=rel, cases=nc, controls=nk,
                         n_eff=n_eff(nc, nk), n_tested=r["n_tested"],
                         n_nominal=r["n_nominal"], n_hits=r["n_hits"],
                         n_genes=len(r["genes"]), n_loci=len(r["loci"]),
                         n_novel_genes=len(r["novel_genes"]),
                         jaccard_vs_R12=jac(r["genes"], ref_main),
                         recovery_main=len(r["genes"] & ref_main) / len(ref_main),
                         recovery_known=rec_known, recovery_novel=rec_novel,
                         diff_gene=rec_known - rec_novel,
                         recovery_known_locus=rec_known_l,
                         recovery_novel_locus=rec_novel_l,
                         diff_locus=rec_known_l - rec_novel_l))
        if mode == "common":
            lists.append(dict(release=rel,
                              genes_FDR05=",".join(sorted(r["genes"])),
                              novel_genes_FDR05=",".join(sorted(r["novel_genes"])) or "-",
                              genes_FDR20=",".join(sorted(r["strat_genes"]))))

out = pd.DataFrame(rows)
out.to_csv(f"{MR}/59a_release_trajectory.tsv", sep="\t", index=False)
pd.DataFrame(lists).to_csv(f"{MR}/59b_release_lists.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
print("\n" + "=" * 100)
print("观测轨迹（主分析 = 共有变异子集）")
print("=" * 100)
cols = ["release", "cases", "n_tested", "n_hits", "n_genes", "n_novel_genes",
        "recovery_main", "recovery_known", "recovery_novel", "diff_gene",
        "recovery_known_locus", "recovery_novel_locus", "diff_locus"]
print(out[out["mode"] == "common"][cols].to_string(
    index=False, float_format=lambda x: f"{x:.3f}"))

# ---------------------------------------------------------------- 阳性对照
print("\n" + "=" * 100)
print("过程阳性对照")
print("=" * 100)
r12 = out[(out["mode"] == "own") & (out.release == "R12")].iloc[0]
p1 = (r12.n_hits == 10) and (r12.n_genes == 6)
print(f"P1 重现校验 R12: {int(r12.n_hits)} 条 / {int(r12.n_genes)} 基因 "
      f"(应为 10 / 6) -> {'PASS' if p1 else 'FAIL'}")
sub = out[out["mode"] == "common"]
rho = np.corrcoef(sub.cases.rank(), sub.n_hits.rank())[0, 1]
print(f"P2 单调性 Spearman rho(cases, hits) = {rho:.3f} -> "
      f"{'PASS' if rho > 0 else 'FAIL'}")
early = lists[0] if lists else None
if early:
    known_big = {"PARP1", "VPS9D1-AS1", "CDK10", "SPATA33", "CHMP1A", "CTU2"}
    hit_early = known_big & set(early["genes_FDR05"].split(","))
    print(f"P3 最早 release({early['release']}) 中的已知大效应基因: "
          f"{sorted(hit_early) if hit_early else '无'} -> "
          f"{'PASS' if hit_early else 'FAIL'}")

# ---------------------------------------------------------------- 预测比对
pred = pd.read_csv(f"{MR}/58a_finngen_reference_predictions.tsv", sep="\t")
pred = pred.drop(columns=[c for c in ("cases", "controls", "n_eff")
                          if c in pred.columns])
chk = sub.merge(pred, on="release", suffixes=("_obs", "_pred"))
chk["hits_in_interval"] = (chk.n_hits >= chk.pred_hits_lo) & (chk.n_hits <= chk.pred_hits_hi)
chk["diff_in_interval"] = (chk.diff_gene >= chk.pred_diff_lo) & (chk.diff_gene <= chk.pred_diff_hi)
chk["recov_in_interval"] = ((chk.recovery_main >= chk.pred_recovery_main_lo) &
                            (chk.recovery_main <= chk.pred_recovery_main_hi))
keep = ["release", "cases", "n_hits", "pred_hits", "pred_hits_lo", "pred_hits_hi",
        "hits_in_interval", "recovery_main", "pred_recovery_main",
        "recov_in_interval", "diff_gene", "pred_known_minus_novel",
        "pred_diff_lo", "pred_diff_hi", "diff_in_interval",
        "diff_locus", "pred_diff_locus"]
chk[keep].to_csv(f"{MR}/59c_prediction_check.tsv", sep="\t", index=False)
print("\n" + "=" * 100)
print("观测 vs 登记预测")
print("=" * 100)
print(chk[keep].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

sub4 = sub[sub.release != "R12"]
n_pos_g = int((sub4.diff_gene > 0).sum())
n_pos_l = int((sub4.diff_locus > 0).sum())
print(f"\n预先指定的合并规则（符号检验，4 个 sub-R12 release）:")
print(f"  基因级差值为正: {n_pos_g}/4"
      f"{'  (单侧 p=0.0625)' if n_pos_g == 4 else ''}")
print(f"  位点级差值为正: {n_pos_l}/4"
      f"{'  (单侧 p=0.0625)' if n_pos_l == 4 else ''}")
print("\n写出 59a / 59b / 59c")
