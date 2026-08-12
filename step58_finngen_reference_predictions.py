#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 58  预测登记：以 FinnGen R12 为参照的降功效模拟

目的：为 FinnGen release 轨迹分析（分析 A）**在查看 R8–R11 真实数据之前**登记定量预测。
设计与判据见 manuscript/PREREG_power_trajectory.md。

为什么不能直接用 53a：53a 的参照是 meta 名单（12,530 例）。release 轨迹的参照必须是
FinnGen R12（5,753 例），否则预测的是另一个量。

*** 预先写死的三件事（勿在看到结果后修改）***
1. 主要名单 = FDR<0.05（研究者实际使用的判决规则），用于整体稳定性
2. 分层对比 = FDR<0.20 的参照名单（因主阈值下新位点类别为空，见 PREREG §2.2-0）
3. 曲线只到观测功效为止，不外推

*** 方向性偏倚（预先声明）***
本模拟模拟的是"同规模的独立研究"，而真实 release 是 R12 的嵌套子集，共享个体使真实
一致性**高于**模拟预测。故：观测 > 预测 = 与嵌套结构一致，不否定模拟；
观测 < 预测 = 真实不稳定性超过模拟，模拟是下界。

输出：58a_finngen_reference_predictions.tsv
"""
import numpy as np
import pandas as pd
from scipy import stats

MR = r"D:/R_ex/MR"
rng = np.random.default_rng(1)
N_REP = 2000                      # 比 53 更多，因为参照名单更小、抖动更大
FDR_MAIN = 0.05
FDR_STRAT = 0.20                  # 分层对比用的放宽阈值（PREREG §2.2-0）
NOVEL = "潜在新位点"

# FinnGen R12 = 参照；各历史 release 的病例/对照来自各自的官方清单
R12_CASE, R12_CTRL = 5753, 378749
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


def wald(bo, seo, be):
    b = bo / be
    se = seo / np.abs(be)
    return 2 * stats.norm.sf(np.abs(b / se))


# ------------------------------------------------------------------ 载入
d = pd.read_csv(f"{MR}/06_locus_annotation.tsv", sep="\t")
d = d[["exposure", "SYMBOL", "category", "SNP", "beta.exposure",
       "beta.outcome", "se.outcome"]].dropna()
d = d[(d["beta.exposure"] != 0) & (d["se.outcome"] > 0)].reset_index(drop=True)
d["novel"] = d["category"].eq(NOVEL)
print(f"FinnGen R12 严格集: {len(d):,} 条 / {d.SYMBOL.nunique():,} 基因")

# ---- 位点 ID：同染色体上相邻工具变量间隔 >1 Mb 即另起一个位点
# 理由：分层参照名单里 CRHR1/KANSL1/KANSL1-AS1 同属 chr17q21.31 MAPT 倒位区，
# 按基因计数会让一个位点算三次。论文既有纪律是按独立位点计数，此处沿用。
d[["_chr", "_pos"]] = d.SNP.str.split(":", expand=True).iloc[:, :2]
d["_pos"] = d["_pos"].astype(int)
d = d.sort_values(["_chr", "_pos"]).reset_index(drop=True)
locus, cur, last_chr, last_pos = [], 0, None, None
for c, p in zip(d._chr, d._pos):
    if c != last_chr or p - last_pos > 1_000_000:
        cur += 1
    locus.append(cur)
    last_chr, last_pos = c, p
d["locus"] = locus
print(f"独立位点数（1 Mb 合并）: {d.locus.nunique():,}")

p_obs = wald(d["beta.outcome"].values, d["se.outcome"].values,
             d["beta.exposure"].values)
fdr_obs = bh(p_obs)

ref_main = set(d.SYMBOL[fdr_obs < FDR_MAIN])
ref_strat = d[fdr_obs < FDR_STRAT]
ref_known = set(ref_strat.SYMBOL[~ref_strat.novel])
ref_novel = set(ref_strat.SYMBOL[ref_strat.novel])
ref_known_loci = set(ref_strat.locus[~ref_strat.novel])
ref_novel_loci = set(ref_strat.locus[ref_strat.novel])
print(f"分层参照（位点级）: 已知 {len(ref_known_loci)} + 新位点 "
      f"{len(ref_novel_loci)} 个独立位点")

print(f"参照名单 FDR<{FDR_MAIN}: {(fdr_obs < FDR_MAIN).sum()} 条 / "
      f"{len(ref_main)} 基因（新位点 "
      f"{d[(fdr_obs < FDR_MAIN) & d.novel].SYMBOL.nunique()} 个）")
print(f"分层参照 FDR<{FDR_STRAT}: {len(ref_strat)} 条 / "
      f"已知 {len(ref_known)} + 新位点 {len(ref_novel)} 基因")

NE_R12 = n_eff(R12_CASE, R12_CTRL)


def simulate(target_ne, reps=N_REP):
    """把 FinnGen R12 降到 target_ne，与 R12 自身的名单比较。"""
    if target_ne >= NE_R12:                      # R12 与自身比较，构造为 1
        return dict(hits=[(fdr_obs < FDR_MAIN).sum()] * reps,
                    jac=[1.0] * reps, rec_main=[1.0] * reps,
                    rec_known=[1.0] * reps, rec_novel=[1.0] * reps,
                    rec_known_loc=[1.0] * reps, rec_novel_loc=[1.0] * reps,
                    novel_hits=[d[(fdr_obs < FDR_MAIN) & d.novel]
                                .SYMBOL.nunique()] * reps)
    se_o = d["se.outcome"].values
    se_s = se_o * np.sqrt(NE_R12 / target_ne)
    extra = np.sqrt(np.maximum(se_s ** 2 - se_o ** 2, 0))
    be = d["beta.exposure"].values
    out = {k: [] for k in
           ("hits", "jac", "rec_main", "rec_known", "rec_novel", "novel_hits",
            "rec_known_loc", "rec_novel_loc")}
    for _ in range(reps):
        bo = d["beta.outcome"].values + rng.normal(0, extra)
        f_main = bh(wald(bo, se_s, be))
        hit_main = f_main < FDR_MAIN
        g = set(d.SYMBOL[hit_main])
        out["hits"].append(int(hit_main.sum()))
        out["novel_hits"].append(int(d[hit_main & d.novel].SYMBOL.nunique()))
        out["jac"].append(len(g & ref_main) / len(g | ref_main)
                          if (g | ref_main) else np.nan)
        out["rec_main"].append(len(g & ref_main) / len(ref_main)
                               if ref_main else np.nan)
        # 分层：用同一次模拟、放宽阈值下的名单，与放宽阈值的参照比
        sub = d[f_main < FDR_STRAT]
        g_s, l_s = set(sub.SYMBOL), set(sub.locus)
        out["rec_known"].append(len(g_s & ref_known) / len(ref_known)
                                if ref_known else np.nan)
        out["rec_novel"].append(len(g_s & ref_novel) / len(ref_novel)
                                if ref_novel else np.nan)
        out["rec_known_loc"].append(len(l_s & ref_known_loci) / len(ref_known_loci)
                                    if ref_known_loci else np.nan)
        out["rec_novel_loc"].append(len(l_s & ref_novel_loci) / len(ref_novel_loci)
                                    if ref_novel_loci else np.nan)
    return out


def q(v, lo=5, hi=95):
    v = np.asarray(v, float)
    return np.nanmean(v), np.nanpercentile(v, lo), np.nanpercentile(v, hi)


rows = []
print("\n" + "=" * 92)
print("预测登记（在查看任何 R8–R11 数据之前）")
print("=" * 92)
hdr = (f"{'rel':5}{'cases':>7}{'N_eff':>9}{'frac':>7}"
       f"{'hits':>16}{'recov(FDR.05)':>18}{'known':>16}{'novel':>16}{'known-novel':>14}")
print(hdr)
for name, nc, nk in RELEASES:
    ne = n_eff(nc, nk)
    s = simulate(ne)
    h = q(s["hits"]); rm = q(s["rec_main"])
    rk = q(s["rec_known"]); rv = q(s["rec_novel"])
    diff = np.asarray(s["rec_known"], float) - np.asarray(s["rec_novel"], float)
    dm, dlo, dhi = q(diff)
    rkl = q(s["rec_known_loc"]); rvl = q(s["rec_novel_loc"])
    dl = np.asarray(s["rec_known_loc"], float) - np.asarray(s["rec_novel_loc"], float)
    dlm, dllo, dlhi = q(dl)
    print(f"{name:5}{nc:7d}{ne:9.0f}{ne/NE_R12:7.3f}"
          f"{h[0]:8.1f} [{h[1]:.0f},{h[2]:.0f}]"
          f"{rm[0]:10.3f} [{rm[1]:.2f},{rm[2]:.2f}]"
          f"{rk[0]:8.3f} [{rk[1]:.2f},{rk[2]:.2f}]"
          f"{rv[0]:8.3f} [{rv[1]:.2f},{rv[2]:.2f}]"
          f"{dm:8.3f} [{dlo:.2f},{dhi:.2f}]")
    rows.append(dict(release=name, cases=nc, controls=nk, n_eff=ne,
                     frac_of_R12=ne / NE_R12,
                     pred_hits=h[0], pred_hits_lo=h[1], pred_hits_hi=h[2],
                     pred_novel_hits=np.mean(s["novel_hits"]),
                     pred_jaccard=q(s["jac"])[0],
                     pred_recovery_main=rm[0], pred_recovery_main_lo=rm[1],
                     pred_recovery_main_hi=rm[2],
                     pred_recovery_known=rk[0], pred_recovery_known_lo=rk[1],
                     pred_recovery_known_hi=rk[2],
                     pred_recovery_novel=rv[0], pred_recovery_novel_lo=rv[1],
                     pred_recovery_novel_hi=rv[2],
                     pred_known_minus_novel=dm, pred_diff_lo=dlo, pred_diff_hi=dhi,
                     pred_recovery_known_locus=rkl[0],
                     pred_recovery_known_locus_lo=rkl[1],
                     pred_recovery_known_locus_hi=rkl[2],
                     pred_recovery_novel_locus=rvl[0],
                     pred_recovery_novel_locus_lo=rvl[1],
                     pred_recovery_novel_locus_hi=rvl[2],
                     pred_diff_locus=dlm, pred_diff_locus_lo=dllo,
                     pred_diff_locus_hi=dlhi))

pd.DataFrame(rows).to_csv(f"{MR}/58a_finngen_reference_predictions.tsv",
                          sep="\t", index=False)
print("\n写出 58a_finngen_reference_predictions.tsv")
print("\n参照名单（FDR<0.05）:", sorted(ref_main))
print("分层参照 已知位点:", sorted(ref_known))
print("分层参照 新位点  :", sorted(ref_novel))
