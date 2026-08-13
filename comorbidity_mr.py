"""
共病 / 多效性 MR：候选基因 eQTL -> 各表型
主检验：候选是否影响「黑色素细胞痣」（色素通路表型）
  - 不影响痣但影响黑色素瘤 → 支持非色素（免疫）机制
  - 阳性对照：MC1R/PARP1 等已知色素位点应对痣有强效应
"""
import csv, gzip, math, os, collections
from statistics import NormalDist

MR = r"D:/R_ex/MR"
PHENO_DIR = os.path.join(MR, "pheno")
ND = NormalDist()

PHENOS = {}
for l in open(os.path.join(MR, "fg_phenos.txt"), encoding="utf-8"):
    n, u, nc, ncon = l.rstrip("\n").split("\t")
    PHENOS[n] = (int(nc), int(ncon))

# ---------- 候选：四重验证基因 + 阳性对照（色素位点） ----------
FOCUS = {"ZFYVE19", "SMC2", "SPSB2", "TPI1", "HLA-C", "KIAA0040"}
CONTROL = {"CDK10", "SPATA33", "CHMP1A", "VPS9D1-AS1", "CTU2", "PARP1", "MDM4"}

# harmonise 后的等位基因方向（effect_allele.outcome 已与暴露效应等位对齐）
alle = {}
for r in csv.DictReader(open(os.path.join(MR, "01_harmonised_all.tsv"),
                             encoding="utf-8"), delimiter="\t"):
    alle[(r["exposure"], r["SNP"])] = (r["effect_allele.outcome"].upper(),
                                       r["other_allele.outcome"].upper())

ann = {}
for r in csv.DictReader(open(os.path.join(MR, "13_meta_locus_annotation.tsv"),
                             encoding="utf-8"), delimiter="\t"):
    a = alle.get((r["exposure"], r["SNP"]))
    if a is None:
        continue
    r["ea_out"], r["oa_out"] = a
    ann[r["exposure"]] = r

# 每个基因取 MR 最强的那条记录作为代表
best = {}
for e, r in ann.items():
    s = r["SYMBOL"]
    if s not in FOCUS and s not in CONTROL:
        continue
    try:
        p = float(r["pval"])
    except ValueError:
        continue
    if s not in best or p < float(best[s]["pval"]):
        best[s] = r

need = {}          # chr:pos -> symbol
for s, r in best.items():
    need[r["SNP"]] = s
print(f"检验基因: {len(best)} 个（候选 {len(FOCUS & set(best))} + 对照 {len(CONTROL & set(best))}）",
      flush=True)

# ---------- 逐表型提取并做 Wald ratio ----------
rows = []
for pheno, (nc, ncon) in PHENOS.items():
    fp = os.path.join(PHENO_DIR, f"{pheno}.gz")
    if not os.path.exists(fp):
        continue
    got = {}
    with gzip.open(fp, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ic, ip = h.index("#chrom"), h.index("pos")
        ir, ia = h.index("ref"), h.index("alt")
        ib, ise, ipv = h.index("beta"), h.index("sebeta"), h.index("pval")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            k = f[ic] + ":" + f[ip]
            if k in need:
                got[k] = (f[ir].upper(), f[ia].upper(),
                          float(f[ib]), float(f[ise]), float(f[ipv]))
    for k, sym in need.items():
        r = best[sym]
        g = got.get(k)
        if g is None:
            continue
        ref, alt, bo, seo, po = g
        ea = r["ea_out"]          # 已与暴露效应等位对齐
        if ea == alt:
            b_out = bo
        elif ea == ref:
            b_out = -bo
        else:
            continue
        be = float(r["beta_exposure"])
        if be == 0 or seo <= 0:
            continue
        b = b_out / be
        se = seo / abs(be)
        z = b / se
        p = 2 * ND.cdf(-abs(z)) if abs(z) < 37 else 0.0
        rows.append(dict(gene=sym, group="候选" if sym in FOCUS else "阳性对照",
                         pheno=pheno, n_case=nc, snp=k,
                         beta=b, se=se, OR=math.exp(b),
                         OR_L=math.exp(b - 1.96 * se), OR_U=math.exp(b + 1.96 * se),
                         pval=p, p_outcome_raw=po))
    print(f"  {pheno:<24} 匹配 {len(got)}/{len(need)}", flush=True)

# BH 校正（按表型内）
by = collections.defaultdict(list)
for r in rows:
    by[r["pheno"]].append(r)
for ph, rs in by.items():
    rs.sort(key=lambda x: x["pval"])
    n = len(rs)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        q = rs[i]["pval"] * n / (i + 1)
        prev = min(prev, q)
        rs[i]["FDR"] = min(prev, 1.0)

out = os.path.join(MR, "20_comorbidity_MR.tsv")
with open(out, "w", newline="", encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# ---------- 汇总 ----------
print("\n" + "=" * 92)
print("★ 主检验：黑色素细胞痣（CD2_BENIGN_MELANOCYTIC, 13357 例）")
print("=" * 92)
nv = sorted([r for r in rows if r["pheno"] == "CD2_BENIGN_MELANOCYTIC"],
            key=lambda x: (x["group"], x["pval"]))
print(f"{'组别':<10}{'基因':<12}{'OR':>7}{'95%CI':>18}{'p':>11}{'FDR':>9}")
print("-" * 70)
for r in nv:
    ci = f"{r['OR_L']:.2f}-{r['OR_U']:.2f}"
    print(f"{r['group']:<10}{r['gene']:<12}{r['OR']:>7.3f}{ci:>18}{r['pval']:>11.2e}{r['FDR']:>9.3f}")

print("\n" + "=" * 92)
print("自身免疫表型（候选基因，p<0.05 的记录）")
print("=" * 92)
for r in sorted([x for x in rows if x["pheno"] != "CD2_BENIGN_MELANOCYTIC"
                 and x["group"] == "候选" and x["pval"] < 0.05],
                key=lambda x: x["pval"]):
    print(f"  {r['gene']:<11}{r['pheno']:<26}OR={r['OR']:.3f}  p={r['pval']:.2e}  FDR={r['FDR']:.3f}  (病例{r['n_case']})")
print(f"\n输出: {out}")
