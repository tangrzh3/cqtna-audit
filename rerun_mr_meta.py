"""
用 meta 结局重跑 Wald ratio MR
复用 01_harmonised_all.tsv 的等位基因对齐（meta 保留了 FinnGen 的 ref/alt 定向）
"""
import csv, gzip, math, os, collections
from statistics import NormalDist

MR = r"D:/R_ex/MR"
META = os.path.join(MR, "meta_melanoma_final.tsv.gz")
ND = NormalDist()

# ---------- 1. meta 结局按 chr:pos 建索引 ----------
need = set()
har = list(csv.DictReader(open(os.path.join(MR, "01_harmonised_all.tsv"),
                               encoding="utf-8"), delimiter="\t"))
for r in har:
    need.add(r["SNP"])
print(f"harmonised 记录: {len(har)}  唯一SNP: {len(need)}", flush=True)

meta = {}
with gzip.open(META, "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ic, ip = h.index("#chrom"), h.index("pos")
    ir, ia = h.index("ref"), h.index("alt")
    ib, ise, ins = h.index("beta"), h.index("sebeta"), h.index("n_studies")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        key = f[ic] + ":" + f[ip]
        if key in need:
            meta[key] = (f[ir].upper(), f[ia].upper(),
                         float(f[ib]), float(f[ise]), f[ins])
print(f"在 meta 中匹配到: {len(meta)} ({len(meta)/len(need):.1%})", flush=True)

# ---------- 2. 替换结局效应并重算 Wald ratio ----------
out = []
skip = 0
for r in har:
    if r["mr_keep"] != "TRUE":
        continue
    m = meta.get(r["SNP"])
    if m is None:
        skip += 1
        continue
    ref, alt, bo, seo, ns = m
    ea = r["effect_allele.outcome"].upper()
    if ea == alt:
        b_out = bo
    elif ea == ref:
        b_out = -bo
    else:
        skip += 1
        continue

    try:
        be, see = float(r["beta.exposure"]), float(r["se.exposure"])
    except ValueError:
        continue
    if be == 0 or seo <= 0:
        continue

    b = b_out / be
    se = seo / abs(be)          # 一阶 delta
    z = b / se
    p = 2 * ND.cdf(-abs(z)) if abs(z) < 37 else \
        2 * math.exp(-z * z / 2) / (abs(z) * math.sqrt(2 * math.pi))

    out.append(dict(
        exposure=r["exposure"], id_exposure=r["id.exposure"], SNP=r["SNP"],
        gene_id=r.get("gene_id", ""), cell_type=r.get("cell_type", ""),
        timepoint=r.get("timepoint", ""), rsid=r.get("rsid", ""),
        F_stat=r.get("F_stat", ""), pval_exposure=r["pval.exposure"],
        beta_exposure=be, se_exposure=see,
        beta_outcome=b_out, se_outcome=seo, n_studies=ns,
        b=b, se=se, pval=p, OR=math.exp(b),
        OR_LCI=math.exp(b - 1.96 * se), OR_UCI=math.exp(b + 1.96 * se)))

print(f"可计算 MR: {len(out)}  跳过: {skip}", flush=True)

# ---------- 3. FDR（严格集内） ----------
def bh(rows, key="pval"):
    idx = sorted(range(len(rows)), key=lambda i: rows[i][key])
    n = len(idx)
    prev = 1.0
    for rank, i in enumerate(reversed(idx), 1):
        q = rows[i][key] * n / (n - rank + 1)
        prev = min(prev, q)
        rows[i]["FDR"] = min(prev, 1.0)

strict = [r for r in out if float(r["pval_exposure"]) < 5e-8
          and float(r["F_stat"] or 0) > 10]
bh(strict)
for r in out:
    r.setdefault("FDR", "")

with open(os.path.join(MR, "12_MR_meta_strict.tsv"), "w", newline="",
          encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(strict[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(strict)

nfdr = sum(1 for r in strict if r["FDR"] < 0.05)
nnom = sum(1 for r in strict if r["pval"] < 0.05)
genes = len(set(r["gene_id"] for r in strict if r["pval"] < 0.05))
print("\n===== meta 结局下的严格集 MR =====")
print(f"检验数        : {len(strict)}")
print(f"FDR<0.05      : {nfdr}   (FinnGen 单研究时: 10)")
print(f"名义 p<0.05   : {nnom}  涉及基因 {genes}   (FinnGen 时: 284 / 132)")
