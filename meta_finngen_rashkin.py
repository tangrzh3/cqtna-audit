"""
FinnGen R12 + Rashkin 2020 黑色素瘤 GWAS 固定效应 IVW meta 分析
Rashkin 只有 OR 和 p 值，SE 由 se = |log(OR)| / |Phi^-1(p/2)| 反推
输出格式与 FinnGen 一致，便于下游脚本复用
"""
import gzip, math, os
from statistics import NormalDist

MR = r"D:/R_ex/MR"
FG = os.path.join(MR, "finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz")
RK = os.path.join(MR, "rashkin2020_melanoma.h.tsv.gz")
OUT = os.path.join(MR, "meta_melanoma_finngen_rashkin.tsv.gz")

ND = NormalDist()
PMIN = 1e-300


def se_from_or_p(orv, p):
    """由 OR 与双侧 p 反推 log-OR 的标准误"""
    if orv is None or orv <= 0:
        return None, None
    b = math.log(orv)
    if abs(b) < 1e-10:
        return None, None
    if p <= 0:
        p = PMIN
    if p >= 1:
        return None, None
    try:
        z = -ND.inv_cdf(max(p, PMIN) / 2.0)
    except Exception:
        return None, None
    if z <= 0:
        return None, None
    return b, abs(b) / z


# ---------- 1. 读 Rashkin 到内存（键：chrom*1e9+pos）----------
print("读取 Rashkin ...", flush=True)
rk = {}
bad = 0
with gzip.open(RK, "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    i_c, i_p = h.index("hm_chrom"), h.index("hm_pos")
    i_ea, i_oa = h.index("hm_effect_allele"), h.index("hm_other_allele")
    i_or, i_pv = h.index("hm_odds_ratio"), h.index("p_value")
    i_af = h.index("hm_effect_allele_frequency")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        try:
            ch = f[i_c]
            if not ch.isdigit():
                continue
            pos = int(f[i_p])
            orv = float(f[i_or]) if f[i_or] not in ("NA", "") else None
            pv = float(f[i_pv]) if f[i_pv] not in ("NA", "") else None
        except (ValueError, IndexError):
            bad += 1
            continue
        if orv is None or pv is None:
            continue
        b, se = se_from_or_p(orv, pv)
        if se is None:
            bad += 1
            continue
        rk[int(ch) * 1_000_000_000 + pos] = (
            f[i_ea].upper(), f[i_oa].upper(), b, se, pv)
print(f"  Rashkin 可用变异: {len(rk):,}  跳过: {bad:,}", flush=True)

# ---------- 2. 流式读 FinnGen 并 meta ----------
print("meta 中 ...", flush=True)
FG_N, RK_N = 384502, 417127
n_out = n_both = n_fg_only = flip = 0
Q_sig = 0

with gzip.open(FG, "rt") as fh, gzip.open(OUT, "wt", newline="") as out:
    h = fh.readline().rstrip("\n").split("\t")
    ic, ip, ir, ia = h.index("#chrom"), h.index("pos"), h.index("ref"), h.index("alt")
    irs, ib, ise, ipv = h.index("rsids"), h.index("beta"), h.index("sebeta"), h.index("pval")
    iaf, ing = h.index("af_alt"), h.index("nearest_genes")
    out.write("#chrom\tpos\tref\talt\trsids\tnearest_genes\tpval\tbeta\tsebeta\t"
              "af_alt\tn_studies\tQ_pval\tbeta_fg\tbeta_rk\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        ch = f[ic]
        if not ch.isdigit():
            continue
        pos = int(f[ip])
        try:
            b1, s1 = float(f[ib]), float(f[ise])
        except ValueError:
            continue
        if s1 <= 0:
            continue
        ref, alt = f[ir].upper(), f[ia].upper()

        rec = rk.get(int(ch) * 1_000_000_000 + pos)
        b2 = s2 = None
        if rec is not None:
            ea2, oa2, b2r, s2r, _ = rec
            if ea2 == alt and oa2 == ref:      # 同向
                b2, s2 = b2r, s2r
            elif ea2 == ref and oa2 == alt:    # 反向，翻转
                b2, s2 = -b2r, s2r
                flip += 1

        if b2 is None:                          # 仅 FinnGen
            bm, sm, ns, qp = b1, s1, 1, ""
            n_fg_only += 1
        else:                                   # IVW 固定效应
            w1, w2 = 1.0 / s1**2, 1.0 / s2**2
            bm = (b1 * w1 + b2 * w2) / (w1 + w2)
            sm = math.sqrt(1.0 / (w1 + w2))
            Q = w1 * (b1 - bm) ** 2 + w2 * (b2 - bm) ** 2
            qp = math.erfc(math.sqrt(Q / 2.0))  # df=1 卡方上尾: P(X>Q)=erfc(sqrt(Q/2))
            if qp < 0.05:
                Q_sig += 1
            qp = f"{qp:.4g}"
            ns = 2
            n_both += 1

        try:
            z = abs(bm) / sm
            pm = 2.0 * ND.cdf(-z) if z < 37 else math.exp(-z * z / 2) / (z * math.sqrt(2 * math.pi)) * 2
        except Exception:
            continue
        out.write(f"{ch}\t{pos}\t{ref}\t{alt}\t{f[irs]}\t{f[ing]}\t{pm:.6g}\t"
                  f"{bm:.6g}\t{sm:.6g}\t{f[iaf]}\t{ns}\t{qp}\t"
                  f"{b1:.6g}\t{'' if b2 is None else f'{b2:.6g}'}\n")
        n_out += 1
        if n_out % 5_000_000 == 0:
            print(f"  已写 {n_out:,}", flush=True)

print(f"\n输出变异: {n_out:,}")
print(f"  两研究均有 (meta): {n_both:,}")
print(f"  仅 FinnGen       : {n_fg_only:,}")
print(f"  等位基因翻转     : {flip:,}")
print(f"  Cochran Q p<0.05 : {Q_sig:,} ({Q_sig/max(n_both,1):.2%})")
print(f"\n合并规模: 12,530 例 / 789,099 对照  N=801,629  s=0.015631")
