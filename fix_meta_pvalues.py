"""
给 meta 结果补 mlogp 列并对 p 设下限，解决 p<1e-308 下溢为 0 的问题。
mlogp 由 z 在对数空间直接算，不经过 p，因此不受下溢影响。
"""
import gzip, math, os

MR = r"D:/R_ex/MR"
SRC = os.path.join(MR, "meta_melanoma_finngen_rashkin.tsv.gz")
DST = os.path.join(MR, "meta_melanoma_final.tsv.gz")
PFLOOR = 1e-300
LOG10 = math.log(10.0)
SQRT2PI = math.sqrt(2.0 * math.pi)


def mlog10p(z):
    """-log10(双侧p)，大 z 用渐近展开避免下溢"""
    z = abs(z)
    if z < 37:
        # p = erfc(z/sqrt2)
        p = math.erfc(z / math.sqrt(2.0))
        if p > 0:
            return -math.log10(p)
    # 渐近: p ≈ 2*phi(z)/z
    ln_p = -z * z / 2.0 - math.log(z * SQRT2PI) + math.log(2.0)
    return -ln_p / LOG10


n = fixed = 0
with gzip.open(SRC, "rt") as fh, gzip.open(DST, "wt", newline="") as out:
    h = fh.readline().rstrip("\n").split("\t")
    ib, ise, ipv = h.index("beta"), h.index("sebeta"), h.index("pval")
    out.write("\t".join(h[:ipv + 1] + ["mlogp"] + h[ipv + 1:]) + "\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        try:
            b, se = float(f[ib]), float(f[ise])
            p = float(f[ipv])
        except ValueError:
            continue
        ml = mlog10p(b / se) if se > 0 else 0.0
        if p <= 0:
            p = PFLOOR
            fixed += 1
        f[ipv] = f"{p:.6g}"
        out.write("\t".join(f[:ipv + 1] + [f"{ml:.4f}"] + f[ipv + 1:]) + "\n")
        n += 1

print(f"处理 {n:,} 行，修正下溢 p=0 共 {fixed:,} 个")
print(f"输出: {DST}")
