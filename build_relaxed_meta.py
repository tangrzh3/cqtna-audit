"""
放宽工具变量集：F>=5 & p<0.05，按 gene x profile 逐个 LD clumping (r2<0.1)
输入：Step 7 生成的 esd/*.esd（cis 全窗口，已带 rsID）
输出：17_relaxed_instruments_meta.tsv
"""
import csv, glob, os, subprocess, collections, sys

MR   = r"D:/R_ex/MR"
ESD  = os.path.join(MR, "esd_meta")
TMP  = os.path.join(MR, "clump_tmp_meta")
PLINK= os.path.join(MR, "bin", "plink2.exe")
BFILE= os.path.join(MR, "ref", "1kg_eur_meta")
os.makedirs(TMP, exist_ok=True)

F_MIN, P_MAX, R2 = 5.0, 0.05, 0.1
KB = 1000          # cis 窗口 ±500kb，用 1000kb 覆盖整个区域

files = sorted(glob.glob(os.path.join(ESD, "*.esd")))
print(f"待处理 gene x profile: {len(files)}", flush=True)

out_rows = []
stats = collections.Counter()
for n, fp in enumerate(files, 1):
    base = os.path.basename(fp)[:-4]
    prof, gene = base.split("__")

    recs = []
    with open(fp, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                b, se, p = float(r["Beta"]), float(r["se"]), float(r["p"])
            except ValueError:
                continue
            if se <= 0:
                continue
            F = (b / se) ** 2
            if F >= F_MIN and p < P_MAX:
                recs.append((r["SNP"], r["Chr"], r["Bp"], r["A1"], r["A2"],
                             r["Freq"], b, se, p, F))
    stats["cand_snps"] += len(recs)
    if len(recs) < 2:
        stats["skip_lt2_before_clump"] += 1
        continue

    # plink --clump 输入
    assoc = os.path.join(TMP, f"{base}.assoc")
    with open(assoc, "w", newline="") as fo:
        fo.write("ID\tP\n")
        for r in recs:
            fo.write(f"{r[0]}\t{r[8]}\n")

    res = subprocess.run(
        [PLINK, "--bfile", BFILE, "--clump", assoc,
         "--clump-p1", str(P_MAX), "--clump-p2", str(P_MAX),
         "--clump-r2", str(R2), "--clump-kb", str(KB),
         "--out", os.path.join(TMP, base), "--threads", "4"],
        capture_output=True, text=True)

    cl = os.path.join(TMP, base + ".clumps")
    if not os.path.exists(cl):
        stats["clump_no_output"] += 1
        continue

    keep = set()
    with open(cl, encoding="utf-8") as fh:
        hdr = fh.readline().split()
        idx = hdr.index("ID") if "ID" in hdr else 2
        for line in fh:
            p = line.split()
            if p:
                keep.add(p[idx])

    sel = [r for r in recs if r[0] in keep]
    if len(sel) < 2:
        stats["skip_lt2_after_clump"] += 1
    stats["kept_snps"] += len(sel)
    if sel:
        stats["exposures_with_snp"] += 1

    for (snp, ch, bp, a1, a2, fr, b, se, p, F) in sel:
        out_rows.append(dict(
            exposure=f"{gene}|{prof}", gene_id=gene, profile=prof,
            SNP=snp, chr=ch, pos=bp,
            effect_allele_exposure=a1, other_allele_exposure=a2,
            eaf_exposure=fr, beta_exposure=b, se_exposure=se,
            pval_exposure=p, F_stat=F))

    if n % 25 == 0:
        print(f"  {n}/{len(files)}  已保留 {len(out_rows)} 个工具变量", flush=True)

with open(os.path.join(MR, "17_relaxed_instruments_meta.tsv"), "w",
          newline="", encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(out_rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(out_rows)

nsnp = collections.Counter()
for r in out_rows:
    nsnp[r["exposure"]] += 1
dist = collections.Counter(nsnp.values())

print("\n===== 汇总 =====")
print("clumping 前候选 SNP :", stats["cand_snps"])
print("clumping 后保留 SNP :", stats["kept_snps"])
print("有工具变量的 exposure:", stats["exposures_with_snp"])
print("clump 前不足2个而跳过:", stats["skip_lt2_before_clump"])
print("clump 后不足2个       :", stats["skip_lt2_after_clump"])
print("\n每个 exposure 的独立 SNP 数分布:")
for k in sorted(dist):
    print(f"  {k} 个SNP: {dist[k]} 个exposure")
print("\n可跑多SNP方法(>=2)的 exposure:", sum(v for k, v in dist.items() if k >= 2))
print("可跑 MR-PRESSO(>=4)的 exposure:", sum(v for k, v in dist.items() if k >= 4))
