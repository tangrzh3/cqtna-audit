"""
构建 SMR BESD 输入：每个候选 gene×profile 一个 .esd，每个 profile 一个 .flist
SNP 主键统一用 rsID（与 gwas_melanoma.ma 和 1000G 参考面板一致）
"""
import csv, gzip, os, sys, collections
import pyarrow.parquet as pq
import pyarrow.compute as pc

MR_DIR      = r"D:/R_ex/MR"
PARQUET_DIR = r"D:/Downloads/CD4_eqtl_step1_clean"
ESD_DIR     = os.path.join(MR_DIR, "esd_meta")
FINNGEN     = os.path.join(MR_DIR, "finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz")
os.makedirs(ESD_DIR, exist_ok=True)

# ---------- 1. 候选清单 ----------
cand = collections.defaultdict(set)          # profile -> {gene_id}
with open(os.path.join(MR_DIR, "13_meta_locus_annotation.tsv"), encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if float(r["pval"]) < 0.05:
            cand[r["exposure"].split("|")[1]].add(r["gene_id"])
print("候选 profile:", len(cand), "| gene×profile:", sum(len(v) for v in cand.values()), flush=True)

# ---------- 2. 读 cis 窗口 ----------
cols = ["gene_id", "chr", "pos", "other_allele", "effect_allele",
        "eaf", "beta", "se", "pval", "tss_distance"]
data = {}                                    # profile -> {gene -> [rows]}
need = set()                                 # (chr, pos)
for prof, genes in cand.items():
    fp = os.path.join(PARQUET_DIR, f"{prof}_step1_clean.parquet")
    t = pq.read_table(fp, columns=cols,
                      filters=[("gene_id", "in", list(genes))])
    d = collections.defaultdict(list)
    for b in t.to_batches():
        py = b.to_pydict()
        for i in range(len(py["gene_id"])):
            row = tuple(py[c][i] for c in cols)
            d[row[0]].append(row)
            need.add((str(row[1]), int(row[2])))
    data[prof] = d
    print(f"  {prof}: {len(d)} 基因, {t.num_rows} SNP", flush=True)
print("需映射的唯一位点:", len(need), flush=True)

# ---------- 3. chr:pos -> rsID ----------
rsmap = {}
with gzip.open(FINNGEN, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ic, ip, irs = hdr.index("#chrom"), hdr.index("pos"), hdr.index("rsids")
    for line in fh:
        f = line.split("\t", irs + 1)
        key = (f[ic], int(f[ip]))
        if key in need and key not in rsmap:
            rs = f[irs].split(",")[0]
            if rs.startswith("rs"):
                rsmap[key] = rs
print(f"映射成功: {len(rsmap)}/{len(need)} ({len(rsmap)/len(need):.1%})", flush=True)

# ---------- 4. 写 .esd + .flist ----------
for prof, d in data.items():
    flist = os.path.join(MR_DIR, f"flistM_{prof}.txt")
    n_probe = 0
    with open(flist, "w", newline="") as fo:
        fo.write("Chr\tProbeID\tGeneticDistance\tProbeBp\tGene\tOrientation\tPathOfEsd\n")
        for gene, rows in d.items():
            recs = []
            for (_g, ch, pos, oa, ea, eaf, beta, se, pval, tssd) in rows:
                rs = rsmap.get((str(ch), int(pos)))
                if rs is None:
                    continue
                recs.append((str(ch), rs, int(pos), ea, oa, eaf, beta, se, pval))
            if len(recs) < 20:               # HEIDI 需要足够 SNP
                continue
            seen, uniq = set(), []
            for r in recs:                   # rsID 去重
                if r[1] in seen:
                    continue
                seen.add(r[1]); uniq.append(r)
            ch = uniq[0][0]
            tss = int(rows[0][2]) - int(rows[0][9])
            path = os.path.join(ESD_DIR, f"{prof}__{gene}.esd")
            with open(path, "w", newline="") as fe:
                fe.write("Chr\tSNP\tBp\tA1\tA2\tFreq\tBeta\tse\tp\n")
                for (c, rs, bp, a1, a2, fr, b, s, p) in uniq:
                    fe.write(f"{c}\t{rs}\t{bp}\t{a1}\t{a2}\t{fr}\t{b}\t{s}\t{p}\n")
            fo.write(f"{ch}\t{gene}\t0\t{tss}\t{gene}\t+\t{path}\n")
            n_probe += 1
    print(f"  {prof}: 写出 {n_probe} 个探针 -> {os.path.basename(flist)}", flush=True)

print("\n完成。下一步用 smr --eqtl-flist ... --make-besd 生成 BESD", flush=True)
