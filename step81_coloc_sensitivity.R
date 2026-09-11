#!/usr/bin/env Rscript
# ============================================================================
# Step 81  coloc 的窗口与先验敏感性分析（回应 R1 "technical failings" 第 5 条）
#
# 审稿意见：「对所有 coloc 结论进行窗口、先验和信号模型敏感性分析」。
# 信号模型部分已在 Step 60 做过（PARP1 排除多信号解释；ZFYVE19/TPI1 无 credible set）。
# **窗口与先验从未变动过**——全部结果都用 coloc.abf 默认先验
# (p1=1e-4, p2=1e-4, p12=1e-5) 与全 cis 窗口。本步补上。
#
# *** 跑之前写死 ***
# 检验对象：正文引用过 PP.H4 的全部位点
#   ④b 的两个案例 PARP1、ZFYVE19；应用示例 TPI1；发现①②的 MC1R 簇代表
#   VPS9D1-AS1、SPATA33、CHMP1A、CDK10、CTU2；meta 轮候选 SMC2、KIAA0040、SPSB2
#
# 变动的两个维度：
#   窗口   ±100 kb / ±250 kb / ±500 kb（默认为全 cis 窗口）
#   p12    1e-6（保守）/ 5e-6 / 1e-5（默认）/ 5e-5（宽松）
#
# 判读（预设）：
#   稳健  = 结论（是否满足 PP.H4/(H3+H4)>0.7 且 H3+H4>0.5）在全部 12 种组合下一致
#   敏感  = 结论随窗口或先验改变 → **须在正文标注该位点的结论依赖参数选择**
# ⚠ 若关键位点（PARP1、ZFYVE19、TPI1）敏感，正文相应结论必须加限定，不得只报默认值。
# ============================================================================

suppressPackageStartupMessages({library(arrow); library(data.table); library(coloc)})
MR <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"; PARQ <- "D:/Downloads/CD4_eqtl_step1_clean"

N_EQTL <- c(CD4_Naive_uns_0h=99, CD4_Naive_stim_16h=99, CD4_Naive_stim_40h=89,
            CD4_Naive_stim_5d=85, CD4_Memory_uns_0h=100, CD4_Memory_stim_16h=95,
            CD4_Memory_stim_40h=89, CD4_Memory_stim_5d=90)
META_N <- 801629; META_S <- 0.015631

TARGETS <- fread(text = "gene,profile
PARP1,CD4_Naive_stim_40h
ZFYVE19,CD4_Naive_stim_5d
TPI1,CD4_Naive_stim_16h
VPS9D1-AS1,CD4_Naive_stim_16h
SPATA33,CD4_Naive_stim_40h
CHMP1A,CD4_Naive_stim_40h
CDK10,CD4_Naive_stim_40h
CTU2,CD4_Naive_stim_40h
SMC2,CD4_Naive_stim_40h
KIAA0040,CD4_Memory_stim_16h
SPSB2,CD4_Memory_uns_0h")

ann <- fread(file.path(MR, "13_meta_locus_annotation.tsv"))
meta <- fread(cmd = paste("gzip -dc", shQuote(file.path(MR, "meta_melanoma_final.tsv.gz"))),
              select = c("#chrom","pos","ref","alt","beta","sebeta"))
setnames(meta, "#chrom", "chrom"); meta[, chrom := as.character(chrom)]
setkey(meta, chrom, pos)

WINDOWS <- c(1e5, 2.5e5, 5e5, NA)          # NA = 全 cis 窗口
P12     <- c(1e-6, 5e-6, 1e-5, 5e-5)

out <- list()
for (k in seq_len(nrow(TARGETS))) {
  g <- TARGETS$gene[k]; prof <- TARGETS$profile[k]
  a <- ann[SYMBOL == g][1]
  if (!nrow(a)) { cat("[skip]", g, "\n"); next }
  ens <- a$exposure; ens <- sub("\\|.*", "", ens)
  gid <- ann[SYMBOL == g, gene_id][1]
  if (is.na(gid)) gid <- sub("_.*", "", ens)
  e <- tryCatch(as.data.table(read_parquet(
        file.path(PARQ, paste0(prof, "_step1_clean.parquet")),
        col_select = c("gene_id","chr","pos","eaf","pval","beta","se"))), error=function(z) NULL)
  if (is.null(e)) { cat("[skip parquet]", g, "\n"); next }
  e <- e[gene_id == gid]
  if (!nrow(e)) { cat("[skip no gene]", g, gid, "\n"); next }
  setnames(e, c("chr","eaf"), c("chrom","maf"))
  e[, chrom := as.character(chrom)]
  e[, maf := pmin(maf, 1 - maf)]
  lead <- e[which.min(pval)]
  for (w in WINDOWS) {
    ee <- if (is.na(w)) e else e[abs(pos - lead$pos) <= w]
    mm <- meta[chrom == ee$chrom[1] & pos %between% c(min(ee$pos), max(ee$pos))]
    d <- merge(ee[, .(pos, maf, beta_e = beta, se_e = se)],
               mm[, .(pos, beta_o = beta, se_o = sebeta)], by = "pos")
    d <- d[is.finite(beta_e) & is.finite(se_e) & is.finite(beta_o) & se_o > 0 & maf > 0 & maf < 1]
    d <- unique(d, by = "pos")
    if (nrow(d) < 50) { cat("   [", g, w, "] nrow(d)=", nrow(d), "跳过
"); next }
    for (p12 in P12) {
      r <- tryCatch(suppressWarnings(coloc.abf(
        list(beta=d$beta_e, varbeta=d$se_e^2, MAF=d$maf, N=N_EQTL[[prof]], type="quant", snp=as.character(d$pos)),
        list(beta=d$beta_o, varbeta=d$se_o^2, MAF=d$maf, N=META_N, s=META_S, type="cc", snp=as.character(d$pos)),
        p12 = p12)), error = function(z) { cat("   [coloc 错误]", g, conditionMessage(z), "
"); NULL })
      if (is.null(r)) next
      s <- as.list(r$summary)
      ratio <- s$PP.H4.abf / (s$PP.H3.abf + s$PP.H4.abf)
      out[[length(out)+1]] <- data.table(gene=g, profile=prof,
        window = ifelse(is.na(w), "full cis", paste0("±", w/1000, " kb")),
        p12 = p12, nsnp = nrow(d), PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf,
        H3H4 = s$PP.H3.abf + s$PP.H4.abf, ratio = ratio,
        pass = (!is.na(ratio) && ratio > 0.7 && (s$PP.H3.abf + s$PP.H4.abf) > 0.5))
    }
  }
  cat("done", g, "\n")
}

res <- rbindlist(out)
fwrite(res, file.path(MR, "81a_coloc_sensitivity.tsv"), sep = "\t")

cat("\n===== 每个位点在 16 种参数组合下的结论一致性 =====\n")
summ <- res[, .(n_combo = .N, n_pass = sum(pass),
                H4_min = round(min(PP.H4),3), H4_max = round(max(PP.H4),3),
                verdict = ifelse(sum(pass) == .N, "稳健(全过)",
                          ifelse(sum(pass) == 0, "稳健(全不过)", "★敏感"))), by = gene]
print(summ[order(verdict, gene)])
cat("\n===== 敏感位点的明细 =====\n")
for (g in summ[verdict == "★敏感", gene])
  print(res[gene == g, .(window, p12, PP.H3 = round(PP.H3,3), PP.H4 = round(PP.H4,3),
                         ratio = round(ratio,3), pass)])
cat("\n写出 81a\n")
