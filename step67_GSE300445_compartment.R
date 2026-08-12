#!/usr/bin/env Rscript
# ============================================================================
# Step 67  GSE300445（Visium 原发皮肤黑色素瘤 4 切片，anti-PD1 语境）
#
# *** 先说明本步不能做什么 ***
# 该数据集与 GSE316760 同为**探针法 Visium**（Human Transcriptome Probe Set，
# 18,085 个靶标）。TPI1 (ENSG00000111669)、GAPDH、PKM、LDHA、ALDOA
# **不在探针集内** —— 厂商设计，非投稿者过滤，两个独立研究完全一致。
# → **TPI1 的区室归属在任何探针法 Visium 数据集上都不可测。**
#   Thrane 2018 之所以能测，是因为它是 poly-A 捕获的全转录组 legacy ST。
#
# *** 本步因此改为一个独立指定的次级分析（在得知 TPI1 不可测之后指定，如实标注）***
#
# 阳性对照（与 Step 19/66 同一套，预设）：
#   PC1 HLA-C 与淋巴区室正相关，且高于与肿瘤区室的相关
#   PC2 SMC2 与增殖模块正相关
#   不通过的切片判为无信息。
#
# 次级问题（方向由 Thrane 预先给定）：
#   缩减糖酵解模块（探针集内可得的 11 个酶，**不含 TPI1**）是否与肿瘤区室正相关、
#   与淋巴区室不正相关。这是通路层面的区室归属，**不是 TPI1 的复现**。
#
# n=4 切片，逐切片报告；不做 meta，不对样本间差异做推断。
# ============================================================================

suppressPackageStartupMessages({library(Matrix); library(data.table)})

DIR <- "D:/Downloads/GSE300445"
OUT <- "D:/R_ex/MR"

MARK <- list(
  Tumour   = c("MLANA","PMEL","TYR","DCT","TYRP1","SOX10","MITF","S100B","PRAME"),
  Lymphoid = c("CD2","CD3D","CD3E","CD3G","CD247","TRAC","LCK","IL7R","CD52",
               "CCL5","CD8A","CD4","SKAP1","CD27"),
  Myeloid  = c("LYZ","CD68","CD14","ITGAX","AIF1","TYROBP","FCER1G"),
  Stromal  = c("COL1A1","COL1A2","DCN","LUM","FBLN1","PDGFRB"),
  # 探针集内可得的糖酵解酶（不含 TPI1/GAPDH/PKM/LDHA/ALDOA —— 探针集不覆盖）
  GlycoReduced = c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
                   "PGK1","PGAM1","ENO1"),
  Prolif   = c("MKI67","TOP2A","CCNB1","CDK1","PCNA")
)
GENES <- c("HLA-C","SMC2","SPSB2")

zsc <- function(x) { s <- sd(x); if (!is.finite(s) || s == 0) return(x * 0); (x - mean(x)) / s }

samples <- unique(sub("_processed.*", "", list.files(DIR, pattern = "matrix.mtx.gz")))
cat("样本:", paste(samples, collapse = ", "), "\n")

res <- list()
for (s in samples) {
  cat("\n========== ", s, " ==========\n", sep = "")
  m  <- readMM(gzfile(file.path(DIR, paste0(s, "_processed_matrix.mtx.gz"))))
  ft <- fread(file.path(DIR, paste0(s, "_processed_features.tsv.gz")), header = FALSE)
  bc <- fread(file.path(DIR, paste0(s, "_processed_barcodes.tsv.gz")), header = FALSE)
  rownames(m) <- make.unique(ft$V2); colnames(m) <- bc$V1
  cat("原始:", nrow(m), "基因 x", ncol(m), "spot\n")

  keep <- Matrix::colSums(m) >= 500 & Matrix::colSums(m > 0) >= 200
  m <- m[, keep, drop = FALSE]
  cat("QC 后:", ncol(m), "spot\n")
  if (ncol(m) < 200) { cat("spot 过少，跳过\n"); next }

  lg <- log1p(Matrix::t(Matrix::t(m) / Matrix::colSums(m)) * 1e4)

  sc <- as.data.table(sapply(MARK, function(gs) {
    g <- intersect(gs, rownames(lg))
    if (length(g) < 3) return(rep(NA_real_, ncol(lg)))
    colMeans(t(apply(as.matrix(lg[g, , drop = FALSE]), 1, zsc)))
  }))
  for (g in GENES)
    set(sc, j = g, value = if (g %in% rownames(lg)) zsc(as.numeric(lg[g, ])) else NA_real_)

  sp <- function(a, b) {
    ok <- is.finite(sc[[a]]) & is.finite(sc[[b]])
    if (sum(ok) < 50) return(c(NA_real_, NA_real_))
    ct <- suppressWarnings(cor.test(sc[[a]][ok], sc[[b]][ok], method = "spearman"))
    c(unname(ct$estimate), ct$p.value)
  }
  rows <- rbindlist(lapply(c("Tumour","Lymphoid","Myeloid","Stromal","Prolif"), function(cmp) {
    a <- sp("GlycoReduced", cmp); b <- sp("HLA-C", cmp); d <- sp("SMC2", cmp)
    data.table(sample = s, compartment = cmp, n_spot = ncol(m),
               Glyco_rho = a[1], Glyco_p = a[2],
               HLAC_rho = b[1], HLAC_p = b[2], SMC2_rho = d[1], SMC2_p = d[2])
  }))
  print(rows[, .(compartment, Glyco_rho = round(Glyco_rho, 3),
                 HLAC_rho = round(HLAC_rho, 3), SMC2_rho = round(SMC2_rho, 3))])

  hl <- rows[compartment == "Lymphoid", HLAC_rho]; ht <- rows[compartment == "Tumour", HLAC_rho]
  s2 <- rows[compartment == "Prolif", SMC2_rho]
  PC1 <- is.finite(hl) && hl > 0 && hl > ht; PC2 <- is.finite(s2) && s2 > 0
  cat(sprintf("PC1 HLA-C 淋巴 %.3f > 肿瘤 %.3f -> %s | PC2 SMC2~增殖 %.3f -> %s | 总判定 %s\n",
              hl, ht, ifelse(PC1, "PASS", "FAIL"), s2, ifelse(PC2, "PASS", "FAIL"),
              ifelse(PC1 && PC2, "可解释", "无信息")))
  rows[, `:=`(PC1 = PC1, PC2 = PC2, usable = PC1 && PC2)]
  res[[s]] <- rows
}

R <- rbindlist(res)
fwrite(R, file.path(OUT, "67a_GSE300445_correlations.tsv"), sep = "\t")

cat("\n\n===== 汇总（仅阳性对照通过的切片）=====\n")
ok <- unique(R[usable == TRUE, sample])
cat("通过阳性对照的切片:", length(ok), "/", length(unique(R$sample)), "\n")
if (length(ok)) {
  S <- R[sample %in% ok]
  for (cmp in c("Tumour", "Lymphoid")) {
    v <- S[compartment == cmp, Glyco_rho]
    cat(sprintf("缩减糖酵解模块 ~ %-9s : %s  (每切片 %s)\n", cmp,
                ifelse(all(v > 0), "全为正", ifelse(all(v <= 0), "全为非正", "不一致")),
                paste(sprintf("%.3f", v), collapse = ", ")))
  }
}
cat("\n写出 67a\n")
