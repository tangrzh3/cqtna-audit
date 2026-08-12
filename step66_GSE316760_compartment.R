#!/usr/bin/env Rscript
# ============================================================================
# Step 66  GSE316760（10x Visium 黑色素瘤 2 张切片）：区室归属的独立复现
#
# 目的：发现⑤（区室归属）目前只建立在 Thrane 一个数据集上（4 患者 8 切片，
#   100 µm spot）。本步用一个**独立数据集**、更细的 spot（Visium 55 µm）
#   原样重跑同一分析。
#
# *** 跑之前写死（勿在看到结果后修改）***
#
# 主检验（方向由 Thrane 结果预先指定）：
#   H1  TPI1 与糖酵解模块（**不含 TPI1**）正相关
#   H2  TPI1 与肿瘤/黑色素细胞区室正相关
#   H3  TPI1 与淋巴区室**不**正相关（Thrane 中为 −0.080）
#
# *** 阳性对照（预设；不通过则该切片判为无信息，不解释 TPI1 结果）***
#   PC1 HLA-C 与淋巴区室正相关，且与肿瘤区室的相关**低于**与淋巴区室的相关
#   PC2 SMC2 与增殖模块正相关
#   —— 这是 Thrane 中的镜像模式，是"阴性可判读"的前提
#
# 统计单元：spot（切片内），每张切片单独报告。**n=2 切片，不做 meta 分析、
#   不对 hot/cold 差异做任何推断**（1 对 1 比较无法与患者间差异区分）。
#
# QC 与处理与 Step 19 逐位一致：spot 总 counts ≥500 且检出基因 ≥200；
#   CP10K + log1p；区室评分 = marker 集在**切片内 z 标准化**后的均值。
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
})

DIR <- "D:/Downloads/GSE316760"
OUT <- "D:/R_ex/MR"

SAMPLES <- c(mel2 = "GSM9459774_mel2_filtered_feature_bc_matrix.h5",
             mel3 = "GSM9459775_mel3_filtered_feature_bc_matrix.h5")

MARK <- list(
  Tumour     = c("MLANA","PMEL","TYR","DCT","TYRP1","SOX10","MITF","S100B","PRAME"),
  Lymphoid   = c("CD2","CD3D","CD3E","CD3G","CD247","TRAC","LCK","IL7R","CD52",
                 "CCL5","CD8A","CD4","SKAP1","CD27"),
  Bcell      = c("MS4A1","CD79A","CD79B","IGHM","BANK1"),
  Myeloid    = c("LYZ","CD68","CD14","ITGAX","AIF1","TYROBP","FCER1G"),
  Stromal    = c("COL1A1","COL1A2","DCN","LUM","FBLN1","PDGFRB"),
  # 糖酵解：**不含 TPI1**（否则与被检验基因循环）
  Glycolysis = c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
                 "ALDOA","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA"),
  Prolif     = c("MKI67","TOP2A","CCNB1","CDK1","PCNA")
)
GENES <- c("TPI1", "HLA-C", "SMC2", "SPSB2")

zsc <- function(x) { s <- sd(x); if (!is.finite(s) || s == 0) return(x * 0); (x - mean(x)) / s }

res <- list(); pc <- list()
for (nm in names(SAMPLES)) {
  cat("\n========== ", nm, " ==========\n", sep = "")
  m <- Read10X_h5(file.path(DIR, SAMPLES[[nm]]))
  cat("原始:", nrow(m), "基因 x", ncol(m), "spot\n")

  keep <- Matrix::colSums(m) >= 500 & Matrix::colSums(m > 0) >= 200
  m <- m[, keep]
  cat("QC 后:", ncol(m), "spot (", sprintf("%.1f%%", 100 * mean(keep)), ")\n")

  cp <- Matrix::t(Matrix::t(m) / Matrix::colSums(m)) * 1e4
  lg <- log1p(cp)

  sc <- sapply(MARK, function(gs) {
    g <- intersect(gs, rownames(lg))
    if (length(g) < 3) return(rep(NA_real_, ncol(lg)))
    colMeans(apply(as.matrix(lg[g, , drop = FALSE]), 1, zsc) |> t())
  })
  sc <- as.data.table(sc)
  for (g in GENES) sc[[g]] <- if (g %in% rownames(lg)) zsc(as.numeric(lg[g, ])) else NA_real_
  sc[, depth := log(Matrix::colSums(m))]
  sc[, sample := nm]

  sp <- function(a, b) {
    ok <- is.finite(sc[[a]]) & is.finite(sc[[b]])
    if (sum(ok) < 50) return(c(NA, NA))
    ct <- suppressWarnings(cor.test(sc[[a]][ok], sc[[b]][ok], method = "spearman"))
    c(unname(ct$estimate), ct$p.value)
  }
  rows <- rbindlist(lapply(c("Glycolysis","Tumour","Lymphoid","Myeloid","Bcell",
                             "Stromal","Prolif"), function(cmp) {
    a <- sp("TPI1", cmp); b <- sp("HLA-C", cmp); d <- sp("SMC2", cmp)
    data.table(sample = nm, compartment = cmp,
               TPI1_rho = a[1], TPI1_p = a[2],
               HLAC_rho = b[1], HLAC_p = b[2],
               SMC2_rho = d[1], SMC2_p = d[2])
  }))
  print(rows[, .(compartment, TPI1_rho = round(TPI1_rho, 3),
                 HLAC_rho = round(HLAC_rho, 3), SMC2_rho = round(SMC2_rho, 3))])

  # ---- 预设阳性对照
  hl <- rows[compartment == "Lymphoid", HLAC_rho]; ht <- rows[compartment == "Tumour", HLAC_rho]
  sp2 <- rows[compartment == "Prolif", SMC2_rho]
  PC1 <- is.finite(hl) && hl > 0 && hl > ht
  PC2 <- is.finite(sp2) && sp2 > 0
  cat(sprintf("PC1 HLA-C 淋巴(%.3f) > 肿瘤(%.3f) 且为正 -> %s\n", hl, ht,
              ifelse(PC1, "PASS", "FAIL")))
  cat(sprintf("PC2 SMC2 ~ 增殖 = %.3f -> %s\n", sp2, ifelse(PC2, "PASS", "FAIL")))
  cat(sprintf("阳性对照总判定: %s\n",
              ifelse(PC1 && PC2, "通过 —— TPI1 结果可解释",
                     "未通过 —— 本切片判为无信息，不解释 TPI1 结果")))
  pc[[nm]] <- data.table(sample = nm, HLAC_lymphoid = hl, HLAC_tumour = ht,
                         SMC2_prolif = sp2, PC1 = PC1, PC2 = PC2,
                         usable = PC1 && PC2)
  res[[nm]] <- rows
}

R <- rbindlist(res); P <- rbindlist(pc)
fwrite(R, file.path(OUT, "66a_GSE316760_correlations.tsv"), sep = "\t")
fwrite(P, file.path(OUT, "66b_GSE316760_positive_controls.tsv"), sep = "\t")

cat("\n\n===== 预设主检验（仅对阳性对照通过的切片）=====\n")
ok <- P[usable == TRUE, sample]
if (!length(ok)) {
  cat("无切片通过阳性对照 -> 本步无信息，不得报告 TPI1 结果\n")
} else {
  S <- R[sample %in% ok]
  for (h in list(c("Glycolysis", "H1 TPI1 ~ 糖酵解(不含TPI1) 为正"),
                 c("Tumour",     "H2 TPI1 ~ 肿瘤区室 为正"),
                 c("Lymphoid",   "H3 TPI1 ~ 淋巴区室 不为正"))) {
    v <- S[compartment == h[1], TPI1_rho]
    cat(sprintf("%-38s %s  (每切片: %s)\n", h[2],
                ifelse(h[1] == "Lymphoid",
                       ifelse(all(v <= 0), "复现", "未复现"),
                       ifelse(all(v > 0), "复现", "未复现")),
                paste(sprintf("%.3f", v), collapse = ", ")))
  }
}
cat("\n写出 66a / 66b\n")
