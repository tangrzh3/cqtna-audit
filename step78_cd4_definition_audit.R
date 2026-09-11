#!/usr/bin/env Rscript
# ============================================================================
# Step 78  审计 Step 77 所用的 CD4 定义
#
# 起因：Step 77 的 CD4 定义是多基因评分的备选方案，且**未做任何纯度检查**——
# 而本文自己的发现⑦正是"按评分切分会切出纯度梯度"。若该定义把 naive CD8
# 大量误判为 CD4，则复现失败可能是定义问题而非生物学问题。
#
# 本步只做诊断，不改结论。三个问题：
#   Q1 被判为 CD4 的细胞，其纯度如何？（CD4/CD8A/CD8B 检出率、CD8 双阴比例）
#   Q2 该定义与"正规注释"相比如何？（若 GC_all_immune_sharing.rds 提供细胞类型）
#   Q3 收紧定义后，Step 77 的结论是否改变？
#
# 参照：发现队列 GSE120575 固定 CD4 群体的纯度为
#   CD4 检出 32.5%、CD8A 检出 52.9%、CD8 双阴 44.6%（HANDOFF 已如实报告）
# ============================================================================

suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
P <- "D:/Downloads/Pozniak"; MR <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"

TCELL <- c("CD3D","CD3E","CD3G","TRAC","TRBC2","CD2")
CD4P  <- c("CD4","IL7R","CD40LG","MAL","LTB","TRAT1","ANXA1","CCR7","AQP3")
CD8P  <- c("CD8A","CD8B","GZMK","NKG7","CCL5","GZMB","PRF1","KLRD1","GNLY")

o <- readRDS(file.path(P, "Entire_TME.rds"))
Xr <- GetAssayData(o, assay = "RNA", layer = "counts")
X <- log1p(Matrix::t(Matrix::t(Xr) / Matrix::colSums(Xr)) * 1e4)
md <- o@meta.data
zs <- function(v) { s <- sd(v); if (!is.finite(s) || s == 0) return(v*0); (v - mean(v))/s }
score <- function(gs, sub = NULL) {
  g <- intersect(gs, rownames(X))
  M <- if (is.null(sub)) X[g, , drop = FALSE] else X[g, sub, drop = FALSE]
  colMeans(t(apply(as.matrix(M), 1, zs)))
}
det <- function(g, idx) if (g %in% rownames(X)) mean(X[g, idx] > 0) else NA_real_

sT <- score(TCELL); s4 <- score(CD4P); s8 <- score(CD8P)
defs <- list(
  "Step77 原定义 (sT>0 & s4>s8)" = sT > 0 & s4 > s8,
  "收紧 A: sT>0.5 & s4>s8"       = sT > 0.5 & s4 > s8,
  "收紧 B: 上 + CD8A/CD8B 未检出" = sT > 0.5 & s4 > s8 & X["CD8A",] == 0 & X["CD8B",] == 0,
  "收紧 C: 上 + CD4 或 CD40LG 检出" = sT > 0.5 & s4 > s8 &
      X["CD8A",] == 0 & X["CD8B",] == 0 & (X["CD4",] > 0 | X["CD40LG",] > 0)
)

cat("===== Q1  各定义下的群体纯度 =====\n")
cat(sprintf("%-34s %7s %8s %8s %8s %9s %9s\n",
            "定义", "n", "CD4检出", "CD8A检出", "CD8B检出", "CD8双阴", "CD3E检出"))
purity <- rbindlist(lapply(names(defs), function(nm) {
  i <- which(defs[[nm]])
  r <- data.table(definition = nm, n = length(i),
                  CD4 = det("CD4", i), CD8A = det("CD8A", i), CD8B = det("CD8B", i),
                  CD8_dneg = mean(X["CD8A", i] == 0 & X["CD8B", i] == 0),
                  CD3E = det("CD3E", i))
  cat(sprintf("%-34s %7d %8.1f%% %8.1f%% %8.1f%% %9.1f%% %9.1f%%\n", nm, r$n,
              100*r$CD4, 100*r$CD8A, 100*r$CD8B, 100*r$CD8_dneg, 100*r$CD3E))
  r
}))
cat("\n参照：发现队列 GSE120575 的固定 CD4 群体 —— CD4 检出 32.5%、CD8A 52.9%、CD8 双阴 44.6%\n")

cat("\n===== Q1b  被原定义判为 CD4 的细胞里，有多少像 naive CD8？ =====\n")
i0 <- which(defs[[1]])
naive_like <- X["CD8A", i0] > 0 | X["CD8B", i0] > 0
cat(sprintf("  原定义 n=%d，其中 CD8A 或 CD8B 有表达: %d (%.1f%%)\n",
            length(i0), sum(naive_like), 100*mean(naive_like)))
cat(sprintf("  其中同时 CCR7 或 MAL 阳性（naive 样）: %.1f%%\n",
            100*mean((X["CCR7", i0] > 0 | X["MAL", i0] > 0) & naive_like)))

cat("\n===== Q2  与簇结构的一致性 =====\n")
cl <- as.character(md$seurat_clusters)
tab <- data.table(cluster = cl, cd4 = defs[[1]], cd4c = defs[[4]])[
  , .(n = .N, frac_orig = mean(cd4), frac_strict = mean(cd4c)), by = cluster][order(-n)]
print(head(tab[frac_orig > .2], 12))

fwrite(purity, file.path(MR, "78a_cd4_definition_purity.tsv"), sep = "\t")
saveRDS(defs, file.path(MR, "78b_cd4_definitions.rds"))
cat("\n写出 78a / 78b\n")
