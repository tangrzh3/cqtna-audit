#!/usr/bin/env Rscript
# ============================================================================
# Step 77  Pozniak (Cell 2024) 队列：患者侧结果的独立复现
#
# 设计、假设方向、判读表、判废条件：manuscript/PREREG_pozniak_replication.md
# **该文件在读入本数据的任何表达值之前已定稿。本脚本只执行，不改设计。**
#
# 数据：KU Leuven RDR doi:10.48804/GSAXBN
#   Entire_TME.rds     59,170 细胞 / 46 样本 / 23 患者，含 BT(治疗前)、OT(治疗中)
#   Malignant_cells.rds  提供 patient → Response(R/NR) 映射（README D01/D02）
#
# 统计与发现队列（Step 74）逐位一致：
#   样本为单位、每样本 ≥20 个 CD4⁺T、16 基因 z 后取均值、
#   置换检验（打乱 Response 标签，保留基因相关结构）+ Wilcoxon 并报、
#   **不用基因方向的二项检验**、去 TPI1 版本并报、重复患者三种处理
# ============================================================================

suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
MR <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"; P <- "D:/Downloads/Pozniak"
set.seed(1); N_PERM <- 20000; MIN_CELLS <- 20L

LOCKED <- c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
            "ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA")
TCELL <- c("CD3D","CD3E","CD3G","TRAC","TRBC2","CD2")
CD4P  <- c("CD4","IL7R","CD40LG","MAL","LTB","TRAT1","ANXA1","CCR7","AQP3")
CD8P  <- c("CD8A","CD8B","GZMK","NKG7","CCL5","GZMB","PRF1","KLRD1","GNLY")

# ---------------------------------------------------------------- 1. 应答映射
mal <- readRDS(file.path(P, "Malignant_cells.rds"))
mm <- mal@meta.data
cat("Malignant_cells meta 列:", paste(names(mm), collapse=", "), "\n")
rc <- grep("^respon", names(mm), ignore.case = TRUE, value = TRUE)[1]
pc <- grep("GC|patient", names(mm), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(rc), !is.na(pc))
map <- unique(data.table(patient = as.character(mm[[pc]]),
                         response = as.character(mm[[rc]])))
map <- map[!is.na(response) & response != ""]
cat("应答映射:", nrow(map), "条 |", paste(capture.output(print(table(map$response))), collapse=" "), "\n")
rm(mal); gc()

# ---------------------------------------------------------------- 2. TME 对象
o <- readRDS(file.path(P, "Entire_TME.rds"))
md <- o@meta.data
md$patient <- as.character(md[["GC number"]])
md$tp <- md[["BT/OT"]]
md$sample <- paste0(md$patient, "_", md$tp)
Xr <- GetAssayData(o, assay = "RNA", layer = "counts")
if (nrow(Xr) == 0) Xr <- GetAssayData(o, assay = "RNA", layer = "data")
cs0 <- Matrix::colSums(Xr[, sample(ncol(Xr), 300)])
cat(sprintf("RNA counts 层: 每细胞和中位 %.0f (CV %.2f), 最大值 %.0f
",
            median(cs0), sd(cs0)/mean(cs0), max(Xr[, 1:200])))
# 预注册 §5：判据是数据本身是否已归一化。此处 CV>>0 且值为整数尺度 → 未归一化
X <- log1p(Matrix::t(Matrix::t(Xr) / Matrix::colSums(Xr)) * 1e4)
cat("→ 判定：未归一化，按库大小归一化 + log1p（与 GSE115978 同一处理）
")

# ---- 归一化性质核查（Step 25 的教训：判据是数据本身，不是与某文件是否一致）
cs <- Matrix::colSums(X[, sample(ncol(X), 300)])
cat(sprintf("\ndata 层每细胞和: 中位 %.1f (CV %.2f) | 最大表达值 %.2f\n",
            median(cs), sd(cs)/mean(cs), max(X[, 1:200])))
cat("→ 判定：已 log 归一化，直接使用，不再二次归一化\n")

# ---------------------------------------------------------------- 3. CD4 定义
zs <- function(v) { s <- sd(v); if (!is.finite(s) || s == 0) return(v*0); (v - mean(v))/s }
score <- function(gs) {
  g <- intersect(gs, rownames(X))
  if (length(g) < 3) return(rep(NA_real_, ncol(X)))
  colMeans(t(apply(as.matrix(X[g, , drop = FALSE]), 1, zs)))
}
sT <- score(TCELL); s4 <- score(CD4P); s8 <- score(CD8P)
is_cd4 <- sT > 0 & s4 > s8
cat(sprintf("\nCD4⁺T（多基因评分）: %d / %d 细胞 (%.1f%%)\n",
            sum(is_cd4), length(is_cd4), 100*mean(is_cd4)))
# cluster 层面的第二定义：簇内 CD4 比例 > 50% 的簇
cl <- as.character(md$seurat_clusters)
frac <- tapply(is_cd4, cl, mean)
cd4_clusters <- names(frac)[frac > .5]
is_cd4_cl <- cl %in% cd4_clusters
cat(sprintf("CD4 主导簇: %s → %d 细胞\n", paste(cd4_clusters, collapse=","), sum(is_cd4_cl)))

# ---------------------------------------------------------------- 4. 检验
gl <- intersect(LOCKED, rownames(X))
cat("signature 可得:", length(gl), "/", length(LOCKED), "\n")

run <- function(keep, defname) {
  d <- data.table(sample = md$sample[keep], patient = md$patient[keep],
                  tp = md$tp[keep])
  A <- t(as.matrix(X[gl, keep, drop = FALSE]))
  agg <- as.data.table(A)[, lapply(.SD, mean), by = .(sample = d$sample)]
  n <- d[, .N, by = sample]
  info <- unique(d[, .(sample, patient, tp)])[n, on = "sample"]
  info <- merge(info, map, by = "patient", all.x = TRUE)
  info <- info[N >= MIN_CELLS & !is.na(response)]
  agg <- agg[match(info$sample, sample)]
  cat(sprintf("\n[%s] 可用样本 %d | 患者 %d\n", defname, nrow(info), uniqueN(info$patient)))
  print(info[, .N, by = .(tp, response)][order(tp, response)])

  out <- list()
  for (t0 in c("BT", "OT")) for (gset in list(list(gl, "locked16"),
                                              list(setdiff(gl,"TPI1"), "locked15_noTPI1"))) {
    sub <- info[tp == t0]
    if (uniqueN(sub$response) < 2 || min(table(sub$response)) < 3) next
    for (mode in c("all", "one_per_patient", "single_sample_patients")) {
      idx <- switch(mode,
        all = seq_len(nrow(sub)),
        one_per_patient = sub[, .I[sample(.N,1)], by = patient]$V1,
        single_sample_patients = which(!sub$patient %in% sub[, .N, by=patient][N>1, patient]))
      s2 <- sub[idx]
      if (uniqueN(s2$response) < 2 || min(table(s2$response)) < 3) next
      M <- as.matrix(agg[match(s2$sample, agg$sample), gset[[1]], with = FALSE])
      sc <- rowMeans(apply(M, 2, zs))
      y <- s2$response
      isNR <- grepl("^Non", y); nr <- sc[isNR]; r <- sc[!isNR]
      dd <- mean(nr) - mean(r)
      w <- suppressWarnings(wilcox.test(nr, r))$p.value
      bs <- replicate(4000, mean(sample(nr,replace=TRUE)) - mean(sample(r,replace=TRUE)))
      pc2 <- sum(replicate(N_PERM, { yp <- sample(y)
        { q <- grepl("^Non", yp); mean(sc[q]) - mean(sc[!q]) } }) >= dd)
      out[[length(out)+1]] <- data.table(definition = defname, timepoint = t0,
        genes = gset[[2]], subset = mode, n_R = sum(!isNR), n_NR = sum(isNR),
        score_diff = dd, ci_lo = quantile(bs,.025), ci_hi = quantile(bs,.975),
        wilcox_p = w, perm_p = (pc2+1)/(N_PERM+1))
    }
  }
  rbindlist(out)
}

res <- rbind(run(is_cd4, "marker score"), run(is_cd4_cl, "CD4-dominant clusters"))

# ---------------------------------------------------------------- 5. 阳性对照
cat("\n===== 预设阳性对照 =====\n")
mal_like <- score(c("MLANA","PMEL","TYR","DCT","TYRP1","SOX10","MITF")) > 0.5 & !is_cd4
pc_tab <- rbindlist(lapply(c("MLANA","PTPRC","TPI1"), function(g) {
  if (!g %in% rownames(X)) return(NULL)
  v <- as.numeric(X[g, ])
  data.table(gene = g, malignant = mean(v[mal_like]), CD4 = mean(v[is_cd4]),
             diff = mean(v[mal_like]) - mean(v[is_cd4]))
}))
print(pc_tab)
PC1 <- pc_tab[gene=="MLANA", diff > 0]; PC2 <- pc_tab[gene=="PTPRC", diff < 0]
PC3 <- pc_tab[gene=="TPI1",  diff > 0]
cat(sprintf("PC1 MLANA %s | PC2 PTPRC %s | PC3 TPI1(区室归属第三队列) %s\n",
            ifelse(PC1,"PASS","FAIL"), ifelse(PC2,"PASS","FAIL"), ifelse(PC3,"PASS","FAIL")))

cat("\n===== 主检验 =====\n")
print(res[, .(definition, timepoint, genes, subset, n_R, n_NR,
              d = round(score_diff,3), lo = round(ci_lo,3), hi = round(ci_hi,3),
              wilcox = signif(wilcox_p,3), perm = signif(perm_p,3))])
fwrite(res, file.path(MR, "77a_pozniak_replication.tsv"), sep = "\t")
fwrite(pc_tab, file.path(MR, "77b_pozniak_positive_controls.tsv"), sep = "\t")
cat("\n写出 77a / 77b\n")
