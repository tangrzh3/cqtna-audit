#!/usr/bin/env Rscript
# ============================================================================
# Step 69  区室归属的两项加强
#   (A) GSE115978：对测序深度做**细胞层面匹配**后重跑（Step 68 的深度混杂）
#   (B) GSE72056（Tirosh）：独立队列复现（TISCH 注释 + TISCH 表达矩阵）
#
# *** 跑之前写死 ***
# (A) 匹配：在**每个样本内**，把 Mal 与 T.CD4 细胞按 log10 深度做最近邻 1:1 匹配
#     （caliper = 0.1 log10 ≈ 26%）。匹配后残余 SMD 必须 |SMD| < 0.1，否则该样本弃用。
#     判据：TPI1 差值方向与 Step 68 一致且样本一致性 >= 12/16 视为稳健。
# (B) 独立复现：同一预设方向（H1 TPI1 在恶性细胞高于 CD4⁺T），
#     同一阳性对照（PC1 MLANA、PC2 PTPRC），同一统计单元（患者）。
#     ⚠ GSE72056 按 CD45 分选且细胞数少，若可配对患者 < 5 则判为功效不足，不解释。
# ============================================================================

suppressPackageStartupMessages({library(data.table); library(hdf5r); library(Matrix)})
OUT <- "D:/R_ex/MR"
MIN_CELLS <- 10L; CALIPER <- 0.1
LOCKED <- c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
            "ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA")

smd <- function(a, b) (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)

# ============================ (A) GSE115978 深度匹配 ========================
cat("############ (A) GSE115978 深度匹配 ############\n")
D <- "D:/数据/黑色素瘤单细胞人/GSE115978"
ann <- fread(file.path(D, "GSE115978_cell.annotations.csv.gz"))
cnt <- fread(file.path(D, "GSE115978_counts.csv.gz"))
g <- cnt[[1]]; cnt[, 1 := NULL]
m <- as.matrix(cnt); rownames(m) <- g; rm(cnt); gc()
ann <- ann[match(colnames(m), cells)]
cs <- colSums(m)
lg <- log2(sweep(m, 2, cs, "/") * 1e5 + 1)
gl <- intersect(LOCKED, rownames(lg))

dt <- data.table(sample = ann$samples, type = ann$cell.types,
                 depth = log10(cs), TPI1 = as.numeric(lg["TPI1", ]),
                 Glyco = colMeans(lg[gl, , drop = FALSE]),
                 GlycoNo = colMeans(lg[setdiff(gl, "TPI1"), , drop = FALSE]))
dt <- dt[type %in% c("Mal", "T.CD4")]

res <- list()
for (s in unique(dt$sample)) {
  A <- dt[sample == s & type == "Mal"]; B <- dt[sample == s & type == "T.CD4"]
  if (nrow(A) < MIN_CELLS || nrow(B) < MIN_CELLS) next
  # 1:1 最近邻匹配（不放回），以细胞数少的一侧为基准
  if (nrow(A) <= nrow(B)) { X <- A; Y <- B; flip <- FALSE } else { X <- B; Y <- A; flip <- TRUE }
  used <- rep(FALSE, nrow(Y)); pick <- integer(nrow(X))
  for (i in seq_len(nrow(X))) {
    d <- abs(Y$depth - X$depth[i]); d[used] <- Inf
    j <- which.min(d)
    pick[i] <- if (is.finite(d[j]) && d[j] <= CALIPER) j else NA_integer_
    if (!is.na(pick[i])) used[j] <- TRUE
  }
  ok <- !is.na(pick)
  if (sum(ok) < MIN_CELLS) next
  Xm <- X[ok]; Ym <- Y[pick[ok]]
  mal <- if (flip) Ym else Xm; cd4 <- if (flip) Xm else Ym
  sm <- smd(mal$depth, cd4$depth)
  if (abs(sm) >= 0.1) { cat(sprintf("  %-8s 匹配后 SMD=%.3f 未达标，弃用\n", s, sm)); next }
  res[[s]] <- data.table(sample = s, n_pair = nrow(mal), smd_depth = sm,
                         TPI1_mal = mean(mal$TPI1), TPI1_cd4 = mean(cd4$TPI1),
                         Glyco_mal = mean(mal$Glyco), Glyco_cd4 = mean(cd4$Glyco),
                         GlycoNo_mal = mean(mal$GlycoNo), GlycoNo_cd4 = mean(cd4$GlycoNo))
}
M <- rbindlist(res)
cat(sprintf("\n可用样本 %d（每样本匹配 %d–%d 对细胞），残余深度 SMD 中位 %.4f\n",
            nrow(M), min(M$n_pair), max(M$n_pair), median(M$smd_depth)))
tst <- function(a, b, lab) {
  d <- M[[a]] - M[[b]]
  p <- suppressWarnings(wilcox.test(M[[a]], M[[b]], paired = TRUE)$p.value)
  cat(sprintf("%-18s Mal %.3f vs CD4 %.3f | 差 %+.3f | %d/%d | p=%.3g\n",
              lab, mean(M[[a]]), mean(M[[b]]), mean(d), sum(d > 0), length(d), p))
  data.table(dataset = "GSE115978_depthmatched", variable = lab,
             mal = mean(M[[a]]), cd4 = mean(M[[b]]), diff = mean(d),
             n_mal_higher = sum(d > 0), n = length(d), p = p)
}
A_out <- rbind(tst("TPI1_mal","TPI1_cd4","TPI1"),
               tst("Glyco_mal","Glyco_cd4","Glyco(16)"),
               tst("GlycoNo_mal","GlycoNo_cd4","Glyco(去TPI1)"))
fwrite(M, file.path(OUT, "69a_GSE115978_depth_matched_samples.tsv"), sep = "\t")
rm(m, lg); gc()

# ============================ (B) GSE72056 独立复现 =========================
cat("\n############ (B) GSE72056（Tirosh）独立复现 ############\n")
meta <- fread("D:/Downloads/SKCM_GSE72056_CellMetainfo_table.tsv")
setnames(meta, c("Celltype (malignancy)", "Celltype (major-lineage)"), c("malig", "lineage"),
         skip_absent = TRUE)
cat("细胞类型分布:\n"); print(head(meta[, .N, by = .(malig, lineage)][order(-N)], 12))

f <- H5File$new("D:/Downloads/SKCM_GSE72056_expression.h5", "r")
gn <- f[["matrix/features/name"]]$read(); bc <- f[["matrix/barcodes"]]$read()
shp <- f[["matrix/shape"]]$read()
X <- sparseMatrix(i = f[["matrix/indices"]]$read() + 1, p = f[["matrix/indptr"]]$read(),
                  x = f[["matrix/data"]]$read(), dims = shp)
f$close_all()
rownames(X) <- make.unique(gn); colnames(X) <- bc
cat("表达矩阵:", nrow(X), "x", ncol(X), "\n")
meta <- meta[match(colnames(X), Cell)]
stopifnot(!any(is.na(meta$Cell)))

# TISCH 矩阵已做 CPM/log 处理；此处不再二次归一化（Step 25 的教训）
cat("每细胞列和中位:", round(median(colSums(X)), 1), "-> ",
    ifelse(median(colSums(X)) > 1e5, "疑似已归一化，直接用", "检查"), "\n")
gl2 <- intersect(LOCKED, rownames(X))
cat("signature 可得:", length(gl2), "/", length(LOCKED), "\n")

d2 <- data.table(sample = meta$Patient, malig = meta$malig, lineage = meta$lineage,
                 TPI1 = as.numeric(X["TPI1", ]),
                 MLANA = as.numeric(X["MLANA", ]), PTPRC = as.numeric(X["PTPRC", ]),
                 Glyco = colMeans(as.matrix(X[gl2, , drop = FALSE])),
                 GlycoNo = colMeans(as.matrix(X[setdiff(gl2, "TPI1"), , drop = FALSE])))
d2[, grp := fifelse(malig == "Malignant cells", "Mal",
             fifelse(lineage == "CD4Tconv", "T.CD4", NA_character_))]
d2 <- d2[!is.na(grp)]
cat("Mal", sum(d2$grp == "Mal"), "| CD4Tconv", sum(d2$grp == "T.CD4"), "\n")

p2 <- d2[, .(n = .N, TPI1 = mean(TPI1), MLANA = mean(MLANA), PTPRC = mean(PTPRC),
             Glyco = mean(Glyco), GlycoNo = mean(GlycoNo)), by = .(sample, grp)]
w2 <- dcast(p2, sample ~ grp, value.var = c("n","TPI1","MLANA","PTPRC","Glyco","GlycoNo"))
w2 <- w2[n_Mal >= MIN_CELLS & n_T.CD4 >= MIN_CELLS]
cat("可配对患者:", nrow(w2), "\n")

if (nrow(w2) < 5) {
  cat("⚠ 可配对患者 < 5 -> 功效不足，按预设不解释结果\n")
  B_out <- data.table()
} else {
  t2 <- function(a, b, lab) {
    d <- w2[[a]] - w2[[b]]
    p <- suppressWarnings(wilcox.test(w2[[a]], w2[[b]], paired = TRUE)$p.value)
    cat(sprintf("%-18s Mal %.3f vs CD4 %.3f | 差 %+.3f | %d/%d | p=%.3g\n",
                lab, mean(w2[[a]]), mean(w2[[b]]), mean(d), sum(d > 0), length(d), p))
    data.table(dataset = "GSE72056", variable = lab, mal = mean(w2[[a]]),
               cd4 = mean(w2[[b]]), diff = mean(d), n_mal_higher = sum(d > 0),
               n = length(d), p = p)
  }
  cat("\n-- 阳性对照 --\n")
  pc <- rbind(t2("MLANA_Mal","MLANA_T.CD4","PC1 MLANA"),
              t2("PTPRC_Mal","PTPRC_T.CD4","PC2 PTPRC"))
  ok1 <- pc[1, diff > 0]; ok2 <- pc[2, diff < 0]
  cat(sprintf("PC1 %s | PC2 %s -> %s\n", ifelse(ok1,"PASS","FAIL"), ifelse(ok2,"PASS","FAIL"),
              ifelse(ok1 && ok2, "可解释", "无信息")))
  cat("\n-- 预设主检验 --\n")
  B_out <- rbind(pc, t2("TPI1_Mal","TPI1_T.CD4","TPI1"),
                 t2("Glyco_Mal","Glyco_T.CD4","Glyco"),
                 t2("GlycoNo_Mal","GlycoNo_T.CD4","Glyco(去TPI1)"))
  fwrite(w2, file.path(OUT, "69b_GSE72056_paired.tsv"), sep = "\t")
}

fwrite(rbind(A_out, B_out, fill = TRUE), file.path(OUT, "69c_tests.tsv"), sep = "\t")
cat("\n写出 69a / 69b / 69c\n")
