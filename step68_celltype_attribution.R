#!/usr/bin/env Rscript
# ============================================================================
# Step 68  区室归属的细胞层面验证（GSE115978 全细胞类型）
#
# 起因：发现⑤（区室归属）目前只有 Thrane 空转一个数据集，而两次尝试用现代
#   空间数据复现都失败——探针法 Visium 的探针集不含 TPI1（Step 66/67）。
#   但空转本来就是**间接**手段：它用 spot 混合物推断哪个区室贡献信号。
#   本数据集直接给出细胞类型标注，可以**直接**回答同一个问题。
#
# 数据：GSE115978（Jerby-Arnon），7,186 细胞，含 Mal 2018 / T.CD4 856 /
#   T.CD8 1759 / B 818 / Macrophage 420 / CAF 106 / Endo 104 / NK 92。
#   ⚠ 该文件为**原始 counts**（Step 27 已确认 CV=1.15，未归一化）→ 必须做库大小归一化。
#
# *** 跑之前写死 ***
# 主假设（方向由 Thrane 空转结果预先给定）：
#   H1  TPI1 在恶性细胞中高于 CD4⁺T
#   H2  锁定糖酵解 signature（16 基因）在恶性细胞中高于 CD4⁺T
# 阳性对照（预设，不通过则本步无信息）：
#   PC1 MLANA 在 Mal 中高于 T.CD4
#   PC2 PTPRC 在 T.CD4 中高于 Mal
# 统计单元：**样本（患者）**，非细胞。取每样本内各细胞类型的均值，配对比较。
#   要求该样本两种细胞类型各 ≥10 个细胞。
#
# ⚠ 本数据按 CD45 分选，细胞比例**不代表组织构成**，故不计算"组织贡献占比"，
#   只比较每细胞表达量。组织水平由谁主导，需结合肿瘤纯度另行说明。
# ============================================================================

suppressPackageStartupMessages({library(data.table); library(Matrix)})

D   <- "D:/数据/黑色素瘤单细胞人/GSE115978"
OUT <- "D:/R_ex/MR"
MIN_CELLS <- 10L

LOCKED <- c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
            "ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA")

ann <- fread(file.path(D, "GSE115978_cell.annotations.csv.gz"))
cat("注释:", nrow(ann), "细胞 |", uniqueN(ann$samples), "样本\n")

cat("读取 counts（约 50 MB gz）...\n")
cnt <- fread(file.path(D, "GSE115978_counts.csv.gz"))
genes <- cnt[[1]]; cnt[, 1 := NULL]
m <- as.matrix(cnt); rownames(m) <- genes; rm(cnt); gc()
cat("counts:", nrow(m), "基因 x", ncol(m), "细胞\n")

stopifnot(all(colnames(m) %in% ann$cells))
ann <- ann[match(colnames(m), cells)]

# 库大小归一化 + log（Step 27 已判定该文件未归一化）
cs <- colSums(m)
lg <- log2(sweep(m, 2, cs, "/") * 1e5 + 1)

zsig <- function(gs) {
  g <- intersect(gs, rownames(lg))
  cat("  signature 可得基因:", length(g), "/", length(gs), "\n")
  colMeans(lg[g, , drop = FALSE])
}
sig <- zsig(LOCKED)
sig_noTPI1 <- zsig(setdiff(LOCKED, "TPI1"))

dt <- data.table(cell = colnames(m), sample = ann$samples, type = ann$cell.types,
                 TPI1 = as.numeric(lg["TPI1", ]),
                 MLANA = as.numeric(lg["MLANA", ]),
                 PTPRC = as.numeric(lg["PTPRC", ]),
                 HLAC  = if ("HLA-C" %in% rownames(lg)) as.numeric(lg["HLA-C", ]) else NA_real_,
                 Glyco = sig, GlycoNoTPI1 = sig_noTPI1,
                 depth = log10(cs))

# ---- 各细胞类型的检出率与均值（描述性）
desc <- dt[, .(n = .N,
               TPI1_detect = mean(TPI1 > 0), TPI1_mean = mean(TPI1),
               Glyco_mean = mean(Glyco), depth = mean(depth)),
           by = type][order(-n)]
cat("\n===== 各细胞类型（全部细胞，描述性）=====\n"); print(desc)

# ---- 配对比较：Mal vs T.CD4，样本为单位
pair <- dt[type %in% c("Mal", "T.CD4"),
           .(n = .N, TPI1 = mean(TPI1), Glyco = mean(Glyco),
             GlycoNoTPI1 = mean(GlycoNoTPI1), MLANA = mean(MLANA),
             PTPRC = mean(PTPRC), depth = mean(depth)),
           by = .(sample, type)]
w <- dcast(pair, sample ~ type,
           value.var = c("n","TPI1","Glyco","GlycoNoTPI1","MLANA","PTPRC","depth"))
w <- w[n_Mal >= MIN_CELLS & n_T.CD4 >= MIN_CELLS]
cat(sprintf("\n可配对样本: %d（每类 >=%d 细胞）\n", nrow(w), MIN_CELLS))

tst <- function(a, b, lab) {
  d <- w[[a]] - w[[b]]
  p <- suppressWarnings(wilcox.test(w[[a]], w[[b]], paired = TRUE)$p.value)
  cat(sprintf("%-16s Mal %.3f vs CD4 %.3f | 差 %+.3f | %d/%d 样本 Mal 更高 | p=%.3g\n",
              lab, mean(w[[a]]), mean(w[[b]]), mean(d), sum(d > 0), length(d), p))
  data.table(variable = lab, mal = mean(w[[a]]), cd4 = mean(w[[b]]),
             diff = mean(d), n_mal_higher = sum(d > 0), n = length(d), p = p)
}
cat("\n===== 阳性对照 =====\n")
pc <- rbind(tst("MLANA_Mal","MLANA_T.CD4","PC1 MLANA"),
            tst("PTPRC_Mal","PTPRC_T.CD4","PC2 PTPRC"))
PC1 <- pc[variable == "PC1 MLANA", diff > 0 & p < 0.05]
PC2 <- pc[variable == "PC2 PTPRC", diff < 0 & p < 0.05]
cat(sprintf("PC1 %s | PC2 %s | 总判定: %s\n", ifelse(PC1,"PASS","FAIL"),
            ifelse(PC2,"PASS","FAIL"),
            ifelse(PC1 && PC2, "通过，主检验可解释", "未通过，本步无信息")))

cat("\n===== 预设主检验 =====\n")
main <- rbind(tst("TPI1_Mal","TPI1_T.CD4","H1 TPI1"),
              tst("Glyco_Mal","Glyco_T.CD4","H2 Glyco(16)"),
              tst("GlycoNoTPI1_Mal","GlycoNoTPI1_T.CD4","   Glyco(去TPI1)"))
cat("\n（测序深度对照）\n"); tst("depth_Mal","depth_T.CD4","depth")

fwrite(desc, file.path(OUT, "68a_celltype_descriptive.tsv"), sep = "\t")
fwrite(w,    file.path(OUT, "68b_paired_sample_level.tsv"), sep = "\t")
fwrite(rbind(pc, main), file.path(OUT, "68c_tests.tsv"), sep = "\t")
cat("\n写出 68a / 68b / 68c\n")
