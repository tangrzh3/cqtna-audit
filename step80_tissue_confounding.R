#!/usr/bin/env Rscript
# ============================================================================
# Step 80  Pozniak 复现：组织部位是否混杂了 BT 的结果
#
# ⚠ 本步是**事后的混杂检查**，不是预注册检验。标注为敏感性分析。
#   起因：该队列的样本来自 Subcutis / Lymph node / Skin 三个部位，
#   而 BT 只有 8 R vs 3 NR。若部位在两组间分布不均，
#   "应答差异"可能是"部位差异"。
#
# 三问（顺序固定，先问混杂是否存在，再谈调整）：
#   Q1 部位在 R/NR 之间是否失衡？
#   Q2 糖酵解 score 本身是否随部位变化？（若否，则 Q1 的失衡不重要）
#   Q3 若两者皆是，调整后 BT 的效应还剩多少？
#
# ⚠ 预先声明：若调整后效应大幅缩小或反向，**以调整后结果为准**，
#   不得因为未调整版本更好看而保留之。
# ============================================================================

suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
MR <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"; P <- "D:/Downloads/Pozniak"; set.seed(1); N_PERM <- 20000
LOCKED <- c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
            "ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA")

o <- readRDS(file.path(P, "GC_all_immune_sharing.rds"))
md <- as.data.table(o@meta.data, keep.rownames = "cell")
Xr <- GetAssayData(o, assay = "RNA", layer = "counts")
X <- log1p(Matrix::t(Matrix::t(Xr) / Matrix::colSums(Xr)) * 1e4)
gl <- intersect(LOCKED, rownames(X))
zs <- function(v){s<-sd(v); if(!is.finite(s)||s==0) return(v*0); (v-mean(v))/s}

cd4 <- md[low_resolution_clusters == "CD4_Tcells"]
A <- t(as.matrix(X[gl, cd4$cell, drop = FALSE]))
agg <- as.data.table(A)[, lapply(.SD, mean), by = .(sample = cd4$sample_ID)]
info <- unique(cd4[, .(sample = sample_ID, patient = patient_ID, tp = Timepoint,
                       response = Response, tissue = Tissue)])
info <- info[cd4[, .N, by = .(sample = sample_ID)], on = "sample"][N >= 20]
agg <- agg[match(info$sample, sample)]
info[, score := rowMeans(apply(as.matrix(agg[, gl, with = FALSE]), 2, zs))]

cat("===== Q1  部位 × 应答（按时点）=====\n")
for (t0 in c("BT", "OT")) {
  cat("\n--", t0, "--\n")
  print(dcast(info[tp == t0, .N, by = .(tissue, response)], tissue ~ response,
              value.var = "N", fill = 0))
}

cat("\n===== Q2  糖酵解 score 是否随部位变化 =====\n")
print(info[, .(n = .N, score = round(mean(score), 3), sd = round(sd(score), 3)),
           by = tissue][order(-n)])
kw <- kruskal.test(score ~ factor(tissue), data = info)
cat(sprintf("  Kruskal-Wallis（全部样本）: p = %.3f\n", kw$p.value))
kw2 <- kruskal.test(score ~ factor(tissue), data = info[tp == "BT"])
cat(sprintf("  仅 BT: p = %.3f\n", kw2$p.value))

cat("\n===== Q3  调整部位后 BT 还剩多少 =====\n")
bt <- info[tp == "BT"]
cat("BT 样本明细（部位 | 应答 | score）:\n")
print(bt[order(response, tissue), .(sample, tissue, response, score = round(score, 3))])

d_raw <- mean(bt[response == "NR", score]) - mean(bt[response == "R", score])
cat(sprintf("\n  未调整 Δ = %+.3f\n", d_raw))

# (a) 部位内中心化后再比较
bt[, score_c := score - mean(score), by = tissue]
d_adj <- mean(bt[response == "NR", score_c]) - mean(bt[response == "R", score_c])
pp <- sum(replicate(N_PERM, {
  y <- sample(bt$response); mean(bt$score_c[y == "NR"]) - mean(bt$score_c[y == "R"])
}) >= d_adj)
cat(sprintf("  部位内中心化后 Δ = %+.3f | 置换 p = %.3f\n", d_adj, (pp + 1) / (N_PERM + 1)))

# (b) 限定到样本最多的单一部位
top <- bt[, .N, by = tissue][order(-N)][1, tissue]
sub <- bt[tissue == top]
cat(sprintf("  限定到 %s（n=%d，R=%d/NR=%d）: ", top, nrow(sub),
            sum(sub$response == "R"), sum(sub$response == "NR")))
if (uniqueN(sub$response) == 2 && min(table(sub$response)) >= 2) {
  d2 <- mean(sub[response == "NR", score]) - mean(sub[response == "R", score])
  cat(sprintf("Δ = %+.3f\n", d2))
} else cat("单一部位内不足以比较\n")

# (c) 秩回归：score ~ response + tissue
bt[, r_score := rank(score)]
fit <- lm(r_score ~ response + tissue, data = bt)
cat("\n  秩回归 score ~ response + tissue:\n")
print(round(summary(fit)$coefficients, 3))

fwrite(info, file.path(MR, "80a_pozniak_sample_level_with_tissue.tsv"), sep = "\t")
cat("\n写出 80a\n")
