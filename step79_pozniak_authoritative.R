#!/usr/bin/env Rscript
# ============================================================================
# Step 79  Pozniak 复现（定论版）——改用作者自带的细胞类型注释
#
# ⚠ 本步取代 Step 77。原因是流程错误，如实记录：
#   预注册 §3 写明"**优先**使用该数据集自带的细胞类型注释，注释粒度不足才用
#   多基因评分"。Step 77 只检查了 Entire_TME.rds（仅有 seurat_clusters）便跳到
#   备选方案，未检查 GC_all_immune_sharing.rds —— 后者含作者注释的
#   low_resolution_clusters（CD4_Tcells、Tregs、CD8_Tcells…）以及
#   Response / Timepoint / patient_ID / Tissue。
#   Step 77 的结果保留在 77a/77b，作为"备选定义"的敏感性分析，不作为主结果。
#
# 主群体（预注册的首选）：low_resolution_clusters == "CD4_Tcells"
# 敏感性：加入 Tregs（发现队列的固定 CD4 群体亦含调节性细胞）
#
# 统计与 Step 74/77 逐位一致：样本为单位、≥20 细胞、16 基因 z 后均值、
# 置换检验 + Wilcoxon、去 TPI1 版本、重复患者三种处理。
# ============================================================================

suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
MR <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"; P <- "D:/Downloads/Pozniak"
set.seed(1); N_PERM <- 20000; MIN_CELLS <- 20L
LOCKED <- c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
            "ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA")

o <- readRDS(file.path(P, "GC_all_immune_sharing.rds"))
md <- as.data.table(o@meta.data, keep.rownames = "cell")
cat("免疫细胞对象:", ncol(o), "细胞 |", nrow(o), "基因\n")
print(md[, .N, by = low_resolution_clusters][order(-N)])

Xr <- GetAssayData(o, assay = "RNA", layer = "counts")
if (!nrow(Xr)) Xr <- GetAssayData(o, assay = "RNA", layer = "data")
cs <- Matrix::colSums(Xr[, sample(ncol(Xr), 300)])
cat(sprintf("\ncounts 层: 每细胞和中位 %.0f (CV %.2f), 最大 %.0f → ", median(cs), sd(cs)/mean(cs), max(Xr[,1:200])))
X <- log1p(Matrix::t(Matrix::t(Xr) / Matrix::colSums(Xr)) * 1e4)
cat("按库大小归一化 + log1p\n")
gl <- intersect(LOCKED, rownames(X)); cat("signature 可得:", length(gl), "/16\n")
zs <- function(v){s<-sd(v); if(!is.finite(s)||s==0) return(v*0); (v-mean(v))/s}

# ---- 纯度诊断（发现⑦的纪律：任何细胞群体都要报纯度）
cat("\n===== 群体纯度 =====\n")
for (nm in c("CD4_Tcells","Tregs","CD8_Tcells")) {
  i <- which(md$low_resolution_clusters == nm)
  if (!length(i)) next
  cat(sprintf("  %-12s n=%5d | CD4 %.1f%% | CD8A %.1f%% | CD8B %.1f%% | CD8双阴 %.1f%% | CD3E %.1f%%\n",
      nm, length(i), 100*mean(X["CD4",i]>0), 100*mean(X["CD8A",i]>0),
      100*mean(X["CD8B",i]>0), 100*mean(X["CD8A",i]==0 & X["CD8B",i]==0),
      100*mean(X["CD3E",i]>0)))
}

run <- function(cells, defname) {
  m <- md[cells]
  A <- t(as.matrix(X[gl, m$cell, drop = FALSE]))
  agg <- as.data.table(A)[, lapply(.SD, mean), by = .(sample = m$sample_ID)]
  info <- unique(m[, .(sample = sample_ID, patient = patient_ID,
                       tp = Timepoint, response = Response)])
  info <- info[m[, .N, by = .(sample = sample_ID)], on = "sample"][N >= MIN_CELLS]
  agg <- agg[match(info$sample, sample)]
  cat(sprintf("\n[%s] 样本 %d | 患者 %d\n", defname, nrow(info), uniqueN(info$patient)))
  print(dcast(info[, .N, by = .(tp, response)], tp ~ response, value.var = "N"))
  out <- list()
  for (t0 in c("BT","OT")) for (gs in list(list(gl,"locked16"), list(setdiff(gl,"TPI1"),"locked15_noTPI1")))
    for (mode in c("all","one_per_patient","single_sample_patients")) {
      sub <- info[tp == t0]
      idx <- switch(mode, all = seq_len(nrow(sub)),
        one_per_patient = sub[, .I[sample(.N,1)], by = patient]$V1,
        single_sample_patients = which(!sub$patient %in% sub[, .N, by=patient][N>1, patient]))
      s2 <- sub[idx]
      if (uniqueN(s2$response) < 2 || min(table(s2$response)) < 3) next
      M <- as.matrix(agg[match(s2$sample, agg$sample), gs[[1]], with = FALSE])
      sc <- rowMeans(apply(M, 2, zs)); y <- s2$response
      nr <- sc[y == "NR"]; r <- sc[y == "R"]; dd <- mean(nr) - mean(r)
      bs <- replicate(4000, mean(sample(nr,replace=TRUE)) - mean(sample(r,replace=TRUE)))
      pp <- sum(replicate(N_PERM, {yp <- sample(y); mean(sc[yp=="NR"]) - mean(sc[yp=="R"])}) >= dd)
      out[[length(out)+1]] <- data.table(definition = defname, timepoint = t0,
        genes = gs[[2]], subset = mode, n_R = length(r), n_NR = length(nr),
        score_diff = dd, ci_lo = quantile(bs,.025), ci_hi = quantile(bs,.975),
        wilcox_p = suppressWarnings(wilcox.test(nr, r))$p.value,
        perm_p = (pp+1)/(N_PERM+1))
    }
  rbindlist(out)
}

res <- rbind(
  run(which(md$low_resolution_clusters == "CD4_Tcells"), "CD4_Tcells (作者注释, 主)"),
  run(which(md$low_resolution_clusters %in% c("CD4_Tcells","Tregs")), "CD4_Tcells + Tregs (敏感性)"))

cat("\n===== 主检验（定论）=====\n")
print(res[, .(definition, timepoint, genes, subset, n_R, n_NR,
              d = round(score_diff,3), lo = round(ci_lo,3), hi = round(ci_hi,3),
              wilcox = signif(wilcox_p,3), perm = signif(perm_p,3))])
cat("\n方向汇总（预设：NR > R）：\n")
for (t0 in c("BT","OT")) {
  s <- res[timepoint == t0]
  cat(sprintf("  %s: %d/%d 为预期方向 | 置换 p<0.05: %d/%d\n",
              t0, sum(s$score_diff > 0), nrow(s), sum(s$perm_p < .05), nrow(s)))
}
fwrite(res, file.path(MR, "79a_pozniak_authoritative.tsv"), sep = "\t")
cat("\n写出 79a\n")
