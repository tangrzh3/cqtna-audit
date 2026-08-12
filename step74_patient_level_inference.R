#!/usr/bin/env Rscript
# ============================================================================
# Step 74  回应审稿 R1-2：患者层面统计单位与相关结构
#
# 审稿意见（成立，本步不辩护，直接重做）：
#   (a) "13/15 基因方向一致" 的二项检验假设各基因方向近似独立。
#       糖酵解基因高度相关，故该 P 值过于乐观。
#   (b) Methods 允许同一患者的多个 patient_tp 样本作为统计单位。
#       实测：Post 时点有 4 个患者各贡献 2 个样本（P1/P3/P5/P23），
#       且 P1 与 P5 的两个样本**应答标签相反**。样本并不独立。
#
# 本步给出四样东西（全部预先指定于本文件，跑之前写死）：
#   1. 基因间相关结构——量化二项检验被违反的程度
#   2. 患者层面复合 signature score（锁定 16 基因；及去 TPI1 的 15 基因）
#      + 组间效应量与 bootstrap 置信区间
#   3. **保留基因相关结构的置换检验**：置换样本的 response 标签，
#      整条 signature 一起动，因此相关结构完全保留
#   4. 同一患者多样本的三种处理：全用 / 每患者取一个 / 排除重复患者
#
# ⚠ 判读纪律：若置换 P 值远大于二项 P 值，则**以置换 P 值为准**，
#   正文与摘要一并改写，不得保留二项数字作为主要证据。
# ============================================================================

suppressPackageStartupMessages({library(Seurat); library(data.table); library(Matrix)})
MR <- "D:/R_ex/MR"; set.seed(1); N_PERM <- 20000

LOCKED <- c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
            "ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA")

obj <- readRDS(file.path(MR, "24_checkpoint_GSE120575_annotated_obj.rds"))
cat("对象:", ncol(obj), "细胞\n")
md <- obj@meta.data
# 细胞集与 Step 26 完全一致：ID 取自旧 rds（该对象**只能取 ID，不可取表达值**），
# 表达值取自 checkpoint 对象。见 HANDOFF 第五节。
fixed_cd4 <- colnames(readRDS(file.path(MR, "GSE120575_CD4_clusters_1_5_6_12_15.rds")))
CD4 <- intersect(fixed_cd4, colnames(obj))
cat("固定 CD4 细胞:", length(CD4), "（应为 3878）
")
stopifnot(length(CD4) == 3878)
X <- GetAssayData(obj, layer = "data")[, CD4, drop = FALSE]
md <- md[CD4, ]

samp <- md[["patient_tp"]]; resp <- md[["response"]]
tp <- ifelse(grepl("^Pre", samp), "Pre", "Post")
gl <- intersect(LOCKED, rownames(X))
cat("signature 可得:", length(gl), "/", length(LOCKED), "\n")

# ---- 样本 × 基因 均值矩阵（统计单位 = 样本）
agg <- t(sapply(split(seq_along(samp), samp), function(i)
  Matrix::rowMeans(X[gl, i, drop = FALSE])))
info <- unique(data.table(sample = samp, tp = tp, response = resp,
                          n = as.integer(table(samp)[samp])))
info <- info[match(rownames(agg), sample)]
info[, patient := sub("^(Pre|Post)_", "", sub("_\\d+$", "", sample))]
info <- info[n >= 20]
agg <- agg[info$sample, , drop = FALSE]
cat("样本:", nrow(info), " | 患者:", uniqueN(info$patient), "\n")
cat("同一患者同一时点多样本:", info[, .N, by = .(tp, patient)][N > 1, .N], "例\n\n")

# ---- 1. 基因间相关结构
cat("===== 1. 基因间相关（样本层，说明二项检验为何过于乐观）=====\n")
for (t0 in c("Pre", "Post")) {
  A <- agg[info$tp == t0, , drop = FALSE]
  R <- cor(A, method = "spearman")
  off <- R[upper.tri(R)]
  cat(sprintf("%-5s 基因间 Spearman: 中位 %.3f, 四分位 [%.3f, %.3f], >0.5 的比例 %.0f%%\n",
              t0, median(off), quantile(off, .25), quantile(off, .75),
              100 * mean(off > .5)))
  # 有效独立基因数（Li & Ji 式特征值法）
  ev <- eigen(cor(A))$values; ev <- ev[ev > 0]
  meff <- 1 + (length(ev) - 1) * (1 - var(ev) / length(ev))
  cat(sprintf("      有效独立基因数 ≈ %.1f / %d  → 二项检验的 n 被高估约 %.1f 倍\n",
              meff, ncol(A), ncol(A) / meff))
}

# ---- 2/3/4
zs <- function(m) apply(m, 2, function(x) (x - mean(x)) / ifelse(sd(x) > 0, sd(x), 1))

run <- function(t0, genes, subset_desc, keep_idx) {
  ii <- which(info$tp == t0)[keep_idx]
  A <- agg[ii, genes, drop = FALSE]; y <- info$response[ii]
  if (length(unique(y)) < 2 || min(table(y)) < 3) return(NULL)
  sc <- rowMeans(zs(A))                       # 复合 signature score
  nr <- sc[y != "Responder"]; r <- sc[y == "Responder"]
  d <- mean(nr) - mean(r)
  w <- suppressWarnings(wilcox.test(nr, r))$p.value
  # bootstrap CI（对样本重抽）
  bs <- replicate(4000, {
    i1 <- sample(seq_along(nr), replace = TRUE); i2 <- sample(seq_along(r), replace = TRUE)
    mean(nr[i1]) - mean(r[i2])
  })
  # 置换：打乱 response 标签，整条 signature 一起动 → 相关结构完全保留
  obs_dir <- sum(colMeans(A[y != "Responder", , drop = FALSE]) >
                 colMeans(A[y == "Responder", , drop = FALSE]))
  pd <- pc <- 0L
  for (k in seq_len(N_PERM)) {
    yp <- sample(y)
    if (length(unique(yp)) < 2) next
    scp <- sc
    dp <- mean(scp[yp != "Responder"]) - mean(scp[yp == "Responder"])
    if (dp >= d) pc <- pc + 1L
    dirp <- sum(colMeans(A[yp != "Responder", , drop = FALSE]) >
                colMeans(A[yp == "Responder", , drop = FALSE]))
    if (dirp >= obs_dir) pd <- pd + 1L
  }
  data.table(timepoint = t0, genes = length(genes), subset = subset_desc,
             n_R = sum(y == "Responder"), n_NR = sum(y != "Responder"),
             score_diff = d, ci_lo = quantile(bs, .025), ci_hi = quantile(bs, .975),
             wilcox_p = w, perm_p_score = (pc + 1) / (N_PERM + 1),
             dir_obs = obs_dir, perm_p_direction = (pd + 1) / (N_PERM + 1))
}

out <- list()
for (t0 in c("Pre", "Post")) {
  n_t <- sum(info$tp == t0)
  sub_info <- info[tp == t0]
  # (i) 全部样本  (ii) 每患者随机取一个  (iii) 排除有重复的患者
  keep_all <- seq_len(n_t)
  keep_one <- sub_info[, .I[sample(.N, 1)], by = patient]$V1
  dup_pat <- sub_info[, .N, by = patient][N > 1, patient]
  keep_nodup <- which(!sub_info$patient %in% dup_pat)
  for (g in list(list(gl, "locked16"), list(setdiff(gl, "TPI1"), "locked15_noTPI1"))) {
    for (s in list(list(keep_all, "all samples"), list(keep_one, "one per patient"),
                   list(keep_nodup, "patients with a single sample"))) {
      r <- run(t0, g[[1]], paste0(g[[2]], " | ", s[[2]]), s[[1]])
      if (!is.null(r)) out[[length(out) + 1]] <- r
    }
  }
}
res <- rbindlist(out)
fwrite(res, file.path(MR, "74a_patient_level_inference.tsv"), sep = "\t")

cat("\n===== 2–4. 患者层面复合 score + 保留相关结构的置换 =====\n")
for (t0 in c("Pre", "Post")) {
  cat("\n--", t0, "--\n")
  r <- res[timepoint == t0]
  for (i in seq_len(nrow(r))) with(r[i], cat(sprintf(
    "%-42s n=%d/%d  Δscore=%+.3f [%.3f, %.3f]  Wilcoxon p=%.4f  置换 p(score)=%.4f  方向 %d/%d 置换 p=%.4f\n",
    subset, n_R, n_NR, score_diff, ci_lo, ci_hi, wilcox_p, perm_p_score,
    dir_obs, genes, perm_p_direction)))
}
cat("\n写出 74a\n")
