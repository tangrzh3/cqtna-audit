#!/usr/bin/env Rscript
# ============================================================================
# Step 65  TPI1 工具变量的 motif disruption
#
# 起因：Step 62/64 已确认该变异所在元件（chr12:6,867,158–6,867,932）是
#   **组成型开放、活化不变（act log2FC=+0.047）、轴不变（mean lfc=−0.013,
#   4 套里 2 套同向, consistent=FALSE）**。
#   故机制不是"元件开关"，若有 TF 参与，只能是组成型开放元件内的**等位特异结合**。
#   本脚本检验该变异是否改变某个 TF 的结合基序。
#
# *** 跑之前写死（勿在看到结果后修改）***
#
# 主假设（方向由已有结果预先指定，非事后挑选）：
#   轴的高端富集 AP-1 家族（23/23 motif 在活化不变 peak 上保留，中位 odds 1.52），
#   TPI1 位于高端。故预先指定的检验是：
#      **该变异对 AP-1 家族 motif 的扰动，是否超过同区域随机变异的背景分布。**
#
# 判读（预设）：
#   命中 = AP-1 家族的最大 |Δscore| 落在背景分布的**前 10%**（经验 p ≤ 0.10）
#   未命中 = 报告为"无证据"，**不得**改换家族、改换阈值或改用其他 motif 重新叙述
#
# 背景（预设）：同一 ±500 kb 区域内的其他双等位 SNP（1000G EUR 面板内），
#   随机取 500 个，种子固定。选同区域是为了匹配序列组成与可及性背景。
#
# ⚠ 本检验只能说明"变异是否改变基序打分"，**不能**证明等位特异结合，
#   后者需 reads 级数据 + 基因型。阴性不得用"静息态/外周血"之类理由解释。
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(TFBSTools); library(JASPAR2020)
  library(motifmatchr); library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings)
})

MR   <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"
GEN  <- BSgenome.Hsapiens.UCSC.hg38
CHR  <- "chr12"; POS <- 6867132L; REF <- "C"; ALT <- "T"
FLANK <- 30L                 # 覆盖 JASPAR 最长基序
N_BG  <- 500L
set.seed(1)

pfm <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))
cat("JASPAR2020 CORE (human) motifs:", length(pfm), "\n")
mnames <- sapply(pfm, name)
AP1 <- grepl("JUN|FOS|BATF|JDP2", mnames, ignore.case = TRUE)
cat("AP-1 家族 motif 数:", sum(AP1), " ->", paste(head(mnames[AP1], 12), collapse = ", "), "\n\n")

# ---------------------------------------------------------------- 背景变异
pv <- fread(file.path(MR, "regions/TPI1_tmp.pvar"), skip = "#CHROM")
setnames(pv, c("#CHROM", "POS", "ID", "REF", "ALT"), c("CHROM", "POS", "ID", "REF", "ALT"),
         skip_absent = TRUE)
pv <- pv[nchar(REF) == 1L & nchar(ALT) == 1L]
TARGET_POS <- POS
bg <- pv[POS != TARGET_POS]
bg <- bg[sample(.N, min(N_BG, .N))]
cat("背景变异:", nrow(bg), "个\n")

vars <- rbind(data.table(POS = POS, REF = REF, ALT = ALT, grp = "target"),
              bg[, .(POS, REF, ALT, grp = "background")])

# ---------------------------------------------------------------- 构建序列
mkseq <- function(pos, allele) {
  s <- as.character(getSeq(GEN, CHR, pos - FLANK, pos + FLANK))
  substr(s, FLANK + 1L, FLANK + 1L) <- allele
  s
}
# 参考碱基核查（必须通过，否则坐标或 build 有误）
ref_base <- as.character(getSeq(GEN, CHR, POS, POS))
cat("基因组参考碱基 @", CHR, POS, "=", ref_base, " | pvar REF =", REF, "
")
stopifnot(identical(ref_base, REF))

ref_seqs <- DNAStringSet(vapply(seq_len(nrow(vars)),
                                function(i) mkseq(vars$POS[i], vars$REF[i]), ""))
alt_seqs <- DNAStringSet(vapply(seq_len(nrow(vars)),
                                function(i) mkseq(vars$POS[i], vars$ALT[i]), ""))
names(ref_seqs) <- names(alt_seqs) <- paste0(vars$grp, "_", vars$POS)

# ---------------------------------------------------------------- 打分
score_of <- function(seqs) {
  m <- matchMotifs(pfm, seqs, out = "scores")
  as.matrix(motifScores(m))          # 行=序列, 列=motif, 值=该序列内最高分
}
cat("打分中（REF）...\n"); S_ref <- score_of(ref_seqs)
cat("打分中（ALT）...\n"); S_alt <- score_of(alt_seqs)
D <- S_ref - S_alt                    # 正 = REF 结合更强 = ALT 破坏基序
colnames(D) <- mnames

i_t <- which(vars$grp == "target")
i_b <- which(vars$grp == "background")

# ---------------------------------------------------------------- 主检验
ap1_target <- max(abs(D[i_t, AP1]))
ap1_bg     <- apply(abs(D[i_b, AP1, drop = FALSE]), 1, max)
emp_p      <- mean(ap1_bg >= ap1_target)
cat(sprintf("\n===== 预设主检验 =====\nAP-1 家族最大 |Δscore|: 目标 %.3f | 背景中位 %.3f\n经验 p = %.3f  ->  %s\n",
            ap1_target, median(ap1_bg), emp_p,
            ifelse(emp_p <= 0.10, "命中（前 10%）", "未命中 -> 报告为无证据")))

# 全 motif 排名（描述性）
dt <- data.table(motif = mnames, delta = as.numeric(D[i_t, ]),
                 ref_score = as.numeric(S_ref[i_t, ]),
                 alt_score = as.numeric(S_alt[i_t, ]),
                 is_AP1 = AP1)
dt[, abs_delta := abs(delta)]
setorder(dt, -abs_delta)
dt[, rank := .I]
# 每个 motif 的经验 p（相对同一 motif 的背景分布）
bgabs <- abs(D[i_b, , drop = FALSE])
dt[, emp_p_motif := sapply(seq_len(.N), function(k) {
  j <- which(mnames == motif[k])[1]
  mean(bgabs[, j] >= abs_delta[k])
})]
cat("\n===== 扰动最大的 12 个 motif（描述性）=====\n")
print(dt[1:12, .(rank, motif, ref_score = round(ref_score, 2),
                 alt_score = round(alt_score, 2), delta = round(delta, 3),
                 is_AP1, emp_p_motif)])
cat("\n===== AP-1 家族成员的排名 =====\n")
print(dt[is_AP1 == TRUE][1:min(10, .N),
        .(rank, motif, delta = round(delta, 3), emp_p_motif)])

fwrite(dt, file.path(MR, "65a_motif_disruption_TPI1.tsv"), sep = "\t")
fwrite(data.table(ap1_target_maxabs = ap1_target,
                  ap1_bg_median = median(ap1_bg),
                  ap1_bg_p90 = quantile(ap1_bg, .90),
                  empirical_p = emp_p, n_background = length(i_b)),
       file.path(MR, "65b_motif_disruption_test.tsv"), sep = "\t")
cat("\n写出 65a / 65b\n")
