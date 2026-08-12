#!/usr/bin/env Rscript
# ============================================================================
# Step 70  TCGA-SKCM：bulk 生存分析作为**区室问题的演示**（不是验证）
#
# *** 跑之前写死，且这是本步存在的全部理由 ***
#
# 本步**不是**为了验证 TPI1 或糖酵解状态与预后有关。Step 68/69 已证明
# TPI1 在恶性细胞中比 CD4⁺T 高 6.7 倍（两队列、16/16 与 11/11 患者），
# 故 bulk 组织中该基因的信号约 98–99% 不来自 CD4⁺T。
#
# 预期：bulk TPI1 / 糖酵解签名与生存**显著相关**（肿瘤糖酵解与预后差是成熟结论）。
# **无论结果显著与否，一律报告为混杂，不得作为对 CD4 假说的支持。**
#
# 演示的逻辑：
#   (1) bulk 中该基因确实与生存相关（预期阳性）
#   (2) 但同一基因在 bulk 中与肿瘤含量/黑色素细胞标记强相关、与免疫标记弱或负相关
#   (3) 且细胞层面已知它主要来自恶性细胞
#   → 因此 (1) 是肿瘤生物学，不是 CD4 生物学。这正是发现⑤要说的。
#
# 对照（预设）：
#   PC  已知的免疫浸润标记（PTPRC/CD3E）应与生存**正**相关（黑色素瘤中免疫浸润
#       与较好预后相关，是成熟结论）。若该对照不成立，说明本队列或处理有问题，
#       本步无信息。
#
# 数据：UCSC Xena TCGA-SKCM HiSeqV2（log2(norm_count+1)）+ TCGA-CDR 生存表
# ============================================================================

suppressPackageStartupMessages({library(data.table); library(survival)})
SP  <- "C:/Users/a1197/AppData/Local/Temp/claude/D--------/63f94216-98ba-4480-9b61-1857f0a83e1b/scratchpad"
OUT <- "D:/R_ex/MR"

LOCKED <- c("SLC2A1","SLC2A3","HK1","HK2","GPI","PFKL","PFKP","PFKFB3",
            "ALDOA","TPI1","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA")
TUMOUR <- c("MLANA","PMEL","TYR","DCT","TYRP1","SOX10","MITF","PRAME")
IMMUNE <- c("PTPRC","CD3D","CD3E","CD2","CD8A","IL7R","LCK")

ex <- fread(file.path(SP, "TCGA.SKCM.sampleMap_HiSeqV2.gz"))
g  <- ex[[1]]; ex[, 1 := NULL]
m  <- as.matrix(ex); rownames(m) <- g; rm(ex); gc()
cat("表达:", nrow(m), "基因 x", ncol(m), "样本\n")

sv <- fread(file.path(SP, "survival_SKCM_survival.txt"))
sv <- sv[!is.na(OS) & !is.na(OS.time) & OS.time > 0]
common <- intersect(colnames(m), sv$sample)
m <- m[, common]; sv <- sv[match(common, sample)]
cat("可分析样本:", length(common), "| 事件数:", sum(sv$OS), "\n")

zs <- function(x) (x - mean(x)) / sd(x)
score <- function(gs) {
  gg <- intersect(gs, rownames(m))
  colMeans(t(apply(m[gg, , drop = FALSE], 1, zs)))
}
d <- data.table(sample = common, OS = sv$OS, time = sv$OS.time,
                TPI1 = zs(as.numeric(m["TPI1", ])),
                Glyco = score(LOCKED), GlycoNo = score(setdiff(LOCKED, "TPI1")),
                Tumour = score(TUMOUR), Immune = score(IMMUNE))

cox <- function(v, lab, extra = NULL) {
  f <- if (is.null(extra)) as.formula(sprintf("Surv(time, OS) ~ %s", v))
       else as.formula(sprintf("Surv(time, OS) ~ %s + %s", v, extra))
  s <- summary(coxph(f, data = d))
  hr <- s$coefficients[v, "exp(coef)"]; p <- s$coefficients[v, "Pr(>|z|)"]
  cat(sprintf("%-34s HR=%.3f  p=%.3g\n", lab, hr, p))
  data.table(variable = lab, HR = hr, p = p, adjusted = ifelse(is.null(extra), "-", extra))
}

cat("\n===== 预设对照：免疫浸润应与较好预后相关 =====\n")
pc <- cox("Immune", "PC  免疫浸润评分")
if (pc$HR < 1 && pc$p < 0.05) cat("  -> PASS，队列与处理正常\n") else
  cat("  -> FAIL，本步无信息，不解释下列结果\n")

cat("\n===== 演示第 (1) 步：bulk 中确实显著（预期阳性，非支持证据）=====\n")
r1 <- rbind(cox("TPI1", "bulk TPI1"),
            cox("Glyco", "bulk 糖酵解签名(16)"),
            cox("GlycoNo", "bulk 糖酵解签名(去TPI1)"))

cat("\n===== 演示第 (2) 步：该信号在 bulk 中跟随哪个区室 =====\n")
for (v in c("Tumour", "Immune")) {
  ct <- cor.test(d$TPI1, d[[v]], method = "spearman")
  cat(sprintf("bulk TPI1 ~ %-7s 评分   rho=%+.3f  p=%.3g\n", v,
              unname(ct$estimate), ct$p.value))
}
ct <- cor.test(d$Glyco, d$Tumour, method = "spearman")
cat(sprintf("bulk 糖酵解 ~ Tumour 评分   rho=%+.3f  p=%.3g\n",
            unname(ct$estimate), ct$p.value))

cat("\n===== 演示第 (3) 步：对肿瘤含量校正后还剩多少 =====\n")
r2 <- rbind(cox("TPI1", "bulk TPI1 | 校正 Tumour", "Tumour"),
            cox("TPI1", "bulk TPI1 | 校正 Tumour+Immune", "Tumour + Immune"),
            cox("Glyco", "bulk 糖酵解 | 校正 Tumour", "Tumour"))

fwrite(rbind(pc, r1, r2), file.path(OUT, "70a_TCGA_cox.tsv"), sep = "\t")
fwrite(d, file.path(OUT, "70b_TCGA_scores.tsv"), sep = "\t")
cat("\n写出 70a / 70b\n")
cat("\n★ 表述纪律：以上任何显著结果**不得**表述为对 CD4 假说的支持。\n")
