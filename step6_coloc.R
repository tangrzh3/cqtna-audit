# ============================================================
# Step 6: 共定位分析 (coloc.abf)
# 候选: 06_locus_annotation.tsv 中名义 p<0.05 的 gene×profile
#   - "潜在新位点" = 主要目标
#   - "已知位点"   = 阳性对照(MC1R/PARP1等,应当跑出高PP.H4)
# ============================================================

library(arrow); library(data.table); library(dplyr); library(coloc)
setwd("D:/R_ex/MR")

PARQUET_DIR <- "D:/Downloads/CD4_eqtl_step1_clean"

# ---- 各profile样本量(由 ma_count/(2*MAF) 反推,8个profile均100%一致) ----
PROFILE_N <- c(
  CD4_Naive_uns_0h   = 99,  CD4_Naive_stim_16h  = 99,
  CD4_Naive_stim_40h = 89,  CD4_Naive_stim_5d   = 85,
  CD4_Memory_uns_0h  = 100, CD4_Memory_stim_16h = 95,
  CD4_Memory_stim_40h= 89,  CD4_Memory_stim_5d  = 90
)

# ---- FinnGen R12 C3_MELANOMA_SKIN_EXALLC (官方manifest) ----
FG_CASES <- 5753; FG_CONTROLS <- 378749
FG_N <- FG_CASES + FG_CONTROLS          # 384502
FG_S <- FG_CASES / FG_N                 # 0.014962

WINDOW <- 500000   # 与parquet的cis窗口一致
MIN_SNP <- 50      # 区域内可用SNP少于此数则跳过

# ============================================================
# 1. 候选清单
# ============================================================
annot <- fread("06_locus_annotation.tsv")
cand <- annot[pval < 0.05]
cand[, profile := sub("^.*\\|", "", exposure)]
cat("候选 gene×profile:", nrow(cand),
    "| 潜在新位点:", sum(cand$category == "潜在新位点"),
    "| 已知位点(阳性对照):", sum(cand$category != "潜在新位点"), "\n")

# ============================================================
# 2. 结局数据(常驻内存,按染色体建key)
# ============================================================
if (!exists("outcome_raw")) {
  outcome_raw <- fread("finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz")
  setnames(outcome_raw, "#chrom", "chrom")
}
outcome_raw[, chrom := as.character(chrom)]
setkey(outcome_raw, chrom, pos)

# ============================================================
# 3. 按profile批量取cis窗口,逐基因跑coloc
# ============================================================
run_one <- function(eq, gwas, n_eqtl) {
  m <- merge(eq, gwas, by = "pos", suffixes = c(".e", ".g"))
  if (nrow(m) == 0) return(NULL)

  # 等位基因一致性:同向保留,反向翻转GWAS的beta
  m[, `:=`(EA = toupper(effect_allele), OA = toupper(other_allele),
           A = toupper(alt), R = toupper(ref))]
  m[, dir := fifelse(OA == R & EA == A,  1,
              fifelse(OA == A & EA == R, -1, NA_real_))]
  m <- m[!is.na(dir)]
  m[, beta_g := beta.g * dir]

  m <- m[!is.na(beta.e) & !is.na(se.e) & se.e > 0 &
         !is.na(beta_g) & !is.na(sebeta) & sebeta > 0 &
         eaf > 0 & eaf < 1]
  m <- unique(m, by = "pos")
  if (nrow(m) < MIN_SNP) return(NULL)

  m[, maf := pmin(eaf, 1 - eaf)]
  m <- m[maf > 0.001]
  if (nrow(m) < MIN_SNP) return(NULL)

  snpid <- as.character(m$pos)
  d1 <- list(type = "quant", beta = m$beta.e, varbeta = m$se.e^2,
             MAF = m$maf, N = n_eqtl, snp = snpid)
  d2 <- list(type = "cc",    beta = m$beta_g, varbeta = m$sebeta^2,
             MAF = m$maf, N = FG_N, s = FG_S, snp = snpid)

  r <- tryCatch(suppressWarnings(coloc.abf(d1, d2)), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  s <- as.list(r$summary)
  data.table(nsnps = s$nsnps, PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf,
             PP.H2 = s$PP.H2.abf, PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf)
}

res_list <- list(); k <- 0
for (prof in unique(cand$profile)) {
  sub <- cand[profile == prof]
  fp  <- file.path(PARQUET_DIR, paste0(prof, "_step1_clean.parquet"))
  if (!file.exists(fp)) { warning("缺文件: ", fp); next }
  cat("\n[", prof, "] 候选基因", uniqueN(sub$gene_id), "个 ... 读取cis窗口\n")

  eq_all <- open_dataset(fp) %>%
    filter(gene_id %in% sub$gene_id) %>%
    select(gene_id, chr, pos, other_allele, effect_allele, eaf, beta, se, pval) %>%
    collect() %>% as.data.table()
  setnames(eq_all, c("beta","se","pval"), c("beta.e","se.e","pval.e"))

  for (i in seq_len(nrow(sub))) {
    g <- sub$gene_id[i]
    eq <- eq_all[gene_id == g]
    if (nrow(eq) < MIN_SNP) next
    ch <- as.character(eq$chr[1])
    lo <- min(eq$pos); hi <- max(eq$pos)
    gwas <- outcome_raw[.(ch)][pos >= lo & pos <= hi,
              .(pos, ref, alt, beta.g = beta, sebeta, pval.g = pval, rsids)]
    if (nrow(gwas) < MIN_SNP) next

    out <- run_one(eq, gwas, PROFILE_N[[prof]])
    if (is.null(out)) next
    k <- k + 1
    res_list[[k]] <- cbind(sub[i, .(gene_id, SYMBOL, exposure, cell_type, timepoint,
                                    category, near_locus, dist_kb,
                                    MR_OR = OR, MR_pval = pval, MR_FDR = FDR)], out)
    if (k %% 25 == 0) cat("  已完成", k, "\n")
  }
}

coloc_res <- rbindlist(res_list)
coloc_res[, PP4_ratio := PP.H4 / (PP.H3 + PP.H4)]
setorder(coloc_res, -PP.H4)
fwrite(coloc_res, "07_coloc_results.tsv", sep = "\t")

# ============================================================
# 4. 汇总
# ============================================================
cat("\n\n================ coloc 汇总 ================\n")
cat("成功跑完:", nrow(coloc_res), "/", nrow(cand), "\n\n")
cat("--- 按位点类别 ---\n")
print(coloc_res[, .(N = .N,
                    PP4_gt_0.7 = sum(PP.H4 > 0.7),
                    ratio_gt_0.7 = sum(PP4_ratio > 0.7),
                    median_PP4 = round(median(PP.H4), 3)), by = category])
cat("\n--- 阳性对照(已知位点)中PP.H4最高的10个 ---\n")
print(head(coloc_res[category != "潜在新位点",
     .(SYMBOL, cell_type, timepoint, near_locus, PP.H4 = round(PP.H4,3),
       PP.H3 = round(PP.H3,3), nsnps)], 10))
cat("\n--- 潜在新位点中PP.H4>0.7的候选 ---\n")
print(coloc_res[category == "潜在新位点" & PP.H4 > 0.7,
     .(SYMBOL, cell_type, timepoint, MR_OR = round(MR_OR,3), MR_pval,
       PP.H4 = round(PP.H4,3), PP4_ratio = round(PP4_ratio,3), nsnps)])
