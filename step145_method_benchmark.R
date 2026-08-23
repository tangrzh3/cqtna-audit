## step145 -- S44：五种提名方法的正面基准。
##
## 规则见 SUPP_method_benchmark_spec.md（S44），冻结于提交 462b27e，
## 先于本脚本与任何基准结果。本脚本不新增任何判据。
##
## 问题：S42 的恒等式只对单变异 Wald 成立。IVW / weighted median /
## coloc / SMR-HEIDI 会不会给出与 outcome GWAS 不同的名单？
##
## 输出：145a_method_benchmark.tsv、145b_console.log

suppressMessages(library(cqtna))

arg <- commandArgs(trailingOnly = TRUE)
root <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
        if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else getwd()
setwd(root)

FDR <- 0.05; LOCUS_KB <- 1000; KNOWN_KB <- 1000; GWS <- 5e-8

fix_gene <- function(g) {
  g <- as.character(g); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad)); g
}

K <- as_cqtna_known(within(read.csv("landi2020_known_loci_grch38.csv",
                                    stringsAsFactors = FALSE),
                           source <- "Landi2020"), build = "GRCh38")

## score(): 给定 被检验集(chr,pos,gene) + 显著标记 + 每条的 outcome P，
## 返回该方法的归属、分解与机会校正量。分区在该方法自己的位置集上算。
score <- function(label, chr, pos, gene, sig, p_out, note = "") {
  keep <- is.finite(as.numeric(pos)) & as.numeric(pos) > 0 & !is.na(sig)
  chr <- as.character(chr)[keep]; pos <- as.numeric(pos)[keep]
  gene <- fix_gene(gene)[keep]; sig <- as.logical(sig)[keep]
  p_out <- as.numeric(p_out)[keep]
  ## as_cqtna_mr 用 p 与 fdr 判显著；这里直接给一个已经体现显著规则的伪 p，
  ## 使 fdr < 0.05 当且仅当该方法判其显著。分区只依赖位置，不受此影响。
  pseudo <- ifelse(sig, 1e-300, 1)
  d <- data.frame(record_id = seq_along(pos), gene = gene, chr = chr,
                  pos = pos, p = pseudo, stringsAsFactors = FALSE)
  mr <- as_cqtna_mr(d, locus_kb = LOCUS_KB, build = "GRCh38",
                    locus_method = "fixed_centre")
  a <- cqtna_attribution(mr, K, known_kb = KNOWN_KB, fdr = FDR,
                         known_from = "any_record")
  kv <- as.logical(a$locus_known); names(kv) <- names(a$locus_known)
  loc <- as.character(mr$locus)
  BT <- a$background_loci; BK <- a$background_known
  sl <- unique(loc[sig])
  fis <- function(SK, ST) if (ST <= 0) NA_real_ else
    stats::fisher.test(matrix(c(SK, ST - SK, BK - SK,
                                (BT - BK) - (ST - SK)), 2),
                       alternative = "greater")$p.value
  ## 每个显著位点的 outcome P = 该位点显著记录中的最小值
  minp <- vapply(sl, function(L) suppressWarnings(min(p_out[sig & loc == L],
                                                      na.rm = TRUE)),
                 numeric(1))
  gws <- is.finite(minp) & minp < GWS
  ST <- length(sl); SK <- sum(kv[sl])
  ST2 <- sum(!gws); SK2 <- sum(kv[sl][!gws])
  p_sig <- if (ST) SK / ST else NA_real_; p_bg <- BK / BT
  data.frame(
    method = label, note = note,
    n_tested_loci = BT, bg_known = BK, bg_share = round(p_bg, 4),
    n_sig_loci = ST, n_known = SK,
    fold = a$fold, fisher_p = a$fisher_p_one_sided,
    A = if (ST) (p_sig - p_bg) / (1 - p_bg) else NA_real_,
    n_gws = sum(gws),
    n_sig_excl_gws = ST2, n_known_excl_gws = SK2,
    fold_excl_gws = if (ST2) (SK2 / ST2) / p_bg else NA_real_,
    p_excl_gws = fis(SK2, ST2),
    stringsAsFactors = FALSE)
}

out <- list()

## ---- 1. 单变异 Wald（主口径） ------------------------------------------------
mel <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
mel$chr <- sub(":.*", "", mel$SNP); mel$pos <- as.numeric(sub(".*:", "", mel$SNP))
mel$fdr_bh <- stats::p.adjust(mel$pval, "BH")
out[[1]] <- score("single-variant Wald", mel$chr, mel$pos, mel$SYMBOL,
                  mel$fdr_bh < FDR, mel$pval,
                  "3,556 records; identity applies")

## ---- 2-4. IVW / IVW-MRE / weighted median（relaxed 集） ----------------------
sens <- read.delim("18_sensitivity_meta.tsv", stringsAsFactors = FALSE)
rel <- read.delim("17_relaxed_instruments_meta.tsv", stringsAsFactors = FALSE)
look <- read.delim("144a_relaxed_outcome_lookup.tsv", stringsAsFactors = FALSE)
rel$key <- paste0(rel$chr, ":", as.integer(rel$pos))
rel$p_out <- look$pval_outcome[match(rel$key, look$key)]
## 代表位置 = 暴露侧 P 最小的工具（S44 §4：不得按 outcome 选）
rel <- rel[order(rel$exposure, rel$pval_exposure), ]
lead <- rel[!duplicated(rel$exposure), ]
## 该暴露的 outcome 最小 P（用于 GWS 分解），取其全部工具中的最小值
minp_by_exp <- tapply(rel$p_out, rel$exposure, function(x)
  suppressWarnings(min(x, na.rm = TRUE)))

for (m in c("Inverse variance weighted",
            "Inverse variance weighted (multiplicative random effects)",
            "Weighted median")) {
  s <- sens[sens$method == m, , drop = FALSE]
  s <- s[!duplicated(s$exposure), , drop = FALSE]
  i <- match(s$exposure, lead$exposure)
  ok <- !is.na(i)
  s <- s[ok, , drop = FALSE]; i <- i[ok]
  fdr <- stats::p.adjust(s$pval, "BH")
  lbl <- if (m == "Inverse variance weighted") "IVW"
         else if (grepl("multiplicative", m)) "IVW (MRE)" else "weighted median"
  out[[length(out) + 1]] <- score(
    lbl, lead$chr[i], lead$pos[i],
    if ("SYMBOL" %in% names(s)) s$SYMBOL else lead$gene_id[i],
    fdr < FDR, unname(minp_by_exp[s$exposure]),
    sprintf("%d exposures, relaxed instruments", nrow(s)))
}

## ---- 5. coloc（条件于 MR p < 0.05） ------------------------------------------
co <- read.delim("14_coloc_meta_results.tsv", stringsAsFactors = FALSE)
## 14 表不带坐标，按 exposure 键回到 13_meta 取位置与 outcome P
i <- match(co$exposure, mel$exposure)
co$chr <- mel$chr[i]; co$pos <- mel$pos[i]; co$p_out <- mel$pval[i]
okc <- !is.na(co$pos)
out[[length(out) + 1]] <- score(
  "coloc (PP.H4 > 0.8)", co$chr[okc], co$pos[okc], co$SYMBOL[okc],
  as.numeric(co$PP.H4)[okc] > 0.8, co$p_out[okc],
  "conditional on MR p < 0.05; downstream filter")
out[[length(out) + 1]] <- score(
  "coloc (PP.H4 > 0.5)", co$chr[okc], co$pos[okc], co$SYMBOL[okc],
  as.numeric(co$PP.H4)[okc] > 0.5, co$p_out[okc],
  "sensitivity threshold")

## ---- 6. SMR / HEIDI（条件于 MR p < 0.05） ------------------------------------
sm <- read.delim("15_SMR_meta_results.tsv", stringsAsFactors = FALSE)
## 15 表的 profile 列只是 cell/timepoint，完整键是 gene|profile
i <- match(paste0(sm$gene, "|", sm$profile), mel$exposure)
sm$chr <- mel$chr[i]; sm$pos <- mel$pos[i]; sm$p_out <- mel$pval[i]
oks <- !is.na(sm$pos) & is.finite(as.numeric(sm$p_SMR))
fdr_smr <- rep(NA_real_, nrow(sm))
fdr_smr[oks] <- stats::p.adjust(as.numeric(sm$p_SMR)[oks], "BH")
out[[length(out) + 1]] <- score(
  "SMR + HEIDI", sm$chr[oks], sm$pos[oks], sm$SYMBOL[oks],
  fdr_smr[oks] < FDR & as.numeric(sm$p_HEIDI)[oks] > 0.05, sm$p_out[oks],
  "conditional on MR p < 0.05; downstream filter")

res <- do.call(rbind, out)
write.table(res, "145a_method_benchmark.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

cat(sprintf("%-22s %7s %7s %6s %6s %9s %7s | %5s %7s %10s\n",
            "method", "tested", "sig", "known", "fold", "P", "A",
            "n_GWS", "fold_ex", "P_excl_GWS"))
for (i in seq_len(nrow(res)))
  cat(sprintf("%-22s %7d %7d %6d %6s %9s %7s | %5d %7s %10s\n",
              res$method[i], res$n_tested_loci[i], res$n_sig_loci[i],
              res$n_known[i],
              ifelse(is.finite(res$fold[i]), sprintf("%.2f", res$fold[i]), "-"),
              ifelse(is.finite(res$fisher_p[i]), sprintf("%.3g", res$fisher_p[i]), "-"),
              ifelse(is.finite(res$A[i]), sprintf("%.3f", res$A[i]), "-"),
              res$n_gws[i],
              ifelse(is.finite(res$fold_excl_gws[i]), sprintf("%.2f", res$fold_excl_gws[i]), "-"),
              ifelse(is.finite(res$p_excl_gws[i]), sprintf("%.3g", res$p_excl_gws[i]), "-")))
cat("\nwrote 145a\n")
