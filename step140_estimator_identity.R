## step140 -- 这个估计量的 z 究竟是不是 outcome 的 z？逐格核验，并量化后果。
##
## Methods 已经写过这件事（"z = beta_out/se_out does not involve the exposure"），
## 但只用来解释 profile 特异性。若要把它接到主结论上，先得确认三件事：
##
##   1. 恒等式在**每一个**格子里都精确成立，不只是 melanoma x CD4；
##   2. BH 阈值换算成 outcome 的 P 是多少 —— 决定"发现"发生在什么显著性水平；
##   3. 显著位点里有多少本身就是 outcome 的全基因组显著信号 ——
##      若多数是，那么"提名"根本没有越过 GWAS 已知的东西。
##
## 输出：140a_identity_check.tsv、140b_threshold_and_gws.tsv、140c_console.log

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

## ---- 1. 恒等式：se == se_out/|beta_exp| 且 p == 2*pnorm(-|z_out|) ----------
## 只有带 beta_outcome/se_outcome/beta_exposure 的表能查，逐个查能查的。
cands <- list(
  list(nm = "melanoma x Soskic_CD4 (meta)", f = "13_meta_locus_annotation.tsv",
       bo = "beta_outcome", so = "se_outcome", bx = "beta_exposure",
       b = "b", se = "se", p = "pval"),
  list(nm = "melanoma x Soskic_CD4 (R12 single)", f = "04_MR_results_strict_all.tsv",
       bo = "beta.outcome", so = "se.outcome", bx = "beta.exposure",
       b = "b", se = "se", p = "pval"),
  list(nm = "melanoma x eQTLGen_blood", f = "92b_mr_eqtlgen_melanoma.tsv",
       bo = "beta.outcome", so = "se.outcome", bx = "beta.exposure",
       b = "b", se = "se", p = "pval"),
  list(nm = "HCC_high x Soskic_CD4", f = "84b_mr_HCC_high_GCST90809296.tsv",
       bo = "beta.outcome", so = "se.outcome", bx = "beta.exposure",
       b = "b", se = "se", p = "pval"),
  list(nm = "HCC_low x Soskic_CD4", f = "84b_mr_HCC_low_FinnGenR12.tsv",
       bo = "beta.outcome", so = "se.outcome", bx = "beta.exposure",
       b = "b", se = "se", p = "pval"))

ident <- list()
cat("1. is the Wald z the outcome z, in every cell that carries both sides?\n")
for (cd in cands) {
  if (!file.exists(cd$f)) {
    cat(sprintf("   %-38s FILE MISSING\n", cd$nm)); next
  }
  d <- read.delim(cd$f, stringsAsFactors = FALSE)
  need <- c(cd$bo, cd$so, cd$bx, cd$b, cd$se, cd$p)
  if (!all(need %in% names(d))) {
    cat(sprintf("   %-38s columns absent: %s\n", cd$nm,
                paste(setdiff(need, names(d)), collapse = ",")))
    next
  }
  bo <- as.numeric(d[[cd$bo]]); so <- as.numeric(d[[cd$so]])
  bx <- as.numeric(d[[cd$bx]]); b <- as.numeric(d[[cd$b]])
  se <- as.numeric(d[[cd$se]]); p <- as.numeric(d[[cd$p]])
  ok <- is.finite(bo) & is.finite(so) & is.finite(bx) & is.finite(se) & is.finite(p)
  d_se <- max(abs(se[ok] - so[ok] / abs(bx[ok])))
  z_out <- bo / so
  ## p 在极小处会下溢，只在可表示范围内比较
  cmp <- ok & p > 1e-290
  d_p <- max(abs(p[cmp] - 2 * pnorm(-abs(z_out[cmp]))))
  d_z <- max(abs((b[ok] / se[ok]) - z_out[ok]))
  ident[[length(ident) + 1]] <- data.frame(
    cell = cd$nm, n = sum(ok),
    max_abs_diff_se = d_se, max_abs_diff_z = d_z, max_abs_diff_p = d_p,
    identity_holds = (d_se < 1e-12 && d_z < 1e-9 && d_p < 1e-12),
    stringsAsFactors = FALSE)
  cat(sprintf("   %-38s n=%5d  |se-se_out/|bx|| %.2e  |z_MR-z_out| %.2e  |p-2F(-|z|)| %.2e  %s\n",
              cd$nm, sum(ok), d_se, d_z, d_p,
              if (d_se < 1e-12 && d_z < 1e-9 && d_p < 1e-12) "IDENTICAL" else "differs"))
}
ID <- do.call(rbind, ident)
write.table(ID, "140a_identity_check.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## ---- 2 & 3. BH 阈值换算成 outcome P；显著位点有多少本就 GWS ---------------
cat("\n2. what outcome P does the BH threshold correspond to, and\n")
cat("3. how many significant loci are already genome-wide significant outcome signals?\n\n")

mel <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
mel$chr <- sub(":.*", "", mel$SNP); mel$pos <- as.numeric(sub(".*:", "", mel$SNP))
d <- data.frame(record_id = seq_len(nrow(mel)), gene = fix_gene(mel$SYMBOL),
                chr = mel$chr, pos = mel$pos, p = mel$pval, stringsAsFactors = FALSE)
keep <- is.finite(d$p) & is.finite(d$pos) & d$pos > 0
d <- d[keep, ]; mel <- mel[keep, ]
K <- as_cqtna_known(within(read.csv("landi2020_known_loci_grch38.csv",
                                    stringsAsFactors = FALSE),
                           source <- "Landi2020"), build = "GRCh38")
mr <- as_cqtna_mr(d, locus_kb = LOCUS_KB, build = "GRCh38",
                  locus_method = "fixed_centre")
a <- cqtna_attribution(mr, K, known_kb = KNOWN_KB, fdr = FDR,
                       known_from = "any_record")
sel <- mr$fdr < FDR
thr <- max(mr$p[sel])            # BH 实际生效的 P 阈值 = outcome 的 P 阈值
kv <- as.logical(a$locus_known); names(kv) <- names(a$locus_known)
sig_loci <- unique(as.character(mr$locus)[sel])

## 每个显著位点：它最强记录的 outcome P，以及是否 < 5e-8
rows <- list()
for (L in sig_loci) {
  ix <- which(as.character(mr$locus) == L & sel)
  best <- ix[which.min(mr$p[ix])]
  rows[[length(rows) + 1]] <- data.frame(
    locus = L, n_sig_records = length(ix),
    gene = paste(head(sort(unique(mr$gene[ix])), 3), collapse = ","),
    min_outcome_p = mr$p[best],
    genome_wide_sig = mr$p[best] < GWS,
    on_known_locus = unname(kv[L]), stringsAsFactors = FALSE)
}
R <- do.call(rbind, rows)
R <- R[order(R$min_outcome_p), ]
cat(sprintf("   BH threshold in outcome-P terms: P <= %.3g  (%d significant records)\n\n",
            thr, sum(sel)))
cat("   locus                       n  gene                 outcome P   GWS?  known?\n")
for (i in seq_len(nrow(R)))
  cat(sprintf("   %-26s %2d  %-18s %10.3g  %-5s %s\n", R$locus[i], R$n_sig_records[i],
              substr(R$gene[i], 1, 18), R$min_outcome_p[i],
              ifelse(R$genome_wide_sig[i], "yes", "no"),
              ifelse(R$on_known_locus[i], "yes", "no")))
cat(sprintf("\n   %d of %d significant loci reach genome-wide significance on their own\n",
            sum(R$genome_wide_sig), nrow(R)))
cat(sprintf("   %d of %d are on a previously reported locus\n",
            sum(R$on_known_locus), nrow(R)))
tb <- table(gws = R$genome_wide_sig, known = R$on_known_locus)
cat("\n   cross-tabulation (genome-wide significant x previously reported):\n")
print(tb)

write.table(R, "140b_threshold_and_gws.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nwrote 140a / 140b\n")
