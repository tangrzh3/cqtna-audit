## step141 -- 4.96 倍的富集是从哪里来的？按 outcome 显著性拆开。
##
## step140 证明了 |z_Wald| 恒等于 |z_outcome|，因此 BH 显著集就是
## "被测变异中 outcome P 最小的那一批"。于是有一个必须问的问题：
## 显著位点里有一部分**本身就是 outcome 的全基因组显著信号**——
## 那些位点当然会落在已发表的名单上，因为已发表的名单正是这么来的。
##
## 本脚本把这批位点逐档剔除，看剩下的部分还有没有富集。
## ⚠ 事后分解，不改动任何注册终点。
##
## 输出：141a_enrichment_decomposition.tsv、141b_console.log

suppressMessages(library(cqtna))
arg <- commandArgs(trailingOnly = TRUE)
setwd(if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
      if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else getwd())
fix_gene <- function(g){g<-as.character(g);b<-is.na(g)|!nzchar(trimws(g));g[b]<-paste0("u",which(b));g}
mel <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors=FALSE)
mel$chr <- sub(":.*","",mel$SNP); mel$pos <- as.numeric(sub(".*:","",mel$SNP))
d <- data.frame(record_id=seq_len(nrow(mel)), gene=fix_gene(mel$SYMBOL),
                chr=mel$chr,pos=mel$pos,p=mel$pval,stringsAsFactors=FALSE)
k <- is.finite(d$p)&is.finite(d$pos)&d$pos>0; d<-d[k,]
K <- as_cqtna_known(within(read.csv("landi2020_known_loci_grch38.csv",stringsAsFactors=FALSE),
                           source<-"Landi2020"), build="GRCh38")
mr <- as_cqtna_mr(d, locus_kb=1000, build="GRCh38", locus_method="fixed_centre")
a  <- cqtna_attribution(mr, K, known_kb=1000, fdr=0.05, known_from="any_record")
kv <- as.logical(a$locus_known); names(kv) <- names(a$locus_known)
sel <- mr$fdr < 0.05
BT <- a$background_loci; BK <- a$background_known
fg <- function(SK,ST) fisher.test(matrix(c(SK,ST-SK,BK-SK,(BT-BK)-(ST-SK)),2),alternative="greater")$p.value

R <- do.call(rbind, lapply(unique(as.character(mr$locus)[sel]), function(L){
  ix <- which(as.character(mr$locus)==L & sel)
  data.frame(locus=L, minp=min(mr$p[ix]), known=unname(kv[L]))}))
cat(sprintf("background: %d/%d known (%.1f%%)\n\n", BK, BT, 100*BK/BT))
out <- list()
for (thr in c(0, 5e-8, 1e-6, 1e-5, 1e-4)) {
  sub <- R[R$minp >= thr, ]           # 去掉本身达到该显著性的位点
  ST <- nrow(sub); SK <- sum(sub$known)
  out[[length(out) + 1]] <- data.frame(
    excluded_below_outcome_p = thr, n_loci = ST, n_known = SK,
    bg_known = BK, bg_loci = BT,
    fold = if (ST) (SK / ST) / (BK / BT) else NA_real_,
    fisher_p = if (ST) fg(SK, ST) else NA_real_, stringsAsFactors = FALSE)
  cat(sprintf("excluding loci with outcome P < %-7g : %d loci, %d known, fold %s, Fisher P %s\n",
      thr, ST, SK,
      if (ST) sprintf("%.2f", (SK / ST) / (BK / BT)) else "-",
      if (ST) sprintf("%.4g", fg(SK, ST)) else "-"))
}
write.table(do.call(rbind, out), "141a_enrichment_decomposition.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nreference: all loci                 : %d loci, %d known, fold %.2f, Fisher P %.4g\n",
            a$significant_loci, a$significant_known, a$fold, a$fisher_p_one_sided))
cat("\nwrote 141a\n")
