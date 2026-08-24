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

## ---------------------------------------------------------------------------
## 评审第三方 2026-08-23 提出的两条，加算于此。两条都不放松任何判据。
##
## R1.3 背景须与显著集同步条件化：只从显著集剔除已达 5e-8 的位点、
##      却把它们留在背景里，问的是一个混合问题。
## R1.4 相邻 bounded locus 未必是两个独立的已发表风险区：
##      16:88860636-89730161 与 16:89871237 落在**完全相同的七个 Landi lead SNP**
##      的 1 Mb 内，即在重复计同一个区域。
## ---------------------------------------------------------------------------
cat("\n", strrep("-", 72), "\n", sep = "")
kn_raw <- read.csv("landi2020_known_loci_grch38.csv", stringsAsFactors = FALSE)
loc <- as.character(mr$locus)
lv <- names(kv)
bg_minp <- tapply(mr$p, loc, min)[lv]
bg_gws  <- is.finite(bg_minp) & bg_minp < 5e-8

fis2 <- function(SK, ST, BK2, BT2)
  stats::fisher.test(matrix(c(SK, ST - SK, BK2 - SK, (BT2 - BK2) - (ST - SK)), 2),
                     alternative = "greater")$p.value

sub <- R[R$minp >= 5e-8, ]
ST <- nrow(sub); SK <- sum(sub$known)
BT2 <- sum(!bg_gws); BK2 <- sum(kv[!bg_gws])
cat(sprintf("R1.3 background also restricted : %d loci, %d known, bg %d/%d, fold %.2f, P %.4g\n",
            ST, SK, BK2, BT2, (SK / ST) / (BK2 / BT2), fis2(SK, ST, BK2, BT2)))

## R1.4：共享任一 Landi lead SNP 的 bounded locus 合并成一个区域，同一规则施于两侧
attr_of <- lapply(lv, function(L) {
  ch <- sub(":.*", "", L); rng <- sub(".*:", "", L)
  lo <- as.numeric(sub("-.*", "", rng)); hi <- as.numeric(sub(".*-", "", rng))
  if (is.na(hi)) hi <- lo
  kn_raw$rsid[kn_raw$chr == ch & kn_raw$pos >= lo - 1e6 & kn_raw$pos <= hi + 1e6]
})
names(attr_of) <- lv
par <- seq_along(lv); names(par) <- lv
find <- function(x) { while (par[[x]] != which(lv == x)) x <- lv[par[[x]]]; x }
for (s in unique(unlist(attr_of))) {
  g <- lv[vapply(attr_of, function(v) s %in% v, logical(1))]
  if (length(g) > 1) for (i in 2:length(g)) {
    a <- find(g[1]); b <- find(g[i])
    if (a != b) par[[a]] <- which(lv == b)
  }
}
reg <- vapply(lv, find, character(1)); names(reg) <- lv
reg_known <- tapply(kv, reg, any)
sig_reg <- unique(reg[unique(loc[mr$fdr < 0.05])])
STm <- length(sig_reg); SKm <- sum(reg_known[sig_reg])
BTm <- length(reg_known); BKm <- sum(reg_known)
cat(sprintf("R1.4 merge loci sharing a lead SNP: %d regions, %d known, bg %d/%d, fold %.2f, P %.4g\n",
            STm, SKm, BKm, BTm, (SKm / STm) / (BKm / BTm), fis2(SKm, STm, BKm, BTm)))
gws_reg <- unique(reg[lv[bg_gws]])
cat(sprintf("     the %d genome-wide-significant bounded loci are %d distinct published regions\n",
            sum(bg_gws), length(gws_reg)))

write.table(data.frame(
  analysis = c("gws-removed, full background", "gws-removed, background restricted",
               "merge loci sharing a lead SNP"),
  n_loci = c(nrow(sub), ST, STm), n_known = c(sum(sub$known), SK, SKm),
  bg_loci = c(BT, BT2, BTm), bg_known = c(BK, BK2, BKm),
  fold = c((sum(sub$known) / nrow(sub)) / (BK / BT), (SK / ST) / (BK2 / BT2),
           (SKm / STm) / (BKm / BTm)),
  fisher_p = c(fg(sum(sub$known), nrow(sub)), fis2(SK, ST, BK2, BT2),
               fis2(SKm, STm, BKm, BTm)),
  n_gws_bounded_loci = c(sum(bg_gws), sum(bg_gws), sum(bg_gws)),
  n_gws_distinct_regions = c(NA, NA, length(gws_reg)),
  stringsAsFactors = FALSE),
  "141b_reviewer_conditioning.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nwrote 141b_reviewer_conditioning.tsv\n")
