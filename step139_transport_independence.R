## step139 -- S41 的结果有多独立？三项事后诊断。
##
## S41 的六行共用同一批冻结的 CD4 工具（3,556 条记录、2,126 个变异）。
## 于是有三个必须回答的问题，本脚本只回答问题，不改判定：
##
##   1. 六行的显著位点彼此重叠多少？重叠越大，六行越不是六个独立证据。
##   2. 六份已知名单彼此重叠多少？若 lung 的名单覆盖了 melanoma 的位点，
##      lung 的"自身名单富集"可能只是 melanoma 的信号。
##   3. 去掉任意单个显著位点，各行的 Fisher P 会怎么动？
##      若某一行靠一个万能多效区（TERT/CDKN2A/HLA 之类）撑着，这里会暴露。
##
## ⚠ 事后诊断。不参与 S41 §7 判定，不得用于升级任何一行。
## 输出：139a_locus_overlap.tsv、139b_list_overlap.tsv、
##       139c_leave_one_out.tsv、139d_console.log

suppressMessages(library(cqtna))

arg <- commandArgs(trailingOnly = TRUE)
root <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
        if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else getwd()
setwd(root)

FDR <- 0.05; LOCUS_KB <- 1000; KNOWN_KB <- 1000

fix_gene <- function(g) {
  g <- as.character(g); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad)); g
}
std <- function(gene, chr, pos, p) {
  d <- data.frame(record_id = seq_along(pos), gene = fix_gene(gene),
                  chr = as.character(chr), pos = as.numeric(pos),
                  p = as.numeric(p), stringsAsFactors = FALSE)
  d[is.finite(d$p) & is.finite(d$pos) & d$pos > 0, , drop = FALSE]
}

list_files <- c(melanoma = "landi2020_known_loci_grch38.csv",
                HCC = "84a_hcc_known_loci_grch38.csv",
                colorectal = "known_loci_colorectal_grch38.csv",
                prostate = "known_loci_prostate_grch38.csv",
                breast = "known_loci_breast_grch38.csv",
                lung = "known_loci_lung_grch38.csv")
K <- list(); raw <- list()
for (nm in names(list_files)) {
  d <- read.csv(list_files[[nm]], stringsAsFactors = FALSE)
  d <- d[!is.na(suppressWarnings(as.numeric(d$pos))), , drop = FALSE]
  raw[[nm]] <- d; d$source <- nm
  K[[nm]] <- as_cqtna_known(d[, c("chr", "pos", "source")], build = "GRCh38")
}

mel <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
mel$chr <- sub(":.*", "", mel$SNP); mel$pos <- as.numeric(sub(".*:", "", mel$SNP))
cc <- read.delim("36a_crosscancer_MR_all.tsv", stringsAsFactors = FALSE)
cc$chr <- sub(":.*", "", cc$SNP); cc$pos <- as.numeric(sub(".*:", "", cc$SNP))

rows <- list(list(nm = "Melanoma", own = "melanoma",
                  d = std(mel$SYMBOL, mel$chr, mel$pos, mel$pval)))
for (cn in c("Lung", "Colorectal", "Breast", "Prostate")) {
  s <- cc[cc$cancer == cn, , drop = FALSE]
  rows[[length(rows) + 1]] <- list(nm = cn, own = tolower(cn),
                                   d = std(s$symbol, s$chr, s$pos, s$pval))
}

## 六行的位点划分由位置决定，而六行位置相同 -> 位点 id 可直接比较。
sig <- list(); genes <- list(); MR <- list()
for (r in rows) {
  m <- as_cqtna_mr(r$d, locus_kb = LOCUS_KB, build = "GRCh38",
                   locus_method = "fixed_centre")
  MR[[r$nm]] <- m
  s <- unique(as.character(m$locus)[m$fdr < FDR])
  sig[[r$nm]] <- s
  genes[[r$nm]] <- sort(unique(m$gene[m$fdr < FDR]))
}

## ---- 1. 显著位点重叠 --------------------------------------------------------
cat("1. significant-locus overlap between rows (Jaccard, and raw intersection)\n")
nms <- names(sig)
ov <- list()
for (i in seq_along(nms)) for (j in seq_along(nms)) if (i < j) {
  a <- sig[[i]]; b <- sig[[j]]
  inter <- length(intersect(a, b)); uni <- length(union(a, b))
  ov[[length(ov) + 1]] <- data.frame(row_a = nms[i], row_b = nms[j],
    n_a = length(a), n_b = length(b), n_shared = inter,
    jaccard = round(inter / uni, 3), stringsAsFactors = FALSE)
}
OV <- do.call(rbind, ov)
print(OV, row.names = FALSE)
write.table(OV, "139a_locus_overlap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

allsig <- table(unlist(sig))
cat(sprintf("\n  loci significant in >1 row: %d of %d distinct\n",
            sum(allsig > 1), length(allsig)))
if (any(allsig > 1)) {
  sh <- names(allsig)[allsig > 1]
  for (L in sh) {
    inrows <- nms[vapply(sig, function(x) L %in% x, logical(1))]
    g <- unique(unlist(lapply(nms, function(n)
      MR[[n]]$gene[as.character(MR[[n]]$locus) == L & MR[[n]]$fdr < FDR])))
    cat(sprintf("    locus %-10s in %s  (%s)\n", L,
                paste(inrows, collapse = "+"), paste(head(g, 4), collapse = ",")))
  }
}

## ---- 2. 名单之间的重叠 ------------------------------------------------------
cat("\n2. known-list overlap, measured as background loci flagged by both\n")
bg <- as.character(MR[["Melanoma"]]$locus)
bgl <- sort(unique(bg))
flag <- list()
for (nm in names(K)) {
  a <- cqtna_attribution(MR[["Melanoma"]], K[[nm]], known_kb = KNOWN_KB,
                         fdr = FDR, known_from = "any_record")
  flag[[nm]] <- names(a$locus_known)[as.logical(a$locus_known)]
}
lo <- list()
for (i in seq_along(names(K))) for (j in seq_along(names(K))) if (i < j) {
  a <- flag[[i]]; b <- flag[[j]]
  lo[[length(lo) + 1]] <- data.frame(list_a = names(K)[i], list_b = names(K)[j],
    n_a = length(a), n_b = length(b), n_shared = length(intersect(a, b)),
    jaccard = round(length(intersect(a, b)) / length(union(a, b)), 3),
    frac_of_a = round(length(intersect(a, b)) / length(a), 3),
    stringsAsFactors = FALSE)
}
LO <- do.call(rbind, lo)
print(LO, row.names = FALSE)
write.table(LO, "139b_list_overlap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## ---- 3. 逐个去掉显著位点 ----------------------------------------------------
cat("\n3. leave-one-significant-locus-out, own list\n")
loo <- list()
for (r in rows) {
  m <- MR[[r$nm]]
  a0 <- cqtna_attribution(m, K[[r$own]], known_kb = KNOWN_KB, fdr = FDR,
                          known_from = "any_record")
  if (a0$significant_loci < 2) next
  kv <- as.logical(a0$locus_known)
  names(kv) <- names(a0$locus_known)
  s <- unique(as.character(m$locus)[m$fdr < FDR])
  BT <- a0$background_loci; BK <- a0$background_known
  fis <- function(ST, SK) stats::fisher.test(
    matrix(c(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)), 2),
    alternative = "greater")$p.value
  worst_p <- -Inf; worst_L <- NA
  for (L in s) {
    ST <- length(s) - 1L; SK <- sum(kv[setdiff(s, L)])
    p <- fis(ST, SK)
    if (p > worst_p) { worst_p <- p; worst_L <- L }
  }
  g <- unique(m$gene[as.character(m$locus) == worst_L & m$fdr < FDR])
  loo[[length(loo) + 1]] <- data.frame(row = r$nm, sig_loci = a0$significant_loci,
    p_full = a0$fisher_p_one_sided, worst_locus = worst_L,
    worst_genes = paste(head(g, 4), collapse = ","),
    p_worst_drop = worst_p,
    still_sig = worst_p < 0.05, stringsAsFactors = FALSE)
  cat(sprintf("  %-11s %2d loci  P %9.3g -> worst single drop %9.3g  (%s: %s)  %s\n",
              r$nm, a0$significant_loci, a0$fisher_p_one_sided, worst_p,
              worst_L, paste(head(g, 3), collapse = ","),
              if (worst_p < 0.05) "still P<0.05" else "LOSES P<0.05"))
}
LOO <- do.call(rbind, loo)
write.table(LOO, "139c_leave_one_out.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nwrote 139a / 139b / 139c\n")
