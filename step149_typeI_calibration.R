## step149 -- S46：归属诊断自身的 I 类错误、功效与 void 率。
##
## 规则见 SUPP_typeI_calibration_spec.md（S46），冻结于本脚本之前。
## 本脚本不新增任何判据，也不挑选名单。
##
## 零假设：染色体内环形平移 outcome 的 |z|，保留空间自相关与位点结构，
##         打断"信号位置 ↔ 已知位点位置"的对应关系。
## 备择：  在平移之后，把落在已知位点 1 Mb 内的记录 |z| 乘以 g。
##
## 输出：149a_calibration.tsv、149b_console.log

suppressMessages(library(cqtna))

arg <- commandArgs(trailingOnly = TRUE)
MR <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
      if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else getwd()
setwd(MR)

FDR <- 0.05; LOCUS_KB <- 1000; KNOWN_KB <- 1000
N_REP <- 10000; SEED <- 1
G_ALT <- c(1.25, 1.5, 2.0)

fix_gene <- function(g) {
  g <- as.character(g); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad)); g
}

lf <- c(melanoma = "landi2020_known_loci_grch38.csv",
        lung = "known_loci_lung_grch38.csv",
        breast = "known_loci_breast_grch38.csv",
        colorectal = "known_loci_colorectal_grch38.csv",
        prostate = "known_loci_prostate_grch38.csv",
        HCC = "84a_hcc_known_loci_grch38.csv")
K <- list()
for (nm in names(lf)) {
  d <- read.csv(lf[[nm]], stringsAsFactors = FALSE)
  d <- d[!is.na(suppressWarnings(as.numeric(d$pos))), , drop = FALSE]
  d$source <- nm
  K[[nm]] <- as_cqtna_known(d[, c("chr", "pos", "source")], build = "GRCh38")
}

## 记录族与位点划分：位置不变，故划分恒定（与 step132 同一理由）
mel <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
mel$chr <- sub(":.*", "", mel$SNP); mel$pos <- as.numeric(sub(".*:", "", mel$SNP))
d0 <- data.frame(record_id = seq_len(nrow(mel)), gene = fix_gene(mel$SYMBOL),
                 chr = mel$chr, pos = mel$pos, p = mel$pval,
                 stringsAsFactors = FALSE)
keep <- is.finite(d0$p) & is.finite(d0$pos) & d0$pos > 0
d0 <- d0[keep, , drop = FALSE]; mel <- mel[keep, , drop = FALSE]
Z <- abs(mel$beta_outcome / mel$se_outcome)

mr0 <- as_cqtna_mr(d0, locus_kb = LOCUS_KB, build = "GRCh38",
                   locus_method = "fixed_centre")
LOC <- as.character(mr0$locus)

## 每份名单的恒定量：位点 known 向量、背景计数、以及"该记录是否在已知位点 1 Mb 内"
KV <- list(); BG <- list(); NEAR <- list()
for (nm in names(K)) {
  a <- cqtna_attribution(mr0, K[[nm]], known_kb = KNOWN_KB, fdr = FDR,
                         known_from = "any_record")
  kv <- as.logical(a$locus_known); names(kv) <- names(a$locus_known)
  KV[[nm]] <- kv
  BG[[nm]] <- c(BT = a$background_loci, BK = a$background_known)
  NEAR[[nm]] <- kv[LOC]          # 记录层：其所在位点是否 known
}

fisher_greater <- function(a, b, c, dd)
  suppressWarnings(stats::fisher.test(matrix(c(a, b, c, dd), 2),
                                      alternative = "greater")$p.value)

## 一次重复：给定 |z|，算某份名单的 fold / P / A
score <- function(zz, nm) {
  p <- 2 * stats::pnorm(-zz)
  sig <- stats::p.adjust(p, "BH") < FDR
  if (!any(sig)) return(c(ST = 0, fold = NA_real_, p = NA_real_, A = NA_real_))
  sl <- unique(LOC[sig])
  kv <- KV[[nm]]; BT <- BG[[nm]][["BT"]]; BK <- BG[[nm]][["BK"]]
  ST <- length(sl); SK <- sum(kv[sl])
  p_bg <- BK / BT; p_sig <- SK / ST
  c(ST = ST, fold = p_sig / p_bg,
    p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)),
    A = (p_sig - p_bg) / (1 - p_bg))
}

## 染色体内环形平移
chrs <- split(seq_along(Z), d0$chr)
shift_z <- function() {
  out <- Z
  for (ix in chrs) {
    n <- length(ix)
    if (n > 1) out[ix] <- Z[ix][((seq_len(n) - 1 + sample.int(n, 1)) %% n) + 1]
  }
  out
}

set.seed(SEED)
res <- list()
settings <- c(list(list(lbl = "null (circular shift)", g = 1)),
              lapply(G_ALT, function(g) list(lbl = sprintf("alt g=%.2f", g), g = g)))

for (st in settings) {
  acc <- lapply(names(K), function(nm)
    list(rej = logical(N_REP), fold = numeric(N_REP), A = numeric(N_REP),
         ST = integer(N_REP)))
  names(acc) <- names(K)
  for (r in seq_len(N_REP)) {
    zz <- shift_z()
    for (nm in names(K)) {
      z2 <- zz
      if (st$g != 1) z2[NEAR[[nm]]] <- z2[NEAR[[nm]]] * st$g
      v <- score(z2, nm)
      acc[[nm]]$rej[r] <- is.finite(v[["p"]]) && v[["p"]] < 0.05
      acc[[nm]]$fold[r] <- v[["fold"]]; acc[[nm]]$A[r] <- v[["A"]]
      acc[[nm]]$ST[r] <- v[["ST"]]
    }
  }
  for (nm in names(K)) {
    a <- acc[[nm]]
    rate <- mean(a$rej)
    verdict <- if (st$g != 1) "-" else
      if (rate > 0.07) "OVER-REJECTS" else
      if (rate < 0.03) "conservative" else "calibrated"
    res[[length(res) + 1]] <- data.frame(
      setting = st$lbl, list_name = nm,
      bg_known = BG[[nm]][["BK"]], bg_loci = BG[[nm]][["BT"]],
      bg_share = BG[[nm]][["BK"]] / BG[[nm]][["BT"]],
      n_rep = N_REP, rejection_rate = rate,
      mcse = sqrt(rate * (1 - rate) / N_REP),
      median_sig_loci = stats::median(a$ST),
      median_fold = stats::median(a$fold, na.rm = TRUE),
      median_A = stats::median(a$A, na.rm = TRUE),
      verdict = verdict, stringsAsFactors = FALSE)
    cat(sprintf("%-22s %-11s bg %.3f | loci %4.0f | fold %5.2f | A %6.3f | reject %.4f +- %.4f  %s\n",
                st$lbl, nm, BG[[nm]][["BK"]] / BG[[nm]][["BT"]],
                stats::median(a$ST), stats::median(a$fold, na.rm = TRUE),
                stats::median(a$A, na.rm = TRUE), rate,
                sqrt(rate * (1 - rate) / N_REP), verdict))
  }
}

R <- do.call(rbind, res)
write.table(R, "149a_calibration.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n", strrep("=", 84), "\n", sep = "")
cat("type I error under the null, by reference-list density\n")
n <- R[R$setting == "null (circular shift)", ]
n <- n[order(n$bg_share), ]
for (i in seq_len(nrow(n)))
  cat(sprintf("  %-11s background %5.1f%%  rejection %.4f +- %.4f  %s\n",
              n$list_name[i], 100 * n$bg_share[i], n$rejection_rate[i],
              n$mcse[i], n$verdict[i]))
cat("\nwrote 149a\n")
