## step133 -- §7 的第二种读法（确定性投影）。
##
## 为什么有这个脚本：§7 规定了 se_ext 怎么投影、规定了三个 effect-retention
## 因子、规定了跑一万次模拟，但**没有规定外部的效应估计怎么抽**。两种自然读法：
##
##   读法 A（step132）：beta_ext ~ Normal(f * beta_disc, se_ext^2)
##       —— 在已经含噪的 beta_disc 之上再加一层抽样噪声。
##   读法 B（本脚本）：beta_ext = f * beta_disc（确定性），se_ext 按 §7 投影，
##       随机性只来自 harmonisation coverage。
##
## 两种读法给出的 N_eff 门槛不同，**这件事本身是要带回给方案作者的**。
## 本脚本不主张哪一种对，只把 B 也算出来。
##
## ⚠ 顺序说明（照实写）：本脚本是在看过 step132 的部分结果之后写的。
## 是那批结果里 f=0.75、R=1.0 竟然给出 28 个显著位点（发现集实得 8 个），
## 才促使我回去重读 §7 并发现它没规定抽样方式。这不是事前想到的。
##
## 输出：133a_projection_reading_b.tsv、133b_console.log

suppressMessages(library(cqtna))

arg <- commandArgs(trailingOnly = TRUE)
root <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
        if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else {
          a <- commandArgs(FALSE)
          f <- sub("^--file=", "", a[grep("^--file=", a)])
          if (length(f)) dirname(normalizePath(f)) else getwd()
        }
setwd(root)

N_EFF_DISC <- 49337
FDR <- 0.05; LOCUS_KB <- 1000; KNOWN_KB <- 1000
N_REP_COV <- 2000; SEED <- 1

fix_gene <- function(g) {
  g <- as.character(g); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad)); g
}
cd4 <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
cd4$chr <- sub(":.*", "", cd4$SNP); cd4$pos <- as.numeric(sub(".*:", "", cd4$SNP))
d <- data.frame(record_id = seq_len(nrow(cd4)), gene = fix_gene(cd4$SYMBOL),
                chr = cd4$chr, pos = cd4$pos, p = cd4$pval, stringsAsFactors = FALSE)
keep <- is.finite(d$p) & is.finite(d$pos) & d$pos > 0
d <- d[keep, , drop = FALSE]; cd4 <- cd4[keep, , drop = FALSE]
mel <- read.csv("landi2020_known_loci_grch38.csv", stringsAsFactors = FALSE)
mel$source <- "Landi2020"; K <- as_cqtna_known(mel, build = "GRCh38")

vkey <- paste0(cd4$chr, ":", cd4$pos)
VAR <- factor(vkey, levels = unique(vkey[order(
  suppressWarnings(as.numeric(sub("X", "23", cd4$chr))), cd4$pos)]))
vidx <- as.integer(VAR); NV <- nlevels(VAR)
Z_DISC <- cd4$beta_outcome / cd4$se_outcome
ZV <- as.numeric(tapply(Z_DISC, vidx, function(z) z[1]))
stopifnot(all(tapply(Z_DISC, vidx, function(z) diff(range(z))) < 1e-9))

run_once <- function(z_ext, keep_var) {
  ok <- keep_var[vidx]
  dd <- d[ok, , drop = FALSE]
  dd$p <- (2 * stats::pnorm(-abs(z_ext))[vidx])[ok]
  mrx <- as_cqtna_mr(dd, locus_kb = LOCUS_KB, build = "GRCh38",
                     locus_method = "fixed_centre")
  ax <- cqtna_attribution(mrx, K, known_kb = KNOWN_KB, fdr = FDR,
                          known_from = "any_record")
  c(ST = ax$significant_loci, SK = ax$significant_known,
    fold = ax$fold, p = ax$fisher_p_one_sided,
    BT = ax$background_loci, BK = ax$background_known)
}
passes <- function(v)
  v[["ST"]] >= 8 && is.finite(v[["fold"]]) && v[["fold"]] > 1 &&
  is.finite(v[["p"]]) && v[["p"]] < 0.05

## 自校准检查：g = f*sqrt(R) = 1 且全覆盖时，读法 B 必须逐项复现发现集。
{
  v <- run_once(ZV, rep(TRUE, NV))
  cat(sprintf("self-calibration g=1: ST %d SK %d fold %.3f P %.5g (expect 8 4 4.962 0.0048407)\n",
              v[["ST"]], v[["SK"]], v[["fold"]], v[["p"]]))
  stopifnot(v[["ST"]] == 8, v[["SK"]] == 4)
}

R_GRID <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16)
F_GRID <- c(0.50, 0.75, 1.00)

set.seed(SEED)
out <- list()
for (cov in c(1.00, 0.95)) for (f in F_GRID) for (R in R_GRID) {
  g <- f * sqrt(R)
  nrep <- if (cov >= 1) 1 else N_REP_COV
  ST <- integer(nrep); PASS <- logical(nrep)
  for (r in seq_len(nrep)) {
    kv <- if (cov >= 1) rep(TRUE, NV) else runif(NV) < cov
    v <- run_once(g * ZV, kv)
    ST[r] <- v[["ST"]]; PASS[r] <- passes(v)
  }
  out[[length(out) + 1]] <- data.frame(
    reading = "B_deterministic", coverage = cov, n_rep = nrep,
    retention = f, R = R, g = g, n_eff_ext = round(N_EFF_DISC * R),
    cases_if_large_control = round(N_EFF_DISC * R / 4),
    median_sig_loci = median(ST), pct_ge8 = 100 * mean(ST >= 8),
    pct_pass = 100 * mean(PASS), qualifies = 100 * mean(PASS) >= 80)
  cat(sprintf("cov %.2f  f %.2f  R %5.1f  g %.3f  N_eff %8.0f  med_loci %5.1f  pass %5.1f%%\n",
              cov, f, R, g, N_EFF_DISC * R, median(ST), 100 * mean(PASS)))
}
res <- do.call(rbind, out)
write.table(res, "133a_projection_reading_b.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nwrote 133a_projection_reading_b.tsv (", nrow(res), " rows)\n", sep = "")
