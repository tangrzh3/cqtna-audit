## step132 -- S40 §7 的前瞻性功效闸门。
##
## 规则见 EXTERNAL_VALIDATION_PROTOCOL_v9_3.md §7（已由 S40 采纳，冻结）。
## 实现细节见 S40_POWER_GATE_SPEC.md。本脚本不新增任何判据。
##
## 输出：132a_power_gate_curve.tsv   每个 (N_eff, f, coverage) 格的通过率
##       132c_console.log            控制台回执
##
## 用法：Rscript step132_external_power_gate.R [dir]   或设 CQTNA_DIR。

suppressMessages(library(cqtna))

## ---- 路径解析（与其余脚本同一约定） ----------------------------------------
arg <- commandArgs(trailingOnly = TRUE)
root <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
        if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else {
          a <- commandArgs(FALSE)
          f <- sub("^--file=", "", a[grep("^--file=", a)])
          if (length(f)) dirname(normalizePath(f)) else getwd()
        }
setwd(root)

N_EFF_DISC <- 49337          # 12,530 / 789,099，正文 Methods
FDR        <- 0.05
LOCUS_KB   <- 1000
KNOWN_KB   <- 1000
N_REP      <- 10000
SEED       <- 1

## ---- 冻结记录族 -------------------------------------------------------------
fix_gene <- function(g) {
  g <- as.character(g); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad)); g
}
cd4 <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
cd4$chr <- sub(":.*", "", cd4$SNP)
cd4$pos <- as.numeric(sub(".*:", "", cd4$SNP))
d <- data.frame(record_id = seq_len(nrow(cd4)), gene = fix_gene(cd4$SYMBOL),
                chr = cd4$chr, pos = cd4$pos, p = cd4$pval,
                stringsAsFactors = FALSE)
keep <- is.finite(d$p) & is.finite(d$pos) & d$pos > 0
d <- d[keep, , drop = FALSE]; cd4 <- cd4[keep, , drop = FALSE]

mel <- read.csv("landi2020_known_loci_grch38.csv", stringsAsFactors = FALSE)
mel$source <- "Landi2020"
K <- as_cqtna_known(mel, build = "GRCh38")

## ---- 恒定量：位点划分、背景、known 向量 -------------------------------------
mr0 <- as_cqtna_mr(d, locus_kb = LOCUS_KB, build = "GRCh38",
                   locus_method = "fixed_centre")
a0 <- cqtna_attribution(mr0, K, known_kb = KNOWN_KB, fdr = FDR,
                        known_from = "any_record")
LOC   <- as.character(mr0$locus)
LEV   <- sort(unique(LOC))            # by_locus 的名字顺序即 tapply 的排序
KNOWN <- a0$locus_known           # tapply(near, locus, any)，恒定
stopifnot(!is.null(names(KNOWN)), !anyNA(KNOWN))
BT <- length(KNOWN); BK <- sum(KNOWN)
stopifnot(BT == a0$background_loci, BK == a0$background_known)

## z 与变异层映射。位点键解析成数值再排（v10 §四-54）。
Z_DISC  <- cd4$beta_outcome / cd4$se_outcome
vkey    <- paste0(sprintf("%02s", cd4$chr), ":", sprintf("%012.0f", cd4$pos))
vord    <- order(suppressWarnings(as.numeric(sub("X", "23", cd4$chr))), cd4$pos)
VAR     <- factor(vkey, levels = unique(vkey[vord]))
vidx    <- as.integer(VAR)
NV      <- nlevels(VAR)
ZV      <- Z_DISC[!duplicated(vidx)][order(unique(vidx))]
ZV      <- tapply(Z_DISC, vidx, function(z) z[1])   # 同一变异 z 相同
stopifnot(all(tapply(Z_DISC, vidx, function(z) diff(range(z))) < 1e-9))
ZV      <- as.numeric(ZV)

fisher_greater <- function(a, b, c, dd)
  suppressWarnings(stats::fisher.test(matrix(c(a, b, c, dd), 2),
                                      alternative = "greater")$p.value)

## ---- 快速路径（仅在 coverage == 1 时有效） ---------------------------------
## 位点划分与 known 向量在全覆盖下恒定；覆盖率 < 1 会重新锚定 fixed_centre 分区，
## 背景本身随之改变，那种情形必须走包的完整路径（下面的 one_rep_exact）。
one_rep_fast <- function(z_ext) {
  pv <- 2 * stats::pnorm(-abs(z_ext))[vidx]
  sig <- stats::p.adjust(pv, "BH") < FDR
  sl <- unique(LOC[sig])
  ST <- length(sl); SK <- sum(KNOWN[sl])
  if (ST == 0) return(c(ST = 0, SK = 0, fold = NA, p = NA, BT = BT, BK = BK))
  c(ST = ST, SK = SK, fold = (SK / ST) / (BK / BT),
    p = fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK)),
    BT = BT, BK = BK)
}

one_rep_exact <- function(z_ext, keep_var) {
  ok <- keep_var[vidx]
  dd <- d[ok, , drop = FALSE]
  dd$p <- (2 * stats::pnorm(-abs(z_ext))[vidx])[ok]
  mrx <- as_cqtna_mr(dd, locus_kb = LOCUS_KB, build = "GRCh38",
                     locus_method = "fixed_centre")
  ax <- cqtna_attribution(mrx, K, known_kb = KNOWN_KB, fdr = FDR,
                          known_from = "any_record")
  c(ST = ax$significant_loci, SK = ax$significant_known,
    fold = if (is.null(ax$fold)) NA_real_ else ax$fold,
    p = if (is.null(ax$fisher_p_one_sided)) NA_real_ else ax$fisher_p_one_sided,
    BT = ax$background_loci, BK = ax$background_known)
}

passes <- function(v)
  v[["ST"]] >= 8 && is.finite(v[["fold"]]) && v[["fold"]] > 1 &&
  is.finite(v[["p"]]) && v[["p"]] < 0.05

## ---- 等价性测试：快速路径 == 包（全覆盖） ----------------------------------
set.seed(20260821)
for (i in 1:25) {
  zt <- rnorm(NV, 0, 2.2)
  fastr <- one_rep_fast(zt)
  exactr <- one_rep_exact(zt, rep(TRUE, NV))
  if (fastr[["ST"]] != exactr[["ST"]] || fastr[["SK"]] != exactr[["SK"]] ||
      fastr[["BT"]] != exactr[["BT"]] || fastr[["BK"]] != exactr[["BK"]] ||
      !isTRUE(all.equal(unname(fastr[["p"]]), unname(exactr[["p"]]),
                        tolerance = 1e-10)))
    stop("fast path disagrees with cqtna at replicate ", i)
}
cat("equivalence test: OK (25 replicates, full coverage)
")

## 行序不变性（v10 §四-53）：把记录族倒序输入，快速路径的判定必须不变。
{
  zt <- rnorm(NV, 0, 2.2)
  a1 <- one_rep_fast(zt)
  ord <- rev(seq_len(nrow(d)))
  d2 <- d[ord, , drop = FALSE]; vidx2 <- vidx[ord]
  dd <- d2; dd$p <- 2 * stats::pnorm(-abs(zt))[vidx2]
  mrr <- as_cqtna_mr(dd, locus_kb = LOCUS_KB, build = "GRCh38",
                     locus_method = "fixed_centre")
  ar <- cqtna_attribution(mrr, K, known_kb = KNOWN_KB, fdr = FDR,
                          known_from = "any_record")
  if (ar$significant_loci != a1[["ST"]] || ar$significant_known != a1[["SK"]] ||
      ar$background_loci != a1[["BT"]])
    stop("row-order invariance failed")
  cat("row-order invariance: OK
")
}

## ---- 主扫描：coverage = 1.00，快速路径 --------------------------------------
R_GRID <- c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16)
F_GRID <- c(0.50, 0.75, 1.00)

set.seed(SEED)
out <- list()
for (f in F_GRID) for (R in R_GRID) {
  ST <- integer(N_REP); PASS <- logical(N_REP)
  ncp <- f * ZV * sqrt(R)
  for (r in seq_len(N_REP)) {
    v <- one_rep_fast(rnorm(NV, ncp, 1))
    ST[r] <- v[["ST"]]; PASS[r] <- passes(v)
  }
  out[[length(out) + 1]] <- data.frame(
    coverage = 1.00, path = "fast", n_rep = N_REP,
    retention = f, R = R, n_eff_ext = round(N_EFF_DISC * R),
    cases_if_large_control = round(N_EFF_DISC * R / 4),
    median_sig_loci = median(ST), pct_ge8 = 100 * mean(ST >= 8),
    pct_pass = 100 * mean(PASS), qualifies = 100 * mean(PASS) >= 80)
  cat(sprintf("cov 1.00  f %.2f  R %5.1f  N_eff %8.0f  med_loci %4.1f  ge8 %5.1f%%  pass %5.1f%%
",
              f, R, round(N_EFF_DISC * R), median(ST),
              100 * mean(ST >= 8), 100 * mean(PASS)))
}

## ---- 覆盖率检查：coverage = 0.95，包的完整路径，仅 f = 0.50 -----------------
## 95% 是 §2.4 的合格下限。这里重复数较少，因为每次重复都要重算分区。
N_REP_COV <- 2000
set.seed(SEED)
for (R in c(4, 6, 8, 10, 12, 16)) {
  ST <- integer(N_REP_COV); PASS <- logical(N_REP_COV)
  ncp <- 0.50 * ZV * sqrt(R)
  for (r in seq_len(N_REP_COV)) {
    v <- one_rep_exact(rnorm(NV, ncp, 1), runif(NV) < 0.95)
    ST[r] <- v[["ST"]]; PASS[r] <- passes(v)
  }
  out[[length(out) + 1]] <- data.frame(
    coverage = 0.95, path = "exact", n_rep = N_REP_COV,
    retention = 0.50, R = R, n_eff_ext = round(N_EFF_DISC * R),
    cases_if_large_control = round(N_EFF_DISC * R / 4),
    median_sig_loci = median(ST), pct_ge8 = 100 * mean(ST >= 8),
    pct_pass = 100 * mean(PASS), qualifies = 100 * mean(PASS) >= 80)
  cat(sprintf("cov 0.95  f 0.50  R %5.1f  N_eff %8.0f  med_loci %4.1f  ge8 %5.1f%%  pass %5.1f%%
",
              R, round(N_EFF_DISC * R), median(ST),
              100 * mean(ST >= 8), 100 * mean(PASS)))
}

res <- do.call(rbind, out)
write.table(res, "132a_power_gate_curve.tsv", sep = "	",
            row.names = FALSE, quote = FALSE)
cat("
wrote 132a_power_gate_curve.tsv (", nrow(res), " rows)
", sep = "")
