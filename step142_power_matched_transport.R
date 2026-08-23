## step142 -- 在匹配功效下比较六个疾病，并换一个不受名单密度支配的统计量。
##
## 为什么需要：S41 的自身 fold 随显著位点数单调下降（8 -> 4.96、10 -> 3.85、
## 23 -> 2.25、28 -> 2.04），所以跨疾病比 fold 是在比功效。
##
## 两个统计量都不合用：
##   fold   = p_sig / p_bg           完全除以背景 —— 背景越密，上限越低（1/p_bg）
##   p_sig  = SK / ST                完全不除背景 —— 我先前误称它为"占上限归一化"，
##                                   它其实恒等于 fold * p_bg，根本没校正密度
## 本脚本采用机会校正量
##   A = (p_sig - p_bg) / (1 - p_bg)
## 读法：**提名把"新位点"那一份额消掉了多大比例**。
## p_sig = p_bg 时 A = 0；全部显著位点都是 known 时 A = 1，与背景密度无关。
##
## 功效匹配：把每个疾病的 |z| 乘以 g（step133 的确定性投影，读法 B），
## 二分搜索使显著 bounded locus 数恰为 k，再在同一个 k 上比较。
##
## ⚠ 事后分析。不改动 S41 §7 的任何判定。
## 输出：142a_power_matched.tsv、142b_console.log

suppressMessages(library(cqtna))

arg <- commandArgs(trailingOnly = TRUE)
root <- if (length(arg) >= 1 && nzchar(arg[1])) arg[1] else
        if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else getwd()
setwd(root)

FDR <- 0.05; LOCUS_KB <- 1000; KNOWN_KB <- 1000
K_TARGETS <- c(8, 10, 15, 20, 25)

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

lf <- c(melanoma = "landi2020_known_loci_grch38.csv",
        colorectal = "known_loci_colorectal_grch38.csv",
        prostate = "known_loci_prostate_grch38.csv",
        breast = "known_loci_breast_grch38.csv",
        lung = "known_loci_lung_grch38.csv")
K <- list()
for (nm in names(lf)) {
  d <- read.csv(lf[[nm]], stringsAsFactors = FALSE)
  d <- d[!is.na(suppressWarnings(as.numeric(d$pos))), , drop = FALSE]
  d$source <- nm
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

## 每行的 |z| 由 p 反解；本投影只用 |z|。
zfp <- function(p) stats::qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)

at_g <- function(r, g) {
  d <- r$d; d$p <- 2 * stats::pnorm(-abs(g * zfp(r$d$p)))
  mr <- as_cqtna_mr(d, locus_kb = LOCUS_KB, build = "GRCh38",
                    locus_method = "fixed_centre")
  a <- cqtna_attribution(mr, K[[r$own]], known_kb = KNOWN_KB, fdr = FDR,
                         known_from = "any_record")
  ST <- a$significant_loci; SK <- a$significant_known
  BT <- a$background_loci;  BK <- a$background_known
  p_sig <- if (ST) SK / ST else NA_real_; p_bg <- BK / BT
  list(ST = ST, SK = SK, BT = BT, BK = BK, p_sig = p_sig, p_bg = p_bg,
       fold = a$fold, A = if (ST) (p_sig - p_bg) / (1 - p_bg) else NA_real_,
       p = a$fisher_p_one_sided)
}

## 二分：位点数关于 g 单调不减，找使 ST == k 的 g
find_g <- function(r, k) {
  lo <- 0.05; hi <- 6
  if (at_g(r, hi)$ST < k) return(NA_real_)
  for (i in 1:44) {
    mid <- (lo + hi) / 2
    if (at_g(r, mid)$ST < k) lo <- mid else hi <- mid
  }
  if (at_g(r, hi)$ST != k) return(NA_real_)   # 位点是离散的，k 可能被跳过
  hi
}

cat("observed, unmatched (S41):\n")
cat(sprintf("%-11s %5s %5s %7s %7s %7s %8s\n",
            "outcome", "loci", "known", "p_sig", "p_bg", "fold", "A"))
out <- list()
for (r in rows) {
  v <- at_g(r, 1)
  out[[length(out) + 1]] <- data.frame(
    outcome = r$nm, matched_k = NA_integer_, g = 1,
    sig_loci = v$ST, sig_known = v$SK, bg_loci = v$BT, bg_known = v$BK,
    p_sig = v$p_sig, p_bg = v$p_bg, fold = v$fold, A_chance_corrected = v$A,
    fisher_p = v$p, stringsAsFactors = FALSE)
  cat(sprintf("%-11s %5d %5d %7.3f %7.3f %7.2f %8.3f\n",
              r$nm, v$ST, v$SK, v$p_sig, v$p_bg, v$fold, v$A))
}

for (k in K_TARGETS) {
  cat(sprintf("\nmatched at %d significant bounded loci:\n", k))
  cat(sprintf("%-11s %6s %5s %7s %7s %7s %8s %10s\n",
              "outcome", "g", "known", "p_sig", "p_bg", "fold", "A", "Fisher P"))
  for (r in rows) {
    g <- find_g(r, k)
    if (!is.finite(g)) {
      cat(sprintf("%-11s %6s  -- cannot be matched to exactly %d loci\n", r$nm, "-", k))
      out[[length(out) + 1]] <- data.frame(
        outcome = r$nm, matched_k = k, g = NA_real_, sig_loci = NA_integer_,
        sig_known = NA_integer_, bg_loci = NA_integer_, bg_known = NA_integer_,
        p_sig = NA_real_, p_bg = NA_real_, fold = NA_real_,
        A_chance_corrected = NA_real_, fisher_p = NA_real_, stringsAsFactors = FALSE)
      next
    }
    v <- at_g(r, g)
    out[[length(out) + 1]] <- data.frame(
      outcome = r$nm, matched_k = k, g = g,
      sig_loci = v$ST, sig_known = v$SK, bg_loci = v$BT, bg_known = v$BK,
      p_sig = v$p_sig, p_bg = v$p_bg, fold = v$fold, A_chance_corrected = v$A,
      fisher_p = v$p, stringsAsFactors = FALSE)
    cat(sprintf("%-11s %6.3f %5d %7.3f %7.3f %7.2f %8.3f %10.3g\n",
                r$nm, g, v$SK, v$p_sig, v$p_bg, v$fold, v$A, v$p))
  }
}

res <- do.call(rbind, out)
write.table(res, "142a_power_matched.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nwrote 142a\n")
