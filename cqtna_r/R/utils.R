## Statistical and interval helpers.
##
## Every function here has an exact counterpart in the Python reference
## implementation (cqtna.py). The oracle test in tests/testthat compares the two
## on the demo data, so any divergence introduced here fails the suite rather
## than propagating silently into a user's report.

#' Benjamini-Hochberg adjusted p-values
#' @param p numeric vector of p-values
#' @return numeric vector of adjusted p-values
#' @keywords internal
#' @noRd
cq_bh <- function(p) stats::p.adjust(p, method = "BH")

#' Simes combination of dependent p-values
#'
#' Used to collapse several records onto one variant, gene or locus. Simes is
#' valid under positive dependence, which is the situation here: records at one
#' locus share an outcome GWAS.
#' @keywords internal
#' @noRd
cq_simes <- function(p) {
  p <- sort(as.numeric(p))
  n <- length(p)
  if (n == 0L) return(NA_real_)
  min(p * n / seq_len(n))
}

#' One-sided Fisher exact test, upper tail
#'
#' P(X >= a) for X hypergeometric. Written against phyper rather than
#' fisher.test so it returns the same tail the Python implementation sums
#' directly, to the last digit.
#' @keywords internal
#' @noRd
cq_fisher_greater <- function(a, b, c, d) {
  if (a + b == 0L) return(NA_real_)
  stats::phyper(a - 1, a + c, b + d, a + b, lower.tail = FALSE)
}

#' Single-linkage clustering of variant positions into independent loci
#'
#' Two variants on one chromosome join the same locus when they are within
#' `window_kb`. Returns an integer locus id per input row, in input order.
#' @keywords internal
#' @noRd
cq_assign_loci <- function(chr, pos, window_kb) {
  chr <- as.character(chr)
  pos <- as.numeric(pos)
  n <- length(chr)
  if (n == 0L) return(integer(0))
  ord <- order(chr, pos)
  out <- integer(n)
  id <- 0L
  lc <- NA_character_
  lp <- NA_real_
  for (i in ord) {
    if (is.na(lc) || chr[i] != lc || (pos[i] - lp) > window_kb * 1000) id <- id + 1L
    out[i] <- id
    lc <- chr[i]
    lp <- pos[i]
  }
  out
}

## Known-locus lookup -----------------------------------------------------
##
## cq_known_index / cq_nearest_bp / cq_is_known share one index so the distance
## and the binary flag can never disagree about whether a variant is inside a
## given window. Keeping them separate is how the Python version first produced
## a distance that contradicted its own flag.

#' @keywords internal
#' @noRd
cq_known_index <- function(known) {
  k <- data.frame(chr = as.character(known$chr), pos = as.numeric(known$pos))
  k <- k[stats::complete.cases(k), , drop = FALSE]
  lapply(split(k$pos, k$chr), sort)
}

#' @keywords internal
#' @noRd
cq_nearest_bp <- function(idx, chr, pos) {
  chr <- as.character(chr)
  pos <- as.numeric(pos)
  out <- rep(Inf, length(pos))
  for (ch in unique(chr)) {
    arr <- idx[[ch]]
    sel <- which(chr == ch)
    if (is.null(arr) || length(arr) == 0L) next
    p <- pos[sel]
    j <- findInterval(p, arr)          # number of arr entries <= p
    lo <- ifelse(j >= 1L, abs(arr[pmax(j, 1L)] - p), Inf)
    hi <- ifelse(j < length(arr), abs(arr[pmin(j + 1L, length(arr))] - p), Inf)
    out[sel] <- pmin(lo, hi)
  }
  out
}

#' @keywords internal
#' @noRd
cq_is_known <- function(idx, chr, pos, window_kb) {
  cq_nearest_bp(idx, chr, pos) <= window_kb * 1000
}

#' Format base pairs as kb for display
#'
#' Rounds; `format = "d"` alone truncates, which turned 73,600 bp into "73 kb"
#' where the Python reference prints "74".
#' @keywords internal
#' @noRd
cq_fmt_kb <- function(bp) formatC(round(bp / 1000), format = "d", big.mark = ",")
