## Statistical and interval helpers.
##
## Every function here has an exact counterpart in the Python reference
## implementation (cqtna.py). The oracle test in tests/testthat compares the two
## on the demo data, so any divergence introduced here fails the suite rather
## than propagating silently into a user's report.

#' Benjamini-Hochberg adjusted p-values
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
#' `window_kb`. Returns a STABLE character key per input row, in input order:
#' `chr:start-end`, the extent of the cluster the row belongs to.
#'
#' The key is a character, not an integer counter, for a reason. Module C
#' compares loci between two outcome tables, and per-table integer ids are not
#' comparable -- id 206 was locus 16:87.7 Mb in one table and 16:89.7 Mb in the
#' other, 2 Mb apart, and intersect() called them the same locus. A key carrying
#' chromosome and extent cannot make that mistake.
#'
#' Note that single-linkage extents depend on which records are present, so keys
#' from two different record sets are still not directly comparable. That is why
#' [cqtna_stability()] re-clusters on the union rather than comparing keys built
#' separately.
#' @keywords internal
#' @noRd
cq_assign_loci <- function(chr, pos, window_kb) {
  chr <- as.character(chr)
  pos <- as.numeric(pos)
  n <- length(chr)
  if (n == 0L) return(character(0))
  ord <- order(chr, pos)
  grp <- integer(n)
  id <- 0L
  lc <- NA_character_
  lp <- NA_real_
  for (i in ord) {
    if (is.na(lc) || chr[i] != lc || (pos[i] - lp) > window_kb * 1000) id <- id + 1L
    grp[i] <- id
    lc <- chr[i]
    lp <- pos[i]
  }
  rng <- vapply(split(pos, grp), function(v) paste0(min(v), "-", max(v)), character(1))
  ch <- vapply(split(chr, grp), function(v) v[1], character(1))
  paste0(ch[as.character(grp)], ":", rng[as.character(grp)])
}

## Known-locus lookup -----------------------------------------------------
##
## cq_known_index / cq_nearest_bp / cq_is_known share one index so the distance
## and the binary flag can never disagree about whether a variant is inside a
## given window.

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
    j <- findInterval(p, arr)
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

#' Locus-level known status -- the single source of truth
#'
#' A locus is "known" if ANY record at that locus lies within `window_kb` of a
#' known lead SNP. Status is a property of the LOCUS, computed once over every
#' record, and everything downstream inherits it: the background count, the
#' significant count, the gene labels, module G's sweeps and the evidence
#' fields.
#'
#' This exists because the first version computed the denominator over all
#' records and the numerator over significant records only, and labelled genes
#' from per-record flags. Those three conventions can disagree: a gene could be
#' listed as novel while the locus it sits on was counted as known. On the demo
#' data they happen to agree, which is exactly why it needed a test rather than
#' an inspection.
#'
#' @return a list with `by_locus` (named logical, one entry per locus) and
#'   `by_record` (logical, one entry per record, the locus status broadcast back)
#' @keywords internal
#' @noRd
cq_locus_known <- function(locus, chr, pos, idx, window_kb) {
  rec <- cq_is_known(idx, chr, pos, window_kb)
  by_locus <- tapply(rec, locus, any)
  list(by_locus = by_locus, by_record = unname(by_locus[as.character(locus)]))
}

#' Format base pairs as kb for display
#'
#' Rounds; `format = "d"` alone truncates, which turned 73,600 bp into "73 kb"
#' where the Python reference prints "74".
#' @keywords internal
#' @noRd
cq_fmt_kb <- function(bp) formatC(round(bp / 1000), format = "d", big.mark = ",")
