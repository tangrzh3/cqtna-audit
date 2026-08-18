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

#' Partition variants into loci
#'
#' Three rules, because the choice is load-bearing and should be visible in the
#' call rather than buried in a helper.
#'
#' \describe{
#'   \item{`"fixed_centre"`}{(default) Non-recursive. Take the first unassigned
#'     variant on a chromosome as a centre, claim every variant within
#'     `locus_kb` of it, repeat. A locus therefore spans at most `locus_kb` and
#'     cannot grow by chaining. Centres are chosen by position, not by
#'     significance, so the partition does not depend on the outcome statistics
#'     -- which matters because the same partition supplies the denominator.}
#'   \item{`"single_linkage"`}{Two variants join when they are within
#'     `locus_kb`, transitively. This is what the source study published, and it
#'     is kept for reproducing those numbers. On a dense resource it chains: in
#'     that study whole-blood loci reached 30.8 Mb at a 1 Mb window, so
#'     "locus-level" counts there were counting blocks.}
#'   \item{`"blocks"`}{Assign by membership of supplied intervals, e.g.
#'     LD blocks. Variants outside every block each form their own locus.}
#' }
#'
#' All three return a stable `chr:start-end` key per row, in input order. The key
#' carries the chromosome, so loci from two tables can never match across
#' chromosomes -- an integer counter allowed exactly that.
#'
#' @param chr,pos variant coordinates.
#' @param window_kb the window. For `"fixed_centre"` it is the radius claimed
#'   from each centre, so the maximum span is `window_kb`; for
#'   `"single_linkage"` it is the joining distance and the span is unbounded.
#' @param method one of the three above.
#' @param blocks for `method = "blocks"`, a data frame with `chr`, `start`, `end`.
#' @keywords internal
#' @noRd
cq_assign_loci <- function(chr, pos, window_kb,
                           method = c("fixed_centre", "single_linkage", "blocks"),
                           blocks = NULL) {
  method <- match.arg(method)
  chr <- as.character(chr)
  pos <- as.numeric(pos)
  n <- length(chr)
  if (n == 0L) return(character(0))
  w <- window_kb * 1000
  grp <- integer(n)

  if (method == "blocks") {
    if (is.null(blocks)) stop("method = \"blocks\" needs a `blocks` table.",
                              call. = FALSE)
    b <- as.data.frame(blocks, stringsAsFactors = FALSE)
    for (cc in c("chr", "start", "end"))
      if (!cc %in% names(b))
        stop("`blocks` is missing column `", cc, "`.", call. = FALSE)
    b$chr <- as.character(b$chr)
    id <- 0L
    key <- character(n)
    hit <- rep(NA_integer_, n)
    for (i in seq_len(nrow(b))) {
      sel <- chr == b$chr[i] & pos >= b$start[i] & pos <= b$end[i] & is.na(hit)
      if (any(sel)) hit[sel] <- i
    }
    for (i in seq_len(n)) {
      key[i] <- if (is.na(hit[i])) paste0(chr[i], ":", pos[i], "-", pos[i])
                else paste0(b$chr[hit[i]], ":", b$start[hit[i]], "-", b$end[hit[i]])
    }
    return(key)
  }

  ord <- order(chr, pos)
  if (method == "single_linkage") {
    id <- 0L; lc <- NA_character_; lp <- NA_real_
    for (i in ord) {
      if (is.na(lc) || chr[i] != lc || (pos[i] - lp) > w) id <- id + 1L
      grp[i] <- id
      lc <- chr[i]; lp <- pos[i]
    }
  } else {
    ## Non-recursive: a centre claims [centre, centre + w] and the next
    ## unclaimed variant starts a new locus. Nothing joins transitively, so the
    ## span is bounded by w whatever the density.
    id <- 0L; lc <- NA_character_; centre <- NA_real_
    for (i in ord) {
      if (is.na(lc) || chr[i] != lc || (pos[i] - centre) > w) {
        id <- id + 1L
        centre <- pos[i]
        lc <- chr[i]
      }
      grp[i] <- id
    }
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
cq_locus_known <- function(locus, chr, pos, fdr, idx, window_kb,
                           known_from = "significant_records", fdr_threshold = 0.05) {
  near <- cq_is_known(idx, chr, pos, window_kb)
  sel <- fdr < fdr_threshold
  if (known_from == "lead_variant") {
    # 每个位点用它最显著的那条记录代表 —— 分子分母同一规则，且不受串联影响
    keep <- !duplicated(locus[order(fdr)])
    ord <- order(fdr)
    rep_near <- stats::setNames(near[ord][keep], as.character(locus[ord][keep]))
    by_locus <- rep_near[sort(names(rep_near))]
    sig_status <- by_locus[unique(as.character(locus[sel]))]
  } else {
    by_locus <- tapply(near, locus, any)          # 背景一律按全部记录
    sig_status <- if (known_from == "any_record")
      by_locus[unique(as.character(locus[sel]))]
    else tapply(near[sel], locus[sel], any)       # 已发表口径：分子只看显著记录
  }
  list(by_locus = by_locus, sig_status = sig_status,
       by_record = unname(by_locus[as.character(locus)]),
       sig_by_record = unname(sig_status[as.character(locus)]),
       known_from = known_from)
}

#' Valid values for `known_from`
#' @keywords internal
#' @noRd
CQ_KNOWN_FROM <- c("significant_records", "any_record", "lead_variant")

#' Format base pairs as kb for display
#'
#' Rounds; `format = "d"` alone truncates, which turned 73,600 bp into "73 kb"
#' where the Python reference prints "74".
#' @keywords internal
#' @noRd
cq_fmt_kb <- function(bp) formatC(round(bp / 1000), format = "d", big.mark = ",")
