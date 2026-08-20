## Density-matched permutation control.
##
## The mismatched-list control answers a narrow question: are these loci dense in
## every disease? It cannot answer the wider one: are they dense in a way that is
## itself disease-specific -- more instruments, wider blocks, more genes.
##
## Everything about this test is a choice, and the choices move the answer. On
## the source study's CD4 cell, holding the seed and the number of permutations
## fixed and changing ONLY the matching tolerance:
##
##     tolerance 0.25 -> 2.33x, empirical P = 0.111
##     tolerance 0.50 -> 4.17x, empirical P = 0.024
##     tolerance 1.00 -> 3.48x, empirical P = 0.039
##
## Significance flips on a parameter nobody has argued for. So this function
## refuses to be used casually: the matching specification is recorded in the
## output, a run that cannot match enough loci returns NA instead of a p-value,
## and cqtna_permutation_sensitivity() sweeps the tolerance so the instability is
## visible rather than latent. Freeze the specification before you look.

#' @keywords internal
#' @noRd
CQ_MATCHABLE <- c("n_records", "span", "n_genes", "chr")
## Canonical locus order: chromosome numerically, then start, then end.
##
## Exists because two different things used to depend on the order rows happened
## to arrive in. The loop over significant loci consumed the RNG stream in that
## order, so reversing the input rows re-cut the same random stream into
## different pieces; and sampling without replacement gave whichever locus came
## first the widest pool. Neither is a property of the data. Sorting here, and
## shuffling the draw order inside each replicate (below), removes both.
##
## String sorting would reintroduce the problem through the collation locale, so
## the key is parsed to numbers.
cq_locus_order <- function(loci) {
  chr <- sub(":.*", "", loci)
  rng <- sub(".*:", "", loci)
  start <- suppressWarnings(as.numeric(sub("-.*", "", rng)))
  end <- suppressWarnings(as.numeric(sub(".*-", "", rng)))
  chr_num <- suppressWarnings(as.numeric(chr))
  chr_num[is.na(chr_num)] <- 1000 + rank(chr[is.na(chr_num)],
                                         ties.method = "min")
  order(chr_num, start, end, loci, method = "radix")
}


#' Density-matched permutation control for locus attribution
#'
#' Draws null sets of loci matched to the observed significant set on the
#' properties that could manufacture an enrichment by themselves, and reports
#' where the observed count falls in that null.
#'
#' @section Why significant_records is refused:
#' Under that convention a locus counts as known when one of ITS SIGNIFICANT
#' records is near a known lead SNP. A background locus has no significant
#' records, so a null locus cannot be scored the way the observed ones were: the
#' null would be built on "any record" while the observation is built on
#' "significant records", and the comparison would be between two different
#' quantities. An earlier version of this function did exactly that. The
#' convention remains available in [cqtna_attribution()] for reproducing
#' published numbers, but it has no valid density-matched null. Use
#' `"any_record"` or `"lead_variant"`.
#'
#' @param mr a `cqtna_mr` object.
#' @param known a `cqtna_known` object for this outcome.
#' @param known_kb,fdr as in [cqtna_attribution()].
#' @param known_from `"any_record"` or `"lead_variant"`; see the section above.
#' @param match_on properties null loci must resemble the observed ones on. Any
#'   of `"n_records"`, `"span"`, `"n_genes"`, `"chr"`. Matching on span is what
#'   separates this from the mismatched-list control; matching on `n_genes` is
#'   what makes "controls for gene density" true rather than aspirational.
#' @param tolerance how close a match must be on the continuous properties, on
#'   the log scale. This parameter can flip significance -- see the note at the
#'   top of this file, and [cqtna_permutation_sensitivity()].
#' @param n_perm number of null sets.
#' @param seed passed to [set.seed()] so an empirical p-value is reproducible.
#' @param min_matched_fraction the run fails, returning `NA` for the p-value, if
#'   fewer than this fraction of significant loci have a matched pool. A p-value
#'   computed while some loci could not be matched is not the test specified.
#' @param min_draw_fraction the run fails, returning `NA` for the p-value, if
#'   fewer than this fraction of draws complete without exhausting a matched
#'   pool. Draws that exhaust are discarded, so a p-value computed from the
#'   survivors of heavy attrition is conditioned on surviving, which is not the
#'   test specified either. Default 0.95.
#' @param replace whether a locus may be drawn twice within one null set. The
#'   default `FALSE` discards a draw that exhausts its pool rather than reusing a
#'   locus, because reuse quietly narrows the null.
#' @return a list carrying the full matching specification, the observed count,
#'   the null distribution, and either an empirical p-value or `NA` together with
#'   `failed_because`.
#' @export
#' @examples
#' mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38")
#' kn <- as_cqtna_known(cqtna_demo("known"), build = "GRCh38")
#' r <- cqtna_permutation_control(mr, kn, n_perm = 200, seed = 1)
#' r$empirical_p
cqtna_permutation_control <- function(mr, known, known_kb = 1000, fdr = 0.05,
                                      known_from = c("any_record", "lead_variant"),
                                      match_on = c("n_records", "span", "n_genes"),
                                      tolerance = 0.25,
                                      n_perm = 1000, seed = NULL,
                                      min_matched_fraction = 1,
                                      min_draw_fraction = 0.95,
                                      replace = FALSE) {
  known_from <- match.arg(known_from)
  match_on <- match.arg(match_on, CQ_MATCHABLE, several.ok = TRUE)
  cq_check_build(mr, known)
  cq_validate_fdr(fdr)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0)
    stop("`tolerance` must be a single positive number.", call. = FALSE)
  if (!is.numeric(min_matched_fraction) || min_matched_fraction < 0 ||
      min_matched_fraction > 1)
    stop("`min_matched_fraction` must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(min_draw_fraction) || min_draw_fraction < 0 ||
      min_draw_fraction > 1)
    stop("`min_draw_fraction` must be between 0 and 1.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  idx <- cq_known_index(known)
  lk <- cq_locus_known(mr$locus, mr$chr, mr$pos, mr$fdr, idx, known_kb,
                       known_from, fdr)
  loci <- names(lk$by_locus)
  sig <- unique(as.character(mr$locus[mr$fdr < fdr]))
  if (!length(sig)) stop("no significant loci at this threshold.", call. = FALSE)
  sig <- sig[cq_locus_order(sig)]      # never the order the rows arrived in

  ## Observed and null are scored by the same rule -- lk$by_locus, which under
  ## both permitted conventions is defined for every locus, significant or not.
  prop <- data.frame(
    locus = loci,
    n_records = as.numeric(tapply(mr$pos, mr$locus, length)[loci]),
    span = as.numeric(tapply(mr$pos, mr$locus,
                             function(v) max(v) - min(v))[loci]) + 1,
    n_genes = as.numeric(tapply(mr$gene, mr$locus,
                                function(g) length(unique(g)))[loci]),
    chr = vapply(strsplit(loci, ":"), `[`, character(1), 1),
    known = as.logical(lk$by_locus[loci]),
    stringsAsFactors = FALSE)
  rownames(prop) <- prop$locus
  obs <- sum(prop[sig, "known"])

  pool <- lapply(sig, function(L) {
    tgt <- prop[L, ]
    ok <- !(prop$locus %in% sig)
    if ("chr" %in% match_on) ok <- ok & prop$chr == tgt$chr
    for (v in intersect(match_on, c("n_records", "span", "n_genes")))
      ok <- ok & abs(log(prop[[v]]) - log(tgt[[v]])) <= tolerance
    cand <- prop$locus[ok]
    cand[cq_locus_order(cand)]         # not the collation locale's order
  })
  names(pool) <- sig
  matched <- vapply(pool, length, integer(1)) > 0L
  matched_fraction <- mean(matched)

  ## Recorded with the result because a permutation p-value is only
  ## reproducible together with the stream that generated it.
  spec <- list(known_from = known_from, match_on = match_on,
               tolerance = tolerance, n_perm = n_perm, seed = seed,
               fdr = fdr, known_kb = known_kb, replace = replace,
               min_matched_fraction = min_matched_fraction,
               min_draw_fraction = min_draw_fraction,
               rng_kind = paste(RNGkind(), collapse = "/"),
               r_version = R.version.string,
               collate = Sys.getlocale("LC_COLLATE"))
  base <- list(observed_known = obs, n_significant_loci = length(sig),
               n_matched_loci = sum(matched),
               matched_fraction = matched_fraction)

  if (matched_fraction < min_matched_fraction)
    return(c(spec, base, list(
      unmatched_loci = names(pool)[!matched],
      empirical_p = NA_real_, fold_vs_null = NA_real_,
      null_distribution = numeric(0),
      failed_because = sprintf(
        paste("only %d of %d significant loci had a matched background pool",
              "(%.0f%%, required %.0f%%). Widen the tolerance, drop a matching",
              "variable, or report that this comparison cannot be made -- but",
              "do not read a p-value off a partial match."),
        sum(matched), length(sig), 100 * matched_fraction,
        100 * min_matched_fraction))))

  ## Sampling without replacement inside a draw means the locus handled first
  ## picks from the widest pool. Fixing that order -- by input, or by any
  ## canonical rule -- hands one locus a systematic advantage. Shuffling it per
  ## replicate makes the order part of the randomisation, so no locus is
  ## privileged and the result does not depend on any data-independent ordering.
  ## Everything else is unchanged: still without replacement, still discarding a
  ## draw whose pool is exhausted (S38 sections 1.5 and 1.6).
  n_sig <- length(sig)
  draws <- vapply(seq_len(n_perm), function(i) {
    drawn <- character(0)
    for (L in sig[sample.int(n_sig)]) {
      p <- pool[[L]]
      avail <- if (replace) p else setdiff(p, drawn)
      if (!length(avail)) return(NA_real_)   # pool exhausted: discard this draw
      drawn <- c(drawn, avail[sample.int(length(avail), 1L)])
    }
    sum(prop[drawn, "known"])
  }, numeric(1))
  n_exhausted <- sum(is.na(draws))
  null <- draws[!is.na(draws)]
  pool_sizes <- vapply(pool, length, integer(1))
  effective_draw_fraction <- length(null) / n_perm

  ## A p-value from the surviving draws is only the stated test if most draws
  ## survived. Below the floor, say so rather than quoting one.
  if (length(null) && effective_draw_fraction < min_draw_fraction)
    return(c(spec, base, list(
      unmatched_loci = character(0),
      n_draws_used = length(null), n_draws_exhausted = n_exhausted,
      effective_draw_fraction = effective_draw_fraction,
      min_pool_size = as.integer(min(pool_sizes)),
      median_pool_size = as.numeric(stats::median(pool_sizes)),
      max_pool_size = as.integer(max(pool_sizes)),
      empirical_p = NA_real_, fold_vs_null = NA_real_,
      null_distribution = numeric(0),
      failed_because = sprintf(
        paste("only %d of %d draws (%.1f%%) completed without exhausting a",
              "matched pool; the required minimum is %.0f%%. The pools overlap",
              "too much for sampling without replacement at this tolerance."),
        length(null), n_perm, 100 * effective_draw_fraction,
        100 * min_draw_fraction))))

  ## Total exhaustion used to return without the diagnostics, so the one case
  ## where a reader most needs the pool sizes was the case that withheld them.
  if (!length(null))
    return(c(spec, base, list(
      unmatched_loci = character(0),
      n_draws_used = 0L, n_draws_exhausted = n_exhausted,
      effective_draw_fraction = 0,
      min_pool_size = as.integer(min(pool_sizes)),
      median_pool_size = as.numeric(stats::median(pool_sizes)),
      max_pool_size = as.integer(max(pool_sizes)),
      empirical_p = NA_real_,
      fold_vs_null = NA_real_, null_distribution = numeric(0),
      failed_because = paste("every draw exhausted its matched pool before a",
                             "full null set could be formed; the pools overlap",
                             "too much for sampling without replacement."))))

  c(spec, base, list(
    unmatched_loci = character(0),
    n_draws_used = length(null), n_draws_exhausted = n_exhausted,
    effective_draw_fraction = effective_draw_fraction,
    min_pool_size = as.integer(min(pool_sizes)),
    median_pool_size = as.numeric(stats::median(pool_sizes)),
    max_pool_size = as.integer(max(pool_sizes)),
    null_mean = mean(null), null_sd = stats::sd(null),
    null_quantiles = stats::quantile(null, c(.5, .95, .99)),
    empirical_p = (1 + sum(null >= obs)) / (1 + length(null)),
    fold_vs_null = if (mean(null) > 0) obs / mean(null) else NA_real_,
    failed_because = NA_character_,
    null_distribution = null))
}

#' Sensitivity of a permutation control to its own matching specification
#'
#' Runs [cqtna_permutation_control()] across a grid of tolerances. If the verdict
#' moves across the grid, the test has not identified an estimand and no single
#' p-value from it should be quoted.
#'
#' @inheritParams cqtna_permutation_control
#' @param tolerances the grid to sweep.
#' @return a data frame, one row per tolerance.
#' @export
cqtna_permutation_sensitivity <- function(mr, known, known_kb = 1000, fdr = 0.05,
                                          known_from = c("any_record", "lead_variant"),
                                          match_on = c("n_records", "span", "n_genes"),
                                          tolerances = c(0.1, 0.25, 0.5, 1, 2),
                                          n_perm = 1000, seed = 1) {
  known_from <- match.arg(known_from)
  rows <- lapply(tolerances, function(tol) {
    r <- cqtna_permutation_control(mr, known, known_kb, fdr, known_from,
                                   match_on, tolerance = tol, n_perm = n_perm,
                                   seed = seed)
    data.frame(tolerance = tol, matched_fraction = r$matched_fraction,
               observed_known = r$observed_known,
               fold_vs_null = if (is.null(r$fold_vs_null)) NA_real_ else r$fold_vs_null,
               empirical_p = r$empirical_p,
               ok = is.na(r$failed_because), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  ok <- out[out$ok & is.finite(out$empirical_p), , drop = FALSE]
  if (nrow(ok) > 1 && any(ok$empirical_p < 0.05) && any(ok$empirical_p >= 0.05))
    warning("the verdict changes across the tolerance grid (p from ",
            sprintf("%.3f to %.3f", min(ok$empirical_p), max(ok$empirical_p)),
            "). Fix the specification in advance and report the sweep, not a ",
            "single p-value.", call. = FALSE)
  out
}
