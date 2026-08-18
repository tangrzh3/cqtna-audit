## Density-matched permutation control.
##
## The mismatched-list control answers a narrow question: are these loci dense in
## every disease? It cannot answer the wider one: are they dense in a way that is
## itself disease-specific -- more instruments, wider blocks, more genes. This
## does, by drawing null locus sets matched on those properties from the study's
## own background.
##
## It is also the control that would have caught the chaining problem without
## needing a second disease's list at all: if wide blocks are what makes a locus
## "known", a null matched on block width says so.

#' Density-matched permutation control for locus attribution
#'
#' Draws `n_perm` null sets of loci from the background, each matched to the
#' observed significant set on the properties that could manufacture an
#' enrichment on their own, and reports where the observed count falls in that
#' null.
#'
#' @param mr a `cqtna_mr` object.
#' @param known a `cqtna_known` object for this outcome.
#' @param known_kb,fdr as in [cqtna_attribution()].
#' @param known_from which convention to test; the null is built under the same
#'   one.
#' @param match_on properties to match null loci on. `"n_records"` is the number
#'   of instruments at the locus, `"span"` its width in bp, `"chr"` its
#'   chromosome. Matching on span is what distinguishes this from the
#'   mismatched-list control.
#' @param n_perm number of null sets.
#' @param seed passed to [set.seed()] so a reported empirical p-value can be
#'   reproduced exactly.
#' @param tolerance how closely a null locus must match on the continuous
#'   properties, as a proportion. Loci are matched within `tolerance` on the log
#'   scale, which is the scale spans actually vary on.
#' @return a list with the observed count, the null distribution, the empirical
#'   one-sided p-value, and how many of the requested matches could be made.
#' @export
#' @examples
#' mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38")
#' kn <- as_cqtna_known(cqtna_demo("known"), build = "GRCh38")
#' p <- cqtna_permutation_control(mr, kn, n_perm = 200, seed = 1)
#' p$empirical_p
cqtna_permutation_control <- function(mr, known, known_kb = 1000, fdr = 0.05,
                                      known_from = c("significant_records",
                                                     "any_record", "lead_variant"),
                                      match_on = c("n_records", "span"),
                                      n_perm = 1000, seed = NULL,
                                      tolerance = 0.25) {
  known_from <- match.arg(known_from)
  match_on <- match.arg(match_on, c("n_records", "span", "chr"), several.ok = TRUE)
  cq_check_build(mr, known)
  cq_validate_fdr(fdr)
  if (!is.null(seed)) set.seed(seed)

  idx <- cq_known_index(known)
  lk <- cq_locus_known(mr$locus, mr$chr, mr$pos, mr$fdr, idx, known_kb,
                       known_from, fdr)
  loci <- names(lk$by_locus)
  sig <- unique(as.character(mr$locus[mr$fdr < fdr]))
  if (!length(sig)) stop("no significant loci at this threshold.", call. = FALSE)

  prop <- data.frame(
    locus = loci,
    n_records = as.numeric(tapply(mr$pos, mr$locus, length)[loci]),
    span = as.numeric(tapply(mr$pos, mr$locus,
                             function(v) max(v) - min(v))[loci]) + 1,
    chr = vapply(strsplit(loci, ":"), `[`, character(1), 1),
    known = as.logical(lk$by_locus[loci]),
    stringsAsFactors = FALSE)
  rownames(prop) <- prop$locus
  obs <- sum(lk$sig_status)

  ## For each significant locus, the pool of background loci resembling it.
  pool <- lapply(sig, function(L) {
    tgt <- prop[L, ]
    ok <- rep(TRUE, nrow(prop))
    if ("chr" %in% match_on) ok <- ok & prop$chr == tgt$chr
    for (v in intersect(match_on, c("n_records", "span")))
      ok <- ok & abs(log(prop[[v]]) - log(tgt[[v]])) <= tolerance
    ok <- ok & !(prop$locus %in% sig)
    p <- prop$locus[ok]
    if (!length(p)) NA_character_ else p
  })
  unmatched <- sum(vapply(pool, function(p) all(is.na(p)), logical(1)))

  null <- vapply(seq_len(n_perm), function(i) {
    drawn <- character(0)
    for (p in pool) {
      if (all(is.na(p))) next
      avail <- setdiff(p, drawn)
      if (!length(avail)) avail <- p
      drawn <- c(drawn, avail[sample.int(length(avail), 1L)])
    }
    sum(prop[drawn, "known"])
  }, numeric(1))

  list(known_from = known_from, matched_on = match_on,
       observed_known = obs, n_significant_loci = length(sig),
       null_mean = mean(null), null_sd = stats::sd(null),
       null_quantiles = stats::quantile(null, c(.5, .95, .99)),
       empirical_p = (1 + sum(null >= obs)) / (1 + n_perm),
       fold_vs_null = if (mean(null) > 0) obs / mean(null) else NA_real_,
       n_perm = n_perm, seed = seed,
       loci_without_a_match = unmatched,
       null_distribution = null)
}
