## The driver, and the evidence tiers.

#' The three diagnostics this tool cannot run
#'
#' Named in every report rather than omitted, because a nomination audited on
#' five of eight diagnostics is audited on five of eight.
#' @export
cqtna_not_automated <- function() {
  list(
    list(name = paste("(ii) colocalisation with explicit multiple-signal",
                      "modelling and an LD reference matched to the outcome cohort"),
         why = "needs regional summary statistics and an LD panel, not a candidate list"),
    list(name = paste("(vii) cell-level matching on lineage composition when",
                      "splitting cells by a score"),
         why = paste("needs the single-cell data and the split itself; a",
                     "cluster-level control can pass while the cell-level one fails")),
    list(name = paste("(viii) code-by-code verification that an endpoint definition",
                      "is unchanged before comparing candidate lists across releases"),
         why = paste("needs the phenotype definitions, which release notes do not",
                     "reliably summarise")))
}

#' Assign an evidence tier to each nominated gene
#'
#' The tier says what *kind* of evidence a gene has. Stability across outcome
#' GWAS is reported beside it rather than folded into it: a gene on a known locus
#' that also fails to replicate is still a gene on a known locus, and collapsing
#' the two produces rows reading "unresolved" over a basis saying "already known".
#'
#' Deliberately conservative. **Nothing reaches `target-supported` from this
#' tool**, because the three diagnostics it cannot run are exactly the ones that
#' would license that word.
#' @keywords internal
#' @noRd
cq_tiers <- function(a, stab, comp) {
  rows <- list()
  for (g in a$known_genes)
    rows[[g]] <- list(tier = "screened",
                      basis = "significant, but on a locus already known for this outcome",
                      stability = "not tested")
  for (g in a$novel_genes)
    rows[[g]] <- list(tier = "unresolved",
                      basis = "significant on a locus not previously reported for this outcome",
                      stability = "not tested")
  if (!is.null(stab)) {
    for (g in names(rows))
      rows[[g]]$stability <- if (g %in% stab$lost) "lost under the second outcome GWAS"
                             else "retained under the second outcome GWAS"
    for (g in stab$gained)
      if (is.null(rows[[g]]))
        rows[[g]] <- list(tier = "unresolved",
                          basis = "significant only under the second outcome GWAS",
                          stability = "gained under the second outcome GWAS")
  }
  if (!is.null(comp) && nrow(comp)) {
    for (i in seq_len(nrow(comp))) {
      g <- comp$gene[i]
      r <- comp$ratio_other_over_target[i]
      if (!is.null(rows[[g]]) && !is.na(r) && r > 2) {
        rows[[g]]$tier <- "state-informative"
        rows[[g]]$basis <- sprintf(
          paste("expression dominated by %s (%sx the target cell type), so",
                "tissue-level data cannot validate a target-cell-specific mechanism"),
          comp$highest_other[i], r)
      }
    }
  }
  rows
}

#' Run the audit
#'
#' Runs modules A-G on a candidate list and returns one object carrying the
#' results, the evidence tiers, and the three diagnostics that were not run.
#'
#' @param mr MR results: a data frame, or a `cqtna_mr` from [as_cqtna_mr()].
#' @param known known loci for **this outcome**, or a `cqtna_known`.
#' @param mismatch a different disease's known-locus list, as the negative
#'   control. Omitting it is allowed and warns: without it you cannot tell
#'   outcome-specific attribution from loci dense in every disease.
#' @param mr_alt MR results from a second outcome GWAS, for module C.
#' @param instruments,expression,peaks optional inputs for modules D, E and F.
#' @param target_cell_type the cell type the exposure was measured in, for module E.
#' @param build genome build; required unless `mr` is already a `cqtna_mr`.
#' @param locus_kb,known_kb,fdr the three conventions. Module G sweeps the first two.
#' @return an object of class `cqtna_audit`.
#' @export
#' @examples
#' au <- cqtna_audit(cqtna_demo("mr"), cqtna_demo("known"),
#'                   mismatch = cqtna_demo("mismatch"), build = "GRCh38")
#' au$A$fold
cqtna_audit <- function(mr, known, mismatch = NULL, mr_alt = NULL,
                        instruments = NULL, expression = NULL, peaks = NULL,
                        target_cell_type = NULL, build = NULL,
                        locus_kb = 1000, known_kb = 1000, fdr = 0.05) {
  if (!inherits(mr, "cqtna_mr")) mr <- as_cqtna_mr(mr, locus_kb, build)
  if (!inherits(known, "cqtna_known")) known <- as_cqtna_known(known, build)

  res <- list()
  res$A <- cqtna_attribution(mr, known, known_kb, fdr)
  if (!is.null(mismatch)) {
    if (!inherits(mismatch, "cqtna_known"))
      mismatch <- suppressWarnings(as_cqtna_known(mismatch, attr(known, "build")))
    res$A_nc <- cqtna_attribution(mr, mismatch, known_kb, fdr,
                                  label = "mismatched (negative control)")
  } else {
    warning("no mismatched-list control supplied; enrichment specific to this ",
            "outcome's genetics is not established.", call. = FALSE)
  }
  res$B <- cqtna_unit_sweep(mr, fdr)
  if (!is.null(mr_alt)) {
    if (!inherits(mr_alt, "cqtna_mr"))
      mr_alt <- as_cqtna_mr(mr_alt, locus_kb, attr(mr, "build"))
    res$C <- cqtna_stability(mr, mr_alt, fdr)
  }
  if (!is.null(instruments)) res$D <- cqtna_ladder(instruments)
  sig_genes <- sort(unique(mr$gene[mr$fdr < fdr]))
  if (!is.null(expression))
    res$E <- cqtna_compartment(expression, sig_genes, target_cell_type)
  if (!is.null(peaks)) res$F <- cqtna_peak_distance(peaks)
  res$G <- cqtna_window_sweep(mr, known, known_kb, locus_kb, fdr)

  res$tiers <- cq_tiers(res$A, res$C, res$E)
  res$not_automated <- cqtna_not_automated()
  res$settings <- list(fdr = fdr, locus_kb = locus_kb, known_kb = known_kb,
                       build = attr(mr, "build"))
  structure(res, class = "cqtna_audit")
}
