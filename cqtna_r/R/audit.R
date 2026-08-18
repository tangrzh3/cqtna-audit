## The driver, and the evidence fields.

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

#' Per-gene evidence, as independent fields
#'
#' NOT a single ordered tier. An earlier version collapsed everything into one
#' label, so a compartment ratio above 2 overwrote the locus evidence and a gene
#' on a known locus came back reading "state-informative" with its known-locus
#' status gone. These are different facts about a gene and none of them ranks
#' above another:
#'
#' \itemize{
#'   \item known_locus_status -- known / novel, inherited from the locus
#'   \item outcome_stability -- retained / lost / gained / not tested
#'   \item compartment_ratio, compartment_flag -- expression outside the target
#'     cell type
#'   \item peak_distance_bp -- eQTL-to-GWAS peak distance
#'   \item manual_diagnostics_completed -- always FALSE from this tool
#'   \item overall_interpretation -- a sentence assembled from the fields above,
#'     never a grade
#' }
#' @keywords internal
#' @noRd
cq_evidence <- function(a, stab, comp, peaks) {
  gained <- if (is.null(stab)) character(0) else stab$gained
  genes <- sort(unique(c(a$known_genes, a$novel_genes, gained)))
  if (!length(genes)) return(data.frame())

  status <- rep("not significant here", length(genes))
  names(status) <- genes
  status[genes %in% a$novel_genes] <- "novel"
  status[genes %in% a$known_genes] <- "known"

  stability <- rep("not tested", length(genes))
  names(stability) <- genes
  if (!is.null(stab)) {
    stability[] <- "retained under the second outcome GWAS"
    stability[genes %in% stab$lost] <- "lost under the second outcome GWAS"
    stability[genes %in% gained] <- "gained under the second outcome GWAS"
    never <- stab$coverage$genes_never_tested_in_2
    stability[genes %in% never] <- "not tested under the second outcome GWAS"
  }

  ratio <- rep(NA_real_, length(genes))
  flag <- rep(NA_character_, length(genes))
  if (!is.null(comp) && nrow(comp)) {
    m <- match(genes, comp$gene)
    ratio <- comp$ratio_other_over_target[m]
    flag <- comp$compartment_flag[m]
  }
  dist <- rep(NA_real_, length(genes))
  if (!is.null(peaks) && nrow(peaks)) dist <- peaks$distance_bp[match(genes, peaks$gene)]

  interp <- vapply(seq_along(genes), function(i) {
    bits <- c(
      if (status[[i]] == "known") "on a locus already known for this outcome"
      else if (status[[i]] == "novel") "on a locus not previously reported for this outcome"
      else NULL,
      if (!is.na(flag[i]) && flag[i] != "comparable across compartments") flag[i] else NULL,
      if (grepl("^lost", stability[[i]])) "not retained under the second outcome GWAS" else NULL)
    if (!length(bits)) "" else paste(bits, collapse = "; ")
  }, character(1))

  data.frame(gene = genes, known_locus_status = unname(status),
             outcome_stability = unname(stability),
             compartment_ratio = unname(ratio), compartment_flag = unname(flag),
             peak_distance_bp = unname(dist),
             manual_diagnostics_completed = FALSE,
             overall_interpretation = interp,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' @keywords internal
#' @noRd
cq_status <- function(x, configured) {
  if (!configured) "not run - no input supplied"
  else if (is.null(x)) "configured but empty"
  else if (is.data.frame(x) && !nrow(x)) "configured but empty"
  else if (is.list(x) && !is.data.frame(x) && !length(x)) "configured but empty"
  else "run"
}

#' Run the audit
#'
#' Runs modules A-G and returns one object carrying the results, the per-gene
#' evidence fields, a run status for every module, and the three diagnostics
#' that were not run.
#'
#' @param mr MR results: a data frame, or a `cqtna_mr` from [as_cqtna_mr()].
#' @param known known loci for **this outcome**, or a `cqtna_known`.
#' @param mismatch a different disease's known-locus list, as the negative
#'   control. Omitting it is allowed and warns. Used for module A and, when
#'   supplied, for a mismatched window sweep reported as `G_nc`.
#' @param mr_alt MR results from a second outcome GWAS, for module C.
#' @param instruments,expression,peaks optional inputs for modules D, E and F.
#' @param target_cell_type the cell type the exposure was measured in, for module E.
#' @param build genome build; required unless `mr` is already a `cqtna_mr`.
#' @param locus_kb,known_kb,fdr the three conventions. If `mr` arrives pre-built
#'   with a different `locus_kb` it is re-clustered, so the window named in the
#'   report is always the window actually used.
#' @param locus_method,blocks how variants are partitioned into loci; passed to
#'   [as_cqtna_mr()]. The default `"fixed_centre"` is non-recursive and bounds a
#'   locus at `locus_kb`; `"single_linkage"` reproduces published numbers but can
#'   chain a chromosome arm into one block.
#' @param known_from how a significant locus inherits its known/novel status;
#'   passed to [cqtna_attribution()] and [cqtna_window_sweep()] so every module
#'   uses one convention. The mismatched-list control is run under the same
#'   convention and the audit warns if it is itself significant, which is the
#'   check that exposes a convention inflated by single-linkage chaining.
#' @return an object of class `cqtna_audit`.
#' @export
#' @examples
#' au <- cqtna_audit(cqtna_demo("mr"), cqtna_demo("known"),
#'                   mismatch = cqtna_demo("mismatch"), build = "GRCh38")
#' au$A$fold
cqtna_audit <- function(mr, known, mismatch = NULL, mr_alt = NULL,
                        instruments = NULL, expression = NULL, peaks = NULL,
                        target_cell_type = NULL, build = NULL,
                        locus_kb = 1000, known_kb = 1000, fdr = 0.05,
                        known_from = c("significant_records", "any_record",
                                       "lead_variant"),
                        locus_method = c("fixed_centre", "single_linkage", "blocks"),
                        blocks = NULL) {
  known_from <- match.arg(known_from)
  locus_method <- match.arg(locus_method)
  cq_validate_window(locus_kb, "locus_kb")
  cq_validate_window(known_kb, "known_kb")
  cq_validate_fdr(fdr)

  reclustered <- FALSE
  if (!inherits(mr, "cqtna_mr")) {
    mr <- as_cqtna_mr(mr, locus_kb, build, locus_method, blocks)
  } else if (!identical(as.numeric(attr(mr, "locus_kb")), as.numeric(locus_kb))) {
    ## The window named in the report has to be the window actually used. An
    ## earlier version took an object clustered at 100 kb, never re-clustered,
    ## and printed 1000 kb in the settings line.
    message("mr was built with locus_kb = ", attr(mr, "locus_kb"),
            " but the audit was asked for ", locus_kb,
            "; re-clustering so the reported window is the one used.")
    mr <- as_cqtna_mr(as.data.frame(mr), locus_kb, attr(mr, "build"),
                      attr(mr, "locus_method"), attr(mr, "blocks"))
    reclustered <- TRUE
  }
  if (!inherits(known, "cqtna_known")) known <- as_cqtna_known(known, build)

  res <- list()
  res$A <- cqtna_attribution(mr, known, known_kb, fdr,
                             known_from = known_from)
  if (!is.null(mismatch)) {
    if (!inherits(mismatch, "cqtna_known"))
      mismatch <- suppressWarnings(as_cqtna_known(mismatch, attr(known, "build")))
    res$A_nc <- cqtna_attribution(mr, mismatch, known_kb, fdr,
                                  label = "mismatched (negative control)",
                                  known_from = known_from)
  } else {
    warning("no mismatched-list control supplied; enrichment specific to this ",
            "outcome's genetics is not established.", call. = FALSE)
  }
  res$B <- cqtna_unit_sweep(mr, fdr)
  if (!is.null(mr_alt)) {
    if (!inherits(mr_alt, "cqtna_mr"))
      mr_alt <- as_cqtna_mr(mr_alt, locus_kb, attr(mr, "build"),
                            attr(mr, "locus_method"), attr(mr, "blocks"))
    res$C <- cqtna_stability(mr, mr_alt, fdr, locus_kb)
  }
  if (!is.null(instruments)) res$D <- cqtna_ladder(instruments)
  sig_genes <- sort(unique(mr$gene[mr$fdr < fdr]))
  if (!is.null(expression))
    res$E <- cqtna_compartment(expression, sig_genes, target_cell_type)
  if (!is.null(peaks)) res$F <- cqtna_peak_distance(peaks)
  res$G <- cqtna_window_sweep(mr, known, known_kb, locus_kb, fdr,
                              known_from = known_from)
  # 位点跨度诊断：单连锁串联会让"位点级"计数变成"区块级"计数
  res$spans <- withCallingHandlers(cqtna_locus_spans(mr, fdr),
                                   warning = function(w) invokeRestart("muffleWarning"))
  res$chaining_warning <- any(res$spans$significant_span_kb > 5 * locus_kb)
  if (!is.null(mismatch))
    res$G_nc <- cqtna_window_sweep(mr, mismatch, known_kb, locus_kb, fdr,
                                   known_from = known_from)

  res$evidence <- cq_evidence(res$A, res$C, res$E, res$F)
  res$not_automated <- cqtna_not_automated()
  res$module_status <- c(
    A = cq_status(res$A, TRUE),
    A_nc = cq_status(res$A_nc, !is.null(mismatch)),
    B = cq_status(res$B, TRUE),
    C = cq_status(res$C, !is.null(mr_alt)),
    D = cq_status(res$D, !is.null(instruments)),
    E = cq_status(res$E, !is.null(expression)),
    F = cq_status(res$F, !is.null(peaks)),
    G = cq_status(res$G, TRUE),
    G_nc = cq_status(res$G_nc, !is.null(mismatch)),
    spans = cq_status(res$spans, TRUE))
  # ★ 阴性对照是裁判。错配名单若也显著，该格作废——这是本工具来源研究
  #   预注册里的原话，也是唯一发现 "any_record" 口径在密集资源上失效的检验。
  res$negative_control_failed <- !is.null(res$A_nc) &&
    is.finite(res$A_nc$fisher_p_one_sided) &&
    res$A_nc$fisher_p_one_sided < 0.05 && isTRUE(res$A_nc$fold > 1)
  if (isTRUE(res$negative_control_failed))
    warning("the mismatched-list control is itself significant (",
            sprintf("%.2fx, P = %.3g", res$A_nc$fold, res$A_nc$fisher_p_one_sided),
            "). Under the criterion this tool's source study pre-registered, ",
            "that voids the cell: the enrichment cannot be attributed to this ",
            "outcome's own genetics. Check cqtna_locus_spans() for chaining, ",
            "and compare known_from = \"lead_variant\".", call. = FALSE)
  res$settings <- list(fdr = fdr, locus_kb = locus_kb, known_kb = known_kb,
                       build = attr(mr, "build"), reclustered = reclustered,
                       known_from = known_from,
                       locus_method = attr(mr, "locus_method"))
  structure(res, class = "cqtna_audit")
}
