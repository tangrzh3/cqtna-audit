## Modules A-G. Each is usable on its own; cqtna_audit() runs them together.
##
## Everything is computed by independent locus wherever a count carries an
## argument, because a significant locus does not name a gene.

#' A. Locus attribution against the outcome's own known loci
#'
#' Fold enrichment of FDR-significant loci on previously reported loci for this
#' outcome, counted by independent locus, with a one-sided Fisher exact test.
#'
#' Run it a second time with a different disease's list as `known` to get the
#' mismatched negative control. Without that control you cannot distinguish
#' attribution specific to your outcome's genetics from loci that are dense in
#' every disease.
#'
#' @param mr a `cqtna_mr` object from [as_cqtna_mr()].
#' @param known a `cqtna_known` object from [as_cqtna_known()].
#' @param known_kb window, in kb, within which a variant counts as landing on a
#'   known locus.
#' @param fdr FDR threshold.
#' @param label what this reference list is, used in the report.
#' @return a list with the counts, `fold`, `fisher_p_one_sided`, and the
#'   significant genes split by whether their locus was already known.
#' @export
#' @examples
#' mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38")
#' kn <- as_cqtna_known(cqtna_demo("known"), build = "GRCh38")
#' cqtna_attribution(mr, kn)$fold
cqtna_attribution <- function(mr, known, known_kb = 1000, fdr = 0.05,
                              label = "known-locus list") {
  cq_check_build(mr, known)
  idx <- cq_known_index(known)
  kn <- cq_is_known(idx, mr$chr, mr$pos, known_kb)
  bg <- tapply(kn, mr$locus, any)
  sel <- mr$fdr < fdr
  sg <- tapply(kn[sel], mr$locus[sel], any)
  BT <- length(bg); BK <- sum(bg)
  ST <- length(sg); SK <- sum(sg)
  fold <- if (ST > 0 && BK > 0) (SK / ST) / (BK / BT) else NA_real_
  p <- if (ST > 0) cq_fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK))
       else NA_real_
  structure(list(
    reference = label,
    background_known = as.integer(BK), background_loci = as.integer(BT),
    background_pct = if (BT) 100 * BK / BT else NA_real_,
    significant_loci = as.integer(ST), significant_known = as.integer(SK),
    pct_known = if (ST) 100 * SK / ST else NA_real_,
    fold = fold,
    fisher_p_one_sided = p,
    known_genes = sort(unique(mr$gene[sel & kn])),
    novel_genes = sort(unique(mr$gene[sel & !kn]))),
    class = "cqtna_attribution")
}

#' B. The same list under four inference units
#'
#' The FDR-significant list recomputed over record, variant, gene and
#' independent locus, collapsing with Simes. A shortest list is not evidence
#' that FDR is controlled under dependence; the point is to state which unit the
#' conclusions are in.
#'
#' @inheritParams cqtna_attribution
#' @return a data frame, one row per unit.
#' @export
cqtna_unit_sweep <- function(mr, fdr = 0.05) {
  out <- data.frame(
    unit = "record", n_tests = nrow(mr),
    n_significant = sum(mr$fdr < fdr),
    n_genes = length(unique(mr$gene[mr$fdr < fdr])),
    n_loci = length(unique(mr$locus[mr$fdr < fdr])),
    stringsAsFactors = FALSE)
  keys <- list(variant = c("chr", "pos"), gene = "gene", locus = "locus")
  for (unit in names(keys)) {
    k <- keys[[unit]]
    grp <- interaction(mr[k], drop = TRUE, sep = "\r")
    pv <- tapply(mr$p, grp, cq_simes)
    q <- cq_bh(as.numeric(pv))
    hit <- names(pv)[q < fdr]
    sub <- mr[grp %in% hit, , drop = FALSE]
    out <- rbind(out, data.frame(
      unit = unit, n_tests = length(pv), n_significant = length(hit),
      n_genes = length(unique(sub$gene)), n_loci = length(unique(sub$locus)),
      stringsAsFactors = FALSE))
  }
  out
}

#' C. List stability against a second outcome GWAS
#'
#' @param mr,mr2 two `cqtna_mr` objects from the same exposure data and
#'   different outcome GWAS.
#' @param fdr FDR threshold.
#' @return a list with Jaccard indices and what moved.
#' @export
cqtna_stability <- function(mr, mr2, fdr = 0.05) {
  a <- unique(mr$gene[mr$fdr < fdr])
  b <- unique(mr2$gene[mr2$fdr < fdr])
  la <- unique(mr$locus[mr$fdr < fdr])
  lb <- unique(mr2$locus[mr2$fdr < fdr])
  jac <- function(x, y) if (length(union(x, y)))
    length(intersect(x, y)) / length(union(x, y)) else NA_real_
  list(genes_outcome1 = length(a), genes_outcome2 = length(b),
       genes_shared = length(intersect(a, b)),
       jaccard_gene = jac(a, b),
       loci_outcome1 = length(la), loci_outcome2 = length(lb),
       jaccard_locus = jac(la, lb),
       lost = sort(setdiff(a, b)), gained = sort(setdiff(b, a)))
}

#' D. Instrument attrition at three levels
#'
#' Instrumentable, analysable against this outcome, nominally associated. These
#' are three different statements and collapsing them attributes an outcome-side
#' limitation to the biology of the pathway.
#'
#' @param instruments data frame with `gene`, `instrumentable`, `analysable`,
#'   `associated` (0/1 or logical).
#' @return a list of the three counts.
#' @export
cqtna_ladder <- function(instruments) {
  d <- as.data.frame(instruments, stringsAsFactors = FALSE)
  need <- c("instrumentable", "analysable", "associated")
  missing <- setdiff(need, names(d))
  if (length(missing))
    stop("instruments is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  list(pathway_genes = nrow(d),
       instrumentable = sum(as.numeric(d$instrumentable) > 0),
       analysable = sum(as.numeric(d$analysable) > 0),
       associated = sum(as.numeric(d$associated) > 0),
       note = paste("report all three; collapsing them attributes an",
                    "outcome-side limit to the biology"))
}

#' E. Compartment attribution
#'
#' For each nominated gene, the ratio of its highest expression in any other
#' cell type to its expression in the cell type the instrument came from. A
#' large ratio means tissue-level data cannot validate a target-cell-specific
#' mechanism.
#'
#' @param expression data frame with `gene`, `cell_type`, `mean_expression`.
#' @param genes character vector of nominated genes.
#' @param target_cell_type the cell type the exposure was measured in.
#' @return a data frame, one row per gene that could be evaluated.
#' @export
cqtna_compartment <- function(expression, genes, target_cell_type) {
  x <- as.data.frame(expression, stringsAsFactors = FALSE)
  need <- c("gene", "cell_type", "mean_expression")
  missing <- setdiff(need, names(x))
  if (length(missing))
    stop("expression is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  rows <- list()
  for (g in genes) {
    s <- x[x$gene == g, , drop = FALSE]
    if (!nrow(s) || is.null(target_cell_type) ||
        !target_cell_type %in% s$cell_type) next
    v <- stats::setNames(as.numeric(s$mean_expression), s$cell_type)
    ref <- v[[target_cell_type]]
    other <- v[names(v) != target_cell_type]
    if (!length(other)) next
    top <- names(other)[which.max(other)]
    rows[[length(rows) + 1L]] <- data.frame(
      gene = g, target_cell_type = target_cell_type,
      target_expression = ref, highest_other = top,
      highest_other_expression = other[[top]],
      ratio_other_over_target = if (ref) other[[top]] / ref else NA_real_,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

#' F. Distance between eQTL and GWAS peaks
#'
#' @param peaks data frame with `gene`, `eqtl_pos`, `gwas_pos`.
#' @return the same rows with `distance_bp`.
#' @export
cqtna_peak_distance <- function(peaks) {
  x <- as.data.frame(peaks, stringsAsFactors = FALSE)
  need <- c("gene", "eqtl_pos", "gwas_pos")
  missing <- setdiff(need, names(x))
  if (length(missing))
    stop("peaks is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  x$distance_bp <- abs(as.numeric(x$eqtl_pos) - as.numeric(x$gwas_pos))
  x[, c("gene", "eqtl_pos", "gwas_pos", "distance_bp")]
}

CQ_SWEEP_KB <- c(100, 250, 500, 1000)

#' G. Window sensitivity, and the distances behind the binary flag
#'
#' Two different conventions get called "the 1 Mb window" and they carry
#' different weight, so they are swept separately: `known_kb` decides what counts
#' as landing on a known locus and moves the fold directly; `locus_kb` decides
#' how independent loci are defined and moves the denominator.
#'
#' A threshold chosen to flatter a result weakens when tightened. If the fold
#' *rises* as the window narrows, the value you reported is the conservative one.
#'
#' The distance for a significant locus is taken over that locus's significant
#' records only, because that is what the binary flag is computed from. Taking it
#' over every record in the locus lets a non-significant variant supply the
#' distance and yields a number that contradicts the flag.
#'
#' @inheritParams cqtna_attribution
#' @param locus_kb the locus-definition window held fixed while `known_kb` is
#'   swept.
#' @param sweep_kb windows to sweep, in kb.
#' @return a list with `known_sweep`, `locus_sweep`, the significant-locus
#'   distances, and whether the chosen threshold falls in a gap.
#' @export
cqtna_window_sweep <- function(mr, known, known_kb = 1000, locus_kb = 1000,
                               fdr = 0.05, sweep_kb = CQ_SWEEP_KB) {
  cq_check_build(mr, known)
  idx <- cq_known_index(known)
  dist_bp <- cq_nearest_bp(idx, mr$chr, mr$pos)
  sel <- mr$fdr < fdr

  enrich <- function(loci, kb) {
    kn <- dist_bp <= kb * 1000
    bg <- tapply(kn, loci, any)
    sg <- tapply(kn[sel], loci[sel], any)
    BT <- length(bg); BK <- sum(bg); ST <- length(sg); SK <- sum(sg)
    fold <- if (ST > 0 && BK > 0) (SK / ST) / (BK / BT) else NA_real_
    data.frame(background_loci = BT, background_known = BK,
               significant_loci = ST, significant_known = SK,
               fold = fold,
               fisher_p_one_sided = if (ST > 0)
                 cq_fisher_greater(SK, ST - SK, BK - SK, (BT - BK) - (ST - SK))
                 else NA_real_)
  }

  known_sweep <- do.call(rbind, lapply(sweep_kb, function(kb)
    cbind(known_window_kb = kb, locus_window_kb = locus_kb,
          enrich(mr$locus, kb))))
  locus_sweep <- do.call(rbind, lapply(sweep_kb, function(kb)
    cbind(locus_window_kb = kb, known_window_kb = known_kb,
          enrich(cq_assign_loci(mr$chr, mr$pos, kb), known_kb))))

  sig_d <- sort(tapply(dist_bp[sel], mr$locus[sel], min))
  sig_d <- sig_d[is.finite(sig_d)]
  bg_d <- tapply(dist_bp, mr$locus, min)
  bg_d <- bg_d[is.finite(bg_d)]

  gap <- NULL
  if (length(sig_d) > 1) {
    below <- sig_d[sig_d <= known_kb * 1000]
    above <- sig_d[sig_d > known_kb * 1000]
    if (length(below) && length(above))
      gap <- list(nearest_below_bp = as.numeric(below[length(below)]),
                  nearest_above_bp = as.numeric(above[1]),
                  threshold_bp = known_kb * 1000,
                  in_gap = as.numeric(above[1]) >=
                    4 * max(as.numeric(below[length(below)]), 1))
  }
  list(known_sweep = known_sweep, locus_sweep = locus_sweep,
       significant_distances_bp = as.numeric(sig_d),
       significant_median_bp = if (length(sig_d)) stats::median(sig_d) else NA_real_,
       background_median_bp = if (length(bg_d)) stats::median(bg_d) else NA_real_,
       threshold_gap = gap)
}
