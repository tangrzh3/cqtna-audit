## Rendering. The refusal and the three un-automated diagnostics are rendered
## here AND carried on the object, because in R a user will reach for
## as.data.frame() and drop everything that is not a column.
##
## Module outputs carry full precision; rounding happens only here, and only
## through sprintf/formatC — never through round().
##
## ⚠ Why that matters. The locus sweep produces 5.325000000000000177. R's
## round(v, 2) returns 5.32 because it corrects toward the decimal literal;
## sprintf("%.2f", v) and Python's round() both return 5.33, rounding the actual
## binary value, which is above the midpoint. 5.33 is the right answer and is
## what the paper reports, so sprintf is what keeps this implementation and the
## Python reference agreeing to the last displayed digit.

#' @export
print.cqtna_attribution <- function(x, ...) {
  cat(sprintf("%s\n", x$reference))
  cat(sprintf("  background : %d/%d loci (%.2f%%)\n",
              x$background_known, x$background_loci, x$background_pct))
  cat(sprintf("  significant: %d/%d (%.1f%%)\n",
              x$significant_known, x$significant_loci, x$pct_known))
  cat(sprintf("  enrichment %.2fx, one-sided Fisher P = %.3g\n",
              x$fold, x$fisher_p_one_sided))
  invisible(x)
}

#' @export
print.cqtna_audit <- function(x, ...) {
  s <- x$settings
  cat(sprintf("CQTNA audit  (FDR < %s, locus window %s kb, known window %s kb, %s)\n\n",
              s$fdr, s$locus_kb, s$known_kb, s$build))

  a <- x$A
  cat(sprintf("A  attribution : %d/%d significant loci known, background %.2f%%  ->  %.2fx, P = %.3g\n",
              a$significant_known, a$significant_loci, a$background_pct,
              a$fold, a$fisher_p_one_sided))
  if (!is.null(x$A_nc)) {
    nc <- x$A_nc
    clean <- is.na(nc$fold) || nc$fold <= 1 || nc$fisher_p_one_sided > 0.05
    cat(sprintf("   mismatched  : %.2fx, P = %.3g  %s\n", nc$fold,
                nc$fisher_p_one_sided,
                if (clean) "(clean)" else
                  "** POSITIVE - attribution may be a density artefact **"))
  } else {
    cat("   mismatched  : NOT SUPPLIED - outcome-specific attribution not established\n")
  }

  cat(sprintf("B  unit       : %s\n", paste(sprintf("%s %d genes/%d loci",
      x$B$unit, x$B$n_genes, x$B$n_loci), collapse = "; ")))
  if (!is.null(x$C))
    cat(sprintf("C  stability  : %d vs %d genes, %d shared, Jaccard %.3f\n",
                x$C$genes_outcome1, x$C$genes_outcome2, x$C$genes_shared,
                x$C$jaccard_gene))
  if (!is.null(x$D))
    cat(sprintf("D  attrition  : %d pathway genes -> %d instrumentable -> %d analysable -> %d associated\n",
                x$D$pathway_genes, x$D$instrumentable, x$D$analysable, x$D$associated))
  if (!is.null(x$G)) {
    g <- x$G
    ks <- g$known_sweep
    cat(sprintf("G  windows    : known %s kb -> fold %s\n",
                paste(ks$known_window_kb, collapse = "/"),
                paste(formatC(ks$fold, format = "f", digits = 2), collapse = "/")))
    if (!is.null(g$threshold_gap))
      cat(sprintf("   threshold  : %s kb sits between %s and %s kb - %s\n",
                  cq_fmt_kb(g$threshold_gap$threshold_bp),
                  cq_fmt_kb(g$threshold_gap$nearest_below_bp),
                  cq_fmt_kb(g$threshold_gap$nearest_above_bp),
                  if (isTRUE(g$threshold_gap$in_gap))
                    "a gap, so the cut is not near your data"
                  else "close to your data, so the result depends on this choice"))
  }

  cat("\nEvidence tiers\n")
  if (length(x$tiers)) {
    tt <- vapply(x$tiers, function(r) r$tier, character(1))
    for (tier in c("screened", "unresolved", "state-informative")) {
      g <- names(tt)[tt == tier]
      if (length(g)) cat(sprintf("  %-18s %s\n", tier, paste(g, collapse = ", ")))
    }
  } else cat("  no significant genes at this threshold\n")

  cat("\nNo gene can reach `target-supported` from this tool. Not run here:\n")
  for (m in x$not_automated) cat(sprintf("  - %s\n", m$name))
  cat("A nomination audited on five of eight diagnostics is audited on five of eight.\n")
  invisible(x)
}

#' @export
summary.cqtna_audit <- function(object, ...) {
  a <- object$A
  out <- data.frame(
    quantity = c("significant loci", "of those, on known loci",
                 "background loci", "of those, known",
                 "fold enrichment", "one-sided Fisher P"),
    value = c(a$significant_loci, a$significant_known, a$background_loci,
              a$background_known, a$fold, signif(a$fisher_p_one_sided, 3)),
    stringsAsFactors = FALSE)
  if (!is.null(object$A_nc))
    out <- rbind(out, data.frame(
      quantity = c("mismatched-list fold", "mismatched-list P"),
      value = c(object$A_nc$fold, signif(object$A_nc$fisher_p_one_sided, 3)),
      stringsAsFactors = FALSE))
  out
}

#' Write the audit as a markdown report
#'
#' Same sections as the Python reference implementation, so the two can be
#' diffed.
#' @param x a `cqtna_audit`.
#' @param file path to write to; `NULL` returns the text.
#' @return the report text, invisibly when written to a file.
#' @export
cqtna_report <- function(x, file = NULL) {
  s <- x$settings
  L <- c(sprintf("# CQTNA report (R %s)", utils::packageVersion("cqtna")), "",
         sprintf("Settings: FDR < %s, locus window %s kb, known-locus window %s kb, build %s.",
                 s$fdr, s$locus_kb, s$known_kb, s$build), "",
         "## A. Locus attribution", "",
         sprintf("- background: %d/%d loci (%s%%) carry a known lead SNP for this outcome",
                 x$A$background_known, x$A$background_loci, x$A$background_pct),
         sprintf("- significant: %d/%d (%s%%)", x$A$significant_known,
                 x$A$significant_loci, x$A$pct_known),
         sprintf("- **enrichment %sx, one-sided Fisher P = %.3g**",
                 x$A$fold, x$A$fisher_p_one_sided))
  if (!is.null(x$A_nc)) {
    clean <- is.na(x$A_nc$fold) || x$A_nc$fold <= 1 || x$A_nc$fisher_p_one_sided > 0.05
    L <- c(L, sprintf("- mismatched-list control: %sx, P = %.3g %s", x$A_nc$fold,
                      x$A_nc$fisher_p_one_sided,
                      if (clean) "(clean)" else
                        "**(POSITIVE -- attribution may be a density artefact)**"))
  } else {
    L <- c(L, paste("-  no mismatched-list control supplied; enrichment specific",
                    "to this outcome's genetics is not established"))
  }
  L <- c(L, "", "## B. Inference unit", "", cq_md_table(x$B), "",
         paste("A shortest list is not evidence that FDR is controlled under",
               "dependence. State which unit the conclusions are in."), "")
  if (!is.null(x$C))
    L <- c(L, "## C. List stability across outcome GWAS", "",
           sprintf("- genes: %d vs %d, %d shared, Jaccard %s", x$C$genes_outcome1,
                   x$C$genes_outcome2, x$C$genes_shared, x$C$jaccard_gene),
           sprintf("- loci: %d vs %d, Jaccard %s", x$C$loci_outcome1,
                   x$C$loci_outcome2, x$C$jaccard_locus),
           sprintf("- lost: %s", if (length(x$C$lost)) paste(x$C$lost, collapse = ", ") else "none"),
           sprintf("- gained: %s", if (length(x$C$gained)) paste(x$C$gained, collapse = ", ") else "none"), "")
  if (!is.null(x$D))
    L <- c(L, "## D. Instrument attrition", "",
           sprintf("- of %d pathway genes: **%d instrumentable, %d analysable against this outcome, %d nominally associated**",
                   x$D$pathway_genes, x$D$instrumentable, x$D$analysable, x$D$associated), "")
  if (!is.null(x$G)) {
    g <- x$G
    L <- c(L, "## G. Window sensitivity and the distances behind the flag", "",
           cq_md_table(cq_round_p(g$known_sweep)), "",
           cq_md_table(cq_round_p(g$locus_sweep)), "")
    if (length(g$significant_distances_bp))
      L <- c(L, sprintf("- distances of significant loci to the nearest known lead SNP (kb): %s",
                        paste(cq_fmt_kb(g$significant_distances_bp), collapse = ", ")),
             sprintf("- median %s kb, against %s kb for all testable loci",
                     cq_fmt_kb(g$significant_median_bp), cq_fmt_kb(g$background_median_bp)))
    if (!is.null(g$threshold_gap))
      L <- c(L, sprintf("- the %s kb threshold falls between %s kb and %s kb - %s",
                        cq_fmt_kb(g$threshold_gap$threshold_bp),
                        cq_fmt_kb(g$threshold_gap$nearest_below_bp),
                        cq_fmt_kb(g$threshold_gap$nearest_above_bp),
                        if (isTRUE(g$threshold_gap$in_gap))
                          "**a gap, so the cut is not near your data**"
                        else "**close to your data, so the result depends on this choice**"))
    L <- c(L, "", paste("A window chosen to flatter a result weakens when tightened.",
                        "If the fold *rises* as the window narrows, the reported",
                        "value is the conservative one."), "")
  }
  L <- c(L, "## Evidence tiers", "")
  if (length(x$tiers)) {
    L <- c(L, "| gene | tier | basis | stability across outcome GWAS |", "|---|---|---|---|")
    ord <- order(vapply(x$tiers, function(r) r$tier, character(1)), names(x$tiers))
    for (i in ord) {
      r <- x$tiers[[i]]
      L <- c(L, sprintf("| %s | **%s** | %s | %s |", names(x$tiers)[i], r$tier,
                        r$basis, r$stability))
    }
  } else L <- c(L, "no significant genes at this threshold")
  L <- c(L, "", paste("**No gene can reach `target-supported` from this tool.**",
                      "The three diagnostics it cannot run are the ones that would",
                      "license that word."), "",
         "## Not run here -- these require manual work", "")
  for (m in x$not_automated) L <- c(L, sprintf("- **%s**  \n  %s", m$name, m$why))
  L <- c(L, "", paste("A nomination audited on five of eight diagnostics is",
                      "audited on five of eight."), "")
  txt <- paste(L, collapse = "\n")
  if (is.null(file)) return(txt)
  writeLines(txt, file, useBytes = TRUE)
  invisible(txt)
}

#' @keywords internal
#' @noRd
cq_round_p <- function(d) {
  d$fold <- formatC(d$fold, format = "f", digits = 2)
  d$fisher_p_one_sided <- signif(d$fisher_p_one_sided, 3)
  d
}

#' @keywords internal
#' @noRd
cq_md_table <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  hdr <- paste0("| ", paste(names(d), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|")
  body <- apply(d, 1, function(r)
    paste0("| ", paste(ifelse(is.na(r), "", as.character(r)), collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}

#' Plot an audit
#'
#' @param x a `cqtna_audit`.
#' @param which `"attribution"` (fold across the known-locus sweep) or
#'   `"distance"` (the continuous distances behind the binary flag).
#' @param ... passed to the underlying plot call.
#' @return `x`, invisibly.
#' @export
plot.cqtna_audit <- function(x, which = c("attribution", "distance"), ...) {
  which <- match.arg(which)
  if (which == "attribution") {
    g <- x$G$known_sweep
    graphics::plot(g$known_window_kb, g$fold, type = "b", pch = 19, log = "x",
                   xlab = "known-locus window (kb)", ylab = "fold enrichment",
                   main = "Fold against the known-locus window", ...)
    graphics::abline(h = 1, lty = 3, col = "grey50")
    graphics::abline(v = x$settings$known_kb, lty = 2, col = "firebrick")
  } else {
    d <- x$G$significant_distances_bp / 1000
    if (!length(d)) {
      graphics::plot.new(); graphics::title("no significant loci"); return(invisible(x))
    }
    graphics::plot(sort(d), seq_along(d), log = "x", pch = 19, type = "b",
                   xlab = "distance to nearest known lead SNP (kb)",
                   ylab = "significant loci (ranked)",
                   main = "Distances behind the binary flag", ...)
    graphics::abline(v = x$settings$known_kb, lty = 2, col = "firebrick")
  }
  invisible(x)
}
