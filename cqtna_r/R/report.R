## Rendering. The refusal, the three un-automated diagnostics and the per-module
## run status are rendered here AND carried on the object, because in R a user
## will reach for as.data.frame() and drop everything that is not a column.
##
## Module outputs carry full precision; rounding happens only here, and only
## through sprintf/formatC -- never through round().
##
## Why that matters. The locus sweep produces 5.325000000000000177. R's
## round(v, 2) returns 5.32 because it corrects toward the decimal literal;
## sprintf("%.2f", v) and Python's round() both return 5.33, rounding the actual
## binary value, which is above the midpoint. 5.33 is the right answer and is
## what the paper reports, so sprintf is what keeps this implementation and the
## Python reference agreeing to the last displayed digit.

#' @keywords internal
#' @noRd
cq_num <- function(v, digits = 2) {
  ifelse(is.na(v), NA_character_, sprintf(paste0("%.", digits, "f"), v))
}

#' @export
print.cqtna_attribution <- function(x, ...) {
  cat(sprintf("%s\n", x$reference))
  cat(sprintf("  background : %d/%d loci (%s%%)\n",
              x$background_known, x$background_loci, cq_num(x$background_pct)))
  cat(sprintf("  significant: %d/%d (%s%%)\n",
              x$significant_known, x$significant_loci, cq_num(x$pct_known, 1)))
  cat(sprintf("  enrichment %sx, one-sided Fisher P = %.3g\n",
              cq_num(x$fold), x$fisher_p_one_sided))
  invisible(x)
}

#' @export
print.cqtna_audit <- function(x, ...) {
  s <- x$settings
  cat(sprintf("CQTNA audit  (FDR < %s, locus window %s kb, known window %s kb, %s)\n",
              s$fdr, s$locus_kb, s$known_kb, s$build))
  if (isTRUE(s$reclustered))
    cat("  note: the supplied object was re-clustered to the requested locus window\n")
  cat("\n")

  a <- x$A
  cat(sprintf("A  attribution : %d/%d significant loci known, background %s%%  ->  %sx, P = %.3g\n",
              a$significant_known, a$significant_loci, cq_num(a$background_pct),
              cq_num(a$fold), a$fisher_p_one_sided))
  if (!is.null(x$A_nc)) {
    nc <- x$A_nc
    clean <- is.na(nc$fold) || nc$fold <= 1 || nc$fisher_p_one_sided > 0.05
    cat(sprintf("   mismatched  : %sx, P = %.3g  %s\n", cq_num(nc$fold),
                nc$fisher_p_one_sided,
                if (clean) "(clean)" else
                  "** POSITIVE - attribution may be a density artefact **"))
  }
  cat(sprintf("B  unit       : %s\n", paste(sprintf("%s %d genes/%d loci",
      x$B$unit, x$B$n_genes, x$B$n_loci), collapse = "; ")))
  if (!is.null(x$C))
    cat(sprintf("C  stability  : %d vs %d genes, %d shared, Jaccard %s; loci %d vs %d, %d shared, Jaccard %s\n",
                x$C$genes_outcome1, x$C$genes_outcome2, x$C$genes_shared,
                cq_num(x$C$jaccard_gene, 3), x$C$loci_outcome1, x$C$loci_outcome2,
                x$C$loci_shared, cq_num(x$C$jaccard_locus, 3)))
  if (!is.null(x$D))
    cat(sprintf("D  attrition  : %d pathway genes -> %d instrumentable -> %d analysable -> %d associated\n",
                x$D$pathway_genes, x$D$instrumentable, x$D$analysable, x$D$associated))
  if (!is.null(x$E) && nrow(x$E))
    cat(sprintf("E  compartment: %s\n", paste(sprintf("%s %sx (%s)", x$E$gene,
        cq_num(x$E$ratio_other_over_target), x$E$highest_other), collapse = "; ")))
  if (!is.null(x$F) && nrow(x$F))
    cat(sprintf("F  peak dist  : %s\n", paste(sprintf("%s %s kb", x$F$gene,
        cq_fmt_kb(x$F$distance_bp)), collapse = "; ")))
  if (!is.null(x$G)) {
    ks <- x$G$known_sweep
    cat(sprintf("G  windows    : known %s kb -> fold %s\n",
                paste(ks$known_window_kb, collapse = "/"),
                paste(cq_num(ks$fold), collapse = "/")))
    n <- x$G$threshold_neighbourhood
    if (is.finite(n$nearest_below_bp) || is.finite(n$nearest_above_bp))
      cat(sprintf("   threshold  : %s kb; nearest significant-locus distances %s kb below, %s kb above\n",
                  cq_fmt_kb(n$threshold_bp),
                  if (is.finite(n$nearest_below_bp)) cq_fmt_kb(n$nearest_below_bp) else "none",
                  if (is.finite(n$nearest_above_bp)) cq_fmt_kb(n$nearest_above_bp) else "none"))
  }

  cat("\nModule status\n")
  for (k in names(x$module_status))
    cat(sprintf("  %-5s %s\n", k, x$module_status[[k]]))

  cat("\nPer-gene evidence: independent fields, not a ranking. See summary(x)$evidence.\n")
  cat("No gene can reach `target-supported` from this tool. Not run here:\n")
  for (m in x$not_automated) cat(sprintf("  - %s\n", m$name))
  cat("A nomination audited on five of eight diagnostics is audited on five of eight.\n")
  invisible(x)
}

#' Summarise an audit
#'
#' Returns a list, not a single table: `attribution` is module A only,
#' `modules` is the run status of every module, and `evidence` is the per-gene
#' fields. An earlier version returned only module A's counts while being named
#' as if it summarised the audit.
#' @param object a `cqtna_audit`.
#' @param ... ignored.
#' @return a list with `attribution`, `modules` and `evidence`.
#' @export
summary.cqtna_audit <- function(object, ...) {
  a <- object$A
  attribution <- data.frame(
    quantity = c("significant loci", "of those, on known loci",
                 "background loci", "of those, known",
                 "fold enrichment", "one-sided Fisher P"),
    value = c(a$significant_loci, a$significant_known, a$background_loci,
              a$background_known, signif(a$fold, 4), signif(a$fisher_p_one_sided, 3)),
    stringsAsFactors = FALSE)
  if (!is.null(object$A_nc))
    attribution <- rbind(attribution, data.frame(
      quantity = c("mismatched-list fold", "mismatched-list P"),
      value = c(signif(object$A_nc$fold, 4),
                signif(object$A_nc$fisher_p_one_sided, 3)),
      stringsAsFactors = FALSE))
  list(attribution = attribution,
       modules = data.frame(module = names(object$module_status),
                            status = unname(object$module_status),
                            stringsAsFactors = FALSE),
       evidence = object$evidence)
}

#' Write the audit as a markdown report
#'
#' Every module gets a section, including the ones that did not run, so nothing
#' is silently skipped.
#' @param x a `cqtna_audit`.
#' @param file path to write to; `NULL` returns the text.
#' @return the report text, invisibly when written to a file.
#' @export
cqtna_report <- function(x, file = NULL) {
  s <- x$settings
  st <- x$module_status
  L <- c(sprintf("# CQTNA report (R %s)", utils::packageVersion("cqtna")), "",
         sprintf("Settings: FDR < %s, locus window %s kb, known-locus window %s kb, build %s.",
                 s$fdr, s$locus_kb, s$known_kb, s$build))
  if (isTRUE(s$reclustered))
    L <- c(L, "", paste("Note: the supplied object had been clustered with a",
                        "different locus window and was re-clustered, so the",
                        "window above is the one actually used."))
  L <- c(L, "", "## A. Locus attribution", "",
         sprintf("- background: %d/%d loci (%s%%) carry a known lead SNP for this outcome",
                 x$A$background_known, x$A$background_loci, cq_num(x$A$background_pct)),
         sprintf("- significant: %d/%d (%s%%)", x$A$significant_known,
                 x$A$significant_loci, cq_num(x$A$pct_known, 1)),
         sprintf("- **enrichment %sx, one-sided Fisher P = %.3g**",
                 cq_num(x$A$fold), x$A$fisher_p_one_sided))
  if (!is.null(x$A_nc)) {
    clean <- is.na(x$A_nc$fold) || x$A_nc$fold <= 1 || x$A_nc$fisher_p_one_sided > 0.05
    L <- c(L, sprintf("- mismatched-list control: %sx, P = %.3g %s",
                      cq_num(x$A_nc$fold), x$A_nc$fisher_p_one_sided,
                      if (clean) "(clean)" else
                        "**(POSITIVE -- attribution may be a density artefact)**"))
  } else {
    L <- c(L, paste("- no mismatched-list control supplied; enrichment specific",
                    "to this outcome's genetics is not established"))
  }
  L <- c(L, "",
         paste("The Fisher test treats independent loci as exchangeable units.",
               "Loci are defined by single-linkage clustering, so they are not",
               "independent tests in the strict sense; read this as a descriptive",
               "enrichment with a sensitivity analysis (module G), not as a",
               "calibrated p-value."), "")

  L <- c(L, "## B. Inference unit", "", cq_md_table(x$B), "",
         paste("A shortest list is not evidence that FDR is controlled under",
               "dependence. Simes-then-BH across units is a sensitivity analysis,",
               "not a proof of FDR control under this correlation structure.",
               "State which unit the conclusions are in."), "")

  L <- c(L, "## C. List stability across outcome GWAS", "")
  if (identical(st[["C"]], "run")) {
    cov <- x$C$coverage
    L <- c(L,
           sprintf("- genes: %d vs %d, %d shared, Jaccard %s", x$C$genes_outcome1,
                   x$C$genes_outcome2, x$C$genes_shared, cq_num(x$C$jaccard_gene, 3)),
           sprintf("- loci: %d vs %d, %d shared, Jaccard %s (clustered jointly at %s kb)",
                   x$C$loci_outcome1, x$C$loci_outcome2, x$C$loci_shared,
                   cq_num(x$C$jaccard_locus, 3), x$C$locus_kb),
           sprintf("- lost: %s", if (length(x$C$lost)) paste(x$C$lost, collapse = ", ") else "none"),
           sprintf("- gained: %s", if (length(x$C$gained)) paste(x$C$gained, collapse = ", ") else "none"),
           sprintf("- coverage: %d gene-variant pairs shared, %d only in the first outcome, %d only in the second",
                   cov$gene_variant_pairs_shared, cov$gene_variant_pairs_only1,
                   cov$gene_variant_pairs_only2))
    if (length(cov$genes_never_tested_in_2))
      L <- c(L, sprintf("- never tested under the second outcome: %s",
                        paste(cov$genes_never_tested_in_2, collapse = ", ")))
    L <- c(L, "")
  } else {
    L <- c(L, sprintf("**%s.** This is the diagnostic that most often changes a conclusion.",
                      st[["C"]]), "")
  }

  for (key in c("D", "E", "F")) {
    ttl <- c(D = "D. Instrument attrition", E = "E. Compartment attribution",
             F = "F. eQTL-to-GWAS peak distance")[[key]]
    L <- c(L, paste0("## ", ttl), "")
    if (!identical(st[[key]], "run")) {
      L <- c(L, sprintf("**%s.**", st[[key]]), "")
      next
    }
    if (key == "D") {
      d <- x$D
      L <- c(L, sprintf("- of %d pathway genes: **%d instrumentable, %d analysable against this outcome, %d nominally associated**",
                        d$pathway_genes, d$instrumentable, d$analysable, d$associated), "")
    } else {
      tab <- if (key == "E") x$E else x$F
      if (key == "E") tab$ratio_other_over_target <- cq_num(tab$ratio_other_over_target)
      # as.character(50000) is "5e+04"; positions and distances must print as integers
      if (key == "F") for (cc in c("eqtl_pos", "gwas_pos", "distance_bp"))
        tab[[cc]] <- formatC(tab[[cc]], format = "d", big.mark = ",")
      L <- c(L, cq_md_table(tab), "")
    }
  }

  for (key in c("G", "G_nc")) {
    ttl <- if (key == "G") "G. Window sensitivity and the distances behind the flag"
           else "G (mismatched list). The same sweep against the wrong disease's loci"
    L <- c(L, paste0("## ", ttl), "")
    if (!identical(st[[key]], "run")) {
      L <- c(L, sprintf("**%s.**", st[[key]]), "")
      next
    }
    g <- x[[key]]
    L <- c(L, cq_md_table(cq_round_p(g$known_sweep)), "",
           cq_md_table(cq_round_p(g$locus_sweep)), "")
    if (length(g$significant_distances_bp)) {
      L <- c(L, sprintf("- distances of significant loci to the nearest known lead SNP (kb): %s",
                        paste(cq_fmt_kb(g$significant_distances_bp), collapse = ", ")),
             sprintf("- median %s kb, against %s kb for all testable loci",
                     cq_fmt_kb(g$significant_median_bp), cq_fmt_kb(g$background_median_bp)))
      sr <- g$significant_distances_sig_records_bp
      if (length(sr) && !isTRUE(all.equal(sr, g$significant_distances_bp)))
        L <- c(L, sprintf(paste("- over each locus's SIGNIFICANT records only (kb): %s",
                                "-- where this differs, the locus is near a known",
                                "lead SNP because of a record that is not itself",
                                "significant"),
                          paste(cq_fmt_kb(sr), collapse = ", ")))
    }
    n <- g$threshold_neighbourhood
    L <- c(L, sprintf("- threshold %s kb; nearest significant-locus distance below it %s kb, above it %s kb%s",
                      cq_fmt_kb(n$threshold_bp),
                      if (is.finite(n$nearest_below_bp)) cq_fmt_kb(n$nearest_below_bp) else "none",
                      if (is.finite(n$nearest_above_bp)) cq_fmt_kb(n$nearest_above_bp) else "none",
                      if (is.finite(n$log10_gap))
                        sprintf(" (log10 ratio %s)", cq_num(n$log10_gap, 2)) else ""),
           "", paste("No verdict is attached to that spacing. Whether the",
                     "threshold is comfortably clear of the data is a judgement",
                     "for the reader."), "")
    if (key == "G")
      L <- c(L, paste("A window chosen to flatter a result weakens when tightened.",
                      "A fold that rises as the window narrows means the reported",
                      "value is the conservative one for this fold -- it does not",
                      "make the wider inference conservative."), "")
  }

  L <- c(L, "## Per-gene evidence", "")
  if (nrow(x$evidence)) {
    e <- x$evidence
    e$compartment_ratio <- cq_num(e$compartment_ratio)
    L <- c(L, cq_md_table(e), "",
           paste("These are independent fields, not an ordered scale. A gene can",
                 "be on a known locus, unstable across outcome GWAS and dominated",
                 "by another compartment at the same time, and each of those",
                 "remains visible."), "")
  } else L <- c(L, "no significant genes at this threshold", "")

  L <- c(L, "## Module status", "", cq_md_table(
    data.frame(module = names(st), status = unname(st), stringsAsFactors = FALSE)), "",
    paste("**No gene can reach `target-supported` from this tool.**",
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
