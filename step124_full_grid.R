## Step 124 -- the complete grid on one statistical unit.
##
## Executes the frozen decision record, manuscript/PREREG_locus_partition.md:
##   main analysis    non-recursive fixed-centre partition, max span 1000 kb,
##                    attribution convention C1 (any_record)
##   sensitivity      the same partition at 500 kb
##   legacy           1 Mb single-linkage with the published convention (C2),
##                    reported for reproduction only, carrying no inference
##
## Every cell is recomputed from a per-record table, so one grid uses one
## definition of its unit. The two cells whose per-record tables were never saved
## -- RA on both resources, and whole blood against HCC -- were re-derived by
## step108 and step94d, which now persist 123b and 123c.
##
## Two columns state each row's standing outright, so that no later reader has
## to reconstruct it from a cell name:
##   region_filter    all_genome | MHC_excluded
##   analysis_status  prereg_primary      registered cell, registered region,
##                                        registered partition
##                    prereg_sensitivity  variation fixed in writing before the
##                                        result was seen (500 kb; HCC-low)
##                    posthoc_sensitivity  chosen after seeing a result
##                                        (the MHC-excluded RA rows)
##                    legacy_reproduction  old algorithm, reproduction only
##
##   Rscript step124_full_grid.R
## Output: 123d_fixed_anchor_full_grid.tsv

suppressMessages(library(cqtna))
MR <- "D:/R_ex/MR"
setwd(MR)

source("step124_cells.R")   # cells, runs, status_of()

rows <- list()
for (cl in cells) for (rn in runs) {
  if (!nrow(cl$d)) next
  mr <- as_cqtna_mr(cl$d, locus_kb = rn$kb, build = "GRCh38",
                    locus_method = rn$method)
  sp <- suppressWarnings(cqtna_locus_spans(mr))
  a <- cqtna_attribution(mr, cl$m, known_from = rn$conv)
  x <- cqtna_attribution(mr, cl$x, known_from = rn$conv)
  rows[[length(rows) + 1L]] <- data.frame(
    cell = cl$name, role = cl$role, analysis = rn$label,
    region_filter = cl$region, analysis_status = status_of(cl, rn),
    partition = rn$method, locus_kb = rn$kb, known_from = rn$conv,
    n_records = nrow(cl$d),
    bg_loci = a$background_loci, bg_known = a$background_known,
    sig_loci = a$significant_loci, sig_known = a$significant_known,
    pct_known = round(a$pct_known, 1),
    fold = round(a$fold, 2), fisher_p = a$fisher_p_one_sided,
    mismatch_fold = round(x$fold, 2), mismatch_p = x$fisher_p_one_sided,
    control = if (is.finite(x$fisher_p_one_sided) && x$fisher_p_one_sided < 0.05 &&
                  isTRUE(x$fold > 1)) "FAILED" else "clean",
    max_sig_span_kb = round(max(sp$significant_span_kb)),
    n_sig_over_window = sum(sp$significant_span_kb > rn$kb),
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
write.table(out, "123d_fixed_anchor_full_grid.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

show <- function(lbl) {
  s <- out[out$analysis == lbl, ]
  cat("\n", strrep("=", 118), "\n", toupper(lbl), " -- ", s$partition[1], " ",
      s$locus_kb[1], " kb, ", s$known_from[1], "\n", sep = "")
  cat(sprintf("%-36s %-14s %-19s %9s %9s %8s %9s   %8s %10s %8s\n",
              "cell", "region", "status", "bg known", "sig known", "fold", "P",
              "mism.", "mism. P", "control"))
  for (i in seq_len(nrow(s))) with(s[i, ], cat(sprintf(
    "%-36s %-14s %-19s %4d/%-4d %4d/%-4d %8.2f %9.3g   %8.2f %10.3g %8s\n",
    cell, region_filter, analysis_status, bg_known, bg_loci, sig_known,
    sig_loci, fold, fisher_p, mismatch_fold, mismatch_p, control)))
  m <- s[s$role == "main", ]
  ok <- m[m$control == "clean", ]
  cat(sprintf("  -> %d registered cells; %d VOID on the mismatched control\n",
              nrow(m), sum(m$control == "FAILED")))
  if (any(m$control == "FAILED"))
    cat(sprintf("     void: %s\n",
                paste(m$cell[m$control == "FAILED"], collapse = ", ")))
  cat(sprintf("  -> of the %d usable cells: %d with fold > 1, %d reach P < 0.05\n",
              nrow(ok), sum(ok$fold > 1, na.rm = TRUE),
              sum(ok$fisher_p < 0.05, na.rm = TRUE)))
  cat(sprintf("  -> widest significant locus %s kb; loci over the window: %d\n",
              format(max(s$max_sig_span_kb), big.mark = ","),
              sum(s$n_sig_over_window)))
}
for (l in c("main", "sensitivity", "legacy")) show(l)
cat("\nwrote 123d_fixed_anchor_full_grid.tsv\n")
