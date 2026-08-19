## Step 128 -- window sensitivity, redone on the frozen partition and extended
## to every cell.
##
## S36 (SUPP_window_sensitivity.md) swept two different "1 Mb" conventions:
##
##   KNOWN_KB  how near a lead SNP a locus must be to count as known
##   LOCUS_KB  how far apart two variants may be and still be one locus
##
## Both sweeps were run under single linkage, and both skipped the two RA cells
## because step108 had not persisted its per-record table. Neither limit still
## holds: the partition is frozen as non-recursive fixed-anchor (S37), and 123b
## now carries RA per record. Redone here so the sensitivity analysis is on the
## same unit as the estimate it is meant to bound.
##
## LOCUS_KB means something different under the two rules, and the difference is
## the point. Under single linkage it is a joining distance and the span is
## unbounded; under fixed anchor it IS the maximum span. So the sweep below is a
## sweep of the maximum locus width, which is the quantity a reader thinks they
## are being shown.
##
##   Rscript step128_window_sensitivity_fixed_anchor.R
## Output: 128a_known_kb_sweep.tsv   KNOWN_KB swept, partition held at 1000 kb
##         128b_locus_kb_sweep.tsv   LOCUS_KB swept, known window held at 1000 kb
##         128c_distances.tsv        per significant locus, distance to the
##                                   nearest lead SNP of its own disease

MR <- "D:/R_ex/MR"
setwd(MR)
source("step124_cells.R")   # cells, runs, status_of()

SWEEP <- c(100, 250, 500, 1000)
FDR   <- 0.05
CONV  <- "any_record"

## The six registered cells plus the two power-sensitivity cells; the post-hoc
## MHC-excluded rows are not swept, because a sensitivity analysis of a
## post-hoc analysis is not something the manuscript can quote.
use <- Filter(function(cl) cl$region == "all_genome", cells)

nearest_bp <- function(known, chr, pos) {
  chr <- as.character(chr); pos <- as.numeric(pos)
  out <- rep(Inf, length(pos))
  by_chr <- split(as.numeric(known$pos), as.character(known$chr))
  for (ch in unique(chr)) {
    arr <- sort(by_chr[[ch]])
    if (!length(arr)) next
    i <- which(chr == ch); p <- pos[i]
    j <- findInterval(p, arr)
    lo <- ifelse(j >= 1L, abs(arr[pmax(j, 1L)] - p), Inf)
    hi <- ifelse(j < length(arr), abs(arr[pmin(j + 1L, length(arr))] - p), Inf)
    out[i] <- pmin(lo, hi)
  }
  out
}

## ------------------------------------------------- A. KNOWN_KB, partition fixed
rows <- list()
for (cl in use) {
  mr <- as_cqtna_mr(cl$d, locus_kb = 1000, build = "GRCh38",
                    locus_method = "fixed_centre")
  for (kb in SWEEP) {
    a <- cqtna_attribution(mr, cl$m, known_kb = kb, known_from = CONV)
    x <- cqtna_attribution(mr, cl$x, known_kb = kb, known_from = CONV)
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cl$name, role = cl$role, known_kb = kb, locus_kb = 1000,
      bg_loci = a$background_loci, bg_known = a$background_known,
      sig_loci = a$significant_loci, sig_known = a$significant_known,
      fold = round(a$fold, 2), fisher_p = a$fisher_p_one_sided,
      mismatch_fold = round(x$fold, 2), mismatch_p = x$fisher_p_one_sided,
      control = if (is.finite(x$fisher_p_one_sided) &&
                    x$fisher_p_one_sided < 0.05 && isTRUE(x$fold > 1))
                  "FAILED" else "clean",
      stringsAsFactors = FALSE)
  }
}
kn <- do.call(rbind, rows)
write.table(kn, "128a_known_kb_sweep.tsv", sep = "\t", row.names = FALSE,
            quote = FALSE)

## ------------------------------------------------- B. LOCUS_KB, known window fixed
rows <- list()
for (cl in use) {
  for (kb in SWEEP) {
    mr <- as_cqtna_mr(cl$d, locus_kb = kb, build = "GRCh38",
                      locus_method = "fixed_centre")
    a <- cqtna_attribution(mr, cl$m, known_kb = 1000, known_from = CONV)
    x <- cqtna_attribution(mr, cl$x, known_kb = 1000, known_from = CONV)
    sp <- suppressWarnings(cqtna_locus_spans(mr))
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cl$name, role = cl$role, known_kb = 1000, locus_kb = kb,
      bg_loci = a$background_loci, bg_known = a$background_known,
      sig_loci = a$significant_loci, sig_known = a$significant_known,
      fold = round(a$fold, 2), fisher_p = a$fisher_p_one_sided,
      mismatch_fold = round(x$fold, 2), mismatch_p = x$fisher_p_one_sided,
      control = if (is.finite(x$fisher_p_one_sided) &&
                    x$fisher_p_one_sided < 0.05 && isTRUE(x$fold > 1))
                  "FAILED" else "clean",
      max_sig_span_kb = round(max(sp$significant_span_kb)),
      stringsAsFactors = FALSE)
  }
}
lk <- do.call(rbind, rows)
write.table(lk, "128b_locus_kb_sweep.tsv", sep = "\t", row.names = FALSE,
            quote = FALSE)

## ------------------------------------------------- C. continuous distances
## Taken over SIGNIFICANT records only. The first version of S36 took the
## minimum over every record at the locus, which reported a 235 kb hit for
## HCC-high that belonged to a non-significant record.
rows <- list()
for (cl in use) {
  mr <- as_cqtna_mr(cl$d, locus_kb = 1000, build = "GRCh38",
                    locus_method = "fixed_centre")
  d <- data.frame(locus = as.character(mr$locus), chr = as.character(mr$chr),
                  pos = as.numeric(mr$pos), gene = as.character(mr$gene),
                  p = as.numeric(mr$p), fdr = as.numeric(mr$fdr),
                  stringsAsFactors = FALSE)
  d$dist <- nearest_bp(cl$m, d$chr, d$pos)
  by <- split(seq_len(nrow(d)), d$locus)
  sig_loci <- names(by)[vapply(by, function(i) any(d$fdr[i] < FDR), logical(1))]
  bg_med <- median(vapply(by, function(i) min(d$dist[i]), numeric(1)))
  for (L in sig_loci) {
    i <- by[[L]]; sg <- i[d$fdr[i] < FDR]
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cl$name, locus = L, chr = d$chr[i][1],
      dist_kb = round(min(d$dist[sg]) / 1000),
      top_gene = d$gene[sg[which.min(d$p[sg])]],
      background_median_kb = round(bg_med / 1000),
      stringsAsFactors = FALSE)
  }
}
ds <- do.call(rbind, rows)
ds <- ds[order(ds$cell, ds$dist_kb), ]
write.table(ds, "128c_distances.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## ------------------------------------------------------------------ report
show <- function(tab, swept, held) {
  cat("\n", strrep("=", 108), "\n", swept, " swept, ", held, " held\n", sep = "")
  for (nm in unique(tab$cell)) {
    s <- tab[tab$cell == nm, ]
    cat(sprintf("\n  %-30s", nm))
    cat(sprintf("\n    %8s %14s %10s %11s %9s\n", swept, "bg known", "sig known",
                "fold", "control"))
    for (i in seq_len(nrow(s))) with(s[i, ], cat(sprintf(
      "    %8d %7d/%-6d %5d/%-4d %6s  P=%-8.3g %8s\n",
      if (swept == "known_kb") known_kb else locus_kb,
      bg_known, bg_loci, sig_known, sig_loci,
      ifelse(is.na(fold), "-", sprintf("%.2f", fold)), fisher_p, control)))
  }
}
show(kn, "known_kb", "locus_kb = 1000")
show(lk, "locus_kb", "known_kb = 1000")

cat("\n", strrep("=", 108), "\n",
    "distance from each significant locus to the nearest lead SNP of its own disease\n",
    sep = "")
for (nm in unique(ds$cell)) {
  s <- ds[ds$cell == nm, ]
  cat(sprintf("  %-30s n=%2d  |  %s  kb   (background median %s kb)\n",
              nm, nrow(s), paste(s$dist_kb, collapse = " "),
              format(s$background_median_kb[1], big.mark = ",")))
}
cat("\nwrote 128a / 128b / 128c\n")
