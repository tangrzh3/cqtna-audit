## Step 125 -- why the mismatched-list control fails, one locus at a time.
##
## The grid reports the mismatched control as a single fold and P per cell. That
## number cannot say WHICH loci carry it, and the two candidate mechanisms leave
## very different traces:
##
##   (a) shared biology / shared density -- a significant record of ours really
##       does sit within 1 Mb of a lead SNP for the wrong disease;
##   (b) C1 spread -- no significant record in the locus is anywhere near the
##       wrong disease's list, and the locus is counted as known only because
##       some NON-significant record elsewhere in the same block is. Convention
##       C1 gives known status to the block, and the block is up to 1000 kb wide.
##
## (b) is a property of the convention, not of the genome, and it is invisible in
## the summary fold. This script separates them.
##
## The distance test here is written from scratch rather than taken from the
## package, so the counts are an independent recomputation; they are checked
## against cqtna_attribution() and the run aborts on any disagreement.
##
##   Rscript step125_mismatch_loci.R
## Output: 125a_mismatch_locus_diagnostics.tsv  (one row per flagged locus)
##         125b_mismatch_spread_summary.tsv     (one row per cell x window)

## Work from wherever this script lives, so the packet runs after extraction.
## Override with:  Rscript <script> /path/to/dir     or  CQTNA_DIR=/path/to/dir
MR <- local({
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) && nzchar(a[1])) return(a[1])
  if (nzchar(Sys.getenv("CQTNA_DIR"))) return(Sys.getenv("CQTNA_DIR"))
  f <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", f[grepl("^--file=", f)])
  if (length(f)) normalizePath(dirname(f[1])) else getwd()
})
setwd(MR)
source("step124_cells.R")   # cells, runs, status_of()

KNOWN_KB <- 1000   # cqtna_attribution()'s default: how near a lead SNP counts
FDR      <- 0.05

## Distance from each record to the nearest lead SNP on its own chromosome.
## Independent of the package: a sort and a binary search, nothing else.
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

wins <- Filter(function(r) r$method == "fixed_centre", runs)

rows <- list(); summ <- list()
for (cl in cells) for (rn in wins) {
  if (!nrow(cl$d)) next
  mr <- as_cqtna_mr(cl$d, locus_kb = rn$kb, build = "GRCh38",
                    locus_method = rn$method)
  d <- data.frame(locus = as.character(mr$locus), chr = as.character(mr$chr),
                  pos = as.numeric(mr$pos), gene = as.character(mr$gene),
                  p = as.numeric(mr$p), fdr = as.numeric(mr$fdr),
                  stringsAsFactors = FALSE)
  d$dist <- nearest_bp(cl$x, d$chr, d$pos)        # x = the MISMATCHED list
  d$near <- d$dist <= KNOWN_KB * 1000
  d$sig  <- d$fdr < FDR

  ## C1: a locus is known if ANY of its records is near. Restrict to loci that
  ## contain at least one significant record -- those are the numerator.
  by <- split(seq_len(nrow(d)), d$locus)
  sig_loci <- names(by)[vapply(by, function(i) any(d$sig[i]), logical(1))]
  known_c1 <- vapply(by[sig_loci], function(i) any(d$near[i]), logical(1))

  ## Cross-check against the package before reporting anything derived from it.
  chk <- cqtna_attribution(mr, cl$x, known_from = rn$conv)
  stopifnot(length(sig_loci) == chk$significant_loci,
            sum(known_c1) == chk$significant_known)

  flagged <- sig_loci[known_c1]
  for (L in flagged) {
    i <- by[[L]]
    sg <- i[d$sig[i]]
    n_sig_near <- sum(d$near[sg])
    ## distance from the locus's own significant records to the wrong list
    dmin_sig <- min(d$dist[sg])
    top <- sg[which.min(d$p[sg])]
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cl$name, region_filter = cl$region,
      analysis_status = status_of(cl, rn), locus_kb = rn$kb, locus = L,
      chr = d$chr[i][1],
      start_mb = round(min(d$pos[i]) / 1e6, 3),
      end_mb   = round(max(d$pos[i]) / 1e6, 3),
      span_kb  = round((max(d$pos[i]) - min(d$pos[i])) / 1000),
      n_records = length(i), n_sig = length(sg),
      n_near = sum(d$near[i]), n_sig_near = n_sig_near,
      min_dist_sig_kb = round(dmin_sig / 1000),
      evidence = if (n_sig_near > 0) "significant_record" else "C1_spread_only",
      top_gene = d$gene[top], top_p = d$p[top],
      stringsAsFactors = FALSE)
  }
  n_spread <- sum(vapply(flagged, function(L) {
    i <- by[[L]]; !any(d$near[i[d$sig[i]]])
  }, logical(1)))
  summ[[length(summ) + 1L]] <- data.frame(
    cell = cl$name, region_filter = cl$region,
    analysis_status = status_of(cl, rn), locus_kb = rn$kb,
    mismatch_fold = round(chk$fold, 2), mismatch_p = chk$fisher_p_one_sided,
    control = if (is.finite(chk$fisher_p_one_sided) &&
                  chk$fisher_p_one_sided < 0.05 && isTRUE(chk$fold > 1))
                "FAILED" else "clean",
    sig_loci = length(sig_loci), mismatch_known = length(flagged),
    spread_only = n_spread,
    pct_spread = if (length(flagged)) round(100 * n_spread / length(flagged), 1)
                 else NA_real_,
    stringsAsFactors = FALSE)
}

det <- do.call(rbind, rows)
det <- det[order(det$locus_kb, det$cell, det$evidence != "C1_spread_only",
                 det$chr, det$start_mb), ]
sm  <- do.call(rbind, summ)
write.table(det, "125a_mismatch_locus_diagnostics.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(sm, "125b_mismatch_spread_summary.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

for (kb in unique(sm$locus_kb)) {
  cat("\n", strrep("=", 104), "\n", kb,
      " kb -- mismatched-known significant loci, split by what makes them known\n",
      sep = "")
  s <- sm[sm$locus_kb == kb, ]
  cat(sprintf("%-36s %-19s %8s %9s %9s %11s %9s\n", "cell", "status",
              "control", "sig loci", "mism.known", "C1 spread", "% spread"))
  for (i in seq_len(nrow(s))) with(s[i, ], cat(sprintf(
    "%-36s %-19s %8s %8d %9d %11d %8s\n", cell, analysis_status, control,
    sig_loci, mismatch_known, spread_only,
    ifelse(is.na(pct_spread), "-", sprintf("%.1f", pct_spread)))))
}

cat("\n", strrep("=", 104), "\n",
    "loci known to the WRONG list with no significant record within 1 Mb of it\n",
    sep = "")
sp <- det[det$evidence == "C1_spread_only", ]
cat(sprintf("%-36s %6s %-5s %16s %8s %7s %8s %14s\n", "cell", "kb", "chr",
            "span (Mb)", "span kb", "n sig", "n near", "nearest sig"))
for (i in seq_len(nrow(sp))) with(sp[i, ], cat(sprintf(
  "%-36s %6d %-5s %7.2f-%-8.2f %8d %7d %8d %11d kb\n",
  cell, locus_kb, chr, start_mb, end_mb, span_kb, n_sig, n_near,
  min_dist_sig_kb)))

cat("\nwrote 125a_mismatch_locus_diagnostics.tsv, 125b_mismatch_spread_summary.tsv\n")
