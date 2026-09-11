## =====================================================================
## Step 50  Motif enrichment on the glycolysis axis, restricted to peaks that
##          do NOT change with activation.
##
## Step 48 found AP-1 enriched in glycolysis-high chromatin -- but with almost
## the same motif ranking as the activation positive control. Step 49 showed why:
## at the chromatin level the axis and the activation programme share 42% of
## their variance (r = 0.647; top-set overlap 3.9x), even though at the RNA level
## glycolysis is nearly orthogonal to the activation module (R2 = 0.024). Matching
## cells on an 11-gene RNA score did not control chromatin activation state.
##
## The AP-1 result therefore cannot be attributed to glycolysis as reported.
## Here the same test is repeated inside the 20,068 peaks whose accessibility does
## not track activation (|activation log2FC| <= 0.553, the median), which retain
## 68% of the axis's spread.
##
## READING THE RESULT
##   AP-1 persists   -> AP-1 is associated with the axis beyond activation
##   AP-1 disappears -> the Step 48 AP-1 signal was activation, and must be
##                      reported as such
## Either way the answer is informative, which is why the test is worth running.
## =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(Biostrings)
  library(TFBSTools); library(JASPAR2020); library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
set.seed(1)
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")
GEN <- BSgenome.Hsapiens.UCSC.hg38
N_FG <- 2000
N_BG_PER_FG <- 5

pfm <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))

peaks_to_gr <- function(x) {
  p <- tstrsplit(x, "[:-]")
  GRanges(p[[1]], IRanges(as.integer(p[[2]]), as.integer(p[[3]])))
}

inv <- fread("49c_activation_invariant_peaks.tsv.gz")
gr <- peaks_to_gr(inv$peak)
ok <- as.logical(seqnames(gr) %in% seqnames(GEN) & width(gr) > 20)
inv <- inv[ok]; gr <- gr[ok]
cat("activation-invariant peaks usable:", nrow(inv), "\n")

gc <- as.numeric(letterFrequency(getSeq(GEN, gr), "GC", as.prob = TRUE))
acc <- inv$total

match_bg <- function(fg_idx) {
  gcb <- cut(gc, quantile(gc, seq(0, 1, .05), na.rm = TRUE),
             include.lowest = TRUE, labels = FALSE)
  acb <- cut(log1p(acc), quantile(log1p(acc), seq(0, 1, .1), na.rm = TRUE),
             include.lowest = TRUE, labels = FALSE)
  key <- paste(gcb, acb)
  pool <- split(seq_along(key), key)
  out <- integer(0)
  for (i in fg_idx) {
    cand <- setdiff(pool[[key[i]]], fg_idx)
    if (length(cand)) out <- c(out, sample(cand, min(N_BG_PER_FG, length(cand))))
  }
  unique(out)
}

enrich <- function(fg_idx, label) {
  bg <- match_bg(fg_idx)
  cat(sprintf("  %s: %d foreground, %d matched background\n",
              label, length(fg_idx), length(bg)))
  idx <- c(fg_idx, bg)
  is_fg <- c(rep(TRUE, length(fg_idx)), rep(FALSE, length(bg)))
  mm <- as.matrix(motifMatches(matchMotifs(pfm, gr[idx], genome = GEN,
                                           out = "matches")))
  res <- rbindlist(lapply(seq_len(ncol(mm)), function(j) {
    a <- sum(mm[is_fg, j]); b <- sum(is_fg) - a
    c_ <- sum(mm[!is_fg, j]); d <- sum(!is_fg) - c_
    ft <- fisher.test(matrix(c(a, b, c_, d), 2), alternative = "greater")
    data.table(motif = name(pfm[[j]]), fg_frac = a / sum(is_fg),
               bg_frac = c_ / sum(!is_fg), odds = unname(ft$estimate),
               p = ft$p.value)
  }))
  res[, FDR := p.adjust(p, "BH")][, contrast := label][order(p)]
}

up <- head(order(-inv$mean_lfc), N_FG)
dn <- head(order(inv$mean_lfc), N_FG)
r_up <- enrich(up, "glyco_high_actInvariant")
r_dn <- enrich(dn, "glyco_low_actInvariant")
fwrite(rbind(r_up, r_dn), "50a_motif_activation_invariant.tsv", sep = "\t")

show <- function(r, ttl) {
  cat("\n", ttl, "\n", sep = "")
  print(head(r[, .(motif, fg_frac = round(fg_frac, 3), bg_frac = round(bg_frac, 3),
                   odds = round(odds, 2), FDR = signif(FDR, 3))], 15))
  cat(sprintf("  motifs at FDR<0.05: %d\n", r[FDR < 0.05, .N]))
}
show(r_up, "GLYCOLYSIS-HIGH, activation-invariant peaks")
show(r_dn, "GLYCOLYSIS-LOW, activation-invariant peaks")

## direct comparison with Step 48
prev <- fread("48b_motif_axis.tsv")
cmp <- merge(prev[contrast == "glyco_high", .(motif, odds_all = odds, FDR_all = FDR)],
             r_up[, .(motif, odds_inv = odds, FDR_inv = FDR)], by = "motif")
ap1 <- cmp[grepl("^(FOS|JUN|BATF|BACH|JDP)", motif)]
cat("\nAP-1 family: all peaks vs activation-invariant peaks\n")
print(head(ap1[order(-odds_all), .(motif, odds_all = round(odds_all, 2),
                                   FDR_all = signif(FDR_all, 2),
                                   odds_inv = round(odds_inv, 2),
                                   FDR_inv = signif(FDR_inv, 2))], 12))
cat(sprintf("\n  AP-1 motifs at FDR<0.05 -- all peaks: %d, activation-invariant: %d\n",
            ap1[FDR_all < 0.05, .N], ap1[FDR_inv < 0.05, .N]))
cat(sprintf("  median AP-1 odds -- all peaks: %.2f, activation-invariant: %.2f\n",
            median(ap1$odds_all), median(ap1$odds_inv)))
fwrite(cmp, "50b_motif_comparison_all_vs_invariant.tsv", sep = "\t")
cat("\nwritten: 50a, 50b\n")
