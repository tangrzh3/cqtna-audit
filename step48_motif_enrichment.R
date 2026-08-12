## =====================================================================
## Step 48  Motif enrichment along the residual glycolysis axis.
##
## The axis has cleared two controls (independent of activation, R2=0.024;
## not a nuclear/cytoplasmic artefact, 4.7% of variance) and its transcriptional
## poles are anabolic growth versus T-cell identity/quiescence. This asks which
## transcription-factor motifs sit in the chromatin that differs along it.
##
## *** POSITIVE CONTROL RUNS FIRST ***
## The same pipeline is applied to rest vs 15 h activation. That contrast must
## recover AP-1 / NFAT / NF-kB opening, which is textbook. If it does not, the
## pipeline is not working and the axis result is not interpretable -- the same
## rule applied in Steps 19, 24 and 43.
##
## Background peaks are matched on GC content and total accessibility, because
## motif content tracks GC and naive backgrounds inflate GC-rich motifs.
## =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(Biostrings)
  library(TFBSTools); library(JASPAR2020); library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
set.seed(1)
setwd("D:/R_ex/MR")
GEN <- BSgenome.Hsapiens.UCSC.hg38
N_FG <- 2000        # peaks per direction
N_BG_PER_FG <- 5    # matched background peaks per foreground peak

pfm <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))
cat("JASPAR CORE vertebrate motifs:", length(pfm), "\n")

peaks_to_gr <- function(x) {
  p <- tstrsplit(x, "[:-]")
  GRanges(p[[1]], IRanges(as.integer(p[[2]]), as.integer(p[[3]])))
}

gc_of <- function(gr) {
  s <- getSeq(GEN, gr)
  as.numeric(letterFrequency(s, "GC", as.prob = TRUE))
}

match_background <- function(fg_idx, gc, acc, n_per) {
  gcb <- cut(gc, quantile(gc, seq(0, 1, .05), na.rm = TRUE),
             include.lowest = TRUE, labels = FALSE)
  acb <- cut(log1p(acc), quantile(log1p(acc), seq(0, 1, .1), na.rm = TRUE),
             include.lowest = TRUE, labels = FALSE)
  key <- paste(gcb, acb)
  pool <- split(seq_along(key), key)
  out <- integer(0)
  for (i in fg_idx) {
    cand <- setdiff(pool[[key[i]]], fg_idx)
    if (!length(cand)) next
    out <- c(out, sample(cand, min(n_per, length(cand))))
  }
  unique(out)
}

enrich <- function(all_gr, gc, acc, fg_idx, label) {
  bg_idx <- match_background(fg_idx, gc, acc, N_BG_PER_FG)
  cat(sprintf("  %s: %d foreground, %d matched background peaks\n",
              label, length(fg_idx), length(bg_idx)))
  idx <- c(fg_idx, bg_idx)
  is_fg <- c(rep(TRUE, length(fg_idx)), rep(FALSE, length(bg_idx)))
  hits <- matchMotifs(pfm, all_gr[idx], genome = GEN, out = "matches")
  mm <- as.matrix(motifMatches(hits))
  res <- rbindlist(lapply(seq_len(ncol(mm)), function(j) {
    a <- sum(mm[is_fg, j]); b <- sum(is_fg) - a
    c_ <- sum(mm[!is_fg, j]); d <- sum(!is_fg) - c_
    ft <- fisher.test(matrix(c(a, b, c_, d), 2), alternative = "greater")
    data.table(motif = name(pfm[[j]]), id = ID(pfm[[j]]),
               fg_frac = a / sum(is_fg), bg_frac = c_ / sum(!is_fg),
               odds = unname(ft$estimate), p = ft$p.value)
  }))
  res[, FDR := p.adjust(p, "BH")]
  res[, contrast := label]
  res[order(p)]
}

## ---------------------------------------------------------------- control
cat("\n=== POSITIVE CONTROL: rest vs 15 h activation ===\n")
ctl <- fread("47b_activation_differential_peaks.tsv.gz")
ctl <- ctl[is.finite(lfc)]
gr_c <- peaks_to_gr(ctl$peak)
ok <- seqnames(gr_c) %in% seqnames(GEN) & width(gr_c) > 20
ctl <- ctl[as.logical(ok)]; gr_c <- gr_c[as.logical(ok)]
gcc <- gc_of(gr_c)
accc <- ctl$rest_cpm + ctl$act_cpm
fg_open <- head(order(-ctl$lfc), N_FG)
ctl_res <- enrich(gr_c, gcc, accc, fg_open, "activation_opening")
fwrite(ctl_res, "48a_motif_activation_control.tsv", sep = "\t")
cat("\n  top 15 motifs in peaks opening on activation:\n")
print(head(ctl_res[, .(motif, fg_frac = round(fg_frac, 3),
                       bg_frac = round(bg_frac, 3),
                       odds = round(odds, 2), FDR = signif(FDR, 3))], 15))
expect <- c("FOS", "JUN", "BATF", "NFATC", "NFKB", "REL", "AP-1")
hit <- ctl_res[FDR < 0.05][grepl(paste(expect, collapse = "|"), motif,
                                 ignore.case = TRUE)]
cat(sprintf("\n  AP-1 / NFAT / NF-kB family motifs at FDR<0.05: %d\n", nrow(hit)))
if (nrow(hit) == 0) {
  cat("  *** POSITIVE CONTROL FAILED -- axis results below are not "  ,
      "interpretable ***\n")
}

## ---------------------------------------------------------------- axis
cat("\n=== AXIS: glycolysis-high vs -low within activated cells ===\n")
ax <- fread("47a_axis_differential_peaks.tsv.gz")
ax <- ax[consistent == TRUE & total >= 100 & is.finite(mean_lfc)]
gr_a <- peaks_to_gr(ax$peak)
ok <- seqnames(gr_a) %in% seqnames(GEN) & width(gr_a) > 20
ax <- ax[as.logical(ok)]; gr_a <- gr_a[as.logical(ok)]
cat("  usable peaks:", nrow(ax), "\n")
gca <- gc_of(gr_a)
acca <- ax$total

up <- head(order(-ax$mean_lfc), N_FG)
dn <- head(order(ax$mean_lfc), N_FG)
r_up <- enrich(gr_a, gca, acca, up, "glyco_high")
r_dn <- enrich(gr_a, gca, acca, dn, "glyco_low")
fwrite(rbind(r_up, r_dn), "48b_motif_axis.tsv", sep = "\t")

cat("\n  top 15 motifs in peaks more accessible in GLYCOLYSIS-HIGH cells:\n")
print(head(r_up[, .(motif, fg_frac = round(fg_frac, 3),
                    bg_frac = round(bg_frac, 3),
                    odds = round(odds, 2), FDR = signif(FDR, 3))], 15))
cat("\n  top 15 motifs in peaks more accessible in GLYCOLYSIS-LOW cells:\n")
print(head(r_dn[, .(motif, fg_frac = round(fg_frac, 3),
                    bg_frac = round(bg_frac, 3),
                    odds = round(odds, 2), FDR = signif(FDR, 3))], 15))

cat(sprintf("\n  motifs at FDR<0.05: high %d, low %d\n",
            r_up[FDR < 0.05, .N], r_dn[FDR < 0.05, .N]))
cat("\nwritten: 48a, 48b\n")
