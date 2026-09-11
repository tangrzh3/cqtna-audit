## =====================================================================
## Step 27  Does GSE115978 carry the same double-normalisation error?
##
## GSE115978_TCD4.rds was built by the same pipeline that produced the faulty
## GSE120575 object (Step 25). Step 18's conclusions rest on it:
##   SMC2 ~ proliferation  rho=0.597 p=0.015 (FDR 0.088)  <- "most robust sc finding"
##   TPI1 ~ exhaustion(post) rho=0.717 p=0.030
##   TPI1 ~ proliferation  rho=0.182 ns  (did not replicate GSE120575)
##   HLA-C ~ proliferation rho=+0.106 ns (direction reversed -> HLA-C withdrawn)
##   detection: TPI1 60.6%, HLA-C 99.9%, SPSB2 5.0%, SMC2 10.0%, ZFYVE19 11.0%
##
## This script (1) diagnoses the scale of both the raw file and the rds, and
## (2) recomputes the Step 18 analyses from the raw file on the correct scale.
##
## GSE115978 has NO response annotation, so only correlations, detection rates
## and treatment-naive vs post comparisons are available.
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat); library(Matrix)})
set.seed(1)
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

CNT <- "D:/数据/黑色素瘤单细胞人/GSE115978/GSE115978_counts.csv.gz"
ANN <- "D:/数据/黑色素瘤单细胞人/GSE115978/GSE115978_cell.annotations.csv.gz"

CAND    <- c("ZFYVE19", "SMC2", "SPSB2", "TPI1", "HLA-C", "KIAA0040")
PROLIF  <- c("MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA")
EXHAUST <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX")
MIN_CELLS <- 20

## =====================================================================
## 1. scale diagnosis
## =====================================================================
cat("=== 1. scale of the raw file ===\n")
raw <- fread(CNT)
genes <- as.character(raw[[1]])
m <- as.matrix(raw[, -1]); rownames(m) <- make.unique(genes)
storage.mode(m) <- "double"
rm(raw); gc()
cat("dim:", nrow(m), "x", ncol(m), "\n")
sub <- sample(ncol(m), min(200, ncol(m)))
cat(sprintf("max = %.2f | raw per-cell sum (median) = %.3g | reverse-log per-cell sum = %.3g\n",
            max(m), median(colSums(m[, sub])), median(colSums(2^m[, sub] - 1))))
is_log <- max(m) <= 25 && median(colSums(2^m[, sub] - 1)) > 3e5
cat("verdict:", if (is_log) "already log2(TPM+1) -- do not re-log" else
                 "raw TPM -- apply log2(x+1)", "\n")
dat <- if (is_log) m else log2(m + 1)
rm(m); gc()

cat("\n=== 2. scale of the stored rds ===\n")
old <- readRDS("GSE115978_TCD4.rds")
do <- GetAssayData(old, assay = "RNA", layer = "data")
cat(sprintf("rds: %d x %d | max = %.2f | reverse-log per-cell sum = %.3g\n",
            nrow(do), ncol(do), max(do),
            median(colSums(2^as.matrix(do[, sample(ncol(do), min(200, ncol(do)))]) - 1))))
common_c <- intersect(colnames(do), colnames(dat))
common_g <- intersect(rownames(do), rownames(dat))
cat("overlap:", length(common_c), "cells,", length(common_g), "genes\n")
if (length(common_c) > 50 && length(common_g) > 500) {
  gs <- sample(common_g, 2000); cs <- sample(common_c, 200)
  vo <- as.numeric(do[gs, cs]); vn <- as.numeric(dat[gs, cs])
  ok <- vo > 0 & vn > 0
  rho <- cor(vo[ok], vn[ok], method = "spearman")
  cat(sprintf("nonzero-value agreement with correctly scaled data: rho = %.4f\n", rho))
  cat(if (rho > 0.999) "  -> rds matches the raw file; Step 18 numbers are unaffected\n"
      else "  -> *** rds does NOT match; it was renormalised. Step 18 must be redone ***\n")
}
rm(old, do); gc()

## =====================================================================
## 3. recompute Step 18 on the correct scale
## =====================================================================
cat("\n=== 3. recomputation on the correct scale ===\n")
ann <- fread(ANN); setnames(ann, 1, "cell")
ann <- ann[cell %in% colnames(dat)]
cat("cell types:\n"); print(table(ann$cell.types))
cd4 <- ann[cell.types == "T.CD4"]
cat("T.CD4 cells:", nrow(cd4), "\n")
d <- dat[, cd4$cell, drop = FALSE]

zsc <- function(genes) {
  g <- intersect(genes, rownames(d))
  x <- t(as.matrix(d[g, , drop = FALSE]))
  z <- scale(x); z[!is.finite(z)] <- 0
  rowMeans(z)
}

cell <- data.table(cell = cd4$cell, sample = cd4$samples,
                   treatment = cd4$treatment.group,
                   ProlifScore = zsc(PROLIF), ExhaustScore = zsc(EXHAUST))
for (g in CAND) cell[[g]] <- if (g %in% rownames(d)) as.numeric(d[g, ]) else NA_real_

detn <- sapply(CAND, function(g) if (g %in% rownames(d)) 100*mean(d[g, ] > 0) else NA)
cat("\ndetection in T.CD4 (%):\n"); print(round(detn, 1))
cat("Step 18 reported: TPI1 60.6, HLA-C 99.9, SPSB2 5.0, SMC2 10.0, ZFYVE19 11.0\n")
fwrite(data.table(gene = CAND, detection_pct = round(detn, 1)),
       "33a_GSE115978_detection_corrected.tsv", sep = "\t")

samp <- cell[, c(.(n = .N, treatment = treatment[1]), lapply(.SD, mean)),
             by = sample,
             .SDcols = c(CAND, "ProlifScore", "ExhaustScore")][n >= MIN_CELLS]
cat("\nsamples with >=", MIN_CELLS, "T.CD4 cells:", nrow(samp), "\n")
print(samp[, .N, by = treatment])
fwrite(samp, "33b_GSE115978_sample_level_corrected.tsv", sep = "\t")

cat("\n--- correlations at sample level ---\n")
cors <- rbindlist(lapply(CAND, function(g) rbindlist(lapply(
  c("ProlifScore", "ExhaustScore"), function(v) {
    ct <- suppressWarnings(cor.test(samp[[g]], samp[[v]], method = "spearman"))
    data.table(gene = g, score = v, n = nrow(samp),
               rho = unname(ct$estimate), p = ct$p.value)
  }))))
cors[, FDR := p.adjust(p, "BH"), by = score]
print(cors[order(score, p)], digits = 3)
fwrite(cors, "33c_GSE115978_correlations_corrected.tsv", sep = "\t")

## post-treatment only, matching the Step 18 "TPI1 ~ exhaustion (post)" claim
post <- samp[grepl("post", treatment, ignore.case = TRUE)]
if (nrow(post) >= 5) {
  cat("\n--- post-treatment samples only (n =", nrow(post), ") ---\n")
  cp <- rbindlist(lapply(CAND, function(g) rbindlist(lapply(
    c("ProlifScore", "ExhaustScore"), function(v) {
      ct <- suppressWarnings(cor.test(post[[g]], post[[v]], method = "spearman"))
      data.table(gene = g, score = v, n = nrow(post),
                 rho = unname(ct$estimate), p = ct$p.value)
    }))))
  print(cp[order(score, p)], digits = 3)
  fwrite(cp, "33d_GSE115978_correlations_post_corrected.tsv", sep = "\t")
}

cat("\n--- Step 18 reference values ---\n")
cat("  SMC2 ~ prolif  rho=0.597 p=0.015 (FDR 0.088)\n")
cat("  TPI1 ~ exhaust (post) rho=0.717 p=0.030\n")
cat("  TPI1 ~ prolif  rho=0.182 ns\n")
cat("  HLA-C ~ prolif rho=+0.106 ns\n")
cat("\nwritten: 33a-33d\n")
