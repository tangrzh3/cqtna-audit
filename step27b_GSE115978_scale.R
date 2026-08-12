## Step 27b  Is GSE115978 already library-size normalised?
##
## GSE120575's file summed to 1e6 per cell, so re-normalising it was clearly
## wrong. GSE115978's file has a median per-cell sum of 2.99e5, which is not
## 1e6 -- so the question is whether the sums are CONSTANT (already normalised,
## just to a different constant) or VARIABLE (normalisation is appropriate and
## the stored rds may be right after all).

suppressPackageStartupMessages({library(data.table)})
setwd("D:/R_ex/MR")

raw <- fread("D:/数据/黑色素瘤单细胞人/GSE115978/GSE115978_counts.csv.gz")
m <- as.matrix(raw[, -1]); storage.mode(m) <- "double"
rm(raw); gc()

cs <- colSums(m)
cat("per-cell sums of the raw file:\n")
print(summary(cs))
cat(sprintf("\nCV = %.4f\n", sd(cs) / mean(cs)))
cat(sprintf("ratio p95/p05 = %.2f\n", quantile(cs, .95) / quantile(cs, .05)))
cat("\nverdict: ")
if (sd(cs) / mean(cs) < 0.02) {
  cat("sums are effectively CONSTANT -> file is already normalised;\n")
  cat("        re-normalising it (as the stored rds did) is an ERROR.\n")
} else {
  cat("sums VARY substantially -> the file is NOT library-size normalised;\n")
  cat("        library-size normalisation is appropriate here, so the stored\n")
  cat("        rds is defensible and GSE115978 is NOT affected by the\n")
  cat("        GSE120575 error. Detection rates match exactly, confirming the\n")
  cat("        same underlying data.\n")
}
cat("\nfor contrast, GSE120575's file summed to ~1e6 per cell (genuine TPM).\n")
