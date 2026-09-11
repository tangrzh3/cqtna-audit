## =====================================================================
## Step 25b  Reconcile the two GSE120575 objects.
##
## Same 3,878 fixed CD4 cells, same test, different answers:
##   Step 24c (old rds GSE120575_CD4_clusters_1_5_6_12_15.rds)
##       TPI1 Post p=0.0022, n=9R/18NR, median_R 3.29 / median_NR 4.60
##   Step 25  (new 24_checkpoint_GSE120575_annotated_obj.rds)
##       TPI1 Post p=0.075,  n=8R/17NR, median_R 4.86 / median_NR 6.90
##
## Two candidate causes, both checked here:
##   (a) different expression scale between the objects
##   (b) different patient/timepoint/response metadata, changing which
##       samples pass the >=20 cell threshold
## Nothing downstream should be trusted until this is resolved.
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat)})
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

old <- readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds")
new <- readRDS("24_checkpoint_GSE120575_annotated_obj.rds")
do  <- GetAssayData(old, assay = "RNA", layer = "data")
dn  <- GetAssayData(new, assay = "RNA", layer = "data")

cat("old:", nrow(old), "x", ncol(old), " | new:", nrow(new), "x", ncol(new), "\n")
common <- intersect(colnames(old), colnames(new))
cat("cell overlap:", length(common), "of", ncol(old), "\n\n")

## ---------------------------------------------------------------- (a) scale
cat("=== (a) expression scale ===\n")
cat(sprintf("old  max=%.2f  TPI1 nonzero median=%.3f  detection=%.1f%%\n",
            max(do), median(do["TPI1", do["TPI1", ] > 0]), 100*mean(do["TPI1", ] > 0)))
cat(sprintf("new  max=%.2f  TPI1 nonzero median=%.3f  detection=%.1f%%\n",
            max(dn), median(dn["TPI1", dn["TPI1", ] > 0]), 100*mean(dn["TPI1", ] > 0)))

sub <- sample(common, min(200, length(common)))
cat(sprintf("\nreverse-log per-cell sum (median):\n  old assuming log2(x+1):   %.3g\n",
            median(colSums(2^as.matrix(do[, sub]) - 1))))
cat(sprintf("  new assuming log2(x+1):   %.3g   <- ~1e6 means log2(TPM+1)\n",
            median(colSums(2^as.matrix(dn[, sub]) - 1))))
cat(sprintf("  old assuming log2(x/10+1): %.3g   <- ~1e6 means log2(TPM/10+1)\n",
            median(colSums((2^as.matrix(do[, sub]) - 1) * 10))))

## direct relation on shared cells/genes
g <- intersect(rownames(do), rownames(dn))
g <- sample(g, min(2000, length(g)))
vo <- as.numeric(do[g, sub]); vn <- as.numeric(dn[g, sub])
ok <- vo > 0 & vn > 0
cat(sprintf("\nnonzero-value agreement: rho=%.4f | median(new-old)=%.3f | median(2^new-1)/(2^old-1)=%.2f\n",
            cor(vo[ok], vn[ok], method = "spearman"),
            median(vn[ok] - vo[ok]),
            median((2^vn[ok] - 1) / (2^vo[ok] - 1))))
cat("  a ratio near 10 means old = log2(TPM/10+1) and new = log2(TPM+1)\n")

## ---------------------------------------------------------------- (b) metadata
cat("\n=== (b) metadata on the shared cells ===\n")
mo <- data.table(cell = common,
                 patient_old = old$patient[match(common, colnames(old))],
                 tp_old = old$timepoint[match(common, colnames(old))])
mn <- data.table(cell = common,
                 patient_new = new$patient[match(common, colnames(new))],
                 tp_new = new$timepoint[match(common, colnames(new))],
                 resp_new = new$response[match(common, colnames(new))])
m <- merge(mo, mn, by = "cell")
cat("patient disagreements:", sum(m$patient_old != m$patient_new), "\n")
cat("timepoint disagreements:", sum(m$tp_old != m$tp_new), "\n")
if (sum(m$patient_old != m$patient_new)) print(head(m[patient_old != patient_new], 10))

resp <- fread("19i_patient_tp_with_response.tsv")[, .(patient, timepoint, response)]
for (lab in c("old", "new")) {
  x <- if (lab == "old") m[, .(cell, patient = patient_old, timepoint = tp_old)]
       else m[, .(cell, patient = patient_new, timepoint = tp_new)]
  s <- x[, .(n = .N), by = .(patient, timepoint)][n >= 20]
  s <- merge(s, resp, by = c("patient", "timepoint"))
  cat(sprintf("\n%s metadata -> samples with >=20 cells: %d  (Pre %d, Post %d)\n",
              lab, nrow(s), sum(s$timepoint == "Pre"), sum(s$timepoint == "Post")))
  print(s[, .N, by = .(timepoint, response)][order(timepoint, response)])
}

## which samples differ
so <- m[, .(n = .N), by = .(patient = patient_old, timepoint = tp_old)][n >= 20]
sn <- m[, .(n = .N), by = .(patient = patient_new, timepoint = tp_new)][n >= 20]
setnames(so, c("patient", "timepoint", "n_old")); setnames(sn, c("patient", "timepoint", "n_new"))
cmp <- merge(so, sn, by = c("patient", "timepoint"), all = TRUE)
cat("\nsamples present in only one version, or with differing cell counts:\n")
print(cmp[is.na(n_old) | is.na(n_new) | n_old != n_new])

cat("\n=== response annotation source check ===\n")
r19 <- fread("19i_patient_tp_with_response.tsv")
cat("19i rows:", nrow(r19), " | Pre:", sum(r19$timepoint == "Pre"),
    " Post:", sum(r19$timepoint == "Post"), "\n")
print(r19[, .N, by = .(timepoint, response)][order(timepoint, response)])
rn <- unique(m[, .(patient_new, tp_new, resp_new)])
chk <- merge(rn, r19, by.x = c("patient_new", "tp_new"),
             by.y = c("patient", "timepoint"), all.x = TRUE)
cat("\nobject response vs 19i response disagreements:",
    sum(!is.na(chk$response) & chk$resp_new != chk$response), "\n")
cat("patient_tp in object but not in 19i:", sum(is.na(chk$response)), "\n")
if (sum(is.na(chk$response))) print(chk[is.na(response)])
