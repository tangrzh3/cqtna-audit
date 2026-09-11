## =====================================================================
## Step 25c  Which is driving the TPI1 Post discrepancy: expression scale,
##           or sample definition?
##
## Step 25b established:
##  (a) the OLD rds (GSE120575_CD4_clusters_1_5_6_12_15.rds) is double-
##      normalised. Reverse-log does not return ~1e6 per cell, and its values
##      correlate with the correctly scaled object at rho=0.874, not 1.0 --
##      so it is not a monotonic rescale but a per-cell library-size
##      renormalisation applied on top of log2(TPM+1).
##      The NEW object reverse-logs to 9.94e5 and matches the raw file.
##  (b) 637 cells carry different patient labels. The new object treats
##      P1_2/P3_2/P5_2/P23_2/P28_2 as separate samples and merges
##      P4_T_enriched into P4; the old object does the opposite.
##      Post samples: 27 (9R/18NR) old vs 25 (8R/17NR) new.
##
## TPI1 Post moved from p=0.0022 to p=0.075. This script holds the expression
## scale fixed at the CORRECT one (new object) and varies only the sample
## definition, so the two causes are separated.
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat)})
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

new  <- readRDS("24_checkpoint_GSE120575_annotated_obj.rds")
d    <- GetAssayData(new, assay = "RNA", layer = "data")
fixed_cd4 <- colnames(readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds"))
resp <- fread("19i_patient_tp_with_response.tsv")[, .(patient, timepoint, response)]

GLYCO <- c("SLC2A3", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
           "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3")
zsc <- function(genes, cells) {
  g <- intersect(genes, rownames(d))
  x <- t(as.matrix(d[g, cells, drop = FALSE]))
  z <- scale(x); z[!is.finite(z)] <- 0
  rowMeans(z)
}

cell <- data.table(
  cell = fixed_cd4,
  patient_new = new$patient[match(fixed_cd4, colnames(new))],
  timepoint   = new$timepoint[match(fixed_cd4, colnames(new))],
  response_obj = new$response[match(fixed_cd4, colnames(new))],
  TPI1 = as.numeric(d["TPI1", fixed_cd4]),
  Glyco = zsc(GLYCO, fixed_cd4))

## sample definition A = as parsed by the new object (replicate biopsies split)
cell[, sample_A := paste0(timepoint, "_", patient_new)]
## sample definition B = old parsing: collapse the "_2" replicate suffix, so a
## patient contributes one Post sample; T_enriched libraries stay with the patient
cell[, patient_B := sub("_2$", "", patient_new)]
cell[, sample_B := paste0(timepoint, "_", patient_B)]

cat("cells:", nrow(cell), "\n")
cat("distinct samples  A:", uniqueN(cell$sample_A), " B:", uniqueN(cell$sample_B), "\n\n")

wtest <- function(x, v) {
  R  <- x[response == "Responder", get(v)]
  NR <- x[response == "Non-responder", get(v)]
  if (length(R) < 3 || length(NR) < 3)
    return(data.table(n_R = length(R), n_NR = length(NR), p = NA_real_))
  w <- suppressWarnings(wilcox.test(R, NR))
  data.table(n_R = length(R), n_NR = length(NR),
             median_R = median(R), median_NR = median(NR),
             diff = median(R) - median(NR), p = w$p.value)
}

run <- function(pat_col, label, min_cells = 20) {
  s <- cell[, .(n_CD4 = .N, TPI1 = mean(TPI1), Glyco = mean(Glyco)),
            by = c("timepoint", pat_col)][n_CD4 >= min_cells]
  setnames(s, pat_col, "patient")
  s <- merge(s, resp, by = c("patient", "timepoint"))
  rbindlist(lapply(c("Pre", "Post"), function(tp)
    rbindlist(lapply(c("TPI1", "Glyco"), function(v)
      cbind(sample_def = label, timepoint = tp, variable = v,
            wtest(s[timepoint == tp], v))), fill = TRUE)), fill = TRUE)
}

cat("=== correct expression scale, sample definition varied ===\n")
res <- rbind(run("patient_new", "A: replicate biopsies separate (new)"),
             run("patient_B",   "B: replicates collapsed (old/19i)"), fill = TRUE)
setorder(res, variable, timepoint, sample_def)
print(res, digits = 3)
fwrite(res, "31c_sample_definition_effect.tsv", sep = "\t")

cat("\nfor reference, published numbers came from the double-normalised object:\n")
cat("  TPI1 Pre p=0.0101 (FDR 0.034) | TPI1 Post p=0.0022 (FDR 0.015)\n")
cat("  Glyco Pre p=0.133             | Glyco Post p=0.0094\n")

## ---- cell-level check, free of both aggregation and normalisation choices ---
cat("\n=== cell-level sanity check (rank test, no patient aggregation) ===\n")
resp_u <- unique(resp, by = c("patient", "timepoint"))
cl <- merge(cell, resp_u, by.x = c("patient_B", "timepoint"),
            by.y = c("patient", "timepoint"))
cell_res <- rbindlist(lapply(c("Pre", "Post"), function(tp)
  rbindlist(lapply(c("TPI1", "Glyco"), function(v)
    cbind(timepoint = tp, variable = v, wtest(cl[timepoint == tp], v))), fill = TRUE)),
  fill = TRUE)
print(cell_res, digits = 3)
cat("NOTE: cell-level p values are anticonservative (cells are not independent);\n")
cat("      shown only to confirm the direction, never as the reported statistic.\n")

## ---- which samples enter/leave, and the two response disagreements ----------
cat("\n=== samples by definition ===\n")
for (pc in c("patient_new", "patient_B")) {
  s <- cell[, .(n = .N), by = c("timepoint", pc)][n >= 20]
  setnames(s, pc, "patient")
  s <- merge(s, resp, by = c("patient", "timepoint"), all.x = TRUE)
  cat("\n", pc, ": ", nrow(s), " samples, ", sum(is.na(s$response)),
      " with no response annotation in 19i\n", sep = "")
  if (any(is.na(s$response))) print(s[is.na(response)])
}

cat("\n=== object response vs 19i response ===\n")
chk <- unique(cell[, .(patient_B, timepoint, response_obj)])
chk <- merge(chk, resp, by.x = c("patient_B", "timepoint"),
             by.y = c("patient", "timepoint"), all.x = TRUE)
print(chk[!is.na(response) & response_obj != response])
cat("\nwritten: 31c\n")
