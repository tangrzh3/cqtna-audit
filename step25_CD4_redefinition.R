## =====================================================================
## Step 25  Is a re-defined CD4 population worth the churn?
##
## THE QUESTION
## The fixed Step 23 CD4 population is impure: CD4 detected 32.5%, CD8A 52.9%,
## CD8-double-negative 44.6%. Step 24c already showed the response results
## survive lineage control, but the strict CD8-negative filter loses power at
## the Pre timepoint (TPI1 Pre p=0.130, n=8 vs 8; was 0.0101 at n=9 vs 10).
##
## SO THE ONLY REAL GAIN from re-defining CD4 is: recover cells that the strict
## single-gene filter throws away because of dropout, and with them recover
## power at Pre. This script measures that gain before anyone re-runs anything.
##
## WHY NOT RE-CLUSTER
## Clustering already produced this population once. At Smart-seq2 dropout the
## CD4/CD8 boundary is exactly where unsupervised clustering is weakest, and a
## second clustering pass gives no guarantee of doing better. Instead we assign
## lineage with MULTI-GENE scores, which are robust to dropout in any single
## gene, and compare three candidate definitions head to head.
##
## Definitions compared
##   D0  fixed Step 23 CD4                      (current, 3,878 cells)
##   D1  D0 minus CD8A/CD8B-positive cells      (strict, used in Step 24c)
##   D2  score-based: CD4score > CD8score among T cells, whole cohort
##   D3  D2 intersected with D0                 (conservative)
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat)})
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

obj  <- readRDS("24_checkpoint_GSE120575_annotated_obj.rds")   # all 16,291 cells
d    <- GetAssayData(obj, assay = "RNA", layer = "data")
fixed_cd4 <- colnames(readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds"))
resp <- fread("19i_patient_tp_with_response.tsv")[, .(patient, timepoint, response)]

cat("full object:", nrow(obj), "x", ncol(obj), "\n")
cat("lineage annotation:\n"); print(table(obj$lineage))
cat("\nfixed Step 23 CD4:", length(fixed_cd4), "cells\n")

GLYCO <- c("SLC2A3", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
           "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3")
## multi-gene lineage panels -- robust to dropout in any single gene
CD4_PANEL <- c("CD4", "IL7R", "CD40LG", "MAL", "LTB", "TRAT1", "ANXA1", "CCR7", "AQP3")
CD8_PANEL <- c("CD8A", "CD8B", "GZMK", "NKG7", "CCL5", "GZMB", "PRF1", "KLRD1", "GNLY")
TCELL     <- c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD2")

zsc <- function(genes, cells = colnames(obj)) {
  g <- intersect(genes, rownames(d))
  x <- t(as.matrix(d[g, cells, drop = FALSE]))
  z <- scale(x); z[!is.finite(z)] <- 0
  rowMeans(z)
}

cell <- data.table(
  cell = colnames(obj), patient = obj$patient, timepoint = obj$timepoint,
  lineage = obj$lineage, nFeature = obj$nFeature_RNA,
  TPI1 = as.numeric(d["TPI1", ]), Glyco = zsc(GLYCO),
  CD4s = zsc(CD4_PANEL), CD8s = zsc(CD8_PANEL), Ts = zsc(TCELL),
  CD8A = as.numeric(d["CD8A", ]), CD8B = as.numeric(d["CD8B", ]),
  CD4g = as.numeric(d["CD4", ]))
cell[, in_fixed := cell %in% fixed_cd4]

## ---------------------------------------------------------------- definitions
is_T <- cell$Ts > 0                       # above-average T-cell score
D0 <- cell$in_fixed
D1 <- D0 & cell$CD8A == 0 & cell$CD8B == 0
D2 <- is_T & (cell$CD4s > cell$CD8s)
D3 <- D2 & D0

defs <- list(D0_fixed = D0, D1_strict_CD8neg = D1, D2_score_cohort = D2, D3_score_x_fixed = D3)

cat("\n=== candidate CD4 definitions ===\n")
summ <- rbindlist(lapply(names(defs), function(nm) {
  k <- defs[[nm]]
  data.table(definition = nm, n_cells = sum(k),
             pct_CD4_detected = 100 * mean(cell$CD4g[k] > 0),
             pct_CD8A_detected = 100 * mean(cell$CD8A[k] > 0),
             pct_CD8B_detected = 100 * mean(cell$CD8B[k] > 0),
             mean_CD4s_minus_CD8s = mean(cell$CD4s[k] - cell$CD8s[k]),
             pct_of_fixed_retained = 100 * mean(k[D0]))
}))
print(summ, digits = 3)
fwrite(summ, "31a_CD4_definition_comparison.tsv", sep = "\t")

## ---------------------------------------------------------------- response test
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

run_def <- function(k, nm, min_cells = 20) {
  s <- cell[k, .(n_CD4 = .N, TPI1 = mean(TPI1), Glyco = mean(Glyco)),
            by = .(patient, timepoint)][n_CD4 >= min_cells]
  s <- merge(s, resp, by = c("patient", "timepoint"))
  rbindlist(lapply(c("Pre", "Post"), function(tp)
    rbindlist(lapply(c("TPI1", "Glyco"), function(v)
      cbind(definition = nm, timepoint = tp, variable = v,
            n_samples = nrow(s[timepoint == tp]),
            wtest(s[timepoint == tp], v))), fill = TRUE)), fill = TRUE)
}

cat("\n=== ICB response stratification under each definition ===\n")
res <- rbindlist(lapply(names(defs), function(nm) run_def(defs[[nm]], nm)), fill = TRUE)
setorder(res, variable, timepoint, definition)
print(res, digits = 3)
fwrite(res, "31b_response_by_CD4_definition.tsv", sep = "\t")

## ---------------------------------------------------------------- overlap
cat("\n=== overlap of D2 (score-based, whole cohort) with the fixed set ===\n")
print(table(fixed = ifelse(D0, "in_fixed", "outside"),
            score = ifelse(D2, "CD4_by_score", "not_CD4")))
cat("\nlineage of D2 cells that were OUTSIDE the fixed CD4 set:\n")
print(sort(table(cell$lineage[D2 & !D0]), decreasing = TRUE))

cat("\nwritten: 31a, 31b\n")
