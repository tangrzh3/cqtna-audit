## =====================================================================
## Step 25d  Definitive re-test using the GEO annotation as deposited.
##
## Step 25c uncovered the decisive fact: P1 and P5 each have TWO post-treatment
## biopsies carrying DIFFERENT response labels (Responder vs Non-responder) --
## these are separate treatment courses, not replicates. The 19i table used for
## every published number collapses them into one sample with a single label,
## so responder and non-responder cells were pooled under one annotation.
##
## Here nothing is re-parsed: sample = patient_tp exactly as deposited,
## response = the per-cell response field of the object, expression = the
## correctly scaled log2(TPM+1) object.
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat)})
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

new <- readRDS("24_checkpoint_GSE120575_annotated_obj.rds")
d   <- GetAssayData(new, assay = "RNA", layer = "data")
fixed_cd4 <- colnames(readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds"))

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
  sample = new$patient_tp[match(fixed_cd4, colnames(new))],
  timepoint = new$timepoint[match(fixed_cd4, colnames(new))],
  response = new$response[match(fixed_cd4, colnames(new))],
  CD8A = as.numeric(d["CD8A", fixed_cd4]),
  CD8B = as.numeric(d["CD8B", fixed_cd4]),
  TPI1 = as.numeric(d["TPI1", fixed_cd4]),
  Glyco = zsc(GLYCO, fixed_cd4))
cell[, cd8_free := CD8A == 0 & CD8B == 0]

## confirm the mixed-label samples
cat("=== samples whose response label differs from a same-patient sibling ===\n")
sib <- cell[, .(n = .N, response = response[1]), by = .(sample, timepoint)]
sib[, base := sub("_2$", "", sub("^(Pre|Post)_", "", sample))]
print(sib[base %in% sib[, .N, by = .(base, timepoint)][N > 1, base]][order(base, sample)])

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

run <- function(dt, label, min_cells = 20) {
  s <- dt[, .(n_CD4 = .N, TPI1 = mean(TPI1), Glyco = mean(Glyco),
              response = response[1]), by = .(sample, timepoint)][n_CD4 >= min_cells]
  rbindlist(lapply(c("Pre", "Post"), function(tp)
    rbindlist(lapply(c("TPI1", "Glyco"), function(v)
      cbind(cells = label, timepoint = tp, variable = v,
            n_samples = nrow(s[timepoint == tp]),
            wtest(s[timepoint == tp], v))), fill = TRUE)), fill = TRUE)
}

cat("\n=== definitive: GEO sample IDs, GEO response, correct expression scale ===\n")
res <- rbind(run(cell, "all fixed CD4"),
             run(cell[cd8_free == TRUE], "CD8-free subset"), fill = TRUE)
setorder(res, variable, timepoint, cells)
print(res, digits = 3)
fwrite(res, "31d_definitive_response.tsv", sep = "\t")

cat("\n--- comparison table ---\n")
cmp <- data.table(
  analysis = c("published (double-normalised, samples collapsed)",
               "correct scale, samples collapsed",
               "correct scale, GEO sample IDs (definitive)"),
  TPI1_Pre  = c(0.0101, 0.00762, res[variable=="TPI1" & timepoint=="Pre" & cells=="all fixed CD4", p]),
  TPI1_Post = c(0.0022, 0.01145, res[variable=="TPI1" & timepoint=="Post" & cells=="all fixed CD4", p]),
  Glyco_Pre = c(0.133,  0.15640, res[variable=="Glyco" & timepoint=="Pre" & cells=="all fixed CD4", p]),
  Glyco_Post= c(0.0094, 0.01340, res[variable=="Glyco" & timepoint=="Post" & cells=="all fixed CD4", p]))
print(cmp, digits = 3)
fwrite(cmp, "31e_result_provenance_comparison.tsv", sep = "\t")
cat("\nwritten: 31d, 31e\n")
