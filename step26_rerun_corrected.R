## =====================================================================
## Step 26  Re-run of the GSE120575 single-cell layer on corrected inputs.
##
## Supersedes: Step 16 (19a-19i), Step 23 (28a-28g), Step 24c (30b).
##
## Two errors fixed (see FINDINGS Step 25):
##   1. expression -- the old rds was double-normalised (log2(TPM+1) fed to
##      NormalizeData). Use 24_checkpoint_GSE120575_annotated_obj.rds, whose
##      reverse-log per-cell sum is 9.94e5 (i.e. genuine TPM).
##   2. sample/response -- 19i collapses Post_P1 (Responder) with Post_P1_2
##      (Non-responder), and Post_P5 with Post_P5_2 in the opposite direction.
##      Use the GEO patient_tp as the sample ID and the GEO response field.
##
## The fixed Step 23 CD4 cell IDs are reused unchanged; only expression and
## annotation are corrected, so this is a re-computation, not a redefinition.
##
## NOTE ON SCORE DEFINITIONS: the proliferation set is the one documented in
## FINDINGS Step 16 (MKI67/TOP2A/CCNB1/CDK1/PCNA). The exhaustion set was not
## recorded, so a standard panel is specified here and used consistently --
## exhaustion numbers are therefore NOT strictly comparable to the old ones.
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat)})
set.seed(1)
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

obj_all   <- readRDS("24_checkpoint_GSE120575_annotated_obj.rds")
fixed_cd4 <- colnames(readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds"))
cd4 <- subset(obj_all, cells = fixed_cd4)
rm(obj_all); gc()
d <- GetAssayData(cd4, assay = "RNA", layer = "data")
cat("CD4 cells:", ncol(cd4), " | reverse-log per-cell sum (median):",
    signif(median(colSums(2^as.matrix(d[, sample(ncol(d), 100)]) - 1)), 3), "\n")

CAND   <- c("ZFYVE19", "SMC2", "SPSB2", "TPI1", "HLA-C", "KIAA0040")
PROLIF <- c("MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA")
EXHAUST<- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX")
GLYCO_ALL <- c("SLC2A1", "SLC2A3", "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM",
               "PFKP", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1",
               "ENO2", "PKM", "LDHA", "PGM1", "PFKFB3", "PFKFB4", "TPI1")

zsc <- function(genes) {
  g <- intersect(genes, rownames(d))
  x <- t(as.matrix(d[g, , drop = FALSE]))
  z <- scale(x); z[!is.finite(z)] <- 0
  rowMeans(z)
}

cell <- data.table(
  cell = colnames(cd4),
  sample = cd4$patient_tp, timepoint = cd4$timepoint, response = cd4$response,
  CD8A = as.numeric(d["CD8A", ]), CD8B = as.numeric(d["CD8B", ]),
  ProlifScore = zsc(PROLIF), ExhaustScore = zsc(EXHAUST))
for (g in CAND) cell[[g]] <- if (g %in% rownames(d)) as.numeric(d[g, ]) else NA_real_
cell[, cd8_free := CD8A == 0 & CD8B == 0]

## glycolysis genes detectable in >=30% of CD4 cells, TPI1 excluded
det <- rowMeans(d[intersect(GLYCO_ALL, rownames(d)), , drop = FALSE] > 0)
glyco_use <- setdiff(names(det)[det >= 0.30], "TPI1")
cat("glycolysis genes >=30% detection (TPI1 excluded):", length(glyco_use), "\n  ",
    paste(glyco_use, collapse = ", "), "\n")
cell[, GlycoScore := zsc(glyco_use)]
fwrite(data.table(gene = names(det), detection = round(det, 3),
                  used_in_module = names(det) %in% glyco_use),
       "32a_glyco_gene_detection_corrected.tsv", sep = "\t")

## sample-level table
samp <- cell[, c(.(n_CD4 = .N, response = response[1]),
                 lapply(.SD, mean)),
             by = .(sample, timepoint),
             .SDcols = c(CAND, "ProlifScore", "ExhaustScore", "GlycoScore")][n_CD4 >= 20]
cat("\nsamples with >=20 CD4 cells:", nrow(samp), "\n")
print(samp[, .N, by = .(timepoint, response)][order(timepoint, response)])
fwrite(samp, "32b_sample_level_corrected.tsv", sep = "\t")

wtest <- function(x, v) {
  R  <- x[response == "Responder", get(v)]
  NR <- x[response == "Non-responder", get(v)]
  if (length(R) < 3 || length(NR) < 3) return(data.table(n_R=length(R), n_NR=length(NR)))
  w <- suppressWarnings(wilcox.test(R, NR))
  data.table(n_R = length(R), n_NR = length(NR),
             median_R = median(R), median_NR = median(NR),
             diff = median(R) - median(NR), p = w$p.value)
}

## =====================================================================
## 1. response stratification -- FULL FAMILY, BH within timepoint
##    (this is what supersedes 19h; FDR is only meaningful family-wise)
## =====================================================================
cat("\n=== 1. responder vs non-responder, full family ===\n")
fam <- c(CAND, "ProlifScore", "ExhaustScore", "GlycoScore")
r1 <- rbindlist(lapply(c("Pre", "Post"), function(tp)
  rbindlist(lapply(fam, function(v)
    cbind(timepoint = tp, variable = v, wtest(samp[timepoint == tp], v))),
    fill = TRUE)), fill = TRUE)
r1[, FDR := p.adjust(p, "BH"), by = timepoint]
setorder(r1, timepoint, p)
print(r1, digits = 3)
fwrite(r1, "32c_response_family_corrected.tsv", sep = "\t")

## CD8-free sensitivity, same family
sf <- cell[cd8_free == TRUE][
  , c(.(n_CD4 = .N, response = response[1]), lapply(.SD, mean)),
  by = .(sample, timepoint), .SDcols = fam][n_CD4 >= 20]
r1s <- rbindlist(lapply(c("Pre", "Post"), function(tp)
  rbindlist(lapply(fam, function(v)
    cbind(timepoint = tp, variable = v, wtest(sf[timepoint == tp], v))),
    fill = TRUE)), fill = TRUE)
r1s[, FDR := p.adjust(p, "BH"), by = timepoint]
cat("\n--- CD8-free sensitivity ---\n")
print(r1s[order(timepoint, p)], digits = 3)
fwrite(r1s, "32d_response_family_CD8free.tsv", sep = "\t")

## =====================================================================
## 2. glycolysis pathway, gene by gene  (supersedes 28d)
## =====================================================================
cat("\n=== 2. glycolysis pathway gene by gene ===\n")
gg_genes <- intersect(GLYCO_ALL, rownames(d))
gcell <- cell[, .(cell, sample, timepoint, response)]
for (g in gg_genes) gcell[[g]] <- as.numeric(d[g, ])
gsamp <- gcell[, c(.(n_CD4 = .N, response = response[1]), lapply(.SD, mean)),
               by = .(sample, timepoint), .SDcols = gg_genes][n_CD4 >= 20]
r2 <- rbindlist(lapply(c("Pre", "Post"), function(tp)
  rbindlist(lapply(gg_genes, function(g)
    cbind(timepoint = tp, gene = g, wtest(gsamp[timepoint == tp], g))),
    fill = TRUE)), fill = TRUE)
r2 <- merge(r2, data.table(gene = names(det), detection = round(100*det, 1)), by = "gene")
r2[, FDR := p.adjust(p, "BH"), by = timepoint]
setorder(r2, timepoint, p)
print(r2, digits = 3)
fwrite(r2, "32e_glyco_gene_by_gene_corrected.tsv", sep = "\t")

cat("\n--- direction consistency across the pathway ---\n")
bt <- r2[!is.na(p), {
  n <- .N; k <- sum(diff < 0)
  .(n_genes = n, n_higher_in_NR = k, n_nominal = sum(p < 0.05),
    binom_p = binom.test(k, n, 0.5, alternative = "greater")$p.value)
}, by = timepoint]
print(bt, digits = 3)
cat("TPI1 rank by p:\n")
print(r2[, .(rank = which(gene == "TPI1"), n_genes = .N), by = timepoint])
fwrite(bt, "32f_glyco_direction_consistency.tsv", sep = "\t")

## =====================================================================
## 3. patient-level correlations  (supersedes 19b)
## =====================================================================
cat("\n=== 3. patient-level correlations with proliferation / exhaustion ===\n")
r3 <- rbindlist(lapply(c("Pre", "Post"), function(tp) {
  s <- samp[timepoint == tp]
  rbindlist(lapply(CAND, function(g) rbindlist(lapply(c("ProlifScore", "ExhaustScore"),
    function(v) {
      ct <- suppressWarnings(cor.test(s[[g]], s[[v]], method = "spearman"))
      data.table(timepoint = tp, gene = g, score = v, n = nrow(s),
                 rho = unname(ct$estimate), p = ct$p.value)
    }))))
}))
r3[, FDR := p.adjust(p, "BH"), by = .(timepoint, score)]
print(r3[order(timepoint, score, p)], digits = 3)
fwrite(r3, "32g_patient_correlations_corrected.tsv", sep = "\t")

## =====================================================================
## 4. side-by-side with the previously published numbers
## =====================================================================
cat("\n=== 4. published vs corrected ===\n")
pub <- data.table(
  timepoint = c("Pre","Pre","Pre","Post","Post","Post"),
  variable = c("TPI1","ProlifScore","ExhaustScore","TPI1","SMC2","KIAA0040"),
  published_p = c(0.0128, 0.0101, 0.0128, 0.0019, 0.029, 0.017),
  published_FDR = c(0.034, 0.034, 0.034, 0.015, 0.077, 0.067))
cmp <- merge(pub, r1[, .(timepoint, variable, corrected_p = p, corrected_FDR = FDR)],
             by = c("timepoint", "variable"), all.x = TRUE)
print(cmp, digits = 3)
fwrite(cmp, "32h_published_vs_corrected.tsv", sep = "\t")

cat("\nwritten: 32a-32h\n")
