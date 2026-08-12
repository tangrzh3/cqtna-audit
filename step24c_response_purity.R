## =====================================================================
## Step 24c  Does the ICB-response signal survive the same lineage control?
##
## Step 24/24b showed the Glyco-hi CD4 "helper/regulatory-leaning, less
## cytotoxic" phenotype is a CD4/CD8 purity artefact: after matching hi/lo on
## CD8ness as well as depth, every module goes null.
##
## The same fixed CD4 population underpins the results the manuscript actually
## depends on:
##   Step 16/19h  TPI1 responder vs non-responder  (Pre FDR 0.034, Post 0.015)
##   Step 23      glycolysis module vs response     (Post p=0.0094, Pre p=0.133)
## Those are BETWEEN-PATIENT comparisons, so the within-patient hi/lo artefact
## does not automatically apply -- but if responders and non-responders differ
## in the CD4/CD8 composition of these clusters, the same confound operates.
##
## Tests
##  1. do R and NR differ in CD4/CD8 composition of the fixed CD4 population?
##  2. recompute TPI1 and the glycolysis module using CD8-free cells only
##  3. recompute using per-cell values residualised on CD8ness
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat)})
setwd("D:/R_ex/MR")

obj  <- readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds")   # 3,878 fixed CD4
d    <- GetAssayData(obj, assay = "RNA", layer = "data")
resp <- fread("19i_patient_tp_with_response.tsv")[, .(patient, timepoint, response)]
cat("fixed CD4 cells:", ncol(obj), "\n")

GLYCO <- c("SLC2A3", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
           "PGAM1", "ENO1", "PKM", "LDHA", "PFKFB3")     # Step 23c module, no TPI1

zsc <- function(genes) {
  g <- intersect(genes, rownames(d))
  x <- t(as.matrix(d[g, , drop = FALSE]))
  z <- scale(x); z[!is.finite(z)] <- 0
  rowMeans(z)
}

cell <- data.table(
  cell = colnames(obj), patient = obj$patient, timepoint = obj$timepoint,
  nFeature = obj$nFeature_RNA,
  TPI1 = as.numeric(d["TPI1", ]),
  Glyco = zsc(GLYCO), CD8ness = zsc(c("CD8A", "CD8B")),
  CD8A = as.numeric(d["CD8A", ]), CD8B = as.numeric(d["CD8B", ]),
  CD4 = as.numeric(d["CD4", ]))
cell[, cd8_free := CD8A == 0 & CD8B == 0]

cat("purity of the fixed CD4 population:\n")
cat(sprintf("  CD4 detected %.1f%% | CD8A detected %.1f%% | CD8B detected %.1f%% | double-negative %.1f%%\n",
            100*mean(cell$CD4 > 0), 100*mean(cell$CD8A > 0),
            100*mean(cell$CD8B > 0), 100*mean(cell$cd8_free)))

wtest <- function(d, v) {
  R  <- d[response == "Responder", get(v)]
  NR <- d[response == "Non-responder", get(v)]
  if (length(R) < 3 || length(NR) < 3) return(data.table(n_R = length(R), n_NR = length(NR)))
  w <- suppressWarnings(wilcox.test(R, NR))
  data.table(n_R = length(R), n_NR = length(NR),
             median_R = median(R), median_NR = median(NR),
             diff = median(R) - median(NR), p = w$p.value)
}

## =====================================================================
## 1. composition difference between responders and non-responders
## =====================================================================
samp <- cell[, .(n_CD4 = .N, pct_CD8A_pos = mean(CD8A > 0),
                 pct_cd8_free = mean(cd8_free), mean_CD8ness = mean(CD8ness),
                 pct_CD4_pos = mean(CD4 > 0), nFeature = mean(nFeature)),
             by = .(patient, timepoint)][n_CD4 >= 20]
samp <- merge(samp, resp, by = c("patient", "timepoint"))

cat("\n=== 1. R vs NR composition of the fixed CD4 population ===\n")
comp <- rbindlist(lapply(c("Pre", "Post"), function(tp)
  rbindlist(lapply(c("pct_CD8A_pos", "pct_cd8_free", "mean_CD8ness", "pct_CD4_pos"),
                   function(v) cbind(timepoint = tp, variable = v,
                                     wtest(samp[timepoint == tp], v))),
            fill = TRUE)), fill = TRUE)
print(comp, digits = 3)
fwrite(comp, "30a_response_composition.tsv", sep = "\t")

## =====================================================================
## 2/3. TPI1 and glycolysis vs response, with and without lineage control
## =====================================================================
cell[, TPI1_resid  := residuals(lm(TPI1  ~ CD8ness))]
cell[, Glyco_resid := residuals(lm(Glyco ~ CD8ness))]

agg <- function(dt, label) {
  s <- dt[, .(n_CD4 = .N, TPI1 = mean(TPI1), Glyco = mean(Glyco),
              TPI1_resid = mean(TPI1_resid), Glyco_resid = mean(Glyco_resid)),
          by = .(patient, timepoint)][n_CD4 >= 20]
  s <- merge(s, resp, by = c("patient", "timepoint"))
  rbindlist(lapply(c("Pre", "Post"), function(tp)
    rbindlist(lapply(c("TPI1", "Glyco", "TPI1_resid", "Glyco_resid"),
                     function(v) cbind(set = label, timepoint = tp, variable = v,
                                       wtest(s[timepoint == tp], v))),
              fill = TRUE)), fill = TRUE)
}

cat("\n=== 2/3. TPI1 and glycolysis vs ICB response ===\n")
out <- rbind(agg(cell, "all fixed CD4"),
             agg(cell[cd8_free == TRUE], "CD8-free cells only"), fill = TRUE)
out[, FDR := p.adjust(p, "BH"), by = .(set, timepoint)]
print(out[order(timepoint, set, variable)], digits = 3)
fwrite(out, "30b_response_TPI1_glyco_lineage_control.tsv", sep = "\t")

cat("\nreference, from Step 16/19h and Step 23c (all fixed CD4, no lineage control):\n")
cat("  TPI1  Pre p=0.0101 (FDR 0.034) | Post p=0.0022 (FDR 0.015)\n")
cat("  Glyco Pre p=0.133             | Post p=0.0094\n")
cat("\nwritten: 30a, 30b\n")
