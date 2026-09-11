## =====================================================================
## Step 23  Is TPI1 marking a whole glycolytic programme in CD4 T cells,
##          or acting as an isolated signal?
##
## This decides how big the biological story can be:
##   TPI1 tracks the module AND the module stratifies ICB response
##        -> the story is "CD4 glycolytic reprogramming", backed by the
##           whole T-cell immunometabolism literature
##   TPI1 alone stratifies, module does not
##        -> TPI1 is an isolated signal; the story stays small
##
## *** THE TRAP HERE IS THE SAME ONE AS IN CellChat ***
## In scRNA-seq every abundant gene correlates with every other abundant gene
## through shared detection depth. A raw TPI1-vs-glycolysis correlation is
## therefore guaranteed to be positive and means nothing on its own. Controls:
##   1. Seurat AddModuleScore (expression-bin-matched control genes)
##   2. partial correlation on nFeature_RNA
##   3. an OXPHOS module as a specificity contrast -- if TPI1 tracks OXPHOS
##      just as well, this is general metabolic activity, not glycolysis
##
## Run with R 4.4.1.  Output: 28a-28e
## =====================================================================

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2)
})
set.seed(1)
MR <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"
setwd(MR)

## core glycolysis, TPI1 deliberately EXCLUDED (it is the variable under test)
GLYCO <- c("SLC2A1", "SLC2A3", "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM",
           "PFKP", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2",
           "PKM", "LDHA", "PGM1", "PFKFB3", "PFKFB4")
## OXPHOS specificity contrast
OXPHOS <- c("NDUFA1", "NDUFA4", "NDUFB1", "NDUFB2", "NDUFS5", "SDHA", "SDHB",
            "UQCRB", "UQCRC1", "UQCRQ", "COX4I1", "COX5A", "COX5B", "COX6C",
            "COX7C", "COX8A", "ATP5F1A", "ATP5F1B", "ATP5MC2", "ATP5ME")

partial_spearman <- function(x, y, z) {
  d <- data.frame(x, y, z); d <- d[complete.cases(d), ]
  if (nrow(d) < 20) return(c(rho = NA, p = NA))
  rx <- residuals(lm(rank(d$x) ~ rank(d$z)))
  ry <- residuals(lm(rank(d$y) ~ rank(d$z)))
  ct <- cor.test(rx, ry, method = "pearson")
  c(rho = unname(ct$estimate), p = ct$p.value)
}

## =====================================================================
## PART 1  GSE120575  (CD45+ sorted -> immune signal cannot be tumour-derived)
## =====================================================================
obj <- readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds")
cat("GSE120575 CD4 cells:", ncol(obj), "\n")

g_use <- intersect(GLYCO, rownames(obj))
o_use <- intersect(OXPHOS, rownames(obj))
cat("glycolysis genes found:", length(g_use), "/", length(GLYCO),
    " | OXPHOS:", length(o_use), "/", length(OXPHOS), "\n")

obj <- AddModuleScore(obj, features = list(g_use), name = "Glyco", seed = 1)
obj <- AddModuleScore(obj, features = list(o_use), name = "Oxphos", seed = 1)
obj$GlycoScore  <- obj$Glyco1
obj$OxphosScore <- obj$Oxphos1

cell <- data.table(
  cell = colnames(obj), patient = obj$patient, timepoint = obj$timepoint,
  nFeature = obj$nFeature_RNA,
  TPI1 = as.numeric(GetAssayData(obj, layer = "data")["TPI1", ]),
  Glyco = obj$GlycoScore, Oxphos = obj$OxphosScore,
  Prolif = obj$ProlifScore, Exhaust = obj$ExhaustScore)

cat("\n--- cell level (n =", nrow(cell), ") ---\n")
r_raw <- cor(cell$TPI1, cell$Glyco, method = "spearman")
r_par <- partial_spearman(cell$TPI1, cell$Glyco, cell$nFeature)
o_raw <- cor(cell$TPI1, cell$Oxphos, method = "spearman")
o_par <- partial_spearman(cell$TPI1, cell$Oxphos, cell$nFeature)
cat(sprintf("TPI1 ~ glycolysis : raw rho=%.3f | partial(nFeature) rho=%.3f p=%.3g\n",
            r_raw, r_par["rho"], r_par["p"]))
cat(sprintf("TPI1 ~ OXPHOS     : raw rho=%.3f | partial(nFeature) rho=%.3f p=%.3g   <- specificity contrast\n",
            o_raw, o_par["rho"], o_par["p"]))

## ---- patient x timepoint aggregation (the statistical unit used throughout)
samp <- cell[, .(n_CD4 = .N, TPI1 = mean(TPI1), Glyco = mean(Glyco),
                 Oxphos = mean(Oxphos), Prolif = mean(Prolif),
                 Exhaust = mean(Exhaust), nFeature = mean(nFeature)),
             by = .(patient, timepoint)][n_CD4 >= 20]

resp <- fread("19i_patient_tp_with_response.tsv")[, .(patient, timepoint, response)]
samp <- merge(samp, resp, by = c("patient", "timepoint"))
fwrite(samp, "28a_GSE120575_glyco_sample_level.tsv", sep = "\t")
cat("\nsamples with >=20 CD4 cells:", nrow(samp), "\n")

cat("\n--- patient level ---\n")
for (v in c("Glyco", "Oxphos", "Prolif", "Exhaust")) {
  ct <- cor.test(samp$TPI1, samp[[v]], method = "spearman")
  cat(sprintf("TPI1 ~ %-8s rho=%+.3f p=%.3g\n", v, ct$estimate, ct$p.value))
}

## ---- THE KEY TEST: does the glycolysis module stratify ICB response? ----
cat("\n=== responder vs non-responder ===\n")
res_rows <- list()
for (tp in c("Pre", "Post")) {
  d <- samp[timepoint == tp]
  for (v in c("TPI1", "Glyco", "Oxphos", "Prolif", "Exhaust")) {
    R  <- d[response == "Responder", get(v)]
    NR <- d[response == "Non-responder", get(v)]
    if (length(R) < 3 || length(NR) < 3) next
    w <- wilcox.test(R, NR)
    res_rows[[length(res_rows) + 1]] <- data.table(
      timepoint = tp, variable = v, n_R = length(R), n_NR = length(NR),
      median_R = median(R), median_NR = median(NR),
      diff = median(R) - median(NR), p = w$p.value)
  }
}
rr <- rbindlist(res_rows)
rr[, FDR := p.adjust(p, "BH"), by = timepoint]
fwrite(rr, "28b_GSE120575_glyco_response.tsv", sep = "\t")
print(rr[order(timepoint, p)], digits = 3)

## ---- is TPI1's response signal independent of the module? --------------
cat("\n=== is TPI1 independent of the glycolysis module? ===\n")
for (tp in c("Pre", "Post")) {
  d <- samp[timepoint == tp]
  d[, y := as.integer(response == "Non-responder")]
  m1 <- glm(y ~ scale(TPI1), data = d, family = binomial)
  m2 <- glm(y ~ scale(TPI1) + scale(Glyco), data = d, family = binomial)
  cat(sprintf("%-5s TPI1 alone      : beta=%+.2f p=%.3g\n", tp,
              coef(summary(m1))[2, 1], coef(summary(m1))[2, 4]))
  cat(sprintf("%-5s TPI1 | Glyco    : beta=%+.2f p=%.3g   (Glyco beta=%+.2f p=%.3g)\n",
              tp, coef(summary(m2))[2, 1], coef(summary(m2))[2, 4],
              coef(summary(m2))[3, 1], coef(summary(m2))[3, 4]))
}
cat("NOTE: n is 19 (Pre) / 27 (Post). These models are underpowered --\n")
cat("      read them as descriptive, not as a formal mediation test.\n")

## =====================================================================
## PART 2  GSE115978 replication of the TPI1 ~ glycolysis link
##         (no response annotation here -- correlation only)
## =====================================================================
cat("\n\n=== GSE115978 replication ===\n")
o2 <- readRDS("GSE115978_TCD4.rds")
cat("T.CD4 cells:", ncol(o2), "\n")
g2 <- intersect(GLYCO, rownames(o2)); x2 <- intersect(OXPHOS, rownames(o2))
o2 <- AddModuleScore(o2, features = list(g2), name = "Glyco", seed = 1)
o2 <- AddModuleScore(o2, features = list(x2), name = "Oxphos", seed = 1)

c2 <- data.table(
  nFeature = o2$nFeature_RNA,
  TPI1 = as.numeric(GetAssayData(o2, layer = "data")["TPI1", ]),
  Glyco = o2$Glyco1, Oxphos = o2$Oxphos1)
p2 <- partial_spearman(c2$TPI1, c2$Glyco, c2$nFeature)
q2 <- partial_spearman(c2$TPI1, c2$Oxphos, c2$nFeature)
cat(sprintf("TPI1 ~ glycolysis : raw rho=%.3f | partial rho=%.3f p=%.3g\n",
            cor(c2$TPI1, c2$Glyco, method = "spearman"), p2["rho"], p2["p"]))
cat(sprintf("TPI1 ~ OXPHOS     : raw rho=%.3f | partial rho=%.3f p=%.3g\n",
            cor(c2$TPI1, c2$Oxphos, method = "spearman"), q2["rho"], q2["p"]))

fwrite(data.table(
  dataset = c(rep("GSE120575", 2), rep("GSE115978", 2)),
  module = rep(c("glycolysis", "OXPHOS"), 2),
  rho_raw = c(r_raw, o_raw, cor(c2$TPI1, c2$Glyco, method = "spearman"),
              cor(c2$TPI1, c2$Oxphos, method = "spearman")),
  rho_partial_nFeature = c(r_par["rho"], o_par["rho"], p2["rho"], q2["rho"]),
  p_partial = c(r_par["p"], o_par["p"], p2["p"], q2["p"])),
  "28c_TPI1_module_correlation_both_cohorts.tsv", sep = "\t")

cat("\nwritten: 28a-28c\n")
