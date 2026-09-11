## =====================================================================
## Step 23c  Dilution confirmed -- rebuild the module from detectable genes
##
## Step 23  : 22-gene glycolysis module separates responders Post (p=0.011)
##            but not Pre (p=0.095)  -> looked like "TPI1 is an isolated signal"
## Step 23b : gene-by-gene, 18/22 (Pre) and 19/22 (Post) glycolytic genes are
##            HIGHER IN NON-RESPONDERS, and 6 / 9 reach nominal p<0.05.
##            The pathway does shift; the module score was diluted by genes
##            detected in <10% of cells (HK3 1.4%, PFKFB4 1.3%, PFKM 4.5%).
##
## Here: (1) formal binomial test of the direction consistency
##       (2) module rebuilt from genes detected in >=30% of CD4 cells,
##           TPI1 excluded, and re-tested against response
##       (3) leave-TPI1-out check -- does the pathway still separate responders
##           when TPI1 itself is removed? (guards against TPI1 driving it)
## Output: 28f, 28g
## =====================================================================

suppressPackageStartupMessages({library(Seurat); library(data.table)})
set.seed(1)
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

MIN_DETECT <- 30   # percent of CD4 cells

gg  <- fread("28d_glycolysis_gene_by_gene_response.tsv")
det <- fread("28e_glycolysis_detection_rate.tsv")

## ---- (1) binomial test of direction consistency ------------------------
cat("=== direction consistency across the glycolytic pathway ===\n")
bt <- gg[, {
  n <- .N; k <- sum(diff < 0)
  .(n_genes = n, n_higher_in_NR = k,
    binom_p = binom.test(k, n, 0.5, alternative = "greater")$p.value)
}, by = timepoint]
print(bt)
cat("\n(diff<0 = higher in non-responders, the TPI1 direction)\n")

## restricted to detectable genes only, TPI1 removed
gg2 <- merge(gg, det, by = "gene", suffixes = c("", ".d"))
gg2 <- gg2[pct_detected >= MIN_DETECT & gene != "TPI1"]
cat(sprintf("\nrestricted to %d genes detected in >=%d%% of cells, TPI1 excluded:\n",
            uniqueN(gg2$gene), MIN_DETECT))
print(gg2[, .(n_genes = .N, n_higher_in_NR = sum(diff < 0),
              n_nominal = sum(p < 0.05),
              binom_p = binom.test(sum(diff < 0), .N, 0.5,
                                   alternative = "greater")$p.value),
          by = timepoint])

## ---- (2) rebuild the module from detectable genes ----------------------
keep <- det[pct_detected >= MIN_DETECT & gene != "TPI1", gene]
cat("\nmodule genes:", paste(keep, collapse = ", "), "\n")

obj  <- readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds")
obj  <- AddModuleScore(obj, features = list(keep), name = "GlycoDet", seed = 1)
resp <- fread("19i_patient_tp_with_response.tsv")[, .(patient, timepoint, response)]

cell <- data.table(patient = obj$patient, timepoint = obj$timepoint,
                   GlycoDet = obj$GlycoDet1,
                   TPI1 = as.numeric(GetAssayData(obj, layer = "data")["TPI1", ]))
samp <- cell[, .(n_CD4 = .N, GlycoDet = mean(GlycoDet), TPI1 = mean(TPI1)),
             by = .(patient, timepoint)][n_CD4 >= 20]
samp <- merge(samp, resp, by = c("patient", "timepoint"))
fwrite(samp, "28f_glyco_detected_module_sample_level.tsv", sep = "\t")

cat("\n=== refined module (TPI1 excluded) vs ICB response ===\n")
out <- list()
for (tp in c("Pre", "Post")) {
  d <- samp[timepoint == tp]
  for (v in c("TPI1", "GlycoDet")) {
    R <- d[response == "Responder", get(v)]; NR <- d[response == "Non-responder", get(v)]
    w <- suppressWarnings(wilcox.test(R, NR))
    out[[length(out) + 1]] <- data.table(
      timepoint = tp, variable = v, n_R = length(R), n_NR = length(NR),
      median_R = median(R), median_NR = median(NR),
      diff = median(R) - median(NR), p = w$p.value)
  }
}
res <- rbindlist(out)
print(res, digits = 3)
fwrite(res, "28g_glyco_detected_module_response.tsv", sep = "\t")

cat("\nKEY: GlycoDet excludes TPI1 entirely. If it separates responders,\n")
cat("     the signal is pathway-wide and not carried by TPI1 alone.\n")

## ---- (3) how much of the module is TPI1? -------------------------------
cat("\n=== TPI1 vs refined module, patient level ===\n")
ct <- cor.test(samp$TPI1, samp$GlycoDet, method = "spearman")
cat(sprintf("rho = %+.3f  p = %.3g\n", ct$estimate, ct$p.value))
cat("\nwritten: 28f, 28g\n")
