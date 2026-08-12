## =====================================================================
## Step 23b  Module dilution, or is TPI1 genuinely special?
##
## Step 23 found: TPI1 separates responders at BOTH timepoints, but the
## 21-gene glycolysis module separates them only Post, not Pre (p=0.095).
## Two readings:
##   (a) the module is diluted -- several glycolytic genes carry signal,
##       averaging 21 genes buries it  -> the "glycolytic reprogramming"
##       story survives
##   (b) TPI1 alone carries the signal                 -> it does not
##
## Test every glycolytic gene individually against ICB response, exactly as
## TPI1 was tested. Where does TPI1 rank among its own pathway?
## Output: 28d, 28e
## =====================================================================

suppressPackageStartupMessages({library(Seurat); library(data.table)})
set.seed(1)
setwd("D:/R_ex/MR")

GLYCO <- c("TPI1", "SLC2A1", "SLC2A3", "HK1", "HK2", "HK3", "GPI", "PFKL",
           "PFKM", "PFKP", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1",
           "ENO2", "PKM", "LDHA", "PGM1", "PFKFB3", "PFKFB4")

obj  <- readRDS("GSE120575_CD4_clusters_1_5_6_12_15.rds")
dat  <- GetAssayData(obj, layer = "data")
use  <- intersect(GLYCO, rownames(obj))
resp <- fread("19i_patient_tp_with_response.tsv")[, .(patient, timepoint, response)]

cell <- data.table(patient = obj$patient, timepoint = obj$timepoint,
                   detect = colSums(dat[use, ] > 0))
for (g in use) cell[[g]] <- as.numeric(dat[g, ])

samp <- cell[, c(.(n_CD4 = .N), lapply(.SD, mean)), by = .(patient, timepoint),
             .SDcols = use][n_CD4 >= 20]
samp <- merge(samp, resp, by = c("patient", "timepoint"))

## detection rate per gene, so low-expressed genes are not read as "no signal"
det <- data.table(gene = use,
                  pct_detected = round(100 * rowMeans(dat[use, ] > 0), 1))

out <- list()
for (tp in c("Pre", "Post")) {
  d <- samp[timepoint == tp]
  for (g in use) {
    R  <- d[response == "Responder", get(g)]
    NR <- d[response == "Non-responder", get(g)]
    if (length(R) < 3 || length(NR) < 3) next
    w <- suppressWarnings(wilcox.test(R, NR))
    out[[length(out) + 1]] <- data.table(
      timepoint = tp, gene = g, n_R = length(R), n_NR = length(NR),
      median_R = median(R), median_NR = median(NR),
      diff = median(R) - median(NR), p = w$p.value)
  }
}
res <- rbindlist(out)
res <- merge(res, det, by = "gene")
res[, FDR := p.adjust(p, "BH"), by = timepoint]
setorder(res, timepoint, p)
fwrite(res, "28d_glycolysis_gene_by_gene_response.tsv", sep = "\t")

for (tp in c("Pre", "Post")) {
  cat("\n=== ", tp, " (responder vs non-responder, each glycolytic gene) ===\n", sep = "")
  print(res[timepoint == tp,
            .(gene, pct_detected, median_R = round(median_R, 2),
              median_NR = round(median_NR, 2), diff = round(diff, 2),
              p = signif(p, 3), FDR = signif(FDR, 3))], nrows = 30)
  n_nom <- res[timepoint == tp & p < 0.05, .N]
  n_fdr <- res[timepoint == tp & FDR < 0.05, .N]
  cat(sprintf("  nominal p<0.05: %d / %d genes | FDR<0.05: %d\n",
              n_nom, res[timepoint == tp, .N], n_fdr))
  cat(sprintf("  TPI1 rank by p: %d of %d\n",
              which(res[timepoint == tp, gene] == "TPI1"), res[timepoint == tp, .N]))
}

## direction consistency: do the others at least point the same way as TPI1?
cat("\n=== direction (negative diff = higher in non-responders, like TPI1) ===\n")
print(res[, .(n_genes = .N, n_higher_in_NR = sum(diff < 0)), by = timepoint])

fwrite(det, "28e_glycolysis_detection_rate.tsv", sep = "\t")
cat("\nwritten: 28d, 28e\n")
