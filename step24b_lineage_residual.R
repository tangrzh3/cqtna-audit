## =====================================================================
## Step 24b  Higher-powered version of the purity control.
##
## Step 24 filtered to CD8A==0 & CD8B==0 pairs and retained only 17.8% of
## pairs (208 cells), so a null there is confounded with loss of power.
## Here the CD8 axis is removed by residualisation instead, keeping all
## 1,168 matched cells.
##
## CD8ness = mean z(CD8A), z(CD8B) at cell level.
## Only the LINEAGE-FREE modules are residualised -- residualising
## CytotoxicScore on CD8A/CD8B would be circular, since those genes are
## in the score.
##   EffectorClean = GZMB, NKG7, PRF1, GNLY, IFNG   (no CD8A/CD8B/CCL5/GZMK)
##   HelperClean   = CD40LG, FOXP3, IL2RA, CTLA4    (no CD4)
##
## Also reported: a matching-based version, where hi/lo are matched on
## CD8ness bin in addition to nFeature bin.
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat)})
setwd("D:/R_ex/MR")

obj <- readRDS("24_checkpoint_Post_Step23CD4_GlycoScore.rds")
d   <- GetAssayData(obj, assay = "RNA", layer = "data")
sg  <- as.data.table(readRDS("24_checkpoint_Glyco_Post_matched.rds"))

EFF_CLEAN    <- c("GZMB", "NKG7", "PRF1", "GNLY", "IFNG")
HELPER_CLEAN <- c("CD40LG", "FOXP3", "IL2RA", "CTLA4")
CYTO_ORIG    <- c("CD8A", "CD8B", "CCL5", "GZMK", "GZMB", "NKG7")
HELPER_ORIG  <- c("CD4", "CD40LG", "FOXP3", "IL2RA", "CTLA4")

zsc <- function(genes) {
  g <- intersect(genes, rownames(d))
  x <- t(as.matrix(d[g, colnames(obj), drop = FALSE]))
  z <- scale(x); z[!is.finite(z)] <- 0
  rowMeans(z)
}
sc <- data.table(cell = colnames(obj),
                 CytotoxicScore = zsc(CYTO_ORIG), HelperRegScore = zsc(HELPER_ORIG),
                 EffectorClean = zsc(EFF_CLEAN), HelperClean = zsc(HELPER_CLEAN),
                 CD8ness = zsc(c("CD8A", "CD8B")))

md <- merge(sg[, .(cell, sample, response, grp, nFeature, Glyco, match_pair)],
            sc, by = "cell")
md <- copy(md)

paired_test <- function(x, col) {
  w <- dcast(x[, .(v = mean(get(col))), by = .(sample, grp)], sample ~ grp,
             value.var = "v")
  w <- w[is.finite(hi) & is.finite(lo)]
  wt <- suppressWarnings(wilcox.test(w$hi, w$lo, paired = TRUE, exact = FALSE))
  data.table(n_samples = nrow(w), mean_delta = mean(w$hi - w$lo),
             n_hi_higher = sum(w$hi > w$lo), n_hi_lower = sum(w$hi < w$lo),
             p = wt$p.value)
}

## =====================================================================
## 1. residualise the lineage-free modules on CD8ness
## =====================================================================
cat("=== CD8ness vs Glyco score (all matched cells) ===\n")
cat(sprintf("  rho = %+.3f\n", cor(md$Glyco, md$CD8ness, method = "spearman")))
cat(sprintf("  CD8ness mean: hi = %+.3f, lo = %+.3f\n",
            md[grp == "hi", mean(CD8ness)], md[grp == "lo", mean(CD8ness)]))

for (col in c("EffectorClean", "HelperClean")) {
  md[[paste0(col, "_resid")]] <- residuals(lm(md[[col]] ~ md$CD8ness))
}

cat("\n=== module tests: raw vs CD8ness-residualised (all 1,168 cells) ===\n")
res <- rbindlist(lapply(
  c("CytotoxicScore", "HelperRegScore", "EffectorClean", "HelperClean",
    "EffectorClean_resid", "HelperClean_resid"),
  function(c0) cbind(module = c0, paired_test(md, c0))), fill = TRUE)
res[, FDR := p.adjust(p, "BH")]
print(res, digits = 4)
fwrite(res, "29f_lineage_residualised_modules.tsv", sep = "\t")

## =====================================================================
## 2. matching on CD8ness bin as well as nFeature bin
## =====================================================================
smd <- function(x, col) {
  hi <- x[grp == "hi", get(col)]; lo <- x[grp == "lo", get(col)]
  (mean(hi) - mean(lo)) / sqrt((var(hi) + var(lo)) / 2)
}
## CD8ness has a large point mass at the double-negative value, so plain
## quantile breaks are not unique. Bin 0 = CD8A and CD8B both undetected;
## bins 1-3 = tertiles of the remaining cells.
md[, cd8_free := CD8ness <= min(CD8ness) + 1e-9]
qs <- unique(quantile(md[cd8_free == FALSE, CD8ness], seq(0, 1, 1/3)))
md[, cd8_bin := 0L]
md[cd8_free == FALSE,
   cd8_bin := as.integer(cut(CD8ness, breaks = qs, include.lowest = TRUE))]
cat("\ncells per CD8ness bin (0 = double-negative):\n"); print(table(md$cd8_bin, md$grp))
md[, depth_bin := floor(nFeature / 200)]

pieces <- list(); k <- 1L
for (s in unique(md$sample)) for (cb in unique(md$cd8_bin)) for (db in unique(md$depth_bin)) {
  z <- md[sample == s & cd8_bin == cb & depth_bin == db]
  hi <- z[grp == "hi"][order(CD8ness, nFeature, cell)]
  lo <- z[grp == "lo"][order(CD8ness, nFeature, cell)]
  n <- min(nrow(hi), nrow(lo)); if (n == 0) next
  pieces[[k]] <- rbind(hi[seq_len(n)], lo[seq_len(n)]); k <- k + 1L
}
mm <- rbindlist(pieces)
cat("\n=== matched on CD8ness bin + nFeature bin ===\n")
cat("cells retained:", nrow(mm), "of", nrow(md),
    sprintf(" (%.1f%%)", 100 * nrow(mm) / nrow(md)),
    "| samples:", uniqueN(mm$sample), "\n")
cat(sprintf("  residual SMD  CD8ness = %+.4f | nFeature = %+.4f\n",
            smd(mm, "CD8ness"), smd(mm, "nFeature")))

res2 <- rbindlist(lapply(
  c("CytotoxicScore", "HelperRegScore", "EffectorClean", "HelperClean"),
  function(c0) cbind(module = c0, paired_test(mm, c0))), fill = TRUE)
res2[, FDR := p.adjust(p, "BH")]
print(res2, digits = 4)
fwrite(res2, "29g_CD8ness_matched_modules.tsv", sep = "\t")

## per-marker in the CD8ness-matched set
mk <- rbindlist(lapply(unique(c(HELPER_ORIG, CYTO_ORIG, "PRF1", "GNLY", "IFNG")),
  function(g) {
    if (!g %in% rownames(d)) return(NULL)
    x <- copy(mm); x[, v := as.numeric(d[g, cell])]
    cbind(gene = g, paired_test(x, "v"))
  }), fill = TRUE)
mk[, FDR := p.adjust(p, "BH")]
cat("\n=== per-marker, CD8ness-matched ===\n")
print(mk[order(p)], digits = 4)
fwrite(mk, "29h_CD8ness_matched_per_marker.tsv", sep = "\t")

cat("\nwritten: 29f-29h\n")
