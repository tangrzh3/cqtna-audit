## =====================================================================
## Step 24  Is the Glyco-hi CD4 phenotype a biological state, or a CD8
##          purity gradient?
##
## WHY THIS IS NEEDED
## 24e_Glyco_hi_lo_CD4_CD8_marker_check.tsv shows a large composition
## difference between the matched groups:
##     CD4  detected 40.1% (hi) vs 19.5% (lo)
##     CD8A detected 47.6% (hi) vs 63.9% (lo)
## Four of the six CytotoxicScore genes (CD8A, CD8B, CCL5, GZMK) are CD8
## lineage/lineage-biased, and CD4 itself is the strongest HelperRegScore
## marker (27/28 samples). So both modules could be reading one axis:
## how much CD8-like contamination sits in each group.
##
## The existing control (24e_Glyco_marker_direction_across_original_clusters)
## shows direction holds within all 5 original CD4 clusters. That excludes
## CLUSTER COMPOSITION as the driver, but not cell-level purity -- contaminating
## CD8-like cells can sit inside every cluster.
##
## DESIGN
##  A. quantify the purity gradient, sample-aware
##  B. keep only matched PAIRS in which BOTH cells are CD8A==0 & CD8B==0,
##     which preserves the 1:1 depth matching, and redo the paired tests
##  C. rescore with lineage-free modules (no CD8A/CD8B/CCL5/GZMK, no CD4)
##  D. ask whether the Glyco score itself tracks CD8-ness
##
## Reproduces the published numbers first (-0.3181 / +0.2347) so that any
## change afterwards is attributable to the control, not to rescoring.
## =====================================================================

suppressPackageStartupMessages({library(data.table); library(Seurat)})
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

OBJ_RDS   <- "24_checkpoint_Post_Step23CD4_GlycoScore.rds"   # Seurat, 2427 Post CD4
GLYCO_RDS <- "24_checkpoint_Glyco_Post_matched.rds"          # 584/584 matched
TPI1_RDS  <- "24_checkpoint_TPI1_Post_matched.rds"

CYTO_ORIG   <- c("CD8A", "CD8B", "CCL5", "GZMK", "GZMB", "NKG7")
HELPER_ORIG <- c("CD4", "CD40LG", "FOXP3", "IL2RA", "CTLA4")
## lineage-free: effector molecules only, no CD8 lineage genes, no CCL5/GZMK
EFF_CLEAN    <- c("GZMB", "NKG7", "PRF1", "GNLY", "IFNG")
HELPER_CLEAN <- c("CD40LG", "FOXP3", "IL2RA", "CTLA4")       # CD4 removed

MIN_PAIRS_PER_SAMPLE <- 5

obj <- readRDS(OBJ_RDS)
d   <- GetAssayData(obj, assay = "RNA", layer = "data")
sg  <- as.data.table(readRDS(GLYCO_RDS))
cat("Post CD4 cells:", ncol(obj), "| matched Glyco cells:", nrow(sg),
    "| samples:", uniqueN(sg$sample), "\n\n")

## ---- module scores, exactly as in the original analysis -----------------
## z-score across ALL post CD4 cells, then rowMeans of each gene subset
mk_scores <- function(sets) {
  genes <- unique(unlist(sets))
  genes <- intersect(genes, rownames(d))
  x <- t(as.matrix(d[genes, colnames(obj), drop = FALSE]))
  z <- scale(x); z[!is.finite(z)] <- 0
  out <- data.table(cell = rownames(z))
  for (nm in names(sets)) {
    g <- intersect(sets[[nm]], colnames(z))
    out[[nm]] <- rowMeans(z[, g, drop = FALSE])
  }
  out
}
sc <- mk_scores(list(CytotoxicScore = CYTO_ORIG, HelperRegScore = HELPER_ORIG,
                     EffectorClean = EFF_CLEAN, HelperClean = HELPER_CLEAN))

md <- merge(sg[, .(cell, sample, response, grp, nFeature, Glyco, TPI1, match_pair)],
            sc, by = "cell")
for (g in c("CD8A", "CD8B", "CD4")) md[[g]] <- as.numeric(d[g, md$cell])

smd <- function(x) {
  hi <- x[grp == "hi", nFeature]; lo <- x[grp == "lo", nFeature]
  (mean(hi) - mean(lo)) / sqrt((var(hi) + var(lo)) / 2)
}
paired_test <- function(x, col) {
  w <- dcast(x[, .(v = mean(get(col))), by = .(sample, grp)], sample ~ grp,
             value.var = "v")
  w <- w[is.finite(hi) & is.finite(lo)]
  if (nrow(w) < 5) return(data.table(n_samples = nrow(w)))
  wt <- suppressWarnings(wilcox.test(w$hi, w$lo, paired = TRUE, exact = FALSE))
  data.table(n_samples = nrow(w), mean_delta = mean(w$hi - w$lo),
             n_hi_higher = sum(w$hi > w$lo), n_hi_lower = sum(w$hi < w$lo),
             p = wt$p.value)
}

## =====================================================================
## 0. reproduce the published numbers
## =====================================================================
cat("=== 0. reproduction check (expect -0.3181 / +0.2347) ===\n")
rep0 <- rbindlist(lapply(c("CytotoxicScore", "HelperRegScore"),
                         function(c0) cbind(module = c0, paired_test(md, c0))))
print(rep0, digits = 4)

## =====================================================================
## A. the purity gradient, sample-aware
## =====================================================================
cat("\n=== A. composition of the matched groups ===\n")
comp <- md[, .(n = .N,
               pct_CD8A_pos = 100 * mean(CD8A > 0),
               pct_CD8B_pos = 100 * mean(CD8B > 0),
               pct_CD4_pos  = 100 * mean(CD4 > 0),
               pct_CD8_free = 100 * mean(CD8A == 0 & CD8B == 0)), by = grp]
print(comp, digits = 3)

pur <- rbindlist(lapply(c("CD8A", "CD8B", "CD4"), function(g) {
  x <- copy(md); x[, v := as.numeric(get(g) > 0)]
  cbind(gene = g, paired_test(x, "v"))
}))
cat("\nsample-aware detection-fraction difference (hi - lo):\n")
print(pur, digits = 4)

## =====================================================================
## D. does the Glyco score itself track CD8-ness?
## =====================================================================
cat("\n=== D. Glyco score vs lineage ===\n")
all_post <- data.table(cell = colnames(obj), Glyco = obj$Glyco1,
                       CD8A = as.numeric(d["CD8A", ]), CD8B = as.numeric(d["CD8B", ]),
                       CD4 = as.numeric(d["CD4", ]))
cat(sprintf("across all 2427 post CD4 cells:  Glyco~CD8A rho=%+.3f | Glyco~CD8B rho=%+.3f | Glyco~CD4 rho=%+.3f\n",
            cor(all_post$Glyco, all_post$CD8A, method = "spearman"),
            cor(all_post$Glyco, all_post$CD8B, method = "spearman"),
            cor(all_post$Glyco, all_post$CD4,  method = "spearman")))

## =====================================================================
## B + C. CD8-free matched pairs, original and lineage-free modules
## =====================================================================
md[, cd8_free := CD8A == 0 & CD8B == 0]
ok_pairs <- md[, .(both_free = all(cd8_free), n = .N), by = match_pair][
  both_free == TRUE & n == 2, match_pair]
pure <- md[match_pair %in% ok_pairs]

cat("\n=== B. CD8-free matched pairs ===\n")
cat("pairs retained:", length(ok_pairs), "of", uniqueN(md$match_pair),
    sprintf(" (%.1f%%)\n", 100 * length(ok_pairs) / uniqueN(md$match_pair)))
cat("cells:", nrow(pure), " | hi/lo:", sum(pure$grp == "hi"), "/", sum(pure$grp == "lo"),
    " | nFeature SMD:", round(smd(pure), 4), "\n")
keep_s <- pure[, .N, by = .(sample, grp)][, .(ok = all(N >= MIN_PAIRS_PER_SAMPLE) &
                                                .N == 2), by = sample][ok == TRUE, sample]
pure_s <- pure[sample %in% keep_s]
cat("samples with >=", MIN_PAIRS_PER_SAMPLE, "cells per group:", length(keep_s),
    "of", uniqueN(md$sample), "\n")

cat("\n--- paired tests ---\n")
res <- rbindlist(lapply(c("CytotoxicScore", "HelperRegScore",
                          "EffectorClean", "HelperClean"), function(c0) {
  rbind(cbind(module = c0, set = "all matched pairs",     paired_test(md, c0)),
        cbind(module = c0, set = "CD8-free pairs",        paired_test(pure, c0)),
        cbind(module = c0, set = "CD8-free, >=5/sample",  paired_test(pure_s, c0)),
        fill = TRUE)
}), fill = TRUE)
res[, FDR := p.adjust(p, "BH"), by = set]
print(res, digits = 4)

fwrite(res,  "29a_purity_control_module_tests.tsv", sep = "\t")
fwrite(comp, "29b_matched_group_composition.tsv", sep = "\t")
fwrite(pur,  "29c_sample_aware_lineage_detection.tsv", sep = "\t")
fwrite(pure[, .(cell, sample, grp, match_pair, nFeature, Glyco)],
       "29d_CD8free_matched_pairs.tsv", sep = "\t")

## =====================================================================
## per-marker breakdown inside CD8-free pairs
## =====================================================================
cat("\n=== per-marker, CD8-free pairs ===\n")
mk <- rbindlist(lapply(c(HELPER_ORIG, CYTO_ORIG, "PRF1", "GNLY", "IFNG"), function(g) {
  if (!g %in% rownames(d)) return(NULL)
  x <- copy(pure); x[, v := as.numeric(d[g, cell])]
  cbind(gene = g, paired_test(x, "v"))
}), fill = TRUE)
mk[, FDR := p.adjust(p, "BH")]
print(mk[order(p)], digits = 4)
fwrite(mk, "29e_purity_control_per_marker.tsv", sep = "\t")

cat("\nwritten: 29a-29e\n")
