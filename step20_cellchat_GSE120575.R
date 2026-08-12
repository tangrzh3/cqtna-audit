## =====================================================================
## Step 20  What do glycolysis-high CD4 T cells talk to?
##          GSE120575 (Sade-Feldman) -- CD45+ sorted, ICB, response-annotated
##
## GROUPING (changed after Step 23):
## The primary split is the GLYCOLYSIS MODULE, not TPI1 alone. Step 23 showed
## the functional phenotype is pathway-wide -- 18/22 (Pre) and 19/22 (Post)
## glycolytic enzymes are higher in non-responders (binomial p=0.0022 /
## 0.00043), and the TPI1-excluded module separates responders post-treatment
## (p=0.0094). A module split is more stable than a single gene, has better
## power, and AddModuleScore's expression-binned control genes blunt the
## depth confound. The TPI1 split is still run, as a sensitivity analysis --
## the two must agree for either to be trusted.
##
## SCOPE: this cohort is CD45+ SORTED. No malignant cells, no fibroblasts, no
## endothelium. Everything here is immune-immune. For the CD4 <-> tumour axis
## use step20b (GSE115978).
##
## DEPTH MATCHING -- still mandatory. Cells scoring high on any expression
## module are, by default, cells with deeper libraries, and CellChat's
## communication probability rises with detection. Unmatched, "high group
## signals more" is guaranteed as an artefact. We 1:1 match on nFeature within
## each sample and print the before/after imbalance.
##
## Run with R 4.4.1.  Outputs: 24*_<tag>_*.tsv/.pdf in D:/R_ex/MR/
## =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat)
  library(CellChat); library(ggplot2); library(patchwork)
})
set.seed(1)
options(stringsAsFactors = FALSE, future.globals.maxSize = 8 * 1024^3)

MR      <- "D:/R_ex/MR"
TPM_GZ  <- "D:/Downloads/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
ANNO_GZ <- "D:/Downloads/GSE120575_patient_ID_single_cells.txt.gz"
setwd(MR)

HI_FRAC <- 0.40
LO_FRAC <- 0.40
MIN_CELLS_PER_GROUP <- 10
MIN_DETECT <- 0.30            # glycolytic gene must be detected in >=30% of CD4 cells

GLYCO <- c("SLC2A1", "SLC2A3", "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM",
           "PFKP", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2",
           "PKM", "LDHA", "PGM1", "PFKFB3", "PFKFB4")   # TPI1 excluded on purpose

## =====================================================================
## 1. load the expression matrix
##
##    Layout quirk: row 1 = cell IDs with a LEADING EMPTY FIELD, row 2 =
##    patient/timepoint, rows 3+ = genes. Gene rows carry ONE extra TRAILING
##    empty field.
##
##    *** DO NOT read the header with fread. *** fread turns the leading empty
##    field into NA, and nzchar(NA) is TRUE, so an "empty field" filter does not
##    remove it. n_cells then becomes 16292 instead of 16291, the value columns
##    are clipped one too far, and EVERY COLUMN LABEL SHIFTS BY ONE POSITION --
##    each cell silently carries its neighbour's ID. readLines + strsplit is
##    used instead, and the stopifnot below is the guard.
## =====================================================================
cat("reading expression matrix (~126 MB gz, takes a few minutes)...\n")
con      <- gzfile(TPM_GZ, "r"); hdr_line <- readLines(con, n = 1); close(con)
cell_id  <- strsplit(hdr_line, "\t", fixed = TRUE)[[1]]
cell_id  <- cell_id[!is.na(cell_id) & nzchar(trimws(cell_id))]
n_cells  <- length(cell_id)
cat("  cells in header:", n_cells, "\n")
stopifnot(n_cells == 16291, !anyNA(cell_id), !any(duplicated(cell_id)))

tpm   <- fread(TPM_GZ, skip = 2, header = FALSE, sep = "\t")
genes <- tpm[[1]]
tpm   <- as.matrix(tpm[, 2:(1 + n_cells), with = FALSE])
rownames(tpm) <- make.unique(as.character(genes))
colnames(tpm) <- cell_id
storage.mode(tpm) <- "double"

## alignment guard: a mis-clipped matrix leaves an all-NA trailing column,
## which shows up as a zero column sum
cs <- colSums(tpm, na.rm = TRUE)
cat("  matrix:", nrow(tpm), "genes x", ncol(tpm), "cells | all-NA columns:",
    sum(colSums(is.na(tpm)) == nrow(tpm)), "\n")
stopifnot(ncol(tpm) == n_cells, !anyNA(colnames(tpm)), min(cs) > 0,
          "TPI1" %in% rownames(tpm))

## scale check: this file is ALREADY log2(TPM + 1), not raw TPM.
## Verified empirically -- colSums(2^x - 1) is ~1e6 per cell (i.e. TPM).
## Re-applying a log here would collapse the dynamic range (4.85 -> 0.57)
## and silently ruin module scores, the hi/lo split and CellChat.
tpm_sum <- median(colSums(2^tpm - 1))
cat("  reverse-log per-cell sum (median):", round(tpm_sum),
    " <- must be ~1e6, confirming log2(TPM+1)\n")
stopifnot(tpm_sum > 5e5, tpm_sum < 2e6)

## =====================================================================
## 2. response / patient annotation (metadata block starts at line 20)
## =====================================================================
anno <- fread(ANNO_GZ, skip = 19, header = TRUE, sep = "\t", fill = TRUE)
setnames(anno, old = names(anno)[1:7],
         new = c("sample_no", "cell", "source", "organism",
                 "patient_tp", "response", "therapy"))
## fill=TRUE turns trailing blank/comment lines into NA rows -- drop them,
## otherwise row.names picks up an NA and AddMetaData errors out
anno <- anno[!is.na(cell) & nzchar(trimws(cell))]
anno[, cell := trimws(cell)]
anno <- unique(anno, by = "cell")
setkey(anno, cell)
cat("  usable annotation rows:", nrow(anno), "\n")

## =====================================================================
## 3. Seurat object.
##    The matrix is ALREADY log2(TPM+1) (verified above), so the data layer
##    is the counts layer unchanged -- no further normalisation.
## =====================================================================
obj <- CreateSeuratObject(counts = as.sparse(tpm), project = "GSE120575",
                          min.cells = 3, min.features = 200)
rm(tpm); gc()

## align on colnames(obj), never on the annotation's own order, so the
## row names of the metadata frame cannot contain NA
md <- anno[J(colnames(obj))]
md[, timepoint := ifelse(grepl("^Pre", patient_tp), "Pre", "Post")]
md[, patient   := sub("^(Pre|Post)_", "", patient_tp)]
cat("  cells with no annotation:", sum(is.na(md$patient_tp)), "/", ncol(obj), "\n")

meta_df <- as.data.frame(md[, .(patient, patient_tp, timepoint, response, therapy)])
rownames(meta_df) <- colnames(obj)
obj <- AddMetaData(obj, meta_df)

keep <- !is.na(obj$patient_tp)
if (any(!keep)) obj <- subset(obj, cells = colnames(obj)[keep])
cat("  final:", ncol(obj), "cells |",
    paste(names(table(obj$response)), table(obj$response), collapse = " | "), "\n")

obj <- SetAssayData(obj, layer = "data",
                    new.data = GetAssayData(obj, layer = "counts"))

## consistency with the object used in Steps 16/19/23: TPI1 detection was 64.9%
d0 <- GetAssayData(obj, layer = "data")
cat("  data max:", round(max(d0), 2),
    "| TPI1 detection:", round(100 * mean(d0["TPI1", ] > 0), 1), "%",
    "(expect ~64.9%) | TPI1 non-zero median:",
    round(median(d0["TPI1", d0["TPI1", ] > 0]), 2), "\n")
rm(d0); gc()

obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 1.0, verbose = FALSE)

## =====================================================================
## 4. lineage annotation by cluster-level marker scoring
## =====================================================================
MARK <- list(
  CD8T      = c("CD8A", "CD8B", "GZMK", "GZMB", "CD3D", "CD3E"),
  CD4Tconv  = c("CD4", "IL7R", "CD40LG", "CD3D", "CD3E", "LTB"),
  Treg      = c("FOXP3", "IL2RA", "IKZF2", "CTLA4", "TNFRSF4"),
  Tprolif   = c("MKI67", "TOP2A", "CCNB1", "CDK1"),
  NK        = c("NKG7", "KLRD1", "FGFBP2", "NCAM1", "KLRF1"),
  B         = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  Plasma    = c("MZB1", "JCHAIN", "IGHG1", "SDC1", "XBP1"),
  MonoMacro = c("LYZ", "CD14", "CD68", "FCGR3A", "AIF1", "TYROBP"),
  DC        = c("CD1C", "FCER1A", "CLEC9A", "LAMP3"),
  pDC       = c("LILRA4", "IRF7", "GZMB", "CLEC4C")
)
dat <- GetAssayData(obj, layer = "data")
sc  <- sapply(MARK, function(g) {
  g <- intersect(g, rownames(obj))
  colMeans(scale(t(as.matrix(dat[g, , drop = FALSE]))))
})
cl       <- Idents(obj)
cl_score <- t(apply(sc, 2, function(v) tapply(v, cl, mean)))
lineage  <- rownames(cl_score)[apply(cl_score, 2, which.max)]
names(lineage) <- colnames(cl_score)
obj$lineage <- lineage[as.character(cl)]

cat("\n--- lineage assignment ---\n"); print(table(obj$lineage))
fwrite(data.table(cluster = colnames(cl_score), lineage = lineage, t(cl_score)),
       "24a_GSE120575_cluster_lineage_scores.tsv", sep = "\t")

## =====================================================================
## 5. CD4 subset -> glycolysis module score -> two candidate splits
##    Module score is computed WITHIN the CD4 subset so that AddModuleScore's
##    control-gene bins come from the same cell population (matches Step 23c).
## =====================================================================
cd4_cells <- colnames(obj)[obj$lineage == "CD4Tconv"]
cat("\nCD4Tconv cells:", length(cd4_cells), "\n")
cd4obj <- subset(obj, cells = cd4_cells)

d4  <- GetAssayData(cd4obj, layer = "data")
det <- rowMeans(d4[intersect(GLYCO, rownames(cd4obj)), , drop = FALSE] > 0)
glyco_use <- names(det)[det >= MIN_DETECT]
cat("glycolytic genes detected in >=", MIN_DETECT * 100, "% of CD4 cells: ",
    length(glyco_use), "\n  ", paste(glyco_use, collapse = ", "), "\n", sep = "")
cat("dropped (too sparse): ",
    paste(sprintf("%s(%.1f%%)", names(det)[det < MIN_DETECT],
                  100 * det[det < MIN_DETECT]), collapse = ", "), "\n")

cd4obj <- AddModuleScore(cd4obj, features = list(glyco_use), name = "Glyco", seed = 1)

md <- data.table(cell = colnames(cd4obj),
                 sample = cd4obj$patient_tp,
                 response = cd4obj$response,
                 nFeature = cd4obj$nFeature_RNA,
                 TPI1 = as.numeric(d4["TPI1", ]),
                 Glyco = cd4obj$Glyco1)
cat("\nTPI1 vs glycolysis module, cell level: rho =",
    round(cor(md$TPI1, md$Glyco, method = "spearman"), 3), "\n")

## greedy 1:1 nearest-neighbour matching on nFeature, within sample
match_on_depth <- function(d) {
  hi <- d[grp == "hi"][order(nFeature)]; lo <- d[grp == "lo"][order(nFeature)]
  if (nrow(hi) == 0 || nrow(lo) == 0) return(d[0])
  keep_hi <- keep_lo <- integer(0); used <- rep(FALSE, nrow(lo))
  for (i in seq_len(nrow(hi))) {
    cand <- which(!used); if (!length(cand)) break
    j <- cand[which.min(abs(lo$nFeature[cand] - hi$nFeature[i]))]
    used[j] <- TRUE; keep_hi <- c(keep_hi, i); keep_lo <- c(keep_lo, j)
  }
  rbind(hi[keep_hi], lo[keep_lo])
}

make_split <- function(var) {
  d <- copy(md)
  d[, v := get(var)]
  d[, r := frank(v, ties.method = "random") / .N, by = sample]
  d[, grp := fifelse(r > 1 - HI_FRAC, "hi", fifelse(r <= LO_FRAC, "lo", NA_character_))]
  d <- d[!is.na(grp)]
  cat(sprintf("\n--- split on %s: depth imbalance BEFORE matching ---\n", var))
  print(d[, .(n = .N, median_nFeature = median(nFeature),
              median_value = round(median(v), 3)), by = grp])
  cat("  Wilcoxon nFeature p =",
      signif(wilcox.test(nFeature ~ grp, data = d)$p.value, 3), "\n")
  dm <- rbindlist(lapply(split(d, d$sample), match_on_depth))
  cat(sprintf("--- split on %s: AFTER matching ---\n", var))
  print(dm[, .(n = .N, median_nFeature = median(nFeature),
               median_value = round(median(v), 3)), by = grp])
  pv <- wilcox.test(nFeature ~ grp, data = dm)$p.value
  cat("  Wilcoxon nFeature p =", signif(pv, 3),
      ifelse(pv < 0.05, "  *** STILL IMBALANCED -- do not interpret ***",
             "  <- balanced, ok to interpret"), "\n")
  dm
}

split_glyco <- make_split("Glyco")
split_tpi1  <- make_split("TPI1")
fwrite(split_glyco, "24b_CD4_split_Glyco.tsv", sep = "\t")
fwrite(split_tpi1,  "24b_CD4_split_TPI1.tsv",  sep = "\t")

## how much do the two splits agree?
ov <- merge(split_glyco[, .(cell, glyco_grp = grp)],
            split_tpi1[, .(cell, tpi1_grp = grp)], by = "cell")
cat("\n--- agreement between the two splits ---\n")
print(table(Glyco = ov$glyco_grp, TPI1 = ov$tpi1_grp))
cat("concordance:", round(100 * mean(ov$glyco_grp == ov$tpi1_grp), 1), "%\n")

## =====================================================================
## 6. CellChat, run once per split
## =====================================================================
run_cellchat <- function(spl, tag) {
  cat("\n==================== CellChat:", tag, "====================\n")
  labels <- setNames(obj$lineage, colnames(obj))
  labels[spl[grp == "hi", cell]] <- paste0("CD4_", tag, "hi")
  labels[spl[grp == "lo", cell]] <- paste0("CD4_", tag, "lo")
  ## CD4 cells not in either extreme are dropped to keep the contrast clean
  labels <- labels[setdiff(colnames(obj), setdiff(cd4_cells, spl$cell))]
  tab <- table(labels); print(tab)
  labels <- labels[labels %in% names(tab)[tab >= MIN_CELLS_PER_GROUP]]
  keep <- names(labels)

  cc <- createCellChat(object = GetAssayData(obj, layer = "data")[, keep],
                       meta = data.frame(labels = factor(labels), row.names = keep),
                       group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = MIN_CELLS_PER_GROUP)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  saveRDS(cc, sprintf("24_cellchat_%s.rds", tag))

  net <- as.data.table(subsetCommunication(cc))
  fwrite(net, sprintf("24c_%s_LR_all.tsv", tag), sep = "\t")

  hi <- paste0("CD4_", tag, "hi"); lo <- paste0("CD4_", tag, "lo")
  ## descriptive contrast: CellChat has no formal test between two cell groups
  ## inside one run. pval is a label-permutation p per link. Rank, do not test.
  cmp <- function(role) {
    other <- setdiff(c("source", "target"), role)
    d <- net[get(role) %in% c(hi, lo)]
    if (!nrow(d)) return(data.table())
    d[, partner := get(other)]
    w <- dcast(d, partner + pathway_name + interaction_name + annotation +
                 ligand + receptor ~ get(role),
               value.var = "prob", fun.aggregate = sum, fill = 0)
    if (!hi %in% names(w)) w[, (hi) := 0]
    if (!lo %in% names(w)) w[, (lo) := 0]
    w[, delta := get(hi) - get(lo)]
    w[order(-delta)]
  }
  out_hi <- cmp("source"); in_hi <- cmp("target")
  fwrite(out_hi, sprintf("24d_%s_outgoing.tsv", tag), sep = "\t")
  fwrite(in_hi,  sprintf("24e_%s_incoming.tsv",  tag), sep = "\t")

  cat("\n=== top OUTGOING enriched in", hi, "===\n")
  print(head(out_hi, 20))
  cat("\n=== top INCOMING enriched in", hi, "===\n")
  print(head(in_hi, 20))

  pdf(sprintf("24_%s_bubble.pdf", tag), width = 11, height = 9)
  print(netVisual_bubble(cc, sources.use = c(hi, lo), remove.isolate = FALSE))
  dev.off()
  ggsave(sprintf("24_%s_signalingRole.pdf", tag),
         netAnalysis_signalingRole_scatter(cc), width = 7, height = 6)

  ## LR pairs to carry into the spatial validation (step21)
  tp <- rbind(head(out_hi[get(hi) > 0], 25)[, .(interaction_name, pathway_name,
                                                annotation, ligand, receptor,
                                                partner, direction = "CD4hi_sends", delta)],
              head(in_hi[get(hi) > 0], 25)[, .(interaction_name, pathway_name,
                                               annotation, ligand, receptor,
                                               partner, direction = "CD4hi_receives", delta)])
  fwrite(unique(tp), sprintf("24f_%s_top_LR_for_spatial.tsv", tag), sep = "\t")
  net
}

net_glyco <- run_cellchat(split_glyco, "Glyco")
net_tpi1  <- run_cellchat(split_tpi1,  "TPI1")

## =====================================================================
## 7. do the two splits agree on which pathways come up?
##    They must, or neither result is trustworthy.
## =====================================================================
top_pw <- function(tag) {
  f <- sprintf("24d_%s_outgoing.tsv", tag)
  if (!file.exists(f)) return(data.table())
  d <- fread(f)
  d[, .(delta = sum(delta)), by = pathway_name][order(-delta)][1:20]
}
g <- top_pw("Glyco"); t <- top_pw("TPI1")
if (nrow(g) && nrow(t)) {
  m <- merge(g, t, by = "pathway_name", suffixes = c("_Glyco", "_TPI1"), all = TRUE)
  fwrite(m, "24g_split_concordance_pathways.tsv", sep = "\t")
  cat("\n=== pathway concordance between the two splits ===\n")
  print(m[order(-delta_Glyco)])
  cat("\noverlap in top-20 outgoing pathways:",
      length(intersect(g$pathway_name, t$pathway_name)), "/ 20\n")
}

cat("\nwritten: 24a-24g, 2 rds, 4 pdf\n")
cat("next: python step21_spatial_LR_coloc.py 24f_Glyco_top_LR_for_spatial.tsv\n")
