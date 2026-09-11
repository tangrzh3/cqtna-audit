## =====================================================================
## Step 20b  The CD4 <-> tumour axis.  GSE115978 (Jerby-Arnon)
##
## WHY THIS SCRIPT EXISTS: GSE120575 is CD45+ sorted, so it cannot show any
## interaction between CD4 T cells and malignant cells. GSE115978 contains
## malignant cells, CAFs and endothelium alongside immune cells -- but it has
## NO ICB response annotation. The two cohorts are complementary:
##   GSE120575 -> immune-immune network, response-annotated
##   GSE115978 -> CD4 <-> tumour axis, no response information
## Do not mix conclusions across them.
##
## Same depth-matching logic as step20 -- see the long comment there for why.
## Extra caution here: SPSB2/SMC2/ZFYVE19 detection in this dataset is 5-11%
## (see FINDINGS Step 18); TPI1 is 60.6%, so only TPI1 is usable for splitting.
## =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat)
  library(CellChat); library(ggplot2)
})
set.seed(1)
options(stringsAsFactors = FALSE, future.globals.maxSize = 8 * 1024^3)

MR      <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR"
CNT_GZ  <- "D:/数据/黑色素瘤单细胞人/GSE115978/GSE115978_counts.csv.gz"
ANN_GZ  <- "D:/数据/黑色素瘤单细胞人/GSE115978/GSE115978_cell.annotations.csv.gz"
setwd(MR)

HI_FRAC <- 0.40
LO_FRAC <- 0.40
MIN_CELLS_PER_GROUP <- 10

## ---------------------------------------------------------------- load
cnt <- fread(CNT_GZ)
genes <- cnt[[1]]
cnt <- as.matrix(cnt[, -1])
rownames(cnt) <- make.unique(as.character(genes))
storage.mode(cnt) <- "double"

ann <- fread(ANN_GZ)
setnames(ann, 1, "cell")
ann <- ann[cell %in% colnames(cnt)]
cnt <- cnt[, ann$cell]
cat("cells:", ncol(cnt), " genes:", nrow(cnt), "\n")
cat("--- cell types provided by the authors ---\n"); print(table(ann$cell.types))
stopifnot("TPI1" %in% rownames(cnt))

obj <- CreateSeuratObject(counts = as.sparse(cnt), project = "GSE115978",
                          min.cells = 3, min.features = 200)
ann <- ann[cell %in% colnames(obj)]
obj <- AddMetaData(obj, as.data.frame(ann[, .(samples, cell.types, treatment.group)],
                                      row.names = ann$cell))
rm(cnt); gc()

## ---------------------------------------------------------------------
## SCALE CHECK -- do not assume. GSE120575 turned out to be already
## log2(TPM+1); re-logging would collapse the dynamic range and silently
## wreck the module score, the hi/lo split and CellChat. This file is a
## different deposit, so decide from the data.
##   reverse-log per-cell sum ~1e6  -> already log2(TPM+1), use as is
##   raw per-cell sum ~1e6          -> raw TPM, apply log2(x+1)
## ---------------------------------------------------------------------
cnt0 <- GetAssayData(obj, layer = "counts")
raw_sum <- median(Matrix::colSums(cnt0))
rev_sum <- median(Matrix::colSums(2^as.matrix(cnt0[1:min(5000, nrow(cnt0)), ]) - 1)) *
           (nrow(cnt0) / min(5000, nrow(cnt0)))
cat(sprintf("max=%.2f | raw per-cell sum median=%.3g | reverse-log (approx)=%.3g\n",
            max(cnt0), raw_sum, rev_sum))

if (max(cnt0) <= 25 && rev_sum > 3e5) {
  cat("-> already log-transformed; using counts layer unchanged\n")
  obj <- SetAssayData(obj, layer = "data", new.data = cnt0)
} else {
  cat("-> raw TPM detected; applying log2(x+1)\n")
  obj <- SetAssayData(obj, layer = "data", new.data = log2(cnt0 + 1))
}
rm(cnt0); gc()

d0 <- GetAssayData(obj, layer = "data")
cat("data max:", round(max(d0), 2),
    "| TPI1 detection:", round(100 * mean(d0["TPI1", ] > 0), 1), "%",
    "(Step 18 recorded 60.6% in T.CD4)\n")
rm(d0); gc()

## ------------------------------------------------- TPI1 split within CD4
md <- data.table(cell = colnames(obj),
                 ctype = obj$cell.types,
                 sample = obj$samples,
                 nFeature = obj$nFeature_RNA,
                 TPI1 = as.numeric(GetAssayData(obj, layer = "data")["TPI1", ]))

CD4_LABEL <- "T.CD4"                     # author label; check the table printed above
stopifnot(CD4_LABEL %in% md$ctype)

cd4 <- md[ctype == CD4_LABEL]
cd4[, r := frank(TPI1, ties.method = "random") / .N, by = sample]
cd4[, grp := fifelse(r > 1 - HI_FRAC, "hi", fifelse(r <= LO_FRAC, "lo", NA_character_))]
cd4 <- cd4[!is.na(grp)]

cat("\n--- depth imbalance BEFORE matching ---\n")
print(cd4[, .(n = .N, median_nFeature = median(nFeature)), by = grp])
cat("  p =", signif(wilcox.test(nFeature ~ grp, data = cd4)$p.value, 3), "\n")

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
cd4m <- rbindlist(lapply(split(cd4, cd4$sample), match_on_depth))
cat("\n--- depth imbalance AFTER matching ---\n")
print(cd4m[, .(n = .N, median_nFeature = median(nFeature)), by = grp])
cat("  p =", signif(wilcox.test(nFeature ~ grp, data = cd4m)$p.value, 3),
    "  <- must be non-significant\n")
fwrite(cd4m, "25a_GSE115978_CD4_TPI1_split.tsv", sep = "\t")

## ---------------------------------------------------------------- CellChat
labels <- setNames(as.character(obj$cell.types), colnames(obj))
labels[cd4m[grp == "hi", cell]] <- "CD4_TPI1hi"
labels[cd4m[grp == "lo", cell]] <- "CD4_TPI1lo"
keep <- setdiff(colnames(obj), setdiff(md[ctype == CD4_LABEL, cell], cd4m$cell))
labels <- labels[keep]
tab <- table(labels); cat("\n--- CellChat groups ---\n"); print(tab)
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
saveRDS(cc, "25_GSE115978_cellchat.rds")

net <- as.data.table(subsetCommunication(cc))
fwrite(net, "25b_GSE115978_LR_all.tsv", sep = "\t")

## focus: everything on the CD4 <-> malignant axis
MAL <- grep("^Mal", unique(net$source), value = TRUE)
ax <- net[(source %in% c("CD4_TPI1hi", "CD4_TPI1lo") & target %in% MAL) |
          (target %in% c("CD4_TPI1hi", "CD4_TPI1lo") & source %in% MAL)]
fwrite(ax, "25c_GSE115978_CD4_tumour_axis.tsv", sep = "\t")
cat("\n=== CD4 <-> malignant links ===\n")
print(ax[order(-prob), .(source, target, pathway_name, interaction_name,
                         prob = signif(prob, 3), pval)][1:min(30, .N)])

pdf("25_GSE115978_bubble_CD4_tumour.pdf", width = 10, height = 9)
print(netVisual_bubble(cc, sources.use = c("CD4_TPI1hi", "CD4_TPI1lo"),
                       remove.isolate = FALSE))
if (length(MAL)) print(netVisual_bubble(cc, sources.use = MAL,
                       targets.use = c("CD4_TPI1hi", "CD4_TPI1lo"), remove.isolate = FALSE))
dev.off()

cat("\nwritten: 25a-25c, 25_GSE115978_cellchat.rds, 1 pdf\n")
