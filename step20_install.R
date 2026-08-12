## Step 20 dependency install -- run ONCE, with R 4.4.1 (D:/R/R-4.4.1)
## Seurat 5.5.1 / data.table / presto / ComplexHeatmap / circlize are already present.
## Missing: CellChat and its hard deps (NMF, BiocNeighbors, Biobase, ggalluvial).

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("BiocNeighbors", "Biobase", "ComplexHeatmap"), ask = FALSE, update = FALSE)

install.packages(c("NMF", "ggalluvial", "svglite", "remotes"))

## CellChat v2 (jinworks fork is the maintained one; sqjin/CellChat is the old v1)
remotes::install_github("jinworks/CellChat", upgrade = "never")

## verify
for (p in c("Seurat", "CellChat", "NMF", "data.table", "ComplexHeatmap")) {
  cat(sprintf("%-16s %s\n", p,
      ifelse(requireNamespace(p, quietly = TRUE),
             as.character(packageVersion(p)), "STILL MISSING")))
}
