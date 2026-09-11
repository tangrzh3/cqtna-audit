## Verify the TPI1 pseudogene Ensembl IDs properly rather than from memory,
## and record their genomic locations -- a cross-chromosome pseudogene cannot
## produce a clean cis peak at TPI1's own promoter, so location matters.
suppressPackageStartupMessages({library(org.Hs.eg.db); library(AnnotationDbi)})
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")

sym <- keys(org.Hs.eg.db, keytype = "SYMBOL")
pg <- grep("^TPI1", sym, value = TRUE)
cat("symbols starting with TPI1:", paste(pg, collapse = ", "), "\n\n")

m <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = pg, keytype = "SYMBOL",
  columns = c("ENSEMBL", "CHR", "GENENAME")))
print(m)
write.table(m, "40b_TPI1_pseudogene_ids.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nwritten: 40b_TPI1_pseudogene_ids.tsv\n")
