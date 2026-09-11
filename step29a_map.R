suppressPackageStartupMessages(library(org.Hs.eg.db))
setwd(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else if (nzchar(Sys.getenv("CQTNA_DIR"))) Sys.getenv("CQTNA_DIR") else "D:/R_ex/MR")
g <- c("SLC2A1","SLC2A3","HK1","HK2","HK3","GPI","PFKL","PFKM","PFKP","ALDOA",
       "ALDOB","ALDOC","TPI1","GAPDH","PGK1","PGAM1","ENO1","ENO2","ENO3","PKM",
       "PKLR","LDHA","LDHB","PGM1","PFKFB3","PFKFB4","SLC16A1","SLC16A3")
m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = g,
                                            keytype = "SYMBOL", columns = "ENSEMBL"))
m <- m[!is.na(m$ENSEMBL), ]
write.table(m, "35a_glyco_ensembl_map.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("mapped", nrow(m), "rows for", length(unique(m$SYMBOL)), "of", length(g), "genes\n")
cat("unmapped:", paste(setdiff(g, m$SYMBOL), collapse = ", "), "\n")
