## Restore the R library recorded by step150.
##
## 150a_environment.lock is not an renv-managed lockfile -- it is an export of
## the library that produced the result tables, written in renv.lock's JSON
## structure precisely so that renv::restore() can read it (step150 header).
## Restoring it therefore reproduces the package VERSIONS, not the machine.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

lock <- Sys.getenv("CQTNA_LOCK", "/repo/150a_environment.lock")
if (!file.exists(lock)) stop("lockfile not found: ", lock)

renv::restore(lockfile = lock, prompt = FALSE)

## The lock records only what was installed on the authoring machine. Anything
## it misses fails loudly here rather than silently at analysis time.
need <- c("susieR", "data.table", "coloc", "Matrix", "arrow", "dplyr",
          "Seurat", "hdf5r", "TFBSTools", "motifmatchr", "JASPAR2020",
          "Biostrings", "BSgenome.Hsapiens.UCSC.hg38", "chromVAR",
          "survival", "org.Hs.eg.db", "testthat")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("not restored from the lock, installing from the snapshot: ",
          paste(missing, collapse = ", "))
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
missing2 <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing2)) stop("still missing: ", paste(missing2, collapse = ", "))
