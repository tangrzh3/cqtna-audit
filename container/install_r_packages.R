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

## ---------------------------------------------------------------------------
## Closure check, written into the image at build time.
##
## WHY: the fallback above installs whatever the snapshot currently holds for
## anything renv failed to restore, and it is NOT version-pinned. That is
## precisely the defect S54 section 9.9 records: the slim acceptance image
## pinned six named packages and let every transitive dependency float, so 42
## of the lock's 226 matched and 11 were at outright different versions --
## discovered only after two acceptance runs had already been reported.
## S54 section 8 makes "dependency resolution does not match 150a" a stop
## condition, so the image must be able to state whether it matches rather than
## leaving a human to check by hand afterwards.
##
## This does not fail the build: whether a mismatch blocks the run is the
## acceptance gate's decision, and step159 phase 1 reads this file to make it.
inst <- as.data.frame(installed.packages()[, c("Package", "Version")],
                      stringsAsFactors = FALSE)
rownames(inst) <- NULL
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
lock_pkgs <- jsonlite::fromJSON(lock)$Packages
rows <- do.call(rbind, lapply(names(lock_pkgs), function(nm) {
  want <- lock_pkgs[[nm]]$Version
  got <- inst$Version[match(nm, inst$Package)]
  data.frame(package = nm, locked = want,
             installed = if (is.na(got)) "<absent>" else got,
             status = if (is.na(got)) "absent"
                      else if (identical(got, want)) "match" else "MISMATCH",
             stringsAsFactors = FALSE)
}))
write.table(rows, "/opt/150a_closure_report.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\n150a closure:", sum(rows$status == "match"), "match,",
    sum(rows$status == "MISMATCH"), "mismatch,",
    sum(rows$status == "absent"), "absent, of", nrow(rows), "\n")
if (any(rows$status == "MISMATCH")) {
  cat("MISMATCHED (locked -> installed):\n")
  bad <- rows[rows$status == "MISMATCH", ]
  cat(paste0("  ", bad$package, " ", bad$locked, " -> ", bad$installed,
             collapse = "\n"), "\n")
}
