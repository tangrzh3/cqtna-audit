## Restore the R library recorded by step150.
##
## 150a_environment.lock is not an renv-managed lockfile -- it is an export of
## the library that produced the result tables, written in renv.lock's JSON
## structure precisely so that renv::restore() can read it (step150 header).
## Restoring it therefore reproduces the package VERSIONS, not the machine.
## Point CRAN at the Posit binary mirror for speed, then let BiocManager
## reinstate the Bioconductor repositories ON TOP of it.
##
## The single line that used to be here set repos to CRAN and nothing else,
## which silently discarded the Bioconductor repos the base image configures.
## 207 of the lock's 226 entries carry Source "Repository" with no repository
## named and the lockfile has no Bioconductor field, so renv finds a package
## only by searching whatever repos are configured -- and every Bioconductor
## package therefore became "failed to find source": S4Vectors, IRanges,
## GenomicRanges, GO.db, TFBSTools, JASPAR2020 and the rest, which then took
## restfulr and the whole restore down with them. BiocManager::repositories()
## returns the Bioc repos for this image's release plus the CRAN already set,
## so both survive. The base is RELEASE_3_19 and the locked versions
## (S4Vectors 0.42.1, GO.db 3.19.1, GenomicRanges 1.56.2) are Bioc 3.19, so
## the release repos carry exactly what the lock asks for.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
if (requireNamespace("BiocManager", quietly = TRUE)) {
  options(repos = BiocManager::repositories())
}
cat("repos in use:\n"); print(getOption("repos"))
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

lock <- Sys.getenv("CQTNA_LOCK", "/repo/150a_environment.lock")
if (!file.exists(lock)) stop("lockfile not found: ", lock)

## ---------------------------------------------------------------------------
## Restore from a DERIVED lockfile, not the deposited one.
##
## 150a_environment.lock records exactly one repository:
##     "Repositories": [ { "Name": "CRAN", "URL": "https://cloud.r-project.org" } ]
## and renv::restore() uses the lockfile's own repository list in preference to
## whatever options(repos=) says. Over forty of the lock's packages are
## Bioconductor, so renv searched CRAN alone and reported "failed to find
## source" for every one of them -- S4Vectors, IRanges, GenomicRanges, GO.db,
## TFBSTools and the rest -- then failed the whole restore. Confirmed against
## the live repositories that this is not a missing-version problem: every
## locked Bioconductor version is still served by the 3.19 repos today
## (S4Vectors 0.42.1, GO.db 3.19.1, GenomicRanges 1.56.2 all match exactly).
## The lock simply never recorded where to look.
##
## ⚠ The deposited lockfile is NOT edited: it is the record of the authoring
## machine and altering it would destroy the thing being compared against.
## What is derived here adds ONLY the repository list and the Bioconductor
## release; every package version comes through untouched, which the closure
## check at the end of this file then verifies against the original.
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
lk <- jsonlite::fromJSON(lock, simplifyVector = FALSE)
rp <- getOption("repos")
lk$R$Repositories <- lapply(names(rp), function(n) list(Name = n, URL = unname(rp[[n]])))
if (requireNamespace("BiocManager", quietly = TRUE)) {
  lk$Bioconductor <- list(Version = as.character(BiocManager::version()))
}
derived <- "/tmp/150a_derived.lock"
writeLines(jsonlite::toJSON(lk, auto_unbox = TRUE, pretty = TRUE), derived)
cat("derived lockfile: ", length(lk$Packages), " packages, ",
    length(lk$R$Repositories), " repositories, Bioc ",
    if (is.null(lk$Bioconductor)) "unset" else lk$Bioconductor$Version, "\n", sep = "")
renv::restore(lockfile = derived, prompt = FALSE)

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
