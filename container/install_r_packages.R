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
## Posit's binary mirror serves the CURRENT build of each package and does not
## expose CRAN's Archive at the path renv looks for, so an older locked version
## is unreachable through it. cluster 2.1.6 is the case that surfaced: the lock
## wants 2.1.6, CRAN now ships 2.1.8.3, and 2.1.6 exists only in the Archive
## (verified reachable, HTTP 200). Appending the canonical CRAN mirror keeps
## the fast binaries as first choice while making archived sources findable at
## all. renv searches every configured repository, so this adds reach without
## changing which version any package resolves to -- the lock still decides.
r <- getOption("repos")
if (!"CRANsrc" %in% names(r)) {
  options(repos = c(r, CRANsrc = "https://cloud.r-project.org"))
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
## Drop what no repository can supply, letting R itself say which those are
## rather than hardcoding a list. Base-priority packages (compiler, graphics,
## grDevices, grid, methods, parallel, splines, stats, stats4, tools, utils)
## ship inside R and cannot be installed; renv tried and failed the restore on
## them. The recommended ones in the lock -- MASS, Matrix, survival, lattice,
## nlme, cluster, codetools, KernSmooth -- are NOT dropped: they are on CRAN,
## so renv installs the locked versions rather than whatever R shipped with,
## which is what keeps the closure honest.
##
## cqtna is dropped too. The lock records it as Source "Repository", but it is
## the local package built from cqtna_r/ in this repository and exists in no
## repository at all -- a second thing step150 recorded inaccurately. It is
## installed from source in a later layer, and the build-time check there does
## library(cqtna), so nothing goes unverified.
base_pkgs <- rownames(installed.packages(priority = "base"))
drop <- intersect(names(lk$Packages), c(base_pkgs, "cqtna"))
lk$Packages <- lk$Packages[setdiff(names(lk$Packages), drop)]
cat("dropped from derived lock (unavailable from any repository): ",
    paste(sort(drop), collapse = ", "), "
", sep = "")

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
## ---------------------------------------------------------------------------
## THREE of those are absent from 150a_environment.lock: chromVAR and
## org.Hs.eg.db, which the Methods names WITH versions, and testthat, which S54
## section 9.1 records as a test-time dependency and which step159 phase 1 then
## runs the cqtna suite with.
##
## What used to be here installed whatever was missing with
## BiocManager::install(), no version. cqtna-audit:0.3.0 happened to resolve
## chromVAR 1.26.0 and org.Hs.eg.db 3.19.1 -- exactly the versions the Methods
## states -- and nothing in the image would have said otherwise had it resolved
## anything else. Correct by luck is not correct by record, and DEPOSIT_GAPS.md
## section 1 called this unfixable when the argument it gave established only
## that 150a must not be REGENERATED. Third-party review, item 3.
##
## 150c_unlocked_packages.tsv is the repair: a deposited supplementary record,
## written by hand and never by a script, naming each package the lock omits
## and the version this image must produce. Bioconductor cannot be pinned by
## package version, so those are pinned by RELEASE and then asserted; CRAN
## packages go through remotes::install_version. Every entry is checked whether
## or not it was missing, because the base image supplying one at some other
## version is exactly the drift this exists to catch.
supp <- Sys.getenv("CQTNA_SUPP", "/repo/150c_unlocked_packages.tsv")
if (!file.exists(supp)) {
  stop("150c_unlocked_packages.tsv not found at ", supp, ". The packages the ",
       "lock omits would then be installed unpinned, which is the defect ",
       "DEPOSIT_GAPS.md section 1 records. Refusing to build that image.")
}
sp <- read.delim(supp, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("package", "version", "repository") %in% names(sp))) {
  stop("150c_unlocked_packages.tsv must carry package/version/repository")
}
unrecorded <- setdiff(need[!need %in% names(jsonlite::fromJSON(lock)$Packages)],
                      sp$package)
if (length(unrecorded)) {
  stop("absent from BOTH 150a and 150c: ", paste(unrecorded, collapse = ", "),
       ". Record the version this image must produce in ",
       "150c_unlocked_packages.tsv rather than letting the resolver pick one.")
}
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
for (i in seq_len(nrow(sp))) {
  nm <- sp$package[i]; want <- sp$version[i]; repo <- sp$repository[i]
  have <- tryCatch(as.character(packageVersion(nm)),
                   error = function(e) NA_character_)
  if (!identical(have, want)) {
    if (grepl("^Bioconductor", repo)) {
      want_rel <- sub("^Bioconductor[[:space:]]+", "", repo)
      got_rel <- as.character(BiocManager::version())
      if (!identical(want_rel, got_rel)) {
        stop(nm, " is pinned to Bioconductor ", want_rel, " but this image is ",
             "on Bioconductor ", got_rel, ". Change the base image or amend ",
             "150c; do not install whatever this release happens to carry.")
      }
      BiocManager::install(nm, ask = FALSE, update = FALSE)
    } else {
      remotes::install_version(nm, version = want, upgrade = "never")
    }
  }
  got <- tryCatch(as.character(packageVersion(nm)), error = function(e) "<absent>")
  if (!identical(got, want)) {
    stop(nm, ": 150c requires ", want, ", this image resolved ", got,
         ". Fix the image; do not proceed on a version no record names.")
  }
  cat(sprintf("pinned from 150c: %-30s %s\n", nm, got))
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
             status = if (identical(nm, "cqtna") && is.na(got)) "later-layer"
                      else if (is.na(got)) "absent"
                      else if (identical(got, want)) "match" else "MISMATCH",
             stringsAsFactors = FALSE)
}))
write.table(rows, "/opt/150a_closure_report.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
## cqtna is installed from source in a later Docker layer, so it is absent
## when this runs and is reported as "later-layer" rather than counted as a
## failure -- the build-time check there loads it. Miscounting it would train
## the reader to ignore this report, which is the one thing it must not be.
## The supplementary record gets its own report, deliberately NOT merged into
## the 226-row one: "226/226" has to keep meaning "everything 150a records is
## present at the recorded version". What it never meant, and what
## DEPOSIT_GAPS.md section 1 spells out, is "everything the paper uses is
## pinned" -- the packages outside the lock were outside the closure check's
## frame too. This is that frame's missing edge, reported beside it.
supp_rows <- do.call(rbind, lapply(seq_len(nrow(sp)), function(i) {
  nm <- sp$package[i]
  got <- tryCatch(as.character(packageVersion(nm)), error = function(e) "<absent>")
  data.frame(package = nm, recorded = sp$version[i], installed = got,
             status = if (identical(got, sp$version[i])) "match" else "MISMATCH",
             source = "150c_unlocked_packages.tsv", stringsAsFactors = FALSE)
}))
write.table(supp_rows, "/opt/150c_closure_report.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\n150c supplementary:", sum(supp_rows$status == "match"), "match,",
    sum(supp_rows$status == "MISMATCH"), "mismatch, of", nrow(supp_rows), "\n")

cat("\n150a closure:", sum(rows$status == "match"), "match,",
    sum(rows$status == "MISMATCH"), "mismatch,",
    sum(rows$status == "absent"), "absent, of", nrow(rows), "\n")
if (any(rows$status == "MISMATCH")) {
  cat("MISMATCHED (locked -> installed):\n")
  bad <- rows[rows$status == "MISMATCH", ]
  cat(paste0("  ", bad$package, " ", bad$locked, " -> ", bad$installed,
             collapse = "\n"), "\n")
}
