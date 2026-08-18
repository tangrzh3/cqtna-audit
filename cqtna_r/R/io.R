## Input coercion and the genome-build guard.

CQ_REQUIRED_MR <- c("record_id", "gene", "chr", "pos", "p")
CQ_REQUIRED_KNOWN <- c("chr", "pos")

#' Coerce MR results into the form CQTNA audits
#'
#' Accepts any data frame with one row per test. Column names follow the
#' `TwoSampleMR` convention where they overlap, and are renamed if a recognised
#' alias is present: `SNP`/`pval` and so on.
#'
#' @param mr data frame with `record_id`, `gene`, `chr`, `pos`, `p`. An
#'   `exposure_profile` column is kept if present.
#' @param locus_kb the locus window in kb. Under `"fixed_centre"` it is the
#'   radius claimed from each centre, bounding a locus at that span; under
#'   `"single_linkage"` it is the joining distance and the span is unbounded.
#' @param locus_method how variants are partitioned into loci. `"fixed_centre"`
#'   is non-recursive and cannot chain, and is the default because a partition
#'   whose blocks can reach tens of megabases is not a partition into loci.
#'   `"single_linkage"` reproduces the source study's published numbers.
#'   `"blocks"` assigns by supplied intervals, e.g. LD blocks.
#' @param blocks for `locus_method = "blocks"`, a data frame with `chr`,
#'   `start`, `end`.
#' @param build genome build label, e.g. `"GRCh38"`. Carried on the object and
#'   checked against the known-locus list. There is no default: a silent build
#'   mismatch produces numbers that look entirely normal.
#' @return a `cqtna_mr` data frame with `locus` and `fdr` added.
#' @export
as_cqtna_mr <- function(mr, locus_kb = 1000, build = NULL,
                        locus_method = c("fixed_centre", "single_linkage", "blocks"),
                        blocks = NULL) {
  locus_method <- match.arg(locus_method)
  if (is.null(build) || !nzchar(build))
    stop("`build` is required. State the genome build of `mr` (e.g. \"GRCh38\").\n",
         "  A known-locus list on a different build gives a plausible-looking ",
         "answer that is wrong throughout.", call. = FALSE)
  d <- as.data.frame(mr, stringsAsFactors = FALSE)
  d <- cq_rename(d, c(SNP = "record_id", pval = "p", p_value = "p",
                      P = "p", chrom = "chr", chromosome = "chr",
                      position = "pos", bp = "pos", gene_name = "gene",
                      symbol = "gene"))
  missing <- setdiff(CQ_REQUIRED_MR, names(d))
  if (length(missing))
    stop("mr is missing required column(s): ", paste(missing, collapse = ", "),
         ".\n  Required: ", paste(CQ_REQUIRED_MR, collapse = ", "), call. = FALSE)
  cq_validate_window(locus_kb, "locus_kb")
  d <- d[!is.na(d$p), , drop = FALSE]
  if (!nrow(d)) stop("mr has no rows with a non-missing p-value.", call. = FALSE)
  if (anyDuplicated(d$record_id))
    warning("record_id is not unique; rows are audited as supplied.", call. = FALSE)
  d$chr <- as.character(d$chr)
  d$pos <- suppressWarnings(as.numeric(d$pos))
  d$p <- suppressWarnings(as.numeric(d$p))
  d$gene <- as.character(d$gene)
  d$record_id <- as.character(d$record_id)
  cq_validate_mr(d)
  d$locus <- cq_assign_loci(d$chr, d$pos, locus_kb, locus_method, blocks)
  d$fdr <- cq_bh(d$p)          # recomputed, never read from the input
  attr(d, "build") <- build
  attr(d, "locus_kb") <- locus_kb
  attr(d, "locus_method") <- locus_method
  attr(d, "blocks") <- blocks
  class(d) <- c("cqtna_mr", "data.frame")
  d
}

#' Coerce a known-locus list
#'
#' @param known data frame of previously reported lead SNPs **for your outcome
#'   trait**, with `chr` and `pos`. A `source` column is strongly recommended;
#'   see the note on circularity in the README.
#' @param build genome build label; must match the MR table.
#' @return a `cqtna_known` data frame.
#' @export
as_cqtna_known <- function(known, build = NULL) {
  if (is.null(build) || !nzchar(build))
    stop("`build` is required for the known-locus list too.", call. = FALSE)
  k <- as.data.frame(known, stringsAsFactors = FALSE)
  k <- cq_rename(k, c(chrom = "chr", chromosome = "chr", position = "pos",
                      bp = "pos"))
  missing <- setdiff(CQ_REQUIRED_KNOWN, names(k))
  if (length(missing))
    stop("known is missing required column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  k$chr <- as.character(k$chr)
  k$pos <- suppressWarnings(as.numeric(k$pos))
  cq_validate_known(k)
  if (!"source" %in% names(k))
    warning("known has no `source` column. Keep the provenance of each lead SNP: ",
            "if any of it comes from the outcome GWAS's own publication, asking ",
            "whether your hits fall on known loci is circular to that extent.",
            call. = FALSE)
  attr(k, "build") <- build
  class(k) <- c("cqtna_known", "data.frame")
  k
}

#' @keywords internal
#' @noRd
cq_rename <- function(d, map) {
  for (from in names(map)) {
    to <- map[[from]]
    if (from %in% names(d) && !(to %in% names(d))) names(d)[names(d) == from] <- to
  }
  d
}

#' @keywords internal
#' @noRd
cq_check_build <- function(mr, known) {
  a <- attr(mr, "build")
  b <- attr(known, "build")
  if (!identical(a, b))
    stop("genome build mismatch: mr is ", a, ", known-locus list is ", b, ".\n",
         "  Lift one over before auditing. A mismatch here does not error ",
         "downstream; it just moves every locus.", call. = FALSE)
  cq_check_chr_style(mr$chr, known$chr)
  invisible(TRUE)
}

#' Read the melanoma demo data shipped with the package
#'
#' The audit's own published data: 3,556 MR records against a 12,530-case
#' melanoma meta-analysis, the 157-locus melanoma known-locus list, and an HCC
#' list as the mismatched negative control.
#'
#' @param what one of `"mr"`, `"mr_alt"`, `"known"`, `"mismatch"`,
#'   `"instruments"`, `"expression"`.
#' @return a data frame
#' @export
#' @examples
#' mr <- cqtna_demo("mr")
#' nrow(mr)
cqtna_demo <- function(what = c("mr", "mr_alt", "known", "mismatch",
                                "instruments", "expression")) {
  what <- match.arg(what)
  f <- c(mr = "mr_meta.tsv", mr_alt = "mr_finngen.tsv",
         known = "known_melanoma.tsv", mismatch = "known_hcc_mismatched.tsv",
         instruments = "instruments.tsv", expression = "expression.tsv")[[what]]
  utils::read.delim(system.file("extdata", f, package = "cqtna"),
                    stringsAsFactors = FALSE)
}
