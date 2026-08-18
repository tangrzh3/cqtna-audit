## Input contract. Every check here rejects something that otherwise produces a
## number rather than an error -- that is the selection criterion for what is
## worth validating.

#' @keywords internal
#' @noRd
cq_bad_rows <- function(ok, max_show = 5L) {
  i <- which(!ok)
  if (!length(i)) return("")
  paste0("\n  offending row(s): ", paste(utils::head(i, max_show), collapse = ", "),
         if (length(i) > max_show) sprintf(" ... and %d more", length(i) - max_show) else "")
}

#' @keywords internal
#' @noRd
cq_stop_col <- function(col, what, ok) {
  stop(sprintf("column `%s`: %s", col, what), cq_bad_rows(ok), call. = FALSE)
}

#' Coerce a 0/1 flag column
#'
#' Accepts logical, numeric and character. Rejects factors outright: a factor of
#' "0"/"1" taken through as.numeric() becomes 1/2, so every "0" counts as TRUE.
#' @keywords internal
#' @noRd
cq_as_flag <- function(v, col) {
  if (is.factor(v))
    stop(sprintf("column `%s` is a factor. Convert it first: as.numeric(as.character(x)) ",
                 col), "or as.logical(x).\n",
         "  A factor of \"0\"/\"1\" read through as.numeric() becomes 1/2, so every ",
         "\"0\" would count as TRUE.", call. = FALSE)
  if (is.logical(v)) return(v)
  if (is.character(v)) {
    t <- trimws(v)
    ok <- is.na(t) | t %in% c("0", "1", "TRUE", "FALSE", "T", "F", "true", "false")
    if (!all(ok)) cq_stop_col(col, "character flag must be 0/1/TRUE/FALSE", ok)
    return(t %in% c("1", "TRUE", "T", "true"))
  }
  v <- as.numeric(v)
  ok <- is.na(v) | v %in% c(0, 1)
  if (!all(ok)) cq_stop_col(col, "numeric flag must be 0 or 1", ok)
  !is.na(v) & v == 1
}

#' @keywords internal
#' @noRd
cq_check_chr_style <- function(a, b, what_a = "mr", what_b = "known") {
  pa <- any(grepl("^chr", a, ignore.case = TRUE))
  pb <- any(grepl("^chr", b, ignore.case = TRUE))
  if (pa != pb)
    stop("chromosome naming differs: ", what_a, " uses ",
         if (pa) "\"chr1\"" else "\"1\"", " and ", what_b, " uses ",
         if (pb) "\"chr1\"" else "\"1\"", " style.\n",
         "  Harmonise them; they will never match otherwise and the audit ",
         "returns a background of zero known loci.", call. = FALSE)
  shared <- intersect(unique(a), unique(b))
  if (!length(shared))
    warning("no chromosome in the MR table appears in the known-locus list. ",
            "Check the naming convention on both.", call. = FALSE)
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
cq_validate_window <- function(v, name) {
  if (is.null(v) || !length(v) || !all(is.finite(v)) || any(v < 0))
    stop(sprintf("`%s` must be finite and non-negative; got %s", name,
                 paste(utils::capture.output(dput(v)), collapse = "")), call. = FALSE)
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
cq_validate_fdr <- function(v) {
  if (length(v) != 1L || !is.finite(v) || v <= 0 || v >= 1)
    stop("`fdr` must be a single number in (0, 1); got ", format(v), call. = FALSE)
  invisible(TRUE)
}

#' Validate the MR table
#'
#' Rejects the inputs that would otherwise flow through and produce a plausible
#' number: non-finite or out-of-range p, missing or non-positive coordinates,
#' blank identifiers.
#' @keywords internal
#' @noRd
cq_validate_mr <- function(d) {
  ok <- is.finite(d$p) & d$p >= 0 & d$p <= 1
  if (!all(ok)) cq_stop_col("p", "must be finite and within [0, 1]", ok)

  ok <- is.finite(d$pos) & d$pos > 0
  if (!all(ok)) cq_stop_col("pos", "must be a finite positive number", ok)

  for (col in c("chr", "gene", "record_id")) {
    v <- as.character(d[[col]])
    ok <- !is.na(v) & nzchar(trimws(v))
    if (!all(ok)) cq_stop_col(col, "must not be missing or blank", ok)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
cq_validate_known <- function(k) {
  ok <- is.finite(k$pos) & k$pos > 0
  if (!all(ok)) cq_stop_col("pos", "must be a finite positive number", ok)
  v <- as.character(k$chr)
  ok <- !is.na(v) & nzchar(trimws(v))
  if (!all(ok)) cq_stop_col("chr", "must not be missing or blank", ok)
  invisible(TRUE)
}
