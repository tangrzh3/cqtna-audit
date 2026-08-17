# The guards are the part of this package that does not exist in the Python
# reference, so they need their own tests. Each one covers a failure that does
# not announce itself: a wrong build, a silently absent negative control, a
# known-locus list of unknown provenance.

mr_df <- cqtna_demo("mr")
kn_df <- cqtna_demo("known")

test_that("a genome build must be declared", {
  expect_error(as_cqtna_mr(mr_df), "build` is required")
  expect_error(as_cqtna_known(kn_df), "build` is required")
})

test_that("a build mismatch is an error, not a warning", {
  mr <- as_cqtna_mr(mr_df, build = "GRCh38")
  kn <- suppressWarnings(as_cqtna_known(kn_df, build = "GRCh37"))
  expect_error(cqtna_attribution(mr, kn), "genome build mismatch")
  expect_error(cqtna_window_sweep(mr, kn), "genome build mismatch")
})

test_that("a known-locus list without provenance warns", {
  expect_warning(as_cqtna_known(kn_df, build = "GRCh38"), "source` column")
  with_source <- kn_df
  with_source$source <- "published"
  expect_silent(as_cqtna_known(with_source, build = "GRCh38"))
})

test_that("omitting the mismatched control warns rather than passing quietly", {
  expect_warning(
    cqtna_audit(mr_df, kn_df, build = "GRCh38"),
    "no mismatched-list control")
})

test_that("missing required columns name themselves", {
  bad <- mr_df[, setdiff(names(mr_df), "gene")]
  expect_error(as_cqtna_mr(bad, build = "GRCh38"), "missing required column\\(s\\): gene")
  expect_error(cqtna_ladder(data.frame(gene = "A")), "missing column")
  expect_error(cqtna_peak_distance(data.frame(gene = "A")), "missing column")
})

test_that("TwoSampleMR-style column names are accepted", {
  alias <- mr_df
  names(alias)[names(alias) == "record_id"] <- "SNP"
  names(alias)[names(alias) == "p"] <- "pval"
  a <- as_cqtna_mr(alias, build = "GRCh38")
  expect_true(all(c("record_id", "p") %in% names(a)))
  expect_equal(nrow(a), nrow(mr_df))
})

test_that("FDR is recomputed rather than trusted", {
  poisoned <- mr_df
  poisoned$fdr <- 0                       # every row significant, if believed
  a <- as_cqtna_mr(poisoned, build = "GRCh38")
  expect_equal(sum(a$fdr < 0.05), 21L)    # the true count, not 3556
})
