## S38-A1: the permutation must not depend on the order rows arrive in.
##
## A third-party reviewer found that reversing the input rows moved the empirical
## p-value. Two causes, both artefacts: the loop over significant loci consumed
## the RNG stream in input order, so the same stream was cut into different
## pieces; and sampling without replacement gave whichever locus went first the
## widest pool. These tests exist so neither can come back silently.

demo_inputs <- function() {
  mr <- cqtna_demo("mr")
  kn <- suppressWarnings(as_cqtna_known(cqtna_demo("known"), build = "GRCh38"))
  list(mr = mr, kn = kn)
}

run_perm <- function(raw, kn, ...) {
  mr <- as_cqtna_mr(raw, locus_kb = 1000, build = "GRCh38",
                    locus_method = "fixed_centre")
  cqtna_permutation_control(mr, kn, known_from = "any_record", tolerance = 1,
                            n_perm = 400, seed = 1, ...)
}

test_that("the empirical p-value does not depend on input row order", {
  d <- demo_inputs()
  raw <- d$mr

  a <- run_perm(raw, d$kn)
  b <- run_perm(raw[rev(seq_len(nrow(raw))), ], d$kn)
  set.seed(99)
  cc <- run_perm(raw[sample(nrow(raw)), ], d$kn)

  expect_equal(a$empirical_p, b$empirical_p)
  expect_equal(a$empirical_p, cc$empirical_p)
  expect_equal(a$fold_vs_null, b$fold_vs_null)
  expect_equal(a$fold_vs_null, cc$fold_vs_null)
  expect_equal(a$observed_known, b$observed_known)
})

test_that("the canonical locus order is numeric, not lexical", {
  loci <- c("2:100-200", "10:100-200", "1:900-1000", "1:100-200", "X:5-6")
  ord <- cqtna:::cq_locus_order(loci)
  expect_equal(loci[ord][1:4],
               c("1:100-200", "1:900-1000", "2:100-200", "10:100-200"))
  # a lexical sort would have put chr10 before chr2
  expect_false(identical(loci[ord], sort(loci)))
})

test_that("diagnostics come back with every successful result", {
  d <- demo_inputs()
  r <- run_perm(d$mr, d$kn)
  for (f in c("n_draws_used", "n_draws_exhausted", "effective_draw_fraction",
              "min_pool_size", "median_pool_size", "max_pool_size",
              "rng_kind", "r_version", "collate"))
    expect_false(is.null(r[[f]]), info = f)
  expect_equal(r$n_draws_used + r$n_draws_exhausted, r$n_perm)
  expect_true(r$effective_draw_fraction <= 1)
})

test_that("heavily overlapping pools are reported, not quietly averaged over", {
  ## Significant loci whose matched pools are a single shared background locus:
  ## after the first draw there is nothing left, so most draws must be discarded
  ## and the result must refuse to give a p-value.
  set.seed(4)
  mk <- function(chr, pos, p, gene) data.frame(
    record_id = seq_along(pos), gene = gene, chr = as.character(chr),
    pos = pos, p = p, stringsAsFactors = FALSE)
  sig_pos <- c(1e6, 2e7, 4e7, 6e7)
  raw <- rbind(
    mk(1, sig_pos, rep(1e-12, 4), paste0("S", seq_along(sig_pos))),
    mk(1, 8e7, 0.9, "B1"))                      # one lone background locus
  raw$record_id <- seq_len(nrow(raw))           # unique, or as_cqtna_mr warns
  kn <- suppressWarnings(as_cqtna_known(
    data.frame(chr = "1", pos = 1e6), build = "GRCh38"))
  mr <- as_cqtna_mr(raw, locus_kb = 1000, build = "GRCh38",
                    locus_method = "fixed_centre")

  r <- cqtna_permutation_control(mr, kn, known_from = "any_record",
                                 tolerance = 10, n_perm = 200, seed = 1)
  ## Four significant loci sharing one background locus: the first draw takes
  ## it and the rest have nothing left, so every draw is discarded.
  expect_true(r$n_draws_exhausted > 0)
  expect_equal(r$n_draws_used, 0L)
  expect_equal(r$effective_draw_fraction, 0)
  expect_true(is.na(r$empirical_p))
  expect_match(r$failed_because, "exhaust")
  expect_false(is.null(r$min_pool_size))
})

test_that("min_draw_fraction is validated", {
  d <- demo_inputs()
  expect_error(run_perm(d$mr, d$kn, min_draw_fraction = 1.5),
               "min_draw_fraction")
})
