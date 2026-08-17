# The regression suite for this package is a published paper's results.
#
# Every number below is reported in the manuscript this tool came out of, and is
# reproduced by the Python reference implementation in demo_out/. If a change to
# the R code moves any of them, the change is wrong until argued otherwise.

# Values come back at full precision; round them the way the report does.
r2 <- function(v) as.numeric(sprintf("%.2f", v))

mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38")
kn <- suppressWarnings(as_cqtna_known(cqtna_demo("known"), build = "GRCh38"))
mm <- suppressWarnings(as_cqtna_known(cqtna_demo("mismatch"), build = "GRCh38"))

test_that("A: locus attribution reproduces the published melanoma numbers", {
  a <- cqtna_attribution(mr, kn)
  expect_equal(a$background_loci, 554L)
  expect_equal(a$background_known, 58L)
  expect_equal(r2(a$background_pct), 10.47)
  expect_equal(a$significant_loci, 7L)
  expect_equal(a$significant_known, 3L)
  expect_equal(as.numeric(sprintf("%.1f", a$pct_known)), 42.9)
  expect_equal(r2(a$fold), 4.09)
  # the published one-sided Fisher P, to full precision
  expect_equal(a$fisher_p_one_sided, 0.028123181205575935, tolerance = 1e-12)
})

test_that("A: the mismatched-list control is clean", {
  nc <- cqtna_attribution(mr, mm, label = "mismatched")
  expect_equal(nc$fold, 0)
  expect_equal(nc$fisher_p_one_sided, 1, tolerance = 1e-12)
})

test_that("B: the same list spans 10 genes by record and 28 by locus", {
  b <- cqtna_unit_sweep(mr)
  expect_equal(b$n_tests[b$unit == "record"], 3556L)
  expect_equal(b$n_significant[b$unit == "record"], 21L)
  expect_equal(b$n_genes[b$unit == "record"], 10L)
  expect_equal(b$n_genes[b$unit == "locus"], 28L)
  expect_equal(b$n_tests[b$unit == "locus"], 554L)
})

test_that("C: the two outcome rounds share 5 of 11 genes", {
  alt <- as_cqtna_mr(cqtna_demo("mr_alt"), build = "GRCh38")
  s <- cqtna_stability(mr, alt)
  expect_equal(s$genes_outcome1, 10L)
  expect_equal(s$genes_outcome2, 6L)
  expect_equal(s$genes_shared, 5L)
  expect_equal(as.numeric(sprintf("%.3f", s$jaccard_gene)), 0.455)
})

test_that("D: the instrument ladder is 28 -> 3 -> 2 -> 1", {
  d <- cqtna_ladder(cqtna_demo("instruments"))
  expect_equal(d$pathway_genes, 28L)
  expect_equal(d$instrumentable, 3L)
  expect_equal(d$analysable, 2L)
  expect_equal(d$associated, 1L)
})

test_that("E: TPI1 is 6.7-fold higher in malignant cells than in CD4 T cells", {
  e <- cqtna_compartment(cqtna_demo("expression"), "TPI1", "CD4_T")
  expect_equal(e$highest_other, "Malignant")
  expect_equal(r2(e$ratio_other_over_target), 6.7)
})

test_that("G: tightening either window raises the fold", {
  g <- cqtna_window_sweep(mr, kn)
  ks <- g$known_sweep
  expect_equal(r2(ks$fold), c(18.26, 8.19, 5.65, 4.09))
  expect_equal(ks$background_known, c(13L, 29L, 42L, 58L))
  expect_equal(ks$fisher_p_one_sided[1], 0.0003361896157234458, tolerance = 1e-12)
  ls <- g$locus_sweep
  # 5.325000000000000177 -- round() would give 5.32 here, sprintf gives 5.33.
  # The paper and the Python reference both report 5.33; see R/report.R.
  expect_equal(r2(ls$fold), c(5.33, 4.63, 4.59, 4.09))
  expect_equal(ls$background_loci, c(1065L, 876L, 728L, 554L))
  # a threshold that flattered the result would weaken when tightened
  expect_true(all(diff(ks$fold) < 0))
})

test_that("G: the 1 Mb cut falls in a gap in the distances", {
  g <- cqtna_window_sweep(mr, kn)
  expect_equal(round(g$significant_distances_bp / 1000),
               c(5, 11, 74, 2167, 7275, 20136, 38494))
  expect_true(g$threshold_gap$in_gap)
  expect_equal(round(g$threshold_gap$nearest_below_bp / 1000), 74)
  expect_equal(round(g$threshold_gap$nearest_above_bp / 1000), 2167)
})

test_that("no gene reaches target-supported, and the tiers hold", {
  au <- cqtna_audit(cqtna_demo("mr"), cqtna_demo("known"),
                    mismatch = cqtna_demo("mismatch"),
                    mr_alt = cqtna_demo("mr_alt"),
                    instruments = cqtna_demo("instruments"),
                    expression = cqtna_demo("expression"),
                    target_cell_type = "CD4_T", build = "GRCh38")
  tiers <- vapply(au$tiers, function(r) r$tier, character(1))
  expect_false(any(tiers == "target-supported"))
  expect_true(all(tiers %in% c("screened", "unresolved", "state-informative")))
  expect_length(au$not_automated, 3L)
  # the refusal must survive into the rendered report, not only the object
  txt <- cqtna_report(au)
  expect_match(txt, "No gene can reach `target-supported`", fixed = TRUE)
  expect_match(txt, "five of eight", fixed = TRUE)
})
