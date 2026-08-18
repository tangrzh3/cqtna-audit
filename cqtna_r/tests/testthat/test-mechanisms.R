# Synthetic tests for the failure mechanisms.
#
# The paper-numbers suite cannot catch these: on the demo data every one of the
# bugs below is latent, because the two conventions happen to agree there. That
# is precisely the argument for testing mechanisms with data built to separate
# them, rather than testing outputs with data that does not.
#
# Each block names the mechanism, then constructs the smallest table that makes
# the wrong answer differ from the right one.

mk_mr <- function(chr, pos, gene, p, build = "GRCh38", locus_kb = 1000,
                  locus_method = "single_linkage") {
  as_cqtna_mr(data.frame(record_id = paste0("r", seq_along(pos)), gene = gene,
                         chr = as.character(chr), pos = pos, p = p,
                         stringsAsFactors = FALSE),
              locus_kb = locus_kb, build = build, locus_method = locus_method)
}
mk_known <- function(chr, pos, build = "GRCh38") {
  suppressWarnings(as_cqtna_known(
    data.frame(chr = as.character(chr), pos = pos, stringsAsFactors = FALSE),
    build = build))
}

# --------------------------------------------------------------------------
# Mechanism 1: several genes at one locus, only some records near a known SNP.
# Locus status is a property of the locus, so every gene there must inherit it.
test_that("all genes at one locus share the locus's known/novel status", {
  # three genes within 1 Mb -> one locus; only GENE_A sits near the known SNP
  mr <- mk_mr(chr = c(1, 1, 1), pos = c(1e6, 1.2e6, 1.4e6),
              gene = c("GENE_A", "GENE_B", "GENE_C"), p = c(1e-8, 1e-8, 1e-8))
  kn <- mk_known(1, 1.05e6)
  a <- cqtna_attribution(mr, kn)

  expect_equal(a$significant_loci, 1L)
  expect_equal(a$significant_known, 1L)
  # the bug: GENE_B and GENE_C would have been labelled novel on a known locus
  expect_setequal(a$known_genes, c("GENE_A", "GENE_B", "GENE_C"))
  expect_length(a$novel_genes, 0L)
  expect_length(intersect(a$known_genes, a$novel_genes), 0L)
})

# --------------------------------------------------------------------------
# Mechanism 2: a NON-significant record is what puts the locus near a known SNP.
# The old code took the denominator over all records and the numerator over
# significant records only, so the same locus could be known in one and novel in
# the other.
test_that("a locus is known even when only a non-significant record is near", {
  mr <- mk_mr(chr = c(1, 1, 2, 2), pos = c(1e6, 1.3e6, 5e6, 5.3e6),
              gene = c("SIG", "NONSIG", "OTHER1", "OTHER2"),
              p = c(1e-9, 0.9, 0.9, 0.9))
  kn <- mk_known(1, 1.35e6)          # within 1 Mb of NONSIG, 350 kb from SIG
  a <- cqtna_attribution(mr, kn, known_from = "any_record")

  expect_equal(a$significant_loci, 1L)
  expect_equal(a$significant_known, 1L)      # inherited from the locus
  expect_setequal(a$known_genes, "SIG")

  # This is the case that separates the conventions, so assert all three.
  # Under "any_record" the locus is known at every window, because NONSIG is
  # 50 kb away. Under the published "significant_records" convention only the
  # 1 Mb window reaches SIG, which is 350 kb away.
  sw <- function(k) cqtna_window_sweep(mr, kn, sweep_kb = c(100, 500, 1000),
                                       known_from = k)$known_sweep$significant_known
  # SIG is 350 kb from the known SNP, NONSIG is 50 kb from it.
  expect_equal(sw("any_record"), c(1L, 1L, 1L))          # NONSIG carries it at 100 kb
  expect_equal(sw("significant_records"), c(0L, 1L, 1L)) # SIG only reaches at 500 kb
  expect_equal(sw("lead_variant"), c(0L, 1L, 1L))        # lead here IS SIG

  # both distance series are always reported, so neither convention can hide
  # which record supplied the proximity
  g <- cqtna_window_sweep(mr, kn, sweep_kb = c(100, 500, 1000))
  expect_equal(round(g$significant_distances_bp / 1000), 50)
  expect_equal(round(g$significant_distances_sig_records_bp / 1000), 350)
})

# --------------------------------------------------------------------------
# Mechanism 3: two outcome tables whose internal locus numbering collides across
# chromosomes. This is the bug that made chr1 match chr16 in the demo data.
test_that("loci on different chromosomes never match across outcome tables", {
  a <- mk_mr(chr = c(1, 5), pos = c(1e6, 9e6), gene = c("G1", "G5"),
             p = c(1e-9, 1e-9))
  b <- mk_mr(chr = c(2, 5), pos = c(1e6, 9e6), gene = c("G2", "G5"),
             p = c(1e-9, 1e-9))
  s <- cqtna_stability(a, b)

  # chr1 and chr2 are each unique; only chr5 is genuinely shared
  expect_equal(s$loci_outcome1, 2L)
  expect_equal(s$loci_outcome2, 2L)
  expect_equal(s$loci_shared, 1L)
  expect_equal(s$jaccard_locus, 1 / 3)
  expect_true(all(grepl("^5:", s$shared_locus_keys)))
})

test_that("far-apart loci on the same chromosome do not match either", {
  a <- mk_mr(chr = c(1, 1), pos = c(1e6, 50e6), gene = c("GA", "GB"),
             p = c(1e-9, 1e-9))
  b <- mk_mr(chr = c(1, 1), pos = c(1e6, 90e6), gene = c("GA", "GC"),
             p = c(1e-9, 1e-9))
  s <- cqtna_stability(a, b)
  expect_equal(s$loci_shared, 1L)            # only the 1 Mb one
  expect_equal(s$jaccard_locus, 1 / 3)
})

# --------------------------------------------------------------------------
# Mechanism 4: coverage differs between the two outcome tables, so "lost" would
# otherwise absorb "never tested".
test_that("genes never tested in the second outcome are reported separately", {
  a <- mk_mr(chr = c(1, 3), pos = c(1e6, 1e6), gene = c("SHARED", "ONLY_IN_1"),
             p = c(1e-9, 1e-9))
  b <- mk_mr(chr = 1, pos = 1e6, gene = "SHARED", p = 1e-9)
  s <- cqtna_stability(a, b)
  expect_equal(s$lost, "ONLY_IN_1")
  expect_equal(s$coverage$genes_never_tested_in_2, "ONLY_IN_1")
  expect_gt(s$coverage$gene_variant_pairs_only1, 0)
})

test_that("two tables clustered at different windows refuse to be compared", {
  a <- mk_mr(1, 1e6, "G", 1e-9, locus_kb = 100)
  b <- mk_mr(1, 1e6, "G", 1e-9, locus_kb = 1000)
  expect_error(cqtna_stability(a, b), "different locus_kb")
  expect_silent(cqtna_stability(a, b, locus_kb = 1000))
})

# --------------------------------------------------------------------------
# Mechanism 5: illegal p-values and coordinates must not reach the arithmetic.
test_that("out-of-range p-values are rejected with the offending rows named", {
  base <- data.frame(record_id = c("a", "b"), gene = c("G1", "G2"),
                     chr = c("1", "1"), pos = c(1e6, 2e6), p = c(0.5, 0.5),
                     stringsAsFactors = FALSE)
  neg <- base; neg$p[2] <- -0.1
  expect_error(as_cqtna_mr(neg, build = "GRCh38"), "within \\[0, 1\\]")
  expect_error(as_cqtna_mr(neg, build = "GRCh38"), "offending row")

  big <- base; big$p[1] <- 1.4
  expect_error(as_cqtna_mr(big, build = "GRCh38"), "within \\[0, 1\\]")

  inf <- base; inf$p[1] <- Inf
  expect_error(as_cqtna_mr(inf, build = "GRCh38"), "finite")
})

test_that("missing or impossible coordinates are rejected", {
  base <- data.frame(record_id = c("a", "b"), gene = c("G1", "G2"),
                     chr = c("1", "1"), pos = c(1e6, 2e6), p = c(0.5, 0.5),
                     stringsAsFactors = FALSE)
  na <- base; na$pos[1] <- NA
  expect_error(as_cqtna_mr(na, build = "GRCh38"), "finite positive")
  zero <- base; zero$pos[2] <- 0
  expect_error(as_cqtna_mr(zero, build = "GRCh38"), "finite positive")
  blank <- base; blank$gene[1] <- ""
  expect_error(as_cqtna_mr(blank, build = "GRCh38"), "must not be missing or blank")
  blankchr <- base; blankchr$chr[2] <- NA
  expect_error(as_cqtna_mr(blankchr, build = "GRCh38"), "must not be missing or blank")
})

test_that("illegal windows and FDR thresholds are rejected", {
  d <- cqtna_demo("mr")
  expect_error(as_cqtna_mr(d, locus_kb = -1, build = "GRCh38"), "non-negative")
  expect_error(as_cqtna_mr(d, locus_kb = NA, build = "GRCh38"), "finite")
  mr <- as_cqtna_mr(d, build = "GRCh38")
  kn <- suppressWarnings(as_cqtna_known(cqtna_demo("known"), build = "GRCh38"))
  expect_error(cqtna_attribution(mr, kn, fdr = 0), "in \\(0, 1\\)")
  expect_error(cqtna_attribution(mr, kn, fdr = 1), "in \\(0, 1\\)")
  expect_error(cqtna_window_sweep(mr, kn, sweep_kb = c(100, -5)), "non-negative")
})

test_that("chromosome naming styles must agree", {
  mr <- mk_mr(1, 1e6, "G", 1e-9)
  kn <- mk_known("chr1", 1e6)
  expect_error(cqtna_attribution(mr, kn), "chromosome naming differs")
})

# --------------------------------------------------------------------------
# Mechanism 6: factor 0/1 columns. as.numeric() on a factor of "0"/"1" gives
# 1/2, so every "0" would count as TRUE.
test_that("factor flag columns are refused rather than silently miscounted", {
  d <- data.frame(gene = c("A", "B", "C"),
                  instrumentable = factor(c("1", "0", "0")),
                  analysable = c(1, 0, 0), associated = c(1, 0, 0))
  expect_error(cqtna_ladder(d), "is a factor")

  d$instrumentable <- c("1", "0", "0")            # character is fine
  expect_equal(cqtna_ladder(d)$instrumentable, 1L)
  d$instrumentable <- c(TRUE, FALSE, FALSE)       # logical is fine
  expect_equal(cqtna_ladder(d)$instrumentable, 1L)
})

test_that("the three instrument levels must be nested", {
  d <- data.frame(gene = c("A", "B"), instrumentable = c(0, 1),
                  analysable = c(1, 1), associated = c(0, 1))
  expect_error(cqtna_ladder(d), "must be nested")
  expect_error(cqtna_ladder(d), "A")              # names the offending gene
})

# --------------------------------------------------------------------------
# Mechanism 7: degenerate shapes.
test_that("a single locus, zero significant and all significant all behave", {
  kn <- mk_known(1, 1e6)
  one <- mk_mr(1, 1e6, "G", 1e-9)
  a <- cqtna_attribution(one, kn)
  expect_equal(a$background_loci, 1L)
  expect_equal(a$fold, 1)

  none <- mk_mr(c(1, 2), c(1e6, 1e6), c("G1", "G2"), c(0.9, 0.9))
  a0 <- cqtna_attribution(none, kn)
  expect_equal(a0$significant_loci, 0L)
  expect_true(is.na(a0$fold))
  expect_length(a0$known_genes, 0L)

  allsig <- mk_mr(c(1, 2), c(1e6, 1e6), c("G1", "G2"), c(1e-9, 1e-9))
  aa <- cqtna_attribution(allsig, kn)
  expect_equal(aa$significant_loci, 2L)
  expect_equal(aa$fold, 1)          # significant set == background set
})

# --------------------------------------------------------------------------
# Mechanism 8: compartment ratio when the target cell type has zero expression.
test_that("zero target expression is Inf, not missing", {
  x <- data.frame(gene = c("G", "G"), cell_type = c("CD4_T", "Tumour"),
                  mean_expression = c(0, 5))
  e <- cqtna_compartment(x, "G", "CD4_T")
  expect_true(is.infinite(e$ratio_other_over_target))
  expect_equal(e$compartment_flag, "absent from target cell type")

  x2 <- data.frame(gene = c("G", "G"), cell_type = c("CD4_T", "Tumour"),
                   mean_expression = c(0, 0))
  e2 <- cqtna_compartment(x2, "G", "CD4_T")
  expect_true(is.na(e2$ratio_other_over_target))
  expect_equal(e2$compartment_flag, "not evaluable")
})

test_that("a peak distance of zero is kept, not dropped", {
  f <- cqtna_peak_distance(data.frame(gene = "G", eqtl_pos = 1e6, gwas_pos = 1e6))
  expect_equal(f$distance_bp, 0)
})

# --------------------------------------------------------------------------
# Mechanism 9: a pre-built object whose window disagrees with the audit call.
test_that("a pre-built object is re-clustered rather than reported wrongly", {
  mr100 <- as_cqtna_mr(cqtna_demo("mr"), locus_kb = 100, build = "GRCh38",
                       locus_method = "single_linkage")
  expect_message(
    au <- suppressWarnings(cqtna_audit(mr100, cqtna_demo("known"),
                                       locus_kb = 1000, build = "GRCh38")),
    "re-clustering")
  expect_true(au$settings$reclustered)
  expect_equal(au$settings$locus_kb, 1000)
  # the reported background must be the 1 Mb one, not the 100 kb one
  expect_equal(au$A$background_loci, 554L)
})

# --------------------------------------------------------------------------
# Mechanism 10: nothing is silently skipped, and non-empty E/F reach the report.
test_that("every module reports a status and E/F appear in the report", {
  pk <- data.frame(gene = "PARP1", eqtl_pos = 1e6, gwas_pos = 1.05e6)
  au <- suppressWarnings(cqtna_audit(
    cqtna_demo("mr"), cqtna_demo("known"), mismatch = cqtna_demo("mismatch"),
    instruments = cqtna_demo("instruments"), expression = cqtna_demo("expression"),
    peaks = pk, target_cell_type = "CD4_T", build = "GRCh38"))

  expect_setequal(names(au$module_status),
                  c("A", "A_nc", "B", "C", "D", "E", "F", "G", "G_nc", "spans"))
  expect_equal(au$module_status[["C"]], "not run - no input supplied")
  expect_equal(au$module_status[["F"]], "run")
  expect_equal(au$module_status[["G_nc"]], "run")

  txt <- cqtna_report(au)
  expect_match(txt, "## E. Compartment attribution", fixed = TRUE)
  expect_match(txt, "## F. eQTL-to-GWAS peak distance", fixed = TRUE)
  expect_match(txt, "50,000", fixed = TRUE)          # peak distance, printed as an integer
  expect_match(txt, "## C. List stability", fixed = TRUE)
  expect_match(txt, "not run - no input supplied", fixed = TRUE)
  expect_match(txt, "G (mismatched list)", fixed = TRUE)
})

# --------------------------------------------------------------------------
# Mechanism 11: evidence fields stay independent -- a compartment ratio must not
# erase the locus evidence.
test_that("a compartment ratio does not overwrite known-locus status", {
  au <- suppressWarnings(cqtna_audit(
    cqtna_demo("mr"), cqtna_demo("known"), mismatch = cqtna_demo("mismatch"),
    mr_alt = cqtna_demo("mr_alt"), expression = cqtna_demo("expression"),
    target_cell_type = "CD4_T", build = "GRCh38"))
  e <- au$evidence
  expect_true(all(c("known_locus_status", "outcome_stability",
                    "compartment_ratio", "compartment_flag",
                    "manual_diagnostics_completed") %in% names(e)))
  expect_false(any(e$manual_diagnostics_completed))
  # every gene keeps a locus status regardless of any other field
  expect_true(all(e$known_locus_status %in%
                    c("known", "novel", "not significant here")))
  expect_true(all(!is.na(e$known_locus_status)))
})

# --------------------------------------------------------------------------
# Mechanism 12: the negative control is the referee. If the wrong disease's
# known-locus list also enriches, the cell is void -- and the audit must say so
# rather than reporting the matched fold as if it stood.
test_that("a failing mismatched control is flagged, not quietly reported", {
  # 20 background loci across 20 chromosomes; the 5 significant ones sit near
  # BOTH lists, so the wrong disease's list enriches just as much as the right one
  sig_chr <- 1:5
  chr <- rep(1:20, each = 1)
  mr <- mk_mr(chr = chr, pos = rep(1e6, 20), gene = paste0("G", 1:20),
              p = ifelse(chr %in% sig_chr, 1e-9, 0.9))
  right <- mk_known(sig_chr, rep(1e6, 5))
  wrong <- mk_known(sig_chr, rep(1.05e6, 5))
  expect_warning(
    au <- cqtna_audit(as.data.frame(mr), as.data.frame(right),
                      mismatch = as.data.frame(wrong), build = "GRCh38"),
    "mismatched-list control is itself significant")
  expect_true(au$negative_control_failed)
})

# --------------------------------------------------------------------------
# Mechanism 13: chaining. Single-linkage joins neighbours and keeps going, so a
# dense table collapses a whole region into one "locus"; locus-level counting
# then measures the block's width.
test_that("single-linkage chaining is detected and reported", {
  pos <- seq(1e6, 21e6, by = 5e5)          # 500 kb apart -> all one 20 Mb locus
  mr <- mk_mr(chr = rep(1, length(pos)), pos = pos,
              gene = paste0("G", seq_along(pos)),
              p = c(1e-9, rep(0.9, length(pos) - 1)))
  sp <- suppressWarnings(cqtna_locus_spans(mr))
  expect_equal(sp$n_loci, 1L)
  expect_gt(sp$span_quantiles_kb[[4]], 19000)
  expect_equal(sp$n_wider_than_5x_window, 1L)
  expect_warning(cqtna_locus_spans(mr), "chaining has merged distinct regions")
})

# --------------------------------------------------------------------------
# Mechanism 14: the density-matched permutation. Everything about this test is a
# choice, so the tests are about what it refuses to do, not about a p-value.
test_that("the permutation refuses a convention it cannot build a null for", {
  mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38",
                    locus_method = "single_linkage")
  kn <- suppressWarnings(as_cqtna_known(cqtna_demo("known"), build = "GRCh38"))
  # significant_records scores the observation on significant records; a
  # background locus has none, so the null would be a different quantity
  expect_error(
    cqtna_permutation_control(mr, kn, known_from = "significant_records"),
    "should be one of")
})

test_that("an incomplete match returns NA and says why, not a p-value", {
  mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38",
                    locus_method = "single_linkage")
  kn <- suppressWarnings(as_cqtna_known(cqtna_demo("known"), build = "GRCh38"))
  tight <- cqtna_permutation_control(mr, kn, tolerance = 0.1, n_perm = 100, seed = 1)
  expect_lt(tight$matched_fraction, 1)
  expect_true(is.na(tight$empirical_p))
  expect_match(tight$failed_because, "matched background pool")
  expect_gt(length(tight$unmatched_loci), 0)

  loose <- cqtna_permutation_control(mr, kn, tolerance = 1, n_perm = 100, seed = 1)
  expect_equal(loose$matched_fraction, 1)
  expect_false(is.na(loose$empirical_p))
})

test_that("the specification is carried on the result and the seed reproduces it", {
  mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38",
                    locus_method = "single_linkage")
  kn <- suppressWarnings(as_cqtna_known(cqtna_demo("known"), build = "GRCh38"))
  a <- cqtna_permutation_control(mr, kn, tolerance = 1, n_perm = 300, seed = 42)
  b <- cqtna_permutation_control(mr, kn, tolerance = 1, n_perm = 300, seed = 42)
  expect_equal(a$empirical_p, b$empirical_p)
  for (k in c("known_from", "match_on", "tolerance", "n_perm", "seed",
              "fdr", "known_kb", "replace", "min_matched_fraction"))
    expect_false(is.null(a[[k]]), info = k)
  expect_true("n_genes" %in% a$match_on)      # gene density is matched by default
})

test_that("draws that exhaust their pool are discarded, not filled by reuse", {
  mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38",
                    locus_method = "single_linkage")
  kn <- suppressWarnings(as_cqtna_known(cqtna_demo("known"), build = "GRCh38"))
  r <- cqtna_permutation_control(mr, kn, tolerance = 1, n_perm = 200, seed = 7)
  expect_equal(r$n_draws_used + r$n_draws_exhausted, 200L)
  expect_length(r$null_distribution, r$n_draws_used)
})

test_that("the tolerance sweep warns when the verdict moves across it", {
  mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38",
                    locus_method = "single_linkage")
  kn <- suppressWarnings(as_cqtna_known(cqtna_demo("known"), build = "GRCh38"))
  s <- cqtna_permutation_sensitivity(mr, kn, n_perm = 300, seed = 1)
  expect_true(all(c("tolerance", "matched_fraction", "empirical_p", "ok") %in% names(s)))
  # tolerances too tight to match everything must be marked, not silently used
  expect_true(any(!s$ok))
  expect_true(all(is.na(s$empirical_p[!s$ok])))
})


# --------------------------------------------------------------------------
# Mechanism 15: the fixed-centre partition. Non-recursive by construction, so a
# locus cannot grow past the window however dense the data. This is the answer
# to chaining, and the reason it is the default.
test_that("fixed-centre loci are bounded by the window, single-linkage are not", {
  pos <- seq(1e6, 21e6, by = 5e5)          # 500 kb apart across 20 Mb
  d <- data.frame(record_id = seq_along(pos), gene = paste0("G", seq_along(pos)),
                  chr = "1", pos = pos, p = c(1e-9, rep(0.9, length(pos) - 1)),
                  stringsAsFactors = FALSE)

  sl <- as_cqtna_mr(d, locus_kb = 1000, build = "GRCh38",
                    locus_method = "single_linkage")
  expect_equal(length(unique(sl$locus)), 1L)              # everything chains
  expect_gt(max(suppressWarnings(cqtna_locus_spans(sl))$span_quantiles_kb), 19000)

  fc <- as_cqtna_mr(d, locus_kb = 1000, build = "GRCh38")  # default
  sp <- cqtna_locus_spans(fc)
  expect_gt(length(unique(fc$locus)), 1L)
  expect_lte(max(sp$span_quantiles_kb), 1000)             # bounded by the window
  expect_equal(sp$n_wider_than_5x_window, 0L)
})

test_that("the partition rule does not depend on the outcome statistics", {
  pos <- c(1e6, 1.2e6, 1.4e6, 3e6)
  base <- data.frame(record_id = 1:4, gene = paste0("G", 1:4), chr = "1",
                     pos = pos, stringsAsFactors = FALSE)
  a <- as_cqtna_mr(cbind(base, p = c(1e-9, 0.9, 0.9, 0.9)), locus_kb = 1000,
                   build = "GRCh38")
  b <- as_cqtna_mr(cbind(base, p = c(0.9, 0.9, 0.9, 1e-9)), locus_kb = 1000,
                   build = "GRCh38")
  # moving which variant is significant must not move the partition, because the
  # same partition supplies the denominator
  expect_equal(a$locus, b$locus)
})

test_that("block assignment uses the supplied intervals", {
  d <- data.frame(record_id = 1:4, gene = paste0("G", 1:4), chr = c("1","1","2","2"),
                  pos = c(1e6, 2e6, 1e6, 9e6), p = rep(1e-9, 4),
                  stringsAsFactors = FALSE)
  blk <- data.frame(chr = c("1", "2"), start = c(5e5, 5e5), end = c(2.5e6, 2e6))
  mr <- as_cqtna_mr(d, build = "GRCh38", locus_method = "blocks", blocks = blk)
  expect_equal(length(unique(mr$locus)), 3L)   # two in block 1, one in block 2, one outside
  # the variant outside every block becomes its own single-position locus
  expect_true("2:9e+06-9e+06" %in% mr$locus || "2:9000000-9000000" %in% mr$locus)
  expect_error(as_cqtna_mr(d, build = "GRCh38", locus_method = "blocks"),
               "needs a `blocks` table")
})

test_that("two tables partitioned by different rules refuse to be compared", {
  d <- data.frame(record_id = 1:2, gene = c("A", "B"), chr = "1",
                  pos = c(1e6, 2e6), p = c(1e-9, 1e-9), stringsAsFactors = FALSE)
  a <- as_cqtna_mr(d, build = "GRCh38", locus_method = "fixed_centre")
  b <- as_cqtna_mr(d, build = "GRCh38", locus_method = "single_linkage")
  expect_error(cqtna_stability(a, b), "different rules")
})
