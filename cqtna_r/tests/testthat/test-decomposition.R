test_that("the decomposition reproduces the paper's melanoma split", {
  # 全部 8 个位点 -> 4.96 倍、P = 0.0048；剔除自身达 5e-8 的 3 个 -> 1.98 倍、P = 0.413
  skip_if_not(file.exists(file.path(Sys.getenv("CQTNA_DIR"), "13_meta_locus_annotation.tsv")),
              "paper inputs not present")
  d <- utils::read.delim(file.path(Sys.getenv("CQTNA_DIR"),
                                   "13_meta_locus_annotation.tsv"))
  d$chr <- sub(":.*", "", d$SNP); d$pos <- as.numeric(sub(".*:", "", d$SNP))
  g <- as.character(d$SYMBOL); bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("u", which(bad))
  mr <- as_cqtna_mr(data.frame(record_id = seq_len(nrow(d)), gene = g,
                               chr = d$chr, pos = d$pos, p = d$pval),
                    locus_kb = 1000, build = "GRCh38",
                    locus_method = "fixed_centre")
  kn <- utils::read.csv(file.path(Sys.getenv("CQTNA_DIR"),
                                  "landi2020_known_loci_grch38.csv"))
  kn$source <- "Landi2020"
  kn <- as_cqtna_known(kn, build = "GRCh38")
  out <- cqtna_decomposition(mr, kn, outcome_p = mr$p, known_from = "any_record")
  expect_equal(out$n_loci, c(8L, 5L, 5L))
  expect_equal(out$n_known, c(4L, 1L, 1L))
  expect_equal(round(out$fold, 2), c(4.96, 1.98, 2.07))
  expect_equal(round(out$fisher_p, 3), c(0.005, 0.413, 0.399))
  expect_equal(attr(out, "n_already_genome_wide"), 3L)
  # 背景里达 5e-8 的位点数；条件化把背景从 66/655 变成 63/652
  expect_equal(attr(out, "n_background_genome_wide"), 3L)
  expect_equal(out$background_loci, c(655L, 655L, 652L))
  expect_equal(out$background_known, c(66L, 66L, 63L))
})

test_that("A has a fixed ceiling where fold's ceiling moves with list density", {
  # 我第一版在这里断言"背景不同但 A 相同"，那是错的：A 是机会校正量，不是密度不变量。
  # 真正成立的不变性是：显著位点全部为 known 时 A = 1，与背景密度无关，
  # 而 fold 在同一情形下等于 1/p_bg，随密度剧烈变化。
  stat <- function(n_bg, n_known_bg, n_sig, n_known_sig) {
    p_bg <- n_known_bg / n_bg; p_sig <- n_known_sig / n_sig
    c(fold = p_sig / p_bg, A = (p_sig - p_bg) / (1 - p_bg))
  }
  sparse <- stat(1000, 100, 10, 10)   # p_bg 0.10, 全部 known
  dense  <- stat(1000, 400, 10, 10)   # p_bg 0.40, 全部 known
  expect_equal(unname(sparse[["A"]]), 1)
  expect_equal(unname(dense[["A"]]), 1)
  expect_equal(unname(sparse[["fold"]]), 10)
  expect_equal(unname(dense[["fold"]]), 2.5)
  expect_gt(sparse[["fold"]] / dense[["fold"]], 3.9)   # fold 差 4 倍，A 一样
})

test_that("outcome_p must be one-to-one with mr", {
  mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38")
  kn <- as_cqtna_known(cqtna_demo("known"), build = "GRCh38")
  expect_error(cqtna_decomposition(mr, kn, outcome_p = mr$p[1:3]),
               "one-to-one")
  expect_error(cqtna_decomposition(mr, kn, outcome_p = mr$p, gws = 0),
               "strictly between 0 and 1")
})

test_that("a locus whose outcome p is missing is not silently called novel", {
  mr <- as_cqtna_mr(cqtna_demo("mr"), build = "GRCh38")
  kn <- as_cqtna_known(cqtna_demo("known"), build = "GRCh38")
  op <- mr$p; op[mr$fdr < 0.05] <- NA_real_
  out <- cqtna_decomposition(mr, kn, outcome_p = op)
  # 全部显著位点的 outcome p 缺失 -> 一个都不能判为"已达全基因组显著"
  expect_equal(attr(out, "n_already_genome_wide"), 0L)
  expect_true(all(out$n_outcome_p_missing[2] == out$n_loci[2]))
})
