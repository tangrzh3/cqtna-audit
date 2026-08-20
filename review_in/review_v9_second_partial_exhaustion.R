library(cqtna)

mk <- function(pos, n, p) {
  data.frame(
    record_id = paste0("r", seq_len(sum(n))),
    gene = unlist(Map(function(k, j) rep(paste0("g", j), k), n, seq_along(n))),
    chr = "1", pos = rep(pos, n), p = rep(p, n),
    stringsAsFactors = FALSE
  )
}

raw <- mk(c(1e6, 4e6, 7e6, 10e6), c(2, 4, 3, 1),
          c(1e-9, 1e-9, 0.8, 0.8))
mr <- as_cqtna_mr(raw, build = "GRCh38", locus_kb = 1000,
                  locus_method = "fixed_centre")
known <- suppressWarnings(as_cqtna_known(
  data.frame(chr = "1", pos = 7e6), build = "GRCh38"
))
r <- cqtna_permutation_control(
  mr, known, match_on = "n_records", tolerance = 0.75,
  n_perm = 10000, seed = 7, min_draw_fraction = 0.95
)
print(unlist(r[c("n_draws_used", "n_draws_exhausted",
                 "effective_draw_fraction", "empirical_p", "fold_vs_null",
                 "min_pool_size", "median_pool_size", "max_pool_size",
                 "failed_because")]))
stopifnot(r$effective_draw_fraction < 0.95,
          is.na(r$empirical_p), is.na(r$fold_vs_null),
          grepl("required minimum", r$failed_because))
