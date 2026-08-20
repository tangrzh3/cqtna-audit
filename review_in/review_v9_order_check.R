library(cqtna)

mk <- function(pos, n, p) {
  data.frame(
    record_id = paste0("r", seq_len(sum(n))),
    gene = unlist(Map(function(k, j) rep(paste0("g", j), k), n, seq_along(n))),
    chr = "1",
    pos = rep(pos, n),
    p = rep(p, n)
  )
}

d <- mk(
  c(1e6, 4e6, 7e6, 10e6, 13e6),
  c(2, 4, 3, 1, 6),
  c(1e-9, 1e-9, 0.8, 0.8, 0.8)
)
mr <- as_cqtna_mr(
  d, build = "GRCh38", locus_kb = 1000, locus_method = "fixed_centre"
)
kn <- suppressWarnings(as_cqtna_known(
  data.frame(chr = "1", pos = c(1e6, 10e6)), build = "GRCh38"
))

run <- function(x) cqtna_permutation_control(
  x, kn, match_on = "n_records", tolerance = 0.75,
  n_perm = 10000, seed = 7
)

a <- run(mr)
z <- run(mr[nrow(mr):1, ])
print(c(
  original_used = a$n_draws_used,
  original_exhausted = a$n_draws_exhausted,
  original_p = a$empirical_p,
  reversed_used = z$n_draws_used,
  reversed_exhausted = z$n_draws_exhausted,
  reversed_p = z$empirical_p
))
