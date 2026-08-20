args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: review_v9_second_actual_order.R <reproduce-dir>")
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)

suppressMessages(library(cqtna))
cd4 <- read.delim(file.path(root, "13_meta_locus_annotation.tsv"),
                  stringsAsFactors = FALSE)
raw <- data.frame(
  record_id = seq_len(nrow(cd4)),
  gene = as.character(cd4$SYMBOL),
  chr = sub(":.*", "", cd4$SNP),
  pos = as.numeric(sub(".*:", "", cd4$SNP)),
  p = as.numeric(cd4$pval),
  stringsAsFactors = FALSE
)
bad <- is.na(raw$gene) | !nzchar(trimws(raw$gene))
raw$gene[bad] <- paste0("unnamed_", which(bad))
raw <- raw[is.finite(raw$p) & is.finite(raw$pos) & raw$pos > 0, ]

known_raw <- read.csv(file.path(root, "landi2020_known_loci_grch38.csv"),
                      stringsAsFactors = FALSE)
known_raw$source <- "Landi2020"
known <- as_cqtna_known(known_raw, build = "GRCh38")

run <- function(d) {
  mr <- as_cqtna_mr(d, locus_kb = 1000, build = "GRCh38",
                    locus_method = "fixed_centre")
  cqtna_permutation_control(
    mr, known, known_from = "any_record", tolerance = 1,
    n_perm = 10000, seed = 1
  )
}

set.seed(99)
ix <- sample.int(nrow(raw))
r <- list(original = run(raw), reversed = run(raw[nrow(raw):1L, ]),
          shuffled = run(raw[ix, ]))
out <- do.call(rbind, lapply(names(r), function(nm) {
  x <- r[[nm]]
  data.frame(
    order = nm, fold = x$fold_vs_null, empirical_p = x$empirical_p,
    used = x$n_draws_used, exhausted = x$n_draws_exhausted,
    effective_fraction = x$effective_draw_fraction,
    min_pool = x$min_pool_size, median_pool = x$median_pool_size,
    max_pool = x$max_pool_size, rng = x$rng_kind,
    r_version = x$r_version, collate = x$collate
  )
}))
print(out, row.names = FALSE)
stopifnot(length(unique(out$fold)) == 1L,
          length(unique(out$empirical_p)) == 1L)
