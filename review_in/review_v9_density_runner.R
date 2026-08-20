args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: review_v9_density_runner.R <step126.R> <workdir> <reverse>")
}

script <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
workdir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
reverse_rows <- identical(tolower(args[[3]]), "true")
code <- readLines(script, warn = FALSE, encoding = "UTF-8")

needle <- 'MR <- "D:/R_ex/MR"'
replacement <- sprintf('MR <- "%s"', workdir)
if (sum(code == needle) != 1L) stop("hard-coded MR assignment not found")
code[code == needle] <- replacement

begin_a <- grep("## ------------------------------------------------ A.", code, fixed = TRUE)[1]
begin_b <- grep("## ------------------------------------------------------- B.", code, fixed = TRUE)[1]
if (!is.finite(begin_a) || !is.finite(begin_b) || begin_b <= begin_a) {
  stop("section markers not found")
}

kept <- c(
  code[seq_len(begin_a - 1L)],
  'pick <- function(d, cell) d[d$cell == cell, , drop = FALSE]',
  code[begin_b:length(code)]
)
target <- "  mr <- as_mr(cl$d)"
if (sum(kept == target) != 1L) stop("permutation MR construction not found")
target_index <- which(kept == target)
kept <- append(
  kept,
  "  if (reverse_rows) mr <- mr[nrow(mr):1L, , drop = FALSE]",
  after = target_index
)

eval(parse(text = kept, srcfile = script), envir = new.env(parent = globalenv()))
