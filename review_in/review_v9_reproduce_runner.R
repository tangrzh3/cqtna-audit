args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: review_v9_reproduce_runner.R <script.R> <isolated_workdir>")
}
script <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
workdir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
code <- readLines(script, warn = FALSE, encoding = "UTF-8")
needle <- 'MR <- "D:/R_ex/MR"'
replacement <- sprintf('MR <- "%s"', workdir)
if (sum(code == needle) != 1L) {
  stop("expected exactly one hard-coded MR assignment in ", basename(script))
}
code[code == needle] <- replacement
eval(parse(text = code, srcfile = script), envir = new.env(parent = globalenv()))
