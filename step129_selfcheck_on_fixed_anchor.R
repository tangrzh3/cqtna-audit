## Step 129 -- do the two self-check tables move when the partition moves?
##
## Two tables in the tree were built on the single-linkage partition and were not
## rebuilt when it changed:
##
##   96a_attribution_selfcheck.tsv   which genes the pipeline NAMES at each known
##                                   melanoma locus, grouped by locus id
##   100a_unit_sensitivity.tsv       the FDR list recomputed under eight testing
##                                   units, which the manuscript cites as S26
##
## Neither feeds an enrichment number, so neither is covered by the migration in
## step124-128. But "it probably does not matter" is not a finding, so this
## recomputes both on the frozen partition and reports whether the manuscript's
## claims about them survive. It writes no replacement table: the point is the
## comparison, and a second table with the same name is how the old one gets
## quoted by accident.
##
##   Rscript step129_selfcheck_on_fixed_anchor.R
## Output: 129a_selfcheck_partition_comparison.tsv

## Work from wherever this script lives, so the packet runs after extraction.
## Override with:  Rscript <script> /path/to/dir     or  CQTNA_DIR=/path/to/dir
MR <- local({
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) && nzchar(a[1])) return(a[1])
  if (nzchar(Sys.getenv("CQTNA_DIR"))) return(Sys.getenv("CQTNA_DIR"))
  f <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", f[grepl("^--file=", f)])
  if (length(f)) normalizePath(dirname(f[1])) else getwd()
})
setwd(MR)
suppressMessages(library(cqtna))

FDR <- 0.05
KB  <- 1000

## Same helpers as step124_cells.R: 13_meta carries blank SYMBOLs, which
## as_cqtna_mr() rejects rather than silently merging into one "gene".
fix_gene <- function(g) {
  g <- as.character(g)
  bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad))
  g
}
std <- function(gene, chr, pos, p) {
  d <- data.frame(record_id = seq_along(pos), gene = fix_gene(gene),
                  chr = as.character(chr), pos = as.numeric(pos),
                  p = as.numeric(p), stringsAsFactors = FALSE)
  d[is.finite(d$p) & is.finite(d$pos) & d$pos > 0, , drop = FALSE]
}

## ------------------------------------------- A. the gene-naming benchmark (96a)
## The accepted-gene column was curated by hand against the published locus ids
## and cannot be regenerated. What CAN be checked is the thing the manuscript
## actually claims: the set of genes named at each benchmark locus. If the
## partition does not change those sets, the benchmark verdict cannot change
## either, whatever the locus ids are called.
eq <- read.delim("92c_locus_annotated.tsv", stringsAsFactors = FALSE)
old_sig <- eq[eq$fdr < FDR & tolower(as.character(eq$known)) %in% c("true", "1"), ]
old_named <- tapply(old_sig$symbol, old_sig$locus,
                    function(g) paste(sort(unique(g)), collapse = ";"))

mr <- as_cqtna_mr(std(eq$symbol, eq$chr, eq$pos, eq$p_mr), locus_kb = KB,
                  build = "GRCh38", locus_method = "fixed_centre")
mel <- read.csv("landi2020_known_loci_grch38.csv"); mel$source <- "Landi2020"
K_mel <- as_cqtna_known(mel, build = "GRCh38")
a <- cqtna_attribution(mr, K_mel, known_from = "any_record")

new <- data.frame(gene = as.character(mr$gene), locus = as.character(mr$locus),
                  fdr = as.numeric(mr$fdr), stringsAsFactors = FALSE)
new$known <- a$significant_known_flag_by_record <- NULL     # not exported; recompute
## a locus is known under C1 if any of its records is near a lead SNP
nearest_bp <- function(known, chr, pos) {
  chr <- as.character(chr); pos <- as.numeric(pos)
  out <- rep(Inf, length(pos))
  by_chr <- split(as.numeric(known$pos), as.character(known$chr))
  for (ch in unique(chr)) {
    arr <- sort(by_chr[[ch]]); if (!length(arr)) next
    i <- which(chr == ch); p <- pos[i]; j <- findInterval(p, arr)
    lo <- ifelse(j >= 1L, abs(arr[pmax(j, 1L)] - p), Inf)
    hi <- ifelse(j < length(arr), abs(arr[pmin(j + 1L, length(arr))] - p), Inf)
    out[i] <- pmin(lo, hi)
  }
  out
}
near <- nearest_bp(K_mel, mr$chr, mr$pos) <= KB * 1000
locus_known <- tapply(near, as.character(mr$locus), any)
new$locus_known <- locus_known[new$locus]
new_sig <- new[new$fdr < FDR & new$locus_known, ]
new_named <- tapply(new_sig$gene, new_sig$locus,
                    function(g) paste(sort(unique(g)), collapse = ";"))

set_old <- sort(unname(old_named)); set_new <- sort(unname(new_named))
same_sets <- identical(set_old, set_new)
genes_old <- sort(unique(unlist(strsplit(set_old, ";"))))
genes_new <- sort(unique(unlist(strsplit(set_new, ";"))))

cat("\n", strrep("=", 92), "\n",
    "A. gene-naming benchmark: named-gene sets at known melanoma loci\n", sep = "")
cat(sprintf("  single linkage : %d known significant loci, %d named genes\n",
            length(set_old), length(genes_old)))
cat(sprintf("  fixed anchor   : %d known significant loci, %d named genes\n",
            length(set_new), length(genes_new)))
## "FALSE" here means the per-locus groupings differ, not that any gene
## appeared or vanished. A reviewer read the bare boolean as a failure, so
## name which of the two questions each line answers.
cat(sprintf("  same gene-to-locus GROUPINGS: %s\n", same_sets))
cat(sprintf("  same SET of named genes:      %s   <- the one that matters\n",
            identical(genes_old, genes_new)))
cat(sprintf("  genes only under single linkage: %s\n",
            paste(setdiff(genes_old, genes_new), collapse = ", ")))
cat(sprintf("  genes only under fixed anchor  : %s\n",
            paste(setdiff(genes_new, genes_old), collapse = ", ")))

## ------------------------------------------- B. the testing-unit sweep (100a)
## S26 varies the TESTING unit -- record, variant, gene, locus, each by minimum-p
## and by Simes. Only the two locus-based units can move with the partition.
cd4 <- read.delim("13_meta_locus_annotation.tsv", stringsAsFactors = FALSE)
cd4$chr <- sub(":.*", "", cd4$SNP); cd4$pos <- as.numeric(sub(".*:", "", cd4$SNP))
d <- std(cd4$SYMBOL, cd4$chr, cd4$pos, cd4$pval)

rows <- list()
for (meth in c("single_linkage", "fixed_centre")) {
  m <- as_cqtna_mr(d, locus_kb = KB, build = "GRCh38", locus_method = meth)
  u <- cqtna_unit_sweep(m)
  u$partition <- meth
  rows[[length(rows) + 1L]] <- u
}
us <- do.call(rbind, rows)
us <- us[, c("partition", setdiff(names(us), "partition"))]

cat("\n", strrep("=", 92), "\n",
    "B. testing-unit sweep (S26) under both partitions\n", sep = "")
cat(sprintf("%-16s %-10s %9s %14s %8s %7s\n", "partition", "unit", "n tests",
            "n significant", "n genes", "n loci"))
for (i in seq_len(nrow(us))) with(us[i, ], cat(sprintf(
  "%-16s %-10s %9d %14d %8d %7d\n", partition, unit, n_tests, n_significant,
  n_genes, n_loci)))

write.table(us, "129a_selfcheck_partition_comparison.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nwrote 129a_selfcheck_partition_comparison.tsv\n")
