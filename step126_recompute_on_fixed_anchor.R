## Step 126 -- everything else the manuscript quotes, recomputed on one unit.
##
## step124 moved the six-cell grid onto the frozen partition (fixed-centre,
## 1000 kb, C1). The manuscript also quotes attribution numbers that live outside
## that grid, and every one of them is a locus-level statistic, so every one of
## them is partition-dependent:
##
##   the R13 transfer test                9.6-fold   (step99, now 99d)
##   the three extra CD4 resources        3.36 / 6.58-fold  (step116, 116b)
##   the matched-background permutation   4.14-4.48-fold    (step92, 92e)
##   the density-matched permutation      3.44-fold  (step103, 103a)
##
## Leaving those on the single-linkage unit while the grid moves is exactly the
## mixing of two statistical units that PREREG_locus_partition.md section 4
## forbids, so they are recomputed here from the per-record tables.
##
## The permutations follow manuscript/PREREG_permutation_estimand.md, frozen
## 2026-08-19: convention any_record, matching on n_records + span + n_genes,
## no chromosome stratification, significant loci excluded from the pool,
## sampling without replacement, 100% matched coverage required, tolerance 1.00
## with the whole tolerance scan reported alongside.
##
##   Rscript step126_recompute_on_fixed_anchor.R
## Output: 126a_offgrid_attribution.tsv   the non-grid attribution cells
##         126b_permutation_scan.tsv      the full tolerance scan, per cell
##         126c_permutation_primary.tsv   tolerance 1.00, the quotable row

MR <- "D:/R_ex/MR"
setwd(MR)
suppressMessages(library(cqtna))

KB   <- 1000            # the frozen main partition
CONV <- "any_record"    # C1
NPERM <- 10000
SEED  <- 1
TOLS  <- c(0.1, 0.25, 0.5, 1, 2)

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
as_mr <- function(d) as_cqtna_mr(d, locus_kb = KB, build = "GRCh38",
                                 locus_method = "fixed_centre")

mel <- read.csv("landi2020_known_loci_grch38.csv"); mel$source <- "Landi2020"
hcc <- read.csv("84a_hcc_known_loci_grch38.csv");   hcc$source <- "step84"
K_mel <- as_cqtna_known(mel, build = "GRCh38")
K_hcc <- as_cqtna_known(hcc, build = "GRCh38")

## ------------------------------------------------ A. the non-grid attributions
r13 <- read.delim(gzfile("99d_r13_records.tsv.gz"), stringsAsFactors = FALSE)
pick <- function(d, cell) d[d$cell == cell, , drop = FALSE]

jobs <- list(
  ## the R13 transfer test and its R12 comparator, computed by the same path
  list(name = "C1 Soskic x melanoma R12", src = "99d", cell = "C1 Soskic x melanoma R12",
       m = K_mel, x = K_hcc, note = "R13 transfer comparator"),
  list(name = "C1 Soskic x melanoma R13", src = "99d", cell = "C1 Soskic x melanoma R13",
       m = K_mel, x = K_hcc, note = "R13 transfer, main cell"),
  list(name = "C2 eQTLGen x melanoma R13", src = "99d", cell = "C2 eQTLGen x melanoma R13",
       m = K_mel, x = K_hcc, note = "R13 transfer"),
  list(name = "C3 Soskic x HCC R13", src = "99d", cell = "C3 Soskic x HCC R13",
       m = K_hcc, x = K_mel, note = "R13 transfer"),
  list(name = "C4 eQTLGen x HCC R13", src = "99d", cell = "C4 eQTLGen x HCC R13",
       m = K_hcc, x = K_mel, note = "R13 transfer"))

## the three additional CD4 resources (S35); 116b holds their per-record tables
for (res in c("Nathan_2022", "Randolph_2021", "Schmiedel_2018"))
  jobs[[length(jobs) + 1L]] <- list(name = paste0(res, " x melanoma"),
                                    src = "116b", cell = res,
                                    m = K_mel, x = K_hcc,
                                    note = "extra CD4 resource (S35)")

rows <- list()
for (j in jobs) {
  d <- if (j$src == "99d") {
    with(pick(r13, j$cell), std(gene, chr, pos, p))
  } else {
    t <- read.delim(sprintf("116b_mr_%s.tsv", j$cell), stringsAsFactors = FALSE)
    std(t$gene_id, t$chr, t$pos, t$p_mr)
  }
  if (!nrow(d)) next
  mr <- as_mr(d)
  a <- cqtna_attribution(mr, j$m, known_from = CONV)
  x <- cqtna_attribution(mr, j$x, known_from = CONV)
  rows[[length(rows) + 1L]] <- data.frame(
    cell = j$name, note = j$note, partition = "fixed_centre", locus_kb = KB,
    known_from = CONV, n_records = nrow(d),
    bg_loci = a$background_loci, bg_known = a$background_known,
    sig_loci = a$significant_loci, sig_known = a$significant_known,
    pct_known = round(a$pct_known, 1),
    fold = round(a$fold, 2), fisher_p = a$fisher_p_one_sided,
    mismatch_fold = round(x$fold, 2), mismatch_p = x$fisher_p_one_sided,
    control = if (is.finite(x$fisher_p_one_sided) && x$fisher_p_one_sided < 0.05 &&
                  isTRUE(x$fold > 1)) "FAILED" else "clean",
    stringsAsFactors = FALSE)
}
off <- do.call(rbind, rows)
write.table(off, "126a_offgrid_attribution.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("\n", strrep("=", 112), "\n",
    "A. non-grid attribution cells, recomputed on fixed_centre ", KB, " kb / ",
    CONV, "\n", sep = "")
cat(sprintf("%-30s %9s %9s %8s %10s %8s %9s %8s\n", "cell", "bg known",
            "sig known", "fold", "P", "mism.", "mism. P", "control"))
for (i in seq_len(nrow(off))) with(off[i, ], cat(sprintf(
  "%-30s %4d/%-4d %4d/%-4d %8.2f %10.3g %8.2f %9.3g %8s\n",
  cell, bg_known, bg_loci, sig_known, sig_loci, fold, fisher_p,
  mismatch_fold, mismatch_p, control)))

## ------------------------------------------------------- B. the permutations
## Only cells whose own attribution is quotable are worth permuting, so the two
## RA cells are included: the mismatched control is the weaker of the two
## controls and RA is where it is least clean, which is precisely where a
## density-matched null is worth having.
ra_rec <- read.delim(gzfile("123b_ra_records.tsv.gz"), stringsAsFactors = FALSE)
eh_rec <- read.delim(gzfile("123c_eqtlgen_hcc_records.tsv.gz"), stringsAsFactors = FALSE)
cd4 <- read.delim("13_meta_locus_annotation.tsv")
cd4$chr <- sub(":.*", "", cd4$SNP); cd4$pos <- as.numeric(sub(".*:", "", cd4$SNP))
eq  <- read.delim("92c_locus_annotated.tsv")
hh  <- read.delim("85a_HCC_high_annotated.tsv")
hl  <- read.delim("85a_HCC_low_annotated.tsv")

ra_known <- read.delim("123b_ra_known_positions.tsv", stringsAsFactors = FALSE)
K_ra <- as_cqtna_known(within(ra_known, source <- "Okada2014"), build = "GRCh38")

perm_cells <- list(
  list(name = "melanoma x Soskic_CD4",
       d = std(cd4$SYMBOL, cd4$chr, cd4$pos, cd4$pval), m = K_mel),
  list(name = "melanoma x eQTLGen_blood",
       d = std(eq$symbol, eq$chr, eq$pos, eq$p_mr), m = K_mel),
  list(name = "HCC_high x Soskic_CD4",
       d = std(hh$gene_id, hh$chr, hh$pos, hh$p_mr), m = K_hcc),
  list(name = "HCC_high x eQTLGen_blood",
       d = with(pick(eh_rec, "HCC_high"), std(gene, chr, pos, p)), m = K_hcc),
  list(name = "HCC_low x Soskic_CD4",
       d = std(hl$gene_id, hl$chr, hl$pos, hl$p_mr), m = K_hcc),
  list(name = "HCC_low x eQTLGen_blood",
       d = with(pick(eh_rec, "HCC_low"), std(gene, chr, pos, p)), m = K_hcc),
  ## RA scored over the whole genome, i.e. the cell S33 registered as primary --
  ## not the MHC-excluded one. RA x eQTLGen is void on its mismatched control and
  ## stays void whatever the permutation says; it is permuted anyway, because
  ## "the density-matched null also fails" and "we never looked" are different
  ## statements and only the first can be written down.
  list(name = "RA x Soskic_CD4",
       d = with(pick(ra_rec, "N1 Soskic x RA"), std(gene, chr, pos, p)), m = K_ra),
  list(name = "RA x eQTLGen_blood (VOID on mismatch)",
       d = with(pick(ra_rec, "N2 eQTLGen x RA"), std(gene, chr, pos, p)), m = K_ra))

scan <- list(); prim <- list()
for (cl in perm_cells) {
  mr <- as_mr(cl$d)
  s <- cqtna_permutation_sensitivity(mr, cl$m, known_from = CONV,
                                     tolerances = TOLS, n_perm = NPERM,
                                     seed = SEED)
  s$cell <- cl$name
  scan[[length(scan) + 1L]] <- s[, c("cell", setdiff(names(s), "cell"))]
  r <- cqtna_permutation_control(mr, cl$m, known_from = CONV, tolerance = 1,
                                 n_perm = NPERM, seed = SEED)
  prim[[length(prim) + 1L]] <- data.frame(
    cell = cl$name, tolerance = 1, n_perm = NPERM, seed = SEED,
    known_from = r$known_from, match_on = paste(r$match_on, collapse = "+"),
    matched_fraction = r$matched_fraction,
    observed_known = r$observed_known,
    fold_vs_null = if (is.null(r$fold_vs_null)) NA_real_ else r$fold_vs_null,
    empirical_p = r$empirical_p,
    failed_because = if (is.na(r$failed_because)) "" else r$failed_because,
    stringsAsFactors = FALSE)
}
sc <- do.call(rbind, scan); pr <- do.call(rbind, prim)
write.table(sc, "126b_permutation_scan.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pr, "126c_permutation_primary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n", strrep("=", 112), "\n",
    "B. density-matched permutation, frozen estimand, ", NPERM, " permutations\n",
    "   tolerance scan -- a row is quotable only where matched coverage is 100%\n",
    sep = "")
for (nm in unique(sc$cell)) {
  s <- sc[sc$cell == nm, ]
  cat("\n  ", nm, "\n", sep = "")
  cat(sprintf("    %9s %10s %9s %11s %8s\n", "tolerance", "coverage",
              "fold", "empirical P", "usable"))
  for (i in seq_len(nrow(s))) with(s[i, ], cat(sprintf(
    "    %9.2f %9.1f%% %9s %11s %8s\n", tolerance, 100 * matched_fraction,
    ifelse(is.na(fold_vs_null), "-", sprintf("%.2f", fold_vs_null)),
    ifelse(is.na(empirical_p), "-", sprintf("%.4f", empirical_p)),
    ifelse(ok, "yes", "NO"))))
}

cat("\n", strrep("=", 112), "\n",
    "   tolerance 1.00 -- the row the manuscript quotes\n", sep = "")
cat(sprintf("%-30s %9s %10s %8s %12s   %s\n", "cell", "coverage", "obs known",
            "fold", "empirical P", "failed because"))
for (i in seq_len(nrow(pr))) with(pr[i, ], cat(sprintf(
  "%-30s %8.1f%% %10d %8s %12s   %s\n", cell, 100 * matched_fraction,
  observed_known, ifelse(is.na(fold_vs_null), "-", sprintf("%.2f", fold_vs_null)),
  ifelse(is.na(empirical_p), "-", sprintf("%.4f", empirical_p)), failed_because)))

cat("\nwrote 126a / 126b / 126c\n")
