## Step 124 -- the complete grid on one statistical unit.
##
## Executes the frozen decision record, manuscript/PREREG_locus_partition.md:
##   main analysis    non-recursive fixed-centre partition, max span 1000 kb,
##                    attribution convention C1 (any_record)
##   sensitivity      the same partition at 500 kb
##   legacy           1 Mb single-linkage with the published convention (C2),
##                    reported for reproduction only, carrying no inference
##
## Every cell is recomputed from a per-record table, so one grid uses one
## definition of its unit. The two cells whose per-record tables were never saved
## -- RA on both resources, and whole blood against HCC -- were re-derived by
## step108 and step94d, which now persist 123b and 123c.
##
##   Rscript step124_full_grid.R
## Output: 123d_fixed_anchor_full_grid.tsv

suppressMessages(library(cqtna))
MR <- "D:/R_ex/MR"
setwd(MR)

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

mel <- read.csv("landi2020_known_loci_grch38.csv"); mel$source <- "Landi2020"
hcc <- read.csv("84a_hcc_known_loci_grch38.csv");   hcc$source <- "step84"
K_mel <- as_cqtna_known(mel, build = "GRCh38")
K_hcc <- as_cqtna_known(hcc, build = "GRCh38")

## Okada RA lead SNPs carry no GRCh38 coordinates of their own; step108 places
## them while streaming the RA sumstats. The placed positions come back inside
## the per-record table's `known` marking, so RA is scored here against the
## positions step108 resolved, read from its own output.
ra_rec <- read.delim(gzfile("123b_ra_records.tsv.gz"), stringsAsFactors = FALSE)
eh_rec <- read.delim(gzfile("123c_eqtlgen_hcc_records.tsv.gz"), stringsAsFactors = FALSE)
ra_known <- read.delim("123b_ra_known_positions.tsv", stringsAsFactors = FALSE)
K_ra <- as_cqtna_known(within(ra_known, source <- "Okada2014"), build = "GRCh38")

cd4 <- read.delim("13_meta_locus_annotation.tsv")
cd4$chr <- sub(":.*", "", cd4$SNP); cd4$pos <- as.numeric(sub(".*:", "", cd4$SNP))
eq  <- read.delim("92c_locus_annotated.tsv")
hh  <- read.delim("85a_HCC_high_annotated.tsv")
hl  <- read.delim("85a_HCC_low_annotated.tsv")

pick <- function(d, cell) d[d$cell == cell, , drop = FALSE]

cells <- list(
  list(name = "melanoma x Soskic_CD4", role = "main",
       d = std(cd4$SYMBOL, cd4$chr, cd4$pos, cd4$pval), m = K_mel, x = K_hcc),
  list(name = "melanoma x eQTLGen_blood", role = "main",
       d = std(eq$symbol, eq$chr, eq$pos, eq$p_mr), m = K_mel, x = K_hcc),
  list(name = "HCC_high x Soskic_CD4", role = "main",
       d = std(hh$gene_id, hh$chr, hh$pos, hh$p_mr), m = K_hcc, x = K_mel),
  list(name = "HCC_high x eQTLGen_blood", role = "main",
       d = with(pick(eh_rec, "HCC_high"), std(gene, chr, pos, p)), m = K_hcc, x = K_mel),
  list(name = "RA x Soskic_CD4", role = "main",
       d = with(pick(ra_rec, "N1 Soskic x RA"), std(gene, chr, pos, p)),
       m = K_ra, x = K_mel),
  list(name = "RA x eQTLGen_blood", role = "main",
       d = with(pick(ra_rec, "N2 eQTLGen x RA"), std(gene, chr, pos, p)),
       m = K_ra, x = K_mel),
  list(name = "HCC_low x Soskic_CD4", role = "power sensitivity",
       d = std(hl$gene_id, hl$chr, hl$pos, hl$p_mr), m = K_hcc, x = K_mel),
  list(name = "HCC_low x eQTLGen_blood", role = "power sensitivity",
       d = with(pick(eh_rec, "HCC_low"), std(gene, chr, pos, p)), m = K_hcc, x = K_mel))

runs <- list(list(method = "fixed_centre",   kb = 1000, conv = "any_record",
                  label = "main"),
             list(method = "fixed_centre",   kb = 500,  conv = "any_record",
                  label = "sensitivity"),
             list(method = "single_linkage", kb = 1000, conv = "significant_records",
                  label = "legacy"))

rows <- list()
for (cl in cells) for (rn in runs) {
  if (!nrow(cl$d)) next
  mr <- as_cqtna_mr(cl$d, locus_kb = rn$kb, build = "GRCh38",
                    locus_method = rn$method)
  sp <- suppressWarnings(cqtna_locus_spans(mr))
  a <- cqtna_attribution(mr, cl$m, known_from = rn$conv)
  x <- cqtna_attribution(mr, cl$x, known_from = rn$conv)
  rows[[length(rows) + 1L]] <- data.frame(
    cell = cl$name, role = cl$role, analysis = rn$label,
    partition = rn$method, locus_kb = rn$kb, known_from = rn$conv,
    n_records = nrow(cl$d),
    bg_loci = a$background_loci, bg_known = a$background_known,
    sig_loci = a$significant_loci, sig_known = a$significant_known,
    pct_known = round(a$pct_known, 1),
    fold = round(a$fold, 2), fisher_p = a$fisher_p_one_sided,
    mismatch_fold = round(x$fold, 2), mismatch_p = x$fisher_p_one_sided,
    control = if (is.finite(x$fisher_p_one_sided) && x$fisher_p_one_sided < 0.05 &&
                  isTRUE(x$fold > 1)) "FAILED" else "clean",
    max_sig_span_kb = round(max(sp$significant_span_kb)),
    n_sig_over_window = sum(sp$significant_span_kb > rn$kb),
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
write.table(out, "123d_fixed_anchor_full_grid.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

show <- function(lbl) {
  s <- out[out$analysis == lbl, ]
  cat("\n", strrep("=", 118), "\n", toupper(lbl), " -- ", s$partition[1], " ",
      s$locus_kb[1], " kb, ", s$known_from[1], "\n", sep = "")
  cat(sprintf("%-28s %-18s %9s %9s %8s %9s   %8s %10s %8s\n", "cell", "role",
              "bg known", "sig known", "fold", "P", "mism.", "mism. P", "control"))
  for (i in seq_len(nrow(s))) with(s[i, ], cat(sprintf(
    "%-28s %-18s %4d/%-4d %4d/%-4d %8.2f %9.3g   %8.2f %10.3g %8s\n",
    cell, role, bg_known, bg_loci, sig_known, sig_loci, fold, fisher_p,
    mismatch_fold, mismatch_p, control)))
  m <- s[s$role == "main", ]
  cat(sprintf("  -> %d/%d main cells with fold > 1; %d reach P < 0.05; controls: %s\n",
              sum(m$fold > 1, na.rm = TRUE), nrow(m),
              sum(m$fisher_p < 0.05, na.rm = TRUE),
              paste(names(table(m$control)), table(m$control), collapse = ", ")))
  cat(sprintf("  -> widest significant locus %s kb; loci over the window: %d\n",
              format(max(s$max_sig_span_kb), big.mark = ","),
              sum(s$n_sig_over_window)))
}
for (l in c("main", "sensitivity", "legacy")) show(l)
cat("\nwrote 123d_fixed_anchor_full_grid.tsv\n")
