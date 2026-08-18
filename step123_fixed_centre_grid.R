## Step 123 -- re-run the attribution grid on the fixed-centre partition.
##
## Decision (2026-08-18): the main analysis uses a NON-RECURSIVE fixed-centre
## window; LD blocks are the key sensitivity analysis; single-linkage is retained
## only to reproduce published numbers.
##
## Why: single-linkage chained a chromosome arm into one block on the dense
## whole-blood resource (30.8 Mb at a 1 Mb window), and on that partition the
## pre-registered mismatched-list control FAILED under the self-consistent
## attribution convention -- 3.90-fold at P = 1.7e-4. A partition whose blocks
## reach tens of megabases is not a partition into loci.
##
##   Rscript step123_fixed_centre_grid.R
## Output: 123a_fixed_centre_grid.tsv

suppressMessages(library(cqtna))
MR <- "D:/R_ex/MR"
setwd(MR)

fix_gene <- function(g) {
  g <- as.character(g)
  bad <- is.na(g) | !nzchar(trimws(g))
  g[bad] <- paste0("unnamed_", which(bad))
  g
}

mel <- read.csv("landi2020_known_loci_grch38.csv"); mel$source <- "Landi2020"
hcc <- read.csv("84a_hcc_known_loci_grch38.csv");   hcc$source <- "step84"
K_mel <- as_cqtna_known(mel, build = "GRCh38")
K_hcc <- as_cqtna_known(hcc, build = "GRCh38")

cd4 <- read.delim("13_meta_locus_annotation.tsv")
cd4$chr <- sub(":.*", "", cd4$SNP); cd4$pos <- as.numeric(sub(".*:", "", cd4$SNP))
eq  <- read.delim("92c_locus_annotated.tsv")
hh  <- read.delim("85a_HCC_high_annotated.tsv")
hl  <- read.delim("85a_HCC_low_annotated.tsv")

cells <- list(
  list(name = "melanoma x Soskic_CD4",
       d = data.frame(record_id = seq_len(nrow(cd4)), gene = fix_gene(cd4$SYMBOL),
                      chr = cd4$chr, pos = cd4$pos, p = cd4$pval),
       matched = K_mel, mismatched = K_hcc),
  list(name = "melanoma x eQTLGen_blood",
       d = data.frame(record_id = seq_len(nrow(eq)), gene = fix_gene(eq$symbol),
                      chr = eq$chr, pos = eq$pos, p = eq$p_mr),
       matched = K_mel, mismatched = K_hcc),
  list(name = "HCC_high x Soskic_CD4",
       d = data.frame(record_id = seq_len(nrow(hh)), gene = fix_gene(hh$gene_id),
                      chr = hh$chr, pos = hh$pos, p = hh$p_mr),
       matched = K_hcc, mismatched = K_mel),
  list(name = "HCC_low x Soskic_CD4",
       d = data.frame(record_id = seq_len(nrow(hl)), gene = fix_gene(hl$gene_id),
                      chr = hl$chr, pos = hl$pos, p = hl$p_mr),
       matched = K_hcc, mismatched = K_mel))

rows <- list()
cat(sprintf("%-26s %-15s %5s %8s %7s %8s   %17s %18s\n",
            "cell", "partition", "kb", "maxspan", "sig", "bg", "matched", "mismatched"))
cat(strrep("-", 122), "\n")
for (cl in cells) for (meth in c("single_linkage", "fixed_centre"))
  for (kb in if (meth == "single_linkage") 1000 else c(500, 1000)) {
    mr <- as_cqtna_mr(cl$d, locus_kb = kb, build = "GRCh38", locus_method = meth)
    sp <- suppressWarnings(cqtna_locus_spans(mr))
    conv <- if (meth == "single_linkage") "significant_records" else "any_record"
    a <- cqtna_attribution(mr, cl$matched,    known_from = conv)
    x <- cqtna_attribution(mr, cl$mismatched, known_from = conv)
    cat(sprintf("%-26s %-15s %5d %8.0f %7d %8d   %6.2fx P=%-8.3g %6.2fx P=%-8.3g\n",
                cl$name, meth, kb, max(sp$significant_span_kb),
                a$significant_loci, a$background_loci,
                a$fold, a$fisher_p_one_sided, x$fold, x$fisher_p_one_sided))
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cl$name, partition = meth, locus_kb = kb, known_from = conv,
      max_sig_span_kb = round(max(sp$significant_span_kb)),
      n_loci_over_window = sum(sp$significant_span_kb > kb),
      sig_loci = a$significant_loci, sig_known = a$significant_known,
      bg_loci = a$background_loci, bg_known = a$background_known,
      fold = a$fold, fisher_p = a$fisher_p_one_sided,
      mismatch_fold = x$fold, mismatch_p = x$fisher_p_one_sided,
      mismatch_control = if (x$fisher_p_one_sided < 0.05 && x$fold > 1)
        "FAILED" else "clean",
      stringsAsFactors = FALSE)
  }

out <- do.call(rbind, rows)
write.table(out, "123a_fixed_centre_grid.tsv", sep = "\t", row.names = FALSE,
            quote = FALSE)
cat("\nwrote 123a_fixed_centre_grid.tsv\n")
cat("\nmismatched-list control by partition:\n")
print(table(out$partition, out$mismatch_control))
