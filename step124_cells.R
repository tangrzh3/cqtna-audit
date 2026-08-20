## Shared cell and run definitions for the fixed-anchor grid.
##
## step124 (the grid) and step125 (the per-locus mismatch diagnostics) both
## source this, so there is exactly one definition of what each cell is and
## which region it covers. Two scripts each carrying their own copy is how the
## MHC-excluded input came to be labelled primary in one place and post-hoc in
## another.
##
## Defines: cells, runs, status_of(), and the loaders they need.
## Expects the working directory to be the MR project root.

suppressMessages(library(cqtna))

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
       region = "all_genome", prereg = TRUE,
       d = std(cd4$SYMBOL, cd4$chr, cd4$pos, cd4$pval), m = K_mel, x = K_hcc),
  list(name = "melanoma x eQTLGen_blood", role = "main",
       region = "all_genome", prereg = TRUE,
       d = std(eq$symbol, eq$chr, eq$pos, eq$p_mr), m = K_mel, x = K_hcc),
  list(name = "HCC_high x Soskic_CD4", role = "main",
       region = "all_genome", prereg = TRUE,
       d = std(hh$gene_id, hh$chr, hh$pos, hh$p_mr), m = K_hcc, x = K_mel),
  list(name = "HCC_high x eQTLGen_blood", role = "main",
       region = "all_genome", prereg = TRUE,
       d = with(pick(eh_rec, "HCC_high"), std(gene, chr, pos, p)), m = K_hcc, x = K_mel),
  ## RA is scored over the whole genome as the primary, because that is what S33
  ## registered: its §4 makes N1 = Soskic x RA the primary cell, and its §5
  ## defines the primary test with no MHC exclusion of any kind. The MHC-excluded
  ## rows appear only in S33 §8, the results-registration table, so they are
  ## post-hoc by construction and are carried below with that label.
  ##
  ## An earlier version of this script had the two positions reversed and
  ## annotated the MHC-excluded input as "pre-registered (S33)". That was wrong,
  ## and it was made after seeing that the whole-genome RA x eQTLGen cell fails
  ## its mismatched-list control -- i.e. a convention changed after the result.
  ## The failure stands: RA x eQTLGen is void at the primary setting, and it is
  ## not rescued by the MHC-excluded or the 500 kb version. See CONVENTION_REVIEW.md.
  list(name = "RA x Soskic_CD4", role = "main",
       region = "all_genome", prereg = TRUE,
       d = with(pick(ra_rec, "N1 Soskic x RA"), std(gene, chr, pos, p)),
       m = K_ra, x = K_mel),
  list(name = "RA x eQTLGen_blood", role = "main",
       region = "all_genome", prereg = TRUE,
       d = with(pick(ra_rec, "N2 eQTLGen x RA"), std(gene, chr, pos, p)),
       m = K_ra, x = K_mel),
  list(name = "RA x Soskic_CD4 (MHC excluded)", role = "MHC sensitivity",
       region = "MHC_excluded", prereg = FALSE,
       d = with(pick(ra_rec, "N1 Soskic x RA, MHC excluded"), std(gene, chr, pos, p)),
       m = K_ra, x = K_mel),
  list(name = "RA x eQTLGen_blood (MHC excluded)", role = "MHC sensitivity",
       region = "MHC_excluded", prereg = FALSE,
       d = with(pick(ra_rec, "N2 eQTLGen x RA, MHC excluded"), std(gene, chr, pos, p)),
       m = K_ra, x = K_mel),
  ## HCC-low is a power sensitivity fixed in PREREG_locus_partition.md §3, which
  ## was frozen before these cells were computed; it is pre-registered but is not
  ## one of the six cells.
  list(name = "HCC_low x Soskic_CD4", role = "power sensitivity",
       region = "all_genome", prereg = TRUE,
       d = std(hl$gene_id, hl$chr, hl$pos, hl$p_mr), m = K_hcc, x = K_mel),
  list(name = "HCC_low x eQTLGen_blood", role = "power sensitivity",
       region = "all_genome", prereg = TRUE,
       d = with(pick(eh_rec, "HCC_low"), std(gene, chr, pos, p)), m = K_hcc, x = K_mel))

runs <- list(list(method = "fixed_centre",   kb = 1000, conv = "any_record",
                  label = "main"),
             list(method = "fixed_centre",   kb = 500,  conv = "any_record",
                  label = "sensitivity"),
             list(method = "single_linkage", kb = 1000, conv = "significant_records",
                  label = "legacy"))

## The status of a row is a property of the row, not of the reader's memory.
## A row is prereg_primary only if the cell is one of the six registered cells,
## scored over the region the registration named, at the registered partition.
status_of <- function(cl, rn) {
  ## A post-hoc cell run through the legacy algorithm is both things at once,
  ## and "legacy_reproduction" alone reads as though it reproduces something
  ## published. The MHC-excluded rows do appear in S33's results table, so they
  ## are reproducible, but they were never a registered analysis.
  if (rn$label == "legacy")
    return(if (isTRUE(cl$prereg)) "legacy_reproduction"
           else "legacy_posthoc_sensitivity")
  if (!isTRUE(cl$prereg)) return("posthoc_sensitivity")
  if (cl$role == "main" && rn$label == "main") return("prereg_primary")
  "prereg_sensitivity"
}
