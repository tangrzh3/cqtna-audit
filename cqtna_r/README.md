# cqtna — Context-specific QTL Target Nomination Audit

An audit for candidate lists produced by cis-eQTL-instrumented Mendelian
randomization. You hand it a list you already have; it reports what kind of
evidence each gene actually carries.

Base R only — no Bioconductor, no tidyverse, no Seurat.

```r
# install.packages("remotes")
remotes::install_github("OWNER/cqtna", subdir = "cqtna_r")
```

```r
library(cqtna)

au <- cqtna_audit(
  mr       = my_mr_results,       # one row per test: record_id, gene, chr, pos, p
  known    = known_loci,          # lead SNPs previously reported for YOUR outcome
  mismatch = another_disease,     # negative control
  build    = "GRCh38"
)

au                          # the report
summary(au)                 # the counts
plot(au, "distance")        # the distances behind the binary flag
cqtna_report(au, "audit.md")
```

## What it runs

| | module | needs |
|---|---|---|
| A | locus attribution against the outcome's own known loci, by independent locus, with a mismatched-list negative control | MR results + known loci |
| B | the significant list recomputed over record, variant, gene and independent locus | MR results |
| C | list stability against a second outcome GWAS | a second MR table |
| D | instrument attrition: instrumentable → analysable → associated | pathway gene table |
| E | compartment attribution: per-cell-type expression ratio | expression table |
| F | distance between eQTL and GWAS peaks | peak table |
| G | window sensitivity: both window conventions swept at 100/250/500/1000 kb, plus the continuous distances behind the binary flag | nothing extra |

**Three of the eight diagnostics cannot be automated** from a candidate list, and
the report names them in every run rather than omitting them:
colocalisation with multiple-signal modelling and matched LD; cell-level lineage
matching; and code-by-code endpoint verification before comparing releases.
`cqtna_not_automated()` returns them, with the failure mode for each.

**Nothing reaches `target-supported` from this tool**, by design. The checks that
would license the word are the three it cannot run.

## Three things to get right before you run it

**Declare your genome build.** `build` is required on both tables and a mismatch
is an error, not a warning. A known-locus list on the wrong build does not fail
loudly; it just moves every locus and returns a plausible number.

**Keep the provenance of your known-locus list.** If any of it comes from the
outcome GWAS's own publication, asking whether your hits fall on "known" loci is
circular to that extent. Keep a `source` column. Note that dropping
outcome-derived loci can only shrink the background, so it can only *raise* the
enrichment — an over-corrected list is not a stricter test and should not be
cited as one.

**Supply the mismatched list.** Without module A's negative control you cannot
tell outcome-specific attribution from loci that are dense in every disease. Even
with it, a mismatched list does not exclude a density that is itself
disease-specific; that needs a permutation matched on instrument and gene
density, which this tool does not do.

## Reading module G

Two different conventions get called "the 1 Mb window", and they carry different
weight, so they are swept separately: `known_kb` decides what counts as landing
on a known locus and moves the fold directly; `locus_kb` decides how independent
loci are defined and moves the denominator.

A threshold chosen to flatter a result weakens when tightened. If your fold
*rises* as the window narrows, the value you reported is the conservative one.
Module G also tells you whether your threshold falls in a gap in the distances or
sits on top of them — the second case means your result depends on that choice.

## Demo

`cqtna_demo()` ships the audit's own published melanoma data. It reproduces the
paper: 3 of 7 significant loci on known melanoma loci against a 10.5% background,
4.09-fold, one-sided P = 0.028, with the mismatched HCC list at 0.00-fold. Module
B shows the same list spanning 10 genes at record level and 28 at locus level,
because a significant locus does not name a gene.

## Testing

The regression suite is a published paper's results — every headline number is
asserted directly. A second suite compares this implementation against the
Python reference (`cqtna.py`) on the same data, because two implementations of
one audit are only safe while something compares them.

```r
testthat::test_local()
```

---

A tool for detecting selective emphasis should be run first against the study
proposing it.
