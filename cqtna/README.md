# CQTNA — Context-specific QTL Target Nomination Audit

A runnable audit for candidate lists produced by cis-eQTL-instrumented Mendelian
randomization. It takes a list you already have and reports what kind of evidence
each gene actually carries.

Dependencies: `numpy`, `pandas`. Nothing else.

```bash
python cqtna.py --demo                    # runs on the audit's own melanoma data
python cqtna.py --config my_audit.json    # runs on yours
```

## What it runs

| | module | needs |
|---|---|---|
| A | locus attribution against the outcome's own known loci, by independent locus, with a mismatched-list negative control | MR results + known-locus list |
| B | the significant list recomputed over record, variant, gene and independent locus | MR results |
| C | list stability against a second outcome GWAS | a second MR results table |
| D | instrument attrition: instrumentable → analysable → associated | pathway gene table |
| E | compartment attribution: per-cell-type expression ratio | expression table |
| F | distance between eQTL and GWAS peaks | peak table |

**Three of the eight diagnostics cannot be automated** and are listed by name in
every report: colocalisation with multiple-signal modelling and matched LD;
cell-level lineage matching; and code-by-code endpoint verification before
comparing releases. See [`cqtna_manual.md`](cqtna_manual.md), which gives the
failure mode for each.

## Input format

Only two files are required. Tab- or comma-separated; the delimiter is sniffed.

**`mr_results.tsv`** — one row per test

| column | | |
|---|---|---|
| `record_id` | required | unique |
| `gene` | required | |
| `chr`, `pos` | required | instrument position, any consistent build |
| `p` | required | MR p-value; FDR is recomputed, not read |
| `exposure_profile` | optional | cell state, timepoint, tissue |

**`known_loci.tsv`** — previously reported lead SNPs **for your outcome trait**

| column | | |
|---|---|---|
| `chr`, `pos` | required | same build as the MR table |
| `rsid`, `source` | optional | `source` is worth keeping — see below |

Optional: `mr_results_alt.tsv` (same schema, second outcome GWAS),
`known_loci_mismatched.tsv` (a different disease's list, as negative control),
`instruments.tsv` (`gene`, `instrumentable`, `analysable`, `associated`),
`expression.tsv` (`gene`, `cell_type`, `mean_expression`),
`peaks.tsv` (`gene`, `eqtl_pos`, `gwas_pos`).

Config is JSON: paths plus `fdr_threshold`, `locus_window_kb`,
`known_window_kb`, `target_cell_type`. See `demo/demo_config.json`.

## Two things to get right before you run it

**Keep the provenance of your known-locus list.** If any of it comes from the
outcome GWAS's own publication, asking whether your hits fall on "known" loci is
circular to that extent. Keep a `source` column, and check each locus for
independent support. Dropping outcome-derived loci can only shrink the
background, so it can only *raise* the enrichment — an over-corrected list is not
a stricter test and should not be cited as one.

**Supply the mismatched list.** Without module A's negative control you cannot
tell outcome-specific attribution from loci that are dense in every disease. Even
with it, a mismatched list does not exclude a density that is itself
disease-specific; that needs a permutation matched on instrument and gene
density, which this tool does not do.

## Reading the output

Tiers are `screened` / `unresolved` / `state-informative`, with stability across
outcome GWAS reported alongside rather than folded in — a gene on a known locus
that also fails to replicate is still a gene on a known locus.

**Nothing reaches `target-supported` from this tool**, by design. The checks that
would license the word are the three it cannot run.

## Demo

`--demo` runs the audit on its own published melanoma data and reproduces the
figures in the paper: 3 of 7 significant loci on known melanoma loci against a
10.5% background, 4.09-fold, one-sided P = 0.028, with the mismatched HCC list at
0.00-fold. Module B shows the same list spanning 10 genes at record level and 28
at locus level, because a significant locus does not name a gene. Module C shows
the two outcome rounds sharing 5 of 11 genes.

A tool for detecting selective emphasis should be run first against the study
proposing it.
