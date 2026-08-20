# cqtna — Context-specific QTL Target Nomination Audit

**cqtna is a research companion implementing the automatable diagnostics used in
one study. It is not a validated general-purpose target-nomination platform.**

An audit for candidate lists produced by cis-eQTL-instrumented Mendelian
randomization. You hand it a list you already have; it reports what kind of
evidence each gene actually carries.

It does **not**: run the MR; perform colocalisation, multiple-signal modelling or
LD matching; do cell-level lineage matching; verify that an endpoint definition is
unchanged across GWAS releases; or produce a `target-supported` conclusion. Its
outputs are for auditing and sensitivity analysis, not for deciding that a gene is
a target.

Base R only — no Bioconductor, no tidyverse, no Seurat.

```r
# install.packages("remotes")
# from the built tarball that ships with the review packet
install.packages("cqtna_0.3.0.tar.gz", repos = NULL, type = "source")
# (no public repository yet; see RELEASE.md before publishing one)
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
| A | locus attribution against the outcome's own known loci, by bounded locus, with a mismatched-list negative control | MR results + known loci |
| B | the significant list recomputed over record, variant, gene and bounded locus | MR results |
| C | list stability against a second outcome GWAS | a second MR table |
| D | instrument attrition: instrumentable → analysable → associated | pathway gene table |
| E | compartment attribution: per-cell-type expression ratio | expression table |
| F | distance between eQTL and GWAS peaks | peak table |
| G | window sensitivity: both window conventions swept at 100/250/500/1000 kb, plus the continuous distances behind the binary flag | nothing extra |

Alongside these: `cqtna_locus_spans()` reports whether single-linkage clustering
has chained distinct regions into one block, and `cqtna_permutation_control()`
draws a null matched on locus size and instrument count — the control the
mismatched list cannot provide.

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
disease-specific; for that, run `cqtna_permutation_control()`, which matches null
loci on size and instrument count.

## What the numbers do and do not mean

Read these before quoting anything from a report.

**The Fisher p-value is descriptive.** It treats loci as exchangeable units, and
they are not independent tests in that sense whichever partition produces them.
The package calls them *bounded* loci rather than *independent* loci for exactly
this reason: `"fixed_centre"` guarantees that two loci share no variant, which is
not the same as statistical independence, and nothing here establishes the
latter. Report the result as an enrichment with its sensitivity analysis, not as
a calibrated p-value.

**Simes-then-BH across inference units is a sensitivity analysis.** It shows how
much the list depends on the unit you chose. It is not a proof that FDR is
controlled under this correlation structure, and module B says so on its face.

**The mismatched-disease list is a weak negative control.** It rules out
enrichment on loci that are dense in every disease. It does not rule out a density
that is itself disease-specific. `cqtna_permutation_control()` addresses that by
matching null loci on size, instrument count and gene count.

⚠ **Do not quote a single p-value from it.** The verdict moves with the matching
tolerance, so the tolerance has to be argued for before you look. On the shipped
demo, same seed and 10,000 permutations, the four tolerances that achieve complete
matching give empirical P between 0.0019 and 0.0096, while the tightest matches
only 75% of loci and correctly returns `NA`. Run
`cqtna_permutation_sensitivity()` and report the whole sweep.

⚠ Two withdrawn figures, kept here so neither reappears. An earlier README quoted
"2.32-fold, P = 0.111", computed while one of seven loci had no matched pool.
An earlier one quoted a 0.022-0.042 range, computed on the single-linkage
partition before the default changed. **Both are withdrawn.** A permutation
p-value is only interpretable together with the partition and the convention it
was computed under, so quote all three or none.

**"The fold rises as the window narrows" means the reported fold is the
conservative one for that fold.** It does not make the wider inference
conservative, and it is not evidence that the attribution is real.

**Locus-level counting assumes loci are loci.** In a dense exposure resource,
1 Mb single-linkage chains a chromosome arm into one block: in the source study
the whole-blood resource produced thirty loci wider than 10 Mb and one significant
"locus" spanning 30.8 Mb across 588 records. Locus-level counts there are counting
blocks. This is why `locus_method` defaults to `"fixed_centre"`, which bounds a
locus at `locus_kb` by construction; `"single_linkage"` remains available and is
what you want only when reproducing a published number computed that way.
`cqtna_locus_spans()` warns when chaining is happening.

A useful diagnostic falls out of this. The three `known_from` conventions
disagree only when the partition chains: on the demo's bounded partition
`"any_record"` and `"significant_records"` return the same 4 of 8 and the same
4.96-fold, whereas on the same data under single linkage they differ. **If your
choice of convention changes the answer, suspect the partition before you argue
about the convention.**

There is no window at which the problem disappears. Tightening to 100 kb still
leaves five significant loci wider than five times the window; what changes is
whether the residual chaining crosses a classification boundary. Where the three
`known_from` conventions agree, that is evidence the chaining did not bite on
those data -- not evidence that the partition is clean.

## Reading module G

Two different conventions get called "the 1 Mb window", and they carry different
weight, so they are swept separately: `known_kb` decides what counts as landing
on a known locus and moves the fold directly; `locus_kb` decides how independent
loci are defined and moves the denominator.

A threshold chosen to flatter a result weakens when tightened. If your fold
*rises* as the window narrows, the value you reported is the conservative one for
that fold. Module G reports the nearest significant-locus distance below and above
your threshold and attaches no verdict: an earlier version printed "falls in a
gap" when their ratio exceeded 4, a cutoff with no basis beyond looking reasonable
on one dataset.

## Demo

`cqtna_demo()` ships the audit's own published melanoma data.

Under the package defaults (`fixed_centre`, 1 Mb) it gives **4 of 8 significant
loci on known melanoma loci against a 10.1% background, 4.96-fold, one-sided
P = 0.0048**, with the mismatched HCC list at 0.00-fold.

Setting `locus_method = "single_linkage"` reproduces the number the source study
published, **3 of 7 against 10.5%, 4.09-fold, P = 0.028** — kept reproducible on
purpose, and carrying no inference. The two differ because the partition differs,
not because the data does; quoting either without saying which partition produced
it is the mistake this package exists to make visible.

Module B shows the same list spanning 10 genes at record level and 28 at locus
level, because a significant locus does not name a gene.

## Testing

Four suites. The first asserts a published paper's headline numbers directly. The
second compares this implementation against the Python reference (`cqtna.py`) on
the same data, because two implementations of one audit are only safe while
something compares them. The third builds synthetic tables that separate failure
mechanisms the real data cannot — several of the bugs this package has fixed were
invisible on the demo fixtures. The fourth checks the input contract.

`inst/validate_synthetic.R` goes further: it plants a known amount of enrichment
in simulated data and checks the tool recovers it, and plants none and checks the
tool does not manufacture it.

```r
testthat::test_local()
```

---

A tool for detecting selective emphasis should be run first against the study
proposing it.
