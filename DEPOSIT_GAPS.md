# Deposit gaps: known, recorded, and not closable today

Things the deposit does not contain that it arguably should. Each is reported
by an audit on every run rather than fixed silently, and each says why it
cannot simply be repaired.

---

## 1. `150a_environment.lock` omits two packages the Methods names

`chromVAR 1.26.0` and `org.Hs.eg.db 3.19.1` appear in the Software section and
in `container/install_r_packages.R`'s `need` list, but **not** among the lock's
226 entries. Confirmed by key against every entry, not by a failed lookup.

In `cqtna-audit:0.3.0` both are present **at exactly the stated versions** —
installed by the fallback that handles anything `renv::restore()` did not
restore. ⚠ **That path is not version-pinned**, so they are correct by luck
rather than by record; a rebuild on another day could get different versions
and nothing would say so.

**Why it is not simply fixed**: regenerating the lock means running `step150`,
which rebuilds it from `installed.packages()` on the **current** machine. That
machine has drifted — numpy 2.5.2 against the recorded 2.0.0 (S54 §9.8) — so
regenerating would overwrite the authoring-machine record with a later one and
destroy the baseline every S54 comparison is made against. **The record being
incomplete is better than the record being wrong.**

**Reported by**: `step127` §1n, as `GAP`, every run.

⚠ **Consequence for the closure check**: `step159`'s gate compares installed
packages against the lock, so a package the lock omits is outside its frame.
226/226 was green while these two went unexamined. The closure figure means
"everything the lock records is present at the recorded version", not
"everything the paper uses is pinned".

## 2. PLINK 2.0.0-a.7.2 is no longer obtainable

`s3://plink2-assets/alpha7/` now holds only 2026 builds;
`plink2_linux_x86_64_20241206.zip`, the version the Methods names, has been
deleted upstream (established by listing the bucket, not by reading a 404).

**Not substituted.** Putting a current build behind a version number the paper
states would make the image quietly wrong rather than openly incomplete.
`/opt/binaries_provenance.txt` in the image records SMR's checksum and PLINK's
absence.

**No analysis on S54's must-run list invokes either binary**, so the acceptance
is unaffected; a reader trying to rerun everything is.

## 3. Four analyses cannot run inside the container

`step158` (S53 visibility overlap) needs externally downloaded bulk RNA-seq;
`step120`, `step121`, `step122` point at authoring-machine paths and at
`landi2020_known_loci_grch38.csv`, which is not distributed. These are data
distribution decisions taken before the container existed, not container
defects, and `step159` records them as `not rerunnable` rather than as
failures.

## 4. Audit coverage of the manuscript's arithmetic is partial

`step160` measures it directly. It is not complete and the uncovered values are
listed in `160a_audit_coverage.tsv`, so a reviewer can see exactly which
numbers no audit checks rather than inferring that a green run means every
number was verified. Two errors have already been found inside that uncovered
set (S54 §9.8's wrong-table 13.53, and the Steiger lower bound), which is the
argument for keeping the list visible rather than quoting a coverage headline.

## 5-6. RESOLVED (2026-09-10): both were my errors, not the paper's

Two entries here previously read "I could not reconcile this" -- the HEIDI
denominator and the winner's-curse medians. **Both now reconcile exactly, and
both failures were mine.**

**HEIDI 253 / 291 = 86.9%.** I had used the wrong tables (`14` and `08` rather
than `16` and `15`) and, more importantly, applied no `PP.H4 < 0.2` filter. The
denominator is not every HEIDI test: it is the records **coloc calls distinct
causal variants**, which is the entire point of the sentence. Encoded in
`step127` section 2h.

**Winner's curse, 0.125 versus 0.202.** Again `14` where it should have been
`16`. The groups are the 244 exposures new to the meta round against the 127
carried over. Encoded in `step127` section 2i.

> ### Worth keeping as a record of how this went wrong
>
> In both cases I reported "the paper is probably right and I cannot confirm
> it" and filed the question for the author -- when the definitions were
> written down in this repository the whole time: the HEIDI one in
> `figures/make_gb_fig4_coloc_heidi.py`, which draws the very same panel, and
> the winner's-curse one in `FINDINGS_step5_pigmentation.md`.
>
> **I did not look.** Deferring to the author was the cautious-LOOKING move and
> it was the wrong one. Asking a person to resolve something the repository
> already answers has a real cost, and here it was entirely avoidable by
> reading the code that had computed the number in the first place.

Both definitions now live in `step127` rather than only in a figure script and
a findings note, so neither number depends on anyone re-reading a file that
nobody opens.

## 7. Two Results numbers with no producing code in the deposit

**`r = 0.937`** — "the two releases' statistics correlate at |z| r = 0.937".
No script in the repository computes a correlation between the R12 and R13
statistics: searched every `step*.py` that mentions R13 for `corrcoef`,
`spearmanr`, `pearsonr` or `.corr(`, and none has one. The R13 outputs that do
exist (`99a`–`99c`) carry counts and attribution, not per-record statistics, so
the number cannot be recomputed from what is deposited either.

⚠ This is the same class of defect `step143` was written to address for the
lung and colorectal locus lists: **a number in the Results whose producing code
is not in the repository.** It is recorded rather than quietly recomputed by a
guessed method, because a reader cannot check the value against anything.

**`0.45%`** — the `--diff-freq` QC removal rate. This comes from SMR's own log
output, documented in `FINDINGS_step5_pigmentation.md`, not from any table
here. Not a defect so much as a category: a figure quoted from an external
tool's console, which no audit over deposited tables can reach.

Both are marked `unauditable` in `step160` with these reasons attached, so
neither drags the coverage figure down as though someone had simply not got to
them.

## 8. Four values whose producing table I could not locate

⚠ **These are NOT classified un-auditable.** They stay in the "auditable, not
done" column, because "I did not find it" and "it cannot be found" are
different claims and only the second belongs in that category. Recording what
was searched so the next attempt does not start over:

| value | claim | searched |
|---|---|---|
| `5.33` / `1.59` | significant-Q rate, correct vs `exp(-Q/2)` form | `11b_heterogeneity` is the only table carrying Cochran's Q and gives 18.25% / 33.46% over 263 records; `12`/`13` (the two-study meta sets) carry no Q column |
| `1.82` | melanoma cell's margin over its best rival BEFORE the prostate comparator was repaired | `130c` holds the post-repair 2.31; the pre-repair value is historical and no table found retains it |
| `2.75` | TPI1 log2 units higher in malignant cells than CD4⁺, 16 of 16 patients | `77b` gives 0.83 for a different dataset (Pozniak); no 16-row per-patient compartment table located |
| `88.86` / `89.73` / `89.87` | the two chr16 bounded-locus windows | `125a` carries locus coordinates but only for the mismatch cells, not melanoma × CD4; derivable by re-running the fixed-centre partition, which is not the same as reading it off a table |

The last one is worth separating: it **is** recomputable, by applying the
partition to melanoma's significant records. I did not, because three of my
"discrepancies" in this exercise turned out to be my own reconstruction being
wrong, and a fourth guess at a definition is a poor trade against three
coordinate digits.
