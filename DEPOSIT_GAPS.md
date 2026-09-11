# Deposit gaps: known, recorded, and not closable today

Things the deposit does not contain that it arguably should. Each is reported
by an audit on every run rather than fixed silently, and each says why it
cannot simply be repaired.

---

## 1. `150a_environment.lock` omits three packages the image needs

`chromVAR 1.26.0` and `org.Hs.eg.db 3.19.1` appear in the Methods Software
paragraph and in `container/install_r_packages.R`'s `need` list, but **not**
among the lock's 226 entries. **`testthat 3.2.1.1` is absent too** — not named
in the Methods, but it is what `step159` phase 1 runs the 271 `cqtna`
assertions with, so the gate's own instrument was the unpinned one. Confirmed
by key against every entry, not by a failed lookup. (S54 §9.1 already noted
testthat as a test-time dependency; what was wrong here was the count.)

### ⚠ What this section said until 2026-09-11, and why it was wrong

It said the gap could not be fixed, and gave as the reason that regenerating
the lock would overwrite the authoring machine's record with a drifted one.

**That reason is correct and still stands** — `step150` rebuilds from
`installed.packages()` on the current machine, which has moved (numpy 2.5.2
against the recorded 2.0.0, S54 §9.8) — **but it establishes only that `150a`
must not be REGENERATED.** It says nothing about whether the gap can be closed
some other way, and no other way had been considered. Third-party review made
the point. The section had substituted one claim for another in the project's
own favour, which is the failure mode this file exists to catch.

### The repair, which does not touch `150a`

`150c_unlocked_packages.tsv` — a supplementary record deposited beside the
lock, written by hand and never by a script, naming each package the lock omits,
the version this image must produce, and where that version comes from.

`container/install_r_packages.R` now reads it and pins. Bioconductor packages
cannot be pinned by package version, so those are pinned by **release** and
then asserted; CRAN packages go through `remotes::install_version`. Every entry
is checked whether or not it was missing, because a base image supplying one at
some other version is exactly the drift this is for. A mismatch stops the build
with both versions named.

The closure report is now two files, deliberately not merged:
`/opt/150a_closure_report.tsv` keeps its 226 rows, so **226/226 goes on meaning
"everything `150a` records is present at the recorded version"**, and
`/opt/150c_closure_report.tsv` carries the three the lock never framed.
`step159` phase 1 reads both.

### ⚠ What is verified, and what is not

**The pin is in the build recipe; the image in hand predates it.**
`cqtna-audit:0.3.0` was built on 2026-09-08. Running the 150c logic against
that image on 2026-09-11 gives `chromVAR 1.26.0`, `org.Hs.eg.db 3.19.1`,
`testthat 3.2.1.1` — all three matching, with the Bioconductor release pin
(3.19) matching the image's own — and nothing absent from both records. So the
record is **corroborated for this image**. What a rebuild adds is that the pin
is **enforced at install time** rather than confirmed after the fact, and
`step159` says exactly which of the two it is looking at rather than letting an
absent report read like a passing one.

⚠ The record remains made **after** the fact, from what the image resolved,
not from the authoring machine. For `chromVAR` and `org.Hs.eg.db` the Methods
states the same versions independently, so those two are corroborated twice
over. For `testthat` **nothing states a version** — 3.2.1.1 is simply what the
image had, now written down so it stops moving. That is pinning, not
provenance, and the file's own `how_the_version_was_established` column says so.

⚠ Also unchanged: S54 §9.1 records the image as carrying `testthat 3.3.2`.
That line describes the **slim** image, whose runs are void under §9.9. The
full image has 3.2.1.1. The line is stale, not wrong about its own subject.

**Reported by**: `step127` §1n, as `GAP`, every run.

⚠ **Consequence for the closure check, unchanged**: a package the lock omits is
outside the 226-row frame. 226/226 was green while these three went unexamined.
The closure figure means "everything the lock records is present at the recorded
version", not "everything the paper uses is pinned" — which is why the second
report exists.

## 1b. ⚠ 102 of 158 analysis scripts named one machine — fixed 2026-09-11

Until 2026-09-11, `MR = r"D:/R_ex/MR"` with no override appeared in **102 of the
158 `step*` scripts**. A reader holding the container and the repository could
not run them: they failed on a path, not on absent data. By number range, 79 of
95 below `step100`, 20 of 21 in 100–119, 3 of 42 from 120 up — the later the
script, the likelier it took an argument.

⚠ **This distorted the coverage figure S54 §8 turns on, and did so invisibly.**
`step159`'s group 7 globbed `step1[2-5][0-9]`, which almost exactly traced the
portability boundary, so widening it to `step*` — the whole point of the
first 2026-09-11 amendment — moved measured coverage by **zero**, 82/191 both
times. The number had been reporting which scripts happen to be portable, not
what the container can reproduce. With the paths fixed the same run reaches
131/191, and group 7 goes from 50 scripts running to 102.

Fixed by giving every one the `argv` / `CQTNA_DIR` / same-default form that
`step101` already used; the default is unchanged in every file, so behaviour on
the authoring machine with no argument is what it always was.

⚠ The fix was made **after** the shortfall was observed and required the
author's explicit authorisation; registered in S54 §10.3 with the ordering, on
the precedent of §9.8. §10.3 also records that the mechanical check first
offered as verification checked the shape of the changed lines and not what
`argv[1]` already meant, and that two scripts had meant something else by it.

## 1c. A permutation that looked seeded and was not (`step29`) — fixed 2026-09-11

`step29_glyco_pathway_genetics.py` line 38 creates `rng =
np.random.default_rng(1)`. The 5000 matched draws at line 176 do not use it:

```python
pool.absz.sample(len(sub), replace=True, random_state=None)
```

So three columns of `35e_pathway_enrichment.tsv` differ on every run. Measured
across two runs of the **same image on the same machine**, which removes the
environment as an explanation:

| column | run 2 | run 3 |
|---|---|---|
| `null_mean` | 0.9384705771575401 | 0.9376783838308392 |
| `null_sd` | 0.26708658312065875 | 0.27100291396109 |
| `p_matched_permutation` | 0.5362927414517097 | 0.5302939412117577 |

`obs_mean_absz` and `p_mannwhitney_unmatched` are identical across both runs,
so the deterministic parts are deterministic.

⚠ **The unused seeded generator on line 38 is worse than no seed at all**: it
tells anyone reading the file that the permutation is reproducible. S54 §1
rules randomness out as an explanation precisely on the grounds that seeded RNGs
are deterministic, and directs that any variation be traced to one of five
listed causes rather than called random fluctuation — here it really is random
fluctuation, because nothing seeded it.

**Scope, checked rather than assumed**: no value in `35e` is among the numbers
the text cites. `step127` does not reconcile it and none of its values appear in
`160a_audit_coverage.tsv`'s 191 rows. No acceptance criterion and no migration
under S54 §5 is affected.

**Fixed 2026-09-11, on the author's authorisation** (S54 §9.8 precedent, since
this is analysis logic changed after results were observed): line 176 now reads
`random_state=rng`, passing the generator that line 38 had been creating and
never using. Verified by running `step29` twice in the container — the two
`35e_pathway_enrichment.tsv` are byte-identical. The values it now fixes on are

| column | value |
|---|---|
| `null_mean` | 0.9411048379462176 |
| `null_sd` | 0.26801494613776417 |
| `p_matched_permutation` | 0.5366926614677064 |

⚠ **This is one arbitrary draw replacing another, not a correction of a wrong
number.** The two pre-fix P values were 0.5363 and 0.5303 against a Monte Carlo
standard error of 0.0071 for 5000 draws — 0.85 SE apart, neither near
significance. Nothing about the finding changed; what changed is that it can now
be reproduced. No §5 migration is registered because no value here reaches the
text.

⚠ **Only running the same container twice finds this.** A cross-environment
comparison cannot: both sides vary, so the difference is charged to the
environment. A single run cannot see it at all.

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

## 3. Analyses that cannot run inside the container

`step158` (S53 visibility overlap) needs externally downloaded bulk RNA-seq;
`step120`, `step121`, `step122` need `landi2020_known_loci_grch38.csv`, which
is not distributed. These are data distribution decisions taken before the
container existed, not container defects, and `step159` records them as
`not rerunnable` rather than as failures.

⚠ **This section used to say those three "point at authoring-machine paths",
folding two different defects into one sentence.** They are different: an
undistributed input is a deposit decision, a hardcoded path is a portability
bug, and the second turned out to affect 102 scripts rather than three — see
§1b. The paths are fixed; the absent input remains absent, and that is the only
reason these three still cannot run.

After the 2026-09-11 fix, group 7 attempts 141 scripts and 102 run. The 43 that
do not divide into: undistributed GSE raw data (single-cell, spatial,
chromatin, `GSE199994`); **no producing script in the deposit at all** (see
§7b); scripts needing a caller-supplied input file (`step91` a PMID list,
`step21` a ligand-receptor table); and one live-service timeout (`step137`,
which queries the GWAS Catalog and exceeded the 3600 s cap after completing in
888 s a run earlier — the duration is not under this project's control).

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

## 7b. Six tables `step127` audits have no producing script in the deposit

Established 2026-09-11 while tracing why widening group 7 changed nothing:

| table | producer |
|---|---|
| `08_SMR_HEIDI_results.tsv` | **none in the repository** |
| `09_steiger_filtering.tsv` | **none** |
| `13_meta_locus_annotation.tsv` | **none** |
| `14_coloc_meta_results.tsv` | **none** |
| `15_SMR_meta_results.tsv` | **none** |
| `16_coloc_meta_results.tsv` | **none** |
| `12_MR_meta_strict.tsv` | `rerun_mr_meta.py` — exists, but is **not** named `step*`, so S54 §3 row 7 never reached it |
| `07_coloc_results.tsv` | `step6_coloc.R` — the table number and the step number differ by one |

Only readers of these files are in the repository. `git log --all
--diff-filter=A` confirms **the only file ever committed named `step0`–`step18`
is `step6_coloc.R`**: steps 7 to 18 were never scripts. They are the early
pigmentation-stage work recorded in `FINDINGS_step5_pigmentation.md`, which
predates the `stepNN_*.py` convention.

⚠ **This is a ceiling on rerunnable coverage that no amount of widening can
lift**, and it is the same class as §7 — a number whose producing code is not
deposited — but at the level of whole tables rather than single values. The
`step127` sections affected are 1k, 1t, 1x, 2h, 2i and 2t.

It is recorded here rather than argued around in the coverage figure: S54
§9.11.1 counts these as not regenerated, which is what they are.

## 8. RESOLVED (2026-09-10): all four were traceable, and I had not looked

This section previously listed four values "whose producing table I could not
locate". **All four are now reconciled**, and every definition was written down
in this repository the whole time.

| value | where it actually was |
|---|---|
| `2.75` | `68c_tests.tsv`, the H1 TPI1 row: diff 2.75, 16 of 16, P = 3.1e-5. Produced by `step68_celltype_attribution.R`, described in `FINDINGS_step5_pigmentation.md`. I had checked `77b`, a different dataset. |
| `5.33` / `1.59` | `meta_finngen_rashkin.py` computes it over the 8.79 M variants present in both studies. I had checked `11b_heterogeneity`, which is a different analysis at 263 records. |
| `1.82` | A **superseded** value: 4.962/2.73 with the pre-repair prostate comparator. Derivation in `findings.md`; the repaired 2.31 is checked in `step127` 1u. |
| `88.86` / `89.73` / `89.87` | Recomputable with `step85e`'s `assign_loci`, whose docstring states it is identical to `cqtna:::cq_assign_loci`. Reusing the deposited implementation removed the risk that had stopped me. |

### ⚠ A rounding trap worth keeping

`5.33%` cannot be reproduced from the deposited meta file, and the reason is
not an error in the paper. `meta_finngen_rashkin.py` counts `qp < 0.05` on the
**unrounded** value, then writes `f"{qp:.4g}"`. Ninety-two records round to
exactly `"0.05"` on the way out, so recomputing from the stored column gives
5.32% while the script's own count gave 5.33%. `step127` 2v reproduces the
number by counting those ninety-two back in, and says why in the code.

**The stored column cannot answer the question that produced it.** Anywhere a
count is made before a rounding write, the file is not a substitute for the
computation.

### What this cost

Three times in this audit I reported "the paper is probably right and I cannot
confirm it" and handed the question back, when the answer was in a script or a
findings note I had not opened. Deferring looks careful and is not: it moves
work to someone who has to re-derive what the repository already records.
**Search the producing code before declaring anything unreconcilable.**
