# Deposit: code and result tables

Companion deposit for *Target nomination from single-variant cis-eQTL Mendelian
randomization re-reads the outcome GWAS: an audit across six diseases and two
exposure resources*.

**⟨repository DOI to be minted at submission⟩** — see "How to mint the DOI" below.

---

## 1. What is here

| | Content |
|---|---|
| `step*.py`, `step*.R`, `step*.sh` | Every numbered analysis script, in execution order (150 files) |
| `figures/*.py` | Every figure script (21 files). Figure titles and in-figure annotations are part of the manuscript's claims and were audited as such (see §4) |
| `*.tsv`, `*.csv` | Numbered result tables. The number prefix matches the step that wrote it |
| `manuscript/` | Manuscript, Methods, every supplementary document (S9-S51) and **every pre-registration and decision record** (`PREREG_*.md`) |
| `FINDINGS_step5_pigmentation.md` | The day-by-day analysis log, appended only, never retro-edited |
| `HANDOFF_v*.md` | Working state, technical pitfalls and the methodological discipline rules the project accumulated, one file per handover; the latest is authoritative |
| `manuscript/assemble.py` | Builds `manuscript/MANUSCRIPT_assembled.md` from the manuscript + Methods sources |
| `container/` | Dockerfile and pinned Python requirements. **Read `CONTAINER.md` first: the container is a runnable environment, not a bitwise reproduction of the machine the reported numbers came from** |

**Deliberately included even though they support no claim** (see Supplementary S12):
the intermediate outputs of analyses we withdrew, and the outputs of tests that
failed their own positive controls. Both are labelled in the log, so the tally of
attempts can be checked against files rather than taken on trust.

## 2. What is NOT here, and where to get it

Third-party data is excluded by `.gitignore` — it must be obtained from its own
source under its own licence. Every input is listed with its accession in
`manuscript/REFERENCES.md`; the large ones are:

| Input | Source |
|---|---|
| Melanoma outcome (FinnGen R12 + Rashkin meta) | FinnGen R12 release; Rashkin et al. 2020 |
| HCC outcomes | GWAS Catalog `GCST90809296`; FinnGen R12 `C3_HEPATOCELLU_CARC_EXALLC` |
| CD4⁺ T dynamic eQTL (exposure) | Soskic et al. |
| Whole-blood cis-eQTL (alternative exposure) | eQTLGen Consortium, 2019-12-11 cis release |
| Single-cell cohorts | GEO `GSE120575`, `GSE115978`, `GSE72056`, `GSE282266`, `GSE199994`, `GSE235863`; KU Leuven RDR `doi:10.48804/GSAXBN` |
| TCGA-SKCM | UCSC Xena |

Large intermediate objects (`*.rds`, `*.h5ad`) are excluded and are regenerable
by the scripts that produced them.

## 3. Reproducing a number

Each numbered result table is written by the step whose number it carries
(`85c_G1_enrichment.tsv` ← `step85_hcc_tests.py`). Scripts read from and write to
the repository root and take no arguments unless documented in their docstring.
Software versions are in the Methods "Software" section; the exact R package
closure is `150a_environment.lock` and the pinned Python set is
`container/requirements.txt`.

## 4. Two audits worth knowing about before you read the code

**Withdrawn claims.** Five claims were withdrawn during revision (listed in the
self-disclosure box of the manuscript). A sweep for their phrasing found residues in three places,
**two of which were inside figure scripts** — text printed onto figures rather
than in prose. Figure scripts now carry comments recording what was changed and
why; do not "restore" those strings.

**Numbering.** Supplementary numbering and section cross-references have been
audited by script after each restructure; both are checked in the assemble step.

## 5. How to mint the DOI (the one step left to the authors)

1. Create a **public GitHub repository** and push this directory
   (`git remote add origin ...` then `git push -u origin main`).
2. In Zenodo, enable the repository under *GitHub -> Repositories*. Zenodo reads
   `cqtna_r/.zenodo.json` for title, licence and creators, so fill the identity
   in **before** cutting the release (step 0 below).
3. Cut a release on GitHub (e.g. `v1.0-submission`). Zenodo mints a DOI
   automatically and archives that snapshot.
4. Write the DOI and the release commit into every file that carries a
   placeholder, in one command rather than by hand:

```bash
python step152_set_identity.py --doi 10.5281/zenodo.XXXXXXX --commit <sha>
```

   Run it with no arguments first; it prints every site and what is still
   unfilled, and exits non-zero while anything remains. Then re-run
   `python manuscript/assemble.py`.

**Step 0, before any of the above** -- the maintainer identity, which the same
script fills:

```bash
python step152_set_identity.py --given <given> --family <family> --email <email> --orcid 0000-0000-0000-0000 --github <owner>/<repo>
```

That touches `cqtna_r/DESCRIPTION`, `LICENSE`, `CITATION.cff`, `.zenodo.json`
and the generated `man/cqtna-package.Rd`. Re-run `R CMD check` afterwards:
`Authors@R` is parsed at build time, so a malformed name fails there and
nowhere else. `cqtna.Rcheck/` is a build product and is deliberately left
alone.

WARNING: Before pushing, confirm no participant-level data is present:

```bash
git ls-files | grep -Ei '\.(rds|h5ad|h5|bam|fastq|vcf)$'
```

That command must return nothing. All analyses here use summary statistics or
publicly deposited de-identified data.

