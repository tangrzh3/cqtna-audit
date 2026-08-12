# Deposit: code and result tables

Companion deposit for *When a pathway can be interrogated genetically, and which
gene gets named: an audit of dynamic-eQTL target nomination in melanoma*.

**⟨repository DOI to be minted at submission⟩** — see "How to mint the DOI" below.

---

## 1. What is here

| | Content |
|---|---|
| `step*.py`, `step*.R` | Every numbered analysis script, in execution order (89 files) |
| `figures/*.py` | Every figure script. Figure titles and in-figure annotations are part of the manuscript's claims and were audited as such (see §4) |
| `*.tsv`, `*.csv` | Numbered result tables. The number prefix matches the step that wrote it |
| `manuscript/` | Manuscript, Methods, all supplementary documents, **all five pre-registrations** |
| `FINDINGS_step5_pigmentation.md` | The day-by-day analysis log, appended only, never retro-edited |
| `HANDOFF_v3.md` | Working state, technical pitfalls, and the methodological discipline rules the project accumulated |
| `assemble.py` | Builds `manuscript/MANUSCRIPT_assembled.md` from the manuscript + Methods sources |

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
Software versions are in Methods §5.19.

## 4. Two audits worth knowing about before you read the code

**Withdrawn claims.** Five claims were withdrawn during revision (listed in Box 1
of the manuscript). A sweep for their phrasing found residues in three places,
**two of which were inside figure scripts** — text printed onto figures rather
than in prose. Figure scripts now carry comments recording what was changed and
why; do not "restore" those strings.

**Numbering.** Supplementary numbering and section cross-references have been
audited by script after each restructure; both are checked in the assemble step.

## 5. How to mint the DOI (the one step left to the authors)

1. Create a **public GitHub repository** and push this directory
   (`git remote add origin …` then `git push -u origin main`).
2. In Zenodo, enable the repository under *GitHub → Repositories*.
3. Cut a release on GitHub (e.g. `v1.0-submission`). Zenodo mints a DOI
   automatically and archives that snapshot.
4. Replace `⟨repository DOI⟩` in Methods §5.20 and at the top of this file with
   the minted DOI, then re-run `python manuscript/assemble.py`.

⚠ Before pushing, confirm no participant-level data is present:

```bash
git ls-files | grep -Ei '\.(rds|h5ad|h5|bam|fastq|vcf)$'
```

That command must return nothing. All analyses here use summary statistics or
publicly deposited de-identified data.
