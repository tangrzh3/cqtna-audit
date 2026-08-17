# CQTNA report (v0.1)

Settings: FDR < 0.05, locus window 1000 kb, known-locus window 1000 kb.

## A. Locus attribution

- background: 58/554 loci (10.47%) carry a known lead SNP for this outcome
- significant: 3/7 (42.9%)
- **enrichment 4.09x, one-sided Fisher P = 0.0281**
- mismatched-list control: 0.0x, P = 1 (clean)

## B. Inference unit

| unit | n_tests | n_significant | n_genes | n_loci |
|---|---|---|---|---|
| record | 3556 | 21 | 10 | 7 |
| variant | 2126 | 15 | 11 | 7 |
| gene | 1255 | 13 | 13 | 9 |
| locus | 554 | 8 | 28 | 8 |

⚠ A shortest list is not evidence that FDR is controlled under dependence. State which unit the conclusions are in.

## C. List stability across outcome GWAS

- genes: 10 vs 6, 5 shared, Jaccard 0.455
- loci: 7 vs 2, Jaccard 0.286
- lost: KANSL1, KIAA0040, MDM4, SMC2, ZFYVE19
- gained: CTU2

## D. Instrument attrition

- of 28 pathway genes: **3 instrumentable, 2 analysable against this outcome, 1 nominally associated**

## E. Compartment attribution

⚠ **configured but empty** — no significant gene appears in the expression table, or no `target_cell_type` was set.

## F. eQTL-to-GWAS peak distance

⚠ **not run** — no input supplied.

## G. Window sensitivity and the distances behind the flag

| known_window_kb | significant_known | significant_loci | background_known | background_loci | fold | fisher_p_one_sided |
|---|---|---|---|---|---|---|
| 100 | 3 | 7 | 13 | 554 | 18.26 | 0.000336 |
| 250 | 3 | 7 | 29 | 554 | 8.19 | 0.00393 |
| 500 | 3 | 7 | 42 | 554 | 5.65 | 0.0115 |
| 1000 | 3 | 7 | 58 | 554 | 4.09 | 0.0281 |

| locus_window_kb | background_loci | significant_loci | fold | fisher_p_one_sided |
|---|---|---|---|---|
| 100 | 1065 | 8 | 5.33 | 0.00382 |
| 250 | 876 | 7 | 4.63 | 0.0203 |
| 500 | 728 | 7 | 4.59 | 0.0208 |
| 1000 | 554 | 7 | 4.09 | 0.0281 |

- distances of significant loci to the nearest known lead SNP (kb): 5, 11, 74, 2,167, 7,275, 20,136, 38,494
- median 2,167 kb, against 10,728 kb for all testable loci
- the 1,000 kb threshold falls between 74 kb and 2,167 kb — **a gap, so the cut is not near your data**

⚠ A window chosen to flatter a result weakens when tightened. If the fold *rises* as the window narrows, the reported value is the conservative one.

## Evidence tiers

| gene | tier | basis | stability across outcome GWAS |
|---|---|---|---|
| CDK10 | **screened** | significant, but on a locus already known for this outcome | retained under the second outcome GWAS |
| CHMP1A | **screened** | significant, but on a locus already known for this outcome | retained under the second outcome GWAS |
| MDM4 | **screened** | significant, but on a locus already known for this outcome | lost under the second outcome GWAS |
| PARP1 | **screened** | significant, but on a locus already known for this outcome | retained under the second outcome GWAS |
| SPATA33 | **screened** | significant, but on a locus already known for this outcome | retained under the second outcome GWAS |
| VPS9D1-AS1 | **screened** | significant, but on a locus already known for this outcome | retained under the second outcome GWAS |
| CTU2 | **unresolved** | significant only under the second outcome GWAS | gained under the second outcome GWAS |
| KANSL1 | **unresolved** | significant on a locus not previously reported for this outcome | lost under the second outcome GWAS |
| KIAA0040 | **unresolved** | significant on a locus not previously reported for this outcome | lost under the second outcome GWAS |
| SMC2 | **unresolved** | significant on a locus not previously reported for this outcome | lost under the second outcome GWAS |
| ZFYVE19 | **unresolved** | significant on a locus not previously reported for this outcome | lost under the second outcome GWAS |

⚠ **No gene can reach `target-supported` from this tool.** The three diagnostics it cannot run are the ones that would license that word.

## Not run here -- these require manual work

- **(ii) colocalisation with explicit multiple-signal modelling and an LD reference matched to the outcome cohort**  
  needs regional summary statistics and an LD panel, not a candidate list
- **(vii) cell-level matching on lineage composition when splitting cells by a score**  
  needs the single-cell data and the split itself; a cluster-level control can pass while the cell-level one fails
- **(viii) code-by-code verification that an endpoint definition is unchanged before comparing candidate lists across releases**  
  needs the phenotype definitions, which release notes do not reliably summarise

A nomination audited on five of eight diagnostics is audited on five of eight.
