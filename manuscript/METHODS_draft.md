# Methods — draft v1 (2026-08-11)

> 结构说明（中文，不进正文）
> 顺序按分析链条，不按数据来源。凡是"跑之前写死"的判据（判废线、窗口、方向）一律写进对应小节，
> 并在 §18 单列纪律小节 —— 那一节是本文方法学定位的一部分，不要删。
> 已知的坑（GRCh38 不 liftover、`nzchar(NA)`、peak 必须按区间对齐、R² 无上界等）写在会影响
> 结果可复现性的位置。⟨⟩ 为待补。

---

## 1. Exposure data: dynamic CD4⁺ T cell cis-eQTLs

Cis-eQTL summary statistics were taken from a published activation time-course
study of primary human CD4⁺ T cells (Soskic et al.), comprising eight profiles:
naive and memory cells at 0 h (resting), 16 h (pre-division), 40 h (after the
first division) and 5 d (effector), from 85–100 donors per profile
(naive 99/99/89/85; memory 100/95/89/90). Per-profile sample sizes were
recovered from `ma_count/(2·MAF)` and were internally consistent for 100% of
variants in all eight profiles. Expression had been inverse-normal transformed,
so exposure units are standard deviations.

Coordinates are GRCh38 and were used as distributed; no liftover was performed
(an early liftover step, based on an incorrect assumption of hg19, was reverted).
All cis-window summary statistics were read from the cleaned per-profile parquet
files, which retain the full cis window required for colocalisation.

## 2. Outcome GWAS and meta-analysis

The primary outcome was a fixed-effect inverse-variance meta-analysis of two
melanoma GWAS with full summary statistics:

- FinnGen R12 `C3_MELANOMA_SKIN_EXALLC`, 5,753 cases / 378,749 controls
  (local release; the OpenGWAS `finn-b-` snapshot was **not** used, as it is an
  R5 extract with 98 cases);
- Rashkin et al. (GCST90011809), 6,777 cases.

The meta-analysis comprises 12,530 cases and 789,099 controls (N = 801,629,
case fraction s = 0.015631). Rashkin reports odds ratios and P values without
standard errors, so standard errors were reconstructed as
se = |log OR| / |Φ⁻¹(P/2)|. Two numerical issues were handled explicitly: for a
two-study meta-analysis Cochran's Q has one degree of freedom, so the upper tail
is erfc(√(Q/2)) and not exp(−Q/2) (the latter, used initially, gave a
significant-Q rate of 1.59%; the correct form gives 5.33% against a 5%
expectation, which shows no gross miscalibration but is not positive
evidence that the reconstruction is calibrated); and 1,740 variants underflowed to P = 0, for which a `mlogp` column
was computed in log space and P was floored at 10⁻³⁰⁰.

Landi et al. (GCST010304, 36,760 cases) provides only 76 loci and could not be
used as an MR outcome; it was used solely as the reference set of known loci
(§7). Five further FinnGen R12 outcomes (lung, colorectal, pancreatic, breast,
prostate) were streamed and filtered to the instrument set without writing the
full files to disk (§9).

## 3. Instrument selection and harmonisation

**Strict set (primary).** Instruments were the top cis-eQTL per gene × profile at
P < 5×10⁻⁸. Variants were matched to the outcome by chr:pos with allele
consistency; 26 variants were discarded for ambiguous indel representation.
After harmonisation, 6,707 records were retained (FinnGen round), of which 3,579
met the strict criterion (2,142 unique variants); the corresponding meta-round
figure is 3,556 exposures over 2,141 unique variants. The minimum F statistic was
22.2, so the conventional F > 10 filter was never binding and the strict set is
defined entirely by P < 5×10⁻⁸. Each exposure carries exactly one instrument, so
the primary estimator is the Wald ratio.

**Relaxed set (sensitivity analyses only).** Candidate instruments were taken at
F ≥ 5 and P < 0.05 from the full cis window and LD-clumped with plink2
(`--clump-r2 0.1 --clump-kb 1000`) against a 1000 Genomes phase 3 GRCh38 EUR
reference panel of 525 unrelated individuals. For the meta round this gave 1,601
independent instruments; all 371 candidate exposures had at least one, 354 (95%)
had at least two, and 210 had at least four. 9.8% of candidate variants were not
present in the reference panel and were dropped.

## 4. Mendelian randomization

Primary analyses used the Wald ratio with first-order delta-method standard
errors. We note the analytic consequence explicitly, because it constrains
interpretation: with b = β_out/β_exp and SE = se_out/|β_exp|, the test statistic
z = β_out/se_out does not involve the exposure, so a given variant returns the
same P value in every profile in which it is the lead eQTL. Profile specificity
is therefore reported as the profile containing the instrument, never as a
difference in P value.

Sensitivity analyses on the relaxed set used inverse-variance weighted (IVW),
IVW with multiplicative random effects, and weighted median estimators.
Weighted mode was evaluated and excluded: with a median of four instruments per
exposure it returned 1 significant result out of 226 testable exposures and has
no discriminating power at this scale. Heterogeneity was assessed by Cochran's Q.
Because clumping at r² < 0.1 retains appreciable correlation between instruments
while IVW assumes independence, IVW significance is inflated (157 exposures at
FDR < 0.05 versus 21 in the strict single-instrument set); we therefore report
sensitivity results as concordance between IVW and weighted median rather than as
independent discovery, and state this limitation in the text.

Multiple testing was controlled by Benjamini–Hochberg FDR within each analysis
family, with the family defined before the analysis was run.

**Provenance of the HCC known-locus reference.** Seventy-three loci were
assembled from a GWAS Catalog query and from the higher-powered outcome
publication's own table, ten from the latter. Because that outcome supplies the
MR p-values, loci taken from it are a circular reference to the extent that they
rest on it alone, and the deposited Catalog file retained no study accession with
which to separate them. Provenance was therefore resolved through the Catalog
REST API: the outcome study's own reported associations were retrieved, and every
locus appearing among them was queried for the number of other studies reporting
it. Seven proved to have substantial independent support and three none. Locus
attribution was recomputed against three references — the published list, the
list with the three outcome-only loci removed (primary), and the list with all
ten paper-sourced loci removed as an over-corrected bound — with the machinery
otherwise unchanged. Note that removing known loci can only shrink the background
and therefore can only raise the enrichment, so the over-corrected bound is not a
stricter test.

**Effect-size matching for the locus-class differential.** The down-sampling
model recomputes Wald ratios and BH-FDR from summary statistics and takes no
class label as input, so recovery is a function of full-power |z| and the global
threshold alone and the known-versus-novel differential can only follow from the
two classes' |z| distributions. To quantify that, each unit's recovery frequency
at half power was regressed on log|z| with and without the class label, and novel
units were matched to known units within a pre-registered caliper of 0.20 on
log|z|, with replacement, nearest first; the residual differential was
bootstrapped over matched pairs. Common support was recorded, since the matched
comparison is only interpretable where the two |z| distributions overlap. The
analysis was run at gene and at independent-locus level, with the locus level
pre-specified as primary. A process control required the new implementation to
reproduce the published recovery curve bitwise before any matched result was
read; this failed on first run because the published script draws from a single
generator consumed across all power points while the reimplementation re-seeded
per point, and was resolved by matching the draw sequence rather than by
loosening the tolerance. Design, criteria and the replacement wording were fixed
in a pre-registration written before the analysis.

**Unit of inference.** The family is gene × profile records, and the record level
is the primary unit throughout: the registered down-sampling predictions, the
sequential-release trajectory, both generalisations and the disease × resource
grid were all computed on it. Records are not independent — 3,556 of them carry
2,126 unique variants and 1,195 unique genes, since one variant can be the lead
cis-eQTL for a gene in several activation profiles and for more than one gene —
so the list was recomputed under every other unit the paper uses as a sensitivity
analysis rather than as an alternative primary. Variant, gene and independent
locus (1 Mb single-linkage) were each collapsed twice, once by taking the minimum
p-value in the group and once by Simes combination, the latter because a minimum
over correlated tests is anticonservative when treated as a single test; a
two-stage hierarchical procedure selecting genes by Simes and then records within
selected genes was run alongside. All seven alternatives return a superset of the
record-level gene list, so the published claims are the most conservative
available under any of these choices, and both the independent-locus count and
the known-locus share are stable across them. Full results are given as a
supplementary table.

Steiger filtering used the recovered per-profile exposure sample sizes,
`units.exposure = "SD"`, `units.outcome = "log odds"`, and outcome prevalence
0.014962. All 6,943 records had the correct direction. TwoSampleMR's R² formula
for SD units, R² = 2β²·EAF·(1−EAF), is unbounded and produced values above 1
(maximum 1.027); we therefore additionally report the bounded form
R² = F/(F + N − 2) (range 0.185–0.933) and confirmed that the direction is
unchanged under it.

## 5. Colocalisation

Colocalisation used `coloc.abf` with default priors over the full cis window,
with the eQTL as a quantitative trait (per-profile N as in §1) and melanoma as a
case-control trait (N and case fraction as in §2), matching variants on a
chr:pos:alleles key and supplying MAF from the eQTL data.

The published criterion PP.H4/(PP.H3+PP.H4) > 0.7 is too permissive here, because
the ratio can be supported by two vanishingly small posteriors; we required in
addition PP.H3+PP.H4 > 0.5. Under this criterion the pigmentation and known-locus
categories fall to zero while novel-locus records are retained, indicating that
the added requirement removes noise rather than signal.

Two diagnostics accompany every colocalisation result: the number of variants in
the region, and the physical distance between the GWAS peak and the eQTL peak,
which we report alongside posterior probabilities because the posterior is
computed over the whole region and can be high when the peaks are tens of
kilobases apart.

For the meta outcome, colocalisation was additionally repeated restricted to
variants present in both contributing studies (`n_studies == 2`), to test whether
mixed per-variant precision distorted the posterior; each region retained
2,400–4,000 variants under this restriction.

`coloc.abf` permits at most one causal variant per trait per region. Because
posterior movement between H3 and H4 can also be produced by an unmodelled second
causal signal, this assumption was tested rather than assumed. Since the exposure
data are identical across outcome rounds, the test was restricted to the outcome
side, where the number of independent signals was assessed in two ways: FinnGen's
published per-endpoint SuSiE fine-mapping, computed with in-sample LD; and
`susie_rss` run by us on 1 Mb windows using an LD matrix from 525 unrelated
1000 Genomes GRCh38 EUR samples (PLINK 2 `--r-unphased square`), with
`estimate_s_rss` reported as an LD-mismatch diagnostic and the permitted number of
signals varied between 10 and 20. The eQTL side was **not** fine-mapped: with
85–100 donors and no in-sample LD, SuSiE results there are not reliable, and this
is stated as a boundary of the analysis rather than worked around.

The MC1R region served as the calibration control, because FinnGen's in-sample
analysis provides an independent answer there. Our proxy-LD procedure recovers
more credible sets than the in-sample analysis (19 versus 3, without saturating as
the permitted number is raised, and with elevated mismatch diagnostics), so its
counts are interpreted in one direction only: exactly one credible set is treated
as conservative evidence of a single signal, while a count of several is not
treated as evidence of multiple signals.

### 5b. Sequential-release power trajectory

Five FinnGen releases carrying an identical melanoma endpoint (R8–R12) were
analysed with the exposure data, instrument set, allele orientation,
harmonisation, estimator, multiple-testing family, locus annotation and 1 Mb
locus-merging rule held fixed, so that the systematic difference between runs is
case number. Summary statistics were streamed and filtered to the instrument
positions without storing the full files. Multi-allelic positions were resolved by
requiring both effect and other alleles to match the instrument, and records
without a match were dropped. The primary analysis was restricted to variants
present in every release, with each release's own coverage reported as a secondary
analysis. Because releases are nested, agreement between them exceeds that
expected between independent studies of equal size, which makes the design
conservative with respect to instability. Predicted recoveries were generated by
down-sampling R12 (2,000 replicates) and registered before any earlier release was
examined; the registration, its kill criteria and its judgment table are provided
as a supplementary document. Because the reference list at FDR < 0.05 contains no
novel-locus gene, the stratified comparison used a pre-specified relaxed reference
(FDR < 0.20), and the absence of novel-locus genes at FDR < 0.05 is reported as a
result in its own right.

### 5c. FinnGen R13 transfer test

The release following R12 replaced `C3_MELANOMA_SKIN_EXALLC` with
`C3_MELANOMA_SKIN_WIDE` (6,226 cases, 372,159 controls) and the hepatocellular
endpoint likewise with `C3_HEPATOCELLU_CARC_WIDE` (1,070 cases). Endpoint
definitions were compared code by code in Risteys rather than read from the
release manifest, whose one-line phenotype description ("including Hilmo") had
led us to the wrong conclusion at first: the two melanoma endpoints select cases
with the same codes (C43, 172, C44 plus melanoma morphology) and differ in the
rule excluding cancers from the controls (`C3_CANCER_WIDE_EXALLC` against
`C3_CANCER_WIDE`). The `_EXALLC` suffix occurs zero times in the R13 manifest, so
the endpoint used in the trajectory has no same-name successor, and the bare
`C3_MELANOMA_SKIN` endpoint was not substituted because R13 defines it more
narrowly still.

R13 was therefore excluded from the sequential-release trajectory of §5b, whose
comparability gate requires verbatim-identical endpoints, and was analysed as a
transfer test with everything except the outcome file held fixed as in §5b. Four
exposure-by-outcome cells were run — the dynamic CD4⁺ T cell and whole-blood
eQTL exposures against each of the two diseases — each paired with an R12
comparator computed through the same code path from the R12 summary statistics
rather than from any previously stored result, because the stored eQTLGen outcome
columns were computed against the melanoma meta-analysis and not against R12.
Locus attribution used each disease's own known-locus list as in §7, with a
mismatched-list negative control. Three process controls were required before any
R13 number was interpreted: exact reproduction of the R12 candidate list through
the new code path, presence of the large-effect MC1R-region genes, and instrument
coverage of at least 95% of the R12 set. Six point predictions with intervals,
the reading table, the direction of every difference between the releases, and
the results register are provided as a supplementary pre-registration document,
registered before any R13 association statistic was read. Instrument-position
extracts from each scanned file were cached so that the analysis can be re-run
without rescanning 3.2 GB of summary statistics.

## 6. SMR and HEIDI

SMR v1.3.1 was run against the same GRCh38 EUR reference panel (525
individuals; 94,979 variants within cis windows) with GWAS summary statistics
trimmed to cis windows and de-duplicated for multi-allelic sites (104,254
variants). BESD files were built per profile, with rsID as the primary key;
chr:pos → rsID mapping succeeded for 90.0% of variants. The `--diff-freq` QC
removed only 0.45% of variants at the default 5% threshold, indicating that the
85–100-donor eQTL sample did not create allele-frequency inconsistency; the
threshold was not relaxed.

The number of variants entering each HEIDI test is reported with the result, as
HEIDI's power depends on it (range 8–20 here).

## 7. Locus annotation and enrichment

The reference set of known loci was the union of lead variants from the three
Landi et al. 2020 GWAS retrieved from the GWAS Catalog REST API and mapped to
GRCh38 via the FinnGen rsID column: 157 loci carrying their own phenotype
labels (melanoma 52, naevus count 43, pigmentation/hair colour 37, remainder
multi-phenotype). Each instrument was annotated with its nearest known locus and
assigned to one of four categories (known melanoma locus, naevus count,
pigmentation/hair, or potential novel locus) using a 1 Mb window.

Enrichment of FDR-significant results in known-locus categories was tested by
Fisher's exact test against the background proportion across all instruments.
Because several genes can share one locus, and records at one locus are not
independent, **the primary analysis counts independent loci**, not records; the
record-level analysis is reported in Supplementary Table S10 only. The
locus-level background proportion is 10.5% (record-level 10.9%). The same test
was applied unchanged to the five other cancers as a control (§9).

## 8. Negative-control phenotype and comorbidity MR

Because vitiligo GWAS with full summary statistics were unavailable (all eleven
catalogued studies lack summary statistics; FinnGen `L12_VITILIGO` has 391
cases), the pigmentation-pathway negative control was FinnGen
`CD2_BENIGN_MELANOCYTIC` (melanocytic naevi, 13,357 cases): a candidate that
affects melanoma but not naevi is unlikely to act through the pigmentation
pathway. Test sensitivity is established by five known pigmentation loci, all of
which do affect naevi. Additional comorbidity outcomes (rheumatoid arthritis,
psoriasis, type 1 diabetes, a composite autoimmune phenotype, and vitiligo with
its case count stated) were analysed identically and are reported as exploratory.

## 9. Cross-cancer analysis

The exposure side was held completely fixed (3,556 exposures, 2,141 instruments)
and only the FinnGen outcome was exchanged, for lung, colorectal, pancreatic,
breast and prostate cancer. Two quantities were compared across outcomes: the
known-locus enrichment of the FDR-significant list (§7, counted by independent
locus), and the effect estimates for the candidate genes.

## 9b. Second-tumour generalisation (hepatocellular carcinoma)

The design, the known-locus reference list, the direction of every test, the
counting unit, the matched-power control and a five-cell reading table covering
success, partial success, uninformativeness and both forms of contradiction were
fixed in a pre-registration written before any HCC outcome record was read
(Supplementary S20). Only §9 of that document, the results register, was
completed afterwards; two deviations from it are logged there and reported below.

*Exposures.* Unchanged from the melanoma analysis in every respect: the same
eight CD4⁺ T cell activation profiles, the same pool of 8,064 exposure records
over 5,437 unique SNPs prior to outcome matching, the same selection thresholds,
harmonisation and estimator. No exposure-side quantity was recomputed. After
matching and harmonisation the analysed strict sets are 3,398 records (2,047
SNPs, 549 independent loci) for HCC-high and 3,686 (2,216 SNPs, 564 loci) for
HCC-low, against 3,556 for melanoma — so the three analyses are of comparable
size on the exposure side.

*Outcomes.* Two power levels. HCC-high: GWAS Catalog GCST90809296, the European
arm of an eleven-cohort meta-analysis (3,748 cases / 1,861,536 controls, GRCh38
harmonised release). HCC-low: FinnGen R12 `C3_HEPATOCELLU_CARC_EXALLC` (947 /
378,749). The published meta-analysis acknowledges FinnGen among its
contributing resources without listing cohorts by analysis arm, so the two
levels are treated as one resource at two power levels rather than as
independent cohorts, and the recovery comparison between them is reported as
descriptive only. Effective sample sizes are 4/(1/N_case + 1/N_control):
14,962 and 3,779 against 49,337 for the melanoma meta-analysis.

*Known-locus reference.* Lead variants for hepatocellular carcinoma
(EFO_0000182) and liver cancer from the GWAS Catalog at p < 5×10⁻⁸, together
with the lead variants tabulated in the outcome publication: 83 rsIDs, of which
73 could be placed on GRCh38 coordinates present in the outcome files. The list
was built and written to disk before the MR step was run. A locus counts as
known if any of its instruments lies within 1 Mb of one of these leads.

*Attribution test.* Independent loci were defined exactly as in §7 — 1 Mb
single-linkage clustering of instrument positions within a chromosome — and the
enrichment of known loci among FDR-significant loci was tested one-sided against
the background proportion over all instrument loci in the strict set. This
matches the melanoma analysis in §7 and makes the two tumours comparable, but it
is not the background named in the pre-registration, which specified matching on
eQTL-p decile and outcome-EAF quintile. Both were therefore run: the matched
version represents each locus by its most significant eQTL record, strata are
deciles crossed with quintiles, and the null is 10,000 stratified resamples of
background loci drawn to the size of the significant set, giving a one-sided
empirical p. Both versions are reported, and the matched version was additionally
applied to melanoma as a cross-check against the value reported in §7.

*Matched-power control.* The down-sampling machinery of §10, calibrated there
against an independently observed lower-power round, was used to reduce the
melanoma outcome to each HCC effective sample size (200 replicates), recomputing
Wald ratios, Benjamini–Hochberg FDR and significant-locus counts under the same
locus definition, and separating known from novel loci. This control was
specified in advance so that a null or weak HCC result could be distinguished
from a power deficit rather than attributed to tumour biology.

*Scope.* The generalisation covers the MR layer only; no colocalisation or
SMR/HEIDI analysis was performed for HCC. It concerns the attribution of
nominated signal, and neither tests nor supports any claim about TPI1 or about
glycolysis, which were excluded from its scope in advance.

## 10. Power–stability simulation

Individual-level outcome data are not available, so power was reduced in summary
space. For a target effective sample size N_sim below the observed N_obs
(N_eff = 4/(1/n_case + 1/n_control)):

    se_sim   = se_obs × √(N_eff,obs / N_eff,sim)
    beta_sim = beta_obs + ε,   ε ~ N(0, se_sim² − se_obs²)

Wald ratios and BH-FDR were recomputed on each of 200 replicates per grid point,
and the resulting FDR < 0.05 gene list was compared with the full-power list by
Jaccard index, sensitivity and precision. BH-FDR was implemented directly (the
local `statsmodels` installation is incompatible with the installed `scipy`).

**Calibration (pre-specified as the condition for reporting the curve at all).**
Reducing the meta outcome to FinnGen's effective size must reproduce the FinnGen
round that had already been analysed independently. It does: 10.2 simulated
versus 10 observed discoveries at FDR < 0.05, and simulated Jaccard 0.52
[0.36, 0.73] containing the observed 0.455.

**Extrapolation above observed power was attempted and abandoned.** Projecting
beyond the observed sample size requires a model of the true effect
distribution. A flat prior (b_true ~ N(b_obs, se_obs²)) overpredicted discoveries
at observed power by 25-fold, and global empirical-Bayes shrinkage underpredicted
them by 100-fold; the true distribution is sparse and no prior we could validate
in these data reproduces the observed count. Both failures were detected by the
pre-specified sanity check that any prior must recover the observed number of
discoveries at observed power. **All curves are therefore reported only up to
observed power and must not be extrapolated.**

Curves were computed for melanoma and for the five other cancers, and for
melanoma additionally stratified by locus class (known pigmentation/naevus versus
novel). Intervals are the 5th–95th percentiles across replicates; because the
reference lists contain 6 and 4 genes respectively, percentages are not
interpreted below the level of the reported intervals.

## 11. Single-cell ICB cohorts

**GSE120575** (Sade-Feldman; 16,291 CD45⁺ cells from 48 melanoma biopsies,
Smart-seq2, anti-CTLA4 and anti-PD1). CD4⁺ T cells were the fixed set of cells in
clusters 1/5/6/12/15 (3,878 cells). Two processing points are essential to
reproduction:

1. *Normalisation.* The distributed expression matrix is already log₂(TPM+1)
   (per-cell sums back-solve to 9.94×10⁵ ≈ 10⁶), so it must be used as
   distributed. An intermediate object in which these values had been passed
   through `NormalizeData()` (i.e. re-normalised) was identified by back-solving
   per-cell sums and by a Spearman correlation of only 0.874 with correctly
   scaled values, and was discarded for all expression-based analyses; it was
   retained only as a source of cell identifiers. The criterion is whether the
   source file is already normalised, not whether a derived object matches it —
   applying the same test to GSE115978 (below) gives the opposite answer.
2. *Samples and response labels.* Patients with two post-treatment biopsies can
   carry different response outcomes for each (e.g. Post_P1 responder,
   Post_P1_2 non-responder), and pre- and post-treatment labels can differ for
   the same patient. The statistical unit is therefore the GEO `patient_tp`
   sample with its own GEO `response` label, requiring ≥20 CD4⁺ T cells; a
   derived table that collapsed these samples was discarded. This yields 9
   responder / 10 non-responder samples pre-treatment and 8 / 20
   post-treatment.

Header parsing of the accompanying annotation file requires care: `nzchar(NA)`
returns TRUE in R, so a leading empty field is not removed by `nzchar` filtering,
which shifts every column label by one. The header is parsed with
`readLines` + `strsplit` with an explicit `!is.na()` filter and a
`stopifnot(n == 16291)` assertion.

**GSE115978** (Jerby-Arnon) was used for replication of correlation structure;
it has no response annotation, so the response comparison cannot be replicated in
it. Its raw matrix is *not* normalised (per-cell sums: median 2.7×10⁵, CV = 1.15,
46-fold range between the 5th and 95th percentiles), so library-size
normalisation followed by log transformation is the correct processing. Only 16
samples (7 treatment-naive, 9 post-treatment) pass the ≥20 CD4⁺ cell threshold.
Detection rates (invariant to any per-cell monotone transform, and therefore
comparable across processing versions) are reported alongside every correlation,
because genes detected in 5–11% of cells (SPSB2, SMC2, ZFYVE19) do not support
reliable per-sample estimates.

**Bulk ICB cohorts** GSE78220 (28 baseline biopsies, 15 responders / 13
non-responders) and GSE91061 (paired pre- and on-treatment) were analysed at the
tissue level, with the interpretive limits set out in the Results.

## 11b. Second-tumour patient stratification (GSE235863)

Design, hypotheses, cell-population definition, score, permutation scheme,
positive controls, reading table and the arm-by-arm minimum attainable P value
were fixed before any expression value was read (Supplementary S21); §0 of that
document lists exactly which metadata had been inspected at the time of writing.

The deposited object contains 191,435 CD45⁺-sorted cells from nine patients
treated with anti-PD-1 plus lenvatinib, in thirty samples spanning paired
pre- and post-treatment tumour and blood. Response is annotated per patient in
the GEO sample records rather than in the object, giving four responders and
five non-responders. Counts were confirmed to be raw integers — per-cell sums
recomputed from the matrix match the deposited total-count field — then
library-size normalised to 10⁴ and log1p-transformed once.

The analysed population is the authors' own `CD4T` annotation excluding the
FOXP3 subcluster, with the Treg-inclusive definition reported alongside because
the direction of the melanoma replication had proved sensitive to that choice;
the score is the locked 16-gene signature, all sixteen genes present, computed as
the mean of per-gene z-scores across samples within an arm, with the
TPI1-removed version reported alongside. Samples required ≥20 CD4⁺ T cells; all
thirty passed with a minimum of 399. Inference is a one-sided exact permutation
of the patient-level response label with the direction taken from the discovery
cohort, evaluated over all label assignments, with the two-sided Wilcoxon test
and the effect size reported alongside. The primary test is the post-treatment
tumour arm; the remaining arms are secondary or exploratory and
Benjamini–Hochberg corrected.

Because the permutation P value has a combinatorial floor of 1/C(n, n_R), each
arm's minimum attainable P was computed from the sample structure in advance:
0.0143 (post-treatment tumour, 4 versus 4), 0.0667 (pre-treatment tumour, 4
versus 2 — unable to reach 0.05 at any effect size), 0.0079 (post-treatment
blood) and 0.0286 (pre-treatment blood). Positive controls, all pre-specified
with direction, were a higher glycolytic score in myeloid cells than in naive
CD4⁺ T cells, PTPRC detection in ≥90% of cells, and the expected polarity of
CD8A/CD8B and CD40LG/IL7R between the annotated CD8 and CD4 compartments. A
leave-one-sample-out check is reported and labelled post-hoc.

## 12. Module scores, response testing and lineage control

Module scores were computed with Seurat `AddModuleScore` (expression-binned
control genes). Gene sets: glycolysis (22 enzymes and transporters: SLC2A1,
SLC2A3, HK1, HK2, HK3, GPI, PFKL, PFKM, PFKP, PFKFB3, PFKFB4, ALDOA, ALDOC,
TPI1, GAPDH, PGK1, PGM1, PGAM1, ENO1, ENO2, PKM, LDHA — analysed with TPI1
excluded wherever the module is used as evidence independent of TPI1); a reduced
version restricted to the 12 enzymes detected in ≥30% of cells is reported where
stated; proliferation (MKI67, TOP2A,
CCNB1, CDK1, PCNA); exhaustion (PDCD1, CTLA4, LAG3, HAVCR2, TIGIT, TOX); and
OXPHOS as a specificity control.

We record one limitation of provenance. The gene set used for the exhaustion
score in an early round of analysis was not preserved, so the recomputed
analyses use the standard set listed above. Exhaustion-related numbers are
therefore not strictly comparable with the pre-correction versions, and we
report only the recomputed values; the proliferation set was recorded and is
unchanged. No conclusion in this paper rests on the exhaustion score alone.

Response comparisons used two-sided Wilcoxon rank-sum tests on **sample-level**
means (donors, not cells, are the statistical unit) with BH-FDR within timepoint
across the full variable family. Correlations between genes and module scores are
reported at the sample level; cell-level correlations are reported as descriptive
only, and always with partial correlation on detected-gene count, because
abundant genes are correlated through shared detection depth.

Because a module score of 22 genes is diluted by enzymes detected in 1–9% of
cells, the pathway was additionally tested enzyme by enzyme with the identical
procedure. Directional consistency across enzymes was originally assessed by a
binomial test against 1/2; that test was withdrawn during revision because the
enzymes are strongly correlated (median pairwise Spearman 0.48, effective number
≈ 11.7 of 16), making it anticonservative by two to three orders of magnitude.
Enzyme-direction counts are therefore reported descriptively, and inference uses
the label permutation that moves the whole signature together and so preserves
the correlation structure exactly.

**Lineage-purity control.** Splitting cells by a score also splits them by
cell-type purity, so any within-population split was controlled at the **cell**
level, not the cluster level. A CD8ness score (CD8A, CD8B, GZMK, NKG7, CCL5,
GZMB, PRF1, KLRD1, GNLY) and a CD4 score (CD4, IL7R, CD40LG, MAL, LTB, TRAT1,
ANXA1, CCR7, AQP3) were used to (i) build lineage-marker-free versions of every
module tested, (ii) residualise scores on CD8ness, and (iii) match high- and
low-score cells on CD8ness and sequencing depth, with residual standardised mean
differences reported for both matching variables. Ties in score-based splitting
were not broken at random (strict > 60th percentile).

## 13. Spatial transcriptomics

Legacy spot-based spatial data (Thrane et al., four patients × two serial
sections = eight sections, 100 µm spots) were used to test compartment
attribution. Spot coordinates are encoded in the column names (`2x9` = x2, y9);
row names are `SYMBOL ENSG…` and were split on whitespace. Spots were retained at
≥500 total counts and ≥200 detected genes (2,317 of 2,345). Counts were
normalised to CP10K and log1p-transformed. Compartment scores (tumour/melanocyte,
lymphoid, B, myeloid, stromal, glycolysis excluding TPI1, proliferation) were
computed as the mean of marker-set genes after **within-section** z-scoring.
Spearman correlations were combined across sections by Fisher z with weights
n − 3, and the number of sections in which the correlation is positive is
reported alongside, because a meta-analytic estimate can be produced by one
section. Per-section joint regressions of the gene on tumour score, lymphoid
score and log sequencing depth are reported with standardised coefficients. Spots
were assigned to a compartment when the within-section z score of one compartment
exceeded the next by ≥0.25, and otherwise labelled mixed.

The pre-specified positive control was HLA-C (an antigen-presentation gene
expected to track the immune compartments) with SMC2 (expected to track
proliferation); the compartment analysis is interpretable only because these
behave as expected in the same sections. A within-lymphoid-subset test was run
and **failed its own positive control** (HLA-C P = 0.23), so that subset is
underpowered and neither its null for TPI1 nor any other result from it is
reported as evidence.

## 14. Purified CD4⁺ T cell multiome (GSE282266)

Single-cell multiome (10x, RNA + ATAC) of purified CD4⁺ T cells stimulated with
anti-CD3/CD28 across rest, 2.5 h, 5 h and 15 h, in four sets (40,495 nuclei for
the 15 h axis analyses). The RNA modality of 10x multiome is snRNA-seq: MALAT1
dominates and GAPDH is detected in 18.7% of nuclei, which is expected; per-cell
detection rates are therefore not compared with scRNA-seq, and RNA-level
comparisons across samples are made on per-sample pseudobulk.

The activation positive control (IL2, IL2RA, IFNG up; LEF1, KLF2, TCF7 down) was
pre-specified with a pass threshold, and passed (7/12 markers by a monotone test;
immediate-early genes peak at 2.5 h and are not captured by a monotone test,
which was stated in advance).

The glycolysis axis was defined by residualising a 16-gene glycolysis module
(SLC2A1, SLC2A3, HK1, HK2, GPI, PFKL, PFKP, PFKFB3, ALDOA, TPI1, GAPDH, PGK1,
PGAM1, ENO1, PKM, LDHA — the enzymes measurable in this snRNA system, a smaller
set than the 22 used in §12) on an eleven-gene activation module (IL2RA, CD69,
TNFRSF9, TNFRSF4, ICOS, BATF, NFKB1, REL, IRF4, MYC, IL2) and on sequencing
depth, with a pre-specified kill criterion: if
R²(glycolysis ~ activation + depth) > 0.70 the axis is not separable from
activation and the analysis stops. Observed R² = 0.024. Cells were then matched
on activation decile × depth quintile before splitting on the residual, and genes
detected in ≥5% of cells in each set were ranked by the split; only genes with
concordant direction in 4/4 sets were carried forward.

Two controls accompany this: (i) a nuclear-retention index (MALAT1 + NEAT1
fraction) with a pre-specified kill criterion (|ρ| > 0.4 with the residual axis,
or > 0.5 of variance removed) to exclude a nuclear/cytoplasmic composition
artefact — observed 0.047 of variance; (ii) module detection rates, checked
**before** interpreting any module-level null. An earlier round of pre-specified
immune modules was withdrawn on this basis: Th2 (median detection 2.4%),
Th17 (4.4%) and FOXP3 (4.1%) cannot be measured in this system, and testing them
is uninformative rather than negative.

## 15. Chromatin accessibility and motif enrichment

Peak sets differ between samples (94,000–134,000 peaks per sample), so
cross-sample comparison must be performed on genomic intervals. Peaks from the
eight samples were merged into a consensus set of 199,740 intervals, and each
sample's counts were mapped onto it with the global index preserved. This is
stated because two earlier implementations failed here: exact peak-name matching
gave an empty intersection, and a subsequent version destroyed the global index
with `reset_index(drop=True)`, mapping every chromosome to chr1. Chromosome
distribution of top-ranked peaks is now a standard diagnostic output.

Differential accessibility was computed for the glycolysis axis and, separately,
for activation (rest versus 15 h). Motif enrichment used JASPAR2020 CORE
(species 9606) with `motifmatchr` against BSgenome.Hsapiens.UCSC.hg38, comparing
foreground peaks with background peaks matched on GC content (5% bins) and
log accessibility (10% bins), by one-sided Fisher's exact test with BH-FDR.

Because axis and activation log₂ fold-changes correlate at r = 0.647 (R² = 0.419)
at the chromatin level — despite near-orthogonality at the RNA level — motif
analysis was repeated restricted to activation-invariant peaks
(|activation log₂FC| ≤ median; 20,068 peaks, retaining 68% of the axis
amplitude). Both versions are reported. Within-family motif similarity means
these analyses identify families, not individual factors, and results are
reported at family level.

## 16. Peripheral blood variance decomposition (GSE199994)

Whether the CD4 glycolysis programme is an individual-level trait was tested
before any biomarker analysis, as an explicit precondition. PBMC multiome from
eight baseline melanoma patients and two healthy donors was processed to 21,213
CD4⁺ T cells (T-cell score positive, CD4 panel score above CD8 panel score).
Intraclass correlation of the per-cell, depth-residualised module score across
donors was compared with (i) an upper-bound control (XIST, RPS4Y1) and (ii) a
floor built from 50 random gene modules matched on set size and detection rate.

Two points are recorded because both changed the answer: gene-level z-scoring
**within** each sample centres every donor at zero and destroys the between-donor
variance that ICC is meant to measure (detected by the upper-bound control:
XIST ICC = 0.0006); standardisation must be performed once across all donors
pooled. And the random-module floor must be recomputed **within the same
subgroup** as the comparison, because the ten-donor floor is not the correct
reference for a patients-only ICC.

The pre-specified downstream analysis conditional on this test is archived
unexecuted (`manuscript/GSE199994_prespecified_design.md`), labelled
"precondition not met".

## 17. Public eQTL resource survey

Replication resources were surveyed across the eQTL Catalogue (758 datasets, of
which 58 are stimulated T cell datasets) and DICE. Two query behaviours matter:
gene-based queries return position-sorted, paginated results, so the minimum P
value on one page is meaningless and queries must be made by rsID; and a gene
absent from a dataset returns HTTP 400 rather than an empty result, which is not
a null result but a quantification decision by the resource. Datasets were also
checked for identity with the exposure data by comparing design and per-timepoint
sample sizes, which excluded one resource as circular.

## 18. Pre-specification, positive controls and stopping rules

The following rules were fixed before the analyses they govern and are reported
because they determined which results are presented:

1. **Every test carries a positive control, and a test whose positive control
   fails is discarded rather than interpreted.** Instances: HLA-C in the spatial
   compartment analysis; AP-1 enrichment in the activation contrast; XIST and a
   matched random-module floor in the variance decomposition. The
   within-lymphoid-subset test (§13) and the pre-specified immune-module analysis
   (§14) were discarded under this rule.
2. **Kill criteria, test direction and window sizes were written into the script
   before it was run** (R² > 0.70 for activation confounding; |ρ| > 0.4 for
   nuclear contamination; ±2 kb as the primary chromatin window, with the ±10 kb
   result explicitly not used because it was post hoc and driven by two genes).
3. **A negative result was not followed by a search for a positive one** in the
   same data.
4. **Simulations require calibration against something not used to build them**:
   the power-stability simulation was required to reproduce an independently
   measured quantity before any curve was reported (§10). We call this *empirical*
   rather than *external* calibration, because FinnGen contributes to the meta
   outcome that is being down-sampled, so the target is a separately observed
   lower-power round rather than a fully independent dataset.
5. **Sanity checks are designed to catch our own errors**, and did: the
   extrapolation attempts (§10), the chromosome distribution of top peaks (§15),
   and the back-solved per-cell sums (§11).
6. **Donors or samples, never cells, are the statistical unit** for inference;
   cell-level statistics are descriptive.

## 19. Software and computing environment

R 4.4.1 with Seurat 5.5.1 (SeuratObject 5.4.0), Matrix 1.7.0, data.table 1.16.0,
coloc 5.2.3, TwoSampleMR 0.7.5, susieR 0.14.2, arrow 25.0.0, hdf5r 1.3.12,
TFBSTools 1.42.0, JASPAR2020 0.99.10, motifmatchr 1.26.0, chromVAR 1.26.0,
BSgenome.Hsapiens.UCSC.hg38 1.4.5, survival 3.8.9, CellChat 2.2.0.9001 and
org.Hs.eg.db 3.19.1. Python 3.12.4 with numpy 2.0.0, pandas 2.2.2, scipy 1.18.0,
pyarrow 25.0.0 and matplotlib 3.11.1. External binaries: SMR v1.3.1
(win-x86_64) and PLINK v2.0.0-a.7.2 (64-bit).

Ordinary least squares and Benjamini–Hochberg FDR were implemented directly in
numpy rather than via `statsmodels`: the installed statsmodels 0.14.4 is
incompatible with scipy 1.18.0 (`scipy._lib._util._lazywhere` has been removed,
so `import statsmodels.api` fails). Analyses are therefore free of that
dependency.

## 20. Data and code availability

All exposure, outcome and functional datasets are public and identified by
accession in the sections above: dynamic CD4⁺ T cell eQTLs; FinnGen releases R8
to R12; the Rashkin pan-cancer GWAS (GCST90011809); Landi melanoma GWAS
(GCST010302–GCST010304); GSE120575, GSE115978, GSE72056, GSE78220, GSE91061,
GSE282266, GSE199994, GSE316760 and GSE300445; the Thrane spatial dataset; and
TCGA-SKCM expression and survival obtained through UCSC Xena.

Analysis code and intermediate result tables will be deposited at
⟨repository DOI⟩. The deposit comprises every numbered analysis script
(`step*.py`, `step*.R`) and figure script, together with the result tables
needed to reproduce each figure and every number reported in the text —
specifically the MR, colocalisation and SMR/HEIDI outputs; the corrected
single-cell patient-level tables that supersede the pre-correction versions;
the cross-cancer and locus-attribution tables; the power-stability and FinnGen
release-trajectory outputs with their registered predictions; the multiome axis
and motif tables; and the compartment-attribution tables. The pre-registration
document and its modification log are included as deposited files rather than
as text in this paper.

Two categories are deliberately included even though they support no claim: the
intermediate outputs of analyses we withdrew, and the outputs of tests that
failed their own positive controls. Both are labelled as such in the deposit, so
that the tally of attempts in Supplementary S12 can be checked against the
files rather than taken on trust.
