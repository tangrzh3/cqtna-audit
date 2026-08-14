<!-- GENERATED FILE - do not edit by hand.
     Sources: MANUSCRIPT_v2_dual_thread.md + METHODS_draft.md
     Regenerate with: python assemble.py   (last built 2026-08-14) -->

# When a pathway can be interrogated genetically, and which gene gets named: an audit of dynamic-eQTL target nomination in melanoma

---

## Abstract

Dynamic expression quantitative trait loci (eQTLs) mapped across immune-cell
activation, combined with Mendelian randomization (MR), are increasingly used to
nominate immune targets in cancer. We applied this framework to melanoma —
CD4⁺ T cell cis-eQTLs from eight activation profiles as exposures, a 12,530-case
meta-analysis as the outcome — and added the checks needed to defend a target
claim.

Nomination proved highly sensitive to the outcome GWAS used. Significant signal is
dominated by loci already known for the outcome (4.1-fold by independent locus);
the pattern holds across five nested power levels of one GWAS resource, where no
novel-locus gene reaches significance at any case number between 2,705 and 5,753,
and it transfers unchanged to that resource's current release, whose test
statistics correlate with the previous ones at 0.94 with standard errors 3.85%
smaller;
it recurs when the exposure side is held fixed and hepatocellular carcinoma is
substituted — three of four significant loci are known HCC loci, including
*PNPLA3*, though at 30% of melanoma's effective sample size that pre-registered
test does not itself reach significance. **The mirror-image test gives the same
answer**: holding the melanoma outcome fixed and replacing the exposure entirely
with a whole-blood eQTL resource 300-fold larger — which changes donors, cell
composition and platform together, and so identifies a resource rather than a
sample size — raises significant loci from 7 to 30 and leaves the attribution
unchanged (66.7% known-locus, 4.44-fold, P = 3.7×10⁻¹¹, against 4.09-fold in
CD4⁺ T cells). Replacing the outcome with a higher-powered meta-analysis — which
changes study composition as well as power, so the two cannot be separated —
increased MR discoveries while *lowering* colocalisation support and replaced the
candidate list entirely: two lists from identical exposure data share no genes. Down-sampling, validated against twelve
pre-registered predictions, shows the loss falls unevenly by effect size and
therefore by locus class — at half
of observed power, known-locus genes are recovered about four times as often as
novel-locus ones. Within this project's fully enumerated
target-substantiation record — eleven tests, not a sample of any framework's
failure rate — seven overturned the claim under examination, three were
inconclusive, and one was significant on discovery but not confirmed on transfer. **In this analysis, and within the power range we could
observe, the reproducible part of the candidate list was the part that did not
constitute a discovery.**

What survives is a constraint on what can be interrogated, and it must be stated
at the level it applies to: of 28 glycolytic genes, **3 are instrumentable in this
exposure resource, 2 remain analysable against this outcome, and 1 yields a
nominal association**. Only the first is a property of the pathway and the
exposure data; the others are joint properties of the pathway and the particular
outcome GWAS, and compressing them would attribute an outcome-side limit to the
biology. Meanwhile the nominated gene is far higher in malignant cells than in
CD4⁺ T cells in all three cohorts examined (6.7-fold, 2.4-fold, and concordant in
direction), so tissue-level validation of it measures tumour rather than T cells.
A patient-level stratification hypothesis nominated by the same programme was
carried from discovery into two pre-registered follow-ups, and no treatment arm
was confirmed in more than one of them.

Dynamic eQTL data are more informative about **what can be interrogated
genetically, and when** than about **which gene should be targeted**.

**Keywords** Mendelian randomization · context-specific eQTL · target nomination ·
colocalisation · statistical power · reproducibility · compartment attribution ·
CD4⁺ T cell metabolism · melanoma

---

## 1. Introduction

Most cis-eQTLs are context dependent. Loci mapped in resting bulk tissue or in
unstimulated primary cells cannot capture regulatory variants that act only while
a cell responds to a stimulus, and much immune-relevant regulation falls into
that category. Dynamic eQTL datasets, in which primary immune cells are profiled
across a stimulation time course, were built to close that gap. Their combination
with Mendelian randomization has become an attractive route to causal target
nomination: germline genotype is fixed before disease, cis-eQTL effects carry a
clear directional prior, and the exposure is measured in the cell type through
which the effect is proposed to act. Several recent studies apply exactly this
design — activation-time-course CD4⁺ T cell eQTLs as instruments, cancer GWAS as
outcomes — and report novel immune targets on that basis.

The inference passes through more joints than its summary statistics reveal. An
instrument is selected in one cell state at one timepoint and usually consists of
a single variant, so the association's P value is supplied entirely by the
outcome GWAS. Whether the eQTL and the disease signal share a causal variant,
rather than lying in linkage disequilibrium with a large neighbouring effect, is
decided by colocalisation, whose resolution depends on outcome power. Which of
two genes at a locus carries the effect is not decided by MR at all. Neither is
the cell type in which the effect operates, because functional replication is
usually available only at tissue level, where a broadly expressed gene reports on
whichever compartment dominates. Each joint has a failure mode that is invisible
unless specifically tested.

Melanoma is an unusually informative setting for examining those failure modes,
for two reasons that pull in opposite directions. Its common-variant architecture
is dominated by pigmentation and naevus loci, several with effects far larger
than anything expected from immune regulation; MC1R alone exceeds P = 10⁻⁵⁹ in
current meta-analyses. That architecture is a built-in end-to-end positive
control — a pipeline that fails to recover it is not working — and simultaneously
a built-in confounder, because long-range LD around such loci can present as a
causal association for any neighbouring gene, including genes with plausible
immune functions. Melanoma is also where CD4⁺ T cell biology has the most direct
clinical relevance, through immune checkpoint blockade, so candidates can be
examined against response-stratified single-cell data rather than annotation
alone.

This paper has two parts, and the division between them is not rhetorical. **Part
I asks what this framework can nominate, and finds that everything depending on
the outcome GWAS is unstable.** Every check we added in order to *strengthen* a
target claim that returned a verdict returned a negative, seven times of eleven
tests, with three inconclusive and one significant on discovery but unconfirmed on
transfer (Box 1). **Part II asks what the same data
still establish once nomination is set aside, and finds that the exposure side
supports claims the outcome side cannot.** Because a single-instrument Wald ratio
has test statistic z = β_out/se_out, the outcome GWAS supplies every MR
significance claim in Part I, and colocalisation depends on it as well as on the
exposure signal, the LD reference and the region and priors chosen. None of the
observations in Part II uses it. The two parts are therefore separated along the
exposure/outcome axis, and Part II does not inherit the *outcome-side* instability
of Part I. It is not, however, independent of the outcome in a stronger sense:
TPI1 and the glycolytic pathway are studied **because** the MR pointed there, so
Part II is outcome-nominated and independently evaluated, not independently
generated. It does inherit the exposure side's own limits —
85–100 donors, winner's curse on lead variants, differential measurability across
timepoints — and, for the patient arm, cohorts of nine to twenty patients.

---

# PART I — The reliability boundary of target nomination

## 2.1 Design, positive controls, and the constraint that organises everything

We instrumented gene expression in CD4⁺ T cells using cis-eQTLs from eight
activation profiles (naive and memory cells at 0 h, 16 h, 40 h and 5 d after
anti-CD3/CD28 stimulation; 85–100 donors per profile) and tested their causal
association with melanoma. After harmonisation and restriction to genome-wide
significant instruments (P < 5×10⁻⁸; minimum F = 22.2, so the F > 10 filter never
bound), 3,556 gene × profile exposures carried exactly one instrument and were
analysed by Wald ratio.

**The single-instrument constraint is stated first because it organises the whole
paper.** For a Wald ratio, b = β_out/β_exp and SE = se_out/|β_exp|, so
z = β_out/se_out does not involve the exposure. The same variant returns an
identical P value in every profile in which it is the lead eQTL — verified for
CDK10, where memory 5 d and 40 h give OR 0.788 versus 0.584 with P = 2.50×10⁻³⁰
in both. Three consequences follow, and they run in different directions:

1. Cell-type and temporal specificity **cannot** be argued from differences in
   MR P value. The admissible statement is which profile contains the instrument.
2. **Given the instrument set, everything in Part I** — which records reach
   significance, what the list contains, how colocalisation resolves — **is a
   function of the outcome GWAS.** The qualifier is not decoration: the exposure
   resource fixes how many instruments exist to be tested at all (§2.2), so the
   outcome determines which of a supplied set is named, not how large that set is.
3. Conversely, **where an instrument appears in the time course is a property of
   the exposure data alone**, and no change in outcome power can alter it. This
   is what Part II rests on.

Of the three studies we identified using this framework, none states this
consequence in its main text; one refers to the single-variant setting without
drawing the inference (Supplementary Table S10).

**Which unit the inference belongs to, and what happens if it is changed.** The
same constraint creates a units problem that has to be settled before any count
below can be read. The 3,556 records carry only **2,126 unique variants and 1,195
unique genes**, because one variant can be the lead eQTL for a gene in several
profiles and for more than one gene; the multiple-testing family is therefore
records, while the candidate list is reported by gene and the attribution result
by independent locus. Records are therefore **not independent**, and
Benjamini–Hochberg over them carries no clean guarantee under that dependence.

We resolve this by separating the two jobs the record level had been doing at
once. **The record-level analysis reproduces the nomination pipeline under
audit** — it is what the studies being examined run, and what the registered
down-sampling predictions, the release trajectory, both generalisations and the
grid were computed on — so it is left exactly as published, in that role.
**The independent locus is the unit in which this audit states its own
conclusions.** That division is not a convenience chosen after the fact: we
recomputed the FDR < 0.05 list under every unit the paper uses — variant, gene
and independent locus, each by both minimum-p and Simes combination, plus a
two-stage hierarchical procedure selecting genes and then records within them
(Supplementary S26) — and **the quantity the audit argues from barely moves**
(7 to 9 independent loci against the 7 reported; 2 under every single unit in the
FinnGen round), while **the attribution conclusion is unit-independent** (the
known-locus share of significant loci ranges 42.9–55.6% against the 42.9%
reported, and is 100% at every unit in the FinnGen round). The audit conclusion
therefore does not rest on the contested unit at all.

The record level is also the most conservative of the eight — no gene is lost
under any other unit and between 1 and 43 are added, so every significance claim
here is a subset of what a looser unit would license. ⚠ **That is an empirical
statement about list length, not a demonstration that FDR is correctly controlled
under this dependence structure, and we do not offer it as one.**

One row of that table is worth reading as a finding rather than a check. Testing
at the locus level leaves the significant-locus count almost unchanged but
inflates the gene list from 10 to 53, because a significant locus does not name a
gene — **22 of those 53 sit in the MHC and a further 7 in the chr17q21.31
inversion**, the two regions this paper already identifies as attribution traps.
That is finding ② arriving through a
different door, and it is the reason we count independent loci rather than gene
records wherever a count carries an argument.

Steiger filtering is near-deterministic at this scale: all 6,943 harmonised
records had the correct direction, with R²_exposure (median 0.152) exceeding
R²_outcome (median 1.5×10⁻⁵) by four orders of magnitude. TwoSampleMR's R²
formula for SD units is unbounded and returned values above 1 (maximum 1.027); we
report the bounded form R² = F/(F + N − 2) (0.185–0.933), under which direction is
unchanged. We therefore treat Steiger as procedure, not evidence; the directional
argument rests on the biology of germline cis-eQTLs.

The outcome combined FinnGen R12 (5,753 cases) with Rashkin et al. (6,777 cases)
by fixed-effect inverse-variance meta-analysis: 12,530 cases and 789,099
controls. It behaves as expected — all 157 known melanoma loci recovered and 136
(86.6%) strengthened; genome-wide significant variants rose from 2,871 to 4,552.
**That recovery is an integrity check on the outcome data and the meta-analysis
step, not a positive control on the nomination pipeline**: it involves no
instrument, no MR and no colocalisation, and a pipeline could pass it while
failing at every joint tested below. The end-to-end controls are the ones in the
next paragraph.
Rashkin standard errors were reconstructed from odds ratios and P values;
Cochran's Q was significant for 5.33% of variants against a 5% expectation, which
shows no gross miscalibration but is not positive evidence that the reconstruction
is calibrated.

Positive controls passed at every layer. PARP1, a known melanoma locus,
colocalised in the FinnGen round (PP.H4 = 0.92 and 0.96) and was confirmed by
SMR/HEIDI, establishing that the LD panel, rsID mapping and BESD construction all
function.

## 2.2 Finding ①: the significant signal sits on loci already known for the outcome

In the FinnGen round, all 10 FDR-significant records fell within 1 Mb of a known
pigmentation or naevus locus (VPS9D1-AS1 50 kb from MC1R; CDK10 48 kb from the
MC1R R151C variant; PARP1 5–14 kb from PARP1). There was no immune signal at
FDR < 0.05 at all.

Under the meta outcome, counting independent loci — records at one locus are not
independent, and counting by record inflates enrichment — 3 of 7 melanoma loci
were known-locus loci, a 4.09-fold enrichment over the 10.5% background
(P = 0.028) (Fig 1, Fig 6a). The record-level figure (12/21, 5.26-fold,
P = 2.7×10⁻⁷) is reported in Supplementary Table S10 only.

Applying the identical exposures to five other cancers, melanoma has the largest
enrichment and the two cancers with no expected relationship to pigmentation sit
near background (breast 1.66×, P = 0.21; colorectal 1.59×, P = 0.49; pancreas no
FDR-significant records). The comparison does not deliver clean specificity —
prostate is also nominally enriched (2.38×, P = 0.031), lung is borderline
(2.86×, P = 0.077), and no outcome survives correction across six tests — so we
report it as qualified support.

This result was then tested twice more, once by holding the disease fixed and
moving power, and once by holding everything fixed and moving the disease.

**Five nested releases of one resource.** FinnGen's sequential releases separate
power from dataset: five carry the identical melanoma endpoint (R8 2,705 cases;
R9 2,993; R10 3,194; R11 3,932; R12 5,753), with exposure data, instruments,
harmonisation, estimator, testing family and locus definition all held fixed.
Predictions were generated by down-sampling R12 and **registered before any
earlier release was examined**; all twelve — discovery counts, overall recovery,
and the known-minus-novel differential at each release — fell inside the
registered intervals (Fig 7, Supplementary S9; the candidate list at each release
is Supplementary Table S11). The known-locus part of the list
is recovered early and then does not move: the same five genes constitute the
FDR-significant list at R9, R10 and R11 without change, four already present at
R8 at 47% of reference power. And **no novel-locus gene reaches FDR < 0.05 in any
release** — not at 2,705 cases, not at 5,753. Because releases are nested,
agreement between them exceeds that between independent studies of equal size, so
this is an upper bound on stability; and because no novel-locus candidate exists
anywhere in the range, the trajectory shows only that they do not appear, not
that they turn over.

**Why the series stops at R12, and what the release after it shows instead.** The
next release, R13, retired the endpoint used above and replaced it with a
successor whose cases are defined by exactly the same codes but whose controls are
screened by a different cancer-exclusion rule (6,226 cases and 372,159 controls
against 5,753 and 378,749). It therefore fails the comparability gate the
trajectory design imposed in advance — the same gate that excluded FinnGen R6 —
and admitting it as a sixth power level would score a change of phenotype
definition as a change in power. We analysed it instead as a **transfer test**,
registered in full before any R13 statistic was read, including the direction of
each difference: the added cases raise power, the control pool contributes only
1.6% of effective sample size, and a broadened cancer exclusion leaves cleaner
controls, so **every difference pushes towards more and stronger signal**. That
asymmetry was recorded precisely so that a successful transfer would count as
partly confounded while a failed one could not be blamed on the definition change
(Supplementary S25).

The list transfers intact. At 6,226 cases the FDR-significant set is the same six
genes at the same two independent loci, still with no novel-locus gene, and locus
attribution is 9.6-fold (P = 0.011); the three remaining exposure-by-disease
pairs move the same way (3.7-, 17.6- and 9.8-fold). **The stability is not an
artefact of reusing data.** The two releases' test statistics correlate at
|z| r = 0.937 rather than being identical; standard errors shrink by a median of
3.85%, against the 3.8% the added cases predict, which is the internal check that
R13 really is the higher-powered computation and not a copy; and 51 of the 275
nominally significant records are replaced. What holds still is the FDR-significant tier;
the tier immediately below it turns over by roughly a fifth in a single release.
Because every difference between the releases points the same way, we do not read
this as an effect of power alone, and because R13 contains R12's participants it
is not an independent replication — it answers only the question a reader
re-running this audit on today's data would be asking.

**A second disease.** If the pattern is a property of the framework rather than
of melanoma, it should recur when only the outcome changes. We repeated the
nomination step against hepatocellular carcinoma with **the exposure side
unaltered to the byte**, at two outcome power levels (a 3,748-case European
meta-analysis; FinnGen R12's 947-case endpoint), with the design, known-locus
reference list, direction of every test and matched-power control registered
before any HCC record was read, together with a reading table that pre-committed
us to reporting a failed or contradictory generalisation as such (Supplementary
S20). Across both levels the FDR-significant set is four independent loci, three
of them carrying a known HCC lead SNP: SAMM50, 26 kb from rs2294915 in the
*PNPLA3* region, and ZNF506 in the chr19 *TM6SF2* region, at both levels. Against
a 5.6% known-locus background the enrichment is 17.6-fold in the lower-powered
outcome (P = 0.0031) and 8.9-fold in the higher-powered one, **where it does not
reach significance (P = 0.110)**. The single novel-locus nomination, SUPV3L1,
appears only at the higher power level and has FDR = 0.98 at the lower one; the
two tumours' novel-locus nominations share no genes.

**One part of that reference list came from the outcome GWAS itself, and we
measured how much (Supplementary S29).** Of its 73 loci, 10 were taken from the
higher-powered outcome publication's own table — asking whether hits computed
against a GWAS fall on loci that GWAS reported is circular to that extent.
Checking each against the GWAS Catalog, **7 carry heavy independent support**
(APOE in 1,399 studies, *TM6SF2* E167K in 1,073, the *PNPLA3* lead in 225) and
**3 rest on the outcome publication alone**. Removing those three changes
nothing at all — 8.85-fold and 17.6-fold to the same decimals — because no
instrument falls within 1 Mb of any of them, so they enter neither the numerator
nor the background. **We report this as a measurement of the circularity (3 of
73, none of them near a significant locus) rather than as a robustness check: by
construction it could not have moved the result.** Two limits belong with it.
Dropping all 10 paper-sourced loci instead *raises* the enrichment (to 10.2- and
20.1-fold), because removing known loci only shrinks the background — so
over-correction flatters this test and should not be cited as a stricter
standard. And the GWAS Catalog carries no cohort-level overlap information, so
this addresses circular *locus attribution* and not sample overlap between the
reference studies and the outcome. The melanoma reference list does not have this
structure at all: it comes from a separate publication rather than from the
outcome GWAS, though sample overlap between that publication and the outcome
meta-analysis cannot be excluded either.

The caution is power, and we measured it rather than asserting it. HCC-high
carries 30.3% of melanoma's effective sample size. Down-sampling melanoma to that
level (§2.4) predicts a median of 2 significant loci [5–95%: 1–6], 1 known and 1
novel; **the observed HCC counts are compatible with that melanoma-derived
prediction interval**, and the non-significant primary test is what a test with
two significant loci is expected to deliver. Compatibility with an interval is
not equivalence: an interval this wide (1 to 6 loci) would also accommodate
outcomes we would have read as a difference, so this establishes that the HCC
result is *not evidence against* the pattern, not that the two tumours behave
alike. Two significant loci cannot distinguish between those readings, and we do
not claim they do. We register this as
**supporting but short of the criterion we set ourselves**: by our own reading
table this is not an established generalisation, though neither is it the
contradictory result the table anticipated. Two further limits were
pre-specified: the larger HCC study cannot be shown to exclude FinnGen, so its
two levels are one resource rather than two cohorts and the recovery comparison
between them is descriptive only; and the exercise covers the MR layer alone.
Nothing in it bears on TPI1 or glycolysis — the three HCC genes are unrelated to
that pathway. What it shows is that a CD4⁺ T cell eQTL panel, held completely
fixed, named *the textbook loci of each of the two diseases we placed on the
outcome side*. Two diseases do not establish that this holds for any disease, and
we do not claim it does; what they establish is that the pattern is not peculiar
to melanoma.

**And a second exposure resource.** Everything above varies the outcome while
holding the exposure fixed, which leaves one alternative open: that the pattern is
a property of *this* eQTL resource — 85–100 donors, CD4⁺ T cells, one laboratory.
We therefore ran the mirror-image test, holding the melanoma outcome byte-for-byte
fixed and replacing the exposure entirely with eQTLGen whole-blood cis-eQTLs
(n = 31,684, a ~300-fold larger eQTL sample, a different cell composition and a
different platform), again pre-registered with its reading table before the
resource was downloaded (Supplementary S22). Instruments rose from 3,556 records
to 12,835 and FDR-significant loci from 7 to 30, as a resource of that size
should. The attribution did not move: **20 of 30 significant loci (66.7%) carry a
known melanoma, naevus or pigmentation lead SNP, a 4.44-fold enrichment over this
resource's own 15.0% background (one-sided P = 3.7×10⁻¹¹)**, against 4.09-fold in
the CD4 analysis — the same fold from two exposure resources that share neither
cell type, sample size, nor platform. The pre-registered matched-background
version agrees (4.14–4.48-fold, empirical P = 10⁻⁴). One qualification is
pre-specified and material: the two resources' novel-locus nominations intersect
in exactly one gene, ZFYVE19 — the boundary of what we predicted (≤1), and
notably the gene that carried the strongest colocalisation in the original
analysis. **"Novel-locus nominations never reproduce" would therefore be too
strong; across exposure resources, one of them did.**

**The two axes crossed.** Varying each axis separately leaves one reading open —
that the whole-blood result is peculiar to melanoma, or that the HCC result is
weak because the CD4 exposure is small. We therefore filled the grid, scoring
every cell against **its own disease's** known-locus list, again pre-registered
(Supplementary S24). All six cells enrich: melanoma 4.09-fold (CD4) and 4.44-fold
(blood); HCC-high 8.85-fold and 5.11-fold; HCC-low 17.6-fold and 8.67-fold. Three
of the six reach P < 0.05 and three do not (P = 0.052–0.110) on numerators as
small as one of two loci, so this is a statement about **direction being
consistent in every cell**, not six independent significant tests. The
pre-registered negative control is what makes it interpretable: scoring the same
cells against the *wrong* disease's list collapses the enrichment (HCC's list on
melanoma, 1.30-fold, P = 0.41; melanoma's list on HCC, 0.00-fold, P = 1.0), so
the effect is specific to each outcome's own genetics.

**That control is narrower than it may appear, and we state its limit rather than
leave it to be inferred.** A mismatched list rules out enrichment on loci that
are indiscriminately dense across diseases; it does not rule out enrichment
driven by a density that is itself disease-specific, since a disease's known loci
and its instrumentable regions can be co-located for reasons unrelated to causal
attribution. The permutation control that does exist here matches significant to
background loci on exposure eQTL-p decile and outcome allele-frequency quintile
(4.36-fold, empirical P = 0.023 for melanoma; 5.45- and 13.7-fold for the two HCC
levels), so instrument strength and allele frequency are controlled — **but not
locus or gene density, which would require a density-matched permutation we have
not run.** The claim we make is therefore that the enrichment is specific to the
outcome's own genetics and is not explained by instrument strength or allele
frequency; **we do not claim that locus density has been excluded.**

One unplanned observation from the grid is worth stating, with a caveat about
what it can attribute. Swapping in the eQTLGen resource multiplied significant
loci in melanoma (7 → 30) but not in HCC (2 → 5, 2 → 3). The natural reading is
that **the exposure resource sets how many instruments exist while outcome power
sets how many of them can reach significance**, which is what z = β_out/se_out
predicts, arrived at here from the data rather than from the algebra.

⚠ **This swap is not a clean manipulation of exposure sample size, and must not
be described as one.** eQTLGen differs from the Soskic data in at least three
respects at once — ~300-fold more donors, whole blood rather than sorted and
stimulated CD4⁺ T cells, and a different platform and processing pipeline. Any of
the three could multiply the instrument count, and this design cannot separate
them. What the comparison does establish is the asymmetry itself: whatever it is
about the exposure resource that changes the instrument count, it changed the
melanoma count fourfold and left the HCC counts almost unmoved, and the only
quantity differing between those two arms is the outcome. **The outcome-side half
of the statement is identified; the exposure-side half names a resource, not a
sample size.**

## 2.3 Finding ②: SMR/HEIDI does not exclude the LD confounding colocalisation identifies

Colocalisation assigned the MC1R-region signals — the strongest MR associations
in the study, reaching P = 4×10⁻³⁷ — to distinct causal variants, with PP.H3
dominant and PP.H4 at or near zero (CHMP1A 0.99; VPS9D1-AS1 0.96; SPATA33
0.92–1.00). That reading is independently supported by FinnGen's own fine-mapping
of this endpoint, computed with in-sample LD, which resolves **three high-purity
credible sets** in the region (log₁₀BF 45.6, 36.6, 19.5). The region genuinely
carries multiple independent causal signals — precisely the configuration in
which a neighbouring gene's eQTL can be tagged by one of them without sharing it.

SMR/HEIDI did not reproduce this. Of 291 records that colocalisation assigned to
distinct causal variants, HEIDI failed to reject homogeneity for 253 (86.9%),
including VPS9D1-AS1 (P_HEIDI = 0.649) and CDK10 (0.086, 0.081) (Fig 2). Only
8–20 variants entered each HEIDI test. Evidence tiers accepting MR + SMR without
colocalisation are therefore unsafe here: they would have reported MC1R LD
spillover as CD4-mediated immune targets.

## 2.4 Findings ③ and ④: power moves the two methods in opposite directions, and the list turns over

Meta-analysis raised MR discoveries from 10 to 21 records while *lowering*
colocalisation support: across the 127 exposures run in both rounds, median
PP.H3+H4 fell from 0.249 to 0.207 and median PP.H4 from 0.124 to 0.095, with only
44.1% improving. Three explanations were tested and rejected: uneven meta
coverage (restricting to variants in both studies moved the median only to
0.202); candidate composition (excluded by the paired design); and winner's-curse
selection (newly entering candidates had *lower* PP.H3+H4, 0.125 versus 0.202).

Posterior support is sensitive to the outcome dataset and its effective
resolution, with movement in both directions (Fig 3). We state it that way rather
than as "power determines colocalisation resolution", because the multiple-signal
sensitivity analysis below could adjudicate only one of the three loci examined.

The posterior movement is not an artefact of the analysis window or the
colocalisation prior. Repeating every reported locus across four windows
(±100 kb to the full cis region) and four values of the shared-causal-variant
prior (p₁₂ from 10⁻⁶ to 5×10⁻⁵) leaves the negative results unchanged in all 16
combinations: the MC1R-region genes and PARP1 fail the criterion under every
setting (PP.H4 ≤ 0.25), and ZFYVE19 passes under every setting (PP.H4 0.86–1.00).
Window size is essentially irrelevant (PP.H4 differs by < 0.01 between ±100 kb and
the full region); what varies is the prior, and it varies only for loci near the
criterion boundary — TPI1 and SPSB2 pass in 8 of 16 combinations, SMC2 and
KIAA0040 in 12 of 16 (Supplementary Table S15). For PARP1, five orders of magnitude more signal
sharpened the peak and moved it 25 kb, to 54 kb from the eQTL peak; PP.H4 fell
from 0.92 to 0.05. CLDN7 (0.62 → 0.05) and ELP5 (0.55 → 0.04) behaved likewise.
In the other direction, the meta peak for ZFYVE19 landed exactly on the eQTL peak
and PP.H4 rose to 0.99.

**The two claims here are of different kinds and are kept separate.** That the
list turns over (④a) is model-free: it is established by comparing lists, and is
reproduced by the down-sampling analysis of §2.5, which recomputes only Wald
ratios and FDR and never invokes colocalisation. That the turnover operates
through a reallocation of posterior mass between H3 and H4 (④b) is a statement
about the colocalisation model, and `coloc.abf` permits at most one causal
variant per trait per region.

We tested ④b directly. Because the exposure data are identical in both rounds,
any posterior change must originate on the outcome side, so the question reduces
to the number of independent outcome signals per region. Our SuSiE procedure,
using an external European LD panel, recovers ≥2 credible sets in the MC1R region
— passing the pre-specified control — but returns 19 where FinnGen's in-sample
analysis returns 3, without saturating as the permitted number of signals is
raised. The procedure over-splits, so its counts read in one direction only: one
credible set is conservative, several are not. Under that constraint, **the PARP1
region carries exactly one credible set at both power levels**, stable when the
permitted number is doubled, with low mismatch diagnostics. The collapse from
0.92 to 0.05 is therefore not an unmodelled second signal. For ZFYVE19 and for the
TPI1/SPSB2 locus, **no credible set is recovered at either power**, so neither can
adjudicate the question and the ZFYVE19 case carries that qualification.

One general fact belongs here. FinnGen fine-maps only 19 regions genome-wide for
this endpoint, all classical pigmentation or melanoma loci, and **none of our
candidate loci is among them**. At this outcome power the candidate list is built
almost entirely on regions the outcome GWAS cannot fine-map.

**The list itself turned over completely.** The FinnGen round produced four
novel-locus genes passing triple validation (PRPSAP2, IMPA1, GCC2, PADI4) and one
passing the full battery (IMPA1). Under the meta outcome none survived — IMPA1's
MR P moved from 9.3×10⁻⁴ to 0.11, PRPSAP2 lost colocalisation support
(PP.H3+H4 = 0.06) — and six different genes took their place (ZFYVE19, SMC2,
KIAA0040, SPSB2, HLA-C, TPI1). **The two lists share no genes**, from identical
exposure data under an identical pipeline.

Two earlier candidate failures belong here because both came from checks added
expecting confirmation. PADI4, supported by an independent breast cancer study
and an ongoing Phase I programme, collapsed under multi-instrument analysis
(OR 1.087, P = 1.7×10⁻³ → IVW 1.033, P = 0.085; weighted median 1.004, P = 0.835
over seven instruments). GDI2 was removed by heterogeneity (Cochran's Q P = 0.01).

**Administering the criticism to ourselves: does the pipeline name the right gene
where the answer is known?** The objection this paper raises against nomination
can be tested on the pipeline that raises it. We took melanoma loci at which a
causal gene is generally accepted and asked whether our own FDR-significant
nomination names it. Of ten such loci reached under the larger exposure resource,
**six name the accepted gene and four do not**, and the four failures are the
informative half. At the MC1R region — the strongest melanoma locus in the genome
— the nomination spans **sixteen genes and MC1R is not among them**. At OCA2/HERC2
it names the pseudogene HERC2P9 rather than HERC2; at *TYR* it names ODF3; at
CDKN2A/MTAP, C9orf66. Where it succeeds it often succeeds cleanly, IRF4 and MX2
each being named alone and correctly.

Two limits on how far this can be pushed. The accepted-gene list is one we
assembled ourselves — fixed before the comparison was run, but not
pre-registered — and ten loci is a small denominator, so the ratio is an
illustration and not an error rate. What it does establish is directional and
does not need precision: a framework that misassigns the gene at the
best-characterised locus in its own disease should not be read as assigning genes
at uncharacterised ones. This is the co-regulation problem of Tambets et al.
observed at the top of the effect-size distribution, and it is why we report
compartment attribution separately from gene attribution — the two fail
independently.

**The reproducible part of the list is the part that is not a discovery.** Finding
④a is an anecdote unless it is shown to be a property of the design, so we
measured candidate-list recovery as a function of outcome power by down-sampling
in summary space, recomputing Wald ratios and FDR, and comparing with the
full-power list (Fig 4). The simulation is calibrated empirically — not
externally, since FinnGen contributes to the meta outcome — against a separately
observed lower-power round: reducing it to FinnGen power reproduces a round
already measured, 10.2 simulated versus 10 real discoveries, with a simulated
Jaccard of 0.52 [0.36, 0.73] containing the observed 0.455. Two attempts to
extrapolate *above* observed power failed, caught by the pre-specified
requirement that any prior reproduce the observed discovery count, so **curves are
reported only up to observed power**. Melanoma recovers 53% of the full-power
list at 5,000 cases and 85% at 10,000, and five other cancers behave the same
way. Stratifying by locus class then dissolves the apparent stability (Fig 4c):
at 10% of observed power, 46.0% of known pigmentation/naevus-locus genes are
recovered against **0.4%** of novel-locus genes; at 50% power, 85.8% versus
22.8%. Reaching 50% recovery requires 2,506 cases for known-locus genes and 8,771
for novel-locus genes — 3.5-fold. Studies using this framework commonly analyse
3,000–8,000 cases, where novel-locus recovery is 6–40%.

**That differential is mostly a statement about effect size, and we tested how
much of it survives conditioning on effect size (Supplementary S28).** Candidates
at novel loci sit against the detection threshold by construction: their
full-power |z| all fall between 3.73 and 4.55, while candidates at known loci
reach 15.99 (medians 4.23 and 11.07). Since recovery under down-sampling is a
monotone function of |z|, part of the gap must follow from that alone — and it
does: **a model using full-power |z| with no class label reproduces 68–80% of it**
(locus and gene level respectively). Matched on |z|, the residual differential is
**+5.2 percentage points [−1.6, +12.0] by independent locus** and +6.4 [+0.9,
+12.0] by gene, against raw gaps of 47.6 and 63.0. The two classes barely overlap
in |z| — only one known locus falls inside the novel range — so the remaining
20–32% cannot be attributed: it is not separable from the failure of matching.
The claim we retain is therefore about threshold proximity rather than about the
category. The practical consequence is unchanged, and does not depend on which
reading is right: **for anyone holding such a candidate list, the novel-locus part
is the part that will not replicate.** What we cannot separate is whether known
loci recover better because they are genuinely larger effects or because novel
loci were selected at the threshold; both are true here, and this design cannot
apportion them.

> In this analysis, and within the range of outcome power we could observe, the
> reproducible part of the candidate list was precisely the part that did not
> constitute a discovery.

Three limits bound that statement: the reference lists contain four and six
genes, so intervals are wide; the full-power list is not ground truth but a
higher-powered list that is itself unstable, so these curves measure agreement
between two imperfect lists rather than recovery of truth; and the novel-list
composition is itself questionable, since KANSL1 lies in the chr17q21.31
inversion and KIAA0040 was later downgraded.

## 2.5 Finding ⑥: MR names the gene that has an instrument

Across 22 glycolytic enzymes measurable in CD4⁺ T cells, 20 are higher in
non-responders post-treatment (8 surviving FDR individually) and 18 of 22
pre-treatment (4 surviving FDR), with **PGAM1** strongest at baseline and TPI1
fourth. These counts are descriptive: because the enzymes are correlated, a
binomial test against one half is anticonservative, and the
correlation-preserving permutation that replaces it returns P = 0.088 and 0.118
respectively (§3.3). Genetically, of 28 glycolytic genes **3 are instrumentable in
this exposure resource, 2 remain analysable after harmonisation with this outcome,
and 1 yields a nominal association** — TPI1 (OR 1.31, P = 1.4×10⁻³), with ENO1
null (P = 0.29) and SLC2A1 lost at harmonisation. PGAM1, GAPDH, PKM, LDHA and PFKP
carry no instrument at all. An instrument-free pathway enrichment test was also
null (mean |z| 0.890 versus 0.937 ± 0.267, P = 0.52), though with nine
gene × profile combinations it is underpowered and we do not claim the pathway is
signal-free.

**Functionally the phenotype is pathway-wide; genetically the framework can see
three of its genes, test two, and name one.** The strongest baseline enzyme is
PGAM1, which carries no instrument; the gene MR names is TPI1. Left undeclared, this invites the reading that the named gene is the
pathway's most important member. A related constraint operates at the resource
level: TPI1 is absent from all eight stimulated T-cell datasets we identified in
the eQTL Catalogue despite TPM 98–194 where quantified, so gene-level filtering
in public resources determines which findings can be tested for replication.

**This finding is the bridge to Part II.** It says the framework's output is
determined by instrument availability rather than biological importance — which
is a limitation for nomination, and simultaneously an invitation to ask what the
instrument's *location in time* tells us.

## 2.6 Finding ⑦: the functional layer has its own failure modes

**A score-based split also splits cell-type purity.** Dividing CD4⁺ T cells by
glycolysis score produced an apparent phenotype — reduced cytotoxic, increased
helper/regulatory signature (Δ −0.318 and +0.235, 24/28 samples, FDR 1.1×10⁻⁴).
The modules contain lineage markers, and the split was confounded: CD4 detected
in 40.1% of glycolysis-high versus 19.5% of glycolysis-low cells, CD8A in 47.6%
versus 63.9%. Removing lineage genes and matching on a CD8ness score in addition
to depth (610 cells, all 28 samples, residual SMD 0.0097 and 0.0000) abolished it
entirely: the lineage-free effector score gave Δ = +0.013 (P = 0.88) with 14/28
concordance — exactly chance — and FOXP3, CTLA4 and GZMB returned P = 0.73, 0.80,
1.00. Effect sizes collapsed 75–100%, so this is not power loss. A cluster-level
control had passed; contaminating cells are distributed within clusters, so the
control must be at the cell level.

**Two processing errors changed every patient-level number.** An object had been
re-normalised (log₂(TPM+1) passed through `NormalizeData()`; Spearman 0.874
against correctly scaled values), and the response table had merged samples with
conflicting outcomes (patients with two post-treatment biopsies can carry
different labels). Both were corrected and all patient-level analyses recomputed
(final: 9 R / 10 NR pre-treatment, 8 R / 20 NR post-treatment). The corrected
numbers strengthened the main claims but eliminated one mechanism we had proposed:
HLA-C's strong negative correlation with proliferation (ρ = −0.612/−0.720)
disappeared (+0.142, P = 0.56; −0.120, P = 0.54). That mechanism had already been
withdrawn when it failed to replicate — replication caught the artefact before we
knew its cause. The same test applied to the second cohort would have been wrong:
its raw file is *not* normalised (CV = 1.15, 46-fold range). The criterion is
whether the source file is already normalised, not whether a derived object
matches it.

## 2.7 Locating an effect when MR cannot (finding ⑤), and TPI1 as a bounded nomination

**Gene attribution.** SPSB2 and TPI1 lie 3–20 kb apart on chr12p13, with the TPI1
instrument bracketed by two SPSB2 instruments; the six candidate genes correspond
to five independent loci. Three orthogonal lines point to TPI1 as the more plausible effector gene, though genetic attribution at this locus remains formally unresolved:
single-cell response stratification is significant for TPI1 (P = 0.0076 and
0.0057) and null for SPSB2 (0.17); TPI1 is detected in 60.6% of CD4⁺ T cells
versus 5.0% for SPSB2 in an independent cohort; and the two genes have entirely
different eQTL dynamics (§3.1). The lines come from different data types, and the
approach generalises.

**Compartment attribution — where the answer is negative.** MR nominates a gene
using instruments measured in one cell type, but validation is almost always
attempted in tissue-level data, where that cell type may contribute a minority of
the signal. Whether it does is directly measurable, and for TPI1 it is decisive.

In melanoma single-cell data with all cell types annotated, TPI1 is far higher in
malignant cells than in CD4⁺ T cells: +2.75 log₂ units (≈6.7-fold) in 16 of 16
patients (P = 3×10⁻⁵), detected in 95.3% of malignant cells versus 60.6% of
CD4⁺ T cells. The difference survives cell-level matching on sequencing depth
(+2.03, 10 of 11 patients; residual standardised mean difference −0.007) and
replicates in an independent cohort (+1.28, **11 of 11** patients, P = 1×10⁻³),
with pre-specified positive controls passing in both (MLANA higher in malignant
cells, PTPRC higher in T cells). The same holds for the locked glycolytic
signature, with and without TPI1 (Fig 5).

The arithmetic consequence is what matters. With a 6.7-fold per-cell ratio and the
cell-type proportions typical of melanoma tissue, CD4⁺ T cells contribute on the
order of one to two per cent of the tissue-level TPI1 signal. **Bulk-tissue TPI1
abundance is dominated by the malignant compartment and therefore cannot validate
a CD4-specific mechanism**, whatever its P value.

The precise form of that statement matters, because the looser version — that a
tissue-level measurement simply *is* a measurement of tumour glycolysis — claims
more than the arithmetic supports. Bulk TPI1 can still covary with clinical
outcome through immune infiltration, tumour purity, or a metabolic state shared
across compartments, and any of those would produce a real association that has
nothing to do with a CD4-specific effect. What the compartment ratio establishes
is that such an association **cannot be attributed** to CD4⁺ T cells, not that
the association is spurious.

Spatial transcriptomics gives the same answer with spatial context (four patients,
eight sections, 2,317 spots): tissue-level TPI1 tracks the glycolytic module
(ρ = +0.174, 8/8 sections) and the tumour compartment (+0.137), and correlates
*negatively* with the lymphoid compartment (−0.080); tumour spots exceed lymphoid
by +0.182 (P = 2.5×10⁻⁷); after regressing out the glycolytic module the tumour
association vanishes (ρ = +0.026) and no hidden lymphoid signal emerges. Its
interpretability rests on a mirror-image positive control in the same sections:
HLA-C tracks lymphoid (+0.187, 8/8) and myeloid (+0.231) and not tumour
(−0.043, 1/8); SMC2 tracks proliferation (+0.105, 8/8). We also record a test that
failed its own control and is therefore uninformative: within the 450
lymphoid-dominant spots neither TPI1 (P = 0.13) nor HLA-C (P = 0.23) behaved as
expected.

The consequence is that the two bulk cohorts in which TPI1 tracked response
(GSE78220 P = 0.045; GSE91061 on-treatment P = 0.030, pre-treatment P = 0.81)
**cannot** be reported as replication of a CD4-specific effect. This is also the
reason Part II uses only CD4-resolved patient data.

The consequence for survival analysis is worth demonstrating, because it is the
form in which such claims are most often made. In TCGA-SKCM (457 patients, 215
events), bulk TPI1 predicts worse overall survival (HR = 1.23, P = 1.6×10⁻³), as
does the glycolytic signature with TPI1 removed (HR = 1.27, P = 0.048). Read
naively, this looks like support. It is not: bulk TPI1 correlates negatively with
immune infiltration (ρ = −0.21, P = 4×10⁻⁶) and positively with tumour content
(ρ = +0.13), immune infiltration predicts better survival (HR = 0.70,
P = 6×10⁻⁷), and adjusting for both removes the TPI1 association
(HR = 1.13, P = 0.06). The bulk association is a statement about tissue
composition.

**This is not in tension with our CD4-resolved result, and the distinction is the
point.** Tissue composition and within-cell-type expression are different
quantities: tumours richer in malignant cells carry more TPI1 signal and fewer
lymphocytes, while separately, among CD4⁺ T cells, non-responders express more
TPI1 than responders. Conflating the two is exactly the error this diagnostic is
meant to catch.

We note the generalisable form, because it costs little to apply: before citing
tissue-level data as validation of a cell-type-specific nomination, compute the
nominated gene's expression ratio across cell types in any annotated single-cell
atlas of that tissue. If the target cell type is a minority contributor, the
tissue-level result is uninformative about the claim regardless of its
significance. This requires no spatial data and takes minutes.

**TPI1 as a nomination, and its exact boundaries.** TPI1 is the only candidate
passing MR, colocalisation, SMR/HEIDI and multi-instrument sensitivity analysis:
OR = 1.305, P = 1.4×10⁻³ but **FDR = 0.119**; colocalisation **PP.H4 = 0.51**, at
the criterion edge; SMR P = 3.0×10⁻³ with HEIDI P = 0.45; over five independent
instruments IVW P = 3.8×10⁻⁴, weighted median P = 2.6×10⁻³, Q P = 0.83. It
entered through the conjunction, and the conjunction carries it because each
component is individually marginal. To that we add the sensitivity result of
§2.4: **SuSiE recovers no credible set in this region at either power**, so the
posterior of 0.51 is computed where the outcome GWAS carries no resolvable
signal. That posterior is also **prior-dependent**: across the window-and-prior
grid of §2.4, TPI1 meets our colocalisation criterion in 8 of 16 combinations,
passing at the default prior and looser (PP.H4 = 0.49 and 0.83) but failing at the
most conservative one (PP.H4 = 0.09). Whether this locus colocalises is therefore
a statement about a prior as much as about the data. Across five other cancers we found **no clear evidence of association**
(P = 0.086–0.69), with point estimates reversed for breast (0.915) and prostate
(0.878), and likewise **no clear evidence of an effect on naevus count**
(OR 0.867, P = 0.079) in a test whose sensitivity is established by five
pigmentation loci that do show one. Neither is evidence that TPI1 acts only on
melanoma: both are non-significant results in modest samples, and we count them as
inconclusive rather than supporting (Box 1).

**What cannot be claimed.** Stated precisely, TPI1 is **a nominal,
prior-sensitive candidate selected by a multi-layer conjunction**: MR
P = 1.4×10⁻³ but FDR = 0.119, PP.H4 = 0.51 passing in 8 of 16 window-and-prior
combinations, in a region the outcome GWAS cannot fine-map, with gene attribution
at chr12p13 unresolved by genetics. It is DepMap-essential with pLI = 0.87, so no
target claim is made, and by finding ④a it may not survive a higher-powered
outcome GWAS. Calling it "the gene MR named" is accurate; calling it a discovery
is not. The genetic side of this nomination is weaker than four passed filters
make it appear, and reporting it as such is the point of the worked example rather
than a caveat attached to it. **Part II does not depend on any of it.**

**A perturbation experiment would not settle it either, and the reason is
measurable.** Whether the nominated gene has immunological consequences is the
question this design cannot answer, and the obvious next experiment — knock it out
and look — is bounded by the same property that makes it DepMap-essential. Across
the **1,471 human CRISPR screens** in BioGRID ORCS in which TPI1 was measured it
scores as a hit in **628 (42.7%)**. That is the profile of a core-essential gene:
GAPDH is 46.6% and PGAM1 47.2% in the same library, while lineage-defining or
trait-specific genes measured in the same screens sit an order of magnitude lower
(IRF4 3.9%, FOXP3 1.5%, MC1R 1.1%, ZFYVE19 1.2%). A knockout of this gene
therefore produces a fitness phenotype in almost any cell type, so **a knockout
cannot isolate a CD4-specific role** — neither a real one nor an in silico
substitute — and the intervention is in any case not comparable to the small,
graded expression shift an eQTL represents. This is why we report no perturbation
experiment rather than a simulated one. It is a boundary on what perturbation
could add here, **not evidence for the nomination**, and it should not be read as
CRISPR data supporting our conclusion. The observation that carries further is
about the pathway rather than the gene: **PGAM1 — the strongest baseline enzyme in
the functional data (§2.5), and the one MR cannot see at all because it carries no
instrument — has the same pan-essential profile.** Essentiality is a property of
glycolysis here, not of the gene MR happened to name, which makes it a second
instance of finding ⑥ rather than an independent observation.

**No candidate is strong on all axes (Fig 8).** Statistical strength, pathway attribution and druggability point at different
genes. ZFYVE19 has the strongest colocalisation (PP.H4 = 0.99, peaks coincident,
zero catalogued pleiotropy, non-essential, LoF-tolerant) but affects naevi
(P = 5.3×10⁻³), shows no functional signal, and is pan-cancer in direction. SMC2
has the largest effect and cleanest pathway attribution but is DepMap-essential
(pLI = 0.99999) and significant in four cancers, identifying a general
proliferative mechanism rather than immune surveillance. HLA-C carries the most
instruments and a vitiligo mirror (OR 0.681, 391 cases) but lies in the most
pleiotropic region of the genome with SJS/TEN safety flags. KIAA0040 gave opposite
directions in two cohorts and was removed. We report the matrix rather than
promoting a target, because the alternative requires exactly the selective
emphasis these diagnostics detect.

---

# PART II — What the same data still establish once nomination is set aside

Part I is longer than Part II because it documents failures, each of which
requires its own control to be interpretable. That asymmetry is not a
statement of relative importance: Part II carries the only claims in this
paper that we expect to survive a larger outcome GWAS, precisely because
they do not depend on one.

> Part II uses only quantities that outcome power cannot change: where in the
> activation time course an instrument exists (exposure side), what programme
> those genes belong to, and how that programme behaves in CD4-resolved patient
> data. No claim in Part II routes through the melanoma MR.
>
> The three datasets use a **single locked gene signature of 16 glycolytic
> genes** — the intersection of the sets used on the eQTL (28), patient (22) and
> multiome (16) sides, which are nested — and every signature-level result is
> repeated with TPI1 removed.

## 3.1 Instrument availability is the binding constraint, and it is severe

Finding ⑥ showed that MR names the gene that carries an instrument. Applied to a
whole pathway, the constraint is more severe than that phrasing conveys, and we
report the full matrix rather than the two genes that happen to be testable
(Supplementary Table S13: 28 glycolytic genes × 8 activation profiles, the
strongest cis-eQTL P value in each cell).

**Three genes of 28 reach P < 5×10⁻⁸ in any profile**: ENO1 (naive 0 h,
P = 5.3×10⁻⁹), SLC2A1 (naive 40 h, P = 3.8×10⁻⁸) and TPI1 (naive 16 h,
P = 1.7×10⁻¹²; memory 16 h, P = 1.2×10⁻⁹). Of these, SLC2A1's variant did not
survive harmonisation with the outcome, so **two are testable by MR**, and one of
those two — ENO1 — is null (OR = 1.02, P = 0.29). The functional data, by
contrast, show a coordinated shift across the whole pathway.

**These three numbers are not the same kind of statement, and we keep them
apart.** *3 of 28 instrumentable* is a property of the pathway and this exposure
resource: it would be unchanged whatever outcome were analysed, and it is the only
one of the three that belongs to Part II's exposure-side thread. *2 of 28
analysable* adds the requirement that the variant harmonise with a particular
outcome GWAS — SLC2A1 fails on this outcome and might not on another, so this is
not a statement that SLC2A1 lacks an instrument. *1 of 28 nominally associated*
adds the outcome's own significance and therefore belongs to Part I, carrying all
the instability Part I documents. Compressing the three into "one usable
instrument" would attribute two outcome-side limitations to the biology of the
pathway. **The exposure-side constraint alone — 3 of 28 — is already the single
most restrictive step in the framework.**

**We do not claim that these instruments are confined to activated states.** An
earlier version of this analysis said so, and the full matrix shows it is wrong:
the three instruments sit at 0 h, 16 h and 40 h respectively, and ENO1's is in
the resting state. Two distinct quantities were being conflated. The timepoint at
which each gene's *strongest* eQTL falls, regardless of significance, does shift
away from rest for this gene class (2 of 24 gene × lineage combinations at 0 h
versus 22.9% of background genes; OR = 0.31), but that comparison does not reach
nominal significance on the locked signature (P = 0.064) and is confounded by
differential measurability — the locked genes are estimated more precisely than
background at every activated timepoint (SE ratio 0.57–0.77) and barely so at rest
(0.86, 0.89), which mechanically moves each gene's strongest eQTL later
(Fig 9a–b). The timepoint at which a *genome-wide significant instrument* exists —
the quantity that actually determines what can be studied — shows no such shift in
these data.

**What survives is the single-gene observation.** TPI1's instrument exists at 16 h
and at no other timepoint, in both lineages: naive P = 8.4×10⁻⁵ (0 h) →
**1.7×10⁻¹²** (16 h, 148 variants below 5×10⁻⁸) → 8.7×10⁻³ (40 h) → 3.5×10⁻³
(5 d); memory 6.2×10⁻³ → **1.2×10⁻⁹** (53 variants) → 3.7×10⁻² → 7.1×10⁻³.
Instrument yield is flat across profiles (51.3–55.0%) and slightly lower at 16 h
than 40 h, so this is not differential detection at the profile level; and the
co-located SPSB2 is constitutive (1.4×10⁻¹¹ at 0 h). Examining effect estimates
rather than P values shows what the window is and is not (Fig 9a): the 16 h
P value is smallest because the standard error is roughly five-fold smaller, and
in naive cells the point estimate is largest **at rest**. Comparing estimates
directly, the 16 h effect exceeds 40 h and 5 d in both lineages (four of four
comparisons, P = 6.0×10⁻⁴ to 0.043) but does not exceed 0 h. The admissible
statement is about measurement, not regulation:

> For this gene, a usable instrument exists only at one timepoint. A study using a
> resting-state resource would find none, and could not distinguish that absence
> from a true null.

Distinguishing activation-dependent *regulation* from activation-dependent
*measurability* requires a genotype × pseudotime interaction test, and that test
has already been run and published for this resource. Its Supplementary Table 8
reports linear and quadratic interaction models for every lead eQTL variant, and
**all three TPI1 lead variants tested return no evidence of interaction**: in
memory cells, `12_6908821_T_C` gives linear P = 0.974 and quadratic P = 0.992; in
naive cells, `12_6868195_A_G` — 1,063 bp from our instrument and in high LD with
it — gives 0.707 and 0.891, and `12_6908616_AAAAACC_A` gives 0.378 and 0.605.
Meanwhile the *expression* of TPI1 is among the most strongly pseudotime-dependent
in the dataset (Moran's I = 0.664, q ≈ 0; Supplementary Table 6).

Three things are therefore separable here, and the published interaction test
tells them apart:

| | TPI1 |
|---|---|
| **Dynamic expression** — does the gene's level change with activation? | **Yes**, strongly (Moran's I = 0.664) |
| **Time-specific eQTL significance** — does the eQTL P value change with time? | **Yes** — genome-wide significant at 16 h only |
| **Dynamic genetic effect** — does the genotype's effect change with time? | **No evidence**, in any of three lead-variant tests |

**We therefore state the conclusion positively rather than as a limitation: the
apparent 16-hour specificity of this instrument reflects time-dependent
detectability, not a demonstrated dynamic genetic effect.** An earlier version of
this manuscript said the interaction test required individual-level data we did
not have; that was wrong — the test is public, we have now used it, and it
supports the measurement reading over the regulation reading.

**Is this gene quantified reliably at all?** The question is forced on us by an
unusual pattern: TPI1 is missing from all eight stimulated T-cell datasets we
identified in the eQTL Catalogue, from 11 of 15 DICE datasets, and from the Visium
probe panels of two independent spatial studies. The common thread is plausibly
its four processed pseudogenes, which make unique quantification — and unique
probe design — difficult. Three arguments say our quantification is nonetheless
sound. First, the cis architecture cannot be manufactured by mismapping: the lead
variant lies 3,816 bp upstream of the TSS, 148 variants exceed 5×10⁻⁸ across
102 kb, and −log₁₀P decays from 11.8 within ±2 kb to 1.0 beyond 200 kb, whereas
TPI1's pseudogenes lie on chromosomes 1, 4, 6 and 7 and would generate trans
associations. Second, for SPSB2 — quantified both by us and by the standardised
processing of the same experiment — our estimates match theirs (four variants,
β within 10%, direction 4/4). Third, TPI1 is well measured wherever it is
measured: detected in 60.6% of CD4⁺ T cells by Smart-seq2 and quantified at
98–194 TPM in the DICE datasets that report it. The absences reflect pipeline and
panel design, not a universal failure of measurability — but they do mean this
finding is unusually hard for others to check, which we regard as part of the
result.

## 3.2 The state: a metabolic programme partially distinct from activation

If instruments point anywhere, they should point to a coherent cell state rather
than a single transcript. In purified CD4⁺ T cells profiled by 10x multiome after
anti-CD3/CD28 stimulation (40,495 nuclei), the locked glycolytic signature is
close to orthogonal to activation intensity at the RNA level (R² = 0.024 for
signature ~ activation + depth, against a pre-specified R² > 0.70 cut-off for
declaring confounding) and is not a nuclear/cytoplasmic composition artefact
(nuclear-retention index removes 0.047 of variance). This answers the obvious
objection that the programme is simply activation — at the RNA level.

Ranking all measurable genes along the axis without pre-specified modules places
ribosomal proteins (12.0% of the top 300, 10.3-fold) and OXPHOS at one end and
T cell identity and quiescence factors (TOX, IKZF2, MYB, TXK, CAMK4, SESN3) at the
other: an anabolic-growth versus T-cell-identity axis.

At the chromatin level the axis shares 42% of its variance with activation
(r = 0.647), so the programme is **partially, not fully, distinct** from
activation and we describe it that way. Restricting to activation-invariant peaks
(retaining 68% of axis amplitude), the high end retains all 23 AP-1 motifs
(median odds 1.81 → 1.52; BATF::JUN 2.22 → 1.62, P = 9.1×10⁻¹⁰) while the low end
is enriched for IRF-family and STAT1::STAT2 interferon-response elements. Motif
similarity within families means these identify families rather than members, and
the restriction is mitigation rather than separation.

An earlier round of pre-specified immune modules was withdrawn on detection
grounds rather than interpreted: Th2 (median detection 2.4%), Th17 (4.4%) and
FOXP3 (4.1%) are not measurable in this system.

**This is convergence across data types, not a mechanistic chain.** Expression and
chromatin move together along one axis. Nothing here shows that the TPI1 eQTL
causes the state, and three attempts to make that link failed: mean accessibility
at the TPI1 promoter is unchanged by activation (log₂FC = +0.047), the peak
containing the instrument does not vary with the axis (mean log₂FC = −0.013,
consistent in 2 of 4 replicates), and the instrument disrupts no AP-1 motif
(maximum |Δscore| = 0, empirical P = 1.0 against 500 regional variants). The
matched test — allele-specific accessibility — needs read-level data with
genotypes that are not available.

## 3.3 In patients: present in every cohort, confirmed in none twice

In CD4⁺ T cells from melanoma patients treated with checkpoint blockade, the
locked signature is higher in non-responders. We report this with the statistic
that respects the data structure, having found that the one we used first does
not, and with an independent cohort, having found that only part of it replicates.

**The statistic.** Glycolytic genes are strongly correlated across samples (median
pairwise Spearman 0.48; effective number of independent genes ≈ 11.7 of 16), so a
binomial test on the count of genes higher in non-responders is anticonservative.
Permuting response labels — which moves the whole signature together and therefore
preserves the correlation structure exactly — gives P = 0.088 (post-treatment) and
0.118 (pre-treatment) for that count, against binomial values three orders of
magnitude smaller. **We therefore report gene-direction counts descriptively and
without P values.**

**The primary estimand and the primary test are fixed once, here, and are the same
in all three cohorts.** The estimand is the difference in the composite signature
score — the mean of the 16 z-scored genes, computed per sample — between
non-responders and responders within one treatment arm. The test is a permutation
of the response label, which moves the whole signature together and so preserves
the correlation structure exactly; it was pre-specified as primary in both
follow-up pre-registrations (Supplementary S18, S21) and is the statistic on which
every significance statement below rests. The Wilcoxon rank test is reported
alongside as a secondary, distribution-free check, and where the two disagree we
say so and treat the disagreement as information about the size of the cohort
rather than as licence to choose. We do not switch statistics between cohorts or
between arms.

**Discovery cohort.** Four patients contribute two post-treatment samples each,
two of them with opposite response labels, so all three ways of handling the
repetition are reported (Supplementary Table S19). Post-treatment: Δ = +0.87 [0.37, 1.31],
permutation P = 0.0005, holding under every subsetting and with TPI1 removed.
Pre-treatment: Δ = +0.74 [0.23, 1.30], permutation **P = 0.009 — significant on
the primary test** — with the secondary Wilcoxon at P = 0.065, on 9 responders and
10 non-responders.

**Independent cohort.** We searched systematically for a replication cohort
(Supplementary S14: five queries across GEO and ArrayExpress, twelve datasets
assessed, inclusion criteria fixed in advance). Exactly one public dataset met all
of them — melanoma, tumour tissue, CD4-resolved, response-annotated, with
pre-treatment samples and without confounding between response and regimen — and
we pre-registered the analysis before reading it (Supplementary S18, including the
reading table against which the outcome below was scored). The results are mixed
and we report both halves.

*Pre-treatment* replicates in magnitude but not in significance. Using the
authors' own cell-type annotation, Δ = +0.56 (CD4 T cells) to +0.65 (including
regulatory T cells), against +0.74 in the discovery cohort, with the predicted
direction in all twelve configurations tested. The effect is not explained by
anatomical site: sites are balanced across response groups, the signature does not
differ by site (P = 0.40), centring within site leaves Δ = +0.63, and restricting
to lymph node — the majority site — gives Δ = +0.89. But with 8 responders and
**3 non-responders**, permutation P = 0.09–0.11. **The effect sizes agree across
the two cohorts; the replication cohort does not reach significance on the primary
test, while the discovery cohort did.** With three non-responders, no configuration
of this cohort could have reached significance at the observed effect size.

*Post-treatment does not replicate.* This was the stronger arm in the discovery
cohort (P = 0.0005). In the independent cohort it is null and its sign depends on
the cell population: Δ = −0.15 to −0.03 for annotated CD4 T cells, +0.14 to +0.30
when regulatory T cells are included. A result whose direction changes with a
population-definition choice is not a result, and we do not present it as one.

**A second tumour type.** The HCC cohort assembled for §2.2 turned out to carry
per-patient response labels — in the GEO sample records rather than in the
deposited object — and to meet every criterion in our replication search except
that the tumour type differs, which was not one of the criteria. It is
CD45⁺-sorted (so it cannot address compartment attribution) but resolves 59,528
CD4⁺ T cells across four responders and five non-responders, with paired
pre- and post-treatment tumour and blood. We pre-registered this analysis too
(Supplementary S21), including the fact — computed from the sample structure
before any expression value was read — that the pre-treatment tumour arm has four
responders against **two** non-responders and therefore *cannot* reach P < 0.05
whatever the effect size, exactly the criterion by which we excluded another
dataset from the replication search.

Post-treatment tumour, the pre-specified primary test, is concordant: Δ = +0.874
against +0.87 in the discovery cohort, one-sided permutation P = 0.043, with all
four positive controls passing. It is stable across the choices that destabilised
the melanoma replication — direction and magnitude are unchanged whether
regulatory T cells are included and whether TPI1 is removed — and a post-hoc
leave-one-out shows Δ between +0.71 and +1.06 with every sample in turn omitted.
Two things must be said with it: the rank test on the same data is not significant
(Wilcoxon P = 0.114, one non-responder falling inside the responder range), and
with eight samples the significance sits at the edge of what the design can
express, seven of eight leave-one-out runs giving P between 0.057 and 0.086.
Pre-treatment tumour shows essentially no effect here (Δ = +0.08) in the arm that
could not have reached significance in any case.

**What we conclude.** On the primary test, both arms were significant in the
discovery cohort (P = 0.0005 post-treatment, P = 0.009 pre-treatment); in the two
follow-ups **no arm reached significance twice**, and the follow-ups fail in mirror
image — pre-treatment matches in magnitude in the same-disease cohort
(Δ = +0.56 versus +0.74) but is underpowered there and absent in HCC, while
post-treatment reverses in the same-disease cohort and is reproduced in magnitude
and significance in HCC (Δ = +0.874, P = 0.043). We therefore describe the
patient-level finding as **suggestive: significant on discovery in both arms,
with one arm reproduced in a second tumour type and neither arm confirmed in a
same-disease cohort** — not as an established stratifier, and not as replicated: the concordant cohort is a different cancer under a different regimen
(anti-PD-1 plus a tyrosine kinase inhibitor), where a negative result would have
been uninterpretable and a positive one therefore cannot carry the weight of a
same-disease replication. We had noted, when only the melanoma replication
existed, that the arm which looked strongest on discovery was the one that failed
— the behaviour of a winner's-curse-inflated estimate. The HCC result does not fit
that reading, since the same arm is reproduced there at full magnitude. We report
the discrepancy rather than a story that accommodates it.

The discovery result is otherwise heavily audited: it survived correction of a
re-normalisation error that changed every patient-level number; redefinition of
samples and response labels from the primary annotation; residualisation on a
lineage score; restriction to CD8-negative cells; restriction to the locked
signature; and removal of TPI1. That auditing is why we can state precisely what
did and did not replicate, rather than only that something did — but it is
internal auditing, and the three cohorts show that surviving it is not the same as
surviving transfer to other data.

Two boundaries belong here. Tissue-level measurements of these genes reflect
tumour rather than CD4 biology (§2.7), which is why only CD4-resolved data are
used. And within the in vitro system the high end of the axis carries *low* TOX,
whereas in patient tumour-infiltrating CD4 cells TPI1 correlates positively with
an exhaustion score (ρ = 0.617–0.726 in two cohorts); early-activation TOX and
exhaustion-associated TOX are not necessarily comparable, but the tension stands
and we do not resolve it.

## 3.4 What Part II licenses, and what it does not

**Licensed.** The glycolytic programme in CD4⁺ T cells is largely inaccessible to
this framework, and the inaccessibility has three separable layers: **3 of 28
genes are instrumentable in this exposure resource, 2 of 28 are analysable against
this outcome, and 1 of 28 yields a nominal association.** Only the first is
licensed by Part II's exposure-side data; the second and third depend on the
outcome and inherit Part I's instability.
The programme corresponds to a definable cell state with coordinated
transcriptional and chromatin structure, partially distinct from activation
intensity. In patients, that state's activity in CD4⁺ T cells is higher in
non-responders in the discovery cohort at both timepoints on the primary test
(P = 0.0005 post-treatment, P = 0.009 pre-treatment); in two pre-registered
follow-ups neither arm was confirmed a second time, the pre-treatment effect
matching in magnitude but not significance in the same-disease cohort and the
post-treatment effect reproducing in the cross-disease one. For TPI1 specifically,
an instrument exists at one timepoint only.

**Not licensed.** We do not claim that instruments for this pathway are confined
to activated states — the full matrix refutes it. We do not claim that genetic
regulation is activation-dependent; measurability is, and the two cannot be
separated by us; the published genotype × pseudotime interaction test settles it
in favour of measurability (§2.7). We do not join the Part I result to the response association: the two are
different estimands, measured in different data, and their directional concordance
is worth stating but is not a causal chain. We
make no target claim: TPI1 is DepMap-essential, and it is the instrumentable
member of a programme rather than its dominant regulator — PGAM1 is stronger at
baseline and carries no instrument. And we claim no independent replication:
the 16 h observation cannot be tested elsewhere because TPI1 is not quantified in
the available resources, and the patient result was carried to the one eligible same-disease cohort and to one
cross-disease cohort, and no single arm was confirmed in both.

> Within this pathway and this resource, dynamic eQTL data are more informative
> about **what can be interrogated genetically, and when** than about **which gene
> should be targeted**.

---

## 4. Discussion

---

> ### Box 1 | The process, in summary
>
> The full record — every attempt, the criterion that decided it, its
> pre-specified status and its date — is Supplementary S12. This box gives only
> the classification, because the classification is what the argument uses.
>
> **Eleven target-level tests were run. They did not divide into successes and
> failures, and we no longer present them that way.**
>
> | Class | n | What it means | Examples |
> |---|---|---|---|
> | **Overturned** | 7 | A positive claim was tested and did not survive, or was explained away by a control | PADI4 under multi-instrument analysis; the entire FinnGen-round candidate list under a higher-powered outcome; HLA-C acting through proliferation; bulk cohorts as CD4-specific evidence; a distinct transcriptional state for glycolysis-high cells; pathway-level MR; cross-cancer specificity |
> | **Inconclusive** | 3 | The test returned no verdict — an underpowered null, which is absence of evidence and not evidence of absence | pathway enrichment of GWAS signal (P = 0.52, n = 9); TPI1 in five other cancers (P = 0.086–0.69); TPI1 and naevus count (P = 0.079) |
> | **Supported, then not confirmed on transfer** | 1 | Significant on discovery, not confirmed in the same arm in either follow-up | TPI1's ICB response stratification (§3.3) |
> | **Not testable** | 3 | The data cannot answer it; reported so the denominator is complete | chromatin test mismatched to the claim; TPI1 absent from all eight stimulated T-cell resources; a precondition test that stopped an analysis (ICC 0.0285 against a matched floor of 0.0403) |
>
> ⚠ **We previously reported this as "eight of eleven overturned, three
> supported".** That tally was wrong in two directions. Two of the three
> "supporting" results were non-significant tests — no association with naevi, no
> association in five other cancers — and a non-significant test in a small sample
> is absence of evidence, not support. The third, the response stratification, was
> significant on discovery but has since not been confirmed in the same arm in
> either follow-up. **The corrected classification is above, and the "8 of 11"
> figure appears nowhere in this paper's conclusions.**
>
> **Five genes were designated lead candidate; the first four designations were
> overturned** — PADI4, ZFYVE19, the ZFYVE19+SMC2 module, SMC2, then TPI1. The
> criterion by which candidates were ranked changed at each step and no ranking
> rule was written down in advance. This is the weakest point in our procedure
> (§4.2).
>
> **Claims withdrawn during revision:** a binomial test on correlated glycolytic
> genes; the statement that instruments for this pathway exist only in activated
> states; Cochran's Q at 5.33% as evidence of calibration; our down-sampling
> calibration described as external rather than empirical; and a post-hoc
> winner's-curse reading of the failed replication, which the second-tumour result
> contradicts. Two processing errors of ours are reported as results rather than in
> Methods (§2.6), because they are instances of finding ⑦.
>
> **Seven of our own methodological claims were pre-registered with failure
> conditions written before the data were read**: the power-trajectory predictions
> (12 of 12 intervals hit); the same-disease replication of the patient
> stratification (partial); the second-tumour generalisation (directionally
> consistent, primary test not significant); the patient stratification in that
> second tumour (primary test passed); the exposure-resource generalisation
> (primary test passed, P = 3.7×10⁻¹¹); the disease × resource grid (all six
> cells concordant in direction, negative control clean); and the transfer test to
> the current release of the outcome resource (all six point predictions inside
> their registered intervals, though the point value of the first was not
> attained). An eighth — extending compartment attribution to a second tumour —
> could not be run at all. **Four passes, two partial results, one reduced in
> scope, one foreclosed.**

---

### 4.1 The seam: why Part II does not inherit Part I's *outcome-side* instability

A paper that spends its first half showing candidate lists are unstable owes an
account of why its second half is exempt. The account is structural, not
special pleading.

With one instrument per exposure, z = β_out/se_out. Given the instrument set,
every quantity in Part I — which genes reach significance, whether colocalisation
calls H3 or H4, what the list contains — is a function of the outcome GWAS, and
we show it changes when that GWAS changes. The size of the instrument set is not:
that is set by the exposure resource, which is why swapping the exposure
multiplies the number of tests without altering where the significant ones land. Every quantity in Part II is either a property of
the exposure data (which timepoint carries the strongest eQTL; whether the
strongest eQTLs of a gene class avoid rest) or of cell-resolved patient data
(whether a state separates response groups). **No amount of outcome-side
instability can move them**, because they never consult the outcome.

This also explains the one place the seam is visibly load-bearing: TPI1 appears
in both parts. In Part I it is a nomination whose genetic evidence is marginal
and, as §2.4 shows, sits in a region without a resolvable outcome signal. In Part
II it is an instance of gene-level instrument availability being state-bound, established by an
exposure-side observation that the Part I result cannot touch. The same gene,
two evidence grades, no borrowing between them.

### 4.2 The central argument of Part I: the failure rate is the finding

A pipeline that passes every positive control and recovers known melanoma biology
at five nested power levels, in two diseases, nevertheless sustained no
target-level claim: of eleven tests, seven overturned the claim under test, three
were inconclusive, and the one that was significant on discovery has not been
confirmed in the same arm in either follow-up (Box 1).

We are explicit about what that classification does and does not license. An
inconclusive test is **not** evidence for the claim it failed to reject: our
earlier framing counted two underpowered nulls — no association with naevi, none
in five other cancers — as support, and that was an absence-of-evidence error,
now corrected. Nor is the count itself a rate to be quoted: eleven tests chosen by
us, on one analysis, is a description of this project's history and not an
estimate of how often the framework fails. What it does license is narrower and,
we think, more useful: **each of the seven overturned claims was overturned by a
check that a study of this kind could have run and typically does not** — which
is why §4.5 states them as procedure rather than as a verdict on anyone's results.
Studies using the same framework that do not report such attempts have not
necessarily avoided these failures; they have not tested for them.

Three of the outcomes in that record are *not* scientific negatives and are not
counted with the eight: a chromatin test mismatched to the claim, since an
activation-specific eQTL is a genotype × timepoint interaction and mean
accessibility is neither necessary nor sufficient for it; unavailable replication,
TPI1 being absent from all eight stimulated T-cell datasets and the one
design-matched resource being the same experiment; and a **precondition test that
stopped an analysis** — before asking whether peripheral-blood CD4 glycolysis
could serve as a correlate, we tested whether it is an individual-level trait at
all, found it indistinguishable from random gene modules within a clinically
homogeneous group (ICC 0.0285 against a matched floor with 95th percentile
0.0403), and cancelled what followed. That is a stopping rule working, at a cost
of about fifteen minutes, on an analysis whose positive and negative results would
both have been uninterpretable.

**How much was chosen, and on what criteria.** A failure tally answers "how many
times did you try to substantiate a target". It does not answer "out of how many
possible targets, pathways and stories did you pick this one". We give that
denominator at four levels (Supplementary S12 §5), because two are clean and two
are not.

*Genes*: across two outcome GWAS the pipeline nominated ten novel-locus genes; six
survived four-way verification and entered functional evaluation. TPI1 was the
fifth designated lead candidate and the first four designations were each
overturned, under a ranking criterion that changed at every step and was never
written down in advance (Box 1). Worse, the single criterion that selected TPI1 —
significance in ICB response stratification at both timepoints — is the one we
have since shown does not hold in either of two further cohorts in the same form
(§3.3). We therefore do not claim TPI1 is the most promising of the six. We claim
only that it is the one MR named, which is a statement about instrument
availability, not about biology.

*Pathways*: the denominator is one, and that is the finding rather than an excuse.
Glycolysis was never selected from a set of pathways; no enrichment scan was run
to choose it. It is the annotation of the gene MR named. The self-refuting detail
is that in the functional data TPI1 is not the strongest glycolytic enzyme —
pre-treatment, PGAM1 is (P = 0.0015, the only enzyme at FDR < 0.05), with TPI1
third (P = 0.010) — and PGAM1 carries no genome-wide significant cis-eQTL here, so
MR cannot see it. A reader who takes the named gene for the important gene is
contradicted by our own functional data.

*Mechanisms*: this is a different denominator from the target-level tests above
and the two must not be pooled. Fifteen mechanistic hypotheses about what the
named gene means were formulated and tested; twelve returned a verdict — eight
overturned, three supported by positive results rather than by non-significant
ones (TPI1 marks the glycolytic programme; the class-level avoidance of the
resting state; compartment attribution), one partial — and three returned none.
*Our own methodological claims*: six were pre-registered with failure conditions
written in advance, and three passed (Box 1).

The consequence for Part II is structural. Because gene-level selection is the
part we cannot defend as pre-specified, Part II's load-bearing conclusions are
placed on compartment attribution and on the pathway's instrumentability rate —
neither of which depends on TPI1 being the right choice among the six — and not on
TPI1 as a candidate target.


### 4.3 What survived, and what Part II adds

Everything that survived shares a property: none depends on a single contrast in a
single dataset under a single processing choice. Locus attribution is a
descriptive fact independent of any model, and recurs at five nested power levels
and in a second disease. The instrument-timing observation is a count, not an
inference, and no change of outcome can move it. Compartment attribution holds in three
cohorts — in 16 of 16 and 11 of 11 patients in the two with patient-level
resolution, and in direction in the third — and survives cell-level depth
matching. Conversely, every claim
we withdrew rested on one comparison in one dataset — and the patient-level
stratification, which survived six audits *within* the discovery cohort, is
precisely the claim that did not survive being taken to other cohorts. Surviving
internal audit is not the same as surviving replication, and the distinction is
worth making explicit because internal audit is the more visible of the two.

What Part II adds is two claims, both modest, neither the discovery of a new
biological state.

First, a **design rule about visibility**. For this pathway, instruments are
scarce and state-bound on the exposure side alone: of 28 genes, **3 carry a
genome-wide significant cis-eQTL in any profile**. Two further attritions — 2 of
28 analysable against this outcome, 1 of 28 nominally associated — are joint
properties of the pathway and the outcome GWAS and are reported separately rather
than folded into the first. Which state carries the instrument differs by gene — the three
sit at three different timepoints, one of them the resting state — so we do not
claim, as an earlier version of this manuscript did, that instruments for this
pathway exist only after stimulation; the full instrument matrix (Supplementary
S13) refutes that. What holds is the per-gene statement: for TPI1 a usable
instrument exists at one timepoint and at no other, driven by measurability rather
than by effect size, since at the same variant its effect estimate is largest at
rest and simply far noisier there. We therefore do not claim that regulation
itself is activation-dependent; separating the two requires individual-level
genotype × pseudotime interaction test — which is public, which we now use, and
which returns no evidence of interaction for any TPI1 lead variant (§2.7). The
design consequence is unaffected by which explanation holds: a study using a resting-state resource would find nothing
here and could not distinguish that from a true null. **The generalisable question
for any pathway is not whether its genes have eQTLs, but in which state they can
be measured well enough to yield one.**

Second, a **stratification hypothesis**, and we are explicit that it is one. T
cell metabolism and dysfunction is a mature field; neither glycolysis in exhausted
T cells nor AP-1 at their regulatory elements is new. What is narrow and specific
here is that the same programme, defined by a signature locked across three data
types and verified after removing the one gene that carries an instrument,
separates checkpoint-blockade responders from non-responders in CD4⁺ T cells —
and that in one discovery cohort and two pre-registered follow-ups, no single
treatment arm is confirmed more than once — the pre-treatment effect matching in
magnitude in the same-disease cohort and vanishing in the cross-disease one, the
post-treatment effect reversing in the first and reproducing at full magnitude in
the second (§3.3). That is a hypothesis whose instability has
been mapped, not an established stratifier.


### 4.4 Limitations

Instruments were single cis-eQTLs, so Steiger filtering is near-deterministic and
reported as procedure. Relaxed instrument sets at r² < 0.1 retain correlation that
inflates IVW significance (157 FDR-significant exposures versus 21 in the strict
set); we report IVW–weighted-median concordance rather than IVW discovery.
Rashkin standard errors were reconstructed; Cochran's Q at 5.33% against a 5%
expectation shows no gross miscalibration but does not establish calibration. Exposure eQTLs come from 85–100 donors, so top-variant effects carry
winner's curse. Single-cell cohorts are small; the CD4 population has limited
purity (CD4 detected in 32.5% of cells), addressed by reporting a CD8-negative
sensitivity analysis rather than redefining the population; four patients
contribute two post-treatment samples each, two of them with opposite response
labels, so all patient-level results are reported under three ways of handling the
repetition. The patient-level analysis rests on three small cohorts (9 versus 10 donors
pre-treatment on discovery; 8 versus 3 in replication), and the replication
cohort's non-responder group is small enough that no configuration could have
reached significance at the observed effect size.
Colocalisation results were checked across four analysis windows and four values
of the shared-variant prior: window size is immaterial, but four of eleven
examined loci — including TPI1 — change criterion status with the prior, so for
those loci the reported classification should be read as conditional on
p₁₂ = 10⁻⁵. Two claims in an earlier version
of this manuscript were withdrawn after checks requested during review: a binomial test on the count of concordant
glycolytic genes, which is anticonservative because those genes are correlated
(effective number ≈ 11.7 of 16), and a statement that instruments for this pathway
exist only in activated states, which the full instrument matrix refutes. The
locked 16-gene signature was fixed as the intersection of the three previously
used sets (membership in Supplementary Table S17), which was done after the
response analyses existed; it is a
post-specified but frozen signature, not a pre-registered one, and every
signature-level result is therefore also reported with the instrumented gene
removed. The compartment-attribution result is
now established in cell-resolved data rather than inferred from mixed spots, so
the earlier concern that a CD4-restricted effect might be diluted below detection
at 100 µm resolution no longer applies to that conclusion; it does still bound
what the spatial data alone can show. We were unable to replicate the spatial
analysis in current probe-based Visium datasets for a structural reason rather
than a biological one: the 10x Visium probe panel does not include TPI1, GAPDH,
PKM, LDHA or ALDOA — nor, in one panel version, any HLA class I gene — so neither
the target genes nor our positive control can be measured there. Replication of
the spatial arm requires a poly-A-capture platform. The power-stability curves
rest on reference lists of six and four genes, are not extrapolated above observed
power, and compare against a full-power list that is itself unstable, so they
measure agreement between two imperfect lists rather than recovery of truth. Both alternative explanations for
the colocalisation results were tested and neither is fully closed: the
multiple-signal explanation is excluded for PARP1 but untestable for ZFYVE19 and
TPI1, and the sequential-release experiment removes cohort heterogeneity for the
known-locus part of the list but cannot address the novel part, which is empty
throughout that range. The second-tumour test in §2.2 is limited by the same
quantity it examines: the larger HCC outcome carries 30.3% of melanoma's effective
sample size, so the pre-registered primary enrichment test did not reach
significance (P = 0.110) and the generalisation is reported as directionally
consistent rather than established; its two power levels are one resource rather
than two cohorts, since the published meta-analysis cannot be shown to exclude
FinnGen; and it covers the MR layer only, without colocalisation or SMR. A CD45⁺
HCC single-cell cohort was obtained to extend the compartment result to a second
tumour and could not serve that purpose — it contains no malignant cells, so the
comparison that defines the result cannot be formed there. The same cohort does
carry response labels and was used for the patient-stratification test in §3.3,
where three further limits apply and are stated there rather than absorbed: its
pre-treatment tumour arm could not have reached significance at any effect size
(four responders against two non-responders), its significant post-treatment
result is not significant on a rank test of the same data, and it is a different
cancer under a different regimen, so it strengthens no claim to the degree a
same-disease replication would. Finally, the
colocalisation supporting our worked nomination sits in a region the outcome GWAS
cannot fine-map — stated rather than absorbed.

Finally, the boundary of what this design can establish is reached, and we name it
rather than leave it to be found. **Whether activation-dependence is regulation or
measurement** is the one question of the three that turned out to be answerable
from public data, and it is answered against the regulation reading: no TPI1 lead
variant shows a genotype × pseudotime interaction, while the gene's expression is
strongly pseudotime-dependent (§2.7). What remains open is narrower — the
interaction test covers lead variants under linear and quadratic models, so a
non-monotonic effect at a non-lead variant is not excluded. **Whether the patient
stratification is real** requires a same-disease, same-regimen, CD4-resolved cohort
beyond the single one our pre-fixed search identified; none exists publicly at the
time of writing. **Whether the state has immunological consequences** requires
perturbation, not observation. These two need
new data rather than new analysis, and we would rather say so than present further
re-analyses of the same data as though they closed the gap.

### 4.5 Eight diagnostics, and what should change in practice

Each finding yields a check that costs little and would have changed what this
analysis reported.

**① Annotate the signal against the outcome's own known loci, counting
independent loci rather than gene records.** Here that is 10 of 10 records in the
first round, 3 of 7 loci under the meta outcome (4.1-fold, P = 0.028), no
novel-locus gene at five nested power levels, and the same pattern in a second
disease. In a search-defined sample of 152 eQTL-instrumented MR
target-nomination papers whose full text we could obtain — 469 records screened,
209 eligible, 154 with an open-access full text (Supplementary S23) — **an
estimated 1.6% [0.3–4.9] compare their significant signal against previously
reported loci for their own outcome trait, and 4.7% [1.3–13.5] report how the
candidate list depends on the outcome GWAS used.**

Those figures are corrected ones, and the correction is the point. Automated
full-text matching returned 7.9% and 35.5%. Manual adjudication of a **random 30%
subsample (46 papers), coded against the extracted passages**, put the precision
of that matching at **0.20 and 0.13**: most hits were spurious, three of five on
the first criterion falling inside reference lists, where a cited GWAS title
containing "novel loci" matches, and thirteen of fifteen on the second being
generic statements about statistical power rather than about dependence on the
outcome GWAS. Applying the measured precision to the full corpus gives the
figures above. **The correction runs in both directions**: our earlier 7.9% for
the first criterion was too generous, but our earlier 0.7% for the second — taken
from the strict pattern set — falls *below* the corrected interval and understated
how often that check is reported. Reporting the generous count for one criterion
and the strict count for the other, as we previously did, was not defensible.

Three limits stay attached. Recall was estimated by probing the negatives for
near-miss wording rather than by reading all 152 papers in full, so it bounds the
matcher's sensitivity from above. **No inter-rater statistic is reported**: Cohen's
κ requires two independent coders and this project has one, and computing it from
one person coding twice, or by treating the regex as a coder, would misrepresent
what was done. Most fundamentally, a full-text audit measures **reporting**, and
whether a check was performed but not written up cannot be recovered from text —
so the claim is about what is reported, and not performing a check is not evidence
that a study's conclusions are wrong.

**② Require colocalisation — preferably with explicit multiple-signal modelling
and an LD reference matched to the outcome cohort — rather than treating
SMR/HEIDI as sufficient, and report how many variants entered each HEIDI test.**
We stop short of calling `coloc.abf` the final arbiter: it permits at most one
causal variant per trait per region, and our own attempt to go beyond that
assumption with a proxy-LD SuSiE procedure over-split the MC1R region (19 credible
sets against 3 from in-sample fine-mapping), so it could constrain PARP1 in one
direction only and could not resolve TPI1 or ZFYVE19 at all. HEIDI passed 253 of 291
records that colocalisation assigned to distinct causal variants, on 8–20 variants
per test, in a region that in-sample fine-mapping shows carries ≥3 independent
signals.

**③ Do not expect power to resolve the disagreement.** It moves MR and
colocalisation in opposite directions, because one depends on a single variant's z
and the other on the shape of the regional signal.

**④ Report the physical distance between GWAS and eQTL peaks alongside
posteriors, and test whether posterior movement survives explicit modelling of
multiple causal signals before calling it a change in resolution.** For PARP1 the
collapse is not an unmodelled second signal; for ZFYVE19 and TPI1 no credible set
exists at either power, so those regions cannot adjudicate it.

**⑤ Compute the nominated gene's cell-type expression ratio in an annotated atlas
before citing tissue-level data as validation.** This costs minutes. Ours is
6.7-fold higher in malignant cells than in CD4⁺ T cells in the first cohort
(16 of 16 patients), 2.4-fold in the second (11 of 11), concordant in direction in
a third, and the difference survives cell-level depth matching — so tissue-level
measurements of it are measurements of tumour.

**⑥ Report instrument availability at each of the three levels it passes through,
and do not compress them.** Of 28 glycolytic genes, **3 carry a genome-wide
significant cis-eQTL in this exposure resource, 2 are analysable against this
outcome after harmonisation, and 1 yields a nominal MR association.** Only the
first is a property of the pathway and the exposure data; the second and third are
joint properties of the pathway and the particular outcome GWAS, and collapsing
them into "one usable instrument" attributes an outcome-side limitation to the
biology.

**⑦ When splitting cells by a score, match on lineage composition as well as
depth, at the cell level.** Splitting by score splits purity; depth matching does
not fix it; our cluster-level control passed while the cell-level control failed.

**⑧ Before comparing candidate lists across releases of the same GWAS resource,
verify code by code that the endpoint definition is the same one.** Release notes
are not sufficient. FinnGen R13 describes the successor to our endpoint as
"including Hilmo", its hospital discharge register, which reads as a newly widened
case definition; the code-level definitions show the cases are identically defined
and that what changed is the rule screening cancers out of the controls. We drew
the wrong conclusion from that one line ourselves before checking the definitions
(Supplementary S12). The cost of not checking is specific: a study that reads the
summary line and compares anyway will score pure phenotype drift as instability of
its candidate list — the very quantity such a comparison exists to measure — and
will do so in the direction that makes its own list look less reproducible than it
is. Where the definitions genuinely differ, the comparison is still worth making,
but it is a transfer test rather than a power point, and the direction of every
difference should be written down before the result is seen. Ours transfers
intact: the same six genes at the same two loci, no novel-locus gene, 3.7- to
17.6-fold attribution across four exposure-by-disease pairs, while not one of
3,434 underlying test statistics is unchanged and 51 of 275 nominally significant
records are replaced.

More generally: a candidate list produced by this framework at current outcome
power should be reported as *the set of genes that passed screening under these
conditions*, not as targets — and because the fragile part of such a list is
precisely its novel-locus part, the headline of a study of this kind is the part
least likely to replicate. What the same data support instead is a statement about
*in which cell state a pathway's regulation can be measured precisely enough to
yield an instrument at all*, and that statement does not degrade with outcome
power because it never depends on it.


---

## 5. Methods

> 结构说明（中文，不进正文）
> 顺序按分析链条，不按数据来源。凡是"跑之前写死"的判据（判废线、窗口、方向）一律写进对应小节，
> 并在 §18 单列纪律小节 —— 那一节是本文方法学定位的一部分，不要删。
> 已知的坑（GRCh38 不 liftover、`nzchar(NA)`、peak 必须按区间对齐、R² 无上界等）写在会影响
> 结果可复现性的位置。⟨⟩ 为待补。

---

### 5.1 Exposure data: dynamic CD4⁺ T cell cis-eQTLs

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

### 5.2 Outcome GWAS and meta-analysis

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

### 5.3 Instrument selection and harmonisation

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

### 5.4 Mendelian randomization

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

### 5.5 Colocalisation

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

#### 5.5b Sequential-release power trajectory

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

#### 5.5c FinnGen R13 transfer test

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

### 5.6 SMR and HEIDI

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

### 5.7 Locus annotation and enrichment

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

### 5.8 Negative-control phenotype and comorbidity MR

Because vitiligo GWAS with full summary statistics were unavailable (all eleven
catalogued studies lack summary statistics; FinnGen `L12_VITILIGO` has 391
cases), the pigmentation-pathway negative control was FinnGen
`CD2_BENIGN_MELANOCYTIC` (melanocytic naevi, 13,357 cases): a candidate that
affects melanoma but not naevi is unlikely to act through the pigmentation
pathway. Test sensitivity is established by five known pigmentation loci, all of
which do affect naevi. Additional comorbidity outcomes (rheumatoid arthritis,
psoriasis, type 1 diabetes, a composite autoimmune phenotype, and vitiligo with
its case count stated) were analysed identically and are reported as exploratory.

### 5.9 Cross-cancer analysis

The exposure side was held completely fixed (3,556 exposures, 2,141 instruments)
and only the FinnGen outcome was exchanged, for lung, colorectal, pancreatic,
breast and prostate cancer. Two quantities were compared across outcomes: the
known-locus enrichment of the FDR-significant list (§7, counted by independent
locus), and the effect estimates for the candidate genes.

### 5.9b Second-tumour generalisation (hepatocellular carcinoma)

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

### 5.10 Power–stability simulation

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

### 5.11 Single-cell ICB cohorts

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

### 5.11b Second-tumour patient stratification (GSE235863)

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

### 5.12 Module scores, response testing and lineage control

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

### 5.13 Spatial transcriptomics

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

### 5.14 Purified CD4⁺ T cell multiome (GSE282266)

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

### 5.15 Chromatin accessibility and motif enrichment

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

### 5.16 Peripheral blood variance decomposition (GSE199994)

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

### 5.17 Public eQTL resource survey

Replication resources were surveyed across the eQTL Catalogue (758 datasets, of
which 58 are stimulated T cell datasets) and DICE. Two query behaviours matter:
gene-based queries return position-sorted, paginated results, so the minimum P
value on one page is meaningless and queries must be made by rsID; and a gene
absent from a dataset returns HTTP 400 rather than an empty result, which is not
a null result but a quantification decision by the resource. Datasets were also
checked for identity with the exposure data by comparing design and per-timepoint
sample sizes, which excluded one resource as circular.

### 5.18 Pre-specification, positive controls and stopping rules

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

### 5.19 Software and computing environment

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

### 5.20 Data and code availability

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

---

## 6. Figure legends

**Fig 1 | Locus attribution.** Manhattan plot of MR P values for all strict-set
exposures under the meta outcome, coloured by locus category. FDR-significant
genes labelled, grouped within 1.5 Mb. Signal concentrates at known pigmentation
and naevus loci.

**Fig 2 | HEIDI does not reject the LD confounding colocalisation identifies.**
Colocalisation PP.H4 against −log₁₀P_HEIDI; points in the lower-left quadrant are
records assigned to distinct causal variants that HEIDI passes. MC1R cluster
labelled.

**Fig 3 | Posterior support moves with the outcome dataset, in both
directions.** Regional plots for PARP1 and ZFYVE19: cis-eQTL, GWAS at FinnGen
power, GWAS at meta power, with peak distances annotated. The two rounds differ
in effective resolution and in cohort composition; the multiple-signal analysis
adjudicates only PARP1, so this figure is not labelled as power determining
resolution (§2.4).

**Fig 4 | Outcome power sets the reproducibility of the candidate list.**
**a**, Down-sampling curve with the separately observed FinnGen round marked;
that round is a component of the meta outcome, so the calibration is empirical
rather than external (§2.5).
**b**, Five cancers. **c**, Stratified by locus class: known loci recover at 46%
where novel loci are at 0.4%; 50% recovery needs 2,506 versus 8,771 cases.
**d**, The same stratification against full-power |z|. The two classes barely
overlap (novel 3.73–4.55; known to 15.99), and a |z|-only model with no class
label reproduces 68–80% of the gap, so panel **c** is read as threshold proximity
rather than as a property of the category (§2.5, Supplementary S28).

**Fig 5 | Compartment attribution: a gene nominated with CD4-specific instruments
is measured, in tissue, as tumour.** **a**, TPI1 by annotated cell type in
melanoma single-cell data (GSE115978); percentages give the fraction of cells with
TPI1 detected. **b**, Malignant versus CD4⁺ T cells paired within patients, in the
full data, after cell-level matching on sequencing depth, and in an independent
cohort; counts give patients with malignant > CD4⁺ T. Note that GSE72056 uses a
different normalisation scale, so comparisons are within cohort. **c**,
Pre-specified positive controls. **d**, Spatial context: meta-analytic Spearman
correlations across eight sections for TPI1 and the HLA-C mirror control.
**e**, TCGA-SKCM overall survival: bulk TPI1 is associated with worse survival
until immune content is adjusted for. Generated by
`figures/make_fig5_compartment.py`.

**Fig 6 | Cross-cancer comparison.** **a**, Known-locus enrichment by independent
locus across six outcomes. **b**, Candidate effects across cancers. **c**, Hit
count against outcome case count.

**Fig 7 | FinnGen sequential releases.** **a**, Discoveries versus case number
with registered prediction intervals. **b**, Known- versus novel-locus recovery
against registered predictions. **c**, Observed versus predicted.

**Fig 8 | Multi-layer evidence matrix.** Six candidates against seven evidence
layers; no candidate is strong on all axes.

**Fig 9 | What the instrument shows, and what it does not.**
**a**, cis-eQTL effect on TPI1 at the same variant across the activation time
course, as β with 95% CI and the P value below each estimate. The 16 h P value is
smallest because the standard error is roughly five-fold smaller; in naive cells
the point estimate is largest at rest. Shading marks the only timepoint carrying a
genome-wide significant instrument for this gene. **b**, Median lead-variant
standard error of the locked 16-gene signature relative to background, by
timepoint: the signature is measured more precisely than background at every
activated timepoint and barely so at rest, which mechanically moves each gene's
strongest eQTL away from the resting state. **c**, Fraction of glycolytic enzymes
higher in non-responders in patient CD4⁺ T cells, for the full 22-gene set, the
locked 16-gene signature, and the signature with TPI1 removed, at both timepoints;
**shown descriptively — inference is by composite signature score and label
permutation (Supplementary Table S19), because these genes are correlated**.

**d**, Family composition of the top of the residualised multiome axis: glycolysis
is enriched 21.3-fold (15 genes) among the axis-defining genes, and none of the
eleven biological modules tested exceeds a composition-matched null, so the axis is
a definable state rather than a restatement of activation intensity. **e**, Motif
enrichment in axis-high versus background peaks, restricted to
activation-invariant peaks. ⚠ Panel **e** is an axis-level chromatin
characterisation and is **not** the claim, withdrawn during revision, that the
instrument acts by disrupting an AP-1 motif — that test returned an empirical
P = 1.0 and is not shown anywhere in this paper.

Panels **d** and **e** were added in revision: an earlier version of this figure
omitted the multiome evidence and said so in its legend, which a reviewer
correctly identified as declaring the gap rather than closing it.


---

## 7. Supplementary

**Figures** (`figures/make_figS1_flow.py`, `figures/make_supp_figures.py`)

| | Content |
|---|---|
| **Fig S1** | Analysis flow with the actual counts at each step, the stopping rule that operated there, and the tests discarded by their own controls |
| **Fig S2** | Steiger filtering: the unbounded R² of the published SD-unit formula (max 1.027, 24 records > 0.8) and the bounded form we report instead |
| **Fig S3** | Multi-instrument estimators, and why weighted mode is reported as inapplicable (1 significant of 226, at a median of 4 instruments) |
| **Fig S4** | Naevi as a pigmentation-pathway negative control; test sensitivity established by the positive controls |
| **Fig S5** | A score-based split also splits lineage purity: effects collapse across three levels of control |
| **Fig S6** | The two processing errors we found in our own completed analyses, with the diagnostics that caught them and every patient-level number before and after |
| **Fig S7** | The metabolic axis: family composition of the top-ranked genes, and motif enrichment before and after restricting to activation-invariant peaks |
| **Fig S8** | The precondition test that closed a planned analysis: ICC against a matched random-module floor, with the XIST upper-bound control |

**Tables and documents**

| | Content | Source |
|---|---|---|
| S9 | **Pre-registration document** (power trajectory), with its modification log | `PREREG_power_trajectory.md` |
| S10 | Record-level locus enrichment (locus-level is reported in the main text) | `36d`, `36e` |
| S11 | Per-release candidate lists across the FinnGen trajectory | `59b` |
| S12 | **Complete timeline of target-substantiation attempts**: inclusion rules, the three attempts that returned support, and the **gene / pathway / mechanism selection denominators** (§5) | `SUPP_attempt_timeline.md` |
| S13 | Instrument availability matrix: 28 glycolytic genes × 8 activation profiles | `75a` |
| S14 | Replication-cohort search: five queries, twelve datasets assessed against pre-fixed criteria | `SUPP_replication_search.md` |
| S15 | Colocalisation sensitivity: 11 loci × 16 window-and-prior combinations | `81a` |
| S16 | Item-by-item comparison with three studies using the same framework | `SUPP_comparator_studies.md` |
| S17 | Locked glycolytic signature: set membership across the three data types | `63a` |
| S18 | **Pre-registration document** (independent-cohort replication), with its reading table | `PREREG_pozniak_replication.md` |
| S19 | Patient-level inference: composite signature score, correlation-preserving label permutation, and all three treatments of repeated patients | `74a` |
| S20 | **Pre-registration document** (second-tumour generalisation to HCC), with its five-cell reading table, results register and two logged deviations; accompanying tables: known-locus reference list, per-outcome MR results, significant-locus attribution, matched-power simulation, and the matched-background version of the enrichment test | `PREREG_hcc_generalisation.md`; `84a`, `84b`, `85a`–`85e` |
| S22 | **Pre-registration document** (exposure-resource generalisation), with its reading table and results register; accompanying tables: eQTLGen instruments, MR records against the fixed melanoma outcome, locus attribution, and the matched-background version | `PREREG_exposure_resource.md`; `92a`–`92e` |
| S24 | **Pre-registration document** (disease × exposure-resource grid), with its reading table, the mismatched-list negative control, and two logged deviations: the grid was reduced from seven diseases to two because per-disease known-locus coordinates could not be resolved for the rest, and the negative control was reformulated | `PREREG_generality_grid.md`; `94d`–`94f` |
| S23 | Search-defined audit of 152 eQTL-MR target-nomination papers: four PubMed queries, eligibility rule, five pre-fixed scoring criteria, a logged scoring bug and a logged post-hoc broadening, the PRISMA-style flow from 469 screened records, and the **validation of the automated coder against manual adjudication of a fixed-seed random 30% subsample** — per-paper codes, the passages each was judged on, measured precision and recall, and the precision-corrected field-level estimates that supersede both the generous and the strict counts | `SUPP_literature_audit.md`; `91a`–`91e`, `104a`–`104c` |
| S21 | **Pre-registration document** (patient stratification in a second tumour type), including the minimum attainable P value of each arm computed from the sample structure before any expression value was read, the positive controls, and the three-cohort comparison; accompanying tables: sample-level scores, per-arm results and positive controls | `PREREG_hcc_part2_generalisation.md`; `87a`–`87c` |
| S29 | **Pre-registration document** (HCC known-locus de-circularisation), with the per-variant GWAS Catalog provenance of all ten outcome-sourced loci, the three references scored, and two items its section 0.3 fixed in advance: that the test could not have changed the result because the circular loci lie far from any significant one, and that over-correction raises the enrichment and must not be cited as a stricter standard | `PREREG_hcc_decircularisation.md`; `102a` |
| S28 | **Pre-registration document** (effect-size matching), registered before any matched analysis was run, including the finding — stated in its section 0 — that the down-sampling model contains no class label and the differential can therefore only follow from the |z| distributions; the pre-committed replacement wording; and two items logged against ourselves: a process control that failed on first run through a random-number-ordering bug, and a design flaw in the pre-registration itself, which permitted the matched comparison to collapse onto a single known anchor | `PREREG_effect_size_matching.md`; `101a`–`101c` |
| S27 | **Self-administered attribution check**: the accepted-causal-gene list assembled before the comparison, the ten loci reached under the eQTLGen exposure, and the named gene at each | `96a` |
| S26 | **Multiple-testing unit sensitivity**: the FDR < 0.05 list recomputed over records, unique variants, genes and independent loci, each by minimum-p and by Simes combination, plus a two-stage hierarchical procedure; reported for both the meta and FinnGen rounds, with the gene list under each unit | `100a`, `100b` |
| S25 | **Pre-registration document** (FinnGen R13 transfer test), with the endpoint mapping that excludes R13 from the power trajectory, the pre-committed direction of every difference between the two releases, six point predictions with intervals, three process controls, the mismatched-list negative control, and its results register — including four items registered against ourselves: a point prediction that did not come true although it fell inside its interval, two identical cells that are arithmetic rather than confirmation, an identical count over non-identical record sets, and a comparator-identity bug of ours that compared a meta-analysis to a release | `PREREG_r13_transfer.md`; `99a`–`99c` |

---

## 8. References

**核实状态**：2026-08-11 通过 CrossRef API 与 NCBI E-utilities 逐条核实。
下列条目的作者、期刊、年、卷、页、DOI **均来自查询返回**，非从记忆填写。

⚠ 核实过程中发现的坑，记下以免复查时重蹈：
- CrossRef 的 `query.bibliographic` 会把 **Faculty Opinions 的推荐短文**、**预印本**、
  **会议摘要**排在前面（Sade-Feldman、Jerby-Arnon、Tirosh、DICE、TCGA-CDR 均如此），
  须过滤 `type == "posted-content"` 与容器名含 "Faculty Opinions"/"Abstract"
- 少数标题检索会返回完全不相干的文献（Landi 2020 曾匹配到蜜蜂抗螨、FinnGen 曾匹配到核桃树病原菌）
- PubMed 检索 FinnGen 时**更正声明（Author Correction）排在原文之前**，须取后者
- 仍待补：**GSE282266 与 GSE199994 的关联论文**（项目中只用了 GEO 号），
  以及 GSE316760 / GSE300445 的著录

---

## A. 同框架对照研究

1. Zheng J, Yang Q, Liu H, et al. Integrating single-cell transcriptome-wide Mendelian
   randomization and differentially expressed gene analyses to prioritize dynamic
   immune-related drug targets for cancers. *Advanced Science* 2025;12:e07451.
   doi:10.1002/advs.202507451

2. Wu X, Ying H, Yang Q, et al. Transcriptome-wide Mendelian randomization during
   CD4⁺ T cell activation reveals immune-related drug targets for cardiometabolic
   diseases. *Nature Communications* 2024;15:9302. doi:10.1038/s41467-024-53621-7

3. Cui K, Zou Q, Qu X, et al. Transcriptome-wide Mendelian randomization and
   single-cell analysis during CD4⁺ T cell activation deciphers immunotherapeutic
   targets for colorectal cancer. *npj Precision Oncology* 2025;10:32.
   doi:10.1038/s41698-025-01236-6

> 逐项对照证据见 `SUPP_comparator_studies.md`；正文对这三篇的主张限定为"正文中未见"。

## B. 暴露与结局

4. **Soskic B, Cano-Gamez K, Smyth DJ, et al. Immune disease risk variants regulate
   gene expression dynamics during CD4⁺ T cell activation. *Nature Genetics*
   2022;54:817–826. doi:10.1038/s41588-022-01066-3** ← 本文暴露数据来源

5. Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a
   well-phenotyped isolated population. *Nature* 2023;613:508–518.
   doi:10.1038/s41586-022-05473-8
   （另见 Author Correction: *Nature* 2023;615:E19. doi:10.1038/s41586-023-05837-8）

6. Rashkin SR, Graff RE, Kachuri L, et al. Pan-cancer study detects genetic risk
   variants and shared genetic basis in two large cohorts. *Nature Communications*
   2020;11:4423. doi:10.1038/s41467-020-18246-6　（GWAS Catalog GCST90011809）

7. Landi MT, Bishop DT, MacGregor S, et al. Genome-wide association meta-analyses
   combining multiple risk phenotypes provide insights into the genetic architecture
   of cutaneous melanoma susceptibility. *Nature Genetics* 2020;52:494–504.
   doi:10.1038/s41588-020-0611-8　（GCST010302 / GCST010303 / GCST010304）

8. Nathan A, Asgari S, Ishigaki K, et al. Single-cell eQTL models reveal dynamic
   T cell state dependence of disease loci. *Nature* 2022;606:120–128.
   doi:10.1038/s41586-022-04713-1　（非对照研究；为 §3.1 核查过的资源之一）

## C. 功能与验证数据集

9. Sade-Feldman M, Yizhak K, Bjorgaard SL, et al. Defining T cell states associated
   with response to checkpoint immunotherapy in melanoma. *Cell* 2018;175:998–1013.e20.
   doi:10.1016/j.cell.2018.10.038　（GSE120575）

10. Jerby-Arnon L, Shah P, Cuoco MS, et al. A cancer cell program promotes T cell
    exclusion and resistance to checkpoint blockade. *Cell* 2018;175:984–997.e24.
    doi:10.1016/j.cell.2018.09.006　（GSE115978）

11. Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of
    metastatic melanoma by single-cell RNA-seq. *Science* 2016;352:189–196.
    doi:10.1126/science.aad0501　（GSE72056；本文独立复现队列）

12. Hugo W, Zaretsky JM, Sun L, et al. Genomic and transcriptomic features of response
    to anti-PD-1 therapy in metastatic melanoma. *Cell* 2016;165:35–44.
    doi:10.1016/j.cell.2016.02.065　（GSE78220）

13. Riaz N, Havel JJ, Makarov V, et al. Tumor and microenvironment evolution during
    immunotherapy with nivolumab. *Cell* 2017;171:934–949.e16.
    doi:10.1016/j.cell.2017.09.028　（GSE91061）

14. Thrane K, Eriksson H, Maaskola J, et al. Spatially resolved transcriptomics enables
    dissection of genetic heterogeneity in stage III cutaneous malignant melanoma.
    *Cancer Research* 2018;78:5970–5979. doi:10.1158/0008-5472.CAN-18-0747

15. Katko A, Potter SJ, et al. Gene regulatory network determinants of rapid recall in
    human memory CD4⁺ T cells. *Cell Reports* 2026;45(4):117103.
    doi:10.1016/j.celrep.2026.117103 — data: NCBI GEO **GSE282266**. [已核实]

16. Boukhaled GM, Gadalla R, et al. Pre-encoded responsiveness to type I interferon in
    the peripheral immune system defines outcome of PD1 blockade therapy.
    *Nature Immunology* 2022;23(8):1273–1283. doi:10.1038/s41590-022-01262-7
    — data: NCBI GEO **GSE199994**. [已核实]

17. Virós A, et al. Spatial transcriptomics of primary cutaneous melanoma.
    NCBI GEO **GSE316760** (submitted 2026-01-16); no associated publication indexed
    at the time of writing. [已核实为数据集，无关联论文]

17b. Pham F, Dufeu M, Benboubker V, et al. Spatial tumour-immune ecosystems shape the
    efficacy of anti-PD1 immunotherapy in primary cutaneous melanoma
    [Spatial Transcriptomics]. NCBI GEO **GSE300445** (submitted 2025-06-23);
    no associated publication indexed at the time of writing.
    [已核实为数据集，无关联论文]

17c. Guo X, et al. Contrasting cytotoxic and regulatory T cell responses underlying
    distinct clinical outcomes to anti-PD-1 plus lenvatinib therapy in hepatocellular
    carcinoma. *Cancer Cell* 2025;43(2):248–268.e9. doi:10.1016/j.ccell.2025.01.001
    — data: NCBI GEO **GSE235863**. [已核实]

17d. Ghouse J, Gellert-Kristensen H, O'Rourke CJ, et al. Genome-wide meta-analysis
    identifies nine loci associated with higher risk of hepatocellular carcinoma.
    *JHEP Reports* 2025;7(9):101485. doi:10.1016/j.jhepr.2025.101485
    — GWAS Catalog **GCST90809296**. [已核实]

17e. Võsa U, Claringbould A, Westra H-J, et al. Large-scale cis- and trans-eQTL analyses
    identify thousands of genetic loci and polygenic scores that regulate blood gene
    expression. *Nature Genetics* 2021;53:1300–1310. doi:10.1038/s41588-021-00913-z
    — eQTLGen Consortium, 2019-12-11 cis-eQTL release (n = 31,684). [已核实]

17f. FinnGen. Release R12 (2024), endpoints `C3_MELANOMA_SKIN_EXALLC` and
    `C3_HEPATOCELLU_CARC_EXALLC`. https://www.finngen.fi/en/access_results
    [数据集，按 FinnGen 引用规范著录]

## D. 参考面板、方法与工具

18. Byrska-Bishop M, Evani US, Zhao X, et al. High-coverage whole-genome sequencing of
    the expanded 1000 Genomes Project cohort including 602 trios. *Cell*
    2022;185:3426–3440.e19. doi:10.1016/j.cell.2022.08.004

19. Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation
    between pairs of genetic association studies using summary statistics.
    *PLoS Genetics* 2014;10:e1004383. doi:10.1371/journal.pgen.1004383

20. Wang G, Sarkar A, Carbonetto P, Stephens M. A simple new approach to variable
    selection in regression, with application to genetic fine mapping. *Journal of the
    Royal Statistical Society Series B* 2020;82:1273–1300. doi:10.1111/rssb.12388

21. Zhu Z, Zhang F, Hu H, et al. Integration of summary data from GWAS and eQTL studies
    predicts complex trait gene targets. *Nature Genetics* 2016;48:481–487.
    doi:10.1038/ng.3538

22. Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic
    causal inference across the human phenome. *eLife* 2018;7:e34408.
    doi:10.7554/eLife.34408

23. Chang CC, Chow CC, Tellier LC, et al. Second-generation PLINK: rising to the
    challenge of larger and richer datasets. *GigaScience* 2015;4:s13742-015-0047-8.
    doi:10.1186/s13742-015-0047-8

24. Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal
    and scalable single-cell analysis. *Nature Biotechnology* 2024;42:293–304.
    doi:10.1038/s41587-023-01767-y

25. Fornes O, Castro-Mondragon JA, Khan A, et al. JASPAR 2020: update of the
    open-access database of transcription factor binding profiles. *Nucleic Acids
    Research* 2020;48:D87–D92. doi:10.1093/nar/gkz1001

26. Schep AN, Wu B, Buenrostro JD, Greenleaf WJ. chromVAR: inferring
    transcription-factor-associated accessibility from single-cell epigenomic data.
    *Nature Methods* 2017;14:975–978. doi:10.1038/nmeth.4401

## E. 资源与数据门户

27. Schmiedel BJ, Singh D, Madrigal A, et al. Impact of genetic polymorphisms on human
    immune cell gene expression. *Cell* 2018;175:1701–1715.e16.
    doi:10.1016/j.cell.2018.10.022　（DICE）

28. Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human
    gene expression and splicing quantitative trait loci. *Nature Genetics*
    2021;53:1290–1299. doi:10.1038/s41588-021-00924-w　（eQTL Catalogue）

29. Sollis E, Mosaku A, Abid A, et al. The NHGRI-EBI GWAS Catalog: knowledgebase and
    deposition resource. *Nucleic Acids Research* 2023;51:D977–D985.
    doi:10.1093/nar/gkac1010

30. Ochoa D, Hercules A, Carmona M, et al. The next-generation Open Targets Platform:
    reimagined, redesigned, rebuilt. *Nucleic Acids Research* 2023;51:D1353–D1359.
    doi:10.1093/nar/gkac1046

31. Sun D, Wang J, Han Y, et al. TISCH: a comprehensive web resource enabling
    interactive single-cell transcriptome visualization of tumor microenvironment.
    *Nucleic Acids Research* 2021;49:D1420–D1430. doi:10.1093/nar/gkaa1020

32. Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics
    data via the Xena platform. *Nature Biotechnology* 2020;38:675–678.
    doi:10.1038/s41587-020-0546-8

33. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data
    resource to drive high-quality survival outcome analytics. *Cell* 2018;173:400–416.e11.
    doi:10.1016/j.cell.2018.02.052

34. Tsherniak A, Vazquez F, Montgomery PG, et al. Defining a cancer dependency map.
    *Cell* 2017;170(3):564–576.e16. doi:10.1016/j.cell.2017.06.010 — portal:
    DepMap, Broad Institute, https://depmap.org (release used: 22Q2 common-essential
    and gene-effect calls). [已核实；⚠ 具体 release 号须与分析脚本核对后定稿]

---

## 投稿前仍须处理

1. ~~条目 15–17、34 的著录~~ ✅ 2026-08-12 完成（含新增 17c–17f：GSE235863、
   GCST90809296、eQTLGen、FinnGen R12）。**唯一残留**：条目 34 的 DepMap release 号
   须与脚本核对
2. 统一到目标期刊的著录格式（本文件为作者-年-卷-页-DOI 的中性格式）
3. 若期刊要求，补全 3 位以上作者的完整列表

## E. 对照研究（2026-08-13 新增，全部经 CrossRef 核实）

35. Reales G, et al. Design and interpretation of eQTL–GWAS colocalisation studies:
    lessons from a large-scale evaluation. *PLOS Genetics* 2026.
    doi:10.1371/journal.pgen.1012141 [已核实]

36. Tambets R, et al. Extensive co-regulation of neighboring genes complicates the
    use of eQTLs in target gene prioritization. *Human Genetics and Genomics
    Advances* 2024;5:100348. doi:10.1016/j.xhgg.2024.100348 [已核实]

37. Rosen JD, et al. Higher eQTL power reveals signals that boost GWAS
    colocalization. *American Journal of Human Genetics* 2026.
    doi:10.1016/j.ajhg.2026.02.009 [已核实]

38. Howe LJ, et al. Evaluating transportability of in vitro cellular models to
    in vivo human phenotypes using gene perturbation data. *Nature Communications*
    2025. doi:10.1038/s41467-025-67199-1 [已核实]

39. Lin Z, Pan W. A robust cis-Mendelian randomization method with application to
    drug target discovery. *Nature Communications* 2024.
    doi:10.1038/s41467-024-50385-y [已核实]

40. Karhunen V, et al. Integrating genetic data with biological insight: a
    practical guide to cis-Mendelian randomization. *American Journal of Human
    Genetics* 2026. doi:10.1016/j.ajhg.2026.03.011 [已核实]

⚠ 正文中标 `[ref]` 的六处引用点对应 35–40，装配前须替换为期刊要求的编号格式。
