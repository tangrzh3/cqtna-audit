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
and it transfers unchanged to that resource's current release although not one of
the 3,434 underlying test statistics does;
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
by independent locus. **We declare the record level primary**, for the sole reason
that it is what every analysis in this paper was built on — the registered
down-sampling predictions, the release trajectory, both generalisations and the
grid. Re-declaring it now, after seeing which unit yields the longest list, is
the selective emphasis this paper exists to detect.

Because that choice cannot be defended by argument alone, we recomputed the
FDR < 0.05 list under every unit the paper uses — variant, gene and independent
locus, each by both minimum-p and Simes combination, plus a two-stage
hierarchical procedure selecting genes and then records within them
(Supplementary S26). Three results matter. First, **the record level is the most
conservative of the seven alternatives**: no gene is lost under any other unit,
and between 1 and 43 are added, so every significance claim in this paper is a
subset of what a looser unit would license. Second, **the quantity the argument
actually runs on barely moves** — 7 to 9 independent loci against the 7 reported,
and 2 under every single unit in the FinnGen round. Third, and decisively for
Part I, **the attribution conclusion is unit-independent**: the known-locus share
of significant loci ranges 42.9–55.6% against the 42.9% reported, and is 100% at
every unit in the FinnGen round.

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
artefact of reusing data.** Not one of the 3,434 test statistics is identical
between the two releases (correlation of |z| = 0.937), standard errors shrink by
3.85% against the 3.8% the added cases predict, and 51 of the 275 nominally
significant records are replaced. What holds still is the FDR-significant tier;
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
the effect is specific to each outcome's own genetics and is not an artefact of
locus density.

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
order of one to two per cent of the tissue-level TPI1 signal. **A tissue-level
measurement of this gene is a measurement of tumour glycolysis**, whatever its
P value.

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
target-nomination papers whose full text we could obtain, **at most 12 (7.9%)
compare their significant signal against previously reported loci for their own
outcome trait, and 1 (0.7%) reports how the candidate list depends on the outcome
GWAS used** (Supplementary S23; three of those studies are also compared item by
item in Supplementary S16). We are explicit that this counts phrases in text, not quality of
practice, and that not performing a check is not evidence that a study's
conclusions are wrong.

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

完整 Methods 见 `manuscript/METHODS_draft.md`（20 节 + §5b）。
定稿前只改该文件再重新组装，避免两处分叉。

§5.1 暴露数据 · §5.2 结局与 meta · §5.3 工具变量 · §5.4 MR 与 Steiger ·
§5.5 共定位（含多信号敏感性与 MC1R 校准对照）· §5.5b **release 功效轨迹（含预注册）** ·
§5.5c **R13 转移测试（含预注册 S25、endpoint 定义核对、提取缓存）** ·
§5.6 SMR/HEIDI · §5.7 位点标注与富集（按独立位点计）· §5.8 阴性对照表型 ·
§5.9 跨癌种 · §5.9b **第二瘤种泛化（HCC，含预注册与两处偏离）** ·
§5.10 功效-稳定性模拟 · §5.11 单细胞队列（归一化判据、样本定义）· §5.11b **第二瘤种患者分层（GSE235863，含预注册）** ·
§5.12 模块评分与谱系对照 · §5.13 空间转录组 · §5.14 纯化 CD4 multiome ·
§5.15 染色质与 motif · §5.16 外周血方差分解 · §5.17 公共 eQTL 资源核查 ·
§5.18 **预设、阳性对照与停止规则** · §5.19 软件 · §5.20 数据与代码

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
| S23 | Search-defined audit of 152 eQTL-MR target-nomination papers: four PubMed queries, eligibility rule, five pre-fixed scoring criteria, a logged scoring bug and a logged post-hoc broadening, and a manual false-negative spot-check | `SUPP_literature_audit.md`; `91a`–`91e` |
| S21 | **Pre-registration document** (patient stratification in a second tumour type), including the minimum attainable P value of each arm computed from the sample structure before any expression value was read, the positive controls, and the three-cohort comparison; accompanying tables: sample-level scores, per-arm results and positive controls | `PREREG_hcc_part2_generalisation.md`; `87a`–`87c` |
| S29 | **Pre-registration document** (HCC known-locus de-circularisation), with the per-variant GWAS Catalog provenance of all ten outcome-sourced loci, the three references scored, and two items its section 0.3 fixed in advance: that the test could not have changed the result because the circular loci lie far from any significant one, and that over-correction raises the enrichment and must not be cited as a stricter standard | `PREREG_hcc_decircularisation.md`; `102a` |
| S28 | **Pre-registration document** (effect-size matching), registered before any matched analysis was run, including the finding — stated in its section 0 — that the down-sampling model contains no class label and the differential can therefore only follow from the |z| distributions; the pre-committed replacement wording; and two items logged against ourselves: a process control that failed on first run through a random-number-ordering bug, and a design flaw in the pre-registration itself, which permitted the matched comparison to collapse onto a single known anchor | `PREREG_effect_size_matching.md`; `101a`–`101c` |
| S27 | **Self-administered attribution check**: the accepted-causal-gene list assembled before the comparison, the ten loci reached under the eQTLGen exposure, and the named gene at each | `96a` |
| S26 | **Multiple-testing unit sensitivity**: the FDR < 0.05 list recomputed over records, unique variants, genes and independent loci, each by minimum-p and by Simes combination, plus a two-stage hierarchical procedure; reported for both the meta and FinnGen rounds, with the gene list under each unit | `100a`, `100b` |
| S25 | **Pre-registration document** (FinnGen R13 transfer test), with the endpoint mapping that excludes R13 from the power trajectory, the pre-committed direction of every difference between the two releases, six point predictions with intervals, three process controls, the mismatched-list negative control, and its results register — including four items registered against ourselves: a point prediction that did not come true although it fell inside its interval, two identical cells that are arithmetic rather than confirmation, an identical count over non-identical record sets, and a comparator-identity bug of ours that compared a meta-analysis to a release | `PREREG_r13_transfer.md`; `99a`–`99c` |

---

## 8. References

References are maintained in `REFERENCES.md` and are assembled into this document
by `assemble.py` at build time. All entries have been verified against CrossRef
and PubMed; no volume, issue or page number has been supplied from memory.

