**⚠ 已被 `MANUSCRIPT_GB.md` 取代（2026-08-11 双主线重构）。勿引用本文件的数字。**

# Dynamic eQTL Mendelian randomization for immune target discovery: an audit in melanoma

> ⚠ **本文件已被 `MANUSCRIPT_v2_dual_thread.md` 取代（2026-08-11 双主线重构）。**
> 保留仅供对照单主线版本的措辞。**勿在此文件上继续改稿。**

**Draft v1 — 2026-08-11**
组装自 `INTRODUCTION_draft.md` / `RESULTS_draft.md` / `METHODS_draft.md` /
`DISCUSSION_draft.md`，数字口径已统一（位点级 4.09×；校准点 10.2 vs 10）。
正文中的 ⟨…⟩ 为待补项，末尾另有清单。

> **备选标题**
> 1. *The reproducible part of a Mendelian randomization candidate list is the part that is not a discovery*
> 2. *Outcome GWAS power, not exposure resolution, sets what dynamic-eQTL Mendelian randomization can find*
> 3. *Seven diagnostics for dynamic-eQTL Mendelian randomization studies of cancer, and what they returned in melanoma*

---

## Abstract

Dynamic expression quantitative trait loci (eQTLs) mapped across immune-cell
activation, combined with Mendelian randomization (MR), are increasingly used to
nominate immune targets in cancer. We applied this framework to melanoma —
CD4⁺ T cell cis-eQTLs from eight activation profiles as exposures, a
12,530-case melanoma meta-analysis as the outcome — and then added checks
intended to strengthen the resulting target claims. Eight attempts to
substantiate a target-level claim were overturned by those checks.

FDR-significant signal is dominated by known pigmentation and naevus loci
(all 10 records in the lower-powered round; 3 of 7 independent loci under the
meta outcome, 4.1-fold enrichment, P = 0.028). SMR/HEIDI failed to reject 253 of
291 records (86.9%) that colocalisation assigned to distinct causal variants,
with only 8–20 variants entering each HEIDI test. Raising outcome power increased
MR discoveries (10 → 21) while *lowering* colocalisation support (median
PP.H3+H4 0.249 → 0.207) and replaced the candidate list entirely: two lists from
identical exposure data share no genes. Down-sampling calibrated against an
independently measured round shows that this is a general property of the design,
and that it is carried entirely by locus class: at half of observed power, 85.8%
of known-locus genes are recovered versus 22.8% of novel-locus genes, and
reaching 50% recovery requires 2,506 versus 8,771 cases. In the 3,000–8,000-case
range typical of studies using this framework, novel-locus recovery is 6–40%.
The functional layers have their own failure modes: splitting cells by a score
also splits lineage purity, and tissue-level replication cannot in principle
support a cell-type-specific claim.

One candidate, TPI1, survives the full battery under sharply bounded claims. We
present seven diagnostics, each grounded in a specific failure, that can be
applied to any study using this design.

**Keywords** ⟨Mendelian randomization; dynamic eQTL; CD4⁺ T cells; melanoma;
colocalisation; reproducibility⟩

---

## 1. Introduction

Most cis-eQTLs are context dependent. Expression quantitative trait loci mapped
in resting bulk tissue or in unstimulated primary cells cannot capture regulatory
variants that act only while a cell is responding to a stimulus, and a
substantial fraction of immune-relevant regulation falls into that category.
Dynamic eQTL datasets, in which primary immune cells are profiled across a
stimulation time course, were designed to close this gap, and their combination
with Mendelian randomization has become an attractive route to causal target
nomination: germline genotype is fixed before disease, cis-eQTL effects have a
clear directional prior, and the exposure is measured in the cell type through
which the effect is proposed to act. Several recent studies have applied exactly
this design — activation-time-course CD4⁺ T cell eQTLs as instruments, cancer
GWAS as outcomes — and have reported novel immune targets on that basis.

The inference, however, passes through more joints than the framework's summary
statistics reveal. An instrument is selected in one cell state at one timepoint
and typically consists of a single variant, so the association's P value is
supplied entirely by the outcome GWAS. Whether the eQTL and the disease signal
share a causal variant, rather than lying in linkage disequilibrium with a large
neighbouring effect, is decided by colocalisation, whose resolution depends on
outcome power. Which of two genes at a locus carries the effect is not decided by
MR at all. Neither is the cell type in which the effect operates, because
functional replication is usually available only at tissue level, where a broadly
expressed gene reports on whichever compartment dominates the tissue. Each of
these joints has a failure mode that is invisible unless it is specifically
tested, and a study that does not test for them will not observe them.

Melanoma is an unusually informative setting in which to examine these failure
modes, for two reasons that operate in opposite directions. Its common-variant
architecture is dominated by pigmentation and naevus loci, several of them with
effects far larger than anything expected from immune regulation; MC1R alone
reaches P < 10⁻⁵⁹ in current meta-analyses. That architecture provides a built-in
end-to-end positive control — a pipeline that fails to recover it is not working
— and simultaneously a built-in confounder, because long-range LD around such
loci can present as a causal association for any gene in the neighbourhood,
including genes with plausible immune functions. Melanoma is also the disease in
which CD4⁺ T cell biology has the most direct clinical relevance, through immune
checkpoint blockade, so that candidate genes can be examined against
response-stratified single-cell data rather than against annotation alone.

We applied the dynamic-eQTL MR framework to melanoma with the intention of
nominating CD4-mediated immune targets, using CD4⁺ T cell eQTLs from eight
activation profiles as exposures and a 12,530-case melanoma meta-analysis as the
outcome, followed by colocalisation, SMR/HEIDI, multi-instrument sensitivity
analysis, a negative-control phenotype, cross-cancer comparison, and functional
work in single-cell, spatial and multiome data. Every additional check was
introduced in the expectation that it would strengthen a target claim. Eight
times it did the opposite. Reporting that tally is the purpose of this paper: we
present seven diagnostics that can be applied to any study using this design,
each grounded in a specific failure we encountered, and we quantify the most
consequential of them — showing that the reproducibility of the candidate list is
set by outcome GWAS power, and that within a list, the reproducible part is
precisely the part that does not constitute a discovery. One candidate, TPI1,
survives the full battery, and we use it to show what a defensible claim from
this framework looks like and where its boundaries fall.

---

## 2. Results

### 2.1 A pipeline that passes its positive controls

We instrumented gene expression in CD4⁺ T cells using cis-eQTLs from eight
activation profiles (naive and memory cells at 0 h, 16 h, 40 h and 5 d after
anti-CD3/CD28 stimulation; n = 85–100 donors per profile) and tested their
causal association with melanoma. After harmonisation against the outcome and
restriction to genome-wide significant instruments (P < 5×10⁻⁸; minimum
F = 22.2, so the F > 10 filter was never binding), 3,556 gene × profile
exposures carried exactly one instrument, and were analysed by Wald ratio.

Two properties of this design constrain what the primary analysis can claim, and
we state them before reporting any result.

**(i) With one instrument per exposure, the MR P value is determined entirely by
the outcome GWAS.** For a Wald ratio, b = β_out/β_exp and SE = se_out/|β_exp|, so
z = β_out/se_out is independent of the exposure. The same variant therefore
returns an identical P value in every cell type and at every timepoint in which
it is the top eQTL — verified here for CDK10, where memory 5 d and 40 h give
OR 0.788 versus 0.584 with P = 2.50×10⁻³⁰ in both. Cell-type and temporal
specificity consequently cannot be argued from differences in P value; the only
admissible statement is which profile contains the instrument, and OR differences
reflect eQTL effect size rather than evidence strength. None of the three
published studies applying this framework makes this distinction.

**(ii) Steiger filtering is near-deterministic at this scale.** All 6,943
harmonised records had the correct causal direction and none were removed;
R²_exposure (median 0.152) exceeds R²_outcome (median 1.5×10⁻⁵) by four orders
of magnitude. We additionally note that TwoSampleMR's R² formula for SD units is
unbounded and returned values above 1 (maximum 1.027, 24 records > 0.8); the
bounded form R² = F/(F + N − 2) gives 0.185–0.933 and leaves the direction
unchanged. We therefore report Steiger as procedure, not as evidence; the
directional argument rests on the biology of germline cis-eQTLs.

For the outcome we combined FinnGen R12 (5,753 cases) with Rashkin et al.
(6,777 cases) by fixed-effect inverse-variance meta-analysis, giving 12,530 cases
and 789,099 controls. The meta-analysis behaves as expected: all 157 known
melanoma loci were recovered and 136 (86.6%) strengthened, and the number of
genome-wide significant variants rose from 2,871 to 4,552 (1.59×). Because
Rashkin reports only odds ratios and P values, standard errors were reconstructed
as se = |log OR| / |Φ⁻¹(P/2)|; Cochran's Q was significant for 5.33% of variants
against a 5% null expectation, which we report as evidence that the
reconstruction is calibrated.

Pipeline positive controls passed at every layer: PARP1, a known melanoma locus,
colocalised with the melanoma signal in the FinnGen round (PP.H4 = 0.92 and 0.96
in two profiles) and was confirmed by SMR/HEIDI, establishing that the LD
reference panel, rsID mapping and BESD construction all function; and the
recovery of known melanoma genetic architecture (below) is itself the strongest
end-to-end control.

### 2.2 Finding ①: the significant signal is pigmentation genetics, not immunology

In the FinnGen round, all 10 FDR-significant records fell within 1 Mb of a known
pigmentation or naevus locus from Landi et al. (VPS9D1-AS1 50 kb from MC1R;
CDK10 48 kb from the MC1R R151C variant rs1805007; PARP1 5–14 kb from PARP1).
There was no immune signal at FDR < 0.05 at all.

Higher outcome power did not change the picture qualitatively. Under the meta
outcome, 12 of 21 FDR-significant records (57.1%) lay in known-locus categories
against a 10.9% background across all instruments. Records at the same locus are
not independent, however, and counting by record inflates enrichment wherever
several genes share one locus, so we repeated the test on independent loci: 3 of
7 melanoma loci were known-locus loci, a 4.09-fold enrichment over the 10.5%
background (Fisher P = 0.028) (Fig 1, Fig 7a). We report the locus-level figure
throughout and the record-level figure in Supplementary Table ⟨S⟩.

We then asked whether this enrichment is a property of melanoma or of the
instrument panel, by applying the identical exposures to five further FinnGen
outcomes (Fig 7a). Melanoma has the largest enrichment, and the two cancers with
no expected relationship to pigmentation biology sit near background (breast
1.66×, P = 0.21; colorectal 1.59×, P = 0.49; pancreas yielded no FDR-significant
records). But the comparison does not deliver a clean specificity result:
prostate is also nominally enriched (2.38×, P = 0.031) and lung is borderline
(2.86×, P = 0.077), and no outcome survives correction across the six tests. We
therefore report this as qualified support — the ranking is as expected and the
two negative controls behave, but the design is underpowered to establish
melanoma specificity, and by record-level counting (which we do not use) lung's
apparent enrichment is driven by two genes at the single chr11 FADS locus.

### 2.3 Finding ②: SMR/HEIDI does not exclude the LD confounding that colocalisation identifies

Colocalisation assigned the MC1R-region signals — the strongest MR associations
in the study, reaching P = 4×10⁻³⁷ — to distinct causal variants, with PP.H3
dominant and PP.H4 at or near zero (CHMP1A PP.H3 = 0.99; VPS9D1-AS1 0.96;
SPATA33 0.92–1.00). This is a statistical statement about the region, and is
stronger than the positional observation that these genes lie near MC1R.

That reading is independently supported by FinnGen's own in-sample fine-mapping,
which resolves three high-purity credible sets in this region (§2.4): the region
genuinely carries multiple independent causal signals, which is precisely the
configuration in which a neighbouring gene's eQTL can be tagged by one of them
without sharing it.

SMR/HEIDI did not reproduce it. Of 291 records that colocalisation assigned to
distinct causal variants, HEIDI failed to reject homogeneity for 253 (86.9%),
including VPS9D1-AS1 (P_HEIDI = 0.649) and CDK10 (0.086, 0.081) (Fig 2). The
reason is visible in the test itself: only 8–20 variants entered each HEIDI test
at this sample size, and the melanoma GWAS signal is weak. Evidence tiers that
accept MR + SMR without colocalisation are therefore unsafe in this setting: they
would have reported the MC1R LD spillover as CD4-mediated immune targets. We use
colocalisation as the primary arbiter and SMR/HEIDI as auxiliary confirmation,
and we report the number of variants entering each HEIDI test.

### 2.4 Findings ③ and ④: outcome power moves MR and colocalisation in opposite directions, and the candidate list turns over completely

Meta-analysis raised MR discoveries from 10 to 21 records at FDR < 0.05 while
*lowering* colocalisation support: in a paired comparison of the 127 exposures
run in both rounds, median PP.H3+H4 fell from 0.249 to 0.207 and median PP.H4
from 0.124 to 0.095, with only 44.1% of exposures improving. The two methods
depend on different information — MR on one variant's z score, colocalisation on
the shape of the regional signal.

We tested and rejected three explanations for the decline. Uneven meta coverage
(Rashkin and FinnGen share 41% of variants; median within-region coverage 41.8%)
would mix two precisions within a region and distort the posterior, but
restricting colocalisation to variants present in both studies moved the median
only from 0.207 to 0.202. A candidate-composition effect is excluded by the
paired design. A winner's-curse selection effect is excluded because newly
entering candidates had *lower*, not higher, PP.H3+H4 (0.125 versus 0.202).

The actual cause is that colocalisation resolution is set by outcome power, and
the resulting errors run in both directions (Fig 3). For PARP1, increasing power
by five orders of magnitude sharpened the melanoma peak and moved it 25 kb, to
54 kb from the eQTL peak; PP.H4 fell from 0.92 to 0.05, revealing the
well-powered answer to be "distinct causal variants" and the FinnGen-era
colocalisation to be a false positive. CLDN7 (0.62 → 0.05) and ELP5
(0.55 → 0.04) behaved likewise. In the opposite direction, the meta peak for
ZFYVE19 landed exactly on the eQTL peak (0 kb) and PP.H4 rose to 0.99, revealing
a false negative.

These two observations are of different kinds, and we separate them because they
carry different model dependence.

That the candidate list **turns over** is an empirical statement about the output
of the standard pipeline. It is established by comparing lists, and it is
reproduced by the down-sampling analysis in §2.5, which recomputes only Wald
ratios and FDR and never invokes colocalisation at any point.

That the turnover operates through a **reallocation of posterior mass between H3
and H4** is a statement about the colocalisation model. `coloc.abf` assumes at
most one causal variant per trait in each region. A region containing two causal
variants — one better tagged at low power, the other emerging at high power —
could produce the same H4 → H3 movement without any change in colocalisation
resolution in the intended sense. The two claims are therefore reported
separately throughout, and we tested the second directly.

Because the exposure data are identical in both rounds, any change in the
posterior must originate on the outcome side, so the question reduces to how many
independent signals the melanoma GWAS carries in each region. FinnGen's own
fine-mapping of this endpoint, computed with in-sample LD, resolves three
high-purity credible sets in the MC1R region (log₁₀BF 45.6, 36.6 and 19.5),
confirming both that regions with multiple independent causal signals exist in
this outcome and that the single-variant assumption is violated there. Applying
SuSiE ourselves with an external European LD panel recovers ≥2 sets in the same
region — passing our pre-specified control — but returns 19 where the in-sample
analysis returns 3, and does not saturate as the permitted number of signals is
raised. The procedure therefore over-splits, and its counts can be read in only
one direction: a finding of exactly one credible set is conservative, a finding
of several is not.

Read under that constraint, the PARP1 region carries **exactly one credible set
at both power levels**, stable when the permitted number of signals is doubled,
with low LD-mismatch diagnostics (0.04 and 0.07). The collapse of PP.H4 from 0.92
to 0.05 is therefore not attributable to an unmodelled second causal signal on
the outcome side; what changed is that the peak sharpened and moved 25 kb. For
ZFYVE19 and for the TPI1/SPSB2 locus, by contrast, **no credible set is recovered
at either power** — the outcome carries no resolvable signal in those regions at
all — so neither locus can adjudicate the question, and the ZFYVE19
"false negative revealed" case must be reported with that qualification.

The wider observation is worth stating on its own. FinnGen fine-maps only 19
regions genome-wide for this endpoint, all of them classical pigmentation or
melanoma loci, and **none of our candidate loci is among them**. A candidate list
produced by this framework at this outcome power is built almost entirely on
regions that the outcome GWAS cannot fine-map.

The consequence for the candidate list is categorical. The FinnGen round yielded
four novel-locus genes passing triple validation (PRPSAP2, IMPA1, GCC2, PADI4)
and one passing the full four-fold battery (IMPA1). Under the meta outcome
**none of them survived** — IMPA1's MR P moved from 9.3×10⁻⁴ to 0.11, PRPSAP2
lost colocalisation support (PP.H3+H4 = 0.06), GCC2 and PADI4 fell below nominal
significance — and six different novel-locus genes took their place (ZFYVE19,
SMC2, KIAA0040, SPSB2, HLA-C, TPI1). **The two lists share no genes.** Both were
produced from identical exposure data by an identical pipeline; only the outcome
dataset changed.

Two earlier candidate-level failures belong here, because both were produced by
checks we added expecting confirmation. PADI4, supported by an independent breast
cancer study and by an ongoing Phase I programme, collapsed under
multi-instrument analysis (single-SNP OR 1.087, P = 1.7×10⁻³ → IVW 1.033,
P = 0.085; weighted median 1.004, P = 0.835 over seven instruments) before the
outcome was even changed. And GDI2 was removed by heterogeneity (Cochran's
Q P = 0.01), which is the sensitivity analysis working as intended.

### 2.5 The reproducible part of the candidate list is the part that is not a discovery

Finding ④ is a single anecdote — two outcome datasets, zero overlap — unless it
can be shown to be a general property of the design. We therefore measured
candidate-list recovery as a function of outcome power by down-sampling in
summary space (se_sim = se_obs × √(N_eff,obs/N_eff,sim), with a matching
inflation of the effect estimate), recomputing Wald ratios and FDR, and comparing
the resulting gene list with the full-power list (Fig 4).

The simulation is externally calibrated rather than self-validating: reducing the
meta outcome to FinnGen power must reproduce a round we had already measured
independently. It does — 10.2 simulated versus 10 real discoveries at FDR < 0.05,
and a simulated gene-list Jaccard of 0.52 [0.36, 0.73] containing the observed
value of 0.455 (Fig 4a). We twice attempted to extrapolate the curve *above*
observed power and twice failed: a flat prior on true effects overpredicted
discoveries at observed power by 25-fold, and global empirical-Bayes shrinkage
underpredicted them by 100-fold, because the true effect distribution is sparse.
Both failures were caught by the pre-specified requirement that any prior
reproduce the observed number of discoveries. **The curves are therefore reported
only up to observed power and must not be extrapolated.**

The melanoma curve shows 53% of the full-power list recovered at 5,000 cases and
85% at 10,000 (Fig 4a), and the five other cancers behave the same way: recovering
half the full-power list requires 40–80% of each study's own observed case count
(Fig 4b).

Stratifying by locus class dissolves melanoma's apparent stability (Fig 4c). At
10% of observed power, 46.0% of known pigmentation/naevus-locus genes are
recovered against **0.4%** of novel-locus genes; at 50% power, 85.8% versus
22.8%. Reaching 50% recovery requires 2,506 cases for known-locus genes and
8,771 for novel-locus genes — a 3.5-fold difference. Studies using this framework
commonly analyse outcomes with 3,000–8,000 cases, a range in which novel-locus
recovery is 6–40%.

> The reproducible part of a candidate list produced by this framework is
> precisely the part that does not constitute a discovery: known large-effect
> pigmentation loci are robust to loss of power, while the novel loci — the
> immune targets such studies set out to find — are the most fragile.

Three limits apply. The novel-locus reference list contains only four genes
(KANSL1, KIAA0040, SMC2, ZFYVE19) and the known-locus list six, so the intervals
are wide and the percentages should not be read to one decimal place. The
full-power list is not ground truth, only a higher-powered list that is itself
unstable — so the true instability is worse than the curves show. And the
composition of the novel list is itself questionable: KANSL1 lies in the
chr17q21.31 MAPT inversion, a known false-positive-prone region, and KIAA0040 was
later downgraded.

### 2.5b A natural experiment: five sequential releases of the same GWAS resource

The comparison in §2.4 changes two things at once — outcome power, and the
dataset itself. To separate them we used FinnGen's sequential releases as a
quasi-longitudinal experiment, holding the exposure data, instrument selection,
harmonisation, MR estimator, multiple-testing family and locus definition fixed
and changing only the release. Five releases carry the identical melanoma
endpoint (R8 2,705 cases; R9 2,993; R10 3,194; R11 3,932; R12 5,753), covering
47–100% of the reference effective sample size. Releases are nested, so
agreement between them is expected to exceed that between independent studies of
the same size; this makes the design conservative with respect to instability.
An earlier release (R6) was excluded after inspection because the endpoint
carries only 143 cases there and is not comparable.

Predictions were generated by down-sampling R12 and were **registered before any
earlier release was examined** ⟨Supplementary; pre-registration document⟩. All
twelve comparisons — discovery counts, overall recovery, and the known-minus-novel
recovery differential at each release — fell inside the registered intervals
(Fig ⟨G⟩). Discoveries rose with case number (Spearman ρ = 0.894), the analysis
reproduced the published R12 result exactly (10 records, 6 genes), and the large
known-locus effects were present in the earliest release.

Two results follow. First, the known-locus part of the list is recovered early
and then does not move: the same five genes constitute the FDR-significant list
at R9, R10 and R11 without change, and four of them are already present at R8, at
47% of reference power. Second, and more consequentially, **no novel-locus gene
reaches FDR < 0.05 in any release** — not at 2,705 cases, and not at 5,753. The
locus-attribution result of §2.2 therefore holds at five independent power levels
rather than in a single dataset.

We are explicit about what this design cannot show. Because no novel-locus
candidate exists anywhere in this power range, the trajectory cannot demonstrate
turnover *among* novel candidates within a single resource; it can only show that
they do not appear. The turnover documented in §2.4 therefore still involves a
step between two datasets, and the contribution of cohort heterogeneity to it is
narrowed but not eliminated. Two further details are worth recording because both
were anticipated in the registration: the apparent rise in novel-locus recovery at
R11 is produced entirely by three genes at a single locus (CRHR1, KANSL1,
KANSL1-AS1, spanning 573 kb of the chr17q21.31 inversion) and disappears when
independent loci rather than genes are counted; and one novel-locus candidate
(MT2A) enters the relaxed list at R11 and is absent at R12 — a candidate visible
only at lower power.

### 2.6 What can locate an effect once MR cannot: orthogonal functional evidence

MR at this resolution cannot say which gene at a locus carries the effect, nor in
which cell type the effect operates. Both questions can be answered, but only by
data of a different kind — and the answer is not always positive.

**Gene attribution.** SPSB2 and TPI1 lie 3–20 kb apart on chr12p13, with the TPI1
instrument bracketed by two SPSB2 instruments; the six candidate genes therefore
correspond to five independent loci. Three orthogonal lines assign the effect to
TPI1: single-cell ICB response stratification is significant for TPI1
(P = 0.0076 pre-treatment, 0.0057 post-treatment) and null for SPSB2 (0.17); TPI1
is detected in 60.6% of CD4⁺ T cells versus 5.0% for SPSB2 in an independent
cohort, so the SPSB2 nulls are not simply dropout in the discovery data; and the
two genes have entirely different eQTL dynamics — SPSB2 is constitutive
(1.4×10⁻¹¹ at 0 h) whereas TPI1 is a transient activation-window pulse
(§2.9). The three lines come from different data types, and this generalises
as a method: use functional data to resolve gene attribution that MR cannot.

**Compartment attribution — and here the answer is negative.** For a broadly
expressed gene, tissue-level replication cannot in principle support a
cell-type-specific causal claim. We tested this directly with spatial
transcriptomics (four patients, eight sections, 2,317 spots) (Fig 5). Tissue-level
TPI1 tracks the glycolytic module (meta ρ = +0.174, 95% CI 0.134–0.214,
P = 3.6×10⁻¹⁷, positive in 8/8 sections) and the tumour/melanocyte compartment
(+0.137), and is *negatively* correlated with the lymphoid compartment (−0.080).
In a per-section joint regression, β_tumour = +0.100 (positive in 8/8 sections)
versus β_lymphoid = −0.024; by compartment assignment, tumour spots exceed
lymphoid spots by +0.182 (P = 2.5×10⁻⁷). After regressing out the glycolytic
module, the residual association with the tumour compartment vanishes (ρ = +0.026,
P = 0.21) and no hidden lymphoid signal emerges (ρ = −0.055).

The negative is interpretable only because the same sections contain a positive
control with the mirror-image behaviour: HLA-C tracks the lymphoid (+0.187, 8/8)
and myeloid (+0.231, 8/8) compartments and not the tumour (−0.043, 1/8), with
tumour < lymphoid at P = 1.9×10⁻⁶; SMC2 tracks proliferation (+0.105, 8/8). The
assay can detect an immune-compartment gene at this power; TPI1 is not one at
tissue level. We also record a test that failed its own control and is therefore
uninformative: within the 450 lymphoid-dominant spots, neither TPI1 (P = 0.13)
nor the HLA-C positive control (P = 0.23) showed the expected relationship, so
that subset is underpowered and its null is not evidence.

The consequence is that the two bulk ICB cohorts in which TPI1 tracked response
(GSE78220 P = 0.045; GSE91061 on-treatment P = 0.030, but pre-treatment P = 0.81)
cannot be reported as replication of a CD4-specific effect. Tissue-level TPI1 is
tumour glycolysis, and tumour glycolysis causing ICB resistance is established
biology. The CD4-specific claim rests on the construction of the instruments and
on a single CD4-resolved cohort.

### 2.7 Finding ⑦: the functional layer has its own failure modes, and two of ours were severe

**A score-based split also splits cell-type purity.** Dividing CD4⁺ T cells by
glycolysis score produced an apparent phenotype — reduced cytotoxic and increased
helper/regulatory signatures (mean Δ −0.318 and +0.235, 24/28 samples, FDR
1.1×10⁻⁴). The two module scores contain lineage markers (CD8A/CD8B and CD4
themselves), and the split was strongly confounded by purity: CD4 was detected in
40.1% of glycolysis-high versus 19.5% of glycolysis-low cells, CD8A in 47.6%
versus 63.9%. Removing lineage genes from the modules, matching on a CD8ness
score in addition to sequencing depth (610 cells, all 28 samples retained,
residual SMD 0.0097 for CD8ness and 0.0000 for depth) abolished the phenotype
entirely: the lineage-free effector score gave Δ = +0.013 (P = 0.88) with 14/28
sample concordance — exactly chance — and FOXP3, CTLA4 and GZMB returned
P = 0.73, 0.80 and 1.00. The effect sizes collapsed by 75–100%, so this is not a
power loss. Critically, a cluster-level control had passed: contaminating cells
are distributed within clusters, so the control must be performed at the cell
level.

**Two data-processing errors changed every patient-level number.** First, a
Seurat object used for the single-cell work had been re-normalised — log₂(TPM+1)
values passed into `NormalizeData()` — giving a Spearman correlation of only
0.874 with correctly scaled values; cell-level rank tests were unaffected but
patient-level means were not. Second, the response annotation table had merged
samples with conflicting outcomes: patients with two post-treatment biopsies can
have different response labels (Post_P1 responder, Post_P1_2 non-responder), and
these had been collapsed. Both were corrected and all patient-level analyses
recomputed from the GEO annotation (final samples: 9 R / 10 NR pre-treatment,
8 R / 20 NR post-treatment).

The corrected numbers strengthened the main claims (TPI1 pre-treatment P = 0.0076,
FDR 0.046; post-treatment P = 0.0057, FDR 0.017; glycolysis module post-treatment
P = 0.0017) but eliminated one mechanism we had proposed: the strong negative
correlation between HLA-C and proliferation (ρ = −0.612/−0.720) disappeared
(ρ = +0.142, P = 0.56; −0.120, P = 0.54). That mechanism had already been
withdrawn when it failed to replicate in an independent cohort — replication
caught the artefact before we knew its cause. We note also that the same test
applied to the second cohort would have been wrong: its raw file is *not*
normalised (CV = 1.15, 46-fold range), so library-size normalisation there is
correct. The criterion is whether the source file is already normalised, not
whether a derived object matches it.

### 2.8 Finding ⑥: MR names the gene that has an instrument

The functional data describe a pathway-level phenotype. Across 22 glycolytic
enzymes measurable in CD4⁺ T cells, 20 are higher in non-responders
post-treatment (binomial P = 6.1×10⁻⁵; 8 surviving FDR < 0.05) and 18 of 22
pre-treatment (P = 2.2×10⁻³; 4 surviving FDR), with PGAM1 the strongest enzyme at
baseline and TPI1 fourth at both timepoints. An earlier module-score analysis had
returned a misleading null at baseline; gene-by-gene testing showed the module
had been diluted by enzymes detected in 1–9% of cells.

Genetically, only 2 of 28 glycolytic genes carry a genome-wide significant
instrument in any profile — TPI1 (OR 1.31, P = 1.4×10⁻³) and ENO1, which is null
(OR 1.02, P = 0.29). PGAM1, GAPDH, PKM, LDHA and PFKP have none. A
pathway-level enrichment test that does not require instruments (top cis-eQTL per
gene, |z| in the melanoma GWAS, compared against eQTL-P- and MAF-matched
background over 5,000 samplings) was also null: mean |z| 0.890 versus
0.937 ± 0.267, empirical P = 0.52, but with only nine gene × profile
combinations it is underpowered and we do not claim the pathway is signal-free.

**The ratio is 1 in 22.** Functionally the phenotype is pathway-wide and the
strongest baseline enzyme is PGAM1; genetically MR can see one gene, and names it
TPI1. This is not a defect of MR, but left undeclared it invites the reading that
the named gene is the pathway's most important member. A related constraint
operates at the resource level: TPI1 is absent from all eight stimulated T-cell
datasets we identified in the eQTL Catalogue despite TPM of 98–194 where it is
quantified, so gene-level filtering in public resources determines which findings
can be tested for replication at all.

One class-level genetic signal did survive, and it was specified in advance by the
TPI1 observation: glycolytic genes' strongest eQTLs avoid the resting state
(10.4% of 48 gene × lineage combinations at 0 h versus 25.3% of 24,077 background
genes; OR = 0.343, P = 0.019). Resting-state eQTL studies systematically miss
genetic regulation of this pathway.

### 2.9 TPI1 as a bounded worked example

TPI1 is the only candidate to pass MR, colocalisation, SMR/HEIDI and
multi-instrument sensitivity analysis, and its evidence structure is the inverse
of the usual one — marginal genetics, strong function. We report it as an
application of the diagnostics above, not as a target.

**Genetics.** MR OR = 1.305, P = 1.4×10⁻³ (FDR = 0.119, i.e. not
FDR-significant); colocalisation PP.H4 = 0.51, at the edge of our criterion
(PP.H4/(PP.H3+PP.H4) > 0.7 and PP.H3+PP.H4 > 0.5); SMR P = 3.0×10⁻³ with HEIDI
P = 0.45; and over five independent instruments IVW P = 3.8×10⁻⁴, weighted median
P = 2.6×10⁻³, Cochran's Q P = 0.83. It entered the list through the conjunction,
not through any single threshold — and the conjunction is what carries it, because
each component is individually marginal. To that we add a further qualification
from the sensitivity analysis in §2.4: SuSiE applied to the melanoma outcome
recovers **no credible set in the TPI1/SPSB2 region at either power level**, so
the colocalisation posterior of 0.51 is computed on a region in which the outcome
GWAS carries no resolvable signal. The genetic side of this example is weaker than
the four passed filters make it appear, and we report it as such.

**A transient activation window (Fig 6a).** The TPI1 cis-eQTL is a 16 h pulse in
both lineages: naive 8.4×10⁻⁵ (0 h) → 1.7×10⁻¹² (16 h, 148 variants at
P < 5×10⁻⁸) → 8.7×10⁻³ (40 h) → 3.5×10⁻³ (5 d); memory 6.2×10⁻³ → 1.2×10⁻⁹
(53 variants) → 3.7×10⁻² → 7.1×10⁻³. This is not an artefact of differential
detection: the fraction of tested genes carrying an instrument is flat across
profiles (51.3–55.0%), and is in fact slightly lower at 16 h than at 40 h. The
co-located SPSB2 shows no such pattern. Two lineages replicate it independently,
and the class-level result above generalises it. We state the limit of this
observation precisely: the absence of instruments at 0 h, 40 h and 5 d means the
eQTL did not reach P < 5×10⁻⁸, not that no causal effect exists at those
timepoints. The admissible claim is that dynamic data supply instruments that
resting-state data do not.

Because this is the study's most original observation, we verified that it is not
an alignment artefact after finding TPI1 absent from public stimulated-T-cell
resources. For SPSB2, quantified in both, our pipeline reproduces the
standardised processing of the same experiment (four variants, β within 10%,
P within an order of magnitude, direction 4/4). The TPI1 signal itself has clean
cis architecture: the lead variant is 3,816 bp upstream of the TSS, 148 variants
exceed 5×10⁻⁸ across 102 kb, and −log₁₀P decays from 11.8 within ±2 kb to 1.0
beyond 200 kb. It is locus-specific (only TPI1 and SPSB2 associate strongly; the
third gene drops to P = 10⁻³), and all TPI1 pseudogenes lie on other chromosomes,
where mismapping would produce trans rather than a cis peak at the TPI1 promoter.

**Function (Fig 6b–c).** In CD4⁺ T cells from ICB-treated melanoma patients, TPI1
is higher in non-responders both before treatment (FDR = 0.046) and after
(FDR = 0.017); the result is unchanged by residualising on lineage score
(P identical to three decimal places) and survives restriction to CD8-negative
cells post-treatment (FDR = 0.005). It sits inside the coordinated pathway-wide
shift described above rather than standing alone.

**Melanoma specificity (Fig 7b).** Across five further cancers the TPI1 effect is
absent (P = 0.086–0.69), with point estimates in the opposite direction for
breast (0.915) and prostate (0.878). The contrast within our own candidate list is
instructive: SMC2 is significant and concordant in four cancers (melanoma 1.415,
lung 1.384, breast 1.304, prostate 1.232), which identifies it as a general
proliferative mechanism rather than melanoma-specific immune surveillance and
retires our earlier "cell-division module" interpretation. HLA-C is significant in
the opposite direction in three other cancers, which we treat as a reason for
caution rather than as mechanism, given the extreme pleiotropy of the MHC region.
TPI1 also affects melanoma without affecting melanocytic naevi (OR 0.867,
P = 0.079) in a test whose sensitivity is established by five known pigmentation
loci that do affect naevi.

**An in vitro axis that is not simply activation.** In purified CD4⁺ T cells
(10x multiome, 40,495 nuclei), the glycolytic programme is close to orthogonal to
activation intensity at the RNA level (R² = 0.024 for glycolysis ~ activation +
depth, against a pre-specified R² > 0.70 cut-off for declaring confounding), and
is not a nuclear/cytoplasmic composition artefact (nuclear-retention index removes
0.047 of variance). Ranking all measurable genes along the axis without
pre-specified modules places ribosomal proteins (12.0% of the top 300, 10.3-fold)
and OXPHOS at one end and T cell identity and quiescence factors (TOX, IKZF2,
MYB, TXK, CAMK4, SESN3) at the other — an anabolic-growth versus
T-cell-identity axis. At the chromatin level the axis shares 42% of its variance
with activation (r = 0.647), so we restricted the analysis to activation-invariant
peaks (|activation log₂FC| ≤ median, retaining 68% of the axis amplitude): the
glycolysis-high end retains all 23 AP-1 motifs (median odds 1.81 → 1.52;
BATF::JUN 2.22 → 1.62, P = 9.1×10⁻¹⁰), while the low end is enriched for IRF
family and STAT1::STAT2 interferon-response elements. Motif similarity within
families means these results identify families, not members, and restriction to
activation-invariant peaks is mitigation rather than separation. We also report a
tension rather than resolving it: in this 15 h system the glycolysis-high end has
*low* TOX, whereas in patient tumour-infiltrating CD4 cells TPI1 correlates
positively with an exhaustion score (ρ = 0.617–0.726 in two cohorts). Early
activation TOX and exhaustion-associated TOX are not necessarily comparable, but
the discrepancy stands.

**What we do not claim.** The MR establishes a causal relationship with melanoma
**risk**; the association with ICB response is observational, from a single
CD4-resolved cohort, and we do not join the two into a causal chain. TPI1 is
DepMap-essential with pLI = 0.87, so no target claim is made. The activation
window has no independent replication: TPI1 is not quantified in any of the eight
stimulated T-cell eQTL datasets available, and the one resource matching our
design is the same experiment as our exposure data, so using it would be
circular. And by finding ④ and Fig 4, TPI1 itself may not survive a
higher-powered outcome GWAS.

### 2.10 No candidate is strong on all axes (Fig 8)

Assembling the layers for the six novel-locus candidates shows that statistical
strength, pathway attribution, and druggability/safety point at different genes.
ZFYVE19 has the strongest colocalisation evidence of any candidate (PP.H4 = 0.99,
GWAS and eQTL peaks coincident, zero catalogued pleiotropy, non-essential,
LoF-tolerant) but affects naevi (OR 0.962, P = 5.3×10⁻³), shows no relationship to
proliferation or ICB response, and is pan-cancer in direction. SMC2 has the
largest effect and the cleanest pathway attribution but is DepMap-essential with
pLI = 0.99999 and is pan-cancer. HLA-C carries the most instruments and an
intriguing mirror in vitiligo (OR 0.681, but from only 391 cases) yet lies in the
most pleiotropic region of the genome and carries SJS/TEN safety flags. KIAA0040
is uninterpretable and gave opposite directions in two bulk cohorts; we removed
it. TPI1 is the only candidate with converging functional evidence, and it is an
essential gene with marginal genetics.

We report this matrix rather than promoting a single target, because the
alternative would require exactly the kind of selective emphasis that the seven
diagnostics above are designed to detect.

---

## 3. Discussion

### 3.1 What the study set out to do and what it became

Dynamic eQTL data from activated CD4⁺ T cells, combined with Mendelian
randomization, has been used to nominate immune targets across several cancers.
We applied that framework to melanoma with the intention of nominating targets,
and instead found that the framework's output is far less stable than its use in
the literature implies. Every check we added in order to *strengthen* a target
claim returned a negative, and it is the accumulation of those negatives —
rather than any single candidate gene — that constitutes our main result.

We report seven methodological findings, each of which is a diagnostic that can
be applied to any dynamic-eQTL-to-disease MR study, together with a worked
example (TPI1) that survives the full battery but only under sharply bounded
claims.

### 3.2 The central argument: the failure rate is the finding

Across the study we made eight attempts to substantiate a target-level claim.
All eight were overturned by checks we ourselves introduced. Critically, these
were not adversarial tests designed to find fault: each was added in the
expectation that it would *support* the claim.

| # | Attempt | Check applied | Outcome | Finding |
|---|---|---|---|---|
| 1 | PADI4 as a novel immune candidate | multi-instrument sensitivity analysis | OR 1.087 (P=1.7×10⁻³) collapsed to IVW 1.033 (P=0.085), weighted median 1.004 (P=0.835) over 7 instruments | single-instrument Wald ratios are not self-validating |
| 2 | FinnGen-round candidate list (PRPSAP2, IMPA1, GCC2, PADI4) | repeat with a higher-powered meta outcome | none survived; **zero overlap** between rounds | ④ |
| 3 | HLA-C acting through CD4 proliferation | replication in an independent scRNA cohort | direction reversed; later traced to a normalisation error (ρ = −0.612/−0.720 → +0.142/−0.120, P ≈ 0.55) | replication is what caught it; the mechanism was an artefact |
| 4 | Bulk ICB cohorts as evidence of a CD4-specific effect | spatial transcriptomics with a cell-type-restricted positive control | tissue-level TPI1 tracks the glycolytic module (ρ = +0.174, 8/8 sections) and is *negatively* correlated with the lymphoid compartment (ρ = −0.080) | ⑤ |
| 5 | A transcriptional state for glycolysis-high CD4 cells | matching on lineage score as well as depth | after matching (610 cells, 28 samples, residual SMD 0.010/0.000) every module went null; the lineage-free effector score gave 14/28 sample concordance, exactly chance | ⑦ |
| 6 | Pathway-level MR of glycolysis | count instruments across the pathway | only 2 of 28 glycolytic genes carry an instrument; the second (ENO1) is null (P=0.29) | ⑥ |
| 7 | Pathway-level enrichment of melanoma GWAS signal at glycolytic eQTLs | expression- and MAF-matched permutation | mean \|z\| 0.890 vs 0.937 in matched background, P=0.52 (n=9, underpowered) | ⑥ |
| 8 | Cross-cancer comparison as proof of melanoma-specific locus attribution | count independent loci rather than gene records | melanoma 4.09× (P=0.028) but prostate 2.38× (P=0.031); neither survives correction across six outcomes | non-independence inflates enrichment estimates |

Two further attempts failed for reasons that are *not* scientific negatives and
should not be counted alongside the eight:

- **A chromatin test of the activation-window eQTL** returned null in the
  pre-specified window, but the test was mismatched to the claim: an
  activation-specific eQTL is a genotype × timepoint interaction, whereas mean
  accessibility is neither a necessary nor a sufficient correlate of one. The
  matched test is allele-specific accessibility, which requires read-level data
  and genotypes that are not publicly available. The null is therefore
  uninformative, and we report it as such.
- **Independent replication in public eQTL resources was not possible.** TPI1 is
  absent from all eight stimulated T-cell datasets we could identify in the eQTL
  Catalogue, despite median TPM of 98–194 where it is quantified; the one
  resource with matching design and timepoints is the same experiment as our
  exposure data, so using it would be circular.

A third category belongs in this discussion for the opposite reason. Before
testing whether peripheral-blood CD4 glycolysis could serve as a correlate of the
tumour phenotype, we tested the precondition — whether it is an individual-level
trait at all. Within a clinically homogeneous patient group it was
indistinguishable from random gene modules (ICC 0.0285 against a matched random
floor with a 95th percentile of 0.0403), so the downstream analysis was cancelled
and archived unexecuted. That is not a ninth failure; it is what a stopping rule
looks like when it works, and it cost about fifteen minutes to avoid an analysis
whose positive and negative results would both have been uninterpretable.

We draw attention to this tally not as a confession but as data. A pipeline that
passes all of its positive controls, and that recovers known melanoma biology,
nevertheless failed to sustain a single target-level claim across eight
independent attempts. Studies using the same framework that do not report such
attempts have not necessarily avoided these failures; they have not tested for
them.

### 3.3 Seven diagnostics, arranged as a progression

The findings are not a checklist. They form a sequence in which each answer
raises the next question.

**① Where does the signal come from?** All 10 FDR-significant associations in the
first round fell within 1 Mb of a known pigmentation or naevus locus (MC1R
rs1805007, 48 kb; PARP1, 5–14 kb); there was no immune signal at FDR<0.05 at all.
Under the higher-powered meta outcome the pattern persisted: counting independent
loci — records at one locus are not independent — 3 of 7 melanoma loci were
known-locus loci, a 4.1-fold enrichment over the 10.5% background (P=0.028;
by record, 12/21 = 57%, 5.3-fold, P=2.7×10⁻⁷, reported in Supplementary). The
same instrument panel applied to five other cancers gives melanoma the largest
enrichment, but prostate is also nominally enriched (2.4-fold, P=0.031) and
neither survives correction across six outcomes, so the cross-cancer comparison
qualifies this finding rather than proving specificity (see attempt 8 above).
None of the three comparable published studies reports this check in any form.

**② Can SMR/HEIDI exclude LD confounding?** Not at this scale. Of 291 records
that colocalisation assigned to distinct causal variants, HEIDI passed 253
(86.9%), with only 8–20 variants entering each HEIDI test. Evidence tiers that
treat MR + SMR without colocalisation as sufficient are therefore unsafe here.

**③ Does more outcome power resolve the problem?** It moves the two methods in
opposite directions. Meta-analysis raised MR discoveries (FDR<0.05: 10 → 21)
while lowering colocalisation support (median PP.H3+H4 0.249 → 0.207 in a paired
comparison), because MR depends on one variant's z score whereas colocalisation
depends on the shape of the regional signal.

**④ How stable is the candidate list?** It is not, and the answer has two parts
that should not be conflated.

*④a — the list turns over (model-free).* Two candidate lists produced from
identical exposure data under two outcome datasets share no genes, and the
down-sampling analysis reproduces this behaviour as a smooth function of outcome
power using MR alone, with colocalisation playing no part in it. This is a
statement about what the standard pipeline outputs, and it does not depend on any
colocalisation model. It is the study's most consequential finding, and the
reason we do not present any single gene as an established target.

*④b — the mechanism inside colocalisation (model-dependent).* The turnover is
accompanied by movement of posterior mass between H3 and H4, in both directions:
PARP1 went from PP.H4 = 0.92 to 0.05 once the GWAS peak sharpened and moved 25 kb
(an apparent false positive removed), while ZFYVE19 rose to 0.99 as the meta peak
landed exactly on the eQTL peak (an apparent false negative revealed). Reading
this as a change in colocalisation *resolution* assumes the model under which it
was computed: `coloc.abf` permits at most one causal variant per trait per
region, so a region with two causal signals — differently tagged at different
power — could generate the same posterior movement without any resolution change
in the intended sense. We tested this. Because the exposure data are identical in
both rounds, the question reduces to the number of independent outcome-side
signals per region. In a procedure calibrated against FinnGen's own in-sample
fine-mapping — and shown by that calibration to *over-split* signals, so that a
count of one is conservative — PARP1 carries exactly one credible set at both
power levels. The posterior collapse there is therefore not an artefact of the
single-variant assumption. For ZFYVE19 and for the TPI1 locus no credible set is
recovered at either power, so those regions cannot adjudicate the question, and we
report the ZFYVE19 case with that qualification. The multiple-signal explanation
is excluded for the case where it can be tested and remains open elsewhere; we do
not claim more. We keep ④a and ④b separate because ④a survives whatever the
answer to ④b turns out to be.

A related fact deserves its own line. FinnGen fine-maps 19 regions genome-wide
for this endpoint, all classical pigmentation or melanoma loci, and none of our
candidate loci is among them. At this outcome power the candidate list is built
almost entirely on regions the outcome GWAS cannot fine-map.

**⑤ What can locate an effect once MR cannot?** Orthogonal functional data, at
two levels. *Gene attribution*: SPSB2 and TPI1 lie 3–20 kb apart and MR cannot
separate them, but three independent lines — single-cell response
stratification, detection rate, and eQTL temporal profile — all point to TPI1.
*Compartment attribution*: for a broadly expressed gene, tissue-level
replication cannot in principle support a cell-type-specific causal claim.
Spatial data adjudicated this and returned a negative for TPI1; the mirror-image
behaviour of HLA-C in the same sections (lymphoid ρ = +0.187, 8/8; tumour −0.043,
1/8) shows the assay had the resolution to detect an immune-compartment gene,
which is what makes the TPI1 negative interpretable.

**⑥ Which gene does MR actually name?** The one with an instrument. Functionally,
20 of 22 glycolytic enzymes are higher in non-responders' CD4 T cells at the
post-treatment timepoint (binomial P=6.1×10⁻⁵; 8 surviving FDR), and the
strongest at baseline is PGAM1. Genetically, only TPI1 and ENO1 carry an
instrument, and ENO1 is null. The ratio is 1 in 22. This asymmetry is not a
defect of MR, but it is routinely left undeclared, and readers will otherwise
read "the gene MR named" as "the most important gene in the pathway". A related
constraint operates at the resource level: gene-level filtering in public eQTL
databases determines which findings can even be tested for replication.

**⑦ Do the functional checks have their own failure modes?** Yes, and one is
severe. Splitting a cell population by a gene or module score also splits it by
cell-type purity. In our data the glycolysis-high and -low groups differed
markedly in composition (CD4 detected in 40.1% vs 19.5% of cells; CD8A in 47.6%
vs 63.9%). Matching on sequencing depth does not fix this. Only after matching
on a lineage score as well did the apparent phenotype disappear. Crucially, a
cluster-level control passed while the cell-level control failed — contaminating
cells are distributed across clusters, so the control must be performed at the
cell level.

### 3.4 What survived, and why

The findings that survived every round of correction share a property: none
depends on a single contrast in a single dataset under a single processing
choice.

- TPI1's association with ICB response survived correction of a normalisation
  error, redefinition of samples, lineage-composition control, and restriction
  to CD8-negative cells, and sits within a pathway-wide pattern.
- The activation-window eQTL replicated across two independent lineages, was
  controlled against global instrument yield (flat at 51.3–55.0% across
  timepoints), and has a class-level counterpart.
- Locus attribution is a descriptive fact that does not depend on any model.

Conversely, every claim we withdrew rested on one comparison in one dataset. We
suggest this as a practical prior: in data of this scale, a conclusion supported
by a single contrast has a low probability of surviving.

### 3.5 TPI1: a bounded example

TPI1 is the only candidate to pass the full battery, and its evidence structure
is the inverse of the usual pattern in MR studies — the genetics are marginal
and the functional data are strong.

*Genetics.* MR P=1.4×10⁻³ but FDR=0.119; colocalisation PP.H4 = 0.51, at the
edge of our criterion. It entered the final list through the conjunction of MR,
colocalisation, SMR/HEIDI and multi-instrument sensitivity analysis, not through
any single threshold. Its cis-eQTL is confined to a narrow activation window
(naive 8.4×10⁻⁵ at 0 h → 1.7×10⁻¹² at 16 h → 8.7×10⁻³ at 40 h; memory
6.2×10⁻³ → 1.2×10⁻⁹ → 3.7×10⁻²), a pattern not attributable to differential
instrument yield, and one that the co-located and constitutively regulated SPSB2
does not share. Glycolytic genes as a class avoid the resting state (10.4% vs
25.3% of genes with their strongest eQTL at 0 h; OR = 0.34, P = 0.019).

*Function.* TPI1 in CD4 T cells is higher in ICB non-responders at both
timepoints (FDR 0.046 pre-treatment, 0.017 post-treatment), robust to lineage
control (CD8-negative subset, post-treatment FDR 0.005), and embedded in a
coordinated pathway-wide shift. Across five other cancers the MR effect is absent
(P = 0.086–0.69), whereas SMC2 is significant in four, indicating that TPI1's
effect is melanoma-specific while SMC2's reflects a general proliferative
mechanism.

*What cannot be claimed.* The MR establishes a causal relationship with melanoma
**risk**; the association with ICB response is observational and comes from a
single CD4-resolved cohort. These are separate claims and we do not join them
into a causal chain. TPI1 is a DepMap-essential gene (pLI 0.87), so no target
claim is made. And by our own finding ④, the current list — including TPI1 — may
not survive a larger outcome GWAS.

### 3.6 Limitations

Beyond those already stated: instruments were single cis-eQTLs, so Steiger
filtering was near-deterministic and is reported as procedure rather than
evidence; relaxed instrument sets at r² < 0.1 retain correlation that inflates
IVW significance (157 FDR-significant exposures versus 21 in the strict
single-instrument set); Rashkin standard errors were reconstructed from odds
ratios and P values, with Cochran's Q significant at 5.33% supporting the
calibration; the single-cell cohorts are small (9 vs 10 donors pre-treatment);
the exposure eQTLs come from 85–100 donors per profile, so top-variant effect
sizes are subject to winner's curse; the CD4 population used for functional work
has limited purity (CD4 detected in 32.5% of cells), which we address by
reporting the CD8-negative sensitivity analysis rather than by redefining the
population; the spatial data comprise four patients at 100 µm resolution, so a
CD4-restricted effect could be diluted below detection — the admissible
conclusion is that tissue-level TPI1 is dominated by tumour glycolysis, not that
TPI1 has no role in CD4 T cells; and the power-stability curves rest on reference
lists of six and four genes, are not extrapolated above observed power, and
compare against a full-power list that is itself unstable, so they understate
rather than overstate instability. Two alternative explanations for the
colocalisation results were addressed directly rather than assumed away, and
neither is fully closed: the multiple-signal explanation is excluded for PARP1,
where it can be tested, but ZFYVE19 and the TPI1 locus carry no resolvable outcome
signal at either power and cannot adjudicate it; and the sequential-release
experiment removes cohort heterogeneity as an explanation for the stability of the
known-locus part of the list, but cannot address the novel part, which is empty
throughout that power range. Finally, the colocalisation supporting our worked
example sits in a region the outcome GWAS cannot fine-map — a limitation we state
rather than absorb.

### 3.7 What should change in practice

For studies applying this framework, we suggest five concrete changes, each
following directly from a finding above: annotate significant signals against
known loci for the trait before interpreting them, counting independent loci
rather than gene records (①); treat colocalisation as the primary arbiter and
SMR/HEIDI as auxiliary, reporting the number of variants entering each HEIDI test
(②); report the physical distance between GWAS and eQTL peaks alongside posterior
probabilities, and test whether posterior movement survives explicit modelling of
multiple causal signals before interpreting it as a change in resolution (④); state explicitly what fraction of a pathway carries
instruments before describing a named gene as the pathway's driver (⑥); and, when
splitting cells by a score, match on lineage composition as well as depth, and
perform the control at the cell level (⑦).

More generally: a candidate list produced by this framework at current outcome
GWAS power should be reported as *the set of genes that passed screening under
these conditions*, not as a set of targets. Our two lists, produced from
identical exposure data under two outcome datasets, share nothing. And because
the fragile part of such a list is precisely its novel-locus part, the headline
finding of a study of this kind is the part least likely to replicate.

---

## 4. Methods

*(完整 Methods 见 `manuscript/METHODS_draft.md`，共 20 节；此处为投稿版正文，
内容与该文件一致。为避免两处正文分叉，定稿前请只改 METHODS_draft.md，
再重新组装本文件。)*

**§4.1 Exposure data** — dynamic CD4⁺ T cell cis-eQTLs, eight profiles, GRCh38 as
distributed (no liftover), per-profile N recovered from `ma_count/(2·MAF)`.
**§4.2 Outcome GWAS and meta-analysis** — FinnGen R12 + Rashkin, SE
reconstruction, df = 1 Cochran's Q, P underflow handling.
**§4.3 Instrument selection and harmonisation** — strict (P < 5×10⁻⁸) and relaxed
(F ≥ 5, clumped at r² < 0.1) sets.
**§4.4 Mendelian randomization** — Wald ratio and its analytic consequence, IVW /
IVW-MRE / weighted median, weighted mode excluded and why, Cochran's Q, BH-FDR
within pre-defined families, Steiger with bounded R².
**§4.5 Colocalisation** — `coloc.abf`, the added PP.H3+PP.H4 > 0.5 requirement,
peak-distance and variant-count diagnostics, `n_studies == 2` re-run.
**§4.6 SMR and HEIDI** — v1.3.1, 525-donor GRCh38 EUR panel, BESD construction,
variant counts reported per test.
**§4.7 Locus annotation and enrichment** — Landi reference set (157 loci), 1 Mb
window, **enrichment counted by independent locus**.
**§4.8 Negative-control phenotype and comorbidity MR** — melanocytic naevi
(13,357 cases) as the pigmentation-pathway control; vitiligo case count stated.
**§4.9 Cross-cancer analysis** — exposure side held fixed, five FinnGen outcomes.
**§4.10 Power–stability simulation** — down-sampling in summary space, 200
replicates, pre-specified external calibration, extrapolation attempted and
abandoned.
**§4.11 Single-cell ICB cohorts** — GSE120575 and GSE115978, the normalisation
criterion, sample and response-label definitions, detection rates reported.
**§4.12 Module scores, response testing and lineage control** — gene sets, sample
as statistical unit, gene-by-gene pathway testing, cell-level CD8ness matching.
**§4.13 Spatial transcriptomics** — Thrane sections, QC, within-section z-scored
compartment scores, Fisher-z meta with per-section concordance, HLA-C/SMC2
positive controls, the discarded within-lymphoid test.
**§4.14 Purified CD4⁺ T cell multiome (GSE282266)** — snRNA caveats, axis
definition with pre-specified kill criteria, nuclear-retention control, module
detection checked before interpretation.
**§4.15 Chromatin accessibility and motif enrichment** — 199,740-interval
consensus peak set, JASPAR2020 CORE with GC- and accessibility-matched
background, activation-invariant restriction.
**§4.16 Peripheral blood variance decomposition (GSE199994)** — ICC with
upper-bound and random-module floor controls; precondition not met, downstream
analysis archived unexecuted.
**§4.17 Public eQTL resource survey** — rsID-based querying, HTTP 400 ≠ null,
circularity check against the exposure experiment.
**§4.18 Pre-specification, positive controls and stopping rules** — six rules
fixed before the analyses they govern, with the tests they discarded.
**§4.19 Software and computing environment** · **§4.20 Data and code
availability**.

---

## 5. Figure legends

**Fig 1 | Locus attribution of CD4⁺ T cell eQTL–melanoma MR signals.**
Manhattan plot of MR P values for all 3,556 strict-set exposures under the meta
outcome, coloured by locus category (novel, naevus count, pigmentation/hair
colour, known melanoma). Triangles mark associations whose P underflowed
floating-point precision. Dashed line, FDR = 0.05; FDR-significant genes are
labelled, grouped where they fall within 1.5 Mb. All FDR-significant signal in
the lower-powered round, and the majority under the meta outcome, falls at known
pigmentation and naevus loci.

**Fig 2 | HEIDI does not reject the LD confounding that colocalisation
identifies.** Each point is one exposure with both colocalisation and SMR/HEIDI
results, plotted as colocalisation PP.H4 (x) against −log₁₀ P_HEIDI (y); colours
as in Fig 1. Shaded region, PP.H4 < 0.5; dashed lines, PP.H4 = 0.7 and
P_HEIDI = 0.05. Points in the lower-left quadrant are records that colocalisation
assigns to distinct causal variants and that HEIDI nevertheless passes; the MC1R
cluster is labelled.

**Fig 3 | Outcome GWAS power determines colocalisation resolution, in both
directions.** Regional plots for PARP1 (top) and ZFYVE19 (bottom): cis-eQTL
(left), melanoma GWAS at FinnGen power (middle), melanoma GWAS at meta power
(right). Dotted vertical line, eQTL peak; circled point, peak of each panel, with
its distance from the eQTL peak annotated. For PARP1 the meta peak sharpens and
moves 25 kb away, and PP.H4 falls from 0.92 to 0.05 (false positive removed); for
ZFYVE19 the meta peak lands on the eQTL peak and PP.H4 reaches 0.99 (false
negative revealed).

**Fig 4 | Outcome GWAS power sets the reproducibility of the candidate list.**
**a**, Melanoma down-sampling curve: mean discoveries at FDR < 0.05 (left axis,
black) and Jaccard index of the recovered gene list against the full-power list
(right axis, green; shaded band, 5th–95th percentile over 200 replicates).
Stars mark the independently measured FinnGen round, which falls inside the
simulated interval — the pre-specified calibration condition. **b**, The same
curve for five cancers, plotted against each study's own observed case count;
dashed line, 50% recovery. **c**, Melanoma recovery stratified by locus class.
Known pigmentation/naevus-locus genes (red, n = 6) are recovered at 46.0% with
10% of observed power, where novel-locus genes (blue, n = 4) are at 0.4%;
reaching 50% recovery requires 2,506 versus 8,771 cases. Grey band, the
3,000–8,000-case range typical of studies using this framework. Curves are not
extrapolated above observed power.

**Fig 5 | Spatial compartment attribution.** ⟨四联图内容需与 `make_fig_spatial.py`
的最终面板对齐后补全⟩ Meta-analytic Spearman correlations of TPI1 with
compartment scores across eight sections; per-section concordance; compartment
contrasts after within-section centring; and the mirror-image behaviour of the
HLA-C positive control. Tissue-level TPI1 tracks the glycolytic module and the
tumour compartment and is negatively correlated with the lymphoid compartment,
while HLA-C shows the opposite pattern in the same sections.

**Fig 6 | TPI1: activation-window eQTL and pathway-level ICB association.**
**a**, Strongest cis-eQTL P value per activation profile for TPI1 and the
co-located SPSB2, in naive and memory lineages: TPI1 is a 16 h pulse in both
lineages, SPSB2 is constitutive. **b**, Glycolytic enzymes in CD4⁺ T cells,
responders versus non-responders, post-treatment (corrected sample definitions;
20 of 22 enzymes higher in non-responders). **c**, TPI1 and the TPI1-excluded
glycolysis module by response at both timepoints, with sample-level points.

**Fig 7 | Cross-cancer comparison.** **a**, Known-locus enrichment counted by
independent locus, for melanoma and five other cancers analysed with identical
exposures. **b**, Candidate gene effects across cancers (odds ratios with
significance markers): TPI1 is melanoma-restricted whereas SMC2 is concordant in
four cancers. **c**, Number of FDR-significant loci against outcome case count
across the six outcomes.

**Fig 8 | Multi-layer evidence matrix for candidate genes.** Six novel-locus
candidates against seven evidence layers (MR, colocalisation, no naevus effect,
no pleiotropy, single-cell ICB response, bulk replication, druggability and
safety). No candidate is strong on all axes; the statistically strongest
candidate and the functionally best-supported candidate are different genes.

---

## 6. Supplementary items (planned)

| Item | Content | Source tables |
|---|---|---|
| S1 | Instrument selection and harmonisation flow; F-statistic distribution | `01`, `04` |
| S2 | Steiger directionality with bounded R² | `09` |
| S3 | Four-estimator sensitivity comparison, including why weighted mode is inapplicable | `18`, `18b` |
| S4 | Naevus negative-control and comorbidity MR | `20` |
| S5 | Lineage-purity control: effect sizes before and after cell-level matching | `29a`–`29h` |
| S6 | Diagnosis of the two processing errors and before/after comparison | `31`–`32h` |
| S7 | Glycolysis axis characterisation and motif enrichment (all peaks vs activation-invariant) | `44a`–`44c`, `50a`–`50b` |
| S8 | Peripheral-blood ICC with upper-bound and random-module floor controls | `52a`–`52d` |
| S Table ⟨S⟩ | Record-level locus enrichment (locus-level reported in main text) | `36d`, `36e` |

---

## 7. References (to be completed)

以下为本项目实际使用、需在投稿前补全著录格式的文献与资源。⟨verify⟩ 表示
项目记录中只有简称，须核对完整著录。

**Primary data sources**

1. Soskic B, et al. Dynamic cis-eQTLs across CD4⁺ T cell activation. ⟨verify: 期刊/年/卷页⟩ — exposure data (8 profiles).
2. FinnGen R12 release. `C3_MELANOMA_SKIN_EXALLC`, `CD2_BENIGN_MELANOCYTIC`, `L12_VITILIGO`, and five cancer endpoints. ⟨verify: release citation⟩
3. Rashkin SR, et al. Pan-cancer GWAS in UK Biobank and Kaiser cohorts, 2020. GWAS Catalog GCST90011809. ⟨verify⟩
4. Landi MT, et al. Genome-wide association meta-analyses combining multiple risk phenotypes provide insights into the genetic architecture of cutaneous melanoma susceptibility. *Nat Genet* 2020. GCST010302/010303/010304.
5. Sade-Feldman M, et al. GSE120575 — single-cell profiling of melanoma under checkpoint blockade. ⟨verify⟩
6. Jerby-Arnon L, et al. GSE115978 — melanoma single-cell atlas. ⟨verify⟩
7. Hugo W, et al. GSE78220 — anti-PD-1 bulk transcriptomes. ⟨verify⟩
8. Riaz N, et al. GSE91061 — nivolumab pre/on-treatment transcriptomes. ⟨verify⟩
9. Thrane K, et al. Spatially resolved transcriptomics of metastatic melanoma, 2018. ⟨verify⟩
10. GSE282266 — purified CD4⁺ T cell 10x multiome activation time course. ⟨verify: associated publication⟩
11. GSE199994 — PBMC multiome, melanoma patients at baseline. ⟨verify: associated publication⟩
12. Byrska-Bishop M, et al. High-coverage 1000 Genomes Project GRCh38 reference panel, 2022. ⟨verify⟩

**Methods and resources**

13. Giambartolomei C, et al. Bayesian test for colocalisation (`coloc`). ⟨verify⟩
14. Zhu Z, et al. Integration of summary data from GWAS and eQTL studies (SMR/HEIDI). *Nat Genet* 2016. ⟨verify⟩
15. Hemani G, et al. The MR-Base platform / TwoSampleMR. ⟨verify⟩
16. Hao Y, et al. Seurat v5. ⟨verify⟩
17. Fornes O, et al. JASPAR2020. ⟨verify⟩
18. Chang CC, et al. PLINK 2. ⟨verify⟩
19. eQTL Catalogue; DICE; GWAS Catalog; Open Targets Platform; DepMap. ⟨verify: 各自的引用格式与访问日期⟩

**Comparator studies using the same framework** ⟨须补全三篇：正文 §1、§2.1(i)、
§2.2、§3.3① 均引用之。项目记录中一篇为 Adv Sci 2025 的乳腺癌研究（PADI4、
MDM4/BI-907828），另两篇见 `D:/文章/孟德尔/pdf_text/`⟩

20. ⟨Comparator study 1 — Adv Sci 2025, dynamic CD4 eQTL MR, breast cancer⟩
21. ⟨Comparator study 2⟩
22. ⟨Comparator study 3⟩

---

## 8. 待补清单（不进正文）

**必须在投稿前补**

1. **三篇对照文献的完整著录**（§7 第 20–22 条）；正文四处引用位置已标好。
2. **参考文献著录格式**：§7 中 ⟨verify⟩ 的 16 条，尤其 GSE282266 与 GSE199994
   的关联论文（本项目只用了 GEO 号）。
3. **Abstract 的关键词**与投稿期刊要求的结构（现为非结构式）。
4. **软件版本逐一核对**（Methods §19）。
5. **耗竭评分基因集的声明**：原始集未记录，现用标准集，须写明与修正前版本不严格可比。
6. **Fig 5 图注**需与 `make_fig_spatial.py` 的最终四联面板对齐。
7. **Supplementary Table ⟨S⟩ 编号**：记录级位点富集表（§2.2、§3.3① 均指向它）。

**可选，待决定**

8. Fig 7a（按独立位点的跨癌种富集）是否并入 Fig 1 —— 二者都在服务发现①。
9. §2.7 的两个数据处理错误是否保留在 Results。**建议保留**：它们是发现⑦的实例，
   移入 Methods 会使 HLA-C 机制的撤回失去解释。
10. 标题选用（见文首三个备选）。
11. 数据与代码发布范围（Methods §20）：建议全部 `step*` 脚本 + `32a`–`32h`、
    `36a`–`36e`、`43a`–`50b`、`53a`–`55a` 结果表随文发布。

**已统一、勿再改动的口径**

- 位点富集一律报**独立位点级 4.09×（P=0.028）**，记录级 5.26% 仅入补充材料
- 功效曲线校准点为 **10.2 模拟 vs 10 真实**（`53b`），不是 9.6
- 单细胞糖酵解模块 = **22 基因（不含 ENO3）**；GSE282266 的轴 = **16 基因**版本
- 单细胞患者层面数字一律取 `32a`–`32h`（修正后），`19*`/`28*` 全部作废
