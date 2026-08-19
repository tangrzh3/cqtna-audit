**⚠ 已被 `MANUSCRIPT_GB.md` 取代。勿引用本文件的数字。**

# Results — draft v1 (2026-08-11)

> 结构说明（中文，不进正文）
> 顺序与 Discussion 的论证骨架一致：先证明管线可靠（阳性对照），再逐层拆掉自己的靶点主张。
> 每个否定紧跟其阳性对照，避免读成"阴性结果论文"。
> 数字一律取定论表（meta 轮；单细胞取 `32a`–`32h`；位点富集按独立位点计）。
> 未定稿处用 ⟨⟩ 标出。

---

## 1. A dynamic CD4⁺ T cell eQTL–melanoma MR pipeline that passes its positive controls

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

## 2. Finding ①: the significant signal is pigmentation genetics, not immunology

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

## 3. Finding ②: SMR/HEIDI does not exclude the LD confounding that colocalisation identifies

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

## 4. Findings ③ and ④: outcome power moves MR and colocalisation in opposite directions, and the candidate list turns over completely

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

## 5. The reproducible part of the candidate list is the part that is not a discovery

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

## 5b. A natural experiment: five sequential releases of the same GWAS resource

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

## 6. What can locate an effect once MR cannot: orthogonal functional evidence

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
(Section 8). The three lines come from different data types, and this generalises
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

## 7. Finding ⑦: the functional layer has its own failure modes, and two of ours were severe

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

## 8. Finding ⑥: MR names the gene that has an instrument

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

## 9. TPI1 as a bounded worked example

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
activation-invariant peaks is mitigation rather than separation.

**What we do not claim.** The MR establishes a causal relationship with melanoma
**risk**; the association with ICB response is observational, from a single
CD4-resolved cohort, and we do not join the two into a causal chain. TPI1 is
DepMap-essential with pLI = 0.87, so no target claim is made. The activation
window has no independent replication: TPI1 is not quantified in any of the eight
stimulated T-cell eQTL datasets available, and the one resource matching our
design is the same experiment as our exposure data, so using it would be
circular. And by finding ④ and Fig 4, TPI1 itself may not survive a
higher-powered outcome GWAS.

## 10. No candidate is strong on all axes (Fig 8)

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

### 待补（写正文前需确认）

1. ⟨Fig 7a 是否并入 Fig 1⟩ —— 二者都在服务发现①，目前分散在两张图。
2. ⟨§2 的记录级 5.26× 是否只放补充材料⟩ —— 与 `FIGURES_plan.md` 的第 2 条一致处理。
3. ⟨三篇参考文献的具体引用位置⟩ —— §1(i)、§2、§3 各需一处。
4. §7 两个数据错误的写法需与期刊风格协调：现按"方法学发现"正面陈述，
   另一选项是移入 Methods/Supplementary。**建议保留在 Results**——
   它们是发现⑦的实例，且删掉会使 HLA-C 机制的撤回失去解释。
