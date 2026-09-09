# Target nomination from single-variant cis-eQTL Mendelian randomization re-reads the outcome GWAS: an audit across six diseases and two exposure resources

⟨author list⟩

⟨affiliations⟩

**Correspondence**: ⟨corresponding author⟩, ⟨corresponding email⟩

## Abstract

**Background.** Context-specific eQTLs with Mendelian randomization are a
recurring design for nominating immune targets in cancer. With one instrument and first-order
standard errors the Wald statistic reduces to |z| = |β_out|/se_out: the exposure
sets the sign, the scale and which variants are eligible, and the outcome GWAS
supplies the whole significance claim. We asked what such a list contains.

**Results.** Instrumenting CD4⁺ T-cell cis-eQTLs from eight activation profiles
against a 12,530-case melanoma meta-analysis, we confirmed the identity
numerically; Benjamini–Hochberg at 0.05 corresponds to an outcome P of
2.1 × 10⁻⁴. Significant signal concentrates on loci already
reported for the outcome, 4.96-fold by bounded locus (P = 0.0048) — but that
enrichment is carried entirely by the three loci the outcome GWAS had already
found unaided at 5 × 10⁻⁸. Removing them leaves five loci, one of them known,
and no locus-level attribution signal we can detect (1.98-fold, P = 0.41).
Naming the gene fails separately: at ten melanoma loci with an accepted causal
gene the pipeline names six, alone at two, and at MC1R it spans sixteen genes
without naming MC1R. Replacing the outcome with a higher-powered meta-analysis
replaced the list entirely. Among 45 ascertainable papers in a two-coder subsample of 46, 10 (22.2%;
Wilson 95% CI 12.5–36.3%) used the exact single-variant Wald design.

**Conclusions.** The reproducible part of such a list is the part that is not a
discovery, and the part that would be a discovery carries no locus-level evidence
we could detect — an absence of evidence, at five loci, not evidence of absence. Eight inexpensive checks follow, nine with the decomposition of the
identity, with an R implementation.

**Keywords** Mendelian randomization · context-specific eQTL · target nomination ·
colocalisation · statistical power · reproducibility · melanoma

---

## Background

Most cis-eQTLs are context dependent, and loci mapped in resting bulk tissue
cannot capture regulatory variants that act only while a cell responds to a
stimulus. Datasets profiling primary immune cells across a stimulation time
course were built to close that gap, and their combination with MR has become an
attractive route to target nomination: genotype is fixed before disease, cis-eQTL
effects carry a directional prior, and the exposure is measured in the cell type
through which the effect is proposed to act. Several recent studies apply exactly
this design — activation-time-course CD4⁺ T cell eQTLs as instruments, cancer
GWAS as outcomes — and report novel immune targets on that basis [1–3].

The inference passes through more joints than its summary statistics reveal.
Because the instrument is usually a single variant, the P value is supplied
entirely by the outcome GWAS. Whether the eQTL and the disease signal share a
causal variant, rather than lying in linkage disequilibrium with a large
neighbouring effect, is decided by colocalisation, whose resolution depends on
outcome power. Which of two genes at a locus carries the effect is not decided by
MR at all; neither is the cell type in which it operates, because functional
replication is usually available only at tissue level, where a broadly expressed
gene reports on whichever compartment dominates. Each joint fails invisibly
unless specifically tested.

Melanoma is an informative setting for two reasons that pull in opposite
directions. Its common-variant architecture is dominated by pigmentation and
naevus loci with effects far larger than anything expected from immune
regulation — MC1R alone exceeds P = 10⁻⁵⁹ — which is simultaneously a built-in
end-to-end positive control and a built-in confounder, since long-range LD around
such loci can present as a causal association for any neighbouring gene. It is
also where CD4⁺ T cell biology has the most direct clinical relevance, so
candidates can be examined against response-stratified single-cell data rather
than annotation alone.

The organising question of this audit is not whether any particular gene is a
target, but what the nomination is a function of.

---

## Results

### Part 1 — Why the nomination is unstable

The first four results establish what a candidate list from this framework is a
function of, and how much of it survives a change in the outcome GWAS.

### Design and positive controls

Eight CD4⁺ T cell activation profiles (naive and memory, 0 h to 5 d) provided
exposures. Instruments were selected at eQTL P < 5×10⁻⁸ with F > 10, one variant
per exposure, giving 3,556 records after harmonisation with the outcome. Every
pre-specified control passed. Recovery of all 157 known melanoma loci, 136 of
them strengthened, with genome-wide significant variants rising from 2,871 to
4,552, is an **integrity check on the outcome data and the meta-analysis step**
rather than a positive control on the nomination pipeline — it involves no
instrument, no MR and no colocalisation. The end-to-end controls are that PARP1
reproduced its published direction and that the MC1R region gave the strongest
associations in the study, at P = 4×10⁻³⁷ (Fig. 1).

One property of the design governs everything that follows. With one instrument
per exposure, the Wald z equals β_out/se_out, so **given the instrument set**, MR
significance is a property of the outcome GWAS; the exposure data determine which
variants are asked about and how many there are, but not which of them answers.
Steiger filtering cannot arbitrate direction here either, since cis-eQTL
and disease R² differ by orders of magnitude and the test passes almost by
construction. We therefore state at the outset what the design can and cannot
decide, and test the consequences empirically rather than asserting them.

### The nomination is the outcome GWAS, read at a lower threshold

With one instrument per exposure the Wald statistic is b = β_out/β_exp with
SE = se_out/|β_exp|, so |z| = |β_out|/se_out. The exposure sets the sign and the
scale of the effect estimate and the set of variants eligible to be tested; it
does not enter the test statistic. We verified this on the two tables that store
both sides: `se` matches se_out/|β_exp| to within 0, and the reported P matches
2Φ(−|z_out|) to within 3 × 10⁻¹⁶, while the sign differs on 51% of records
(Supplementary S42). **A candidate list from this design is therefore the outcome
GWAS restricted to variants that happen to be lead cis-eQTLs, thresholded by the
multiple-testing burden rather than by 5 × 10⁻⁸.** Everything below follows from
that, and we report it as the mechanism rather than as a caveat.

The threshold is the first consequence. Benjamini–Hochberg at 0.05 over the 3,556
records corresponds to an outcome P of 2.1 × 10⁻⁴ — three and a half orders of
magnitude more permissive than genome-wide significance.

The second is what the resulting list contains. In the FinnGen round, all ten
FDR-significant records fell within 1 Mb of a known
pigmentation or naevus locus — VPS9D1-AS1 50 kb from MC1R, CDK10 48 kb from the
MC1R R151C variant, PARP1 5–14 kb from PARP1 — with no immune signal at
FDR < 0.05 at all. Under the meta outcome, counting bounded loci rather than
gene records, 4 of 8 significant loci were known-locus loci: a 4.96-fold
enrichment over the 10.1% background (one-sided P = 0.0048; Fig. 1). The same
enrichment counted at the record level, which is the unit the source literature
reports, is Supplementary S10.

**That enrichment is carried entirely by loci the outcome GWAS had already found
without any exposure data.** Three of the eight reach 5 × 10⁻⁸ unaided — two at
16q24.3 over MC1R and one at PARP1 — and all three are on the published list,
which is how they came to be published. The cross-tabulation is completely
separated: every genome-wide significant locus is a known locus and no
genome-wide significant locus is novel. Removing those three leaves five loci,
one of them known: 1.98-fold, one-sided P = 0.41, and 2.07-fold, P = 0.40 when
the background is restricted to sub-threshold loci as well, which is the
conditional comparison the question implies (Supplementary S42). **What this
design adds beyond genome-wide significance shows no locus-level attribution
signal at all** — though at five loci it is also underpowered to show one, so this
is an absence of evidence and we do not read it as evidence of absence.

Two of the eight bounded loci are not two published regions. The windows at
16:88.86–89.73 Mb and 16:89.87 Mb fall within 1 Mb of an identical set of seven
Landi lead variants, so the fixed-anchor partition counts the MC1R region twice;
the three genome-wide-significant loci are two distinct published regions, not
three. Merging bounded loci that share an attributed lead variant, and applying
the same rule to the background, gives seven regions of which three are known
against a background of 52 of 641: 5.28-fold, P = 0.014 (Supplementary S42).
The direction is unchanged and the evidence is weaker than the unmerged count
suggests, which is the honest reading of a partition that bounds window width
without guaranteeing that adjacent windows are independent. The same
decomposition is consistent with the colocalisation result reported below, where
2 of 284 records reach PP.H4 > 0.8. It does not imply it: a variant that drives
both expression and disease can colocalise well however the significance arose,
so the low rate is an empirical finding compatible with outcome-dominated
selection rather than a consequence of it.

Loci here
are the non-recursive fixed-anchor partition fixed in Supplementary S37 — a
1 Mb window claimed from each anchor, so a locus spans at most 1 Mb whatever the
variant density; the single-linkage rule used in the source literature chains a
dense resource into blocks of tens of megabases and is reported only where we
reproduce published numbers (Methods). Applied unchanged to five other cancers,
each scored against its own known-locus list rather than against melanoma's
(Supplementary S41), the attribution recurs: own-list enrichment survives Holm
correction across all five outcomes that have a usable list — melanoma 4.96-fold
(P = 0.0097 corrected), lung 3.85-fold (0.0138), colorectal 3.15-fold (0.0029),
breast 2.25-fold (0.0011) and prostate 2.04-fold (2.3 × 10⁻⁵); pancreas has no
list that clears the size and build-rate gates and is not scored. The five
outcomes share almost no significant loci with one another (melanoma and lung
share none; 5 of 69 distinct loci are significant in more than one outcome). Two
qualifications bound this. The first is that folds are not comparable across
diseases at all. Fold divides by the background known share, so its ceiling is
the reciprocal of that share, and across our lists the share runs from 0.101 for
melanoma to 0.403 for prostate — ceilings of 9.92 and 2.48. Matching the
diseases on power by scaling each until it yields the same number of significant
loci, and scoring them with the chance-corrected share
A = (p_sig − p_bg)/(1 − p_bg), the two statistics rank them almost in reverse: at
ten loci, fold gives melanoma 5.95, lung 3.85, breast 3.10 and colorectal 2.52,
while A gives breast 0.859, colorectal 0.707, melanoma 0.555 and lung 0.330
(Supplementary S43). **Melanoma's apparent lead is a property of its reference
list being sparse, not of melanoma.** The second is that only melanoma exceeds
every mismatched list by a margin (2.30-fold against
the best rival, against 1.29 for lung and below 1 for the rest), while lung's
result does not survive dropping its single strongest locus (P = 0.0138 to
0.0569 without the chr11 FADS1/TMEM258 cluster, the same cluster identified as
the driver of lung's enrichment in our earlier melanoma-list comparison). We
therefore report transport as recurring but read it as confirmatory only for
melanoma.

We then tested the same attribution three further ways, changing one thing at a
time (Fig. 2).

#### Five nested power levels of one resource

FinnGen's sequential releases carry
an identical melanoma endpoint at 2,705, 2,993, 3,194, 3,932 and 5,753 cases,
with exposure data, instruments, harmonisation, estimator, testing family and
locus definition held fixed. Predictions were generated by down-sampling the
largest release and registered before any earlier release was examined; all twelve
— discovery counts, overall recovery, and the known-minus-novel differential at
each release — fell inside the registered intervals. The known-locus part of the
list is recovered early and then does not move: the same five genes constitute the
significant list at three consecutive releases without change, four already
present at 47% of reference power. And **no novel-locus gene reaches FDR < 0.05 in
any release**; the per-release lists are in Supplementary S11. Because releases are nested, agreement between them exceeds that
between independent studies of equal size, so this is an upper bound on stability.

That series stops where it does for a reason. The release after it, R13, retired
the endpoint we had used and replaced it with one whose cases are defined by the
same codes but whose controls are screened by a different exclusion rule, so it
cannot be a sixth power level; treating it as one would score a change of
phenotype definition as a change in power. Analysed instead as a transfer test —
pre-registered, with the direction of every difference recorded before the result
was seen — the list survives it unchanged: at 6,226 cases the significant list is
the same six genes at the same two loci, still with no novel-locus gene, and
locus attribution is 10.0-fold (P = 9.1×10⁻⁵), all four significant loci known. The stability is not an artefact of
reusing data: the two releases' statistics correlate at |z| r = 0.937 rather than
being identical, standard errors shrink by a median 3.85% against the 3.8% the
added cases predict — the internal check that R13 is the higher-powered
computation and not a copy — and 51 of the 275 nominally significant records are
replaced. What holds still is the
significant tier; the tier below it churns. Because the extra cases and the
broadened control exclusion both push towards more signal, a successful transfer
here is partly confounded and we do not read it as power alone.

#### A second disease

Repeating the nomination against hepatocellular carcinoma
with the exposure side unaltered, at two outcome power levels (3,748 and 947
cases), gives four significant loci across both levels, three carrying a known
HCC lead SNP: SAMM50, 26 kb from rs2294915 in the *PNPLA3* region, and ZNF506 in
the chr19 *TM6SF2* region at both levels. Enrichment is 15.5-fold in the
lower-powered outcome (P = 0.0041) and 7.9-fold in the higher-powered one, where
it does not reach significance (P = 0.123). The single novel nomination, SUPV3L1,
appears only at high power and has FDR = 0.98 at low power.

Ten of that reference list's 73 loci were taken from the higher-powered outcome
publication's own table, which is circular to that extent. Checked individually
against the GWAS Catalog, 7 carry heavy independent support (APOE in 1,399
studies, *TM6SF2* E167K in 1,073) and **3 rest on the outcome publication alone**;
removing those three leaves both enrichments unchanged to the same decimals,
because no instrument lies within 1 Mb of any of them. We report that as a
measurement of the circularity — 3 of 73, none near a significant locus — not as
a robustness check, since by construction it could not have moved the result.
Dropping all 10 instead *raises* the enrichment, as removing known loci only
shrinks the background, so over-correction flatters this test; and the Catalog
carries no cohort-level information, so this addresses circular locus attribution
and not sample overlap (Supplementary S29). The melanoma reference list comes
from a separate publication and does not have this structure, though sample
overlap with the outcome meta-analysis cannot be excluded.

The higher-powered
HCC study carries 30.3% of melanoma's effective sample size, so we down-sampled
melanoma to match: the predicted median is 2 significant loci [5–95%: 1–6], one
known and one novel, and **the observed HCC counts are compatible with that
melanoma-derived prediction interval**. The non-significant primary test is what
a test with two significant loci is expected to deliver. Compatibility with an
interval spanning 1 to 6 loci is not equivalence — that interval would also
accommodate outcomes we would have read as a difference — so this shows the HCC
result is not evidence against the pattern, not that the two diseases behave
alike once power is equalised.

#### A second exposure resource

Holding the melanoma outcome fixed byte-for-byte
and replacing the exposure entirely with eQTLGen whole-blood cis-eQTLs
(n = 31,684, roughly 300-fold larger, different cells, different platform) raises
instruments from 3,556 to 12,835 records and significant loci from 7 to 30. The
attribution does not move: **23 of 34 significant loci (67.6%) carry a known
melanoma, naevus or pigmentation lead SNP, 6.31-fold over this resource's own
10.7% background (P = 1.3×10⁻¹⁵)**, against 4.96-fold in CD4⁺ T cells. A
matched-background version agrees, drawing each significant locus a background
locus of the same eQTL-p decile and allele-frequency quintile (6.55-fold,
empirical P ≤ 1×10⁻⁴, the floor at 10,000 permutations; 5.66-fold on the decile alone). The two
resources' novel nominations intersect in exactly one gene, ZFYVE19, so the
stronger claim that novel nominations never reproduce is not supported.

#### Three further exposure resources, and why they settle nothing

We tried three more CD4⁺ resources — Nathan et al. [4], Schmiedel et al. [5]
and the influenza-stimulation panel of Randolph et al. [6] (eQTL Catalogue
QTD000588). None settles the question, and the identity
says why: the exposure fixes which variants are eligible for testing, not which
of them are significant, so a resource with more instruments enlarges the tested
set without changing what makes an entry cross the threshold. Two returned too
few harmonised records to score, and the third fails its mismatched control at
the registered setting. The per-resource counts and gates are in Supplementary
S35.

#### The whole-blood rheumatoid arthritis cell is void

Scored on the same records, it reaches
3.94-fold (P = 3.3×10⁻¹⁷), but its pre-registered mismatched-list control also
enriches — melanoma's known loci give 2.05-fold on the RA nominations,
P = 2.9×10⁻⁴ — and S33 §5 voids a cell whose mismatched list enriches. It stays
void: excluding the MHC does not clean it (1.83-fold, P = 0.0052), and the
500 kb partition is a pre-registered sensitivity, not a second chance. **The
non-cancer generalisation therefore rests on the CD4⁺ cell alone.**

**The control that clears this cell is clean at exactly one setting.** We swept
both 1 Mb conventions — the window within which a locus counts as known, and the
width of a locus — over 100, 250, 500 and 1000 kb (Supplementary S36). Across the
seven distinct settings, the RA CD4⁺ cell's mismatched-list control fails at six
and is clean only at the registered 1 Mb main analysis. It fails in a direction
that matters: as the known-locus window narrows, the wrong disease's list
enriches **more**, monotonically — 1.71-, 2.38-, 3.03- and 4.17-fold at 1000,
500, 250 and 100 kb. So the control does not pass at 1 Mb because attribution
there is specific; it passes because 1 Mb is coarse enough that both reference
lists hit a great deal and the ratio is diluted towards one. **The cell clears
its negative control at the setting where that control has least resolution.**

Two observations locate the cause. First, the contrast is within one exposure
resource: melanoma and RA are scored on the same Soskic instruments, the same
partition and the same convention, and the melanoma cell's mismatch enrichment
is 0.00-fold at all eight settings while the RA cell's fails at six. The
difference is the disease pair, not the machinery — RA and melanoma share immune
loci, melanoma and hepatocellular carcinoma do not. We note against ourselves
that the fragility is not confined to RA: melanoma on whole blood traces almost
the same mismatch curve (4.11-, 3.44-, 2.55-, 1.59-fold as the window widens) and
also fails at the tightest setting, so a narrow window on a dense resource
strains this control generally. What is specific to RA is how far the failure
extends — six settings against one. Second, a per-locus breakdown
rules out the convention: of the eight loci the wrong list claims at 500 kb, none
is claimed through a non-significant record; every one has a significant
nomination genuinely within 1 Mb of a melanoma lead. Consistent with that, RA's
hit distances are spread across the intermediate range rather than bimodal — 6 of
41 within 50 kb against 15 of 34 for melanoma on whole blood — so RA's hits sit
precisely in the band that a narrowing window removes.

**A mismatched list is therefore not an orthogonal negative control for this
outcome pair.** The enrichment itself is not the fragile part: it is flat across
locus widths (3.70- to 3.96-fold) and rises as the known-locus window narrows
(3.91- to 7.30-fold), the direction every other cell moves in. It is the control
that is scale-bound, not the estimate.

#### Which disease is nominated as the mismatch changes the verdict

A single
mismatched list cannot separate specific attribution from a fortunate pick, so we
scored every cell against every reference list in the study with at least 30
placeable lead SNPs belonging to a different disease — six per cell
(Supplementary S39). **This analysis was designed after the results it comments
on and its criterion is weaker than the registered one, so we report it as a
post-hoc diagnostic and draw no confirmatory conclusion from it.** Its one
informative result is about the design rather than about RA: of the six
unrelated lists, only one leaves the RA CD4⁺ cell standing — melanoma, the list
this study designated, at 1.71-fold (P = 0.10) — while **hepatocellular
carcinoma's at 3.05-fold (P = 0.0028), lung's at 3.32-fold (P = 1.5 × 10⁻⁵),
colorectal's at 1.71-fold (P = 0.0018), breast's at 1.69-fold (P = 0.0042) and
prostate's at 1.39-fold (P = 0.026) would each have voided it.** Whether this cell survives
its negative control therefore depends on which comparator was nominated, which
is a property of single-list controls in general and one we would apply to any
study using them.

The ranking that falls out is descriptive. Taking each cell's own fold against
the largest fold from any of its six comparators, the six cancer cells span
2.31- to 5.52-fold while **the two RA cells sit at 1.18- and 1.32-fold, the
narrowest in the grid**; RA on CD4⁺ T cells also inverts at the narrowest locus
width, the only cell at any of 64 settings that does.

**We therefore treat the non-cancer cell as exploratory.** Its enrichment is a
real observation on its own reference list, and the numbers and the full sweeps
are given here and in Supplementary S36 and S39; but a cell that clears its
registered negative control at one of seven settings, whose survival depends on
which comparator was named, and whose margin over unrelated diseases is the
smallest present, does not carry the same weight as the tumour cells. **The
confirmatory claim in this paper is that the attribution appears in melanoma and
in hepatocellular carcinoma. Whether it extends beyond cancer is open.**

Two further limits. Okada 2014 predates a decade of later RA loci, which count as
novel here and bias the 3.91-fold downward; for the same reason enrichment
magnitudes are **not comparable across diseases**, each reference list differing
in age and completeness. The strongest supported statement is correspondingly
narrow. Within cancer, enrichment on known outcome loci is not confined to one
tumour or one exposure resource, and its negative controls are clean at every
setting we swept. Outside cancer we can say that the same concentration appears
in an autoimmune outcome whose pathogenic compartment is the exposure cell type,
but **not** that it is specific to that outcome's own genetics, because the
control which would establish that does not hold away from the registered
window. Generality beyond cancer is therefore open, not established.

#### Both axes crossed

The main grid is the three diseases at their higher-powered outcome crossed with
the two exposure resources — six cells — with the lower-powered hepatocellular
outcome held back as a power sensitivity. Scored against each disease's own known
loci, one cell is void on its mismatched control and the other five enrich:

| Cell | Own fold | One-sided P |
|---|---:|---:|
| melanoma × CD4⁺ | 4.96 | 0.0048 |
| melanoma × whole blood | 6.31 | 1.3 × 10⁻¹⁵ |
| HCC-high × CD4⁺ | 7.88 | 0.123 |
| HCC-high × whole blood | 10.87 | 0.0015 |
| RA × CD4⁺ | 3.91 | 2.3 × 10⁻⁷ |
| RA × whole blood | — | void on its mismatched control |

Four of the five reach nominal significance and remain so after
Benjamini–Hochberg across the five. **That is not the relevant caveat**: the
cells share exposure data and reference lists, so they are consistency of
direction rather than five independent replications. The full cell-by-cell
account, the registration history and the power sensitivity are in
Supplementary S48.

Two permutations close the gap from different sides, because a locus dense in
records is not the same thing as a locus dense in a way that is itself
disease-specific. Matching background loci on eQTL-p decile and
allele-frequency quintile, which controls instrument strength and frequency,
gives 6.66-fold (empirical P = 0.0016) for melanoma on CD4⁺ T cells. Matching on
density — record count, span and gene count, with complete matched coverage
required and the whole tolerance scan reported, under rules fixed in
Supplementary S38 before the run — gives **4.87-fold (P = 0.0050)** there,
6.14-fold on whole blood, 8.24-fold for HCC-high on whole blood, 12.85- and
13.53-fold for the two HCC-low cells, and **4.52-fold for RA on CD4⁺ T cells**.
Values quoted at ≤ 1 × 10⁻⁴ sit at the resolution floor of 10,000 permutations
and are not point estimates. On bounded loci, density matching therefore barely
attenuates the melanoma cell, 4.96-fold unmatched against 4.87-fold matched.
**That is a statement about one specified null, not about density in general**
— the 3.44-fold reported earlier came from quantile-stratified matching, a
different test we do not treat as the same quantity (Supplementary S38 §4).

### MR and colocalisation moved in opposite directions when the outcome was replaced

Replacing the FinnGen outcome with the meta-analysis changed outcome power and
study composition together, so what follows is an observation about that
replacement and not a measurement of power. Meta-analysis raised MR discoveries
from 10 to 21 records while *lowering* colocalisation support: across the 127 exposures run in both rounds, median
PP.H3+H4 fell from 0.249 to 0.207 and median PP.H4 from 0.124 to 0.095, with only
44.1% improving. Three explanations were tested and rejected, including
winner's-curse selection, for which newly entering candidates had *lower*
PP.H3+H4 (0.125 versus 0.202).

The movement is not an artefact of analysis window or prior. Across four windows
(±100 kb to the full cis region) and four values of the shared-causal-variant
prior, the negative results are unchanged in all 16 combinations; window size is
essentially irrelevant, and what varies is the prior, and only for loci near the
criterion boundary — TPI1 and SPSB2 pass in 8 of 16 combinations, SMC2 and
KIAA0040 in 12 of 16. For PARP1, five orders of magnitude more signal sharpened
the peak and moved it 25 kb, to 54 kb from the eQTL peak, and PP.H4 fell from 0.92
to 0.05; in the other direction, the meta peak for ZFYVE19 landed exactly on the
eQTL peak and PP.H4 rose to 0.99. One general fact belongs with this: FinnGen
fine-maps only 19 regions genome-wide for this endpoint, all classical
pigmentation or melanoma loci, and none of our candidate loci is among them. At
this outcome power the candidate list is built almost entirely on regions the
outcome GWAS cannot fine-map.

**The candidate list turned over completely.** The FinnGen round produced four
novel-locus genes passing triple validation (PRPSAP2, IMPA1, GCC2, PADI4) and one
passing the full battery. Under the meta outcome none survived — IMPA1's MR P moved
from 9.3×10⁻⁴ to 0.11 — and six different genes took their place. **The two lists
share no genes**, from identical exposure data under an identical pipeline.

### The reproducible part of the list is the part that is not a discovery

We measured candidate-list recovery as a function of outcome power by
down-sampling in summary space, recomputing Wald ratios and FDR, and comparing
with the full-power list (Fig. 3). The simulation is calibrated empirically — not
externally, since FinnGen contributes to the meta outcome — against a separately
observed lower-power round, which it reproduces (10.2 simulated versus 10 real
discoveries; simulated Jaccard 0.52 [0.36, 0.73] containing the observed 0.455).
Curves are reported only up to observed power.

Melanoma recovers 53% of the full-power list at 5,000 cases and 85% at 10,000, and
five other cancers behave the same way. Stratifying by locus class dissolves that
apparent stability: at 10% of observed power, 46.0% of known pigmentation and
naevus-locus genes are recovered against **0.4%** of novel-locus genes; at 50%
power, 85.8% versus 22.8%. Reaching 50% recovery requires 2,506 cases for
known-locus genes and 8,771 for novel-locus genes — 3.5-fold. Studies using this
framework commonly analyse 3,000–8,000 cases, where novel-locus recovery is 6–40%.

That differential is mostly a statement about effect size, and we tested how much
of it survives conditioning on effect size. Candidates at novel loci sit against
the detection threshold by construction — their full-power |z| all lie between
3.73 and 4.55, while known-locus candidates reach 15.99 (medians 4.23 and 11.07) —
and recovery under down-sampling is a monotone function of |z|. **A model using
|z| with no class label reproduces 68–80% of the gap**; matched on |z|, the
residual is **+5.2 percentage points [−1.6, +12.0] by bounded locus** and
+6.4 [+0.9, +12.0] by gene, against raw gaps of 47.6 and 63.0. Only one known
locus falls inside the novel |z| range, so the remaining 20–32% is not separable
from the failure of matching and we do not attribute it. What we retain is
threshold proximity rather than a property of the category — and the practical
consequence is the same either way: **for anyone holding such a list, the
novel-locus part is the part that will not replicate** (Supplementary S28).

> In this analysis, and within the range of outcome power we could observe, the
> reproducible part of the candidate list was precisely the part that did not
> constitute a discovery.

### Part 2 — Why a significant locus is not a target gene

Even where the signal is stable, the step from a significant locus to a named
gene fails independently, and fails in ways the significance test cannot
register.

### Colocalisation and SMR/HEIDI disagree, and the disagreement is consequential

Colocalisation assigned the MC1R-region signals — the strongest MR associations in
the study — to distinct causal variants, with PP.H3 dominant and PP.H4 at or near
zero (CHMP1A 0.99; VPS9D1-AS1 0.96; SPATA33 0.92–1.00). FinnGen's own in-sample
fine-mapping of this endpoint independently resolves three high-purity credible
sets in the region (log₁₀BF 45.6, 36.6, 19.5), so it genuinely carries multiple
independent causal signals: precisely the configuration in which a neighbouring
gene's eQTL can be tagged by one of them without sharing it.

SMR/HEIDI did not reproduce this. Of 291 records that colocalisation assigned to
distinct causal variants, HEIDI failed to reject homogeneity for 253 (86.9%),
including VPS9D1-AS1 (P_HEIDI = 0.649) and CDK10 (0.086, 0.081), on 8–20 variants
per test (Fig. 4). Evidence tiers that accept MR plus SMR without colocalisation
would have reported MC1R linkage spillover as CD4-mediated immune targets.

### Does the pipeline attribute the right gene where the answer is known?

The criticism this paper makes of nomination applies to this paper, and can be
administered to it. At melanoma loci where a causal gene is generally accepted, we
asked whether the pipeline's own significant nomination names that gene. Of ten
such loci reached under the larger exposure resource, **the accepted gene is
among those named at six and absent at four — and is named alone at only two**,
and the failures are the informative ones. At the MC1R region —
the strongest melanoma locus in the genome — the nomination spans **sixteen genes
and MC1R is not among them**. At the OCA2/HERC2 locus it names the pseudogene
HERC2P9 rather than HERC2. At *TYR* it names ODF3, and at the CDKN2A/MTAP locus
C9orf66. Where it succeeds it often does so cleanly: IRF4 and MX2 are each named
alone and correctly.

Two limits bound how far this can be pushed. The accepted-gene list is one we
assembled ourselves — fixed before the comparison was run, but **not
pre-registered** — and ten loci is a small denominator, so **six of ten is an
illustration and not an error rate**. What the check does establish is
directional and does not need precision: a framework that misassigns the gene at
the best-characterised locus in its own disease should not be read as assigning
genes at uncharacterised ones.

Turning it into a rate would need an external, pre-defined benchmark, and we
pre-registered the attempt (Supplementary S34) with a floor of 30 evaluable loci
below which no rate would be reported. It could not be met: scoring our 30
significant loci against 2,435 curated GWAS locus-to-gene assignments, restricted
in advance to evidence classes that do not themselves rest on molecular-QTL
colocalisation, leaves **16 evaluable loci and only 2 when the benchmark is
restricted to melanoma and related traits**. The binding constraint is coverage,
not circularity — the class restriction costs two loci. For a disease whose
architecture is dominated by pigmentation and naevus loci, curated locus-to-gene
truth barely exists, and three-quarters of what exists here comes from
drug-target pairs, which follow therapeutic attention rather than allelic
architecture.

This is the co-regulation problem of Tambets et al. [7] observed at the top of
the effect-size distribution, and it is why we report compartment attribution
separately from gene attribution — the two fail independently.

### Which unit the significance belongs to

Recomputing the significant list over record, variant, gene and bounded locus
changes its length: a shorter list under a stricter unit is not evidence that
FDR was controlled under dependence, and the point is to state which unit a
conclusion is in rather than to pick the flattering one. The four counts, the
Simes collapsing and the consequences for every number quoted here are in
Supplementary S50.

### Part 3 — One nomination, followed as far as it goes

The identity says a nomination list is the outcome GWAS at a lower threshold. It
does not say what happens when someone takes a single entry from such a list and
tries to substantiate it, which is what the field actually does. We did that
once, for the one glycolytic gene this pipeline nominates, and report it as a
worked example subordinate to the audit rather than as a second finding.

**It fails at every tier, and it fails in the way the identity predicts.** The
gene is nominal, not resolved genetically; its apparent time specificity is a
statement about precision rather than amplitude; the axis it marks is not
specific to it; the tissue in which the signal is largest is not the tissue the
instrument came from; the functional contrasts dissolve into lineage composition
once that is matched; and the patient-level result appears in each cohort and
replicates in none twice.

The six rungs are recorded in full in Supplementary S47 and summarised here.
Of 28 glycolytic genes, 3 are instrumentable in this exposure resource, 2 remain
analysable against this outcome after harmonisation, and 1 yields a nominal
association; the full gene-by-profile availability matrix is Supplementary S13
(Fig. 5). For that gene, TPI1, a usable instrument exists at 16 h and at no
other timepoint, but the effect estimate is largest at rest and merely noisier
there, so the window is a statement about precision rather than amplitude (Fig. 6). The
axis the gene marks is enriched 21.3-fold for glycolysis and is not specific to
TPI1, which is what a leave-one-out reconstruction with the gene removed from
both the defining score and the enrichment family shows (Supplementary S31;
Fig. 7). In melanoma single-cell data TPI1 is 2.75 log2 units higher in malignant
cells than in the CD4+ T cells the instrument came from, in 16 of 16 patients; the
compartment attribution across cell-type, within-patient paired, spatial and
bulk data is Fig. 8.
Splitting cells by glycolytic score splits them by lineage purity, and the
helper and regulatory differences collapse once lineage composition is matched
at cell level. In patients, both arms are significant in the discovery cohort
and the result appears in each of the three cohorts without replicating in any
two of them; the locked signature is Supplementary S17 and the inference,
including all three treatments of repeated patients, Supplementary S19 (Fig. 9).

**Nothing on that ladder makes TPI1 a target, and we do not present it as one.**
Its value here is the price it puts on each rung of a ladder the identity
predicts is not climbable, and every one of those prices was paid before we
understood why.

### Part 4 — These checks are rarely reported

From 469 screened records, 209 were eligible and 154 had an open-access full
text, of which 152 were scorable. In that sample, **not one paper in a random 30%
subsample compared its significant signal against previously reported loci for
its own outcome trait (0% [0–3.4]), and an estimated 7.1% [2.5–16.1] report how
the candidate list depends on the outcome GWAS used**. About 59% performed
colocalisation and 47% used SMR or HEIDI. Three studies using this framework
closely enough to compare item by item are set beside ours in Supplementary
S16.

Those are corrected figures. Automated matching returned 7.9% and 35.5%; **two
coders independently scored a random 30% subsample (46 papers, 92 judgements),
blind to each other and to the automated score**, and its precision proved to be
0.00 and 0.20. Most hits were spurious — three of five on the first criterion lay
inside reference lists, where a cited GWAS title containing "novel loci" matches,
and twelve of fifteen on the second were generic statements about statistical
power. **Neither of our published numbers survived**: 7.9% was too generous, and
after adjudication no paper in the subsample performed that check at all, while
0.7% for the second criterion was too strict and falls below the corrected
interval. We had reported the permissive count for one criterion and the strict
count for the other.

Agreement was 97.8% on each criterion with one disagreement each; Cohen's κ was
0.79 for list stability and **0.00 for locus attribution — a prevalence artefact,
since with one positive in 46 the expected agreement equals the observed**, which
is why the raw agreement is given beside it. Both disagreements were adjudicated
jointly and both resolved against the first coder; κ is computed before that.
Recall was bounded from above rather than measured, and a full-text audit
measures reporting, not practice (Methods).

How far the identity reaches into this literature is a question about the
literature rather than about the algebra, and our first attempt to answer it
failed. Scanning the 154 cached full texts — two more than the 152 scorable ones,
because scorability required the five pre-fixed criteria to be assessable and the
scan needs only the text — 71 contain a phrase asserting a single-variant
instrument and 130 name a multi-instrument estimator. **Neither count survives
inspection.** 23 of the 71 never mention a cis-eQTL anywhere; requiring the
phrase to fall within 500 characters of one leaves 13, and within the same
sentence, 10. **A count that falls from 71 to 10 under a minimal proximity
requirement is a phrase-hit tally of unknown direction, not a bound in either
direction, and we withdraw it as one.**

We therefore coded the quantity the question needs — the estimator that produced
each paper's primary cis nomination list — by hand, as a sixth criterion on the
same fixed-seed 46-paper subsample used for the others, by two coders
independently under rules fixed before any paper was read. Agreement was 71.7%
and Cohen's κ 0.534, both before adjudication; the thirteen disagreements were
resolved jointly, each with a written justification citing the evidence.
**Those justifications were destroyed by a defect in our own scoring script
and rewritten by the two coders afterwards; the codes were preserved
throughout and none changed, so the distribution is unaffected, but the
justifications are reconstructions and we report them as such.** **Of 45
ascertainable papers, 10 rest on a single-variant Wald ratio, 28 on a
multi-instrument estimator and 7 on mixed pipelines; one could not be
ascertained and is not redistributed.** The identity therefore governs
**10 of the 45 ascertainable papers in this subsample (22.2%; Wilson 95% CI
12.5–36.3%)**, or at most 37.8% if every mixed pipeline is counted as partly
governed. **This is a statement about these 45 papers and the sampling frame
that produced them, not about current practice as a whole**; within that
frame the design we analyse is a minority, and we state it as the limit on
the claim rather than leaving the scope open (Supplementary S51). The scan itself is reported as a corpus-discovery tool
rather than as a result.

---

## Discussion

The results have a single mechanism, and it is algebraic rather than empirical.
With one instrument and first-order standard errors the Wald statistic is
|β_out|/se_out, so a candidate list from this design is the outcome GWAS
restricted to lead cis-eQTLs and thresholded at whatever the multiple-testing
burden allows — here 2.1 × 10⁻⁴. One consequence is algebraic and the rest are not, and the
distinction matters. What follows algebraically is narrow: signal must land where
the outcome GWAS has signal, and profile specificity cannot appear in a P value,
because the P value contains no exposure term. Everything else below is an
empirical observation that the identity makes likely rather than necessary.
Whether the outcome's strongest signals among lead cis-eQTLs sit on previously
published loci depends on where cis-eQTLs are and on the disease's genetic
architecture, not on the algebra — here three of our eight significant loci reach
5 × 10⁻⁸ unaided, all three are on the reference list, and removing them leaves
1.98-fold, P = 0.41 (2.07-fold, P = 0.40 with the background restricted to
sub-threshold loci as well). The attribution recurred in every further disease we
could test, which the identity did not require; we nonetheless claim it as
observed in two tumours rather than in tumours generally, for the reasons in
the limitations. Colocalisation was poor in this
application, 2 of 284 above PP.H4 0.8; the identity does not require that either,
since a variant that drives both expression and disease can colocalise well.

We think this reframing is the useful contribution, and we state plainly what it
is not. The reduction itself is elementary and not new. It does not apply to
multi-instrument estimators, where combining instruments breaks the identity —
our relaxed IVW and weighted-median analyses are outside it. How much of the
literature that exempts is now measured rather than assumed: of 45 ascertainable
papers in a two-coder subsample, 10 (22.2%; Wilson 95% CI 12.5–36.3%) produce
their primary cis nomination list from a single-variant Wald ratio
(Supplementary S51). **Within that sample the design this paper analyses is a
minority**, and we do not extrapolate the fraction beyond it. Whether the argument weakens gracefully as instruments
are added is not something we tested, and we do not assert it. It does
not make MR effect estimates wrong; the sign and the scale genuinely come from the exposure,
and only the significance claim does not. And it does not show that the
sub-threshold part of a nomination list is spurious: at five loci we are
underpowered to say anything about it, which is a different statement from
showing there is nothing there.

Everything in this framework that depends on the outcome GWAS proved unstable,
and the instability has a shape: signal lands on loci the outcome already knows
about, the novel part of the list turns over completely between outcomes, and a
larger, compositionally different meta-analysis moved MR and colocalisation in
opposite directions. The pattern held when we changed the disease, when we
changed the exposure resource by two orders of magnitude, and in five of the six
crossed cells; in the sixth the mismatched-locus control itself enriched and we
report the cell as void rather than as a result. That control is the load-bearing
part of the sentence and it is scale-dependent, most so exactly where the outcome
pair shares biology — so the claim it licenses is scale-dependent too. Scoring
every cell against six unrelated lists rather than one, a post-hoc diagnostic
and not a registered test, turns that into a design lesson for any study using
this control: which disease is nominated as the mismatch determines the verdict,
and a control with one arbitrary comparator measures the comparator. The full
argument for both is in the Results, where the void cell and the comparator
sweep are reported at length. What did not depend on the outcome was a smaller
and more durable statement: in which cell state a pathway's regulation can be
measured precisely enough to yield an instrument at all.

Of eleven target-level tests that returned a verdict, seven overturned the claim
under examination, three were inconclusive and one was significant on discovery
but not confirmed on transfer; five genes were successively designated lead
candidate under a ranking criterion that was never fixed in advance, and the
first four designations were overturned. That is this project's own history and
not an estimate of how often the framework fails, and the itemised record — every
attempt, stopping rule and withdrawal, the candidate-selection timeline, the
seventeen pre-registrations and the technical account of two processing errors of
ours — is Supplementary S12.

Eight checks follow from the failures above. Each is stated in full, with the
numbers behind it, in Supplementary S49, and all nine — these eight plus the
decomposition of S42 — are implemented in the `cqtna` package.

**(i)** Annotate the signal against the outcome's own known loci, counting.
**(ii)** Require colocalisation, preferably with explicit multiple-signal modelling.
**(iii)** Do not assume that more outcome power will resolve the disagreement.
**(iv)** Report the physical distance between GWAS and eQTL peaks alongside.
**(v)** Compute the nominated gene's cell-type expression ratio in an annotated
single-cell dataset for the tissue, before any functional interpretation.
**(vi)** Report instrument availability at each of the three levels it passes.
**(vii)** When splitting cells by a score, match on lineage composition as well as.
**(viii)** Before comparing candidate lists across releases of the same GWAS.

The worked example is best read as a ladder rather than a verdict, the same
layered reporting Howe et al. use when a cellular model proves only partly
transportable [8]. Of seven claims about TPI1, two are supported — dynamic
*expression* across activation, and membership of a definable CD4⁺ metabolic
state — and neither is a claim that the gene is a target; the other five are
not supported, unstable, or answered no. **What MR supplies here is a
genetically tractable entry point into a pathway, and the audit establishes
that the gene marks an expression state — not that the state is genetically
driven, and not that the gene is a target.** Instrument availability and
biological importance are separable: the strongest enzyme in the functional
data carries no instrument, and the gene that does carry one is principally
expressed by tumour. The claim-by-claim table with its numbers is
Supplementary S47; the complete record of which claims were raised, tested and
withdrawn, with the selection denominators, is Supplementary S12, and the
seventeen pre-registration documents are enumerated in Methods.

**Relation to existing guidance.** These eight add to, and do not replace,
current cis-MR practice, which already stresses that cis analyses must be
tailored to local biology, that an expression biomarker is not an intervention,
that co-regulation of neighbouring genes can act as horizontal pleiotropy, and
that MR is one strand of triangulation [9,10]. Tambets et al. showed with
approximate ground truth that neighbouring-gene co-regulation routinely produces
colocalisation evidence for several genes at once, so eQTL data suit candidate
*generation* and evidence *composition* better than gene *attribution* [7]; our
MC1R region is the same phenomenon at higher effect size, and TPI1 at chr12p13
is the same problem unresolved. Reales et al., auditing over a million
colocalisation tests, likewise found gene assignment shifting with platform
coverage, cell context and resource size [11]. What we add is specific to the
nomination step and to the outcome side: locus attribution against the outcome's
own known loci, sensitivity to the outcome release, candidate-list stability,
layered instrument visibility, and cell-compartment attribution.

**Limitations.**

*Limitations — design.* Three diseases and two exposure resources, the second
hepatocellular power level drawn from the same resource as the first, show the
pattern is not peculiar to melanoma or to one small eQTL dataset — not that it
holds for any disease or resource. One of the six cells is void on its own
mismatched control, one of the remaining five reaches neither P < 0.05 nor the
density-matched null, and the surviving non-cancer cell does so at the registered
window only, its control failing at six of seven settings and harder as the
window narrows. We claim the attribution as observed in two tumours, not in
tumours generally, and treat the autoimmune cell as exploratory. Three diseases
cannot estimate generality, and the six comparator lists are all cancers:
unrelated to rheumatoid arthritis, but not a random sample of diseases.

*Limitations — data.* Both resources carrying the attribution finding are
European-ancestry and the second is whole blood, so a positive result there
cannot be read back as CD4-specific; the three further resources we tried, which
include Peruvian and African-ancestry donors, were underpowered or failed their
control. **That limit is structural.** In the eQTL Catalogue [12],
multi-condition stimulation resources at adequate donor numbers exist only for
myeloid cells; for CD4⁺ T cells the only one is used here. **The dynamic axis of
the framework this paper audits rests on a single resource**, which bounds any
audit of it, ours included. The patient analyses rest on three small cohorts of
which the strongest is the discovery cohort, and the genotype × pseudotime test
covers lead variants under linear and quadratic models, so a non-monotonic effect
at a non-lead variant is not excluded.

*Limitations — estimator.* The identity governs the single-variant Wald design
and nothing else. Our relaxed IVW and weighted-median analyses lie outside it and
are reported as concordance rather than discovery, and the seven-way benchmark is
a boundary check rather than seven independent tests. How much of this literature
runs the design the identity governs is measured, not assumed, and it is a
minority (Supplementary S51).

*Limitations — diagnostic.* The mismatched-list control is informative only
inside a bounded power range: as the significant set grows towards the background
the fold approaches 1 and unrelated lists begin to enrich too. Scored against
their own lists, the own-list fold falls from 4.96 at eight significant loci to
2.04 at twenty-eight, and at the upper end comparators enrich together — for
prostate and breast, whose lists share 43% of their flagged loci,
own-versus-mismatch cannot discriminate by construction. A clean control at low
power and a failed one at high power are not comparable, and melanoma's clean
controls are partly a property of its having only eight significant loci. The
power-recovery curves rest on reference lists of four and six genes and compare
against a full-power list that is itself unstable, so they measure agreement
between two imperfect lists rather than recovery of truth.

*Limitations — reproducibility, and what needs new data.* The external validation
protocol is adopted and frozen (Supplementary S40) but not registered, and no
qualifying cohort has been secured, so it carries no realised power and stands as
a commitment. The analyses ran on one machine; the deposited container fixes the
software versions and not the machine, and which numbers survive that change is
recorded there rather than assumed. Two questions cannot be settled here at all:
whether the patient stratification is real, which needs a same-disease,
same-regimen, CD4-resolved cohort that does not exist publicly, and whether
the state has immunological consequences, which needs perturbation. A
genome-scale CD4⁺ T perturbation map has since been published, and we
pre-registered a lookup of the nominated gene in it: the gene is not
ascertainable there, because its expression is outside the assay's probe
panel, so the question stands open rather than answered (Supplementary S52).

**Conclusions.** A candidate list produced by a single-variant cis-eQTL MR
design at current outcome power should be reported as the set of genes that passed screening under these
conditions, not as targets — and because the fragile part of such a list is
precisely its novel-locus part, the headline of a study of this kind is the part
least likely to replicate. The eight checks above cost little and would have
changed what this analysis reported at nearly every stage.

---

## Figures

**Fig. 1 | Locus attribution under both outcomes, by bounded locus.**
**a**, One point per bounded locus, both outcomes: a mirrored Manhattan plot in
which each locus contributes its strongest record, FinnGen R12 above and the
meta-analysis below, coloured by locus class. The MC1R region gives the
strongest association in the study, P = 4×10⁻³⁷, and PARP1 reproduces its
published direction; these are the end-to-end positive controls. **b**,
Attribution, by locus: for each outcome, the share of significant loci that are
known loci against that outcome's own background share, 10.1% for the meta
outcome. Loci are the non-recursive fixed-anchor partition frozen in
Supplementary S37 — 1 Mb claimed from each anchor — not the single-linkage rule
used in the source literature, which chains a dense resource into blocks of tens
of megabases. Counting is by bounded locus throughout; the record-level version
is Supplementary S10.

**Fig. 2 | Generality: the attribution under one change at a time.**
**a**, One resource, five nested power levels, and a transfer test: FinnGen R8
to R12 carry an identical melanoma endpoint at rising case numbers, with
everything else held fixed. R13 is separated by a rule and coloured apart
because it is not a sixth power level — the endpoint definition and the control
exclusions both change — and must not be read as a continuation of the curve.
**b**, Both axes crossed: three diseases by two exposure resources. The second
disease and the non-cancer outcome are two columns of this grid and the second
exposure resource is its second row. The lower-powered hepatocellular level sits
below a separating rule because it is a power sensitivity, not a seventh and
eighth cell. **c**, The same cells scored against the wrong disease: each cell
rescored with another disease's known-locus list. Without this panel, **b** can
only show that significant loci fall on some loci, not that they fall on the
loci of their own disease. All panels are computed on the frozen partition
(Supplementary S37).

**Fig. 3 | Power and list stability by locus class, and the same stratification
against full-power |z|.**
**a**, Down-sampling reproduces a separately observed round, which is the
positive control for the whole simulation. **b**, The same curve in five
cancers. **c**, The reproducible part of the candidate list is the part that is
not a discovery: recovery against power, split by locus class. **d**,
Novel-locus candidates sit against the detection threshold by construction —
their full-power |z| all lie between 3.73 and 4.55, while known-locus candidates
reach 15.99. **e**, Conditioning on effect size absorbs most of the class gap:
raw gaps of 47.6 percentage points by bounded locus and 63.0 by gene fall to
+5.2 and +6.4 after matching on |z|. The locus-level interval crosses zero, so
**the residual is not claimed as an effect of the class label**; only one known
locus lies inside the novel |z| range, and that failure of common support is why.

**Fig. 4 | Colocalisation versus SMR/HEIDI, with in-sample fine-mapping.**
**a**, HEIDI fails to reject the LD-confounded signals colocalisation
identifies: of the 291 records colocalisation assigns to distinct causal
variants, HEIDI does not reject homogeneity for 253 (86.9%). **b**, In-sample
fine-mapping, and what proxy LD does to it. FinnGen's published in-sample LD
gives three high-purity credible sets at MC1R, so the region genuinely carries
several independent signals — the configuration in which a neighbouring gene's
eQTL is tagged by one of them without sharing it. The same region run on
1000 Genomes EUR proxy LD returns nineteen, which is over-splitting.
**Proxy-LD credible-set counts are therefore one-directional: exactly one set is
conservative evidence of a single signal, but several sets are not evidence of
several.** In-sample counts exist for MC1R only and are not extrapolated to the
other regions.

**Fig. 5 | Instrument availability across the glycolytic pathway, at three
levels.**
**a**, Strongest cis-eQTL per gene and activation profile: 28 glycolytic genes
by 8 profiles, with the genome-wide threshold marked and genes absent from the
resource shown separately. **b**, The same pathway counted at three levels,
naming which gene drops out where and why. The three counts are not the same
kind of statement and the panel keeps them apart: the first depends only on the
pathway and the exposure data and is unchanged by the outcome; the second adds
harmonisation against one specific outcome GWAS, so a gene lost there is not a
gene without an instrument; only the third involves the outcome's own
significance.

**Fig. 6 | The activation window in which TPI1 is instrumentable: effect size
against precision across the eight profiles.**
**a**, The same variant in all eight profiles, β with 95% confidence intervals.
**b**, Effect size against precision. On this plane genome-wide significance is
a ray through the origin, |β| = 5.45 × SE. The resting profile sits highest on
the effect-size axis and still falls below the ray: **the largest effect and the
only usable instrument are not the same timepoint.** This is a statement about
where the effect is measurable, not about where it is regulated, and it neither
contradicts nor repeats the published genotype-by-pseudotime interaction test,
which asks a different question and is negative at all three lead variants [4].

**Fig. 7 | The CD4⁺ metabolic axis and its chromatin signature.**
**a**, What the axis is made of, with no pre-specified modules: family
composition at the top of the residualised axis, glycolysis enriched 21.3-fold.
**b**, Eleven biological modules against a composition-matched null, the
threshold being the 95th percentile of |SMD| across fifty random modules. The
activation module remained within the matched-null distribution. Two modules
exceeded the threshold: OXPHOS, and the proliferation module, which was
specified in advance as a sensitivity control (Methods; Supplementary S47).
**c**, A concordant chromatin signature: motif enrichment in axis-high versus
background peaks. **d**, AP-1 enrichment persists when peaks
are restricted to those that do not change with activation. Panels **c** and
**d** characterise the axis, and are not the withdrawn claim that the instrument
acts by disrupting an AP-1 motif, which was tested separately and failed.

**Fig. 8 | Compartment attribution.**
**a**, TPI1 by cell type, with detection rates. **b**, Paired within patients:
malignant versus CD4⁺ cells, unmatched, depth-matched, and in an independent
cohort side by side. This is the primary evidence in the figure. **c**,
Pre-specified positive controls. **d**, Spatial context, as a cross-section
meta-analysis. **e**, Bulk survival is composition, in TCGA-SKCM as distributed by Xena
[13,14]. Spatial and bulk data are
context here rather than evidence; the cell-level comparison in **b** is what
the attribution rests on.

**Fig. 9 | Patients across three cohorts.**
**a**, The same pre-specified test in three cohorts: effect size and one-sided
permutation P for both arms. Both arms are significant in the discovery cohort
and neither is confirmed twice, and the two failures are mirror images — the
pre-treatment arm matches in magnitude in the same-disease cohort but is
underpowered there and absent in the cross-disease one, while the post-treatment
arm inverts in the same-disease cohort and returns at full magnitude in the
cross-disease one. **b**, What the cross-disease design could have shown at all.
Its pre-treatment arm compares four responders with two non-responders, so no
effect size whatever can reach P < 0.05. That floor was computed and registered
before any expression value was read (Supplementary S21), and a test that could
not have been significant is not the same thing as a test that was run and was
not significant.

Figure numbering follows this version's reading order and is independent of the
numbering used in the full reference document. The self-administered attribution
check is reported in the text and tabulated in Supplementary S27; it carries no
main figure.

## Methods

### Exposure data

Cis-eQTL summary statistics came from a published activation time-course of
primary human CD4⁺ T cells [15], comprising eight profiles — naive and memory
cells at 0 h, 16 h, 40 h and 5 d — from 85–100 donors each (naive 99/99/89/85;
memory 100/95/89/90). Per-profile sample sizes were recovered from
`ma_count/(2·MAF)` and were internally consistent for 100% of variants in all
eight profiles. Expression had been inverse-normal transformed, so exposure units
are standard deviations. Coordinates are GRCh38 and were used as distributed; no
liftover was performed. The second exposure resource was the eQTLGen Consortium
whole-blood cis-eQTL release (n = 31,684) [16].

### Outcome GWAS

The primary outcome was a fixed-effect inverse-variance meta-analysis of FinnGen
R12 `C3_MELANOMA_SKIN_EXALLC` (5,753 cases / 378,749 controls; the OpenGWAS
`finn-b-` snapshot was not used, being an R5 extract with 98 cases) and Rashkin
et al. (GCST90011809, 6,777 cases) [17,18], totalling 12,530 cases and 789,099
controls (N = 801,629, case fraction 0.015631). Rashkin reports odds ratios and P
values without standard errors, so standard errors were reconstructed as
se = |log OR|/|Φ⁻¹(P/2)|. For a two-study meta-analysis Cochran's Q has one
degree of freedom, so the upper tail is erfc(√(Q/2)); the exp(−Q/2) form used
initially gave a significant-Q rate of 1.59% against 5.33% for the correct form.
1,740 variants underflowed to P = 0, for which log-space `mlogp` was computed and
P floored at 10⁻³⁰⁰. Landi et al. (GCST010304) supplied only 76 loci and was used
solely as the known-locus reference [19]. Further outcomes were FinnGen R12 lung,
colorectal, pancreatic, breast and prostate cancer; FinnGen R12
`C3_HEPATOCELLU_CARC_EXALLC` (947 cases) [20] and GCST90809296, the European arm of an
eleven-cohort hepatocellular meta-analysis (3,748 cases / 1,861,536 controls)
[21]; and FinnGen R13 `M13_RHEUMA` (16,775 cases / 308,847 controls). Effective
sample sizes, 4/(1/N_case + 1/N_control), are 49,337 for the melanoma meta,
14,962 and 3,779 for the two hepatocellular levels, and 63,643 for rheumatoid
arthritis.

### Instrument selection, harmonisation and MR

Instruments were the top cis-eQTL per gene × profile at P < 5×10⁻⁸, matched to
the outcome by chr:pos with allele consistency; 26 variants were discarded for
ambiguous indel representation. The strict set is 3,556 records over 2,126 unique
variants for the meta round, and 3,579 over 2,142 for the single-outcome round.
The minimum F statistic was 36.1 in either strict set and 22.2 across the
harmonised set before instrument selection, so the F > 10 filter was never
binding at any stage. Each exposure carries one instrument, so the primary
estimator is the Wald ratio with first-order delta-method standard errors. The
analytic consequence is stated because it constrains interpretation: with
b = β_out/β_exp and SE = se_out/|β_exp|, the statistic z = β_out/se_out does not
involve the exposure, so a variant returns the same P value in every profile in
which it is the lead eQTL. Profile specificity is therefore reported as the
profile containing the instrument, never as a difference in P value.

A relaxed set (F ≥ 5, P < 0.05, LD-clumped with PLINK 2 at
`--clump-r2 0.1 --clump-kb 1000` against 525 unrelated 1000 Genomes GRCh38 EUR
samples [22,23]) supported IVW, multiplicative-random-effects IVW and weighted
median sensitivity analyses [24]. Weighted mode was evaluated and excluded: at a
median of four instruments per exposure it returned 1 significant result of 226
testable exposures. Because clumping at r² < 0.1 retains correlation while IVW
assumes independence, IVW significance is inflated (157 exposures at FDR < 0.05
against 21 in the strict set), so sensitivity results are reported as concordance
between IVW and weighted median rather than as independent discovery. Multiple
testing used Benjamini–Hochberg FDR within each analysis family, the family fixed
before the analysis ran. Steiger filtering used the recovered per-profile sample
sizes, SD exposure units, log-odds outcome units and outcome prevalence 0.014962;
all records had the correct direction. TwoSampleMR's R² formula for SD units is
unbounded and produced values above 1 (maximum 1.027), so the bounded form
R² = F/(F + N − 2) (range 0.185–0.933) is additionally reported.

**Unit of inference.** The testing family is gene × profile records, and the
record level is primary throughout, because it reproduces the pipeline under
audit: the registered down-sampling predictions, the release trajectory, both
generalisations and the grid were computed on it. Records are not independent —
3,556 carry 2,126 unique variants and 1,195 unique genes — so the list was
recomputed under seven alternatives as a sensitivity analysis: variant, gene and
bounded locus, each collapsed by minimum-p and by Simes
combination, plus a two-stage hierarchical procedure selecting genes by Simes and
then records within selected genes. All seven return a superset of the
record-level gene list. The audit's own conclusions are stated in independent
loci.

### Colocalisation, SMR and HEIDI

Colocalisation used `coloc.abf` with default priors over the full cis window
[25], the eQTL as a quantitative trait and the disease as case-control, matching
on a chr:pos:alleles key. The published criterion PP.H4/(PP.H3+PP.H4) > 0.7 is
too permissive here, since the ratio can rest on two vanishingly small
posteriors, so PP.H3+PP.H4 > 0.5 was additionally required. Every result carries
two diagnostics: the number of variants in the region, and the physical distance
between the GWAS and eQTL peaks. Sensitivity across four windows (±100 kb to the
full cis region) and four values of the shared-causal-variant prior (p₁₂ from
10⁻⁶ to 5×10⁻⁵) is reported for every locus. Because `coloc.abf` permits at most
one causal variant per trait per region, the assumption was tested on the outcome
side using FinnGen's published in-sample SuSiE fine-mapping and `susie_rss` run
on 1 Mb windows with a proxy LD matrix from the same 525-sample panel [22,26], with
`estimate_s_rss` as a mismatch diagnostic. The eQTL side was not fine-mapped: at
85–100 donors without in-sample LD, SuSiE is not reliable there, and this is a
stated boundary. Our proxy-LD procedure recovers 19 credible sets at MC1R against
3 from in-sample fine-mapping, so its counts are read in one direction only —
exactly one credible set is conservative evidence of a single signal, several is
not evidence of multiple signals.

SMR v1.3.1 [27] was run against the same reference panel (94,979 within-cis
variants) with GWAS statistics trimmed to cis windows and de-duplicated for
multi-allelic sites. The `--diff-freq` QC removed 0.45% of variants at the
default threshold, which was not relaxed. The number of variants entering each
HEIDI test is reported with the result, as HEIDI's power depends on it (range
8–20 here).

### Locus annotation and enrichment

The melanoma reference set was the union of lead variants from the three Landi
et al. 2020 GWAS [19] retrieved from the GWAS Catalog REST API [28] and mapped to
GRCh38: 157 loci. Hepatocellular loci were lead variants for EFO_0000182 and
liver cancer at P < 5×10⁻⁸ together with those tabulated in the outcome
publication — 83 rsIDs, 73 placeable on GRCh38. The rheumatoid arthritis
reference was the 87 lead rsIDs of Okada et al. 2014 [29] (GCST002318), published
before FinnGen existed and from cohorts that exclude it; 83 were placeable.
Loci are a **non-recursive fixed-anchor partition** of instrument
positions within a chromosome — bounded rather than independent, since the
partition guarantees only that two loci share no variant: the first unassigned variant becomes an anchor and
claims every variant within 1 Mb of it, then the next unassigned variant becomes
the next anchor. Anchors are taken by position and never by significance, because
the same partition supplies the numerator and the denominator, and a locus
therefore spans at most 1 Mb whatever the variant density. A locus counts as
known if **any** of its instruments lies within 1 Mb of a reference lead; the
same convention scores background and significant loci, since a background locus
has no significant record to be scored by. The partition rule, the 1 Mb main
window, the 500 kb sensitivity and this convention were fixed together in
Supplementary S37 before the last four cells of the grid were computed.

The single-linkage rule used in the source literature — variants joining
transitively at 1 Mb — is retained only where we reproduce published numbers, and
carries no inference here. It is not interchangeable with the above: on the
whole-blood resource it chains a chromosome arm into a single 30.8 Mb "locus" of
588 records, and a block that wide contains a known lead SNP for almost any
disease. Under that rule the pre-registered mismatched-list control fails on the
whole-blood melanoma cell, which is what prompted the change; under the
fixed-anchor rule the widest significant locus is bounded at 1 Mb by construction
and that control is clean (Supplementary S37 §2). Enrichment among
FDR-significant loci was tested one-sided by Fisher's exact test against the
background proportion over all instrument loci. Every reference list was built
and written to disk before the MR step it scores was run.

Three controls bound this test. A **mismatched-list** control scores each
disease's significant loci against another disease's reference. Because the
verdict of that control depends on which disease is nominated, it is accompanied
by a **multi-list** diagnostic (Supplementary S39): every cell is scored against
every reference list in the study that has at least 30 placeable lead SNPs and
belongs to a different disease, six per cell, and we report how many also enrich
together with the ratio of the outcome's own fold to the largest fold from any
comparator. **This diagnostic was specified after the results it comments on, and
its criterion — that the outcome's own list exceed every comparator — is weaker
than the registered rule that voids a cell whose mismatched list enriches at all.
It is reported as exploratory and is not used to upgrade any cell.** The 30-lead
floor is the one already used in Supplementary S34 and no new lists were
obtained; the breast and lung lists were, however, rebuilt from the same
GWAS Catalog rsIDs after we found that the original coordinate lookup had
silently dropped whole request batches, leaving 6 and 0 placeable leads where
740 and 271 were available (Supplementary S41). Under the repaired panel the
melanoma CD4⁺ cell acquires one nominally enriching comparator, breast at
2.15-fold (P = 0.0495), so its multi-list verdict moves from no enriching
comparator to one; the cell's margin over its best rival widens from 1.82- to
2.31-fold, because the published prostate comparator had been computed from
35 of its 1,160 leads. Under Holm correction across the six comparators the
breast result does not reach significance, so this verdict depends on which
multiplicity rule is applied and we report both. A **strength-and-
frequency** permutation matches significant to background loci on eQTL-p decile
crossed with outcome allele-frequency quintile, with 10,000 stratified resamples
and a one-sided empirical P. A **density** permutation follows the estimand
in Supplementary S38, whose matching rule, coverage requirement and sampling
scheme were fixed in a document first committed before the run; the choice of
1.00 as the reported tolerance was recorded in the same commit as the results and
is therefore an analysis choice rather than an independently timestamped one.
Background loci are matched to
significant loci on record count, physical span and unique-gene count
simultaneously, drawn without replacement from a pool that excludes every
significant locus, over 10,000 replicates with a one-sided empirical P. The
matching tolerance is fixed at 1.00 and the full tolerance scan is reported
alongside; a cell whose matched coverage falls short of 100% returns no P value
rather than a partially matched one, and four cells do so at the tightest
tolerance. An earlier quantile-stratified version (Supplementary S30) is retained
but is a different test and is not pooled with this one. All three controls are
computed on the same locus partition as the estimates they bound.

**Circularity of the hepatocellular reference.** Ten of its 73 loci came from the
higher-powered outcome publication's own table. Provenance was resolved through
the Catalog REST API by retrieving that study's reported associations and
querying each locus for the number of other studies reporting it: seven have
substantial independent support and three rest on that publication alone.
Attribution was recomputed against the published list, the list without those
three (primary), and the list without all ten as an over-corrected bound. Note
that removing known loci shrinks the background and can only raise the
enrichment, so the over-corrected version is not a stricter test.

### Power trajectory, transfer test and down-sampling

Five FinnGen releases carrying an identical melanoma endpoint (R8–R12) were
analysed with exposure data, instrument set, allele orientation, harmonisation,
estimator, testing family and locus-merging rule held fixed, so the systematic
difference between runs is case number. The primary analysis was restricted to
variants present in every release. Predictions were generated by down-sampling
R12 over 2,000 replicates and registered before any earlier release was examined.
Because releases are nested, agreement between them exceeds that expected between
independent studies of equal size, so the design is conservative with respect to
instability.

R13 replaced `C3_MELANOMA_SKIN_EXALLC` with `C3_MELANOMA_SKIN_WIDE` (6,226
cases). Endpoint definitions were compared code by code in Risteys rather than
read from the release manifest, whose one-line description ("including Hilmo")
had led us to the wrong conclusion at first: the two endpoints select cases with
the same codes and differ in the rule excluding cancers from the controls. R13
was therefore excluded from the trajectory, whose comparability gate requires
verbatim-identical endpoints, and analysed as a transfer test with everything but
the outcome file held fixed. Three process controls were required before any R13
number was interpreted: exact reproduction of the R12 candidate list through the
new code path, presence of the large-effect MC1R-region genes, and instrument
coverage of at least 95% of the R12 set.

Power was reduced in summary space, individual-level data being unavailable: for
a target effective size, se_sim = se_obs·√(N_obs/N_sim) and
beta_sim = beta_obs + ε with ε ~ N(0, se_sim² − se_obs²). Wald ratios and BH-FDR
were recomputed over 200 replicates per grid point and the FDR < 0.05 list
compared with the full-power list by Jaccard index, sensitivity and precision.
Calibration was pre-specified as the condition for reporting any curve: reducing
the meta outcome to FinnGen's effective size reproduces the independently
observed FinnGen round (10.2 simulated against 10 observed discoveries; simulated
Jaccard 0.52 [0.36, 0.73] containing the observed 0.455). We call this empirical
rather than external calibration, because FinnGen contributes to the outcome
being down-sampled. Two attempts to extrapolate above observed power were
abandoned — a flat prior overpredicted discoveries at observed power 25-fold and
global empirical-Bayes shrinkage underpredicted them 100-fold — both caught by
the pre-specified requirement that any prior reproduce the observed discovery
count. **All curves are reported only up to observed power and must not be
extrapolated.**

**Effect-size matching.** The down-sampling model takes no class label as input,
so the known-versus-novel differential can only follow from the two classes' |z|
distributions. Each unit's recovery frequency at half power was regressed on
log|z| with and without the class label, and novel units were matched to known
units within a pre-registered caliper of 0.20 on log|z|, with replacement,
nearest first, the residual bootstrapped over matched pairs. Common support was
recorded, since the matched comparison is interpretable only where the |z|
distributions overlap. The locus level was pre-specified as primary. A process
control required the reimplementation to reproduce the published recovery curve
bitwise before any matched result was read; this failed on first run because the
published script draws from one generator consumed across all power points while
the reimplementation re-seeded per point, and was resolved by matching the draw
sequence rather than by loosening the tolerance.

### Self-administered attribution check

Of the 30 loci reaching FDR < 0.05 under the whole-blood exposure, 20 are known
loci, and those 20 are the domain of this check. A locus's nomination is the set
of FDR-significant genes there whose own instrument lies within 1 Mb of a known
lead. Melanoma loci at which a causal gene is generally accepted were listed, and
the nomination at each was compared against that gene under two rules fixed
before the comparison: *lenient*, the accepted gene appears anywhere in the
nomination, and *strict*, the nomination is that gene alone. Ten of the 20 had an
accepted gene. The list was assembled by us and fixed before the comparison was
run, but was **not** pre-registered, and ten loci is a small denominator; the
resulting ratios are reported as an illustration and not as an error rate.

### Single-cell, spatial and chromatin analyses

Full processing details for the ICB cohorts (GSE120575, GSE115978, GSE72056,
GSE78220, GSE91061, GSE235863) [30–35], the three spatial datasets [36–38], the purified
CD4⁺ multiome (GSE282266) [39] and the peripheral-blood variance decomposition
(GSE199994) [40] are given in Supplementary S12; the points that changed a result
are summarised here. Donors or samples, never cells, are the statistical unit for
inference. Module scores used Seurat `AddModuleScore` [41] with expression-binned
control genes, and every module used as evidence independent of the nominated
gene was recomputed with that gene removed. Because splitting cells by a score
also splits them by lineage purity, every within-population split was controlled
at the **cell** level — lineage-marker-free module versions, residualisation on a
CD8ness score, and matching of high- and low-score cells on CD8ness and
sequencing depth, with residual standardised mean differences reported for both.

The metabolic axis was defined by residualising a 16-gene glycolysis module on an
eleven-gene activation module and on sequencing depth, with a pre-specified kill
criterion — if R²(glycolysis ~ activation + depth) > 0.70 the axis is not
separable from activation and the analysis stops. Observed R² = 0.024. Cells were
matched on activation decile × depth quintile before splitting on the residual,
and only genes concordant in direction across all four sets were carried forward.
A nuclear-retention index (MALAT1 + NEAT1 fraction) carried its own kill
criterion, and module detection rates were checked before any module-level null
was interpreted; three pre-specified immune modules were withdrawn on that basis
as unmeasurable rather than negative. The axis was originally circular with
respect to the nominated gene, which entered both the defining score and the
enrichment family, and was rebuilt with that gene removed from both, everything
else held fixed, leaving its rank a held-out quantity.

**Replication of the axis.** GSE166188 (DOGMA-seq, STIM arm, 16 h) was analysed
with the construction held bitwise identical to the discovery run and only the
dataset changed, pre-registered before it was run (Supplementary S32). CD4 cells
were gated on CLR-normalised surface protein (CD4 above median, CD8 at or below),
not on transcript. Its two arms are two lysis buffers on one batch of material —
technical replicates, not two donors — so this replicates across study, platform
and lineage calling and not across donors. Two dimension conventions differ
between the protein and RNA matrices in this resource, and dimension assertions
were added after an earlier script indexed the protein matrix by the RNA
convention; the protein file also retains unfiltered empty droplets, so gating
requires intersecting with the RNA barcodes first.

Chromatin peaks differ between samples, so cross-sample comparison was performed
on merged genomic intervals (199,740 consensus intervals) with the global index
preserved. Motif enrichment used JASPAR2020 CORE with `motifmatchr` against
hg38 [42,43], comparing foreground with background peaks matched on GC content
and log accessibility. Because axis and activation log₂ fold-changes correlate at
r = 0.647 at the chromatin level despite near-orthogonality at the RNA level, the
analysis was repeated restricted to activation-invariant peaks; both versions are
reported, at motif-family rather than individual-factor level.

### Literature audit

Four PubMed queries returned 469 records, of which 209 were eligible and 154 had
an open-access full text, 152 of them scorable (Supplementary S23). Five scoring
criteria were fixed before scoring. Because automated matching is imprecise on
this task, a fixed-seed random 30% subsample (46 papers, 92 judgements) was
scored by **two coders independently, each blind to the other and to the
automated score**: the handover file carried neither the first coder's decisions
nor the automated label, and rows were shuffled under a fixed seed so that
automated positives were not clustered. Measured precision was 0.00 and 0.20 on
the two criteria, and the field-level estimates reported in the text are
precision-corrected and supersede both the generous and the strict automated
counts. Cohen's κ is reported beside the raw agreement, because at one positive
in 46 the expected agreement equals the observed and κ collapses to 0 while
agreement is 97.8%. Disagreements were adjudicated jointly, never by the first
coder alone, and κ is computed before adjudication. Recall was estimated by
probing negatives for near-miss wording rather than by reading all 152 in full,
so it bounds sensitivity from above; and a full-text audit measures reporting,
not practice.

On the same corpus we counted, by phrase matching, whether each paper names a
Wald ratio, IVW, weighted median, MR-Egger or mode-based estimator and whether
it states anywhere that its cis instrument is a single variant. Those counts are
reported as corpus description and not as an estimate of anything: the
single-variant phrase is matched anywhere in the text, so it is not tied to the
cis analysis, to the primary candidate list or to the estimator used, and the
count collapses from 71 to 10 when the phrase is required to appear in the same
sentence as a cis-eQTL mention. The interpretable quantity is which estimator
produced each paper's primary cis nomination list, and it was coded by hand as a
sixth criterion, C6, on the fixed-seed 46-paper subsample already used for the
other criteria. Two coders scored it independently, each blind to the other and
to the automated scan, under rules fixed before any paper was read, recording the
evidence passage and its location for every judgement; a judgement without a
recorded quote was to be treated as not ascertainable. Four codes were fixed in
advance — single-variant Wald, multi-instrument, mixed, and not ascertainable —
with a decision order placing an explicit per-gene instrument count above any
estimator named elsewhere in the paper, because naming IVW does not establish
that the cis list came from it. Cohen's κ is reported before adjudication and
beside the raw agreement, for the prevalence reason already given for C1;
disagreements were then resolved jointly, never by one coder alone. Papers coded
not ascertainable are reported as such and never redistributed (Supplementary
S51).

### Pre-specification and stopping rules

The following were fixed before the analyses they govern and are reported because
they determined which results are presented. Every test carries a positive
control, and a test whose positive control fails is discarded rather than
interpreted — the within-lymphoid-subset spatial test (HLA-C P = 0.23) and the
pre-specified immune-module analysis were discarded under this rule. Kill
criteria, test direction and window sizes were written into each script before it
ran. A negative result was not followed by a search for a positive one in the
same data. Simulations require calibration against a quantity not used to build
them. Seventeen pre-registration documents (Supplementary S9, S18, S20, S21, S22,
S24, S25, S28, S29, S30, S32, S33, S34, S35, S52, S53, S54) each carry a
reading
table written
before the run, and only their results register was completed afterwards;
deviations are logged in the document they belong to. Five further documents
sit beside that series and are not pre-registrations, which each of them says
in its own header: a frozen decision record on the locus partition (S37); an
estimand frozen except for one tolerance that remained an analysis choice
(S38); a post-hoc exploratory control (S39); an external-validation protocol
adopted but not registered, on three counts of its own — no qualifying cohort
exists, no registry DOI was obtained, and no environment was locked beforehand
(S40); and a transport test whose primary endpoint was fixed in advance but
which we do not call pre-registered, because it carries no registry DOI (S41).

### Software

R 4.4.1 with Seurat 5.5.1, Matrix 1.7.0, data.table 1.16.0, coloc 5.2.3,
TwoSampleMR 0.7.5, susieR 0.14.2, arrow 25.0.0, hdf5r 1.3.12, TFBSTools 1.42.0,
JASPAR2020 0.99.10, motifmatchr 1.26.0, chromVAR 1.26.0,
BSgenome.Hsapiens.UCSC.hg38 1.4.5, survival 3.8.9 and org.Hs.eg.db 3.19.1.
Python 3.12.4 with numpy 2.0.0, pandas 2.2.2, scipy 1.18.0, pyarrow 25.0.0 and
matplotlib 3.11.1. External binaries: SMR v1.3.1 and PLINK v2.0.0-a.7.2. Ordinary
least squares and Benjamini–Hochberg FDR were implemented directly in numpy
rather than through `statsmodels`, whose installed version is incompatible with
the installed scipy.

### Data and code availability

All datasets are public and identified by accession above. Analysis code,
intermediate result tables and the seventeen pre-registration documents will be
deposited at ⟨repository DOI⟩, comprising every numbered analysis and figure
script together with the tables needed to reproduce each figure and every number
reported in the text. CQTNA, a runnable implementation of the diagnostics in the
Discussion, is included in the deposit, as is a container definition pinning
R, Python and the two external binaries to the versions used here. The
analyses were run on Windows and the container is Linux, so it fixes the
software versions and not the machine; the deposit states which numbers were
checked across the two and which were not.

## Declarations

### Ethics approval and consent to participate

This study analysed only publicly available data: published or consortium-released
genome-wide association summary statistics, published cis-eQTL summary statistics,
and de-identified single-cell, spatial and chromatin datasets deposited in public
repositories. No new human data were collected and no individual-level
participant data were accessed. Ethical approval and informed consent were
obtained by the original studies and are described in their primary publications
and repository records, which are cited by accession in the Methods.

### Consent for publication

Not applicable.

### Availability of data and materials

All datasets are public and identified by accession in the Methods. Analysis
code, intermediate result tables, the pre-registration documents and the
container definition are in the deposit described under Data and code
availability.

### Competing interests

⟨competing interests statement⟩

### Funding

⟨funding statement, including grant numbers and the role of each funder in
design, analysis, interpretation and writing⟩

### Authors' contributions

⟨author contributions, by initials⟩

### Acknowledgements

⟨acknowledgements⟩

## Supplementary information

S9–S54, including the seventeen pre-registration documents with their reading tables
and results registers; **S41**, the transport grid scoring six diseases against
their own reference lists; **S42**, the estimator identity, its numerical
verification and the decomposition of the attribution by outcome significance;
**S43**, the power-matched cross-disease comparison and the chance-corrected
attribution statistic; **S44**, the seven-way method benchmark, run as a boundary check on the
identity rather than as seven independent tests; **S45**, the
fine-mapped decomposition and two failure modes of proxy-LD fine-mapping;
**S46**, the diagnostic's own type I error, power and void rate; **S47**, the
worked example in full; **S48**, the crossed grid cell by cell with its
registration history; **S49**, the eight checks in full; **S50**, which
inference unit the significance belongs to; **S51**, how far the identity reaches into the audited
literature, the phrase-matching attempt we withdrew and the hand coding that
replaced it; **S52**, the pre-registered lookup of the nominated gene in a
genome-scale CD4⁺ T perturbation map; **S53**, whether cis-eQTL-instrumentable and
perturbation-readable genes are the same genes; and **S54**, the rules under
which the deposited container becomes the canonical environment, fixed before
the image was built.
Earlier items include the literature-audit corpus with every PMID, both scoring
passes, the two-coder subsample and the coding rules (S23); the multiple-testing
unit sensitivity analysis, in which the significant list is recomputed over four
units by minimum-p and by Simes (S26); the self-administered attribution check
(S27); the complete record of target-substantiation attempts with the selection
denominators and every stopping-rule instance, and the technical account of the
two processing errors (S12); the replication-cohort search (S14); and the
colocalisation window-and-prior sensitivity analysis (S15).

## References

1. Zheng J, Yang Q, Liu H, et al. Integrating single-cell transcriptome-wide Mendelian randomization and differentially expressed gene analyses to prioritize dynamic immune-related drug targets for cancers. *Adv Sci* 2025;12:e07451. doi:10.1002/advs.202507451
2. Wu X, Ying H, Yang Q, et al. Transcriptome-wide Mendelian randomization during CD4⁺ T cell activation reveals immune-related drug targets for cardiometabolic diseases. *Nat Commun* 2024;15:9302. doi:10.1038/s41467-024-53621-7
3. Cui K, Zou Q, Qu X, et al. Transcriptome-wide Mendelian randomization and single-cell analysis during CD4⁺ T cell activation deciphers immunotherapeutic targets for colorectal cancer. *npj Precis Oncol* 2025;10:32. doi:10.1038/s41698-025-01236-6
4. Nathan A, Asgari S, Ishigaki K, et al. Single-cell eQTL models reveal dynamic T cell state dependence of disease loci. *Nature* 2022;606:120–128. doi:10.1038/s41586-022-04713-1
5. Schmiedel BJ, Singh D, Madrigal A, et al. Impact of genetic polymorphisms on human immune cell gene expression. *Cell* 2018;175:1701–1715.e16. doi:10.1016/j.cell.2018.10.022 (DICE)
6. Randolph HE, Fiege JK, Thielen BK, et al. Genetic ancestry effects on the response to viral infection are pervasive but cell type specific. *Science* 2021;374:1127–1133. doi:10.1126/science.abg0928
7. Tambets R, et al. Extensive co-regulation of neighboring genes complicates the use of eQTLs in target gene prioritization. *Hum Genet Genomics Adv* 2024;5:100348. doi:10.1016/j.xhgg.2024.100348
8. Howe LJ, et al. Evaluating transportability of in vitro cellular models to in vivo human phenotypes using gene perturbation data. *Nat Commun* 2025. doi:10.1038/s41467-025-67199-1
9. Lin Z, Pan W. A robust cis-Mendelian randomization method with application to drug target discovery. *Nat Commun* 2024. doi:10.1038/s41467-024-50385-y
10. Karhunen V, et al. Integrating genetic data with biological insight: a practical guide to cis-Mendelian randomization. *Am J Hum Genet* 2026. doi:10.1016/j.ajhg.2026.03.011
11. Reales G, et al. Design and interpretation of eQTL–GWAS colocalisation studies: lessons from a large-scale evaluation. *PLoS Genet* 2026. doi:10.1371/journal.pgen.1012141
12. Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human gene expression and splicing quantitative trait loci. *Nat Genet* 2021;53:1290–1299. doi:10.1038/s41588-021-00924-w (eQTL Catalogue)
13. Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nat Biotechnol* 2020;38:675–678. doi:10.1038/s41587-020-0546-8
14. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell* 2018;173:400–416.e11. doi:10.1016/j.cell.2018.02.052
15. Soskic B, Cano-Gamez K, Smyth DJ, et al. Immune disease risk variants regulate gene expression dynamics during CD4⁺ T cell activation. *Nat Genet* 2022;54:817–826. doi:10.1038/s41588-022-01066-3
16. Võsa U, Claringbould A, Westra H-J, et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. *Nat Genet* 2021;53:1300–1310. doi:10.1038/s41588-021-00913-z (eQTLGen, 2019-12-11 cis-eQTL release, n = 31,684)
17. Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature* 2023;613:508–518. doi:10.1038/s41586-022-05473-8 (Author Correction: *Nature* 2023;615:E19. doi:10.1038/s41586-023-05837-8)
18. Rashkin SR, Graff RE, Kachuri L, et al. Pan-cancer study detects genetic risk variants and shared genetic basis in two large cohorts. *Nat Commun* 2020;11:4423. doi:10.1038/s41467-020-18246-6 (GCST90011809)
19. Landi MT, Bishop DT, MacGregor S, et al. Genome-wide association meta-analyses combining multiple risk phenotypes provide insights into the genetic architecture of cutaneous melanoma susceptibility. *Nat Genet* 2020;52:494–504. doi:10.1038/s41588-020-0611-8 (GCST010302–GCST010304)
20. FinnGen. Release R12 (2024), endpoints `C3_MELANOMA_SKIN_EXALLC` and `C3_HEPATOCELLU_CARC_EXALLC`; Release R13, endpoints `C3_MELANOMA_SKIN_WIDE` and `M13_RHEUMA`. https://www.finngen.fi/en/access_results
21. Ghouse J, Gellert-Kristensen H, O'Rourke CJ, et al. Genome-wide meta-analysis identifies nine loci associated with higher risk of hepatocellular carcinoma. *JHEP Rep* 2025;7:101485. doi:10.1016/j.jhepr.2025.101485 (GCST90809296)
22. Byrska-Bishop M, Evani US, Zhao X, et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell* 2022;185:3426–3440.e19. doi:10.1016/j.cell.2022.08.004
23. Chang CC, Chow CC, Tellier LC, et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 2015;4:s13742-015-0047-8. doi:10.1186/s13742-015-0047-8
24. Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal inference across the human phenome. *eLife* 2018;7:e34408. doi:10.7554/eLife.34408
25. Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet* 2014;10:e1004383. doi:10.1371/journal.pgen.1004383
26. Wang G, Sarkar A, Carbonetto P, Stephens M. A simple new approach to variable selection in regression, with application to genetic fine mapping. *J R Stat Soc Series B* 2020;82:1273–1300. doi:10.1111/rssb.12388
27. Zhu Z, Zhang F, Hu H, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. *Nat Genet* 2016;48:481–487. doi:10.1038/ng.3538
28. Sollis E, Mosaku A, Abid A, et al. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. *Nucleic Acids Res* 2023;51:D977–D985. doi:10.1093/nar/gkac1010
29. Okada Y, Wu D, Trynka G, et al. Genetics of rheumatoid arthritis contributes to biology and drug discovery. *Nature* 2014;506(7488):376–381. doi:10.1038/nature12873 (GWAS Catalog GCST002318)
30. Sade-Feldman M, Yizhak K, Bjorgaard SL, et al. Defining T cell states associated with response to checkpoint immunotherapy in melanoma. *Cell* 2018;175:998–1013.e20. doi:10.1016/j.cell.2018.10.038 (GSE120575)
31. Jerby-Arnon L, Shah P, Cuoco MS, et al. A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. *Cell* 2018;175:984–997.e24. doi:10.1016/j.cell.2018.09.006 (GSE115978)
32. Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 2016;352:189–196. doi:10.1126/science.aad0501 (GSE72056)
33. Hugo W, Zaretsky JM, Sun L, et al. Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. *Cell* 2016;165:35–44. doi:10.1016/j.cell.2016.02.065 (GSE78220)
34. Riaz N, Havel JJ, Makarov V, et al. Tumor and microenvironment evolution during immunotherapy with nivolumab. *Cell* 2017;171:934–949.e16. doi:10.1016/j.cell.2017.09.028 (GSE91061)
35. Guo X, et al. Contrasting cytotoxic and regulatory T cell responses underlying distinct clinical outcomes to anti-PD-1 plus lenvatinib therapy in hepatocellular carcinoma. *Cancer Cell* 2025;43:248–268.e9. doi:10.1016/j.ccell.2025.01.001 (GSE235863)
36. Thrane K, Eriksson H, Maaskola J, et al. Spatially resolved transcriptomics enables dissection of genetic heterogeneity in stage III cutaneous malignant melanoma. *Cancer Res* 2018;78:5970–5979. doi:10.1158/0008-5472.CAN-18-0747
37. Virós A, et al. Spatial transcriptomics of primary cutaneous melanoma. NCBI GEO GSE316760 (submitted 2026-01-16); no associated publication indexed at the time of writing.
38. Pham F, Dufeu M, Benboubker V, et al. Spatial tumour-immune ecosystems shape the efficacy of anti-PD1 immunotherapy in primary cutaneous melanoma. NCBI GEO GSE300445 (submitted 2025-06-23); no associated publication indexed at the time of writing.
39. Katko A, Potter SJ, et al. Gene regulatory network determinants of rapid recall in human memory CD4⁺ T cells. *Cell Rep* 2026;45:117103. doi:10.1016/j.celrep.2026.117103 (GSE282266)
40. Boukhaled GM, Gadalla R, et al. Pre-encoded responsiveness to type I interferon in the peripheral immune system defines outcome of PD1 blockade therapy. *Nat Immunol* 2022;23:1273–1283. doi:10.1038/s41590-022-01262-7 (GSE199994)
41. Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat Biotechnol* 2024;42:293–304. doi:10.1038/s41587-023-01767-y
42. Fornes O, Castro-Mondragon JA, Khan A, et al. JASPAR 2020: update of the open-access database of transcription factor binding profiles. *Nucleic Acids Res* 2020;48:D87–D92. doi:10.1093/nar/gkz1001
43. Schep AN, Wu B, Buenrostro JD, Greenleaf WJ. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nat Methods* 2017;14:975–978. doi:10.1038/nmeth.4401
