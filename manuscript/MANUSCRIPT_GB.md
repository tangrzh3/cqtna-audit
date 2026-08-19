# Outcome GWAS architecture shapes target nomination from context-specific eQTLs: an audit across three diseases and two exposure resources

## Abstract

**Background.** Context-specific expression quantitative trait loci (eQTLs) with
Mendelian randomization (MR) are widely used to nominate immune targets in
cancer. The instrument is usually a single variant, so the Wald statistic reduces
to z = β_out/se_out and the outcome GWAS supplies the significance claim. We audited
this, instrumenting CD4⁺ T cell cis-eQTLs from eight activation profiles against
a 12,530-case melanoma meta-analysis, then varying outcome and exposure resource
independently.

**Results.** Significant signal concentrates on loci already known for the
outcome: 4.96-fold by independent locus. It recurs across five nested power levels, where no
novel-locus gene reaches significance at any case number; in a second disease;
and when the exposure is replaced by a whole-blood eQTL dataset 300-fold larger
(6.31-fold, P = 1.3×10⁻¹⁵). Of six cells of three diseases by two resources, one
is void on its pre-registered mismatched-locus control; the remaining five all
enrich, four significantly. Outside cancer the picture is weaker and we report it
as such: rheumatoid arthritis enriches 3.91-fold (P = 2.3×10⁻⁷) on CD4⁺ T cells,
but its mismatched-list control is clean only at the registered window and fails
at every tighter one, so the non-cancer cell is descriptive rather than a
demonstration of generality. Naming the
gene at a significant locus then fails separately. Of 291 records that colocalisation assigned to distinct causal variants, SMR/HEIDI
failed to reject homogeneity for 253, so tiers accepting MR plus SMR alone would
have reported MC1R linkage spillover as CD4-mediated targets. At ten melanoma loci
with a generally accepted causal gene, our pipeline names it at six and fails at
four; at MC1R it spans sixteen genes without naming MC1R. Replacing the
outcome with a higher-powered, compositionally different meta-analysis increased
MR discoveries, lowered colocalisation support and replaced the candidate list
entirely: two lists from identical exposure data share no genes. Recovery under
down-sampling falls by effect size and hence by locus class — at half power
known-locus genes recover 85.8% against 22.8%, 68–80% of that gap following from
|z| alone. Of 28 glycolytic genes, 3 are instrumentable, 2 analysable and 1
nominally associated. In a two-coder subsample of 152 comparable studies, no
paper reported the locus-attribution check.

**Conclusions.** Within the power range we could observe, the reproducible part
of a candidate list produced by this framework was the part that did not
constitute a discovery. **Which loci a nomination lands on tracks the outcome
GWAS; how many instruments exist tracks the exposure resource** — the second
naming a resource, not a sample size. What survives is a statement about the
state in which a pathway's regulation can be measured precisely enough to
instrument at all. Eight inexpensive checks follow.

**Keywords** Mendelian randomization · context-specific eQTL · target nomination ·
colocalisation · statistical power · reproducibility · compartment attribution ·
melanoma

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

### Significant signal sits on loci already known for the outcome

In the FinnGen round, all ten FDR-significant records fell within 1 Mb of a known
pigmentation or naevus locus — VPS9D1-AS1 50 kb from MC1R, CDK10 48 kb from the
MC1R R151C variant, PARP1 5–14 kb from PARP1 — with no immune signal at
FDR < 0.05 at all. Under the meta outcome, counting independent loci rather than
gene records, 4 of 8 significant loci were known-locus loci: a 4.96-fold
enrichment over the 10.1% background (one-sided P = 0.0048; Fig. 1). Loci here
are the non-recursive fixed-anchor partition fixed in Supplementary S37 — a
1 Mb window claimed from each anchor, so a locus spans at most 1 Mb whatever the
variant density; the single-linkage rule used in the source literature chains a
dense resource into blocks of tens of megabases and is reported only where we
reproduce published numbers (Methods). Applied
unchanged to five other cancers, melanoma showed the largest enrichment and the
two cancers with no expected relationship to pigmentation sat near background,
though the comparison does not deliver clean specificity and no outcome survives
correction across six tests.

We then tested the same attribution three further ways, changing one thing at a
time (Fig. 2).

**Five nested power levels of one resource.** FinnGen's sequential releases carry
an identical melanoma endpoint at 2,705, 2,993, 3,194, 3,932 and 5,753 cases,
with exposure data, instruments, harmonisation, estimator, testing family and
locus definition held fixed. Predictions were generated by down-sampling the
largest release and registered before any earlier release was examined; all twelve
— discovery counts, overall recovery, and the known-minus-novel differential at
each release — fell inside the registered intervals. The known-locus part of the
list is recovered early and then does not move: the same five genes constitute the
significant list at three consecutive releases without change, four already
present at 47% of reference power. And **no novel-locus gene reaches FDR < 0.05 in
any release**. Because releases are nested, agreement between them exceeds that
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

**A second disease.** Repeating the nomination against hepatocellular carcinoma
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

**A second exposure resource.** Holding the melanoma outcome fixed byte-for-byte
and replacing the exposure entirely with eQTLGen whole-blood cis-eQTLs
(n = 31,684, roughly 300-fold larger, different cells, different platform) raises
instruments from 3,556 to 12,835 records and significant loci from 7 to 30. The
attribution does not move: **23 of 34 significant loci (67.6%) carry a known
melanoma, naevus or pigmentation lead SNP, 6.31-fold over this resource's own
10.7% background (P = 1.3×10⁻¹⁵)**, against 4.96-fold in CD4⁺ T cells. A
matched-background version agrees, drawing each significant locus a background
locus of the same eQTL-p decile and allele-frequency quintile (6.49-fold,
empirical P = 10⁻⁴; 5.66-fold on the decile alone). The two
resources' novel nominations intersect in exactly one gene, ZFYVE19, so the
stronger claim that novel nominations never reproduce is not supported.

**Three further exposure resources, and why they settle nothing.** eQTLGen is a
large whole-blood dataset, not a second stimulation time course, so we
pre-registered an attempt to add context-specific CD4⁺ resources with the outcome
held byte-identical (Supplementary S35): anti-CD3/CD28 at 4 h, influenza
infection at 6 h, and an activation-state resource in unstimulated memory cells.
None supports the attribution claim. Two fall below the registered floors — 2
significant loci and none at all — and are reported as underpowered rather than
as nulls, since a locus-level test on such numerators has no power to be either.
The third clears both floors at 3.41-fold, but **its pre-registered
mismatched-list control fails: the wrong disease's known loci give 5.79-fold
(P = 0.042) against 3.41-fold (P = 0.11) for melanoma's own**, each list hitting
two of six significant loci and not the same two — and a mismatched list that
also enriches voids the cell. Why these resources cannot decide is itself the
point: at 89 to 147 donors they yield 112 to 860 usable instruments and 0 to 6
significant loci, against 3,556 and 7 for the CD4⁺ time course and 12,835 and 30
for whole blood. **The finding rests on two exposure resources and we could not
extend it to a third.**

**Outside cancer, where the exposure cell type is part of the causal
compartment.** In neither melanoma nor HCC is the CD4⁺ T cell an accepted causal
cell type, so the attribution finding could be what happens when the exposure
tissue is irrelevant to the disease. Rheumatoid arthritis removes that objection:
**CD4⁺ T cells constitute a well-established pathogenic immune compartment in
RA** — one among several that also include B cells, myeloid cells and synovial
fibroblasts — and FinnGen's RA endpoint has 16,775 cases against the melanoma
meta's 12,530, an effective sample size about 29% larger, so **a null could not
readily be attributed to a smaller outcome effective sample size**. Scored
against Okada et al.'s 2014 lead SNPs [47], published before FinnGen existed and from
cohorts that exclude it, registered in advance (Supplementary S33): **16 of 41
FDR-significant loci carry a known RA lead SNP against a 10.0% background,
3.91-fold, P = 2.3×10⁻⁷**. S33 §4 makes the CD4⁺ cell the primary and §5 defines
the test over the whole genome, so that is the cell reported here; the MHC
positive control passes, and a post-hoc version excluding the MHC region leaves
3.48-fold (P = 2.1×10⁻⁴). We label that version post-hoc rather than
pre-registered because it appears only in S33's results table, written after the
RA result was seen — the same after-the-fact change of convention this paper
criticises elsewhere, and one we made and corrected (Supplementary S37).

**The whole-blood RA cell is void.** Scored on the same records, it reaches
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
loci, melanoma and hepatocellular carcinoma do not. Second, a per-locus breakdown
rules out the convention: of the eight loci the wrong list claims at 500 kb, none
is claimed through a non-significant record; every one has a significant
nomination genuinely within 1 Mb of a melanoma lead. Consistent with that, RA's
hit distances are spread across the intermediate range rather than bimodal — 6 of
41 within 50 kb against 15 of 34 for melanoma on whole blood — so RA's hits sit
precisely in the band that a narrowing window removes.

**A mismatched list is therefore not an orthogonal negative control for this
outcome pair, and we do not read it as one.** By the frozen reading table
(Supplementary S37 §4) the wording is downgraded accordingly: the enrichment
itself is stable across every setting (3.70- to 3.96-fold), but **we report the
non-cancer cell as descriptive attribution rather than as a demonstration that
the phenomenon generalises beyond cancer**, and give the full sweep in
Supplementary S36 rather than the single setting that passes.

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

**Both axes crossed.** The main grid is the three diseases at their higher-powered
outcome crossed with the two exposure resources — six cells — with the
lower-powered HCC outcome held back as a power sensitivity, since the two HCC
levels are drawn from one resource and are not independent (Supplementary S20 §9).
Scoring every cell against its own disease's known loci, **one of the six is void
and the other five all enrich**: 4.96 and 6.31 for melanoma, 7.88 and 10.87 for
HCC-high, 3.91 for RA on CD4⁺ T cells, with RA on whole blood void on its
mismatched control as above. **Four of the five reach P < 0.05**; the one that
does not is HCC-high on CD4⁺ T cells (P = 0.123), on a numerator of one of two
loci, so this is consistency of direction rather than five independent
significant tests. Note that the two non-significant cells are no longer the same
disease: HCC-high on whole blood now reaches P = 0.0015, because the fixed-anchor
partition splits that resource's 562 chained blocks into 1,666 bounded loci and
the background known rate falls from 7.8% to 5.5%. The lower-powered HCC outcome
enriches on both resources as well (15.5-fold, P = 0.0041, and 13.5-fold,
P = 0.00063), which is the sensitivity and not a seventh and eighth independent
cell. A pre-registered mismatched-locus control makes the rest interpretable:
scoring the same combinations against the wrong disease's list collapses the
enrichment in every surviving cell — HCC's list on melanoma gives 1.59-fold
(P = 0.29) on whole blood and 0.00-fold on CD4⁺ T cells, melanoma's list on HCC
gives 0.00-fold (P = 1.0) — so in those cells the effect is specific to each
outcome's own genetics. A mismatched list is a narrower control than it
looks — it rules out enrichment on loci indiscriminately dense across diseases,
not a density that is itself disease-specific — so two permutations close the gap
from different sides. Matching on eQTL-p decile and allele-frequency quintile
controls instrument strength and frequency (6.66-fold, empirical P = 0.0016 for
melanoma on CD4⁺ T cells). Matching on **density** — the estimand frozen in
Supplementary S38 before it was run, drawing background loci matched on record
count, span and gene count, requiring 100% matched coverage and reporting the
whole tolerance scan — gives **4.95-fold (empirical P = 0.0042)** for melanoma on
CD4⁺ T cells, 6.16-fold (P = 10⁻⁴) on whole blood, 8.47-fold (P = 0.0041) for
HCC-high on whole blood, 12.95- and 12.46-fold (P = 0.0049 and 0.0010) for
HCC-low, and **4.50-fold (P = 10⁻⁴) for RA on CD4⁺ T cells**. Every one of those
holds across all tolerances that reach full matched coverage; the tightest
tolerance fails coverage on four cells and returns no P value rather than a
partially matched one (Supplementary S38 §1.7).

Two readings follow. First, on bounded loci density matching no longer
attenuates: melanoma on CD4⁺ T cells is 4.96-fold unmatched and 4.95-fold
density-matched. This is not evidence that the attenuation reported earlier was
spurious — the earlier 3.44-fold came from quantile-stratified matching, a
different test that we do not treat as the same quantity (Supplementary S38 §4) —
but under the frozen estimand and the corrected partition, locus density does not
account for the enrichment. Second, **HCC-high on CD4⁺ T cells does not clear the
density-matched null** (P = 0.105, and P ≥ 0.05 at every usable tolerance), so by
S38's reading table its attribution claim is descriptive only; that is the same
cell and the same power limit that leaves it at P = 0.123 on the Fisher test, so
the two controls agree rather than conflict. The density-matched control now
covers all eight cells including RA, so the RA cell no longer rests on the
mismatched-list control alone — which matters, because that control is the less
clean one precisely in RA.

One observation from the grid places this work against the complementary
literature. Swapping in the eQTLGen resource multiplied significant loci in
melanoma (7 → 30) but not in HCC (2 → 5 and 2 → 3), suggesting that **the exposure
resource sets how many instruments exist while outcome power sets how many of
them can reach significance.** The swap is not a clean manipulation of exposure
sample size and we do not present it as one: donors, cell composition and
platform all change together, and this design cannot separate them. What is
identified is the asymmetry — the same swap moved melanoma fourfold and HCC
hardly at all, and the outcome is the only thing differing between those arms —
so the outcome-side half of that statement is supported while the exposure-side
half names a resource rather than a sample size. Rosen et al. have shown that raising *eQTL* sample
size uncovers additional independent regulatory signals and closes part of the
gap between eQTL and GWAS colocalisation [43]. Our design asks the mirror
question — with the exposure held fixed, does a *higher-powered outcome* make
target nomination more reliable? — and the answer here is that it does not: the
meta-analysis produced more MR discoveries, lower colocalisation support, and a
candidate list with no genes in common with the previous one. That outcome
differs from its predecessor in study composition as well as power, so this is
not a power experiment either; only the nested-release series isolates power, and
there the significant list is stable rather than replaced. The two axes are not
interchangeable, and a study underpowered on one cannot be rescued by the other.
Reporting which axis a claim rests on should be routine.

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
residual is **+5.2 percentage points [−1.6, +12.0] by independent locus** and
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

This is the co-regulation problem of Tambets et al. [42] observed at the top of
the effect-size distribution, and it is why we report compartment attribution
separately from gene attribution — the two fail independently.

### Which unit the significance belongs to

A units problem has to be settled before any of the counts above is read. Those
3,556 records carry only **2,126 unique variants and 1,195
unique genes**, because one variant can be the lead eQTL for a gene in several
profiles and for more than one gene, so the testing family is records while the
list is reported by gene and the attribution by independent locus. Records are
therefore not independent, and a Benjamini–Hochberg procedure over them does not
have a clean guarantee under that dependence.

We resolve this by separating the two things the record level is being asked to
do. **The record-level analysis reproduces the nomination pipeline being
audited** — it is what the studies under examination run, and what every
registered prediction, the release trajectory, both generalisations and the grid
were computed on, so it stays exactly as published. **The independent locus is
the unit the audit's own conclusions are stated in.** That is not a convenience:
the FDR < 0.05 list recomputed under all seven alternatives — variant, gene and
locus, each by minimum-p and by Simes, plus a two-stage hierarchical procedure —
leaves the locus count stable at 6 to 9, against 8 in the
main analysis, and 2 under every unit in the FinnGen round (Supplementary S26,
recomputed on the frozen partition). That sweep varies the **testing unit** while
holding the partition fixed, so it shows the conclusion survives the choice of
testing unit; the partition itself is varied separately, over four widths in
Supplementary S36. The audit conclusion therefore does not rest on the disputed
unit at all. The gene-naming benchmark below is likewise unaffected: changing the
partition redistributes the same 56 named genes over 23 loci instead of 20,
without adding or removing a single gene, so a benchmark scored on which gene is
named cannot move. The record level is also the most conservative of the eight — no
gene is lost under any other unit and up to 43 are added — but that is an
empirical observation about list length, **not evidence that FDR is correctly
controlled under this dependence**, and we do not use it as such.

One row reads as a finding rather than a check: testing at locus level leaves the
locus count almost unchanged but inflates the gene list from 10 to 53, of which 22
are in the MHC and 7 in the chr17q21.31 inversion, because a significant locus
does not name a gene.

### Part 3 — What survives the audit

What is left once nomination is set aside is a narrower statement about where a
pathway can be interrogated at all, and it is reported as a ladder rather than a
verdict.

### What can be interrogated, stated at the level it applies to

Across 22 glycolytic enzymes measurable in CD4⁺ T cells, 20 are higher in
non-responders post-treatment and 18 of 22 pre-treatment, with PGAM1 strongest at
baseline and TPI1 fourth. Genetically, of 28 glycolytic genes **3 are
instrumentable in this exposure resource** (a genome-wide significant cis-eQTL in
some activation profile), **2 remain analysable against this outcome** after
harmonisation, and **1 yields a nominal association** (Fig. 5).

These three numbers are not the same kind of statement and we keep them apart.
The first is a property of the pathway and the exposure data and would be
unchanged whatever outcome were analysed. The second adds the requirement that the
variant harmonise with a particular outcome GWAS — SLC2A1 fails on this outcome
and might not on another, so it is not a statement that SLC2A1 lacks an
instrument. The third adds the outcome's own significance. Compressing them into
"one usable instrument" would attribute two outcome-side limitations to the
biology of the pathway. The exposure-side constraint alone — 3 of 28 — is the
single most restrictive step in the framework, and functionally the strongest
baseline enzyme, PGAM1, is among the 25 that carry no instrument at all.

### Time-specific significance is not a dynamic genetic effect

For the one instrumented gene, TPI1, a usable instrument exists at 16 h and at no
other timepoint. The effect estimate is largest at rest and simply far noisier
there, so the window is a statement about precision, not about amplitude
(Fig. 6a,b). The published genotype × pseudotime interaction test [8] settles what
that window is not: **none of the three TPI1 lead variants shows an
interaction** — memory cells, linear P = 0.974 and quadratic P = 0.992; naive
cells, 0.707 and 0.891 for the variant 1,063 bp from our instrument, 0.378 and
0.605 for the other — while the gene's *expression* is among the most strongly
pseudotime-dependent in the dataset (Moran's I = 0.664). Three things are
separable, and the interaction test tells them apart: dynamic expression (yes,
strongly), time-specific eQTL significance (yes, at 16 h only) and dynamic
genetic effect (no evidence, in any of three tests). The design consequence holds
regardless: a study using a resting-state resource would find nothing here and
could not distinguish that from a true null.

### What the instrumented gene marks

Because MR names one gene of a pathway, we asked what that gene marks. In
purified CD4⁺ T cell multiome data, residualising expression on activation
intensity leaves an axis enriched 21.3-fold for glycolysis and 10.3-fold for
ribosomal proteins, while **the activation module itself sits at a
composition-matched null (SMD = −0.003), and of the eleven biological modules
tested only OXPHOS and the proliferation positive control exceed that null** —
a definable anabolic state rather than a restatement of activation strength,
with a concordant chromatin signature whose AP-1 motif enrichment persists in
activation-invariant peaks (Fig. 7).

As first built the analysis was circular with respect to TPI1, which entered both
the defining score and the enrichment family. Rebuilt without it, everything else
held fixed, the enrichment is essentially unmoved (**21.0-fold on 14 genes
against 21.3 on 15**) and TPI1 becomes a held-out test: on an axis derived
without it, **TPI1 ranks 19th of 7,653 detected genes**, positive in all four
sets. Pre-registered replication in a second dataset (Supplementary S32) —
DOGMA-seq from a different study and platform, where CD4 is called by surface
antibody rather than by transcript, the calling that finding ⑦ shows to be
unreliable — places it **18th of 8,224 and 28th of 7,165** in two arms, with
glycolysis enriched 27.4- and 23.9-fold. Those two arms are two lysis buffers on
one batch of material, so this replicates across study, platform and lineage
calling and **not across donors**.

The same ranking restates the instrument-visibility point from the other side:
**PGAM1 ranks 8th, above TPI1**, and 3rd in both replication arms, with LDHA,
GAPDH, PKM and ENO1 at or above it — the enzymes MR cannot instrument sit as high
on the axis as the one it can. All of this characterises what the nominated gene
co-varies with, not that its variant causes the state; a direct test of the
latter, whether the instrument disrupts an AP-1 motif, returned an empirical
P = 1.0.

### Where the effect is, when MR cannot say

In melanoma single-cell data with all cell types annotated, TPI1 is far higher in
malignant cells than in CD4⁺ T cells: +2.75 log₂ units (≈6.7-fold) in 16 of 16
patients (P = 3×10⁻⁵), detected in 95.3% of malignant cells against 60.6% of CD4⁺
T cells. The difference survives cell-level matching on sequencing depth (+2.03,
10 of 11 patients), replicates in an independent cohort (+1.28, 11 of 11,
P = 1×10⁻³) with pre-specified positive controls passing in both, and is
concordant in direction in a third. The same holds for the locked glycolytic
signature with and without TPI1 (Fig. 8).

The arithmetic consequence is what matters. With those ratios and the cell-type
proportions typical of melanoma tissue, CD4⁺ T cells contribute on the order of
one to two per cent of the tissue-level TPI1 signal, so **bulk-tissue TPI1
abundance is dominated by the malignant compartment and cannot validate a
CD4-specific mechanism**, whatever its P value. The looser claim — that such a
measurement simply *is* a measurement of tumour glycolysis — would overstate it:
bulk TPI1 can still covary with outcome through immune infiltration, tumour
purity or a metabolic state shared across compartments. What the ratio
establishes is that such an association cannot be *attributed* to CD4⁺ T cells,
not that it is spurious — and the check costs minutes.

### The functional layer has its own failure modes

Splitting cells by a score also splits them by lineage purity. When cells were
divided by glycolytic score, apparent differences in helper and regulatory
phenotype collapsed after matching on lineage composition at the cell level, and
14 of 28 module comparisons were concordant — exactly chance. A cluster-level
control passed while the cell-level control failed, so whether a purity control
works can depend on the level at which it is applied. Two processing errors of
our own, detected by these controls, are described in Supplementary S12; they are
instances of the same failure mode.

The perturbation layer is bounded in the same way. TPI1 is a hit in 628 of the
**1,471 human CRISPR screens** in BioGRID ORCS [40] that measured it (**42.7%**), a
core-essential profile alongside GAPDH and PGAM1 (46.6% and 47.2%) and an order
of magnitude above lineage-defining genes in the same screens (IRF4 3.9%, FOXP3
1.5%, MC1R 1.1%). Since the knockout is lethal in almost any cell type, **a
knockout cannot isolate a CD4-specific role** — a boundary on what perturbation
could add here, not evidence for the nomination. PGAM1, the strongest enzyme in
the functional data and the one MR cannot see, has the same profile, so the
property belongs to the pathway rather than to the named gene.

### In patients: present in each cohort, confirmed in none twice

Inference throughout is a composite signature score with a permutation of the
response label, fixed as primary in advance and applied identically in all three
cohorts. Both arms are significant in the discovery cohort (post-treatment
Δ = +0.87, P = 0.0005; pre-treatment Δ = +0.74, P = 0.009). Two follow-ups were
pre-registered — the one qualifying public replication cohort a fixed-criteria
search returned, and a cross-disease cohort that became available during revision.
**No arm was confirmed a second time**, and the two fail in mirror image:
pre-treatment matches in magnitude in the same-disease cohort but is underpowered
there and absent in the cross-disease one, while post-treatment reverses in the
same-disease cohort and reproduces at full magnitude in the cross-disease one
(Δ = +0.874, P = 0.043) (Fig. 9). The per-arm estimates, the three treatments of
repeated patients and the sign dependence on regulatory T cells are Supplementary
S19 and S21.

### These checks are rarely reported

From 469 screened records, 209 were eligible and 154 had an open-access full
text, of which 152 were scorable. In that sample, **not one paper in a random 30%
subsample compared its significant signal against previously reported loci for
its own outcome trait (0% [0–3.4]), and an estimated 7.1% [2.5–16.1] report how
the candidate list depends on the outcome GWAS used**. About 59% performed
colocalisation and 47% used SMR or HEIDI.

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
Recall was estimated by probing negatives for near-miss wording rather than by
reading all 152 in full, so it bounds sensitivity from above; and a full-text
audit measures reporting, not practice.

---

## Discussion

Everything in this framework that depends on the outcome GWAS proved unstable, and
the instability has a specific shape: signal lands on loci the outcome already
knows about, the novel part of the list turns over completely between outcomes,
and replacing the outcome with a larger and compositionally different
meta-analysis moved MR and colocalisation in opposite directions. The pattern
held when we changed the disease, when we changed the exposure resource by two
orders of magnitude, and in five of the six cells of the two crossed — three
diseases by two resources, four of them individually significant. In the sixth,
rheumatoid arthritis on whole blood, the mismatched-locus control itself
enriched, and we report the cell as void rather than as a result. The control is
the load-bearing part of that sentence, and it does not carry equal weight
everywhere: in the four cancer cells it is clean at every window we swept, while
in rheumatoid arthritis on CD4⁺ T cells it is clean only at the registered
window and fails at every tighter one. Where a negative control is scale-
dependent, so is the claim it licenses. What did not depend on the outcome
was a smaller and more durable statement: in which cell state a pathway's
regulation can be measured precisely enough to yield an instrument at all.

Of eleven target-level tests that returned a verdict, seven overturned the claim
under examination, three were inconclusive and one was significant on discovery
but not confirmed on transfer; five genes were successively designated lead
candidate under a ranking criterion that was never fixed in advance, and the
first four designations were overturned. That is this project's own history and
not an estimate of how often the framework fails, and the itemised record — every
attempt, stopping rule and withdrawal, the candidate-selection timeline, the
fourteen pre-registrations and the technical account of two processing errors of
ours — is Supplementary S12.

Eight checks follow directly, each cheap and each capable of changing what a study
of this kind reports.

**(i) Annotate the signal against the outcome's own known loci, counting
independent loci rather than gene records.** Every result in this paper's first
half rests on that distinction, and no paper in a two-coder subsample reported it.

**(ii) Require colocalisation, preferably with explicit multiple-signal modelling
and an LD reference matched to the outcome cohort, rather than treating SMR/HEIDI
as sufficient — and report how many variants entered each HEIDI test.** HEIDI
passed 253 of 291 records that colocalisation assigned to distinct causal
variants, on 8–20 variants per test. We stop short of naming any single
colocalisation method the final arbiter: ours permits at most one causal variant
per region, and our attempt to go beyond that assumption over-split MC1R relative
to in-sample fine-mapping, so it could constrain one locus in one direction only.

**(iii) Do not assume that more outcome power will resolve the disagreement.** In
the one comparison available here the two moved in opposite directions, but that
outcome differed from its predecessor in study composition as well as power, so
the divergence is not attributable to power alone. The two quantities are
sensitive to different things by construction — one to a single variant's z, the
other to the shape of the regional signal — which is reason enough not to treat a
larger outcome GWAS as the remedy for their disagreement.

**(iv) Report the physical distance between GWAS and eQTL peaks alongside
posteriors, and test whether posterior movement survives explicit modelling of
multiple causal signals before calling it a change in resolution.** For PARP1 the
collapse is not an unmodelled second signal; for two other loci no credible set
exists at either power, so those regions cannot adjudicate it.

**(v) Compute the nominated gene's cell-type expression ratio in an annotated
atlas before citing tissue-level data as validation.** This costs minutes, and
here it shows that a tissue-level measurement of the nominated gene is dominated
by the malignant compartment and therefore cannot be attributed to CD4⁺ T cells —
which is a bound on what such a measurement can establish, not a claim that the
association is spurious.

**(vi) Report instrument availability at each of the three levels it passes
through — instrumentable in the exposure resource, analysable against this
outcome, nominally associated — and do not compress them.** Here, 3, 2 and 1 of 28.

**(vii) When splitting cells by a score, match on lineage composition as well as
depth, at the cell level.** A cluster-level control passed while the cell-level
control failed.

**(viii) Before comparing candidate lists across releases of the same GWAS
resource, verify code-by-code that the endpoint definition is the same one.**
Release notes are not sufficient: FinnGen R13 describes our endpoint's successor
as "including Hilmo", which reads as newly added hospital-register cases, whereas
the code-level definitions show the cases are identically defined and what
changed is the control-exclusion rule. A study that reads the summary line and
compares anyway will attribute pure phenotype drift to instability of its
candidate list — the very quantity such a comparison is meant to measure. Where
the definition does differ, the comparison is still worth making but is a
transfer test rather than a power point, and the direction of each difference
should be recorded before the result is seen.

The worked example should be read in that light, and is best stated as a ladder
rather than a verdict — the same layered reporting Howe et al. use when a cellular
model proves only partly transportable [44]:

| Claim about TPI1 | Status here |
|---|---|
| Causal target for melanoma | **Not supported** — FDR = 0.119, prior-dependent colocalisation, region not fine-mappable |
| Dynamic *genetic* effect across activation | **Not supported** — no genotype × pseudotime interaction at any of three lead variants |
| Dynamic *expression* across activation | **Supported** — Moran's I = 0.664 |
| Membership of a definable CD4⁺ metabolic state | **Supported** — ranks 19th of 7,653 on an axis built without it, and 18th and 28th in a second dataset with protein-based lineage calling; axis enrichment 21.0–27.4-fold |
| Predicts checkpoint-blockade response | **Unstable** — significant on discovery, no arm confirmed twice |
| Functional consequence isolable by perturbation | **No** — a hit in 628 of 1,471 human CRISPR screens (42.7%), a core-essential profile that cannot isolate a CD4-specific role |
| Recognised as a target for this disease elsewhere | **No** — in Open Targets [35], TPI1's strongest disease associations are triosephosphate isomerase deficiency and neurodegenerative disease; melanoma is not among its leading associations |

Stated flatly: TPI1 is a nominal, prior-sensitive candidate selected by a
multi-layer conjunction — MR P = 1.4×10⁻³ but FDR = 0.119, PP.H4 = 0.51 passing
in 8 of 16 window-and-prior combinations, in a region the outcome GWAS cannot
fine-map, with gene attribution at chr12p13 unresolved by genetics. **What MR
supplies here is a genetically tractable entry point into a pathway, and the
audit establishes that the gene marks an expression state — not that the state is
genetically driven, and not that the gene is a target.** Its value is as a
demonstration that instrument availability and biological importance are
separable: the strongest enzyme in the functional data carries no instrument, and
the gene that does carry one is principally expressed by tumour. The complete
record of which claims about it were raised, tested and withdrawn, together with
the selection denominators at gene, pathway and mechanism level and the twelve
pre-registration documents, is Supplementary S12 and S18–S38.

**Relation to existing guidance.** These eight are additions to, not a
replacement for, current cis-MR practice. Existing guidance already stresses that
cis analyses must be tailored to local biology, that an expression biomarker is
not an intervention, that co-regulation of neighbouring genes can act as
horizontal pleiotropy, and that MR is one strand of triangulation [45,46]. Tambets
et al. showed with approximate ground truth that neighbouring-gene co-regulation
routinely produces colocalisation evidence for several genes at once, so eQTL
data suit candidate *generation* and evidence *composition* better than gene
*attribution* [42]; our MC1R region is the same phenomenon at higher effect size,
and TPI1 at chr12p13 is the same problem unresolved. Reales et al., auditing over
a million colocalisation tests across immune diseases and cell types, likewise
found that gene assignment shifts with platform coverage, cell context and
resource size [41]. What we add is specific to the nomination step and to the
outcome side: locus attribution against the outcome's own known loci, sensitivity
to the outcome release, candidate-list stability, layered instrument visibility,
and cell-compartment attribution.

**Limitations.** The grid covers three diseases and two exposure resources, and
the second HCC power level is drawn from the same resource as the first, so it
shows the pattern is not peculiar to melanoma or to one small eQTL dataset — not
that it holds for any disease or any resource. One of its six cells is void on
its own mismatched-list control, one of the remaining five does not reach
P < 0.05, and that same cell does not clear the density-matched null either. The
one non-cancer cell that survives does so only at the registered window: across
seven window settings its mismatched-list control fails at six, and fails harder
as the window narrows. We therefore claim generality across tumours, and only
descriptive attribution outside them. A stronger non-cancer claim would need a
disease that does not share immune architecture with the comparator, or several
pre-selected unrelated reference lists in place of a single mismatched one. Three diseases is also not a sample from which generality can be
estimated: they were chosen for the contrasts they provide, not at random. The
power-recovery curves rest on reference lists of four and six genes and compare
against a full-power list that is itself unstable, so they measure agreement
between two imperfect lists rather than recovery of truth, and are not
extrapolated above observed power. The two exposure resources that carry the
attribution finding are both European-ancestry, and the second is whole blood, so
a positive result there cannot be read back as CD4-specific; the three additional
resources we tried, which include Peruvian and African-ancestry donors, were
underpowered or failed their control, so they neither extend the finding nor
broaden its ancestry. That limit is structural rather than incidental. Across the
whole eQTL Catalogue [33], stimulation resources profiling several conditions at
adequate donor numbers exist only for myeloid cells — one macrophage resource
covers 24 conditions at n = 199 — while for CD4⁺ T cells the only multi-condition
stimulation resource is the one used here, and every other stimulated CD4⁺
dataset carries a single condition. **The dynamic axis of the framework this
paper audits is currently supported by one resource**, which bounds how far any
audit of it, including ours, can go. The patient analyses rest on three small cohorts and the strongest
of them is the discovery cohort. The genotype × pseudotime test covers lead
variants under linear and quadratic models, so a non-monotonic effect at a
non-lead variant is not excluded. Two questions this design cannot settle need new
data rather than new analysis: whether the patient stratification is real, which
requires a same-disease, same-regimen, CD4-resolved cohort that does not exist
publicly; and whether the state has immunological consequences, which requires
perturbation.

**Conclusions.** A candidate list produced by this framework at current outcome
power should be reported as the set of genes that passed screening under these
conditions, not as targets — and because the fragile part of such a list is
precisely its novel-locus part, the headline of a study of this kind is the part
least likely to replicate. The eight checks above cost little and would have
changed what this analysis reported at nearly every stage.

---

## Figures

Fig. 1 Locus attribution under both outcomes, by independent locus ·
Fig. 2 Generality: five nested releases, the transfer test, the second disease,
the non-cancer outcome, the second exposure resource, and the crossed grid with
its mismatched-locus control · Fig. 3 Power and list stability by locus class, and
the same stratification against full-power |z| · Fig. 4 Colocalisation versus
SMR/HEIDI, with in-sample fine-mapping · Fig. 5 Instrument availability across the
pathway at three levels · Fig. 6 The activation window in which TPI1 is
instrumentable: effect size against precision across the eight profiles ·
Fig. 7 The CD4⁺ metabolic axis and its chromatin signature ·
Fig. 8 Compartment attribution · Fig. 9 Patients across three cohorts.

Figure numbering follows this version's reading order and is independent of the
numbering used in the full reference document. The self-administered attribution
check is reported in the text and tabulated in Supplementary S27; it carries no
main figure.

## Methods

### Exposure data

Cis-eQTL summary statistics came from a published activation time-course of
primary human CD4⁺ T cells [4], comprising eight profiles — naive and memory
cells at 0 h, 16 h, 40 h and 5 d — from 85–100 donors each (naive 99/99/89/85;
memory 100/95/89/90). Per-profile sample sizes were recovered from
`ma_count/(2·MAF)` and were internally consistent for 100% of variants in all
eight profiles. Expression had been inverse-normal transformed, so exposure units
are standard deviations. Coordinates are GRCh38 and were used as distributed; no
liftover was performed. The second exposure resource was the eQTLGen Consortium
whole-blood cis-eQTL release (n = 31,684) [21].

### Outcome GWAS

The primary outcome was a fixed-effect inverse-variance meta-analysis of FinnGen
R12 `C3_MELANOMA_SKIN_EXALLC` (5,753 cases / 378,749 controls; the OpenGWAS
`finn-b-` snapshot was not used, being an R5 extract with 98 cases) and Rashkin
et al. (GCST90011809, 6,777 cases) [5,6], totalling 12,530 cases and 789,099
controls (N = 801,629, case fraction 0.015631). Rashkin reports odds ratios and P
values without standard errors, so standard errors were reconstructed as
se = |log OR|/|Φ⁻¹(P/2)|. For a two-study meta-analysis Cochran's Q has one
degree of freedom, so the upper tail is erfc(√(Q/2)); the exp(−Q/2) form used
initially gave a significant-Q rate of 1.59% against 5.33% for the correct form.
1,740 variants underflowed to P = 0, for which log-space `mlogp` was computed and
P floored at 10⁻³⁰⁰. Landi et al. (GCST010304) supplied only 76 loci and was used
solely as the known-locus reference [7]. Further outcomes were FinnGen R12 lung,
colorectal, pancreatic, breast and prostate cancer; FinnGen R12
`C3_HEPATOCELLU_CARC_EXALLC` (947 cases) [22] and GCST90809296, the European arm of an
eleven-cohort hepatocellular meta-analysis (3,748 cases / 1,861,536 controls)
[20]; and FinnGen R13 `M13_RHEUMA` (16,775 cases / 308,847 controls). Effective
sample sizes, 4/(1/N_case + 1/N_control), are 49,337 for the melanoma meta,
14,962 and 3,779 for the two hepatocellular levels, and 63,643 for rheumatoid
arthritis.

### Instrument selection, harmonisation and MR

Instruments were the top cis-eQTL per gene × profile at P < 5×10⁻⁸, matched to
the outcome by chr:pos with allele consistency; 26 variants were discarded for
ambiguous indel representation. The strict set is 3,556 records over 2,141 unique
variants for the meta round. The minimum F statistic was 22.2, so the F > 10
filter was never binding. Each exposure carries one instrument, so the primary
estimator is the Wald ratio with first-order delta-method standard errors. The
analytic consequence is stated because it constrains interpretation: with
b = β_out/β_exp and SE = se_out/|β_exp|, the statistic z = β_out/se_out does not
involve the exposure, so a variant returns the same P value in every profile in
which it is the lead eQTL. Profile specificity is therefore reported as the
profile containing the instrument, never as a difference in P value.

A relaxed set (F ≥ 5, P < 0.05, LD-clumped with PLINK 2 at
`--clump-r2 0.1 --clump-kb 1000` against 525 unrelated 1000 Genomes GRCh38 EUR
samples [23,28]) supported IVW, multiplicative-random-effects IVW and weighted
median sensitivity analyses [27]. Weighted mode was evaluated and excluded: at a
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
independent locus (1 Mb single-linkage), each collapsed by minimum-p and by Simes
combination, plus a two-stage hierarchical procedure selecting genes by Simes and
then records within selected genes. All seven return a superset of the
record-level gene list. The audit's own conclusions are stated in independent
loci.

### Colocalisation, SMR and HEIDI

Colocalisation used `coloc.abf` with default priors over the full cis window
[24], the eQTL as a quantitative trait and the disease as case-control, matching
on a chr:pos:alleles key. The published criterion PP.H4/(PP.H3+PP.H4) > 0.7 is
too permissive here, since the ratio can rest on two vanishingly small
posteriors, so PP.H3+PP.H4 > 0.5 was additionally required. Every result carries
two diagnostics: the number of variants in the region, and the physical distance
between the GWAS and eQTL peaks. Sensitivity across four windows (±100 kb to the
full cis region) and four values of the shared-causal-variant prior (p₁₂ from
10⁻⁶ to 5×10⁻⁵) is reported for every locus. Because `coloc.abf` permits at most
one causal variant per trait per region, the assumption was tested on the outcome
side using FinnGen's published in-sample SuSiE fine-mapping and `susie_rss` run
on 1 Mb windows with a proxy LD matrix from the same 525-sample panel [23,25], with
`estimate_s_rss` as a mismatch diagnostic. The eQTL side was not fine-mapped: at
85–100 donors without in-sample LD, SuSiE is not reliable there, and this is a
stated boundary. Our proxy-LD procedure recovers 19 credible sets at MC1R against
3 from in-sample fine-mapping, so its counts are read in one direction only —
exactly one credible set is conservative evidence of a single signal, several is
not evidence of multiple signals.

SMR v1.3.1 [26] was run against the same reference panel (94,979 within-cis
variants) with GWAS statistics trimmed to cis windows and de-duplicated for
multi-allelic sites. The `--diff-freq` QC removed 0.45% of variants at the
default threshold, which was not relaxed. The number of variants entering each
HEIDI test is reported with the result, as HEIDI's power depends on it (range
8–20 here).

### Locus annotation and enrichment

The melanoma reference set was the union of lead variants from the three Landi
et al. 2020 GWAS [7] retrieved from the GWAS Catalog REST API [34] and mapped to
GRCh38: 157 loci. Hepatocellular loci were lead variants for EFO_0000182 and
liver cancer at P < 5×10⁻⁸ together with those tabulated in the outcome
publication — 83 rsIDs, 73 placeable on GRCh38. The rheumatoid arthritis
reference was the 87 lead rsIDs of Okada et al. 2014 [47] (GCST002318), published
before FinnGen existed and from cohorts that exclude it; 83 were placeable.
Independent loci are a **non-recursive fixed-anchor partition** of instrument
positions within a chromosome: the first unassigned variant becomes an anchor and
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
disease's significant loci against another disease's reference. A **strength-and-
frequency** permutation matches significant to background loci on eQTL-p decile
crossed with outcome allele-frequency quintile, with 10,000 stratified resamples
and a one-sided empirical P. A **density** permutation follows the estimand
frozen in Supplementary S38 before it was run: background loci are matched to
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
GSE78220, GSE91061, GSE235863) [9–13,19], the spatial dataset [14], the purified
CD4⁺ multiome (GSE282266) [15] and the peripheral-blood variance decomposition
(GSE199994) [16] are given in Supplementary S12; the points that changed a result
are summarised here. Donors or samples, never cells, are the statistical unit for
inference. Module scores used Seurat `AddModuleScore` [29] with expression-binned
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
hg38 [30,31], comparing foreground with background peaks matched on GC content
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

### Pre-specification and stopping rules

The following were fixed before the analyses they govern and are reported because
they determined which results are presented. Every test carries a positive
control, and a test whose positive control fails is discarded rather than
interpreted — the within-lymphoid-subset spatial test (HLA-C P = 0.23) and the
pre-specified immune-module analysis were discarded under this rule. Kill
criteria, test direction and window sizes were written into each script before it
ran. A negative result was not followed by a search for a positive one in the
same data. Simulations require calibration against a quantity not used to build
them. Twelve pre-registration documents (Supplementary S9, S18, S20, S21, S22,
S24, S25, S28, S29, S30, S32, S33) each carry a reading table written before the
run, and only their results register was completed afterwards; deviations are
logged in the document they belong to.

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
intermediate result tables and the fourteen pre-registration documents will be
deposited at ⟨repository DOI⟩, comprising every numbered analysis and figure
script together with the tables needed to reproduce each figure and every number
reported in the text. CQTNA, a runnable implementation of the diagnostics in the
Discussion, is included in the deposit.

## Supplementary information

S9–S38, including the fourteen pre-registration documents with their reading tables
and results registers; the multiple-testing unit sensitivity analysis; the
self-administered attribution check; the complete record of target-substantiation
attempts with the selection denominators and every stopping-rule instance; the
technical account of the two processing errors; the replication-cohort search;
the colocalisation window-and-prior sensitivity analysis; and the literature
audit.

## References

1. Zheng J, Yang Q, Liu H, et al. Integrating single-cell transcriptome-wide Mendelian randomization and differentially expressed gene analyses to prioritize dynamic immune-related drug targets for cancers. *Adv Sci* 2025;12:e07451. doi:10.1002/advs.202507451
2. Wu X, Ying H, Yang Q, et al. Transcriptome-wide Mendelian randomization during CD4⁺ T cell activation reveals immune-related drug targets for cardiometabolic diseases. *Nat Commun* 2024;15:9302. doi:10.1038/s41467-024-53621-7
3. Cui K, Zou Q, Qu X, et al. Transcriptome-wide Mendelian randomization and single-cell analysis during CD4⁺ T cell activation deciphers immunotherapeutic targets for colorectal cancer. *npj Precis Oncol* 2025;10:32. doi:10.1038/s41698-025-01236-6
4. Soskic B, Cano-Gamez K, Smyth DJ, et al. Immune disease risk variants regulate gene expression dynamics during CD4⁺ T cell activation. *Nat Genet* 2022;54:817–826. doi:10.1038/s41588-022-01066-3
5. Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature* 2023;613:508–518. doi:10.1038/s41586-022-05473-8 (Author Correction: *Nature* 2023;615:E19. doi:10.1038/s41586-023-05837-8)
6. Rashkin SR, Graff RE, Kachuri L, et al. Pan-cancer study detects genetic risk variants and shared genetic basis in two large cohorts. *Nat Commun* 2020;11:4423. doi:10.1038/s41467-020-18246-6 (GCST90011809)
7. Landi MT, Bishop DT, MacGregor S, et al. Genome-wide association meta-analyses combining multiple risk phenotypes provide insights into the genetic architecture of cutaneous melanoma susceptibility. *Nat Genet* 2020;52:494–504. doi:10.1038/s41588-020-0611-8 (GCST010302–GCST010304)
8. Nathan A, Asgari S, Ishigaki K, et al. Single-cell eQTL models reveal dynamic T cell state dependence of disease loci. *Nature* 2022;606:120–128. doi:10.1038/s41586-022-04713-1
9. Sade-Feldman M, Yizhak K, Bjorgaard SL, et al. Defining T cell states associated with response to checkpoint immunotherapy in melanoma. *Cell* 2018;175:998–1013.e20. doi:10.1016/j.cell.2018.10.038 (GSE120575)
10. Jerby-Arnon L, Shah P, Cuoco MS, et al. A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. *Cell* 2018;175:984–997.e24. doi:10.1016/j.cell.2018.09.006 (GSE115978)
11. Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 2016;352:189–196. doi:10.1126/science.aad0501 (GSE72056)
12. Hugo W, Zaretsky JM, Sun L, et al. Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. *Cell* 2016;165:35–44. doi:10.1016/j.cell.2016.02.065 (GSE78220)
13. Riaz N, Havel JJ, Makarov V, et al. Tumor and microenvironment evolution during immunotherapy with nivolumab. *Cell* 2017;171:934–949.e16. doi:10.1016/j.cell.2017.09.028 (GSE91061)
14. Thrane K, Eriksson H, Maaskola J, et al. Spatially resolved transcriptomics enables dissection of genetic heterogeneity in stage III cutaneous malignant melanoma. *Cancer Res* 2018;78:5970–5979. doi:10.1158/0008-5472.CAN-18-0747
15. Katko A, Potter SJ, et al. Gene regulatory network determinants of rapid recall in human memory CD4⁺ T cells. *Cell Rep* 2026;45:117103. doi:10.1016/j.celrep.2026.117103 (GSE282266)
16. Boukhaled GM, Gadalla R, et al. Pre-encoded responsiveness to type I interferon in the peripheral immune system defines outcome of PD1 blockade therapy. *Nat Immunol* 2022;23:1273–1283. doi:10.1038/s41590-022-01262-7 (GSE199994)
17. Virós A, et al. Spatial transcriptomics of primary cutaneous melanoma. NCBI GEO GSE316760 (submitted 2026-01-16); no associated publication indexed at the time of writing.
18. Pham F, Dufeu M, Benboubker V, et al. Spatial tumour-immune ecosystems shape the efficacy of anti-PD1 immunotherapy in primary cutaneous melanoma. NCBI GEO GSE300445 (submitted 2025-06-23); no associated publication indexed at the time of writing.
19. Guo X, et al. Contrasting cytotoxic and regulatory T cell responses underlying distinct clinical outcomes to anti-PD-1 plus lenvatinib therapy in hepatocellular carcinoma. *Cancer Cell* 2025;43:248–268.e9. doi:10.1016/j.ccell.2025.01.001 (GSE235863)
20. Ghouse J, Gellert-Kristensen H, O'Rourke CJ, et al. Genome-wide meta-analysis identifies nine loci associated with higher risk of hepatocellular carcinoma. *JHEP Rep* 2025;7:101485. doi:10.1016/j.jhepr.2025.101485 (GCST90809296)
21. Võsa U, Claringbould A, Westra H-J, et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. *Nat Genet* 2021;53:1300–1310. doi:10.1038/s41588-021-00913-z (eQTLGen, 2019-12-11 cis-eQTL release, n = 31,684)
22. FinnGen. Release R12 (2024), endpoints `C3_MELANOMA_SKIN_EXALLC` and `C3_HEPATOCELLU_CARC_EXALLC`; Release R13, endpoints `C3_MELANOMA_SKIN_WIDE` and `M13_RHEUMA`. https://www.finngen.fi/en/access_results
23. Byrska-Bishop M, Evani US, Zhao X, et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell* 2022;185:3426–3440.e19. doi:10.1016/j.cell.2022.08.004
24. Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genet* 2014;10:e1004383. doi:10.1371/journal.pgen.1004383
25. Wang G, Sarkar A, Carbonetto P, Stephens M. A simple new approach to variable selection in regression, with application to genetic fine mapping. *J R Stat Soc Series B* 2020;82:1273–1300. doi:10.1111/rssb.12388
26. Zhu Z, Zhang F, Hu H, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. *Nat Genet* 2016;48:481–487. doi:10.1038/ng.3538
27. Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal inference across the human phenome. *eLife* 2018;7:e34408. doi:10.7554/eLife.34408
28. Chang CC, Chow CC, Tellier LC, et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 2015;4:s13742-015-0047-8. doi:10.1186/s13742-015-0047-8
29. Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat Biotechnol* 2024;42:293–304. doi:10.1038/s41587-023-01767-y
30. Fornes O, Castro-Mondragon JA, Khan A, et al. JASPAR 2020: update of the open-access database of transcription factor binding profiles. *Nucleic Acids Res* 2020;48:D87–D92. doi:10.1093/nar/gkz1001
31. Schep AN, Wu B, Buenrostro JD, Greenleaf WJ. chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nat Methods* 2017;14:975–978. doi:10.1038/nmeth.4401
32. Schmiedel BJ, Singh D, Madrigal A, et al. Impact of genetic polymorphisms on human immune cell gene expression. *Cell* 2018;175:1701–1715.e16. doi:10.1016/j.cell.2018.10.022 (DICE)
33. Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human gene expression and splicing quantitative trait loci. *Nat Genet* 2021;53:1290–1299. doi:10.1038/s41588-021-00924-w (eQTL Catalogue)
34. Sollis E, Mosaku A, Abid A, et al. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. *Nucleic Acids Res* 2023;51:D977–D985. doi:10.1093/nar/gkac1010
35. Ochoa D, Hercules A, Carmona M, et al. The next-generation Open Targets Platform: reimagined, redesigned, rebuilt. *Nucleic Acids Res* 2023;51:D1353–D1359. doi:10.1093/nar/gkac1046
36. Sun D, Wang J, Han Y, et al. TISCH: a comprehensive web resource enabling interactive single-cell transcriptome visualization of tumor microenvironment. *Nucleic Acids Res* 2021;49:D1420–D1430. doi:10.1093/nar/gkaa1020
37. Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nat Biotechnol* 2020;38:675–678. doi:10.1038/s41587-020-0546-8
38. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell* 2018;173:400–416.e11. doi:10.1016/j.cell.2018.02.052
39. Tsherniak A, Vazquez F, Montgomery PG, et al. Defining a cancer dependency map. *Cell* 2017;170:564–576.e16. doi:10.1016/j.cell.2017.06.010 (DepMap, Broad Institute)
40. Oughtred R, Rust J, Chang C, et al. The BioGRID database: a comprehensive biomedical resource of curated protein, genetic, and chemical interactions. *Protein Sci* 2021;30:187–200. doi:10.1002/pro.3978 (BioGRID ORCS)
41. Reales G, et al. Design and interpretation of eQTL–GWAS colocalisation studies: lessons from a large-scale evaluation. *PLoS Genet* 2026. doi:10.1371/journal.pgen.1012141
42. Tambets R, et al. Extensive co-regulation of neighboring genes complicates the use of eQTLs in target gene prioritization. *Hum Genet Genomics Adv* 2024;5:100348. doi:10.1016/j.xhgg.2024.100348
43. Rosen JD, et al. Higher eQTL power reveals signals that boost GWAS colocalization. *Am J Hum Genet* 2026. doi:10.1016/j.ajhg.2026.02.009
44. Howe LJ, et al. Evaluating transportability of in vitro cellular models to in vivo human phenotypes using gene perturbation data. *Nat Commun* 2025. doi:10.1038/s41467-025-67199-1
45. Lin Z, Pan W. A robust cis-Mendelian randomization method with application to drug target discovery. *Nat Commun* 2024. doi:10.1038/s41467-024-50385-y
46. Karhunen V, et al. Integrating genetic data with biological insight: a practical guide to cis-Mendelian randomization. *Am J Hum Genet* 2026. doi:10.1016/j.ajhg.2026.03.011
47. Okada Y, Wu D, Trynka G, et al. Genetics of rheumatoid arthritis contributes to biology and drug discovery. *Nature* 2014;506(7488):376–381. doi:10.1038/nature12873 (GWAS Catalog GCST002318)
