# Outcome GWAS architecture shapes target nomination from context-specific eQTLs: an audit across two diseases and two exposure resources

<!-- Genome Biology format: structured abstract (~350 w); main text 6,000-8,000.
     The failure history and revision record, every stopping-rule instance, the
     technical account of our own processing errors, the candidate-selection
     timeline and repeated claim boundaries are in Supplementary, not here. The
     eight diagnostics are the backbone. Every number is identical to
     MANUSCRIPT_v2_dual_thread.md, the single reference document; this version is
     a pure derivative with nothing of its own. -->

## Abstract

**Background.** Context-specific expression quantitative trait loci (eQTLs)
combined with Mendelian randomization (MR) are widely used to nominate immune
targets in cancer. An instrument selected in one cell state at one timepoint is
usually a single variant, so the Wald ratio has test statistic z = β_out/se_out
and the outcome GWAS supplies every MR significance claim. We audited what that
implies in practice, using CD4⁺ T cell cis-eQTLs from eight activation profiles
against a 12,530-case melanoma meta-analysis, then varied outcome and exposure
resource independently.

**Results.** Significant signal concentrates on loci already known for the
outcome: 4.09-fold by independent locus. The pattern recurs across five nested
power levels of one GWAS resource, where no novel-locus gene reaches significance
at any case number between 2,705 and 5,753, and it transfers unchanged to that
resource's current release although not one of the 3,434 test statistics does; in
a second disease, where three of four significant loci are known hepatocellular
carcinoma loci including *PNPLA3*; and when the exposure resource is replaced by
a whole-blood eQTL dataset 300-fold larger, which raises significant loci from 7
to 30 and leaves enrichment at 4.44-fold (P = 3.7×10⁻¹¹). All six
disease-by-resource combinations enrich, against a mismatched-locus control that
does not. Replacing the outcome with a higher-powered meta-analysis — which
changes study composition as well as power — increased MR discoveries while
lowering colocalisation support and replaced the candidate list entirely: two
lists from identical exposure data share no genes. Down-sampling shows the loss
falls unevenly by effect size and hence by locus class: at half power,
known-locus genes recover 85.8% of the time against 22.8%, of which 68–80% is
reproduced by a model using full-power |z| alone. Of 28 glycolytic genes, 3 are
instrumentable, 2 analysable against this outcome and 1 yields a nominal
association; that gene is 6.7-fold higher in malignant cells than in CD4⁺ T
cells, so tissue-level validation of it measures tumour. At most 7.9% of 152
comparable studies perform the locus-attribution check.

**Conclusions.** In this analysis, and within the power range we could observe,
the reproducible part of a candidate list produced by this framework was the part
that did not constitute a discovery. **Which loci a nomination lands on is set by
the outcome GWAS; how many instruments exist to be tested is set by the exposure
resource.** What survives both is a statement about the state in which a pathway's
regulation can be measured precisely enough to yield an instrument at all, and
eight inexpensive checks follow directly.

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
GWAS as outcomes — and report novel immune targets on that basis.

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

### Design and positive controls

Eight CD4⁺ T cell activation profiles (naive and memory, 0 h to 5 d) provided
exposures. Instruments were selected at eQTL P < 5×10⁻⁸ with F > 10, one variant
per exposure, giving 3,556 records after harmonisation with the outcome. Every
pre-specified positive control passed: all 157 known melanoma loci were recovered
and 136 strengthened under the meta outcome; genome-wide significant variants rose
from 2,871 to 4,552; PARP1 reproduced its published direction; and the MC1R region
gave the strongest associations in the study, at P = 4×10⁻³⁷ (Fig. 1).

One property of the design governs everything that follows. With one instrument
per exposure, the Wald z equals β_out/se_out, so **given the instrument set**, MR
significance is a property of the outcome GWAS; the exposure data determine which
variants are asked about and how many there are, but not which of them answers. Steiger filtering cannot arbitrate direction here either, since cis-eQTL
and disease R² differ by orders of magnitude and the test passes almost by
construction. We therefore state at the outset what the design can and cannot
decide, and test the consequences empirically rather than asserting them.

The same constraint creates a units problem that must be settled before any count
below is read. Those 3,556 records carry only **2,126 unique variants and 1,195
unique genes**, because one variant can be the lead eQTL for a gene in several
profiles and for more than one gene, so the testing family is records while the
list is reported by gene and the attribution by independent locus. **The record
level is primary**, for the single reason that every analysis here was built on
it — the registered predictions, the release trajectory, both generalisations,
the grid. Choosing a unit now, after seeing which one yields the longest list, is
the selective emphasis this paper exists to detect. We instead recomputed the
FDR < 0.05 list under all seven alternatives — variant, gene and locus, each by
minimum-p and by Simes, plus a two-stage hierarchical procedure (Supplementary
S26). **The record level is the most conservative of them**: no gene is lost
under any other unit and up to 43 are added, so every claim here is a subset of
what a looser unit would license. The quantities the argument runs on are stable
— 7 to 9 independent loci against the 7 reported, 2 under every unit in the
FinnGen round — and the attribution result does not depend on the unit at all
(known-locus share 42.9–55.6% against 42.9% reported; 100% at every unit in
FinnGen). One row reads as a finding rather than a check: testing at locus level
leaves the locus count almost unchanged but inflates the gene list from 10 to 53,
of which 22 are in the MHC and 7 in the chr17q21.31 inversion, because a
significant locus does not name a gene.

### Significant signal sits on loci already known for the outcome

In the FinnGen round, all ten FDR-significant records fell within 1 Mb of a known
pigmentation or naevus locus — VPS9D1-AS1 50 kb from MC1R, CDK10 48 kb from the
MC1R R151C variant, PARP1 5–14 kb from PARP1 — with no immune signal at
FDR < 0.05 at all. Under the meta outcome, counting independent loci rather than
gene records, 3 of 7 significant loci were known-locus loci: a 4.09-fold
enrichment over the 10.5% background (one-sided P = 0.028; Fig. 1). Applied
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
locus attribution is 9.6-fold (P = 0.011). The stability is not an artefact of
reusing data. Not one of the 3,434 test statistics is unchanged between the two
releases, standard errors shrink by 3.85% against the 3.8% the added cases
predict, and 51 of the 275 nominally significant records are replaced. What holds still is the
significant tier; the tier below it churns. Because the extra cases and the
broadened control exclusion both push towards more signal, a successful transfer
here is partly confounded and we do not read it as power alone.

**A second disease.** Repeating the nomination against hepatocellular carcinoma
with the exposure side unaltered, at two outcome power levels (3,748 and 947
cases), gives four significant loci across both levels, three carrying a known
HCC lead SNP: SAMM50, 26 kb from rs2294915 in the *PNPLA3* region, and ZNF506 in
the chr19 *TM6SF2* region at both levels. Enrichment is 17.6-fold in the
lower-powered outcome (P = 0.0031) and 8.9-fold in the higher-powered one, where
it does not reach significance (P = 0.110). The single novel nomination, SUPV3L1,
appears only at high power and has FDR = 0.98 at low power. The higher-powered
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
attribution does not move: **20 of 30 significant loci (66.7%) carry a known
melanoma, naevus or pigmentation lead SNP, 4.44-fold over this resource's own
15.0% background (P = 3.7×10⁻¹¹)**, against 4.09-fold in CD4⁺ T cells. A
matched-background version agrees (4.14–4.48-fold, empirical P = 10⁻⁴). The two
resources' novel nominations intersect in exactly one gene, ZFYVE19, so the
stronger claim that novel nominations never reproduce is not supported.

**Both axes crossed.** Scoring every disease-by-resource combination against its
own disease's known loci, all six enrich: 4.09 and 4.44 for melanoma, 8.85 and
5.11 for HCC-high, 17.6 and 8.67 for HCC-low. Three of six reach P < 0.05, on
numerators as small as one of two loci, so this is consistency of direction rather
than six independent significant tests. A pre-registered mismatched-locus control
makes it interpretable: scoring the same combinations against the wrong disease's
list collapses the enrichment — HCC's list on melanoma gives 1.30-fold (P = 0.41),
melanoma's list on HCC gives 0.00-fold (P = 1.0) — so the effect is specific to
each outcome's own genetics and is not an artefact of locus density.

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
gap between eQTL and GWAS colocalisation [37]. Our design asks the mirror
question — with the exposure held fixed, does a *higher-powered outcome* make
target nomination more reliable? — and the answer here is that it does not: the
meta-analysis produced more MR discoveries, lower colocalisation support, and a
candidate list with no genes in common with the previous one. That outcome
differs from its predecessor in study composition as well as power, so this is
not a power experiment either; only the nested-release series isolates power, and
there the significant list is stable rather than replaced. The two axes are not
interchangeable, and a study underpowered on one cannot be rescued by the other.
Reporting which axis a claim rests on should be routine.

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
per test (Fig. 3). Evidence tiers that accept MR plus SMR without colocalisation
would have reported MC1R linkage spillover as CD4-mediated immune targets.

### More outcome power moves the two methods in opposite directions

Meta-analysis raised MR discoveries from 10 to 21 records while *lowering*
colocalisation support: across the 127 exposures run in both rounds, median
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
with the full-power list (Fig. 4). The simulation is calibrated empirically — not
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
there, so the window is a statement about precision rather than amplitude
(Fig. 6a,b). The published genotype × pseudotime interaction test settles the
question directly: **none of the three TPI1 lead variants shows an interaction** —
in memory cells, linear P = 0.974 and quadratic P = 0.992; in naive cells, 0.707
and 0.891 for the variant 1,063 bp from our instrument, and 0.378 and 0.605 for
the other — while the gene's *expression* is among the most strongly
pseudotime-dependent in the dataset (Moran's I = 0.664).

Three things are separable here, and the interaction test tells them apart:
dynamic expression (yes, strongly), time-specific eQTL significance (yes,
genome-wide significant at 16 h only), and dynamic genetic effect (no evidence, in
any of three tests). The design consequence holds regardless: a study using a
resting-state resource would find nothing here and could not distinguish that
from a true null.

### The programme the instrumented gene belongs to is a definable cell state

Because MR names one gene of a pathway, we asked what that gene marks. In
purified CD4⁺ T cell multiome data, residualising expression on activation
intensity leaves an axis enriched 21.3-fold for glycolysis and 10.3-fold for
ribosomal proteins, while **none of the eleven biological modules tested exceeds
a composition-matched null**, so the axis is a definable anabolic state rather
than a restatement of activation strength; its chromatin signature is concordant,
with AP-1 motif enrichment that persists in activation-invariant peaks (Fig. 9).
This characterises what the nominated gene co-varies with and not that its
variant causes the state — a separate test of that, whether the instrument
disrupts an AP-1 motif, returned an empirical P = 1.0.

### Where the effect is, when MR cannot say

In melanoma single-cell data with all cell types annotated, TPI1 is far higher in
malignant cells than in CD4⁺ T cells: +2.75 log₂ units (≈6.7-fold) in 16 of 16
patients (P = 3×10⁻⁵), detected in 95.3% of malignant cells versus 60.6% of CD4⁺
T cells. The difference survives cell-level matching on sequencing depth (+2.03,
10 of 11 patients) and replicates in an independent cohort (+1.28, 11 of 11,
P = 1×10⁻³), with pre-specified positive controls passing in both, and is
concordant in direction in a third cohort. The same holds for the locked
glycolytic signature with and without TPI1 (Fig. 7).

The arithmetic consequence is what matters. With those ratios and the cell-type
proportions typical of melanoma tissue, CD4⁺ T cells contribute on the order of
one to two per cent of the tissue-level TPI1 signal. **A tissue-level measurement
of this gene is a measurement of tumour glycolysis**, whatever its P value. This
bounds what bulk validation of such a nomination can establish, and the check that
would have caught it costs minutes.

### Does the pipeline attribute the right gene where the answer is known?

The criticism this paper makes of nomination applies to this paper, and can be
administered to it. At melanoma loci where a causal gene is generally accepted, we
asked whether the pipeline's own significant nomination names that gene. Of ten
such loci reached under the larger exposure resource, **six name the accepted gene
and four do not**, and the failures are the informative ones. At the MC1R region —
the strongest melanoma locus in the genome — the nomination spans **sixteen genes
and MC1R is not among them**. At the OCA2/HERC2 locus it names the pseudogene
HERC2P9 rather than HERC2. At *TYR* it names ODF3, and at the CDKN2A/MTAP locus
C9orf66. Where it succeeds it often does so cleanly: IRF4 and MX2 are each named
alone and correctly.

This is the co-regulation problem of Tambets et al. [36] observed at the top of
the effect-size distribution, and it bounds the whole design: a framework that
misassigns the gene at the best-characterised locus in its own disease should not
be read as assigning genes at uncharacterised ones. It is also why we report
compartment attribution separately from gene attribution — the two fail
independently.

### The functional layer has its own failure modes

Splitting cells by a score also splits them by lineage purity. When cells were
divided by glycolytic score, apparent differences in helper and regulatory
phenotype collapsed after matching on lineage composition at the cell level, and
14 of 28 module comparisons were concordant — exactly chance. A cluster-level
control passed while the cell-level control failed, so the level at which a
purity control is applied determines whether it works. Two processing errors of
our own, detected by these controls, are described in Supplementary S12; they are
instances of the same failure mode.

The perturbation layer is bounded in the same way. TPI1 is a hit in 628 of the
**1,471 human CRISPR screens** in BioGRID ORCS that measured it (**42.7%**), a
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
(Δ = +0.874, P = 0.043) (Fig. 8). The per-arm estimates, the three treatments of
repeated patients and the sign dependence on regulatory T cells are Supplementary
S19 and S21.

### These checks are rare in current practice

In a search-defined sample of 152 eQTL-instrumented MR target-nomination papers
whose full text we could obtain, **at most 12 (7.9%) compare their significant
signal against previously reported loci for their own outcome trait, and 1 (0.7%)
reports how the candidate list depends on the outcome GWAS used**. About 59%
performed colocalisation and 47% used SMR or HEIDI. This counts phrases in full
text rather than quality of practice, and not performing a check is not evidence
that a study's conclusions are wrong.

---

## Discussion

Everything in this framework that depends on the outcome GWAS proved unstable, and
the instability has a specific shape: signal lands on loci the outcome already
knows about, the novel part of the list turns over completely between outcomes,
and more power moves MR and colocalisation in opposite directions. The pattern
held when we changed the disease, when we changed the exposure resource by two
orders of magnitude, and across all six combinations of the two — against a
mismatched-locus control that showed nothing. What did not depend on the outcome
was a smaller and more durable statement: in which cell state a pathway's
regulation can be measured precisely enough to yield an instrument at all.

Of eleven target-level tests that returned a verdict, seven overturned the claim
under examination, three were inconclusive and one was significant on discovery
but not confirmed on transfer; five genes were successively designated lead
candidate under a ranking criterion that was never fixed in advance, and the
first four designations were overturned. That is this project's own history and
not an estimate of how often the framework fails, and the itemised record — every
attempt, stopping rule and withdrawal, the candidate-selection timeline, the
eight pre-registrations and the technical account of two processing errors of
ours — is Supplementary S12.

Eight checks follow directly, each cheap and each capable of changing what a study
of this kind reports.

**(i) Annotate the signal against the outcome's own known loci, counting
independent loci rather than gene records.** Every result in this paper's first
half rests on that distinction, and at most 7.9% of comparable studies apply it.

**(ii) Require colocalisation, preferably with explicit multiple-signal modelling
and an LD reference matched to the outcome cohort, rather than treating SMR/HEIDI
as sufficient — and report how many variants entered each HEIDI test.** HEIDI
passed 253 of 291 records that colocalisation assigned to distinct causal
variants, on 8–20 variants per test. We stop short of naming any single
colocalisation method the final arbiter: ours permits at most one causal variant
per region, and our attempt to go beyond that assumption over-split MC1R relative
to in-sample fine-mapping, so it could constrain one locus in one direction only.

**(iii) Do not expect more outcome power to resolve the disagreement.** It moves
MR and colocalisation in opposite directions, because one depends on a single
variant's z and the other on the shape of the regional signal.

**(iv) Report the physical distance between GWAS and eQTL peaks alongside
posteriors, and test whether posterior movement survives explicit modelling of
multiple causal signals before calling it a change in resolution.** For PARP1 the
collapse is not an unmodelled second signal; for two other loci no credible set
exists at either power, so those regions cannot adjudicate it.

**(v) Compute the nominated gene's cell-type expression ratio in an annotated
atlas before citing tissue-level data as validation.** This costs minutes, and
here it shows that tissue-level measurements of the nominated gene are
measurements of tumour.

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
model proves only partly transportable [38]:

| Claim about TPI1 | Status here |
|---|---|
| Causal target for melanoma | **Not supported** — FDR = 0.119, prior-dependent colocalisation, region not fine-mappable |
| Dynamic *genetic* effect across activation | **Not supported** — no genotype × pseudotime interaction at any of three lead variants |
| Dynamic *expression* across activation | **Supported** — Moran's I = 0.664 |
| Membership of a definable CD4⁺ metabolic state | **Supported** — axis enrichment 21.3-fold, chromatin concordant |
| Predicts checkpoint-blockade response | **Unstable** — significant on discovery, no arm confirmed twice |
| Functional consequence isolable by perturbation | **No** — a hit in 628 of 1,471 human CRISPR screens (42.7%), a core-essential profile that cannot isolate a CD4-specific role |
| Recognised as a target for this disease elsewhere | **No** — in Open Targets, TPI1's strongest disease associations are triosephosphate isomerase deficiency and neurodegenerative disease; melanoma is not among its leading associations |

Stated flatly: TPI1 is a nominal, prior-sensitive candidate selected by a
multi-layer conjunction: MR P = 1.4×10⁻³
but FDR = 0.119, PP.H4 = 0.51 passing in 8 of 16 window-and-prior combinations, in
a region the outcome GWAS cannot fine-map, with gene attribution at chr12p13
unresolved by genetics. Calling it the gene MR named is accurate; calling it a
discovery is not. Its value is as a demonstration that instrument availability and
biological importance are separable — the strongest enzyme in the functional data
carries no instrument, and the gene that does carry one is principally expressed
by tumour. The complete record of which claims about it were raised, tested and
withdrawn, together with the selection denominators at gene, pathway and mechanism
level and the eight pre-registration documents, is Supplementary S12 and S18–S28.

**Relation to existing guidance.** These eight are additions to, not a
replacement for, current cis-MR practice. Existing guidance already stresses that
cis analyses must be tailored to local biology, that an expression biomarker is
not an intervention, that co-regulation of neighbouring genes can act as
horizontal pleiotropy, and that MR is one strand of triangulation [39,40]. Tambets
et al. showed with approximate ground truth that neighbouring-gene co-regulation
routinely produces colocalisation evidence for several genes at once, so eQTL
data suit candidate *generation* and evidence *composition* better than gene
*attribution* [36]; our MC1R region is the same phenomenon at higher effect size,
and TPI1 at chr12p13 is the same problem unresolved. Reales et al., auditing over
a million colocalisation tests across immune diseases and cell types, likewise
found that gene assignment shifts with platform coverage, cell context and
resource size [35]. What we add is specific to the nomination step and to the
outcome side: locus attribution against the outcome's own known loci, sensitivity
to the outcome release, candidate-list stability, layered instrument visibility,
and cell-compartment attribution.

**Limitations.** The grid covers two diseases and two exposure resources, and the
two HCC power levels are drawn from one resource, so it shows the pattern is not
peculiar to melanoma or to one small eQTL dataset — not that it holds for any
disease or any resource; three of its six cells do not reach P < 0.05. The
power-recovery curves rest on reference lists of four and six genes and compare
against a full-power list that is itself unstable, so they measure agreement
between two imperfect lists rather than recovery of truth, and are not
extrapolated above observed power. Both exposure resources are European-ancestry,
and the second is whole blood, so a positive result there cannot be read back as
CD4-specific. The patient analyses rest on three small cohorts and the strongest
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

## Methods

Full Methods accompany this manuscript.

## Supplementary information

S9–S28, including the eight pre-registration documents with their reading tables
and results registers; the multiple-testing unit sensitivity analysis; the
self-administered attribution check; the complete record of target-substantiation
attempts with the selection denominators and every stopping-rule instance; the
technical account of the two processing errors; the replication-cohort search;
the colocalisation window-and-prior sensitivity analysis; and the literature
audit.
