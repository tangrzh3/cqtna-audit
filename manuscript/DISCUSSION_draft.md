**⚠ 已被 `MANUSCRIPT_v2_dual_thread.md` §4 取代。勿引用本文件的数字。**

# Discussion — draft v1

> 结构说明（中文，不进正文）
> 全文论证骨架是「失败率本身就是结论」。七个方法学发现按递进排列，而不是并列清单。
> TPI1 放在解法演示的位置，且明确标注它靠功能证据的冗余存活、而非遗传统计强度。
> 每个否定都配一个仍然成立的阳性对照，避免读成"阴性结果论文"。

---

## 1. Opening — what the study set out to do and what it actually became

Dynamic eQTL data from activated CD4+ T cells, combined with Mendelian
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

---

## 2. The central argument: the failure rate is the finding

Across the study we made eight attempts to substantiate a target-level claim.
All eight were overturned by checks we ourselves introduced. Critically, these
were not adversarial tests designed to find fault: each was added in the
expectation that it would *support* the claim.

| # | Attempt | Check applied | Outcome | Methodological finding |
|---|---|---|---|---|
| 1 | PADI4 as a novel immune candidate | multi-instrument sensitivity analysis | OR 1.087 (P=1.7×10⁻³) collapsed to IVW 1.033 (P=0.085), weighted median 1.004 (P=0.835) over 7 instruments | single-instrument Wald ratios are not self-validating |
| 2 | FinnGen-round candidate list (PRPSAP2, IMPA1, GCC2, PADI4) | repeat with a higher-powered meta outcome | none survived; **zero overlap** between rounds | ④ |
| 3 | HLA-C acting through CD4 proliferation | replication in an independent scRNA cohort | direction reversed; later traced to a normalisation error (ρ = −0.612/−0.720 → +0.142/−0.120, P ≈ 0.55) | replication is what caught it; the mechanism was an artefact |
| 4 | Bulk ICB cohorts as evidence of a CD4-specific effect | spatial transcriptomics with a cell-type-restricted positive control | tissue-level TPI1 tracks the glycolytic module (ρ = +0.174, 8/8 sections) and is *negatively* correlated with the lymphoid compartment (ρ = −0.080) | ⑤ compartment attribution |
| 5 | A transcriptional state for glycolysis-high CD4 cells | matching on lineage score as well as depth | after matching (610 cells, 28 samples, residual SMD 0.010/0.000) every module went null; the lineage-free effector score gave 14/14 sample concordance, exactly chance | ⑦ |
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

We draw attention to this tally not as a confession but as data. A pipeline that
passes all of its positive controls, and that recovers known melanoma biology
with 5.3-fold enrichment, nevertheless failed to sustain a single target-level
claim across eight independent attempts. Studies using the same framework that
do not report such attempts have not necessarily avoided these failures; they
have not tested for them.

---

## 3. Seven diagnostics, arranged as a progression

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
reason we do not present any single gene as an established target. Five
sequential releases of a single GWAS resource, analysed with everything but the
release held fixed, matched the registered predictions at all twelve comparisons
and showed the known-locus part of the list recovering early and then not moving,
with no novel-locus gene reaching significance at any power in that range. That
experiment strengthens the known-locus half of the claim with real longitudinal
data; it cannot demonstrate turnover among novel candidates, because in that power
range there are none.

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

---

## 4. What survived, and why

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

---

## 5. TPI1: a bounded example

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

---

## 6. Limitations

Beyond those already stated: instruments were single cis-eQTLs, so Steiger
filtering was near-deterministic and is reported as procedure rather than
evidence; relaxed instrument sets at r² < 0.1 retain correlation that inflates
IVW significance (157 FDR-significant exposures versus 21 in the strict
single-instrument set); Rashkin standard errors were reconstructed from odds
ratios and P values, with Cochran's Q significant at 5.33% supporting the
calibration; the single-cell cohorts are small (9 vs 10 donors pre-treatment);
and the CD4 population used for functional work has limited purity (CD4 detected
in 32.5% of cells), which we address by reporting the CD8-negative sensitivity
analysis rather than by redefining the population. Two alternative explanations
for the colocalisation results were addressed directly rather than assumed away,
and neither is fully closed: the multiple-signal explanation is excluded for
PARP1, where it can be tested, but ZFYVE19 and the TPI1 locus carry no resolvable
outcome signal at either power and cannot adjudicate it; and the sequential-release
experiment removes cohort heterogeneity as an explanation for the stability of the
known-locus part of the list, but cannot address the novel part, which is empty
throughout that power range. Finally, the colocalisation supporting our worked
example sits in a region the outcome GWAS cannot fine-map — a limitation we state
rather than absorb.

---

## 7. What should change in practice

For studies applying this framework, we suggest five concrete changes, each
following directly from a finding above: annotate significant signals against
known loci for the trait before interpreting them (①); treat colocalisation as
the primary arbiter and SMR/HEIDI as auxiliary, reporting the number of variants
entering each HEIDI test (②); report the physical distance between GWAS and eQTL
peaks alongside posterior probabilities, and test whether posterior movement
survives explicit modelling of multiple causal signals before interpreting it as
a change in resolution (④); state explicitly what fraction of a
pathway carries instruments before describing a named gene as the pathway's
driver (⑥); and, when splitting cells by a score, match on lineage composition
as well as depth, and perform the control at the cell level (⑦).

More generally: a candidate list produced by this framework at current outcome
GWAS power should be reported as *the set of genes that passed screening under
these conditions*, not as a set of targets. Our two lists, produced from
identical exposure data under two outcome datasets, share nothing.
