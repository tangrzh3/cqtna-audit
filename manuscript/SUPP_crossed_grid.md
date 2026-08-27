# S48：双轴交叉网格的完整叙述

**编号**：S48　**日期**：2026-08-24
**来源**：本节原为正文 Part 1 的 `Both axes crossed`，按第三方评审的压缩意见移出。
**内容逐字保留，未作删改。** 正文保留判定表与两句读法。

---

#### Both axes crossed

The main grid is the three diseases at their higher-powered
outcome crossed with the two exposure resources — six cells — with the
lower-powered HCC outcome held back as a power sensitivity, since the two HCC
levels are drawn from one resource and are not independent (Supplementary S20 §9).
Scoring every cell against its own disease's known loci, **one of the six is void
and the other five all enrich**: 4.96 and 6.31 for melanoma, 7.88 and 10.87 for
HCC-high, 3.91 for RA on CD4⁺ T cells, with RA on whole blood void on its
mismatched control as above. **Four of the five reach nominal one-sided Fisher
P < 0.05** — and remain so after Benjamini–Hochberg across the five, though that
is not the relevant caveat: the cells share exposure data and reference lists, so
they are consistency of direction rather than five independent replications. The
one that does not reach nominal significance is HCC-high on CD4⁺ T cells
(P = 0.123), on a numerator of one of two loci. Of the four that do, the RA cell
is reported as exploratory for the reasons given above, so **the confirmatory
count is three tumour cells**. Note that the two non-significant cells are no longer the same
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
melanoma on CD4⁺ T cells). Matching on **density** — background loci matched on record
count, span and gene count, requiring 100% matched coverage and reporting the
whole tolerance scan; the matching rules were fixed in Supplementary S38 before
the run, the tolerance reported was chosen alongside the results (S38 §2) —
gives **4.87-fold (empirical P = 0.0050)** for melanoma on
CD4⁺ T cells, 6.14-fold (P ≤ 1×10⁻⁴) on whole blood, 8.24-fold (P = 0.0045) for
HCC-high on whole blood, 12.85- and 12.09-fold (P = 0.0065 and 0.0012) for
HCC-low, and **4.52-fold (P ≤ 1×10⁻⁴) for RA on CD4⁺ T cells**. Values given as ≤ 1×10⁻⁴ are at the resolution floor of 10,000 permutations, not point estimates. Every one of those
holds across all tolerances that reach full matched coverage; the tightest
tolerance fails coverage on four cells and returns no P value rather than a
partially matched one (Supplementary S38 §1.7). These values follow a
correction: an external reviewer found the sampler's result moved with the order
of the input rows, an artefact with no statistical content, and the repaired
version is order-invariant by construction and by test (Supplementary S38-A1).
Every cell's verdict is unchanged by the repair. The values themselves moved by
up to 0.37-fold (hepatocellular carcinoma at its lower power on whole blood,
12.46 to 12.09), six of the eight downwards, with the melanoma CD4⁺ cell moving
0.08-fold. We report the sampler as randomized-greedy matching: it fills a
matched set greedily in a random order redrawn for each replicate, which is a
well-defined randomised null but is not uniform sampling from the set of all
feasible complete matchings.

Two readings follow. First, on bounded loci density matching barely
attenuates: melanoma on CD4⁺ T cells is 4.96-fold unmatched and 4.87-fold
density-matched. This is not evidence that the attenuation reported earlier was
spurious — the earlier 3.44-fold came from quantile-stratified matching, a
different test that we do not treat as the same quantity (Supplementary S38 §4) —
but under this matching rule — record count, span and gene count, at the
tolerances that achieve complete matching — the enrichment persists against the
density-matched null. That is a statement about one specified null, not a
demonstration that density plays no part. Second, **HCC-high on CD4⁺ T cells does not clear the
density-matched null** (P = 0.107, and P ≥ 0.05 at every usable tolerance), so by
S38's reading table its attribution claim is descriptive only; that is the same
cell and the same power limit that leaves it at P = 0.123 on the Fisher test, so
the two controls agree rather than conflict. The density-matched control
covers all eight cells including RA, but **it does not relieve the RA cell of its
dependence on the mismatched-list control**: the two ask different questions —
whether the enrichment follows from locus density, and whether it is specific to
this outcome's own genetics — and passing the first says nothing about the
second. RA's difficulty is entirely with the second.

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
gap between eQTL and GWAS colocalisation (Rosen et al. 2026). Our design asks the mirror
question — with the exposure held fixed, does a *higher-powered outcome* make
target nomination more reliable? — and the answer here is that it does not: the
meta-analysis produced more MR discoveries, lower colocalisation support, and a
candidate list with no genes in common with the previous one. That outcome
differs from its predecessor in study composition as well as power, so this is
not a power experiment either; only the nested-release series isolates power, and
there the significant list is stable rather than replaced. The two axes are not
interchangeable, and a study underpowered on one cannot be rescued by the other.
Reporting which axis a claim rests on should be routine.
