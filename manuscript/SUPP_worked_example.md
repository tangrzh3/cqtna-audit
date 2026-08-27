# S47：一次提名，追到底——工作示例的完整记录

**编号**：S47　**日期**：2026-08-23
**来源**：本节原为正文 Part 3 的六个小节，按第三方评审意见移出，
使正文保持为一条审计线。**内容逐字保留，未作删改。**

⚠ 这不是本文的第二个发现。恒等式说提名名单是结局 GWAS 在更宽阈值下的重读，
但没有说"有人拿其中一条去坐实"会怎样——而那正是这个领域在做的事。
做一次，就给那把梯子的每一级标了价。**TPI1 不是靶点，本文也不如此主张。**

---

#### What can be interrogated, stated at the level it applies to

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

#### Time specificity is precision, not amplitude

For the one instrumented gene, TPI1, a usable instrument exists at 16 h and at no
other timepoint. The effect estimate is largest at rest and simply far noisier
there, so the window is a statement about precision, not about amplitude
(Fig. 6a,b). The published genotype × pseudotime interaction test (Nathan et al. 2022) settles what
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

#### What the instrumented gene marks

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

#### Where the effect is, when MR cannot say

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

#### The functional layer has its own failure modes

Splitting cells by a score also splits them by lineage purity. When cells were
divided by glycolytic score, apparent differences in helper and regulatory
phenotype collapsed after matching on lineage composition at the cell level, and
14 of 28 module comparisons were concordant — exactly chance. A cluster-level
control passed while the cell-level control failed, so whether a purity control
works can depend on the level at which it is applied. Two processing errors of
our own, detected by these controls, are described in Supplementary S12; they are
instances of the same failure mode.

The perturbation layer is bounded in the same way. TPI1 is a hit in 628 of the
**1,471 human CRISPR screens** in BioGRID ORCS (Oughtred et al. 2021) that measured it (**42.7%**), a
core-essential profile alongside GAPDH and PGAM1 (46.6% and 47.2%) and an order
of magnitude above lineage-defining genes in the same screens (IRF4 3.9%, FOXP3
1.5%, MC1R 1.1%). Since the knockout is lethal in almost any cell type, **a
knockout cannot isolate a CD4-specific role** — a boundary on what perturbation
could add here, not evidence for the nomination. PGAM1, the strongest enzyme in
the functional data and the one MR cannot see, has the same profile, so the
property belongs to the pathway rather than to the named gene.

#### In patients: present in each cohort, confirmed in none twice

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

---

## 逐层判定表（2026-08-26 由正文 Discussion 移入）

**来源**：正文 Discussion 的 TPI1 状态表与其后一段。
**内容逐字保留，未作删改。** 正文保留两句读法与指针。

---

The worked example should be read in that light, and is best stated as a ladder
rather than a verdict — the same layered reporting Howe et al. use when a cellular
model proves only partly transportable (Howe et al. 2022):

| Claim about TPI1 | Status here |
|---|---|
| Causal target for melanoma | **Not supported** — FDR = 0.119, prior-dependent colocalisation, region not fine-mappable |
| Dynamic *genetic* effect across activation | **Not supported** — no genotype × pseudotime interaction at any of three lead variants |
| Dynamic *expression* across activation | **Supported** — Moran's I = 0.664 |
| Membership of a definable CD4⁺ metabolic state | **Supported** — ranks 19th of 7,653 on an axis built without it, and 18th and 28th in a second dataset with protein-based lineage calling; axis enrichment 21.0–27.4-fold |
| Predicts checkpoint-blockade response | **Unstable** — significant on discovery, no arm confirmed twice |
| Functional consequence isolable by perturbation | **No** — a hit in 628 of 1,471 human CRISPR screens (42.7%), a core-essential profile that cannot isolate a CD4-specific role |
| Recognised as a target for this disease elsewhere | **No** — in Open Targets (Ochoa et al. 2023), TPI1's strongest disease associations are triosephosphate isomerase deficiency and neurodegenerative disease; melanoma is not among its leading associations |

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
the selection denominators at gene, pathway and mechanism level, is
Supplementary S12; the fourteen pre-registration documents are enumerated in
Methods.
