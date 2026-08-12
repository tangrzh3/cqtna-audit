**⚠ 已被 `MANUSCRIPT_v2_dual_thread.md` §1 取代。**

# Introduction — draft v1 (2026-08-11)

> 结构说明（中文，不进正文）
> 四段：框架为何有吸引力 → 推断链的关节在哪 → 为何选黑色素瘤（自带阳性对照与自带混杂）
> → 本文做了什么、结论是什么。最后一段直接把"失败率是结论"摆出来，不留悬念。
> 不在 Introduction 里预告 TPI1 的生物学，只说明它是示例。

---

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
