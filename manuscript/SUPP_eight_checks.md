# S49：八项检查的完整说明

**编号**：S49　**日期**：2026-08-24
**来源**：本节原在正文 Discussion，按第三方评审的压缩意见移出。
**内容逐字保留，未作删改。** 正文保留八条的编号与标题句。

⚠ 第九项检查（`cqtna_decomposition`，把提名名单按结局是否已达全基因组显著分层）
见 S42；R 实现见 `cqtna` 包。

---

**(i) Annotate the signal against the outcome's own known loci, counting
bounded loci rather than gene records.** Every result in this paper's first
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
