# S50：显著性属于哪个推断单位

**编号**：S50　**日期**：2026-08-24
**来源**：本节原为正文 Part 2 的一节，按第三方评审的压缩意见移出。
**内容逐字保留，未作删改。**

---

### Which unit the significance belongs to

A units problem has to be settled before any of the counts above is read. Those
3,556 records carry only **2,126 unique variants and 1,195
unique genes**, because one variant can be the lead eQTL for a gene in several
profiles and for more than one gene, so the testing family is records while the
list is reported by gene and the attribution by bounded locus. Records are
therefore not independent, and a Benjamini–Hochberg procedure over them does not
have a clean guarantee under that dependence.

We resolve this by separating the two things the record level is being asked to
do. **The record-level analysis reproduces the nomination pipeline being
audited** — it is what the studies under examination run, and what every
registered prediction, the release trajectory, both generalisations and the grid
were computed on, so it stays exactly as published. **The bounded locus is
the unit the audit's own conclusions are stated in.** We call these loci bounded
rather than independent deliberately: a non-overlapping partition guarantees that
two loci share no variant, not that they are statistically independent, and the
Fisher p-values below are nominal for that reason. That is not a convenience:
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
