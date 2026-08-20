> ## ⚠ 来源与地位（项目方 2026-08-20 加注）
>
> **本文件由第三方复核者起草**，基于 `review_packet_v9.zip`
> （SHA-256 `1B27C900…2F3968`）。**不是本项目自己写的预注册，项目也尚未采纳它。**
>
> 加这段抬头是因为：本项目的全部纪律建立在"谁在什么时间写下什么"之上，
> 而一份来自外部、写得很完整的方案，如果不标出处，
> 三个月后会被当成我们自己冻结的设计——**那正是本项目反复栽过的那类错误**。
>
> **采纳它需要三件事，都还没发生**：
>
> 1. **一个合格的外部队列**。方案 §2 要求与 FinnGen R12/R13、Rashkin
>    GCST90011809 及 Landi 参照名单的贡献队列**完全无participant 重叠**，
>    并需数据保管方签署的队列谱系表。**目前没有这样的数据集在手**，
>    所以本方案**没有任何已实现的统计功效**。
> 2. **作者决定是否接受它划定的范围**——尤其 §1 把确证性问题**收窄到
>    皮肤黑色素瘤的 outcome 侧一条**，明确排除 RA、CD4 特异性、
>    患者疗效分层与因果基因归属。这与 v9.3 的降级判断一致，但它是一个科学承诺。
> 3. **公开注册**（OSF / Registered Report）并拿到 DOI，之后才谈得上"事前"。
>
> **★ 更新（2026-08-20）：本项目已采纳本方案全文，记录见 `manuscript/PREREG_external_validation.md`（S40）。**
> 本文件**保持第三方起草时的样子，未经改写**——采纳记录另起一份，
> 因为"谁写的"和"谁采纳的"是两件必须能分开看的事。
>
> ⚠ 采纳**不等于**预注册。上列三项前提（合格外部队列、作者接受范围划定、公开注册取得 DOI）
> 中，第二项已完成，**第一、三项仍未完成**，因此本方案目前没有任何已实现的统计功效，
> 对外也不得称 pre-registered。详见 S40 §3。

---

# External validation protocol for review packet v9.3

**Design lock basis:** `review_packet_v9.zip`, SHA-256 `1B27C9001BC789CC1B0839076D5E1227BD567C1FB786889D8235C50B5E2F3968` (2026-08-20).

## One-sentence argument

In an outcome-independent cohort, we will test whether a candidate list generated from the frozen CD4⁺ T-cell cis-eQTL instruments again concentrates on melanoma loci known before outcome access, while the novel-locus component shows weaker transport, using the same bounded-locus convention and prospectively fixed validity controls.

## Terminology ledger

| Canonical term | Definition used here | Terms not to substitute |
|---|---|---|
| independent external validation | No participant or cohort contributes to FinnGen R12/R13, Rashkin GCST90011809, or the GWAS used to define the primary known-locus list | nested release, transfer test |
| bounded locus | `fixed_centre`, `locus_kb = 1000` as implemented in CQTNA; non-recursive assignment with a maximum 1-Mb span | independent locus, LD block |
| C1 attribution | A significant bounded locus is known if any qualifying record lies within the frozen known-locus window | causal-locus attribution |
| primary known-locus list | The packet's frozen Landi-derived list, used only if the validation cohort contributed no participants to it | updated Catalog list |
| updated reference sensitivity | A time-stamped list compiled before outcome access, excluding the validation cohort and its publication | primary reference |
| void | The own-list result is not interpreted because a prospectively specified negative-control list also enriches | negative result |
| inconclusive | The cohort fails a prespecified information, coverage, or valid-draw gate | null result |
| randomized-greedy matching | Canonically sorted, per-replicate randomized greedy matching without replacement | uniform sampling of all feasible matchings |

---

## Draft preregistration (English)

### 1. Confirmatory question and scope

The confirmatory question is restricted to the outcome-side claim in cutaneous melanoma. It does not test whether stimulated CD4⁺ T-cell eQTLs generalise to other cell types, ancestries, autoimmune diseases, patient response, or causal gene assignment. Rheumatoid arthritis is excluded from confirmatory validation because the v9.3 analysis classifies that evidence as exploratory and window-sensitive.

### 2. Validation dataset eligibility

Before receiving any association statistic, the external data custodian will supply phenotype, ancestry, case/control, imputation and cohort-lineage metadata. A dataset is eligible only if all conditions below are met.

1. Participants are absent from FinnGen R12/R13, Rashkin GCST90011809 and every cohort contributing to the Landi-derived primary reference list. A cohort-level overlap table and a signed confirmation from the data custodian are required.
2. The primary phenotype is invasive cutaneous melanoma or a definition demonstrably compatible with it. In-situ-only, acral-only, uveal and mucosal melanoma analyses are excluded from the primary endpoint and may be reported separately.
3. The primary analysis is restricted to a genetically defined European-ancestry stratum because the frozen exposure eQTLs and main LD framework are European-ancestry resources. Other ancestry strata are prespecified stress tests, not substitutes for the primary analysis.
4. At least 95% of the frozen exposure records have an unambiguous outcome allele, effect estimate and standard error after deterministic harmonisation. The absolute coverage difference between primary-known and primary-novel background loci must be no more than 5 percentage points.
5. The dataset passes the information gate in Section 7 before outcome statistics are released to the analysis team.

Preferred acquisition is either (i) a de novo GWAS in an independent biobank or (ii) a consortium-generated leave-cohort-out meta-analysis excluding all overlapping discovery/reference cohorts. A new melanoma meta-analysis may not be used as a single external dataset merely because its publication is newer; cohort composition, rather than publication date, determines independence.

### 3. Roles, blinding and lock

The authors will deposit the protocol, executable container, `renv.lock`, CQTNA source tarball, frozen instrument table, known-locus lists, comparator lists, test-family manifest and SHA-256 manifest before external association results are transferred.

The external analyst will receive only the variant request and harmonisation specification, without known/novel labels, expected directions or candidate ranks. The external analyst will return allele, beta, standard error, allele frequency, effective sample size and QC fields. Regional summary statistics will be returned only for the prespecified colocalisation sensitivity. The analysis repository will be unblinded only after the registration DOI and commit hash are public.

### 4. Frozen primary analysis

The primary exposure is the Soskic CD4⁺ T-cell activation time-course resource used in v9.3. No instrument is reselected using the validation outcome.

- Testing family: all frozen records that pass deterministic allele harmonisation.
- Record statistic: the same single-variant Wald test used in v9.3.
- Multiple testing: Benjamini–Hochberg FDR across the frozen harmonised record family; threshold 0.05.
- Locus partition: `fixed_centre`, `locus_kb = 1000` (maximum span 1 Mb), unchanged from v9.3.
- Attribution: C1 (`any_record`) against the frozen Landi-derived melanoma list with a 1,000-kb known-locus window.
- Background: every testable bounded locus in the same harmonised family.
- Primary statistic: one-sided Fisher exact enrichment of primary-known loci among FDR-significant bounded loci relative to all testable bounded loci, reported as fold enrichment, odds ratio, exact P value and 95% confidence interval.

The primary analysis is not interpreted when fewer than eight FDR-significant bounded loci are observed; this is an information failure, not evidence against enrichment.

### 5. Prospectively fixed validity gates

All gates are evaluated before declaring primary success.

1. **Specificity gate.** The same significant and background loci are scored against four frozen comparator lists (HCC, lung, colorectal and prostate). The melanoma fold must exceed the largest comparator fold. Comparator P values are adjusted by Holm's method; any adjusted one-sided P < 0.05 voids the primary cell.
2. **Density gate.** The v9.3 randomized-greedy density-matched null is run with 10,000 draws, seed 1 and tolerance 1.00. At least 95% of draws must complete. The full tolerance grid (0.10, 0.25, 0.50, 1.00 and 2.00) is reported. A coverage-failing row is NA, never evidence.
3. **Provenance gate.** No locus added from the external outcome's own publication may enter either the primary or updated known-locus list.
4. **Harmonisation gate.** Palindromic or allele-inconsistent variants are removed by the frozen rule; neither effect direction nor P value may be used to resolve an ambiguity.

### 6. Decision rule

The external validation is declared successful only if:

1. all dataset, coverage, information and provenance gates pass;
2. at least eight FDR-significant bounded loci are available;
3. the primary fold enrichment is greater than 1 with one-sided Fisher P < 0.05;
4. the melanoma fold exceeds every prospectively frozen comparator and no comparator survives Holm correction; and
5. the density-matched primary row is valid and directionally agrees (fold > 1).

If the primary Fisher test fails after all information gates pass, the external validation is negative. If a comparator enriches, the result is void. If an information or coverage gate fails, the result is inconclusive. These categories may not be collapsed.

### 7. Prospective information and power gate

Power is evaluated before outcome transfer using metadata and the frozen discovery estimates only. Effective sample size is

`N_eff = 4 / (1/N_case + 1/N_control)`.

For each candidate validation cohort, external standard errors are projected as

`se_ext = se_disc × sqrt(N_eff,disc / N_eff,ext)`.

Ten thousand complete-pipeline simulations are run under effect-retention factors of 0.50, 0.75 and 1.00, including harmonisation coverage, BH testing, locus formation and all primary gates. A cohort qualifies for a confirmatory run only if, under the conservative 0.50 scenario, at least 80% of simulations produce eight or more significant bounded loci and pass the primary Fisher endpoint. Otherwise the cohort may be analysed only as a labelled pilot. The threshold is not revised after seeing outcome statistics.

### 8. Key secondary analyses

Secondary analyses do not rescue a failed or void primary endpoint.

1. **Discovery-list transport.** The eight v9.3 melanoma bounded loci and their representative records are frozen before transfer. Report sign concordance, Bonferroni replication and effect-size calibration, separately for the four primary-known and four primary-novel loci. Because each stratum contains only four loci, this is descriptive unless a future protocol supplies a justified larger frozen set.
2. **Continuous all-locus evidence.** Compare external locus-level evidence between primary-known and primary-novel background loci using a cluster-aware model or matched permutation prespecified for record count, locus span, unique-gene count, allele frequency and eQTL strength. This preserves information when few loci cross FDR.
3. **Candidate-list stability.** Report locus- and gene-level Jaccard index, sensitivity and precision between v9.3 and external FDR lists, split by known/novel status and compared with the prospectively simulated interval at the external effective sample size.
4. **Second exposure resource.** Repeat the locked analysis with eQTLGen. This tests resource robustness but is not an independent replication because it shares the external outcome and reference lists with the primary cell.
5. **Updated knowledge sensitivity.** Repeat attribution using a date-stamped melanoma list compiled before unblinding and excluding the external cohort. The original Landi-derived list remains primary to preserve exact reproducibility.
6. **Regional colocalisation.** At frozen loci only, repeat single- and multi-signal colocalisation using validation-ancestry LD. Report PP.H3, PP.H4 and prior/window sensitivities. This does not convert locus replication into causal-gene validation.
7. **Ancestry stress test.** Apply the frozen pipeline to other ancestry strata only when instrument coverage and ancestry-matched LD meet separate preregistered gates. Report it as transportability, not as a replacement for the European-ancestry primary analysis.

### 9. Optional high-impact extension

A second, independently powered HCC validation may be added as a separate preregistration. It must exclude FinnGen and all cohorts in GCST90809296, use an HCC known-locus list that excludes both the current and external outcome publications, and pass its own eight-locus information gate. It is not pooled with melanoma. Rheumatoid arthritis should remain exploratory unless a new protocol prospectively fixes non-oncology and non-immunology comparator lists and reproduces the result across more than one window.

An orthogonal causal-gene module may be added only with an adjudication panel blinded to the MR nominations and a benchmark restricted to evidence not derived from molecular-QTL colocalisation. Fewer than 20 evaluable loci should be reported descriptively rather than as an error-rate estimate.

### 10. Outputs

The validation release will contain:

- cohort-lineage and overlap audit;
- preregistration DOI, Git commit and input-hash manifest;
- harmonisation flow diagram and exclusions by known/novel status;
- primary 2×2 table plus comparator and density-control tables;
- discovery-to-external effect plot and known/novel transport plot;
- complete null, void and inconclusive results;
- executable container and machine-readable result JSON/TSV.

---

## 中文设计说明

### 为什么以独立 melanoma outcome 为主，而不是再做一个疾病网格

v9.3 已把可确认结论限制在 melanoma/HCC，并把 RA 降为探索性。当前最缺的不是更多相关格子，而是一个与 FinnGen、Rashkin 及已知位点来源均无样本重叠的 outcome。先在同一疾病中复现，能够把“疾病差异、祖源差异、表型差异”与真正的外部重复分开。

### 这个方案验证什么、不能验证什么

它直接验证两条：冻结暴露工具后，显著名单是否仍富集于 outcome 已知位点；已知位点部分是否比 novel-locus 部分更可运输。它不验证 CD4 特异性、患者疗效分层或某个基因就是因果靶点。后几项需要不同数据和不同实验，不能捆成一个“外部验证成功”。

### 为什么不能直接用 FinnGen DF13 或最新大型 meta

DF13 与 R12 共享参与者，是 transfer/temporal update，不是独立外部队列。新的大型 meta 只有在提供排除 FinnGen、Rashkin 和 Landi 贡献队列后的 new-only/leave-cohort-out 结果时才合格；“论文更新”不等于“样本独立”。

### 最小版与冲击 10 分期刊版

- **最小可信版**：一个真正独立的 melanoma outcome，盲法运行以上主终点和四名单控制。
- **更强版**：在最小版基础上增加独立 HCC outcome，并加入盲法因果基因基准；两项分别判读，不合并 P 值。
- **不够的版本**：只跑 DF13、只换 500/1,000-kb 窗口、只增加另一个 eQTL 数据库、或用包含原队列的新 meta。这些属于敏感性/transfer，不是外部验证。

## Assumptions or missing inputs

1. No eligible external cohort, exact phenotype definition, case/control count or cohort-overlap certificate has yet been supplied.
2. The 50% effect-retention power rule is deliberately conservative and must be simulated on the frozen record table before registration; it is not a claimed achieved power.
3. Exact comparator-list filenames and checksums should be inserted into the final hash manifest.
4. If the external cohort contributed to the Landi study, the current primary reference cannot be used without circularity; a different fully independent cohort is required.

## Claim-evidence map

| Claim | Evidence planned | Status before execution |
|---|---|---|
| Outcome-known loci dominate the external nomination list | independent melanoma GWAS; frozen C1 Fisher endpoint and controls | needs new data |
| Novel-locus candidates transport less well | eight-locus transport plus continuous all-locus analysis | needs new data; binary comparison underpowered |
| Result is not an artefact of one exposure resource | locked eQTLGen secondary analysis | needs new data; not an independent replication |
| Finding generalises beyond melanoma | independent HCC extension | optional; not supported by the melanoma validation alone |
| Nominated genes are causal genes | blinded non-QTL benchmark/perturbation evidence | outside the primary validation |

## Why this structure

- It tests the strongest v9.3 claim with one independent variable changed: the outcome cohort.
- It preserves the frozen 1,000-kb `fixed_centre`/C1 convention and turns the previously post-hoc multi-list comparison into a prospective validity gate.
- It distinguishes a negative finding from a void control and an underpowered/inconclusive run.
- It prevents a newer but overlapping meta-analysis from being presented as external replication.
