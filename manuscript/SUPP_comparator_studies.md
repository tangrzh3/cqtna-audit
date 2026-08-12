# Supplementary — Item-by-item comparison with three studies using the same framework

回应外部审稿："'None of the three published studies…' 属于强主张，需要逐项证据表。"

## 1. 被比较的三项研究

| # | 研究 | 期刊/年 | 暴露 | 结局 | 本文档中的简称 |
|---|---|---|---|---|---|
| C1 | Zheng J, Yang Q, Liu H, et al. *Integrating single-cell transcriptome-wide Mendelian randomization and differentially expressed gene analyses to prioritize dynamic immune-related drug targets for cancers* | Adv Sci 2025 (e07451) | 动态 CD4⁺T eQTL | 多种癌症 | Adv Sci |
| C2 | Wu X, Ying H, Yang Q, et al. *Transcriptome-wide Mendelian randomization during CD4⁺ T cell activation reveals immune-related drug targets for cardiometabolic diseases* | Nat Commun 2024 (53621-7) | 同上 | 心血管代谢疾病 | Nat Commun |
| C3 | Cui K, Zou Q, Qu X, et al. *Transcriptome-wide Mendelian randomization and single-cell analysis during CD4⁺ T cell activation deciphers immunotherapeutic targets for colorectal cancer* | npj Precis Oncol 2025 (01236-6) | 同上 | 结直肠癌 | npj PO |

⚠ `pdf_text/s41586-022-04713-1.txt`（Nathan et al., *Nature* 2022，单细胞 eQTL 的
动态 T 细胞状态依赖）**不是**同框架的对照研究，而是一个 eQTL 资源论文，
且是我们在 Step 39 核查过的数据集之一。不列入本表。

## 2. 逐项对照

判定方式：对每份全文做关键词检索（不区分大小写），并人工核对命中上下文。
"—" = 全文无相关表述。计数为词频，仅作定位用，结论以上下文核对为准。

| 诊断项 | 本文 | Adv Sci | Nat Commun | npj PO |
|---|---|---|---|---|
| **① 显著信号是否按该性状的已知位点做归因** | 是（按独立位点计，4.09×，P=0.028）| — | — | — |
| ↳ 关键词 "known loci/locus" 命中 | — | **0** | **0** | **0** |
| ↳ 关键词 "GWAS Catalog" 命中 | — | 0 | 0 | 3（用于注释，非归因）|
| ↳ 关键词 "previously reported" 命中 | — | 0 | 1 | 0 |
| **② coloc 与 SMR/HEIDI 的关系** | coloc 为主判据，SMR/HEIDI 为辅，报告每次 HEIDI 的 SNP 数 | coloc（49）+ SMR（3），无 HEIDI | coloc（36），**无 SMR、无 HEIDI** | coloc（30）+ SMR（25）+ HEIDI（6），并列 |
| **③④ 是否检验结局功效对名单的影响** | 是（两轮 + 五 release + 模拟）| — | — | — |
| ↳ 关键词 "power" 命中 | — | 14（多为统计功效一般性表述）| 3 | 1 |
| **④b 是否报告 GWAS 峰与 eQTL 峰距离 / 多信号** | 是 | — | — | — |
| ↳ "credible set / fine-map" 命中 | — | **0** | **0** | **0** |
| **⑥ 是否声明通路中有工具变量的基因比例** | 是（1/22；2/28）| — | — | — |
| **⑦ 是否做细胞层面的谱系纯度对照** | 是 | — | — | — |
| ↳ "purity / contamination" 命中 | — | **0** | **0** | **0** |
| **单工具变量下 p 由结局唯一决定，是否点明** | 是（正文首节）| — | 1 处提及单 SNP，未展开推论 | — |
| 敏感性分析（多工具变量）| 是（IVW/IVW-MRE/WMed；weighted mode 因不适用而排除）| MR-Egger 1、WMed 2 | MR-Egger 2、WMed 3、异质性 8 | MR-Egger 1、**WMed 0** |
| Steiger | 是（并说明其近乎必然通过）| 7 | 7 | 4 |
| 复现/独立队列 | 尝试并报告不可得 | "replicat" 0 | 3 | 5 |

## 3. 可以写进正文的表述（已按证据收紧）

**可以写**：
> 在我们检索到的三项使用同一框架的研究中，均未见对显著信号按该性状已知位点做
> 归因的分析，未见细胞层面的谱系纯度对照，也未见 credible set 或精细定位层面的
> 讨论；三项中仅一项使用 HEIDI，且与 coloc 并列而非分层。

**不可以写**（原稿曾出现，已改）：
> ~~"None of the three published studies reports any of these seven checks."~~
> 七项中有若干项（如敏感性分析、Steiger）三项研究**都做了**，只是做法与我们不同。
> 强主张只在有逐项证据的四项上成立。

## 4. 本表的局限

1. 判定基于 PDF 抽取的正文文本，**未包含补充材料**。若某项检查写在补充材料中，
   本表会漏判。正文表述已相应限定为"未见于正文"。
2. PDF 文本抽取存在断词（如 "Tc e l la c t i v a t i o n"），故多词短语检索不可靠；
   本表只使用单词级检索 + 人工核对。
3. 三项研究的结局性状不同（癌症/心血管代谢/结直肠癌），**①位点归因这一项对
   非色素相关性状的具体形式会不同**，但"是否对照该性状的已知位点"这一问题本身适用。

⚠ 投稿前须补：三项研究的完整著录信息、以及是否需向作者核实补充材料内容。
