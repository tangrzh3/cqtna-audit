# References

**核实状态**：2026-08-11 通过 CrossRef API 与 NCBI E-utilities 逐条核实。
下列条目的作者、期刊、年、卷、页、DOI **均来自查询返回**，非从记忆填写。

⚠ 核实过程中发现的坑，记下以免复查时重蹈：
- CrossRef 的 `query.bibliographic` 会把 **Faculty Opinions 的推荐短文**、**预印本**、
  **会议摘要**排在前面（Sade-Feldman、Jerby-Arnon、Tirosh、DICE、TCGA-CDR 均如此），
  须过滤 `type == "posted-content"` 与容器名含 "Faculty Opinions"/"Abstract"
- 少数标题检索会返回完全不相干的文献（Landi 2020 曾匹配到蜜蜂抗螨、FinnGen 曾匹配到核桃树病原菌）
- PubMed 检索 FinnGen 时**更正声明（Author Correction）排在原文之前**，须取后者
- 仍待补：**GSE282266 与 GSE199994 的关联论文**（项目中只用了 GEO 号），
  以及 GSE316760 / GSE300445 的著录

---

## A. 同框架对照研究

1. Zheng J, Yang Q, Liu H, et al. Integrating single-cell transcriptome-wide Mendelian
   randomization and differentially expressed gene analyses to prioritize dynamic
   immune-related drug targets for cancers. *Advanced Science* 2025;12:e07451.
   doi:10.1002/advs.202507451

2. Wu X, Ying H, Yang Q, et al. Transcriptome-wide Mendelian randomization during
   CD4⁺ T cell activation reveals immune-related drug targets for cardiometabolic
   diseases. *Nature Communications* 2024;15:9302. doi:10.1038/s41467-024-53621-7

3. Cui K, Zou Q, Qu X, et al. Transcriptome-wide Mendelian randomization and
   single-cell analysis during CD4⁺ T cell activation deciphers immunotherapeutic
   targets for colorectal cancer. *npj Precision Oncology* 2025;10:32.
   doi:10.1038/s41698-025-01236-6

> 逐项对照证据见 `SUPP_comparator_studies.md`；正文对这三篇的主张限定为"正文中未见"。

## B. 暴露与结局

4. **Soskic B, Cano-Gamez K, Smyth DJ, et al. Immune disease risk variants regulate
   gene expression dynamics during CD4⁺ T cell activation. *Nature Genetics*
   2022;54:817–826. doi:10.1038/s41588-022-01066-3** ← 本文暴露数据来源

5. Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a
   well-phenotyped isolated population. *Nature* 2023;613:508–518.
   doi:10.1038/s41586-022-05473-8
   （另见 Author Correction: *Nature* 2023;615:E19. doi:10.1038/s41586-023-05837-8）

6. Rashkin SR, Graff RE, Kachuri L, et al. Pan-cancer study detects genetic risk
   variants and shared genetic basis in two large cohorts. *Nature Communications*
   2020;11:4423. doi:10.1038/s41467-020-18246-6　（GWAS Catalog GCST90011809）

7. Landi MT, Bishop DT, MacGregor S, et al. Genome-wide association meta-analyses
   combining multiple risk phenotypes provide insights into the genetic architecture
   of cutaneous melanoma susceptibility. *Nature Genetics* 2020;52:494–504.
   doi:10.1038/s41588-020-0611-8　（GCST010302 / GCST010303 / GCST010304）

8. Nathan A, Asgari S, Ishigaki K, et al. Single-cell eQTL models reveal dynamic
   T cell state dependence of disease loci. *Nature* 2022;606:120–128.
   doi:10.1038/s41586-022-04713-1　（非对照研究；为 §3.1 核查过的资源之一）

## C. 功能与验证数据集

9. Sade-Feldman M, Yizhak K, Bjorgaard SL, et al. Defining T cell states associated
   with response to checkpoint immunotherapy in melanoma. *Cell* 2018;175:998–1013.e20.
   doi:10.1016/j.cell.2018.10.038　（GSE120575）

10. Jerby-Arnon L, Shah P, Cuoco MS, et al. A cancer cell program promotes T cell
    exclusion and resistance to checkpoint blockade. *Cell* 2018;175:984–997.e24.
    doi:10.1016/j.cell.2018.09.006　（GSE115978）

11. Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of
    metastatic melanoma by single-cell RNA-seq. *Science* 2016;352:189–196.
    doi:10.1126/science.aad0501　（GSE72056；本文独立复现队列）

12. Hugo W, Zaretsky JM, Sun L, et al. Genomic and transcriptomic features of response
    to anti-PD-1 therapy in metastatic melanoma. *Cell* 2016;165:35–44.
    doi:10.1016/j.cell.2016.02.065　（GSE78220）

13. Riaz N, Havel JJ, Makarov V, et al. Tumor and microenvironment evolution during
    immunotherapy with nivolumab. *Cell* 2017;171:934–949.e16.
    doi:10.1016/j.cell.2017.09.028　（GSE91061）

14. Thrane K, Eriksson H, Maaskola J, et al. Spatially resolved transcriptomics enables
    dissection of genetic heterogeneity in stage III cutaneous malignant melanoma.
    *Cancer Research* 2018;78:5970–5979. doi:10.1158/0008-5472.CAN-18-0747

15. Katko A, Potter SJ, et al. Gene regulatory network determinants of rapid recall in
    human memory CD4⁺ T cells. *Cell Reports* 2026;45(4):117103.
    doi:10.1016/j.celrep.2026.117103 — data: NCBI GEO **GSE282266**. [已核实]

16. Boukhaled GM, Gadalla R, et al. Pre-encoded responsiveness to type I interferon in
    the peripheral immune system defines outcome of PD1 blockade therapy.
    *Nature Immunology* 2022;23(8):1273–1283. doi:10.1038/s41590-022-01262-7
    — data: NCBI GEO **GSE199994**. [已核实]

17. Virós A, et al. Spatial transcriptomics of primary cutaneous melanoma.
    NCBI GEO **GSE316760** (submitted 2026-01-16); no associated publication indexed
    at the time of writing. [已核实为数据集，无关联论文]

17b. Pham F, Dufeu M, Benboubker V, et al. Spatial tumour-immune ecosystems shape the
    efficacy of anti-PD1 immunotherapy in primary cutaneous melanoma
    [Spatial Transcriptomics]. NCBI GEO **GSE300445** (submitted 2025-06-23);
    no associated publication indexed at the time of writing.
    [已核实为数据集，无关联论文]

17c. Guo X, et al. Contrasting cytotoxic and regulatory T cell responses underlying
    distinct clinical outcomes to anti-PD-1 plus lenvatinib therapy in hepatocellular
    carcinoma. *Cancer Cell* 2025;43(2):248–268.e9. doi:10.1016/j.ccell.2025.01.001
    — data: NCBI GEO **GSE235863**. [已核实]

17d. Ghouse J, Gellert-Kristensen H, O'Rourke CJ, et al. Genome-wide meta-analysis
    identifies nine loci associated with higher risk of hepatocellular carcinoma.
    *JHEP Reports* 2025;7(9):101485. doi:10.1016/j.jhepr.2025.101485
    — GWAS Catalog **GCST90809296**. [已核实]

17e. Võsa U, Claringbould A, Westra H-J, et al. Large-scale cis- and trans-eQTL analyses
    identify thousands of genetic loci and polygenic scores that regulate blood gene
    expression. *Nature Genetics* 2021;53:1300–1310. doi:10.1038/s41588-021-00913-z
    — eQTLGen Consortium, 2019-12-11 cis-eQTL release (n = 31,684). [已核实]

17f. FinnGen. Release R12 (2024), endpoints `C3_MELANOMA_SKIN_EXALLC` and
    `C3_HEPATOCELLU_CARC_EXALLC`. https://www.finngen.fi/en/access_results
    [数据集，按 FinnGen 引用规范著录]

## D. 参考面板、方法与工具

18. Byrska-Bishop M, Evani US, Zhao X, et al. High-coverage whole-genome sequencing of
    the expanded 1000 Genomes Project cohort including 602 trios. *Cell*
    2022;185:3426–3440.e19. doi:10.1016/j.cell.2022.08.004

19. Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation
    between pairs of genetic association studies using summary statistics.
    *PLoS Genetics* 2014;10:e1004383. doi:10.1371/journal.pgen.1004383

20. Wang G, Sarkar A, Carbonetto P, Stephens M. A simple new approach to variable
    selection in regression, with application to genetic fine mapping. *Journal of the
    Royal Statistical Society Series B* 2020;82:1273–1300. doi:10.1111/rssb.12388

21. Zhu Z, Zhang F, Hu H, et al. Integration of summary data from GWAS and eQTL studies
    predicts complex trait gene targets. *Nature Genetics* 2016;48:481–487.
    doi:10.1038/ng.3538

22. Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic
    causal inference across the human phenome. *eLife* 2018;7:e34408.
    doi:10.7554/eLife.34408

23. Chang CC, Chow CC, Tellier LC, et al. Second-generation PLINK: rising to the
    challenge of larger and richer datasets. *GigaScience* 2015;4:s13742-015-0047-8.
    doi:10.1186/s13742-015-0047-8

24. Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal
    and scalable single-cell analysis. *Nature Biotechnology* 2024;42:293–304.
    doi:10.1038/s41587-023-01767-y

25. Fornes O, Castro-Mondragon JA, Khan A, et al. JASPAR 2020: update of the
    open-access database of transcription factor binding profiles. *Nucleic Acids
    Research* 2020;48:D87–D92. doi:10.1093/nar/gkz1001

26. Schep AN, Wu B, Buenrostro JD, Greenleaf WJ. chromVAR: inferring
    transcription-factor-associated accessibility from single-cell epigenomic data.
    *Nature Methods* 2017;14:975–978. doi:10.1038/nmeth.4401

## E. 资源与数据门户

27. Schmiedel BJ, Singh D, Madrigal A, et al. Impact of genetic polymorphisms on human
    immune cell gene expression. *Cell* 2018;175:1701–1715.e16.
    doi:10.1016/j.cell.2018.10.022　（DICE）

28. Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human
    gene expression and splicing quantitative trait loci. *Nature Genetics*
    2021;53:1290–1299. doi:10.1038/s41588-021-00924-w　（eQTL Catalogue）

29. Sollis E, Mosaku A, Abid A, et al. The NHGRI-EBI GWAS Catalog: knowledgebase and
    deposition resource. *Nucleic Acids Research* 2023;51:D977–D985.
    doi:10.1093/nar/gkac1010

30. Ochoa D, Hercules A, Carmona M, et al. The next-generation Open Targets Platform:
    reimagined, redesigned, rebuilt. *Nucleic Acids Research* 2023;51:D1353–D1359.
    doi:10.1093/nar/gkac1046

31. Sun D, Wang J, Han Y, et al. TISCH: a comprehensive web resource enabling
    interactive single-cell transcriptome visualization of tumor microenvironment.
    *Nucleic Acids Research* 2021;49:D1420–D1430. doi:10.1093/nar/gkaa1020

32. Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics
    data via the Xena platform. *Nature Biotechnology* 2020;38:675–678.
    doi:10.1038/s41587-020-0546-8

33. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data
    resource to drive high-quality survival outcome analytics. *Cell* 2018;173:400–416.e11.
    doi:10.1016/j.cell.2018.02.052

34. Tsherniak A, Vazquez F, Montgomery PG, et al. Defining a cancer dependency map.
    *Cell* 2017;170(3):564–576.e16. doi:10.1016/j.cell.2017.06.010 — portal:
    DepMap, Broad Institute, https://depmap.org (release used: 22Q2 common-essential
    and gene-effect calls). [已核实；⚠ 具体 release 号须与分析脚本核对后定稿]

---

## 投稿前仍须处理

1. ~~条目 15–17、34 的著录~~ ✅ 2026-08-12 完成（含新增 17c–17f：GSE235863、
   GCST90809296、eQTLGen、FinnGen R12）。**唯一残留**：条目 34 的 DepMap release 号
   须与脚本核对
2. 统一到目标期刊的著录格式（本文件为作者-年-卷-页-DOI 的中性格式）
3. 若期刊要求，补全 3 位以上作者的完整列表

## E. 对照研究（2026-08-13 新增，全部经 CrossRef 核实）

35. Reales G, et al. Design and interpretation of eQTL–GWAS colocalisation studies:
    lessons from a large-scale evaluation. *PLOS Genetics* 2026.
    doi:10.1371/journal.pgen.1012141 [已核实]

36. Tambets R, et al. Extensive co-regulation of neighboring genes complicates the
    use of eQTLs in target gene prioritization. *Human Genetics and Genomics
    Advances* 2024;5:100348. doi:10.1016/j.xhgg.2024.100348 [已核实]

37. Rosen JD, et al. Higher eQTL power reveals signals that boost GWAS
    colocalization. *American Journal of Human Genetics* 2026.
    doi:10.1016/j.ajhg.2026.02.009 [已核实]

38. Howe LJ, et al. Evaluating transportability of in vitro cellular models to
    in vivo human phenotypes using gene perturbation data. *Nature Communications*
    2025. doi:10.1038/s41467-025-67199-1 [已核实]

39. Lin Z, Pan W. A robust cis-Mendelian randomization method with application to
    drug target discovery. *Nature Communications* 2024.
    doi:10.1038/s41467-024-50385-y [已核实]

40. Karhunen V, et al. Integrating genetic data with biological insight: a
    practical guide to cis-Mendelian randomization. *American Journal of Human
    Genetics* 2026. doi:10.1016/j.ajhg.2026.03.011 [已核实]

⚠ 正文中标 `[ref]` 的六处引用点对应 35–40，装配前须替换为期刊要求的编号格式。

---

## ⚠ 2026-08-14：投稿版重新编号（GB 已采用）

本文件的 17b–17f 是字母编号，投稿稿不能用。**GB 的 References 节已改为连续 1–45**，
映射如下（本文件保持原编号不动，作为溯源记录）：

| 本文件 | 投稿版 |
|---|---|
| 1–17 | 1–17（不变）|
| 17b, 17c, 17d, 17e, 17f | 18, 19, 20, 21, 22 |
| 18–40 | **+5**（如 19 coloc → 24；21 SMR → 26；35–40 → 40–45）|

GB 正文的六处旧引用已同步改为 40–45。**NC 尚未加参考文献列表**，加时须用同一套 1–45 编号。
⚠ 仍未解决：条目 39（DepMap）的 release 号须与脚本核对；条目 22 已把 FinnGen R13 的
两个 endpoint（`C3_MELANOMA_SKIN_WIDE`、`M13_RHEUMA`）并入同一条，若期刊要求按 release
分列须拆开。

### ✅ 2026-08-14 已补：BioGRID ORCS（投稿版编号 **40**）

**Oughtred R, Rust J, Chang C, et al. The BioGRID database: a comprehensive
biomedical resource of curated protein, genetic, and chemical interactions.
*Protein Sci* 2021;30:187–200. doi:10.1002/pro.3978**

核实方式：PubMed 全库检索 `"BioGRID ORCS"` **只返回这一条**（PMID 33070389）；
著录字段（15 位作者、卷 30、期 1、页 187–200、2021）取自 **CrossRef API 返回**，
非从记忆填写。**未找到 ORCS 的独立专文**，故引这篇 BioGRID 主文。
✅ **2026-08-17 已核对**：release = **BioGRID ORCS 2.0.18**。
下载归档 `orcs/human_screens.tar.gz` 里 1,953 个 human screen 文件的文件名
（`BIOGRID-ORCS-SCREEN_<id>-2.0.18.screen.tab.txt`）**全部**带同一个版本串，
无混版。TPI1 出现在其中 **1,471** 个 screen（`orcs/tpi1_all_screens.tsv` 恰 1,471 行，
无表头），命中 628 / 未命中 843 = **42.7%**，与正文一致。
著录时应写明 release：**BioGRID ORCS release 2.0.18**。
（各基因的 screen 数不同——FOXP3 1,050、MC1R 1,404、GAPDH 1,437 等——
1,471 是 **TPI1 自己的**分母，不是全库 screen 总数 1,953。）

**编号后果**：新条目插为 40，原 40–45（Reales…Karhunen）顺延为 **41–46**，
GB 与 NC 的正文引用已同步。下表的"投稿版"映射相应更新：18–40 → **+6**。

---

### ~~⚠ 缺一条：BioGRID ORCS~~（已解决，保留原文备查）

正文（GB 与全文源）报"TPI1 在 BioGRID ORCS 的 1,471 项人类 CRISPR 筛选中命中 628 项
（42.7%）"，**但本文件没有 BioGRID ORCS 的著录条目**，GB 该处因此无引用。
候选原文为 Oughtred R 等的 BioGRID 论文（*Protein Sci* 2021）或 ORCS 专文
（*Nucleic Acids Res* 2019/2021），**但尚未经 CrossRef 核实，故不写入编号列表**——
本项目的纪律是不从记忆填著录。**投稿前必须核实并插入**；插入位置建议为 39 之后
（新 40），其后 40–45 顺延为 41–46，GB 与 NC 的六处引用须同步。
⚠ 同时须核实所用的 ORCS release 号，与 1,471 这个分母对应。

### ★ 新增条目 47：Okada 2014（RA 已知位点参照系）—— 2026-08-17

**这条是漏的，不是可选的。** GB 与 NC 都写 "Okada et al. 2014 lead SNPs"
（GB 正文 1 处、GB Methods 1 处并注明 GCST002318、NC 正文 1 处），
而 1–46 号编号表里**根本没有这篇**。RA 是本文"非癌结局"这一格的全部依据，
其已知位点名单的来源必须可引。

核实方式（与 ORCS 同一套纪律，**不从记忆填**）：
1. GWAS Catalog REST `studies/GCST002318` → trait = Rheumatoid arthritis，
   `pubmedId` = **24390342**，author = Okada Y，journal = Nature；
2. NCBI E-utilities esummary（PMID 24390342）→
   *Nature* **2014 Feb 20**（epub 2013 Dec 25），**卷 506，期 7488，页 376–81**，
   **97 位作者**，DOI **10.1038/nature12873**。

> **47. Okada Y, Wu D, Trynka G, et al. Genetics of rheumatoid arthritis
> contributes to biology and drug discovery. *Nature* 2014;506(7488):376–381.
> doi:10.1038/nature12873 (GWAS Catalog GCST002318)**

**编号处置：追加为 47，不插队。** 1–46 的编号在上一窗口刚因插入 ORCS 而整体动过一次
（18–40 → +6），再插一次风险高于收益。若期刊要求按出现顺序排号，
**留到排版阶段一次性重排**，不要在投稿稿里再动。
GB 的 2 处与 NC 的 1 处 Okada 提及均已加上 `[47]`。

---

### 未被 GB 正文引用的条目 —— 2026-08-17 重新逐条核过

**原记录（7 条）不准确，实为 8 条**，且其中一条属于漏引而非多余。
核法：从 GB 正文（`## References` 之前的全部内容）提取所有 `[n]`、`[n,m]`、`[n–m]`
形式的引用（**范围用 en dash，早先的核对漏了这种写法**），与编号表 1–46 求差。

| 条目 | 内容 | GB 正文是否提到该主题 | 处置 |
|---|---|---|---|
| **33** | Kerimov 等，**eQTL Catalogue** | **提到**（"Across the whole eQTL Catalogue…"） | ✅ **已补引 `[33]`**，属漏引，非多余 |
| 17 | Virós 等，空间转录组 GSE31676 | 否 | 待定 |
| 18 | Pham 等，空间免疫生态 GSE300445 | 否 | 待定 |
| 32 | Schmiedel 等，**DICE** | 否（GB 只说"第三个资源"，未点名） | 待定 |
| 36 | Sun 等，**TISCH** | 否 | 待定 |
| 37 | Goldman 等，**UCSC Xena** | 否 | 待定 |
| 38 | Liu 等，**TCGA-CDR** | 否 | 待定 |
| 39 | Tsherniak 等，**DepMap** | 否（"DepMap" 一词在 GB 正文中 0 次） | 待定，且见下 |

余下 7 条对应的分析都只在全文源／补充材料里，GB 压缩版没有承载它们的句子。
**三条路，须作者定，不宜代拍**：
(a) 从 GB 的 References 中删去这 7 条——**会引起 17→46 段的重新编号**，
    而编号刚在上一窗口因插入 ORCS 而动过一次，再动风险高；
(b) 保留，并在 GB 的 Methods 已有的"…are given in Supplementary S12"这类指针处
    补上引用（GB 现在就是这么处理 `[9–13,19]` 的）——**但必须先逐条核实
    该数据集确实用于对应的补充分析**，不能为了安置引用而挂靠；
(c) 保留，另立 "Supplementary references" 一节。
⚠ GB 正文已贴着 8,000 词上限，(b) 若要加句子须同时减内容。

### ⚠ 条目 39（DepMap）的 release 号：**核不了，因为没有脚本**

HANDOFF v5/v6 都记着"DepMap release 号须与脚本核对"。2026-08-17 查证结果：
**仓库里既没有使用 DepMap 的脚本，也没有任何 DepMap 数据文件**
（`*.py` 全文检索 `DepMap|depmap|CRISPRGeneEffect|Achilles` 无命中；
全盘 `find` 同名文件无命中；`D:/Downloads` 亦无）。

即正文的 "TPI1 is DepMap-essential"（全文源 3 处）与著录里写的
"release used: 22Q2 common-essential and gene-effect calls"，
**来自一次未留痕的门户手查**，无法与脚本核对——**这项待办按原样是做不成的**。

可行的替代（与 ORCS/CrossRef 用的是同一套纪律）：
1. 上 DepMap 门户重查，把 **release 号、查询日期、TPI1 的 gene-effect 数值与
   common-essential 判定**记成一张小表落盘，著录按该 release 定稿；或
2. 把措辞改成不依赖具体 release 的说法；或
3. 删去该主张（它在 GB 正文里本来就不存在，只在全文源里）。

⚠ 顺带：**GB 从头到尾没出现过 "DepMap"**，所以对 GB 而言条目 39 属上表的
"未被引用"一类；此事只影响全文源与 assembled。

---

## NC 稿的引用密度 —— 2026-08-17 处理结果

**处理前：全文只有 2 处引用（`[40]` ORCS、`[42]` Tambets）。处理后 16 条被引。**

新加的 14 条都挂在 **NC 本来就在陈述的论断**上，没有为了安置引用而新造句子
（只多了 35 词，全是方括号标记）：

| 位置 | 加了 | 依据 |
|---|---|---|
| Introduction「a common route to causal target nomination」 | `[1–3]` | **本文审计的正是这三篇**，此前整个 Introduction 一条引用都没有 |
| Results「Eight CD4⁺ T cell activation profiles」 | `[4]` | Soskic 时程 |
| Abstract「a 12,530-case meta-analysis」 | `[6,7]` | 结局 GWAS |
| 「FinnGen's sequential releases」 | `[5,22]` | FinnGen 与 R12 endpoint |
| 「eQTLGen whole-blood cis-eQTLs」 | `[21]` | Võsa |
| 「Okada et al. 2014 lead SNPs」 | `[47]` | 见上节新增条目 |
| 「Colocalisation assigned…」 | `[24]` | coloc |
| 「…on 8–20 variants per test」 | `[26]` | SMR/HEIDI |
| 「The published genotype × pseudotime interaction test」 | `[8]` | Nathan |
| 「Raising eQTL sample size is known to…colocalisation gap」 | `[43]` | Rosen |

### ⚠ 剩下 31 条未被 NC 引用，原因是**结构性的，不是漏标**

未引清单：9–20、23、25、27–39、41、44–46。这些全是**数据集与软件类**引用
（ICB 队列、空转、multiome、1000G、susieR、TwoSampleMR、PLINK、Seurat、
JASPAR、chromVAR、GWAS Catalog、Open Targets、TISCH、Xena、TCGA-CDR、DepMap…）。

**GB 把它们全部放在 Methods 里；而 NC 稿目前没有 Methods 章节**
（章节只有 Abstract · Introduction · Results · Discussion · Display items · References）。

→ **正确的修法不是继续往正文塞引用，而是补上 NC 的 Methods。**
   Nature Communications 本来就要求 Methods，投稿时无论如何得写；
   一旦写了，这 31 条自然各归其位。在此之前**不应**把数据集引用硬塞进 Results，
   那会既伤可读性又不合格式。

⚠ 另：GB 的 Discussion 引了 `[41]`（Reales）、`[44]`（Howe）、`[45,46]`（Lin、Karhunen）
来做文献对位，NC 的 Discussion 是纯自证结构、没有对应句子。
**要不要在 NC 补这几句对位讨论，是取舍问题，留给作者定**——
NC 正文有字数上限，加了就得减别的。

- Randolph HE, Fiege JK, Thielen BK, Mickelson CK, Shiratori M, Barroso-Batista J,
  Langlois RA, Barreiro LB. Genetic ancestry effects on the response to viral
  infection are pervasive but cell type specific. *Science* 2021;374:1127–1133.
  doi:10.1126/science.abg0928 · PMID 34822289 · PMC8957271
  （核实：2026-08-26，NCBI E-utilities esummary；用于正文“three further exposure
  resources”中的 eQTL Catalogue QTD000588）
