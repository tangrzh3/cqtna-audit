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
