# Supplementary — Systematic search for an independent replication cohort

回应审稿：患者侧结果需要在独立的、带应答标注的、CD4 可解析的队列中复现。
本文件给出**完整检索过程与逐个数据集的排除理由**，使 Limitations 中
"we identified one suitable cohort" 这句话可被核查，而不是只能相信。

检索日期：2026-08-11。工具：NCBI E-utilities（db=gds）、EBI BioStudies API、
KU Leuven RDR Dataverse API。

---

## 1. 纳入标准（检索前固定）

| # | 要求 | 理由 |
|---|---|---|
| 1 | 接受免疫检查点抑制剂治疗的患者样本 | 与发现队列同一临床情境 |
| 2 | 单细胞分辨率，**CD4⁺T 可解析** | 组织水平测的是肿瘤，见正文 §2.7 |
| 3 | **有 responder / non-responder 标注** | 主检验所需 |
| 4 | **肿瘤浸润**（非仅外周血） | 发现队列的结果在肿瘤浸润 CD4 中 |
| 5 | 治疗前样本 | Pre 是本文最弱的一支，最需要复现 |
| 6 | 应答与治疗方案不混杂 | 否则组间差异不可归因 |

## 2. 检索式

```
# GEO (db=gds)，五组
A "Homo sapiens"[Organism] AND melanoma AND (single cell OR scRNA)
    AND (anti-PD-1 OR pembrolizumab OR nivolumab OR immunotherapy OR checkpoint)
    AND "expression profiling by high throughput sequencing"[DataSet Type]
B melanoma AND (responder OR "response to therapy") AND "single cell"
C ("checkpoint blockade" OR "immune checkpoint") AND single-cell AND ("T cell" OR CD4)
    AND (responder OR "non-responder")
D melanoma AND (tumor-infiltrating OR TIL) AND scRNA AND (anti-PD-1 OR ...)
    AND 20[n_samples]:500[n_samples]
E melanoma AND (neoadjuvant OR "pathologic response") AND single-cell     → 0 条
```
命中后逐个拉取 `!Sample_characteristics_ch1`，检查是否存在应答字段。

## 3. 逐个数据集的判定

| 数据集 | 癌种 | 组织 | 应答标注 | 判定 |
|---|---|---|---|---|
| **Pozniak et al. 2024**（RDR doi:10.48804/GSAXBN）| 黑色素瘤 | 转移灶 | **BT/OT × R/NR，在 RDS meta_data 中** | ✅ **纳入** |
| E-MTAB-13770 | 黑色素瘤 | 肿瘤浸润 + PBMC | 有（R/NR_adj/NR_nadj/TF）| ❌ **应答与方案混杂**：6 个应答者全部为 Ipi+Nivo 双药，非应答者多为单药或含 BRAF/MEK 抑制剂；治疗前 CD45+ 仅 3 个应答者；部位跨脑/淋巴结/皮肤/鼻腔/肺 |
| GSE294273 | 黑色素瘤 | **淋巴结** | 4 R / 4 resistant / 4 untreated | ❌ 组织不符 + 4 对 4（Wilcoxon 最小可能 P = 0.029）|
| GSE295942 | 黑色素瘤 | 血 | GEO metadata 中无 | ❌ 组织不符 + 无标注 |
| GSE174401 | 黑色素瘤 | 脑/软脑膜转移 | 无（仅 site）| ❌ 无标注 |
| GSE120575 | 黑色素瘤 | 转移灶 | 有 | ——**本文发现队列** |
| GSE115978 / GSE72056 | 黑色素瘤 | 转移灶 | **无** | ❌ 无应答标注（本文用于区室归属，不用于应答分层）|
| GSE235863 | 肝细胞癌 | 肝肿瘤 + 血 | 有（CR/PR vs SD/PD）| △ 跨癌种**泛化**，非复现；联合治疗（anti-PD-1+lenvatinib）；CD45+ 文件覆盖九个患者。**2026-08-12 已下载并逐项核实，见下注①** |
| GSE243572 | 肝细胞癌 | **仅 PBMC** | 有（长期应答 vs 早期进展）| △ 同上，且仅外周血（组织不符，标准 4 不满足）|
| GSE270235 | 三阴乳腺癌 | 仅 PBMC | 无 | ❌ |
| GSE301720 | 头颈鳞癌 | — | — | ❌ 样本量 8 |
| GSE324815 / GSE325923 | 黑色素瘤 | 肿瘤 | 联合 CD40 激动剂 | ❌ 干预不可比 |

## 4. 结论

**符合全部六项标准的公开数据集有且仅有一个**（Pozniak et al., *Cell* 2024）。
其预注册分析见 `PREREG_pozniak_replication.md`。

两个跨癌种数据集（GSE235863、GSE243572）保留作为**泛化**检验，
其性质与复现不同，不得混为一谈。

### 注① GSE235863 的核实结果（2026-08-12，Step 86）

下载并检查后，本表的判定须补充如下事实（对本文有两个方向相反的后果）：

| 项 | 值 | 对应标准 |
|---|---|---|
| 规模 | 191,435 细胞 × 28,671 基因，30 样本，9 患者 | — |
| 应答 | **4 应答（P5/P11/P18/P27）vs 5 无应答（P1/P15/P26/P51/P52）**；标注在 GEO 的 GSM 标题/characteristics，**不在 h5ad 的 `obs` 内** | 标准 3 ✅ |
| 组织 | 肝肿瘤 86,252 细胞 / 外周血 105,183 细胞，配对 pre/post | 标准 4、5 ✅ |
| CD4 | **59,528 细胞**，9 个亚簇（含 FOXP3、CXCL13、CTLA4）| 标准 2 ✅ |
| 方案 | 全部患者同一方案（抗 PD-1 + 仑伐替尼），应答与方案不混杂 | 标准 6 ✅ |
| 细胞类型 | **CD45⁺ 分选，7 个免疫簇，无恶性细胞** | — |

**后果一（否定）**：因无恶性细胞，该数据集**不能**用作发现⑤（区室归属）的第四个队列。
`PREREG_hcc_generalisation.md` §7 中"若含恶性细胞则可做"的前提不成立。

**后果二（须明说）**：该数据集**满足本节事先固定的全部六项标准**，
唯一不符的是"癌种与发现队列相同"——而这一条**不在**六项标准之内。
因此 §4"符合全部六项标准的公开数据集有且仅有一个"这一表述，
准确的说法是：**在黑色素瘤内有且仅有一个**；若允许跨癌种，则有两个。
本文对此的处理是保守的：跨癌种数据集**不计入复现**，
因为在另一癌种、另一治疗方案（含 TKI）中的阴性结果无法与"原结论不成立"区分。

**已执行（Step 87）**：据 GSE235863 做了 Part II 的分层检验，
另写预注册（S21）后再跑，并按**泛化**而非复现报告。
结果：治疗后肿瘤臂 Δ=+0.874、单侧置换 p=0.043（发现队列 +0.87），
治疗前肿瘤臂 Δ≈0 且该臂 4 对 2、**事先算出最小可达 p=0.0667，不可能显著**——
即本节排除 GSE294273 所用的那把尺子，这里对自己也用了一次。
→ 与 Pozniak 合看：**两个独立队列各支持一支、各否定另一支**，详见正文 §3.3 与 S21 §10.2。

## 5. 这一节本身与论文主线的关系

本文已两次记录"结论难以被独立检验"的具体形式：TPI1 在 8/8 刺激态 eQTL 数据集中
未被定量；探针法 Visium 的探针集不含目标基因与阳性对照。本次检索是第三个实例——
**在公开数据中，同时满足"ICB 治疗 + 肿瘤浸润 + CD4 可解析 + 有应答标注 + 方案不混杂"
的黑色素瘤单细胞队列，只有一个。**

这不是抱怨，而是对该类结论**可核验性**的一个计量：一项研究若报告此类发现，
其读者事实上只有一次独立检验的机会。
