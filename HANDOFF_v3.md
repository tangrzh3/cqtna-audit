# 项目交接文档 v3

**项目**：CD4⁺T 细胞活化动态 eQTL → 黑色素瘤的孟德尔随机化研究
**最后更新**：2026-08-13（**Step 93 后**；Step 83 之后新增 84–93，见文末第九至十二节）
**取代**：`HANDOFF_v2.md`（其第七节"方法论纪律"与第八节"技术坑"仍有效，本文件补充）

> **新会话阅读顺序**
> 1. 本文件
> 2. `manuscript/MANUSCRIPT_v2_dual_thread.md` ← **唯一正文源**
> 3. `FINDINGS_step5_pigmentation.md` 的 Step 56 起（本轮全部记录）
> 4. **五份预注册**：`PREREG_power_trajectory.md`（S9）、`PREREG_pozniak_replication.md`（S18）、
>    `PREREG_hcc_generalisation.md`（S20）、`PREREG_hcc_part2_generalisation.md`（S21）、
>    **`PREREG_exposure_resource.md`（S22）**；
>    以及 `SUPP_attempt_timeline.md`（完整分母，含新增 §5.4b 层 D）

---

## 一、当前定位（一句话）

**双主线方法学论文**：Part I 证明 dynamic-eQTL MR 的靶点提名高度依赖结局 GWAS；
Part II 说明同一批数据在放下提名之后还能确立什么。

**接缝定义在暴露侧/结局侧**：单工具变量下 `z = β_out/se_out`，故 Part I 的一切
significance 由结局供给；Part II 只用暴露侧与细胞解析的患者数据。
⚠ 但 Part II **不是**独立生成的——TPI1 与糖酵解是被 MR 点名的，
正文已明写 "outcome-nominated and independently evaluated, not independently generated"。

---

## 二、经过一轮外部审稿后，各主张的当前强度

### 最强（可放心写）

| 主张 | 证据 |
|---|---|
| **① 显著信号落在已知色素/痣位点** | meta 轮按独立位点 4.09×（P=0.028）；**FinnGen R8–R12 五个嵌套功效水平上，无任何新位点基因达 FDR<0.05**；跨癌种排序符合预期 |
| **⑤ 区室归属** | TPI1 在恶性细胞比 CD4 高 **6.7×**；**三个独立队列**（GSE115978 16/16 患者、GSE72056 11/11、Pozniak PC3 +0.83）；深度匹配后仍在；TCGA bulk 生存关联可被免疫含量解释掉 |
| **② SMR/HEIDI 不能排除 LD 混淆** | 291 条中放行 253（86.9%）；FinnGen 样本内精细定位证实 MC1R 区确有 ≥3 个独立信号；**coloc 的阴性结论对窗口与先验全稳健（0/16 通过）** |
| **★★ 功效-稳定性收口** | 已知位点 vs 新位点恢复率 85.8% vs 22.8%（半功效）；2,506 vs 8,771 例（3.5×）；**12/12 预注册预测命中** |
| **⑥ 通路可工具化率** | 28 个糖酵解基因中 **3 个**达 p<5e-8、**2 个**可测、**1 个**有效应 |

### 中等（须带限定）

| 主张 | 限定 |
|---|---|
| ④a 名单随结局翻转 | 两轮交集为零 + 模拟曲线。但**power 与 dataset heterogeneity 未分离**；FinnGen 轨迹只能证明已知位点稳定，新位点在该区间根本不出现 |
| **①的跨瘤种泛化（HCC，Step 84–86 新增）** | 4 个显著位点中 3 个是已知 HCC 位点（含 PNPLA3）、两瘤种新位点交集=0、匹配功效下同量级；**但 G1 主检验 P=0.110 未达自设判据**。只能写"方向一致的支持性证据"，**不得**写"泛化已确立" |
| ④b H3/H4 后验重分配 | PARP1 排除多信号解释；**ZFYVE19 与 TPI1 无 credible set，无法判定** |
| ⑦ 功能层的纯度陷阱 | 稳健，但只在发现队列做过 |
| TPI1 的 16h 窗口 | 两谱系一致，但**同一实验**；无独立复现（8/8 数据集不定量 TPI1）；效应量在 0h 最大、16h 只是 SE 小 5 倍 |

### 弱（已下调，勿再拔高）

| 主张 | 现状 |
|---|---|
| **患者糖酵解应答分层** | 发现队列 OT 置换 p=0.0005、BT p=0.009（Wilcoxon 0.065 不显著）。**Pozniak（同瘤种）：BT 同量级（+0.56 vs +0.74）但 p≈0.09（NR 仅 3 例）；OT 未复现且方向随 Treg 翻转。** **★ Step 87 新增 HCC（跨瘤种）：post Δ=+0.874、置换 p=0.043（同量级且显著）、方向不随 Treg 翻转；pre Δ≈0（该臂 4 对 2，事先不可能显著）。** → 正文定性改为 **suggestive，第二瘤种上 post 支同向且显著，但没有任何一支拿到两个队列的一致支持**；**仍不得写 replicated** |
| TPI1 作为候选 | FDR=0.119 不过 · PP.H4=0.51 且**依赖先验**（8/16 组合通过）· 区域内**无 credible set** · DepMap 必需基因 · 无独立复现 |

**★ Step 83 新增的、必须遵守的表述约束**：
TPI1 是项目先后指定的**第 5 个**主角（前四个 PADI4 / ZFYVE19 / ZFYVE19+SMC2 / SMC2 全被推翻），
且四次更替之间**排序判据逐次更换，无任何跑前写死的排序规则**；
确立 TPI1 的唯一判据（Step 16 应答分层）现已被本文自己下调为 suggestive and unreplicated。
→ **不得**称 TPI1 为"最佳候选"或"证据链最完整的候选"，只能称"**被 MR 点名的那一个**"。
Part II 承重明确落在**区室归属 + 通路可工具化率**，不落在 TPI1 这个基因上。

---

## 三、本轮被撤回或更正的（勿恢复）

1. **基因方向计数的二项检验** —— 糖酵解基因高度相关（有效独立基因 ≈11.7/16），
   二项 P 值被高估 **20–300 倍**。置换检验给出 0.088/0.118。
   **正文只把方向计数作描述，不带 P 值。**
2. **"instruments exist only in activated states"** —— 完整工具变量矩阵推翻：
   **ENO1 的工具变量在 0h 静息态**。三个工具变量分别在 0h/16h/40h。
3. **"Nomination is contingent on outcome power"** —— 改为 "highly sensitive to the
   outcome GWAS used"，并加 "power is one demonstrated contributor"。
4. **Rashkin 的 Cochran Q = 5.33%** —— 由 "calibration evidence" 改为
   "shows no gross miscalibration but is not positive evidence"。
5. **Methods §18 "external calibration"** —— 改为 **empirical**，
   因 FinnGen 是被降采样的 meta 的组成部分。
6. **Step 77 的 Pozniak 结果** —— 用了备选 CD4 定义（流程错误），
   被 Step 79 取代（作者注释）。备选定义把 BT 效应稀释 3 倍。
7. **"效应量衰减 5 倍 / winner's curse"** —— 那是 Step 77 定义错误造成的，**不成立**。

---

## 四、还欠什么

### B 组：需要决策（不是技术问题）

| # | 事项 | 说明 |
|---|---|---|
| ~~**B1**~~ | ~~跨疾病泛化 vs 限定为 case study~~ | ✅ **Step 84–86 完成**：选了泛化。见"九、Step 84–86 增补" |
| B2 | R3 的结构压缩 | ✅ **第一阶段完成（Step 89）**：Part I 10 节→7 节、Discussion 8 节→5 节、新增 **Box 1** 收纳全部过程/错误/撤回；交叉引用已核。⚠ **第二阶段未做**：正文仍约 **11,200 词**，Abstract 已由 633→423 词。距 7,000 词还差约 4,000，且砍的是限定条件，**须先定目标期刊**（见第十一节）|
| B3 | Fig 9 是否放入 multiome 状态证据 | 我当时只改了图注声明"本图不含"，**未满足 R3 的实际要求** |
| ~~B4~~ | ~~R1-3 的假设选择分母~~ | ✅ **Step 83 完成**。见 `SUPP_attempt_timeline.md` §5 + 正文新增 §4.2b。三层分母：基因 10 个曾被点名 / 通路 1（无筛选）/ 机制假设 15 个检验、12 个有判决 |

### C 组：数据限制，只能写 Limitations

- 同一资源内观察新位点候选的进入/退出（FinnGen 全区间无新位点候选）
- ZFYVE19/TPI1 的多信号判定（两区在两个功效下均无 credible set）
- 患者结果的充分功效复现（唯一合格队列 NR 仅 3 例）
- **时间线的外部可核验性**（项目无 git；FINDINGS 是逐日追加但不构成版本记录）

---

## 五、文件地图（本轮新增/变更）

### 正文（`manuscript/`）
| 文件 | 状态 |
|---|---|
| **`MANUSCRIPT_v2_dual_thread.md`** | **唯一正文源** |
| `METHODS_draft.md` | 完整 Methods（20 节），正文只放摘要。**改这里再重新组装** |
| `REFERENCES.md` | 34 条已联网核实（CrossRef+PubMed），仅 4 条数据集著录待补 |
| `PREREG_power_trajectory.md` / `PREREG_pozniak_replication.md` / **`PREREG_hcc_generalisation.md`** / **`PREREG_hcc_part2_generalisation.md`** | **四份预注册**，各含结果登记区与修改记录（后两份为 Step 84–87 新增）|
| `SUPP_attempt_timeline.md` | **完整分母：11 次判决 = 8 推翻 + 3 支持**，五类排除，4 条事后决定 |
| `SUPP_comparator_studies.md` | 三篇对照研究逐项证据表 |
| `SUPP_replication_search.md` | 复现队列检索：5 组检索式、12 个数据集逐个排除 |
| `MANUSCRIPT_full_v1.md`、`RESULTS_draft.md`、`DISCUSSION_draft.md`、`INTRODUCTION_draft.md` | **已标注被取代，勿引用其数字** |

### 关键结果表（本轮）
| 文件 | 内容 |
|---|---|
| `58a` / `59a`–`59c` | FinnGen release 轨迹的登记预测与观测 |
| `63a`–`63c` | 锁定 16 基因 signature（三集嵌套）+ leave-TPI1-out |
| `65a`/`65b` | motif disruption（AP-1 假设否定，经验 p=1.0） |
| `68a`–`68c` / `69a`–`69c` | 区室归属：细胞层面 + 深度匹配 + 独立队列 |
| `70a`/`70b` | TCGA-SKCM 区室演示 |
| **`74a`** | **患者层面：复合 score + 保留相关结构的置换 + 重复患者三处理** |
| **`75a`** | **28 基因 × 8 profile 工具变量矩阵（S13）** |
| `77a`/`77b` | Pozniak（备选定义，已被 79 取代，保留为敏感性） |
| **`79a`** | **Pozniak 定论（作者注释）** |
| `80a` | Pozniak 组织混杂检查 |
| **`81a`** | **coloc 窗口 × 先验敏感性，11 位点 × 16 组合（S15）** |

### 图
`Fig1`–`Fig4`（原有）· **`Fig5_compartment`（重做，五面板）** · `Fig6` · `Fig7` ·
`Fig8` · **`Fig9_part2`（新）** · **`FigS1_flow`（新）**。
`FigE_spatial_compartment` 已退役。

### 新数据（本轮下载）
| 路径 | 内容 |
|---|---|
| `D:/Downloads/Pozniak/` | `Entire_TME.rds` 2.6 GB、`Malignant_cells.rds`、`GC_all_immune_sharing.rds`（**含作者细胞注释 + Response/Timepoint，用这个**）|
| `D:/Downloads/GSE316760/`、`GSE300445/` | 探针法 Visium，**目标基因不在探针集内，不可用** |
| `D:/R_ex/MR/release_extracts/` | FinnGen R8–R11 的工具变量提取 |
| scratchpad | TCGA-SKCM（Xena）表达 + 生存 |

---

## 六、技术坑（v2 第八节之外的新增）

14. **`coloc.abf` 遇到同一 pos 的多行会静默失败返回空** —— 不报错。
    必须先 `unique(d, by = "pos")`。**本项目第二次被重复位置坑**
    （第一次是 Step 59 的等位匹配）。
15. **FinnGen 跨代命名不一致** —— 清单 R8–R10 为 `summary_stats/R{N}_manifest.tsv`，
    R11–R12 为根目录 `finngen_R{N}_pheno_n.tsv`；精细定位目录 R8/R9 为 `finemapping/`，
    R10 起为 `finemap/full/{finemap,susie}/`。
16. **CrossRef 标题检索会把 Faculty Opinions 短文、预印本、会议摘要排在原文之前**
    （5 篇中招），且偶尔返回完全不相干文献。须过滤 `type=="posted-content"`
    与容器名含 Faculty Opinions/Abstract，并按标题相似度择优。
    **PubMed 查 FinnGen 时 Author Correction 排在原文之前。**
17. **探针法 Visium（Visium Human Transcriptome Probe Set）不含** TPI1/GAPDH/PKM/
    LDHA/ALDOA，某些版本连整个 HLA class I 都没有。**判别方法：特征数 ~18,000
    = 探针法；~32,000–36,000 = poly-A 全转录组。**
18. **Pozniak `Entire_TME.rds` 的 RNA `data` 层不是归一化值**（max=1225）；
    且**细胞类型注释不在这个对象里**，在 `GC_all_immune_sharing.rds`。
19. Seurat 对象取表达值前**务必确认 assay 与 layer**——同一对象的 SCT/RNA
    两个 assay 的 `data` 层性质完全不同。

---

## 七、方法论纪律（v2 第七节的补充，本轮验证有效）

7. **预注册的"优先/备选"顺序必须真的执行** —— Step 77 因为只看了一个对象就跳到
   备选定义，导致效应量被稀释 3 倍，结论一度写错。**跳到备选之前，先穷尽首选。**
8. **审稿人要求的检验要真跑，不要只改措辞** —— 本轮两次"跑出来推翻自己"：
   二项检验（Step 74）、工具变量矩阵（Step 75）。两次都是先改措辞会漏掉的。
9. **不得用平台差异、群体定义差异或样本量解释掉阴性** —— 样本量作为限定同时报告，
   但不能作为开脱。
10. **区分"n 次独立确认"与"n 个嵌套配置"** —— Pozniak 的 24/24 方向一致
    实际只有约 2 个半独立比较，正文已写明。
11. **跨文件的机械改动（编号、引用、清单）改完必须全局核一遍** —— Step 82 记为
    "清单补入 S13/S14/S15"，实际只改了正文引用、文末清单没同步，造成 S10/S13/S14
    各指两样东西，直到 Step 83 才发现。核法：把正文全部 `Supplementary S\d+` 抓出来
    去重，与文末清单逐条对照，两边都不能有孤儿。
12. **供给纪律的文件，自己要先被遵守** —— `SUPP_attempt_timeline.md` §3 明写
    "八次不得单独出现"，而正文从 Step 61 起一直单独出现，无人核。
    **凡是补充材料里写下的"正文必须如何"，都要回正文验一遍。**
13. **撤回一条主张时必须全文搜其同义表述**（Step 87–88，**共抓到 5 处**）——
    ① §4.5 仍在正面陈述已撤回的 "instruments exist only in activated states"；
    ② §2.7 仍在报已撤回的 binomial P = 6.1×10⁻⁵ / 2.2×10⁻³；
    ③ **Fig 6 面板 b 把 `binomial P = 6.1e-5` 印在图上**；
    ④ **Fig 3 图注与 suptitle 写着 "power determines colocalisation resolution"**，
       而 §2.4 明写"我们不这样表述"；
    ⑤ **Fig 4 图注与出图脚本把 FinnGen 轮称作 "independent measurement"**，
       而它是 meta 结局的组成部分（撤回 5 的同一问题）。
    → **检索范围必须包含 `figures/*.py` 的全部字符串字面量**（`suptitle`/`set_title`/
    `ax.text`/`annotate`/`set_xlabel`）**与正文图注**。**5 处中 3 处在图里，图是盲区。**
    第 11 条覆盖编号层，本条覆盖论断层。**Step 88 已把 8 条全部跑完并复核干净。**

---

## 八、下一步建议

~~**优先级 1**：决定 B1~~ ✅ Step 84–86 完成（做了泛化）。
~~**优先级 2**：B4~~ ✅ Step 83 完成。

~~**优先级 1**：是否用 `GSE235863` 做 Part II 患者分层检验~~ ✅ **Step 87 完成**（用户决定：做，
先写预注册）。见"十、Step 87"。

**优先级 1**：B2/B3（结构压缩、Fig 9）——纯写作与制图，是目前仅剩的结构性工作。

~~**优先级 1.5**：7 条撤回的全文检索~~ ✅ **Step 88 完成**（8 条全部跑完，抓到 5 处残留并修完，
已复核干净）。⚠ 两条后果要记住：**图脚本是盲区**（5 处中 3 处在图里）；
**Fig 3/Fig 4/Fig 6 三张图需重出**（脚本已改，图未重新生成）。

**投稿前必做**：4 条数据集著录（GSE282266、GSE199994、两个 Visium）
**+ 新增 3 条**（GCST90809296 / PMID 40823170、FinnGen R12 HCC、GSE235863）·
期刊格式统一 · Methods 的 ⟨repository DOI⟩。

---

## 九、Step 84–86 增补（2026-08-12，B1 已执行）

### 1. 做了什么

**先写预注册**（`manuscript/PREREG_hcc_generalisation.md`，在读入任何 HCC 结局 SNP 之前，
含五格判读表，事先接受"泛化失败"与"结果相反"两种反证），**再跑**。
暴露侧一个字节没改，只换结局：HCC-high = GCST90809296（3,748/1,861,536，GRCh38）、
HCC-low = FinnGen R12（947/378,749）。

### 2. 结果（全部数字见 FINDINGS Step 84–85，或预注册 §9）

| | 结果 |
|---|---|
| 显著位点 | 两个功效层合计 **4 个独立位点，3 个是已知 HCC 位点**：SAMM50（距 PNPLA3 lead 26 kb）、ZNF506（TM6SF2 区，两层都在）、SUPV3L1（新位点，只在高功效层，低功效层 FDR=0.98）|
| G1 位点归属 | HCC-low **17.6×，P=0.0031**；**HCC-high（主检验）8.85×，P=0.110 —— 未过** |
| G2 交集 | **0**（melanoma 4 个新位点基因 vs HCC 1 个），命中事先预测 |
| 匹配功效对照 | melanoma 降到 HCC-high 的 Neff（30.3%）后中位 2 个位点［1–6］= 1 已知/1 新，**与 HCC 观测完全同量级** |
| 判读 | **不落在任何单一格，最接近 A 但不满足 A 的字面条件**。正文只能写"方向一致的支持性证据"|

**两处必须一并报告的偏离**（已登记于预注册 §10）：
① G1 的背景口径与预注册字面不同（用了与 melanoma 可比的未匹配背景）——
**已补跑预注册字面版**（`85e`），两种口径结论相同，且匹配版跑黑色素瘤得 4.36×/p=0.023，
与正文 4.09×/P=0.028 互证；
② 预注册 §7 说两个 HCC 单细胞数据"无应答标注"，**这是错的**（见下）。

### 3. 已写入正文的位置

`MANUSCRIPT_v2_dual_thread.md`：新增 **§2.6b**（放在 FinnGen 轨迹之后）·
Abstract 加一段 · §4.2b 加"我们自己的方法学主张"一段 · §4.6 Limitations 加三条 ·
S20 加入清单（S9–S20 已核，无孤儿）；
`METHODS_draft.md` 新增 **§9b**；
`SUPP_attempt_timeline.md` 新增 **§5.4b（层 D：本文自身方法学主张的分母 = 3）**；
`SUPP_replication_search.md` 加"注①"。

### 4. ⚠ Step 86 发现的、需要决策的事

`GSE235863` 已下载（`D:/R_ex/MR/hcc/GSE235863/`，h5ad 已解压 724 MB，191,435 细胞）。

**否定的一半**：CD45⁺ 分选，**无恶性细胞** → **不能**做发现⑤的第四个队列，该项确定不做。

**肯定的一半**：**它有应答标注**（在 GEO 的 GSM 标题里，不在 h5ad 的 `obs` 里），
可拆 **4 应答 vs 5 无应答**，配对 pre/post，肿瘤 86,252 + 血 105,183 细胞，
CD4T 59,528 细胞含 9 个亚簇，TPI1/PGAM1 等目标基因全在。
→ **它满足 `SUPP_replication_search.md` 事先固定的全部六项标准**，唯一不符的是"癌种不同"，
而那一条不在六项之内。而正文目前把患者糖酵解应答分层定为
*suggestive and unreplicated*（唯一合格队列 Pozniak 只有 3 例 NR）。

**决策**：做还是不做。做则**必须先另写一份预注册**（本次的预注册不授权，
且已在 §11 写明"另写"），写完之前不得读该对象的任何表达值。
截至目前只看过：`obs` 列名与取值域、各簇细胞数、tissue/patient/sample 计数、
`var` 中若干基因是否存在、GEO 的 GSM 标题。**未读任何表达数值，未做任何连接。**

⚠ 反面理由同样要写进新预注册：换癌种 + 换方案（含 TKI 仑伐替尼），
阴性结果无法与"原结论不成立"区分——这正是当初把它排除在复现之外的理由。

### 5. 本轮新增文件

| 文件 | 内容 |
|---|---|
| `step84_hcc_generalisation.py` / `84a`、`84b` | 已知位点参照系（73 条）+ 两个结局的 MR |
| `step85_hcc_tests.py` / `85a`–`85d` | G1/G2/G3 + 匹配功效对照 |
| `step85e_matched_background.py` / `85e` | 预注册字面版背景（含 melanoma 交叉验证）|
| `manuscript/PREREG_hcc_generalisation.md` | 预注册，§9 结果登记区已填满，§10 两处偏离，§11 新机会，§12 正文表述约束 |
| `hcc/`（约 1.5 GB）| 两个结局的 sumstats + `GSE235863` 的 h5ad |

### 6. 技术坑（续第六节）

20. **GEO 的 h5ad 里没有的临床标注，往往在 GSM 层的 `!Sample_title` / `characteristics` 里**。
    `GSE235863` 的应答标注就只在那儿。取法：
    `curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE...&targ=gsm&form=text&view=brief"`。
    **只看 h5ad 的 `obs` 就断言"无标注"会出错——本项目已经犯过一次（预注册 §7）。**
21. **NCBI FTP 大文件容易断在中途**（本次 484 MB / 758 MB 处 exit 56）。
    用 `curl -L -C - --retry 10 --retry-all-errors` 续传即可，不必重下。

---

## 十、Step 87（2026-08-12）：Part II 患者分层在 HCC 中的检验 —— 主检验通过

**用户决定**：做，且先写预注册。已执行：`manuscript/PREREG_hcc_part2_generalisation.md`
（含 §0"写这份文件时已经看过什么"的据实声明、§6 跑前算出的每臂最小可达 p、§8 五格判读表），
写完之后才读表达值。脚本 `step87_hcc_part2.py` → `87a`–`87c`。

### 1. 结果

阳性对照 4/4 通过（PC1 糖酵解 Myeloid>naive CD4 +0.210 · PC2 PTPRC 91.9% ·
PC3a +1.190 · PC3b +0.714）。原始整数计数确认，库大小归一化 + log1p，未二次归一化。

| 臂 | Δ(NR−R) | 置换 p 单侧 | Wilcoxon | 最小可达 p | 对照发现队列 |
|---|---|---|---|---|---|
| **H1 肿瘤 post（主）** | **+0.874** | **0.043** | 0.114 | 0.0143 | +0.87 |
| H2 肿瘤 pre | +0.077 | 0.267 | 0.533 | **0.0667（不可能显著）** | +0.742 |
| H3 血 post | −0.051 | 0.571 | 0.905 | 0.0079 | — |
| H4 血 pre | +0.746 | 0.057 | 0.114 | 0.0286 | — |

判读落 **A**。Treg 在/不在、TPI1 在/不在四种配置结论一致（但那是同一批 8 样本，非四次独立确认）。

### 2. ⚠ 必须与"通过"一并写的三条

1. 置换 p=0.043 但**秩检验 p=0.114**（一个 NR 落进 R 群内，未完全分离）
2. **事后**留一检验：Δ 保持 +0.71~+1.06（效应量不由单样本驱动），
   但 8 次中 7 次 p 落 0.057–0.086 —— **显著性在 8 个样本能表达的边缘**
3. 跨瘤种 + 跨方案（含 TKI），**不得写 replicated**

### 3. ★ 三队列合看（正文与 S21 §10.2 的核心表）

| 时点 | 发现 GSE120575 | Pozniak（同瘤种）| HCC（跨瘤种）|
|---|---|---|---|
| Pre / BT | +0.742, p=0.009 | **+0.555, p≈0.09**（同向不显著）| **+0.077 ≈ 0** |
| Post / OT | +0.87, p=0.0005 | **方向相反**（随 Treg 翻转）| **+0.874, p=0.043** |

> **两个独立队列各支持一支、各否定另一支，互为镜像；没有任何一支拿到两次一致支持。**

⚠ 这也**推翻了**此前写在正文里的一句事后解释："最强的那一支失败，符合 winner's curse"——
HCC 里恰恰是那一支被同量级复现。正文 §3.3 已改为报告这个矛盾，不再给出容纳它的说法。

### 4. 已写入的位置

正文 §3.3（大改）· Abstract Part II 段 · §4.2b（层 D 改为四次）· §4.6 Limitations ·
S21 入清单（S9–S21 已核，无孤儿）· `METHODS_draft.md` 新增 §11b ·
`SUPP_attempt_timeline.md` §5.4b 第 4 行 + §5.4 的 #15 + §5.5 层 D ·
`SUPP_replication_search.md` 注① 补"已执行"。

### 5. 口径（勿再改动）

- post 支：**"在第二个瘤种上同向、同量级、显著"**；pre 支：**"该臂事先不可能显著"**
- **不得**写 patient-level 结论"已被复现"；**不得**用 H4（血 pre，p=0.057）单独叙述
- 四种配置一致 ≠ 四次独立确认

---

## 十一、Step 89 之后的长度问题（新增，投稿前的第一道门槛）

**首次量化**：正文叙述约 **11,200 词**（结构压缩后），Abstract **423 词**。

| 期刊 | Abstract | 正文 | 现状 |
|---|---|---|---|
| Genome Biology / Genome Medicine | ~350 | 无硬限，惯例 6–8k | Abstract 略超；正文超 40–60% |
| Nature Communications | **150** | **5,000** | 两项都远超 |

→ **目标期刊必须先定**，否则第二阶段压缩没有靶子。按 Step 88 后的评估，
推荐 Genome Biology / Genome Medicine 路线；若走 Nat Commun，正文要砍掉一半以上，
本文最值钱的东西（每条结论后面的限定）会被削掉，**不建议**。

**第二阶段压缩的候选靶点**（按字数排）：
§2.7（1,154）· §3.3（1,100）· §2.4（1,050）· §3.1（812）· §4.4 Limitations（777）· Box 1（725）· Intro（604）。

⚠ 压缩原则（定死，防止越压越弱）：
**优先降级到补充材料，而不是删除**。每删一条限定，都要确认它没有在别处被引用为"我们已声明过"。

---

## 十二、Step 91–93（2026-08-12/13）：第二轮审稿的"距离 8 分"与"距离 9 分"

### 1. 距离 8 分：七条全做完（Step 90 + 93）

见 FINDINGS Step 90 与 93。要点：Part II 三层口径拆分 · 11 项重分类为四分法
（**"8/11" 已从摘要与结论全部撤下**）· 患者分析定死唯一 primary test（置换检验）·
R3 的 10 处版本矛盾清零 · TPI1/coloc/whichever disease 三处收紧 ·
中文笔记移出正文（`MANUSCRIPT_worknotes.md`）· Methods 由 `assemble.py` 合成 ·
参考文献 ⟨⟩ 归零 · 代码仓库 402 文件已提交 · Fig 9 补 d/e 两面板。

### 2. 距离 9 分：四条里做掉两条

| 审稿人给的四选一 | 状态 |
|---|---|
| **换一套 eQTL 资源重复位点归因** | ✅ **Step 92 完成，判读落 A**（见下）|
| **系统性文献审计** | ✅ **Step 91 完成，152 篇全文**（见下）|
| 同疾病同治疗的独立患者复制 | ❌ 已穷尽：S14 检索只有一个合格队列，已用 |
| in-sample LD 解 TPI1/ZFYVE19 coloc | ❌ **被堵死**：FinnGen 全基因组只精细定位 19 个区，本文候选位点无一在内（正文已写明）；meta 结局无 in-sample LD |

**Step 92 结论（可进摘要）**：结局固定、暴露换成 eQTLGen 全血（n=31,684，约 300 倍），
FDR<0.05 位点由 7 增至 30，已知位点占比由 42.9% 升至 **66.7%**，
富集 **4.44× vs 4.09×** —— **两套毫无共同点的暴露资源给出几乎相同的倍数**。
⚠ 但新位点基因交集 = 1（ZFYVE19），故**"新位点提名全部不可复现"是过强说法**，正文已改。

**Step 91 结论**：152 篇全文中，**至多 12 篇（7.9%）**做过与自身结局已知位点的比较，
**1 篇（0.7%）**报告过名单对结局 GWAS 的依赖。已写入 §4.5 诊断①。

### 3. 新增纪律第 14 条

**撤回一个统计量时，产出它的结果表列名也要一并标记。**
Step 93 给 Fig 9 加面板时，又在 panel c 上抓到被撤回的二项 P 值——
而 Step 88 明明已按第 13 条把图脚本扫了一遍。原因是
`63c_patient_locked.tsv` 里**数据列就叫 `binom_p`**：
只要有人再拿这张表画图，撤回的量就会从**数据侧**复活。
→ 已在图脚本写死 "Do not restore the P values" 及理由。

### 4. 仍未做

- **正文长度**：约 11,600 词（又加了两块新结果），Abstract 约 450 词。
  第二阶段压缩仍待定期刊。
- **repository DOI**：需作者推 GitHub + Zenodo 发 release（`DEPOSIT.md` §5 有四步说明）
- DepMap 的具体 release 号须与脚本核对
- 距离 10 分的四条（个体级 genotype×time、功能扰动等）均需新实验，非本轮范围
