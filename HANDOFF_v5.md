# 项目交接文档 v5

**项目**：CD4⁺T 细胞活化动态 eQTL → 靶点提名的审计
**最后更新**：2026-08-14（Step 108 后，第四轮审稿意见已到）
**取代**：`HANDOFF_v4.md`（其第七节技术坑、第八节方法论纪律仍有效，本文件继承并扩充）

> **新会话阅读顺序**
> 1. 本文件（尤其 §三：审稿要求，那是下一轮的全部工作）
> 2. `manuscript/MANUSCRIPT_GB.md`（**投稿母稿**）
> 3. 全文源 `manuscript/MANUSCRIPT_v2_dual_thread.md`（**唯一正文源**）
> 4. **十二份预注册**（S9/S18/S20/S21/S22/S24/S25/S28/S29/S30/S32/S33）

---

## 一、当前定位（一句话）

**单变异 eQTL-MR 产出的是"在这份数据条件下通过筛选的候选"，不是稳定靶点；
其可复现的部分由已知结局位点、阈值邻近性与基因/区室错归属决定。**

审稿人认为这个结论**目前证据已足以支撑**。而以下三条**目前不足以支撑**，
正文里任何接近它们的表述都要收窄：

1. ❌ outcome power 本身导致 colocalisation 下降
2. ❌ dynamic eQTL 普遍能揭示遗传调控的免疫状态
3. ❌ TPI1 代表一个遗传驱动的 CD4 代谢状态

**评分**：当前 8.3 → 只做文字与一致性修复可到 **8–8.5**；
要到 9.5–10 必须**增加一项能改变证据级别的外部验证**（§三-C）。

---

## 二、本窗口做了什么（Step 99–108）

| Step | 内容 | 结果 |
|---|---|---|
| 99 | **S25 R13 转移测试** | 判读 A：名单逐字不变而 3,434 条 z 无一相同 |
| 100 | **S26 推断单位敏感性** | record 级最保守；归属结论不依赖单位 |
| 101 | **S28 效应量匹配** | \|z\| 解释 68–80%；残差不可归因（两类几乎无共同支撑）|
| 102 | **S29 HCC 去循环** | 循环性 3/73，且该检验在设计上无力（如实写明）|
| 103 | **S30 密度匹配置换** | 4.09→**3.44×**（P=0.041）存活；衰减一并报 |
| 104–105 | **文献审计验证 + κ** | **两个已发表数字都错**：7.9%→**0%**，0.7%→**7.1%** |
| 106 | **leave-TPI1-out 轴** | TPI1 在不含它的轴上排 19/7,653；**PGAM1 排第 8，更高** |
| 107 | **S32 第二数据集复制** | DOGMA-seq 两臂 top 0.22%/0.39%；PGAM1 均排第 3 |
| 108 | **S33 非癌症结局（RA）** | **3.46×**（P=1.65e-5）；排除 MHC 仍 3.03× |
| — | **CQTNA 工具** | 可运行审计，自动化 8 条诊断中的 5 条 |

**三幕重组**已完成（GB Results 分 Part 1/2/3）。

---

## 三、★ 第四轮审稿要求（下一轮的全部工作）

### A. 投稿前必须做（四项）

#### A1. ~~修复 step108 脚本~~ —— **已完成，审稿人读的是旧 commit**

⚠ **这一条不要再做。** 审稿人引用的 `step108:199` 是 **`60cf4b9`（修复前）** 的行号，
那一行确实是 `want_rs = okada | set(landi.rsid...)`。
**修复已在 `3e89415` 提交**，当前文件第 225 行是 `ok_pos = placed["okada_RA"]`。

本窗口已做闭环验证：**重跑脚本，输出与已提交的 `108a` 逐位一致**
（3.459× / P=1.65e-05 / 83 of 87）。`108_console.log` 是修复后的日志。

→ 下一轮只需在 cover letter 里说明这一点即可。

#### A2. 网格口径不闭环 —— **真问题，必须修**

摘要说 "all six disease-by-resource combinations"，但实际是：

| | melanoma | HCC-high | HCC-low | RA |
|---|---|---|---|---|
| Soskic | ✓ | ✓ | ✓ | ✓ |
| eQTLGen | ✓ | ✓ | ✓ | ✓ |

= **4 个 outcome dataset × 2 资源 = 8 格，跨 3 个疾病**。
标题写 "three diseases and two exposure resources"（=6 格），
摘要写 six combinations，正文的六格不含 RA —— **四处互相矛盾**。

**二选一，全文统一**：
- (a) 主网格 = melanoma / HCC-high / RA × 2 资源 = **6 格**，HCC-low 作功效敏感性；
- (b) 明写 **"eight outcome-dataset-by-resource cells across three diseases"**。

**建议 (a)**：它让标题的"三疾病×两资源"真正成立，HCC 两层不独立本来就已登记（S20 §9）。

#### A3. "功效导致 coloc 下降" 的因果措辞

正文自己承认 melanoma meta 同时改了功效与组成（GB ~284 行），
但 Discussion 第 iii 条仍写成 "more outcome power moves MR and colocalisation
in opposite directions"。

**二选一**：
- 降级为 **"Do not assume that more outcome power will resolve the disagreement"**，
  并把观察限定为"在功效更高但组成不同的 meta 结局中"；
- 或在 **FinnGen 嵌套 releases 上跑同一套区域 coloc 轨迹**（那才隔离了功效）。

⚠ 这是本轮**第四次**同类问题（前三次：eQTLGen 轴、R13、FinnGen→meta 摘要）。
**建议下一轮开工第一件事就是全文搜一遍所有因果动词**，别再逐处发现。

#### A4. 单因果 coloc 的软肋

至少要做到三条：
1. **明确标出哪些结论不依赖 coloc**（④a 与 Fig 4 本来就不依赖，S9 §1 已拆过）；
2. "PP.H4 随功效下降" **限定为单因果 coloc 下的观察**；
3. 核心区域用**条件分析 / masking coloc / FinnGen credible-set 条件化**。

### B. 提交成熟度（本窗口已核实，全部属实）

| 项 | 状态 |
|---|---|
| GB 顶部内部 HTML 注释 | 仍在，须删 |
| GB 第 505、761 行的 **⚠** 标记 | **确认存在**（2 处），投稿版不该有 |
| **Fig. 6 错位** | **确认**：正文 464 行用 `Fig. 6a,b` 指 TPI1 时序，图目录 756 行称它是 self-attribution check。**三幕重排后图目录没跟着重算** |
| Methods | 仍是 "Full Methods accompany this manuscript." 一句 |
| References | **GB/NC 都没有参考文献列表** |
| **GB 缺 accepted-gene 限定** | ⚠ **确认缺失**（全文源有，GB 无）：清单由我们自行整理、比较前固定但**未预注册**、分母只有 10 个位点，**6/10 不是错误率**。**必须恢复**，否则最强结果之一显得过度确定 |

### C. 决定能否接近 10 分的三项（任选其一即可显著提级）

1. **第二个真正的 dynamic / context-specific eQTL 资源**
   （eQTLGen 是全血大样本，**不是**第二个 dynamic 资源——这是最大证据缺口）。
   或跨多个公开细胞特异 eQTL 资源做简化版 known-locus attribution benchmark。
2. **可信的多信号 / 条件化 coloc**，稳固"MR 与 coloc 背离"这一核心发现。
3. **独立、盲法、预定义的 gene-attribution benchmark**（不是作者整理的 10 个位点）。

### D. 应该做但不紧急

- **跨格综合统计**：现在靠"所有倍数>1"+若干单格 Fisher，各格互相依赖。
  应做预定义的 omnibus / hierarchical 分析，回答"跨疾病与资源，
  显著提名落入已知结局位点的总体概率提高多少"。
- **1 Mb 窗口敏感性**：100/250/500 kb 各跑一遍；LD-block 或 credible-set 版本；
  报告**连续距离分布**而不只是二分类。方向若稳定，归属现象会从"有趣统计现象"升级。
- **摘要重排**：GB 摘要 **435 词 vs 目标 350**，且塞进了 TPI1 6.7 倍这类次要结果，
  却没突出两个最有杀伤力的：**HEIDI 在 253/291 条 coloc-H3 记录上未拒绝**、
  **10 个答案明确的位点里 pipeline 有 4 个没命名公认基因**。
- **TPI1 压缩约三分之一**，收口为：
  *MR 提供的是一个 genetically tractable entry point；审计证明它标记一个表达状态，
  但不能升级为遗传驱动状态或治疗靶点。*
  ⚠ genotype×pseudotime 已为阴性，**不得**再把 16h 工具可用性解释成
  activation-dependent genetic regulation。

---

## 四、审稿人认可的强项（不要在压缩中砍掉）

- **MC1R 自我校准**：框架在最明确的位点都可能不命名 MC1R —— 全文最有说服力的结果之一
- 候选名单完全翻转、known-vs-novel 恢复率差异
- 预注册 + 错误披露 + 多层反证
- RA 排除 MHC 后仍稳健

---

## 五、技术坑（v4 第七节之外的新增）

25. **ADT 矩阵与 RNA 矩阵的转置约定相反**（DOGMA-seq）：行是条形码、列是蛋白。
    按 RNA 约定索引会**静默返回某条形码的谱**。已加维度断言。
    ⚠ `step42` 有同一 bug，但其输出为空且未被引用，**无已发表数字受影响**。
26. **ADT 文件含未过滤空液滴**（711k）而 RNA 已过滤（7k）。
    门控前必须**先与 RNA 条形码取交集再算中位数**。
27. **多参照系扫描必须按名分别收集坐标**。合并成一个集合再整体赋值，
    会让参照系混入它本该被对照的那份名单（S33 §8.1）。
28. **`git add -A` 会扫进 `node_modules`**（本窗口踩过，提交超时并留下 index.lock）。
    `.gitignore` 已补 `node_modules/`、`105a_artifact_work/`、`.idea/`、`outputs/`。
29. **`efoTraits/{id}/studies` 端点 404**；可用的是
    `associations/search/findByEfoTrait?efoTrait=<名称>`，
    以及 `studies/{GCST}/associations`（按研究取关联，一次请求）。

## 六、方法论纪律（v4 第八节之外的新增）

21. **判读表的 D 分支（"失败则不解释"）真的会救人。** 本窗口 S32 两次错误运行
    都停在 D，没有产出任何排名，因此不存在"看过错误结果再改"的问题。
22. **低基线率下 κ 会塌到 0 而一致率仍是 97.8%。** 必须**并排报原始一致率**，
    否则 κ=0 会被读成"两人不一致"。
23. **分歧不得由第一编码者单方面裁定**；κ 报**调和之前**的值。
24. **外部复核指出的 bug 要独立复算，不直接采信数字**——但也**不预设自己是对的**。
    本窗口两次：S33 的 bug 属实（我错），step108 的"未修"不属实（审稿人读了旧 commit）。
25. **重排章节后必须重算图号与图目录。** 三幕重组把 Fig 6 弄错位了，
    而重排当时只检查了正文引用顺序、没检查图目录与正文的对应。

---

## 七、文件地图（本窗口新增）

**预注册**：`PREREG_r13_transfer.md`(S25) · `PREREG_effect_size_matching.md`(S28) ·
`PREREG_hcc_decircularisation.md`(S29) · `PREREG_density_matched.md`(S30) ·
`PREREG_axis_replication.md`(S32) · `PREREG_noncancer_outcome.md`(S33)
**脚本**：`step99`–`step108` · `step105`（盲编码与 κ）· `cqtna/`
**结果表**：`99a`–`99c` · `100a`–`100b` · `101a`–`101c` · `102a` · `103a` ·
`104a`–`104c` · `105b`–`105c` · `106a`–`106c` · `107a` · `108a`
**数据**：`ra/`（R13 RA，805 MB，gzip 已校验）· `r13_extracts/`（提取缓存）
**工具**：`cqtna/cqtna.py` + `cqtna/demo/` + `cqtna_manual.md`（三条不可自动化诊断）

### git

本窗口 commit：`69af5ad` → `3e89415`，工作树干净。
⚠ 我用 `git add -A` 时扫进过 `105a_research_cache/`（45 篇全文，7.5 MB）
与 `105a_research_fetch.py`——**不是我写的，也没审阅就提交了**。如需撤下说一声。

---

## 八、当前字数

| | 词数 | 摘要 | 状态 |
|---|---|---|---|
| 全文源 | 19,680 | — | 唯一正文源 |
| assembled | 27,580 | — | 由 `assemble.py` 生成，**永不手改** |
| **GB** | **7,823** | **435** | 正文在 6,000–8,000 内；**摘要超目标 85 词** |
| NC | 3,947 | **150** | 均达标 ✅ |

---

## 九、下一轮建议顺序

1. **全文搜因果动词**（A3 的根因，已犯四次）
2. **统一网格口径**（A2，建议方案 a）
3. **恢复 GB 的 accepted-gene 限定**（B，最强结果缺限定）
4. **修图号 + 删 ⚠ 与内部注释 + 补 Methods/References**（B）
5. **摘要重排 + TPI1 压缩**（D）
6. 然后再决定 C 里做哪一项——**那是唯一能把分数推过 9 的**
