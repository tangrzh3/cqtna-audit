# C1 侦察：有没有第二个真正的 dynamic / context-specific eQTL 资源

**写于 2026-08-14。**
**状态**：⚠ **基于 2026-08-10 的本地普查结果，本次未能刷新。**
本会话中 eQTL Catalogue API（`https://www.ebi.ac.uk/eqtl/api/v2/*`）
**全端点返回 HTTP 500**（能连上服务器，非沙箱阻断），
WebSearch/WebFetch 亦因模型配置错误不可用。
**因此下表须在动工前重跑 `step39_find_stim_datasets.py` 复核。**

来源：`39b_stim_dataset_coverage.tsv`（Step 39），
循环性判定来源：`step39b_circularity_check.py`（Step 39b）。

---

## 一、eQTL Catalogue 里的"刺激态 T 细胞"数据集（8 个）

| 研究 | 数据集 | 细胞 / 条件 | n | 能否用作第二个 dynamic 暴露资源 |
|---|---|---|---|---|
| **Nathan_2022** | QTD000666 | CD4⁺ activated | **248** | ★ **最佳候选** |
| Nathan_2022 | QTD000684 | CD8⁺ activated | 232 | 非 CD4，可作跨细胞类型对照 |
| **Schmiedel_2018 (DICE)** | QTD000484 | CD4 T anti-CD3/CD28 **4h** | 89 | ○ 可用但只有单一时点 |
| Schmiedel_2018 (DICE) | QTD000494 | CD8 T anti-CD3/CD28 4h | 88 | 非 CD4 |
| Cytoimmgen | QTD000693 | CD4 Naive STIM **16h** | 99 | ❌ **循环** |
| Cytoimmgen | QTD000690 | CD4 Naive STIM 40h | 94 | ❌ **循环** |
| Cytoimmgen | QTD000691 | CD4 Memory STIM 40h | 94 | ❌ **循环** |
| Cytoimmgen | QTD000692 | CD4 Memory STIM 5D | 93 | ❌ **循环** |

### Cytoimmgen 为什么判为循环

设计与每时点样本量与本文暴露（Soskic）几乎逐格吻合：
Soskic 为 Naive 0h/16h/40h/5d = 99/99/89/85，Memory = 100/95/89/90；
Cytoimmgen 的 16h 为 99（Naive）与 95（Memory），40h 为 94/94，5d 为 93。
→ 高度疑为**同一实验的标准化再处理**。
**不得**用它做"独立复制"，此判定 Step 39b 已登记，Methods §17 已写入。

---

## 二、判断

**Nathan_2022（QTD000666，CD4⁺ activated，n = 248）是唯一真正合格的候选。**

有利：
- **n = 248**，是本文暴露（85–100/时点）的 2.5 倍以上；
- **独立研究、独立平台**（单细胞、以细胞状态而非固定时点定义 context）；
- 已在 eQTL Catalogue 标准化处理，全 summary stats 可经 API 取得；
- 本文**已经引用它**（参考文献 [8]），且**已经用过它**——
  genotype × pseudotime 交互检验就来自这份数据。

⚠ 三条必须事先写下的问题：
1. **它已被用于 TPI1 的交互检验。** 再拿它做暴露资源不构成循环
   （交互检验用的是 TPI1 单基因的动态遗传效应，归属检验用的是全基因组提名落点），
   但**两次使用必须在正文中分别说明**，否则读者会以为是同一次分析。
2. **"activated" 的 condition_label 在表中显示为 `naive`**，字段语义需核实——
   这正是 §一 表格必须重跑复核的原因之一。
3. **TPI1 在其中未被量化**（与 DICE、Cytoimmgen 相同）。
   对**归属基准**无影响（不依赖 TPI1），但**对 TPI1 那条线索无帮助**。

**DICE（Schmiedel_2018）是次优候选**：n=89 与本文相当，但只有 4h 单一时点，
无法支撑"dynamic"这一层，只能支撑"context-specific"。
可作为**多资源简化归属基准**的第三格。

---

## 三、建议的做法

审稿意见的 C1 有两个变体，本地证据支持**后者**更划算：

| 变体 | 可行性 | 判断 |
|---|---|---|
| 单独把 Nathan_2022 当第二个 dynamic 暴露资源，重跑全套归属 | 可行 | 一格，n 更大，但仍只是"第二个资源" |
| **跨多个细胞特异 eQTL 资源做简化版 known-locus attribution benchmark** | 可行 | ★ 用 Nathan + DICE + eQTLGen（+ 本文暴露）四格，直接回答"这是不是该框架的一般性质"，而不是"再加一个数据点" |

⚠ 无论选哪个，**动工前须**：
1. 重跑 Step 39 复核上表（API 恢复后）；
2. 对每个新资源重跑 Step 39b 式的循环性检查
   （比对设计、每时点 n、以及是否与 Soskic 或 eQTLGen 同源）；
3. 写第十四份预注册——**本文件不是预注册**，只是侦察记录。
