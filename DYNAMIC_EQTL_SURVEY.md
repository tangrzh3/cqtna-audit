# C1 侦察：有没有第二个真正的 dynamic / context-specific eQTL 资源

**首版 2026-08-14（基于 Step 39，2026-08-10 的 API 读数）**
**本版 2026-08-14 修订（`step114`）：改用静态元数据表重跑，结论有实质变化。**

## 数据来源的变更（重要）

eQTL Catalogue 的 **API 全端点持续返回 HTTP 500**（`/api/v2/*` 与不带版本的
`/api/*` 皆然；api-docs 页与 FTP 树可达，故是 API 后端故障，不是整站）。

改用**静态元数据表**，即 API `datasets` 端点背后的同一份清单：

    github.com/eQTL-Catalogue/eQTL-Catalogue-resources/data_tables/dataset_metadata_r8.tsv

**1,205 条记录，466 个基因表达（`ge`）数据集。** 此路可用，C1 不再被 API 阻塞。

---

## ⚠ 一处与旧记录的冲突，未能判定

| 数据集 | Step 39（API，2026-08-10）| `dataset_metadata_r7` | `dataset_metadata_r8` |
|---|---|---|---|
| Nathan_2022 CD4+_activated (QTD000666) | **248** | **147** | **147** |
| Nathan_2022 CD8+_activated (QTD000684) | **232** | **42** | **42** |

r7 与 r8 互相一致，与 API 旧读数不一致。**无法向 API 复核（其仍 500）。**
按两份独立静态表一致，**暂以 147 为准**，但此冲突须记录而非抹去。

→ **本文件首版所写"n = 248，是本文暴露的 2.5 倍以上"是错的，据此更正为 n = 147，约 1.5 倍。**
Nathan 仍是最强候选，但优势没有原先记录的那么大。

---

## 一、刺激态 / 非静息 T 细胞数据集（本次 64 + 11 条）

### ★ 合格候选

| 研究 | 数据集 | 细胞 / 条件 | n | 判断 |
|---|---|---|---|---|
| **Nathan_2022** | QTD000666 | CD4⁺ activated | **147** | ★ 最强候选；独立研究、单细胞、以细胞状态定义 context |
| **Randolph_2021** | QTD000588 | CD4 T，**流感感染 6 h** | **89** | ★ **本次新发现**，首版遗漏。真实刺激语境，独立研究 |
| Randolph_2021 | QTD000590 | CD8 T，流感 6 h | 88 | 非 CD4，可作跨细胞类型对照 |
| Randolph_2021_reannotated | QTD000852–862 | 六个 T 亚群，流感 6 h | 59–89 | 同一实验的再注释，**不得与上条并列计数** |
| Schmiedel_2018 (DICE) | QTD000484 | CD4 T anti-CD3/CD28 **4 h** | 89 | ○ 单一时点，只撑 "context-specific"，撑不起 "dynamic" |

⚠ **Randolph_2021 是首版遗漏的**。它落在 pass 2 而非 pass 1，因为其
`condition_label` 是 `Influenza_6h`，不含 "stim"/"activ" 等关键词。
**教训：刺激语境的命名不一定含刺激类词汇，关键词表须包含具体刺激物名称。**

### ❌ 排除：Cytoimmgen（循环）

**本次从 4 个数据集扩到约 50 个**，覆盖 TEM / TCM / TN1 / TN2 / nTreg /
TN_IFN / TM_cycling / T_NFKB / ER-stress 等多个亚群，时点 16H / 40H / 5D。

循环性判定（Step 39b 规则）**依然成立**：其 CD4 臂的样本量与本文暴露（Soskic）
逐格吻合——匹配到 89、90、95、99 四个值，且 naive/memory × 16h/40h/5d 的设计一致。
**判为同一实验的标准化再处理，不得用作独立复制。**
⚠ `Cytoimmgen_reannotated`（UNS_16H 等）同源，**一并排除**。

### △ 语境特异但非动态：IBDverse

11 个数据集，`condition_label` 全为 `naive`（未刺激），按**细胞亚群**分层，
n 达 79–339（CD8⁺ TRM 339、CD4⁺ memory 318、CD4⁺ PASK 286）。

→ **不是 dynamic 资源**，但样本量远大于本文暴露，
适合作"跨资源简化归属基准"（C1 的第二个变体）里的一格，
用以区分"细胞亚群特异"与"刺激动态"两种 context。

---

## 二、修订后的判断

C1 的两个变体，本次证据更偏向**后者**：

| 变体 | 可行性 | 判断 |
|---|---|---|
| 单独把 Nathan_2022 当第二个 dynamic 暴露资源 | 可行 | n=147 而非 248，优势缩水；仍是一格 |
| **跨资源简化 known-locus attribution benchmark** | 可行且更强 | ★ 现在有 **Nathan（活化态）+ Randolph（流感 6 h）+ DICE（anti-CD3/CD28 4 h）+ IBDverse（亚群，大 n）+ eQTLGen（全血）+ 本文暴露** 六格，能同时区分"动态刺激"与"细胞亚群"两类 context |

**推荐后者。** 它直接回答"位点归属是不是该框架的一般性质"，
而不是再添一个数据点；且不依赖任何单一资源的样本量。

---

## 三、动工前仍须做的

1. **Randolph_2021 与 Nathan_2022 的循环性逐一手查**——
   `step114` 的旗标只是提示（两个以上样本量巧合），不是判定。
   两者与 Soskic 的匹配分别为 [89] 与 []，均未触发旗标，但仍须核对研究设计。
2. **确认各资源的 TPI1 量化状况**（首版已知 TPI1 在 DICE / Cytoimmgen / Nathan 中均未量化）。
   ⚠ 对**归属基准无影响**（不依赖 TPI1），仅影响 TPI1 那条线索。
3. **Nathan 的双重使用须在正文分别交代**：
   它已被用于 TPI1 的 genotype × pseudotime 检验（参考文献 [8]）。
   再作暴露资源不构成循环（一个是单基因动态遗传效应，一个是全基因组提名落点），
   但**两次使用必须分别说明**。
4. **写第十四份预注册**——本文件不是预注册，只是侦察记录。
