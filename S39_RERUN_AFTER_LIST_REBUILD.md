# S39 用修复后的名单重跑：三格判定改变

**日期**：2026-08-22
**触发规则**：`manuscript/PREREG_transport_grid.md`（S41）§5.2，**写于重建之前**（提交 `bd5cdfb`）
**输入变化**：`137a_list_rebuild_audit.tsv`（名单重建），已发表版本留存为 `130*.pre137.tsv`

---

## 1. 发生了什么

S39 的错配对照面板按设计应有六份名单。已发表版本**实际只有四份**——
`breast` 与 `lung` 被 `MIN_LEAD = 30` 的门槛排除，因为它们的 GRCh38 名单
分别只建出 **6** 个和 **0** 个位点。

**那不是这两个疾病的性质，是一次静默的工具失败。**
`step94b2_rsid_positions.py` 在任何一批 rsID 查询抛异常时执行 `continue`，
整批丢弃且不留记录。重建后（`step137`，含重试、逐个回退与失败审计）：

| 名单 | 已发表 | 重建后 | 建成率 |
|---|---:|---:|---:|
| lung | **0** | **271** | 85.8% |
| breast | **6** | **740** | 96.6% |
| prostate | **35** | **1,160** | 98.7% |
| colorectal | 642 | 828 | 99.4% |
| pancreas | 0 | 0 | 无源 rsID |

91 个未解析 rsID 逐条记于 `137b_unresolved_rsids.tsv`（90 个 Ensembl 未返回、
1 个无 1–22/X 的 GRCh38 映射）。

---

## 2. 判定对照（`130c` 新 vs `130c.pre137`）

| 格 | 已发表 | 重建后 | k（富集对照数）| margin_all |
|---|---|---|---|---|
| **melanoma × Soskic_CD4** | **A control clean** | **B not orthogonal** | 0 → **1** | 1.82 → **2.31** |
| melanoma × eQTLGen_blood | B not orthogonal | B not orthogonal | 1 → 4 | 2.82 → 2.83 |
| HCC_high × Soskic_CD4 | A control clean | A control clean | 0 → 0 | 3.79 → 4.50 |
| **HCC_high × eQTLGen_blood** | **A control clean** | **B not orthogonal** | 0 → 1 | 2.09 → 4.03 |
| RA × Soskic_CD4 | B not orthogonal | B not orthogonal | 2 → **5** | 1.28 → **1.18** |
| RA × eQTLGen_blood | B not orthogonal | B not orthogonal | 3 → **6（全部）** | 1.32 → 1.32 |
| HCC_low × Soskic_CD4 | A control clean | A control clean | 0 → 0 | 7.63 → 4.84 |
| **HCC_low × eQTLGen_blood** | **B not orthogonal** | **A control clean** | 1 → 0 | 4.24 → 5.52 |

**三格改变，其中一格是正文的头号格子。**

---

## 3. ⚠ 头号格子要怎么说，才不算说过头

`melanoma × Soskic_CD4` 从 "A control clean" 降为 "B not orthogonal"，
触发者是**新进入面板的 breast**：fold 2.15、**P = 0.0495**。

**必须同时说出的三件事**：

1. **P = 0.0495 距阈值不到千分之一。** S39 的判定规则用的是**名义 P**，无多重校正
   （`step130` 第 127 行）。若按 S40 协议 §5.1 的 **Holm** 校正（六个对照），
   0.0495 × 6 ≈ 0.30，**不会触发任何降级**。**这一格的翻转完全取决于用哪一条规则。**
2. **自身名单的优势反而变大了**：margin_all 从 1.82 升到 **2.31**。
   因为已发表版里 `prostate` 的 2.73 倍是**用 3% 完整度的名单**算出来的；
   换成完整名单后 prostate 掉到 1.24 倍，最强对手变成 breast 的 2.15 倍。
   **自身 fold 4.96 未变。**
3. **主网格的错配对照不受影响。** `123d` 的 `mismatch_fold`/`mismatch_p`
   用的是单一错配名单（melanoma 对 HCC），两份名单本来就完好。
   **受影响的是 S39 这个补充性多名单诊断，不是注册的主对照。**

⚠ **所以正确的表述是"降级"而不是"作废"**，且必须写明它对判定规则敏感。

---

## 4. RA 的结论变强了

正文现有措辞（[MANUSCRIPT_GB.md:312](manuscript/MANUSCRIPT_GB.md)）：

> "of the four unrelated lists, **two leave the RA CD4⁺ cell standing** —
> prostate at 1.04-fold (P = 0.59) and melanoma at 1.71-fold (P = 0.10)"

重建后，六个对照里**只剩 melanoma 一个**（1.71 倍，P = 0.102）没有富集：

| 对照 | fold | P | |
|---|---:|---:|---|
| melanoma | 1.71 | 0.102 | 未富集 |
| lung | 3.32 | 1.5×10⁻⁵ | 富集 |
| HCC | 3.05 | 0.0028 | 富集 |
| colorectal | 1.71 | 0.0018 | 富集 |
| breast | 1.69 | 0.0042 | 富集 |
| prostate | 1.39 | 0.0264 | 富集 |

⚠ **prostate 从"留它站着"（1.04 倍，P = 0.59）翻成"富集"（1.39 倍，P = 0.026）——
那个 P = 0.59 是 3% 完整度的名单给出的。**
`RA × eQTLGen` 更彻底：**六个对照全部富集**。

**➜ 把 RA 降为探索性的判断不但不受影响，证据反而更硬。**

---

## 5. 正文要改的地方

| 位置 | 现状 | 应改为 |
|---|---|---|
| L308 | "four per cell" | **six per cell** |
| L311–315 | "of the four unrelated lists, two leave the RA CD4⁺ cell standing — prostate … and melanoma" | 六个里只剩 melanoma；prostate 现为 1.39 倍 P = 0.026 |
| L320 | "the two RA cells sit at 1.28- and 1.32-fold" | **1.18** 与 1.32 |
| L934 | "The four available comparator lists are all cancers" | **six**（仍全为癌症，该论证的性质不变）|
| L1136 | "four per cell" | **six per cell** |
| L1142 | "no new lists were obtained" | 仍然成立——**没有取任何新名单，只是把原有名单正确地建了出来**。建议加一句说明 |
| **新增** | — | melanoma × CD4 的 S39 判定为 B，触发者 breast（2.15 倍，P = 0.0495），并注明对判定规则敏感（§3）|

⚠ **"at least 30 placeable lead SNPs" 这个措辞本身是准确的**，
但读者会理解成"这些疾病没有足够的已定位位点"——
而 breast 实际有 740 个、lung 有 271 个。**"placeable"承担不了它现在承担的意思。**

---

## 6. 方向

⚠ 这次修复**同时往两个方向推**：
让一格降级（melanoma × CD4：A → B），也让一格升级（HCC_low × eQTLGen：B → A），
还让头号格子的 margin 变大（1.82 → 2.31），并让 RA 的负面结论更硬。
**不是单向的。** 这一点值得记下来，因为前几轮的自查里错误方向一直是单向有利的。
