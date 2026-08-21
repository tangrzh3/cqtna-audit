# S41 迁移网格的结果，以及 S39 降级对文章结论的影响

**日期**：2026-08-22
**产物**：`138a_transport_main.tsv`、`138b_transport_verdict.tsv`、`138c_console.log`
**规则**：`manuscript/PREREG_transport_grid.md`（S41），冻结于 `cf59123` / `bd5cdfb`，
**均先于本结果**

---

# 甲、S41：机制迁移了

## 1. 主表

| 结局 | 显著位点 | 自身 fold | P | **P (Holm×6)** | 置换 P | 判定 |
|---|---:|---:|---:|---:|---:|---|
| **Melanoma** | 8 | 4.96 | 0.00484 | **0.00968** | 0.005 | **supportive** |
| **Lung** | 10 | 3.85 | 0.0138 | **0.0138** | 0.0107 | **supportive** |
| Colorectal | 6 | 3.15 | 0.000976 | **0.00293** | — | inconclusive（功效：6 < 8）|
| Breast | 23 | 2.25 | 0.000265 | **0.00106** | 0.0009 | **void** |
| Prostate | 28 | 2.04 | 4.66×10⁻⁶ | **2.33×10⁻⁵** | 1.0×10⁻⁴ | **void** |
| Pancreas | — | — | — | — | — | inconclusive（无名单）|

**五个有名单的疾病，自身名单富集全部在跨六行 Holm 校正后仍然显著。**
参照行核对：melanoma 返回 8 个位点、4.96 倍、P = 0.00484，与主分析逐位相同。

## 2. Lung 是这一轮真正的新结果

10 个显著位点、自身 fold 3.85、置换 P = 0.0107、
**五个错配对照经 Holm 后无一显著**（最强者 melanoma 2.98 倍，P = 0.07）。

⚠ **而 lung 的名单昨天还不存在**——它的 GRCh38 名单原本建出 **0** 个位点。
这个正面结果是名单修复的直接产物。

## 3. Breast 与 Prostate 为什么 void——**不是自身检验失败**

两者的自身富集都很强且有置换支持（2.25 倍 P = 2.7×10⁻⁴；2.04 倍 P = 4.7×10⁻⁶）。
它们 void 是因为**错配对照也富集，且自身 fold 不是最高的**：

- Prostate：lung 2.75（P = 0.0051）、HCC 2.18、melanoma 2.13 都高过自身的 2.04；
- Breast：HCC 2.65、lung 2.09 高过自身的 2.25 中的前者，colorectal 经 Holm 后显著。

**这正是这个终点在高功效下的行为**，`step133` 用单一疾病缩放已经证明过：
显著集越接近背景，fold 必然趋向 1，于是**任何名单看起来都像富集**。
六行的观测复现了这个方向：

    8 个位点 → 4.96 倍 ｜ 10 → 3.85 ｜ 23 → 2.25 ｜ 28 → 2.04

## 4. ★ 但原始 fold 跨疾病不可比——密度上限差 4 倍

⚠ **以下为事后描述性观察，不属于 S41 §7 的判定，不得用于升级任何一行。**

fold = (自身 known 占显著集比) ÷ (背景 known 占比)。
背景 known 占比越高，**fold 的理论上限越低**：

| 结局 | 背景 known | fold 上限 | 实得 | **占上限** | 自身/最强对照 |
|---|---:|---:|---:|---:|---:|
| Melanoma | 10.1% | 9.92 | 4.96 | 50% | **2.30** |
| Lung | 10.4% | 9.63 | 3.85 | 40% | 1.29 |
| Colorectal | 31.8% | 3.15 | 3.15 | **100%** | 0.98 |
| Breast | 29.0% | 3.45 | 2.25 | 65% | 0.85 |
| Prostate | 40.3% | 2.48 | 2.04 | **82%** | 0.74 |

**Prostate 的 fold 最高只可能是 2.48，它拿到 2.04。Colorectal 的 6 个显著位点全部是 known，
即打满上限。** 按占上限归一化，排序**完全颠倒**：
colorectal 100% > prostate 82% > breast 65% > melanoma 50% > lung 40%。

**➜ "melanoma 的富集最强"这句话在原始 fold 上成立，在归一化后不成立。
差别几乎全部由已知位点名单的密度造成，而不是由疾病造成。**

## 5. melanoma 在哪一项上确实是特殊的

不是 fold，是**与错配对照的间距**：

    melanoma 2.30 ｜ lung 1.29 ｜ colorectal 0.98 ｜ breast 0.85 ｜ prostate 0.74

**只有 melanoma 和 lung 的自身名单高过全部对照**，而 melanoma 的余量是 lung 的 1.8 倍。
这一项不受密度上限影响（分子分母同背景）。

## 6. S41 不能说什么

- **不是外部验证。** 六行全部来自 FinnGen（melanoma 行还含 Rashkin），
  同一批参与者、同一条管线。**任何一行都不得称 external validation。**
- ⚠ **melanoma 行的结局与其余五行不同**：它是 FinnGen R12 + Rashkin 的 meta，
  其余五行是 FinnGen 单一 endpoint。**跨行比较带着这个混杂**，
  melanoma 的功效高于同名 FinnGen endpoint。
- 不检验 CD4 特异性、因果基因归属、患者疗效。

---

# 乙、S39 降级对文章结论的影响

## 1. 结论：**核心主张不变，两处措辞必须改，加一条新的限制**

### 1.1 注册的主对照没有被触及

正文的注册负对照是**单一错配名单**：`123d` 里 melanoma × CD4 的
`mismatch_fold = 0`、`mismatch_p = 1`（用 HCC 名单）。
melanoma 与 HCC 两份名单本来就完好，**不在本次修复范围内**。

被改变的是 **S39 那个补充性多名单诊断**，而正文自己写明它
"designed after the results it comments on… **reported as exploratory and is
not used to upgrade any cell**"。**一个明示为探索性的诊断变化，不能降级一个注册结果。**

### 1.2 这次降级完全取决于用哪条规则

触发者 breast：fold 2.15、**P = 0.0495**。

- S39 用**名义 P、无校正**（`step130` 第 127 行）→ 触发，判定 A → B；
- S40 协议 §5.1 用 **Holm** → 0.0495 × 6 ≈ 0.30 → **不触发**；
- `step138`（S41）对对照族用 Holm → melanoma 行 **supportive，全部对照干净**。

**➜ 同一份数据，两条规则给出相反的标签。正文必须写明这一点，
而不是只报其中一个。**

### 1.3 头号格子的余量反而变大

`margin_all` 1.82 → **2.31**。因为已发表版里 prostate 的 2.73 倍是
**用 3% 完整度的名单**算出来的；补全后 prostate 掉到 1.24 倍。
**自身 fold 4.96 未动。**

### 1.4 RA 的负面结论变硬

六个对照里只剩 melanoma（1.71 倍，P = 0.102）不富集；
`RA × eQTLGen` 六个对照**全部富集**。prostate 从 P = 0.59 翻成 P = 0.026。
**把 RA 降为探索性的判断证据更足。**

## 2. ★ 一条必须新加的限制

S41 显示：**在高功效下这个终点会饱和**——显著集逼近背景，fold 趋向 1，
错配对照跟着一起富集（prostate 28 个位点时，四个对照有三个 Holm 后显著）。

**所以 melanoma 那一格之所以"对照干净"，部分是因为它只有 8 个显著位点。**
这不是说结果是假的——melanoma 的余量 2.30 是全网格最高，
是唯一一个自身名单明显高过全部对照的——
**但"对照干净"这件事本身对功效敏感，不能当作与功效无关的特异性证据。**

**➜ 建议在 Discussion 的限制段加一句**：

> The mismatched-list control is informative only within a bounded power range.
> As the significant set grows towards the background the fold necessarily
> approaches 1 and unrelated lists begin to enrich as well: across six outcomes
> the own-list fold falls monotonically from 4.96 at eight significant loci to
> 2.04 at twenty-eight. A clean mismatched control at low power and a failed one
> at high power are therefore not directly comparable.

## 3. 正文要改的地方（合并本轮全部）

| 位置 | 现状 | 应改为 |
|---|---|---|
| L308 · L1136 | "four per cell" | **six per cell** |
| L311–315 | "of the four unrelated lists, **two** leave the RA CD4⁺ cell standing — prostate at 1.04-fold (P = 0.59) and melanoma at 1.71-fold (P = 0.10)" | 六个里**只剩 melanoma**（1.71 倍，P = 0.102）；prostate 现为 1.39 倍 P = 0.026 |
| L320 | "the two RA cells sit at **1.28**- and 1.32-fold" | **1.18** 与 1.32 |
| L934 | "The **four** available comparator lists are all cancers" | **six**（仍全为癌症，该论证性质不变）|
| L1142 | "no new lists were obtained" | 仍成立——**没有取任何新名单，只是把原有名单正确地建了出来**；建议加一句说明 |
| L133–137 | "**Applied unchanged** to five other cancers, melanoma showed the largest enrichment…" | **必须重写**：那句依据的 `36e` 用的是**单连锁**约定、且五个癌种打的是 **melanoma 名单**（错配对照）。现在有了 S41 的真迁移结果，应当直接换成 S41 |
| 新增 | — | melanoma × CD4 的 S39 判定为 B，触发者 breast（2.15 倍，P = 0.0495），**并注明它对判定规则敏感**（§1.2）|
| 新增 | — | §2 的功效饱和限制 |

⚠ "at least 30 **placeable** lead SNPs" 字面准确，但读者会理解成
"这些疾病没有足够的已定位位点"——breast 实际有 740 个、lung 有 271 个。
**"placeable"承担不了它现在承担的意思，建议改写。**

## 4. 净影响

**加分**：多了一个真正的跨疾病迁移结果（lung，全对照干净），
五个疾病自身名单富集全部经 Holm 存活；RA 的降级更有据。

**减分**：头号格子的一个**探索性**诊断从 A 降到 B，且该降级对规则敏感；
新增一条功效饱和的限制。

**不变**：注册的主对照、主网格的六格判定、8 个显著位点、4.96 倍、P = 0.0048。
