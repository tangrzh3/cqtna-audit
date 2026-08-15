# GB 九图 ↔ 图件文件 对照与缺口审计

日期：2026-08-15（HANDOFF_v6 §三-A1 的落实）
适用：`manuscript/MANUSCRIPT_GB.md` §Figures 所列的 **Fig 1–9**

> `FIGURES_plan.md` 记的是全文源 `MANUSCRIPT_v2_dual_thread.md` 的**八图**方案，
> 与 GB 的九图**不是同一套编号，也不是同一套分组**。
> GB 重排后，若干面板从原来的 `Fig9_part2` / `FigS7` 里被拆出来独立成图。
> **两份方案都有效，不要互相覆盖**；本文件只管 GB 这一套。

---

## 一、对照表

| GB | 图注要求 | 现有资产 | 数据源 | 状态 |
|---|---|---|---|---|
| **1** | 两个结局下的位点归属，**按独立位点** | — | `06` · `13` | ✅ **已建**：`Fig1_locus_attribution`。镜像 Manhattan，每位点一点。重算得 meta 3/7 = 4.09× P=0.028、FinnGen 轮 2/2 = 9.57× P=0.011，与 `36e`/`99c` 一致。⚠ 截断改为 8（`fig_A` 用的 30 会把按位点收拢后的图压到中线） |
| **2** | 泛化：五个嵌套 release · 迁移检验 · 第二疾病 · 非癌结局 · 第二暴露资源 · 双轴交叉网格 + 错配位点对照 | — | 见 §二 | ✅ **已建**：`Fig2_generality`（`make_gb_fig2_generality.py`），三面板。**只读 `119a`/`119c`，不读 `94d`/`94e`**（§三） |
| **3** | 按位点类别的功效-稳定性，**以及同一分层对全功效 \|z\|** | — | `53a`–`55a` · `101a` · `101b` | ✅ **已建**：`Fig3_power_stability`。复用 `FigF` 的 a/b/c，补 d（两类候选的 \|z\| 区间）与 e（原始差对 \|z\| 匹配后残余）。⚠ 正文的"仅用 \|z\| 重现 68–80% 差距"**无落盘的表**，图上未画 |
| **4** | 共定位 vs SMR/HEIDI，**含 in-sample 精细定位** | — | `15` · `16` · `60a` | ✅ **已建**：`Fig4_coloc_heidi`。a = 253/291 = 86.9%；b = MC1R 样本内 3 个 credible set 对代理 LD 的 10+9=19，并在图上写明单向解读规则。样本内计数只有 MC1R 有，其余三区标注为无 |
| **5** | 通路上三个层级的工具变量可得性 | — | `75a` `34a` `34b` | ✅ **已建**：`Fig5_instrument_ladder`（`make_gb_fig5_instrument_ladder.py`）。28 → 3 → 2 → 1，与正文逐字吻合 |
| **6** | TPI1 可工具化的活化窗口：**八个 profile 上效应量对精度** | — | `120a`（step120 新建） | ✅ **已建**：`Fig6_tpi1_window`。a = β±95%CI；b = β–SE 平面加 P=5e-8 射线。⚠ 原 `Fig9_part2` 面板 a 的八组数是**硬编码**的（注称 Step 62，无脚本无表），`step120` 从 parquet 重算，八个全部逐位重现 |
| **7** | CD4⁺ 代谢轴及其染色质签名 | — | `43b` · `43c` · `44b` · `48b` · `50b` | ✅ **已建**：`Fig7_axis_chromatin`，四面板。⚠⚠ **panel b 与正文现有措辞冲突，按表画**，见 §五 |
| **8** | 区室归属 | `Fig5_compartment`（a–e，五面板） | `22c` `23e` `26c` `69b` `70a` | ✅ 内容吻合，**只需改名**。（`FigE_spatial_compartment` 是更窄的纯空间版，不是这张） |
| **9** | 三个队列的患者 | — | `74a` · `79a` · `87b` | ✅ **已建**：`Fig9_patients`。a = 三队列×两 arm 森林图；b = 跨病种治疗前 arm 的 **P 下限 0.067**（4R/2NR，设计上不可能达 0.05，S21 预注册）。同病种队列是 `79a`（Pozniak）。⚠ 见 §五第 2 条 |

**小结**：真正"没有文件"的只有 Fig 2 与 Fig 5；但 Fig 1/3/4 各缺一个面板，
Fig 6/7/9 需要从现有合成图里拆分重组，Fig 8 只需改名。
HANDOFF_v6 的判断（"不是编号问题，是图件本身没按九图口径导出"）成立。

**九图脚本一览**（全部 `figures/make_gb_figN_*.py`，各自输出 `FigN_*.pdf` + `.png`）：

| GB | 脚本 | 输出文件名 |
|---|---|---|
| 1 | `make_gb_fig1_locus_attribution.py` | `Fig1_locus_attribution` |
| 2 | `make_gb_fig2_generality.py` | `Fig2_generality` |
| 3 | `make_gb_fig3_power_stability.py` | `Fig3_power_stability` |
| 4 | `make_gb_fig4_coloc_heidi.py` | `Fig4_coloc_heidi` |
| 5 | `make_gb_fig5_instrument_ladder.py` | `Fig5_instrument_ladder` |
| 6 | `make_gb_fig6_tpi1_window.py` | `Fig6_tpi1_window` |
| 7 | `make_gb_fig7_axis_chromatin.py` | `Fig7_axis_chromatin` |
| 8 | `make_gb_fig8_compartment.py` | `Fig8_compartment` |
| 9 | `make_gb_fig9_patients.py` | `Fig9_patients` |

⚠ 旧的 `FigA`–`FigG`、`FigS1`–`FigS8`、`Fig5_compartment`、`Fig6_glycolysis_TPI1`、
`Fig7_crosscancer`、`Fig9_part2` **全部保留不动**——它们是全文源八图方案的图件，
两套并存。新脚本一律不覆盖旧文件名。

---

## 二、Fig 2 的六个组成部分与各自的权威数据源

| 组成 | 数字 | 表 |
|---|---|---|
| 五个嵌套 release | R8–R12，2,705→5,753 例；命中 7→10；已知位点恢复早且不动 | `59a_release_trajectory.tsv`（`own` 与 `common` 两口径）· `58a`（登记预测区间） |
| 迁移检验 R13 | 6,226 例，同六基因同两位点；归属 9.58×，P=0.0107 | `99c_r13_attribution.tsv` |
| 第二疾病 HCC | high 8.85×/P=0.110 · low 17.62×/P=0.0031 | `85c_G1_enrichment.tsv` |
| 非癌结局 RA | Soskic 3.459×/P=1.65e-5 · eQTLGen 3.566×/P=7.70e-14；去 MHC 后 3.026/3.319 | `108a_ra_attribution.tsv` |
| 第二暴露资源 eQTLGen | melanoma 4.44×/P=3.72e-11（vs Soskic 4.09×） | `92c` / `92d` / `92e`（匹配背景 4.14–4.48×，emp P=1e-4） |
| 交叉网格 + 错配对照 | 见下 | ⚠ **无单一权威表**，见 §三 |

### 主网格（三疾病 × 二资源 = 六格，四格 P<0.05）

| | Soskic CD4 | eQTLGen 全血 |
|---|---|---|
| melanoma | 4.09×，P = 0.028 ✔ | 4.44×，P = 3.7e-11 ✔ |
| HCC-high | 8.85×，P = 0.110 | 5.11×，P = 0.052 |
| RA | 3.46×，P = 1.7e-5 ✔ | 3.57×，P = 7.7e-14 ✔ |

功效敏感性（**不是第七、八格**）：HCC-low × Soskic 17.62×/P=0.0031 · HCC-low × eQTLGen 8.67×/P=0.017。

### 错配位点对照

| 组合 | 倍数 | P | 出处 |
|---|---|---|---|
| eQTLGen × melanoma 用 **HCC 名单** | 1.30× | 0.411 | ⚠ 仅见 `PREREG_generality_grid.md` §8 表，**无 TSV** |
| Soskic × HCC-low 用 **melanoma 名单** | 0.00× | 1.0 | ⚠ 同上，**无 TSV** |
| Soskic × melanoma R13 用 **HCC 名单** | 0.00× | 1.0 | `99c_r13_attribution.tsv` ✔ |
| Soskic × RA 用 **melanoma 名单** | 1.596× | 0.161 | `108a_ra_attribution.tsv` ✔ |

---

## 三、⚠ 画 Fig 2 之前必须先修的三处（新发现，不是审稿人提的）

> **状态：已由 `step119_grid_rebuild.py` 处理**，产出 `119a_grid_main.tsv`、
> `119b_grid_summary.tsv`、`119c_mismatch_controls.tsv`。
> 原表 `94d`/`94e` **保留不动**，但**任何新脚本一律不得再读它们**。

### 1. `94d_grid.tsv` 里混着两行**已作废、且无任何作废标记**的结果

`PREREG_generality_grid.md` §8.2「偏离 1」写明：五个 FinnGen 癌种的已知位点坐标映射有 bug
（`step94c` 把已知位点 rsID 拿到**已被过滤成工具变量位置的**结局提取里查坐标，
致背景已知位点只有 2/559 与 5/559，而黑色素瘤为 58/554），
故 **lung（0/8，fold 0.00）与 colorectal（3/22，fold 15.25）"不得引用"**。

但这两行**原样躺在 `94d_grid.tsv` 里**，与有效行并列，没有 `void` 列、没有注释。
任何脚本（包括图件脚本）读这张表都会把它们当成有效格子。

**处置**：`step119` 重建时给每行加 `status` 列（`main` / `sensitivity` / `void_bad_coords`），
并在表头注明作废理由与出处。原表保留不动。

### 2. `94e_grid_summary.tsv` 是**旧网格**的汇总，与现行口径不符

现内容 `n_testable=6, n_fold_gt1=5, pct=83.3`。这 6 格是
Soskic×{melanoma, HCC-high, HCC-low} + eQTLGen×{melanoma, **lung**, **colorectal**}
——**含两行作废数据，且不含 RA**。而现行主网格是 melanoma/HCC-high/RA × 2 资源，
预注册 §8 记的 U1 是 **6/6 = 100%**，不是 5/6 = 83.3%。

**处置**：`step119` 重算并另存 `119b_grid_summary.tsv`；`94e` 保留 + 加日期批注。

### 3. 两个错配对照数字没有生成脚本

`1.30×/P=0.411` 与 `0.00×/P=1.0` 只存在于预注册 §8 的记录里，
仓库中**没有任何 TSV 承载**，也没有能重跑出它们的脚本
（`step94c` 里编码的 NC 是"melanoma 名单套胰腺癌"，因胰腺癌名单为空而被放弃，
改成两名单互错配的那次运行**没有落盘**）。

这与上一窗口发现的 `96a` 无生成脚本属同一类问题，且**可以廉价补上**：
两份名单都在仓库里（`landi2020_known_loci_grch38.csv` 157 个位点 ·
`84a_hcc_known_loci_grch38.csv` 73 个位点），逐格重打分即可。

**处置**：`step119` 重算两个错配对照，落 `119c_mismatch_controls.tsv`，
并与预注册 §8 记录的数字逐位比对；**若对不上，以重算为准并披露差异**。

**结果：两个都逐位重现，正文数字成立。**

| 对照 | 重算 | 预注册 §8 记录 |
|---|---|---|
| eQTLGen × melanoma 用 HCC 名单 | 3/30 显著位点已知，背景 43/559 = 7.7% → **1.30×，P = 0.4109** | 1.30×，P = 0.411 ✔ |
| Soskic × HCC-low 用 melanoma 名单 | 0/2 显著位点已知，背景 59/564 = 10.5% → **0.00×，P = 1.0** | 0.00×，P = 1.0 ✔ |

重算同时确认主网格 **六格 fold>1 = 6/6**（预注册 U1 判据 ≥80%，达标）、
**四格 P<0.05**（与 2026-08-14 的更正一致）。

---

## 四、执行顺序

1. ~~`step119`：重建网格主表 + 汇总 + 错配对照~~ ✅
2. ~~**Fig 2** 合成图~~ ✅
3. ~~Fig 1 补 FinnGen 轮面板 + 改为位点级~~ ✅
4. ~~Fig 3 补 \|z\| 分层面板 · Fig 4 补 in-sample 精细定位面板~~ ✅
5. ~~Fig 6 / 7 / 9 从 `Fig9_part2`、`FigS7` 拆分重组~~ ✅（Fig 6 另需 `step120`）
6. ~~统一导出 `Fig1`–`Fig9`（PDF + PNG）~~ ✅

⚠ 文件名一律在 `save()` 传参处改，**不要手动重命名 PDF**（脚本之间互相引用）。

---

## 五、⚠ 画图过程中查出的正文问题（新发现，须在投稿前处理）

### 1. ★ "十一个生物学模块无一超过组成匹配零模型"——**这句话是错的**

出现在：`MANUSCRIPT_GB.md` 第 517 行正文 ·
`MANUSCRIPT_v2_dual_thread.md` 第 1756 行 Fig9 图注（`MANUSCRIPT_assembled.md` 同源）。

`43c_residual_axis_verdict.tsv` 说的是：11 个生物学模块里 **有 2 个超过**零模型
（阈值 = 50 个随机模块 |SMD| 的 95 分位 = **0.0720**，已按 `43b` 重算逐行核对一致）：

| 模块 | \|SMD\| | 超过零模型 |
|---|---|---|
| **OXPHOS** | 0.2125 | **是** |
| **Proliferation** | 0.1180 | **是** |
| IFN_response | 0.0628 | 否 |
| Exhaustion | 0.0619 | 否 |
| Th1 / Tfh / Cytotoxic / Th17 / Treg / Th2 | ≤ 0.046 | 否 |
| **Activation** | **0.0032** | 否 |

而且 **Proliferation 是 `step43_glyco_vs_activation.py` 里预先写死的阳性对照**，
脚本明写"若阳性对照不超过零模型，则该检验功效不足、其零模型无信息"
——**它超过是检验成立的必要条件，不是反例**。

⚠ 现有的 `make_fig9_part2.py` 面板 d 注释同样错，且自相矛盾：
它渲染出的字面是"all 11 biological modules within the matched null (**9/11**)"。

**正文真正要说的那件事仍然成立，且改对后更强**：
`Activation` 模块自己就贴在零模型上（SMD = −0.003），
即这条轴不是"活化强度"的换个说法；超过零模型的两个是**阳性对照**与 **OXPHOS**，
后者与"合成代谢状态"的读法**一致而非相左**（糖酵解 + OXPHOS + 核糖体蛋白）。

**建议改法**：把"none of the eleven … exceeds a composition-matched null"改为
"the Activation module itself sits at the null (SMD = −0.003), and of the eleven
biological modules only the Proliferation positive control and OXPHOS exceed it"。

`Fig7_axis_chromatin` 的 panel b **已按表画全部 11 个模块**，不按现有正文措辞画。

### 2. 同病种复制队列的 +0.56 与发现队列的 +0.74 **不是同一个 signature**

正文："pre-treatment matches in magnitude in the same-disease cohort
(Δ = +0.56 versus +0.74)"。

- `79a` 同病种 BT：`locked16` = **0.5548**（→ +0.55）· `locked15_noTPI1` = **0.5648**（→ +0.56）
- `74a` 发现 Pre：`locked16` = **0.7417**（→ +0.74）· `locked15_noTPI1` = 0.7099（→ +0.71）

即 **+0.56 取自去 TPI1 的 15 基因版，+0.74 取自 16 基因版**，两个数不同口径。
统一到 `locked16` 应写 **+0.55 对 +0.74**。`Fig9_patients` 三个队列一律用 `locked16`。

（跨病种那格没有这个问题：正文的 +0.874 与 +0.08 都取自 `87b` 的 **nonTreg** 行，
与 `79a` 的主定义"CD4_Tcells（作者注释，主）"口径一致。）

### 3. "仅用 \|z\| 的模型重现 68–80% 的差距"没有落盘的表

`step101_effect_size_matching.py` 不输出这个量，`101a`–`101c` 里也没有。
`Fig3_power_stability` 因此不画它，只画 `101a`/`101b` 里有的量
（\|z\| 区间、原始差、匹配后残余）。投稿前须补表或在正文注明其出处。
